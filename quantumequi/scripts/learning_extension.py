"""Matched supervised H2 distance-message networks; frozen quantum labels only.

H2 has one distinct pair distance: this graph experiment is radial regression,
not evidence of chemical-space transfer or a calibrated uncertainty model.
"""
from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import math
from pathlib import Path
import platform
import time

import numpy as np
from scipy.stats import spearmanr
import torch
from torch import nn

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "quantumequi/results/electronic/h2_reference.csv"
SPLIT = ROOT / "quantumequi/results/electronic/surrogate_split.json"
BASELINES = ROOT / "quantumequi/results/electronic/surrogate_metrics.csv"
OBJECTIVES = ("energy_only", "energy_gradient")
MAIN_SEEDS = (7301, 7302, 7303)
MAIN_EPOCHS = 800


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def state_digest(state):
    digest = hashlib.sha256()
    for key in sorted(state):
        value = state[key].detach().cpu().contiguous()
        digest.update(key.encode())
        digest.update(str(value.dtype).encode())
        digest.update(np.asarray(value.shape, dtype=np.int64).tobytes())
        digest.update(value.numpy().tobytes())
    return digest.hexdigest()


def dump(path, obj):
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")


def rows_write(path, rows):
    if not rows:
        return
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def read_reference(path=DATA, split_path=SPLIT):
    with Path(path).open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    split = json.loads(Path(split_path).read_text(encoding="utf-8"))["split_point_ids"]
    if len({r["point_id"] for r in rows}) != len(rows):
        raise ValueError("Duplicate reference point ID")
    for row in rows:
        for key in ("R_A", "energy_total_Hartree", "dE_dR_Hartree_A"):
            row[key] = float(row[key])
            if not np.isfinite(row[key]):
                raise ValueError("Nonfinite reference")
        if row["R_A"] <= 0 or row["SCF_converged"].lower() != "true":
            raise ValueError("Invalid or unconverged reference")
        if row["point_id"] not in split.get(row["split"], []):
            raise ValueError("Reference split differs from frozen partition")
    for group, ids in split.items():
        if set(ids) != {r["point_id"] for r in rows if r["split"] == group}:
            raise ValueError("Reference has missing partition members")
    return rows


def fit_preprocessing(rows):
    train = [row for row in rows if row["split"] == "train"]
    if len(train) < 2:
        raise ValueError("At least two training points required")
    result = {"training_point_ids": [row["point_id"] for row in train]}
    for name, field in [("distance", "R_A"), ("energy", "energy_total_Hartree"), ("gradient", "dE_dR_Hartree_A")]:
        values = np.array([row[field] for row in train], dtype=float)
        result[name + "_mean"] = float(values.mean())
        result[name + "_std"] = float(values.std(ddof=0))
        if result[name + "_std"] <= 0 or not np.isfinite(values).all():
            raise ValueError("Training scale must be finite and positive")
    return result


def coordinates_from_distances(distances):
    r = torch.as_tensor(distances, dtype=torch.float64)
    if r.ndim != 1 or not torch.isfinite(r).all() or (r <= 0).any():
        raise ValueError("Positive finite 1D distances required")
    zeros = torch.zeros_like(r)
    first = torch.stack([zeros, zeros, -r/2], dim=-1)
    second = torch.stack([zeros, zeros, r/2], dim=-1)
    return torch.stack([first, second], dim=1)


class HydrogenDistanceMPNN(nn.Module):
    """Hydrogen-only scalar message network with summed atomic readout.

    Two directed edges for H2; permutation-invariant sum and distance-only
    messages. No coordinate head, cutoff, pretrained parameters or QM calls.
    """
    def __init__(self, preprocessing, width=16, layers=2):
        super().__init__()
        if width < 2 or layers < 1:
            raise ValueError("Invalid architecture")
        self.width, self.layers_count = width, layers
        self.hydrogen_embedding = nn.Parameter(torch.randn(width) * 0.1)
        self.messages = nn.ModuleList([nn.Sequential(nn.Linear(2*width+1, width), nn.SiLU(), nn.Linear(width, width), nn.SiLU()) for _ in range(layers)])
        self.updates = nn.ModuleList([nn.Sequential(nn.Linear(2*width, width), nn.SiLU(), nn.Linear(width, width)) for _ in range(layers)])
        self.readout = nn.Sequential(nn.Linear(width, width), nn.SiLU(), nn.Linear(width, 1))
        for key in ("distance_mean", "distance_std", "energy_mean", "energy_std"):
            self.register_buffer(key, torch.tensor(preprocessing[key], dtype=torch.float64))
        self.forward_batches = 0
        self.forward_structures = 0
        self.double()

    def forward(self, coords):
        if coords.ndim != 3 or coords.shape[1:] != (2, 3):
            raise ValueError("This benchmark supports batches of two hydrogen atoms only")
        if not torch.isfinite(coords).all():
            raise ValueError("Nonfinite coordinates")
        batch = len(coords)
        self.forward_batches += 1
        self.forward_structures += batch
        # Edge i <- j with j=1-i. Exclude self edges before taking norms.
        distances = torch.linalg.vector_norm(coords - coords.flip(1), dim=-1, keepdim=True)
        if (distances <= 0).any():
            raise ValueError("Coincident atoms are outside this radial model's differentiable domain")
        edge = (distances - self.distance_mean) / self.distance_std
        node = self.hydrogen_embedding.view(1, 1, -1).expand(batch, 2, -1)
        for message, update in zip(self.messages, self.updates):
            received = message(torch.cat([node, node.flip(1), edge], dim=-1))
            node = node + update(torch.cat([node, received], dim=-1))
        # The fixed reference offset is split equally between the two H atoms.
        atom_energy = self.energy_std * self.readout(node).squeeze(-1) + self.energy_mean / 2
        return atom_energy.sum(dim=1)


def energy_forces(model, coords, create_graph=False):
    x = coords.detach().clone().requires_grad_(True)
    energy = model(x)
    force = -torch.autograd.grad(energy.sum(), x, create_graph=create_graph)[0]
    return energy, force


def radial_predictions(model, distances, create_graph=False):
    energy, forces = energy_forces(model, coordinates_from_distances(distances), create_graph)
    gradient = (forces[:, 0, 2] - forces[:, 1, 2]) / 2
    return energy, gradient


def validation_score(energy_rmse, gradient_rmse, preprocessing):
    return energy_rmse / preprocessing["energy_std"] + gradient_rmse / preprocessing["gradient_std"]


def fit_one(rows, preprocessing, seed, objective, epochs, out):
    if objective not in OBJECTIVES or epochs < 1:
        raise ValueError("Invalid objective or epoch count")
    torch.manual_seed(seed)
    model = HydrogenDistanceMPNN(preprocessing)
    initial = copy.deepcopy(model.state_dict())
    optimizer = torch.optim.Adam(model.parameters(), lr=0.003)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=epochs, eta_min=0.00015)
    train = [row for row in rows if row["split"] == "train"]
    val = [row for row in rows if row["split"] == "validation"]
    train_coords = coordinates_from_distances([row["R_A"] for row in train])
    targets_e = torch.tensor([row["energy_total_Hartree"] for row in train], dtype=torch.float64)
    targets_g = torch.tensor([row["dE_dR_Hartree_A"] for row in train], dtype=torch.float64)
    val_e = np.array([row["energy_total_Hartree"] for row in val])
    val_g = np.array([row["dE_dR_Hartree_A"] for row in val])
    curve, snapshots, best = [], {}, None
    start = time.perf_counter()
    for epoch in range(1, epochs + 1):
        optimizer.zero_grad(set_to_none=True)
        if objective == "energy_gradient":
            predicted_e, forces = energy_forces(model, train_coords, create_graph=True)
            predicted_g = (forces[:, 0, 2] - forces[:, 1, 2]) / 2
            gradient_loss = torch.mean(((predicted_g-targets_g)/preprocessing["gradient_std"])**2)
        else:
            predicted_e = model(train_coords)
            gradient_loss = predicted_e.sum() * 0
        energy_loss = torch.mean(((predicted_e-targets_e)/preprocessing["energy_std"])**2)
        loss = energy_loss + gradient_loss
        if not torch.isfinite(loss):
            raise RuntimeError("Nonfinite training objective")
        loss.backward()
        grad_norm = torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
        optimizer.step()
        selected_lr = optimizer.param_groups[0]["lr"]
        scheduler.step()
        ep, gp = radial_predictions(model, [row["R_A"] for row in val])
        ermse = float(np.sqrt(np.mean((ep.detach().numpy()-val_e)**2)))
        grmse = float(np.sqrt(np.mean((gp.detach().numpy()-val_g)**2)))
        score = validation_score(ermse, grmse, preprocessing)
        improvement = best is None or score < best["score"]
        if improvement:
            best = {"epoch": epoch, "score": score, "state_dict": copy.deepcopy(model.state_dict()),
                    "validation_energy_RMSE_Hartree": ermse, "validation_gradient_RMSE_Hartree_A": grmse}
        if epoch % 200 == 0 or epoch == epochs:
            snapshots[str(epoch)] = copy.deepcopy(model.state_dict())
        curve.append({"objective": objective, "seed": seed, "epoch": epoch,
                      "pre_update_training_energy_normalized_MSE": float(energy_loss.detach()),
                      "pre_update_training_gradient_normalized_MSE": float(gradient_loss.detach()),
                      "pre_update_total_loss": float(loss.detach()), "gradient_norm_before_clip": float(grad_norm),
                      "learning_rate": selected_lr, "validation_energy_RMSE_Hartree": ermse,
                      "validation_gradient_RMSE_Hartree_A": grmse, "selection_score": score,
                      "new_best_checkpoint": improvement})
    training_counts = {"forward_batches": model.forward_batches, "molecular_energy_evaluations": model.forward_structures,
                       "training_batches": epochs, "validation_batches": epochs,
                       "training_geometries_per_batch": len(train), "validation_geometries_per_batch": len(val),
                       "coordinate_gradient_batches": epochs * (2 if objective == "energy_gradient" else 1),
                       "optimizer_steps": epochs}
    checkpoint = {"architecture": {"name": "HydrogenDistanceMPNN", "width": 16, "layers": 2},
                  "preprocessing": preprocessing, "seed": seed, "objective": objective, "epochs": epochs,
                  "initial_state": initial, "selected_state": best["state_dict"], "final_state": copy.deepcopy(model.state_dict()),
                  "periodic_states": snapshots, "selected_epoch": best["epoch"],
                  "source_data_sha256": sha(DATA), "code_sha256": sha(__file__),
                  "selection_rule": "validation energy RMSE/train energy std + validation gradient RMSE/train gradient std"}
    torch.save(checkpoint, out / "checkpoints" / f"{objective}_seed{seed}.pt")
    model.load_state_dict(best["state_dict"])
    record = {"objective": objective, "seed": seed, "epochs": epochs, "selected_epoch": best["epoch"],
              "selected_validation_score": best["score"], "parameter_count": sum(p.numel() for p in model.parameters()),
              "initial_state_sha256": state_digest(initial), "selected_state_sha256": state_digest(best["state_dict"]),
              "fit_seconds": time.perf_counter()-start, "training_counts": training_counts}
    return model, record, curve


def evaluate_rows(model, rows, objective, seed):
    ep, gp = radial_predictions(model, [row["R_A"] for row in rows])
    ep, gp = ep.detach().numpy(), gp.detach().numpy()
    predictions = [{"objective": objective, "seed": seed, "point_id": row["point_id"], "R_A": row["R_A"], "split": row["split"],
                    "reference_energy_Hartree": row["energy_total_Hartree"], "predicted_energy_Hartree": float(e),
                    "reference_gradient_Hartree_A": row["dE_dR_Hartree_A"], "predicted_gradient_Hartree_A": float(g),
                    "energy_error_Hartree": float(e-row["energy_total_Hartree"]),
                    "gradient_error_Hartree_A": float(g-row["dE_dR_Hartree_A"])} for row, e, g in zip(rows, ep, gp)]
    metrics = []
    for split in dict.fromkeys(row["split"] for row in rows):
        part = [row for row in predictions if row["split"] == split]
        er, gr = np.array([row["energy_error_Hartree"] for row in part]), np.array([row["gradient_error_Hartree_A"] for row in part])
        metrics.append({"objective": objective, "seed": seed, "split": split, "points": len(part),
                        "energy_RMSE_Hartree": float(np.sqrt(np.mean(er**2))), "gradient_RMSE_Hartree_A": float(np.sqrt(np.mean(gr**2))),
                        "energy_MAE_Hartree": float(np.mean(np.abs(er))), "gradient_MAE_Hartree_A": float(np.mean(np.abs(gr)))})
    return predictions, metrics


def mathematical_audit(model, objective, seed):
    covariance, differences = [], []
    rng = np.random.default_rng(1907)
    for distance in [0.74, 1.14, 2.30]:
        x = coordinates_from_distances([distance])
        energy, force = energy_forces(model, x)
        e0, f0 = float(energy.detach()[0]), force.detach().numpy()[0]
        q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
        for reflected in [False, True]:
            transform = q.copy()
            if reflected:
                transform[:, 0] *= -1
            moved = (x.numpy()[0] @ transform.T + [2.3, -1.7, 0.8])[::-1].copy()
            ep, fp = energy_forces(model, torch.tensor(moved[None], dtype=torch.float64))
            expected_f = (f0 @ transform.T)[::-1]
            covariance.append({"objective": objective, "seed": seed, "R_A": distance, "determinant": float(np.linalg.det(transform)),
                               "energy_invariance_error_Hartree": abs(float(ep.detach()[0])-e0),
                               "force_covariance_error_Hartree_A": float(np.max(np.abs(fp.detach().numpy()[0]-expected_f))),
                               "net_force_Hartree_A": float(np.max(np.abs(f0.sum(axis=0)))),
                               "torque_Hartree": float(np.max(np.abs(np.cross(x.numpy()[0], f0).sum(axis=0))))})
        for step in [0.01, 0.001, 0.0001, 0.00001]:
            displaced = x.repeat(12, 1, 1).numpy()
            for i in range(6):
                displaced[2*i].ravel()[i] += step
                displaced[2*i+1].ravel()[i] -= step
            with torch.no_grad():
                ep = model(torch.tensor(displaced, dtype=torch.float64)).numpy()
            fd = -(ep[0::2]-ep[1::2])/(2*step)
            differences.append({"objective": objective, "seed": seed, "R_A": distance, "step_A": step,
                                "energy_geometries_evaluated": 12,
                                "maximum_force_FD_error_Hartree_A": float(np.max(np.abs(fd-f0.ravel())))})
    return covariance, differences


def ensemble_tables(predictions):
    grouped = {}
    for row in predictions:
        grouped.setdefault((row["objective"], row["point_id"]), []).append(row)
    per_point = []
    for (objective, point), rows in grouped.items():
        if len(rows) != 3 or len({row["seed"] for row in rows}) != 3:
            raise ValueError("Exactly three distinct seeds required for ensemble diagnostics")
        ref = rows[0]
        e, g = np.array([r["predicted_energy_Hartree"] for r in rows]), np.array([r["predicted_gradient_Hartree_A"] for r in rows])
        per_point.append({"objective": objective, "point_id": point, "R_A": ref["R_A"], "split": ref["split"],
                          "reference_energy_Hartree": ref["reference_energy_Hartree"], "reference_gradient_Hartree_A": ref["reference_gradient_Hartree_A"],
                          "ensemble_energy_mean_Hartree": float(e.mean()), "ensemble_energy_std_Hartree": float(e.std(ddof=1)),
                          "ensemble_gradient_mean_Hartree_A": float(g.mean()), "ensemble_gradient_std_Hartree_A": float(g.std(ddof=1)),
                          "ensemble_energy_abs_error_Hartree": abs(float(e.mean())-ref["reference_energy_Hartree"]),
                          "ensemble_gradient_abs_error_Hartree_A": abs(float(g.mean())-ref["reference_gradient_Hartree_A"]),
                          "n_seeds": 3})
    summaries = []
    for objective in OBJECTIVES:
        for split in ("train", "validation", "test", "ood_stretch"):
            rows = [r for r in per_point if r["objective"] == objective and r["split"] == split]
            row = {"objective": objective, "split": split, "points": len(rows), "spread_scope": "sample SD across three initializations; uncalibrated"}
            for label, unit in [("energy", "Hartree"), ("gradient", "Hartree_A")]:
                errors = np.array([r[f"ensemble_{label}_abs_error_{unit}"] for r in rows])
                spread = np.array([r[f"ensemble_{label}_std_{unit}"] for r in rows])
                correlation = spearmanr(spread, errors).statistic if np.ptp(spread) and np.ptp(errors) else None
                row.update({f"{label}_ensemble_RMSE_{unit}": float(np.sqrt(np.mean(errors**2))),
                            f"{label}_mean_spread_{unit}": float(spread.mean()),
                            f"{label}_spread_error_spearman": None if correlation is None else float(correlation),
                            f"{label}_fraction_abs_error_le_2std": float(np.mean(errors <= 2*spread))})
            summaries.append(row)
    return per_point, summaries


def run(pilot=False):
    torch.set_num_threads(1)
    torch.use_deterministic_algorithms(True)
    start = time.perf_counter()
    out = ROOT / "quantumequi/results/extensions/learning"
    if pilot:
        out /= "pilot"
    out.mkdir(parents=True, exist_ok=True)
    if (out / "summary.json").exists():
        raise RuntimeError("Existing study is frozen; use a separate checkout for a deliberate rerun")
    (out / "checkpoints").mkdir(exist_ok=True)
    (out / "executed_code.py.txt").write_bytes(Path(__file__).read_bytes())
    rows = read_reference()
    preprocessing = fit_preprocessing(rows)
    dump(out / "preprocessing.json", preprocessing)
    epochs, seeds = (200, [7301]) if pilot else (MAIN_EPOCHS, MAIN_SEEDS)
    # Pilot never scores test or OOD labels, and skips post-selection OOD audits.
    scored_rows = [r for r in rows if r["split"] in ("train", "validation")] if pilot else rows
    curves, records, predictions, metrics, covariance, fd = [], [], [], [], [], []
    for seed in seeds:
        for objective in OBJECTIVES:
            model, record, curve = fit_one(rows, preprocessing, seed, objective, epochs, out)
            p, met = evaluate_rows(model, scored_rows, objective, seed)
            if not pilot:
                c, d = mathematical_audit(model, objective, seed)
                covariance.extend(c)
                fd.extend(d)
            record["total_model_forward_batches_including_postselection_checks"] = model.forward_batches
            record["total_molecular_energy_evaluations_including_postselection_checks"] = model.forward_structures
            records.append(record)
            curves.extend(curve)
            predictions.extend(p)
            metrics.extend(met)
            print(json.dumps({"objective": objective, "seed": seed, "selected_epoch": record["selected_epoch"], "fit_seconds": record["fit_seconds"]}), flush=True)
    for seed in seeds:
        pair = [r for r in records if r["seed"] == seed]
        if pair[0]["initial_state_sha256"] != pair[1]["initial_state_sha256"]:
            raise AssertionError("Matched objectives must share initial parameters")
    rows_write(out / "learning_curves.csv", curves)
    rows_write(out / "predictions.csv", predictions)
    rows_write(out / "metrics.csv", metrics)
    rows_write(out / "covariance.csv", covariance)
    rows_write(out / "force_finite_differences.csv", fd)
    if not pilot:
        ensemble, spread = ensemble_tables(predictions)
        rows_write(out / "ensemble_predictions.csv", ensemble)
        rows_write(out / "ensemble_diagnostics.csv", spread)
        with BASELINES.open(encoding="utf-8", newline="") as handle:
            baseline_rows = list(csv.DictReader(handle))
        rows_write(out / "frozen_baseline_metrics.csv", baseline_rows)
    summary = {"scope": "trained scalar distance-message H2 models; one-molecule radial interpolation/extrapolation only",
               "pilot": pilot, "architecture": {"name": "HydrogenDistanceMPNN", "width": 16, "layers": 2,
                                                  "parameter_count": records[0]["parameter_count"], "dtype": "float64"},
               "training": {"seeds": list(seeds), "epochs_per_run": epochs, "batch_size": 17,
                            "optimizer": "Adam", "initial_learning_rate": 0.003,
                            "schedule": "cosine to 0.00015", "gradient_norm_clip": 5.0,
                            "weight_decay": 0.0, "joint_gradient_loss_weight": 1.0,
                            "objective": "MSE(E/sigma_E) plus optional MSE(dE/dR/sigma_gradient); scales from train only",
                            "checkpoint_selection": "joint validation energy and gradient normalized RMSE sum for BOTH objectives",
                            "checkpoint_retention": "initial, best validation, final and every 200 epochs in each run archive"},
               "data_access": {"training_only_preprocessing": True, "validation_gradient_used_for_both_objectives": True,
                               "test_OOD_used_for_selection": False, "pilot_test_OOD_evaluated": False,
                               "split_counts": {s: sum(r["split"] == s for r in rows) for s in ("train", "validation", "test", "ood_stretch")}},
               "runs": records,
               "counts": {"training_runs": len(records), "optimizer_steps": len(records)*epochs,
                          "validation_epochs": len(curves), "main_metric_rows": len(metrics), "prediction_rows": len(predictions),
                          "covariance_rows": len(covariance), "full_Cartesian_FD_rows": len(fd),
                          "molecular_energy_evaluations": sum(r["total_molecular_energy_evaluations_including_postselection_checks"] for r in records),
                          "model_forward_batches": sum(r["total_model_forward_batches_including_postselection_checks"] for r in records),
                          "new_quantum_jobs": 0, "experiments": 0},
               "counts_scope": "Only this run; pilot, original study and unit tests excluded from main counts. Batched geometry evaluations are distinct from Python forward calls.",
               "limits": ["H2 one-distance regression does not establish molecular GNN transfer", "RHF/STO-3G labels retain physical approximation error",
                          "Three-seed spread is not calibrated predictive uncertainty", "No FCI/DFT/Cu/electrode training labels",
                          "Energy-only fitting still uses validation gradient labels and train gradient scale for selection"],
               "code_sha256": sha(__file__),
               "inputs_sha256": {"frozen_H2_reference": sha(DATA), "frozen_split": sha(SPLIT), "frozen_baseline_metrics": sha(BASELINES)},
               "environment": {"python": platform.python_version(), "numpy": np.__version__, "torch": torch.__version__, "threads": torch.get_num_threads()},
               "sources": ["https://proceedings.mlr.press/v70/gilmer17a.html",
                           "https://proceedings.neurips.cc/paper/2017/hash/303ed4c69846ab36c2904d3ba8573050-Abstract.html",
                           "https://docs.pytorch.org/docs/stable/generated/torch.autograd.grad.html"],
               "license_note": "Original compact implementation; method attribution only, not a SchNet reproduction. No model or data download.",
               "wall_seconds": time.perf_counter()-start}
    summary["output_sha256"] = {p.relative_to(out).as_posix(): sha(p) for p in sorted(out.rglob("*"))
                                if p.is_file() and "pilot" not in p.relative_to(out).parts and p.name != "summary.json"}
    dump(out / "summary.json", summary)
    print(json.dumps({"output": str(out), "counts": summary["counts"], "seconds": summary["wall_seconds"]}), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pilot", action="store_true")
    args = parser.parse_args()
    run(args.pilot)
