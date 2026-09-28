"""Validate frozen SynthaPore evidence without training or simulation reruns.

Checks provenance, saved-data arithmetic, publication structure, and QA records.
A passing result establishes artifact consistency, not chemical validation.
"""
from collections import Counter, defaultdict
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import sys
from urllib.parse import unquote
import xml.etree.ElementTree as ET

import numpy as np
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parent
RESULT = ROOT / "results/publication_validation.json"
STEMS = ("Fig1_Equivariance_Denoising", "Fig2_MD_Validation",
         "Fig3_NEB_Validation", "Fig4_Pore_Adsorption")


def main(require_cross_review=False):
    checks = []

    def check(name, ok, detail=None):
        item = {"check": name, "passed": bool(ok)}
        if detail is not None:
            value = str(detail)
            for path, label in ((ROOT, "synthapore"), (REPO, "repository")):
                value = value.replace(str(path).replace("\\", "\\\\"), label).replace(str(path), label)
            item["detail"] = value
        checks.append(item)

    def guarded(name, function):
        try:
            function()
        except Exception as exc:
            check(name + " readable and valid schema", False, f"{type(exc).__name__}: {exc}")

    def read(name):
        return json.loads((ROOT / name).read_text(encoding="utf-8-sig"))

    def rows(name, expected=None):
        with (ROOT / name).open(encoding="utf-8-sig", newline="") as handle:
            values = list(csv.DictReader(handle))
        if expected is not None:
            check(str(name) + " row count", len(values) == expected)
        return values

    def digest(path):
        return hashlib.sha256(Path(path).read_bytes()).hexdigest()

    def hashes(mapping, base, label):
        check(label + " nonempty hash inventory", bool(mapping))
        for name, expected in mapping.items():
            path = base / name
            check(label + " " + name, path.is_file() and digest(path) == expected)

    def near(a, b, tol=1e-10):
        return math.isclose(float(a), float(b), rel_tol=tol, abs_tol=tol)

    def num(row, key):
        return float(row[key])

    def grouped(values, keys):
        output = defaultdict(list)
        for row in values:
            output[tuple(row[k] for k in keys)].append(row)
        return output

    def module_provenance():
        for module in ("equivariant", "dynamics", "pore"):
            for prefix in ("", "pilot/"):
                name = f"results/{module}/{prefix}summary.json"

                def validate_summary(name=name, module=module, prefix=prefix):
                    summary = read(name)
                    base = ROOT / f"results/{module}/{prefix}"
                    check(name + " pilot label", summary["pilot"] == bool(prefix))
                    mapping = dict(summary["source_sha256"])
                    if module == "equivariant" and prefix:
                        snapshot = summary["archived_source_snapshot"]
                        old = mapping.pop("synthapore/scripts/equivariant_reviewed.py")
                        path = base / snapshot["file"]
                        check(name + " historical source preserved without rewriting its hash",
                              old == snapshot["sha256"] and path.is_file() and digest(path) == old
                              and bool(snapshot["reason"]))
                    hashes(mapping, REPO, name + " source")
                    hashes(summary["output_sha256"], base, name + " output")
                guarded(name, validate_summary)
        supplement = read("results/equivariant/captured_source_audit.json")
        hashes(supplement["source_sha256"], REPO, "captured model supplementary audit source")
        hashes(supplement["output_sha256"], ROOT / "results/equivariant", "captured model supplementary output")

    def source():
        record = read("source/source_record.json")
        original = read("results/original/execution.json")
        compat = read("results/compatibility/execution.json")
        audit = read("results/source_audit.json")
        check("archived specification bytes", digest(ROOT / "source/specification.md") == record["specification_sha256"])
        check("original source and execution bytes", digest(ROOT / "source/synthapore_engine.py") ==
              record["original_python_sha256"] == original["source_sha256"])
        check("compatibility source execution and audit bytes", digest(ROOT / "source/synthapore_engine_compat.py") ==
              record["compatibility_python_sha256"] == compat["source_sha256"] == audit["source_sha256"])
        for mode, execution in (("original", original), ("compatibility", compat)):
            hashes(execution["output_sha256"], ROOT / "results" / mode, mode + " execution output")
        check("original failure and repaired completion retained", original["exit_code"] == 1 and compat["exit_code"] == 0)
        check("execution uses no physical instrument", not any(r["physical_instrument_connected"] for r in (record, original, compat)))
        first = (ROOT / "source/synthapore_engine.py").read_bytes()
        second = (ROOT / "source/synthapore_engine_compat.py").read_bytes()
        before = b"mlip = NeuralInteratomicPotential(num_species=10, hidden_dim=48, num_layers=3)"
        after = before.replace(b"num_species=10", b"num_species=30")
        check("exactly one species-table execution repair", len(record["repairs"]) == 1 and
              first.count(before) == 1 and first.replace(before, after) == second)
        check("embedding repair RNG caveat explicit", "random" in record["random_initialization_caveat"].lower())
        check("source audit captured result and weight identity", digest(ROOT / "results/compatibility/captured_results.json") ==
              audit["captured_sha256"] and digest(ROOT / "results/compatibility/untrained_source_weights.pt") == audit["weights_sha256"])
        check("source audit script bytes", digest(ROOT / "scripts/audit_source.py") == audit["script_sha256"])
        hashes(audit["output_sha256"], ROOT, "source audit output")
        c = audit["counts"]
        check("source audit bounded costs", c == {"source_NEB_replays": 1, "NEB_iterations": 45,
              "NEB_images": 7, "NEB_replay_energy_force_calls": 315, "final_band_energy_force_calls": 7,
              "saved_MD_frames_audited": 30, "new_MD_trajectories": 0})
        md = rows("results/source_audit/MD_telemetry_audit.csv", 30)
        bands = rows("results/source_audit/NEB_final_images.csv", 7)
        replay = rows("results/source_audit/NEB_replay_evaluations.csv", 315)
        check("source replay 45 iterations by seven images", Counter(int(r["iteration"]) for r in replay) == {i: 7 for i in range(45)})
        check("source MD timestamp offset retained", all(near(num(r, "actual_integration_time_fs"), num(r, "source_time_label_fs") + .5) for r in md)
              and near(md[0]["actual_integration_time_fs"], audit["MD"]["first_actual_frame_fs"])
              and near(md[-1]["actual_integration_time_fs"], audit["MD"]["last_actual_frame_fs"]))
        check("source temperature mean from 30 saved samples", near(np.mean([num(r, "source_temperature_before_scaling_K") for r in md]),
              audit["MD"]["recomputed_mean_K"]) and near(round(audit["MD"]["recomputed_mean_K"], 2), audit["MD"]["reported_mean_K"]))
        check("source stale band energies retained", near(max(abs(num(r, "energy_change_nominal_eV")) for r in bands),
              audit["NEB"]["max_stale_energy_difference_nominal_eV"]) and all(near(num(r, "energy_change_nominal_eV"),
              num(r, "reevaluated_postupdate_energy_nominal_eV") - num(r, "reported_preupdate_energy_nominal_eV")) for r in bands))
        check("source final force and pair-distance diagnostics retained", near(min(num(r, "min_pair_distance_A") for r in bands),
              audit["NEB"]["minimum_final_band_pair_distance_A"]) and near(max(num(r, "projected_max_atomic_force_nominal_eV_A")
              for r in bands if r["projected_max_atomic_force_nominal_eV_A"]), audit["NEB"]["final_max_projected_atomic_force_nominal_eV_A"]))
        check("source NEB convergence and free-energy limitations explicit", audit["NEB"]["exact_source_return_matches_capture"]
              and not any(audit["NEB"][k] for k in ("source_has_force_tolerance", "source_has_hessian_check", "barrier_is_free_energy", "source_NEB_uses_MD_external_field")))
        check("untrained source has no diffusion or experimental evidence", not any(audit["source_claims"][k] for k in (
              "trained_parameters", "diffusion_noise_schedule", "reverse_diffusion_sampler", "atomic_framework_coordinates", "DFT_or_experiments")))
        check("source atom inventory and graph declared", audit["input"]["atomic_number_counts"] == {"6": 4, "7": 1, "1": 2, "29": 1}
              and audit["input"]["directed_edges"] == 56 and not audit["input"]["charge_spin_specified"])

    def pilot_accounting():
        base = "results/equivariant/pilot/"
        summary = read(base + "summary.json")
        c = summary["counts"]
        check("pilot denoising cost partition", [c[k] for k in ("neural_training_runs", "training_epochs", "synthetic_shapes",
              "noisy_nodes", "metric_rows", "per_shape_metric_rows", "predicted_coordinate_rows")] == [1, 2, 368, 2944, 16, 1472, 11776])
        for file, count in (("training_losses.csv", 2), ("denoising_metrics.csv", 16), ("per_shape_metrics.csv", 1472),
                            ("predicted_coordinates.csv", 11776), ("equivariance_probes.csv", 22), ("force_finite_differences.csv", 6)):
            rows(base + file, count)
        check("pilot and main use the identical fixed denoising dataset", digest(ROOT / base / "synthetic_shapes.json") ==
              digest(ROOT / "results/equivariant/synthetic_shapes.json"))
        base = "results/dynamics/pilot/"
        c = read(base + "summary.json")["counts"]
        nve = rows(base + "nve_timestep_summary.csv", 1)
        thermostats = rows(base + "thermostat_replicas.csv", 3)
        neb = rows(base + "neb_summary.csv", 1)
        source_neb = rows(base + "source_neb_summary.csv", 2)
        check("pilot dynamics costs separate from main", sum(int(r["integration_steps"]) for r in nve+thermostats) == c["MD_integration_steps"] == 700
              and sum(int(r["energy_force_evaluations"]) for r in nve+thermostats) == c["MD_energy_force_evaluations"] == 704
              and sum(int(r["iterations_executed"]) for r in neb) == c["reviewed_NEB_iterations"] == 135
              and sum(int(r["energy_force_point_evaluations"]) for r in neb) == c["reviewed_NEB_energy_force_point_evaluations"] == 952
              and sum(int(r["iterations_executed"]) for r in source_neb) == c["source_NEB_iterations"] == 46
              and sum(int(r["source_point_energy_force_calls"]) for r in source_neb) == c["source_NEB_point_energy_force_calls"] == 322)
        base = "results/pore/pilot/"
        c = read(base + "summary.json")["counts"]
        geometry = rows(base + "geometry_grid_convergence.csv", 36)
        fits = rows(base + "BET_window_fits.csv", 24)
        noise_fits = rows(base + "noise_BET_fits.csv", 16)
        observations = rows(base + "noise_observations.csv", 48)
        fields = {(r["geometry"], int(r["resolution"])) for r in geometry}
        check("pilot pore costs separate from main", len(fields) == c["clearance_fields"] == 6
              and sum(n*n for _, n in fields) == c["clearance_point_evaluations"] == 24000
              and c["total_BET_fits"] == len(fits)+len(noise_fits) == 40
              and len({r["seed_index"] for r in observations}) == c["noise_seeds"] == 8)

    def equivariant():
        base = "results/equivariant/"
        summary, pilot = read(base + "summary.json"), read(base + "pilot/summary.json")
        c = summary["counts"]
        check("main denoising accounting", [c[k] for k in ("neural_training_runs", "training_epochs", "synthetic_shapes",
              "noisy_nodes", "metric_rows", "per_shape_metric_rows", "predicted_coordinate_rows", "DFT", "new_experiments")] ==
              [3, 180, 368, 2944, 24, 2208, 17664, 0, 0])
        check("pilot plus main four fits and 182 epochs", c["neural_training_runs"] + pilot["counts"]["neural_training_runs"] == 4
              and c["training_epochs"] + pilot["counts"]["training_epochs"] == 182)
        data = read(base + "synthetic_shapes.json")
        by_id = {r["shape_id"]: r for r in data}
        check("368 unique shape IDs with disjoint partitions", len(data) == len(by_id) == 368 and
              Counter(r["split"] for r in data) == summary["data"]["split_sizes"] == {"train": 192, "validation": 48, "test": 64, "ood_warped": 64})
        clean_hashes = [hashlib.sha256(np.asarray(r["clean"], dtype=np.float64).tobytes()).hexdigest() for r in data]
        check("no exact clean-coordinate duplicates across or within partitions", len(set(clean_hashes)) == len(data))
        check("four independent data seeds and three noise levels", len(set(summary["data"]["seeds"].values())) == 4 and
              set(r["noise_sigma"] for r in data) == {.08, .16, .24})
        check("centered noise and finite eight-node coordinates", all(np.asarray(r["clean"]).shape == np.asarray(r["noisy"]).shape == (8, 3)
              and np.isfinite(np.asarray(r["clean"]) + np.asarray(r["noisy"])).all()
              and np.max(np.abs(np.asarray(r["clean"]).mean(0) - np.asarray(r["noisy"]).mean(0))) < 1e-12 for r in data))
        coeff = read(base + "baseline_coefficients.json")
        for sigma, alpha in coeff["noise_level_to_alpha"].items():
            group = [r for r in data if r["split"] == "train" and near(r["noise_sigma"], sigma)]
            x = np.array([r["noisy"] for r in group])
            y = np.array([r["clean"] for r in group])
            lap = (np.roll(x, 1, axis=1) + np.roll(x, -1, axis=1)) / 2 - x
            expected = np.sum(lap * (y - x)) / np.sum(lap**2)
            check("smoother coefficient fitted only on train at sigma " + sigma, coeff["fit_partition"] == "train" and near(alpha, expected))
        losses = rows(base + "training_losses.csv", 180)
        check("three fixed training seeds and 60 recorded epochs", Counter(r["seed"] for r in losses) == {str(i): 60 for i in (4441, 4442, 4443)}
              and all({int(r["epoch"]) for r in group} == set(range(1, 61)) for group in grouped(losses, ["seed"]).values()))
        for fit in summary["training"]["fits"]:
            group = [r for r in losses if int(r["seed"]) == fit["seed"]]
            best = min(group, key=lambda r: num(r, "validation_coordinate_MSE"))
            check(f"validation-only checkpoint selection seed {fit['seed']}", int(best["epoch"]) == fit["selected_epoch"] and
                  near(best["validation_coordinate_MSE"], fit["validation_coordinate_MSE"]) and fit["epochs"] == 60 and fit["parameters"] == 13176)
        metrics = rows(base + "denoising_metrics.csv", 24)
        per_shape = rows(base + "per_shape_metrics.csv", 2208)
        prediction_rows = rows(base + "predicted_coordinates.csv", 17664)
        preds = grouped(prediction_rows, ["model", "shape_id"])
        model_names = {"identity", "PCA_plane_projection", "training_fitted_cycle_smoother"} | {f"EGNN_seed{i}" for i in (4441, 4442, 4443)}
        check("six models each predict every unique shape once", len(preds) == 2208 and {k[0] for k in preds} == model_names
              and all(len(v) == 8 and {int(r["node"]) for r in v} == set(range(8)) for v in preds.values())
              and all({k[1] for k in preds if k[0] == model} == set(by_id) for model in model_names))
        recomputed = {}
        prediction_arrays = {}
        upper = np.triu_indices(8, k=1)
        for (model, sid), group in preds.items():
            p = np.array([[num(r, axis) for axis in ("x", "y", "z")] for r in sorted(group, key=lambda r: int(r["node"]))])
            y = np.asarray(by_id[sid]["clean"])
            # Stored neural coordinates were produced as float32; reproduce the
            # original pair-distance arithmetic, whose norms precede promotion.
            if model.startswith("EGNN_"):
                p = p.astype(np.float32)
            d_pred = np.linalg.norm(p[:, None] - p[None, :], axis=-1)[upper]
            d_true = np.linalg.norm(y[:, None] - y[None, :], axis=-1)[upper]
            recomputed[model, sid] = (np.sqrt(np.mean((p-y)**2)), np.sqrt(np.mean((d_pred-d_true)**2)), np.sqrt(np.mean((p-y).mean(0)**2)))
            prediction_arrays[model, sid] = p
        metric_keys = ("coordinate_RMSE", "pair_distance_RMSE", "centroid_RMS_error")
        check("all 2208 per-shape metrics recompute from saved coordinates", all(r["split"] == by_id[r["shape_id"]]["split"]
              and near(r["noise_sigma"], by_id[r["shape_id"]]["noise_sigma"])
              and all(near(r[k], x) for k, x in zip(metric_keys, recomputed[r["model"], r["shape_id"]])) for r in per_shape))
        for row in metrics:
            group = [values for (model, sid), values in recomputed.items() if model == row["model"] and by_id[sid]["split"] == row["split"]]
            check("aggregate RMS metrics " + row["model"] + " " + row["split"], len(group) == int(row["shapes"])
                  and all(near(row[k], v) for k, v in zip(metric_keys, np.sqrt(np.mean(np.asarray(group)**2, axis=0)))))
        baseline_ok = {"identity": True, "PCA_plane_projection": True, "training_fitted_cycle_smoother": True}
        for sid, datum in by_id.items():
            noisy = np.asarray(datum["noisy"])
            center = noisy.mean(0)
            _, _, axes = np.linalg.svd(noisy-center, full_matrices=False)
            expected = noisy - ((noisy-center) @ axes[-1])[:, None] * axes[-1]
            alpha = coeff["noise_level_to_alpha"][str(datum["noise_sigma"])]
            smoother = noisy + alpha * ((np.roll(noisy, 1, axis=0) + np.roll(noisy, -1, axis=0))/2-noisy)
            for model, target in (("identity", noisy), ("PCA_plane_projection", expected), ("training_fitted_cycle_smoother", smoother)):
                baseline_ok[model] &= np.allclose(prediction_arrays[model, sid], target, rtol=1e-10, atol=1e-10)
        for model, ok in baseline_ok.items():
            check("saved baseline independently reconstructs " + model, ok)
        indexed = {(r["model"], r["split"]): num(r, "coordinate_RMSE") for r in metrics}
        check("negative OOD outcomes retained", indexed["PCA_plane_projection", "ood_warped"] > indexed["identity", "ood_warped"]
              and indexed["EGNN_seed4441", "ood_warped"] > indexed["training_fitted_cycle_smoother", "ood_warped"])
        audit = read(base + "source_and_reviewed_audit.json")
        probes = rows(base + "equivariance_probes.csv", 22)
        fd = rows(base + "force_finite_differences.csv", 6)
        check("small audit exact evaluation accounting", audit["counts"] == c["energy_audit"] and
              c["energy_audit"]["potential_energy_force_evaluations"] == {"source": 123, "reviewed": 121})
        for model, reference in audit["references"].items():
            shift = -reference["charge_sum"] * np.dot(audit["field_vector"], audit["translation_vector"])
            check(model + " charge and energy translation accounting", near(sum(reference["charges"]), reference["charge_sum"])
                  and near(reference["translation_expected_energy_shift_minus_Q_E_dot_t"], shift)
                  and near(reference["translated_energy"] - reference["energy"], reference["observed_energy_shift"]))
        check("reviewed fixed neutral charge and force covariance", abs(audit["references"]["reviewed"]["charge_sum"]) < 1e-12
              and all(num(r, "force_max_error") < 1e-12 and num(r, "energy_error") < 1e-12 for r in probes if r["model"] == "reviewed"))
        check("source cutoff leak and double-conversion failure retained", audit["cutoff"]["source_max_scalar_change_from_edges_at_distance_8"] > .01
              and audit["cutoff"]["reviewed_max_scalar_change_from_edges_at_distance_8"] == 0
              and audit["cutoff"]["reviewed_max_coordinate_change_from_edges_at_distance_8"] == 0 and not audit["dtype"]["source_double_succeeded"])
        check("reviewed central finite difference agreement", all(num(r, "force_max_abs_error") < 1e-6 for r in fd if r["model"] == "reviewed"))
        check("trained double-precision geometric covariance", len(summary["trained_covariance_checks"]) == 3 and
              all(r["reflection_rotation_translation_max_error"] < 1e-12 for r in summary["trained_covariance_checks"]))
        exact = read(base + "captured_source_audit.json")
        rows(base + "captured_source_transformations.csv", 11)
        rows(base + "captured_source_finite_differences.csv", 3)
        check("supplement audits actual 55394-parameter captured model", (exact["num_species"], exact["hidden_dim"], exact["layers"], exact["model_parameter_count"])
              == (30, 48, 3, 55394) and exact["counts"]["energy_force_forwards"] == 157
              and exact["counts"]["training_runs"] == exact["counts"]["MD_trajectories"] == exact["counts"]["NEB_runs"] == 0)
        check("captured float32 charge and field translation identity", near(sum(exact["per_atom_charges"]), exact["unconstrained_total_charge"], 1e-7)
              and near(-exact["unconstrained_total_charge"] * np.dot(exact["field_vector"], exact["translation"]), exact["expected_shift_minus_Q_E_dot_t"], 1e-7)
              and near(exact["observed_energy_shift"], exact["expected_shift_minus_Q_E_dot_t"], 1e-6)
              and exact["origin_force_identity_max_error"] < 1e-7 and exact["origin_force_change_max"] > .001)

    def dynamics():
        base = "results/dynamics/"
        summary = read(base + "summary.json")
        c = summary["counts"]
        nve = rows(base + "nve_timestep_summary.csv", 20)
        thermostats = rows(base + "thermostat_replicas.csv", 24)
        therm_summary = rows(base + "thermostat_summary.csv", 3)
        orders = rows(base + "nve_convergence_orders.csv", 16)
        all_md = nve + thermostats
        check("MD steps and force evaluation accounting", sum(int(r["integration_steps"]) for r in all_md) == c["MD_integration_steps"] == 350000
              and sum(int(r["energy_force_evaluations"]) for r in all_md) == c["MD_energy_force_evaluations"] == 350044
              and c["NVE_trajectories"] == 20 and c["thermostat_trajectories"] == 24)
        kb = summary["units"]["kB_eV_K"]
        check("harmonic frequency and mass conversion", near(summary["units"]["mass_conversion_eV_fs2_A2_per_amu"], 103.64269652680505)
              and near(summary["harmonic_reference"]["normal_mode_frequency_fs_1"], math.sqrt(2/(12.011*103.64269652680505))))
        check("MD DOF, burn-in and canonical reference balances", all(int(r["degrees_of_freedom"]) == 21
              and near(num(r, "integration_steps") * num(r, "timestep_fs"), r["duration_fs"])
              and int(r["energy_force_evaluations"]) == int(r["integration_steps"]) + 1
              and near(num(r, "production_samples_correlated") * num(r, "timestep_fs"), num(r, "duration_fs") - num(r, "burn_in_fs"))
              and near(r["canonical_reference_mean_kinetic_eV"], 21*kb*300/2)
              and near(r["canonical_reference_variance_kinetic_eV2"], 21*(kb*300)**2/2)
              and near(r["canonical_reference_variance_temperature_K2"], 2*300**2/21) for r in all_md))
        by_seed_dt = {(r["seed"], float(r["timestep_fs"])): r for r in nve}
        states = read(base + "nve_initial_final_states.json")
        check("20 NVE coordinate and velocity snapshots", len(states) == 20)
        for state in states:
            s = state["summary"]
            q0, v0, q1, v1 = [np.asarray(state[k]) for k in ("initial_position_A", "initial_velocity_A_fs", "final_position_A", "final_velocity_A_fs")]
            omega, time = s["omega_fs_1"], s["duration_fs"]
            phase = omega*time
            exact_q = q0*np.cos(phase) + v0/omega*np.sin(phase)
            exact_v = -q0*omega*np.sin(phase) + v0*np.cos(phase)
            mass = s["mass_amu"]*summary["units"]["mass_conversion_eV_fs2_A2_per_amu"]
            energy = lambda q, v: s["spring_eV_A2"]*np.sum((q-q.mean(0))**2)/2 + mass*np.sum(v*v)/2
            check("independent exact harmonic solution " + state["case"], near(np.sqrt(np.mean((q1-exact_q)**2)), s["final_position_error_A_vs_exact_NVE"])
                  and near(np.sqrt(np.mean((v1-exact_v)**2)), s["final_velocity_error_A_fs_vs_exact_NVE"])
                  and near(energy(q0, v0), s["initial_total_eV"]) and near(energy(q1, v1), s["final_total_eV"]))
        for order in orders:
            fine = by_seed_dt[order["seed"], float(order["fine_dt_fs"])]
            coarse = by_seed_dt[order["seed"], float(order["coarse_dt_fs"])]
            expected_energy = math.log(num(coarse, "max_abs_total_energy_change_eV")/num(fine, "max_abs_total_energy_change_eV"), 2)
            expected_position = math.log(num(coarse, "final_position_error_A_vs_exact_NVE")/num(fine, "final_position_error_A_vs_exact_NVE"), 2)
            check("saved NVE convergence order " + order["seed"] + " " + order["fine_dt_fs"], near(order["energy_error_order"], expected_energy)
                  and near(order["position_error_order"], expected_position) and 1.9 < expected_position < 2.1)
        for file in sorted((ROOT / base).glob("*_trajectory.csv")):
            trajectory = rows(file.relative_to(ROOT))
            check(file.name + " saved energy and temperature balances", all(near(num(r, "potential_eV") + num(r, "kinetic_eV"), r["total_eV"])
                  and near(r["temperature_COM_dof_K"], 2*num(r, "kinetic_eV")/(21*kb))
                  and near(r["temperature_3N_K"], 2*num(r, "kinetic_eV")/(24*kb)) for r in trajectory)
                  and all(num(b, "time_fs") > num(a, "time_fs") for a, b in zip(trajectory, trajectory[1:])))
        for item in therm_summary:
            group = [r for r in thermostats if r["thermostat"] == item["thermostat"]]
            mappings = {"mean_temperature_K": "replicate_mean_temperature_K", "variance_temperature_K2": "mean_within_trajectory_temperature_variance_K2",
                        "mean_kinetic_eV": "replicate_mean_kinetic_eV", "mean_potential_eV": "replicate_mean_potential_eV"}
            check("thermostat replica summaries " + item["thermostat"], len(group) == int(item["replicas"]) == 8
                  and all(near(np.mean([num(r, source) for r in group]), item[target]) for source, target in mappings.items())
                  and near(item["variance_ratio_to_canonical"], num(item, "mean_within_trajectory_temperature_variance_K2")/(2*300**2/21)))
        therm_index = {r["thermostat"]: r for r in therm_summary}
        check("Berendsen suppressed fluctuations and wrong-DOF temperature retained", num(therm_index["berendsen"], "variance_ratio_to_canonical") < 1e-5
              and num(therm_index["berendsen_source_dof"], "replicate_mean_temperature_K") > 340)
        neb = rows(base + "neb_summary.csv", 18)
        bands = grouped(rows(base + "neb_final_bands.csv", 210), ["case"])
        history = grouped(rows(base + "neb_force_history.csv", 3147), ["case"])
        check("reviewed NEB counts and point costs", sum(int(r["iterations_executed"]) for r in neb) == c["reviewed_NEB_iterations"] == 3129
              and sum(int(r["energy_force_point_evaluations"]) for r in neb) == c["reviewed_NEB_energy_force_point_evaluations"] == 37337
              and c["reviewed_NEB_cases"] == 18 and c["reviewed_NEB_endpoint_check_point_evaluations"] == 36)
        for row in neb:
            name = row["case"]
            group = sorted(bands[name,], key=lambda r: int(r["image_index"]))
            trace = history[name,]
            energies = [(num(r, "x_A")**2-1)**2 + 4*(num(r, "y_A")-.65*(1-num(r, "x_A")**2))**2 for r in group]
            ci = int(row["climbing_index"])
            x, y = num(group[ci], "x_A"), num(group[ci], "y_A")
            t = y-.65*(1-x*x)
            grad = np.array([4*x*(x*x-1)+10.4*x*t, 8*t])
            hessian = np.array([[12*x*x-4+10.4*t+13.52*x*x, 10.4*x], [10.4*x, 8]])
            eig = np.linalg.eigvalsh(hessian)
            check("analytic final-band energies " + name, len(group) == int(row["n_images"])
                  and all(near(value, r["energy_eV"]) for value, r in zip(energies, group))
                  and near(energies[ci]-energies[0], row["barrier_eV"]) and near(num(row, "barrier_eV")-1, row["barrier_error_eV"]))
            check("actual final-force and saddle convergence " + name, row["converged"] == "True" and row["stop_reason"] == "force_tolerance"
                  and num(row, "max_NEB_force_eV_A") < num(row, "tolerance_eV_A")
                  and near(np.linalg.norm(grad), row["climbing_true_force_eV_A"])
                  and near(np.linalg.norm([x, y-.65]), row["saddle_coordinate_error_A"])
                  and np.allclose(eig, json.loads(row["saddle_hessian_eigenvalues_eV_A2"]), rtol=1e-10, atol=1e-10)
                  and sum(eig < 0) == int(row["saddle_negative_hessian_eigenvalues"]) == 1)
            check("NEB trace final-state consistency " + name, len(trace) == int(row["iterations_executed"])+1
                  and [int(r["iteration"]) for r in trace] == list(range(len(trace)))
                  and near(trace[-1]["max_NEB_force_eV_A"], row["max_NEB_force_eV_A"])
                  and near(trace[-1]["barrier_eV"], row["barrier_eV"]) and row["potential_energy_not_free_energy"] == "True")
        source_neb = rows(base + "source_neb_summary.csv", 12)
        check("source NEB costs and false convergence counted", sum(int(r["iterations_executed"]) for r in source_neb) == c["source_NEB_iterations"] == 618
              and sum(int(r["source_point_energy_force_calls"]) for r in source_neb) == c["source_NEB_point_energy_force_calls"] == 7210
              and sum(int(r["n_images"]) for r in source_neb) == c["source_NEB_postupdate_audit_point_evaluations"] == 140
              and sum(r["passes_independent_1e_5_force_check"] == "False" for r in source_neb) == summary["checks"]["source_false_convergence_declarations"] == 10)
        endpoints = read(base + "endpoint_minimization.json")
        check("stationary endpoint optimization accounting", len(endpoints) == c["endpoint_optimizations"] == 2
              and sum(r["energy_force_point_evaluations"] for r in endpoints) == c["endpoint_energy_force_point_evaluations"] == 778
              and all(r["converged"] and r["force_norm_eV_A"] < 1e-10 and min(r["hessian_eigenvalues_eV_A2"]) > 0 for r in endpoints))
        check("analytic dynamics has no DFT or experiment", c["DFT_calculations"] == c["experimental_runs"] == 0)

    def pore():
        base = "results/pore/"
        summary = read(base + "summary.json")
        c, adsorption = summary["counts"], summary["adsorption"]
        geometry = rows(base + "geometry_grid_convergence.csv", 90)
        cdf = rows(base + "clearance_CDF.csv", 123)
        hist = rows(base + "clearance_histogram.csv", 120)
        bias = rows(base + "endpoint_bias.csv", 20)
        iso = rows(base + "source_isotherm_points.csv", 20)
        fits = rows(base + "BET_window_fits.csv", 24)
        noise = rows(base + "noise_observations.csv", 384)
        noise_fits = rows(base + "noise_BET_fits.csv", 128)
        noise_summary = rows(base + "noise_summary.csv", 6)
        fields = {(r["geometry"], int(r["resolution"])) for r in geometry}
        check("pore grid and inverse-fit accounting", len(fields) == c["clearance_fields"] == 15
              and sum(n*n for _, n in fields) == c["clearance_point_evaluations"] == 1636800
              and c["probe_masks"] == len(geometry) and c["total_BET_fits"] == len(fits)+len(noise_fits) == 152
              and c["noise_seeds"] == len({r["seed_index"] for r in noise}) == 64)
        geometry_ok = []
        for row in geometry:
            radius, probe, cell = (num(row, k) for k in ("channel_inradius_or_disk_radius_A", "probe_A", "cell_A"))
            if row["geometry"] == "circle":
                exact = math.pi*max(radius-probe, 0)**2
            elif row["geometry"] == "hexagon":
                exact = 2*math.sqrt(3)*max(radius-probe, 0)**2
            else:
                exact = cell**2 - math.pi*(radius+probe)**2
            estimate = num(row, "accessible_points")/num(row, "grid_points")*cell**2
            geometry_ok.append(near(row["analytic_area_A2"], exact) and near(row["midpoint_area_A2"], estimate)
                and near(row["grid_accessible_fraction"], estimate/cell**2) and near(row["analytic_accessible_fraction"], exact/cell**2)
                and near(row["abs_fraction_error"], abs(estimate-exact)/cell**2)
                and int(row["grid_points"]) == int(row["resolution"])**2)
        check("90 accessibility rows reconstruct analytic areas and midpoint counts", all(geometry_ok))
        for key, group in grouped(geometry, ["geometry", "resolution"]).items():
            values = sorted(group, key=lambda r: num(r, "probe_A"))
            check("probe accessibility monotone " + " ".join(key), all(num(a, "accessible_points") >= num(b, "accessible_points") for a, b in zip(values, values[1:])))
        check("clearance CDF bounds and errors", all(0 <= num(r, "empirical_area_weighted_CDF") <= 1
              and 0 <= num(r, "analytic_CDF") <= 1 and near(abs(num(r, "empirical_area_weighted_CDF")-num(r, "analytic_CDF")), r["absolute_CDF_error"]) for r in cdf))
        check("clearance histograms retain normalized point probabilities", all(near(num(r, "count")/num(r, "grid_void_points"), r["empirical_probability"]) for r in hist)
              and all(near(sum(num(r, "empirical_probability") for r in group), 1) for group in grouped(hist, ["geometry"]).values()))
        check("endpoint bias percentage-point arithmetic", all(near((num(r, "estimated_circle_area_fraction")-num(r, "analytic_circle_area_fraction"))*100,
              r["signed_error_percentage_points"]) for r in bias))
        factor = adsorption["nitrogen_cross_section_m2"]*adsorption["Avogadro_mol_1"]/adsorption["STP_molar_volume_cm3_mol"]
        check("BET area unit conversion", near(factor, adsorption["derived_area_factor_m2_per_cm3_STP"]))
        fit_arithmetic = []
        for row in fits + noise_fits:
            slope, intercept = num(row, "slope"), num(row, "intercept")
            vm, constant = 1/(slope+intercept), 1+slope/intercept
            physical = vm > 0 and constant > 0
            area_ok = (near(row["apparent_area_m2_g"], vm*factor) and
                       near(row["source_conversion_area_m2_g"], vm*adsorption["source_area_factor"])) if physical else (
                           row["apparent_area_m2_g"] == row["source_conversion_area_m2_g"] == "")
            fit_arithmetic.append(near(row["fitted_Vm_cm3_g"], vm) and near(row["fitted_C"], constant)
                and area_ok and (row["positive_parameters"] == "True") == physical)
        check("152 BET slope/intercept and area conversions independently recompute", all(fit_arithmetic))
        def regression(points, row):
            selected = [(p, v) for p, v in points if num(row, "lower_requested") <= p <= num(row, "upper_requested")]
            p, v = np.array(selected).T
            y = p/(v*(1-p))
            design = np.column_stack([np.ones(len(p)), p])
            intercept, slope = np.linalg.lstsq(design, y, rcond=None)[0]
            residual = y-(intercept+slope*p)
            return len(p) == int(row["points"]) and near(intercept, row["intercept"]) and near(slope, row["slope"]) and near(
                np.sqrt(np.mean(residual**2)), row["transformed_RMSE"]) and near(1-np.sum(residual**2)/np.sum((y-y.mean())**2), row["transformed_R2"])
        names = {"source_rounded": "source_rounded_volume_cm3_g", "source_unrounded": "source_unrounded_volume_cm3_g", "BET_exact_control": "BET_only_control_cm3_g"}
        check("24 deterministic regressions recompute from saved isotherm points", all(regression([(num(r, "relative_pressure"), num(r, names[row["data"]])) for r in iso], row) for row in fits))
        noise_groups = grouped(noise, ["seed_index"])
        check("128 noisy regressions use recorded observations", all(regression([(num(r, "relative_pressure"), num(r, "noisy_volume_cm3_g"))
              for r in noise_groups[row["seed_index"],]], row) for row in noise_fits))
        for row in noise_summary:
            selection = [r for r in noise_fits if near(r["lower_requested"], row["lower"]) and near(r["upper_requested"], row["upper"])]
            values = np.array([num(r, row["parameter"]) for r in selection])
            check("noise replication quantiles " + row["parameter"] + " " + row["lower"], len(values) == int(row["total_fits"]) == 64
                  and sum(r["positive_parameters"] == "True" for r in selection) == int(row["positive_fits"])
                  and near(values.mean(), row["mean"]) and near(values.std(ddof=1), row["sample_sd"])
                  and all(near(x, row[k]) for x, k in zip(np.quantile(values, [.025, .5, .975]), ("p025", "median", "p975"))))
        check("unphysical fitted BET constants retained", any(num(r, "fitted_C") < 0 and r["positive_parameters"] == "False" for r in fits))
        audit = read(base + "source_pore_audit.json")
        check("source probe argument does not change outputs", audit["class_only_execution_count"] == len(audit["calls"]) == c["source_class_calls"] == 8
              and not audit["probe_argument_read_in_body"] and all(len({r["source_reported_porosity_pct"] for r in group}) == 1
              for group in grouped(audit["calls"], ["cell_A"]).values()))
        check("piecewise synthetic isotherm discontinuity preserved", near(audit["empirical_right_value_cm3_g"]-audit["BET_left_limit_cm3_g"], audit["jump_cm3_g"])
              and audit["jump_cm3_g"] > 200 and not audit["source_type_IV_label_validated"])
        check("pore study has no atomistic adsorption or experiment", c["GCMC_runs"] == c["atomistic_frameworks"] == c["experimental_measurements"] == 0
              and not audit["PSD_computed"] and not audit["atomic_framework_coordinates_supplied"])

    def figures():
        manifest = read("results/figure_manifest.json")["figures"]
        expected = {f"{stem}_{language}.{fmt}" for stem in STEMS for language in ("english", "chinese") for fmt in ("png", "svg")}
        check("16 reviewed figures with specified names", len(manifest) == 16 and {Path(r["file"]).name for r in manifest} == expected
              and Counter(r["format"] for r in manifest) == {"png": 8, "svg": 8})
        for figure in manifest:
            path = ROOT / figure["file"]
            if not path.is_file():
                path = ROOT / "reports/figures" / figure["file"]
            name = Path(figure["file"]).name
            check("figure current bytes " + name, path.is_file() and digest(path) == figure["sha256"])
            hashes(figure["sources"], ROOT, name + " source")
            if figure["format"] == "png":
                with Image.open(path) as im:
                    check("PNG pixels and resolution " + name, list(im.size) == figure["pixels"] == [2400, 1380]
                          and all(abs(value-300) < .1 for value in im.info.get("dpi", (0, 0))))
            else:
                tree = ET.parse(path)
                ns = {"s": "http://www.w3.org/2000/svg"}
                nodes = tree.findall(".//s:text", ns)
                labels = " ".join("".join(node.itertext()) for node in nodes)
                check("editable SVG without embedded raster " + name, len(nodes) > 10 and not tree.findall(".//s:image", ns)
                      and figure["editable_text"] and not figure["embedded_raster"])
                check("visible figure evidence boundary " + name, bool(re.search(
                    r"UNTRAINED|SYNTHETIC|ANALYTIC|ASSUMED|MODEL|NUMERICAL|AUDIT|未训练|合成|解析|假定|假设|模型|数值|审计", labels, re.I)))
        qa = read("results/figure_qa.json")
        actual = {Path(r["file"]).name: r["sha256"] for r in manifest}
        recorded = {Path(name).name: value for name, value in qa["reviewed_figure_hashes"].items()}
        check("visual QA covers all current eight PNG and eight SVG", actual == recorded and qa["png_files_visually_reviewed"] == 8
              and qa["svg_files_rendered_and_visually_reviewed"] == 8)
        check("visual QA identifies current generator", qa["generator_sha256"] == digest(ROOT / "scripts/plot_reviewed.py"))

    def cross_review():
        path = ROOT / "results/cross_review.json"
        if not path.is_file():
            if require_cross_review:
                check("required peer cross-review record present", False)
            return
        record = read("results/cross_review.json")
        check("cross-review records have no unresolved blocking findings", bool(record["reviews"])
              and all(r["blocking_findings"] == [] for r in record["reviews"]))
        for group, mapping in record["reviewed_artifact_sha256"].items():
            hashes(mapping, ROOT, "cross-review current evidence " + group)

    def reports():
        texts = {}
        for language in ("english", "chinese"):
            path = ROOT / "reports" / f"synthapore_report_{language}.md"
            text = path.read_text(encoding="utf-8")
            texts[language] = text
            check(language + " six ordered main sections", re.findall(r"^## (\d+)\.", text, re.M) == list("123456"))
            subsections = re.findall(r"^### (\d+)\.(\d+) ", text, re.M)
            check(language + " subsection numbering follows main sections", all(1 <= int(a) <= 6 and int(b) >= 1 for a, b in subsections)
                  and len(subsections) == len(set(subsections)))
            for number, stem in enumerate(STEMS, 1):
                check(language + f" figure {number} caption and PNG SVG links", ("**" + ("Figure " if language == "english" else "图 ") + str(number)) in text
                      and all(stem+"_"+language+"."+fmt in text for fmt in ("png", "svg")))
            check(language + " no placeholders and balanced math or code fences", not re.search(r"@@[A-Z_]+@@", text)
                  and text.count(chr(96)*3) % 2 == 0 and text.count("\\[") == text.count("\\]") and text.count("\\(") == text.count("\\)"))
        numbers = lambda text: [re.findall(r"(?<![A-Za-z0-9.])[-−]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?", line)
                                for line in text.splitlines() if line.startswith("|")]
        check("separate English report contains no Chinese prose", not re.search(r"[\u4e00-\u9fff]", texts["english"]))
        check("English Chinese numerical table parity", numbers(texts["english"]) == numbers(texts["chinese"]))

    def local_links():
        check("module README present", (ROOT / "README.md").is_file())
        paths = list(ROOT.glob("*.md")) + list((ROOT / "reports").rglob("*.md"))
        check("authored documentation available", bool(paths))
        for path in sorted(paths):
            text = path.read_text(encoding="utf-8")
            for link in re.findall(r"\]\(([^)]+)\)", text):
                if link.startswith(("https://", "http://", "mailto:", "#", "data:")):
                    continue
                name = unquote(link.split("#", 1)[0].strip("<>"))
                target = path.parent / name
                check(path.relative_to(ROOT).as_posix() + " local link " + link, target.exists() or target.resolve() == RESULT.resolve())

    for name, function in (("module provenance", module_provenance), ("source", source), ("pilot accounting", pilot_accounting), ("equivariant", equivariant),
                           ("dynamics", dynamics), ("pore", pore), ("figures", figures), ("cross-review", cross_review),
                           ("reports", reports), ("local links", local_links)):
        guarded(name, function)
    result = {"passed": all(r["passed"] for r in checks), "checks": len(checks), "failed": [r for r in checks if not r["passed"]],
              "scope": "Saved provenance, independent arithmetic from stored coordinates/trajectories/fits, negative outcomes, bilingual structure/table parity, local links and recorded visual QA. No training, MD, NEB, pore simulation or physical validation performed by this validator.",
              "cross_review_required": require_cross_review, "validator_sha256": digest(Path(__file__)), "items": checks}
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(result, indent=2, ensure_ascii=False)+"\n", encoding="utf-8", newline="\n")
    print(json.dumps({k: result[k] for k in ("passed", "checks", "failed")}, indent=2, ensure_ascii=False))
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--require-cross-review", action="store_true", help="Fail if the final peer-review record is absent")
    sys.exit(main(parser.parse_args().require_cross_review))
