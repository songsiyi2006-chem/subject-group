"""Finite-pool benchmark extension: cached oracle calls and held-out surrogate checks.

No new transport solve, experimental run, instrument connection or test-set tuning
occurs here. The immutable 81-point simulation pool is the entire target domain.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
from importlib.metadata import version
import json
import math
import os
from pathlib import Path
import platform
import sys
import time

import numpy as np
from scipy.stats import t as student_t
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, Matern

sys.path.insert(0, str(Path(__file__).resolve().parent))
from metrics_control import hypervolume_2d, sequential_campaign

ROOT = Path(__file__).resolve().parents[2]
CONTROL = ROOT / "electratwin/results/control"
OUT = ROOT / "electratwin/results/benchmark_extension"
TARGETS = ["STY_div_50000", "negative_SEC_kwh_kg"]
METHODS = ["gp_mc_ehvi", "random", "maximin_spacefill"]
MODELS = ["fixed_matern_gp", "training_mean", "quadratic_regression"]
REFERENCE = [0.0, -2.0]


def read_csv(path):
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(path, rows):
    if not rows:
        raise ValueError("Cannot write empty CSV")
    with Path(path).open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_pool():
    rows = read_csv(CONTROL / "candidate_pool.csv")
    indices = [int(r["candidate_index"]) for r in rows]
    if indices != list(range(81)):
        raise ValueError("Expected immutable ordered 81-point pool")
    x = np.array([[float(r["flow_rate_uL_min"]), float(r["overpotential_V"])] for r in rows])
    y = np.array([[float(r["STY_assumed_product_kg_m3_day"]) / 50000,
                   -float(r["SEC_assumed_product_kwh_kg"])] for r in rows])
    prior = json.loads((CONTROL / "optimization_summary.json").read_text(encoding="utf-8"))
    if sha(CONTROL / "candidate_pool.csv") != prior["candidate_pool_sha256"]:
        raise ValueError("Pool no longer matches frozen control provenance")
    if not np.isfinite(x).all() or not np.isfinite(y).all() or len(np.unique(x, axis=0)) != 81:
        raise ValueError("Invalid pool")
    return x, y


class BudgetedCachedOracle:
    """Expose only requested cached rows, and enforce a strict per-run call budget."""
    def __init__(self, objectives, budget):
        self._objectives = np.array(objectives, dtype=float, copy=True)
        self.budget = int(budget)
        self.calls = []

    def __call__(self, index):
        index = int(index)
        if len(self.calls) >= self.budget:
            raise RuntimeError("Oracle call would exceed budget")
        if index in self.calls:
            raise RuntimeError("Repeated oracle call within campaign")
        if not 0 <= index < len(self._objectives):
            raise ValueError("Unknown candidate")
        self.calls.append(index)
        return self._objectives[index].copy()


def maximin_campaign(candidate_x, evaluate, *, budget=25, initial_count=5, seed=0):
    """Sequential maximin distance uses inputs and selected indices, never outcomes."""
    x = np.asarray(candidate_x, dtype=float)
    if x.ndim != 2 or x.shape[1] != 2 or not np.isfinite(x).all():
        raise ValueError("Expected finite (n,2) inputs")
    if len(np.unique(x, axis=0)) != len(x) or not 2 <= initial_count <= budget <= len(x):
        raise ValueError("Duplicate inputs or invalid budget")
    scale = np.ptp(x, axis=0)
    if np.any(scale <= 0):
        raise ValueError("Both input dimensions must vary")
    # The finite campaign domain is known before any objective observations.
    xn = (x - x.min(axis=0)) / scale
    rng = np.random.default_rng(seed)
    initial = rng.choice(len(x), size=initial_count, replace=False).tolist()
    selected, observed, records = [], [], []
    for step in range(budget):
        distance = None
        if step < initial_count:
            index = initial[step]
        else:
            remaining = np.array([i for i in range(len(x)) if i not in selected])
            delta = xn[remaining, None, :] - xn[np.asarray(selected), :][None, :, :]
            min_distances = np.sqrt(np.sum(delta * delta, axis=-1)).min(axis=1)
            best = int(np.argmax(min_distances))  # Stable tie: lowest candidate index.
            index, distance = int(remaining[best]), float(min_distances[best])
        value = np.asarray(evaluate(index), dtype=float)
        if value.shape != (2,) or not np.isfinite(value).all():
            raise ValueError("Invalid oracle response")
        selected.append(index)
        observed.append(value.tolist())
        records.append({"seed": seed, "method": "maximin_spacefill", "evaluation": step + 1,
            "candidate_index": index, "flow_rate_uL_min": float(x[index, 0]), "overpotential_V": float(x[index, 1]),
            "objective_STY_div_50000": float(value[0]), "objective_negative_SEC": float(value[1]),
            "acquisition_MC_EHVI": None, "hypervolume": hypervolume_2d(observed, REFERENCE),
            "phase": "initial_shared" if step < initial_count else "sequential",
            "acquisition_input_maximin_distance": distance})
    return records


def descriptive_t(values):
    values = np.asarray(values, dtype=float)
    if values.ndim != 1 or not len(values) or not np.isfinite(values).all():
        raise ValueError("Expected nonempty finite vector")
    mean = float(values.mean())
    sd = float(values.std(ddof=1)) if len(values) > 1 else None
    margin = float(student_t.ppf(0.975, len(values) - 1) * sd / np.sqrt(len(values))) if sd is not None else None
    return {"n": len(values), "mean": mean, "sample_sd": sd,
            "lower_95pct_t": mean - margin if margin is not None else None,
            "upper_95pct_t": mean + margin if margin is not None else None,
            "wins": int(np.sum(values > 1e-12)), "losses": int(np.sum(values < -1e-12)),
            "ties": int(np.sum(np.abs(values) <= 1e-12))}


def compare_saved_prefix(records):
    """Compare the original eight seeds at budget 15 without rerunning that study."""
    previous = read_csv(CONTROL / "sequential_campaigns.csv")
    lookup = {(int(r["seed"]), r["method"], int(r["evaluation"])): r for r in records}
    comparisons, missing = [], 0
    for old in previous:
        key = (int(old["seed"]), old["method"], int(old["evaluation"]))
        if key not in lookup:
            missing += 1
            continue
        new = lookup[key]
        comparisons.append({"seed": key[0], "method": key[1], "evaluation": key[2],
            "same_candidate": int(new["candidate_index"]) == int(old["candidate_index"]),
            "hypervolume_absolute_difference": abs(float(new["hypervolume"]) - float(old["hypervolume"]))})
    passed = bool(comparisons) and all(r["same_candidate"] and r["hypervolume_absolute_difference"] < 1e-12 for r in comparisons)
    return {"prior_records": len(previous), "compared_records": len(comparisons), "missing_records": missing,
            "all_available_prefixes_match": passed,
            "maximum_HV_absolute_difference": max((r["hypervolume_absolute_difference"] for r in comparisons), default=None),
            "all_original_240_records_reproduced": passed and missing == 0}


def run_campaign_extension(x, y, *, seed_count=64, budget=25, mc_draws=256):
    full_hv = hypervolume_2d(y, REFERENCE)
    records, budget_checks = [], []
    for seed in range(seed_count):
        initial_ids = None
        for method in METHODS:
            oracle = BudgetedCachedOracle(y, budget)
            if method == "maximin_spacefill":
                run = maximin_campaign(x, oracle, budget=budget, initial_count=5, seed=seed)
            else:
                run = sequential_campaign(x, oracle, budget=budget, initial_count=5,
                                          seed=seed, method=method, mc_draws=mc_draws, reference=REFERENCE)
                for row in run:
                    row["acquisition_input_maximin_distance"] = None
            ids = [int(r["candidate_index"]) for r in run]
            if initial_ids is None:
                initial_ids = ids[:5]
            if len(run) != budget or oracle.calls != ids or len(set(ids)) != budget or ids[:5] != initial_ids:
                raise RuntimeError("Campaign budget, uniqueness or shared-initial audit failed")
            for row in run:
                row["full_pool_HV_fraction"] = row["hypervolume"] / full_hv
                row["evidence_role"] = "cached_deterministic_PDE_pool_evaluation_use"
            records.extend(run)
            budget_checks.append({"seed": seed, "method": method, "requested_budget": budget,
                "oracle_calls": len(oracle.calls), "unique_calls": len(set(oracle.calls)),
                "shared_initial_ids": initial_ids.copy(), "passed": True})
    budgets = [b for b in [5, 10, 15, 20, 25] if b <= budget]
    prefixes = [r for r in records if r["evaluation"] in budgets]
    summary, pairs = [], []
    for b in budgets:
        at_budget = [r for r in prefixes if r["evaluation"] == b]
        for method in METHODS:
            f = np.array([r["full_pool_HV_fraction"] for r in at_budget if r["method"] == method])
            summary.append({"budget": b, "method": method, "seed_count": len(f),
                "mean_full_pool_HV_fraction": float(f.mean()), "sample_sd_full_pool_HV_fraction": float(f.std(ddof=1)),
                "min_full_pool_HV_fraction": float(f.min()), "max_full_pool_HV_fraction": float(f.max())})
        indexed = {(r["seed"], r["method"]): float(r["hypervolume"]) for r in at_budget}
        for left, right in [("gp_mc_ehvi", "random"), ("gp_mc_ehvi", "maximin_spacefill"), ("maximin_spacefill", "random")]:
            diff = [indexed[(seed, left)] - indexed[(seed, right)] for seed in range(seed_count)]
            stats = descriptive_t(diff)
            pairs.append({"budget": b, "left_method": left, "right_method": right,
                "difference_direction": "left_minus_right_raw_HV", **stats,
                "interval_scope": "seed_sensitivity_only_fixed_pool_not_physical_uncertainty_or_general_superiority"})
    return {"records": records, "prefixes": prefixes, "learning_summary": summary,
            "paired_summary": pairs, "budget_checks": budget_checks,
            "saved_prefix_replication": compare_saved_prefix(records), "full_pool_HV": full_hv}


def make_holdout_splits(x, random_count=20):
    x = np.asarray(x, dtype=float)
    if x.shape != (81, 2):
        raise ValueError("The preregistered split design requires the 81-point pool")
    splits = []
    all_ids = np.arange(len(x))
    for seed in range(2000, 2000 + random_count):
        permutation = np.random.default_rng(seed).permutation(len(x))
        splits.append({"split_id": f"random_{seed}", "kind": "random_54_27", "seed": seed,
                       "train_ids": np.sort(permutation[:54]).tolist(), "test_ids": np.sort(permutation[54:]).tolist()})
    for axis, name in [(0, "flow"), (1, "eta")]:
        levels = np.unique(x[:, axis])
        if len(levels) != 9:
            raise ValueError("Expected 9 levels per input")
        for block, excluded_levels in enumerate(np.array_split(levels, 3)):
            test = all_ids[np.isin(x[:, axis], excluded_levels)]
            train = all_ids[~np.isin(x[:, axis], excluded_levels)]
            splits.append({"split_id": f"block_{name}_{block}", "kind": f"block_{name}", "seed": None,
                           "excluded_axis_levels": excluded_levels.tolist(), "train_ids": train.tolist(), "test_ids": test.tolist()})
    for split in splits:
        train, test = split["train_ids"], split["test_ids"]
        if len(train) != 54 or len(test) != 27 or set(train) & set(test) or set(train) | set(test) != set(all_ids):
            raise RuntimeError("Holdout split leakage or coverage failure")
    return splits


def quadratic_features(x):
    x = np.asarray(x, dtype=float)
    return np.column_stack([np.ones(len(x)), x[:, 0], x[:, 1], x[:, 0] ** 2,
                            x[:, 0] * x[:, 1], x[:, 1] ** 2])


def fit_predict_models(train_x, train_y, test_x):
    """No held-out target argument exists; every learned transform uses training rows."""
    tx, ty, vx = map(lambda a: np.asarray(a, dtype=float), [train_x, train_y, test_x])
    if tx.ndim != 2 or tx.shape[1] != 2 or ty.shape != (len(tx), 2) or vx.ndim != 2 or vx.shape[1] != 2:
        raise ValueError("Invalid training/test shapes")
    if not all(np.isfinite(a).all() for a in [tx, ty, vx]):
        raise ValueError("Nonfinite model inputs")
    lower, upper = tx.min(axis=0), tx.max(axis=0)
    scale = upper - lower
    if np.any(scale <= 0):
        raise ValueError("Training inputs must vary on both axes")
    xn, vn = (tx - lower) / scale, (vx - lower) / scale
    means, stds, fitted_y_means, fitted_y_stds = [], [], [], []
    for objective in range(2):
        kernel = ConstantKernel(1.0, constant_value_bounds="fixed") * Matern(
            length_scale=[0.3, 0.3], length_scale_bounds="fixed", nu=2.5)
        gp = GaussianProcessRegressor(kernel=kernel, alpha=1e-8, optimizer=None, normalize_y=True)
        gp.fit(xn, ty[:, objective])
        mu, sd = gp.predict(vn, return_std=True)
        means.append(mu)
        stds.append(sd)
        fitted_y_means.append(float(gp._y_train_mean))
        fitted_y_stds.append(float(gp._y_train_std))
    features = quadratic_features(xn)
    coefficients, _, rank, singular = np.linalg.lstsq(features, ty, rcond=None)
    predictions = {
        "fixed_matern_gp": {"mean": np.column_stack(means), "latent_sd": np.column_stack(stds)},
        "training_mean": {"mean": np.broadcast_to(ty.mean(axis=0), (len(vx), 2)).copy(), "latent_sd": None},
        "quadratic_regression": {"mean": quadratic_features(vn) @ coefficients, "latent_sd": None},
    }
    metadata = {"input_min_from_training": lower.tolist(), "input_max_from_training": upper.tolist(),
                "input_scale_from_training": scale.tolist(), "GP_outcome_mean_from_training": fitted_y_means,
                "GP_outcome_scale_from_training": fitted_y_stds,
                "quadratic_coefficients_training_only": coefficients.tolist(), "quadratic_design_rank": int(rank),
                "quadratic_design_condition_number": float(singular[0] / singular[-1]),
                "heldout_targets_seen_by_fit": False,
                "test_outside_training_axis_box": np.any((vx < lower) | (vx > upper), axis=1).tolist()}
    return predictions, metadata


def prediction_statistics(truth, mean, latent_sd=None):
    truth, mean = np.asarray(truth, dtype=float), np.asarray(mean, dtype=float)
    if truth.ndim != 1 or truth.shape != mean.shape or not len(truth) or not np.isfinite(truth).all() or not np.isfinite(mean).all():
        raise ValueError("Expected equal nonempty finite one-dimensional targets")
    error = mean - truth
    stats = {"n_test": len(truth), "MAE": float(np.mean(np.abs(error))),
             "RMSE": float(np.sqrt(np.mean(error * error))), "mean_prediction_minus_truth": float(error.mean()),
             "latent_interval_coverage_95pct": None, "mean_latent_sd": None,
             "standardized_residual_mean_truth_minus_prediction": None,
             "standardized_residual_RMS": None, "max_abs_standardized_residual": None,
             "zero_latent_sd_count": None}
    if latent_sd is not None:
        sd = np.asarray(latent_sd, dtype=float)
        if sd.shape != truth.shape or not np.isfinite(sd).all() or np.any(sd < 0):
            raise ValueError("Invalid latent standard deviations")
        nonzero = sd > 0
        z = (truth[nonzero] - mean[nonzero]) / sd[nonzero]
        stats.update({"latent_interval_coverage_95pct": float(np.mean(np.abs(error) <= 1.96 * sd)),
                      "mean_latent_sd": float(sd.mean()),
                      "standardized_residual_mean_truth_minus_prediction": float(z.mean()) if len(z) else None,
                      "standardized_residual_RMS": float(np.sqrt(np.mean(z * z))) if len(z) else None,
                      "max_abs_standardized_residual": float(np.max(np.abs(z))) if len(z) else None,
                      "zero_latent_sd_count": int(np.sum(~nonzero))})
    return stats


def run_holdouts(x, y, random_count=20):
    splits = make_holdout_splits(x, random_count)
    prediction_rows, statistics_rows, split_metadata = [], [], []
    for split in splits:
        train_ids, test_ids = np.array(split["train_ids"]), np.array(split["test_ids"])
        predictions, meta = fit_predict_models(x[train_ids], y[train_ids], x[test_ids])
        split_metadata.append({**split, **meta})
        for model in MODELS:
            value = predictions[model]
            for target, name in enumerate(TARGETS):
                mu = value["mean"][:, target]
                sd = value["latent_sd"][:, target] if value["latent_sd"] is not None else None
                statistics_rows.append({"split_id": split["split_id"], "split_kind": split["kind"],
                    "model": model, "target": name, **prediction_statistics(y[test_ids, target], mu, sd)})
                for position, index in enumerate(test_ids):
                    sigma = float(sd[position]) if sd is not None else None
                    truth = float(y[index, target])
                    mean = float(mu[position])
                    prediction_rows.append({"split_id": split["split_id"], "split_kind": split["kind"],
                        "model": model, "target": name, "candidate_index": int(index),
                        "flow_rate_uL_min": float(x[index, 0]), "overpotential_V": float(x[index, 1]),
                        "truth_cached_PDE_objective": truth, "prediction": mean,
                        "prediction_minus_truth": mean - truth, "GP_latent_sd": sigma,
                        "GP_latent_95pct_lower": mean - 1.96 * sigma if sigma is not None else None,
                        "GP_latent_95pct_upper": mean + 1.96 * sigma if sigma is not None else None,
                        "standardized_residual_truth_minus_prediction": (truth - mean) / sigma if sigma and sigma > 0 else None,
                        "GP_latent_interval_covers_cached_truth": abs(truth - mean) <= 1.96 * sigma if sigma is not None else None,
                        "outside_training_axis_box": meta["test_outside_training_axis_box"][position],
                        "evidence_role": "held_out_deterministic_PDE_prediction_not_physical_uncertainty"})
    grouped = []
    for kind in sorted({s["kind"] for s in splits}):
        for model in MODELS:
            for target in TARGETS:
                group = [r for r in prediction_rows if r["split_kind"] == kind and r["model"] == model and r["target"] == target]
                sd = [r["GP_latent_sd"] for r in group] if model == "fixed_matern_gp" else None
                grouped.append({"split_kind": kind, "model": model, "target": target,
                    "split_count": len({r["split_id"] for r in group}),
                    "aggregation": "descriptive_pooled_heldout_uses_overlapping_splits_not_independent_new_data",
                    **prediction_statistics([r["truth_cached_PDE_objective"] for r in group], [r["prediction"] for r in group], sd)})
    return {"predictions": prediction_rows, "statistics": statistics_rows, "metadata": split_metadata,
            "grouped_statistics": grouped}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pilot", action="store_true", help="Bounded 2-seed, 10-budget, 8-split check in pilot/; no full-run overwrite")
    args = parser.parse_args()
    destination = OUT / "pilot" if args.pilot else OUT
    destination.mkdir(parents=True, exist_ok=True)
    seed_count, budget, draws, random_count = (2, 10, 256, 2) if args.pilot else (64, 25, 256, 20)
    started = datetime.now(timezone.utc).isoformat()
    start = time.perf_counter()
    x, y = load_pool()
    source_before = {"pool": sha(CONTROL / "candidate_pool.csv"), "old_control": sha(Path(__file__).parent / "metrics_control.py"),
                     "old_campaign": sha(CONTROL / "sequential_campaigns.csv"), "extension": sha(Path(__file__))}
    campaigns = run_campaign_extension(x, y, seed_count=seed_count, budget=budget, mc_draws=draws)
    holdouts = run_holdouts(x, y, random_count=random_count)
    if source_before["pool"] != sha(CONTROL / "candidate_pool.csv") or source_before["old_control"] != sha(Path(__file__).parent / "metrics_control.py"):
        raise RuntimeError("Frozen input changed during run")
    if not args.pilot and not campaigns["saved_prefix_replication"]["all_original_240_records_reproduced"]:
        raise RuntimeError("Original early-budget prefixes were not reproduced")
    write_csv(destination / "campaign_evaluations.csv", campaigns["records"])
    write_csv(destination / "budget_prefixes.csv", campaigns["prefixes"])
    write_csv(destination / "learning_curve_summary.csv", campaigns["learning_summary"])
    write_csv(destination / "paired_HV_differences.csv", campaigns["paired_summary"])
    write_json(destination / "oracle_budget_checks.json", campaigns["budget_checks"])
    write_csv(destination / "holdout_predictions.csv", holdouts["predictions"])
    write_csv(destination / "holdout_split_metrics.csv", holdouts["statistics"])
    write_csv(destination / "holdout_grouped_metrics.csv", holdouts["grouped_statistics"])
    write_json(destination / "holdout_split_metadata.json", holdouts["metadata"])
    summary = {"schema_version": "1.0", "mode": "bounded_pilot" if args.pilot else "full_64_seed_extension",
        "started_UTC": started, "completed_UTC": datetime.now(timezone.utc).isoformat(),
        "wall_seconds": time.perf_counter() - start,
        "evidence_role": "executed_cached_numerical_benchmark_not_new_PDE_or_experiment",
        "new_PDE_solves": 0, "new_experimental_runs": 0, "hardware_connections": 0,
        "pool_size": len(x), "pool_grid": [9, 9], "frozen_pool_PDE_grid": [80, 32],
        "campaign": {"seeds": list(range(seed_count)), "methods": METHODS,
            "budget_including_initial": budget, "shared_initial_count": 5, "MC_draws_per_candidate": draws,
            "total_runs": seed_count * len(METHODS), "cached_oracle_evaluation_uses": len(campaigns["records"]),
            "budget_prefixes": sorted({r["evaluation"] for r in campaigns["prefixes"]}),
            "objectives_maximized": TARGETS, "reference": REFERENCE, "full_pool_HV": campaigns["full_pool_HV"],
            "all_budget_and_uniqueness_checks_passed": all(c["passed"] for c in campaigns["budget_checks"]),
            "all_methods_have_shared_initial_per_seed": True,
            "outcomes_available_to_maximin_selection": False,
            "maximin_ties": "lowest remaining candidate index",
            "saved_original_prefix_replication": campaigns["saved_prefix_replication"]},
        "holdout": {"split_count": len(holdouts["metadata"]), "random_seeds": list(range(2000, 2000 + random_count)),
            "flow_blocks": 3, "eta_blocks": 3, "train_count_per_split": 54, "test_count_per_split": 27,
            "models": MODELS, "targets": TARGETS, "prediction_rows": len(holdouts["predictions"]),
            "split_metric_rows": len(holdouts["statistics"]),
            "preprocessing": "input min/max and GP outcome normalization fit using training rows only",
            "GP": "independent normalized Matern 5/2, fixed length scales [0.3,0.3], amplitude 1, alpha 1e-8, optimizer None",
            "GP_interval": "mean +/- 1.96 latent posterior sd; numerical interpolation diagnostic, not physical uncertainty",
            "polynomial": "OLS training-only [1,x1,x2,x1^2,x1*x2,x2^2] on train-normalized coordinates",
            "test_set_hyperparameter_tuning": False},
        "limitations": ["Entire domain is the frozen deterministic 81-point simulation pool with assumed product mass and voltage law",
            "No additional information about transport-model discrepancy, chemistry, instrumentation or industrial performance",
            "Paired t intervals describe seed sensitivity under approximate independent normal seed-effect assumptions",
            "No p-value, universal method superiority or industrial optimum claimed",
            "Learning-curve prefix budgets and paired comparisons are correlated descriptive views",
            "Overlapping holdout uses are not independent new observations",
            "Blocked holdout extrapolates only beyond the training subset; all points remain inside the declared pool domain",
            "GP latent sd can be miscalibrated and is not PDE discretization or physical measurement uncertainty"],
        "source_sha256": source_before,
        "runtime_versions": {"python": platform.python_version(), "numpy": np.__version__, "scipy": version("scipy"),
                             "scikit_learn": version("scikit-learn")},
        "thread_environment": {name: os.environ.get(name, "unspecified") for name in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"]},
        "output_sha256": {p.name: sha(p) for p in sorted(destination.iterdir()) if p.is_file() and p.name != "summary.json"}}
    write_json(destination / "summary.json", summary)
    print(json.dumps({"mode": summary["mode"], "wall_seconds": summary["wall_seconds"],
        "campaign_evaluation_uses": len(campaigns["records"]), "holdout_splits": len(holdouts["metadata"]),
        "saved_prefix_replication": campaigns["saved_prefix_replication"],
        "last_budget_learning_curve": [r for r in campaigns["learning_summary"] if r["budget"] == budget],
        "last_budget_paired_comparisons": [r for r in campaigns["paired_summary"] if r["budget"] == budget]}, indent=2))


if __name__ == "__main__":
    main()
