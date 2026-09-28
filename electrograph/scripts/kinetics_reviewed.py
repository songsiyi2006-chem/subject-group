"""Exact count-based SSA and CTMC checks for an ASSUMED independent-site cycle.

No spatial interactions, periodic boundary operation, fitted electrochemical
parameters, or experimental validation are present. Counts are a sufficient
state for independent sites of each fixed rate class, so this aggregation does
not approximate the corresponding site-resolved Markov process.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform
import random
import time
from pathlib import Path

import numpy as np
import scipy
from scipy.linalg import expm
from scipy.stats import t as student_t

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parent
STATES = ("empty", "substrate", "radical", "intermediate", "product")
REACTIONS = ("adsorption", "SET", "coupling", "PCET", "desorption")
ELECTRON_INVENTORY = (0, 0, 1, 1, 2)
ELEMENTARY_CHARGE_C = 1.602176634e-19
SOURCE_URLS = [
    "https://pubs.acs.org/doi/10.1021/j100540a008",
    "https://www.cl.cam.ac.uk/teaching/2526/Bioinfo/papers/DanielGillespie1.pdf",
    "https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.expm.html",
]


def rate_constants(eta_V: float = 0.48) -> list[float]:
    """Source's assumed rates in s^-1, retaining its F, R, and T constants."""
    if not math.isfinite(eta_V) or not -1 <= eta_V <= 1:
        raise ValueError("eta_V must be finite and within the supported [-1, 1] V range")
    factor = math.exp(0.5 * 96485.33 * eta_V / (8.314 * 298.15))
    return [45.0, 120.0 * factor, 350.0, 200.0 * factor, 80.0]


def _rates(values):
    a = np.asarray(values, dtype=float)
    if a.shape != (5,) or not np.all(np.isfinite(a)) or np.any(a < 0):
        raise ValueError("exactly five finite nonnegative first-order rates are required")
    return a


def generator_matrix(rates) -> np.ndarray:
    """Row generator Q; row probabilities obey dp/dt = p Q."""
    rates = _rates(rates)
    q = np.diag(-rates)
    for i in range(5):
        q[i, (i + 1) % 5] = rates[i]
    return q


def stationary_solution(rates) -> dict:
    rates = _rates(rates)
    if np.any(rates == 0):
        raise ValueError("a unique cyclic steady state requires strictly positive rates")
    residence = 1.0 / rates
    mean = float(residence.sum())
    return {
        "occupancy": (residence / mean).tolist(),
        "TOF_s_1": 1.0 / mean,
        "mean_cycle_wait_s": mean,
        "cycle_wait_variance_s2": float(np.sum(residence**2)),
        "rate_log_elasticity": (residence / mean).tolist(),
        "renewal_long_time_product_Fano_factor": float(np.sum(residence**2) / mean**2),
    }


def ctmc_window(rates, horizon_s=1.0, burn_in_s=0.2, initial_probability=None) -> dict:
    """Finite-window mean from an augmented matrix exponential, without inversion."""
    rates = _rates(rates)
    if not math.isfinite(horizon_s) or not math.isfinite(burn_in_s) or not 0 <= burn_in_s < horizon_s:
        raise ValueError("require 0 <= burn_in_s < finite horizon_s")
    p0 = np.array([1., 0, 0, 0, 0] if initial_probability is None else initial_probability, dtype=float)
    if p0.shape != (5,) or not np.all(np.isfinite(p0)) or np.any(p0 < 0) or not np.isclose(p0.sum(), 1):
        raise ValueError("initial_probability must be a five-state probability vector")
    block = np.zeros((10, 10))
    block[:5, :5] = generator_matrix(rates).T
    block[5:, :5] = np.eye(5)
    initial = np.r_[p0, np.zeros(5)]
    stop = expm(block * horizon_s) @ initial
    start = expm(block * burn_in_s) @ initial
    integrated = stop[5:] - start[5:]
    duration = horizon_s - burn_in_s
    extents = rates * integrated
    return {
        "final_probability": stop[:5].tolist(),
        "burn_in_probability": start[:5].tolist(),
        "time_average_probability": (integrated / duration).tolist(),
        "expected_reaction_counts_per_site": extents.tolist(),
        "expected_TOF_s_1": float(extents[4] / duration),
        "expected_electrons_per_site_s": float((extents[1] + extents[3]) / duration),
        "probability_sum_error": float(abs(stop[:5].sum() - 1)),
        "window_occupancy_sum_error": float(abs(integrated.sum() / duration - 1)),
    }


def cycle_wait_cdf(rates, time_s):
    """CDF for one complete cycle: sum of five independent exponential waits."""
    rates = _rates(rates)
    if np.any(rates == 0) or not math.isfinite(time_s) or time_s < 0:
        raise ValueError("positive rates and nonnegative finite time required")
    transient = np.diag(-rates)
    for i in range(4):
        transient[i, i + 1] = rates[i]
    return float(1 - (expm(transient * time_s)[0]).sum())


def simulate(rates=None, n_sites=256, horizon_s=1.0, burn_in_s=0.2, seed=20260928,
             *, rate_classes=None, initial_counts=None, record_every=0,
             max_events=None):
    """Direct SSA, censoring the first event beyond horizon_s.

    rate_classes is [{"sites": integer, "rates": five rates}, ...]. Sites in a
    class share fixed rates; classes never exchange sites. max_events is a
    diagnostic stopping rule ONLY: if reached, results use the actual shorter
    horizon, explicitly flag event_limit, and are not a fixed-time estimate.
    """
    if not math.isfinite(horizon_s) or not math.isfinite(burn_in_s) or not 0 <= burn_in_s < horizon_s:
        raise ValueError("require 0 <= burn_in_s < finite horizon_s")
    if not isinstance(seed, int) or isinstance(seed, bool):
        raise ValueError("integer seed required")
    if not isinstance(record_every, int) or record_every < 0:
        raise ValueError("record_every must be a nonnegative integer")
    if max_events is not None and (not isinstance(max_events, int) or max_events < 1):
        raise ValueError("max_events must be a positive integer or None")
    if rate_classes is None:
        rate_classes = [{"sites": n_sites, "rates": rate_constants() if rates is None else rates}]
    if not rate_classes:
        raise ValueError("at least one rate class is required")
    class_sizes, flat_rates = [], []
    for cls in rate_classes:
        size = cls["sites"]
        if not isinstance(size, int) or isinstance(size, bool) or size < 1:
            raise ValueError("each class must contain a positive integer number of sites")
        class_sizes.append(size)
        flat_rates.extend(_rates(cls["rates"]).tolist())
    n_sites = sum(class_sizes)
    counts = [x for size in class_sizes for x in [size, 0, 0, 0, 0]]
    if initial_counts is not None:
        initial = np.asarray(initial_counts)
        if initial.shape != (len(class_sizes), 5) or not np.all(np.isfinite(initial)):
            raise ValueError("initial_counts must have shape (classes, 5)")
        if np.any(initial < 0) or not np.all(initial == np.floor(initial)) or not np.array_equal(initial.sum(axis=1), class_sizes):
            raise ValueError("initial_counts must be nonnegative integers preserving each class size")
        counts = initial.astype(int).ravel().tolist()
    initial = counts.copy()
    burn_counts = counts.copy() if burn_in_s == 0 else None
    channel_count = len(counts)
    next_channel = [5 * (i // 5) + (i + 1) % 5 for i in range(channel_count)]
    integrals = [0.] * channel_count
    events, observed_events = [0] * channel_count, [0] * channel_count
    propensities = [count * rate for count, rate in zip(counts, flat_rates)]
    rng = random.Random(seed)
    draw = rng.random
    trajectory = []

    def snapshot(t, event_channel=-1):
        row = {"time_s": t, "event_channel": event_channel,
               "cumulative_events": sum(events), "cumulative_products": sum(events[4::5])}
        row.update({f"n_{state}": sum(counts[i::5]) for i, state in enumerate(STATES)})
        row.update({f"events_{reaction}": sum(events[i::5]) for i, reaction in enumerate(REACTIONS)})
        trajectory.append(row)

    if record_every:
        snapshot(0.)
    now, steps, censored, stop_reason = 0., 0, False, "horizon"
    while now < horizon_s:
        total = sum(propensities)
        scheduled = now - math.log1p(-draw()) / total if total > 0 else math.inf
        end = min(scheduled, horizon_s)
        if burn_counts is None and end >= burn_in_s:
            burn_counts = counts.copy()
        observed_dt = max(0., end - max(now, burn_in_s))
        if observed_dt:
            for index in range(channel_count):
                integrals[index] += counts[index] * observed_dt
        now = end
        if scheduled > horizon_s or total == 0:
            censored = total > 0
            stop_reason = "horizon" if total > 0 else "absorbing"
            break
        threshold = draw() * total
        accumulated, chosen = 0., channel_count - 1
        for index, propensity in enumerate(propensities):
            accumulated += propensity
            if threshold < accumulated:
                chosen = index
                break
        successor = next_channel[chosen]
        counts[chosen] -= 1
        counts[successor] += 1
        events[chosen] += 1
        if now > burn_in_s:
            observed_events[chosen] += 1
        elif now == burn_in_s:
            burn_counts = counts.copy()
        propensities[chosen] = counts[chosen] * flat_rates[chosen]
        propensities[successor] = counts[successor] * flat_rates[successor]
        steps += 1
        if record_every and steps % record_every == 0:
            snapshot(now, chosen)
        if max_events is not None and steps >= max_events:
            stop_reason = "event_limit"
            break
    if now <= burn_in_s or burn_counts is None:
        raise ValueError("event cap stopped before observation window; lower burn-in or increase cap")
    if record_every and (not trajectory or trajectory[-1]["time_s"] != now):
        snapshot(now)
    duration = now - burn_in_s
    aggregate_counts = lambda values: [sum(values[i::5]) for i in range(5)]
    final_pop, burn_pop = aggregate_counts(counts), aggregate_counts(burn_counts)
    event_totals, obs = aggregate_counts(events), aggregate_counts(observed_events)
    coverage = [x / (n_sites * duration) for x in aggregate_counts(integrals)]
    state_balance = [final_pop[i] - burn_pop[i] - obs[(i - 1) % 5] + obs[i] for i in range(5)]
    electrons = obs[1] + obs[3]
    charge_inventory_change = sum((final_pop[i] - burn_pop[i]) * ELECTRON_INVENTORY[i] for i in range(5))
    return {
        "seed": seed, "sites": n_sites, "rate_classes": rate_classes,
        "requested_horizon_s": horizon_s, "simulated_time_s": now,
        "burn_in_s": burn_in_s, "observation_duration_s": duration,
        "stop_reason": stop_reason, "censored_final_event": censored,
        "total_events": steps, "reaction_counts_full": event_totals,
        "reaction_counts_observed": obs, "products_observed": obs[4],
        "electrons_observed": electrons,
        "TOF_s_1": obs[4] / (n_sites * duration),
        "electrons_per_site_s": electrons / (n_sites * duration),
        "current_for_explicit_sites_A": electrons * ELEMENTARY_CHARGE_C / duration,
        "current_scope": "literal modeled site count; not macroscopic electrode current",
        "initial_counts_by_class": np.array(initial).reshape(-1, 5).tolist(),
        "final_counts_by_class": np.array(counts).reshape(-1, 5).tolist(),
        "burn_in_counts_by_class": np.array(burn_counts).reshape(-1, 5).tolist(),
        "reaction_counts_observed_by_class": np.array(observed_events).reshape(-1, 5).tolist(),
        "final_counts": final_pop, "burn_in_counts": burn_pop,
        "time_average_coverage": coverage,
        "time_integrated_counts_by_class": np.array(integrals).reshape(-1, 5).tolist(),
        "state_balance_residuals": state_balance,
        "electron_product_inventory_residual": electrons - 2 * obs[4] - charge_inventory_change,
        "occupancy_sum_error": abs(sum(coverage) - 1),
        "trajectory": trajectory,
    }


def _write_csv(path, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _record(result, family, eta, case):
    row = {"family": family, "case": case, "eta_V": eta}
    for name in ("seed", "sites", "requested_horizon_s", "simulated_time_s", "burn_in_s",
                 "observation_duration_s", "stop_reason", "censored_final_event", "total_events",
                 "products_observed", "electrons_observed", "TOF_s_1", "electrons_per_site_s",
                 "current_for_explicit_sites_A", "electron_product_inventory_residual", "occupancy_sum_error"):
        row[name] = result[name]
    for name in ("final_counts_by_class", "burn_in_counts_by_class", "reaction_counts_observed_by_class", "time_integrated_counts_by_class"):
        row[name] = json.dumps(result[name], separators=(",", ":"))
    for i, state in enumerate(STATES):
        row[f"final_{state}"] = result["final_counts"][i]
        row[f"burn_in_{state}"] = result["burn_in_counts"][i]
        row[f"mean_{state}"] = result["time_average_coverage"][i]
        row[f"state_balance_{state}"] = result["state_balance_residuals"][i]
        row[f"events_{REACTIONS[i]}"] = result["reaction_counts_observed"][i]
    return row


def _group(rows, analytic_tof, steady_tof):
    values = np.array([r["TOF_s_1"] for r in rows])
    n = len(values)
    sd = float(values.std(ddof=1))
    half = float(student_t.ppf(.975, n - 1) * sd / math.sqrt(n))
    return {"trajectories": n, "mean_TOF_s_1": float(values.mean()), "sd_TOF_s_1": sd,
            "mean_MC_lower_95pct_t": float(values.mean() - half),
            "mean_MC_upper_95pct_t": float(values.mean() + half),
            "finite_window_CTMC_TOF_s_1": analytic_tof,
            "steady_TOF_s_1": steady_tof,
            "relative_mean_error_vs_finite_CTMC": float((values.mean() - analytic_tof) / analytic_tof),
            "total_events": sum(r["total_events"] for r in rows)}


def run_study(output_dir=None, pilot=False):
    out = Path(output_dir) if output_dir else ROOT / "results" / "kinetics" / ("pilot" if pilot else "")
    out.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    rows, scans, scale, hetero, legacy = [], [], [], [], []
    eta_values = [.2, .3, .4, .5, .6, .7] if not pilot else [.2, .7]
    n_seeds = 32 if not pilot else 2
    for eta_index, eta in enumerate(eta_values):
        rates = rate_constants(eta)
        finite, steady = ctmc_window(rates), stationary_solution(rates)
        group = []
        for offset in range(n_seeds):
            r = simulate(rates, seed=20260928 + 1000 * eta_index + offset)
            row = _record(r, "potential_scan", eta, f"eta_{eta:.1f}")
            rows.append(row)
            group.append(row)
        scans.append({"eta_V": eta, **_group(group, finite["expected_TOF_s_1"], steady["TOF_s_1"]),
                      "CTMC_probability_sum_error": finite["probability_sum_error"],
                      "CTMC_window_occupancy_sum_error": finite["window_occupancy_sum_error"],
                      **{f"steady_{s}": steady["occupancy"][i] for i, s in enumerate(STATES)}})
        print(f"potential eta={eta:.1f} complete: {n_seeds} trajectories", flush=True)
    baseline = simulate(rate_constants(.48), seed=7101, record_every=1)
    rows.append(_record(baseline, "baseline", .48, "all_event_trajectory"))
    _write_csv(out / "baseline_trajectory.csv", baseline.pop("trajectory"))
    (out / "baseline.json").write_text(json.dumps(baseline, indent=2) + "\n", encoding="utf-8")
    finite48, steady48 = ctmc_window(rate_constants(.48)), stationary_solution(rate_constants(.48))
    if not pilot:
        for sites in [64, 256, 1024]:
            group = []
            for offset in range(16):
                r = simulate(rate_constants(.48), n_sites=sites, seed=50000 + sites + offset)
                row = _record(r, "site_scaling", .48, f"sites_{sites}")
                rows.append(row)
                group.append(row)
            scale.append({"sites": sites, **_group(group, finite48["expected_TOF_s_1"], steady48["TOF_s_1"]),
                          "asymptotic_TOF_sd_renewal": math.sqrt(steady48["renewal_long_time_product_Fano_factor"] * steady48["TOF_s_1"] / (sites * .8))})
        rates_a = rate_constants(.48)
        rates_b = rates_a.copy()
        rates_b[0] *= .25
        rates_b[4] *= .5
        classes = [{"sites": 192, "rates": rates_a}, {"sites": 64, "rates": rates_b}]
        group = []
        for offset in range(16):
            r = simulate(rate_classes=classes, seed=61000 + offset)
            row = _record(r, "two_rate_classes", .48, "75pct_reference_25pct_slow")
            rows.append(row)
            group.append(row)
        finite_b, steady_b = ctmc_window(rates_b), stationary_solution(rates_b)
        mixture_finite = .75 * finite48["expected_TOF_s_1"] + .25 * finite_b["expected_TOF_s_1"]
        mixture_steady = .75 * steady48["TOF_s_1"] + .25 * steady_b["TOF_s_1"]
        hetero.append({"case": "75pct_reference_25pct_slow", **_group(group, mixture_finite, mixture_steady),
                       "classes": classes, "spatial_interactions": False})
        for offset in range(16):
            r = simulate(rate_constants(.48), burn_in_s=0., seed=72000 + offset, max_events=3500)
            row = _record(r, "fixed_event_transient", .48, "3500_events_all_empty_start")
            rows.append(row)
            # A CTMC at each realized stopping time is descriptive only: the stopping time
            # depends on the path, so it is NOT an unbiased reference for this ensemble.
            legacy.append({"seed": r["seed"], "simulated_time_s": r["simulated_time_s"],
                           "TOF_s_1": r["TOF_s_1"], "total_events": r["total_events"],
                           "steady_TOF_s_1": steady48["TOF_s_1"],
                           "interpretation": "path-dependent transient stopping; no steady-state guarantee"})
    analytic_rows = []
    for eta in [.2, .3, .4, .48, .5, .6, .7]:
        rates = rate_constants(eta)
        steady = stationary_solution(rates)
        for index, step in enumerate(REACTIONS):
            for factor in [.5, 1., 2.]:
                modified = rates.copy()
                modified[index] *= factor
                result = stationary_solution(modified)
                analytic_rows.append({"eta_V": eta, "step": step, "rate_multiplier": factor,
                                      "baseline_rate_s_1": rates[index], "TOF_s_1": result["TOF_s_1"],
                                      "baseline_TOF_s_1": steady["TOF_s_1"],
                                      "baseline_log_elasticity": steady["rate_log_elasticity"][index]})
    transient_rows = []
    for t in np.linspace(0, .2, 101):
        p = expm(generator_matrix(rate_constants(.48)).T * t) @ np.array([1., 0, 0, 0, 0])
        transient_rows.append({"time_s": float(t), **{s: float(p[i]) for i, s in enumerate(STATES)}})
    wait_rows = [{"time_s": float(t), "cycle_wait_CDF": cycle_wait_cdf(rate_constants(.48), t)} for t in np.linspace(0, .3, 151)]
    files = {"trajectories.csv": rows, "potential_summary.csv": scans,
             "analytic_rate_perturbations.csv": analytic_rows,
             "analytic_transient.csv": transient_rows, "cycle_wait_distribution.csv": wait_rows}
    if not pilot:
        files.update({"site_scaling.csv": scale, "fixed_event_transients.csv": legacy})
    for name, data in files.items():
        _write_csv(out / name, data)
    sources = [Path(__file__), ROOT / "source" / "electrograph_kmc_core.py", REPO / "tests" / "test_electrograph_kinetics.py"]
    summary = {
        "schema_version": 1, "mode": "pilot" if pilot else "full",
        "evidence": "ASSUMED independent-site stochastic model; numerical verification only",
        "source_corrections": ["post-event coverage logging", "actual executed event count", "fixed-time censoring",
                               "time-weighted occupancy", "computed potential sweep replaces hard-coded TOFs"],
        "model": {"states": STATES, "steps": REACTIONS, "reference_eta_V": .48,
                  "reference_rates_s_1": rate_constants(.48), "horizon_s": 1., "burn_in_s": .2,
                  "electrons_per_complete_cycle": 2, "inter_site_interactions": False,
                  "spatial_periodic_boundary_operations": False, "fitted_rates": False},
        "counts": {"study_SSA_trajectories": len(rows), "potential_trajectories": len(eta_values) * n_seeds,
                   "baseline_trajectories": 1, "site_scaling_trajectories": 0 if pilot else 48,
                   "two_class_trajectories": 0 if pilot else 16,
                   "fixed_event_diagnostic_trajectories": 0 if pilot else 16,
                   "executed_SSA_events": sum(row["total_events"] for row in rows),
                   "analytic_rate_perturbation_cases": len(analytic_rows),
                   "analytic_transient_times": len(transient_rows), "analytic_cycle_CDF_times": len(wait_rows),
                   "pilot_excluded_from_full_counts": True, "unit_tests_excluded_from_study_counts": True},
        "potential_summary": scans, "site_scaling": scale, "two_class_summary": hetero,
        "reference_steady": steady48, "reference_finite_window": finite48,
        "fixed_event_diagnostic": {"count": len(legacy), "mean_TOF_s_1": float(np.mean([r["TOF_s_1"] for r in legacy])) if legacy else None,
                                   "min_time_s": min((r["simulated_time_s"] for r in legacy), default=None),
                                   "max_time_s": max((r["simulated_time_s"] for r in legacy), default=None)},
        "rate_limited_TOF_upper_bound_s_1": 1. / (1. / 45 + 1. / 350 + 1. / 80),
        "source_hardcoded_TOFs_not_used": [4.2, 12.8, 38.4, 95.1, 184.6, 260.2],
        "numerical_checks": {"max_occupancy_sum_error": max(r["occupancy_sum_error"] for r in rows),
                             "all_site_and_electron_balances_exact": all(r["electron_product_inventory_residual"] == 0 and all(r[f"state_balance_{s}"] == 0 for s in STATES) for r in rows),
                             "all_fixed_time_horizons_exact": all(r["simulated_time_s"] == 1. for r in rows if r["family"] != "fixed_event_transient"),
                             "reference_CTMC_probability_sum_error": finite48["probability_sum_error"],
                             "reference_CTMC_window_sum_error": finite48["window_occupancy_sum_error"]},
        "interval_scope": "Student-t intervals for mean across independent PRNG seeds conditional on fixed assumed rates; not parameter/physical confidence; simultaneous coverage not adjusted",
        "limitations": ["Five irreversible pseudo-first-order reactions; no thermodynamic reverse rates or electrolyte transport",
                        "Two classes represent assigned static rate variation only; no neighbors, hopping, geometry or periodic operations",
                        "No experimental or DFT rate calibration; labels do not establish a POP-SAC mechanism",
                        "High-potential saturation follows the assumed cycle; not a measured electrocatalytic trend",
                        "Fixed event stopping is path dependent and its apparent TOF is a transient diagnostic",
                        "Current refers only to the explicitly enumerated sites and is not electrode-scale current"],
        "runtime_seconds": time.perf_counter() - started,
        "runtime_versions": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
        "source_sha256": {p.relative_to(REPO).as_posix(): _sha(p) for p in sources if p.exists()},
        "output_sha256": {p.name: _sha(p) for p in sorted(out.glob("*.csv"))} | {"baseline.json": _sha(out / "baseline.json")},
        "sources": SOURCE_URLS,
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"counts": summary["counts"], "numerical_checks": summary["numerical_checks"],
                      "runtime_seconds": summary["runtime_seconds"]}, indent=2), flush=True)
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pilot", action="store_true")
    parser.add_argument("--output-dir", type=Path)
    arguments = parser.parse_args()
    run_study(arguments.output_dir, arguments.pilot)
