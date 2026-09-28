"""Scenario propagation and un-clipped Jansen Sobol diagnostics for ElectraTwin.

The independent parameter distributions are declared assumptions, not calibrated
posteriors. Paired-row bootstrap intervals diagnose resampling sensitivity of a
single scrambled QMC design; they are not rigorous QMC confidence guarantees.
No physical equipment, network connection or experimental observation is used.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import platform
import time

import numpy as np
import scipy
from scipy.stats import qmc

try:
    from .transport_reviewed import solve_transport
except ImportError:
    from transport_reviewed import solve_transport

BASE = Path(__file__).resolve().parents[1]
PARAMETERS = (
    {"name": "flow_rate_uL_min", "low": 405., "high": 495., "distribution": "uniform", "unit": "uL/min"},
    {"name": "overpotential_V", "low": .46, "high": .50, "distribution": "uniform", "unit": "V; imposed model index"},
    {"name": "diffusivity_m2_s", "low": .77e-9, "high": 1.43e-9, "distribution": "uniform", "unit": "m2/s"},
    {"name": "current_scale_A_m2", "low": .04, "high": .16, "distribution": "log_uniform", "unit": "A/m2"},
    {"name": "height_m", "low": .00027, "high": .00033, "distribution": "uniform", "unit": "m"},
)
OUTPUTS = ("conversion_pct", "current_A", "STY_assumed_product_kg_m3_day", "SEC_assumed_product_kwh_kg")
SCENARIO_SEED = 20260928
BOOTSTRAP_SEED = 20260929
# Chosen in unit-quantile space before any output was examined.
GRID_UNIT_POINTS = np.array([
    [0, 0, 0, 0, 0], [1, 1, 1, 1, 1], [0, 1, 0, 1, 0], [1, 0, 1, 0, 1],
    [.5, .5, .5, .5, .5], [.25, .75, .25, .75, .25],
    [.75, .25, .75, .25, .75], [.5, .5, 0, 1, 1],
], dtype=float)


def transform_unit(unit):
    """Map unit quantiles to the five independent declared physical scenarios."""
    u = np.asarray(unit, dtype=float)
    if u.ndim < 1 or u.shape[-1] != len(PARAMETERS) or not np.isfinite(u).all() or np.any((u < 0) | (u > 1)):
        raise ValueError("Expected finite unit coordinates with final dimension five and values in [0,1]")
    x = np.empty_like(u)
    for i, parameter in enumerate(PARAMETERS):
        lo, hi = parameter["low"], parameter["high"]
        x[..., i] = lo * np.exp(u[..., i] * np.log(hi / lo)) if parameter["distribution"] == "log_uniform" else lo + u[..., i] * (hi - lo)
    return x


def make_design(power=8, seed=SCENARIO_SEED):
    if isinstance(power, bool) or int(power) != power or not 1 <= power <= 16:
        raise ValueError("Require an integer power from 1 to 16")
    engine = qmc.Sobol(d=10, scramble=True, rng=np.random.default_rng(seed), optimization=None)
    ab = engine.random_base2(m=int(power))
    a, b = ab[:, :5].copy(), ab[:, 5:].copy()
    hybrids = np.repeat(a[None, :, :], 5, axis=0)
    for i in range(5):
        hybrids[i, :, i] = b[:, i]
    return a, b, hybrids


def jansen_indices(y_a, y_b, y_ab):
    """Jansen estimators for AB_i = A with column i from B; no clipping."""
    a, b, ab = [np.asarray(v, dtype=float) for v in (y_a, y_b, y_ab)]
    if a.ndim == 1:
        a = a[:, None]
    if b.ndim == 1:
        b = b[:, None]
    if ab.ndim == 2:
        ab = ab[:, :, None]
    if a.ndim != 2 or b.shape != a.shape or ab.ndim != 3 or ab.shape[1:] != a.shape or len(a) < 2:
        raise ValueError("Expected YA/YB=(N,outputs), YAB=(parameters,N,outputs), N>=2")
    if not all(np.isfinite(v).all() for v in (a, b, ab)):
        raise ValueError("Jansen input must be finite")
    combined = np.concatenate([a, b], axis=0)
    variance = np.var(combined, axis=0, ddof=1)
    scale = np.max(np.abs(combined), axis=0)
    if np.any(variance <= np.finfo(float).eps * np.maximum(scale**2, np.finfo(float).tiny)):
        raise ValueError("Sensitivity indices are undefined for constant or numerically unresolved output variance")
    first = 1 - np.mean((b[None, :, :] - ab)**2, axis=1) / (2 * variance)
    total = np.mean((a[None, :, :] - ab)**2, axis=1) / (2 * variance)
    return {"first_order": first, "total_order": total, "variance": variance}


def paired_bootstrap(y_a, y_b, y_ab, replicates=500, seed=BOOTSTRAP_SEED):
    """Resample matched rows across A, B and every hybrid; QMC diagnostic only."""
    if isinstance(replicates, bool) or int(replicates) != replicates or replicates < 1:
        raise ValueError("Positive integer bootstrap count required")
    a, b, ab = [np.asarray(v, dtype=float) for v in (y_a, y_b, y_ab)]
    jansen_indices(a, b, ab)
    rng = np.random.default_rng(seed)
    results = []
    for _ in range(int(replicates)):
        indices = rng.integers(0, len(a), size=len(a))
        results.append(jansen_indices(a[indices], b[indices], ab[:, indices]))
    return {key: np.stack([item[key] for item in results]) for key in ["first_order", "total_order"]}


def metrics_from_summary(summary, eta):
    """Use assumed product mass and voltage, not a verified chemical identity."""
    if not summary["converged"] or not summary["nonnegative"]:
        raise ValueError("Transport result failed numerical feasibility")
    for key in ["linear_residual_relative", "material_balance_relative_error", "charge_balance_relative_error"]:
        if not math.isfinite(summary[key]) or summary[key] >= 1e-8:
            raise ValueError("Transport numerical balance/residual check failed")
    reaction = summary["reacted_mol_s"]
    current, volume = summary["current_A"], summary["reactor_volume_m3"]
    voltage = 1.85 + eta + 18.5 * current
    if min(reaction, current, volume, voltage) <= 0 or not all(math.isfinite(v) for v in [reaction, current, volume, voltage]):
        raise ValueError("Positive finite forward-product rate, current, volume and voltage required")
    mass_kg_s = reaction * .33143
    values = {
        "conversion_pct": summary["conversion_pct"], "current_A": current,
        "STY_assumed_product_kg_m3_day": mass_kg_s * 86400 / volume,
        "SEC_assumed_product_kwh_kg": current * voltage / (mass_kg_s * 3.6e6),
        "cell_voltage_assumed_V": voltage,
    }
    if not all(math.isfinite(values[key]) and values[key] > 0 for key in OUTPUTS):
        raise ValueError("Scenario outputs must be positive finite values")
    return values


def evaluate_point(unit, *, nx=60, ny=24):
    unit = np.asarray(unit, dtype=float)
    if unit.shape != (5,):
        raise ValueError("One five-coordinate point required")
    physical = transform_unit(unit)
    kwargs = {p["name"]: float(v) for p, v in zip(PARAMETERS, physical)}
    result = solve_transport(nx=nx, ny=ny, include_field=False, **kwargs)
    s = result["summary"]
    values = metrics_from_summary(s, kwargs["overpotential_V"])
    return {
        **{f"u_{p['name']}": float(v) for p, v in zip(PARAMETERS, unit)}, **kwargs,
        "nx": nx, "ny": ny, **values,
        "reacted_mol_s": s["reacted_mol_s"], "reactor_volume_m3": s["reactor_volume_m3"],
        "linear_residual_relative": s["linear_residual_relative"],
        "material_balance_relative_error": s["material_balance_relative_error"],
        "charge_balance_relative_error": s["charge_balance_relative_error"],
        "minimum_species_mM": s["minimum_species_mM"], "converged": s["converged"],
        "evidence_role": "executed_independent_parameter_scenario_not_calibrated_uncertainty",
    }


def local_design():
    points = [("center", None, None, None, np.full(5, .5))]
    for step in [.02, .01]:
        for i in range(5):
            for sign in [-1, 1]:
                u = np.full(5, .5)
                u[i] += sign * step
                points.append((f"step{step}_p{i}_{sign}", i, step, sign, u))
    return points


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save_json(path, data):
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8", newline="\n")


def save_csv(path, rows):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def run_analysis(output_dir=None, power=8, bootstrap_replicates=500):
    if not 6 <= power <= 9 or not 20 <= bootstrap_replicates <= 5000:
        raise ValueError("Bounded run requires power 6..9 and bootstrap count 20..5000")
    output = Path(output_dir) if output_dir else BASE / "results" / "uncertainty"
    output.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    executed_utc = datetime.now(timezone.utc).isoformat()
    a, b, hybrids = make_design(power)
    n = len(a)
    records, block_outputs = [], []
    for block, u in [("A", a), ("B", b)] + [(f"AB_{i}", hybrids[i]) for i in range(5)]:
        current = []
        for row_index, point in enumerate(u):
            row = {"evaluation_index": len(records), "block": block, "paired_row_index": row_index, **evaluate_point(point)}
            records.append(row)
            current.append([row[key] for key in OUTPUTS])
        block_outputs.append(np.asarray(current))
        print(f"Completed {block}: {len(records)} / {7*n} main transport evaluations", flush=True)
    save_csv(output / "design_evaluations.csv", records)
    ya, yb, yab = block_outputs[0], block_outputs[1], np.stack(block_outputs[2:])
    prefix_rows = []
    prefix_sizes = [v for v in [64, 128, 256, 512] if v <= n]
    for prefix in prefix_sizes:
        indices = jansen_indices(ya[:prefix], yb[:prefix], yab[:, :prefix])
        for i, p in enumerate(PARAMETERS):
            for j, metric in enumerate(OUTPUTS):
                prefix_rows.append({"base_N": prefix, "cached_PDE_evaluation_uses": prefix*7,
                    "parameter": p["name"], "output": metric,
                    "first_order_S": float(indices["first_order"][i,j]), "total_order_ST": float(indices["total_order"][i,j]),
                    "pooled_AB_sample_variance": float(indices["variance"][j]), "indices_clipped": False})
    save_csv(output / "sobol_prefix_estimates.csv", prefix_rows)
    final_indices = jansen_indices(ya, yb, yab)
    bootstrap = paired_bootstrap(ya, yb, yab, bootstrap_replicates)
    bootstrap_rows, interval_rows = [], []
    for i, p in enumerate(PARAMETERS):
        for j, metric in enumerate(OUTPUTS):
            first = bootstrap["first_order"][:, i, j]
            total = bootstrap["total_order"][:, i, j]
            interval_rows.append({"parameter": p["name"], "output": metric,
                "first_order_S": float(final_indices["first_order"][i,j]),
                "S_bootstrap_p025": float(np.quantile(first, .025)), "S_bootstrap_p975": float(np.quantile(first, .975)),
                "total_order_ST": float(final_indices["total_order"][i,j]),
                "ST_bootstrap_p025": float(np.quantile(total, .025)), "ST_bootstrap_p975": float(np.quantile(total, .975)),
                "interval_role": "paired_row_resampling_stability_not_QMC_confidence_guarantee"})
            for replicate in range(bootstrap_replicates):
                bootstrap_rows.append({"replicate": replicate, "parameter": p["name"], "output": metric,
                                       "first_order_S": float(first[replicate]), "total_order_ST": float(total[replicate])})
    save_csv(output / "sobol_bootstrap_replicates.csv", bootstrap_rows)
    save_csv(output / "sobol_resampling_intervals.csv", interval_rows)
    base = np.concatenate([ya, yb])
    quantiles = []
    for j, key in enumerate(OUTPUTS):
        values = base[:,j]
        quantiles.append({"output": key, "base_AB_count": 2*n,
            "minimum": float(values.min()), "p025": float(np.quantile(values,.025)),
            "median": float(np.quantile(values,.5)), "p975": float(np.quantile(values,.975)),
            "maximum": float(values.max()), "mean": float(values.mean()), "sample_sd": float(values.std(ddof=1)),
            "role": "scenario_distribution_not_experimental_confidence_interval"})
    save_csv(output / "base_AB_quantiles.csv", quantiles)
    grid_rows, grid_differences = [], []
    for index, point in enumerate(GRID_UNIT_POINTS):
        coarse, fine = evaluate_point(point), evaluate_point(point, nx=180, ny=72)
        grid_rows.extend([{"preselected_point": index, "grid_role": label, **row} for label,row in [("coarse",coarse),("refined",fine)]])
        for key in OUTPUTS:
            grid_differences.append({"preselected_point": index, "output": key,
                "coarse": coarse[key], "refined": fine[key], "coarse_minus_refined": coarse[key]-fine[key],
                "relative_difference_to_refined": (coarse[key]-fine[key])/fine[key]})
    save_csv(output / "grid_check_evaluations.csv", grid_rows)
    save_csv(output / "grid_check_differences.csv", grid_differences)
    local_rows, local_results = [], {}
    for label,i,step,sign,point in local_design():
        row = evaluate_point(point)
        local_rows.append({"label":label,"parameter_index":i,"unit_quantile_step":step,"sign":sign,**row})
        local_results[(i,step,sign)] = row
    save_csv(output / "local_difference_evaluations.csv", local_rows)
    center = local_results[(None,None,None)]
    physical_center = transform_unit(np.full(5,.5))
    local_derivatives = []
    for step in [.02,.01]:
        for i,p in enumerate(PARAMETERS):
            derivative_x_to_u = physical_center[i]*math.log(p["high"]/p["low"]) if p["distribution"] == "log_uniform" else p["high"]-p["low"]
            for key in OUTPUTS:
                du = (local_results[(i,step,1)][key]-local_results[(i,step,-1)][key])/(2*step)
                dx = du/derivative_x_to_u
                local_derivatives.append({"parameter":p["name"],"output":key,"unit_quantile_step":step,
                    "center_parameter":float(physical_center[i]),"center_output":center[key],
                    "derivative_dY_du":du,"derivative_dY_dphysical_parameter":float(dx),
                    "local_dimensionless_elasticity":float(dx*physical_center[i]/center[key]),
                    "role":"central_difference_in_unit_quantile_space_transformed_by_center_chain_rule"})
    save_csv(output / "local_derivatives.csv", local_derivatives)
    all_rows = records + grid_rows + local_rows
    indices_summary = {metric:{p["name"]:{"S":float(final_indices["first_order"][i,j]),"ST":float(final_indices["total_order"][i,j])}
                             for i,p in enumerate(PARAMETERS)} for j,metric in enumerate(OUTPUTS)}
    summary = {
        "schema_version":"electratwin.uncertainty.v1", "executed_utc":executed_utc,
        "evidence_role":"executed_scenario_propagation_not_calibrated_uncertainty_or_experiment",
        "parameters":list(PARAMETERS),"independence":"Independent selected input distributions; no empirical covariance or fitted posterior.",
        "fixed_assumptions":{"substrate_inlet_mM":50.,"product_inlet_mM":0.,"temperature_K":298.15,
            "length_m":.06,"width_m":.012,"n_electrons":2,"product_MW_g_mol":331.43,
            "voltage_law":"1.85 + eta + 18.5 * current_A","side_reactions":False},
        "sampling":{"method":"SciPy scrambled Sobol, LMS plus digital shift, random_base2",
            "dimension":10,"base_N":n,"base_power":power,"rng_seed":SCENARIO_SEED,
            "construction":"First five coordinates form A; last five B; AB_i copies A and replaces only column i by B_i.",
            "optimization":None,"skipped_points":0,"thinned_points":0,"coarse_mesh":[60,24],
            "main_PDE_solves":7*n,"base_AB_quantile_count":2*n,"cached_prefix_sizes":prefix_sizes},
        "jansen":{"S":"1 - mean((YB - YAB_i)^2)/(2*Var)","ST":"mean((YA - YAB_i)^2)/(2*Var)",
            "variance":"Sample variance ddof=1 of concatenated YA and YB for each output and prefix.",
            "negative_or_above_one_estimates_clipped":False,"constant_output_policy":"Raise ValueError; no fabricated indices."},
        "bootstrap":{"replicates":bootstrap_replicates,"rng_seed":BOOTSTRAP_SEED,
            "method":"Resample matched row indices simultaneously across A, B and every AB_i; reestimate pooled variance each time.",
            "reported_quantiles":[.025,.975],"interpretation":"Exploratory paired-row resampling stability only. Rows of a scrambled QMC net are not iid; this bootstrap destroys its balance and is not a rigorous QMC confidence interval. No independent scrambling replicates were run."},
        "sobol_at_full_N":indices_summary,"base_AB_quantiles":quantiles,
        "grid_check":{"preselected_unit_points":GRID_UNIT_POINTS.tolist(),"coarse_mesh":[60,24],"refined_mesh":[180,72],
            "PDE_solves":len(grid_rows),"point_selection":"Declared before examining outputs, not selected for agreement.",
            "maximum_absolute_conversion_difference_pp":max(abs(r["coarse_minus_refined"]) for r in grid_differences if r["output"]=="conversion_pct"),
            "maximum_absolute_relative_difference_by_output":{key:max(abs(r["relative_difference_to_refined"]) for r in grid_differences if r["output"]==key) for key in OUTPUTS},
            "scope":"Eight points only; no domain-wide discretization error bound or refined-grid Sobol indices."},
        "local_differences":{"PDE_solves":len(local_rows),"unit_quantile_center":[.5]*5,"unit_quantile_steps":[.02,.01],
            "interpretation":"Local slopes on the coarse model, not global Sobol effects. Physical derivatives apply the transform Jacobian at the center."},
        "counts":{"main_PDE_solves":7*n,"grid_check_PDE_solves":len(grid_rows),"local_difference_PDE_solves":len(local_rows),
            "total_PDE_solves":len(all_rows),"bootstrap_output_rows":len(bootstrap_rows),"additional_PDE_solves_for_prefix_or_bootstrap":0},
        "numerical_checks":{"all_converged":all(r["converged"] for r in all_rows),
            "all_reported_outputs_positive":all(all(r[k]>0 for k in OUTPUTS) for r in all_rows),
            "maximum_material_balance_relative_error":max(r["material_balance_relative_error"] for r in all_rows),
            "maximum_charge_balance_relative_error":max(r["charge_balance_relative_error"] for r in all_rows),
            "maximum_linear_residual_relative":max(r["linear_residual_relative"] for r in all_rows)},
        "limitations":["The distribution ranges and independence are selected scenarios, not measurement error bars or calibrated Bayesian posteriors.",
            "Sobol rankings depend on those distributions; eta and current scale need not be physically independent.",
            "A/P molar closure does not establish elemental mass closure for a real coupling reaction.",
            "Only one scrambling realization and modest base N; prefix changes and resampling diagnostics do not guarantee convergence.",
            "No negative sensitivity estimate is hidden by clipping; finite-sample S can be negative and ST can exceed one.",
            "Coarse-grid numerical error, model discrepancy and omitted side reactions/potential/heat/gas/migration are not integrated into the scenario distribution.",
            "Product mass and voltage law are assumed; STY and SEC are proxies and no instrument is connected."],
        "sources":[{"title":"SciPy 1.18 Sobol documentation","url":"https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.qmc.Sobol.html",
            "verified_date":"2026-09-28","access_scope":"Official HTML opened and read; power-of-two sampling, LMS+shift scrambling, rng API and loss of balance when skipping or thinning verified.",
            "supports":"QMC construction and sampling limits; does not endorse the selected parameter distributions or iid bootstrap confidence claims."}],
        "source_sha256":{"electratwin/scripts/uncertainty_analysis.py":sha256(Path(__file__)),
            "electratwin/scripts/transport_reviewed.py":sha256(BASE/"scripts/transport_reviewed.py")},
        "output_sha256":{path.name:sha256(path) for path in sorted(output.glob("*.csv"))},
        "runtime":{"python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__,
            "requested_numerical_threads":1,"elapsed_seconds":time.perf_counter()-started},
    }
    save_json(output/"summary.json",summary)
    print(json.dumps({"total_PDE_solves":summary["counts"]["total_PDE_solves"],"elapsed_seconds":summary["runtime"]["elapsed_seconds"],
                      "numerical_checks":summary["numerical_checks"],"grid_max_conversion_difference_pp":summary["grid_check"]["maximum_absolute_conversion_difference_pp"]}),flush=True)
    return summary


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output-dir",type=Path)
    parser.add_argument("--power",type=int,default=8)
    parser.add_argument("--bootstrap-replicates",type=int,default=500)
    args=parser.parse_args()
    run_analysis(args.output_dir,args.power,args.bootstrap_replicates)


if __name__=="__main__":
    main()
