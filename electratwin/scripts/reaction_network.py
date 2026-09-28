"""Explicit four-species wall-network demonstrator with conservative transport.

A -> P (desired), A -> B (parallel), P -> D (overoxidation). These labels and
rate/electron parameters have no identified molecular mechanism or calibration.
No hardware access. Run from the repository root with numerical threads = 1.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import platform
from pathlib import Path
import sys
import time
import warnings

import numpy as np
import scipy
from scipy.sparse import coo_matrix, eye, kron
from scipy.sparse.linalg import MatrixRankWarning, spsolve

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from electratwin.scripts.transport_reviewed import solve_transport

F = 96485.33212
R = 8.31446261815324
SPECIES = ("A", "P", "B", "D")
CHANNELS = ("AP", "AB", "PD")
DEFAULT_REFERENCE_RATES = (5e-5, 3e-6, 8e-6)
DEFAULT_TRANSFER_FACTORS = (.50, .65, .85)
OUTPUT = ROOT / "electratwin/results/reaction_network"


def _finite(name, value, *, positive=False, nonnegative=False):
    if isinstance(value, (bool, str)) or not np.isfinite(value) or (positive and value <= 0) or (nonnegative and value < 0):
        raise ValueError(f"Invalid {name}")


def network_rates(potential_index_V=.48, *, temperature_K=298.15,
                  reference_potential_V=.45, reference_rates_m_s=DEFAULT_REFERENCE_RATES,
                  transfer_factors=DEFAULT_TRANSFER_FACTORS):
    """Illustrative irreversible exponential kinetics, never fitted to data."""
    _finite("potential index", potential_index_V, nonnegative=True)
    _finite("reference potential", reference_potential_V, nonnegative=True)
    _finite("temperature", temperature_K, positive=True)
    if len(reference_rates_m_s) != 3 or len(transfer_factors) != 3:
        raise ValueError("Exactly three rate constants and transfer factors required")
    rates = []
    for base, beta in zip(reference_rates_m_s, transfer_factors):
        _finite("reference rate", base, nonnegative=True)
        _finite("transfer factor", beta, nonnegative=True)
        exponent = beta * F * (potential_index_V - reference_potential_V) / (R * temperature_K)
        if abs(exponent) > 100:
            raise ValueError("Potential outside bounded illustrative kinetics")
        rates.append(float(base * math.exp(exponent)))
    return dict(zip(CHANNELS, rates))


def _validate_rates(rates):
    if not isinstance(rates, dict) or set(rates) != set(CHANNELS):
        raise ValueError("Rates must contain exactly AP, AB, PD; reverse pathways are not supported")
    for name, value in rates.items():
        _finite(name + " rate", value, nonnegative=True)
    return np.array([rates[k] for k in CHANNELS], dtype=float)


def loss_matrix(rates):
    k1, k2, k3 = _validate_rates(rates)
    return np.array([[k1 + k2, 0, 0, 0], [-k1, k3, 0, 0],
                     [-k2, 0, 0, 0], [0, -k3, 0, 0]], dtype=float)


def analytic_plug_outlet(inlet_mM, rates, area_over_flow_s_m):
    """Independent mixed-cross-section PFR limit, in surface exposure A_wall/Q.

    Closed form for parallel plus consecutive first-order reactions. It is only
    the D_y -> infinity, D_x=0, plug-velocity limit of the spatial model.
    """
    feed = np.asarray(inlet_mM, dtype=float)
    if feed.shape != (4,) or not np.isfinite(feed).all() or np.any(feed < 0) or feed.sum() <= 0:
        raise ValueError("Require four nonnegative finite concentrations and a nonzero total")
    _finite("surface exposure", area_over_flow_s_m, nonnegative=True)
    k1, k2, k3 = _validate_rates(rates)
    exposure = area_over_flow_s_m
    a, b = k1 + k2, k3
    ca = feed[0] * math.exp(-a * exposure)
    if a == b:
        made_p = k1 * feed[0] * exposure * math.exp(-a * exposure)
    else:
        gap = abs(b - a)
        made_p = k1 * feed[0] * math.exp(-min(a, b) * exposure) * -math.expm1(-gap * exposure) / gap
    cp = feed[1] * math.exp(-b * exposure) + made_p
    cb = feed[2] + (k2 / a * feed[0] * -math.expm1(-a * exposure) if a else 0)
    cd = feed.sum() - ca - cp - cb
    return dict(zip(SPECIES, (float(ca), float(cp), float(cb), float(cd))))


def solve_network(flow_rate_uL_min=450.0, potential_index_V=.48, *, nx=80, ny=32,
                  inlet_mM=(50., 0., 0., 0.), length_m=.06, height_m=.0003, width_m=.012,
                  diffusivity_m2_s=1.1e-9, axial_diffusivity_m2_s=None,
                  temperature_K=298.15, rate_constants_m_s=None, electrons=(2, 2, 2),
                  velocity_profile="poiseuille", include_field=True):
    """Solve four species explicitly; no species is reconstructed by subtraction.

    div(u C_s-D grad C_s)=0; Danckwerts total-flux inlet, convective outlet,
    no-flux upper wall; lower-wall outward flux J=L C_wall. Robin elimination:
    C_wall=(I+dy L/(2D))^-1 C_firstcell. Exactly shared sinks/sources and reaction
    currents preserve species balances and total molecular count, not elements.
    """
    started = time.perf_counter()
    for name, value in [("Q", flow_rate_uL_min), ("L", length_m), ("H", height_m),
                        ("W", width_m), ("temperature", temperature_K)]:
        _finite(name, value, positive=True)
    _finite("D", diffusivity_m2_s, nonnegative=True)
    _finite("potential index", potential_index_V, nonnegative=True)
    for value in [nx, ny]:
        _finite("grid count", value, positive=True)
        if int(value) != value or value < 2:
            raise ValueError("Require integer grid counts >=2")
    nx, ny = int(nx), int(ny)
    if nx * ny * 4 > 50000:
        raise ValueError("Bounded CPU pilot: at most 50000 total species unknowns")
    feed = np.asarray(inlet_mM, dtype=float)
    if feed.shape != (4,) or not np.isfinite(feed).all() or np.any(feed < 0) or feed.sum() <= 0:
        raise ValueError("Require four nonnegative finite feed concentrations and nonzero total")
    if len(electrons) != 3:
        raise ValueError("Require three electron counts")
    for value in electrons:
        _finite("electron count", value, positive=True)
        if int(value) != value:
            raise ValueError("Electron counts must be positive integers")
    n_e = np.array(electrons, dtype=int)
    if velocity_profile not in {"poiseuille", "plug"}:
        raise ValueError("Unknown velocity profile")
    rates = network_rates(potential_index_V, temperature_K=temperature_K) if rate_constants_m_s is None else rate_constants_m_s
    k = _validate_rates(rates)
    L = loss_matrix(rates)
    dax = diffusivity_m2_s if axial_diffusivity_m2_s is None else axial_diffusivity_m2_s
    _finite("axial D", dax, nonnegative=True)
    dx, dy = length_m / nx, height_m / ny
    x, y = (np.arange(nx) + .5) * dx, (np.arange(ny) + .5) * dy
    q = flow_rate_uL_min * 1e-9 / 60
    u_avg = q / (width_m * height_m)
    edges = np.arange(ny + 1) / ny
    u = u_avg * ny * np.diff(3 * edges**2 - 2 * edges**3) if velocity_profile == "poiseuille" else np.full(ny, u_avg)
    adv = u * dy * width_m
    gx, gy = dax * dy * width_m / dx, diffusivity_m2_s * dx * width_m / dy
    # D=0 is an uncoupled wall-access limit. Wall values are undefined there.
    wall_map = np.linalg.solve(np.eye(4) + dy / (2 * diffusivity_m2_s) * L, np.eye(4)) if diffusivity_m2_s else None
    robin = L @ wall_map if wall_map is not None else np.zeros((4, 4))
    rows, cols, values = [], [], []
    for i in range(nx):
        for j in range(ny):
            p = i * ny + j
            diagonal = adv[j]
            if i > 0:
                rows.append(p); cols.append(p - ny); values.append(-adv[j] - gx)
                diagonal += gx
            if i < nx - 1:
                rows.append(p); cols.append(p + ny); values.append(-gx)
                diagonal += gx
            if j > 0:
                rows.append(p); cols.append(p - 1); values.append(-gy)
                diagonal += gy
            if j < ny - 1:
                rows.append(p); cols.append(p + 1); values.append(-gy)
                diagonal += gy
            rows.append(p); cols.append(p); values.append(diagonal)
    transport = coo_matrix((values, (rows, cols)), shape=(nx * ny, nx * ny)).tocsr()
    matrix = kron(transport, eye(4), format="csr")
    wr, wc, wv = [], [], []
    area = dx * width_m
    for i in range(nx):
        cell = i * ny * 4
        for s in range(4):
            for t in range(4):
                if robin[s, t]:
                    wr.append(cell + s); wc.append(cell + t); wv.append(robin[s, t] * area)
    matrix += coo_matrix((wv, (wr, wc)), shape=matrix.shape).tocsr()
    rhs = np.zeros((nx, ny, 4))
    rhs[0, :, :] = adv[:, None] * feed[None, :]
    with warnings.catch_warnings():
        warnings.simplefilter("error", MatrixRankWarning)
        c = spsolve(matrix, rhs.ravel()).reshape(nx, ny, 4)
    wall_c = c[:, 0, :] @ wall_map.T if wall_map is not None else None
    wall_flux = wall_c @ L.T if wall_c is not None else np.zeros((nx, 4))
    reaction_profiles = np.column_stack([k[0] * wall_c[:, 0], k[1] * wall_c[:, 0], k[2] * wall_c[:, 1]]) if wall_c is not None else np.zeros((nx, 3))
    extents = reaction_profiles.sum(axis=0) * area
    channel_current = extents * n_e * F
    current = float(channel_current.sum())
    inlet_flows, outlet_flows = q * feed, adv @ c[-1, :, :]
    net_gain = outlet_flows - inlet_flows
    wall_integrals = wall_flux.sum(axis=0) * area
    species_residual = inlet_flows - outlet_flows - wall_integrals
    mass_scale = float(q * feed.sum())
    max_species_error = float(np.max(np.abs(species_residual)) / mass_scale)
    total_error = float(abs(outlet_flows.sum() - inlet_flows.sum()) / mass_scale)
    extent_bulk = np.array([net_gain[1] + net_gain[3], net_gain[2], net_gain[3]])
    charge_bulk = float(F * np.dot(n_e, extent_bulk))
    charge_error = abs(current - charge_bulk) / max(F * max(n_e) * mass_scale, abs(current))
    residual = float(np.max(np.abs(matrix @ c.ravel() - rhs.ravel())))
    relative_residual = residual / float(np.max(np.abs(rhs)))
    minimum = min(float(c.min()), float(wall_c.min()) if wall_c is not None else float(c.min()))
    nonnegative = minimum >= -1e-10 * float(feed.sum())
    converged = bool(np.isfinite(c).all() and relative_residual < 1e-8 and max_species_error < 1e-8
                     and total_error < 1e-8 and charge_error < 1e-8 and nonnegative)
    total_extent = float(extents[0] + extents[1])
    net_fe = float(n_e[0] * F * net_gain[1] / current * 100) if current > 0 else None
    summary = {
        "converged": converged, "linear_residual_relative": relative_residual,
        "max_species_balance_relative_error": max_species_error,
        "total_molar_balance_relative_error": total_error,
        "charge_balance_relative_error": charge_error,
        "species_balance_residual_mol_s": dict(zip(SPECIES, map(float, species_residual))),
        "outlet_mM": dict(zip(SPECIES, map(float, outlet_flows / q))),
        "inlet_molar_flow_mol_s": dict(zip(SPECIES, map(float, inlet_flows))),
        "outlet_molar_flow_mol_s": dict(zip(SPECIES, map(float, outlet_flows))),
        "net_outlet_gain_mol_s": dict(zip(SPECIES, map(float, net_gain))),
        "wall_net_outward_flux_mol_s": dict(zip(SPECIES, map(float, wall_integrals))),
        "extent_mol_s": dict(zip(CHANNELS, map(float, extents))),
        "channel_current_A": dict(zip(CHANNELS, map(float, channel_current))),
        "current_A": current, "current_from_bulk_balance_A": charge_bulk,
        "conversion_A_pct": float((1 - outlet_flows[0] / inlet_flows[0]) * 100) if inlet_flows[0] > 0 else None,
        "net_P_yield_pct": float(net_gain[1] / inlet_flows[0] * 100) if inlet_flows[0] > 0 else None,
        "net_P_selectivity_pct": float(net_gain[1] / total_extent * 100) if total_extent > 0 else None,
        "net_P_faradaic_efficiency_pct": net_fe,
        "gross_AP_charge_fraction_pct": float(channel_current[0] / current * 100) if current > 0 else None,
        "nonnegative": bool(nonnegative), "minimum_concentration_mM": minimum,
        "maximum_total_concentration_error_mM": float(np.max(np.abs(c.sum(axis=2) - feed.sum()))),
        "discrete_flow_m3_s": float(adv.sum()), "prescribed_flow_m3_s": q,
        "reactor_volume_m3": length_m * height_m * width_m,
        "residence_time_s": length_m * height_m * width_m / q,
        "total_unknowns": nx * ny * 4, "solve_elapsed_s": time.perf_counter() - started,
    }
    result = {
        "schema_version": "electratwin.wall_network.v1",
        "evidence_class": "executed_hypothetical_network_not_experimental_validation",
        "species": list(SPECIES), "reaction_channels": ["A -> P", "A -> B", "P -> D"],
        "assumptions": ["Abstract species labels and uncalibrated irreversible first-order wall kinetics",
                        "One molecule enters and one molecule leaves each step: total molar count conserved, not a claim of elemental balance",
                        "Equal constant diffusivity, prescribed flow and temperature; no potential, migration, thermal or gas field",
                        "Three electron counts default to 2 each by assumption, not a molecular mechanism",
                        "Net desired FE is charge-equivalent net P outlet gain divided by all three integrated reaction currents",
                        "With P in the inlet, net desired FE may be negative and is not an absolute product-generation FE",
                        "First-order axial upwind convection, central diffusion and exact matrix Robin elimination"],
        "inputs": {"flow_rate_uL_min": flow_rate_uL_min, "potential_index_V": potential_index_V,
                   "nx": nx, "ny": ny, "inlet_mM": feed.tolist(), "length_m": length_m,
                   "height_m": height_m, "width_m": width_m, "diffusivity_m2_s": diffusivity_m2_s,
                   "axial_diffusivity_m2_s": dax, "temperature_K": temperature_K,
                   "rate_constants_m_s": dict(zip(CHANNELS, map(float, k))),
                   "explicit_rate_override": rate_constants_m_s is not None, "electrons": n_e.tolist(),
                   "reference_rates_m_s": dict(zip(CHANNELS, DEFAULT_REFERENCE_RATES)),
                   "transfer_factors": dict(zip(CHANNELS, DEFAULT_TRANSFER_FACTORS)),
                   "reference_potential_V": .45, "velocity_profile": velocity_profile, "F_C_mol": F},
        "summary": summary,
    }
    if include_field:
        result["field"] = {"x_mm": (x * 1000).tolist(), "y_um": (y * 1e6).tolist(),
                           "velocity_m_s": u.tolist(),
                           "concentration_mM": {name: c[:, :, s].tolist() for s, name in enumerate(SPECIES)},
                           "wall_concentration_mM": {name: wall_c[:, s].tolist() for s, name in enumerate(SPECIES)} if wall_c is not None else None,
                           "wall_outward_flux_mol_m2_s": {name: wall_flux[:, s].tolist() for s, name in enumerate(SPECIES)},
                           "reaction_flux_mol_m2_s": {name: reaction_profiles[:, t].tolist() for t, name in enumerate(CHANNELS)}}
    return result


def _json(path, data):
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")


def _csv(path, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def _row(result):
    s, p = result["summary"], result["inputs"]
    scalar = {key: value for key, value in s.items() if not isinstance(value, dict)}
    return {"flow_rate_uL_min": p["flow_rate_uL_min"], "potential_index_V": p["potential_index_V"],
            "nx": p["nx"], "ny": p["ny"], **scalar,
            **{f"outlet_{key}_mM": value for key, value in s["outlet_mM"].items()},
            **{f"current_{key}_A": value for key, value in s["channel_current_A"].items()},
            **{f"extent_{key}_mol_s": value for key, value in s["extent_mol_s"].items()},
            **{f"balance_{key}_mol_s": value for key, value in s["species_balance_residual_mol_s"].items()}}


def run_study(output=OUTPUT):
    output = Path(output); output.mkdir(parents=True, exist_ok=True)
    registry, summaries = [], []
    def evaluate(label, **kwargs):
        r = solve_network(**kwargs)
        registry.append({"network_solve_id": len(registry) + 1, "case": label, **_row(r)})
        summaries.append(r["summary"])
        return r
    baseline = evaluate("baseline")
    _json(output / "baseline_solution.json", baseline)
    pool = []
    for q in np.linspace(100, 1500, 7):
        for eta in np.linspace(.2, .75, 7):
            pool.append(_row(evaluate("parameter_pool", flow_rate_uL_min=float(q), potential_index_V=float(eta), include_field=False)))
    _csv(output / "parameter_pool.csv", pool)
    grid = [_row(evaluate("joint_grid", nx=nx, ny=ny, include_field=False)) for nx, ny in [(40, 16), (80, 32), (160, 64)]]
    fine = grid[-1]
    for row in grid:
        row["P_yield_delta_from_finest_pp"] = row["net_P_yield_pct"] - fine["net_P_yield_pct"]
        row["FE_delta_from_finest_pp"] = row["net_P_faradaic_efficiency_pct"] - fine["net_P_faradaic_efficiency_pct"]
    _csv(output / "grid_convergence.csv", grid)
    limits = []
    explicit_k = 5e-5
    cases = [("no_side_reactions", dict(rate_constants_m_s={"AP": explicit_k, "AB": 0, "PD": 0})),
             ("parallel_only", dict(rate_constants_m_s={"AP": 2e-5, "AB": 1e-5, "PD": 0})),
             ("zero_reaction", dict(rate_constants_m_s={"AP": 0, "AB": 0, "PD": 0})),
             ("zero_diffusion", dict(diffusivity_m2_s=0)),
             ("P_feed_overoxidation", dict(inlet_mM=(0, 50, 0, 0), rate_constants_m_s={"AP": 0, "AB": 0, "PD": 2e-5}))]
    for name, args in cases:
        r = evaluate(name, include_field=False, **args)
        limits.append({"case": name, **_row(r)})
    _csv(output / "limiting_cases.csv", limits)
    legacy = solve_transport(nx=80, ny=32, forward_rate_m_s=explicit_k, reverse_rate_m_s=0, include_field=False)
    old_comparison = {"network_conversion_A_pct": limits[0]["conversion_A_pct"],
                      "frozen_solver_conversion_pct": legacy["summary"]["conversion_pct"],
                      "conversion_difference_pp": limits[0]["conversion_A_pct"] - legacy["summary"]["conversion_pct"],
                      "current_difference_A": limits[0]["current_A"] - legacy["summary"]["current_A"],
                      "matched_conditions": "80x32; k_AP=5e-5 m/s, k_AB=k_PD=0, frozen solver reverse rate=0",
                      "frozen_solver_additional_PDE_solves": 1}
    _json(output / "frozen_solver_comparison.json", old_comparison)
    analytic = []
    rates = {"AP": 2e-6, "AB": 0., "PD": 4e-6}
    exact = analytic_plug_outlet((50, 0, 0, 0), rates, .06 * .012 / (450e-9 / 60))
    for nx in [40, 80, 160, 320]:
        r = evaluate("analytic_sequential_plug_limit", nx=nx, ny=8, diffusivity_m2_s=1e-4,
                     axial_diffusivity_m2_s=0, velocity_profile="plug", rate_constants_m_s=rates,
                     include_field=False)
        analytic.append({"analytic_P_mM": exact["P"], "P_error_mM": r["summary"]["outlet_mM"]["P"] - exact["P"],
                         "max_species_error_mM": max(abs(r["summary"]["outlet_mM"][s] - exact[s]) for s in SPECIES), **_row(r)})
    _csv(output / "analytic_sequential_limit.csv", analytic)
    invalid = []
    for name, args in [("reverse_flow", {"flow_rate_uL_min": -10}), ("negative_feed", {"inlet_mM": (50, -1, 0, 0)}),
                       ("zero_feed", {"inlet_mM": (0, 0, 0, 0)}), ("reverse_rate", {"rate_constants_m_s": {"AP": -1, "AB": 0, "PD": 0}}),
                       ("unsupported_reverse_path", {"rate_constants_m_s": {"AP": 1, "AB": 0, "PD": 0, "PA": 1}}),
                       ("unknown_limit", {"nx": 200, "ny": 100}), ("nonfinite_feed", {"inlet_mM": (float("nan"), 0, 0, 0)}),
                       ("negative_potential_index", {"potential_index_V": -.2})]:
        try:
            solve_network(**args)
            invalid.append({"case": name, "rejected_before_PDE_solve": False, "error": None})
        except ValueError as exc:
            invalid.append({"case": name, "rejected_before_PDE_solve": True, "error": str(exc)})
    _json(output / "input_rejections.json", invalid)
    _csv(output / "solve_registry.csv", registry)
    sources = [Path(__file__), ROOT / "electratwin/scripts/transport_reviewed.py", ROOT / "tests/test_electratwin_network.py"]
    source_hashes = {str(path.relative_to(ROOT)).replace("\\", "/"): hashlib.sha256(path.read_bytes()).hexdigest() for path in sources if path.exists()}
    low_fe = min(pool, key=lambda r: r["net_P_faradaic_efficiency_pct"])
    high_fe = max(pool, key=lambda r: r["net_P_faradaic_efficiency_pct"])
    summary = {
        "schema_version": "electratwin.network.study.v1", "evidence_class": "numerical_verification_of_hypothetical_network",
        "baseline": baseline["summary"], "baseline_inputs": baseline["inputs"], "source_sha256": source_hashes,
        "runtime": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
        "PDE_counts": {"network_solves": len(registry), "frozen_comparison_solves": 1, "total_PDE_solves": len(registry) + 1,
                       "baseline": 1, "pool": 49, "joint_grid": 3, "limits": 5, "analytic_limits": 4,
                       "invalid_inputs_rejected_without_solving": len(invalid), "unit_test_solves_excluded": True},
        "all_network_cases_converged": all(s["converged"] for s in summaries),
        "all_invalid_inputs_rejected": all(r["rejected_before_PDE_solve"] for r in invalid),
        "maximum_species_balance_relative_error": max(s["max_species_balance_relative_error"] for s in summaries),
        "maximum_total_molar_balance_relative_error": max(s["total_molar_balance_relative_error"] for s in summaries),
        "maximum_charge_balance_relative_error": max(s["charge_balance_relative_error"] for s in summaries),
        "maximum_total_unknowns": max(s["total_unknowns"] for s in summaries),
        "joint_grid_P_yield_pct": [r["net_P_yield_pct"] for r in grid],
        "joint_grid_FE_pct": [r["net_P_faradaic_efficiency_pct"] for r in grid],
        "last_joint_refinement_P_yield_delta_pp": abs(grid[-1]["net_P_yield_pct"] - grid[-2]["net_P_yield_pct"]),
        "last_joint_refinement_FE_delta_pp": abs(grid[-1]["net_P_faradaic_efficiency_pct"] - grid[-2]["net_P_faradaic_efficiency_pct"]),
        "pool_min_FE": low_fe, "pool_max_FE": high_fe,
        "analytic_limit_outlet_mM": exact, "finest_analytic_max_species_error_mM": analytic[-1]["max_species_error_mM"],
        "frozen_solver_comparison": old_comparison,
        "limits": ["All kinetics, potential sensitivities, abstract species and electron stoichiometries are assumptions",
                   "Total molecular count is conserved; no atom mapping or elemental reaction balance has been established",
                   "49 conditions are model calculations, not observations or a validated operating window",
                   "Three grids are a numerical sensitivity study; first-order upwind local-gradient errors remain",
                   "Negative net P FE in product-feed degradation is retained; no clipping or claimed measured selectivity",
                   "No hardware, thermal field, electrostatic field, migration, gas, catalyst degradation, or safety validation"]
    }
    _json(output / "study_summary.json", summary)
    print(json.dumps({key: summary[key] for key in ["PDE_counts", "all_network_cases_converged",
                      "maximum_species_balance_relative_error", "maximum_total_molar_balance_relative_error",
                      "maximum_charge_balance_relative_error", "joint_grid_P_yield_pct", "joint_grid_FE_pct",
                      "finest_analytic_max_species_error_mM"]}, indent=2))
    return summary


if __name__ == "__main__":
    run_study()
