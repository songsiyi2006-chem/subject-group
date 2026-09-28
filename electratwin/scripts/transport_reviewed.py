"""Bounded, conservative 2D transport demonstrator; no validated reactor claim.

Run from the repository root: python electratwin/scripts/transport_reviewed.py
All concentrations use mol/m^3, numerically equal to mM. No hardware access.
"""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path
import platform
import time
import warnings

import numpy as np
import scipy
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import MatrixRankWarning, spsolve

F = 96485.33212
R = 8.31446261815324
ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "electratwin/results/transport"


def _finite(name, value, *, positive=False, nonnegative=False):
    if not np.isfinite(value) or (positive and value <= 0) or (nonnegative and value < 0):
        raise ValueError(f"Invalid {name}: require finite " + ("positive" if positive else "nonnegative" if nonnegative else "value"))


def effective_rates(overpotential_V, current_scale_A_m2=0.08,
                    reference_concentration_mM=50.0, temperature_K=298.15,
                    alpha=0.5, n_electrons=2):
    """Phenomenological rates, not a mechanistic multi-electron BV derivation.

    eta is an imposed index relative to a fixed reference state, not a solved
    interfacial potential. The exponent uses F (effective transfer parameter);
    nF converts net molecular turnover to current. The current scale equals j0
    only under the model's specified equal reference concentrations.
    """
    for name, value in [("eta", overpotential_V), ("alpha", alpha)]:
        _finite(name, value)
    for name, value in [("reference concentration", reference_concentration_mM), ("temperature", temperature_K)]:
        _finite(name, value, positive=True)
    _finite("current scale", current_scale_A_m2, nonnegative=True)
    if not 0 < alpha < 1 or isinstance(n_electrons, bool) or int(n_electrons) != n_electrons or n_electrons < 1:
        raise ValueError("Require 0 < alpha < 1 and positive integer electron count")
    anodic = alpha * F * overpotential_V / (R * temperature_K)
    cathodic = -(1 - alpha) * F * overpotential_V / (R * temperature_K)
    if max(abs(anodic), abs(cathodic)) > 100:
        raise ValueError("Potential outside the bounded illustrative rate model")
    k0 = current_scale_A_m2 / (n_electrons * F * reference_concentration_mM)
    return k0 * math.exp(anodic), k0 * math.exp(cathodic)


def solve_transport(flow_rate_uL_min=450.0, overpotential_V=0.48, *,
                    nx=100, ny=48, inlet_A_mM=50.0, inlet_P_mM=0.0,
                    length_m=0.06, height_m=0.0003, width_m=0.012,
                    diffusivity_m2_s=1.1e-9, temperature_K=298.15,
                    current_scale_A_m2=0.08, reference_concentration_mM=50.0,
                    alpha=0.5, n_electrons=2, forward_rate_m_s=None,
                    reverse_rate_m_s=None, axial_diffusivity_m2_s=None,
                    velocity_profile="poiseuille", include_field=True):
    """Solve div(u C - D grad C)=0 with conservative cell-centered volumes.

    Inlet: prescribed total flux u*C_in (Danckwerts). Outlet: zero diffusive
    flux, outgoing upwind convection. Top wall: no flux. Bottom: A <-> P,
    J=kf*C_A,s-kr*C_P,s. Equal species D and 1:1 molecular balance imply
    C_A+C_P=C_total. Half-cell diffusion resistance determines wall C.

    Default axial D equals transverse D. A distinct axial D and plug velocity
    are provided only for explicit diagnostic limits. No concentration, yield,
    current, or convergence flag is clipped or forced successful.
    """
    started = time.perf_counter()
    for name, value in [("flow", flow_rate_uL_min), ("length", length_m),
                        ("height", height_m), ("width", width_m), ("temperature", temperature_K)]:
        _finite(name, value, positive=True)
    for name, value in [("inlet A", inlet_A_mM), ("inlet P", inlet_P_mM), ("D", diffusivity_m2_s)]:
        _finite(name, value, nonnegative=True)
    c_total = inlet_A_mM + inlet_P_mM
    if c_total <= 0:
        raise ValueError("At least one inlet species concentration must be positive")
    if any(isinstance(v, bool) or int(v) != v or v < 2 for v in (nx, ny)) or nx * ny > 50000:
        raise ValueError("Require integer grid counts >=2 and <=50000 cells")
    nx, ny = int(nx), int(ny)
    if velocity_profile not in {"poiseuille", "plug"}:
        raise ValueError("velocity_profile must be poiseuille or plug")
    kf, kr = effective_rates(overpotential_V, current_scale_A_m2,
                             reference_concentration_mM, temperature_K, alpha, n_electrons)
    if (forward_rate_m_s is None) != (reverse_rate_m_s is None):
        raise ValueError("Supply both explicit forward and reverse rates or neither")
    if forward_rate_m_s is not None:
        _finite("forward rate", forward_rate_m_s, nonnegative=True)
        _finite("reverse rate", reverse_rate_m_s, nonnegative=True)
        kf, kr = float(forward_rate_m_s), float(reverse_rate_m_s)
    dax = diffusivity_m2_s if axial_diffusivity_m2_s is None else axial_diffusivity_m2_s
    _finite("axial D", dax, nonnegative=True)
    dx, dy = length_m / nx, height_m / ny
    x = (np.arange(nx) + 0.5) * dx
    y = (np.arange(ny) + 0.5) * dy
    q = flow_rate_uL_min * 1e-9 / 60
    u_avg = q / (width_m * height_m)
    # Exact strip-averaged Poiseuille velocity conserves imposed Q on any grid.
    eta_edges = np.arange(ny + 1) / ny
    primitive = 3 * eta_edges**2 - 2 * eta_edges**3
    u = u_avg * ny * np.diff(primitive) if velocity_profile == "poiseuille" else np.full(ny, u_avg)
    f = u * dy * width_m
    gx = dax * dy * width_m / dx
    gy = diffusivity_m2_s * dx * width_m / dy
    if diffusivity_m2_s == 0:
        wall_k, wall_b = 0.0, 0.0
    else:
        resistance = dy / (2 * diffusivity_m2_s)
        denominator = 1 + (kf + kr) * resistance
        wall_k = (kf + kr) / denominator
        wall_b = kr * c_total / denominator
    wall_area = dx * width_m
    rows, cols, values = [], [], []
    rhs = np.zeros(nx * ny)
    for i in range(nx):
        for j in range(ny):
            p = i * ny + j
            diagonal = f[j]  # outflow east; inlet supplied on RHS at i=0
            if i == 0:
                rhs[p] = f[j] * inlet_A_mM
            else:
                rows.append(p); cols.append(p - ny); values.append(-f[j])
            if i > 0:
                diagonal += gx
                rows.append(p); cols.append(p - ny); values.append(-gx)
            if i < nx - 1:
                diagonal += gx
                rows.append(p); cols.append(p + ny); values.append(-gx)
            if j > 0:
                diagonal += gy
                rows.append(p); cols.append(p - 1); values.append(-gy)
            if j < ny - 1:
                diagonal += gy
                rows.append(p); cols.append(p + 1); values.append(-gy)
            if j == 0:
                diagonal += wall_k * wall_area
                rhs[p] += wall_b * wall_area
            rows.append(p); cols.append(p); values.append(diagonal)
    matrix = coo_matrix((values, (rows, cols)), shape=(nx * ny, nx * ny)).tocsr()
    with warnings.catch_warnings():
        warnings.simplefilter("error", MatrixRankWarning)
        concentration = spsolve(matrix, rhs).reshape(nx, ny)
    product = c_total - concentration
    wall_flux = wall_k * concentration[:, 0] - wall_b
    if diffusivity_m2_s == 0:
        # Uncoupled zero-D wall does not transmit species from the fluid.
        wall_A = np.full(nx, kr * c_total / (kf + kr)) if kf + kr else concentration[:, 0].copy()
    else:
        wall_A = concentration[:, 0] - wall_flux * dy / (2 * diffusivity_m2_s)
    wall_P = c_total - wall_A
    wall_current = n_electrons * F * wall_flux
    reacted = float(np.sum(wall_flux) * wall_area)
    current = float(np.sum(wall_current) * wall_area)
    inlet_A_flow = q * inlet_A_mM
    outlet_A_flow = float(f @ concentration[-1])
    outlet_P_flow = float(f @ product[-1])
    balance = inlet_A_flow - outlet_A_flow - reacted
    mass_scale = max(q * c_total, abs(reacted), np.finfo(float).tiny)
    charge_from_bulk = n_electrons * F * (inlet_A_flow - outlet_A_flow)
    charge_scale = max(n_electrons * F * q * c_total, abs(current), np.finfo(float).tiny)
    linear_residual = float(np.linalg.norm(matrix @ concentration.ravel() - rhs, ord=np.inf))
    relative_linear_residual = linear_residual / max(float(np.linalg.norm(rhs, ord=np.inf)), np.finfo(float).tiny)
    min_species = min(float(concentration.min()), float(product.min()), float(wall_A.min()), float(wall_P.min()))
    nonnegative = min_species >= -1e-10 * c_total
    converged = bool(np.isfinite(concentration).all() and relative_linear_residual < 1e-9
                     and abs(balance) / mass_scale < 1e-8 and nonnegative)
    summary = {
        "converged": converged, "solver": "scipy.sparse.linalg.spsolve",
        "linear_residual_mol_s": linear_residual,
        "linear_residual_relative": relative_linear_residual,
        "material_balance_residual_mol_s": balance,
        "material_balance_relative_error": abs(balance) / mass_scale,
        "charge_balance_relative_error": abs(current - charge_from_bulk) / charge_scale,
        "charge_from_bulk_A": charge_from_bulk, "current_A": current,
        "reacted_mol_s": reacted,
        "outlet_product_gain_mol_s": outlet_P_flow - q * inlet_P_mM,
        "inlet_A_mol_s": inlet_A_flow, "outlet_A_mol_s": outlet_A_flow,
        "outlet_A_mM": outlet_A_flow / q, "outlet_P_mM": outlet_P_flow / q,
        "conversion_pct": (1 - outlet_A_flow / inlet_A_flow) * 100 if inlet_A_flow else None,
        "net_product_from_total_feed_pct": reacted / (q * c_total) * 100,
        "nonnegative": nonnegative, "minimum_species_mM": min_species,
        "minimum_bulk_A_mM": float(concentration.min()), "maximum_bulk_A_mM": float(concentration.max()),
        "minimum_bulk_P_mM": float(product.min()), "maximum_bulk_P_mM": float(product.max()),
        "maximum_total_species_deviation_mM": float(np.max(np.abs(concentration + product - c_total))),
        "discrete_flow_rate_m3_s": float(f.sum()), "prescribed_flow_rate_m3_s": q,
        "average_velocity_m_s": u_avg, "reactor_volume_m3": length_m * height_m * width_m,
        "residence_time_s": length_m * height_m * width_m / q,
        "forward_rate_m_s": kf, "reverse_rate_m_s": kr,
        "solve_elapsed_s": time.perf_counter() - started,
    }
    inputs = dict(flow_rate_uL_min=float(flow_rate_uL_min), overpotential_V=float(overpotential_V),
                  nx=nx, ny=ny, inlet_A_mM=float(inlet_A_mM), inlet_P_mM=float(inlet_P_mM),
                  length_m=length_m, height_m=height_m, width_m=width_m,
                  diffusivity_m2_s=diffusivity_m2_s, axial_diffusivity_m2_s=dax,
                  temperature_K=temperature_K, current_scale_A_m2=current_scale_A_m2,
                  reference_concentration_mM=reference_concentration_mM, alpha=alpha,
                  n_electrons=n_electrons, F_C_mol=F, R_J_mol_K=R,
                  velocity_profile=velocity_profile, explicit_rates=forward_rate_m_s is not None)
    result = {
        "schema_version": "electratwin.transport.v1", "evidence_class": "executed_unvalidated_model",
        "model": "Steady 2D finite-volume equal-diffusivity A/P transport with effective reversible Robin wall",
        "assumptions": ["1:1 A/P molecular bookkeeping, constant density and equal diffusion coefficients",
                        "Prescribed isothermal velocity and potential index; no potential, heat, migration, or gas equation",
                        "Current is nF times net wall molecular flux; no parasitic current or selectivity model",
                        "Danckwerts total-flux inlet; zero axial diffusive outlet flux; upper wall insulating",
                        "First-order upwind convection; central diffusion; half-cell wall resistance",
                        "Effective rate constants are uncalibrated; no first-principles or industrial validation"],
        "inputs": inputs, "summary": summary,
    }
    if include_field:
        result["field"] = {"x_mm": (x * 1000).tolist(), "y_um": (y * 1e6).tolist(),
                           "A_mol_m3": concentration.tolist(), "P_mol_m3": product.tolist(),
                           "velocity_m_s": u.tolist(), "wall_A_mol_m3": wall_A.tolist(),
                           "wall_P_mol_m3": wall_P.tolist(), "wall_flux_mol_m2_s": wall_flux.tolist(),
                           "wall_current_A_m2": wall_current.tolist()}
    return result


def _write_json(path, data):
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def _write_csv(path, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def run_verification(output=OUTPUT):
    output = Path(output); output.mkdir(parents=True, exist_ok=True)
    baseline = solve_transport()
    _write_json(output / "baseline_solution.json", baseline)
    grid_results = []
    for nx, ny in [(50, 12), (100, 24), (200, 48), (400, 96), (200, 96), (400, 48), (80, 32)]:
        item = solve_transport(nx=nx, ny=ny, include_field=False)
        grid_results.append({"nx": nx, "ny": ny, "cells": nx * ny, **item["summary"]})
    fine = grid_results[3]["conversion_pct"]
    for row in grid_results:
        row["conversion_delta_from_finest_pp"] = row["conversion_pct"] - fine
    _write_csv(output / "grid_convergence.csv", grid_results)
    sweep = []
    for q in [150.0, 450.0, 800.0, 1200.0]:
        for eta in [0.25, 0.45, 0.65]:
            item = solve_transport(q, eta, nx=80, ny=32, include_field=False)
            sweep.append({"flow_rate_uL_min": q, "overpotential_V": eta, **item["summary"]})
    _write_csv(output / "parameter_sweep.csv", sweep)
    screening_grid_audit, screening_solutions = [], []
    for q, eta in [(100, .2), (100, .75), (1500, .2), (1500, .75), (800, .475)]:
        coarse = solve_transport(q, eta, nx=80, ny=32, include_field=False)["summary"]
        refined = solve_transport(q, eta, nx=400, ny=96, include_field=False)["summary"]
        screening_solutions.extend([coarse, refined])
        screening_grid_audit.append({
            "flow_rate_uL_min": q, "overpotential_V": eta,
            "screening_nx": 80, "screening_ny": 32, "reference_nx": 400, "reference_ny": 96,
            "screening_conversion_pct": coarse["conversion_pct"],
            "reference_conversion_pct": refined["conversion_pct"],
            "screening_minus_reference_pp": coarse["conversion_pct"] - refined["conversion_pct"],
            "screening_current_A": coarse["current_A"], "reference_current_A": refined["current_A"],
            "reference_converged": refined["converged"],
            "reference_material_balance_relative_error": refined["material_balance_relative_error"],
        })
    _write_csv(output / "screening_domain_grid_check.csv", screening_grid_audit)
    limits = []
    cases = [
        ("zero_reaction", dict(forward_rate_m_s=0, reverse_rate_m_s=0)),
        ("zero_diffusivity", dict(diffusivity_m2_s=0)),
        ("equal_feed_equilibrium", dict(inlet_A_mM=25, inlet_P_mM=25, overpotential_V=0)),
        ("reverse_only_feed", dict(inlet_A_mM=0, inlet_P_mM=50, overpotential_V=-0.48)),
        ("no_axial_diffusion_diagnostic", dict(axial_diffusivity_m2_s=0)),
    ]
    for label, params in cases:
        item = solve_transport(nx=80, ny=32, include_field=False, **params)
        limits.append({"case": label, **item["summary"]})
    _write_csv(output / "limiting_cases.csv", limits)
    # Independent analytic plug-flow limit: rapid transverse mixing, D_x=0,
    # irreversible wall loss. The continuum solution is C_out/C_in=exp(-k A/Q).
    plug = []
    kf, q = 1e-6, 450e-9 / 60
    exact = (1 - math.exp(-kf * 0.06 * 0.012 / q)) * 100
    for nx in [40, 80, 160, 320]:
        item = solve_transport(nx=nx, ny=8, velocity_profile="plug", diffusivity_m2_s=1e-4,
                               axial_diffusivity_m2_s=0, forward_rate_m_s=kf, reverse_rate_m_s=0,
                               include_field=False)
        plug.append({"nx": nx, "ny": 8, "analytic_conversion_pct": exact,
                     "error_vs_analytic_pp": item["summary"]["conversion_pct"] - exact, **item["summary"]})
    _write_csv(output / "analytic_plug_limit.csv", plug)
    all_summaries = [baseline["summary"]] + grid_results + sweep + limits + plug + screening_solutions
    verification = {
        "schema_version": "electratwin.transport.verification.v1",
        "evidence_class": "numerical_verification_not_experimental_validation",
        "runtime": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__,
                    "hardware": "CPU; one numerical thread requested by launch environment"},
        "baseline": baseline["summary"], "grid_cases": len(grid_results), "parameter_sweep_cases": len(sweep),
        "limiting_cases": len(limits), "analytic_limit_cases": len(plug),
        "screening_domain_cases": len(screening_grid_audit),
        "screening_domain_maximum_abs_conversion_delta_pp": max(abs(r["screening_minus_reference_pp"]) for r in screening_grid_audit),
        "screening_domain_all_reference_checks_passed": all(r["reference_converged"] for r in screening_grid_audit),
        "total_transport_solves": len(all_summaries),
        "all_finite_volume_solutions_passed_numerical_checks": all(r["converged"] for r in all_summaries),
        "maximum_material_balance_relative_error": max(r["material_balance_relative_error"] for r in all_summaries),
        "maximum_charge_balance_relative_error": max(r["charge_balance_relative_error"] for r in all_summaries),
        "finest_grid_conversion_pct": fine,
        "joint_grid_conversion_pct": [r["conversion_pct"] for r in grid_results[:4]],
        "last_joint_refinement_delta_pp": abs(grid_results[3]["conversion_pct"] - grid_results[2]["conversion_pct"]),
        "screening_grid_error_vs_finest_pp": grid_results[-1]["conversion_delta_from_finest_pp"],
        "analytic_plug_limit_conversion_pct": exact, "finest_plug_error_pp": plug[-1]["error_vs_analytic_pp"],
        "important_limits": ["Conservative discretization does not validate kinetics or chemistry",
                             "First-order axial upwind numerical diffusion remains and is quantified by refinement",
                             "FE is structurally 100% for a sole forward reaction, so STY-FE optimization degenerates",
                             "Negative current in the reverse-feed limit is retained, not clipped",
                             "No thermal, migration, potential, side-reaction, transient, or hardware simulation"]
    }
    _write_json(output / "verification.json", verification)
    print(json.dumps(verification, ensure_ascii=False, indent=2, allow_nan=False))
    return verification


if __name__ == "__main__":
    run_verification()
