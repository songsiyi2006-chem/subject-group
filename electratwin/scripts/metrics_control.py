"""Audited metrics, finite-pool sequential MC-EHVI and simulation-only SCPI.

All evaluations in this module are numerical demonstrations. No hardware APIs,
network connections, experimental selectivity or product identities are inferred.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import platform
import sys
from importlib.metadata import version

import numpy as np
from scipy.stats import t as student_t
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, Matern

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "electratwin/results/control"
FARADAY = 96485.33212


def _finite(name, value, minimum=None, strictly_positive=False):
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    if minimum is not None and value < minimum:
        raise ValueError(f"{name} must be >= {minimum}")
    if strictly_positive and value <= 0:
        raise ValueError(f"{name} must be positive")
    return value


def green_metrics(*, substrate_mw, coupling_partner_mw, product_mw,
                  flow_rate_mL_min, reactor_volume_mL, inlet_conc_M,
                  conversion_pct, chemoselectivity_pct, current_A,
                  cell_voltage_V, solvent_mass_flow_g_min,
                  electrons_transferred=2, partner_equivalents=1.0,
                  recovery_fraction=1.0, additional_inputs_g_min=None,
                  full_boundary_verified=False, balanced_reaction_verified=False):
    """Preserve unrounded computed values; flag impossible FE and open boundaries.

    Product rate assumes 1 product per converted substrate. Additional inputs
    must be incremental external mass flows, not double-counted solution mass.
    Formal atom economy is withheld without an independently balanced reaction.
    """
    positive = {"substrate_mw": substrate_mw, "coupling_partner_mw": coupling_partner_mw,
                "product_mw": product_mw, "flow_rate_mL_min": flow_rate_mL_min,
                "reactor_volume_mL": reactor_volume_mL, "inlet_conc_M": inlet_conc_M}
    for key, value in positive.items():
        _finite(key, value, strictly_positive=True)
    for key, value in {"conversion_pct": conversion_pct, "chemoselectivity_pct": chemoselectivity_pct,
                       "current_A": current_A, "cell_voltage_V": cell_voltage_V,
                       "solvent_mass_flow_g_min": solvent_mass_flow_g_min,
                       "partner_equivalents": partner_equivalents}.items():
        _finite(key, value, minimum=0)
    if conversion_pct > 100 or chemoselectivity_pct > 100:
        raise ValueError("Conversion and chemoselectivity must be in [0, 100]")
    _finite("recovery_fraction", recovery_fraction, minimum=0)
    if recovery_fraction > 1:
        raise ValueError("recovery_fraction must be <= 1")
    if isinstance(electrons_transferred, bool) or int(electrons_transferred) != electrons_transferred or electrons_transferred <= 0:
        raise ValueError("electrons_transferred must be a positive integer")
    extra = dict(additional_inputs_g_min or {})
    for key, value in extra.items():
        _finite(f"additional_inputs_g_min.{key}", value, minimum=0)
    substrate_mol_min = inlet_conc_M * flow_rate_mL_min / 1000
    product_mol_min = substrate_mol_min * conversion_pct / 100 * chemoselectivity_pct / 100
    product_g_min = product_mol_min * product_mw
    recovered_g_min = product_g_min * recovery_fraction
    partial_inputs = substrate_mol_min * (substrate_mw + partner_equivalents * coupling_partner_mw) + solvent_mass_flow_g_min
    all_inputs = partial_inputs + sum(extra.values())
    flags = ["product_identity_and_molecular_weights_are_user_assumptions"]
    if not balanced_reaction_verified:
        flags.append("atom_economy_unverified_balanced_reaction_missing")
    if not full_boundary_verified:
        flags.append("partial_mass_boundary_not_full_process_PMI")
    fe = product_mol_min / 60 * electrons_transferred * FARADAY / current_A * 100 if current_A else None
    if fe is None:
        flags.append("faradaic_efficiency_undefined_zero_current")
    elif fe > 100 + 1e-6:
        flags.append("faradaic_efficiency_exceeds_100_retained_not_clipped")
    pmi = all_inputs / recovered_g_min if recovered_g_min else None
    sec = cell_voltage_V * current_A / (product_g_min * 60) if product_g_min else None
    if recovered_g_min == 0:
        flags.append("mass_metrics_undefined_zero_recovered_product")
    if product_g_min == 0:
        flags.append("specific_energy_undefined_zero_product")
    if pmi is not None and pmi < 1:
        flags.append("mass_intensity_below_1_check_stoichiometry_or_boundary")
    assumed_ae = 100 * product_mw / (substrate_mw + coupling_partner_mw)
    if assumed_ae > 100:
        flags.append("assumed_1_to_1_atom_economy_exceeds_100")
    return {
        "evidence_role": "calculated_from_declared_inputs_not_experiment",
        "inputs": {**positive, "conversion_pct": conversion_pct, "chemoselectivity_pct": chemoselectivity_pct,
                   "current_A": current_A, "cell_voltage_V": cell_voltage_V,
                   "solvent_mass_flow_g_min": solvent_mass_flow_g_min, "electrons_transferred": electrons_transferred,
                   "partner_equivalents": partner_equivalents, "recovery_fraction": recovery_fraction,
                   "additional_inputs_g_min": extra, "full_boundary_verified": full_boundary_verified,
                   "balanced_reaction_verified": balanced_reaction_verified},
        "product_mol_min": product_mol_min, "product_g_min": product_g_min,
        "recovered_product_g_min": recovered_g_min,
        "assumed_1_to_1_atom_economy_pct": assumed_ae,
        "atom_economy_pct": assumed_ae if balanced_reaction_verified else None,
        "space_time_yield_kg_m3_day": product_g_min * 1.44 / (reactor_volume_mL * 1e-6),
        "partial_input_mass_g_min": partial_inputs, "declared_input_mass_g_min": all_inputs,
        "declared_boundary_mass_intensity": pmi,
        "full_process_PMI": pmi if full_boundary_verified else None,
        "declared_boundary_waste_ratio_assuming_all_nonproduct_is_waste": pmi - 1 if pmi is not None else None,
        "faradaic_efficiency_pct": fe, "electrical_energy_kwh_kg_reactor_product": sec,
        "electrical_energy_kwh_kg_recovered_product": sec / recovery_fraction if sec is not None and recovery_fraction else None,
        "flags": flags,
    }


def pareto_mask(points):
    """Nondomination for finite, maximized objectives; identical ties are retained."""
    p = np.asarray(points, dtype=float)
    if p.ndim != 2 or p.shape[1] != 2 or not np.isfinite(p).all():
        raise ValueError("Expected finite (n, 2) objectives")
    return np.array([not np.any(np.all(p >= row, axis=1) & np.any(p > row, axis=1)) for row in p])


def hypervolume_2d(points, reference=(0.0, -2.0)):
    """Exact union of rectangles above a fixed reference, both objectives maximized."""
    p = np.asarray(points, dtype=float).reshape(-1, 2)
    ref = np.asarray(reference, dtype=float)
    if ref.shape != (2,) or not np.isfinite(ref).all() or not np.isfinite(p).all():
        raise ValueError("Nonfinite points/reference")
    p = p[np.all(p > ref, axis=1)]
    if not len(p):
        return 0.0
    p = p[np.argsort(-p[:, 0], kind="stable")]
    hv, current_y = 0.0, ref[1]
    for x, y in p:
        if y > current_y:
            hv += (x - ref[0]) * (y - current_y)
            current_y = y
    return float(hv)


def hypervolume_improvement(samples, observed, reference=(0.0, -2.0)):
    """Vectorized exact HVI for individual candidate draws against fixed observations."""
    z = np.asarray(samples, dtype=float)
    obs = np.asarray(observed, dtype=float).reshape(-1, 2)
    ref = np.asarray(reference, dtype=float)
    if z.shape[-1] != 2 or not np.isfinite(z).all() or not np.isfinite(obs).all():
        raise ValueError("Expected finite objectives")
    width = np.maximum(z[..., 0] - ref[0], 0)
    gain = width * np.maximum(z[..., 1] - ref[1], 0)
    obs = obs[np.all(obs > ref, axis=1)]
    obs = obs[np.argsort(-obs[:, 0], kind="stable")]
    ylo = ref[1]
    for x, y in obs:
        if y > ylo:
            gain -= np.minimum(width, x - ref[0]) * np.clip(z[..., 1] - ylo, 0, y - ylo)
            ylo = y
    return np.maximum(gain, 0)  # Roundoff guard only for geometrical area.


def sequential_campaign(candidate_x, evaluate, *, budget=15, initial_count=5, seed=0,
                        method="gp_mc_ehvi", mc_draws=256, reference=(0.0, -2.0)):
    """Oracle is called only for selected points; no unseen objectives reach the GP.

    Independent fixed-kernel Matérn-5/2 GPs; empirical observation normalization
    is refit on observed data only. Finite candidate-pool, single-point MC-EHVI.
    This is not noisy/batch EHVI or an experimentally validated controller.
    """
    x = np.asarray(candidate_x, dtype=float)
    if x.ndim != 2 or x.shape[1] != 2 or not np.isfinite(x).all():
        raise ValueError("candidate_x must be finite (n, 2)")
    if len(np.unique(x, axis=0)) != len(x):
        raise ValueError("Duplicate candidate inputs are disallowed")
    if not (2 <= initial_count <= budget <= len(x)):
        raise ValueError("Require 2 <= initial_count <= budget <= pool size")
    if method not in {"gp_mc_ehvi", "random"} or mc_draws < 1:
        raise ValueError("Unsupported method or invalid MC count")
    scale = np.ptp(x, axis=0)
    if np.any(scale <= 0):
        raise ValueError("Both input dimensions must vary")
    normalized_x = (x - x.min(axis=0)) / scale
    rng = np.random.default_rng(seed)
    initial = rng.choice(len(x), size=initial_count, replace=False).tolist()
    acquisition_rng = np.random.default_rng(seed + 100000)
    selected, observed, records = [], [], []
    for step in range(budget):
        remaining = np.array([i for i in range(len(x)) if i not in selected])
        score = None
        if step < initial_count:
            index = initial[step]
        elif method == "random":
            index = int(rng.choice(remaining))
        else:
            means, stds = [], []
            for objective in range(2):
                kernel = ConstantKernel(1.0, constant_value_bounds="fixed") * Matern(
                    length_scale=[0.3, 0.3], length_scale_bounds="fixed", nu=2.5)
                gp = GaussianProcessRegressor(kernel=kernel, alpha=1e-8, optimizer=None, normalize_y=True)
                gp.fit(normalized_x[selected], np.asarray(observed)[:, objective])
                mu, sd = gp.predict(normalized_x[remaining], return_std=True)
                means.append(mu)
                stds.append(sd)
            # Common random normal draws across candidates reduce comparison noise.
            normal = acquisition_rng.standard_normal((1, mc_draws, 2))
            draws = np.stack(means, axis=-1)[:, None, :] + np.stack(stds, axis=-1)[:, None, :] * normal
            improvement = hypervolume_improvement(draws, observed, reference)
            # Physical-domain mask on posterior draws: nonnegative STY, positive SEC.
            physical = (draws[..., 0] >= 0) & (draws[..., 1] < 0)
            expected = np.mean(np.where(physical, improvement, 0.0), axis=1)
            best = int(np.argmax(expected))
            index, score = int(remaining[best]), float(expected[best])
        value = np.asarray(evaluate(index), dtype=float)
        if value.shape != (2,) or not np.isfinite(value).all():
            raise ValueError("Oracle returned invalid objective vector")
        selected.append(index)
        observed.append(value.tolist())
        records.append({"seed": seed, "method": method, "evaluation": step + 1,
                        "candidate_index": index, "flow_rate_uL_min": float(x[index, 0]),
                        "overpotential_V": float(x[index, 1]), "objective_STY_div_50000": float(value[0]),
                        "objective_negative_SEC": float(value[1]), "acquisition_MC_EHVI": score,
                        "hypervolume": hypervolume_2d(observed, reference),
                        "phase": "initial_shared" if step < initial_count else "sequential"})
    return records


class SimulationOnlySCPI:
    """SCPI-like toy state machine. No vendor compatibility or actual IO claims."""
    def __init__(self, resource="SIMULATOR::ELECTRATWIN"):
        if resource != "SIMULATOR::ELECTRATWIN":
            raise ValueError("Only the explicit local simulator resource is accepted")
        self.resource = resource
        self.connected = False
        self.output_enabled = False
        self.current_A = 0.0
        self.voltage_limit_V = 10.0
        self.fault = None
        self.transcript = []

    def connect(self):
        self.connected = True
        return {"status": "SIMULATOR_CONNECTED", "physical_connection": False}

    def disconnect(self):
        self.output_enabled = False
        self.current_A = 0.0
        self.connected = False

    def _reply(self, cmd, status, value=None):
        result = {"command": cmd, "status": status, "value": value,
                  "connected": self.connected, "output_enabled": self.output_enabled,
                  "current_setpoint_A": self.current_A, "voltage_limit_V": self.voltage_limit_V,
                  "fault": self.fault, "evidence_role": "simulation_only"}
        self.transcript.append(result)
        return result

    def _error(self, cmd, code):
        self.output_enabled = False
        self.fault = code
        return self._reply(cmd, "ERROR", code)

    def _compliance(self, cmd):
        if self.output_enabled and 1.95 + self.current_A * 22.4 > self.voltage_limit_V:
            return self._error(cmd, "SIMULATED_VOLTAGE_COMPLIANCE_TRIP")
        return self._reply(cmd, "OK")

    def send_command(self, command):
        cmd = str(command).strip().upper()
        if not self.connected:
            return self._error(cmd, "NOT_CONNECTED")
        if cmd == "*IDN?":
            return self._reply(cmd, "OK", "ELECTRATWIN,SIMULATOR,NO-PHYSICAL-DEVICE,1.0")
        if cmd == ":OUTP OFF":
            self.output_enabled = False
            return self._reply(cmd, "OK")
        if cmd == "*CLS":
            if self.output_enabled or self.current_A != 0:
                return self._error(cmd, "RESET_REQUIRES_OUTPUT_OFF_AND_ZERO_CURRENT")
            self.fault = None
            return self._reply(cmd, "OK")
        if cmd == ":SYST:ERR?":
            return self._reply(cmd, "OK", self.fault or "NO_ERROR")
        if cmd == ":SOUR:FUNC CURR":
            return self._reply(cmd, "OK", "CURRENT")
        if cmd == ":OUTP ON":
            if self.fault is not None:
                return self._error(cmd, "FAULT_LATCHED_RESET_REQUIRED")
            self.output_enabled = True
            return self._compliance(cmd)
        if cmd in {":MEAS:VOLT?", ":MEAS:CURR?"}:
            if not self.output_enabled:
                return self._error(cmd, "OUTPUT_DISABLED_NO_MEASUREMENT")
            value = 1.95 + 22.4 * self.current_A if cmd == ":MEAS:VOLT?" else self.current_A
            return self._reply(cmd, "OK", value)
        parts = cmd.split()
        if len(parts) == 2 and parts[0] in {":SOUR:CURR", ":SENS:VOLT:PROT"}:
            try:
                value = float(parts[1])
            except ValueError:
                return self._error(cmd, "INVALID_NUMERIC_PARAMETER")
            low, high = (0.0, 0.1) if parts[0] == ":SOUR:CURR" else (0.01, 10.0)
            if not math.isfinite(value) or not low <= value <= high:
                return self._error(cmd, "PARAMETER_OUT_OF_SIMULATOR_RANGE")
            if parts[0] == ":SOUR:CURR":
                self.current_A = value
            else:
                self.voltage_limit_V = value
            return self._compliance(cmd)
        return self._error(cmd, "UNKNOWN_OR_UNSUPPORTED_COMMAND")


def _save_json(name, data):
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / name).write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def _save_csv(name, rows):
    OUT.mkdir(parents=True, exist_ok=True)
    with (OUT / name).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def metric_sensitivity_demo():
    base = dict(substrate_mw=223.27, coupling_partner_mw=110.18, product_mw=331.43,
                flow_rate_mL_min=0.45, reactor_volume_mL=0.216, inlet_conc_M=0.05,
                conversion_pct=20, chemoselectivity_pct=100, current_A=0.014472799818,
                cell_voltage_V=2.5, solvent_mass_flow_g_min=0.45 * 0.786)
    cases = {"hypothetical_single_reaction": green_metrics(**base),
             "one_electron_current_two_electron_reporting": green_metrics(**{**base, "current_A": base["current_A"] / 2}),
             "assumed_added_inputs_and_80pct_recovery": green_metrics(**base, recovery_fraction=0.8,
                 additional_inputs_g_min={"electrolyte": 0.005, "workup_solvent": 0.5, "wash_water": 0.2}),
             "zero_product": green_metrics(**{**base, "conversion_pct": 0})}
    _save_json("metric_boundary_sensitivity.json", {
        "evidence_role": "hypothetical_arithmetic_sensitivity_not_measurement",
        "cases": cases,
        "source_corrections": ["No FE clipping", "No fabricated 9999/999 metric sentinels",
                               "No ACS conformance asserted for an incomplete inventory",
                               "No formal atom economy without a verified reaction equation",
                               "Electrical SEC excludes pumps, thermal control and downstream processing"]})
    return cases


def scpi_demo():
    sim = SimulationOnlySCPI()
    sim.send_command(":OUTP ON")
    sim.connect()
    for command in ["*IDN?", "*CLS", ":SOUR:FUNC CURR", ":SOUR:CURR 0.02", ":OUTP ON",
                    ":MEAS:VOLT?", ":SENS:VOLT:PROT 2.0", ":OUTP ON", ":SOUR:CURR 0", "*CLS",
                    ":SENS:VOLT:PROT 10", ":SOUR:CURR NaN", ":SOUR:CURR 0", "*CLS",
                    ":SOUR:VOLT 3", ":SOUR:CURR 0", "*CLS", ":SOUR:CURR 0.01", ":OUTP ON",
                    ":MEAS:CURR?", ":OUTP OFF"]:
        sim.send_command(command)
    sim.disconnect()
    _save_json("scpi_simulator_trace.json", {"physical_connections": 0, "simulation_only": True,
        "vendor_compatibility_verified": False, "scope": "SCPI-like state-machine tests, not a real instrument driver",
        "trace": sim.transcript, "final_state": {"connected": sim.connected, "output_enabled": sim.output_enabled,
                                                 "current_A": sim.current_A}})


def evaluate_pool(nx=80, ny=32):
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from transport_reviewed import solve_transport
    rows = []
    for flow in np.linspace(100, 1500, 9):
        for eta in np.linspace(0.2, 0.75, 9):
            result = solve_transport(flow_rate_uL_min=float(flow), overpotential_V=float(eta), nx=nx, ny=ny,
                                     include_field=False)
            s = result["summary"]
            if not s["converged"] or not s["nonnegative"]:
                raise RuntimeError("Numerical transport evaluation failed feasibility checks")
            current = float(s["current_A"])
            voltage = 1.85 + eta + 18.5 * current  # Source assumptions, no fitted voltage data.
            metrics = green_metrics(substrate_mw=223.27, coupling_partner_mw=110.18, product_mw=331.43,
                flow_rate_mL_min=float(flow / 1000), reactor_volume_mL=float(s["reactor_volume_m3"] * 1e6), inlet_conc_M=0.05,
                conversion_pct=float(s["conversion_pct"]), chemoselectivity_pct=100,
                current_A=current, cell_voltage_V=float(voltage), solvent_mass_flow_g_min=float(flow / 1000 * 0.786))
            rows.append({"candidate_index": len(rows), "flow_rate_uL_min": float(flow), "overpotential_V": float(eta),
                "conversion_pct": float(s["conversion_pct"]), "current_A": current, "cell_voltage_assumed_V": float(voltage),
                "STY_assumed_product_kg_m3_day": metrics["space_time_yield_kg_m3_day"],
                "SEC_assumed_product_kwh_kg": metrics["electrical_energy_kwh_kg_reactor_product"],
                "FE_pct": metrics["faradaic_efficiency_pct"],
                "partial_mass_intensity": metrics["declared_boundary_mass_intensity"],
                "material_balance_relative_error": float(s["material_balance_relative_error"]),
                "charge_balance_relative_error": float(s["charge_balance_relative_error"]),
                "nx": nx, "ny": ny, "evidence_role": "uncalibrated_transport_simulation_with_assumed_product_mass"})
    objectives = np.array([[r["STY_assumed_product_kg_m3_day"] / 50000, -r["SEC_assumed_product_kwh_kg"]] for r in rows])
    pm = pareto_mask(objectives)
    # Round only for diagnosing physical degeneracy, never for reporting FE or optimization.
    sty_fe = np.array([[r["STY_assumed_product_kg_m3_day"], 100.0] for r in rows])
    for i, row in enumerate(rows):
        row["pareto_STY_negative_SEC"] = bool(pm[i])
        row["pareto_STY_constant_FE_model_identity"] = bool(pareto_mask(sty_fe)[i])
    _save_csv("candidate_pool.csv", rows)
    return rows, objectives


def baseline_metrics():
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from transport_reviewed import solve_transport
    result = solve_transport(flow_rate_uL_min=450.0, overpotential_V=0.48,
                             nx=100, ny=48, include_field=False)
    s = result["summary"]
    current = float(s["current_A"])
    voltage = 1.85 + 0.48 + 18.5 * current
    metrics = green_metrics(substrate_mw=223.27, coupling_partner_mw=110.18, product_mw=331.43,
        flow_rate_mL_min=0.45, reactor_volume_mL=float(s["reactor_volume_m3"] * 1e6), inlet_conc_M=0.05,
        conversion_pct=float(s["conversion_pct"]), chemoselectivity_pct=100, current_A=current,
        cell_voltage_V=voltage, solvent_mass_flow_g_min=0.45 * 0.786)
    _save_json("baseline_metrics.json", {"evidence_role": "uncalibrated_model_with_assumed_product_identity",
        "additional_PDE_evaluations_outside_81_point_pool": 1,
        "conditions": {"flow_rate_uL_min": 450.0, "overpotential_V": 0.48, "nx": 100, "ny": 48},
        "transport_summary": s, "assumed_voltage_formula": "1.85 + eta + 18.5 * current_A",
        "cell_voltage_assumed_V": voltage, "metrics": metrics})


def run_benchmarks(*, seeds=8, budget=15, mc_draws=256, nx=80, ny=32):
    if seeds < 1:
        raise ValueError("seeds must be positive")
    rows, objectives = evaluate_pool(nx, ny)
    candidate_x = np.array([[r["flow_rate_uL_min"], r["overpotential_V"]] for r in rows])
    full_hv = hypervolume_2d(objectives)
    records, summaries, campaign_fronts = [], [], []
    for seed in range(seeds):
        pair = {}
        for method in ["gp_mc_ehvi", "random"]:
            run = sequential_campaign(candidate_x, lambda i: objectives[i], budget=budget, initial_count=5,
                                      seed=seed, method=method, mc_draws=mc_draws)
            records.extend(run)
            keep = pareto_mask([[r["objective_STY_div_50000"], r["objective_negative_SEC"]] for r in run])
            campaign_fronts.extend([r for r, include in zip(run, keep) if include])
            pair[method] = run[-1]["hypervolume"]
        summaries.append({"seed": seed, "budget_per_method": budget, "initial_count": 5,
                          "EHVI_final_HV": pair["gp_mc_ehvi"], "random_final_HV": pair["random"],
                          "paired_difference": pair["gp_mc_ehvi"] - pair["random"],
                          "EHVI_full_pool_HV_fraction": pair["gp_mc_ehvi"] / full_hv,
                          "random_full_pool_HV_fraction": pair["random"] / full_hv})
    _save_csv("sequential_campaigns.csv", records)
    _save_csv("campaign_pareto_records.csv", campaign_fronts)
    _save_csv("paired_budget_comparison.csv", summaries)
    d = np.array([r["paired_difference"] for r in summaries])
    eh = np.array([r["EHVI_full_pool_HV_fraction"] for r in summaries])
    ra = np.array([r["random_full_pool_HV_fraction"] for r in summaries])
    fe = np.array([r["FE_pct"] for r in rows])
    seed_interval = None
    if seeds > 1:
        margin = float(student_t.ppf(0.975, seeds - 1) * d.std(ddof=1) / np.sqrt(seeds))
        seed_interval = [float(d.mean() - margin), float(d.mean() + margin)]
    _save_json("optimization_summary.json", {
        "evidence_role": "executed_numerical_benchmark_not_wet_lab_or_industrial_control",
        "pool_design": {"flow_uL_min": [100, 1500], "overpotential_V": [0.2, 0.75], "shape": [9, 9],
                        "unique_PDE_evaluations": len(rows), "nx": nx, "ny": ny},
        "budget": {"seeds": seeds, "methods": 2, "per_method_including_initial": budget,
                   "shared_initial_per_seed": 5, "sequential_per_run": budget - 5,
                   "total_recorded_evaluation_uses": len(records),
                   "offline_precomputation_for_reference_excluded_from_per_method_budget": len(rows)},
        "objectives_maximized": ["STY / 50000", "-electrical_SEC_kWh_per_kg"],
        "fixed_reference": [0.0, -2.0], "full_pool_hypervolume": full_hv,
        "full_pool_nondominated_count": int(pareto_mask(objectives).sum()),
        "FE_range_pct": [float(fe.min()), float(fe.max())],
        "FE_is_fixed_by_single_reaction_charge_identity": True,
        "STY_FE_model_identity_front_count": sum(r["pareto_STY_constant_FE_model_identity"] for r in rows),
        "method": {"name": "single-candidate finite-pool GP MC-EHVI", "posterior_draws_per_candidate": mc_draws,
                   "kernel": "independent normalized Matérn-5/2, fixed length scales [0.3,0.3], amplitude 1",
                   "numerical_nugget": 1e-8, "hyperparameter_optimization": False,
                   "observed_only_outcome_normalization": True, "unseen_objectives_available_to_GP": False,
                   "draw_constraints": "nonnegative STY and positive SEC; no learned process-safety constraints",
                   "experimental_noise_model": None},
        "comparison": {"EHVI_mean_pool_HV_fraction": float(eh.mean()), "random_mean_pool_HV_fraction": float(ra.mean()),
                       "EHVI_min_max_fraction": [float(eh.min()), float(eh.max())],
                       "random_min_max_fraction": [float(ra.min()), float(ra.max())],
                       "mean_paired_HV_difference": float(d.mean()), "paired_HV_difference_sample_sd": float(d.std(ddof=1)) if seeds > 1 else None,
                       "paired_seed_mean_difference_95pct_t_interval": seed_interval,
                       "interval_scope": "Descriptive seed sensitivity under independent approximately normal seed effects; not chemical uncertainty or generalized method superiority",
                       "EHVI_wins": int(np.sum(d > 1e-12)), "random_wins": int(np.sum(d < -1e-12)),
                       "ties": int(np.sum(np.abs(d) <= 1e-12))},
        "selection": "All nondominated evaluated candidates exported in campaign_pareto_records.csv; no arbitrary median declared industrial optimum",
        "limitations": ["FE is structurally 100% because side reactions are absent; it is not a measured selectivity claim",
            "331.43 g/mol product and 223.27/110.18 g/mol reactants are attachment assumptions without verified identities",
            "Cell voltage and solvent density are assumed, not coupled to a potential or thermal PDE",
            "Independent GP posteriors ignore physical objective correlations",
            "Only finite 81-point pool and eight default seeds; superiority cannot be generalized",
            "No real experimental data, hardware connection, DFT, uncertainty calibration or industrial validation",
            "Reference point and STY scale are declared a priori and are not tuned using benchmark outcomes"],
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "transport_source_sha256": hashlib.sha256((Path(__file__).parent / "transport_reviewed.py").read_bytes()).hexdigest(),
        "candidate_pool_sha256": hashlib.sha256((OUT / "candidate_pool.csv").read_bytes()).hexdigest(),
        "runtime_versions": {"python": platform.python_version(), "numpy": np.__version__,
                             "scipy": version("scipy"), "scikit_learn": version("scikit-learn")}})
    return {"unique_PDE_evaluations": len(rows), "campaign_records": len(records), "comparison": {
        "EHVI_mean_HV_fraction": float(eh.mean()), "random_mean_HV_fraction": float(ra.mean()),
        "EHVI_wins": int(np.sum(d > 1e-12)), "random_wins": int(np.sum(d < -1e-12))}}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--skip-optimization", action="store_true")
    parser.add_argument("--seeds", type=int, default=8)
    parser.add_argument("--budget", type=int, default=15)
    parser.add_argument("--mc-draws", type=int, default=256)
    parser.add_argument("--nx", type=int, default=80)
    parser.add_argument("--ny", type=int, default=32)
    args = parser.parse_args()
    metric_sensitivity_demo()
    scpi_demo()
    if not args.skip_optimization:
        baseline_metrics()
        print(json.dumps(run_benchmarks(seeds=args.seeds, budget=args.budget, mc_draws=args.mc_draws,
                                        nx=args.nx, ny=args.ny), indent=2))
    else:
        print("Metric sensitivity and local-only simulator trace written.")


if __name__ == "__main__":
    main()
