"""Run the immutable submitted program, then an explicitly recorded API repair.

No real instrument is accessed: the submitted driver only changes in-memory state.
The child captures model arrays after runpy returns; it does not alter source logic.
"""
from __future__ import annotations
import argparse
import ast
from datetime import datetime, timezone
import difflib
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import re
import runpy
import subprocess
import sys
import time

BASE = Path(__file__).resolve().parents[1]
REPO = BASE.parent


def write_json(path, data):
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sanitize(text):
    for path, replacement in [(str(REPO), "<repository>"), (str(Path(sys.executable).parent.parent), "<chem-ai4s-runtime>")]:
        text = text.replace(path, replacement).replace(path.replace("\\", "/"), replacement)
    text = re.sub(r"[A-Za-z]:[\\/]Users[\\/][^\r\n\"']+", "<local-runtime-path>", text)
    return text


def child(mode):
    source = BASE / "source" / ("electratwin_core.py" if mode == "original" else "electratwin_core_numpy_compat.py")
    values = runpy.run_path(str(source), run_name="__main__")
    write_json(Path("captured_model_arrays.json"), {
        "evidence_role": "executed_source_model_not_experiment",
        "capture_method": "runpy return globals after the unchanged main program completed",
        "pde_result": values["pde_result"],
        "campaign_results": values["campaign_results"],
        "campaign_history": values["controller"].history,
        "rdkit_available": values["HAS_RDKIT"],
        "molecule_structures_processed": 0,
        "physical_instrument_connected": False,
        "driver_output_enabled_after_program": values["driver"].output_enabled,
    })


def audit_saved_results():
    """Recompute conservation/FE diagnostics without rerunning the expensive source PDE."""
    import contextlib
    import io
    import numpy as np

    capture_path = BASE / "results" / "compatibility" / "captured_model_arrays.json"
    captured = json.loads(capture_path.read_text(encoding="utf-8"))
    pde = captured["pde_result"]
    c = np.asarray(pde["spatial_concentration_2D"])
    j = np.asarray(pde["anode_current_profile_x"])
    length, height, width, diffusivity, faraday = 0.06, 0.0003, 0.012, 1.1e-9, 96485.33
    nx, ny = c.shape
    dx, dy = length / (nx - 1), height / (ny - 1)
    q = pde["flow_rate_uL_min"] * 1e-9 / 60
    cin = pde["inlet_conc_mM"]
    y = np.linspace(0, height, ny)
    u = 6 * q / (width * height) * y / height * (1 - y / height)
    q_discrete = width * np.sum(u * dy)
    removed = width * np.sum(u * (cin - c[-1]) * dy)
    current = float(np.trapezoid(j, dx=dx) * width)
    conv = np.maximum(u[1:-1], 1e-7) / dx
    diff = diffusivity / dy**2
    equation_residual = conv[None, :] * (c[1:, 1:-1] - c[:-1, 1:-1]) - diff * (c[1:, 2:] - 2*c[1:, 1:-1] + c[1:, :-2])
    normalized_equation_residual = abs(equation_residual) / ((conv[None, :] + 2*diff)*cin)
    outlet = float(np.sum(u * c[-1]) / np.sum(u))
    residual_boundary = diffusivity * (c[1:, 1] - c[1:, 0]) / dy - j[1:] / faraday
    rows = []
    for record in captured["campaign_history"]:
        mol_product = .05 * record["flow_rate_uL_min"] * 1e-6 * record["conversion_pct"] / 100 * (94 - 5*record["overpotential_V"]) / 100 / 60
        raw_fe = 2 * faraday * mol_product / record["total_current_A"] * 100
        rows.append({**record, "reconstructed_unclipped_FE_pct": raw_fe,
                     "raw_fe_exceeds_100": bool(raw_fe > 100), "selectivity_assumed_pct": 94 - 5*record["overpotential_V"]})
    source = (BASE / "source" / "electratwin_core.py").read_text(encoding="utf-8")
    tree = ast.parse(source)
    imports = sorted({n.names[0].name for n in ast.walk(tree) if isinstance(n, ast.Import)} | {n.module for n in ast.walk(tree) if isinstance(n, ast.ImportFrom)})
    model = runpy.run_path(str(BASE / "source" / "electratwin_core_numpy_compat.py"), run_name="source_audit_module")
    driver = model["SCPIHardwareDriver"]("SIMULATED_RESOURCE_ONLY")
    driver.connect()
    unknown_response = driver.send_command("UNKNOWN_COMMAND")
    excessive_response = driver.send_command(":SOUR:CURR 1000000")
    negative_response = driver.send_command(":SOUR:CURR -100")
    driver.send_command(":OUTP OFF")

    # Replace only the solver instance with a deterministic stub to inspect planner budget semantics.
    class StubSolver:
        L, H, W = length, height, width
        def solve(self, flow, eta):
            return {"total_anodic_current_A": .01, "conversion_pct": 50.}
    budget_runs = {}
    for budget in [0, 1, 7, 12]:
        controller = model["AutonomousFlowController"](StubSolver())
        with contextlib.redirect_stdout(io.StringIO()):
            plan = controller.run_optimization_campaign(budget_iterations=budget)
        budget_runs[str(budget)] = plan["total_campaign_runs"]
    hv = model["HypervolumeCalculator"].compute_2d_hypervolume(np.array([[1., 3.], [2., 2.], [3., 1.], [1., 1.]]), np.array([0., 0.]))
    result = {
        "evidence_role": "independent_arithmetic_and_code_audit_of_source_model",
        "source_capture_sha256": digest(capture_path),
        "source_role": "unvalidated continuum toy model; not experimental data or validated digital twin",
        "physical_instrument_connected": False,
        "single_case": {
            "nx": nx, "ny": ny, "source_reported_conversion_pct": pde["conversion_pct"],
            "unrounded_mixed_cup_conversion_pct": (1-outlet/cin)*100,
            "source_reported_current_A": pde["total_anodic_current_A"], "integrated_current_A_unrounded": current,
            "source_iterations": pde["iterations"], "source_update_norm": pde["final_residual"],
            "discrete_equation_residual_max_mol_m3_s": float(abs(equation_residual).max()),
            "discrete_equation_residual_max_normalized": float(normalized_equation_residual.max()),
            "anode_boundary_residual_max_mol_m2_s_excluding_inlet": float(abs(residual_boundary).max()),
            "imposed_Q_m3_s": q, "discretely_integrated_Q_m3_s": float(q_discrete),
            "flow_quadrature_error_relative": float(q_discrete/q-1),
            "substrate_removed_by_mixed_cup_flux_mol_s": float(removed),
            "substrate_removed_from_current_assuming_source_n1_mol_s": current / faraday,
            "relative_current_mass_flux_mismatch_n1": float((current/faraday-removed)/removed),
            "boundary_stoichiometric_n": 1, "green_metrics_stoichiometric_n": 2,
            "first_inlet_node_contributes_current_A": float(width*dx*j[0]/2),
            "first_inlet_node_current_fraction": float(width*dx*j[0]/2/current),
            "inlet_corner_anode_bc_applied": False,
            "Pe_length_uavg_L_over_D": q/(width*height)*length/diffusivity,
            "Pe_gap_uavg_H_over_D": q/(width*height)*height/diffusivity,
            "Pe_dx_uavg_dx_over_D": q/(width*height)*dx/diffusivity,
            "volume_mL": length*height*width*1e6,
            "residence_time_s": length*height*width/q,
            "notes": ["Source update norm is successive-iterate change, not a reported equation residual.",
                      "Independent equation residual uses the implemented axial-upwind/transverse-diffusion stencil, not the full documented 2D PDE.",
                      "Mass mismatch compares the trapezoidal source current with the source mixed-cup discretization; no clipping is used in this audit.",
                      "Large average axial Pe does not justify dropping axial diffusion everywhere near no-slip walls or inlet corners."]},
        "campaign": {"rows": rows, "evaluations": len(rows), "clipped_fe_rows": sum(r["raw_fe_exceeds_100"] for r in rows),
                     "minimum_unclipped_FE_pct": min(r["reconstructed_unclipped_FE_pct"] for r in rows),
                     "maximum_unclipped_FE_pct": max(r["reconstructed_unclipped_FE_pct"] for r in rows),
                     "fe_reconstruction_basis": "Rounded conversion and rounded current stored by the original controller; matches the inputs used by its green-metrics routine.",
                     "source_pareto_count": captured["campaign_results"]["pareto_solutions_found"],
                     "source_selected_compromise": captured["campaign_results"]["best_industrial_compromise"],
                     "source_initial_hypervolume": captured["campaign_results"]["initial_hypervolume"],
                     "budget_parameter_probe_using_stub_solver": budget_runs,
                     "gp_surrogate_fitted": False, "posterior_uncertainty_computed": False, "ehvi_computed": False, "constraint_model_fitted": False,
                     "design": "Fixed 4 x 3 Cartesian grid; budget value only appears in the printed message.",
                     "hypervolume_hand_case_actual": hv, "hypervolume_hand_case_expected": 6.0,
                     "hypervolume_note": "The tested finite 2D maximization rectangle union is correct; it is not expected hypervolume improvement or evidence of Bayesian optimization."},
        "static_findings": [
            {"id": "PDE_AXIAL_DIFFUSION_ABSENT", "finding": "Docstring includes D*C_xx, but interior code only implements axial upwind convection and transverse diffusion."},
            {"id": "CONVERGENCE_FLAG_UNCONDITIONAL", "finding": "Return dict sets converged=True even if max_iter=2500 is exhausted; residual is iterate change."},
            {"id": "CLIPPED_OBSERVABLES", "finding": "Source clips conversion into [0,99.9] percent and FE into [0,100] percent. Clipping hides rather than diagnoses invalid observables or charge/material inconsistency."},
            {"id": "KINETICS_LINEAR_AT_FIXED_ETA", "finding": "For prescribed eta the boundary law is affine in substrate concentration. No electrostatic potential, product balance or kinetic parameter inference is solved."},
            {"id": "KINETIC_PARAMETERS_UNCALIBRATED", "finding": "j0=0.08 A/m2, alpha=0.5 and prescribed eta are assumed constants without voltammetry or kinetic calibration; this is not a first-principles molecular calculation."},
            {"id": "MULTIPHYSICS_NOT_IMPLEMENTED", "finding": "Temperature, diffusivity and velocity profile are fixed; no Navier-Stokes, energy, electrolyte potential, migration, gas, pressure or fouling equations are solved."},
            {"id": "PARTIAL_MASS_BOUNDARY", "finding": "PMI numerator includes substrate, equimolar coupling partner and feed solvent only; electrolyte, workup, purification, water, cleaning, solvent recovery and losses are absent."},
            {"id": "REACTION_IDENTITY_UNVERIFIED", "finding": "Fixed molecular weights and partner comment do not identify atom-balanced substrate/product structures; atom economy is an assumed 1:1 ratio."},
            {"id": "MOLECULAR_GRAPH_PROFILER_ABSENT", "finding": "RDKit is imported, but no molecular graph, SMILES, descriptor, conformer or reaction mapping is actually processed."},
            {"id": "VOLTAGE_SELECTIVITY_ASSUMPTIONS", "finding": "Ucell=1.85+eta+18.5I and selectivity=94-5eta are assumed algebraic laws without calibration data."},
            {"id": "INDUSTRIAL_STATUS_UNVERIFIED", "finding": "No facility trial, product assay, process qualification, validated instrument I/O or industrial acceptance is supplied."}],
        "driver": {"transport_imports_present": [n for n in imports if n in ["socket", "pyvisa", "serial", "requests"]],
                   "identity_response_role": "hardcoded simulated identity; not queried from equipment",
                   "unknown_command_response": unknown_response, "excessive_setpoint_response": excessive_response,
                   "negative_setpoint_response": negative_response, "voltage_limit_attribute_enforced": False,
                   "real_transport": False, "calibrated_measurement": False,
                   "model_laws_inconsistent": "Controller uses 1.85+eta+18.5I; emulator uses 1.95+22.4I+Gaussian noise."},
    }
    write_json(BASE / "results" / "source_audit.json", result)
    print(json.dumps({"source_audit": "written", "FE_clipped_rows": result["campaign"]["clipped_fe_rows"], "mass_flux_mismatch": result["single_case"]["relative_current_mass_flux_mismatch_n1"]}))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--child", choices=["original", "compatibility"])
    parser.add_argument("--timeout-seconds", type=int, default=600)
    parser.add_argument("--audit-only", action="store_true")
    args = parser.parse_args()
    if args.child:
        child(args.child)
        return
    if args.audit_only:
        audit_saved_results()
        return
    original = BASE / "source" / "electratwin_core.py"
    compatibility = BASE / "source" / "electratwin_core_numpy_compat.py"
    source_text = original.read_text(encoding="utf-8")
    assert source_text.count("np.trapz(") == 1
    repaired_text = source_text.replace("np.trapz(", "np.trapezoid(")
    compatibility.write_text(repaired_text, encoding="utf-8", newline="\n")
    patch = "".join(difflib.unified_diff(source_text.splitlines(True), repaired_text.splitlines(True), fromfile="electratwin_core.py", tofile="electratwin_core_numpy_compat.py"))
    (BASE / "source" / "numpy_compatibility.patch").write_text(patch, encoding="utf-8", newline="\n")
    record = {
        "specification_sha256": digest(BASE / "source" / "specification.md"),
        "original_python_sha256": digest(original),
        "compatibility_python_sha256": digest(compatibility),
        "extraction": "First outer Python fence; line endings normalized to LF and final newline appended. No source logic edits.",
        "compatibility_edit": "One call-name replacement: np.trapz to np.trapezoid; no numerical or scientific model corrections.",
        "numpy_version": importlib.metadata.version("numpy"),
        "python_version": sys.version.split()[0],
        "author_and_laboratory_claims": "Copied source statements; affiliation, authorship and industrial validation are not verified or endorsed.",
        "real_hardware_access": False,
        "real_experiment": False,
    }
    write_json(BASE / "source" / "source_record.json", record)
    for mode in ["original", "compatibility"]:
        target = BASE / "results" / mode
        target.mkdir(parents=True, exist_ok=True)
        env = os.environ.copy()
        env.update({"MPLBACKEND": "Agg", "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "PYTHONIOENCODING": "utf-8"})
        start = time.perf_counter()
        utc = datetime.now(timezone.utc).isoformat()
        timed_out = False
        try:
            proc = subprocess.run([sys.executable, str(Path(__file__).resolve()), "--child", mode], cwd=target, env=env, capture_output=True, timeout=args.timeout_seconds)
            stdout, stderr, code = proc.stdout, proc.stderr, proc.returncode
        except subprocess.TimeoutExpired as exc:
            stdout, stderr, code = exc.stdout or b"", exc.stderr or b"", None
            timed_out = True
        elapsed = time.perf_counter() - start
        (target / "stdout.log").write_text(sanitize(stdout.decode("utf-8", "replace")), encoding="utf-8", newline="\n")
        (target / "stderr.log").write_text(sanitize(stderr.decode("utf-8", "replace")), encoding="utf-8", newline="\n")
        result = {"mode": mode, "started_utc": utc, "elapsed_seconds": elapsed, "exit_code": code, "timeout_seconds": args.timeout_seconds, "timed_out": timed_out,
                  "source_sha256": digest(original if mode == "original" else compatibility), "physical_instrument_connected": False,
                  "source_output_claims_validated": False, "output_scope": "Source failure record or compatibility-only source-model execution; no corrected-physics claim."}
        write_json(target / "execution.json", result)
        print(json.dumps(result))
    if (BASE / "results" / "compatibility" / "captured_model_arrays.json").exists():
        audit_saved_results()


if __name__ == "__main__":
    main()
