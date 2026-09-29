"""Read-only, fixed-commit history extraction and arithmetic; never launches chemistry.

Run with an existing Git object store:
    python reproduce_audit.py --repository PATH_TO_PINCER_REPOSITORY
The default output is this file's parent directory. Only selected public values,
small CSV files, and short native-log excerpts are retained. Full original-byte
hashes and JSON/CSV/line selectors bind each extract to its upstream source.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import io
import json
import re
import subprocess
from collections import Counter
from pathlib import Path

COMMIT = "2edffb123791bd61acbfdeff763f603ee50c2287"
REPOSITORY = "https://github.com/songsiyi2006-chem/pincer-catmech-ai"


def digest(data):
    return hashlib.sha256(data).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    entries, cache, values = {}, {}, {}

    def git(*parts):
        return subprocess.check_output(
            ["git", "-c", "core.longpaths=true", "-C", str(args.repository), *parts],
            timeout=90,
        )

    assert git("rev-parse", COMMIT).decode().strip() == COMMIT

    def source(sid, path):
        if path not in cache:
            cache[path] = git("show", f"{COMMIT}:{path}")
        raw = cache[path]
        entries[sid] = {
            "source_id": sid, "repository": REPOSITORY, "commit": COMMIT,
            "path": path, "sha256": digest(raw), "bytes": len(raw),
            "git_blob_sha1": git("rev-parse", f"{COMMIT}:{path}").decode().strip(),
            "url": f"{REPOSITORY}/blob/{COMMIT}/{path}",
        }
        return raw

    def selected(sid, path, keys, list_fields=None):
        obj = json.loads(source(sid, path))
        data = {"/" + k: obj[k] for k in keys}
        for key, fields in (list_fields or {}).items():
            for i, row in enumerate(obj[key]):
                for field in fields:
                    if field in row:
                        data[f"/{key}/{i}/{field}"] = row[field]
        target = out / f"{sid}_excerpt.json"
        target.write_text(json.dumps({"source": entries[sid], "selected_values": data}, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        entries[sid]["local_excerpt"] = target.name
        entries[sid]["excerpt_sha256"] = digest(target.read_bytes())
        entries[sid]["selectors"] = list(data)
        values[sid] = obj
        return obj

    def csv_source(sid, path):
        raw = source(sid, path)
        target = out / f"{sid}_{Path(path).name}"
        target.write_bytes(raw)
        entries[sid]["local_exact_copy"] = target.name
        rows = list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig"))))
        entries[sid]["data_rows"] = len(rows)
        entries[sid]["columns"] = list(rows[0]) if rows else []
        values[sid] = rows
        return rows

    def evidence(sid, selector):
        return {k: entries[sid][k] for k in ["source_id", "commit", "path", "sha256"]} | {"selector": selector}

    def table(name, rows):
        target = out / name
        with target.open("w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0]))
            writer.writeheader()
            for row in rows:
                writer.writerow({k: json.dumps(v, separators=(",", ":")) if isinstance(v, (dict, list)) else v for k, v in row.items()})
        return {"file": name, "rows": len(rows), "sha256": digest(target.read_bytes())}

    status = selected("C01", "data/phase3/PHASE3_STATUS.json", ["scientific_status", "design_count", "temperature_base_grid_per_design", "readiness_rows", "accepted_physical_18_step_kinetic_predictions", "conditions", "interpretation"])
    spin = selected("C02", "data/phase3/spin/summary.json", ["run_complete", "planned_target_states", "physical_state_attempts", "states_converged", "states_failed_or_missing", "missing_geometry_slots", "s2_and_identity_eligible_states", "computed_vertical_gaps", "eligible_diagnostic_gaps", "mecp_accepted", "validated_catalyst_gaps", "interpretation"], {"target_protocol_runs": ["method", "basis", "grid_radial", "grid_spherical", "scf_algorithm", "states_converged", "target_states", "sha256"]})
    matrix = selected("C03", "data/phase3/spin/pincer_vertical_diagnostic_003/summary.json", ["method", "basis", "grid_radial", "grid_spherical", "scf_algorithm", "e_convergence", "d_convergence", "target_states", "states_converged", "states_failed_or_missing"], {"state_records": ["catalyst_id", "metal", "multiplicity", "charge", "status", "energy_hartree", "spin_squared", "spin_contamination", "spin_contamination_flag", "intended_catalyst_identity_pass", "wall_seconds", "peak_rss_mib", "geometry_sha256", "source_xyz_sha256", "evidence"]})
    states = csv_source("C04", "data/phase3/spin/pincer_vertical_diagnostic_003/spin_states.csv")
    gaps = csv_source("C05", "data/phase3/spin/pincer_vertical_diagnostic_003/vertical_gaps.csv")
    mecp = selected("C06", "data/phase3/spin/pincer_mecp_from_matrix_001/summary.json", ["status", "quantum_pair_evaluations", "quantum_state_attempts", "first_order_crossing", "minimum_verified", "validated_catalyst_label", "evidence"])
    egnn = selected("C07", "data/phase3/egnn/summary.json", ["status", "torch_version", "device", "threads", "configuration", "parameter_count", "elapsed_seconds", "best_validation_epoch", "completed_or_partial_epochs", "split_policy", "metrics", "heads", "eligible_structures", "rejected_records", "auxiliary_pairs", "composition_groups", "families", "symmetry", "predictive_assessment", "split", "scope", "limitations", "model_sha256"])
    selected("C08", "data/phase3/egnn/dataset_manifest.json", ["eligible_structures", "rejected_records", "auxiliary_pairs", "composition_groups", "reference_policy", "scientific_scope", "unsupported_targets", "dataset_files_sha256"])
    assoc = csv_source("C09", "data/phase3/solvation/association_energies.csv")
    thermo = csv_source("C10", "data/phase3/solvation/association_thermochemistry.csv")
    many = csv_source("C11", "data/phase3/solvation/many_body_nonadditivity.csv")
    selected("C12", "data/phase3/solvation/summary.json", ["protocol", "scope", "n_solvent", "orientation_seeds_per_size", "program_version", "converged_optimizations", "limitation"])
    wire = csv_source("C13", "data/phase3/proton_wire/attempts.csv")
    selected("C14", "data/phase3/proton_wire/summary.json", ["status", "method", "solvent", "thermodynamic_baseline_temperature_K", "electronic_smearing_temperature_K", "tBuOK_total_base_equivalent_baseline", "total_base_speciation_resolved", "potassium_explicit", "tert_butoxide_explicit", "scope", "accepted_TS_count", "activation_free_energy_barrier_count", "elapsed_seconds", "continuation_count", "distinct_seed_attempt_count"], {"attempts": ["status", "accepted_TS", "failure_reasons", "actual_NEB_evaluation_reached", "last_reported_NEB_step", "neb_converged", "accepted_endpoint_minima", "imaginary_frequency_count", "free_energy_barrier_computed"]})
    kin = selected("C15", "data/phase3/kinetics/software_verification/summary.json", ["evidence_kind", "grid_count", "timepoints_per_grid", "temperature_K", "initial_free_base_equiv", "maximum_jacobian_central_difference_error", "maximum_metal_error_M", "maximum_base_equivalent_error_M", "maximum_metal_sensitivity_residual_M", "physical_catalyst_predictions", "limitations"])
    phase4 = selected("C16", "data/phase4/PHASE4_STATUS.json", ["scientific_status", "actual_local_quantum_attempts", "scf_converged", "spin_quality_numeric_pass", "production_spin_gaps_validated", "new_TS_accepted", "new_MECP_accepted", "physical_kinetic_predictions", "thermochemical_reanalysis_rows", "baseline_rows_reconstructed", "baseline_max_difference_kcal_mol", "public_crystal_independent_molecules", "public_crystal_atoms_each", "hpc_inputs_prepared", "hpc_jobs_submitted", "quantum_attempt_count_definition", "psi4_diagnostic_attempts", "xtb_precursor_native_calls", "xtb_precursor_minimum_checked_endpoints", "xtb_precursor_distinct_basins_established", "xtb_precursor_mapped_all_atom_RMSD_A", "xtb_precheck_only_run_native_calls", "new_experimental_observations", "synthetic_math_audit"])
    p4jobs = selected("C17", "data/phase4/local_pilot_001/summary.json", ["catalyst_id", "scope", "geometry_optimized_by_this_run", "production_spin_gap_validated", "accepted_MECP_count", "scf_converged_count", "spin_quality_pass_count"], {"jobs": ["name", "basis", "method", "multiplicity", "status", "accepted_chemical_label", "energy_hartree", "spin_squared", "expected_spin_squared", "spin_contamination", "spin_quality_pass", "max_gradient_hartree_bohr", "gradient_kind", "scf_converged", "worker_wall_seconds", "peak_rss_mib"]})
    sensitivity = selected("C18", "data/phase4/solvation_sensitivity_001/summary.json", ["rows", "new_quantum_jobs", "baseline_rows_checked", "baseline_max_error_kcal_mol", "baseline_temperature_K", "cutoff_sensitivity_at_baseline", "scientific_status", "assumptions", "physical_kinetics_validated", "bulk_speciation_validated"])
    selected("C19", "data/phase4/public_precursor_relaxation_002/summary.json", ["method", "solvent", "solvation_state", "ground_spin_determined", "validated_solution_mechanism", "new_transition_state_count", "new_MECP_count", "scope", "optimizations_converged", "mapped_structures_retained", "verified_model_minima", "native_calculations_started", "source_cif_sha256"], {"structures": ["id", "status", "optimization_converged", "minimum_certified", "accepted_chemical_label", "native_calculations_started", "fresh_max_force_eV_A", "electronic_energy_eV", "full_internal_mode_count", "imaginary_count", "lowest_frequency_cm1", "hessian_minimum_pass"]})
    xtb = selected("C20", "projects/cu-np-electroreduction/results/pilot_summary.json", ["scope", "native_calls", "completed_native_calls", "initial_minimum_rejected", "fixed_geometry_charging_energies", "hamiltonian_spread_eV", "predeclared_spread_threshold_eV", "method_sensitivity_stop", "interpretation"], {"rows": ["id", "settings", "energy_hartree", "status", "frequency_count", "minimum_vibrational_frequency_cm1", "significant_imaginary_modes", "minimum_accepted"]})
    dft = selected("C21", "projects/cu-np-operando/data/dft_pilot/summary.json", ["scope", "hartree_to_eV", "runs", "energy_pairs", "completed_quantum_runs", "failed_or_terminated_runs", "preflight_rejections_not_counted_as_quantum_runs", "method_spread_eV", "limitations"])
    selected("C22", "projects/cu-np-operando/data/dft_pilot/geometry_provenance.json", ["molecule_name", "formula", "smiles", "atom_count", "neutral_electron_count", "neutral_charge_multiplicity", "anion_charge_multiplicity", "coordinate_sha256", "source_repository_commit", "source_relative_path", "prior_generator", "prior_minimum_validation", "coordinate_origin", "new_protocol"])
    selected("C23", "projects/cu-np-operando/data/dft_pilot/memory_preflight_rejection.json", ["requested_job", "status", "native_psi4_process_launched", "reason"])
    license_path = "data/phase4/public_structure/import_v002/THIRD_PARTY_LICENSE.md"
    license_raw = source("C24", license_path)
    (out / "C24_THIRD_PARTY_LICENSE.md").write_bytes(license_raw)
    entries["C24"]["local_exact_copy"] = "C24_THIRD_PARTY_LICENSE.md"

    status_counts = dict(Counter(r["status"] for r in states))
    assert status_counts == {"converged": 8, "timeout": 40, "missing_source_geometry": 6}
    eligible = [r for r in matrix["state_records"] if r["status"] == "converged" and not r.get("spin_contamination_flag", True) and r.get("intended_catalyst_identity_pass", False)]
    assert len(eligible) == spin["s2_and_identity_eligible_states"] == 6
    assert sum(bool(r["gap_high_minus_low_hartree"]) for r in gaps) == 1
    checks = {"spin_status_counts_recomputed": status_counts, "eligible_states_recomputed": len(eligible), "phase3_mecp_not_launched": mecp["quantum_state_attempts"] == 0, "egnn_split_sum": sum(v["pairs"] for v in egnn["split"].values()), "solvation_csv_row_count": len(assoc), "association_thermochemistry_row_count": len(thermo), "wire_records": len(wire), "phase4_driver_count_identity": phase4["actual_local_quantum_attempts"] == phase4["psi4_diagnostic_attempts"] + phase4["xtb_precursor_native_calls"], "new_quantum_calculations_in_this_audit": 0}
    tables = []
    gate_rows = []
    for label, key in [("Target slots", "planned_target_states"), ("Physical state attempts", "physical_state_attempts"), ("SCF converged", "states_converged"), ("Spin and identity eligible", "s2_and_identity_eligible_states"), ("Raw computed vertical gaps", "computed_vertical_gaps"), ("Eligible diagnostic gaps", "eligible_diagnostic_gaps"), ("Accepted MECP", "mecp_accepted")]:
        gate_rows.append({"stage": label, "count": spin[key], "evidence": [evidence("C02", "/" + key)]})
    tables.append(table("posthoc_spin_gates.csv", gate_rows))
    sol_rows = []
    for n in (1, 2, 3):
        a = min((r for r in assoc if int(r["n_solvent"]) == n and r["status"] == "converged" and r["connectivity_retained"] == "True"), key=lambda r: float(r["energy_eV"]))
        t = {float(r["temperature_K"]): r for r in thermo if r["name"] == a["name"]}
        b = next(r for r in many if r["name"] == a["name"])
        s = next(r for r in sensitivity["cutoff_sensitivity_at_baseline"] if r["n_alcohol"] == n)
        sol_rows.append({"n_tBuOH": n, "selected_seed": int(a["seed"]), "association_energy_kcal_mol": float(a["association_energy_kcal_mol"]), "association_G_298_15_kcal_mol": float(t[298.15]["delta_G_qRRHO_1M_kcal_mol"]), "association_G_383_15_kcal_mol": float(t[383.15]["delta_G_qRRHO_1M_kcal_mol"]), "beyond_pair_kcal_mol": float(b["beyond_pair_nonadditivity_kcal_mol"]), "cutoff_sensitivity_range_kcal_mol": s["cutoff_range_kcal_mol"], "evidence": [evidence("C09", f"CSV[name={a['name']}]; energy, association and seed columns"), evidence("C10", f"CSV[name={a['name']},temperature_K in {{298.15,383.15}}].delta_G_qRRHO_1M_kcal_mol"), evidence("C11", f"CSV[name={a['name']}].beyond_pair_nonadditivity_kcal_mol"), evidence("C18", f"/cutoff_sensitivity_at_baseline/{n-1}")]})
    tables.append(table("posthoc_solvation_selected.csv", sol_rows))
    wire_rows = []
    for r in wire:
        wire_rows.append({"attempt": r["attempt"], "last_step": int(r["last_NEB_step"]), "residual_eV_A": float(r["last_NEB_fmax_eV_A"]), "target_eV_A": float(r["NEB_fmax_target_eV_A"]), "residual_over_target": float(r["last_NEB_fmax_eV_A"]) / float(r["NEB_fmax_target_eV_A"]), "band_maximum_eV": float(r["last_band_max_above_reactant_eV"]), "accepted_TS": r["accepted_TS"], "evidence": [evidence("C13", f"CSV[attempt={r['attempt']}]; ratio=last_NEB_fmax_eV_A/NEB_fmax_target_eV_A")]})
    tables.append(table("posthoc_proton_wire.csv", wire_rows))
    ml_rows = []
    for split in ("train", "validation", "test"):
        v = egnn["metrics"][split]
        ml_rows.append({"split": split, "pairs": v["samples"], "model_MAE_eV": v["mae_eV"], "model_RMSE_eV": v["rmse_eV"], "equal_reference_MAE_eV": v["equal_reference_energy_baseline_mae_eV"], "training_mean_MAE_eV": v["training_mean_baseline_mae_eV"], "improvement_over_equal_eV": v["equal_reference_energy_baseline_mae_eV"] - v["mae_eV"], "evidence": [evidence("C07", f"/metrics/{split}; improvement=equal_reference_energy_baseline_mae_eV-mae_eV")]})
    tables.append(table("posthoc_egnn_metrics.csv", ml_rows))
    p4_rows = []
    for i, v in enumerate(p4jobs["jobs"]):
        p4_rows.append({"job": v["name"], "basis": v["basis"], "multiplicity": v["multiplicity"], "status": v["status"], "energy_hartree": v.get("energy_hartree"), "S2": v.get("spin_squared"), "max_gradient_hartree_bohr": v.get("max_gradient_hartree_bohr"), "accepted_chemical_label": v["accepted_chemical_label"], "evidence": [evidence("C17", f"/jobs/{i}")]})
    tables.append(table("posthoc_phase4_jobs.csv", p4_rows))
    dft_rows = []
    for i, v in enumerate(dft["runs"]):
        if v["status"] == "resource_preflight_not_launched":
            continue
        dft_rows.append({"run_id": v["run_id"], "method": v["method"], "basis": v["basis"], "charge": v["charge"], "status": v["status"], "energy_hartree": v.get("energy_hartree"), "elapsed_seconds": v["elapsed_seconds"], "scf_stability_analysis_performed": v["scf_stability_analysis_performed"], "evidence": [evidence("C21", f"/runs/{i}")]})
    tables.append(table("posthoc_dft_jobs.csv", dft_rows))

    # Native-text corroboration is a read-only check, with line-numbered excerpts.
    native_checks = []
    native_targets = []
    for v in matrix["state_records"]:
        if v["status"] != "converged":
            continue
        candidate = next(x for x in v["native_evidence"] if x["path"].endswith("psi4.out"))
        path = candidate["path"].replace("\\", "/").split("pincer-catmech-ai/", 1)[1]
        native_targets.append((v["catalyst_id"] + "_M" + str(v["multiplicity"]), path, v["energy_hartree"], candidate["sha256"]))
    for v in dft["runs"]:
        if v["status"] == "completed":
            path = "projects/cu-np-operando/data/dft_pilot/runs/" + v["run_id"] + "/psi4.out"
            native_targets.append((v["run_id"], path, v["energy_hartree"], None))
    for i, (name, path, expected, recorded_hash) in enumerate(native_targets):
        sid = f"N{i+1:02d}"
        try:
            raw = source(sid, path)
            lines = raw.decode("utf-8", errors="replace").splitlines()
            final = [(n+1, line, float(match.group(1))) for n, line in enumerate(lines) if (match := re.search(r"Final Energy:\s*(-?\d+\.\d+)", line))]
            keep = [{"line": n+1, "text": line} for n, line in enumerate(lines) if re.search(r"Final Energy:|SCF converged|Psi4 exiting|Iterations converged|Psi4: An Open", line)]
            rec = {"name": name, "source": entries[sid], "selected_native_lines": keep, "expected_summary_energy_hartree": expected, "native_final_energy_hartree": final[-1][2] if final else None, "energy_absolute_difference_hartree": abs(final[-1][2]-expected) if final else None, "recorded_native_hash_matches": digest(raw) == recorded_hash if recorded_hash else None, "passed": bool(final) and abs(final[-1][2]-expected) < 1e-10}
        except (subprocess.SubprocessError, OSError) as error:
            rec = {"name": name, "path": path, "commit": COMMIT, "passed": False, "read_failure": type(error).__name__, "interpretation": "Native text could not be checked; do not treat the summary as independently corroborated by this audit."}
        native_checks.append(rec)
    (out / "native_energy_checks.json").write_text(json.dumps(native_checks, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    # Commit messages establish historical sequence, not scientific correctness.
    raw_log = git("log", "--reverse", "--format=%H%x09%aI%x09%s", COMMIT).decode().splitlines()
    history = [dict(zip(["commit", "author_date", "subject"], line.split("\t", 2))) for line in raw_log]
    (out / "commit_history.json").write_text(json.dumps(history, indent=2) + "\n", encoding="utf-8")
    checks["native_energy_checks_attempted"] = len(native_checks)
    checks["native_energy_checks_passed"] = sum(v["passed"] for v in native_checks)
    package = {"fixed_commit": COMMIT, "repository": REPOSITORY, "source_entries": entries, "posthoc_tables": tables, "checks": checks, "history": history, "license_note": "Third-party crystal-derived material retains the archived CC BY-NC 4.0 notice C24. No whole-repository LICENSE file or pyproject license field was located at this commit; do not infer a permissive license from prose. These bounded source extracts are retained for the repository owner's requested audit; dependencies and third-party material keep their own terms.", "read_only_scientific_scope": "Git extraction, CSV/JSON selection, status counts, arithmetic ratios and native-final-energy comparisons only; zero new quantum calls, NEB runs or neural training runs."}
    (out / "catalog.json").write_text(json.dumps(package, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"source_files": len(entries), "posthoc_tables": len(tables), "checks": checks}))


if __name__ == "__main__":
    main()
