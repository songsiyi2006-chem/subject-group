"""Independent saved-matrix energy and SCF-record checks; no engine calls."""
import csv
import hashlib
import json
from pathlib import Path
import re
import numpy as np

BASE = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    records = []
    inputs = {}
    for prefix, expected in (("", 55), ("pilot/", 3)):
        out = BASE / prefix
        read = lambda name: json.loads((out/name).read_text(encoding="utf-8"))
        with (out/"quantum_job_accounting.csv").open(encoding="utf-8", newline="") as handle:
            jobs = {row["job_id"]: row for row in csv.DictReader(handle)}
        matrices = read("ao_matrices.json")
        log_records = read("psi4_log_provenance.json")
        errors, electron_errors, energy_errors = [], [], []
        log_ok = len(jobs) == len(matrices) == len(log_records) == expected
        for matrix in matrices:
            job = jobs[matrix["job_id"]]
            S, F, C, eps, D, h = [np.asarray(matrix[k]) for k in ("S", "Fock_Hartree", "C_alpha", "epsilon_alpha_Hartree", "density_alpha", "Hcore_Hartree")]
            electronic = float(np.sum(D*(h+F))) * (1. if job["reference"] == "RHF" else .5)
            errors.append(abs(electronic-float(job["electronic_energy_Hartree"])))
            electron_errors.append(abs(np.trace(D@S)-1))
            errors.append(float(np.max(np.abs(F@C-(S@C)*eps))))
            path = out/"psi4_outputs"/(matrix["job_id"]+".out")
            text = path.read_text(encoding="utf-8")
            log_ok &= "Energy and wave function converged." in text and "[LOCAL_HOST]" in text
            value = re.search(r"@(?:RHF|UHF) Final Energy:\s+([-0-9.]+)", text)
            energy_errors.append(abs(float(value[1])-float(job["energy_Hartree"])))
        for log in log_records:
            log_ok &= sha(out/log["file"]) == log["public_sha256"]
        records.append({"scope": "main" if not prefix else "pilot", "jobs": expected,
                        "max_density_energy_or_Fock_residual_Hartree": max(errors),
                        "max_alpha_electron_trace_error": max(electron_errors),
                        "max_log_vs_saved_energy_error_Hartree": max(energy_errors),
                        "all_recorded_SCF_converged_and_log_hashes_match": bool(log_ok),
                        "passed": bool(log_ok and max(errors) < 1e-11 and max(electron_errors) < 1e-12 and max(energy_errors) < 1e-12)})
        for name in ("ao_matrices.json", "quantum_job_accounting.csv", "psi4_log_provenance.json"):
            inputs[prefix+name] = sha(out/name)
    result = {"passed": all(r["passed"] for r in records), "records": records,
              "scope": "Independent stored-density energy identity, alpha-electron traces, Fock residuals, log convergence and hashes. No quantum recalculation.",
              "inputs_sha256": inputs, "hash_base": "results/electronic", "script_sha256": sha(Path(__file__))}
    (BASE/"saved_energy_identity_check.json").write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8", newline="\n")
    print(json.dumps(result, indent=2))
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
