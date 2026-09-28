"""Auditable generalized eigenproblems and a bounded H2 RHF/STO-3G study.

Run --stage all with a Python environment containing Psi4, NumPy and SciPy.
The source overlap audit, ab initio reference, and distance surrogate have
separate evidence meanings. No Cu catalyst quantum calculation is performed.
"""
import argparse
import ast
from collections import Counter
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import re
import shutil
import time
from typing import Dict, List

import numpy as np
import scipy
import scipy.linalg as la
from scipy.interpolate import CubicSpline

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "source/quantum_egnn_neb_engine.py"
SOURCES = ["https://psicode.org/psi4manual/master/scf.html",
           "https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.eigh.html",
           "https://github.com/psi4/psi4/blob/master/samples/fd-gradient/input.dat"]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False)+"\n", encoding="utf-8", newline="\n")


def write_rows(path, rows):
    if not rows:
        raise ValueError("CSV output must have at least one row")
    with Path(path).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def matrix_pair(H, S):
    H, S = np.asarray(H, dtype=float), np.asarray(S, dtype=float)
    if H.ndim != 2 or H.shape[0] != H.shape[1] or H.size == 0 or S.shape != H.shape:
        raise ValueError("H and S must be finite square matrices of the same nonzero shape")
    if not np.isfinite(H).all() or not np.isfinite(S).all():
        raise ValueError("Matrices must be finite")
    if not np.allclose(H, H.T, rtol=1e-12, atol=1e-12) or not np.allclose(S, S.T, rtol=1e-12, atol=1e-12):
        raise ValueError("Real symmetric matrices required; asymmetric input is not silently repaired")
    return (H+H.T)/2, (S+S.T)/2


def generalized_eigh(H, S, mode="strict", relative_cutoff=1e-10, absolute_cutoff=1e-12, occupied=None):
    """Strict SPD solve or explicit canonical positive-subspace Ritz solve.

    Canonical truncation returns projected and FULL-space residuals separately.
    It cannot repair an inconsistent Hamiltonian on the nullspace of S.
    Negative eigenvalues larger than floating-point tolerance are rejected.
    """
    H, S = matrix_pair(H, S)
    if mode not in ("strict", "canonical") or not np.isfinite([relative_cutoff, absolute_cutoff]).all() or relative_cutoff < 0 or absolute_cutoff <= 0:
        raise ValueError("Invalid mode or positive cutoff controls")
    values, vectors = la.eigh(S)
    scale = max(1., float(np.max(np.abs(values))))
    negative_tolerance = 64*np.finfo(float).eps*scale*len(values)
    if values[0] < -negative_tolerance:
        raise ValueError("Indefinite overlap rejected")
    cutoff = max(absolute_cutoff, relative_cutoff*max(0., float(values[-1])))
    keep = values > cutoff
    rank = int(np.count_nonzero(keep))
    if rank == 0 or (mode == "strict" and rank != len(values)):
        raise ValueError("Overlap is singular or below the strict SPD cutoff")
    if occupied is not None and (not isinstance(occupied, (int, np.integer)) or occupied < 0 or occupied > rank):
        raise ValueError("Occupied orbital count exceeds retained metric rank")
    X = vectors[:, keep] / np.sqrt(values[keep])[None, :]
    energies, transformed = la.eigh(X.T @ H @ X)
    C = X @ transformed
    residual = H @ C - (S @ C)*energies[None, :]
    normalization = max(1., la.norm(H, 2)*la.norm(C, 2), la.norm(S, 2)*la.norm(C*energies[None, :], 2))
    null = vectors[:, ~keep]
    return {"energies": energies, "coefficients": C, "overlap_eigenvalues": values,
            "rank": rank, "discarded": len(values)-rank, "cutoff": cutoff,
            "orthonormality_max_abs": float(np.max(np.abs(C.T@S@C-np.eye(rank)))),
            "full_residual_max_abs": float(np.max(np.abs(residual))),
            "full_residual_relative_2norm": float(la.norm(residual, 2)/normalization),
            "projected_residual_max_abs": float(np.max(np.abs(X.T@residual))),
            "discarded_H_coupling_2norm": float(la.norm(null.T@H@X, 2)) if null.size else 0.,
            "interpretation": "full-space SPD eigensolution" if mode == "strict" else "Ritz solution restricted to retained positive-overlap subspace"}


def closed_shell_occupied(electrons, rank):
    if not isinstance(electrons, (int, np.integer)) or electrons < 0 or electrons % 2:
        raise ValueError("Closed-shell occupancy requires a nonnegative even integer electron count")
    if electrons//2 > rank:
        raise ValueError("Electron count cannot fit the retained orbital space")
    return electrons//2


def clipped_source_solution(H, S, floor):
    H, S = matrix_pair(H, S)
    if floor <= 0 or not math.isfinite(floor):
        raise ValueError("Positive floor required")
    values, U = la.eigh(S)
    X = (U / np.sqrt(np.maximum(values, floor))) @ U.T
    energies, coefficients = la.eigh(X@H@X)
    C = X@coefficients
    original_residual = H@C-(S@C)*energies[None, :]
    modified = (U*np.maximum(values, floor))@U.T
    modified_residual = H@C-(modified@C)*energies[None, :]
    return {"floor": floor, "energies": energies, "coefficients": C,
            "original_metric_orthonormality_max_abs": float(np.max(np.abs(C.T@S@C-np.eye(len(S))))),
            "original_generalized_residual_max_abs": float(np.max(np.abs(original_residual))),
            "modified_generalized_residual_max_abs": float(np.max(np.abs(modified_residual))),
            "modified_metric_distance_2norm": float(la.norm(modified-S, 2))}


def source_class():
    tree = ast.parse(SOURCE.read_text(encoding="utf-8"))
    cls = next(node for node in tree.body if isinstance(node, ast.ClassDef) and node.name == "ExtendedHuckelQuantumSolver")
    namespace = {"np": np, "la": la, "List": List, "Dict": Dict}
    exec(compile(ast.Module(body=[cls], type_ignores=[]), "archived_source_EHT_class", "exec"), namespace)
    return namespace[cls.name]


def source_matrices(elements, coordinates):
    solver = source_class()(elements, np.asarray(coordinates, dtype=float))
    n = solver.n_basis
    S = np.array([[solver._compute_overlap(i, j) for j in range(n)] for i in range(n)])
    diagonals = np.array([orbital["h_ii"] for orbital in solver.basis_orbitals])
    H = .5*1.75*S*(diagonals[:, None]+diagonals[None, :])
    np.fill_diagonal(H, diagonals)
    return H, S, solver.basis_orbitals


def audit_source(out):
    started = time.perf_counter()
    capture_path = ROOT / "results/original/captured_results.json"
    capture = json.loads(capture_path.read_text(encoding="utf-8"))
    spectra, floors, cases, matrices = [], [], [], {}
    for name, field in (("reactant", "R_reactant"), ("product", "R_product")):
        H, S, orbitals = source_matrices(capture["elements"], capture[field])
        values = la.eigvalsh(S)
        canonical = generalized_eigh(H, S, mode="canonical")
        strict_error, occupancy_error = "", ""
        try:
            generalized_eigh(H, S)
        except ValueError as exc:
            strict_error = str(exc)
        try:
            closed_shell_occupied(40, canonical["rank"])
        except ValueError as exc:
            occupancy_error = str(exc)
        same_atom = [(i, j) for i in range(len(S)) for j in range(i+1, len(S)) if orbitals[i]["atom_idx"] == orbitals[j]["atom_idx"]]
        rows_equal = max(float(np.max(np.abs(S[i]-S[j]))) for i, j in same_atom)
        cases.append({"case": name, "basis_dimension": len(S), "atoms": len(capture["elements"]), "valence_electrons": 40,
                      "source_occupied_orbitals": 20, "positive_overlap_rank": canonical["rank"],
                      "discarded_null_directions": canonical["discarded"], "minimum_overlap_eigenvalue": float(values[0]),
                      "same_atom_offdiagonal_pairs": len(same_atom), "same_atom_overlaps_all_one": all(S[i, j] == 1 for i, j in same_atom),
                      "same_atom_overlap_row_max_difference": rows_equal,
                      "strict_solver_rejection": strict_error, "closed_shell_capacity_rejection": occupancy_error,
                      "canonical_projected_residual_max_abs": canonical["projected_residual_max_abs"],
                      "canonical_full_residual_max_abs": canonical["full_residual_max_abs"],
                      "canonical_metric_orthonormality_max_abs": canonical["orthonormality_max_abs"],
                      "canonical_discarded_H_coupling_2norm": canonical["discarded_H_coupling_2norm"],
                      "canonical_status": "numerical projected control only; cannot host 40 electrons or repair the supplied EHT model"})
        for index, value in enumerate(values):
            spectra.append({"case": name, "index": index, "overlap_eigenvalue": float(value), "retained_at_1e_10_relative": bool(value > canonical["cutoff"])})
        for floor in (1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8):
            solution = clipped_source_solution(H, S, floor)
            e = solution["energies"]
            floors.append({"case": name, "clipping_floor": floor, "HOMO_eV": float(e[19]), "LUMO_eV": float(e[20]),
                          "occupied_eigenvalue_sum_times_two_eV": float(2*sum(e[:20])),
                          **{k: v for k, v in solution.items() if k not in ("floor", "energies", "coefficients")}})
        matrices[name] = {"H_eV": H.tolist(), "S_dimensionless": S.tolist(), "orbitals": orbitals,
                          "coordinates_A": capture[field], "canonical_energies_eV": canonical["energies"].tolist(),
                          "canonical_coefficients": canonical["coefficients"].tolist()}
    np_error = next(r for r in floors if r["case"] == "reactant" and r["clipping_floor"] == 1e-5)
    returned = capture["eht_results"]
    source_match = abs(round(np_error["occupied_eigenvalue_sum_times_two_eV"]*23.060541945329, 2)-returned["electronic_energy_kcal_mol"]) < .03
    controls = []
    control_matrices = []
    for seed in range(901, 907):
        rng = np.random.default_rng(seed)
        A, B = rng.normal(size=(6, 6)), rng.normal(size=(6, 6))
        S, H = A.T@A+.2*np.eye(6), (B+B.T)/2
        ours = generalized_eigh(H, S)
        direct, _ = la.eigh(H, S)
        controls.append({"case": f"SPD_seed{seed}", "retained_rank": ours["rank"], "energy_max_abs_vs_scipy": float(np.max(np.abs(ours["energies"]-direct))),
                         "orthonormality_max_abs": ours["orthonormality_max_abs"], "full_residual_max_abs": ours["full_residual_max_abs"],
                         "projected_residual_max_abs": ours["projected_residual_max_abs"]})
        control_matrices.append({"case": f"SPD_seed{seed}", "H": H.tolist(), "S": S.tolist(), "energies": ours["energies"].tolist(), "coefficients": ours["coefficients"].tolist()})
    for name, H in (("compatible_PSD", np.diag([2., 3., 0.])), ("incompatible_PSD", np.array([[2., 0., .2], [0., 3., .4], [.2, .4, 1.]]))):
        S = np.diag([1., 2., 0.])
        ours = generalized_eigh(H, S, mode="canonical")
        controls.append({"case": name, "retained_rank": ours["rank"], "energy_max_abs_vs_scipy": None,
                         "orthonormality_max_abs": ours["orthonormality_max_abs"], "full_residual_max_abs": ours["full_residual_max_abs"],
                         "projected_residual_max_abs": ours["projected_residual_max_abs"]})
        control_matrices.append({"case": name, "H": H.tolist(), "S": S.tolist(), "energies": ours["energies"].tolist(), "coefficients": ours["coefficients"].tolist()})
    write_rows(out / "source_overlap_spectrum.csv", spectra)
    write_rows(out / "source_clipping_sensitivity.csv", floors)
    write_rows(out / "solver_controls.csv", controls)
    dump(out / "source_matrices.json", matrices)
    dump(out / "solver_control_matrices.json", control_matrices)
    audit = {"cases": cases, "reactant_rounded_capture_matches_reconstruction": source_match,
             "captured_eht_results": returned, "source_orbital_sum_is_not_ab_initio_total_energy": True,
             "same_atom_defect": "All distinct same-atom s/p/d functions have overlap 1 and identical overlap rows; angular orbital orthogonality is absent",
             "clipping_interpretation": "Solves a modified metric with added nullspace eigenvalues; not the original secular problem",
             "canonical_interpretation": "Controls rank and conditioning but does not supply a physical orbital overlap or a 40-electron solution",
             "counts": {"source_matrix_assemblies": 2, "clipped_metric_controls": 12, "source_canonical_controls": 2,
                        "source_strict_rejections": 2, "source_electron_capacity_rejections": 2, "independent_SPD_controls": 6,
                        "independent_PSD_controls": 2, "ab_initio_jobs": 0, "new_neural_forwards": 0},
             "runtime_seconds": time.perf_counter()-started,
             "inputs_sha256": {SOURCE.relative_to(ROOT).as_posix(): sha(SOURCE), capture_path.relative_to(ROOT).as_posix(): sha(capture_path)}}
    dump(out / "source_audit.json", audit)
    return audit


def geometry_plan(pilot=False):
    if pilot:
        return [{"point_id": "pilot_R0.740", "R_A": .74, "split": "pilot"}]
    points = []
    for i in range(33):
        split = "train" if i % 2 == 0 else "validation" if i % 4 == 1 else "test"
        points.append({"point_id": f"scan_{i:03d}", "R_A": round(.50+.04*i, 8), "split": split})
    points.extend({"point_id": f"ood_{i:03d}", "R_A": round(1.90+.10*i, 8), "split": "ood_stretch"} for i in range(9))
    return points


def sanitize_log(path):
    text = path.read_text(encoding="utf-8", errors="replace")
    raw_digest = sha(path)
    # Keep numerical and SCF evidence; remove host installation/scratch identity.
    text = re.sub(r"[A-Za-z]:[\\/][^\n\r\"<>]*", "[LOCAL_PATH]", text)
    text = re.sub(r"(?m)(called on|Host:)\s+[^\r\n]+", r"\1 [LOCAL_HOST]", text)
    path.write_text(text, encoding="utf-8", newline="\n")
    return raw_digest


def quantum_reference(out, pilot=False):
    import psi4
    started = time.perf_counter()
    psi4.set_num_threads(1)
    psi4.set_memory("500 MB")
    raw = out / "psi4_outputs"
    scratch = out / "scratch"
    raw.mkdir(exist_ok=True)
    scratch.mkdir(exist_ok=True)
    psi4.core.IOManager.shared_object().set_default_path(str(scratch))
    bohr_A = float(psi4.constants.bohr2angstroms)
    logs, points, matrices, checks, jobs = [], [], [], [], []

    def calculate(job_id, distance=None, gradient=True):
        log = raw / (job_id+".out")
        psi4.core.clean()
        psi4.core.clean_options()
        psi4.core.set_output_file(str(log), False)
        if distance is None:
            geometry = "0 2\nH 0 0 0\nunits angstrom\nsymmetry c1\nno_reorient\nno_com\n"
            method = "UHF"
        else:
            geometry = f"0 1\nH 0 0 {-distance/2:.12f}\nH 0 0 {distance/2:.12f}\nunits angstrom\nsymmetry c1\nno_reorient\nno_com\n"
            method = "RHF"
        mol = psi4.geometry(geometry)
        psi4.set_options({"basis": "sto-3g", "reference": method, "scf_type": "pk", "guess": "core", "e_convergence": 1e-12,
                          "d_convergence": 1e-12, "maxiter": 100, "s_orthogonalization": "symmetric"})
        began = time.perf_counter()
        if gradient:
            grad, wfn = psi4.gradient("scf", molecule=mol, return_wfn=True)
            grad = np.asarray(grad).copy()/bohr_A
            energy = float(wfn.energy())
        else:
            energy, wfn = psi4.energy("scf", molecule=mol, return_wfn=True)
            energy = float(energy)
            grad = None
        mints = psi4.core.MintsHelper(wfn.basisset())
        matrix = {"job_id": job_id, "S": np.asarray(mints.ao_overlap()).tolist(), "Fock_Hartree": np.asarray(wfn.Fa()).tolist(),
                  "C_alpha": np.asarray(wfn.Ca()).tolist(), "epsilon_alpha_Hartree": np.asarray(wfn.epsilon_a()).tolist(),
                  "density_alpha": np.asarray(wfn.Da()).tolist(), "Hcore_Hartree": (np.asarray(mints.ao_kinetic())+np.asarray(mints.ao_potential())).tolist(),
                  "gradient_Hartree_A": None if grad is None else grad.tolist()}
        nuclear = float(mol.nuclear_repulsion_energy())
        psi4.core.close_outfile()
        raw_digest = sanitize_log(log)
        jobs.append({"job_id": job_id, "molecule": "H" if distance is None else "H2", "R_A": distance,
                     "reference": method, "basis": "STO-3G", "charge": 0, "multiplicity": 2 if distance is None else 1,
                     "energy_Hartree": energy, "nuclear_repulsion_Hartree": nuclear, "electronic_energy_Hartree": energy-nuclear,
                     "analytic_gradient_requested": gradient, "elapsed_seconds": time.perf_counter()-began, "SCF_completed": True})
        logs.append({"file": log.relative_to(out).as_posix(), "raw_before_path_redaction_sha256": raw_digest,
                     "public_sha256": sha(log), "redaction": "machine paths and host identity only; numerical SCF output retained"})
        matrices.append(matrix)
        psi4.core.clean()
        return energy, grad, matrix, nuclear

    for point in geometry_plan(pilot):
        energy, grad, matrix, nuclear = calculate(point["point_id"], point["R_A"])
        eigen = generalized_eigh(matrix["Fock_Hartree"], matrix["S"])
        C, S, F = [np.asarray(matrix[k]) for k in ("C_alpha", "S", "Fock_Hartree")]
        eps = np.asarray(matrix["epsilon_alpha_Hartree"])
        bond_gradient = .5*(grad[1, 2]-grad[0, 2])
        points.append({**point, "energy_total_Hartree": energy, "energy_electronic_Hartree": energy-nuclear,
                       "nuclear_repulsion_Hartree": nuclear, "dE_dR_Hartree_A": float(bond_gradient),
                       "force_atom2_z_Hartree_A": float(-grad[1, 2]), "net_force_max_Hartree_A": float(np.max(np.abs(grad.sum(0)))),
                       "Fock_generalized_residual_max": float(np.max(np.abs(F@C-(S@C)*eps))),
                       "MO_orthonormality_max": float(np.max(np.abs(C.T@S@C-np.eye(len(C))))),
                       "eigenvalue_difference_reviewed_solver": float(np.max(np.abs(eigen["energies"]-eps))),
                       "orbital_sum_times_two_Hartree": float(2*eps[0]), "SCF_converged": True})
    # Independent displaced-energy checks; these jobs are not training examples.
    anchors = [.74] if pilot else [.70, 1.14, 1.78]
    steps = [.001] if pilot else [.001, .0005]
    for anchor in anchors:
        match = next(p for p in points if abs(p["R_A"]-anchor) < 1e-9)
        for step in steps:
            plus = calculate(f"fd_R{anchor:.3f}_h{step:.4f}_plus", anchor+step, gradient=False)[0]
            minus = calculate(f"fd_R{anchor:.3f}_h{step:.4f}_minus", anchor-step, gradient=False)[0]
            fd = (plus-minus)/(2*step)
            checks.append({"R_A": anchor, "step_A": step, "plus_energy_Hartree": plus, "minus_energy_Hartree": minus,
                           "finite_difference_dE_dR_Hartree_A": fd, "analytic_dE_dR_Hartree_A": match["dE_dR_Hartree_A"],
                           "absolute_error_Hartree_A": abs(fd-match["dE_dR_Hartree_A"])})
    atomic = None
    if not pilot:
        energy, _, _, _ = calculate("isolated_H_doublet", None, gradient=False)
        atomic = {"method": "UHF/STO-3G", "charge": 0, "multiplicity": 2, "energy_Hartree": energy,
                  "twice_atomic_H_energy_Hartree": 2*energy,
                  "scope": "Same minimal basis separated-atom reference; restricted H2 does not recover correct dissociation"}
    write_rows(out / "h2_reference.csv", points)
    write_rows(out / "h2_gradient_finite_differences.csv", checks)
    write_rows(out / "quantum_job_accounting.csv", jobs)
    dump(out / "ao_matrices.json", matrices)
    dump(out / "psi4_log_provenance.json", logs)
    summary = {"evidence": "Executed ab initio nonrelativistic RHF/STO-3G Born-Oppenheimer total energies and analytic gradients for H2; UHF isolated H reference",
               "pilot": pilot, "engine": "Psi4", "engine_version": psi4.__version__, "python": platform.python_version(),
               "numpy": np.__version__, "threads": 1, "memory_MB_requested": 500, "SCF_algorithm": "PK", "charge": 0, "H2_multiplicity": 1,
               "SCF_energy_and_density_thresholds": 1e-12, "bohr_in_angstrom": bohr_A,
               "split_sizes": dict(Counter(p["split"] for p in points)), "split_rule": "Fixed R grid before labels: alternating training points; intervening odd indices alternate validation/test; longer bonds are OOD",
               "partition_scope": "Geometry interpolation/extrapolation of ONE molecule, not molecular generalization", "isolated_H_reference": atomic,
               "counts": {"H2_scan_SCF_jobs": len(points), "H2_analytic_gradients": len(points), "displaced_energy_SCF_jobs": 2*len(checks),
                          "atomic_H_SCF_jobs": int(not pilot), "total_SCF_jobs": len(jobs), "failed_SCF_jobs": 0,
                          "saved_AO_matrix_sets": len(matrices), "finite_difference_comparisons": len(checks), "Cu_cluster_ab_initio_jobs": 0, "DFT_jobs": 0},
               "maximum_force_balance_error_Hartree_A": max(p["net_force_max_Hartree_A"] for p in points),
               "maximum_generalized_residual": max(p["Fock_generalized_residual_max"] for p in points),
               "maximum_MO_orthonormality_error": max(p["MO_orthonormality_max"] for p in points),
               "max_gradient_FD_error_Hartree_A": max(p["absolute_error_Hartree_A"] for p in checks),
               "scan_minimum": min(points, key=lambda p: p["energy_total_Hartree"]),
               "limitations": ["Minimal basis and RHF neglect electron correlation; stretched RHF H2 is not an accurate dissociation curve",
                               "Grid minimum is not a geometry optimization", "No thermochemistry, DFT, Cu-N4, indole or electrochemical environment in this quantum reference"],
               "runtime_seconds": time.perf_counter()-started, "sources": SOURCES}
    dump(out / "reference_summary.json", summary)
    return summary


def radial_features(r, centers, lengthscale):
    r = np.asarray(r, dtype=float).reshape(-1)
    centers = np.asarray(centers, dtype=float).reshape(-1)
    if not np.isfinite(r).all() or not np.isfinite(centers).all() or lengthscale <= 0:
        raise ValueError("Finite distances and positive radial lengthscale required")
    delta = r[:, None]-centers[None, :]
    phi = np.exp(-.5*(delta/lengthscale)**2)
    return np.column_stack([np.ones(len(r)), r, phi]), np.column_stack([np.zeros(len(r)), np.ones(len(r)), -delta/lengthscale**2*phi])


def fit_radial(r, energy, gradient, lengthscale, ridge, use_gradients):
    r, energy, gradient = np.asarray(r), np.asarray(energy), np.asarray(gradient)
    if len(r) < 3 or len(np.unique(r)) != len(r) or len(energy) != len(r) or len(gradient) != len(r) or ridge <= 0:
        raise ValueError("Distinct aligned training distances and positive ridge required")
    features, derivatives = radial_features(r, r, lengthscale)
    offset, e_scale, g_scale = float(energy.mean()), float(energy.std()), float(gradient.std())
    if min(e_scale, g_scale) <= 0:
        raise ValueError("Nonconstant training data required")
    design, target = features/e_scale, (energy-offset)/e_scale
    if use_gradients:
        design = np.vstack([design, derivatives/g_scale])
        target = np.concatenate([target, gradient/g_scale])
    design = np.vstack([design, math.sqrt(ridge)*np.eye(design.shape[1])])
    target = np.concatenate([target, np.zeros(design.shape[1])])
    coefficients = np.linalg.lstsq(design, target, rcond=None)[0]
    return {"centers_A": r.tolist(), "coefficients_Hartree": coefficients.tolist(), "energy_offset_Hartree": offset,
            "lengthscale_A": lengthscale, "ridge": ridge, "uses_training_gradients": use_gradients,
            "training_energy_scale_Hartree": e_scale, "training_gradient_scale_Hartree_A": g_scale}


def predict_radial(model, distances):
    features, derivatives = radial_features(distances, model["centers_A"], model["lengthscale_A"])
    coefficients = np.asarray(model["coefficients_Hartree"])
    return model["energy_offset_Hartree"] + features@coefficients, derivatives@coefficients


def fit_surrogates(out):
    started = time.perf_counter()
    with (out / "h2_reference.csv").open(newline="", encoding="utf-8") as handle:
        points = list(csv.DictReader(handle))
    train = [p for p in points if p["split"] == "train"]
    validation = [p for p in points if p["split"] == "validation"]
    if len(train) != 17 or len(validation) != 8:
        raise ValueError("Surrogate fitting expects frozen 17-train/8-validation main scan")
    array = lambda group, key: np.array([float(p[key]) for p in group])
    r, e, g = [array(train, k) for k in ("R_A", "energy_total_Hartree", "dE_dR_Hartree_A")]
    rv, ev, gv = [array(validation, k) for k in ("R_A", "energy_total_Hartree", "dE_dR_Hartree_A")]
    candidates, models = [], {}
    for use_gradient in (False, True):
        name = "RBF_energy_force" if use_gradient else "RBF_energy_only"
        choices = []
        for lengthscale in (.06, .12, .24, .48):
            for ridge in (1e-10, 1e-6):
                model = fit_radial(r, e, g, lengthscale, ridge, use_gradient)
                pe, pg = predict_radial(model, rv)
                ermse, grmse = float(np.sqrt(np.mean((pe-ev)**2))), float(np.sqrt(np.mean((pg-gv)**2)))
                score = ermse/e.std() + grmse/g.std()
                row = {"model": name, "lengthscale_A": lengthscale, "ridge": ridge, "validation_energy_RMSE_Hartree": ermse,
                       "validation_gradient_RMSE_Hartree_A": grmse, "validation_dimensionless_joint_score": float(score)}
                candidates.append(row)
                choices.append((score, model))
        score, model = min(choices, key=lambda p: p[0])
        model["selected_validation_score"] = float(score)
        models[name] = model
    all_r = array(points, "R_A")
    pred = {name: predict_radial(model, all_r) for name, model in models.items()}
    spline = CubicSpline(r, e, bc_type="natural", extrapolate=True)
    pred["cubic_spline_energy_only"] = (spline(all_r), spline(all_r, 1))
    bins = np.clip(np.searchsorted(r, all_r, side="right")-1, 0, len(r)-2)
    slopes = np.diff(e)/np.diff(r)
    pred["linear_interpolation_energy_only"] = (e[bins] + slopes[bins]*(all_r-r[bins]), slopes[bins])
    pred["training_mean_zero_force"] = (np.full(len(points), e.mean()), np.zeros(len(points)))
    predictions, metrics = [], []
    true_e, true_g = array(points, "energy_total_Hartree"), array(points, "dE_dR_Hartree_A")
    for name, (pe, pg) in pred.items():
        for i, point in enumerate(points):
            predictions.append({"model": name, "point_id": point["point_id"], "split": point["split"], "R_A": float(point["R_A"]),
                                "reference_energy_Hartree": float(true_e[i]), "predicted_energy_Hartree": float(pe[i]),
                                "reference_dE_dR_Hartree_A": float(true_g[i]), "predicted_dE_dR_Hartree_A": float(pg[i]),
                                "predicted_force_atom2_z_Hartree_A": float(-pg[i])})
        for split in ("train", "validation", "test", "ood_stretch"):
            mask = np.array([p["split"] == split for p in points])
            metrics.append({"model": name, "split": split, "points": int(mask.sum()),
                            "energy_RMSE_Hartree": float(np.sqrt(np.mean((pe[mask]-true_e[mask])**2))),
                            "gradient_RMSE_Hartree_A": float(np.sqrt(np.mean((pg[mask]-true_g[mask])**2))),
                            "energy_MAE_Hartree": float(np.mean(np.abs(pe[mask]-true_e[mask]))),
                            "gradient_MAE_Hartree_A": float(np.mean(np.abs(pg[mask]-true_g[mask])))})
    write_rows(out / "surrogate_hyperparameter_validation.csv", candidates)
    write_rows(out / "surrogate_predictions.csv", predictions)
    write_rows(out / "surrogate_metrics.csv", metrics)
    dump(out / "radial_surrogate_models.json", models)
    dump(out / "surrogate_split.json", {"split_point_ids": {split: [p["point_id"] for p in points if p["split"] == split] for split in ("train", "validation", "test", "ood_stretch")},
          "training_distances_A": r.tolist(), "hyperparameter_selection": "Joint validation RMSE scaled by TRAIN energy and gradient standard deviations; no test or OOD selection or refit",
          "seed": None, "determinism": "Fixed distances, deterministic linear least squares; no pseudorandom learning seed"})
    summary = {"evidence": "Distance-only supervised energy/gradient surrogate of actual H2 RHF/STO-3G reference, not an EGNN or Cu catalyst potential",
               "counts": {"training_points": 17, "validation_points": 8, "test_points": 8, "OOD_points": 9,
                          "RBF_linear_fits": 16, "selected_RBF_models": 2, "baseline_models": 3, "prediction_rows": 210, "metric_rows": 20,
                          "neural_training_runs": 0, "new_quantum_jobs_during_fitting": 0},
               "selected_hyperparameters": {name: {k: model[k] for k in ("lengthscale_A", "ridge", "selected_validation_score")} for name, model in models.items()},
               "energy_gradient_consistency": "RBF and cubic derivatives are analytic derivatives of the fitted energy; linear derivative is piecewise constant with knots not differentiable",
               "limitations": ["Same H2 molecule at interleaved distances, not graph or chemical-space generalization", "OOD model error is against imperfect restricted HF",
                               "No neural network, uncertainty calibration, molecular dissociation guarantee or independent external dataset"],
               "inputs_sha256": {"h2_reference.csv": sha(out / "h2_reference.csv")}, "input_hash_base": "this electronic result directory",
               "runtime_seconds": time.perf_counter()-started}
    dump(out / "surrogate_summary.json", summary)
    return summary


def assemble_summary(out, pilot):
    sections = {}
    for name, file in (("source_audit", "source_audit.json"), ("quantum_reference", "reference_summary.json"), ("surrogate", "surrogate_summary.json")):
        if (out/file).is_file():
            sections[name] = json.loads((out/file).read_text(encoding="utf-8"))
    outputs = {p.relative_to(out).as_posix(): sha(p) for p in sorted(out.rglob("*")) if p.is_file()
               and p.name not in ("summary.json", "timer.dat") and "scratch" not in p.parts and "environment" not in p.parts and "pilot" not in p.relative_to(out).parts}
    summary = {"schema_version": 1, "pilot": pilot, "evidence": "Separated source-matrix diagnostic, actual ab initio H2 reference, and deterministic distance surrogate",
               "completed_sections": list(sections), "sections": {k: {"counts": v["counts"]} for k, v in sections.items()},
               "inputs_sha256": {"scripts/electronic_reviewed.py": sha(out / "executed_script.py.txt"), "source/quantum_egnn_neb_engine.py": sha(SOURCE)},
               "input_hash_base": "quantumequi module root", "outputs_sha256": outputs, "output_hash_base": "this result directory",
               "executed_source_snapshot": {"file": "executed_script.py.txt", "sha256": sha(out / "executed_script.py.txt"),
                                            "reason": "Exact code executed for this run; pilot may precede log-host redaction improvements"},
               "versions": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
               "sources": SOURCES, "no_physical_experiment": True}
    dump(out / "summary.json", summary)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=("all", "audit", "quantum", "fit", "summarize"), default="all")
    parser.add_argument("--pilot", action="store_true")
    args = parser.parse_args()
    out = ROOT / "results/electronic" / ("pilot" if args.pilot else "")
    out.mkdir(parents=True, exist_ok=True)
    # Psi4 may create timer.dat in cwd; confine all engine products here.
    os.chdir(out)
    if args.stage != "summarize":
        snapshot = out / "executed_script.py.txt"
        if snapshot.exists() and sha(snapshot) != sha(Path(__file__)):
            raise RuntimeError("Code differs from frozen executed snapshot; preserve prior run in a separate directory before rerunning")
        shutil.copyfile(Path(__file__), snapshot)
    if args.stage in ("all", "audit"):
        audit_source(out)
    if args.stage in ("all", "quantum"):
        quantum_reference(out, args.pilot)
    if args.stage in ("all", "fit") and not args.pilot:
        fit_surrogates(out)
    summary = assemble_summary(out, args.pilot)
    print(json.dumps({"pilot": args.pilot, "completed_sections": summary["completed_sections"], "counts": summary["sections"]}, indent=2))


if __name__ == "__main__":
    main()
