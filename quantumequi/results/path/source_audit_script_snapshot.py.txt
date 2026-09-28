"""Source NEB replay and callback-based CI-NEB numerical benchmarks.

The inertial optimizer/improved tangent are adapted from the repository's
SynthaPore numerical implementation, not a new method. All analytic units are
dimensionless. Source neural outputs retain nominal, uncalibrated source units.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import platform
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "quantumequi"
OUT = BASE / "results/path"
SOURCE = BASE / "source/quantum_egnn_neb_engine.py"
PREDECESSOR = ROOT / "synthapore/scripts/dynamics_reviewed.py"
REFERENCES = [
    "https://doi.org/10.1063/1.1323224",
    "https://doi.org/10.1063/1.1329672",
    "https://doi.org/10.1103/PhysRevLett.97.170201",
]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump(path, data):
    Path(path).write_text(json.dumps(data, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def table(path, rows):
    with Path(path).open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


class AnalyticSurface:
    """A designed reaction coordinate confined to a curved 2D or 3D valley."""
    def __init__(self, name):
        if name not in ("tilted_sine_2d", "periodic_curve_3d"):
            raise ValueError("Unknown analytic surface")
        self.name = name
        self.dimension = 2 if name == "tilted_sine_2d" else 3
        self.stiffness = np.array([5.]) if self.dimension == 2 else np.array([4., 6.])

    def terms(self, x):
        x = np.asarray(x)
        if self.dimension == 2:
            u, up, upp = (x*x-1)**2+.35*x, 4*x*(x*x-1)+.35, 12*x*x-4
            f = (.55*np.sin(1.7*x))[..., None]
            fp = (.935*np.cos(1.7*x))[..., None]
            fpp = (-1.5895*np.sin(1.7*x))[..., None]
        else:
            u, up, upp = 2*(1-np.cos(x)), 2*np.sin(x), 2*np.cos(x)
            f = np.stack([.45*np.sin(x), .30*(1-np.cos(x))], axis=-1)
            fp = np.stack([.45*np.cos(x), .30*np.sin(x)], axis=-1)
            fpp = np.stack([-.45*np.sin(x), .30*np.cos(x)], axis=-1)
        return u, up, upp, f, fp, fpp

    def energy_force(self, points):
        q = np.asarray(points, dtype=float)
        if q.shape[-1] != self.dimension or not np.isfinite(q).all():
            raise ValueError("Wrong dimension or nonfinite coordinates")
        u, up, _, f, fp, _ = self.terms(q[..., 0])
        delta = q[..., 1:] - f
        energy = u + .5*np.sum(self.stiffness*delta**2, axis=-1)
        gx = up-np.sum(self.stiffness*delta*fp, axis=-1)
        gradient = np.concatenate([gx[..., None], self.stiffness*delta], axis=-1)
        return energy, -gradient

    def hessian(self, point):
        q = np.asarray(point, float)
        if q.shape != (self.dimension,) or not np.isfinite(q).all():
            raise ValueError("One finite point required")
        _, _, upp, f, fp, fpp = self.terms(q[0])
        delta = q[1:]-f
        h = np.diag(np.concatenate([[0.], self.stiffness]))
        h[0, 0] = upp+np.sum(self.stiffness*(fp*fp-delta*fpp))
        h[0, 1:] = h[1:, 0] = -self.stiffness*fp
        return h

    def references(self):
        if self.dimension == 2:
            roots = np.sort(np.roots([4., 0., -4., .35]).real)
            xs = roots[[0, 2, 1]]
        else:
            xs = np.array([0., 2*np.pi, np.pi])
        points = np.array([np.r_[x, self.terms(x)[3]] for x in xs])
        energies, forces = self.energy_force(points)
        return {"left": points[0].tolist(), "right": points[1].tolist(),
                "saddle": points[2].tolist(), "energies": energies.tolist(),
                "forward_barrier": float(energies[2]-energies[0]),
                "reverse_barrier": float(energies[2]-energies[1]),
                "reaction_energy": float(energies[1]-energies[0]),
                "stationary_force_max": float(np.linalg.norm(forces, axis=1).max()),
                "hessian_eigenvalues": [np.linalg.eigvalsh(self.hessian(q)).tolist() for q in points],
                "valley_not_assumed_to_be_MEP": True}


def improved_tangent(previous, current, following, ep, ec, en):
    """Energy-weighted tangent, generalized to an arbitrary coordinate shape."""
    ahead, behind = following-current, current-previous
    if en > ec > ep:
        tangent = ahead
    elif en < ec < ep:
        tangent = behind
    else:
        big, small = max(abs(en-ec), abs(ep-ec)), min(abs(en-ec), abs(ep-ec))
        tangent = ahead*big+behind*small if en > ep else ahead*small+behind*big
        if np.linalg.norm(tangent) < 1e-14:
            tangent = ahead+behind
    norm = np.linalg.norm(tangent)
    if norm <= 1e-14:
        raise ValueError("Degenerate neighboring images")
    return tangent/norm


def source_tangent(previous, current, following):
    ahead, behind = following-current, current-previous
    na, nb = np.linalg.norm(ahead), np.linalg.norm(behind)
    if min(na, nb) <= 1e-14:
        raise ValueError("Coincident source images")
    tangent = ahead/na+behind/nb
    norm = np.linalg.norm(tangent)
    if norm <= 1e-14:
        raise ValueError("Antiparallel source tangent")
    return tangent/norm


def projected_forces(band, energies, true, spring=2., climb=True, source=False):
    q, e, f = np.asarray(band), np.asarray(energies), np.asarray(true)
    if q.shape != f.shape or e.shape != (len(q),) or len(q) < 3:
        raise ValueError("Mismatched band, energy, force arrays")
    out = np.zeros_like(q, dtype=float)
    ci = 1+int(np.argmax(e[1:-1]))
    for i in range(1, len(q)-1):
        tangent = source_tangent(q[i-1], q[i], q[i+1]) if source else improved_tangent(
            q[i-1], q[i], q[i+1], *e[i-1:i+2])
        parallel = np.sum(f[i]*tangent)*tangent
        if climb and i == ci:
            out[i] = f[i]-2*parallel
        else:
            out[i] = f[i]-parallel+spring*(np.linalg.norm(q[i+1]-q[i])-np.linalg.norm(q[i]-q[i-1]))*tangent
    return out, ci


def source_update(band, energies, true, spring=1.8, climb=False, sequential=True):
    """Diagnostic algebra reproduces source update; no neural evaluation."""
    old = np.asarray(band, float)
    q = old.copy()
    ci = 1+int(np.argmax(energies[1:-1]))
    for i in range(1, len(q)-1):
        geometry = q if sequential else old
        t = source_tangent(geometry[i-1], geometry[i], geometry[i+1])
        parallel = np.sum(true[i]*t)*t
        force = true[i]-2*parallel if climb and i == ci else true[i]-parallel+spring*(
            np.linalg.norm(geometry[i+1]-geometry[i])-np.linalg.norm(geometry[i]-geometry[i-1]))*t
        q[i] = old[i]+.04*force
    return q


def optimize_neb(energy_force, endpoints, n_images=11, tolerance=1e-5,
                 max_iterations=5000, seed=11, spring=2., climb_after=50,
                 initial_band=None, perturbation=.06):
    """Snapshot-based CI-NEB with explicit force stop and endpoint checks.

    Callback consumes a batch of arbitrary shaped image coordinates and returns
    one scalar energy per image and a force array of matching shape.
    """
    if not isinstance(n_images, int) or n_images < 5 or not np.isfinite([tolerance, spring, perturbation]).all():
        raise ValueError("Invalid NEB parameters")
    if tolerance <= 0 or spring <= 0 or perturbation < 0 or max_iterations < 1 or climb_after < 0:
        raise ValueError("Invalid NEB parameters")
    ends = np.asarray(endpoints, float)
    if ends.ndim < 2 or ends.shape[0] != 2 or not np.isfinite(ends).all():
        raise ValueError("Two finite endpoints required")
    _, endpoint_force = energy_force(ends)
    endpoint_max = float(np.linalg.norm(np.reshape(endpoint_force, (2, -1)), axis=1).max())
    if endpoint_max > tolerance:
        raise ValueError("Endpoints are not stationary within requested tolerance")
    if initial_band is None:
        q = np.linspace(ends[0], ends[1], n_images)
        rng = np.random.default_rng(seed)
        scale = np.sin(np.linspace(0, np.pi, n_images))[1:-1].reshape((-1,)+(1,)*(q.ndim-1))
        q[1:-1] += perturbation*rng.normal(size=q[1:-1].shape)*scale
    else:
        q = np.array(initial_band, dtype=float, copy=True)
        if q.shape != (n_images,)+ends.shape[1:] or not np.isfinite(q).all() or not np.array_equal(q[[0, -1]], ends):
            raise ValueError("Initial band does not match fixed endpoints")
    velocity = np.zeros_like(q)
    dt, alpha, positive, calls = .025, .1, 0, 0
    history, converged = [], False
    for iteration in range(max_iterations+1):
        energies, true = energy_force(q); calls += n_images
        if np.asarray(energies).shape != (n_images,) or np.asarray(true).shape != q.shape or not np.isfinite(energies).all() or not np.isfinite(true).all():
            raise ValueError("Potential callback returned invalid energies or forces")
        climbing = iteration >= climb_after
        force, ci = projected_forces(q, energies, true, spring, climbing)
        fmax = float(np.linalg.norm(force[1:-1].reshape(n_images-2, -1), axis=1).max())
        fci = float(np.linalg.norm(true[ci]))
        history.append({"iteration": iteration, "climbing": climbing, "climbing_index": ci,
                        "max_NEB_force": fmax, "climbing_true_force": fci,
                        "forward_barrier": float(energies[ci]-energies[0]), "optimizer_step_parameter": dt})
        if climbing and fmax <= tolerance and fci <= tolerance:
            converged = True
            break
        if iteration == max_iterations:
            break
        if iteration == climb_after:
            velocity[:] = 0
            dt, alpha, positive = .025, .1, 0
        velocity += dt*force
        power = float(np.sum(velocity*force))
        if power > 0:
            velocity = (1-alpha)*velocity+alpha*np.linalg.norm(velocity)*force/max(np.linalg.norm(force), 1e-30)
            positive += 1
            if positive > 5:
                dt = min(dt*1.1, .15)
                alpha *= .99
        else:
            velocity[:] = 0
            dt *= .5
            alpha, positive = .1, 0
        displacement = dt*velocity
        max_move = float(np.linalg.norm(displacement[1:-1].reshape(n_images-2, -1), axis=1).max())
        if max_move > .08:
            displacement *= .08/max_move
        q[1:-1] += displacement[1:-1]
    return {"summary": {"n_images": n_images, "seed": seed, "force_tolerance": tolerance,
            "iterations": iteration, "iteration_cap": max_iterations, "converged": converged,
            "stop_reason": "force_tolerance" if converged else "iteration_limit",
            "band_energy_force_point_evaluations": calls, "endpoint_energy_force_point_evaluations": 2,
            "climbing_index": ci, "max_NEB_force": fmax, "climbing_true_force": fci,
            "endpoint_max_force": endpoint_max, "forward_barrier": float(energies[ci]-energies[0]),
            "reverse_barrier": float(energies[ci]-energies[-1]), "reaction_energy": float(energies[-1]-energies[0])},
            "band": q.tolist(), "energies": np.asarray(energies).tolist(), "history": history}


def load_source():
    spec = importlib.util.spec_from_file_location("quantumequi_path_original", SOURCE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def source_replay_audit(output=OUT):
    """Exactly replay the unchanged class in both endpoint directions."""
    import torch
    output = Path(output); output.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    torch.set_num_threads(1)
    module = load_source()
    capture_path = BASE/"results/original/captured_results.json"
    weights_path = BASE/"results/original/untrained_source_weights.pt"
    capture = json.loads(capture_path.read_text(encoding="utf-8"))
    z = torch.tensor(capture["z_atomic_numbers"], dtype=torch.long)
    model = module.EquivariantPotentialModel(num_species=35, hidden_dim=32, num_layers=3)
    model.load_state_dict(torch.load(weights_path, map_location="cpu", weights_only=True))
    model.eval()
    method = module.ClimbingImageNEB.run_neb
    ends = np.array([capture["R_reactant"], capture["R_product"]], float)
    cases, all_rows, saved = [], [], {}
    for direction in ("forward", "reverse"):
        initial, final = ends if direction == "forward" else ends[::-1]
        evaluations, frame_data, snapshots, e_history, f_history = [], {}, [], [], []

        def potential(numbers, coords):
            index = len(evaluations) % 7
            caller = sys._getframe(1)
            if index == 0:
                assert caller.f_code is method.__code__
                snapshots.append(caller.f_locals["images"].copy())
                e_history.append([]); f_history.append([])
            energy, force = model(numbers, coords)
            e_history[-1].append(float(energy.detach()))
            f_history[-1].append(force.detach().numpy().copy())
            evaluations.append({"direction": direction, "iteration": len(evaluations)//7,
                                "image": index, "energy_nominal": float(energy.detach()),
                                "true_max_atomic_force_nominal": float(force.detach().norm(dim=1).max())})
            return energy, force

        def profile(frame, event, arg):
            if event == "return" and frame.f_code is method.__code__:
                frame_data["final_band"] = frame.f_locals["images"].copy()
                frame_data["preupdate_energies"] = frame.f_locals["final_energies"].copy()
                frame_data["climbing_active"] = frame.f_locals["climbing_active"]

        previous = sys.getprofile()
        try:
            sys.setprofile(profile)
            result = module.ClimbingImageNEB(potential, z, k_spring=1.8, num_images=7).run_neb(
                initial, final, max_iterations=30)
        finally:
            sys.setprofile(previous)
        band = frame_data["final_band"]
        fresh_e, fresh_f = [], []
        for coords in band:
            e, f = model(z, torch.tensor(coords, dtype=torch.float32))
            fresh_e.append(float(e.detach())); fresh_f.append(f.detach().numpy())
        fresh_e, fresh_f = np.array(fresh_e), np.array(fresh_f)
        residual, ci = projected_forces(band, fresh_e, fresh_f, spring=1.8, climb=True, source=True)
        actual_ci = result["transition_state_image_idx"]
        sequential_errors, snapshot_differences = [], []
        for step, before in enumerate(snapshots):
            expected = source_update(before, e_history[step], f_history[step], climb=step >= 15)
            next_band = snapshots[step+1] if step < 29 else band
            sequential_errors.append(float(np.max(abs(expected-next_band))))
            simultaneous = source_update(before, e_history[step], f_history[step], climb=step >= 15, sequential=False)
            snapshot_differences.append(float(np.max(abs(expected-simultaneous))))
        minimum_distances = []
        for i, coords in enumerate(band):
            distance = np.linalg.norm(coords[:, None]-coords[None], axis=-1)
            np.fill_diagonal(distance, np.inf)
            minimum_distances.append(float(distance.min()))
            all_rows.append({"direction": direction, "image": i,
                "preupdate_energy_nominal": float(frame_data["preupdate_energies"][i]),
                "postupdate_energy_nominal": float(fresh_e[i]),
                "stale_energy_difference_nominal": float(fresh_e[i]-frame_data["preupdate_energies"][i]),
                "true_max_atomic_force_nominal": float(np.linalg.norm(fresh_f[i], axis=1).max()),
                "projected_max_atomic_force_nominal": float(np.linalg.norm(residual[i], axis=1).max()) if 0 < i < 6 else None,
                "minimum_pair_distance_A": minimum_distances[-1]})
        record = {"direction": direction, "source_return": result,
            "source_energy_force_calls": len(evaluations), "final_band_energy_force_calls": 7,
            "iterations": len(snapshots), "source_reported_candidate_index": actual_ci,
            "fresh_highest_interior_index": ci,
            "source_candidate_true_max_atomic_force_nominal": float(np.linalg.norm(fresh_f[actual_ci], axis=1).max()),
            "final_max_projected_atomic_force_nominal": float(np.linalg.norm(residual[1:-1], axis=2).max()),
            "endpoint_max_atomic_forces_nominal": [float(np.linalg.norm(fresh_f[i], axis=1).max()) for i in (0, 6)],
            "max_stale_energy_difference_nominal": float(np.max(abs(fresh_e-frame_data["preupdate_energies"]))),
            "fresh_forward_barrier_candidate_nominal": float(fresh_e[ci]-fresh_e[0]),
            "fresh_max_all_images_relative_nominal": float(fresh_e.max()-fresh_e[0]),
            "minimum_final_band_pair_distance_A": min(minimum_distances),
            "sequential_update_reconstruction_max_error_A": max(sequential_errors),
            "same_step_snapshot_vs_sequential_max_difference_A": max(snapshot_differences),
            "capture_match": result == capture["neb_results"] if direction == "forward" else None}
        cases.append(record)
        saved[direction] = {"atomic_numbers": z.tolist(), "final_band_A": band.tolist(),
            "candidate_index": actual_ci, "candidate_A": band[actual_ci].tolist(),
            "final_energies_nominal": fresh_e.tolist(), "final_forces_nominal": fresh_f.tolist()}
        np.savez_compressed(output/f"source_{direction}_history.npz",
            preupdate_bands_A=np.array(snapshots), preupdate_energies_nominal=np.array(e_history),
            preupdate_forces_nominal=np.array(f_history), final_band_A=band)
        table(output/f"source_{direction}_evaluations.csv", evaluations)
    forward = np.array(saved["forward"]["final_band_A"])
    reversed_back = np.array(saved["reverse"]["final_band_A"])[::-1]
    dump(output/"source_final_coordinates.json", saved)
    table(output/"source_final_images.csv", all_rows)
    snapshot = output/"source_audit_script_snapshot.py.txt"
    snapshot.write_bytes(Path(__file__).read_bytes())
    report = {"scope": "Exact unchanged random-weight source NEB replay; no chemical potential validation",
        "cases": cases, "reversal_mapped_max_coordinate_difference_A": float(np.max(abs(forward-reversed_back))),
        "reversal_mapped_RMS_coordinate_difference_A": float(np.sqrt(np.mean((forward-reversed_back)**2))),
        "source_has_force_stop": False, "source_is_inplace_sequential": True,
        "source_forces_evaluated_on_preupdate_snapshot": True,
        "source_tangent_uses_already_updated_left_neighbor": True,
        "source_output_energy_one_update_stale": True,
        "source_energies_are_nominal_kcal_mol_not_calibrated": True,
        "counts": {"unchanged_source_path_runs": 2, "source_iterations": 60,
            "source_energy_force_calls": sum(r["source_energy_force_calls"] for r in cases),
            "final_band_energy_force_calls": 14, "new_training_runs": 0, "DFT_calculations": 0},
        "versions": {"python": platform.python_version(), "numpy": np.__version__, "torch": torch.__version__},
        "runtime_seconds": time.perf_counter()-started,
        "source_sha256": {p.relative_to(ROOT).as_posix(): digest(p) for p in (SOURCE, capture_path, weights_path, snapshot)},
        "output_sha256": {p.name: digest(p) for p in sorted(output.glob("source_*")) if p.is_file() and p.name not in ("source_audit.json", snapshot.name)}}
    dump(output/"source_audit.json", report)
    print(json.dumps({"source_counts": report["counts"], "source_cases": cases,
                      "reversal_max_A": report["reversal_mapped_max_coordinate_difference_A"]}, indent=2))
    return report


def run_benchmark(pilot=False):
    started = time.perf_counter()
    output = OUT/"pilot" if pilot else OUT
    output.mkdir(parents=True, exist_ok=True)
    rows, history_rows, image_rows, references = [], [], [], {}
    for name in ("tilted_sine_2d", "periodic_curve_3d"):
        surface = AnalyticSurface(name); ref = surface.references(); references[name] = ref
        for n in ([9] if pilot else [7, 11, 17]):
            for tolerance in ([1e-4] if pilot else [1e-3, 1e-5]):
                for seed in ([11] if pilot else [11, 22, 33]):
                    case = f"{name}_n{n}_tol{tolerance:g}_seed{seed}"
                    result = optimize_neb(surface.energy_force, [ref["left"], ref["right"]],
                                          n_images=n, tolerance=tolerance, seed=seed)
                    s = result["summary"]; q = np.array(result["band"])
                    saddle = q[s["climbing_index"]]
                    eigen = np.linalg.eigvalsh(surface.hessian(saddle))
                    s.update({"case": case, "surface": name, "dimension": surface.dimension,
                        "forward_barrier_error": s["forward_barrier"]-ref["forward_barrier"],
                        "reverse_barrier_error": s["reverse_barrier"]-ref["reverse_barrier"],
                        "saddle_coordinate_error": float(np.linalg.norm(saddle-ref["saddle"])),
                        "negative_hessian_modes": int(np.count_nonzero(eigen < -1e-7)),
                        "min_hessian_eigenvalue": float(eigen.min()), "max_hessian_eigenvalue": float(eigen.max())})
                    rows.append(s)
                    history_rows.extend(dict(case=case, **h) for h in result["history"])
                    for i, point in enumerate(q):
                        image_rows.append({"case": case, "surface": name, "image": i,
                            "x": point[0], "y": point[1], "z": point[2] if surface.dimension == 3 else None,
                            "energy": result["energies"][i]})
    # Reversal uses an explicitly reversed identical initial band; no random
    # perturbation asymmetry or different starting conditions are introduced.
    reversal = []
    for name in references:
        sf = AnalyticSurface(name); ref = references[name]
        initial = np.linspace(ref["left"], ref["right"], 11)
        for direction in ("forward", "reverse"):
            q = initial if direction == "forward" else initial[::-1]
            r = optimize_neb(sf.energy_force, q[[0, -1]], initial_band=q,
                             n_images=11, tolerance=1e-5)
            reversal.append({"surface": name, "direction": direction, **r["summary"], "band": r["band"]})
    pairs = []
    for name in references:
        a, b = [r for r in reversal if r["surface"] == name]
        pairs.append({"surface": name,
            "mapped_max_coordinate_difference": float(np.max(abs(np.array(a["band"])-np.array(b["band"])[::-1]))),
            "forward_minus_reversed_reverse_barrier": a["forward_barrier"]-b["reverse_barrier"],
            "both_converged": a["converged"] and b["converged"]})
    table(output/"analytic_cases.csv", rows)
    table(output/"analytic_force_history.csv", history_rows)
    table(output/"analytic_final_images.csv", image_rows)
    dump(output/"analytic_references.json", references)
    dump(output/"analytic_reversal.json", {"runs": reversal, "comparisons": pairs})
    summary = {"pilot": pilot, "scope": "Two designed dimensionless analytic surfaces; numerical convergence only",
        "method_provenance": "Callback-based generalization of SynthaPore improved-tangent/FIRE-style solver; no new optimizer claimed",
        "surfaces": {
            "tilted_sine_2d": "V=(x*x-1)^2+0.35*x+2.5*(y-0.55*sin(1.7*x))^2",
            "periodic_curve_3d": "V=2*(1-cos(x))+2*(y-0.45*sin(x))^2+3*(z-0.30*(1-cos(x)))^2"},
        "units": "dimensionless coordinates, energy, force; optimizer step parameter is not time",
        "counts": {"sweep_cases": len(rows), "sweep_iterations": sum(r["iterations"] for r in rows),
            "sweep_band_point_evaluations": sum(r["band_energy_force_point_evaluations"] for r in rows),
            "sweep_endpoint_point_evaluations": 2*len(rows), "reversal_cases": len(reversal),
            "reversal_iterations": sum(r["iterations"] for r in reversal),
            "reversal_band_point_evaluations": sum(r["band_energy_force_point_evaluations"] for r in reversal),
            "reversal_endpoint_point_evaluations": 2*len(reversal),
            "sweep_history_rows": len(history_rows), "sweep_final_image_rows": len(image_rows),
            "analytic_Hessians_at_sweep_candidates": len(rows), "reference_Hessians": 6,
            "DFT_calculations": 0, "new_training_runs": 0},
        "checks": {"all_sweep_cases_converged": all(r["converged"] for r in rows),
            "all_sweep_candidates_index_one": all(r["negative_hessian_modes"] == 1 for r in rows),
            "maximum_absolute_forward_barrier_error": max(abs(r["forward_barrier_error"]) for r in rows),
            "reversal": pairs},
        "limits": ["Two low-dimensional designed surfaces do not validate a chemical potential or reaction mechanism",
            "Stationary endpoints are algebraic references; molecular endpoint relaxation is not performed",
            "No free energies, atomistic electrochemistry, quantum barriers, MD, or new trained weights",
            "Finite image discretization affects the band; a converged CI is separately checked against the known saddle"],
        "versions": {"python": platform.python_version(), "numpy": np.__version__},
        "references": REFERENCES, "runtime_seconds": time.perf_counter()-started,
        "source_sha256": {p.relative_to(ROOT).as_posix(): digest(p) for p in (Path(__file__), PREDECESSOR)},
        "output_sha256": {p.name: digest(p) for p in sorted(output.glob("analytic_*")) if p.is_file()}}
    dump(output/"summary.json", summary)
    print(json.dumps({"counts": summary["counts"], "checks": summary["checks"]}, indent=2))
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-audit", action="store_true")
    parser.add_argument("--pilot", action="store_true")
    args = parser.parse_args()
    source_replay_audit() if args.source_audit else run_benchmark(args.pilot)
