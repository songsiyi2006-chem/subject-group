"""Numerical MD and CI-NEB benchmarks with explicit units and analytic references.

The spring and curved double-well potentials are designed mathematical models,
not trained interatomic potentials or a chemical reaction potential surface.
"""
from __future__ import annotations
import argparse
import ast
import csv
import hashlib
import importlib.metadata
import json
import math
from pathlib import Path
import platform
import sys
import time
from typing import Dict

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parent
EV_J = 1.602176634e-19
AMU_KG = 1.66053906660e-27
MASS_EV_FS2_A2 = AMU_KG * 1e10 / EV_J
KB_EV_K = 1.380649e-23 / EV_J
EV_KCAL_MOL = EV_J * 6.02214076e23 / 4184
SOURCES = [
    "https://journals.aps.org/pr/abstract/10.1103/PhysRev.159.98",
    "https://pure.rug.nl/ws/files/64380902/1.448118.pdf",
    "https://arxiv.org/abs/1203.5428",
    "https://henkelmanlab.org/pubs/henkelman00_9978.pdf",
    "https://doi.org/10.1063/1.1329672",
]


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
                          encoding="utf-8", newline="\n")


def write_csv(path, rows):
    with Path(path).open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def remove_com_velocity(velocities, masses):
    v, m = np.asarray(velocities, float), np.asarray(masses, float)
    if v.ndim != 2 or v.shape[1] != 3 or m.shape != (len(v),) or len(m) < 2:
        raise ValueError("Require N>=2 masses and an (N,3) velocity array")
    if not np.all(np.isfinite(v)) or not np.all(np.isfinite(m)) or np.any(m <= 0):
        raise ValueError("Finite velocities and positive finite masses required")
    return v - np.sum(v * m[:, None], axis=0) / m.sum()


def kinetic_energy(velocities, masses):
    return float(.5 * MASS_EV_FS2_A2 * np.sum(np.asarray(masses)[:, None] * np.asarray(velocities)**2))


def temperature(kinetic_eV, n_atoms, remove_com=True):
    dof = 3 * n_atoms - (3 if remove_com else 0)
    if dof <= 0 or not math.isfinite(kinetic_eV) or kinetic_eV < 0:
        raise ValueError("Positive degrees of freedom and finite nonnegative kinetic energy required")
    return 2 * kinetic_eV / (dof * KB_EV_K)


def harmonic_energy_force(position, spring_eV_A2=2.):
    """Translation-invariant isotropic spring cluster; masses are equal in studies."""
    q = np.asarray(position, float)
    relative = q - q.mean(axis=0)
    return float(.5 * spring_eV_A2 * np.sum(relative**2)), -spring_eV_A2 * relative


def verlet_step(position, velocity, force, masses, timestep_fs, spring_eV_A2=2.):
    mass = np.asarray(masses)[:, None] * MASS_EV_FS2_A2
    next_position = position + velocity * timestep_fs + .5 * force / mass * timestep_fs**2
    energy, next_force = harmonic_energy_force(next_position, spring_eV_A2)
    next_velocity = velocity + .5 * (force + next_force) / mass * timestep_fs
    return next_position, next_velocity, next_force, energy


def simulate_md(*, seed=20260928, timestep_fs=.5, duration_fs=2000., n_atoms=8,
                mass_amu=12.011, spring_eV_A2=2., target_K=300.,
                thermostat="nve", relaxation_fs=50., burn_in_fs=0., sample_every=10):
    if thermostat not in ("nve", "langevin_baoab", "berendsen", "berendsen_source_dof"):
        raise ValueError("Unknown ensemble algorithm")
    if not isinstance(n_atoms, int) or n_atoms < 2 or not isinstance(seed, int):
        raise ValueError("Integer seed and N>=2 required")
    vals = [timestep_fs, duration_fs, mass_amu, spring_eV_A2, target_K, relaxation_fs]
    if not all(math.isfinite(v) and v > 0 for v in vals) or not 0 <= burn_in_fs < duration_fs:
        raise ValueError("Positive finite integration parameters and valid burn-in required")
    steps = int(round(duration_fs / timestep_fs))
    if steps < 1 or not math.isclose(steps * timestep_fs, duration_fs, abs_tol=1e-10) or sample_every < 1:
        raise ValueError("Duration must be an integer number of timesteps; positive stride required")
    omega = math.sqrt(spring_eV_A2 / (mass_amu * MASS_EV_FS2_A2))
    if omega * timestep_fs >= 2:
        raise ValueError("Harmonic Verlet stability condition omega*dt < 2 violated")
    if thermostat.startswith("berendsen") and timestep_fs > relaxation_fs:
        raise ValueError("Berendsen coupling requires dt <= relaxation time")
    rng = np.random.default_rng(seed)
    masses = np.full(n_atoms, mass_amu)
    position = rng.normal(size=(n_atoms, 3)) * math.sqrt(KB_EV_K * target_K / spring_eV_A2)
    position -= position.mean(axis=0)
    velocity = remove_com_velocity(rng.normal(size=(n_atoms, 3)) *
                                  math.sqrt(KB_EV_K * target_K / (mass_amu * MASS_EV_FS2_A2)), masses)
    initial_q, initial_v = position.copy(), velocity.copy()
    potential, force = harmonic_energy_force(position, spring_eV_A2)
    initial_total = potential + kinetic_energy(velocity, masses)
    max_error, max_com, telemetry, production_T, production_K, production_U = 0., 0., [], [], [], []
    c = math.exp(-timestep_fs / relaxation_fs)
    noise_scale = math.sqrt((1-c*c) * KB_EV_K * target_K / (mass_amu * MASS_EV_FS2_A2))
    for step in range(steps + 1):
        if step:
            if thermostat == "langevin_baoab":
                mass = masses[:, None] * MASS_EV_FS2_A2
                velocity += .5 * timestep_fs * force / mass
                position += .5 * timestep_fs * velocity
                noise = rng.normal(size=velocity.shape)
                noise -= noise.mean(axis=0)
                velocity = c * velocity + noise_scale * noise
                position += .5 * timestep_fs * velocity
                potential, force = harmonic_energy_force(position, spring_eV_A2)
                velocity += .5 * timestep_fs * force / mass
            else:
                position, velocity, force, potential = verlet_step(position, velocity, force, masses, timestep_fs, spring_eV_A2)
                if thermostat.startswith("berendsen"):
                    current_T = temperature(kinetic_energy(velocity, masses), n_atoms, thermostat != "berendsen_source_dof")
                    factor2 = 1 + timestep_fs / relaxation_fs * (target_K / current_T - 1)
                    if not factor2 > 0:
                        raise ArithmeticError("Nonpositive thermostat scale")
                    velocity *= math.sqrt(factor2)
        kinetic = kinetic_energy(velocity, masses)
        T = temperature(kinetic, n_atoms, True)
        total = kinetic + potential
        max_error = max(max_error, abs(total-initial_total))
        com_speed = float(np.linalg.norm(np.sum(masses[:, None] * velocity, axis=0) / masses.sum()))
        max_com = max(max_com, com_speed)
        actual_time = step * timestep_fs
        if actual_time > burn_in_fs:
            production_T.append(T)
            production_K.append(kinetic)
            production_U.append(potential)
        if step % sample_every == 0 or step == steps:
            telemetry.append({"step": step, "time_fs": actual_time, "potential_eV": potential,
                              "kinetic_eV": kinetic, "total_eV": total, "temperature_COM_dof_K": T,
                              "temperature_3N_K": temperature(kinetic, n_atoms, False),
                              "COM_speed_A_fs": com_speed})
    exact_q = initial_q * math.cos(omega*duration_fs) + initial_v * math.sin(omega*duration_fs) / omega
    exact_v = initial_v * math.cos(omega*duration_fs) - initial_q * omega * math.sin(omega*duration_fs)
    f = 3*n_atoms-3
    summary = {
        "seed": seed, "thermostat": thermostat, "n_atoms": n_atoms, "degrees_of_freedom": f,
        "timestep_fs": timestep_fs, "duration_fs": duration_fs, "burn_in_fs": burn_in_fs,
        "integration_steps": steps, "energy_force_evaluations": steps+1, "saved_samples": len(telemetry),
        "production_samples_correlated": len(production_T),
        "mass_amu": mass_amu, "spring_eV_A2": spring_eV_A2, "omega_fs_1": omega,
        "target_temperature_K": target_K, "initial_total_eV": initial_total,
        "final_total_eV": total, "max_abs_total_energy_change_eV": max_error,
        "max_relative_total_energy_change": max_error/initial_total,
        "max_COM_speed_A_fs": max_com,
        "final_position_error_A_vs_exact_NVE": float(np.sqrt(np.mean((position-exact_q)**2))) if thermostat == "nve" else None,
        "final_velocity_error_A_fs_vs_exact_NVE": float(np.sqrt(np.mean((velocity-exact_v)**2))) if thermostat == "nve" else None,
        "mean_temperature_K": float(np.mean(production_T)), "variance_temperature_K2": float(np.var(production_T, ddof=1)),
        "mean_kinetic_eV": float(np.mean(production_K)), "variance_kinetic_eV2": float(np.var(production_K, ddof=1)),
        "mean_potential_eV": float(np.mean(production_U)), "variance_potential_eV2": float(np.var(production_U, ddof=1)),
        "canonical_reference_mean_kinetic_eV": .5*f*KB_EV_K*target_K,
        "canonical_reference_variance_kinetic_eV2": .5*f*(KB_EV_K*target_K)**2,
        "canonical_reference_variance_temperature_K2": 2*target_K**2/f,
        "mean_scope": "Finite correlated production trajectory; not independent samples or chemical validation",
    }
    return {"summary": summary, "telemetry": telemetry, "initial_position_A": initial_q.tolist(),
            "initial_velocity_A_fs": initial_v.tolist(), "final_position_A": position.tolist(),
            "final_velocity_A_fs": velocity.tolist()}


def curved_pes(points, curvature=4., bend=.65):
    """V=(x^2-1)^2 + curvature*(y-bend*(1-x^2))^2, in chosen eV/A units."""
    p = np.asarray(points, float)
    if p.shape[-1] != 2 or not np.all(np.isfinite(p)) or curvature <= 0:
        raise ValueError("Finite two-dimensional coordinates and positive curvature required")
    x, y = p[..., 0], p[..., 1]
    delta = y - bend * (1-x*x)
    energy = (x*x-1)**2 + curvature * delta**2
    gradient = np.stack([4*x*(x*x-1) + 4*curvature*bend*x*delta, 2*curvature*delta], axis=-1)
    return energy, -gradient


def pes_hessian(point, curvature=4., bend=.65):
    x, y = np.asarray(point, float)
    delta = y-bend*(1-x*x)
    return np.array([[12*x*x-4+4*curvature*bend*(delta+2*bend*x*x), 4*curvature*bend*x],
                     [4*curvature*bend*x, 2*curvature]])


def minimize_endpoint(guess, tolerance=1e-11, max_iterations=2000):
    point = np.asarray(guess, float).copy()
    evaluations, history = 0, []
    for iteration in range(max_iterations + 1):
        energy, force = curved_pes(point); evaluations += 1
        norm = float(np.linalg.norm(force))
        history.append({"iteration": iteration, "energy_eV": float(energy), "force_norm_eV_A": norm})
        if norm < tolerance:
            break
        if iteration == max_iterations:
            break
        scale = .05
        while scale > 1e-14:
            trial = point + scale*force
            candidate, _ = curved_pes(trial); evaluations += 1
            if candidate <= energy - 1e-4*scale*norm**2:
                point = trial
                break
            scale *= .5
        else:
            break
    eigen = np.linalg.eigvalsh(pes_hessian(point))
    return {"point_A": point.tolist(), "energy_eV": float(energy), "force_norm_eV_A": norm,
            "hessian_eigenvalues_eV_A2": eigen.tolist(), "converged": bool(norm < tolerance and min(eigen) > 0),
            "iterations": iteration, "energy_force_point_evaluations": evaluations, "history": history}


def improved_tangent(previous, current, following, ep, ec, en):
    ahead, behind = following-current, current-previous
    if en > ec > ep:
        tangent = ahead
    elif en < ec < ep:
        tangent = behind
    else:
        big, small = max(abs(en-ec), abs(ep-ec)), min(abs(en-ec), abs(ep-ec))
        tangent = ahead*big + behind*small if en > ep else ahead*small + behind*big
        if np.linalg.norm(tangent) < 1e-14:
            tangent = ahead+behind
    norm = np.linalg.norm(tangent)
    if norm <= 1e-14:
        raise ValueError("Degenerate neighboring images cannot define a tangent")
    return tangent/norm


def neb_forces(band, energies, true_forces, spring=4.5, climb=False, legacy=False):
    forces = np.zeros_like(band)
    ci = 1 + int(np.argmax(energies[1:-1]))
    for i in range(1, len(band)-1):
        if legacy:
            tangent = band[i+1]-band[i-1]
            tangent /= np.linalg.norm(tangent)+1e-8
        else:
            tangent = improved_tangent(band[i-1], band[i], band[i+1], *energies[i-1:i+2])
        parallel = np.dot(true_forces[i], tangent)*tangent
        if climb and i == ci:
            forces[i] = true_forces[i]-2*parallel
        else:
            spacing = np.linalg.norm(band[i+1]-band[i])-np.linalg.norm(band[i]-band[i-1])
            forces[i] = true_forces[i]-parallel+spring*spacing*tangent
    return forces, ci


def optimize_neb(n_images=11, tolerance=1e-5, max_iterations=6000, seed=1,
                 spring=4.5, climb_after=50, endpoints=None, perturbation=.08):
    if not isinstance(n_images, int) or n_images < 5 or tolerance <= 0 or max_iterations < 1 or spring <= 0 or climb_after < 0:
        raise ValueError("Invalid NEB parameters")
    ends = np.array([[-1., 0.], [1., 0.]] if endpoints is None else endpoints, float)
    if ends.shape != (2, 2) or not np.all(np.isfinite(ends)):
        raise ValueError("Two finite 2D endpoints required")
    endpoint_energy, endpoint_force = curved_pes(ends)
    endpoint_norm = float(np.max(np.linalg.norm(endpoint_force, axis=1)))
    if endpoint_norm > tolerance:
        raise ValueError("Endpoints must be minimized before optimizing the band")
    band = np.linspace(ends[0], ends[1], n_images)
    rng = np.random.default_rng(seed)
    band[1:-1] += perturbation*rng.normal(size=band[1:-1].shape)*np.sin(np.linspace(0, math.pi, n_images))[1:-1, None]
    velocity = np.zeros_like(band)
    dt, alpha, positive = .025, .1, 0
    history, evaluations, converged = [], 0, False
    for iteration in range(max_iterations+1):
        energies, true = curved_pes(band); evaluations += n_images
        climb = iteration >= climb_after
        force, ci = neb_forces(band, energies, true, spring, climb)
        max_force = float(np.max(np.linalg.norm(force[1:-1], axis=1)))
        ci_force = float(np.linalg.norm(true[ci]))
        history.append({"iteration": iteration, "climbing": climb, "climbing_index": ci,
                        "barrier_eV": float(energies[ci]-energies[0]),
                        "max_NEB_force_eV_A": max_force, "climbing_true_force_eV_A": ci_force,
                        "optimizer_step_parameter": dt})
        if climb and max_force <= tolerance and ci_force <= tolerance:
            converged = True
            break
        if iteration == max_iterations:
            break
        if iteration == climb_after:
            velocity[:] = 0; dt, alpha, positive = .025, .1, 0
        velocity += dt*force
        power = float(np.sum(velocity*force))
        if power > 0:
            vnorm, fnorm = np.linalg.norm(velocity), np.linalg.norm(force)
            velocity = (1-alpha)*velocity + alpha*vnorm*force/max(fnorm, 1e-30)
            positive += 1
            if positive > 5:
                dt = min(dt*1.1, .15)
                alpha *= .99
        else:
            velocity[:] = 0
            dt *= .5
            alpha, positive = .1, 0
        displacement = dt*velocity
        max_move = float(np.max(np.linalg.norm(displacement[1:-1], axis=1)))
        if max_move > .08:
            displacement *= .08/max_move
        band[1:-1] += displacement[1:-1]
    eigen = np.linalg.eigvalsh(pes_hessian(band[ci]))
    return {"summary": {
        "n_images": n_images, "seed": seed, "tolerance_eV_A": tolerance,
        "iteration_cap": max_iterations, "iterations_executed": iteration,
        "converged": converged, "stop_reason": "force_tolerance" if converged else "iteration_limit",
        "energy_force_point_evaluations": evaluations, "climbing_index": ci,
        "max_NEB_force_eV_A": max_force, "climbing_true_force_eV_A": ci_force,
        "endpoint_max_force_eV_A": endpoint_norm, "barrier_eV": float(energies[ci]-energies[0]),
        "barrier_error_eV": float(energies[ci]-energies[0]-1),
        "saddle_coordinate_error_A": float(np.linalg.norm(band[ci]-np.array([0., .65]))),
        "saddle_negative_hessian_eigenvalues": int(np.count_nonzero(eigen < -1e-6)),
        "saddle_hessian_eigenvalues_eV_A2": eigen.tolist(), "potential_energy_not_free_energy": True},
        "band_A": band.tolist(), "energies_eV": energies.tolist(), "history": history}


def source_neb_probe(n_images=7, iterations=45):
    """Execute the unchanged source class against the same analytic PES.

    A return-frame trace captures its post-update band without editing the class.
    Its final energy variable is also retained to diagnose one-update staleness.
    """
    import torch
    source_path = ROOT / "source/synthapore_engine.py"
    tree = ast.parse(source_path.read_text(encoding="utf-8"))
    cls = next(node for node in tree.body if isinstance(node, ast.ClassDef) and node.name == "ClimbingImageNEB")
    env = {"np": np, "torch": torch, "Dict": Dict, "NeuralInteratomicPotential": object}
    exec(compile(ast.Module(body=[cls], type_ignores=[]), "archived_source_ClimbingImageNEB", "exec"), env)

    class AnalyticPotential:
        calls = 0
        def __call__(self, atomic_numbers, position, edge_index):
            self.calls += 1
            x, y = position[0, 0], position[0, 1]
            delta = y-.65*(1-x*x)
            energy = (x*x-1)**2+4*delta**2
            fx = -(4*x*(x*x-1)+16*.65*x*delta)
            fy = -8*delta
            return energy, torch.stack([fx, fy, x*0]).reshape(1, 3)

    capture = {}
    method = env["ClimbingImageNEB"].optimize_pathway
    def trace(frame, event, arg):
        if frame.f_code is method.__code__ and event == "return":
            capture["band"] = np.stack([x.detach().numpy()[0, :2] for x in frame.f_locals["band_coords"]])
            capture["last_preupdate_energies"] = np.array(frame.f_locals["current_energies"])
        return trace
    potential = AnalyticPotential()
    old_trace = sys.gettrace()
    try:
        sys.settrace(trace)
        output = env["ClimbingImageNEB"](potential, num_images=n_images, spring_k=4.5).optimize_pathway(
            torch.tensor([6]), torch.tensor([[-1., 0., 0.]], dtype=torch.float64),
            torch.tensor([[1., 0., 0.]], dtype=torch.float64), torch.empty((2, 0), dtype=torch.long),
            max_iterations=iterations)
    finally:
        sys.settrace(old_trace)
    band = capture["band"]
    final_energy, true = curved_pes(band)
    residual, ci = neb_forces(band, final_energy, true, 4.5, iterations-1 > 10, legacy=True)
    return {"n_images": n_images, "iterations_executed": iterations, "source_reported_status": output["status"],
            "source_point_energy_force_calls": potential.calls,
            "source_rounded_barrier_kcal_mol": output["calculated_activation_barrier_kcal_mol"],
            "source_preupdate_barrier_eV": float(max(capture["last_preupdate_energies"])-capture["last_preupdate_energies"][0]),
            "fresh_postupdate_barrier_eV": float(max(final_energy)-final_energy[0]),
            "maximum_stale_energy_difference_eV": float(np.max(abs(final_energy-capture["last_preupdate_energies"]))),
            "fresh_max_NEB_force_eV_A": float(np.max(np.linalg.norm(residual[1:-1], axis=1))),
            "fresh_climbing_true_force_eV_A": float(np.linalg.norm(true[ci])),
            "passes_independent_1e_5_force_check": bool(np.max(np.linalg.norm(residual[1:-1], axis=1)) < 1e-5 and np.linalg.norm(true[ci]) < 1e-5),
            "band_A": band.tolist(), "fresh_energies_eV": final_energy.tolist(),
            "original_energy_conversion_eV_to_kcal_mol": 23.0605,
            "scope": "Unchanged source class on analytic benchmark; no neural potential or chemical barrier"}


def run_study(output_dir=None, pilot=False):
    out = Path(output_dir) if output_dir else ROOT / "results/dynamics" / ("pilot" if pilot else "")
    out.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    nve, thermostat_rows, md_details, neb_rows, histories, bands, source_rows = [], [], [], [], [], [], []
    seed_values = [101] if pilot else [101, 102, 103, 104]
    dt_values = [1.] if pilot else [.25, .5, 1., 2., 4.]
    for seed in seed_values:
        for dt in dt_values:
            result = simulate_md(seed=seed, timestep_fs=dt, duration_fs=100. if pilot else 2000.)
            summary = result["summary"]
            nve.append(summary)
            md_details.append({"case": f"nve_s{seed}_dt{dt}", **{k: v for k, v in result.items() if k != "telemetry"}})
            if seed == seed_values[0]:
                write_csv(out / f"nve_dt{dt:g}_trajectory.csv", result["telemetry"])
    print(f"NVE trajectories complete: {len(nve)}", flush=True)
    for method in ["langevin_baoab", "berendsen", "berendsen_source_dof"]:
        for seed in ([201] if pilot else range(201, 209)):
            result = simulate_md(seed=seed, thermostat=method, duration_fs=100. if pilot else 6000.,
                                 burn_in_fs=20. if pilot else 1000.)
            thermostat_rows.append(result["summary"])
            if seed == 201:
                write_csv(out / f"{method}_trajectory.csv", result["telemetry"])
    print(f"Thermostat trajectories complete: {len(thermostat_rows)}", flush=True)
    endpoints = [minimize_endpoint(guess) for guess in [(-1.2, .15), (1.15, -.18)]]
    dump(out / "endpoint_minimization.json", endpoints)
    if not all(e["converged"] for e in endpoints):
        raise ArithmeticError("Endpoint minimization failed")
    ends = [e["point_A"] for e in endpoints]
    for images in ([7] if pilot else [7, 11, 17]):
        for tol in ([1e-4] if pilot else [1e-3, 1e-5]):
            for seed in ([1] if pilot else [1, 2, 3]):
                result = optimize_neb(n_images=images, tolerance=tol, seed=seed, endpoints=ends)
                case = f"images{images}_tol{tol:g}_seed{seed}"
                neb_rows.append({"case": case, **result["summary"]})
                histories.extend({"case": case, **r} for r in result["history"])
                for index, (point, energy) in enumerate(zip(result["band_A"], result["energies_eV"])):
                    bands.append({"case": case, "image_index": index, "x_A": point[0], "y_A": point[1], "energy_eV": energy})
    print(f"Reviewed NEB cases complete: {len(neb_rows)}", flush=True)
    for images in ([7] if pilot else [7, 11, 17]):
        for iterations in ([1, 45] if pilot else [1, 10, 45, 150]):
            source_rows.append(source_neb_probe(images, iterations))
    dump(out / "source_neb_probes.json", source_rows)
    source_table = [{k: v for k, v in r.items() if k not in ("band_A", "fresh_energies_eV")} for r in source_rows]
    write_csv(out / "nve_timestep_summary.csv", nve)
    write_csv(out / "thermostat_replicas.csv", thermostat_rows)
    write_csv(out / "neb_summary.csv", [{**r, "saddle_hessian_eigenvalues_eV_A2": json.dumps(r["saddle_hessian_eigenvalues_eV_A2"])} for r in neb_rows])
    write_csv(out / "neb_force_history.csv", histories)
    write_csv(out / "neb_final_bands.csv", bands)
    write_csv(out / "source_neb_summary.csv", source_table)
    dump(out / "nve_initial_final_states.json", md_details)
    convergence = []
    if not pilot:
        for seed in seed_values:
            group = sorted([r for r in nve if r["seed"] == seed], key=lambda r: r["timestep_fs"])
            for fine, coarse in zip(group, group[1:]):
                convergence.append({"seed": seed, "fine_dt_fs": fine["timestep_fs"], "coarse_dt_fs": coarse["timestep_fs"],
                    "energy_error_order": math.log(coarse["max_abs_total_energy_change_eV"]/fine["max_abs_total_energy_change_eV"], 2),
                    "position_error_order": math.log(coarse["final_position_error_A_vs_exact_NVE"]/fine["final_position_error_A_vs_exact_NVE"], 2)})
        write_csv(out / "nve_convergence_orders.csv", convergence)
    thermostat_grouped = []
    for method in ["langevin_baoab", "berendsen", "berendsen_source_dof"]:
        group = [r for r in thermostat_rows if r["thermostat"] == method]
        thermostat_grouped.append({"thermostat": method, "replicas": len(group),
            "replicate_mean_temperature_K": float(np.mean([r["mean_temperature_K"] for r in group])),
            "mean_within_trajectory_temperature_variance_K2": float(np.mean([r["variance_temperature_K2"] for r in group])),
            "canonical_reference_variance_temperature_K2": group[0]["canonical_reference_variance_temperature_K2"],
            "variance_ratio_to_canonical": float(np.mean([r["variance_temperature_K2"] for r in group])) / group[0]["canonical_reference_variance_temperature_K2"],
            "replicate_mean_kinetic_eV": float(np.mean([r["mean_kinetic_eV"] for r in group])),
            "replicate_mean_potential_eV": float(np.mean([r["mean_potential_eV"] for r in group])),
            "uncertainty_scope": "Correlated finite trajectories; replica summaries are diagnostic, no calibrated confidence interval"})
    write_csv(out / "thermostat_summary.csv", thermostat_grouped)
    summary = {
        "schema_version": 1, "pilot": pilot,
        "evidence": "Numerical algorithms on designed analytic potentials; no chemical reaction barrier or trained MLIP",
        "units": {"length": "angstrom", "time": "fs", "energy": "eV", "force": "eV/angstrom", "mass": "amu",
                  "mass_conversion_eV_fs2_A2_per_amu": MASS_EV_FS2_A2, "kB_eV_K": KB_EV_K,
                  "eV_to_kcal_mol": EV_KCAL_MOL, "mass_constant_kg_per_amu": AMU_KG},
        "counts": {"NVE_trajectories": len(nve), "thermostat_trajectories": len(thermostat_rows),
                   "MD_integration_steps": sum(r["integration_steps"] for r in nve+thermostat_rows),
                   "MD_energy_force_evaluations": sum(r["energy_force_evaluations"] for r in nve+thermostat_rows),
                   "reviewed_NEB_cases": len(neb_rows), "reviewed_NEB_iterations": sum(r["iterations_executed"] for r in neb_rows),
                   "reviewed_NEB_energy_force_point_evaluations": sum(r["energy_force_point_evaluations"] for r in neb_rows),
                   "reviewed_NEB_endpoint_check_point_evaluations": 2*len(neb_rows),
                   "source_NEB_cases": len(source_rows), "source_NEB_iterations": sum(r["iterations_executed"] for r in source_rows),
                   "source_NEB_point_energy_force_calls": sum(r["source_point_energy_force_calls"] for r in source_rows),
                   "source_NEB_postupdate_audit_point_evaluations": sum(r["n_images"] for r in source_rows),
                   "endpoint_optimizations": 2, "endpoint_energy_force_point_evaluations": sum(r["energy_force_point_evaluations"] for r in endpoints),
                   "pilot_excluded_from_main_counts": True, "tests_excluded_from_study_counts": True, "DFT_calculations": 0, "experimental_runs": 0},
        "harmonic_reference": {"potential": "0.5*k*sum_i|r_i-mean(r)|^2", "k_eV_A2": 2., "equal_mass_amu": 12.011,
                               "atoms": 8, "removed_dof": "3 center-of-mass translations; no rotational or bond constraints",
                               "normal_mode_frequency_fs_1": nve[0]["omega_fs_1"], "NVE_seeds": seed_values},
        "NVE_convergence_orders": convergence,
        "thermostat_summary": thermostat_grouped,
        "NEB_reference": {"potential": "(x^2-1)^2 + 4*(y-0.65*(1-x^2))^2",
                          "known_minima_A": [[-1., 0.], [1., 0.]], "known_saddle_A": [0., .65],
                          "known_barrier_eV": 1., "saddle_Hessian_eigenvalues_eV_A2": [-4., 8.],
                          "valley_is_not_assumed_to_be_MEP": True, "optimizer": "FIRE-style inertial relaxation of NEB projected force",
                          "tangent": "Henkelman-Jonsson energy-weighted tangent", "climb_after_iterations": 50,
                          "spring_eV_A2": 4.5, "fixed_endpoint_policy": "minimized first; reject endpoints above force tolerance"},
        "NEB_cases": neb_rows,
        "checks": {"all_reviewed_NEB_converged": all(r["converged"] for r in neb_rows),
                   "all_reviewed_saddles_have_one_negative_mode": all(r["saddle_negative_hessian_eigenvalues"] == 1 for r in neb_rows),
                   "all_endpoints_converged": all(r["converged"] for r in endpoints),
                   "max_MD_COM_speed_A_fs": max(r["max_COM_speed_A_fs"] for r in nve+thermostat_rows),
                   "source_false_convergence_declarations": sum(not r["passes_independent_1e_5_force_check"] for r in source_rows)},
        "source_MD_audit": ["Source conversion eV/A to N is dimensionally correct, but the random MLIP energy scale is untrained",
            "Source removes COM initially but uses laboratory-frame kinetic energy divided by 3N; internal thermal temperature needs COM-subtracted kinetic energy and 3N-3, especially if an external net force restores drift",
            "Source records updated coordinates at step*dt rather than (step+1)*dt",
            "Source records pre-rescaling kinetic energy and temperature although velocity has already been rescaled",
            "Berendsen weak coupling does not establish canonical fluctuations or equilibrium"],
        "limitations": ["Equal-mass harmonic modes are a numerical benchmark, not a molecular force field",
            "Langevin sampling has finite timestep bias and serial correlations; variance comparisons are descriptive",
            "No full constrained dynamics, rigid rotation removal, electrochemical bath, or trained neural forces",
            "Two-dimensional NEB convergence is not a molecular transition-state or free-energy validation",
            "Image refinement tests this designed PES only; unknown higher-dimensional saddles may be missed"],
        "sources": SOURCES, "runtime_seconds": time.perf_counter()-started,
        "versions": {"python": platform.python_version(), "numpy": np.__version__,
                     "torch_for_unchanged_source_class": importlib.metadata.version("torch")},
        "source_sha256": {p.relative_to(REPO).as_posix(): sha(p) for p in
            [Path(__file__), ROOT/"source/synthapore_engine.py", REPO/"tests/test_synthapore_dynamics.py"] if p.exists()},
        "output_sha256": {p.name: sha(p) for p in sorted(out.iterdir()) if p.suffix in (".csv", ".json") and p.name != "summary.json"},
    }
    dump(out/"summary.json", summary)
    print(json.dumps({"counts": summary["counts"], "checks": summary["checks"], "runtime_seconds": summary["runtime_seconds"]}, indent=2))
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pilot", action="store_true")
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    run_study(args.output_dir, args.pilot)
