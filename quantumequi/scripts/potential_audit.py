"""Audit the exact captured, untrained EGNN; no chemical accuracy labels.

Float64 probes convert the saved float32 weights, holding their represented
values fixed. They do not recover lost precision or constitute training.
"""
import copy
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import time
import numpy as np
import torch

BASE = Path(__file__).resolve().parents[1]
OUT = BASE/'results/potential'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dump(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n', encoding='utf-8', newline='\n')


def table(path, rows):
    with path.open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def source_module():
    path = BASE/'source/quantum_egnn_neb_engine.py'
    spec = importlib.util.spec_from_file_location('quantumequi_source_potential_audit', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_model(dtype=torch.float64):
    model = source_module().EquivariantPotentialModel(num_species=35, hidden_dim=32, num_layers=3)
    state = torch.load(BASE/'results/original/untrained_source_weights.pt', map_location='cpu', weights_only=True)
    model.load_state_dict(state)
    return model.to(dtype=dtype).eval()


def evaluate(model, numbers, coordinates):
    dtype = next(model.parameters()).dtype
    coords = torch.tensor(np.asarray(coordinates), dtype=dtype)
    z = torch.tensor(np.asarray(numbers), dtype=torch.long)
    e, f = model(z, coords)
    return float(e.detach()), f.detach().numpy().astype(float)


def orthogonal(seed, reflection=False):
    q, r = np.linalg.qr(np.random.default_rng(seed).normal(size=(3, 3)))
    q = q @ np.diag(np.sign(np.diag(r)))
    if np.linalg.det(q) < 0:
        q[:, 0] *= -1
    if reflection:
        q[:, 0] *= -1
    return q


def finite_difference(model, numbers, coords, step):
    values = np.zeros_like(coords, dtype=float)
    for j in range(coords.size):
        plus, minus = coords.copy(), coords.copy()
        plus.flat[j] += step
        minus.flat[j] -= step
        ep, _ = evaluate(model, numbers, plus)
        em, _ = evaluate(model, numbers, minus)
        values.flat[j] = -(ep-em)/(2*step)
    return values


def main():
    tic = time.perf_counter()
    torch.set_num_threads(1)
    OUT.mkdir(parents=True, exist_ok=True)
    source = json.loads((BASE/'results/original/captured_results.json').read_text(encoding='utf-8'))
    z = np.array(source['z_atomic_numbers'])
    geometries = {'reactant': np.array(source['R_reactant']),
                  'product': np.array(source['R_product']), 'source_candidate': np.array(source['ts_coords'])}
    rows, differences, baselines, coordinates = [], [], [], []
    evaluations = 0
    models = {}
    for label, dtype in [('float32', torch.float32), ('float64', torch.float64)]:
        model = models[label] = load_model(dtype)
        for geometry, coords in geometries.items():
            energy, force = evaluate(model, z, coords)
            evaluations += 1
            baselines.append(dict(dtype=label, geometry=geometry, nominal_energy=energy,
                                 max_force_component=float(abs(force).max()),
                                 net_force_norm=float(np.linalg.norm(force.sum(0))),
                                 torque_norm=float(np.linalg.norm(np.cross(coords-coords.mean(0), force).sum(0)))))
            for i in range(len(z)):
                coordinates.append(dict(dtype=label, geometry=geometry, atom=i, atomic_number=int(z[i]),
                    x_A=coords[i,0], y_A=coords[i,1], z_A=coords[i,2],
                    force_x_nominal=force[i,0], force_y_nominal=force[i,1], force_z_nominal=force[i,2]))
            for seed in (101, 202, 303, 404):
                for reflection in (False, True):
                    q = orthogonal(seed, reflection)
                    shift = np.array([3.1, -2.2, .7])
                    permutation = np.random.default_rng(seed+17).permutation(len(z))
                    transformed = (coords @ q.T + shift)[permutation]
                    changed_e, changed_f = evaluate(model, z[permutation], transformed)
                    evaluations += 1
                    expected_f = (force @ q.T)[permutation]
                    rows.append(dict(dtype=label, geometry=geometry, seed=seed, reflection=reflection,
                        transform='orthogonal_plus_translation_plus_permutation',
                        energy_abs_error=abs(changed_e-energy), force_max_abs_error=float(abs(changed_f-expected_f).max())))
            for step in (1e-2, 1e-3, 1e-4, 1e-5):
                fd = finite_difference(model, z, coords, step)
                evaluations += 2*coords.size
                differences.append(dict(dtype=label, geometry=geometry, displacement_A=step,
                    force_max_abs_error=float(abs(fd-force).max()), force_RMSE=float(np.sqrt(np.mean((fd-force)**2)))))
    model = models['float64']
    coords = torch.tensor(geometries['reactant'], dtype=torch.float64)
    energy, _ = model(torch.tensor(z), coords)
    evaluations += 1
    gradients = torch.autograd.grad(energy, tuple(model.parameters()), allow_unused=True)
    unused = [dict(name=name, count=p.numel()) for (name,p),g in zip(model.named_parameters(), gradients) if g is None]
    modified = copy.deepcopy(model)
    with torch.no_grad():
        for p in modified.layers[-1].phi_x.parameters():
            p.add_(7.0)
    e0, f0 = evaluate(model, z, geometries['reactant'])
    e1, f1 = evaluate(modified, z, geometries['reactant'])
    evaluations += 2
    table(OUT/'transformations.csv', rows)
    table(OUT/'force_finite_differences.csv', differences)
    table(OUT/'baselines.csv', baselines)
    table(OUT/'coordinates_forces.csv', coordinates)
    double = [r for r in rows if r['dtype']=='float64']
    single = [r for r in rows if r['dtype']=='float32']
    saved_reactant = next(r for r in baselines if r['dtype']=='float32' and r['geometry']=='reactant')
    result = {
        'scope': 'Differentiation and symmetry of exact captured untrained weights, not quantum-force or chemical validation',
        'model_parameter_count': sum(p.numel() for p in model.parameters()),
        'unused_energy_parameters': unused,
        'unused_energy_parameter_count': sum(p['count'] for p in unused),
        'last_coordinate_head_ablation': {'added_to_all_parameters': 7.0,
            'energy_difference': e1-e0, 'force_max_difference': float(abs(f1-f0).max()),
            'interpretation': 'The last updated coordinates are not consumed by the energy head; those coordinate-head parameters are disconnected from energy and force.'},
        'float32_capture_energy_matches': saved_reactant['nominal_energy'] == source['e_val'],
        'float32_max_energy_covariance_error': max(r['energy_abs_error'] for r in single),
        'float32_max_force_covariance_error': max(r['force_max_abs_error'] for r in single),
        'float64_max_energy_covariance_error': max(r['energy_abs_error'] for r in double),
        'float64_max_force_covariance_error': max(r['force_max_abs_error'] for r in double),
        'float64_weight_policy': 'Exact saved float32 parameter values converted to float64; no reinitialization or training',
        'energy_units': 'Nominal kcal/mol labels assigned by source; no calibration',
        'force_units': 'Nominal kcal/(mol Angstrom) labels assigned by source; no calibration',
        'counts': {'geometries': 3, 'dtypes': 2, 'transformation_probes': len(rows),
                   'full_Cartesian_finite_difference_checks': len(differences),
                   'energy_force_model_evaluations': evaluations, 'training_epochs': 0},
        'elapsed_seconds': time.perf_counter()-tic,
        'inputs_sha256': {p.relative_to(BASE).as_posix(): sha(p) for p in [
            BASE/'source/quantum_egnn_neb_engine.py', BASE/'results/original/captured_results.json',
            BASE/'results/original/untrained_source_weights.pt', Path(__file__)]},
        'outputs_sha256': {p.name: sha(p) for p in sorted(OUT.glob('*.csv'))},
        'versions': {'numpy': np.__version__, 'torch': torch.__version__},
    }
    dump(OUT/'summary.json', result)
    print(json.dumps({k:result[k] for k in ['model_parameter_count','unused_energy_parameter_count','counts',
        'float64_max_energy_covariance_error','float64_max_force_covariance_error']}, indent=2))


if __name__ == '__main__':
    main()
