"""Replay the separate audit of EXACT captured source weights, without training.

Kept beside its results to preserve the already frozen denoising-study source.
"""
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import numpy as np
import torch

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
SOURCE = ROOT / 'synthapore/source/synthapore_engine.py'
WEIGHTS = ROOT / 'synthapore/results/compatibility/untrained_source_weights.pt'
CAPTURE = ROOT / 'synthapore/results/compatibility/captured_results.json'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    torch.set_num_threads(1)
    spec = importlib.util.spec_from_file_location('captured_source_classes', SOURCE)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    captured = json.loads(CAPTURE.read_text(encoding='utf-8'))
    saved = torch.load(WEIGHTS, weights_only=True, map_location='cpu')
    model = module.NeuralInteratomicPotential(num_species=30, hidden_dim=48, num_layers=3)
    model.load_state_dict(saved); model.eval()
    z = torch.tensor(captured['atomic_numbers'], dtype=torch.long)
    x = torch.tensor(captured['r_init'], dtype=torch.float32)
    edges = torch.tensor(captured['edge_index'], dtype=torch.long)
    field = torch.tensor(captured['e_field'], dtype=torch.float32)
    shift = torch.tensor([4., -3., 7.])
    count = 0

    def evaluate(xx, zz=z, ee=field):
        nonlocal count
        count += 1
        energy, forces = model(zz, xx.detach().clone().requires_grad_(True), edges, ee)
        return float(energy.detach()), forces.detach().double().numpy()

    base, f0 = evaluate(x)
    results = []
    rng = np.random.default_rng(6221)
    for reflected in [False, True]:
        for index in range(4):
            r, _ = np.linalg.qr(rng.normal(size=(3, 3)))
            if np.linalg.det(r) < 0:
                r[:, 0] *= -1
            if reflected:
                r[:, 0] *= -1
            R = torch.tensor(r, dtype=torch.float32)
            e, f = evaluate(x @ R.T, ee=R @ field)
            results.append({'probe': 'reflection' if reflected else 'rotation', 'replica': index,
                            'energy_abs_error': abs(e - base), 'force_max_abs_error': float(np.max(abs(f - f0 @ r.T)))})
    permutation = torch.tensor([7, 2, 4, 0, 1, 5, 6, 3])
    e, f = evaluate(x[permutation], z[permutation])
    results.append({'probe': 'permutation', 'replica': 0, 'energy_abs_error': abs(e - base),
                    'force_max_abs_error': float(np.max(abs(f - f0[permutation])))})
    translated_energy, translated_force = evaluate(x + shift)
    results.append({'probe': 'translation_with_field', 'replica': 0,
                    'energy_abs_error': abs(translated_energy - base),
                    'force_max_abs_error': float(np.max(abs(translated_force - f0)))})
    nofield, nf = evaluate(x, ee=None); translated_nofield, ntf = evaluate(x + shift, ee=None)
    results.append({'probe': 'translation_without_field', 'replica': 0,
                    'energy_abs_error': abs(translated_nofield - nofield),
                    'force_max_abs_error': float(np.max(abs(ntf - nf)))})
    finite = []
    for step in [.01, .001, .0001]:
        estimate = np.zeros((len(z), 3))
        for i in range(len(z)):
            for axis in range(3):
                delta = torch.zeros_like(x); delta[i, axis] = step
                plus, _ = evaluate(x + delta); minus, _ = evaluate(x - delta)
                estimate[i, axis] = -(plus - minus) / (2 * step)
        finite.append({'step': step, 'force_max_abs_error': float(np.max(abs(estimate - f0))),
                       'force_RMS_error': float(np.sqrt(np.mean((estimate - f0) ** 2)))})
    xx = x.detach().clone().requires_grad_(True)
    h, current = model.embedding(z), xx
    for layer in model.layers:
        h, current = layer(h, current, edges)
    charges = model.charge_head(h).reshape(-1)
    gradient = torch.autograd.grad(charges.sum(), xx)[0].detach().double().numpy()
    identity_error = float(np.max(abs(translated_force - f0 - float((field * shift).sum()) * gradient)))
    pairh = model.embedding(z[[0, 7]])
    far = torch.tensor([[0., 0., 0.], [8., 0., 0.]])
    pair_edges = torch.tensor([[0, 1], [1, 0]], dtype=torch.long)
    connected_h, connected_x = model.layers[0](pairh, far, pair_edges)
    isolated_h, isolated_x = model.layers[0](pairh, far, torch.empty((2, 0), dtype=torch.long))
    try:
        model.double()(z, x.double().requires_grad_(True), edges, field.double())
        dtype_failure = None
    except RuntimeError as error:
        dtype_failure = str(error)
    for filename, rows in [('captured_source_transformations.csv', results), ('captured_source_finite_differences.csv', finite)]:
        with (OUT / filename).open('w', encoding='utf-8', newline='') as handle:
            writer = csv.DictWriter(handle, list(rows[0]), lineterminator='\n'); writer.writeheader(); writer.writerows(rows)
    summary = {'scope': 'EXACT captured untrained source weights and source initial typed coordinate array; no MD/NEB rerun',
               'num_species': 30, 'hidden_dim': 48, 'layers': 3,
               'model_parameter_count': sum(p.numel() for p in model.parameters()),
               'embedding_repair_note': 'Changing embedding size from 10 to 30 changes initialization RNG consumption. This audit loads the actual captured repaired-run weights, not a separately reseeded approximation.',
               'base_energy_arbitrary_units': base, 'field_vector': field.tolist(), 'translation': shift.tolist(),
               'unconstrained_total_charge': float(charges.detach().sum()),
               'per_atom_charges': charges.detach().tolist(), 'total_charge_gradient': gradient.tolist(),
               'translated_energy': translated_energy, 'observed_energy_shift': translated_energy - base,
               'expected_shift_minus_Q_E_dot_t': -float(charges.detach().sum()) * float((field * shift).sum()),
               'origin_force_identity_max_error': identity_error,
               'origin_force_change_max': float(np.max(abs(translated_force - f0))),
               'cutoff8_scalar_change': float((connected_h - isolated_h).detach().abs().max()),
               'cutoff8_coordinate_change': float((connected_x - isolated_x).detach().abs().max()),
               'double_failure': dtype_failure,
               'counts': {'energy_force_forwards': count, 'charge_only_forwards': 1, 'extra_layer_forwards': 2,
                          'expected_double_failure_probe': 1, 'transform_rows': len(results), 'finite_difference_rows': len(finite),
                          'training_runs': 0, 'MD_trajectories': 0, 'NEB_runs': 0},
               'source_sha256': {p.relative_to(ROOT).as_posix(): digest(p) for p in [SOURCE, WEIGHTS, CAPTURE, Path(__file__)]},
               'output_sha256': {name: digest(OUT / name) for name in ['captured_source_transformations.csv', 'captured_source_finite_differences.csv']}}
    (OUT / 'captured_source_audit.json').write_text(json.dumps(summary, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == '__main__':
    main()
