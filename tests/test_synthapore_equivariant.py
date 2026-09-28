"""Small numerical invariants and archived-study checks, without training."""
import csv
import hashlib
import json
from pathlib import Path
import unittest

import numpy as np
import torch

from synthapore.scripts.equivariant_reviewed import (
    OUT, ChargeConstrainedPotential, EquivariantLayer, GeometryDenoiser,
    SmoothRadial, cycle_laplacian, orthogonal, pca_plane, synthetic_dataset,
)


class TestSynthaPoreEquivariant(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        torch.set_num_threads(1)

    def setUp(self):
        torch.manual_seed(715)
        self.x = torch.randn(5, 3, dtype=torch.float64)
        self.z = torch.tensor([6, 1, 8, 7, 29])
        self.field = torch.tensor([.2, -.1, .04], dtype=torch.float64)

    def test_radial_cutoff_value_and_first_derivative(self):
        layer = SmoothRadial().double()
        d = torch.tensor([6., 6.001, 9.], dtype=torch.float64, requires_grad=True)
        _, gate = layer(d)
        derivative = torch.autograd.grad(gate.sum(), d)[0]
        torch.testing.assert_close(gate, torch.zeros_like(gate), atol=0, rtol=0)
        torch.testing.assert_close(derivative, torch.zeros_like(gate), atol=0, rtol=0)

    def test_full_message_and_coordinate_cutoff(self):
        layer = EquivariantLayer(8).double()
        h = torch.randn(1, 2, 8, dtype=torch.float64)
        x = torch.tensor([[[0., 0, 0], [8., 0, 0]]], dtype=torch.float64)
        a, xa = layer(h, x)
        b, xb = layer(h, x, torch.zeros(2, 2, dtype=torch.float64))
        torch.testing.assert_close(a, b, atol=0, rtol=0)
        torch.testing.assert_close(xa, xb, atol=0, rtol=0)

    def test_double_dtype_and_finite_training_gradients(self):
        model = GeometryDenoiser().double()
        x = torch.randn(2, 8, 3, dtype=torch.float64)
        pred = model(x, torch.tensor([.1, .2], dtype=torch.float64))
        self.assertEqual(pred.dtype, torch.float64)
        (pred ** 2).mean().backward()
        for layer in model.layers:
            grad = layer.coordinate[0].weight.grad
            self.assertTrue(torch.isfinite(grad).all())
            self.assertGreater(float(grad.abs().sum()), 0)
            self.assertEqual(layer.radial.beta.dtype, torch.float64)

    def test_layer_rotations_reflections_translations(self):
        layer = EquivariantLayer(8).double()
        h = torch.randn(1, 5, 8, dtype=torch.float64)
        x = self.x[None]
        R = torch.tensor(orthogonal(np.random.default_rng(7), True), dtype=torch.float64)
        shift = torch.tensor([3., -1, 2], dtype=torch.float64)
        a, xa = layer(h, x)
        b, xb = layer(h, x @ R.T + shift)
        torch.testing.assert_close(a, b, rtol=1e-12, atol=1e-12)
        torch.testing.assert_close(xb, xa @ R.T + shift, rtol=1e-12, atol=1e-12)

    def test_layer_general_graph_permutation(self):
        layer = EquivariantLayer(8).double()
        h = torch.randn(1, 5, 8, dtype=torch.float64)
        adjacency = torch.tensor([[0, 1, 0, 1, 0], [1, 0, 1, 0, 0], [0, 1, 0, 1, 1],
                                  [1, 0, 1, 0, 0], [0, 0, 1, 0, 0]], dtype=torch.float64)
        order = torch.tensor([3, 0, 4, 1, 2])
        a, xa = layer(h, self.x[None], adjacency)
        b, xb = layer(h[:, order], self.x[order][None], adjacency[order][:, order])
        torch.testing.assert_close(b, a[:, order], atol=1e-12, rtol=1e-12)
        torch.testing.assert_close(xb, xa[:, order], atol=1e-12, rtol=1e-12)

    def test_denoiser_permutation_with_relabelled_adjacency(self):
        model = GeometryDenoiser().double()
        x = torch.randn(2, 8, 3, dtype=torch.float64)
        sigma = torch.tensor([.08, .24], dtype=torch.float64)
        order = torch.tensor([3, 1, 5, 7, 0, 2, 4, 6])
        a = np.zeros((8, 8))
        for i in range(8):
            a[i, (i - 1) % 8] = a[i, (i + 1) % 8] = 1
        adjacency = torch.tensor(a, dtype=torch.float64)
        pred = model(x, sigma)
        permuted = model(x[:, order], sigma, adjacency[order][:, order])
        torch.testing.assert_close(permuted, pred[:, order], atol=1e-12, rtol=1e-12)

    def test_denoiser_preserves_observed_centroid(self):
        model = GeometryDenoiser().double()
        x = torch.randn(3, 8, 3, dtype=torch.float64)
        y = model(x, torch.tensor([.08, .16, .24], dtype=torch.float64))
        torch.testing.assert_close(y.mean(1), x.mean(1), atol=1e-12, rtol=1e-12)

    def test_neutral_charge_and_origin_invariance(self):
        model = ChargeConstrainedPotential().double()
        e, f, q = model(self.z, self.x, self.field)
        et, ft, qt = model(self.z, self.x + torch.tensor([8., -2., 7.]), self.field)
        self.assertLess(abs(float(q.sum())), 1e-12)
        self.assertLess(abs(float(qt.sum())), 1e-12)
        torch.testing.assert_close(e, et, atol=1e-12, rtol=1e-12)
        torch.testing.assert_close(f, ft, atol=1e-12, rtol=1e-12)

    def test_fixed_charged_system_translation_law_and_net_force(self):
        Q = 1.25
        model = ChargeConstrainedPotential(total_charge=Q).double()
        shift = torch.tensor([8., -2., 7.], dtype=torch.float64)
        e, f, q = model(self.z, self.x, self.field)
        et, ft, qt = model(self.z, self.x + shift, self.field)
        self.assertAlmostEqual(float(q.sum()), Q, places=12)
        torch.testing.assert_close(et - e, -Q * (shift * self.field).sum(), atol=1e-12, rtol=1e-12)
        torch.testing.assert_close(f, ft, atol=1e-12, rtol=1e-12)
        torch.testing.assert_close(f.sum(0), Q * self.field, atol=1e-12, rtol=1e-12)

    def test_potential_force_joint_field_covariance(self):
        model = ChargeConstrainedPotential().double()
        R = torch.tensor(orthogonal(np.random.default_rng(15), True), dtype=torch.float64)
        e, f, _ = model(self.z, self.x, self.field)
        er, fr, _ = model(self.z, self.x @ R.T, R @ self.field)
        torch.testing.assert_close(e, er, atol=1e-12, rtol=1e-12)
        torch.testing.assert_close(fr, f @ R.T, atol=1e-12, rtol=1e-12)

    def test_force_central_difference(self):
        model = ChargeConstrainedPotential().double()
        _, forces, _ = model(self.z, self.x, self.field)
        displacement = torch.zeros_like(self.x); displacement[2, 1] = 1e-4
        eplus, _ = model.energy_charges(self.z, self.x + displacement, self.field)
        eminus, _ = model.energy_charges(self.z, self.x - displacement, self.field)
        difference = -(eplus - eminus) / 2e-4
        self.assertLess(abs(float(difference.detach() - forces[2, 1])), 1e-8)

    def test_invalid_species_shapes_fields_and_noise(self):
        model = ChargeConstrainedPotential().double()
        for z in [torch.tensor([0, 1, 8, 7, 29]), torch.tensor([119, 1, 8, 7, 29])]:
            with self.assertRaises(ValueError):
                model(z, self.x, self.field)
        with self.assertRaises(ValueError):
            model(self.z, self.x, torch.ones(4))
        denoiser = GeometryDenoiser()
        with self.assertRaises(ValueError):
            denoiser(torch.zeros(2, 8, 3), torch.tensor([0., .2]))
        with self.assertRaises(ValueError):
            denoiser(torch.full((2, 8, 3), float('nan')), torch.tensor([.1, .2]))
        with self.assertRaises(ValueError):
            denoiser(torch.zeros(2, 8, 3), torch.tensor([.1, .2]), torch.full((8, 8), -1.))

    def test_dataset_disjoint_shapes_and_centered_noise(self):
        records = synthetic_dataset()
        self.assertEqual(len(records), 368)
        self.assertEqual(len({r['shape_id'] for r in records}), 368)
        fingerprints = set()
        for r in records:
            clean, noisy = np.array(r['clean']), np.array(r['noisy'])
            np.testing.assert_allclose(noisy.mean(0), clean.mean(0), atol=1e-14)
            distances = np.sort(np.linalg.norm(clean[:, None] - clean[None], axis=-1).ravel())
            fingerprints.add(tuple(np.round(distances, 9)))
        self.assertEqual(len(fingerprints), 368)

    def test_PCA_projection_reproduces_already_planar_points(self):
        x = np.array([[[0., 0., 2.], [1., 0., 2.], [0., 1., 2.], [-1., 1., 2.]]])
        np.testing.assert_allclose(pca_plane(x), x, atol=1e-12)
        # A rotation does not select a laboratory-axis plane.
        R = orthogonal(np.random.default_rng(99), True)
        np.testing.assert_allclose(pca_plane(x @ R.T), pca_plane(x) @ R.T, atol=1e-12)

    def test_saved_smoother_fits_training_only(self):
        data = synthetic_dataset()
        saved = json.loads((OUT / 'baseline_coefficients.json').read_text())
        for sigma in [.08, .16, .24]:
            rows = [r for r in data if r['split'] == 'train' and r['noise_sigma'] == sigma]
            noisy = np.array([r['noisy'] for r in rows]); clean = np.array([r['clean'] for r in rows])
            lap = cycle_laplacian(noisy)
            alpha = float(np.sum(lap * (clean - noisy)) / np.sum(lap ** 2))
            self.assertAlmostEqual(alpha, saved['noise_level_to_alpha'][str(sigma)], places=14)

    def test_checkpoint_reload_matches_saved_coordinates(self):
        checkpoint = torch.load(OUT / 'denoiser_seed4441.pt', map_location='cpu', weights_only=True)
        model = GeometryDenoiser().eval(); model.load_state_dict(checkpoint['state_dict'])
        records = [r for r in synthetic_dataset() if r['split'] == 'test'][:2]
        with torch.no_grad():
            actual = model(torch.tensor([r['noisy'] for r in records], dtype=torch.float32),
                           torch.tensor([r['noise_sigma'] for r in records], dtype=torch.float32)).numpy()
        with (OUT / 'predicted_coordinates.csv').open() as stream:
            rows = [r for r in csv.DictReader(stream) if r['model'] == 'EGNN_seed4441' and r['shape_id'] in {v['shape_id'] for v in records}]
        saved = np.array([[float(r[a]) for a in ['x', 'y', 'z']] for r in rows]).reshape(2, 8, 3)
        np.testing.assert_allclose(actual, saved, atol=2e-6, rtol=1e-6)

    def test_archived_source_defects_retained(self):
        audit = json.loads((OUT / 'source_and_reviewed_audit.json').read_text())
        self.assertFalse(audit['dtype']['source_double_succeeded'])
        self.assertGreater(audit['cutoff']['source_max_coordinate_change_from_edges_at_distance_8'], 1e-3)
        self.assertEqual(audit['cutoff']['reviewed_max_coordinate_change_from_edges_at_distance_8'], 0)
        self.assertGreater(abs(audit['references']['source']['charge_sum']), .1)
        self.assertLess(abs(audit['references']['reviewed']['charge_sum']), 1e-12)
        captured = json.loads((OUT / 'captured_source_audit.json').read_text())
        self.assertEqual(captured['model_parameter_count'], 55394)
        self.assertEqual(captured['counts']['energy_force_forwards'], 157)
        self.assertGreater(captured['origin_force_change_max'], 1e-4)
        self.assertIsNotNone(captured['double_failure'])

    def test_study_accounting_and_output_hashes(self):
        summary = json.loads((OUT / 'summary.json').read_text())
        self.assertFalse(summary['pilot'])
        self.assertEqual(summary['counts']['neural_training_runs'], 3)
        self.assertEqual(summary['counts']['training_epochs'], 180)
        self.assertEqual(summary['counts']['predicted_coordinate_rows'], 368 * 8 * 6)
        self.assertEqual(summary['counts']['DFT'], 0)
        for filename, expected in summary['output_sha256'].items():
            self.assertEqual(hashlib.sha256((OUT / filename).read_bytes()).hexdigest(), expected)


if __name__ == '__main__':
    unittest.main()
