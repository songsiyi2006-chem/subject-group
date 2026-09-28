"""Fast invariant and bookkeeping tests; no neural training is run here."""
import csv
import json
from pathlib import Path
import unittest

import numpy as np
import torch

from electrograph.scripts.graph_learning import (
    DATA, OUT, EDGE_DIM, NODE_DIM, Graph, HydrationMPNN, assign_splits,
    batch_graphs, campaign, choose_ucb, descriptors, gp_posterior,
    group_key, load_records, matern52, molecular_graph, structure_key,
)


class TestElectrographLearning(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        torch.set_num_threads(1)

    def test_graph_shapes_and_reciprocal_edges(self):
        graph = molecular_graph('CC(=O)Oc1ccccc1')
        self.assertEqual(graph.x.shape, (10, NODE_DIM))
        self.assertEqual(graph.edge_attr.shape[1], EDGE_DIM)
        pairs = list(zip(*graph.edge_index.tolist()))
        for i, pair in enumerate(pairs):
            j = pairs.index(pair[::-1])
            torch.testing.assert_close(graph.edge_attr[i], graph.edge_attr[j])

    def test_isolated_atom_and_unknown_element(self):
        graph = molecular_graph('[Xe]')
        self.assertEqual(graph.edge_index.shape, (2, 0))
        self.assertEqual(graph.edge_attr.shape, (0, EDGE_DIM))
        self.assertEqual(float(graph.x[0, 12]), 1.0)
        model = HydrationMPNN()
        pred, latent = model(batch_graphs([graph]))
        self.assertTrue(torch.isfinite(pred).all())
        self.assertEqual(latent.shape, (1, 48))
        pred.sum().backward()
        self.assertGreater(float(model.embedding[0].weight.grad.abs().sum()), 0)

    def test_invalid_and_empty_input_rejected(self):
        for smiles in ('', 'C1CC', 'not-a-molecule'):
            with self.assertRaises(ValueError):
                molecular_graph(smiles)
        with self.assertRaises(ValueError):
            batch_graphs([])

    def test_atom_permutation_invariance(self):
        torch.manual_seed(812)
        graph = molecular_graph('CC(=O)Oc1ccccc1C(=O)O')
        permutation = torch.randperm(len(graph.x))
        inverse = torch.argsort(permutation)
        changed = Graph(graph.x[permutation], inverse[graph.edge_index], graph.edge_attr.clone())
        model = HydrationMPNN().eval()
        a, za = model(batch_graphs([graph]))
        b, zb = model(batch_graphs([changed]))
        torch.testing.assert_close(a, b, rtol=1e-5, atol=2e-6)
        torch.testing.assert_close(za, zb, rtol=1e-5, atol=2e-6)

    def test_batch_separation_and_order(self):
        torch.manual_seed(81)
        graphs = [molecular_graph(s) for s in ('CCO', 'C', 'c1ccncc1', '[Na+]')]
        model = HydrationMPNN().eval()
        joint, _ = model(batch_graphs(graphs))
        separate = torch.cat([model(batch_graphs([g]))[0] for g in graphs])
        reversed_pred, _ = model(batch_graphs(graphs[::-1]))
        torch.testing.assert_close(joint, separate, rtol=1e-5, atol=2e-6)
        torch.testing.assert_close(joint, reversed_pred.flip(0), rtol=1e-5, atol=2e-6)

    def test_equivalent_chiral_smiles_have_same_prediction(self):
        left, right = 'N[C@@H](C)C(=O)O', 'C[C@H](N)C(=O)O'
        self.assertEqual(structure_key(left), structure_key(right))
        torch.manual_seed(812)
        model = HydrationMPNN().eval()
        a, _ = model(batch_graphs([molecular_graph(left)]))
        b, _ = model(batch_graphs([molecular_graph(right)]))
        torch.testing.assert_close(a, b, rtol=1e-5, atol=2e-6)

    def test_edge_network_gets_finite_gradient(self):
        torch.manual_seed(6)
        model = HydrationMPNN()
        pred, _ = model(batch_graphs([molecular_graph('CCN'), molecular_graph('C=C')]))
        ((pred - torch.tensor([1., -2.])) ** 2).mean().backward()
        for layer in model.layers:
            grad = layer.edge_net[0].weight.grad
            self.assertTrue(torch.isfinite(grad).all())
            self.assertGreater(float(grad.abs().sum()), 0)

    def test_structure_and_group_identity(self):
        self.assertEqual(structure_key('OCC'), structure_key('CCO'))
        self.assertEqual(group_key('Cc1ccccc1'), group_key('Oc1ccccc1'))
        self.assertNotEqual(group_key('CCC'), group_key('CCCC'))
        self.assertTrue(group_key('CCC').startswith('acyclic:'))

    def test_saved_splits_have_no_group_or_structure_leakage(self):
        records = load_records()
        self.assertEqual(len(records), 256)
        self.assertEqual(len({r['smiles'] for r in records}), len(records))
        groups = {name: {r['group_key'] for r in records if r['split'] == name}
                  for name in ('train', 'validation', 'test')}
        self.assertTrue(groups['train'].isdisjoint(groups['validation']))
        self.assertTrue(groups['train'].isdisjoint(groups['test']))
        self.assertTrue(groups['validation'].isdisjoint(groups['test']))
        replay = assign_splits([{**r, 'expt_kcal_mol': 1e6 - r['expt_kcal_mol']} for r in records])
        self.assertEqual([r['split'] for r in replay], [r['split'] for r in records])

    def test_finite_descriptors_and_recorded_target(self):
        records = load_records()
        X = np.stack([descriptors(r['smiles']) for r in records])
        self.assertEqual(X.shape, (256, 10))
        self.assertTrue(np.isfinite(X).all())
        p = json.loads((DATA / 'provenance.json').read_text(encoding='utf-8'))
        self.assertEqual(p['target_unit'], 'kcal/mol')
        self.assertIn('hydration', p['dataset'])

    def test_gp_matches_independent_dense_conditioning(self):
        x = np.array([[0., 1.], [1., 0.], [2., 3.], [-1., .2]])
        ids, y = [0, 2], np.array([-2., 4.])
        mean, sd = gp_posterior(x, ids, y)
        K = matern52(x[ids], x[ids]) + 1e-5 * np.eye(2)
        k = matern52(x, x[ids])
        reference_mean = y.mean() + y.std() * k @ np.linalg.solve(K, (y - y.mean()) / y.std())
        reference_sd = y.std() * np.sqrt(np.maximum(1 - np.diag(k @ np.linalg.solve(K, k.T)), 0))
        np.testing.assert_allclose(mean, reference_mean, atol=1e-10)
        np.testing.assert_allclose(sd, reference_sd, atol=1e-10)
        self.assertLess(float(sd[ids].max()), .02)

    def test_gp_positive_affine_output_equivariance(self):
        x = np.arange(12.).reshape(6, 2)
        a, sa = gp_posterior(x, [1, 4], [1., 3.])
        b, sb = gp_posterior(x, [1, 4], [7., 15.])
        np.testing.assert_allclose(b, 4 * a + 3, atol=1e-10)
        np.testing.assert_allclose(sb, 4 * sa, atol=1e-10)

    def test_gp_rejects_duplicates_nonfinite_and_exhausted_pool(self):
        x = np.arange(6.).reshape(3, 2)
        for ids, values in [([0, 0], [1., 2.]), ([0], [np.nan]), ([3], [1.]), ([], [])]:
            with self.assertRaises(ValueError):
                gp_posterior(x, ids, values)
        with self.assertRaises(ValueError):
            choose_ucb(x, [0, 1, 2], [1., 2., 3.])

    def test_campaign_budget_shared_initial_and_monotonic_best(self):
        x = np.linspace(-1, 1, 40).reshape(20, 2)
        y = -(x[:, 0] - .2) ** 2
        initial = []
        for method in ('latent_gp_ucb', 'descriptor_gp_ucb', 'random'):
            rows = campaign(x, y, seed=7, budget=12, method=method)
            ids = [r['candidate_index'] for r in rows]
            self.assertEqual(len(set(ids)), 12)
            self.assertEqual([r['observed_count_before_selection'] for r in rows], list(range(12)))
            self.assertTrue(np.all(np.diff([r['best_observed_objective'] for r in rows]) >= 0))
            initial.append(ids[:4])
        self.assertEqual(initial[0], initial[1])
        self.assertEqual(initial[0], initial[2])
        with self.assertRaises(ValueError):
            campaign(x, y, seed=7, budget=21)

    def test_first_adaptive_choice_cannot_read_unobserved_outcomes(self):
        x = np.linspace(-1, 1, 40).reshape(20, 2)
        y = np.sin(x[:, 0])
        initial = np.random.default_rng(21).choice(20, size=4, replace=False)
        changed = y.copy()
        changed[np.setdiff1d(np.arange(20), initial)] += 10000
        first = campaign(x, y, seed=21, budget=5)
        second = campaign(x, changed, seed=21, budget=5)
        self.assertEqual(first[4]['candidate_index'], second[4]['candidate_index'])
        self.assertEqual(first[4]['posterior_sd_selected'], second[4]['posterior_sd_selected'])

    def test_saved_model_reload_and_output_accounting(self):
        summary = json.loads((OUT / 'summary.json').read_text(encoding='utf-8'))
        self.assertFalse(summary['pilot'])
        self.assertEqual(summary['counts']['cached_oracle_calls'], 8 * 3 * 16)
        self.assertEqual(summary['counts']['new_experiments'], 0)
        checkpoint = torch.load(OUT / 'mpnn_seed20260928.pt', map_location='cpu', weights_only=True)
        model = HydrationMPNN().eval()
        model.load_state_dict(checkpoint['state_dict'])
        records = load_records()[:8]
        with torch.no_grad():
            pred, _ = model(batch_graphs([molecular_graph(r['smiles']) for r in records]))
        pred = pred.numpy() * checkpoint['target_scale'] + checkpoint['target_mean']
        with (OUT / 'predictions.csv').open(encoding='utf-8') as handle:
            saved = {r['molecule_id']: float(r['predicted_hydration_kcal_mol']) for r in csv.DictReader(handle)
                     if r['model'] == 'mpnn_seed20260928'}
        np.testing.assert_allclose(pred, [saved[r['molecule_id']] for r in records], atol=5e-5, rtol=1e-5)


if __name__ == '__main__':
    unittest.main()
