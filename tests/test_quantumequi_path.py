"""Bounded invariants and saved exact-source replay checks; no model retraining."""
import csv
import hashlib
import json
from pathlib import Path
import unittest

import numpy as np

from quantumequi.scripts.path_reviewed import (
    AnalyticSurface, OUT, ROOT, improved_tangent, optimize_neb, projected_forces,
    source_tangent, source_update,
)


class TestQuantumEquiPath(unittest.TestCase):
    def test_both_analytic_forces_central_difference(self):
        for name in ("tilted_sine_2d", "periodic_curve_3d"):
            sf = AnalyticSurface(name)
            q = np.array([.27, -.31, .19])[:sf.dimension]
            _, force = sf.energy_force(q)
            for j in range(sf.dimension):
                d = np.eye(sf.dimension)[j]*1e-6
                fd = -(sf.energy_force(q+d)[0]-sf.energy_force(q-d)[0])/2e-6
                self.assertAlmostEqual(float(fd), force[j], places=8)

    def test_both_hessians_force_difference(self):
        for name in ("tilted_sine_2d", "periodic_curve_3d"):
            sf = AnalyticSurface(name)
            q = np.array([.27, -.31, .19])[:sf.dimension]
            h = np.empty((sf.dimension, sf.dimension))
            for j in range(sf.dimension):
                d = np.eye(sf.dimension)[j]*1e-5
                h[:, j] = -(sf.energy_force(q+d)[1]-sf.energy_force(q-d)[1])/2e-5
            np.testing.assert_allclose(h, sf.hessian(q), atol=2e-8)

    def test_stationary_reference_indices(self):
        for name in ("tilted_sine_2d", "periodic_curve_3d"):
            ref = AnalyticSurface(name).references()
            self.assertLess(ref["stationary_force_max"], 1e-12)
            for eigen in ref["hessian_eigenvalues"][:2]:
                self.assertGreater(min(eigen), 0)
            self.assertEqual(sum(e < 0 for e in ref["hessian_eigenvalues"][2]), 1)

    def test_asymmetric_forward_reverse_energetics(self):
        r = AnalyticSurface("tilted_sine_2d").references()
        self.assertGreater(abs(r["reaction_energy"]), .5)
        self.assertAlmostEqual(r["forward_barrier"]-r["reverse_barrier"], r["reaction_energy"], places=13)
        self.assertGreater(r["forward_barrier"], r["reverse_barrier"])

    def test_second_surface_is_three_dimensional(self):
        sf = AnalyticSurface("periodic_curve_3d")
        points = np.array([np.r_[x, sf.terms(x)[3]] for x in [0., .8, 2., 4., 2*np.pi]])
        self.assertEqual(np.linalg.matrix_rank(points-points.mean(0)), 3)
        self.assertAlmostEqual(sf.references()["forward_barrier"], 4.)

    def test_monotonic_tangent_selects_uphill_segment(self):
        a, b, c = np.array([0., 0.]), np.array([1., .1]), np.array([1.5, 1.])
        np.testing.assert_allclose(improved_tangent(a, b, c, 0, 1, 2), (c-b)/np.linalg.norm(c-b))
        np.testing.assert_allclose(improved_tangent(a, b, c, 2, 1, 0), (b-a)/np.linalg.norm(b-a))

    def test_energy_weighted_tangent_path_reversal(self):
        a, b, c = np.array([0., .2]), np.array([1., .3]), np.array([2., 1.])
        t = improved_tangent(a, b, c, .3, 2, 1)
        tr = improved_tangent(c, b, a, 1, 2, .3)
        np.testing.assert_allclose(t, -tr, atol=1e-15)

    def test_nonclimbing_force_perpendicular_component(self):
        band = np.array([[0., 0.], [1., 0.], [2., 0.]])
        force = np.array([[0., 0.], [3., 4.], [0., 0.]])
        projected, _ = projected_forces(band, [0., 1., 2.], force, climb=False)
        np.testing.assert_allclose(projected[1], [0., 4.])
        np.testing.assert_array_equal(projected[[0, -1]], 0)

    def test_climbing_force_independent_of_spring(self):
        band = np.array([[0., 0.], [1., .2], [2., 0.]])
        true = np.array([[0., 0.], [3., 4.], [0., 0.]])
        a, _ = projected_forces(band, [0., 1., 0.], true, spring=.1)
        b, _ = projected_forces(band, [0., 1., 0.], true, spring=100.)
        np.testing.assert_allclose(a, b, atol=0)

    def test_fixed_endpoints_and_fresh_energy(self):
        sf = AnalyticSurface("tilted_sine_2d"); ref = sf.references()
        result = optimize_neb(sf.energy_force, [ref["left"], ref["right"]], n_images=7, tolerance=1e-4)
        q = np.array(result["band"]); s = result["summary"]
        self.assertTrue(s["converged"])
        self.assertLessEqual(s["max_NEB_force"], 1e-4)
        self.assertLessEqual(s["climbing_true_force"], 1e-4)
        np.testing.assert_array_equal(q[[0, -1]], [ref["left"], ref["right"]])
        np.testing.assert_array_equal(sf.energy_force(q)[0], result["energies"])
        self.assertEqual(s["band_energy_force_point_evaluations"], 7*(s["iterations"]+1))

    def test_iteration_limit_is_not_convergence(self):
        sf = AnalyticSurface("tilted_sine_2d"); r = sf.references()
        result = optimize_neb(sf.energy_force, [r["left"], r["right"]], max_iterations=1)
        self.assertFalse(result["summary"]["converged"])
        self.assertEqual(result["summary"]["stop_reason"], "iteration_limit")
        self.assertEqual(len(result["history"]), 2)

    def test_nonstationary_endpoints_rejected(self):
        sf = AnalyticSurface("tilted_sine_2d")
        with self.assertRaises(ValueError):
            optimize_neb(sf.energy_force, [[-1., 0.], [1., 0.]])

    def test_degenerate_tangents_rejected(self):
        x = np.array([0., 0.])
        with self.assertRaises(ValueError):
            improved_tangent(x, x, x, 0., 0., 0.)
        with self.assertRaises(ValueError):
            source_tangent(x, x, np.ones(2))
        with self.assertRaises(ValueError):
            source_tangent(x, np.ones(2), x)

    def test_reviewed_reversal_identical_initial_band(self):
        sf = AnalyticSurface("tilted_sine_2d"); r = sf.references()
        band = np.linspace(r["left"], r["right"], 7)
        a = optimize_neb(sf.energy_force, band[[0, -1]], n_images=7, initial_band=band, tolerance=1e-4)
        b = optimize_neb(sf.energy_force, band[[-1, 0]], n_images=7, initial_band=band[::-1], tolerance=1e-4)
        np.testing.assert_allclose(a["band"], b["band"][::-1], atol=2e-12)

    def test_callback_supports_atom_coordinate_shape(self):
        sf = AnalyticSurface("tilted_sine_2d"); ref = sf.references()
        def callback(q):
            e, f = sf.energy_force(q[:, 0, :])
            return e, f[:, None, :]
        r = optimize_neb(callback, np.array([ref["left"], ref["right"]])[:, None, :],
                         n_images=7, tolerance=1e-4)
        self.assertEqual(np.array(r["band"]).shape, (7, 1, 2))
        self.assertTrue(r["summary"]["converged"])

    def test_sequential_update_has_direction_bias(self):
        q = np.array([[0., 0.], [.9, .15], [2.2, -.1], [3., 0.]])
        e = np.array([0., 1., .7, .2]); f = np.array([[0., 0.], [1., 2.], [-1., 1.], [0., 0.]])
        forward = source_update(q, e, f)
        reverse = source_update(q[::-1], e[::-1], f[::-1])[::-1]
        self.assertGreater(np.max(abs(forward-reverse)), 1e-5)
        a = source_update(q, e, f, sequential=False)
        b = source_update(q[::-1], e[::-1], f[::-1], sequential=False)[::-1]
        np.testing.assert_allclose(a, b, atol=1e-14)

    def test_seeded_solver_reproducibility(self):
        sf = AnalyticSurface("periodic_curve_3d"); r = sf.references()
        a = optimize_neb(sf.energy_force, [r["left"], r["right"]], max_iterations=3, seed=99)
        b = optimize_neb(sf.energy_force, [r["left"], r["right"]], max_iterations=3, seed=99)
        self.assertEqual(a, b)

    def test_invalid_parameters_and_band(self):
        sf = AnalyticSurface("periodic_curve_3d"); ref = sf.references(); ends = [ref["left"], ref["right"]]
        for kw in [{"n_images": 2}, {"tolerance": 0}, {"spring": -1}, {"perturbation": float("nan")},
                   {"max_iterations": 0}, {"initial_band": np.zeros((11, 3))}]:
            with self.assertRaises(ValueError):
                optimize_neb(sf.energy_force, ends, **kw)

    def test_archived_source_replay_provenance_and_candidate(self):
        a = json.loads((OUT/"source_audit.json").read_text())
        for p, h in a["source_sha256"].items():
            self.assertEqual(hashlib.sha256((ROOT/p).read_bytes()).hexdigest(), h)
        for p, h in a["output_sha256"].items():
            self.assertEqual(hashlib.sha256((OUT/p).read_bytes()).hexdigest(), h)
        self.assertTrue(a["cases"][0]["capture_match"])
        self.assertEqual(a["counts"]["source_energy_force_calls"], 420)
        self.assertEqual(a["counts"]["final_band_energy_force_calls"], 14)
        self.assertGreater(a["cases"][0]["source_candidate_true_max_atomic_force_nominal"], .01)
        self.assertGreater(a["reversal_mapped_max_coordinate_difference_A"], 1e-4)

    def test_archived_inplace_steps_reconstruct_exactly(self):
        for direction in ("forward", "reverse"):
            archive = np.load(OUT/f"source_{direction}_history.npz")
            bands = archive["preupdate_bands_A"]
            self.assertEqual(bands.shape, (30, 7, 8, 3))
            for step in range(30):
                actual = bands[step+1] if step < 29 else archive["final_band_A"]
                expected = source_update(bands[step], archive["preupdate_energies_nominal"][step],
                    archive["preupdate_forces_nominal"][step], climb=step >= 15)
                np.testing.assert_array_equal(expected, actual)


if __name__ == "__main__":
    unittest.main()
