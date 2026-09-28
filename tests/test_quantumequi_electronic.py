"""Numerical invariants; no quantum jobs or main-study fitting in this suite."""
import importlib.util
from pathlib import Path
import unittest

import numpy as np
import scipy.linalg as la

PATH = Path(__file__).resolve().parents[1] / "quantumequi/scripts/electronic_reviewed.py"
SPEC = importlib.util.spec_from_file_location("qe_electronic_reviewed_tests", PATH)
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)


class GeneralizedEigenTests(unittest.TestCase):
    def setUp(self):
        self.H = np.array([[-2., .2, -.1], [.2, 1., .3], [-.1, .3, 2.]])
        self.S = np.array([[1., .15, 0.], [.15, 1., .1], [0., .1, 1.]])

    def test_spd_matches_independent_lapack(self):
        answer = module.generalized_eigh(self.H, self.S)
        np.testing.assert_allclose(answer["energies"], la.eigh(self.H, self.S)[0], atol=1e-12)
        self.assertLess(answer["full_residual_max_abs"], 1e-12)

    def test_metric_orthonormality(self):
        C = module.generalized_eigh(self.H, self.S)["coefficients"]
        np.testing.assert_allclose(C.T@self.S@C, np.eye(3), atol=1e-12)

    def test_permuted_basis_preserves_eigenvalues(self):
        order = [2, 0, 1]
        changed = module.generalized_eigh(self.H[np.ix_(order, order)], self.S[np.ix_(order, order)])
        np.testing.assert_allclose(changed["energies"], module.generalized_eigh(self.H, self.S)["energies"], atol=1e-12)

    def test_nonsingular_congruence_preserves_spectrum(self):
        A = np.array([[1., .2, 0.], [0., 2., .1], [0., 0., .5]])
        transformed = module.generalized_eigh(A.T@self.H@A, A.T@self.S@A)
        np.testing.assert_allclose(transformed["energies"], module.generalized_eigh(self.H, self.S)["energies"], atol=1e-12)

    def test_rejects_asymmetric_hamiltonian(self):
        self.H[0, 1] = .4
        with self.assertRaises(ValueError):
            module.generalized_eigh(self.H, self.S)

    def test_rejects_nonfinite(self):
        self.S[0, 0] = np.nan
        with self.assertRaises(ValueError):
            module.generalized_eigh(self.H, self.S)

    def test_rejects_indefinite_even_in_canonical_mode(self):
        with self.assertRaisesRegex(ValueError, "Indefinite"):
            module.generalized_eigh(np.eye(2), np.diag([1., -.01]), mode="canonical")

    def test_strict_rejects_singular_metric(self):
        with self.assertRaisesRegex(ValueError, "singular"):
            module.generalized_eigh(np.eye(2), np.diag([1., 0.]))

    def test_consistent_psd_has_full_residual_zero(self):
        result = module.generalized_eigh(np.diag([2., 3., 0.]), np.diag([1., 2., 0.]), mode="canonical")
        self.assertEqual(result["rank"], 2)
        self.assertLess(result["full_residual_max_abs"], 1e-12)

    def test_inconsistent_psd_projected_solution_is_not_full_solution(self):
        H = np.array([[2., .4], [.4, 1.]])
        result = module.generalized_eigh(H, np.diag([1., 0.]), mode="canonical")
        self.assertLess(result["projected_residual_max_abs"], 1e-12)
        self.assertGreater(result["full_residual_max_abs"], .3)

    def test_occupancy_capacity_enforced(self):
        with self.assertRaises(ValueError):
            module.generalized_eigh(np.eye(2), np.diag([1., 0.]), mode="canonical", occupied=2)

    def test_odd_electrons_not_silently_floored(self):
        with self.assertRaises(ValueError):
            module.closed_shell_occupied(3, 3)
        self.assertEqual(module.closed_shell_occupied(4, 3), 2)

    def test_rank_threshold_control(self):
        S = np.diag([1., 1e-9])
        H = np.diag([2., 3e-9])
        self.assertEqual(module.generalized_eigh(H, S, mode="canonical", relative_cutoff=1e-8)["rank"], 1)
        self.assertEqual(module.generalized_eigh(H, S, relative_cutoff=1e-10)["rank"], 2)

    def test_all_null_overlap_rejected(self):
        with self.assertRaises(ValueError):
            module.generalized_eigh(np.eye(2), np.zeros((2, 2)), mode="canonical")

    def test_wrong_shapes_rejected(self):
        with self.assertRaises(ValueError):
            module.generalized_eigh(np.eye(2), np.eye(3))

    def test_source_same_atom_orbital_degeneracy(self):
        H, S, orbitals = module.source_matrices(["C"], [[0., 0., 0.]])
        self.assertEqual(len(orbitals), 4)
        np.testing.assert_array_equal(S, np.ones((4, 4)))
        self.assertEqual(np.linalg.matrix_rank(S), 1)
        with self.assertRaises(ValueError):
            module.closed_shell_occupied(4, 1)

    def test_clipping_solves_modified_metric(self):
        result = module.clipped_source_solution(np.eye(2), np.diag([1., 0.]), 1e-4)
        self.assertGreater(result["original_metric_orthonormality_max_abs"], .9)
        self.assertGreater(result["original_generalized_residual_max_abs"], 99)
        self.assertLess(result["modified_generalized_residual_max_abs"], 1e-8)


class SurrogateTests(unittest.TestCase):
    def test_frozen_partition_disjoint_and_ood_outside_training(self):
        plan = module.geometry_plan()
        self.assertEqual(len({p["point_id"] for p in plan}), 42)
        self.assertEqual(len({p["R_A"] for p in plan}), 42)
        self.assertEqual([sum(p["split"] == s for p in plan) for s in ("train", "validation", "test", "ood_stretch")], [17, 8, 8, 9])
        self.assertGreater(min(p["R_A"] for p in plan if p["split"] == "ood_stretch"), max(p["R_A"] for p in plan if p["split"] == "train"))

    def test_radial_feature_derivatives(self):
        r, centers, h = np.array([.7, 1.2]), np.array([.6, 1., 1.4]), 1e-6
        _, derivative = module.radial_features(r, centers, .2)
        fd = (module.radial_features(r+h, centers, .2)[0]-module.radial_features(r-h, centers, .2)[0])/(2*h)
        np.testing.assert_allclose(derivative, fd, atol=1e-9)

    def test_fitted_gradient_is_energy_derivative(self):
        r = np.linspace(.5, 1.5, 9)
        model = module.fit_radial(r, np.exp(-r), -np.exp(-r), .3, 1e-9, True)
        query, h = np.array([.71, 1.18]), 1e-5
        derivative = module.predict_radial(model, query)[1]
        fd = (module.predict_radial(model, query+h)[0]-module.predict_radial(model, query-h)[0])/(2*h)
        np.testing.assert_allclose(derivative, fd, atol=1e-8)


if __name__ == "__main__":
    unittest.main()
