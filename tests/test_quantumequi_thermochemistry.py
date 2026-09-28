"""Independent identities and failure gates for the thermochemistry API."""
import copy
import math
import unittest

import numpy as np

from quantumequi.scripts import thermochemistry_reviewed as tc


class ThermochemistryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cases = {}
        for name in ("linear_minimum", "linear_saddle", "nonlinear_minimum", "nonlinear_saddle"):
            x, m, fn, h = tc.control_case(name)
            e, g, ad = tc.torch_derivatives(fn, x)
            cls.cases[name] = (x, m, fn, h, tc.harmonic_analysis(ad, x, m, g))

    def test_nominal_frequency_factor_from_SI(self):
        # Independent Newton/metre force constant, 1 amu oscillator.
        expected = math.sqrt((4184 / 6.02214076e23) * 1e20 / 1.66053906892e-27) / (2 * math.pi * 299792458 * 100)
        self.assertAlmostEqual(tc.frequency_factor(), expected, places=10)
        self.assertAlmostEqual(tc.frequency_factor(), 108.591358535, places=8)
        self.assertGreater(tc.frequency_factor("hartree", "bohr"), 5000)

    def test_unit_scaling_equivalence(self):
        x, m, _, h, a = self.cases["nonlinear_minimum"]
        ev_per_kcal = tc.KCAL_MOL_J_PARTICLE / tc.EV_J
        other = tc.harmonic_analysis(h * ev_per_kcal, x, m, np.zeros_like(x), energy_unit="eV")
        np.testing.assert_allclose(a["projected_signed_frequencies_cm"], other["projected_signed_frequencies_cm"], rtol=1e-13)

    def test_mass_weighted_rigid_rank(self):
        for name, expected in [("linear_minimum", 5), ("nonlinear_minimum", 6)]:
            x, m, _, _, _ = self.cases[name]
            sub = tc.rigid_subspaces(x, m)
            self.assertEqual(sub["rigid_rank"], expected)
            q, b = sub["rigid_basis"], sub["internal_basis"]
            np.testing.assert_allclose(q.T @ b, 0, atol=1e-14)
            np.testing.assert_allclose(b.T @ b, np.eye(x.size - expected), atol=1e-14)

    def test_translation_rotation_frequency_invariance(self):
        x, m, _, h, a = self.cases["nonlinear_minimum"]
        angle = 0.713
        rotation = np.array([[math.cos(angle), -math.sin(angle), 0], [math.sin(angle), math.cos(angle), 0], [0, 0, 1.]])
        transform = np.kron(np.eye(len(x)), rotation)
        b = tc.harmonic_analysis(transform @ h @ transform.T, x @ rotation.T + [12, -3, 8], m, np.zeros_like(x))
        np.testing.assert_allclose(a["projected_signed_frequencies_cm"], b["projected_signed_frequencies_cm"], rtol=1e-13)
        np.testing.assert_allclose(a["inertia_amu_A2"], b["inertia_amu_A2"], rtol=1e-13)

    def test_isotope_frequency_scaling(self):
        x, m, _, h, a = self.cases["linear_minimum"]
        b = tc.harmonic_analysis(h, x, 4*m, np.zeros_like(x))
        np.testing.assert_allclose(b["projected_signed_frequencies_cm"], np.array(a["projected_signed_frequencies_cm"])/2)

    def test_diatomic_analytic_reduced_mass(self):
        _, m, _, _, a = self.cases["linear_minimum"]
        expected = tc.frequency_factor() * math.sqrt(70 * (1/m[0] + 1/m[1]))
        self.assertAlmostEqual(a["projected_signed_frequencies_cm"][0], expected, places=10)

    def test_minimum_saddle_classification(self):
        for name, (_, _, _, _, a) in self.cases.items():
            self.assertTrue(a["stationary"])
            self.assertEqual(a["negative_internal_modes"], int(name.endswith("saddle")))
            self.assertEqual(a["stationary_point_classification"], "first_order_saddle" if name.endswith("saddle") else "minimum")

    def test_blind_first_six_deletion_erases_saddle(self):
        a = self.cases["nonlinear_saddle"][-1]
        self.assertEqual(a["negative_internal_modes"], 1)
        self.assertLess(a["raw_signed_frequencies_cm"][0], -1)
        self.assertFalse(any(f < -0.01 for f in a["raw_signed_frequencies_cm"][6:]))

    def test_nonstationary_and_wrong_index_rejected(self):
        a = copy.deepcopy(self.cases["nonlinear_minimum"][-1])
        a["stationary"] = False
        with self.assertRaisesRegex(ValueError, "nonstationary"):
            tc.ideal_gas_rrho(a, [12, 16, 14])
        with self.assertRaisesRegex(ValueError, "classification"):
            tc.ideal_gas_rrho(self.cases["nonlinear_saddle"][-1], [12, 16, 14])

    def test_zero_internal_mode_rejected(self):
        x, m, _, h, _ = self.cases["linear_minimum"]
        a = tc.harmonic_analysis(h*0, x, m, np.zeros_like(x))
        with self.assertRaises(ValueError):
            tc.ideal_gas_rrho(a, m)

    def test_fd_sign_and_quadratic_exactness(self):
        x, _, _, h, _ = self.cases["nonlinear_minimum"]
        def force(y):
            return -(h @ (y-x).ravel()).reshape(x.shape)
        got, asymmetric = tc.finite_difference_hessian(force, x, 0.0001)
        np.testing.assert_allclose(got, h, atol=1e-10)
        self.assertLess(asymmetric, 1e-10)

    def test_partition_low_high_frequency_numerical_stability(self):
        terms = tc.vibrational_terms([1e-10, 1e7], 298.15)
        self.assertTrue(all(np.isfinite(v) for v in terms.values()))
        self.assertGreater(terms["vibrational_entropy_cal_mol_K"], 40)

    def test_vibration_classical_and_frozen_limits(self):
        classical = tc.vibrational_terms([1e-6], 300)
        frozen = tc.vibrational_terms([1e6], 300)
        self.assertAlmostEqual(classical["vibrational_Cv_cal_mol_K"], tc.R_CAL, places=8)
        self.assertAlmostEqual(classical["vibrational_thermal_enthalpy_kcal_mol"], tc.R_J*300/4184, places=8)
        self.assertEqual(frozen["vibrational_entropy_cal_mol_K"], 0)
        self.assertEqual(frozen["vibrational_thermal_enthalpy_kcal_mol"], 0)

    def test_pressure_gibbs_identity(self):
        _, m, _, _, a = self.cases["nonlinear_minimum"]
        t = 298.15
        low, high = [tc.ideal_gas_rrho(a, m, t, pressure_pa=p) for p in [1e5, 1e6]]
        self.assertAlmostEqual(high["G_correction_kcal_mol"] - low["G_correction_kcal_mol"], tc.R_J*t/4184*math.log(10), places=12)

    def test_symmetry_and_electronic_degeneracy(self):
        _, m, _, _, a = self.cases["nonlinear_minimum"]
        base = tc.ideal_gas_rrho(a, m)
        symmetry = tc.ideal_gas_rrho(a, m, symmetry_number=2)
        doublet = tc.ideal_gas_rrho(a, m, electronic_degeneracy=2)
        delta = tc.R_J*298.15/4184*math.log(2)
        self.assertAlmostEqual(symmetry["G_correction_kcal_mol"] - base["G_correction_kcal_mol"], delta, places=12)
        self.assertAlmostEqual(doublet["G_correction_kcal_mol"] - base["G_correction_kcal_mol"], -delta, places=12)

    def test_enthalpy_includes_translation_rotation_and_vibration(self):
        _, m, _, _, a = self.cases["nonlinear_minimum"]
        row = tc.ideal_gas_rrho(a, m)
        extra = row["thermal_enthalpy_excluding_ZPVE_kcal_mol"] - row["vibrational_thermal_enthalpy_kcal_mol"]
        self.assertAlmostEqual(extra, 4*tc.R_J*298.15/4184, places=12)
        self.assertAlmostEqual(row["G_correction_kcal_mol"], row["zpve_kcal_mol"] + row["thermal_enthalpy_excluding_ZPVE_kcal_mol"] - 298.15*row["total_entropy_cal_mol_K"]/1000, places=12)

    def test_first_order_saddle_partition_excludes_exactly_one(self):
        _, m, _, _, a = self.cases["nonlinear_saddle"]
        row = tc.ideal_gas_rrho(a, m, expected_stationary_point="first_order_saddle")
        self.assertEqual(row["unstable_modes_excluded"], 1)
        self.assertEqual(row["positive_modes_used"], 2)

    def test_atom_partition_has_no_rotation(self):
        x, m = [[0, 0, 0]], [40]
        a = tc.harmonic_analysis(np.zeros((3, 3)), x, m, [[0, 0, 0]])
        row = tc.ideal_gas_rrho(a, m)
        self.assertEqual(a["rigid_rank"], 3)
        self.assertEqual(row["rotation_entropy_cal_mol_K"], 0)
        self.assertEqual(row["positive_modes_used"], 0)

    def test_frequency_floor_is_explicit_and_counted(self):
        _, m, _, _, original = self.cases["nonlinear_minimum"]
        a = copy.deepcopy(original)
        a["projected_signed_frequencies_cm"] = [1, 500, 1000]
        base, floor = tc.ideal_gas_rrho(a, m), tc.ideal_gas_rrho(a, m, frequency_floor_cm=100)
        self.assertEqual(floor["modes_raised_by_floor"], 1)
        self.assertGreater(base["vibrational_entropy_cal_mol_K"], floor["vibrational_entropy_cal_mol_K"])
        self.assertGreater(floor["zpve_kcal_mol"], base["zpve_kcal_mol"])

    def test_invalid_inputs_are_rejected(self):
        with self.assertRaises(ValueError):
            tc.frequency_factor("unknown")
        with self.assertRaises(ValueError):
            tc.rigid_subspaces([[0, 0, 0]], [0])
        with self.assertRaises(ValueError):
            tc.rigid_subspaces([[0, 0, 0], [0, 0, 0]], [1, 1])
        with self.assertRaises(ValueError):
            tc.vibrational_terms([0], 300)
        with self.assertRaises(ValueError):
            tc.vibrational_terms([-20], 300)
        with self.assertRaises(ValueError):
            tc.ideal_gas_rrho(self.cases["linear_minimum"][-1], [12, 16], pressure_pa=0)


if __name__ == "__main__":
    unittest.main()
