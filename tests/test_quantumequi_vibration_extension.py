"""Bounded spectral and units tests; no electronic-structure calculation."""
import csv
import json
import math
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import numpy as np

from quantumequi.scripts.vibration_extension import (
    AMU_KG, HARTREE_J, HBAR, NUCLEAR_MASS_AMU, bound_partition,
    exact_morse_levels, fci_parameterization, find_equilibrium_records,
    harmonic_partition, harmonic_spacing, kinetic_coefficient,
    morse_potential, reduced_mass, solve_radial_fd,
)


class TestQuantumEquiVibration(unittest.TestCase):
    def setUp(self):
        self.mu = NUCLEAR_MASS_AMU["H2"]/2
        self.potential = lambda r: morse_potential(r, .2, 1.8, 1.)

    def test_reduced_mass_and_SI_kinetic_units(self):
        self.assertEqual(reduced_mass(2., 2.), 1.)
        self.assertAlmostEqual(reduced_mass(1., 2.), 2/3)
        c = HBAR**2/(2*self.mu*AMU_KG)/HARTREE_J/1e-20
        self.assertAlmostEqual(kinetic_coefficient(self.mu)/c, 1., places=14)
        self.assertAlmostEqual(kinetic_coefficient(2*self.mu)/c, .5, places=14)

    def test_Morse_minimum_and_curvature(self):
        self.assertEqual(float(self.potential(1.)), 0.)
        h = 1e-5
        curvature = (self.potential(1+h)+self.potential(1-h))/h**2
        self.assertAlmostEqual(float(curvature), 2*.2*1.8**2, places=8)
        self.assertAlmostEqual(float(self.potential(50.)), .2, places=14)

    def test_exact_bound_spectrum_cutoff_and_spacing(self):
        q = exact_morse_levels(.2, 1.8, self.mu)
        energy = q["energies_Hartree"]; n = q["quantum_numbers"]; lam = q["lambda"]
        self.assertTrue(np.all(energy < .2))
        self.assertTrue(np.all(np.diff(energy) > 0))
        self.assertTrue(np.all(n+.5 < lam))
        self.assertGreaterEqual(n[-1]+1.5, lam)
        self.assertLess(energy[0], .5*q["harmonic_spacing_Hartree"])
        self.assertLess(energy[1]-energy[0], q["harmonic_spacing_Hartree"])

    def test_weak_well_has_no_bound_full_line_state(self):
        exact = exact_morse_levels(1e-7, 3., self.mu)
        self.assertEqual(len(exact["energies_Hartree"]), 0)
        numerical = solve_radial_fd(lambda r: morse_potential(r, 1e-7, 3., 1.), self.mu, 0., 8., 200, 1e-7)
        self.assertEqual(numerical["bound_states_returned"], 0)

    def test_free_particle_box_matches_discrete_formula(self):
        n = 200
        result = solve_radial_fd(lambda r: np.zeros_like(r), self.mu, 0., 1., n, .05)
        j = np.arange(1, len(result["energies_Hartree"])+1)
        discrete = 4*kinetic_coefficient(self.mu)*n*n*np.sin(j*np.pi/(2*n))**2
        np.testing.assert_allclose(result["energies_Hartree"], discrete, atol=2e-12)

    def test_second_order_Morse_grid_convergence(self):
        exact = exact_morse_levels(.2, 1.8, self.mu)["energies_Hartree"][0]
        a = solve_radial_fd(self.potential, self.mu, 0., 8., 400, .2)["energies_Hartree"][0]
        b = solve_radial_fd(self.potential, self.mu, 0., 8., 800, .2)["energies_Hartree"][0]
        self.assertGreater(abs(a-exact)/abs(b-exact), 3.95)
        self.assertLess(abs(a-exact)/abs(b-exact), 4.05)

    def test_wavefunction_quadrature_and_eigen_residual(self):
        result = solve_radial_fd(self.potential, self.mu, 0., 8., 400, .2)
        self.assertLess(result["orthonormality_max"], 1e-12)
        self.assertLess(result["eigen_residual_max_Hartree"], 1e-10)
        self.assertEqual(result["interior_nodes"], 399)
        self.assertTrue(np.all(result["r_A"] > 0))
        self.assertTrue(np.all(result["r_A"] < 8))

    def test_hard_wall_is_not_full_line_Morse(self):
        exact = exact_morse_levels(.2, 1.8, self.mu)["energies_Hartree"][0]
        q = solve_radial_fd(lambda r: morse_potential(r, .2, 1.8, .12), self.mu, 0., 8., 800, .2)
        self.assertGreater(q["energies_Hartree"][0]-exact, .001)

    def test_short_domain_removes_near_dissociation_states(self):
        short = solve_radial_fd(self.potential, self.mu, 0., 2., 300, .2)
        long = solve_radial_fd(self.potential, self.mu, 0., 8., 1200, .2)
        self.assertLess(short["bound_states_returned"], long["bound_states_returned"])

    def test_isotope_harmonic_and_anharmonic_trends(self):
        h = exact_morse_levels(.2, 1.8, NUCLEAR_MASS_AMU["H2"]/2)
        d = exact_morse_levels(.2, 1.8, NUCLEAR_MASS_AMU["D2"]/2)
        expected = math.sqrt(NUCLEAR_MASS_AMU["D2"]/NUCLEAR_MASS_AMU["H2"])
        self.assertAlmostEqual(h["harmonic_spacing_Hartree"]/d["harmonic_spacing_Hartree"], expected, places=13)
        self.assertGreater(len(d["energies_Hartree"]), len(h["energies_Hartree"]))
        self.assertLess(d["energies_Hartree"][0], h["energies_Hartree"][0])
        gap_ratio = np.diff(h["energies_Hartree"][:2])[0]/np.diff(d["energies_Hartree"][:2])[0]
        self.assertGreater(gap_ratio, 1.)
        self.assertLess(gap_ratio, expected)

    def test_bound_partition_zero_temperature_limit(self):
        q = bound_partition([.01, .03, .04], 1.)
        self.assertAlmostEqual(q["U_vib_Hartree"], .01, places=14)
        self.assertAlmostEqual(q["F_vib_bound_only_Hartree"], .01, places=14)
        self.assertAlmostEqual(q["S_vib_J_mol_K"], 0., places=10)
        self.assertTrue(q["continuum_excluded"])

    def test_finite_state_partition_high_temperature_limit(self):
        q = bound_partition([0., .01], 1e12)
        self.assertAlmostEqual(q["log_q_relative_to_v0"], math.log(2), places=8)
        self.assertAlmostEqual(q["U_vib_Hartree"], .005, places=10)

    def test_harmonic_partition_matches_many_exact_levels(self):
        spacing = .008
        levels = spacing*(np.arange(100)+.5)
        a, b = bound_partition(levels, 2000.), harmonic_partition(spacing, 2000.)
        self.assertAlmostEqual(a["F_vib_bound_only_Hartree"], b["F_vib_harmonic_Hartree"], places=13)
        self.assertAlmostEqual(a["U_vib_Hartree"], b["U_vib_harmonic_Hartree"], places=13)

    def test_bad_units_domains_and_temperatures_rejected(self):
        for value in (0., -1., float("inf")):
            with self.assertRaises(ValueError):
                kinetic_coefficient(value)
            with self.assertRaises(ValueError):
                bound_partition([.01], value)
        for left, right, n in [(-1., 8., 400), (1., 0., 400), (0., 8., 2)]:
            with self.assertRaises(ValueError):
                solve_radial_fd(self.potential, self.mu, left, right, n, .2)
        for energies in ([], [.2, .1], [float("nan")], [-.1]):
            with self.assertRaises(ValueError):
                bound_partition(energies, 300.)

    def test_parameterization_recovers_fixture_not_curve_fit(self):
        # Synthetic fixture tests the schema and algebra only, never saved as FCI evidence.
        record = {"basis": "TEST_FIXTURE", "method": "FCI", "R_e_A": 1., "E_min_Hartree": -1.,
                  "D_e_Hartree": .2, "curvature_Hartree_A2": 2*.2*1.8**2, "E_infinity_Hartree": -.8}
        with TemporaryDirectory() as temp:
            root = Path(temp); a = root/"equilibrium.json"; b = root/"curve.csv"
            a.write_text(json.dumps([record]))
            rows = [{"basis": "TEST_FIXTURE", "method": "FCI", "R_A": r,
                     "energy_total_Hartree": -1+float(self.potential(r)), "point_id": str(i)}
                    for i, r in enumerate([.7, 1., 1.5, 2., 3.])]
            with b.open("w", newline="") as f:
                w = csv.DictWriter(f, list(rows[0])); w.writeheader(); w.writerows(rows)
            model, residual = fci_parameterization(a, b)
        self.assertAlmostEqual(model[0]["a_A_inverse"], 1.8, places=14)
        self.assertFalse(model[0]["is_nonlinear_fit"])
        self.assertLess(max(abs(r["Morse_minus_FCI_Hartree"]) for r in residual), 1e-14)

    def test_missing_or_duplicate_equilibrium_schema_rejected(self):
        with self.assertRaises(ValueError):
            find_equilibrium_records({"incomplete": []})
        row = {"basis": "x", "R_e_A": 1., "E_min_Hartree": -1., "D_e_Hartree": .2, "curvature_Hartree_A2": 1.}
        with self.assertRaises(ValueError):
            find_equilibrium_records([row, row])

    def test_failed_FCI_equilibrium_rejected_before_curve_use(self):
        record = {"basis": "TEST_FIXTURE", "method": "FCI", "R_e_A": 1., "E_min_Hartree": -1.,
                  "D_e_Hartree": .2, "curvature_Hartree_A2": 1., "E_infinity_Hartree": -.8,
                  "optimization_success": False}
        with TemporaryDirectory() as temp:
            a, b = Path(temp)/"equilibrium.json", Path(temp)/"curve.csv"
            a.write_text(json.dumps([record])); b.write_text("basis,method,R_A,energy_total_Hartree,point_id\n")
            with self.assertRaisesRegex(ValueError, "optimization"):
                fci_parameterization(a, b)

    def test_failed_FCI_curve_row_rejected(self):
        record = {"basis": "TEST_FIXTURE", "method": "FCI", "R_e_A": 1., "E_min_Hartree": -1.,
                  "D_e_Hartree": .2, "curvature_Hartree_A2": 1., "E_infinity_Hartree": -.8}
        with TemporaryDirectory() as temp:
            a, b = Path(temp)/"equilibrium.json", Path(temp)/"curve.csv"
            a.write_text(json.dumps([record]))
            rows = [{"basis": "TEST_FIXTURE", "method": "FCI", "R_A": r,
                     "energy_total_Hartree": -1., "point_id": str(i),
                     "status": "failed" if i == 2 else "converged"}
                    for i, r in enumerate([.7, 1., 1.5, 2., 3.])]
            with b.open("w", newline="") as f:
                w = csv.DictWriter(f, list(rows[0])); w.writeheader(); w.writerows(rows)
            with self.assertRaisesRegex(ValueError, "Failed FCI"):
                fci_parameterization(a, b)


if __name__ == "__main__":
    unittest.main()
