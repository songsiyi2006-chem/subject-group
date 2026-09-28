"""Independent conservation, analytic-limit, and input-contract tests."""
import math
import unittest
from unittest.mock import patch

import numpy as np

from electratwin.scripts.transport_reviewed import F, effective_rates, solve_transport


class TestElectraTwinTransport(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = solve_transport(nx=60, ny=24)

    def test_boundary_integrals_close_material_and_charge(self):
        r = self.result
        p, s, f = r["inputs"], r["summary"], r["field"]
        molar_in = p["flow_rate_uL_min"] * 1e-9 / 60 * p["inlet_A_mM"]
        velocity = np.array(f["velocity_m_s"])
        outlet = np.array(f["A_mol_m3"])[-1]
        molar_out = np.sum(velocity * outlet) * p["height_m"] / p["ny"] * p["width_m"]
        wall_rate = np.sum(f["wall_flux_mol_m2_s"]) * p["length_m"] / p["nx"] * p["width_m"]
        self.assertLess(abs(molar_in - molar_out - wall_rate) / molar_in, 1e-10)
        self.assertAlmostEqual(s["current_A"], 2 * F * (molar_in - molar_out), places=11)
        self.assertAlmostEqual(s["outlet_product_gain_mol_s"], wall_rate, places=16)

    def test_wall_kinetics_and_diffusion_use_same_flux(self):
        r = self.result
        p, s, f = r["inputs"], r["summary"], r["field"]
        ca, cp = np.array(f["wall_A_mol_m3"]), np.array(f["wall_P_mol_m3"])
        reaction = s["forward_rate_m_s"] * ca - s["reverse_rate_m_s"] * cp
        gradient_flux = p["diffusivity_m2_s"] * (np.array(f["A_mol_m3"])[:, 0] - ca) / (p["height_m"] / p["ny"] / 2)
        np.testing.assert_allclose(reaction, f["wall_flux_mol_m2_s"], rtol=1e-12, atol=1e-16)
        np.testing.assert_allclose(reaction, gradient_flux, rtol=1e-12, atol=1e-16)

    def test_zero_reaction_preserves_uniform_feed(self):
        r = solve_transport(nx=20, ny=12, inlet_A_mM=37, inlet_P_mM=13,
                            forward_rate_m_s=0, reverse_rate_m_s=0)
        np.testing.assert_allclose(r["field"]["A_mol_m3"], 37, atol=1e-10)
        np.testing.assert_allclose(r["field"]["P_mol_m3"], 13, atol=1e-10)
        self.assertEqual(r["summary"]["current_A"], 0)

    def test_zero_diffusivity_blocks_transport_to_wall(self):
        r = solve_transport(nx=20, ny=12, diffusivity_m2_s=0)
        np.testing.assert_allclose(r["field"]["A_mol_m3"], 50, atol=1e-10)
        self.assertEqual(r["summary"]["current_A"], 0)
        self.assertLess(abs(r["summary"]["conversion_pct"]), 1e-10)

    def test_equilibrium_equal_forward_reverse_rates(self):
        r = solve_transport(nx=20, ny=12, inlet_A_mM=25, inlet_P_mM=25, overpotential_V=0)
        self.assertLess(abs(r["summary"]["current_A"]), 1e-12)
        np.testing.assert_allclose(r["field"]["A_mol_m3"], 25, atol=1e-10)

    def test_reverse_reaction_is_negative_not_clipped(self):
        forward = self.result
        reverse = solve_transport(nx=60, ny=24, inlet_A_mM=0, inlet_P_mM=50, overpotential_V=-0.48)
        self.assertLess(reverse["summary"]["current_A"], 0)
        self.assertIsNone(reverse["summary"]["conversion_pct"])
        self.assertAlmostEqual(reverse["summary"]["current_A"], -forward["summary"]["current_A"], places=11)
        np.testing.assert_allclose(reverse["field"]["A_mol_m3"], forward["field"]["P_mol_m3"], atol=1e-10)

    def test_electron_stoichiometry_changes_current_not_fixed_rate_conversion(self):
        args = dict(nx=25, ny=12, forward_rate_m_s=1e-5, reverse_rate_m_s=1e-8, include_field=False)
        one = solve_transport(n_electrons=1, **args)["summary"]
        two = solve_transport(n_electrons=2, **args)["summary"]
        self.assertAlmostEqual(one["conversion_pct"], two["conversion_pct"], places=12)
        self.assertAlmostEqual(two["current_A"], 2 * one["current_A"], places=13)

    def test_channel_numbering_up_scale_preserves_conversion(self):
        args = dict(nx=25, ny=12, include_field=False)
        one = solve_transport(width_m=.012, flow_rate_uL_min=450, **args)["summary"]
        two = solve_transport(width_m=.024, flow_rate_uL_min=900, **args)["summary"]
        self.assertAlmostEqual(one["conversion_pct"], two["conversion_pct"], places=10)
        self.assertAlmostEqual(two["current_A"], 2 * one["current_A"], places=11)

    def test_exact_flow_integral_on_coarse_grid(self):
        for ny in [2, 3, 7, 32]:
            s = solve_transport(nx=5, ny=ny, include_field=False)["summary"]
            self.assertAlmostEqual(s["discrete_flow_rate_m3_s"] / (450e-9 / 60), 1, places=14)

    def test_analytic_mixed_transverse_plug_flow_limit(self):
        k, q, area = 1e-6, 450e-9 / 60, .06 * .012
        exact = 100 * (1 - math.exp(-k * area / q))
        errors = []
        for nx in [40, 80, 160]:
            s = solve_transport(nx=nx, ny=8, diffusivity_m2_s=1e-4, axial_diffusivity_m2_s=0,
                                forward_rate_m_s=k, reverse_rate_m_s=0, velocity_profile="plug",
                                include_field=False)["summary"]
            errors.append(abs(s["conversion_pct"] - exact))
        self.assertTrue(errors[2] < errors[1] < errors[0])
        self.assertLess(errors[-1], .003)
        self.assertGreater(errors[0] / errors[1], 1.8)

    def test_nonnegative_species_and_one_to_one_molecular_balance(self):
        a = np.array(self.result["field"]["A_mol_m3"])
        p = np.array(self.result["field"]["P_mol_m3"])
        self.assertGreater(a.min(), 0)
        self.assertGreater(p.min(), 0)
        np.testing.assert_allclose(a + p, 50, atol=1e-12)

    def test_axial_diffusion_option_has_effect(self):
        with_axial = solve_transport(nx=30, ny=12, diffusivity_m2_s=1e-6, include_field=False)
        without_axial = solve_transport(nx=30, ny=12, diffusivity_m2_s=1e-6,
                                        axial_diffusivity_m2_s=0, include_field=False)
        self.assertGreater(abs(with_axial["summary"]["conversion_pct"] - without_axial["summary"]["conversion_pct"]), .01)

    def test_corrupt_solver_return_does_not_force_converged(self):
        with patch("electratwin.scripts.transport_reviewed.spsolve", return_value=np.full(6 * 4, -2.0)):
            r = solve_transport(nx=6, ny=4)
        self.assertFalse(r["summary"]["converged"])
        self.assertFalse(r["summary"]["nonnegative"])
        self.assertEqual(r["summary"]["minimum_bulk_A_mM"], -2)
        self.assertGreater(r["summary"]["conversion_pct"], 100)

    def test_bad_inputs_are_rejected(self):
        for bad in [dict(flow_rate_uL_min=0), dict(diffusivity_m2_s=-1), dict(nx=1), dict(nx=1000, ny=100),
                    dict(inlet_A_mM=0, inlet_P_mM=0), dict(n_electrons=0), dict(alpha=1),
                    dict(overpotential_V=float("nan")), dict(forward_rate_m_s=1), dict(velocity_profile="unknown")]:
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                solve_transport(**bad)

    def test_effective_rates_have_explicit_concentration_scale(self):
        a, b = effective_rates(.2, reference_concentration_mM=50)
        doubled_a, doubled_b = effective_rates(.2, reference_concentration_mM=100)
        self.assertAlmostEqual(a, 2 * doubled_a, places=15)
        self.assertAlmostEqual(b, 2 * doubled_b, places=15)


if __name__ == "__main__":
    unittest.main()
