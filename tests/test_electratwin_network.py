"""Network verification against analytic limits and independent balance laws."""
import math
import unittest
from unittest.mock import patch

import numpy as np

from electratwin.scripts.reaction_network import F, analytic_plug_outlet, network_rates, solve_network
from electratwin.scripts.transport_reviewed import solve_transport


class TestElectraTwinNetwork(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.baseline = solve_network(nx=30, ny=12)

    def test_all_four_species_are_present_and_nonnegative(self):
        fields = self.baseline["field"]["concentration_mM"]
        self.assertEqual(set(fields), {"A", "P", "B", "D"})
        for name, values in fields.items():
            self.assertEqual(np.shape(values), (30, 12), name)
            self.assertGreaterEqual(np.min(values), -1e-10, name)

    def test_total_molecular_count_conserved_spatially(self):
        total = sum(np.array(v) for v in self.baseline["field"]["concentration_mM"].values())
        np.testing.assert_allclose(total, 50, atol=1e-10)

    def test_species_inlet_outlet_and_wall_flux_independently_close(self):
        r = self.baseline
        p, f = r["inputs"], r["field"]
        q = p["flow_rate_uL_min"] * 1e-9 / 60
        dy, dx = p["height_m"] / p["ny"], p["length_m"] / p["nx"]
        for index, species in enumerate(("A", "P", "B", "D")):
            inward = q * p["inlet_mM"][index]
            outward = np.dot(f["velocity_m_s"], np.array(f["concentration_mM"][species])[-1]) * dy * p["width_m"]
            reacted = np.sum(f["wall_outward_flux_mol_m2_s"][species]) * dx * p["width_m"]
            self.assertLess(abs(inward - outward - reacted) / (q * 50), 1e-11)

    def test_robin_wall_satisfies_kinetics_and_half_cell_gradient(self):
        r = self.baseline
        p, f = r["inputs"], r["field"]
        w = {s: np.array(v) for s, v in f["wall_concentration_mM"].items()}
        k = p["rate_constants_m_s"]
        expected = {"A": (k["AP"] + k["AB"]) * w["A"], "P": k["PD"] * w["P"] - k["AP"] * w["A"],
                    "B": -k["AB"] * w["A"], "D": -k["PD"] * w["P"]}
        for s in ("A", "P", "B", "D"):
            diffusive = p["diffusivity_m2_s"] * (np.array(f["concentration_mM"][s])[:, 0] - w[s]) / (p["height_m"] / p["ny"] / 2)
            np.testing.assert_allclose(f["wall_outward_flux_mol_m2_s"][s], expected[s], atol=1e-15, rtol=1e-12)
            np.testing.assert_allclose(diffusive, expected[s], atol=1e-15, rtol=1e-12)

    def test_current_accounting_from_products_counts_overoxidation_twice(self):
        s = self.baseline["summary"]
        gain = s["net_outlet_gain_mol_s"]
        from_bulk = 2 * F * (gain["P"] + gain["B"] + 2 * gain["D"])
        self.assertAlmostEqual(s["current_A"], from_bulk, places=11)
        self.assertAlmostEqual(s["current_A"], sum(s["channel_current_A"].values()), places=14)
        self.assertAlmostEqual(s["net_P_faradaic_efficiency_pct"], 200 * F * gain["P"] / from_bulk, places=9)
        self.assertLess(s["net_P_faradaic_efficiency_pct"], s["gross_AP_charge_fraction_pct"])

    def test_no_side_reactions_matches_frozen_single_species_solver(self):
        rates = {"AP": 5e-5, "AB": 0., "PD": 0.}
        network = solve_network(nx=30, ny=12, rate_constants_m_s=rates)
        frozen = solve_transport(nx=30, ny=12, forward_rate_m_s=5e-5, reverse_rate_m_s=0)
        np.testing.assert_allclose(network["field"]["concentration_mM"]["A"], frozen["field"]["A_mol_m3"], atol=1e-10)
        np.testing.assert_allclose(network["field"]["concentration_mM"]["P"], frozen["field"]["P_mol_m3"], atol=1e-10)
        self.assertAlmostEqual(network["summary"]["current_A"], frozen["summary"]["current_A"], places=11)
        self.assertAlmostEqual(network["summary"]["net_P_faradaic_efficiency_pct"], 100, places=9)

    def test_parallel_branch_ratio_is_exact_without_overoxidation(self):
        r = solve_network(nx=30, ny=12, rate_constants_m_s={"AP": 2e-5, "AB": 1e-5, "PD": 0})
        c = r["field"]["concentration_mM"]
        np.testing.assert_allclose(np.array(c["P"]), 2 * np.array(c["B"]), atol=1e-10)
        self.assertAlmostEqual(r["summary"]["net_P_faradaic_efficiency_pct"], 100 * 2 / 3, places=9)

    def test_zero_reaction_preserves_arbitrary_inlet_mixture(self):
        r = solve_network(nx=10, ny=6, inlet_mM=(20, 10, 5, 15), rate_constants_m_s={"AP": 0, "AB": 0, "PD": 0})
        for name, feed in zip(("A", "P", "B", "D"), (20, 10, 5, 15)):
            np.testing.assert_allclose(r["field"]["concentration_mM"][name], feed, atol=1e-10)
        self.assertEqual(r["summary"]["current_A"], 0)
        self.assertIsNone(r["summary"]["net_P_faradaic_efficiency_pct"])

    def test_zero_diffusion_has_no_wall_access(self):
        r = solve_network(nx=10, ny=6, diffusivity_m2_s=0)
        self.assertEqual(r["summary"]["current_A"], 0)
        self.assertIsNone(r["field"]["wall_concentration_mM"])
        np.testing.assert_allclose(r["field"]["concentration_mM"]["A"], 50, atol=1e-10)

    def test_product_feed_degradation_retains_negative_net_FE(self):
        r = solve_network(nx=15, ny=8, inlet_mM=(0, 50, 0, 0), rate_constants_m_s={"AP": 0, "AB": 0, "PD": 2e-5})
        self.assertAlmostEqual(r["summary"]["net_P_faradaic_efficiency_pct"], -100, places=9)
        self.assertIsNone(r["summary"]["conversion_A_pct"])
        self.assertLess(r["summary"]["net_outlet_gain_mol_s"]["P"], 0)

    def test_analytic_parallel_formula(self):
        r = analytic_plug_outlet((50, 0, 0, 0), {"AP": 2e-6, "AB": 1e-6, "PD": 0}, 100000)
        expected_A = 50 * math.exp(-.3)
        self.assertAlmostEqual(r["A"], expected_A, places=12)
        self.assertAlmostEqual(r["P"], (50 - expected_A) * 2 / 3, places=12)
        self.assertAlmostEqual(r["B"], (50 - expected_A) / 3, places=12)
        self.assertAlmostEqual(r["D"], 0, places=12)

    def test_analytic_equal_sequential_rates_and_near_equal_limit(self):
        exact = analytic_plug_outlet((50, 0, 0, 0), {"AP": 1e-6, "AB": 0, "PD": 1e-6}, 100000)
        near = analytic_plug_outlet((50, 0, 0, 0), {"AP": 1e-6, "AB": 0, "PD": 1e-6 * (1 + 1e-12)}, 100000)
        self.assertAlmostEqual(exact["P"], 5 * math.exp(-.1), places=12)
        self.assertAlmostEqual(near["P"], exact["P"], places=10)

    def test_independent_sequential_plug_formula_and_grid_refinement(self):
        # dA/dtheta=-k1 A, dP/dtheta=k1 A-k2 P; k2=2k1.
        k1, k2, exposure = 2e-6, 4e-6, .06 * .012 / (450e-9 / 60)
        expected_P = 50 * k1 / (k2 - k1) * (math.exp(-k1 * exposure) - math.exp(-k2 * exposure))
        errors = []
        for nx in (40, 80, 160):
            r = solve_network(nx=nx, ny=8, velocity_profile="plug", diffusivity_m2_s=1e-4,
                              axial_diffusivity_m2_s=0, rate_constants_m_s={"AP": k1, "AB": 0, "PD": k2}, include_field=False)
            errors.append(abs(r["summary"]["outlet_mM"]["P"] - expected_P))
        self.assertTrue(errors[2] < errors[1] < errors[0])
        self.assertGreater(errors[0] / errors[1], 1.8)
        self.assertLess(errors[2], .012)

    def test_electron_counts_change_current_not_concentrations(self):
        one = solve_network(nx=10, ny=6, electrons=(1, 1, 1), include_field=False)["summary"]
        two = solve_network(nx=10, ny=6, electrons=(2, 2, 2), include_field=False)["summary"]
        self.assertEqual(one["outlet_mM"], two["outlet_mM"])
        self.assertAlmostEqual(two["current_A"], 2 * one["current_A"], places=13)

    def test_rate_potential_dependence_is_explicit_and_distinct(self):
        lower, upper = network_rates(.45), network_rates(.55)
        self.assertEqual(lower, {"AP": 5e-5, "AB": 3e-6, "PD": 8e-6})
        factors = [upper[name] / lower[name] for name in ("AP", "AB", "PD")]
        self.assertTrue(1 < factors[0] < factors[1] < factors[2])

    def test_invalid_inputs_reject_before_linear_solver(self):
        cases = [dict(flow_rate_uL_min=-1), dict(inlet_mM=(50, -1, 0, 0)), dict(inlet_mM=(0, 0, 0, 0)),
                 dict(rate_constants_m_s={"AP": -1, "AB": 0, "PD": 0}), dict(electrons=(2, 0, 2)),
                 dict(electrons=(2, 2.5, 2)), dict(nx=200, ny=100), dict(potential_index_V=-.1),
                 dict(rate_constants_m_s={"AP": 1, "AB": 0, "PD": 0, "PA": 1}), dict(inlet_mM=(float("nan"), 0, 0, 0))]
        with patch("electratwin.scripts.reaction_network.spsolve") as linear_solver:
            for kwargs in cases:
                with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                    solve_network(**kwargs)
            linear_solver.assert_not_called()

    def test_bad_numeric_solution_is_flagged_without_clipping(self):
        with patch("electratwin.scripts.reaction_network.spsolve", return_value=np.full(6 * 4 * 4, -1.0)):
            result = solve_network(nx=6, ny=4)
        self.assertFalse(result["summary"]["converged"])
        self.assertFalse(result["summary"]["nonnegative"])
        self.assertGreater(result["summary"]["conversion_A_pct"], 100)


if __name__ == "__main__":
    unittest.main()
