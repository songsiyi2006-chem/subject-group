"""Numerical and conservation checks for the assumed independent-site SSA."""
import unittest

import numpy as np
from scipy.linalg import expm

from electrograph.scripts.kinetics_reviewed import (
    STATES, ELEMENTARY_CHARGE_C, ctmc_window, cycle_wait_cdf,
    generator_matrix, rate_constants, simulate, stationary_solution,
)


class TestElectroGraphKinetics(unittest.TestCase):
    def test_generator_conserves_probability_and_has_only_cyclic_edges(self):
        q = generator_matrix([1, 2, 3, 4, 5])
        np.testing.assert_allclose(q.sum(axis=1), 0)
        self.assertEqual(np.count_nonzero(q), 10)
        self.assertTrue(np.all(np.diag(q) < 0))

    def test_stationary_distribution_annihilates_generator(self):
        rates = rate_constants()
        s = stationary_solution(rates)
        np.testing.assert_allclose(np.array(s["occupancy"]) @ generator_matrix(rates), 0, atol=1e-12)
        self.assertAlmostEqual(sum(s["occupancy"]), 1.)
        self.assertLess(s["TOF_s_1"], 1 / (1 / 45 + 1 / 350 + 1 / 80))

    def test_equal_rate_cycle_has_uniform_stationary_and_erlang_cdf(self):
        s = stationary_solution([2] * 5)
        np.testing.assert_allclose(s["occupancy"], .2)
        self.assertAlmostEqual(s["TOF_s_1"], .4)
        import math
        expected = 1 - math.exp(-2) * sum(2**j / math.factorial(j) for j in range(5))
        self.assertAlmostEqual(cycle_wait_cdf([2] * 5, 1), expected, places=13)

    def test_analytic_window_with_stationary_start(self):
        rates = [1, 2, 3, 4, 5]
        s = stationary_solution(rates)
        w = ctmc_window(rates, 3, .7, s["occupancy"])
        np.testing.assert_allclose(w["time_average_probability"], s["occupancy"], atol=1e-12)
        self.assertAlmostEqual(w["expected_TOF_s_1"], s["TOF_s_1"], places=12)

    def test_finite_ctmc_state_flow_balances(self):
        w = ctmc_window([4, 7, 2, 8, 5], 1, .2)
        e = np.array(w["expected_reaction_counts_per_site"])
        delta = np.array(w["final_probability"]) - np.array(w["burn_in_probability"])
        np.testing.assert_allclose(delta, np.roll(e, 1) - e, atol=1e-12)

    def test_absorbing_zero_rates_preserve_state_for_complete_horizon(self):
        r = simulate([0] * 5, n_sites=7, seed=1)
        self.assertEqual(r["total_events"], 0)
        self.assertEqual(r["simulated_time_s"], 1.)
        self.assertEqual(r["time_average_coverage"], [1, 0, 0, 0, 0])
        self.assertEqual(r["stop_reason"], "absorbing")

    def test_censoring_does_not_execute_event_beyond_horizon(self):
        r = simulate([1e-12] * 5, n_sites=1, seed=7, horizon_s=.01, burn_in_s=0, record_every=1)
        self.assertTrue(r["censored_final_event"])
        self.assertEqual(r["total_events"], 0)
        self.assertEqual(r["trajectory"][-1]["time_s"], .01)
        self.assertEqual(r["final_counts"], [1, 0, 0, 0, 0])

    def test_time_weighted_occupancy_reconstructs_from_full_trajectory(self):
        r = simulate([4, 7, 2, 8, 5], n_sites=9, seed=22, burn_in_s=.17, record_every=1)
        accumulated = np.zeros(5)
        for left, right in zip(r["trajectory"][:-1], r["trajectory"][1:]):
            dt = max(0., right["time_s"] - max(left["time_s"], .17))
            accumulated += dt * np.array([left[f"n_{state}"] for state in STATES])
        np.testing.assert_allclose(accumulated / (9 * .83), r["time_average_coverage"], atol=1e-14)

    def test_event_counts_site_inventory_and_charge_are_exact(self):
        r = simulate(n_sites=12, seed=77, horizon_s=.1, burn_in_s=.02)
        self.assertEqual(sum(r["final_counts"]), 12)
        self.assertEqual(r["state_balance_residuals"], [0] * 5)
        self.assertEqual(r["electron_product_inventory_residual"], 0)
        self.assertEqual(sum(r["reaction_counts_full"]), r["total_events"])
        self.assertAlmostEqual(r["current_for_explicit_sites_A"], r["electrons_observed"] * ELEMENTARY_CHARGE_C / .08)

    def test_seed_reproducibility_and_no_global_rng_dependency(self):
        first = simulate(n_sites=8, seed=12, horizon_s=.05, burn_in_s=0)
        np.random.seed(188)
        second = simulate(n_sites=8, seed=12, horizon_s=.05, burn_in_s=0)
        self.assertEqual(first, second)

    def test_heterogeneous_classes_conserve_their_own_sizes(self):
        r = simulate(rate_classes=[{"sites": 3, "rates": [1, 2, 3, 4, 5]},
                                   {"sites": 7, "rates": [5, 4, 3, 2, 1]}], seed=88)
        self.assertEqual([sum(x) for x in r["final_counts_by_class"]], [3, 7])
        self.assertEqual(r["state_balance_residuals"], [0] * 5)

    def test_fixed_event_diagnostic_returns_actual_post_event_state(self):
        r = simulate([1] * 5, n_sites=1, seed=1, horizon_s=100, burn_in_s=0, max_events=3, record_every=1)
        self.assertEqual(r["stop_reason"], "event_limit")
        self.assertEqual(r["total_events"], 3)
        self.assertEqual(r["final_counts"], [0, 0, 0, 1, 0])
        self.assertEqual(r["trajectory"][-1]["n_intermediate"], 1)

    def test_invalid_inputs_rejected(self):
        cases = [{"rates": [-1, 2, 3, 4, 5]}, {"rates": [1, 2]}, {"n_sites": 0},
                 {"horizon_s": .2, "burn_in_s": .2}, {"seed": 1.2}, {"max_events": 0},
                 {"initial_counts": [[1, 0, 0, 0, 0]]}, {"horizon_s": float("inf")}]
        for kwargs in cases:
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                simulate(**kwargs)

    def test_cycle_cdf_monotone_normalized(self):
        values = [cycle_wait_cdf(rate_constants(), t) for t in [0, .01, .03, .1, 1]]
        self.assertEqual(values[0], 0.)
        self.assertTrue(np.all(np.diff(values) >= 0))
        self.assertAlmostEqual(values[-1], 1., places=12)

    def test_analytic_rate_elasticity_matches_small_perturbation(self):
        rates = rate_constants()
        baseline = stationary_solution(rates)
        for i in range(5):
            plus, minus = rates.copy(), rates.copy()
            plus[i] *= np.exp(1e-4)
            minus[i] *= np.exp(-1e-4)
            derivative = (np.log(stationary_solution(plus)["TOF_s_1"]) - np.log(stationary_solution(minus)["TOF_s_1"])) / 2e-4
            self.assertAlmostEqual(derivative, baseline["rate_log_elasticity"][i], places=8)


if __name__ == "__main__":
    unittest.main()
