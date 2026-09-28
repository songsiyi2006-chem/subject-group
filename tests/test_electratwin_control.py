"""Metric identities, independent rectangle areas and local simulator state gates."""
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "electratwin/scripts"))
import metrics_control as m


def example(**changes):
    values = dict(substrate_mw=100, coupling_partner_mw=20, product_mw=118,
                  flow_rate_mL_min=1, reactor_volume_mL=1, inlet_conc_M=1,
                  conversion_pct=50, chemoselectivity_pct=80, current_A=2 * m.FARADAY * 0.0004 / 60,
                  cell_voltage_V=3, solvent_mass_flow_g_min=1)
    values.update(changes)
    return m.green_metrics(**values)


class ElectraTwinMetricsTests(unittest.TestCase):
    def test_known_amount_units_and_charge_energy_identity(self):
        r = example()
        self.assertAlmostEqual(r["product_mol_min"], 0.0004)
        self.assertAlmostEqual(r["product_g_min"], 0.0472)
        self.assertAlmostEqual(r["space_time_yield_kg_m3_day"], 67968)
        self.assertAlmostEqual(r["faradaic_efficiency_pct"], 100)
        self.assertAlmostEqual(r["electrical_energy_kwh_kg_reactor_product"], 2 * m.FARADAY * 3 / (3600 * 118))

    def test_raw_impossible_fe_is_retained_and_flagged(self):
        r = example(current_A=m.FARADAY * 0.0004 / 60)
        self.assertAlmostEqual(r["faradaic_efficiency_pct"], 200)
        self.assertIn("faradaic_efficiency_exceeds_100_retained_not_clipped", r["flags"])

    def test_open_inventory_and_unidentified_reaction_withhold_formal_metrics(self):
        r = example()
        self.assertIsNone(r["full_process_PMI"])
        self.assertIsNone(r["atom_economy_pct"])
        self.assertAlmostEqual(r["assumed_1_to_1_atom_economy_pct"], 118 / 120 * 100)
        self.assertAlmostEqual(r["declared_boundary_mass_intensity"], 1.12 / 0.0472)

    def test_more_inputs_and_lower_recovery_raise_mass_intensity(self):
        r = example(additional_inputs_g_min={"water": 0.5, "electrolyte": 0.1}, recovery_fraction=0.5)
        self.assertAlmostEqual(r["declared_boundary_mass_intensity"], 1.72 / (0.0472 * 0.5))
        self.assertAlmostEqual(r["electrical_energy_kwh_kg_recovered_product"], 2 * r["electrical_energy_kwh_kg_reactor_product"])

    def test_zero_denominators_are_null_not_fabricated_large_values(self):
        r = example(conversion_pct=0, current_A=0)
        self.assertIsNone(r["faradaic_efficiency_pct"])
        self.assertIsNone(r["declared_boundary_mass_intensity"])
        self.assertIsNone(r["electrical_energy_kwh_kg_reactor_product"])
        self.assertEqual(r["space_time_yield_kg_m3_day"], 0)

    def test_bad_units_and_nonfinite_values_reject(self):
        for changes in [{"reactor_volume_mL": 0}, {"conversion_pct": 101}, {"current_A": -1},
                        {"product_mw": float("nan")}, {"recovery_fraction": 1.01},
                        {"electrons_transferred": 1.5}, {"additional_inputs_g_min": {"water": -1}}]:
            with self.assertRaises(ValueError):
                example(**changes)


class ElectraTwinHypervolumeTests(unittest.TestCase):
    def test_rectangle_union_known_geometry_duplicates_and_dominated_points(self):
        # Rectangles 3x1 and 1x3 overlap by 1x1: union area 5.
        p = [[3, 1], [1, 3], [1, 1], [3, 1], [-1, 10]]
        self.assertAlmostEqual(m.hypervolume_2d(p, [0, 0]), 5)
        self.assertAlmostEqual(m.hypervolume_2d(p[::-1], [0, 0]), 5)
        self.assertEqual(m.hypervolume_2d([], [0, 0]), 0)

    def test_reference_translation_and_positive_axis_scaling(self):
        p = np.array([[3, 1], [1, 3]])
        self.assertAlmostEqual(m.hypervolume_2d(p + [2, -5], [2, -5]), 5)
        self.assertAlmostEqual(m.hypervolume_2d(p * [2, 4], [0, 0]), 40)

    def test_vectorized_improvement_matches_independent_union_difference(self):
        rng = np.random.default_rng(111)
        observed = rng.uniform(0, 2, size=(8, 2))
        candidates = rng.uniform(-0.5, 3, size=(40, 2))
        predicted = m.hypervolume_improvement(candidates, observed, [0, 0])
        baseline = m.hypervolume_2d(observed, [0, 0])
        expected = [m.hypervolume_2d(np.vstack([observed, c]), [0, 0]) - baseline for c in candidates]
        np.testing.assert_allclose(predicted, expected, atol=1e-12)

    def test_pareto_directions_and_ties(self):
        self.assertEqual(m.pareto_mask([[3, 1], [1, 3], [1, 1], [3, 1]]).tolist(), [True, True, False, True])
        with self.assertRaises(ValueError):
            m.pareto_mask([[float("nan"), 1]])

    def test_budget_controls_actual_oracle_calls_and_shared_initial(self):
        x = np.array([[a, b] for a in range(4) for b in range(4)], dtype=float)
        runs = []
        for method in ["gp_mc_ehvi", "random"]:
            calls = []
            def oracle(i):
                self.assertNotIn(i, calls)
                calls.append(i)
                return [1 + x[i, 0] / 3, -(0.1 + x[i, 1] / 3)]
            run = m.sequential_campaign(x, oracle, budget=9, initial_count=3, seed=7, method=method, mc_draws=32)
            self.assertEqual(len(calls), 9)
            self.assertEqual(calls, [r["candidate_index"] for r in run])
            self.assertTrue(np.all(np.diff([r["hypervolume"] for r in run]) >= -1e-12))
            runs.append(run)
        self.assertEqual([r["candidate_index"] for r in runs[0][:3]], [r["candidate_index"] for r in runs[1][:3]])
        self.assertTrue(all(r["acquisition_MC_EHVI"] is not None for r in runs[0][3:]))

    def test_invalid_budget_or_duplicate_pool_rejects(self):
        x = np.array([[0, 0], [0, 1], [1, 0], [1, 1]])
        for budget in [1, 6]:
            with self.assertRaises(ValueError):
                m.sequential_campaign(x, lambda i: [1, -1], budget=budget, initial_count=2)
        with self.assertRaises(ValueError):
            m.sequential_campaign(np.vstack([x, x[0]]), lambda i: [1, -1], budget=4, initial_count=2)


class ElectraTwinSimulatorTests(unittest.TestCase):
    def setUp(self):
        self.s = m.SimulationOnlySCPI()

    def test_real_resource_names_are_rejected(self):
        for resource in ["TCPIP0::192.0.2.1::inst0::INSTR", "GPIB0::12::INSTR", "USB0::INSTR"]:
            with self.assertRaises(ValueError):
                m.SimulationOnlySCPI(resource)

    def test_disconnected_output_cannot_enable(self):
        self.assertEqual(self.s.send_command(":OUTP ON")["status"], "ERROR")
        self.assertFalse(self.s.output_enabled)

    def test_identity_explicitly_simulated_and_unknown_commands_fail(self):
        self.s.connect()
        self.assertIn("SIMULATOR", self.s.send_command("*IDN?")["value"])
        self.assertEqual(self.s.send_command(":UNSUPPORTED ACTION")["status"], "ERROR")
        self.assertFalse(self.s.output_enabled)

    def test_parameters_require_valid_finite_values(self):
        self.s.connect()
        for value in ["NaN", "Inf", "-0.1", "0.101", "abc"]:
            self.assertEqual(self.s.send_command(":SOUR:CURR " + value)["status"], "ERROR")
            self.assertFalse(self.s.output_enabled)
        self.assertEqual(self.s.current_A, 0)

    def test_compliance_trip_latches_fault_and_requires_safe_reset(self):
        self.s.connect()
        self.s.send_command(":SOUR:CURR 0.02")
        self.s.send_command(":OUTP ON")
        self.assertTrue(self.s.output_enabled)
        self.assertAlmostEqual(self.s.send_command(":MEAS:VOLT?")["value"], 2.398)
        self.assertEqual(self.s.send_command(":SENS:VOLT:PROT 2")["status"], "ERROR")
        self.assertFalse(self.s.output_enabled)
        self.assertEqual(self.s.send_command(":OUTP ON")["status"], "ERROR")
        self.assertEqual(self.s.send_command("*CLS")["status"], "ERROR")
        self.s.send_command(":SOUR:CURR 0")
        self.assertEqual(self.s.send_command("*CLS")["status"], "OK")
        self.assertIsNone(self.s.fault)

    def test_output_off_measurement_rejected_and_disconnect_is_safe(self):
        self.s.connect()
        self.s.send_command(":SOUR:CURR 0.01")
        self.s.send_command(":OUTP ON")
        self.s.send_command(":OUTP OFF")
        self.assertEqual(self.s.send_command(":MEAS:VOLT?")["status"], "ERROR")
        self.s.disconnect()
        self.assertFalse(self.s.connected)
        self.assertFalse(self.s.output_enabled)
        self.assertEqual(self.s.current_A, 0)


if __name__ == "__main__":
    unittest.main()
