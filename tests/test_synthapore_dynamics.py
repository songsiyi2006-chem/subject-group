"""Small analytic invariants; the publication study is not rerun by these tests."""
import math
import unittest
import numpy as np
from synthapore.scripts.dynamics_reviewed import (
    AMU_KG, EV_J, MASS_EV_FS2_A2, KB_EV_K, curved_pes, pes_hessian,
    harmonic_energy_force, kinetic_energy, temperature, remove_com_velocity,
    verlet_step, simulate_md, minimize_endpoint, improved_tangent, neb_forces,
    optimize_neb, source_neb_probe,
)


class TestSynthaPoreDynamics(unittest.TestCase):
    def test_mass_and_force_conversion_match_SI(self):
        force, mass = 1.7, 12.011
        via_SI = (force*EV_J*1e10)/(mass*AMU_KG)*1e-20
        self.assertAlmostEqual(force/(mass*MASS_EV_FS2_A2), via_SI, places=15)

    def test_mass_weighted_COM_removal(self):
        masses = np.array([1., 12., 16.])
        velocity = remove_com_velocity(np.arange(9.).reshape(3, 3), masses)
        np.testing.assert_allclose((velocity*masses[:, None]).sum(axis=0), 0, atol=1e-13)

    def test_COM_temperature_uses_3N_minus_3(self):
        K = .5*21*KB_EV_K*300
        self.assertAlmostEqual(temperature(K, 8), 300.)
        self.assertAlmostEqual(temperature(K, 8, False), 262.5)

    def test_kinetic_energy_matches_SI(self):
        v, mass = np.array([[.02, -.01, .03], [.01, .02, .01]]), np.array([12., 16.])
        si = .5*np.sum(mass[:, None]*AMU_KG*(v*1e5)**2)/EV_J
        self.assertAlmostEqual(kinetic_energy(v, mass), si, places=14)

    def test_harmonic_translation_invariance_and_zero_net_force(self):
        q = np.arange(12.).reshape(4, 3)/10
        energy, force = harmonic_energy_force(q)
        e2, f2 = harmonic_energy_force(q+np.array([8., -5., 4.]))
        self.assertAlmostEqual(energy, e2, places=12)
        np.testing.assert_allclose(force, f2, atol=1e-13)
        np.testing.assert_allclose(force.sum(axis=0), 0, atol=1e-13)

    def test_harmonic_force_matches_finite_energy_difference(self):
        q = np.arange(12.).reshape(4, 3)/10
        _, force = harmonic_energy_force(q)
        h = 1e-6
        for i in range(4):
            for j in range(3):
                plus, minus = q.copy(), q.copy()
                plus[i,j] += h; minus[i,j] -= h
                self.assertAlmostEqual(force[i,j], -(harmonic_energy_force(plus)[0]-harmonic_energy_force(minus)[0])/(2*h), places=8)

    def test_verlet_time_reversal(self):
        q = np.array([[-.1, .02, 0.], [.1, -.02, 0.]])
        v, mass = np.array([[.001, .002, 0.], [-.001, -.002, 0.]]), np.array([12., 12.])
        _, force = harmonic_energy_force(q)
        q1, v1, f1, _ = verlet_step(q, v, force, mass, .5)
        q2, v2, _, _ = verlet_step(q1, -v1, f1, mass, .5)
        np.testing.assert_allclose(q2, q, atol=1e-15)
        np.testing.assert_allclose(v2, -v, atol=1e-15)

    def test_verlet_second_order_energy_envelope(self):
        a = simulate_md(seed=17, timestep_fs=1, duration_fs=200)["summary"]
        b = simulate_md(seed=17, timestep_fs=.5, duration_fs=200)["summary"]
        order = math.log(a["max_abs_total_energy_change_eV"]/b["max_abs_total_energy_change_eV"], 2)
        self.assertGreater(order, 1.95)
        self.assertLess(order, 2.05)

    def test_recording_times_include_initial_and_final_states(self):
        r = simulate_md(duration_fs=20., timestep_fs=.5, sample_every=10)
        self.assertEqual(r["telemetry"][0]["time_fs"], 0.)
        self.assertEqual(r["telemetry"][-1]["time_fs"], 20.)
        self.assertEqual(r["summary"]["energy_force_evaluations"], 41)

    def test_thermostat_records_post_rescaling_energy_and_temperature(self):
        r = simulate_md(duration_fs=20., thermostat="berendsen")
        last = r["telemetry"][-1]
        K = kinetic_energy(np.array(r["final_velocity_A_fs"]), np.full(8, 12.011))
        self.assertAlmostEqual(last["kinetic_eV"], K, places=14)
        self.assertAlmostEqual(last["temperature_COM_dof_K"], temperature(K, 8), places=12)

    def test_langevin_explicit_seed_and_COM_constraint(self):
        a = simulate_md(duration_fs=20., thermostat="langevin_baoab", seed=81)
        b = simulate_md(duration_fs=20., thermostat="langevin_baoab", seed=81)
        self.assertEqual(a, b)
        self.assertLess(a["summary"]["max_COM_speed_A_fs"], 1e-14)

    def test_invalid_md_inputs_and_unstable_timestep_rejected(self):
        for kwargs in ({"n_atoms": 1}, {"timestep_fs": 100}, {"duration_fs": 10.1},
                       {"mass_amu": -1}, {"thermostat": "unknown"}, {"burn_in_fs": 2000}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                simulate_md(**kwargs)

    def test_curved_pes_gradient_and_hessian(self):
        p, h = np.array([.33, -.12]), 1e-5
        _, force = curved_pes(p)
        hessian = pes_hessian(p)
        for j in range(2):
            delta = np.eye(2)[j]*h
            ep, fp = curved_pes(p+delta); em, fm = curved_pes(p-delta)
            self.assertAlmostEqual(force[j], -(ep-em)/(2*h), places=7)
            np.testing.assert_allclose(hessian[:, j], -(fp-fm)/(2*h), atol=1e-8)

    def test_analytic_minima_and_saddle_index(self):
        e, f = curved_pes(np.array([[-1., 0.], [0., .65], [1., 0.]]))
        np.testing.assert_allclose(e, [0, 1, 0], atol=1e-14)
        np.testing.assert_allclose(f, 0, atol=1e-14)
        np.testing.assert_allclose(np.linalg.eigvalsh(pes_hessian([0., .65])), [-4., 8.], atol=1e-14)

    def test_endpoint_minimization_is_convergence_checked(self):
        r = minimize_endpoint([-1.2, .15])
        self.assertTrue(r["converged"])
        np.testing.assert_allclose(r["point_A"], [-1, 0], atol=1e-9)
        self.assertLess(r["force_norm_eV_A"], 1e-11)

    def test_improved_tangent_unit_and_degenerate_rejection(self):
        t = improved_tangent(np.array([0.,0.]), np.array([1.,0.]), np.array([2.,1.]), 0., 1., 3.)
        np.testing.assert_allclose(t, np.array([1.,1.])/math.sqrt(2))
        with self.assertRaises(ValueError):
            improved_tangent(np.zeros(2), np.zeros(2), np.zeros(2), 0., 0., 0.)

    def test_climbing_image_has_no_spring_component(self):
        band = np.array([[-1.,0.], [-.3,.3], [0.,.4], [.3,.3], [1.,0.]])
        energy, true = curved_pes(band)
        f1, ci = neb_forces(band, energy, true, spring=1., climb=True)
        f2, _ = neb_forces(band, energy, true, spring=50., climb=True)
        np.testing.assert_allclose(f1[ci], f2[ci], atol=1e-14)
        np.testing.assert_allclose(f1[[0,-1]], 0, atol=1e-14)

    def test_reviewed_neb_converges_and_recomputes_final_energy(self):
        r = optimize_neb(n_images=7, tolerance=1e-4, seed=1)
        self.assertTrue(r["summary"]["converged"])
        self.assertEqual(r["summary"]["saddle_negative_hessian_eigenvalues"], 1)
        self.assertLess(r["summary"]["saddle_coordinate_error_A"], 1e-3)
        energy, _ = curved_pes(np.array(r["band_A"]))
        np.testing.assert_allclose(r["energies_eV"], energy, atol=1e-14)
        np.testing.assert_allclose(np.array(r["band_A"])[[0,-1]], [[-1,0],[1,0]], atol=1e-14)

    def test_iteration_limit_is_not_convergence_and_bad_endpoints_rejected(self):
        r = optimize_neb(max_iterations=1)
        self.assertFalse(r["summary"]["converged"])
        self.assertEqual(r["summary"]["stop_reason"], "iteration_limit")
        with self.assertRaises(ValueError):
            optimize_neb(endpoints=[[-1.2,0],[1.1,.2]])

    def test_original_class_unconditional_status_and_stale_energy(self):
        r = source_neb_probe(n_images=7, iterations=1)
        self.assertEqual(r["source_reported_status"], "converged")
        self.assertFalse(r["passes_independent_1e_5_force_check"])
        self.assertGreater(r["maximum_stale_energy_difference_eV"], .1)
        self.assertEqual(r["source_point_energy_force_calls"], 7)


if __name__ == "__main__":
    unittest.main()
