"""Analytic-geometry and inverse-model checks, not material validation."""
import unittest

import numpy as np

from synthapore.scripts.pore_reviewed import (
    AREA_PER_MONOLAYER_CM3, bet_volume, channel_area, channel_clearance,
    clearance_cdf, disk_square_intersection_area, fit_bet, grid_points,
    minimum_image, periodic_accessible_area, periodic_disk_clearance,
    source_class_audit, source_synthetic_volume,
)


class PoreBenchmarkTests(unittest.TestCase):
    def test_circle_erosion_analytic_and_closure(self):
        self.assertAlmostEqual(channel_area("circle", 10, 2), np.pi*64)
        self.assertEqual(channel_area("circle", 10, 10), 0)
        self.assertEqual(channel_area("circle", 10, 11), 0)

    def test_hexagon_uses_apothem_not_circumradius(self):
        self.assertAlmostEqual(channel_area("hexagon", 10, 2), 2*np.sqrt(3)*64)
        vertices = np.array([[10, 10/np.sqrt(3)], [0, 20/np.sqrt(3)]])
        np.testing.assert_allclose(channel_clearance(vertices, "hexagon", 10), 0, atol=1e-13)
        self.assertGreater(channel_clearance([[0, 0]], "hexagon", 10)[0], 0)

    def test_probe_increase_cannot_increase_grid_access(self):
        points = grid_points(26, 100)
        for shape in ("circle", "hexagon"):
            clearance = channel_clearance(points, shape)
            counts = [np.count_nonzero(clearance > p) for p in (0, 1.82, 6, 9.8, 10)]
            self.assertEqual(counts, sorted(counts, reverse=True))
            self.assertEqual(counts[-1], 0)

    def test_midpoint_geometry_agrees_with_exact_area(self):
        points = grid_points(26, 320)
        for shape in ("circle", "hexagon"):
            fraction = np.mean(channel_clearance(points, shape) > 1.82)
            exact = channel_area(shape, 9.8, 1.82)/26**2
            self.assertLess(abs(fraction-exact), .002)

    def test_channel_area_scale_covariance(self):
        points = np.array([[0, 0], [2, 4], [-3, 1]])
        for shape in ("circle", "hexagon"):
            self.assertAlmostEqual(channel_area(shape, 20, 4), 4*channel_area(shape, 10, 2))
            np.testing.assert_allclose(channel_clearance(2*points, shape, 20), 2*channel_clearance(points, shape, 10))

    def test_grid_rule_midpoints_excludes_duplicate_endpoints(self):
        points = grid_points(10, 4)
        self.assertEqual(points.shape, (16, 2))
        self.assertEqual(np.min(points), -3.75)
        self.assertEqual(np.max(points), 3.75)
        self.assertEqual(np.max(grid_points(10, 4, "source_endpoints")), 5)

    def test_geometry_invalid_parameters_rejected(self):
        for call in (lambda: channel_area("circle", 10, -1), lambda: channel_area("invalid", 10),
                     lambda: grid_points(10, 1), lambda: grid_points(10, 3.5),
                     lambda: minimum_image([np.nan], 10), lambda: periodic_disk_clearance([[0, 0]], radius_A=14),
                     lambda: channel_clearance([[1, 2, 3]], "circle")):
            with self.assertRaises(ValueError):
                call()

    def test_minimum_image_far_integer_translation_invariance(self):
        delta = np.array([[12.7, -13.2], [1.2, 5.0], [-25, 28]])
        shifted = delta + np.array([[52, -78], [-26, 104], [130, -52]])
        np.testing.assert_allclose(minimum_image(delta, 26), minimum_image(shifted, 26), atol=1e-13)

    def test_minimum_image_matches_explicit_neighbor_cells(self):
        points = np.random.default_rng(72).uniform(-13, 13, (50, 2))
        shifts = np.array([(26*i, 26*j) for i in (-1, 0, 1) for j in (-1, 0, 1)])
        centers = np.array([12.4, 2.5])+shifts
        brute = np.min(np.linalg.norm(points[:, None]-centers[None], axis=2), axis=1)-3
        np.testing.assert_allclose(periodic_disk_clearance(points), brute, atol=1e-13)
        # Across the cell boundary this point is only 0.7 A from the disk center.
        self.assertAlmostEqual(periodic_disk_clearance([[-12.9, 2.5]])[0], -2.3)

    def test_periodic_area_accounts_for_overlapping_exclusion_disks(self):
        self.assertAlmostEqual(periodic_accessible_area(26, 3, 1.82), 26**2-np.pi*4.82**2)
        self.assertAlmostEqual(disk_square_intersection_area(13, 26), np.pi*13**2)
        self.assertAlmostEqual(disk_square_intersection_area(13*np.sqrt(2), 26), 26**2)
        self.assertEqual(periodic_accessible_area(26, 3, 20), 0)
        near_left = disk_square_intersection_area(13-1e-8, 26)
        near_right = disk_square_intersection_area(13+1e-8, 26)
        self.assertLess(abs(near_right-near_left), 2e-6)

    def test_clearance_cdf_normalized_monotone_and_not_diameter(self):
        for shape in ("circle", "hexagon", "periodic_disk"):
            cdf = clearance_cdf(np.linspace(-1, 20, 100), shape)
            self.assertEqual(cdf[0], 0)
            self.assertEqual(cdf[-1], 1)
            self.assertTrue(np.all(np.diff(cdf) >= -1e-14))
        self.assertAlmostEqual(clearance_cdf([4.9], "circle")[0], .75)

    def test_exact_bet_inverse_recovers_prescribed_parameters(self):
        p = np.linspace(.01, .95, 20)
        fit = fit_bet(p, bet_volume(p), .01, .30)
        self.assertAlmostEqual(fit["fitted_Vm_cm3_g"], 145, places=10)
        self.assertAlmostEqual(fit["fitted_C"], 115, places=9)
        self.assertTrue(fit["positive_parameters"])
        self.assertLess(fit["volume_RMSE_cm3_g"], 1e-10)
        self.assertTrue(fit["v_times_one_minus_p_nondecreasing"])

    def test_invalid_bet_inputs_rejected(self):
        for call in (lambda: bet_volume([0, .1]), lambda: bet_volume([.1], C=-2),
                     lambda: fit_bet([.1, .2], [10, 20]),
                     lambda: fit_bet([.1, .1, .2], [10, 11, 20]),
                     lambda: fit_bet([.1, .2, .3], [10, -5, 20])):
            with self.assertRaises(ValueError):
                call()

    def test_unphysical_fit_retained_with_rejection_flags(self):
        p = np.array([.1, .2, .3, .4])
        # A positive transformed response with a negative intercept gives negative C.
        y = -.01+.2*p
        v = p/(y*(1-p))
        fit = fit_bet(p, v, .05, .45)
        self.assertFalse(fit["positive_parameters"])
        self.assertLess(fit["fitted_C"], 0)
        self.assertIsNone(fit["apparent_area_m2_g"])

    def test_source_branch_discontinuity_and_probe_audit(self):
        audit = source_class_audit()
        self.assertFalse(audit["probe_argument_read_in_body"])
        self.assertGreater(audit["jump_cm3_g"], 100)
        values = [r["source_reported_porosity_pct"] for r in audit["calls"] if r["cell_A"] == 26]
        self.assertEqual(len(set(values)), 1)
        self.assertGreater(float(source_synthetic_volume(.35))-float(source_synthetic_volume(.35-1e-10)), 100)
        self.assertFalse(audit["PSD_computed"])

    def test_fit_window_failure_not_hidden_by_high_r_squared(self):
        p = np.linspace(.01, .95, 20)
        low = fit_bet(p, source_synthetic_volume(p), .01, .30)
        wide = fit_bet(p, source_synthetic_volume(p), .01, .95)
        self.assertLess(abs(low["fitted_Vm_cm3_g"]-145), 1e-10)
        self.assertGreater(abs(wide["fitted_Vm_cm3_g"]-145), 20)
        self.assertAlmostEqual(AREA_PER_MONOLAYER_CM3, 4.352578, places=6)


if __name__ == "__main__":
    unittest.main()
