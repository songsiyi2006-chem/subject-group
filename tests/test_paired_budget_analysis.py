"""Hand-solvable synthetic cases, not chemical validation."""
import copy
import csv
import math
import tempfile
import unittest
from pathlib import Path
from scientific_analysis import paired_budget as p


class PairedBudgetTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.path = Path(self.temp.name)/'synthetic.csv'
        self.rows = [dict(seed=s, method=m, evaluation=b, full_pool_HV_fraction=v,
                          evidence_role=p.ROLE)
                     for m in p.METHODS for s in [4, 9]
                     for b, v in [(5, .2), (10, .5), (25, .9)]]
        self.save(self.rows)

    def save(self, rows):
        with self.path.open('w', newline='', encoding='utf-8') as f:
            w = csv.DictWriter(f, fieldnames=list(self.rows[0]))
            w.writeheader(); w.writerows(rows)

    def test_exact_sign_hand_cases(self):
        self.assertEqual(p.exact_sign_p(0, 0), 1)
        self.assertEqual(p.exact_sign_p(8, 0), 1/128)
        self.assertEqual(p.exact_sign_p(0, 8), 1/128)
        self.assertEqual(p.exact_sign_p(4, 4), 1)
        for a, b in [(-1, 2), (True, 1), (1.5, 2)]:
            with self.assertRaises(ValueError): p.exact_sign_p(a, b)

    def test_holm_hand_cases(self):
        self.assertEqual(p.holm_adjust([.01, .04, .03]), [.03, .06, .06])
        self.assertEqual(p.holm_adjust([1, 1]), [1, 1])
        for x in [[], [float('nan')], [-.01], [1.1]]:
            with self.assertRaises(ValueError): p.holm_adjust(x)

    def test_percentile_and_auc_hand_values(self):
        self.assertEqual(p.percentile([0, 10], .25), 2.5)
        self.assertEqual(p.normalized_auc([5, 10, 25], [0, .5, 1]), .625)
        for budgets, values in [([1, 1], [0, 1]), ([1], [1]), ([1, 2], [0])]:
            with self.assertRaises(ValueError): p.normalized_auc(budgets, values)

    def test_constant_paired_difference(self):
        r = p.paired_summary({1: 2, 2: 4}, {1: 1, 2: 3}, 100)
        self.assertEqual(r['mean_difference'], 1)
        self.assertEqual(r['bootstrap_mean_ci95'], [1, 1])
        self.assertEqual(r['sd_difference'], 0)
        self.assertEqual(r['max_absolute_leave_one_out_shift'], 0)

    def test_ties_and_seed_order(self):
        a, b = {2: 3, 1: 1}, {1: 1, 2: 3}
        r = p.paired_summary(a, b, 100)
        self.assertEqual(r['ties'], 2)
        self.assertEqual(r['sign_p_raw'], 1)
        self.assertEqual(r, p.paired_summary(dict(reversed(list(a.items()))), b, 100))

    def test_seed_alignment_and_nonfinite_rejected(self):
        for a, b in [({1: 1}, {1: 2}), ({1: 1, 2: 2}, {1: 1, 3: 2}),
                     ({1: math.inf, 2: 1}, {1: 1, 2: 2})]:
            with self.assertRaises(ValueError): p.paired_summary(a, b, 100)
        with self.assertRaises(ValueError): p.paired_summary({1: 1, 2: 2}, {1: 1, 2: 2}, 99)

    def test_leave_one_out_outlier(self):
        r = p.paired_summary(dict(enumerate([0, 0, 0, 4])), dict(enumerate([0]*4)), 100)
        self.assertEqual(r['mean_difference'], 1)
        self.assertEqual(r['leave_one_out_mean_min'], 0)
        self.assertAlmostEqual(r['leave_one_out_mean_max'], 4/3)

    def test_panel_permutation_invariance(self):
        first = p.load_panel(self.path)
        self.save(list(reversed(self.rows)))
        self.assertEqual(p.load_panel(self.path), first)

    def test_duplicate_and_unbalanced_panel_rejected(self):
        for rows in [self.rows+[self.rows[0]], self.rows[1:], self.rows[3:]]:
            self.save(rows)
            with self.assertRaises(ValueError): p.load_panel(self.path)

    def test_invalid_role_values_and_monotonicity(self):
        for field, value in [('evidence_role', 'experiment'), ('method', 'unknown'),
                             ('full_pool_HV_fraction', 'nan'), ('full_pool_HV_fraction', 1.1),
                             ('evaluation', -1), ('seed', -1)]:
            rows = copy.deepcopy(self.rows); rows[0][field] = value; self.save(rows)
            with self.assertRaises(ValueError): p.load_panel(self.path)
        rows = copy.deepcopy(self.rows); rows[1]['full_pool_HV_fraction'] = .1; self.save(rows)
        with self.assertRaises(ValueError): p.load_panel(self.path)

    def test_unreached_is_right_censored(self):
        report = p.analyze(self.path, 100, threshold=.99)
        self.assertEqual(len(report['comparisons']), 8)
        for row in report['curve_metrics']:
            self.assertEqual(row['crossing_status'], 'right_censored')
            self.assertEqual(row['crossing_lower_bound'], 25)
            self.assertIsNone(row['crossing_upper_bound'])
            self.assertIsNone(row['first_recorded_budget_reaching_threshold'])

    def test_interval_and_left_censoring(self):
        r = p.analyze(self.path, 100, threshold=.5)
        for row in r['curve_metrics']:
            self.assertEqual((row['crossing_lower_bound'], row['crossing_upper_bound']), (5, 10))
            self.assertEqual(row['crossing_status'], 'interval_censored')
        r = p.analyze(self.path, 100, threshold=.1)
        self.assertTrue(all(x['crossing_status'] == 'left_censored' for x in r['curve_metrics']))

    def test_reproducible_report_and_refuse_overwrite(self):
        a = p.analyze(self.path, 100)
        self.assertEqual(a, p.analyze(self.path, 100))
        out = Path(self.temp.name)/'report'; p.write_report(out, a)
        with self.assertRaises(FileExistsError): p.write_report(out, a)
        self.assertTrue((out/'comparisons.csv').exists())


if __name__ == '__main__':
    unittest.main()
