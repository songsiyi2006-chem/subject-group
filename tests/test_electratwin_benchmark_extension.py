"""Budget, leakage and independent error-statistic checks for the cached extension."""
import csv
import json
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "electratwin/scripts"))
import benchmark_extension as b


class BenchmarkBudgetTests(unittest.TestCase):
    def setUp(self):
        self.x = np.array([[q, eta] for q in range(4) for eta in range(4)], dtype=float)
        self.y = np.column_stack([0.5 + self.x[:, 0] / 10, -0.3 - self.x[:, 1] / 10])

    def test_cached_oracle_forbids_duplicates_and_budget_overrun(self):
        oracle = b.BudgetedCachedOracle(self.y, 2)
        oracle(1)
        with self.assertRaises(RuntimeError):
            oracle(1)
        oracle(2)
        with self.assertRaises(RuntimeError):
            oracle(3)
        self.assertEqual(oracle.calls, [1, 2])

    def test_three_methods_share_initial_and_use_exact_budget_without_repeats(self):
        initial = []
        for method in b.METHODS:
            oracle = b.BudgetedCachedOracle(self.y, 9)
            if method == "maximin_spacefill":
                rows = b.maximin_campaign(self.x, oracle, budget=9, initial_count=5, seed=7)
            else:
                rows = b.sequential_campaign(self.x, oracle, budget=9, initial_count=5, seed=7, method=method, mc_draws=16)
            ids = [r["candidate_index"] for r in rows]
            self.assertEqual(ids, oracle.calls)
            self.assertEqual(len(set(ids)), 9)
            self.assertTrue(np.all(np.diff([r["hypervolume"] for r in rows]) >= -1e-12))
            initial.append(ids[:5])
        self.assertEqual(initial[0], initial[1])
        self.assertEqual(initial[1], initial[2])

    def test_maximin_selection_ignores_permuted_and_extreme_objectives(self):
        altered = self.y[::-1] * [100, 0.01]
        first = b.maximin_campaign(self.x, b.BudgetedCachedOracle(self.y, 12), budget=12, initial_count=5, seed=13)
        second = b.maximin_campaign(self.x, b.BudgetedCachedOracle(altered, 12), budget=12, initial_count=5, seed=13)
        self.assertEqual([r["candidate_index"] for r in first], [r["candidate_index"] for r in second])
        self.assertEqual([r["acquisition_input_maximin_distance"] for r in first], [r["acquisition_input_maximin_distance"] for r in second])

    def test_budget_prefixes_do_not_change_selection(self):
        short = b.maximin_campaign(self.x, b.BudgetedCachedOracle(self.y, 8), budget=8, initial_count=5, seed=9)
        long = b.maximin_campaign(self.x, b.BudgetedCachedOracle(self.y, 12), budget=12, initial_count=5, seed=9)
        self.assertEqual(short, long[:8])

    def test_descriptive_t_handles_identical_initial_differences(self):
        stats = b.descriptive_t(np.zeros(64))
        self.assertEqual(stats["mean"], 0)
        self.assertEqual(stats["lower_95pct_t"], 0)
        self.assertEqual(stats["upper_95pct_t"], 0)
        self.assertEqual(stats["ties"], 64)


class BenchmarkHoldoutTests(unittest.TestCase):
    def test_all_26_splits_disjoint_complete_and_blocked(self):
        x, _ = b.load_pool()
        splits = b.make_holdout_splits(x)
        self.assertEqual(len(splits), 26)
        self.assertEqual(len({s["split_id"] for s in splits}), 26)
        for split in splits:
            train, test = set(split["train_ids"]), set(split["test_ids"])
            self.assertEqual(len(train), 54)
            self.assertEqual(len(test), 27)
            self.assertFalse(train & test)
            self.assertEqual(train | test, set(range(81)))
            if split["kind"] in {"block_flow", "block_eta"}:
                axis = 0 if split["kind"] == "block_flow" else 1
                self.assertFalse(set(x[list(train), axis]) & set(x[list(test), axis]))

    def test_error_statistics_match_hand_computation(self):
        result = b.prediction_statistics([1, 2, 3], [1, 3, 1], [1, 0.5, 1])
        self.assertAlmostEqual(result["MAE"], 1)
        self.assertAlmostEqual(result["RMSE"], np.sqrt(5 / 3))
        self.assertAlmostEqual(result["mean_prediction_minus_truth"], -1 / 3)
        self.assertAlmostEqual(result["latent_interval_coverage_95pct"], 1 / 3)
        self.assertAlmostEqual(result["standardized_residual_mean_truth_minus_prediction"], 0)
        self.assertAlmostEqual(result["standardized_residual_RMS"], np.sqrt(8 / 3))
        self.assertAlmostEqual(result["max_abs_standardized_residual"], 2)

    def test_quadratic_baseline_reconstructs_independent_linear_example(self):
        train = np.array([[a, c] for a in [-2, -1, 0, 1, 2] for c in [-2, -1, 0, 1, 2]], dtype=float)
        test = np.array([[2.2, 0.5], [-1.3, 0.7], [0.2, -2.3]])
        def truth(x):
            return np.column_stack([2 + x[:, 0] - 0.5 * x[:, 1], -3 + 0.1 * x[:, 0] + 0.2 * x[:, 1]])
        models, meta = b.fit_predict_models(train, truth(train), test)
        np.testing.assert_allclose(models["quadratic_regression"]["mean"], truth(test), atol=1e-12)
        np.testing.assert_allclose(models["training_mean"]["mean"], np.tile([2, -3], (3, 1)), atol=1e-12)
        self.assertEqual(meta["quadratic_design_rank"], 6)

    def test_test_inputs_do_not_change_preprocessing_or_training_coefficients(self):
        train = np.array([[a, c] for a in range(3) for c in range(3)], dtype=float)
        y = np.column_stack([train[:, 0] ** 2 + train[:, 1], -2 + 0.1 * train[:, 0]])
        test = np.array([[0.5, 1.5], [1.3, 0.7]])
        first, first_meta = b.fit_predict_models(train, y, test)
        second, second_meta = b.fit_predict_models(train, y, np.vstack([test, [10000, -10000]]))
        for key in ["input_min_from_training", "input_max_from_training", "input_scale_from_training",
                    "GP_outcome_mean_from_training", "GP_outcome_scale_from_training", "quadratic_coefficients_training_only"]:
            self.assertEqual(first_meta[key], second_meta[key])
        for name in b.MODELS:
            np.testing.assert_allclose(first[name]["mean"], second[name]["mean"][:2], atol=1e-12)
        self.assertFalse(second_meta["heldout_targets_seen_by_fit"])

    def test_invalid_error_vectors_and_standard_deviations_reject(self):
        with self.assertRaises(ValueError):
            b.prediction_statistics([1, 2], [1])
        with self.assertRaises(ValueError):
            b.prediction_statistics([1, 2], [1, 2], [1, -1])
        zero = b.prediction_statistics([1], [1], [0])
        self.assertEqual(zero["latent_interval_coverage_95pct"], 1)
        self.assertEqual(zero["zero_latent_sd_count"], 1)
        self.assertIsNone(zero["standardized_residual_RMS"])


class BenchmarkSavedResultsTests(unittest.TestCase):
    @unittest.skipUnless((b.OUT / "campaign_evaluations.csv").exists(), "Full benchmark has not been run yet")
    def test_saved_4800_calls_and_original_240_early_prefixes(self):
        rows = b.read_csv(b.OUT / "campaign_evaluations.csv")
        self.assertEqual(len(rows), 4800)
        groups = {}
        for row in rows:
            groups.setdefault((row["seed"], row["method"]), []).append(row)
        self.assertEqual(len(groups), 192)
        for run in groups.values():
            self.assertEqual([int(r["evaluation"]) for r in run], list(range(1, 26)))
            self.assertEqual(len({r["candidate_index"] for r in run}), 25)
        check = b.compare_saved_prefix(rows)
        self.assertTrue(check["all_original_240_records_reproduced"])

    @unittest.skipUnless((b.OUT / "holdout_predictions.csv").exists(), "Full holdouts have not been run yet")
    def test_saved_predictions_belong_only_to_each_split_test_membership(self):
        metadata = json.loads((b.OUT / "holdout_split_metadata.json").read_text(encoding="utf-8"))
        test_sets = {r["split_id"]: set(r["test_ids"]) for r in metadata}
        predictions = b.read_csv(b.OUT / "holdout_predictions.csv")
        self.assertEqual(len(predictions), 4212)
        for row in predictions:
            self.assertIn(int(row["candidate_index"]), test_sets[row["split_id"]])


if __name__ == "__main__":
    unittest.main()
