import unittest
import pathlib
import tempfile

import numpy as np
import pandas as pd
from scipy.integrate import trapezoid
from sklearn.metrics import average_precision_score, precision_recall_curve, roc_auc_score

from funmirbench.subreport_metrics import (
    coverage_frontier, masks, metrics, percentile_ensemble, threshold_curves, top_weights, weighted_cdf,
)
from funmirbench.subreports import component, key_results
from funmirbench.benchmark import rank_last, run_benchmark, score_metrics


class DetailedMetricTests(unittest.TestCase):
    def test_report_key_results_use_primary_fps_and_explain_unavailable_metrics(self):
        frame = pd.DataFrame([
            {"universe": "fps", "experiment_id": "E1", "predictor": "a", "average_precision": .5,
             "auroc": .7, "spearman": .2, "classification_note": "", "spearman_na_reason": ""},
            {"universe": "fps", "experiment_id": "E1", "predictor": "b", "average_precision": .6,
             "auroc": .8, "spearman": .1, "classification_note": "", "spearman_na_reason": ""},
            {"universe": "ips", "experiment_id": "E1", "predictor": "b", "average_precision": .99,
             "auroc": .99, "spearman": .99, "classification_note": "", "spearman_na_reason": ""},
        ])
        lines = key_results(frame)
        self.assertIn("Highest mean FPS AP: b (0.600).", lines)
        self.assertIn("Highest mean FPS AUROC: b (0.800).", lines)
        self.assertIn("Highest mean FPS AP: Tool B (0.600).", key_results(frame, {"b": "Tool B"}))
        unavailable = frame.iloc[[0]].assign(average_precision=np.nan, auroc=np.nan, spearman=np.nan,
                                              classification_note="Needs both labelled states",
                                              spearman_na_reason="Needs 3 finite genes")
        lines = key_results(unavailable)
        self.assertTrue(any("Classification unavailable" in line for line in lines))
        self.assertTrue(any("Continuous metrics unavailable" in line for line in lines))

    def test_curves_match_sklearn_with_ties_and_have_no_gene_order_effect(self):
        rng = np.random.default_rng(7)
        scores, pos = rng.integers(-2, 4, 150), rng.random(150) < .25
        curve = threshold_curves(scores, pos)
        ap = np.sum(np.diff(np.r_[0, curve.recall]) * curve.precision)
        self.assertAlmostEqual(ap, average_precision_score(pos, scores))
        self.assertAlmostEqual(trapezoid(np.r_[0, curve.recall], np.r_[0, curve.fpr]), roc_auc_score(pos, scores))
        reverse = threshold_curves(scores[::-1], pos[::-1])
        pd.testing.assert_frame_equal(curve, reverse)
        constant = threshold_curves(np.zeros(100), np.arange(100) < 10)
        self.assertEqual(constant.running_es.iloc[0], 0)
        self.assertEqual(constant.precision.iloc[0], .1)

    def test_ap_and_trapezoidal_pr_auc_are_reported_as_distinct_summaries(self):
        table = pd.DataFrame({"state_d05": ["responsive", "no_response_above_margin",
                                             "responsive", "no_response_above_margin"],
                              "eligible": [True] * 4, "oriented_logFC_shrunk": [2., 0., 1., -1.]})
        raw = pd.DataFrame({"tool": [3., 2., 2., 1.]})
        rec, _ = metrics(table, raw, raw.tool, "tool", "fps", 1)
        self.assertAlmostEqual(rec["average_precision"], 5 / 6)
        self.assertAlmostEqual(rec["pr_auc"], 11 / 12)
        self.assertAlmostEqual(rec["prevalence"], .5)

    def test_fractional_top_selection_preserves_budget_and_weighted_distribution(self):
        rank = pd.Series([2, 1, 1, 1, 0], index=list("abcde"))
        w = top_weights(rank, 2)
        self.assertAlmostEqual(w.sum(), 2)
        self.assertEqual(w["a"], 1)
        self.assertAlmostEqual(w["b"], 1 / 3)
        self.assertAlmostEqual(w["c"], 1 / 3)
        c = weighted_cdf(np.array([1., 2., 3., np.nan, 5.]), w.to_numpy())
        self.assertAlmostEqual(c.cdf.iloc[-1], 1)
        self.assertAlmostEqual(c.weight.sum(), 1 + 2 / 3)
        self.assertEqual(top_weights(rank, 100).sum(), 5)

    def test_empty_single_state_and_nonfinite_effects(self):
        self.assertTrue(threshold_curves(np.array([1, 2]), np.array([True, True])).empty)
        self.assertTrue(weighted_cdf(np.array([np.nan])).empty)
        self.assertTrue(top_weights(pd.Series(dtype=float), 10).empty)

    def test_universes_and_unresolved_are_separate_from_negatives(self):
        table = pd.DataFrame({"state_d05": ["responsive", "no_response_above_margin", "unresolved", "responsive"],
                              "eligible": [True, True, True, False], "oriented_logFC_shrunk": [2., 0., 1., np.nan]})
        raw = pd.DataFrame({"a": [3., 2., 1., np.nan], "b": [3., np.nan, 1., np.nan]})
        ranking = raw.a.fillna(0)
        rec, curve = metrics(table, raw, ranking, "a", "fps", 1)
        self.assertEqual(rec["n_unresolved"], 1)
        self.assertEqual(rec["n_usable"], 3)
        self.assertEqual(rec["n_no_response"], 1)
        self.assertEqual(int(curve.selected.iloc[-1]), 3)
        self.assertEqual(masks(raw, "a", "ips").sum(), 2)
        self.assertEqual(masks(raw, "a", "own_scored").sum(), 3)
        rec, curve = metrics(table, raw, ranking, "a", "ips", 1)
        self.assertTrue(curve.empty)
        self.assertTrue(np.isnan(rec["auroc"]))
        self.assertIn("Needs at least", rec["classification_note"])
        self.assertIn("at least 3", rec["spearman_na_reason"])
        self.assertEqual(rec["pearson_na_reason"], rec["spearman_na_reason"])

    def test_fps_metrics_match_main_with_missing_scores_ties_and_unresolved_genes(self):
        genes = pd.Index([f"g{i}" for i in range(8)])
        table = pd.DataFrame({
            "state_d05": ["responsive", "responsive", "no_response_above_margin", "no_response_above_margin",
                          "unresolved", "unresolved", "responsive", "no_response_above_margin"],
            "eligible": [True, True, True, True, True, True, False, True],
            "oriented_logFC_shrunk": [2., 1.5, 0., -.2, .6, .3, np.nan, -.5],
            "expression_tenth": [1, 1, 2, 2, 3, 3, 4, 4],
        }, index=genes)
        raw = pd.DataFrame({"tool": [4., 4., 1., np.nan, 2., np.nan, 8., 0.]}, index=genes)
        ranked = rank_last(raw.tool, genes)
        main = score_metrics(ranked, table, min_per_state=1)
        detailed, curve = metrics(table, raw, ranked, "tool", "fps", 1)
        self.assertAlmostEqual(detailed["auroc"], main["auroc_d05"])
        self.assertAlmostEqual(detailed["spearman"], main["spearman_shrunk"])
        self.assertEqual(detailed["n_unresolved"], 2)
        self.assertEqual(detailed["n_responsive"], 3)
        self.assertEqual(detailed["n_no_response"], 3)
        self.assertEqual(curve.threshold.iloc[-1], -1.)

    def test_metric_helpers_are_invariant_to_gene_order_with_ties(self):
        from funmirbench.metric_helpers import auroc, spearman

        scores = np.array([3., 3., 2., 1., 1., -2.])
        positive = np.array([True, False, True, False, True, False])
        effects = np.array([4., 0., 2., -1., 1., -3.])
        order = np.array([5, 2, 0, 4, 1, 3])
        self.assertAlmostEqual(auroc(scores, positive), auroc(scores[order], positive[order]))
        self.assertAlmostEqual(spearman(scores, effects), spearman(scores[order], effects[order]))

    def test_combination_percentiles_leave_missing_pairs_for_the_fps_floor(self):
        raw = pd.DataFrame({"a": [3., 1., np.nan, np.nan], "b": [np.nan, 2., 1., np.nan]})
        percentiles = raw.rank(method="average", pct=True)
        combined = percentile_ensemble(percentiles, ["a", "b"])
        self.assertTrue(np.isnan(combined.iloc[3]))
        ranked = rank_last(combined, raw.index)
        self.assertLess(ranked.iloc[3], ranked.iloc[:3].min())
        all_missing = percentile_ensemble(percentiles[["a"]].iloc[[3]], ["a"])
        self.assertTrue(all_missing.isna().all())
        self.assertTrue((rank_last(all_missing, all_missing.index) == 0).all())

    def test_safe_paths(self):
        for name in ["../outside", "/absolute", "..", "", "a" * 200]:
            safe = component(name)
            self.assertNotIn("/", safe)
            self.assertNotIn(safe, ["", ".", ".."])
        self.assertNotEqual(component("a/b"), component("a_b"))

    def test_frontier_discards_dominated_points_at_the_same_coverage(self):
        points = pd.DataFrame({"coverage": [1., 1., .8, .7], "average_precision": [.2, .6, .5, .8]})
        expected = {(1., .6), (.7, .8)}
        for data in [points, points.iloc[::-1]]:
            frontier = coverage_frontier(data)
            self.assertEqual(set(zip(frontier.coverage, frontier.average_precision)), expected)

    def test_invalid_report_settings_fail_before_writing_outputs(self):
        import yaml
        with tempfile.TemporaryDirectory() as tmp:
            path = pathlib.Path(tmp) / "config.yaml"
            for settings in [False, {"detailed": "false"}, {"top_n": 0}, {"top_n": 1.5}, {"top_n": True}, {"top_n": float("inf")}]:
                with self.subTest(settings=settings):
                    path.write_text(yaml.safe_dump({"reports": settings}))
                    with self.assertRaisesRegex(ValueError, "reports"):
                        run_benchmark(path)
            self.assertEqual(list(path.parent.iterdir()), [path])
