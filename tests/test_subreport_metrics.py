import unittest
import pathlib
import tempfile

import numpy as np
import pandas as pd
from scipy.integrate import trapezoid
from sklearn.metrics import average_precision_score, precision_recall_curve, roc_auc_score

from funmirbench.subreport_metrics import coverage_frontier, masks, metrics, threshold_curves, top_weights, weighted_cdf
from funmirbench.subreports import component
from funmirbench.benchmark import run_benchmark


class DetailedMetricTests(unittest.TestCase):
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
