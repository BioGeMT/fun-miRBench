import unittest
import pathlib
import tempfile
from unittest.mock import patch

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.integrate import trapezoid
from sklearn.metrics import average_precision_score, precision_recall_curve, roc_auc_score

from funmirbench.subreport_metrics import (
    coverage_frontier, masks, metrics, percentile_ensemble, threshold_curves, top_weights, weighted_cdf,
)
from funmirbench.subreports import component, key_results, pdf_report, responsive_gene_rank_data, summarize_universes
from funmirbench.subreport_plots import classification, responsive_gene_heatmap
from funmirbench.benchmark import group_bootstrap, rank_last, run_benchmark, score_metrics


class DetailedMetricTests(unittest.TestCase):
    def test_report_moves_tall_plots_to_separate_pages_without_duplicate_captions(self):
        from funmirbench.report_pdf import _save_page

        frame = pd.DataFrame([{"experiment_id": "E1", "predictor": "tool", "universe": "fps",
                               "n_responsive": 40, "n_no_response": 100, "average_precision": .4,
                               "auroc": .7, "spearman": .2}])
        pages = []
        def capture(pdf, fig):
            pages.append([text.get_text() for ax in fig.axes for text in ax.texts])
            _save_page(pdf, fig)

        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            image = root / "tall.png"
            plt.imsave(image, np.ones((1000, 1000, 3)))
            with patch("funmirbench.subreports._save_page", side_effect=capture):
                pdf_report(root / "REPORT.pdf", "Test report", ["Experiment: E1"],
                           [(image, "First plot", ["First caption"]), (image, "Second plot", ["Second caption"])], frame)
        self.assertEqual(len(pages), 4)
        self.assertIn("First plot", pages[1])
        self.assertNotIn("Second plot", pages[1])
        self.assertIn("Second plot", pages[2])
        self.assertNotIn("First plot", pages[2])

    def test_ten_predictor_curve_legend_has_no_overlapping_labels(self):
        tools = {f"tool_{i}": f"Synthetic predictor {i}" for i in range(10)}
        curves = {tool: threshold_curves(np.array([3., 2., 1., 0.]), np.array([True, False, True, False]))
                  for tool in tools}
        with patch("funmirbench.subreport_plots.save") as save:
            classification(curves, {tool: {"prevalence": .5} for tool in tools}, tools,
                           pathlib.Path("unused"), "classification")
        fig = save.call_args.args[0]
        try:
            fig.canvas.draw()
            renderer = fig.canvas.get_renderer()
            boxes = [text.get_window_extent(renderer) for text in fig.legends[0].get_texts()]
            self.assertEqual(len(boxes), 10)
            for i, box in enumerate(boxes):
                self.assertFalse(any(box.overlaps(other) for other in boxes[i + 1:]))
                self.assertGreaterEqual(box.y0, fig.bbox.y0)
            for ax in fig.axes:
                xlabel = ax.xaxis.label.get_window_extent(renderer)
                self.assertFalse(any(xlabel.overlaps(box) for box in boxes))
        finally:
            plt.close(fig)

    def test_summary_intervals_match_main_and_ignore_unrelated_predictors_and_row_order(self):
        values = [.01, .03, .04, .15, .3, .5, .7, .9]
        frame = pd.DataFrame({"predictor": "b", "universe": "fps",
                              "action": ["keep"] * 6 + ["flag"] * 2,
                              "dependency_group": list("aabbcdef"),
                              **{metric: values for metric in
                                 ["average_precision", "pr_auc", "auroc", "spearman", "enrichment_score"]}})
        frame.loc[2, "spearman"] = np.nan
        baseline = summarize_universes(frame, bootstrap=200, seed=7)
        extra = frame.assign(predictor="a", universe="ips")
        expanded = summarize_universes(pd.concat([extra, frame]).sample(frac=1, random_state=1),
                                       bootstrap=200, seed=7)
        pd.testing.assert_frame_equal(baseline.reset_index(drop=True),
                                      expanded[expanded.predictor == "b"].reset_index(drop=True))
        for row in baseline.itertuples(index=False):
            selected = frame if row.subset == "all" else frame[frame.action == "keep"]
            expected = group_bootstrap(selected[row.metric], selected.dependency_group, n=200,
                                       rng=np.random.default_rng(8))
            np.testing.assert_allclose([row.ci95_low, row.ci95_high], expected)

    def test_responsive_heatmap_labels_fit_for_one_five_and_ten_predictors(self):
        names = ["TargetScan v8", "miRDB", "microT-CNN", "miRBind2-3UTR", "miRAW",
                 "Synthetic excellent", "Synthetic medium", "Synthetic poor", "Synthetic random", "Synthetic sparse"]
        for count in [1, 5, 10]:
            with self.subTest(count=count):
                tools = {f"tool_{i}": name for i, name in enumerate(names[:count])}
                matrix = pd.DataFrame(.5, index=[f"ENSG{i:011d}" for i in range(20)], columns=tools)
                with patch("funmirbench.subreport_plots.save") as save:
                    responsive_gene_heatmap(matrix, tools, pathlib.Path("unused"))
                fig = save.call_args.args[0]
                try:
                    fig.canvas.draw()
                    renderer = fig.canvas.get_renderer()
                    boxes = [label.get_window_extent(renderer) for label in fig.axes[0].get_xticklabels()]
                    self.assertFalse(any(a.overlaps(b) for a, b in zip(boxes, boxes[1:])))
                    for box in boxes:
                        self.assertGreaterEqual(box.y0, fig.bbox.y0)
                        self.assertLessEqual(box.x1, fig.bbox.x1)
                    xlabel = fig.axes[0].xaxis.label.get_window_extent(renderer)
                    self.assertGreaterEqual(xlabel.y0, fig.bbox.y0)
                finally:
                    plt.close(fig)

    def test_responsive_gene_heatmap_selects_effect_ordered_truth_positives(self):
        table = pd.DataFrame({
            "state_d05": ["responsive", "responsive", "responsive", "no_response_above_margin", "responsive"],
            "eligible": [True, True, True, True, False],
            "oriented_logFC_shrunk": [2., 3., 2., 9., 20.],
        }, index=["g1", "g2", "g3", "g4", "g5"])
        percentiles = pd.DataFrame({"tool_a": [.1, .9, .5, .8, 1.], "tool_b": [.8, .4, .7, .9, 1.]},
                                   index=table.index)

        matrix, effects = responsive_gene_rank_data(table, percentiles, limit=3)

        self.assertEqual(list(matrix.index), ["g2", "g1", "g3"])
        self.assertEqual(list(matrix.columns), ["tool_a", "tool_b"])
        self.assertEqual(effects.to_dict(), {"g2": 3., "g1": 2., "g3": 2.})
        self.assertEqual(matrix.loc["g2", "tool_a"], .9)

    def test_responsive_gene_heatmap_handles_empty_and_invalid_limit(self):
        table = pd.DataFrame({"state_d05": ["responsive"], "eligible": [False],
                              "oriented_logFC_shrunk": [1.]}, index=["g1"])
        percentiles = pd.DataFrame({"tool": [.5]}, index=table.index)
        matrix, effects = responsive_gene_rank_data(table, percentiles)
        self.assertTrue(matrix.empty)
        self.assertTrue(effects.empty)
        with self.assertRaisesRegex(ValueError, "limit"):
            responsive_gene_rank_data(table, percentiles, limit=0)

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
        self.assertNotIn("pearson", rec)

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

