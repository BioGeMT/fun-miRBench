import math
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from funmirbench.evaluate import evaluate_joined_dataframe
from funmirbench.evaluate_common import (
    REFERENCE_TOOL_IDS,
    _is_publication_tool,
    _reference_scores,
)


GENE_IDS = [f"ENSG00000000{index:03d}" for index in range(1, 11)]
REFERENCE_IDS = {
    "reference_random",
    "reference_expression_high",
    "reference_expression_low",
}


def _joined(*, with_control_counts=True):
    frame = pd.DataFrame({
        "dataset_id": "dataset",
        "gene_id": GENE_IDS,
        "logFC": [-2.0, -1.8, -1.5, 0.1, 0.2, -0.1, 0.0, 0.3, -0.2, 0.4],
        "FDR": [0.001, 0.01, 0.02, 0.5, 0.6, 0.7, 0.9, 0.8, 0.4, 0.3],
        "perturbation": "Overexpression",
        "score_tool": [0.9, 0.2, 0.7, 0.1, 0.8, 0.3, 0.4, 0.5, 0.6, 0.05],
    })
    if with_control_counts:
        # GT positives are the three most highly expressed genes.
        frame["control_mean_normalized_count"] = [1000.0, 900.0, 800.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0]
    return frame


def _evaluate(joined, **kwargs):
    with tempfile.TemporaryDirectory() as temporary_directory:
        root = Path(temporary_directory)
        return evaluate_joined_dataframe(
            joined,
            plots_dir=root / "plots",
            reports_dir=root / "reports",
            fdr_threshold=0.05,
            effect_threshold=1.0,
            predictor_top_fraction=0.1,
            dataset_id="dataset",
            perturbation="Overexpression",
            **kwargs,
        )


class ReferenceBaselineTests(unittest.TestCase):
    def test_reference_tools_are_not_publication_tools(self):
        self.assertEqual(set(REFERENCE_TOOL_IDS), REFERENCE_IDS)
        for tool_id in REFERENCE_IDS:
            with self.subTest(tool_id=tool_id):
                self.assertFalse(_is_publication_tool(tool_id))
        self.assertTrue(_is_publication_tool("targetscan"))

    def test_random_reference_is_seeded_and_row_order_independent(self):
        joined = _joined()
        first = _reference_scores(joined, dataset_id="dataset", seed=7)
        second = _reference_scores(joined, dataset_id="dataset", seed=7)
        other_seed = _reference_scores(joined, dataset_id="dataset", seed=8)
        shuffled = joined.sample(frac=1.0, random_state=3)
        reordered = _reference_scores(shuffled, dataset_id="dataset", seed=7)

        column = "score_reference_random"
        pd.testing.assert_series_equal(first[column], second[column])
        self.assertFalse(first[column].equals(other_seed[column]))
        by_gene = pd.Series(first[column].to_numpy(), index=joined["gene_id"])
        reordered_by_gene = pd.Series(reordered[column].to_numpy(), index=shuffled["gene_id"])
        pd.testing.assert_series_equal(by_gene, reordered_by_gene.loc[by_gene.index])

    def test_expression_references_score_control_expression_in_both_directions(self):
        scores = _reference_scores(_joined(), dataset_id="dataset", seed=0)

        self.assertEqual(
            scores["score_reference_expression_high"].tolist(),
            _joined()["control_mean_normalized_count"].tolist(),
        )
        self.assertEqual(
            scores["score_reference_expression_low"].tolist(),
            (-_joined()["control_mean_normalized_count"]).tolist(),
        )

    def test_evaluator_reports_reference_rows_separately(self):
        result = _evaluate(_joined())

        self.assertEqual([row["tool_id"] for row in result["metric_rows"]], ["tool"])
        reference_rows = {
            row["tool_id"]: row
            for row in result["reference_metric_rows"]
            if row["population"] == "all_usable_genes"
        }
        self.assertEqual(set(reference_rows), REFERENCE_IDS)
        self.assertTrue(all(row["is_reference"] for row in reference_rows.values()))
        self.assertAlmostEqual(reference_rows["reference_expression_high"]["auroc"], 1.0)
        self.assertAlmostEqual(reference_rows["reference_expression_low"]["auroc"], 0.0)
        self.assertEqual(reference_rows["reference_random"]["coverage"], 1.0)

    def test_reference_rows_do_not_change_publication_metrics(self):
        without = _evaluate(_joined(), reference_baselines=False)
        with_references = _evaluate(_joined(), reference_baselines=True)

        self.assertEqual(without["reference_metric_rows"], [])
        pd.testing.assert_frame_equal(
            pd.DataFrame(without["metric_rows"]),
            pd.DataFrame(with_references["metric_rows"]),
        )

    def test_expression_references_need_control_counts(self):
        result = _evaluate(_joined(with_control_counts=False))

        self.assertEqual(
            sorted((row["population"], row["tool_id"]) for row in result["reference_metric_rows"]),
            [("all_usable_genes", "reference_random"), ("matched_to_tool", "reference_random")],
        )

    def test_global_reference_rows_are_labelled_all_usable_genes(self):
        result = _evaluate(_joined())

        global_rows = [row for row in result["reference_metric_rows"] if not row["matched_tool_id"]]
        self.assertEqual({row["tool_id"] for row in global_rows}, REFERENCE_IDS)
        self.assertEqual({row["population"] for row in global_rows}, {"all_usable_genes"})
        self.assertEqual({row["rows_scored"] for row in global_rows}, {10})

    def test_matched_reference_rows_use_each_predictor_scored_genes(self):
        joined = _joined()
        # The sparse predictor scores one GT positive and one background gene only.
        joined["score_sparse"] = [0.4, None, None, 0.6, None, None, None, None, None, None]
        result = _evaluate(joined)

        self.assertEqual(sorted(row["tool_id"] for row in result["metric_rows"]), ["sparse", "tool"])
        matched = [
            row
            for row in result["reference_metric_rows"]
            if row["population"] == "matched_to_sparse"
        ]
        self.assertEqual({row["tool_id"] for row in matched}, REFERENCE_IDS)
        self.assertEqual({row["matched_tool_id"] for row in matched}, {"sparse"})
        self.assertEqual({row["rows_scored"] for row in matched}, {2})
        self.assertEqual({row["positives_scored"] for row in matched}, {1})
        by_tool = {row["tool_id"]: row for row in matched}
        self.assertAlmostEqual(by_tool["reference_expression_high"]["auroc"], 1.0)
        self.assertAlmostEqual(by_tool["reference_expression_low"]["auroc"], 0.0)
        sparse_row = [row for row in result["metric_rows"] if row["tool_id"] == "sparse"][0]
        self.assertEqual(sparse_row["rows_scored"], 2)
        dense_matched = [
            row
            for row in result["reference_metric_rows"]
            if row["population"] == "matched_to_tool"
        ]
        self.assertEqual({row["rows_scored"] for row in dense_matched}, {10})

    def test_matched_reference_is_unavailable_when_expression_is_missing_on_predictor_genes(self):
        joined = _joined()
        # The sparse predictor scores genes 1, 3 and 4; control expression is missing for gene 3.
        joined["score_sparse"] = [0.4, None, 0.5, 0.6, None, None, None, None, None, None]
        joined.loc[2, "control_mean_normalized_count"] = None
        result = _evaluate(joined)

        matched = {
            row["tool_id"]: row
            for row in result["reference_metric_rows"]
            if row["population"] == "matched_to_sparse"
        }
        self.assertEqual(set(matched), REFERENCE_IDS)
        for tool_id in ("reference_expression_high", "reference_expression_low"):
            with self.subTest(tool_id=tool_id):
                row = matched[tool_id]
                self.assertEqual(
                    row["reference_unavailable_reason"],
                    "reference_missing_on_1_of_3_predictor_genes",
                )
                self.assertEqual(row["matched_genes"], 3)
                self.assertEqual(row["reference_missing_genes"], 1)
                for metric_name in ("aps", "pr_auc", "auroc", "spearman", "pearson"):
                    self.assertTrue(math.isnan(row[metric_name]), metric_name)
        random_row = matched["reference_random"]
        self.assertEqual(random_row["reference_unavailable_reason"], "")
        self.assertEqual(random_row["rows_scored"], 3)
        self.assertEqual(random_row["matched_genes"], 3)
        self.assertEqual(random_row["reference_missing_genes"], 0)
        self.assertFalse(math.isnan(random_row["auroc"]))


if __name__ == "__main__":
    unittest.main()
