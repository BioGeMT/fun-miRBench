import unittest

import pandas as pd

from funmirbench.reference_audit import (
    assign_expression_bins,
    prepare_de_frame,
    summarize_expression_bins,
    summarize_gene_response_propensity,
)


class ReferenceAuditTests(unittest.TestCase):
    def test_overexpression_positives_use_negative_logfc(self):
        df = pd.DataFrame({
            "gene_id": ["ENSG1", "ENSG2", "ENSG3"],
            "logFC": [-1.5, 1.8, -0.5],
            "FDR": [0.01, 0.01, 0.01],
            "control_mean": [10.0, 20.0, 30.0],
        })

        prepared, summary = prepare_de_frame(
            df,
            dataset_id="dataset",
            mirna="hsa-miR-1",
            perturbation="Overexpression",
            fdr_threshold=0.05,
            effect_threshold=1.0,
        )

        self.assertEqual(summary["positives"], 1)
        self.assertTrue(prepared.loc[prepared["gene_id"].eq("ENSG1"), "is_current_positive"].iloc[0])
        self.assertFalse(prepared.loc[prepared["gene_id"].eq("ENSG2"), "is_current_positive"].iloc[0])

    def test_knockout_positives_use_positive_logfc(self):
        df = pd.DataFrame({
            "gene_id": ["ENSG1", "ENSG2"],
            "logFC": [-1.5, 1.8],
            "FDR": [0.01, 0.01],
        })

        prepared, summary = prepare_de_frame(
            df,
            dataset_id="dataset",
            mirna="hsa-miR-1",
            perturbation="Knockout",
            fdr_threshold=0.05,
            effect_threshold=1.0,
        )

        self.assertEqual(summary["positives"], 1)
        self.assertFalse(prepared.loc[prepared["gene_id"].eq("ENSG1"), "is_current_positive"].iloc[0])
        self.assertTrue(prepared.loc[prepared["gene_id"].eq("ENSG2"), "is_current_positive"].iloc[0])

    def test_normalized_control_mean_is_preferred_for_baseline_expression(self):
        df = pd.DataFrame({
            "gene_id": ["ENSG1"],
            "logFC": [-1.5],
            "FDR": [0.01],
            "raw_control_mean": [10.0],
            "normalized_control_mean": [12.5],
        })

        prepared, summary = prepare_de_frame(
            df,
            dataset_id="dataset",
            mirna="hsa-miR-1",
            perturbation="Overexpression",
            fdr_threshold=0.05,
            effect_threshold=1.0,
        )

        self.assertEqual(summary["baseline_expression_column"], "normalized_control_mean")
        self.assertEqual(prepared["baseline_expression"].iloc[0], 12.5)

    def test_expression_bins_keep_zero_separate(self):
        bins = assign_expression_bins(
            pd.Series([0.0, 1.0, 2.0, 100.0, None]),
            bins=2,
        )

        self.assertEqual(bins.iloc[0], "zero")
        self.assertEqual(bins.iloc[4], "missing")
        self.assertIn(bins.iloc[1], {"q1_positive", "q2_positive"})
        self.assertIn(bins.iloc[3], {"q1_positive", "q2_positive"})

    def test_expression_bin_summary_counts_current_positives(self):
        df = pd.DataFrame({
            "gene_id": ["ENSG1", "ENSG2", "ENSG3"],
            "logFC": [-2.0, -0.2, 1.5],
            "FDR": [0.01, 0.9, 0.01],
            "control_mean": [0.0, 10.0, 100.0],
        })
        prepared, _ = prepare_de_frame(
            df,
            dataset_id="dataset",
            mirna="hsa-miR-1",
            perturbation="Overexpression",
            fdr_threshold=0.05,
            effect_threshold=1.0,
        )

        summary = summarize_expression_bins(prepared, bins=2)

        self.assertEqual(int(summary["rows"].sum()), 3)
        self.assertEqual(int(summary["positives"].sum()), 1)

    def test_gene_response_propensity_uses_assessable_rows(self):
        frame = pd.DataFrame({
            "gene_id": ["ENSG1", "ENSG1", "ENSG2"],
            "dataset_id": ["d1", "d2", "d1"],
            "is_assessable": [True, True, False],
            "is_current_positive": [True, False, False],
            "expected_effect": [2.0, 0.5, 3.0],
            "logFC": [-2.0, -0.5, -3.0],
        })

        propensity = summarize_gene_response_propensity(frame)

        self.assertEqual(propensity["gene_id"].tolist(), ["ENSG1"])
        self.assertEqual(propensity["assessable_rows"].iloc[0], 2)
        self.assertEqual(propensity["positive_fraction"].iloc[0], 0.5)


if __name__ == "__main__":
    unittest.main()
