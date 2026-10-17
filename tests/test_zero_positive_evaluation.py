import math
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from funmirbench.evaluate import evaluate_joined_dataframe
from funmirbench.validate_experiments import validate_experiments


GENE_IDS = [f"ENSG00000000{index:03d}" for index in range(1, 9)]


def _sparse_joined():
    # Genes 1-2 are GT positives; the sparse tool only scores background genes.
    return pd.DataFrame({
        "gene_id": GENE_IDS,
        "logFC": [-2.0, -2.2, 0.5, -0.4, 0.8, 0.1, -0.6, 0.3],
        "FDR": [0.01, 0.01, 0.5, 0.4, 0.3, 0.6, 0.2, 0.7],
        "perturbation": "Overexpression",
        "score_dense": [0.9, 0.8, 0.1, 0.5, 0.2, 0.3, 0.6, 0.4],
        "score_sparse": [np.nan, np.nan, 0.1, 0.5, 0.2, 0.3, 0.6, 0.4],
    })


def _zero_positive_joined():
    return pd.DataFrame({
        "gene_id": GENE_IDS,
        "logFC": [-0.9, -0.7, 0.5, -0.4, 0.8, 0.1, -0.6, 0.3],
        "FDR": [0.01, 0.01, 0.5, 0.4, 0.3, 0.6, 0.2, 0.7],
        "perturbation": "Overexpression",
        "score_tool": [0.9, 0.8, 0.1, 0.5, 0.2, 0.3, 0.6, 0.4],
    })


def _expected_spearman(joined, score_col):
    scored = joined[joined[score_col].notna()]
    return float(scored[score_col].corr(-scored["logFC"], method="spearman"))


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


def _continuous_row(result, tool_id):
    rows = [row for row in result["continuous_metric_rows"] if row["tool_id"] == tool_id]
    if len(rows) != 1:
        raise AssertionError(f"expected one continuous row for {tool_id}, got {rows}")
    return rows[0]


class ZeroPositiveEvaluationTests(unittest.TestCase):
    def test_sparse_tool_without_positives_still_gets_continuous_metric(self):
        joined = _sparse_joined()
        result = _evaluate(joined)

        # Legacy class-metric outputs are unchanged by default.
        self.assertEqual([row["tool_id"] for row in result["metric_rows"]], ["dense"])
        skipped = result["skipped_tool_rows"][0]
        self.assertEqual(skipped["tool_id"], "sparse")
        self.assertTrue(skipped["skip_reason"].startswith("No positives remain"))
        self.assertAlmostEqual(skipped["spearman"], _expected_spearman(joined, "score_sparse"))

        sparse = _continuous_row(result, "sparse")
        self.assertAlmostEqual(sparse["spearman"], _expected_spearman(joined, "score_sparse"))
        self.assertEqual(sparse["positives_scored"], 0)
        self.assertIn("No positives remain", sparse["class_metrics_na_reason"])
        dense = _continuous_row(result, "dense")
        self.assertEqual(dense["class_metrics_na_reason"], "")
        self.assertAlmostEqual(dense["spearman"], result["metric_rows"][0]["spearman"])

    def test_fully_scored_zero_positive_dataset_does_not_raise(self):
        joined = _zero_positive_joined()
        result = _evaluate(joined)

        self.assertEqual(result["metric_rows"], [])
        self.assertIn("No selected predictor", result["dataset_skip_reason"])
        row = _continuous_row(result, "tool")
        self.assertAlmostEqual(row["spearman"], _expected_spearman(joined, "score_tool"))

    def test_continuous_only_policy_keeps_zero_positive_metric_rows(self):
        joined = _zero_positive_joined()
        result = _evaluate(joined, zero_positive_policy="continuous_only")

        self.assertIsNone(result.get("dataset_skip_reason"))
        self.assertEqual(len(result["metric_rows"]), 1)
        row = result["metric_rows"][0]
        self.assertEqual(row["tool_id"], "tool")
        self.assertAlmostEqual(row["spearman"], _expected_spearman(joined, "score_tool"))
        for metric_name in ("aps", "pr_auc", "auroc", "gsea_es"):
            with self.subTest(metric=metric_name):
                self.assertTrue(math.isnan(row[metric_name]))
        self.assertIn("No positives remain", row["class_metrics_na_reason"])

    def test_unknown_zero_positive_policy_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "zero_positive_policy"):
            _evaluate(_zero_positive_joined(), zero_positive_policy="keep")

    def test_validation_exclusion_of_zero_positive_datasets_is_configurable(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            pd.DataFrame({
                "gene_id": ["ENSG000001", "ENSG000002"],
                "logFC": [0.0, 0.2],
                "control_mean_normalized_count": [10.0, 20.0],
                "FDR": [0.5, 0.5],
            }).to_csv(root / "de.tsv", sep="\t", index=False)
            pd.DataFrame({
                "id": ["dataset"],
                "mirna_name": ["hsa-miR-1"],
                "organism": ["Homo sapiens"],
                "experiment_type": ["Overexpression"],
                "de_table_path": ["de.tsv"],
            }).to_csv(root / "metadata.tsv", sep="\t", index=False)

            legacy = validate_experiments(root / "metadata.tsv", root=root)
            kept = validate_experiments(
                root / "metadata.tsv",
                root=root,
                zero_positive_policy="continuous_only",
            )

        self.assertEqual(legacy.excluded_dataset_ids, {"dataset"})
        self.assertEqual(kept.excluded_dataset_ids, set())
        self.assertEqual(kept.issues, [])
        self.assertEqual(len(kept.notices), 1)
        self.assertIn("no positive genes", kept.notices[0].message)


if __name__ == "__main__":
    unittest.main()
