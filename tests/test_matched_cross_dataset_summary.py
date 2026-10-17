import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score

from funmirbench.cross_dataset import write_cross_dataset_summaries, write_reference_baseline_tables


GENE_IDS = [f"ENSG00000000{index:03d}" for index in range(1, 7)]


CONTROL_COUNTS = [600.0, 500.0, 400.0, 300.0, 200.0, 100.0]


def _joined(dataset_id, *, score_a, score_b=None):
    frame = pd.DataFrame({
        "dataset_id": dataset_id,
        "mirna": "hsa-miR-1",
        "perturbation": "Overexpression",
        "gene_id": GENE_IDS,
        "logFC": [-2.0, -1.5, 0.1, -2.5, 0.3, 0.2],
        "FDR": [0.01, 0.01, 0.5, 0.01, 0.6, 0.7],
        "control_mean_normalized_count": CONTROL_COUNTS,
        "score_a": score_a,
    })
    if score_b is not None:
        frame["score_b"] = score_b
    return frame


def _metric_row(dataset_id, tool_id, aps):
    return {
        "dataset_id": dataset_id,
        "mirna": "hsa-miR-1",
        "cell_line": "",
        "perturbation": "Overexpression",
        "geo_accession": "",
        "tool_id": tool_id,
        "coverage": 1.0,
        "positive_coverage": 1.0,
        "aps": aps,
        "pr_auc": aps,
        "spearman": 0.5,
        "auroc": aps,
    }


class MatchedCrossDatasetSummaryTests(unittest.TestCase):
    def setUp(self):
        # Tool b only scores genes 1-3 and 5; dataset ds3 was evaluated for tool a only.
        self.joined_frames = [
            _joined(
                "ds1",
                score_a=[0.9, 0.2, 0.3, 0.8, 0.1, 0.4],
                score_b=[0.7, 0.6, 0.1, np.nan, 0.2, np.nan],
            ),
            _joined(
                "ds2",
                score_a=[0.6, 0.5, 0.9, 0.7, 0.1, 0.2],
                score_b=[0.9, 0.1, 0.8, np.nan, 0.3, np.nan],
            ),
            _joined("ds3", score_a=[0.9, 0.8, 0.1, 0.7, 0.2, 0.3]),
        ]
        self.metric_rows = [
            _metric_row("ds1", "a", 0.9),
            _metric_row("ds1", "b", 0.8),
            _metric_row("ds2", "a", 0.7),
            _metric_row("ds2", "b", 0.6),
            _metric_row("ds3", "a", 1.0),
        ]

    def _write(self, root, **kwargs):
        return write_cross_dataset_summaries(
            self.metric_rows,
            root / "tables",
            root / "plots",
            joined_frames=self.joined_frames,
            fdr_threshold=0.05,
            effect_threshold=1.0,
            **kwargs,
        )

    def test_existing_own_population_summary_is_kept(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            outputs = self._write(Path(temporary_directory))
            summary = pd.read_csv(
                outputs["tables"]["cross_dataset_predictor_summary"], sep="\t"
            ).set_index("tool_id")

        self.assertEqual(summary.loc["a", "aps_count"], 3)
        self.assertAlmostEqual(summary.loc["a", "aps_mean"], (0.9 + 0.7 + 1.0) / 3.0)

    def test_matched_summary_uses_common_genes_and_datasets(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            outputs = self._write(Path(temporary_directory))
            summary = pd.read_csv(
                outputs["tables"]["cross_dataset_predictor_summary_matched"], sep="\t"
            ).set_index("tool_id")
            per_dataset = pd.read_csv(
                outputs["tables"]["matched_metrics_per_experiment"], sep="\t"
            )

        self.assertEqual(sorted(summary.index), ["a", "b"])
        self.assertEqual(set(summary["population"]), {"matched"})
        self.assertEqual(set(summary["matched_tool_ids"]), {"a,b"})
        self.assertEqual(set(summary["matched_dataset_count"]), {2})
        self.assertEqual(summary.loc["a", "aps_comparisons"], 2)
        self.assertEqual(sorted(per_dataset["dataset_id"].unique()), ["ds1", "ds2"])

        common = self.joined_frames[0].dropna(subset=["score_a", "score_b"])
        y_true = ((-common["logFC"] > 1.0) & (common["FDR"] < 0.05)).astype(int)
        expected_aps = average_precision_score(y_true, common["score_a"])
        ds1_a = per_dataset[(per_dataset["dataset_id"] == "ds1") & (per_dataset["tool_id"] == "a")]
        self.assertEqual(int(ds1_a["common_rows"].iloc[0]), 4)
        self.assertAlmostEqual(float(ds1_a["aps"].iloc[0]), expected_aps)


def _expected_effect_spearman(frame, score_col):
    return float(frame[score_col].corr(-frame["logFC"], method="spearman"))


class MatchedCohortTests(unittest.TestCase):
    """Every metric is summarized over one dataset cohort shared by all compared tools."""

    def setUp(self):
        # Same genes in both datasets; tool b is constant in ds2, so its correlations
        # are undefined there while its class metrics are still defined.
        self.joined_frames = [
            _joined(
                "ds1",
                score_a=[0.9, 0.2, 0.3, 0.8, 0.1, 0.4],
                score_b=[0.7, 0.6, 0.1, 0.5, 0.2, 0.3],
            ),
            _joined(
                "ds2",
                score_a=[0.6, 0.5, 0.9, 0.7, 0.1, 0.2],
                score_b=[0.5, 0.5, 0.5, 0.5, 0.5, 0.5],
            ),
        ]
        self.metric_rows = [
            _metric_row("ds1", "a", 0.9),
            _metric_row("ds1", "b", 0.8),
            _metric_row("ds2", "a", 0.7),
            _metric_row("ds2", "b", 0.6),
        ]

    def _write(self, root, **kwargs):
        return write_cross_dataset_summaries(
            self.metric_rows,
            root / "tables",
            root / "plots",
            joined_frames=self.joined_frames,
            fdr_threshold=0.05,
            effect_threshold=1.0,
            **kwargs,
        )

    def _read(self, outputs, key):
        return pd.read_csv(outputs["tables"][key], sep="\t")

    def test_each_metric_uses_one_cohort_for_every_tool(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            outputs = self._write(Path(temporary_directory))
            summary = self._read(outputs, "cross_dataset_predictor_summary_matched").set_index("tool_id")
            per_dataset = self._read(outputs, "matched_metrics_per_experiment")

        self.assertEqual(summary.loc["a", "spearman_comparisons"], 1)
        self.assertEqual(summary.loc["b", "spearman_comparisons"], 1)
        self.assertEqual(summary.loc["a", "aps_comparisons"], 2)
        self.assertEqual(summary.loc["b", "aps_comparisons"], 2)
        self.assertAlmostEqual(
            summary.loc["a", "spearman_mean"],
            _expected_effect_spearman(self.joined_frames[0], "score_a"),
        )
        # The undefined correlation is kept as missing, never replaced by 0.
        ds2_b = per_dataset[(per_dataset["dataset_id"] == "ds2") & (per_dataset["tool_id"] == "b")]
        self.assertTrue(pd.isna(ds2_b["spearman"].iloc[0]))

    def test_roster_records_configured_included_and_excluded_tools(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            outputs = self._write(Path(temporary_directory), tool_ids=["a", "b", "c"])
            summary = self._read(outputs, "cross_dataset_predictor_summary_matched")
            roster = self._read(outputs, "matched_tool_roster").set_index("tool_id")

        self.assertEqual(set(summary["configured_tool_ids"]), {"a,b,c"})
        self.assertEqual(set(summary["matched_tool_ids"]), {"a,b"})
        self.assertEqual(set(summary["excluded_tool_ids"]), {"c"})
        self.assertEqual(roster.loc[["a", "b"], "included"].tolist(), [True, True])
        self.assertFalse(roster.loc["c", "included"])
        self.assertIn("no metric rows", roster.loc["c", "exclusion_reason"])

    def test_references_are_scored_on_the_common_population(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            outputs = self._write(Path(temporary_directory), reference_baselines=True, reference_seed=0)
            summary = self._read(outputs, "cross_dataset_predictor_summary_matched").set_index("tool_id")
            per_dataset = self._read(outputs, "matched_metrics_per_experiment")
            roster = self._read(outputs, "matched_tool_roster").set_index("tool_id")

        reference_ids = ["reference_expression_high", "reference_expression_low", "reference_random"]
        self.assertTrue(set(reference_ids) <= set(summary.index))
        self.assertTrue(summary.loc[reference_ids, "is_reference"].all())
        self.assertFalse(summary.loc[["a", "b"], "is_reference"].any())
        self.assertEqual(set(summary["matched_tool_ids"]), {"a,b"})
        self.assertEqual(roster.loc["reference_random", "tool_type"], "reference")
        ds1 = per_dataset[per_dataset["dataset_id"] == "ds1"].set_index("tool_id")
        self.assertEqual(set(ds1["common_rows"]), {6})
        common = self.joined_frames[0]
        y_true = ((-common["logFC"] > 1.0) & (common["FDR"] < 0.05)).astype(int)
        self.assertAlmostEqual(
            ds1.loc["reference_expression_high", "aps"],
            average_precision_score(y_true, common["control_mean_normalized_count"]),
        )

    def test_reference_with_missing_expression_on_common_genes_is_unavailable(self):
        self.joined_frames[0].loc[2, "control_mean_normalized_count"] = np.nan
        with tempfile.TemporaryDirectory() as temporary_directory:
            outputs = self._write(Path(temporary_directory), reference_baselines=True, reference_seed=0)
            per_dataset = self._read(outputs, "matched_metrics_per_experiment")

        ds1 = per_dataset[per_dataset["dataset_id"] == "ds1"].set_index("tool_id")
        for tool_id in ("reference_expression_high", "reference_expression_low"):
            with self.subTest(tool_id=tool_id):
                self.assertEqual(
                    ds1.loc[tool_id, "metrics_na_reason"],
                    "reference_missing_on_1_of_6_common_genes",
                )
                self.assertTrue(pd.isna(ds1.loc[tool_id, "aps"]))
                self.assertTrue(pd.isna(ds1.loc[tool_id, "spearman"]))
        self.assertFalse(pd.isna(ds1.loc["reference_random", "aps"]))


def _reference_row(dataset_id, tool_id, matched_tool_id, **metrics):
    return {
        "dataset_id": dataset_id,
        "tool_id": tool_id,
        "is_reference": True,
        "population": f"matched_to_{matched_tool_id}" if matched_tool_id else "all_usable_genes",
        "matched_tool_id": matched_tool_id,
        "coverage": 1.0,
        "positive_coverage": 1.0,
        "aps": metrics.get("aps", 0.5),
        "pr_auc": metrics.get("pr_auc", 0.5),
        "auroc": metrics.get("auroc", 0.5),
        "spearman": metrics.get("spearman", 0.1),
        "pearson": metrics.get("pearson", 0.1),
    }


class PairedReferenceSummaryTests(unittest.TestCase):
    """Tool-vs-reference comparisons use one per-metric cohort for each pair."""

    def setUp(self):
        self.metric_rows = [
            {**_metric_row("ds1", "b", 0.8), "spearman": 0.6, "pearson": 0.5},
            # Tool b has constant scores in ds2, so its correlations are undefined there.
            {**_metric_row("ds2", "b", 0.6), "spearman": float("nan"), "pearson": float("nan")},
        ]
        self.reference_rows = [
            _reference_row("ds1", "reference_random", "b", aps=0.3, spearman=0.1),
            _reference_row("ds2", "reference_random", "b", aps=0.2, spearman=0.4),
            _reference_row("ds1", "reference_random", "", aps=0.25, spearman=0.0),
            _reference_row("ds2", "reference_random", "", aps=0.25, spearman=0.0),
        ]

    def _write(self, root):
        return write_reference_baseline_tables(
            self.reference_rows,
            root / "per_experiment",
            root / "combined",
            metric_rows=self.metric_rows,
        )

    def test_paired_summary_uses_metric_specific_common_cohort(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            outputs = self._write(Path(temporary_directory))
            paired = pd.read_csv(outputs["paired_summary"], sep="\t")

        paired = paired.set_index(["tool_id", "reference_tool_id", "metric"])
        spearman = paired.loc[("b", "reference_random", "spearman")]
        self.assertEqual(spearman["comparisons"], 1)
        self.assertEqual(spearman["population"], "matched_to_b")
        self.assertEqual(spearman["summary_type"], "paired_common_cohort")
        self.assertAlmostEqual(spearman["tool_mean"], 0.6)
        self.assertAlmostEqual(spearman["reference_mean"], 0.1)
        self.assertAlmostEqual(spearman["mean_difference_reference_minus_tool"], 0.1 - 0.6)
        aps = paired.loc[("b", "reference_random", "aps")]
        self.assertEqual(aps["comparisons"], 2)
        self.assertAlmostEqual(aps["mean_difference_reference_minus_tool"], ((0.3 - 0.8) + (0.2 - 0.6)) / 2.0)

    def test_per_reference_aggregate_is_labelled_descriptive(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            outputs = self._write(Path(temporary_directory))
            summary = pd.read_csv(outputs["summary"], sep="\t", keep_default_na=False)

        self.assertEqual(set(summary["summary_type"]), {"descriptive_unpaired"})
        matched_b = summary[summary["population"] == "matched_to_b"].iloc[0]
        self.assertEqual(int(matched_b["spearman_count"]), 2)


if __name__ == "__main__":
    unittest.main()
