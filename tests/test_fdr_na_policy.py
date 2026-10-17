import math
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from funmirbench import DatasetMeta
from funmirbench.de_table import add_fdr_derived_columns
from funmirbench.evaluate import evaluate_joined_dataframe
from funmirbench.join import load_experiment_table
from funmirbench.validate_experiments import validate_experiments


def _joined():
    # The last gene has a missing FDR but a reported PValue (e.g. DESeq2 independent
    # filtering). It is strongly down-regulated and top-scored by the tool.
    return pd.DataFrame({
        "gene_id": [f"ENSG00000000{index:03d}" for index in range(1, 9)],
        "logFC": [-2.0, -1.8, -1.5, 0.1, 0.2, -0.1, 0.0, -2.5],
        "FDR": [0.001, 0.01, 0.02, 0.5, 0.6, 0.7, 0.9, np.nan],
        "PValue": [1e-4, 1e-3, 2e-3, 0.3, 0.4, 0.5, 0.8, 1e-3],
        "perturbation": "Overexpression",
        "score_tool": [0.9, 0.8, 0.7, 0.1, 0.2, 0.3, 0.4, 0.95],
    })


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


class FdrNaPolicyTests(unittest.TestCase):
    def test_default_policy_keeps_legacy_negative_label(self):
        out = add_fdr_derived_columns(_joined())

        self.assertEqual(out["benchmark_FDR"].iloc[-1], 1.0)

    def test_exclude_policy_leaves_benchmark_fdr_missing(self):
        out = add_fdr_derived_columns(_joined(), fdr_na_policy="exclude")

        self.assertTrue(math.isnan(out["benchmark_FDR"].iloc[-1]))
        self.assertEqual(out["plot_FDR"].iloc[-1], 1.0)
        self.assertEqual(out["FDR"].isna().sum(), 1)

    def test_unknown_policy_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "fdr_na_policy"):
            add_fdr_derived_columns(_joined(), fdr_na_policy="positive")

    def test_evaluator_default_counts_missing_fdr_row_as_negative(self):
        row = _evaluate(_joined())["metric_rows"][0]

        self.assertEqual(row["rows_total"], 8)
        self.assertAlmostEqual(row["auroc"], 12.0 / 15.0)
        self.assertEqual(row["rows_excluded_fdr_na"], 0)

    def test_evaluator_exclude_policy_drops_row_from_both_classes(self):
        row = _evaluate(_joined(), fdr_na_policy="exclude")["metric_rows"][0]

        self.assertEqual(row["rows_total"], 7)
        self.assertEqual(row["positives_total"], 3)
        self.assertAlmostEqual(row["auroc"], 1.0)
        self.assertEqual(row["rows_excluded_fdr_na"], 1)
        self.assertEqual(row["rows_excluded_fdr_na_scored"], 1)

    def test_join_applies_policy_when_loading_de_table(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            _joined().drop(columns=["perturbation", "score_tool"]).to_csv(
                root / "de.tsv", sep="\t", index=False
            )
            meta = DatasetMeta(
                id="dataset",
                miRNA="hsa-miR-1",
                cell_line="",
                tissue="",
                perturbation="Overexpression",
                organism="Homo sapiens",
                geo_accession="",
                data_path="de.tsv",
                root=root,
            )

            legacy = load_experiment_table(meta)
            excluded = load_experiment_table(meta, fdr_na_policy="exclude")

        self.assertEqual(legacy["benchmark_FDR"].iloc[-1], 1.0)
        self.assertTrue(math.isnan(excluded["benchmark_FDR"].iloc[-1]))

    def test_validation_applies_policy_to_ground_truth_classes(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            # Without the missing-FDR row there are no negatives at all.
            pd.DataFrame({
                "gene_id": ["ENSG000001", "ENSG000002", "ENSG000003"],
                "logFC": [-2.0, -2.0, 0.0],
                "control_mean_normalized_count": [10.0, 20.0, 30.0],
                "PValue": [1e-4, 1e-4, 0.01],
                "FDR": [0.001, 0.001, np.nan],
            }).to_csv(root / "de.tsv", sep="\t", index=False)
            pd.DataFrame({
                "id": ["dataset"],
                "mirna_name": ["hsa-miR-1"],
                "organism": ["Homo sapiens"],
                "experiment_type": ["Overexpression"],
                "de_table_path": ["de.tsv"],
            }).to_csv(root / "metadata.tsv", sep="\t", index=False)

            legacy = validate_experiments(root / "metadata.tsv", root=root)
            excluded = validate_experiments(
                root / "metadata.tsv",
                root=root,
                fdr_na_policy="exclude",
            )

        self.assertEqual(legacy.issues, [])
        self.assertEqual([issue.check for issue in excluded.issues], ["ground_truth_classes"])
        self.assertIn("no negative genes", excluded.issues[0].message)


if __name__ == "__main__":
    unittest.main()
