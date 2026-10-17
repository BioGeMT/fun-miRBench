import tempfile
import unittest
from pathlib import Path

import pandas as pd

from funmirbench import DatasetMeta
from funmirbench.evaluate import evaluate_joined_dataframe
from funmirbench.evaluate_common import _prepare_scored_frame
from funmirbench.join import load_experiment_table


COVARIATES = [
    "baseMean",
    "control_mean_normalized_count",
    "treated_mean_normalized_count",
    "lfcSE",
    "stat",
    "analysis_status",
]


def _meta(root, data_path):
    return DatasetMeta(
        id="dataset",
        miRNA="hsa-miR-1",
        cell_line="",
        tissue="",
        perturbation="Overexpression",
        organism="Homo sapiens",
        geo_accession="",
        data_path=data_path,
        root=root,
    )


def _de_with_covariates():
    return pd.DataFrame({
        "gene_id": ["ENSG000001", "ENSG000002", "ENSG000003", "ENSG000004"],
        "baseMean": [100.0, 200.0, 300.0, 400.0],
        "logFC": [-2.0, -1.5, 0.1, 0.2],
        "control_mean_normalized_count": [110.0, 210.0, 310.0, 410.0],
        "treated_mean_normalized_count": [27.0, 74.0, 330.0, 470.0],
        "lfcSE": [0.2, 0.3, 0.4, 0.5],
        "stat": [-10.0, -5.0, 0.25, 0.4],
        "PValue": [1e-5, 1e-4, 0.8, 0.7],
        "FDR": [1e-4, 1e-3, 0.9, 0.8],
        "analysis_status": ["tested", "tested", "tested", "tested"],
        "score_shadow": [9.0, 9.0, 9.0, 9.0],
    })


class DeCovariateTests(unittest.TestCase):
    def test_join_keeps_available_de_covariates(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            _de_with_covariates().to_csv(root / "de.tsv", sep="\t", index=False)

            out = load_experiment_table(_meta(root, "de.tsv"))

        for column in COVARIATES:
            with self.subTest(column=column):
                self.assertIn(column, out.columns)
        self.assertEqual(out["lfcSE"].tolist(), [0.2, 0.3, 0.4, 0.5])
        self.assertNotIn("score_shadow", out.columns)

    def test_join_tolerates_missing_covariates(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            _de_with_covariates()[["gene_id", "logFC", "FDR"]].to_csv(
                root / "de.tsv", sep="\t", index=False
            )

            out = load_experiment_table(_meta(root, "de.tsv"))

        self.assertEqual(
            list(out.columns),
            ["dataset_id", "mirna", "perturbation", "gene_id", "logFC", "FDR", "plot_FDR", "benchmark_FDR"],
        )

    def test_scored_frame_keeps_covariates_for_diagnostics(self):
        joined = _de_with_covariates().drop(columns=["score_shadow"])
        joined["score_tool"] = [0.9, 0.8, 0.1, 0.2]

        scored, _ = _prepare_scored_frame(
            joined,
            score_col="score_tool",
            fdr_threshold=0.05,
            effect_threshold=1.0,
            perturbation="Overexpression",
        )

        for column in COVARIATES:
            with self.subTest(column=column):
                self.assertIn(column, scored.columns)

    def test_covariates_are_never_scored_as_predictors(self):
        joined = _de_with_covariates().drop(columns=["score_shadow"])
        joined["score_tool"] = [0.9, 0.8, 0.1, 0.2]

        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            result = evaluate_joined_dataframe(
                joined,
                plots_dir=root / "plots",
                reports_dir=root / "reports",
                fdr_threshold=0.05,
                effect_threshold=1.0,
                predictor_top_fraction=0.1,
                dataset_id="dataset",
                perturbation="Overexpression",
            )

        self.assertEqual([row["tool_id"] for row in result["metric_rows"]], ["tool"])


if __name__ == "__main__":
    unittest.main()
