import json
import math
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from funmirbench.benchmark import run_benchmark


GENE_IDS = [f"ENSG00000000{index:03d}" for index in range(1, 13)]
REFERENCE_IDS = {
    "reference_random",
    "reference_expression_high",
    "reference_expression_low",
}


def _write_fixture(root, *, evaluation_lines=(), unscored_tool=False):
    pd.DataFrame({
        "gene_id": GENE_IDS,
        "logFC": [-2.0, -1.8, -1.5, -2.2, 0.1, 0.2, -0.1, 0.0, 0.3, -0.2, 0.4, -2.5],
        "control_mean_normalized_count": [100.0 * (index + 1) for index in range(12)],
        "PValue": [1e-4, 1e-3, 2e-3, 1e-4, 0.3, 0.4, 0.5, 0.8, 0.6, 0.2, 0.3, 1e-3],
        # The last gene has a missing FDR but a reported PValue.
        "FDR": [0.001, 0.01, 0.02, 0.001, 0.5, 0.6, 0.7, 0.9, 0.8, 0.4, 0.5, None],
    }).to_csv(root / "de_pos.tsv", sep="\t", index=False)
    pd.DataFrame({
        "gene_id": GENE_IDS,
        "logFC": [-0.9, -0.8, -0.5, -0.2, 0.1, 0.2, -0.1, 0.0, 0.3, -0.2, 0.4, 0.5],
        "control_mean_normalized_count": [100.0 * (index + 1) for index in range(12)],
        "PValue": [0.01] * 12,
        "FDR": [0.04] * 12,
    }).to_csv(root / "de_zero.tsv", sep="\t", index=False)
    predictor_rows = []
    for index, gene_id in enumerate(GENE_IDS):
        predictor_rows.append(("tool_a", gene_id, 1.0 - index / 20.0))
        if index % 2 == 0:
            predictor_rows.append(("tool_b", gene_id, float(index)))
    for tool_id in ("tool_a", "tool_b"):
        pd.DataFrame(
            [
                {"Ensembl_ID": gene_id, "miRNA_Name": "hsa-miR-1", "Score": score}
                for row_tool, gene_id, score in predictor_rows
                if row_tool == tool_id
            ]
        ).to_csv(root / f"{tool_id}.tsv", sep="\t", index=False)
    predictions = (
        "tool_id\tofficial_name\tscore_direction\tpredictor_output_path\n"
        "tool_a\tTool A\thigher_is_stronger\ttool_a.tsv\n"
        "tool_b\tTool B\thigher_is_stronger\ttool_b.tsv\n"
    )
    if unscored_tool:
        # Configured, but it only scores a miRNA that no selected experiment uses.
        pd.DataFrame(
            {"Ensembl_ID": GENE_IDS[:3], "miRNA_Name": "hsa-miR-999", "Score": [1.0, 2.0, 3.0]}
        ).to_csv(root / "tool_c.tsv", sep="\t", index=False)
        predictions += "tool_c\tTool C\thigher_is_stronger\ttool_c.tsv\n"
    (root / "predictions.tsv").write_text(predictions, encoding="utf-8")
    (root / "experiments.tsv").write_text(
        "id\tmirna_name\torganism\texperiment_type\tde_table_path\n"
        "ds_pos\thsa-miR-1\tHomo sapiens\tOE\tde_pos.tsv\n"
        "ds_zero\thsa-miR-1\tHomo sapiens\tOverexpression\tde_zero.tsv\n",
        encoding="utf-8",
    )
    config_path = root / "benchmark.yaml"
    config_path.write_text(
        "experiments_tsv: experiments.tsv\n"
        "predictions_tsv: predictions.tsv\n"
        "evaluation:\n"
        "  protein_coding_only: false\n"
        + "".join(f"  {line}\n" for line in evaluation_lines),
        encoding="utf-8",
    )
    return config_path


def _run(config_path, run_name):
    with (
        patch("requests.get", side_effect=AssertionError("unexpected network request")),
        patch("funmirbench.benchmark.build_run_dir_name", return_value=run_name),
    ):
        return run_benchmark(config_path)


def _read_table(out_dir, relative_path):
    return pd.read_csv(out_dir / relative_path, sep="\t")


class BenchmarkOptionTests(unittest.TestCase):
    def test_default_options_keep_legacy_selection_and_add_reference_outputs(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            out_dir = _run(_write_fixture(root), "defaults")
            summary = json.loads((out_dir / "summary.json").read_text(encoding="utf-8"))
            aps = _read_table(out_dir, "tables/per_experiment/aps_per_experiment.tsv")
            references = _read_table(out_dir, "tables/per_experiment/reference_baselines_per_experiment.tsv")
            continuous = _read_table(out_dir, "tables/per_experiment/continuous_metrics_per_experiment.tsv")
            reference_summary_exists = (
                out_dir / "tables/combined/cross_dataset_reference_summary.tsv"
            ).is_file()
            matched_summary_exists = (
                out_dir / "tables/combined/cross_dataset_predictor_summary_matched.tsv"
            ).is_file()
            readme = (out_dir / "README.md").read_text(encoding="utf-8")

        self.assertEqual(summary["dataset_ids"], ["ds_pos"])
        self.assertEqual([item["dataset_id"] for item in summary["skipped_datasets"]], ["ds_zero"])
        self.assertEqual(
            summary["evaluation_options"],
            {
                "fdr_na_policy": "negative",
                "zero_positive_policy": "exclude",
                "reference_baselines": True,
                "reference_seed": 0,
            },
        )
        self.assertFalse(REFERENCE_IDS & set(aps.columns))
        self.assertEqual(set(references["tool_id"]), REFERENCE_IDS)
        self.assertTrue(references["is_reference"].all())
        self.assertTrue(reference_summary_exists)
        self.assertTrue(matched_summary_exists)
        self.assertEqual(set(continuous["tool_id"]), {"tool_a", "tool_b"})
        self.assertEqual(set(continuous.loc[continuous["dataset_id"] == "ds_pos", "rows_total"]), {12})
        self.assertIn("reference", readme.lower())

    def test_opt_in_options_are_plumbed_from_config(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            config_path = _write_fixture(
                root,
                evaluation_lines=(
                    "fdr_na_policy: exclude",
                    "zero_positive_policy: continuous_only",
                    "reference_baselines: false",
                ),
            )
            out_dir = _run(config_path, "opt-in")
            summary = json.loads((out_dir / "summary.json").read_text(encoding="utf-8"))
            spearman = _read_table(out_dir, "tables/per_experiment/spearman_per_experiment.tsv")
            aps = _read_table(out_dir, "tables/per_experiment/aps_per_experiment.tsv")
            continuous = _read_table(out_dir, "tables/per_experiment/continuous_metrics_per_experiment.tsv")
            reference_table_exists = (
                out_dir / "tables/per_experiment/reference_baselines_per_experiment.tsv"
            ).exists()

        self.assertEqual(summary["dataset_ids"], ["ds_pos", "ds_zero"])
        self.assertEqual(summary["skipped_datasets"], [])
        self.assertEqual(summary["evaluation_options"]["fdr_na_policy"], "exclude")
        self.assertEqual(summary["evaluation_options"]["zero_positive_policy"], "continuous_only")
        self.assertFalse(reference_table_exists)
        zero_spearman = spearman.set_index("dataset_id").loc["ds_zero", "tool_a"]
        self.assertTrue(math.isfinite(float(zero_spearman)))
        zero_aps = aps.set_index("dataset_id")
        self.assertTrue("ds_zero" not in zero_aps.index or pd.isna(zero_aps.loc["ds_zero", "tool_a"]))
        self.assertEqual(set(continuous.loc[continuous["dataset_id"] == "ds_pos", "rows_total"]), {11})

    def test_matched_roster_and_reference_populations_are_reported(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            out_dir = _run(_write_fixture(root, unscored_tool=True), "roster")
            roster = _read_table(out_dir, "tables/combined/matched_tool_roster.tsv").set_index("tool_id")
            references = _read_table(out_dir, "tables/per_experiment/reference_baselines_per_experiment.tsv")
            matched = _read_table(out_dir, "tables/combined/cross_dataset_predictor_summary_matched.tsv")
            readme = (out_dir / "README.md").read_text(encoding="utf-8")

        self.assertFalse(roster.loc["tool_c", "included"])
        self.assertIn("no metric rows", roster.loc["tool_c", "exclusion_reason"])
        self.assertEqual(set(matched["configured_tool_ids"]), {"tool_a,tool_b,tool_c"})
        self.assertEqual(set(matched["matched_tool_ids"]), {"tool_a,tool_b"})
        self.assertTrue(REFERENCE_IDS <= set(matched["tool_id"]))
        self.assertEqual(
            set(references["population"]),
            {"all_usable_genes", "matched_to_tool_a", "matched_to_tool_b"},
        )
        matched_b = references[references["population"] == "matched_to_tool_b"]
        self.assertEqual(set(matched_b["rows_scored"]), {6})
        self.assertIn("`tool_c`", readme)
        self.assertIn("no metric rows", readme)
        self.assertIn("all_usable_genes", readme)
        self.assertIn("matched_to_<tool_id>", readme)

    def test_reference_summaries_are_labelled_descriptive_or_paired(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            out_dir = _run(_write_fixture(root), "paired")
            summary = json.loads((out_dir / "summary.json").read_text(encoding="utf-8"))
            descriptive = _read_table(out_dir, "tables/combined/cross_dataset_reference_summary.tsv")
            paired = _read_table(out_dir, "tables/combined/cross_dataset_tool_vs_reference_paired.tsv")
            readme = (out_dir / "README.md").read_text(encoding="utf-8")

        self.assertIn(
            "cross_dataset_tool_vs_reference_paired",
            summary["cross_dataset_outputs"]["tables"],
        )
        self.assertEqual(set(descriptive["summary_type"]), {"descriptive_unpaired"})
        self.assertEqual(set(paired["summary_type"]), {"paired_common_cohort"})
        self.assertEqual(set(paired["tool_id"]), {"tool_a", "tool_b"})
        self.assertEqual(set(paired["reference_tool_id"]), REFERENCE_IDS)
        self.assertIn("mean_difference_reference_minus_tool", paired.columns)
        self.assertIn("descriptive_unpaired", readme)
        self.assertIn("cross_dataset_tool_vs_reference_paired.tsv", readme)


if __name__ == "__main__":
    unittest.main()
