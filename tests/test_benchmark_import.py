import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from matplotlib import pyplot as plt
import pandas as pd

from funmirbench.benchmark import run_benchmark
from funmirbench.evaluate import evaluate_joined_dataframe
from funmirbench.evaluate_common import _add_figure_heading


class BenchmarkImportTests(unittest.TestCase):
    def test_benchmark_module_imports_with_metadata_schema_changes(self):
        import funmirbench.benchmark  # noqa: F401

    def test_missing_experiment_is_local_validation_error(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            (root / "predictor.tsv").write_text("Score\n1\n", encoding="utf-8")
            (root / "predictions.tsv").write_text(
                "tool_id\tofficial_name\tscore_direction\tpredictor_output_path\n"
                "tool\tTool\thigher_is_stronger\tpredictor.tsv\n",
                encoding="utf-8",
            )
            (root / "experiments.tsv").write_text(
                "id\tmirna_name\torganism\texperiment_type\tde_table_path\n"
                "dataset\thsa-miR-1\tHomo sapiens\tOverexpression\tmissing.tsv\n",
                encoding="utf-8",
            )
            config_path = root / "benchmark.yaml"
            config_path.write_text(
                "experiments_tsv: experiments.tsv\n"
                "predictions_tsv: predictions.tsv\n"
                "evaluation:\n"
                "  protein_coding_only: false\n",
                encoding="utf-8",
            )

            with (
                patch("requests.get", side_effect=AssertionError("unexpected network request")),
                self.assertRaisesRegex(ValueError, "funmirbench-download-data"),
            ):
                run_benchmark(config_path)

    def test_run_directory_builder_context_and_collision_suffix(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            (root / "predictor.tsv").write_text("Score\n1\n", encoding="utf-8")
            (root / "predictions.tsv").write_text(
                "tool_id\tofficial_name\tscore_direction\tpredictor_output_path\n"
                "tool\tTool\thigher_is_stronger\tpredictor.tsv\n",
                encoding="utf-8",
            )
            (root / "de.tsv").write_text(
                "gene_id\tlogFC\tcontrol_mean_normalized_count\tFDR\n"
                "ENSG000001\t-2.0\t10.0\t0.01\n"
                "ENSG000002\t0.0\t20.0\t0.5\n",
                encoding="utf-8",
            )
            (root / "experiments.tsv").write_text(
                "id\tmirna_name\torganism\texperiment_type\tde_table_path\n"
                "dataset\thsa-miR-1\tHomo sapiens\tOverexpression\tde.tsv\n",
                encoding="utf-8",
            )
            config_path = root / "benchmark.yaml"
            config_path.write_text(
                "experiments_tsv: experiments.tsv\n"
                "predictions_tsv: predictions.tsv\n"
                "tags: [smoke]\n"
                "evaluation:\n"
                "  protein_coding_only: false\n",
                encoding="utf-8",
            )
            joined = pd.DataFrame({
                "gene_id": ["ENSG000001", "ENSG000002"],
                "logFC": [-2.0, 0.0],
                "FDR": [0.01, 0.5],
                "score_tool": [float("nan"), float("nan")],
            })
            first_run = root / "results" / "run-id"
            second_run = root / "results" / "run-id__r2"
            first_run.mkdir(parents=True)
            second_run.mkdir()
            marker = first_run / "existing.txt"
            marker.write_text("keep", encoding="utf-8")

            with (
                patch(
                    "funmirbench.benchmark.build_run_dir_name",
                    return_value="run-id",
                ) as build_name,
                patch(
                    "funmirbench.benchmark.load_predictor_score_cache",
                    return_value={},
                ),
                patch(
                    "funmirbench.benchmark.build_joined",
                    return_value=(joined, {"tool": str(root / "predictor.tsv")}),
                ),
            ):
                out_dir = run_benchmark(config_path)

            self.assertEqual(out_dir, root / "results" / "run-id__r3")
            self.assertEqual(marker.read_text(encoding="utf-8"), "keep")
            build_kwargs = build_name.call_args.kwargs
            self.assertEqual([item.id for item in build_kwargs["experiments"]], ["dataset"])
            self.assertEqual(build_kwargs["tool_ids"], ["tool"])
            self.assertEqual(build_kwargs["eval_cfg"], {"protein_coding_only": False})
            self.assertEqual(build_kwargs["tags"], ["smoke"])

            summary = json.loads((out_dir / "summary.json").read_text(encoding="utf-8"))
            self.assertEqual(summary["tags"], ["smoke"])
            self.assertEqual(summary["dataset_ids"], [])
            self.assertEqual(summary["skipped_datasets"][0]["dataset_id"], "dataset")
            self.assertTrue((out_dir / "README.md").is_file())
            self.assertTrue((out_dir / "REPORT.pdf").is_file())

    def test_dataset_is_skipped_when_no_predictor_rows_are_evaluable(self):
        joined = pd.DataFrame({
            "gene_id": ["ENSG000001", "ENSG000002"],
            "logFC": [-2.0, 0.0],
            "FDR": [0.01, 0.5],
            "score_left": [0.9, float("nan")],
            "score_right": [float("nan"), 0.8],
        })

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

        self.assertEqual(result["metric_rows"], [])
        self.assertEqual(len(result["skipped_tool_rows"]), 2)
        self.assertIn("No selected predictor", result["dataset_skip_reason"])

    def test_publication_figure_heading_remains_disabled(self):
        fig = plt.figure()
        try:
            _add_figure_heading(fig, title="Title", subtitle="Subtitle")
            self.assertEqual(fig.texts, [])
        finally:
            plt.close(fig)


if __name__ == "__main__":
    unittest.main()
