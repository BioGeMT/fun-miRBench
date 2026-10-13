import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from funmirbench.benchmark import run_benchmark


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

    def test_run_directory_builder_receives_benchmark_context(self):
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
                patch("funmirbench.benchmark._finalize_run_bundle") as finalize,
            ):
                out_dir = run_benchmark(config_path)

            self.assertEqual(out_dir, root / "results" / "run-id")
            build_kwargs = build_name.call_args.kwargs
            self.assertEqual([item.id for item in build_kwargs["experiments"]], ["dataset"])
            self.assertEqual(build_kwargs["tool_ids"], ["tool"])
            self.assertEqual(build_kwargs["eval_cfg"], {"protein_coding_only": False})
            self.assertEqual(build_kwargs["tags"], ["smoke"])
            self.assertEqual(finalize.call_args.kwargs["skipped_datasets"][0]["dataset_id"], "dataset")


if __name__ == "__main__":
    unittest.main()
