import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

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


if __name__ == "__main__":
    unittest.main()
