import gzip
import json
import tempfile
import unittest
import zipfile
from pathlib import Path
from unittest.mock import patch

from funmirbench.experiment_store import ensure_zenodo_experiment_cached
from funmirbench.predictor_store import sync_zenodo_predictors
from funmirbench.zenodo_store import compute_md5


class FakeResponse:
    def __init__(self, payload):
        self.payload = payload

    def raise_for_status(self):
        return None

    def iter_content(self, chunk_size):
        for index in range(0, len(self.payload), chunk_size):
            yield self.payload[index:index + chunk_size]


class ZenodoArtifactTests(unittest.TestCase):
    def test_experiment_download_accepts_flat_gzipped_artifact(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            repo = Path(temporary_directory)
            payload = b"gene_id\tlogFC\tFDR\nENSG1\t1.2\t0.01\n"
            compressed = gzip.compress(payload)
            compressed_path = repo / "remote.tsv.gz"
            compressed_path.write_bytes(compressed)
            destination = repo / "data/experiments/processed/21671476/example.tsv"
            registry = {
                "example.tsv.gz": {
                    "checksum": f"md5:{compute_md5(compressed_path)}",
                    "url": "https://example.test/example.tsv.gz",
                }
            }

            with patch("funmirbench.experiment_store.requests.get", return_value=FakeResponse(compressed)):
                cached = ensure_zenodo_experiment_cached(
                    "data/experiments/processed/21671476/example.tsv",
                    repo=repo,
                    registry=registry,
                )

            self.assertEqual(cached, destination)
            self.assertEqual(destination.read_bytes(), payload)
            manifest_path = destination.with_name("example.tsv.funmirbench.json")
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            self.assertEqual(manifest["remote_filename"], "example.tsv.gz")
            self.assertEqual(manifest["remote_checksum"], registry["example.tsv.gz"]["checksum"])

            destination.write_bytes(b"corrupt local cache\n")
            with self.assertRaises(ValueError):
                ensure_zenodo_experiment_cached(
                    "data/experiments/processed/21671476/example.tsv",
                    repo=repo,
                    registry=registry,
                )

    def test_predictor_sync_extracts_nested_gzipped_members(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            repo = Path(temporary_directory)
            archive_path = repo / "predictions.zip"
            payload = b"Ensembl_ID\tGene_Name\tmiRNA_ID\tmiRNA_Name\tScore\nENSG1\tGENE1\tMIMAT1\thsa-miR-1\t0.5\n"
            with zipfile.ZipFile(archive_path, "w") as archive:
                archive.writestr(
                    "published/predictors/tool/tool_standardized.tsv.gz",
                    gzip.compress(payload),
                )
            predictions = {
                "tool": {
                    "predictor_output_path": "data/predictions/tool/tool_standardized.tsv",
                }
            }
            registry = {
                "fun-miRBench_standardized_predictions.zip": {
                    "filename": "fun-miRBench_standardized_predictions.zip",
                    "checksum": "",
                    "url": "https://example.test/predictions.zip",
                }
            }

            with patch("funmirbench.predictor_store._download_archive", return_value=archive_path):
                paths = sync_zenodo_predictors(predictions, repo=repo, registry=registry)

            destination = repo / "data/predictions/tool/tool_standardized.tsv"
            self.assertEqual(paths, [destination])
            self.assertEqual(destination.read_bytes(), payload)


if __name__ == "__main__":
    unittest.main()
