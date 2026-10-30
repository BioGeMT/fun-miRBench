import gzip
import hashlib
import json
import os
import tempfile
import unittest
import zipfile
from pathlib import Path
from unittest.mock import patch

import requests

from funmirbench.experiment_store import (
    ensure_zenodo_experiment_cached,
    sync_zenodo_experiments,
)
from funmirbench.benchmark_data_store import sync_family_expanded_targetscan, sync_benchmark_data
from funmirbench.predictor_store import _download_archive, sync_zenodo_predictors
from funmirbench.zenodo_store import compute_md5


class FakeResponse:
    def __init__(self, payload):
        self.payload = payload

    def raise_for_status(self):
        return None

    def iter_content(self, chunk_size):
        for index in range(0, len(self.payload), chunk_size):
            yield self.payload[index:index + chunk_size]


class FailingResponse(FakeResponse):
    def iter_content(self, chunk_size):
        yield b"partial archive"
        raise requests.ConnectionError("download interrupted")


class ZenodoArtifactTests(unittest.TestCase):
    def test_experiment_sync_leaves_pipeline_file_untouched(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            repo = Path(temporary_directory)
            destination = repo / "data/experiments/processed/example.tsv"
            destination.parent.mkdir(parents=True)
            destination.write_bytes(b"pipeline output\n")
            registry = {
                "example.tsv": {
                    "checksum": "md5:00000000000000000000000000000000",
                    "url": "https://example.test/example.tsv",
                }
            }

            with patch("funmirbench.experiment_store.requests.get") as request:
                paths = sync_zenodo_experiments(
                    ["data/experiments/processed/example.tsv"],
                    repo=repo,
                    registry=registry,
                    force=True,
                )

            self.assertEqual(paths, [])
            self.assertEqual(destination.read_bytes(), b"pipeline output\n")
            request.assert_not_called()

    def test_experiment_sync_validates_existing_cache(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            repo = Path(temporary_directory)
            payload = b"gene_id\tlogFC\tFDR\nENSG1\t1.2\t0.01\n"
            destination = repo / "data/experiments/processed/21671476/example.tsv"
            destination.parent.mkdir(parents=True)
            destination.write_bytes(payload)
            registry = {
                "example.tsv": {
                    "checksum": f"md5:{compute_md5(destination)}",
                    "url": "https://example.test/example.tsv",
                }
            }

            with patch("funmirbench.experiment_store.requests.get") as request:
                paths = sync_zenodo_experiments(
                    ["data/experiments/processed/21671476/example.tsv"],
                    repo=repo,
                    registry=registry,
                )

            self.assertEqual(paths, [destination])
            request.assert_not_called()

            destination.write_bytes(b"corrupt local cache\n")
            with self.assertRaisesRegex(ValueError, "failed checksum verification"):
                sync_zenodo_experiments(
                    ["data/experiments/processed/21671476/example.tsv"],
                    repo=repo,
                    registry=registry,
                )

    def test_compressed_experiment_cache_requires_manifest(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            repo = Path(temporary_directory)
            payload = b"gene_id\tlogFC\tFDR\nENSG1\t1.2\t0.01\n"
            compressed = gzip.compress(payload)
            compressed_path = repo / "remote.tsv.gz"
            compressed_path.write_bytes(compressed)
            destination = repo / "data/experiments/processed/21671476/example.tsv"
            destination.parent.mkdir(parents=True)
            destination.write_bytes(payload)
            registry = {
                "example.tsv.gz": {
                    "checksum": f"md5:{compute_md5(compressed_path)}",
                    "url": "https://example.test/example.tsv.gz",
                }
            }

            with self.assertRaisesRegex(ValueError, "Rerun with --force"):
                sync_zenodo_experiments(
                    ["data/experiments/processed/21671476/example.tsv"],
                    repo=repo,
                    registry=registry,
                )

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
                    "published/predictors/targetscan/targetscan_standardized.tsv.gz",
                    gzip.compress(payload),
                )
            predictions = {
                "targetscan": {
                    "predictor_output_path": "data/predictions/targetscan/targetscan_standardized.tsv",
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

            destination = repo / "data/predictions/targetscan/targetscan_standardized.tsv"
            self.assertEqual(paths, [destination])
            self.assertEqual(destination.read_bytes(), payload)

    def test_missing_user_predictor_does_not_download_archive(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            repo = Path(temporary_directory)
            published_destination = (
                repo / "data/predictions/targetscan/targetscan_standardized.tsv"
            )
            published_destination.parent.mkdir(parents=True)
            published_destination.write_bytes(b"published predictor\n")
            predictions = {
                "targetscan": {
                    "predictor_output_path": (
                        "data/predictions/targetscan/targetscan_standardized.tsv"
                    ),
                },
                "custom": {
                    "predictor_output_path": "data/predictions/custom/missing.tsv",
                },
            }

            with patch(
                "funmirbench.predictor_store.fetch_zenodo_file_registry"
            ) as fetch_registry, patch(
                "funmirbench.predictor_store._download_archive"
            ) as download_archive:
                paths = sync_zenodo_predictors(predictions, repo=repo)

            self.assertEqual(paths, [published_destination])
            fetch_registry.assert_not_called()
            download_archive.assert_not_called()

    def test_predictor_sync_skips_user_member_and_supports_expanded_paths(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            repo = Path(temporary_directory)
            user_dir = repo / "user-predictors"
            user_dir.mkdir()
            user_destination = user_dir / "custom.tsv"
            user_destination.write_bytes(b"user predictor\n")
            published_destination = (
                repo.parent
                / f"{repo.name}-external"
                / "targetscan_standardized.tsv"
            )
            archive_path = repo / "predictions.zip"
            payload = b"published predictor\n"
            with zipfile.ZipFile(archive_path, "w") as archive:
                archive.writestr("published/targetscan_standardized.tsv", payload)
            predictions = {
                "custom": {
                    "predictor_output_path": "$FUNMIRBENCH_USER_PREDICTORS/custom.tsv",
                },
                "targetscan": {
                    "predictor_output_path": str(published_destination),
                },
            }
            registry = {
                "fun-miRBench_standardized_predictions.zip": {
                    "filename": "fun-miRBench_standardized_predictions.zip",
                    "checksum": "",
                    "url": "https://example.test/predictions.zip",
                }
            }

            try:
                with patch.dict(
                    os.environ,
                    {"FUNMIRBENCH_USER_PREDICTORS": str(user_dir)},
                ), patch(
                    "funmirbench.predictor_store._download_archive",
                    return_value=archive_path,
                ):
                    sync_zenodo_predictors(
                        predictions,
                        repo=repo,
                        registry=registry,
                        force=True,
                    )

                self.assertEqual(user_destination.read_bytes(), b"user predictor\n")
                self.assertEqual(published_destination.read_bytes(), payload)
                self.assertFalse((repo / "$FUNMIRBENCH_USER_PREDICTORS").exists())
            finally:
                published_destination.unlink(missing_ok=True)
                published_destination.parent.rmdir()

    def test_failed_predictor_download_removes_temporary_archive(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            meta = {
                "filename": "predictions.zip",
                "size": 100,
                "checksum": "",
                "url": "https://example.test/predictions.zip",
            }

            with patch(
                "funmirbench.predictor_store.requests.get",
                return_value=FailingResponse(b""),
            ), self.assertRaisesRegex(requests.ConnectionError, "interrupted"):
                _download_archive(meta, directory=directory, timeout=10)

            self.assertEqual(
                list(directory.glob(".funmirbench_predictors.*.zip")),
                [],
            )

    def test_benchmark_data_is_unpacked_once_and_checked(self):
        with tempfile.TemporaryDirectory() as tmp:
            repo = Path(tmp)
            archive = repo / "benchmark_data.zip"
            with zipfile.ZipFile(archive, "w") as handle:
                handle.writestr("benchmark_data/comparisons.tsv", "experiment_id\n")
                handle.writestr("benchmark_data/genes/E1.tsv.gz", "table")
            payload = archive.read_bytes()
            archive.unlink()
            meta = {"filename": "benchmark_data.zip", "url": "https://example.org/benchmark_data.zip", "size": len(payload)}
            registry = {"benchmark_data.zip": {**meta, "checksum": "md5:" + hashlib.md5(payload).hexdigest()}}
            with patch("funmirbench.predictor_store.requests.get", return_value=FakeResponse(payload)) as get:
                final = sync_benchmark_data(repo=repo, registry=registry)
                self.assertEqual(sync_benchmark_data(repo=repo, registry=registry), final)
            self.assertEqual(get.call_count, 1)
            self.assertEqual(final, (repo / "data/benchmark_data").resolve())
            self.assertEqual((final / "genes/E1.tsv.gz").read_text(), "table")
            self.assertEqual([p.name for p in final.parent.iterdir()], ["benchmark_data"])
            # An existing folder, such as tables built from reads, is only replaced on request.
            (final / "own.tsv").write_text("x")
            with patch("funmirbench.predictor_store.requests.get", return_value=FakeResponse(payload)):
                sync_benchmark_data(repo=repo, registry=registry)
                self.assertTrue((final / "own.tsv").exists())
                sync_benchmark_data(repo=repo, registry=registry, force=True)
            self.assertFalse((final / "own.tsv").exists())

    def test_family_expanded_targetscan_is_decompressed_to_its_registered_path(self):
        with tempfile.TemporaryDirectory() as tmp:
            repo = Path(tmp)
            (repo / "metadata").mkdir()
            (repo / "metadata/predictions_info.tsv").write_text(
                "tool_id\tpredictor_output_path\ntargetscan_family_expanded\tdata/predictions/t/expanded.tsv\n")
            payload = gzip.compress(b"scores")
            name = "targetscan_family_expanded_standardized.tsv.gz"
            registry = {name: {"filename": name, "url": "https://example.org/x", "size": len(payload),
                               "checksum": "md5:" + hashlib.md5(payload).hexdigest()}}
            with patch("funmirbench.predictor_store.requests.get", return_value=FakeResponse(payload)):
                path = sync_family_expanded_targetscan(repo=repo, registry=registry)
            self.assertEqual(path.read_bytes(), b"scores")
            self.assertEqual([p.name for p in path.parent.iterdir()], ["expanded.tsv"])


if __name__ == "__main__":
    unittest.main()
