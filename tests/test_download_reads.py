import hashlib
import importlib.util
import pathlib
import tempfile
import unittest
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import requests

SPEC = importlib.util.spec_from_file_location("download_reads", pathlib.Path(__file__).resolve().parents[1] /
                                            "pipelines/processing/src/download_reads.py")
downloader = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(downloader)


class DownloadReadsTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = pathlib.Path(tmp.name)
        repo = patch.object(downloader, "REPO", self.root)
        repo.start()
        self.addCleanup(repo.stop)
        sleep = patch.object(downloader.time, "sleep")
        self.sleep = sleep.start()
        self.addCleanup(sleep.stop)
        self.payload = b"verified FASTQ contents"
        self.row = SimpleNamespace(path="reads/run.fastq.gz", source_url="https://example.test/run.fastq.gz",
                                   expected_bytes=len(self.payload),
                                   expected_md5=hashlib.md5(self.payload).hexdigest())
        self.dest = self.root / self.row.path
        self.part = self.dest.with_name(self.dest.name + ".part")

    def response(self, chunks=None):
        response = MagicMock()
        response.__enter__.return_value = response
        response.iter_content.return_value = [self.payload] if chunks is None else chunks
        return response

    def test_verified_existing_file_does_not_download(self):
        self.dest.parent.mkdir()
        self.dest.write_bytes(self.payload)
        with patch.object(downloader.requests, "get") as get:
            self.assertTrue(downloader.fetch(self.row).startswith("ok"))
        get.assert_not_called()

    def test_proxy_failure_retries_then_promotes_verified_file(self):
        with patch.object(downloader.requests, "get", side_effect=[
                requests.exceptions.ProxyError("temporary proxy 503"), self.response()]) as get:
            self.assertTrue(downloader.fetch(self.row).startswith("fetched"))
        self.assertEqual(get.call_count, 2)
        self.sleep.assert_called_once_with(2)
        self.assertEqual(self.dest.read_bytes(), self.payload)
        self.assertFalse(self.part.exists())

    def test_bad_checksum_retries_even_when_size_matches(self):
        with patch.object(downloader.requests, "get", side_effect=[
                self.response([b"x" * len(self.payload)]), self.response()]):
            downloader.fetch(self.row)
        self.assertEqual(self.dest.read_bytes(), self.payload)
        self.assertFalse(self.part.exists())

    def test_interrupted_stream_restarts_file_without_appending_partial_bytes(self):
        def broken_stream():
            yield b"partial transfer"
            raise requests.exceptions.ChunkedEncodingError("connection lost")
        first = self.response(broken_stream())
        with patch.object(downloader.requests, "get", side_effect=[first, self.response()]):
            downloader.fetch(self.row)
        self.assertEqual(self.dest.read_bytes(), self.payload)
        self.assertFalse(self.part.exists())

    def test_retry_exhaustion_cleans_partial_and_preserves_existing_destination(self):
        self.dest.parent.mkdir()
        self.dest.write_bytes(b"old invalid file")
        with patch.object(downloader.requests, "get", side_effect=lambda *a, **kw: self.response([b"bad"])) as get:
            with self.assertRaisesRegex(downloader.DownloadError, "after 3 attempt"):
                downloader.fetch(self.row, attempts=3)
        self.assertEqual(get.call_count, 3)
        self.assertEqual(self.dest.read_bytes(), b"old invalid file")
        self.assertFalse(self.part.exists())

    def test_permanent_http_failure_is_not_retried(self):
        response = self.response()
        http_response = requests.Response()
        http_response.status_code = 404
        response.raise_for_status.side_effect = requests.HTTPError(response=http_response)
        with patch.object(downloader.requests, "get", return_value=response) as get:
            with self.assertRaisesRegex(downloader.DownloadError, "HTTP 404 after 1 attempt"):
                downloader.fetch(self.row)
        self.assertEqual(get.call_count, 1)
        self.sleep.assert_not_called()
        self.assertFalse(self.part.exists())


if __name__ == "__main__":
    unittest.main()
