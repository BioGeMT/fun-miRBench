import gzip
import hashlib
import io
import json
import pathlib
import tempfile
import unittest
import zipfile
from unittest.mock import patch

import numpy as np
import pandas as pd
import requests
import yaml

from funmirbench.benchmark import load_predictions
from funmirbench.response_vectors import alignment, build_vectors, control_expression_tenths
from funmirbench.vector_engine import PredictionPanel, batch_alignment, load_matrix, response_view
from funmirbench.vector_sources import download, registry, unpack_tables
from funmirbench.vector_explore import replot, run_exploration


class RankEngineTests(unittest.TestCase):
    def test_batch_matches_reference_with_ties_missing_scores_and_response_masks(self):
        rng = np.random.default_rng(24)
        for n in (0, 1, 2, 3, 19, 40, 300, 1000):
            with self.subTest(n=n):
                genes = pd.Index([f"g{i}" for i in range(n)])
                response = rng.integers(-3, 4, n).astype(float)
                response[rng.random(n) < .1] = np.nan
                expression = rng.integers(0, 4, n).astype(float)
                expression[rng.random(n) < .1] = np.nan
                table = pd.DataFrame({"eligible": rng.random(n) > .1, "oriented_logFC_shrunk": response,
                                      "control_mean_normalized_count": expression}, index=genes)
                table["expression_tenth"] = control_expression_tenths(table.control_mean_normalized_count)
                matrix = rng.integers(-2, 3, (7, n)).astype(float)
                matrix[rng.random(matrix.shape) < .3] = np.nan
                matrix[0] = np.nan
                matrix[1] = 1.
                if n == 0:
                    plain, within, *_ = batch_alignment(PredictionPanel(matrix), response_view(table))
                    self.assertTrue(np.isnan(plain).all())
                    self.assertTrue(np.isnan(within).all())
                    continue
                plain, within, support, counts, _ = batch_alignment(PredictionPanel(matrix), response_view(table))
                for i in range(len(matrix)):
                    reference = alignment(build_vectors(table, pd.Series(matrix[i], index=genes)), table.expression_tenth)
                    np.testing.assert_allclose(plain[i], reference["spearman"], atol=1e-14, equal_nan=True)
                    np.testing.assert_allclose(within[i], reference["spearman_within_tenths"], atol=1e-14, equal_nan=True)
                    self.assertEqual(support[i], reference["within_tenths_bins"])
                    self.assertEqual(counts[i], reference["n_within_tenths_genes"])

    def test_compact_loader_matches_reference_both_directions_zip_gzip_and_cache_invalidation(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            payload = ("Ensembl_ID\tmiRNA_Name\tScore\n"
                       "g1.1\tm1\t-5\ng1.2\tm1\t2\ng2\tm1\t0\ng3\tm2\t7\n"
                       "g1\tm1\tinf\ng1\tm1\t-inf\ng2\tm2\tbad\noutside\tm1\t99\n")
            path = root / "scores.tsv.gz"
            path.write_bytes(gzip.compress(payload.encode()))
            archive = root / "scores.zip"
            with zipfile.ZipFile(archive, "w") as z:
                z.writestr("nested/scores.tsv.gz", path.read_bytes())
            genes, mirnas = pd.Index(["g1", "g2", "g3"]), pd.Index(["m1", "m2", "missing"])
            for direction in ("higher_is_stronger", "lower_is_stronger"):
                expected = load_predictions(path, direction, set(mirnas), set(genes))
                for source in ({"path": str(path)}, {"path": str(archive), "member": "nested/scores.tsv.gz"}):
                    matrix = load_matrix(source, direction, mirnas, genes, root)
                    for i, name in enumerate(mirnas):
                        reference = expected.get(name, pd.Series(dtype=float)).reindex(genes).to_numpy(float)
                        np.testing.assert_allclose(matrix[i], reference, equal_nan=True)
                    with patch("funmirbench.vector_engine.prediction_stream", side_effect=AssertionError("reparsed")):
                        np.testing.assert_allclose(load_matrix(source, direction, mirnas, genes, root), matrix, equal_nan=True)
            path.write_bytes(gzip.compress(payload.replace("m2\t7", "m2\t8").encode()))
            updated = load_matrix({"path": str(path)}, "higher_is_stronger", mirnas, genes, root)
            self.assertEqual(updated[1, 2], 8.)

    def test_rank_cache_respects_memory_limit_and_does_not_change_results(self):
        matrix = np.arange(100., dtype=float).reshape(5, 20)
        panel = PredictionPanel(matrix, max_cache_bytes=500)
        for cutoff in range(3, 21):
            mask = np.arange(20) < cutoff
            result = panel.ranked(mask)
            self.assertLessEqual(panel.cache_bytes, 500)
            np.testing.assert_allclose(result, PredictionPanel(matrix).ranked(mask))


class Response:
    def __init__(self, payload, status=200, headers=None, broken=False):
        self.payload, self.status_code, self.headers, self.broken = payload, status, headers or {}, broken
    def __enter__(self):
        return self
    def __exit__(self, *args):
        pass
    def raise_for_status(self):
        pass
    def iter_content(self, size):
        yield self.payload
        if self.broken:
            raise requests.ConnectionError("interrupted")


class SourceTests(unittest.TestCase):
    def http_error(self, status):
        response = requests.Response()
        response.status_code = status
        return requests.HTTPError(f"HTTP {status}", response=response)

    def test_metadata_retries_transient_errors_then_caches_success(self):
        payload = {"example.zip": {"checksum": "md5:example"}}
        with tempfile.TemporaryDirectory() as tmp, \
                patch("funmirbench.vector_sources.fetch_zenodo_file_registry",
                      side_effect=[self.http_error(500), requests.Timeout(), payload]) as fetch, \
                patch("funmirbench.vector_sources.time.sleep") as sleep:
            root = pathlib.Path(tmp)
            self.assertEqual(registry("23839843", root), payload)
            self.assertEqual(registry("23839843", root), payload)
            self.assertEqual(fetch.call_count, 3)
            self.assertEqual([c.args[0] for c in sleep.call_args_list], [2, 4])

    def test_metadata_outage_uses_pinned_public_files_with_verified_checksums(self):
        for record, name, checksum in (
                ("23839843", "benchmark_data.zip", "42bd622a8760671ca80c477e7fa81b87"),
                ("21671476", "fun-miRBench_standardized_predictions.zip", "9e619a17668ea780fa3f4de85a7f5476")):
            with self.subTest(record=record), tempfile.TemporaryDirectory() as tmp, \
                    patch("funmirbench.vector_sources.fetch_zenodo_file_registry",
                          side_effect=self.http_error(500)) as fetch, \
                    patch("funmirbench.vector_sources.time.sleep"):
                root = pathlib.Path(tmp)
                files = registry(record, root)
                self.assertEqual(files[name]["checksum"], "md5:" + checksum)
                self.assertIn(f"/records/{record}/files/{name}?download=1", files[name]["url"])
                self.assertEqual(registry(record, root), files)
                self.assertEqual(fetch.call_count, 3)

    def test_metadata_missing_record_fails_without_retry_or_substitution(self):
        with tempfile.TemporaryDirectory() as tmp, \
                patch("funmirbench.vector_sources.fetch_zenodo_file_registry",
                      side_effect=self.http_error(404)) as fetch, \
                patch("funmirbench.vector_sources.time.sleep") as sleep:
            with self.assertRaises(requests.HTTPError):
                registry("23839843", pathlib.Path(tmp))
            self.assertEqual(fetch.call_count, 1)
            sleep.assert_not_called()

    def test_metadata_unknown_record_has_no_pinned_fallback(self):
        with tempfile.TemporaryDirectory() as tmp, \
                patch("funmirbench.vector_sources.fetch_zenodo_file_registry", side_effect=self.http_error(503)), \
                patch("funmirbench.vector_sources.time.sleep"):
            with self.assertRaises(requests.HTTPError):
                registry("99999999", pathlib.Path(tmp))
            self.assertFalse((pathlib.Path(tmp) / "zenodo/99999999/registry.json").exists())

    def meta(self, payload):
        return {"filename": "data.zip", "size": len(payload),
                "checksum": "md5:" + hashlib.md5(payload).hexdigest(), "url": "https://example.test/data.zip"}

    def test_resume_and_cached_verification(self):
        payload = b"abcdefgh"
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            responses = [Response(payload[:3], broken=True),
                         Response(payload[3:], 206, {"Content-Range": "bytes 3-7/8"})]
            with patch("funmirbench.vector_sources.requests.get", side_effect=responses) as get, \
                    patch("funmirbench.vector_sources.time.sleep"):
                path = download(self.meta(payload), root)
                self.assertEqual(get.call_args_list[1].kwargs["headers"], {"Range": "bytes=3-"})
            self.assertEqual(path.read_bytes(), payload)
            with patch("funmirbench.vector_sources.requests.get", side_effect=AssertionError("downloaded again")):
                self.assertEqual(download(self.meta(payload), root), path)

    def test_ignored_range_restarts_instead_of_appending(self):
        payload = b"abcdefgh"
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            (root / "data.zip.part").write_bytes(payload[:3])
            with patch("funmirbench.vector_sources.requests.get", return_value=Response(payload)):
                self.assertEqual(download(self.meta(payload), root).read_bytes(), payload)

    def test_bad_checksum_is_never_promoted(self):
        payload = b"abcdefgh"
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            with patch("funmirbench.vector_sources.requests.get", side_effect=lambda *a, **kw: Response(b"xxxxxxxx")), \
                    patch("funmirbench.vector_sources.time.sleep"):
                with self.assertRaisesRegex(ValueError, "checksum"):
                    download(self.meta(payload), root)
            self.assertFalse((root / "data.zip").exists())

    def test_archive_traversal_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            archive = root / "tables.zip"
            with zipfile.ZipFile(archive, "w") as z:
                z.writestr("../escape.txt", "bad")
            with self.assertRaisesRegex(ValueError, "unexpected"):
                unpack_tables(archive, root)
            self.assertFalse((root.parent / "escape.txt").exists())

    def test_changed_extracted_files_are_restored_from_verified_archive(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            archive = root / "tables.zip"
            with zipfile.ZipFile(archive, "w") as z:
                z.writestr("benchmark_data/comparisons.tsv", "original")
                z.writestr("benchmark_data/experiments.tsv", "original")
            tables = unpack_tables(archive, root)
            (tables / "comparisons.tsv").write_text("changed")
            self.assertEqual((unpack_tables(archive, root) / "comparisons.tsv").read_text(), "original")


class ExploreFlowTests(unittest.TestCase):
    def test_reports_warm_cache_and_replot_without_network_or_prediction_parsing(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            tables = root / "tables"
            tables.mkdir()
            genes = [f"g{i}" for i in range(300)]
            for exp, response in (("e1", np.linspace(-1, 1, 300)), ("e2", np.zeros(300))):
                pd.DataFrame({"gene_id": genes, "eligible": True, "oriented_logFC_shrunk": response,
                              "control_mean_normalized_count": np.arange(300)}).to_csv(
                                  tables / f"{exp}.tsv.gz", sep="\t", index=False)
            pd.DataFrame({"experiment_id": ["e1", "e2"], "table": ["e1.tsv.gz", "e2.tsv.gz"]}).to_csv(
                tables / "experiments.tsv", sep="\t", index=False)
            pd.DataFrame({"experiment_id": ["e1", "e2"], "mirna_name": "m1", "perturbation": "OE",
                          "action": ["keep", "flag"], "dependency_group": "shared"}).to_csv(
                              tables / "comparisons.tsv", sep="\t", index=False)
            pd.DataFrame({"mirna_name": ["m1", "m2", "m3"],
                          "mirna_sequence": ["UAAAAAAA", "UCCCCCCC", "UGGGGGGG"]}).to_csv(
                              root / "mirnas.tsv", sep="\t", index=False)
            pd.DataFrame({"tool_id": ["normal", "inverse", "absent"], "official_name": ["Normal", "Inverse", "Absent"],
                          "score_direction": ["higher_is_stronger", "lower_is_stronger", "higher_is_stronger"],
                          "predictor_output_path": "unused"}).to_csv(root / "predictors.tsv", sep="\t", index=False)
            scores = pd.concat([pd.DataFrame({"Ensembl_ID": genes, "miRNA_Name": m, "Score": v}) for m, v in
                                (("m1", np.arange(300)), ("m2", np.arange(300)[::-1]), ("m3", np.arange(300) % 2))])
            scores.to_csv(root / "scores.tsv.gz", sep="\t", index=False)
            scores.assign(miRNA_Name="outside").to_csv(root / "absent.tsv.gz", sep="\t", index=False)
            sources = {t: {"path": str(root / ("absent.tsv.gz" if t == "absent" else "scores.tsv.gz")),
                            "origin": "local fixture"} for t in ("normal", "inverse", "absent")}
            provenance = {"benchmark_record": "23839843", "prediction_sources": sources}
            config = root / "explore.yaml"
            config.write_text(yaml.safe_dump({"mirnas_tsv": "mirnas.tsv", "predictions_tsv": "predictors.tsv",
                                              "cache_dir": "cache", "out_dir": "results",
                                              "exploration": {"detail_experiments": 0, "dpi": 70}}))
            with patch("funmirbench.vector_explore.prepare_inputs", return_value=(tables, sources, provenance)):
                out = run_exploration(config)
                values = pd.read_csv(out / "specificity.tsv", sep="\t")
                normal = values[(values.predictor == "normal") & (values.metric == "spearman") &
                                (values.experiment_id == "e1")].iloc[0]
                inverse = values[(values.predictor == "inverse") & (values.metric == "spearman") &
                                 (values.experiment_id == "e1")].iloc[0]
                self.assertAlmostEqual(normal.alignment, 1.)
                self.assertAlmostEqual(inverse.alignment, -1.)
                self.assertTrue(values[values.predictor == "absent"].alignment.isna().all())
                self.assertTrue(values[values.experiment_id == "e2"].alignment.isna().all())
                self.assertTrue((out / "REPORT.pdf").read_bytes().rstrip().endswith(b"%%EOF"))
                report = (out / "REPORT.html").read_text()
                self.assertIn("data:image/png;base64,", report)
                self.assertIn("Checks still pending", report)
                with patch("funmirbench.vector_engine.load_matrix", side_effect=AssertionError("recomputed")), \
                        patch("funmirbench.vector_report.figures", side_effect=AssertionError("redrawn")):
                    warm = run_exploration(config)
                self.assertNotEqual(warm, out)
                self.assertTrue(json.loads((warm / "timings.json").read_text())["analysis_cache_hit"])
            with patch("funmirbench.vector_sources.requests.get", side_effect=AssertionError("network")), \
                    patch("funmirbench.vector_engine.load_matrix", side_effect=AssertionError("parsed predictions")):
                replot(out, details=0, dpi=70)
                with patch("funmirbench.vector_report.figures", side_effect=AssertionError("redrawn")):
                    replot(out, details=0, dpi=70)
                # Editing saved numerical results must invalidate report-only figures.
                from funmirbench.vector_report import figures
                changed = pd.read_csv(out / "specificity.tsv", sep="\t")
                changed.loc[0, "alignment"] = .25
                changed.to_csv(out / "specificity.tsv", sep="\t", index=False, na_rep="NA")
                with patch("funmirbench.vector_report.figures", wraps=figures) as redraw:
                    replot(out, details=0, dpi=70)
                    redraw.assert_called_once()
            self.assertTrue((out / "REPORT.pdf").read_bytes().rstrip().endswith(b"%%EOF"))


if __name__ == "__main__":
    unittest.main()
