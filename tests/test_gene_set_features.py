import gzip
import pathlib
import sys
import tempfile
import types
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

from funmirbench.gene_set_features import (
    GeneSetFeatureSettings,
    assign_coverage_groups,
    build_utr3_feature_table,
    ensure_feature_resources,
    format_pvalue,
    intersection_gene_set,
    parse_settings,
    summarize_features,
)
from funmirbench.plots import plot_gene_set_features
from funmirbench.benchmark import run_benchmark


class GeneSetFeatureTests(unittest.TestCase):
    def test_longest_protein_coding_transcript_supplies_length_and_conservation(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            gtf = root / "genes.gtf.gz"
            lines = [
                'chr1\ttest\tthree_prime_utr\t1\t5\t.\t+\t.\tgene_id "ENSG1.2"; transcript_id "ENST1.1"; transcript_biotype "protein_coding";',
                'chr1\ttest\tthree_prime_utr\t4\t10\t.\t+\t.\tgene_id "ENSG1.2"; transcript_id "ENST1.1"; transcript_biotype "protein_coding";',
                'chr1\ttest\tthree_prime_utr\t50\t61\t.\t+\t.\tgene_id "ENSG1.2"; transcript_id "ENST2.1"; transcript_biotype "protein_coding";',
                'chr1\ttest\tthree_prime_utr\t100\t500\t.\t+\t.\tgene_id "ENSG1.2"; transcript_id "ENST3.1"; transcript_biotype "lncRNA";',
            ]
            with gzip.open(gtf, "wt", encoding="utf-8") as handle:
                handle.write("\n".join(lines) + "\n")

            class BigWig:
                def chroms(self):
                    return {"chr1": 1000}

                def values(self, chrom, start, end):
                    return [0.25] * (end - start)

                def close(self):
                    pass

            fake_module = types.SimpleNamespace(open=lambda _path: BigWig())
            with patch.dict(sys.modules, {"pyBigWig": fake_module}):
                features = build_utr3_feature_table(gtf, root / "track.bw", pd.Index(["ENSG1", "ENSG2"]))

        first = features.set_index("gene_id").loc["ENSG1"]
        self.assertEqual(first.transcript_id, "ENST2")
        self.assertEqual(first.utr3_length_nt, 12)
        self.assertAlmostEqual(first.utr3_mean_phastcons, 0.25)
        self.assertEqual(first.utr3_scored_bases, 12)
        self.assertTrue(pd.isna(features.set_index("gene_id").loc["ENSG2", "utr3_length_nt"]))

    def test_intersection_uses_only_evaluated_mirnas_and_benchmark_genes(self):
        predictions = {
            "a": {"m1": pd.Series([1, 2], index=["g1", "g2"]),
                  "m2": pd.Series([1], index=["g3"])},
            "b": {"m1": pd.Series([1, 2], index=["g1", "g3"]),
                  "m2": pd.Series([1], index=["g2"])},
        }
        result = intersection_gene_set(predictions, {"m1"}, pd.Index(["g1", "g2", "g3"]))
        self.assertEqual(result, {"g1"})

    def test_summary_reports_matched_counts_medians_and_two_sided_test(self):
        frame = pd.DataFrame({
            "gene_id": ["g1", "g2", "g3", "g4"],
            "coverage_group": ["IGS", "IGS", "non-IGS", "non-IGS"],
            "utr3_length_nt": [100., 200., 10., np.nan],
            "utr3_mean_phastcons": [.8, .9, .1, .2],
        })
        summary = summarize_features(frame).set_index("feature_column")
        length = summary.loc["utr3_length_nt"]
        self.assertEqual(length.n_igs_with_value, 2)
        self.assertEqual(length.n_non_igs_with_value, 1)
        self.assertEqual(length.median_igs, 150)
        self.assertEqual(length.median_non_igs, 10)
        self.assertTrue(0 <= length.pvalue_two_sided <= 1)
        self.assertAlmostEqual(length.log10_pvalue_two_sided, np.log10(length.pvalue_two_sided))

    def test_underflowed_pvalue_is_preserved_and_formatted_from_log10(self):
        frame = pd.DataFrame({"coverage_group": ["IGS"] * 1000 + ["non-IGS"] * 1000,
                              "utr3_length_nt": [2.] * 1000 + [1.] * 1000,
                              "utr3_mean_phastcons": [.9] * 1000 + [.1] * 1000})
        summary = summarize_features(frame)
        length = summary.iloc[0]
        self.assertTrue(pd.isna(length.pvalue_two_sided))
        self.assertLess(length.log10_pvalue_two_sided, -300)
        self.assertEqual(length.pvalue_two_sided_display,
                         format_pvalue(length.pvalue_two_sided, length.log10_pvalue_two_sided))
        self.assertRegex(length.pvalue_two_sided_display, r"\de-\d+")

    def test_feature_groups_cover_the_full_gene_universe(self):
        features = pd.DataFrame({"gene_id": ["g1"], "transcript_id": ["tx1"],
                                 "utr3_length_nt": [100], "utr3_mean_phastcons": [.5],
                                 "utr3_scored_bases": [100]})
        rows = assign_coverage_groups(features, pd.Index(["g1", "g2"]), {"g1"})
        self.assertEqual(rows.coverage_group.tolist(), ["IGS", "non-IGS"])
        self.assertTrue(pd.isna(rows.loc[1, "utr3_mean_phastcons"]))

    def test_settings_default_off_and_rejects_non_boolean_enabled(self):
        settings = parse_settings(None, pathlib.Path("/repo"))
        self.assertFalse(settings.enabled)
        self.assertEqual(settings.gtf_path, pathlib.Path("/repo/data/experiments/raw/refs/ensembl_v115/Homo_sapiens.GRCh38.115.gtf.gz"))
        with self.assertRaisesRegex(ValueError, "enabled must be a boolean"):
            parse_settings({"enabled": "true"}, pathlib.Path("/repo"))

    def test_missing_annotation_resources_are_downloaded_and_cached(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            settings = GeneSetFeatureSettings(True, root / "annotation.gtf.gz", root / "track.bw")

            class Response:
                def __init__(self, payload):
                    self.payload = payload
                    self.headers = {"content-length": str(len(payload))}

                def __enter__(self):
                    return self

                def __exit__(self, *_args):
                    return None

                def raise_for_status(self):
                    return None

                def iter_content(self, chunk_size):
                    return [self.payload]

            responses = [Response(b"\x1f\x8bgtf"), Response(b"\x26\xfc\x8f\x88bigwig")]
            with patch.dict(sys.modules, {"pyBigWig": types.SimpleNamespace()}):
                with patch("funmirbench.gene_set_features.requests.get", side_effect=responses) as get:
                    ensure_feature_resources(settings)
                    ensure_feature_resources(settings)

            self.assertEqual(get.call_count, 2)
            self.assertEqual(settings.gtf_path.read_bytes(), b"\x1f\x8bgtf")
            self.assertEqual(settings.phastcons_path.read_bytes(), b"\x26\xfc\x8f\x88bigwig")

    def test_plot_contains_both_feature_distributions(self):
        frame = pd.DataFrame({"coverage_group": ["IGS", "IGS", "non-IGS", "non-IGS"],
                              "utr3_length_nt": [100., 200., 50., 75.],
                              "utr3_mean_phastcons": [.8, .9, .2, .3]})
        summary = summarize_features(frame)
        with patch("funmirbench.plots.save") as save:
            plot_gene_set_features(frame, summary, pathlib.Path("unused"))
        fig = save.call_args.args[0]
        try:
            self.assertEqual(len(fig.axes), 2)
            self.assertIn("3′UTR length", fig.axes[0].get_title())
            self.assertIn("phastCons100way", fig.axes[1].get_xlabel())
        finally:
            import matplotlib.pyplot as plt
            plt.close(fig)

    def test_enabled_feature_comparison_is_written_into_the_main_report(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            data = root / "benchmark_data"
            (data / "genes").mkdir(parents=True)
            genes = [f"ENSG{i:011d}" for i in range(50)]
            states = ["responsive"] * 15 + ["no_response_above_margin"] * 20 + ["unresolved"] * 15
            table_path = data / "genes/E1.tsv.gz"
            pd.DataFrame({"gene_id": genes, "eligible": True,
                          "oriented_logFC_shrunk": np.linspace(2, -1, len(genes)),
                          "control_mean_normalized_count": np.arange(len(genes)) + 1,
                          "state_d05": states}).to_csv(table_path, sep="\t", index=False)
            pd.DataFrame([{"experiment_id": "E1", "table": "genes/E1.tsv.gz"}]).to_csv(
                data / "experiments.tsv", sep="\t", index=False)
            pd.DataFrame([{"experiment_id": "E1", "mirna_name": "miR-1", "perturbation": "OE",
                           "action": "keep", "dependency_group": "g1", "flags": "", "cell_line": "HeLa"}]).to_csv(
                data / "comparisons.tsv", sep="\t", index=False)
            pd.DataFrame([{"mirna_name": "miR-1", "mirna_sequence": "UAGCAGCAC"}]).to_csv(
                root / "mirnas.tsv", sep="\t", index=False)
            predictions = []
            for name, count in (("a", 40), ("b", 30)):
                path = root / f"{name}.tsv"
                pd.DataFrame({"Ensembl_ID": genes[:count], "miRNA_Name": "miR-1",
                              "Score": np.arange(count, dtype=float)}).to_csv(path, sep="\t", index=False)
                predictions.append({"tool_id": name, "official_name": name, "score_direction": "higher_is_stronger",
                                   "predictor_output_path": str(path)})
            pd.DataFrame(predictions).to_csv(root / "predictors.tsv", sep="\t", index=False)
            (root / "genes.gtf.gz").touch()
            (root / "track.bw").touch()
            config_path = root / "config.yaml"
            config_path.write_text(
                "input_dir: benchmark_data\nmirnas_tsv: mirnas.tsv\npredictions_tsv: predictors.tsv\n"
                "out_dir: results\nevaluation: {min_per_state: 1, bootstrap: 2}\n"
                "reports:\n  detailed: false\n  gene_set_features:\n    enabled: true\n"
                "    gtf_path: genes.gtf.gz\n    phastcons_path: track.bw\n",
                encoding="utf-8")
            features = pd.DataFrame({"gene_id": genes, "transcript_id": ["tx"] * len(genes),
                                     "utr3_length_nt": np.arange(len(genes)) + 100,
                                     "utr3_mean_phastcons": np.linspace(0, 1, len(genes)),
                                     "utr3_scored_bases": [100] * len(genes)})
            with patch("funmirbench.benchmark.build_utr3_feature_table", return_value=features):
                out = run_benchmark(config_path)
            self.assertTrue((out / "plots/gene_set_features.png").is_file())
            self.assertTrue((out / "gene_set_features.tsv").is_file())
            summary = pd.read_csv(out / "gene_set_feature_summary.tsv", sep="\t")
            self.assertEqual(summary.n_igs.iloc[0], 30)
            self.assertEqual(summary.n_non_igs.iloc[0], 20)
            self.assertEqual((out / "REPORT.pdf").read_bytes()[:4], b"%PDF")


if __name__ == "__main__":
    unittest.main()
