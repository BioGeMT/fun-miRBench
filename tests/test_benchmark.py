import pathlib
import tempfile
import unittest

import numpy as np
import pandas as pd
import yaml
from sklearn.metrics import roc_auc_score

from funmirbench import benchmark, plots, report_pdf

MIR = "UAGCAGCACGUAAAUAUUGGCG"  # miR-16-5p, seed nt 2-8 AGCAGCA


class MetricTests(unittest.TestCase):
    def test_auroc_matches_sklearn_with_ties(self):
        rng = np.random.default_rng(1)
        score = rng.integers(0, 5, 200).astype(float)
        pos = rng.random(200) < 0.3
        self.assertAlmostEqual(benchmark.auroc(score, pos), roc_auc_score(pos, score))

    def test_binned_auroc_ignores_between_bin_signal(self):
        # The score separates the classes only through the bin, so within-bin AUROC is 0.5.
        bins = np.repeat([1, 2], 100)
        pos = np.r_[np.zeros(80, bool), np.ones(20, bool), np.zeros(20, bool), np.ones(80, bool)]
        score = bins.astype(float)
        self.assertGreater(benchmark.auroc(score, pos), 0.7)
        self.assertAlmostEqual(benchmark.binned_auroc(score, pos, bins), 0.5)

    def test_rank_last_puts_unscored_genes_below(self):
        genes = pd.Index(["a", "b", "c"])
        s = benchmark.rank_last(pd.Series({"a": -5.0, "b": 2.0}), genes)
        self.assertLess(s["c"], s["a"])
        self.assertTrue((benchmark.rank_last(pd.Series(dtype=float), genes) == 0).all())

    def test_group_bootstrap(self):
        rng = np.random.default_rng(0)
        v = pd.Series([1.0, 1.0, 3.0, 3.0])
        self.assertTrue(np.isnan(benchmark.group_bootstrap(v, pd.Series(["g"] * 4), n=50, rng=rng)[0]))
        lo, hi = benchmark.group_bootstrap(v, pd.Series(["a", "a", "b", "b"]), n=500, rng=rng)
        self.assertGreaterEqual(lo, 1.0)
        self.assertLessEqual(hi, 3.0)


class RecoveryTests(unittest.TestCase):

    def test_recall_at_small_budget_counts_only_the_top_genes(self):
        n = 40
        score = np.arange(n, 0, -1, dtype=float)
        positive = np.zeros(n, dtype=bool)
        positive[[0, 30]] = True  # one positive in the top 10, one below
        recalls = benchmark.recall_at_budgets(score, positive, np.array([f"g{i:02d}" for i in range(n)]))
        self.assertEqual(recalls[10], 0.5)
        self.assertEqual(recalls[50], 1.0)

    def test_precision_at_recall_admits_tied_genes_together(self):
        score = np.array([3.0, 2.0, 2.0, 2.0, 1.0, 1.0])
        positive = np.array([False, True, False, False, True, False])
        row = benchmark.precision_at_recall(score, positive)  # half of 2 positives needs the block at 2.0: 1 of 4 genes
        self.assertEqual(row["precision"], 0.25)
        self.assertEqual(row["false_positives_per_true_positive"], 3.0)
        self.assertAlmostEqual(row["precision_over_random"], 0.25 / (2 / 6))


class PlotTests(unittest.TestCase):
    def test_experiments_won_skips_missing_values_and_keeps_tool_order_on_ties(self):
        per = pd.DataFrame({"experiment_id": ["E1", "E1", "E2", "E2", "E3", "E3"],
                            "predictor": ["a", "b"] * 3,
                            "auroc_d05": [0.6, 0.6, 0.5, 0.7, np.nan, np.nan]})
        self.assertEqual(plots.experiments_won(per, ["b", "a"], "auroc_d05").to_dict(), {"b": 2, "a": 0})


class ReportTests(unittest.TestCase):
    def test_tool_means_use_all_comparisons_and_the_drawn_tools(self):
        tools = pd.DataFrame({"tool_id": ["targetscan_family_expanded", "miraw"],
                              "official_name": ["TargetScan v8", "miRAW"]})
        summary = pd.DataFrame([
            {"subset": subset, "predictor": tool, "metric": metric, "mean": value}
            for subset, value in (("all", 0.5), ("unflagged", 0.9))
            for tool in ("targetscan_family_expanded", "miraw", "random")
            for metric in ("coverage", "auroc_d05", "spearman_shrunk")])
        means = report_pdf.tool_means(summary, tools)
        self.assertEqual(list(means.index), ["TargetScan v8", "miRAW"])
        self.assertTrue((means == 0.5).all(axis=None))


class EndToEndTests(unittest.TestCase):
    def test_small_run(self):
        rng = np.random.default_rng(3)
        genes = [f"ENSG{i:011d}" for i in range(400)]
        with tempfile.TemporaryDirectory() as tmp:
            tmp = pathlib.Path(tmp)
            final = tmp / "final"
            (final / "genes").mkdir(parents=True)
            target = rng.random(400) < 0.3
            rows = []
            for exp in ("E1", "E2", "E3"):
                o = np.where(target, 1.2, 0.0) + rng.normal(0, 0.3, 400)
                state = np.where(o > 1.0, "responsive", np.where(o < 0.3, "no_response_above_margin", "unresolved"))
                pd.DataFrame({"gene_id": genes, "eligible": True, "oriented_logFC_shrunk": o,
                              "control_mean_normalized_count": rng.gamma(2, 50, 400), "state_d05": state}).to_csv(
                    final / "genes" / f"{exp}.tsv.gz", sep="\t", index=False)
                rows.append({"experiment_id": exp, "table": f"genes/{exp}.tsv.gz"})
            pd.DataFrame(rows).to_csv(final / "experiments.tsv", sep="\t", index=False)
            pd.DataFrame({"experiment_id": ["E1", "E2", "E3", "E4"], "mirna_name": ["hsa-miR-16-5p"] * 4,
                          "cell_line": ["HeLa"] * 4, "perturbation": ["OE", "OE", "KO", "OE"],
                          "action": ["keep", "flag", "keep", "set_aside"], "flags": ["", "", "", "global_shift"],
                          "dependency_group": ["g1", "g2", "g3", "g4"]}).to_csv(final / "comparisons.tsv", sep="\t", index=False)
            pd.DataFrame({"mirna_name": ["hsa-miR-16-5p", "hsa-miR-15b-5p", "hsa-miR-1-3p"],
                          "mirna_sequence": [MIR, "UAGCAGCACAUCAUGGUUUACA", "UGGAAUGUAAAGAAGUAUGUAU"]}).to_csv(
                tmp / "mirnas.tsv", sep="\t", index=False)
            pred = pd.concat([
                pd.DataFrame({"Ensembl_ID": genes[:300], "miRNA_Name": "hsa-miR-16-5p", "Score": target[:300] * 1.0}),
                pd.DataFrame({"Ensembl_ID": genes, "miRNA_Name": "hsa-miR-1-3p", "Score": rng.random(400)}),
                pd.DataFrame({"Ensembl_ID": genes, "miRNA_Name": "hsa-miR-15b-5p", "Score": target * 1.0})])
            pred_dir = tmp / "data/predictions/tool"
            pred_dir.mkdir(parents=True)
            pred.to_csv(pred_dir / "tool.tsv", sep="\t", index=False)
            pd.DataFrame({"tool_id": ["tool"], "official_name": ["Tool"], "score_direction": ["higher_is_stronger"],
                          "predictor_output_path": [str(pred_dir / "tool.tsv")]}).to_csv(
                tmp / "predictions.tsv", sep="\t", index=False)
            (tmp / "benchmark.yaml").write_text(yaml.safe_dump({
                "input_dir": "final", "mirnas_tsv": "mirnas.tsv",
                "predictions_tsv": "predictions.tsv", "evaluation": {"bootstrap": 100}, "out_dir": "results"}))
            with self.assertLogs("funmirbench.benchmark", "WARNING") as logs:
                out = benchmark.run_benchmark(tmp / "benchmark.yaml")
            self.assertIn("skipping 1 others: E3", logs.output[0])
            self.assertEqual(out.parent, (tmp / "results").resolve())
            self.assertTrue((out / "benchmark_config.yaml").exists())
            per = pd.read_csv(out / "per_comparison.tsv", sep="\t")
            tool = per[per.predictor == "tool"].set_index("experiment_id")
            self.assertAlmostEqual(tool.loc["E1", "coverage"], 0.75)
            # Unscored genes rank last, so a tool that covers 75% of the genes scores below full coverage.
            self.assertGreater(tool.loc["E1", "auroc_d05"], 0.6)
            self.assertLess(tool.loc["E1", "auroc_d05"], 0.9)
            # miR-15b-5p shares the nt 2-8 seed, so it is not a decoy; only miR-1-3p is, and a z needs two decoys.
            self.assertTrue(per.decoy_z_shrunk.isna().all())
            summary = pd.read_csv(out / "summary.tsv", sep="\t")
            self.assertEqual(set(summary.subset), {"all", "unflagged"})
            spearman = summary[summary.metric == "spearman_shrunk"].set_index(["subset", "predictor"]).n_comparisons
            self.assertEqual(spearman["all", "tool"], 2)  # the KO comparison E3 is not scored
            self.assertEqual(spearman["unflagged", "tool"], 1)
            self.assertEqual(set(per.experiment_id), {"E1", "E2"})
            for stem in ("coverage", "oe_all_genes", "recovery"):
                self.assertTrue((out / f"plots/{stem}.png").exists())
                self.assertTrue((out / f"plots/{stem}.svg").exists())
            coverage = pd.read_csv(out / "coverage.tsv", sep="\t").set_index("predictor")
            self.assertEqual(coverage.loc["tool", "genes_covered"], 300)
            self.assertEqual(coverage.loc["tool", "genes_measured"], 400)
            self.assertEqual(coverage.loc["tool", "experiments_retained"], 2)
            self.assertEqual((out / "REPORT.pdf").read_bytes()[:4], b"%PDF")

    def test_missing_final_tables_name_the_download_command(self):
        with tempfile.TemporaryDirectory() as tmp:
            config = pathlib.Path(tmp) / "benchmark.yaml"
            config.write_text(yaml.safe_dump({"input_dir": "final"}))
            with self.assertRaisesRegex(FileNotFoundError, "funmirbench-download-data"):
                benchmark.run_benchmark(config)


if __name__ == "__main__":
    unittest.main()
