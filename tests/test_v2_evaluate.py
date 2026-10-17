import argparse
import pathlib
import tempfile
import unittest

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

from funmirbench import v2_evaluate as v2

MIR = "UAGCAGCACGUAAAUAUUGGCG"  # miR-16-5p, seed nt 2-8 AGCAGCA


class SeedSiteTests(unittest.TestCase):
    def test_site_hierarchy(self):
        s6 = MIR[1:7].translate(v2.COMPLEMENT)[::-1]  # pairs with nt 2-7
        m8 = MIR[7].translate(v2.COMPLEMENT)
        self.assertEqual(v2.seed_site_counts("GG" + m8 + s6 + "AGG", MIR), (1, 0, 0, 0))
        self.assertEqual(v2.seed_site_counts("GG" + m8 + s6 + "CGG", MIR), (0, 1, 0, 0))
        self.assertEqual(v2.seed_site_counts("GGC" + s6 + "AGG", MIR), (0, 0, 1, 0))
        self.assertEqual(v2.seed_site_counts("GGC" + s6 + "CGG", MIR), (0, 0, 0, 1))
        self.assertEqual(v2.seed_site_counts("CCCCCCCC", MIR), (0, 0, 0, 0))

    def test_seed_score_weights_and_missing_utr(self):
        s6 = MIR[1:7].translate(v2.COMPLEMENT)[::-1]
        m8 = MIR[7].translate(v2.COMPLEMENT)
        utrs = {"G1": m8 + s6 + "A" + "CC" + "C" + s6 + "C"}
        score = v2.seed_score(utrs, pd.Index(["G1", "G2"]), MIR)
        self.assertEqual(score["G1"], 3.5)
        self.assertEqual(score["G2"], 0.0)


class MetricTests(unittest.TestCase):
    def test_auroc_matches_sklearn_with_ties(self):
        rng = np.random.default_rng(1)
        score = rng.integers(0, 5, 200).astype(float)
        pos = rng.random(200) < 0.3
        self.assertAlmostEqual(v2.auroc(score, pos), roc_auc_score(pos, score))

    def test_binned_auroc_ignores_between_bin_signal(self):
        # The score separates the classes only through the bin, so within-bin AUROC is 0.5.
        bins = np.repeat([1, 2], 100)
        pos = np.r_[np.zeros(80, bool), np.ones(20, bool), np.zeros(20, bool), np.ones(80, bool)]
        score = bins.astype(float)
        self.assertGreater(v2.auroc(score, pos), 0.7)
        self.assertAlmostEqual(v2.binned_auroc(score, pos, bins), 0.5)

    def test_rank_last_puts_unscored_genes_below(self):
        genes = pd.Index(["a", "b", "c"])
        s = v2.rank_last(pd.Series({"a": -5.0, "b": 2.0}), genes)
        self.assertLess(s["c"], s["a"])
        self.assertTrue((v2.rank_last(pd.Series(dtype=float), genes) == 0).all())

    def test_group_bootstrap(self):
        rng = np.random.default_rng(0)
        v = pd.Series([1.0, 1.0, 3.0, 3.0])
        self.assertTrue(np.isnan(v2.group_bootstrap(v, pd.Series(["g"] * 4), n=50, rng=rng)[0]))
        lo, hi = v2.group_bootstrap(v, pd.Series(["a", "a", "b", "b"]), n=500, rng=rng)
        self.assertGreaterEqual(lo, 1.0)
        self.assertLessEqual(hi, 3.0)


class EndToEndTests(unittest.TestCase):
    def test_small_run(self):
        rng = np.random.default_rng(3)
        genes = [f"ENSG{i:011d}" for i in range(400)]
        s6 = MIR[1:7].translate(v2.COMPLEMENT)[::-1]
        with tempfile.TemporaryDirectory() as tmp:
            tmp = pathlib.Path(tmp)
            final = tmp / "final"
            (final / "genes").mkdir(parents=True)
            target = rng.random(400) < 0.3
            rows = []
            for exp, kind in (("E1", "OE"), ("E2", "KO")):
                o = np.where(target, 1.2, 0.0) + rng.normal(0, 0.3, 400)
                state = np.where(o > 1.0, "responsive", np.where(o < 0.3, "no_response_above_margin", "unresolved"))
                pd.DataFrame({"gene_id": genes, "eligible": True, "oriented_logFC": o, "oriented_logFC_shrunk": o,
                              "control_mean_normalized_count": rng.gamma(2, 50, 400), "state_d05": state,
                              "state_d10": state}).to_csv(final / "genes" / f"{exp}.tsv.gz", sep="\t", index=False)
                rows.append({"experiment_id": exp, "table": f"genes/{exp}.tsv.gz"})
            pd.DataFrame(rows).to_csv(final / "experiments.tsv", sep="\t", index=False)
            pd.DataFrame({"experiment_id": ["E1", "E2"], "mirna_name": ["hsa-miR-16-5p"] * 2,
                          "perturbation": ["OE", "KO"], "action": ["keep", "flag"],
                          "dependency_group": ["g1", "g2"]}).to_csv(final / "comparisons.tsv", sep="\t", index=False)
            with open(tmp / "utr.fa", "w") as fh:
                for g, t in zip(genes, target):
                    fh.write(f">{g}|transcript=T\n{'CCCC' + (s6 + 'A' if t else 'GGGGGGG') + 'CCCC'}\n")
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
            pd.DataFrame({"tool_id": ["tool"], "score_direction": ["higher_is_stronger"],
                          "predictor_output_path": [str(pred_dir / "tool.tsv")]}).to_csv(
                tmp / "predictions.tsv", sep="\t", index=False)
            args = argparse.Namespace(final_dir=final, utr_fasta=tmp / "utr.fa", out=tmp / "out", comparisons=None,
                                      mirnas=tmp / "mirnas.tsv", predictions=tmp / "predictions.tsv", tools=None,
                                      min_per_state=30, bootstrap=100, seed=0)
            v2.evaluate(args)
            per = pd.read_csv(tmp / "out/per_comparison.tsv", sep="\t")
            decoy = pd.read_csv(tmp / "out/decoy.tsv", sep="\t")
            tool = per[per.predictor == "tool"].set_index("experiment_id")
            self.assertAlmostEqual(tool.loc["E1", "coverage"], 0.75)
            seed = per[per.predictor == "seed_site"].set_index("experiment_id")
            self.assertGreater(seed.loc["E1", "auroc_d05"], 0.8)
            # Unscored genes rank last, so a tool that covers 75% of the genes scores below full coverage.
            self.assertGreater(tool.loc["E1", "auroc_d05"], 0.6)
            self.assertLess(tool.loc["E1", "auroc_d05"], seed.loc["E1", "auroc_d05"])
            # miR-15b-5p shares the nt 2-8 seed, so it is not a decoy; only miR-1-3p is.
            self.assertTrue((decoy.n_decoys == 1).all())
            summary = pd.read_csv(tmp / "out/summary.tsv", sep="\t")
            self.assertEqual(set(summary.subset), {"all", "unflagged"})
            self.assertTrue((tmp / "out/paired_vs_reference.tsv").exists())


if __name__ == "__main__":
    unittest.main()
