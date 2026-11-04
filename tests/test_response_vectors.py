import pathlib
import tempfile
import unittest

import numpy as np
import pandas as pd
import yaml

from funmirbench.response_vectors import alignment, build_vectors, run_analysis, specificity


class VectorTests(unittest.TestCase):
    def table(self):
        return pd.DataFrame({"eligible": [True, True, True, False, True, True],
                             "oriented_logFC_shrunk": [3., 2., 1., 9., np.nan, -1.],
                             "control_mean_normalized_count": [10.] * 6},
                            index=pd.Index(["a", "b", "c", "d", "e", "f"], name="gene_id"))

    def test_alignment_uses_gene_ids_and_preserves_both_missingness_masks(self):
        table = self.table()
        vectors = build_vectors(table, pd.Series({"b": 8., "a": 8., "c": 2., "outside": 999.}))
        self.assertEqual(list(vectors.index), list(table.index))
        self.assertTrue(vectors.loc[["d", "e"], "observed_response"].isna().all())
        self.assertFalse(vectors.loc["f", "predictor_scored"])
        self.assertTrue(np.isnan(vectors.loc["f", "prediction_score"]))
        self.assertLess(vectors.loc["f", "prediction_ordinal"], vectors.loc["c", "prediction_ordinal"])
        usable = vectors[vectors.response_usable]
        expected = alignment(vectors, pd.Series(1, index=table.index))["spearman"]
        self.assertAlmostEqual(np.dot(usable.response_rank_unit, usable.prediction_rank_unit), expected)
        self.assertAlmostEqual(np.linalg.norm(usable.response_rank_unit), 1.)
        # Gene order and a monotonic score transformation do not change alignment.
        shuffled = build_vectors(table.iloc[::-1], pd.Series({"a": 80., "b": 80., "c": 20.}))
        self.assertAlmostEqual(alignment(shuffled, pd.Series(1, index=table.index))["spearman"], expected)

    def test_missing_predictions_are_undefined_not_a_zero_alignment(self):
        vectors = build_vectors(self.table(), pd.Series(dtype=float))
        result = alignment(vectors, pd.Series(1, index=vectors.index))
        self.assertTrue(np.isnan(result["spearman"]))
        self.assertEqual(result["alignment_na_reason"], "constant_prediction")
        self.assertTrue(vectors.prediction_rank_unit.isna().all())

    def test_invalid_masks_and_duplicate_genes_are_rejected(self):
        table = self.table()
        table["eligible"] = "False"
        with self.assertRaisesRegex(ValueError, "boolean"):
            build_vectors(table, pd.Series(dtype=float))
        table = self.table()
        table.index = ["a"] * len(table)
        with self.assertRaisesRegex(ValueError, "unique"):
            build_vectors(table, pd.Series(dtype=float))

    def test_specificity_excludes_same_seed_and_missing_and_counts_ties(self):
        candidates = pd.DataFrame({"candidate_mirna": ["own", "sibling", "d1", "d2", "missing"],
                                   "candidate_seed": ["AAAAAAA", "AAAAAAA", "CCCCCCC", "GGGGGGG", "UUUUUUU"],
                                   "spearman": [0.5, 1., 0.5, 0.1, np.nan]})
        result = specificity(candidates, "own", "spearman")
        self.assertEqual(result["n_decoys"], 2)
        self.assertEqual(result["cognate_rank"], 1.5)
        self.assertEqual(result["fraction_decoys_better"], 0.25)
        candidates.loc[candidates.candidate_mirna == "d2", "spearman"] = 0.5
        result = specificity(candidates, "own", "spearman")
        self.assertTrue(np.isnan(result["decoy_z"]))
        self.assertEqual(result["specificity_na_reason"], "no_decoy_variation")
        candidates.loc[candidates.candidate_mirna == "own", "spearman"] = np.nan
        self.assertEqual(specificity(candidates, "own", "spearman")["specificity_na_reason"],
                         "cognate_alignment_unavailable")


class AnalysisTests(unittest.TestCase):
    def test_end_to_end_without_discrete_labels_and_with_missing_predictor(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            (root / "tables").mkdir()
            genes = [f"ENSG{i:011d}" for i in range(300)]
            response = np.linspace(-2, 2, len(genes))
            pd.DataFrame({"gene_id": genes, "eligible": True, "oriented_logFC_shrunk": response,
                          "control_mean_normalized_count": np.arange(len(genes)) + 1.}).to_csv(
                              root / "tables/e.tsv.gz", sep="\t", index=False)
            pd.DataFrame({"experiment_id": ["e"], "table": ["e.tsv.gz"]}).to_csv(
                root / "tables/experiments.tsv", sep="\t", index=False)
            pd.DataFrame({"experiment_id": ["e", "ko", "aside"], "mirna_name": ["own"] * 3,
                          "perturbation": ["OE", "KO", "OE"], "action": ["flag", "keep", "set_aside"],
                          "dependency_group": ["g1", "g2", "g3"]}).to_csv(
                              root / "tables/comparisons.tsv", sep="\t", index=False)
            pd.DataFrame({"mirna_name": ["own", "sibling", "opposite", "unscored"],
                          "mirna_sequence": ["UAAAAAAA", "CAAAAAAA", "UCCCCCCC", "UGGGGGGG"]}).to_csv(
                              root / "mirnas.tsv", sep="\t", index=False)
            # lower_is_stronger must flip the sign; reverse input row order to check alignment by ID.
            pd.concat([pd.DataFrame({"Ensembl_ID": genes, "miRNA_Name": mirna, "Score": scores})
                       for mirna, scores in (("own", -response), ("sibling", -response),
                                             ("opposite", response))]).iloc[::-1].to_csv(
                                                 root / "scores.tsv", sep="\t", index=False)
            pd.DataFrame({"Ensembl_ID": genes, "miRNA_Name": "outside_registry", "Score": 1.}).to_csv(
                root / "absent.tsv", sep="\t", index=False)
            pd.DataFrame({"tool_id": ["tool", "absent"], "score_direction": ["lower_is_stronger"] * 2,
                          "predictor_output_path": ["scores.tsv", "absent.tsv"]}).to_csv(
                              root / "predictors.tsv", sep="\t", index=False)
            config = root / "config.yaml"
            config.write_text(yaml.safe_dump({"input_dir": "tables", "mirnas_tsv": "mirnas.tsv",
                                              "predictions_tsv": "predictors.tsv", "out_dir": "results"}))
            out = run_analysis(config)
            self.assertEqual(out.parent.parent, root / "results")
            candidates = pd.read_csv(out / "candidate_alignment.tsv", sep="\t")
            self.assertEqual(set(candidates.experiment_id), {"e"})
            self.assertEqual(len(candidates), 8)
            own = candidates[(candidates.predictor == "tool") & candidates.is_cognate].iloc[0]
            self.assertAlmostEqual(own.spearman, 1.)
            self.assertAlmostEqual(own.spearman_within_tenths, 1.)
            summary = pd.read_csv(out / "specificity.tsv", sep="\t")
            row = summary[(summary.predictor == "tool") & (summary.metric == "spearman")].iloc[0]
            self.assertEqual(row.n_decoys, 1)  # missing predictions and same-seed sibling excluded
            self.assertEqual(row.cognate_rank, 1.)
            self.assertEqual(row.specificity_na_reason, "fewer_than_two_usable_decoys")
            self.assertTrue(summary[summary.predictor == "absent"].alignment.isna().all())
            files = pd.read_csv(out / "vector_files.tsv", sep="\t")
            vectors = pd.read_csv(out / files.iloc[0].file, sep="\t")
            self.assertEqual(vectors.gene_id.tolist(), genes)
            self.assertAlmostEqual(np.dot(vectors.response_rank_unit, vectors.prediction_rank_unit), 1.)
            self.assertEqual(set(pd.read_csv(out / "skipped_comparisons.tsv", sep="\t").experiment_id),
                             {"ko", "aside"})
            self.assertTrue((out / "analysis_config.yaml").is_file())


if __name__ == "__main__":
    unittest.main()
