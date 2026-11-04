import pathlib
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd
import yaml

from funmirbench.response_vectors import (
    alignment, build_vectors, control_expression_tenths, create_output_dir, read_response_table,
    run_analysis, specificity,
)
from funmirbench.benchmark import load_predictions


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

    def test_tied_expression_does_not_make_alignment_depend_on_gene_order(self):
        n = 400
        genes = pd.Index([f"g{i:03}" for i in range(n)])
        table = pd.DataFrame({"eligible": True, "oriented_logFC_shrunk": np.repeat([0., 1.], n // 2),
                              "control_mean_normalized_count": 1.}, index=genes)
        score = pd.Series(np.tile(np.arange(n // 2), 2), index=genes)
        shuffled = table.iloc[np.ravel(np.column_stack([np.arange(n // 2), np.arange(n // 2, n)]))]
        for data in (table, shuffled):
            bins = control_expression_tenths(data.control_mean_normalized_count)
            self.assertEqual(bins.nunique(), 1)
            result = alignment(build_vectors(data, score), bins)
            self.assertAlmostEqual(result["spearman_within_tenths"], 0.)
            self.assertEqual(result["n_within_tenths_genes"], n)

    def test_unknown_baseline_is_not_zero_and_does_not_enter_a_bin(self):
        values = pd.Series([0., 0., 1., np.nan, np.inf])
        bins = control_expression_tenths(values)
        self.assertEqual(bins[0], bins[1])
        self.assertTrue(bins.loc[[3, 4]].isna().all())
        self.assertGreater(bins[2], bins[0])
        self.assertTrue(control_expression_tenths(pd.Series([np.nan] * 3)).isna().all())
        with self.assertRaisesRegex(ValueError, "nonnegative"):
            control_expression_tenths(pd.Series([-1., 0., 1.]))

    def test_specificity_requires_the_same_bin_support(self):
        candidates = pd.DataFrame({"candidate_mirna": ["own", "d1", "d2"],
                                   "candidate_seed": ["AAAAAAA", "CCCCCCC", "GGGGGGG"],
                                   "spearman_within_tenths": [.8, .9, .2],
                                   "within_tenths_bins": ["1", "1,2", "1"]})
        result = specificity(candidates, "own", "spearman_within_tenths")
        self.assertEqual(result["n_decoys"], 1)
        self.assertEqual(result["n_decoys_excluded_support"], 1)
        self.assertEqual(result["cognate_rank"], 1.)
        candidates.loc[2, "within_tenths_bins"] = "2"
        result = specificity(candidates, "own", "spearman_within_tenths")
        self.assertTrue(np.isnan(result["cognate_rank"]))
        self.assertEqual(result["specificity_na_reason"], "no_comparable_different_seed_decoys")

    def test_too_few_or_nonfinite_responses_are_reported_as_na(self):
        table = self.table()
        table["oriented_logFC_shrunk"] = pd.Series(
            [np.inf, np.nan, 1., 1., -np.inf, 2.], index=table.index, dtype="Float64")
        vectors = build_vectors(table, pd.Series(1., index=table.index))
        result = alignment(vectors, pd.Series(1, index=table.index))
        self.assertEqual(result["n_usable"], 2)
        self.assertEqual(result["alignment_na_reason"], "fewer_than_three_usable_genes")
        self.assertTrue(np.isnan(result["spearman"]))

    def test_prediction_loader_discards_nonfinite_before_aggregating(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = pathlib.Path(tmp) / "scores.tsv.gz"
            pd.DataFrame({"Ensembl_ID": ["a.1", "a.2", "a.3", "b", "c", "d"],
                          "miRNA_Name": "m", "Score": [1., np.inf, -np.inf, np.nan, 2., "bad"]}).to_csv(
                              path, sep="\t", index=False)
            for direction, expected in (("higher_is_stronger", 1.), ("lower_is_stronger", -1.)):
                loaded = load_predictions(path, direction, {"m"}, {"a", "b", "c", "d"})["m"]
                self.assertEqual(set(loaded.index), {"a", "c"})
                self.assertEqual(loaded["a"], expected)

    def test_versioned_response_ids_normalize_and_collisions_fail(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = pathlib.Path(tmp) / "response.tsv"
            table = self.table().reset_index()
            table["gene_id"] += ".1"
            table.to_csv(path, sep="\t", index=False)
            loaded = read_response_table(path)
            self.assertEqual(list(loaded.index), list(self.table().index))
            table.loc[1, "gene_id"] = "a.2"
            table.to_csv(path, sep="\t", index=False)
            with self.assertRaisesRegex(ValueError, "unique"):
                read_response_table(path)
            table.loc[1, "gene_id"] = "b.1"
            table.loc[0, "gene_id"] = ".1"
            table.to_csv(path, sep="\t", index=False)
            with self.assertRaisesRegex(ValueError, "gene_id"):
                read_response_table(path)

    def test_same_second_runs_get_distinct_directories(self):
        with tempfile.TemporaryDirectory() as tmp, patch(
                "funmirbench.response_vectors.build_run_dir_name", return_value="20261102_120000"):
            first = create_output_dir(pathlib.Path(tmp))
            marker = first / "keep.txt"
            marker.write_text("original")
            second = create_output_dir(pathlib.Path(tmp))
            self.assertNotEqual(first, second)
            self.assertEqual(marker.read_text(), "original")


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

    def test_reordered_experiments_and_versioned_gene_ids_align_to_one_universe(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            genes = [f"ENSG{i:011d}" for i in range(300)]
            response = np.linspace(-1, 1, len(genes))
            for name in ("e1", "e2"):
                table = pd.DataFrame({"gene_id": [g + ".3" for g in genes], "eligible": True,
                                      "oriented_logFC_shrunk": response,
                                      "control_mean_normalized_count": 0.})
                if name == "e2":
                    table = table.iloc[::-1]
                table.to_csv(root / f"{name}.tsv", sep="\t", index=False)
            pd.DataFrame({"experiment_id": ["e1", "e2"], "table": ["e1.tsv", "e2.tsv"]}).to_csv(
                root / "experiments.tsv", sep="\t", index=False)
            pd.DataFrame({"experiment_id": ["e1", "e2"], "mirna_name": "own", "perturbation": "OE",
                          "action": "keep", "dependency_group": "g"}).to_csv(root / "comparisons.tsv", sep="\t", index=False)
            pd.DataFrame({"mirna_name": ["own"], "mirna_sequence": ["  TAAAAAAA  "]}).to_csv(
                root / "mirnas.tsv", sep="\t", index=False)
            pd.DataFrame({"Ensembl_ID": genes, "miRNA_Name": "own", "Score": response}).to_csv(
                root / "scores.tsv", sep="\t", index=False)
            pd.DataFrame({"tool_id": ["tool"], "score_direction": ["higher_is_stronger"],
                          "predictor_output_path": ["scores.tsv"]}).to_csv(root / "predictors.tsv", sep="\t", index=False)
            config = root / "config.yaml"
            config.write_text(yaml.safe_dump({"input_dir": ".", "mirnas_tsv": "mirnas.tsv",
                                              "predictions_tsv": "predictors.tsv", "out_dir": "results"}))
            out = run_analysis(config)
            candidates = pd.read_csv(out / "candidate_alignment.tsv", sep="\t")
            np.testing.assert_allclose(candidates.spearman, 1.)
            np.testing.assert_allclose(candidates.spearman_within_tenths, 1.)
            files = pd.read_csv(out / "vector_files.tsv", sep="\t")
            for path in files.file:
                self.assertEqual(pd.read_csv(out / path, sep="\t").gene_id.tolist(), genes)
            # A genuinely different universe must still fail before creating outputs.
            table = pd.read_csv(root / "e2.tsv", sep="\t").iloc[:-1]
            table.to_csv(root / "e2.tsv", sep="\t", index=False)
            with self.assertRaisesRegex(ValueError, "same gene list"):
                run_analysis(config)
            self.assertEqual(len(list((root / "results").glob("*/response_vectors"))), 1)

    def test_empty_yaml_has_an_actionable_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = pathlib.Path(tmp) / "empty.yaml"
            path.write_text("")
            with self.assertRaisesRegex(ValueError, "YAML mapping"):
                run_analysis(path)


if __name__ == "__main__":
    unittest.main()
