import pathlib
import tempfile
import unittest

import numpy as np
import pandas as pd

from funmirbench import single_cell_eval

SELECTED = pd.DataFrame({"comparison_id": ["main_one", "pair", "screen__DICER1_KO"], "group": ["main", "supplement", "machinery"]})


def final_row(comparison_id, mirna_or_gene, dataset="DATA", status="fitted"):
    return {"comparison_id": comparison_id, "mirna_or_gene": mirna_or_gene, "orientation": "OE", "dataset": dataset, "design": "~ condition",
            "replicate_type": "technical", "final_dir": dataset, "table": f"genes/{comparison_id}.tsv.gz", "status": status}


class ReadComparisonsTests(unittest.TestCase):
    def write(self, directory, rows):
        directory = pathlib.Path(directory)
        SELECTED.to_csv(directory / "selected.tsv", sep="\t", index=False)
        pd.DataFrame(rows).to_csv(directory / "comparisons.tsv", sep="\t", index=False)
        return directory

    def test_keeps_only_selected_comparisons_with_their_group(self):
        with tempfile.TemporaryDirectory() as directory:
            rows = [final_row("main_one", "hsa-miR-1-3p"), final_row("pair", "hsa-miR-2-3p"), final_row("screen__DICER1_KO", "DICER1"), final_row("other", "hsa-miR-4-3p")]
            root = self.write(directory, rows)
            kept = single_cell_eval.read_comparisons(root, root / "selected.tsv")
        self.assertEqual(sorted(kept.comparison_id), ["main_one", "pair", "screen__DICER1_KO"])
        self.assertEqual(kept.set_index("comparison_id").group["main_one"], "main")

    def test_a_selected_comparison_missing_from_the_final_table_stops_the_run(self):
        with tempfile.TemporaryDirectory() as directory:
            root = self.write(directory, [final_row("main_one", "hsa-miR-1-3p")])
            with self.assertRaisesRegex(ValueError, "run run.sh for all groups"):
                single_cell_eval.read_comparisons(root, root / "selected.tsv")

    def test_a_selected_comparison_that_is_not_fitted_stops_the_run(self):
        with tempfile.TemporaryDirectory() as directory:
            root = self.write(directory, [final_row("main_one", "hsa-miR-1-3p", status="not_fitted"), final_row("pair", "hsa-miR-2-3p"), final_row("screen__DICER1_KO", "DICER1")])
            with self.assertRaisesRegex(ValueError, "not fitted"):
                single_cell_eval.read_comparisons(root, root / "selected.tsv")


class EvaluatorInputTests(unittest.TestCase):
    def test_one_row_per_named_mirna_and_machinery_left_out(self):
        rows = pd.DataFrame([final_row("main_one", "hsa-miR-1-3p"), final_row("pair", "hsa-miR-2-3p;hsa-miR-3-3p"), final_row("screen__DICER1_KO", "DICER1")])
        rows["group"] = ["main", "supplement", "machinery"]
        rows.loc[0, "orientation"] = "KO"
        mirnas = pd.DataFrame({"mirna_name": ["hsa-miR-1-3p"], "mirna_sequence": ["UGGAAUGUAAAGAAGUAUGUAU"], "seed_2_8": ["GGAAUGU"], "seed_family": ["a"]})
        families = pd.DataFrame({"mirna_name": ["hsa-miR-2-3p", "hsa-miR-3-3p"], "mirna_sequence": ["UAACACUGUCUGGUAAAGAUGG", "UCACAUUGCCAGGGAUUUCC"], "seed": ["AACACUG", "CACAUUG"]})
        tools = pd.DataFrame({"tool_id": ["tool"], "official_name": ["Tool"], "score_direction": ["higher_is_stronger"], "predictor_output_path": ["/tool.tsv"]})
        with tempfile.TemporaryDirectory() as directory:
            out = pathlib.Path(directory) / "input"
            comparisons = single_cell_eval.evaluator_input(rows, pathlib.Path("/final"), mirnas, families, tools, {}, out)
            written = pd.read_csv(out / "mirnas.tsv", sep="\t")
        self.assertEqual(comparisons.experiment_id.tolist(), ["main_one", "pair__hsa-miR-2-3p", "pair__hsa-miR-3-3p"])
        self.assertEqual(comparisons.perturbation.unique().tolist(), ["OE"])
        self.assertEqual(comparisons.orientation.tolist(), ["KO", "OE", "OE"])  # the evaluator's OE does not hide the real direction
        self.assertEqual(comparisons.action.tolist(), ["keep", "supplement", "supplement"])
        self.assertEqual(sorted(written.mirna_name), ["hsa-miR-1-3p", "hsa-miR-2-3p", "hsa-miR-3-3p"])


class FamilyTests(unittest.TestCase):
    def test_sets_follow_the_share_rules(self):
        expressed = [("HIGH1", 0.6), ("HIGH2", 0.2)]
        low = [(f"LOW{i}", 0.00001) for i in range(5)]
        middle = [("MIDDLE", 0.001)]
        profiles = pd.DataFrame([{"line": "A", "seed": seed, "share": share} for seed, share in expressed + low + middle])
        sets = single_cell_eval.family_sets(profiles)["A"]
        self.assertEqual(sets["expressed"], ["HIGH1", "HIGH2"])
        np.testing.assert_allclose(sets["weights"], [0.75, 0.25])
        self.assertEqual(len(sets["decoys"]), single_cell_eval.N_DECOYS)
        for decoy in sets["decoys"]:
            self.assertEqual(len(set(decoy)), 2)
            self.assertTrue(set(decoy) <= {seed for seed, _ in low})  # MIDDLE is neither expressed nor a decoy

    def test_too_few_low_families_is_an_error(self):
        profiles = pd.DataFrame({"line": "A", "seed": ["HIGH1", "HIGH2", "LOW"], "share": [0.5, 0.4, 0.00001]})
        with self.assertRaisesRegex(ValueError, "low families"):
            single_cell_eval.family_sets(profiles)

    def test_a_family_with_no_scored_mirna_does_not_count(self):
        ranks = {"a1": np.array([1.0, 0.0]), "a2": np.array([0.0, 0.0]), "b1": np.array([0.0, 1.0])}
        members = {"A": ["a1", "a2"], "B": ["b1"], "C": ["unscored"]}
        combined = single_cell_eval.combine_families(["A", "B", "C"], np.array([0.5, 0.25, 0.25]), members, ranks, 2)
        np.testing.assert_allclose(combined, [(0.5 * 0.5) / 0.75, 0.25 / 0.75])


class FigureMainTests(unittest.TestCase):
    def test_only_main_comparisons_whose_seed_effect_interval_is_above_0_are_drawn(self):
        rows = pd.DataFrame({"comparison_id": ["acted", "zero_in_interval", "no_seed_row", "machinery_row"], "group": ["main", "main", "main", "machinery"]})
        seed = pd.DataFrame({"comparison_id": ["acted", "zero_in_interval", "machinery_row"], "ci_low": [0.04, -0.01, 0.5]})
        self.assertEqual(single_cell_eval.figure_main(rows, seed).comparison_id.tolist(), ["acted"])


class FigureMachineryTests(unittest.TestCase):
    def test_only_figure_screens_with_a_clear_transcript_drop_and_a_passed_seed_check_are_drawn(self):
        rows = pd.DataFrame({"comparison_id": ["clear", "unclear", "other_screen", "seed_failed", "main_row"], "group": ["machinery"] * 4 + ["main"],
                             "dataset": ["Replogle2022_K562_GW", "Replogle2022_K562_GW", "XAtlasOrion_HCT116", "Replogle2022_K562_GW", "DATA"]})
        drop = pd.DataFrame({"comparison_id": ["clear", "unclear", "other_screen", "seed_failed"], "target_ci95_high": [-0.5, 0.1, -2.0, -1.0]})
        machinery_seed = pd.DataFrame({"comparison_id": ["clear", "unclear", "other_screen", "seed_failed"], "passes": [True, True, True, False]})
        self.assertEqual(single_cell_eval.figure_machinery(rows, drop, machinery_seed).comparison_id.tolist(), ["clear"])

    def test_a_gene_median_pools_its_screens_and_keeps_gene_types_apart(self):
        record = lambda comparison, screen, gene, gene_type, value: {"comparison_id": comparison, "screen": screen, "gene": gene, "gene_type": gene_type,
                                                                     "predictor": "tool", "auroc_d05": np.nan, "spearman_shrunk": value, "decoy_z_shrunk": value}
        machinery = pd.DataFrame([record("a", "Replogle2022_K562_GW", "XPO5", "miRNA-making", 1.0), record("b", "Nadig2025_HepG2", "XPO5", "miRNA-making", 3.0),
                                  record("c", "Replogle2022_K562_GW", "DICER1", "miRNA-making", 4.0), record("d", "Replogle2022_K562_GW", "AGO2", "repression", 5.0),
                                  record("left_out", "Nadig2025_HepG2", "AGO2", "repression", 100.0)])
        kept = pd.DataFrame({"comparison_id": ["a", "b", "c", "d"]})
        tables, n_making = single_cell_eval.machinery_tables(machinery, kept, pd.DataFrame({"tool_id": ["tool"]}))
        self.assertEqual(n_making, 2)
        self.assertEqual(tables["spearman_shrunk"].tool.to_dict(), {"DICER1": 4.0, "XPO5": 2.0, "AGO2": 5.0})
        self.assertNotIn("auroc_d05", tables)

    def test_two_experiments_of_one_gene_in_one_screen_are_refused(self):
        record = lambda comparison: {"comparison_id": comparison, "screen": "Replogle2022_K562_GW", "gene": "XPO5", "gene_type": "miRNA-making",
                                     "predictor": "tool", "spearman_shrunk": 1.0, "decoy_z_shrunk": 1.0}
        with self.assertRaises(ValueError):
            single_cell_eval.machinery_tables(pd.DataFrame([record("a"), record("b")]), pd.DataFrame({"comparison_id": ["a", "b"]}), pd.DataFrame({"tool_id": ["tool"]}))

if __name__ == "__main__":
    unittest.main()
