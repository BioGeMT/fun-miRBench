import json
import pathlib
import tempfile
import unittest
from unittest.mock import patch

import pandas as pd

from funmirbench import add_experiment as add

CONFIG = {"dataset_id": "my_experiment", "mirna_name": "hsa-miR-new-5p", "experiment_type": "Overexpression",
          "source": {"mode": "count_matrix"}, "metadata": {"tested_cell_line": "HeLa"},
          "comparison": {"control_columns": ["c1", "c2"], "treated_columns": ["t1", "t2"]}}


class AddExperimentTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = pathlib.Path(tmp.name)
        self.data = self.root / "benchmark_data"
        self.data.mkdir()
        pd.DataFrame({"experiment_id": ["E1"], "n_eligible": [10]}).to_csv(self.data / "experiments.tsv", sep="\t", index=False)
        # a published column the new row does not have stays empty for it
        pd.DataFrame({"experiment_id": ["E1"], "flags": [""], "gc_rho": [0.1]}).to_csv(self.data / "comparisons.tsv", sep="\t", index=False)
        (self.root / "metadata").mkdir()
        pd.DataFrame({"mirna_name": ["hsa-miR-16-5p"], "mirna_sequence": ["UAGCAGCACGUAAAUAUUGGCG"], "seed_2_8": ["AGCAGCA"],
                      "seed_family": ["fam1"]}).to_csv(self.root / add.MIRNAS, sep="\t", index=False)
        fasta = self.root / add.MATURE_FASTA
        fasta.parent.mkdir(parents=True)
        fasta.write_text(">hsa-miR-new-5p MIMAT1 Homo sapiens miR-new-5p\nUGGAAUGUAAAGAAGUAUGUAU\n")
        (self.root / add.GENES).touch()
        (self.data / "manifest.json").write_text(json.dumps({"experiments": 1, "tables": {"E1": "x"},
                                                             "size_factors": add.SIZE_FACTORS}))

    def run_add(self, config=CONFIG, eligible=7):
        def label(row, *, genes, out, rscript):
            (out / "genes" / f"{row.experiment_id}.tsv.gz").write_text("table")
            return {"experiment_id": row.experiment_id, "n_eligible": eligible, "table": f"genes/{row.experiment_id}.tsv.gz",
                    "table_sha256": "abc", "d05_responsive": 1, "d05_no_response_above_margin": 2, "d05_unresolved": 4}

        path = self.root / "config.yaml"
        with patch.object(add, "load_yaml", return_value=config), \
                patch.object(add, "run_ingestion_config", return_value={"run_dir": str(self.root)}) as self.ingest, \
                patch.object(add, "run_one", side_effect=label), patch.object(add.shutil, "which", return_value="Rscript"):
            return add.add_experiment(path, self.data, self.root)

    def test_new_experiment_gets_its_rows_and_mirna(self):
        self.assertTrue(self.run_add().is_file())
        comps = pd.read_csv(self.data / "comparisons.tsv", sep="\t", keep_default_na=False).set_index("experiment_id")
        self.assertEqual(comps.loc["my_experiment", "cell_line"], "HeLa")
        self.assertEqual(comps.loc["my_experiment", "source_mode"], "count_matrix")
        self.assertEqual(comps.loc["my_experiment", "action"], "keep")
        self.assertEqual(str(comps.loc["my_experiment", "gc_rho"]), "")
        self.assertEqual(float(comps.loc["E1", "gc_rho"]), 0.1)
        self.assertEqual(list(pd.read_csv(self.data / "experiments.tsv", sep="\t").experiment_id), ["E1", "my_experiment"])
        mirnas = pd.read_csv(self.root / add.MIRNAS, sep="\t").set_index("mirna_name")
        self.assertEqual(mirnas.loc["hsa-miR-new-5p", "seed_2_8"], "GGAAUGU")

    def test_an_id_already_in_the_tables_is_refused(self):
        with self.assertRaisesRegex(ValueError, "already in"):
            self.run_add({**CONFIG, "dataset_id": "E1"})

    def test_tables_built_with_other_size_factors_are_refused(self):
        manifest = json.loads((self.data / "manifest.json").read_text())
        for value in (None, "median of ratios to all samples"):
            with self.subTest(value):
                (self.data / "manifest.json").write_text(json.dumps({**manifest, "size_factors": value}))
                with self.assertRaisesRegex(ValueError, "not built with"):
                    self.run_add()
                self.ingest.assert_not_called()
                self.assertEqual(list(pd.read_csv(self.data / "experiments.tsv", sep="\t").experiment_id), ["E1"])
        with self.subTest("interrupted first addition made the old way"):
            (self.data / "comparisons.tsv").unlink()
            (self.data / "manifest.json").write_text(json.dumps({**manifest, "size_factors": None}))
            with self.assertRaisesRegex(ValueError, "not built with"):
                self.run_add()
            self.ingest.assert_not_called()

    def test_counts_without_a_usable_gene_add_nothing(self):
        with self.assertRaisesRegex(ValueError, "Ensembl gene IDs"):
            self.run_add(eligible=0)
        self.assertEqual(list(pd.read_csv(self.data / "experiments.tsv", sep="\t").experiment_id), ["E1"])
        self.assertFalse((self.data / "genes/my_experiment.tsv.gz").exists())

    def fail_replacing(self, name):
        real = add.os.replace

        def fail(src, dst):
            if pathlib.Path(dst).name == name:
                raise OSError("disk full")
            real(src, dst)

        return patch.object(add.os, "replace", side_effect=fail)

    def test_a_failed_write_can_be_retried_without_a_duplicate(self):
        for name in ("metadata", "experiments.tsv", "manifest.json", "comparisons.tsv"):
            target = {"metadata": add.MIRNAS.name}.get(name, name)
            with self.subTest(name), self.fail_replacing(target), self.assertRaises(OSError):
                self.run_add()
            self.assertEqual(list(pd.read_csv(self.data / "comparisons.tsv", sep="\t").experiment_id), ["E1"])
            self.assertEqual(list(pd.read_csv(self.root / add.MIRNAS, sep="\t").mirna_name), ["hsa-miR-16-5p"])
        json.loads((self.data / "manifest.json").read_text())  # still valid after the failed writes
        self.run_add()
        self.assertEqual(list(pd.read_csv(self.data / "experiments.tsv", sep="\t").experiment_id), ["E1", "my_experiment"])
        self.assertEqual(list(pd.read_csv(self.data / "comparisons.tsv", sep="\t").experiment_id), ["E1", "my_experiment"])
        self.assertEqual(sorted(json.loads((self.data / "manifest.json").read_text())["tables"]), ["E1", "my_experiment"])

    def test_a_failed_first_addition_to_a_new_folder_can_be_retried(self):
        for name in ("experiments.tsv", "comparisons.tsv", "manifest.json"):
            (self.data / name).unlink()
        with self.fail_replacing("manifest.json"), self.assertRaises(OSError):
            self.run_add()
        self.assertFalse((self.data / "comparisons.tsv").exists())
        self.run_add()
        self.assertEqual(list(pd.read_csv(self.data / "comparisons.tsv", sep="\t").experiment_id), ["my_experiment"])
        self.assertEqual(list(pd.read_csv(self.data / "experiments.tsv", sep="\t").experiment_id), ["my_experiment"])

    def test_without_published_tables_a_new_folder_is_started(self):
        for name in ("experiments.tsv", "comparisons.tsv", "manifest.json"):
            (self.data / name).unlink()
        self.data.rmdir()
        self.run_add()
        self.assertEqual(list(pd.read_csv(self.data / "comparisons.tsv", sep="\t").experiment_id), ["my_experiment"])
        manifest = json.loads((self.data / "manifest.json").read_text())
        self.assertEqual(manifest["tables"], {"my_experiment": "abc"})
        self.assertEqual(manifest["comparisons_sha256"], add.sha256(self.data / "comparisons.tsv"))

    def test_only_oe_with_replicates_is_accepted(self):
        with self.assertRaisesRegex(ValueError, "overexpression"):
            add.check_config({**CONFIG, "experiment_type": "Knockout"})
        with self.assertRaisesRegex(ValueError, "at least 2 control"):
            add.check_config({**CONFIG, "comparison": {"control_columns": ["c1"], "treated_columns": ["t1", "t2"]}})


if __name__ == "__main__":
    unittest.main()
