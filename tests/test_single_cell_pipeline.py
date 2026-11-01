import gzip
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import io as scipy_io
from scipy import sparse

from funmirbench.single_cell_pipeline import experiments_config, load_samples, pseudobulk, symbol_map_from_gtf

GENES = [("ENSG00000000001.5", "GENEA"), ("ENSG00000000002", "GENEB"), ("ENSG00000198804", "MT-CO1")]
GENES += [(f"ENSG{n:011d}", f"G{n}") for n in range(10, 40)] + [("EGFP_IR_Puro", "EGFP_IR_Puro")]


def write_sample(root: Path, name: str, counts: np.ndarray) -> dict:
    """10x files of one sample; counts is genes x cells for the first three genes, the others (30 genes and a
    reporter transgene) have no counts."""
    counts = np.vstack([counts, np.zeros((len(GENES) - counts.shape[0], counts.shape[1]), dtype=int)])
    files = {key: root / f"{name}_{key}.{ext}" for key, ext in
             (("matrix", "mtx.gz"), ("features", "tsv.gz"), ("barcodes", "tsv.gz"))}
    with gzip.open(files["matrix"], "wb") as handle:
        scipy_io.mmwrite(handle, sparse.coo_matrix(counts))
    with gzip.open(files["features"], "wt") as handle:
        handle.writelines(f"{gene_id}\t{symbol}\tGene Expression\n" for gene_id, symbol in GENES)
    with gzip.open(files["barcodes"], "wt") as handle:
        handle.writelines(f"{name}-cell{i}\n" for i in range(counts.shape[1]))
    return {key: str(path) for key, path in files.items()}


class SingleCellPseudobulkTests(unittest.TestCase):
    def test_passing_cells_are_summed_per_sample(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            # cells: kept, kept, too few genes (one detected), mitochondrial fraction 0.5
            control = write_sample(root, "ctrl", np.array([[5, 2, 9, 1], [3, 4, 0, 0], [0, 1, 0, 1]]))
            treated = write_sample(root, "mimic", np.array([[1, 7], [6, 2], [1, 0]]))
            config = {"dataset_id": "TEST", "samples": [
                {"sample_id": "ctrl_d2", "group": "control", "batch": "d2", **control},
                {"sample_id": "mimic_d2", "group": "treated", "batch": "d2", **treated}]}
            samples = load_samples(config, config_path=root / "config.yaml", repo=root, force=False)

            matrix, cells, summary, _ = pseudobulk(samples, {"min_genes": 2, "max_mito_fraction": 0.2}, "MT-",
                                                   min_cells=2)

            first = ["ENSG00000000001", "ENSG00000000002", "ENSG00000198804"]
            self.assertEqual(matrix.shape[0], len(GENES) - 1)
            self.assertNotIn("EGFP_IR_Puro", matrix.index)
            self.assertEqual(matrix.loc[first, "ctrl_d2"].tolist(), [7, 7, 1])
            self.assertEqual(matrix.loc[first, "mimic_d2"].tolist(), [8, 8, 1])
            self.assertEqual(summary.features_left_out.tolist(), [1, 1])
            self.assertEqual(summary.features_left_out_examples.tolist(), ["EGFP_IR_Puro", "EGFP_IR_Puro"])
            self.assertEqual(cells.passed.tolist(), [True, True, False, False, True, True])
            self.assertEqual(summary.cells_passed.tolist(), [2, 2])

    def test_hashed_arms_of_one_lane_become_separate_samples(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            lane = write_sample(root, "lane", np.array([[5, 2, 9, 4], [3, 4, 1, 2], [0, 1, 0, 0]]))
            calls = root / "calls.csv"
            calls.write_text("Barcode,Call\nlane-cell0,A955_x\nlane-cell1,A965_x\nlane-cell2,A955_x\n"
                             "lane-cell3,doublet_x\n")
            def arm(tag):
                return {"path": str(calls), "barcode_column": "Barcode", "label_column": "Call", "label": tag,
                        "label_regex": "^([^_]+)_"}
            config = {"dataset_id": "TEST", "samples": [
                {"sample_id": "saline", "group": "control", "cell_labels": arm("A955"), **lane},
                {"sample_id": "mimic", "group": "treated", "cell_labels": arm("A965"), **lane}]}
            samples = load_samples(config, config_path=root / "config.yaml", repo=root, force=False)

            matrix, cells, summary, _ = pseudobulk(samples, {}, "MT-", min_cells=1)

            first = ["ENSG00000000001", "ENSG00000000002", "ENSG00000198804"]
            self.assertEqual(matrix.loc[first, "saline"].tolist(), [14, 4, 0])
            self.assertEqual(matrix.loc[first, "mimic"].tolist(), [2, 4, 1])
            self.assertEqual(summary.cells_selected.tolist(), [2, 1])
            self.assertEqual(cells.barcode.tolist(), ["lane-cell0", "lane-cell2", "lane-cell1"])

    def test_gene_symbols_map_to_genes_with_a_unique_symbol(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            gtf = Path(temporary_directory) / "genes.gtf"
            record = '1\tensembl\tgene\t1\t9\t.\t+\t.\tgene_id "{}"; gene_name "{}";\n'
            gtf.write_text(record.format("ENSG00000000001", "GENEA") + record.format("ENSG00000000002", "DUP")
                           + record.format("ENSG00000000003", "DUP"))

            self.assertEqual(symbol_map_from_gtf(gtf), {"GENEA": "ENSG00000000001"})

    def test_genes_missing_from_one_annotation_are_left_out(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            control = write_sample(root, "ctrl", np.array([[5, 2], [3, 4], [0, 1]]))
            treated = write_sample(root, "mimic", np.array([[1, 7], [6, 2], [1, 0]]))
            lines = gzip.open(treated["features"], "rt").readlines()
            with gzip.open(treated["features"], "wt") as handle:  # GENEB under another ID in the treated file
                handle.writelines([lines[0], lines[1].replace("ENSG00000000002", "ENSG00000000099"), *lines[2:]])
            config = {"dataset_id": "TEST", "samples": [
                {"sample_id": "ctrl", "group": "control", **control},
                {"sample_id": "mimic", "group": "treated", **treated}]}
            samples = load_samples(config, config_path=root / "config.yaml", repo=root, force=False)

            matrix, _, summary, _ = pseudobulk(samples, {}, "MT-", min_cells=1)

            self.assertNotIn("ENSG00000000002", matrix.index)
            self.assertNotIn("ENSG00000000099", matrix.index)
            self.assertEqual(summary.genes_not_in_every_file.tolist(), [1, 1])

    def test_too_few_cells_is_an_error(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            control = write_sample(root, "ctrl", np.array([[5, 0], [3, 0], [0, 4]]))
            treated = write_sample(root, "mimic", np.array([[1, 7], [6, 2], [1, 0]]))
            config = {"dataset_id": "TEST", "samples": [
                {"sample_id": "ctrl", "group": "control", **control},
                {"sample_id": "mimic", "group": "treated", **treated}]}
            samples = load_samples(config, config_path=root / "config.yaml", repo=root, force=False)

            with self.assertRaisesRegex(ValueError, "ctrl: 1 cells pass QC"):
                pseudobulk(samples, {"max_mito_fraction": 0.2}, "MT-", min_cells=2)

    def test_pseudo_replicates_split_the_passing_cells_of_each_sample(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            control = write_sample(root, "ctrl", np.array([[5, 2, 9, 1], [3, 4, 2, 6], [0, 1, 0, 0]]))
            treated = write_sample(root, "mimic", np.array([[1, 7, 3], [6, 2, 2], [1, 0, 0]]))
            config = {"dataset_id": "TEST", "samples": [
                {"sample_id": "ctrl", "group": "control", **control},
                {"sample_id": "mimic", "group": "treated", **treated}]}
            samples = load_samples(config, config_path=root / "config.yaml", repo=root, force=False)

            matrix, cells, _, design = pseudobulk(samples, {}, "MT-", min_cells=1, pseudo_replicates=2)

            self.assertEqual(design.column.tolist(), ["ctrl_pr1", "ctrl_pr2", "mimic_pr1", "mimic_pr2"])
            first = ["ENSG00000000001", "ENSG00000000002", "ENSG00000198804"]
            self.assertEqual((matrix.ctrl_pr1 + matrix.ctrl_pr2).loc[first].tolist(), [17, 15, 1])
            self.assertEqual(sorted(cells.pseudo_replicate[cells.sample_id == "ctrl"].tolist()), [1, 1, 2, 2])

    def test_batches_reach_the_deseq2_config_only_when_they_differ(self):
        samples = pd.DataFrame([{"column": "c2", "group": "control", "batch": "d2"},
                                {"column": "t2", "group": "treated", "batch": "d2"},
                                {"column": "c3", "group": "control", "batch": "d3"},
                                {"column": "t3", "group": "treated", "batch": "d3"}])
        config = {"dataset_id": "TEST", "mirna_name": "hsa-miR-98-5p", "experiment_type": "Overexpression"}

        paired = experiments_config(config, samples, Path("/x/counts.tsv.gz"))
        unpaired = experiments_config(config, samples.assign(batch=""), Path("/x/counts.tsv.gz"))

        self.assertEqual(paired["comparison"]["control_columns"], ["c2", "c3"])
        self.assertEqual(paired["comparison"]["batch_by_sample"], {"c2": "d2", "t2": "d2", "c3": "d3", "t3": "d3"})
        self.assertNotIn("batch_by_sample", unpaired["comparison"])
        self.assertEqual(yaml.safe_load(yaml.safe_dump(paired))["source"]["mode"], "count_matrix")


if __name__ == "__main__":
    unittest.main()
