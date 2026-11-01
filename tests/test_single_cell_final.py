import gzip
import os
import shutil
import stat
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse, stats

from funmirbench.single_cell_cooccurrence import run_cooccurrence
from funmirbench.single_cell_final import CellSums, ratio_table, read_comparison_list, select_units
from funmirbench.single_cell_relabel import check_only_labels_changed, run_relabel
from funmirbench.single_cell_store import (
    assign_lane_sets,
    map_features_to_gene_ids,
    pseudo_replicate_units,
    raw_counts_from_lognorm,
    read_10x_h5,
    read_cell_store,
    write_cell_store,
)

R_ENV = Path("/SCRATCH/dtzim01/fun-experiment-data/tools/envs/deseq2-shrinkage/bin")
REPO = Path(__file__).resolve().parents[1]
N_GENES = 1000
GENE_IDS = [f"ENSG{n:011d}" for n in range(1, N_GENES + 1)]
COMPARISON_DEFAULTS = {"dataset": "TEST", "category": "raise", "mirna_or_gene": "hsa-miR-1-3p", "orientation": "OE",
                       "treated_arm": "treated", "control_arm": "ctrl", "min_cells": "5", "caveats": ""}


def write_gtf(root: Path) -> Path:
    gtf = root / "genes.gtf"
    gtf.write_text("".join(f'1\tensembl\tgene\t1\t9\t.\t+\t.\tgene_id "{gene}"; gene_version "1"; gene_name "G{i}";\n'
                           for i, gene in enumerate(GENE_IDS)))
    return gtf


def obs_table(arms, batches, units, *, replicate_type="biological", guides=None, passed=None) -> pd.DataFrame:
    arms = list(arms)
    is_control = np.array([arm in ("ctrl", "no_guide") for arm in arms])
    return pd.DataFrame({
        "dataset": "TEST", "sample": list(units), "batch": list(batches), "replicate_unit": list(units),
        "replicate_type": replicate_type, "arm": arms, "is_control": is_control,
        "control_type": ["non_targeting" if control else "" for control in is_control],
        "guide": arms, "n_guides": np.ones(len(arms), dtype=int) if guides is None else np.asarray(guides),
        "passed_qc": np.ones(len(arms), dtype=bool) if passed is None else np.asarray(passed),
    }, index=[f"cell{i}" for i in range(len(arms))])


def write_store(root: Path, counts, obs: pd.DataFrame, name: str = "store.h5ad") -> Path:
    path = root / name
    write_cell_store(path, counts, obs, pd.Series(GENE_IDS), gtf=root / "genes.gtf", source="test")
    return path


def comparison_row(**values):
    return pd.Series({**COMPARISON_DEFAULTS, "comparison_id": "c", "other_pathway_arms": "", **values})


def simulate_cells(fold_change: np.ndarray, *, cells_per_unit: int = 60, seed: int = 1):
    """Three donors, each one batch with control and treated cells: Poisson counts around a shared gene profile with
    a depth per cell and a small factor per donor and arm. Treated cells have the means times fold_change."""
    rng = np.random.default_rng(seed)
    means = rng.lognormal(mean=0.2, sigma=1.2, size=N_GENES)
    rows, arms, donors = [], [], []
    for donor in ("d1", "d2", "d3"):
        for arm, change in (("ctrl", 1.0), ("treated", fold_change)):
            unit_means = means * change * rng.lognormal(0, 0.05, N_GENES)
            depth = rng.lognormal(0, 0.3, (cells_per_unit, 1))
            rows.append(rng.poisson(unit_means * depth))
            arms += [arm] * cells_per_unit
            donors += [donor] * cells_per_unit
    return sparse.csr_matrix(np.vstack(rows)), arms, donors


class CellStoreTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.gtf = write_gtf(self.root)

    def tearDown(self):
        self.temporary.cleanup()

    def test_store_round_trip_keeps_counts_and_maps_versioned_and_unknown_features(self):
        features = pd.Series([f"{GENE_IDS[0]}.3", *GENE_IDS[1:-1], "ENSG99999999999"])
        counts = sparse.random(4, N_GENES, density=0.1, random_state=0, dtype=np.int32, data_rvs=lambda n: np.arange(1, n + 1))
        obs = obs_table(["ctrl", "ctrl", "treated", "treated"], "bbbb", "uvuv")

        write_cell_store(self.root / "s.h5ad", counts, obs, features, gtf=self.gtf, source="authors")
        store = read_cell_store(self.root / "s.h5ad")

        self.assertEqual((store.X != counts.tocsr()).nnz, 0)
        self.assertEqual(store.var.gene_id.iloc[0], GENE_IDS[0])
        self.assertEqual(store.var.gene_id.iloc[-1], "")  # not in the GTF
        self.assertEqual(store.var_names[0], f"{GENE_IDS[0]}.3")  # index stays as deposited
        self.assertEqual(store.uns["source"], "authors")

    def test_symbols_that_two_features_share_are_left_without_a_gene_id(self):
        mapped = map_features_to_gene_ids(pd.Series(["G0", "G1", "G1", "G2", "unknown"]), self.gtf)
        self.assertEqual(mapped.tolist(), [GENE_IDS[0], "", "", GENE_IDS[2], ""])

    def test_fractional_counts_are_rejected(self):
        counts = sparse.csr_matrix(np.full((2, N_GENES), 1.5))
        with self.assertRaisesRegex(ValueError, "raw integers"):
            write_store(self.root, counts, obs_table(["ctrl", "treated"], "bb", "uv"))

    def test_store_without_a_required_column_is_rejected(self):
        counts = sparse.csr_matrix(np.ones((2, N_GENES), dtype=int))
        with self.assertRaisesRegex(ValueError, "missing columns.*replicate_unit"):
            write_store(self.root, counts, obs_table(["ctrl", "treated"], "bb", "uv").drop(columns="replicate_unit"))

    def test_lognormalised_values_give_back_the_raw_counts(self):
        rng = np.random.default_rng(0)
        raw = sparse.csr_matrix(rng.poisson(2, (50, 200)) * rng.integers(1, 40, (50, 1)))
        library = np.asarray(raw.sum(axis=1)).ravel()
        normalised = raw.astype(np.float64).multiply(1e4 / library[:, None]).tocsr()
        normalised.data = np.log1p(normalised.data).astype(np.float32)  # HKU files hold float32

        rebuilt = raw_counts_from_lognorm(normalised, library)

        self.assertEqual((rebuilt != raw).nnz, 0)
        with self.assertRaisesRegex(ValueError, "not integer"):
            raw_counts_from_lognorm(normalised, library * 1.37)

    def test_pseudo_replicates_are_balanced_reproducible_and_skip_failing_cells(self):
        samples = np.array(["s1"] * 30 + ["s2"] * 30, dtype=object)
        arms = np.array(["a", "b"] * 30, dtype=object)
        passed = np.ones(60, dtype=bool)
        passed[0] = False

        units = pseudo_replicate_units(samples, arms, passed, seed=7)

        self.assertEqual(units[0], "s1")
        self.assertEqual(pd.Series(units[(samples == "s2") & (arms == "a")]).value_counts().tolist(), [5, 5, 5])
        self.assertTrue((units == pseudo_replicate_units(samples, arms, passed, seed=7)).all())
        self.assertFalse((units == pseudo_replicate_units(samples, arms, passed, seed=8)).all())

    def test_last_incomplete_lane_set_joins_the_one_before(self):
        sets = assign_lane_sets(np.arange(1, 21))
        self.assertEqual(pd.Series(sets).value_counts().sort_index().tolist(), [8, 12])
        self.assertEqual(sets[7], "lane_set_1")
        self.assertEqual(sets[8], "lane_set_2")
        with self.assertRaisesRegex(ValueError, "integers"):
            assign_lane_sets(["lane1", "lane2"])

    def test_10x_h5_of_both_layouts(self):
        matrix = sparse.csc_matrix(np.array([[1, 0], [0, 2], [3, 4]]))
        fields = {"data": matrix.data, "indices": matrix.indices, "indptr": matrix.indptr, "shape": matrix.shape,
                  "barcodes": [b"A-1", b"B-1"]}
        with h5py.File(self.root / "v3.h5", "w") as handle:
            group = handle.create_group("matrix")
            for key, value in fields.items():
                group[key] = value
            group["features/id"] = [b"ENSG00000000001", b"ENSG00000000002", b"HTO1"]
            group["features/name"] = [b"A", b"B", b"HTO1"]
            group["features/feature_type"] = [b"Gene Expression", b"Gene Expression", b"Antibody Capture"]
        with h5py.File(self.root / "v2.h5", "w") as handle:
            group = handle.create_group("GRCh38")
            for key, value in fields.items():
                group[key] = value
            group["genes"] = [b"ENSG00000000001", b"ENSG00000000002", b"ENSG00000000003"]
            group["gene_names"] = [b"A", b"B", b"C"]

        counts, features, barcodes = read_10x_h5(self.root / "v3.h5")
        self.assertEqual(counts.toarray().tolist(), [[1, 0], [0, 2]])  # the antibody row is dropped
        self.assertEqual(features[0].tolist(), ["ENSG00000000001", "ENSG00000000002"])
        self.assertEqual(barcodes.tolist(), ["A-1", "B-1"])
        counts, features, _ = read_10x_h5(self.root / "v2.h5")
        self.assertEqual(counts.shape, (3, 2))
        self.assertEqual(features[1].tolist(), ["A", "B", "C"])


class UnitTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        write_gtf(self.root)

    def tearDown(self):
        self.temporary.cleanup()

    def sums(self, counts, obs) -> CellSums:
        return CellSums(read_cell_store(write_store(self.root, counts, obs)))

    def fractional_sums(self, counts, obs) -> CellSums:
        path = self.root / "fractional.h5ad"
        write_cell_store(path, sparse.csr_matrix(counts), obs, pd.Series(GENE_IDS), gtf=self.root / "genes.gtf",
                         source="test", fractional=True)
        return CellSums(read_cell_store(path))

    def test_fractional_counts_are_rounded_after_summing(self):
        counts = np.zeros((4, N_GENES))
        counts[:, 0] = [0.4, 0.4, 0.4, 0.5]
        counts[:, 1] = [0.2, 0.2, 0.0, 0.0]
        obs = obs_table(["treated", "treated", "treated", "ctrl"], "aaaa", ["u1", "u1", "u1", "u1"])

        treated, control = self.fractional_sums(counts, obs).sides(comparison_row())

        self.assertTrue(np.issubdtype(treated.sums.dtype, np.integer))
        self.assertEqual(treated.sums.toarray()[0, :2].tolist(), [1, 0])
        self.assertEqual(control.sums.toarray()[0, :2].tolist(), [0, 0])

    def test_fractional_counts_are_rounded_after_removing_guide_carriers_at_high_guide_load(self):
        counts = np.zeros((4, N_GENES))
        counts[:, 0] = [0.4, 0.4, 0.4, 0.4]
        obs = obs_table(["treated;other", "other;x", "other;y", "x;y"], "aaaa", ["u1", "u1", "u1", "u1"], guides=[2, 2, 2, 2])

        treated, control = self.fractional_sums(counts, obs).sides(comparison_row())

        self.assertEqual(treated.sums.toarray()[0, 0], 0)  # 0.4
        self.assertEqual(control.sums.toarray()[0, 0], 1)  # 0.4 x 3, not round(1.6) - round(0.4) = 2

    def test_passing_cells_are_summed_per_unit_and_arm_within_the_batch(self):
        counts = np.zeros((6, N_GENES), dtype=int)
        counts[:, 0] = [1, 2, 4, 8, 16, 32]
        obs = obs_table(["ctrl", "ctrl", "treated", "treated", "ctrl", "treated"], "aaaabb", ["u1", "u1", "u1", "u2", "u3", "u3"],
                        passed=[True, True, True, True, False, True])

        treated, control = self.sums(counts, obs).sides(comparison_row())

        self.assertEqual(control.table[["batch", "replicate_unit", "n_cells"]].values.tolist(), [["a", "u1", 2]])
        self.assertEqual(control.sums[:, 0].toarray().ravel().tolist(), [3])  # the failing control cell is out
        self.assertEqual(treated.table.replicate_unit.tolist(), ["u1", "u2", "u3"])
        self.assertEqual(treated.sums[:, 0].toarray().ravel().tolist(), [4, 8, 32])

    def test_units_with_too_few_cells_or_no_shared_batch_are_not_fitted(self):
        arms = ["ctrl"] * 10 + ["treated"] * 10 + ["ctrl"] * 10 + ["treated"] * 2
        batches = ["a"] * 20 + ["b"] * 12
        units = ["a1"] * 10 + ["a1"] * 10 + ["b1"] * 10 + ["b1"] * 2
        sums = self.sums(sparse.csr_matrix(np.ones((32, N_GENES), dtype=int)), obs_table(arms, batches, units))
        treated, control = sums.sides(comparison_row())

        kept_treated, kept_control, reason = select_units(treated, control, min_cells=5)
        self.assertEqual(kept_treated.table.batch.tolist(), ["a"])  # batch b's treated unit is below min_cells
        self.assertEqual(kept_control.table.batch.tolist(), ["a"])  # so its control unit goes too
        self.assertIn("1 treated and 1 control units", reason)

        _, _, reason = select_units(*sums.sides(comparison_row(treated_arm="absent")), min_cells=5)
        self.assertEqual(reason, "the treated arm has no passing cells")

        arms = ["ctrl"] * 10 + ["treated"] * 10
        separate = self.sums(sparse.csr_matrix(np.ones((20, N_GENES), dtype=int)),
                             obs_table(arms, ["a"] * 10 + ["b"] * 10, ["a1"] * 10 + ["b1"] * 10))
        _, _, reason = select_units(*separate.sides(comparison_row()), min_cells=5)
        self.assertIn("no batch holds both arms", reason)

    def test_two_units_per_side_in_shared_batches_can_be_fitted(self):
        arms = (["ctrl"] * 5 + ["treated"] * 5) * 2
        units = ["a"] * 10 + ["b"] * 10
        sums = self.sums(sparse.csr_matrix(np.ones((20, N_GENES), dtype=int)), obs_table(arms, units, units))
        treated, control, reason = select_units(*sums.sides(comparison_row()), min_cells=5)
        self.assertEqual(reason, "")
        self.assertEqual(len(treated.table), 2)
        self.assertEqual(len(control.table), 2)

    def test_high_guide_load_control_has_no_guide_of_the_target_or_the_other_pathway_genes(self):
        arms = ["DICER1", "DICER1;DROSHA", "DROSHA", "OTHER;THIRD", "THIRD;FOURTH", "no_guide", "DICER1;THIRD", "OTHER"]
        batches = "aaaaaabb"
        counts = np.zeros((8, N_GENES), dtype=int)
        counts[:, 0] = [1, 2, 4, 8, 16, 32, 64, 128]
        sums = self.sums(counts, obs_table(arms, batches, batches, guides=[2, 2, 2, 2, 2, 0, 2, 2]))
        self.assertTrue(sums.high_guide_load)

        treated, control = sums.sides(comparison_row(treated_arm="DICER1", other_pathway_arms="DROSHA"))

        self.assertEqual(treated.table[["batch", "n_cells"]].values.tolist(), [["a", 2], ["b", 1]])
        self.assertEqual(treated.sums[:, 0].toarray().ravel().tolist(), [3, 64])
        # batch a: OTHER;THIRD, THIRD;FOURTH and the no-guide cell; batch b: OTHER
        self.assertEqual(control.table[["batch", "n_cells"]].values.tolist(), [["a", 3], ["b", 1]])
        self.assertEqual(control.sums[:, 0].toarray().ravel().tolist(), [56, 128])

    def test_low_guide_load_uses_the_listed_control_arm(self):
        counts = sparse.csr_matrix(np.ones((3, N_GENES), dtype=int))
        sums = self.sums(counts, obs_table(["ctrl", "treated", "DROSHA"], "aaa", "uuu"))
        self.assertFalse(sums.high_guide_load)
        _, control = sums.sides(comparison_row())
        self.assertEqual(control.table.n_cells.tolist(), [1])

    def test_ratio_of_the_sums_recovers_a_known_fold_change_after_library_normalisation(self):
        change = np.ones(N_GENES)
        change[:100] = 0.25
        counts, arms, donors = simulate_cells(change, cells_per_unit=200)
        sums = self.sums(counts, obs_table(arms, donors, donors))
        genes = pd.DataFrame({"gene_id": [*GENE_IDS, "ENSG99999999999"], "gene_name": "x"})

        table = ratio_table(genes, sums.gene_ids, *sums.sides(comparison_row()))

        high = table[table.control_sum.astype(float) > 5000]
        changed = high.gene_id.isin(GENE_IDS[:100])
        # library normalisation moves the ratio up by the share of lost counts, about 0.1
        self.assertAlmostEqual(high.log2_ratio_of_sums[changed].median(), -2, delta=0.2)
        self.assertAlmostEqual(high.log2_ratio_of_sums[~changed].median(), 0, delta=0.2)
        self.assertTrue(pd.isna(table.log2_ratio_of_sums.iloc[-1]))  # not in the store
        self.assertEqual(table.analysis_status.unique().tolist(), ["not_fitted"])

    def test_comparison_list_errors(self):
        path = self.root / "list.tsv"
        def write(**values):
            pd.DataFrame([comparison_row(**values)]).drop(columns="other_pathway_arms").to_csv(path, sep="\t", index=False)

        write()
        self.assertEqual(read_comparison_list(path, "TEST").min_cells.tolist(), [5])
        write(orientation="down")
        with self.assertRaisesRegex(ValueError, "orientation"):
            read_comparison_list(path, "TEST")
        write(category="cooccurrence")
        with self.assertRaisesRegex(ValueError, "category"):
            read_comparison_list(path, "TEST")
        write(dataset="OTHER")
        with self.assertRaisesRegex(ValueError, "dataset"):
            read_comparison_list(path, "TEST")
        write(min_cells="0")
        with self.assertRaisesRegex(ValueError, "min_cells"):
            read_comparison_list(path, "TEST")


class CommandTests(unittest.TestCase):
    """The console entrypoint, with a stub Rscript whose DESeq2 step fails as the case needs."""

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        write_gtf(self.root)
        self.genes = self.root / "genes.tsv"
        pd.DataFrame({"gene_id": GENE_IDS, "gene_name": [f"G{i}" for i in range(N_GENES)]}).to_csv(self.genes, sep="\t", index=False)
        arms = (["ctrl"] * 5 + ["treated"] * 5) * 2 + ["alone"] * 5
        units = ["a"] * 10 + ["b"] * 10 + ["a"] * 5
        rng = np.random.default_rng(0)
        self.store = write_store(self.root, sparse.csr_matrix(rng.poisson(1, (25, N_GENES))), obs_table(arms, units, units))
        rows = [comparison_row(comparison_id="two_units"),
                comparison_row(comparison_id="one_unit", treated_arm="alone")]
        self.comparisons = self.root / "list.tsv"
        pd.DataFrame(rows).to_csv(self.comparisons, sep="\t", index=False)

    def tearDown(self):
        self.temporary.cleanup()

    def run_command(self, deseq2_step: str) -> subprocess.CompletedProcess:
        bin_dir = self.root / "bin"
        bin_dir.mkdir()
        stub = bin_dir / "Rscript"
        stub.write_text(f'#!/bin/bash\nif [ "$1" = "-e" ]; then exit 0; fi\n{deseq2_step}\n')
        stub.chmod(stub.stat().st_mode | stat.S_IEXEC)
        env = {**os.environ, "PATH": f"{bin_dir}:{os.environ['PATH']}"}
        return subprocess.run([sys.executable, "-m", "funmirbench.single_cell_final", "--store", str(self.store),
                               "--comparisons", str(self.comparisons), "--genes", str(self.genes), "--out", str(self.root / "out"),
                               "--workers", "2"], env=env, cwd=REPO, capture_output=True, text=True, check=False)

    def test_failed_fit_is_a_result_and_stops_nothing(self):
        done = self.run_command('echo "Error: every gene contains at least one zero" >&2; exit 1')

        self.assertEqual(done.returncode, 0, done.stderr)
        table = pd.read_csv(self.root / "out" / "comparisons.tsv", sep="\t").set_index("comparison_id")
        self.assertEqual(table.loc["two_units", "status"], "fit_failed")
        self.assertIn("every gene contains at least one zero", table.loc["two_units", "reason"])
        self.assertEqual(table.loc["one_unit", "status"], "not_fitted")
        self.assertEqual(table.loc["one_unit", "n_treated_cells"], 5)
        self.assertEqual(table.loc["one_unit", "source"], "test")
        for comparison_id in table.index:
            ratios = pd.read_csv(self.root / "out" / "genes" / f"{comparison_id}.tsv.gz", sep="\t")
            self.assertEqual(len(ratios), N_GENES)
            self.assertGreater(ratios.log2_ratio_of_sums.notna().sum(), 0)
        self.assertIn("one_unit", done.stdout)
        units = (self.root / "out" / "units.tsv").read_text()
        self.assertIn("two_units\tskipped", units)

    def test_killed_fit_is_an_infrastructure_error_and_exits_non_zero(self):
        done = self.run_command("kill -9 $$")

        self.assertEqual(done.returncode, 1)
        table = pd.read_csv(self.root / "out" / "comparisons.tsv", sep="\t").set_index("comparison_id")
        self.assertEqual(table.loc["two_units", "status"], "error")
        self.assertEqual(table.loc["one_unit", "status"], "not_fitted")
        self.assertIn("two_units\tfailed", (self.root / "out" / "units.tsv").read_text())

    def test_high_guide_load_store_goes_through_the_command_with_its_caveat(self):
        arms = (["DICER1"] * 5 + ["DICER1;DROSHA"] * 5 + ["DROSHA"] * 5 + ["OTHER;THIRD"] * 5) * 2
        units = ["a"] * 20 + ["b"] * 20
        rng = np.random.default_rng(0)
        self.store = write_store(self.root, sparse.csr_matrix(rng.poisson(1, (40, N_GENES))),
                                 obs_table(arms, units, units, guides=np.full(40, 2)))
        pd.DataFrame([comparison_row(comparison_id="dicer", treated_arm="DICER1", control_arm="OTHER",
                                     other_pathway_arms="DROSHA", caveats="listed")]).to_csv(self.comparisons, sep="\t", index=False)

        done = self.run_command("exit 1")

        self.assertEqual(done.returncode, 0, done.stderr)
        row = pd.read_csv(self.root / "out" / "comparisons.tsv", sep="\t").iloc[0]
        self.assertEqual((row.n_treated_cells, row.n_control_cells), (20, 10))
        self.assertTrue(row.caveats.startswith("listed; high guide load: median 2 guides per cell"))
        self.assertEqual(row.control_selection, "cells_of_the_batch_without_these_guides")
        design = pd.read_csv(self.root / "out" / "comparisons" / "dicer" / "design.tsv", sep="\t")
        self.assertEqual(design.batch.tolist(), ["a", "b", "a", "b"])


def r_environment_ready() -> bool:
    return (R_ENV / "Rscript").exists()


class RelabelTests(unittest.TestCase):
    def test_a_comparison_id_in_two_datasets_stops_the_relabel_before_any_work(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            for dataset in ("A", "B"):
                (root / "finished" / dataset).mkdir(parents=True)
                pd.DataFrame({"comparison_id": ["same"], "status": ["not_fitted"]}).to_csv(
                    root / "finished" / dataset / "comparisons.tsv", sep="\t", index=False)

            with self.assertRaisesRegex(ValueError, "duplicated comparison_id"):
                run_relabel(root / "finished", root / "out", root / "genes.tsv", workers=1)
            self.assertFalse((root / "out").exists())


@unittest.skipUnless(r_environment_ready(), f"no R with DESeq2 and apeglm at {R_ENV}")
class EndToEndTests(unittest.TestCase):
    def test_relabel_keeps_all_but_the_labels_and_only_the_predicted_direction_is_responsive(self):
        with tempfile.TemporaryDirectory() as temporary_directory, mock.patch.dict(
                os.environ, {"PATH": f"{R_ENV}:{os.environ['PATH']}", "OMP_NUM_THREADS": "1"}):
            from funmirbench.single_cell_final import run_final

            root = Path(temporary_directory)
            write_gtf(root)
            change = np.ones(N_GENES)
            change[:100] = 0.5
            change[100:200] = 2.0
            counts, arms, donors = simulate_cells(change)
            store = write_store(root, counts, obs_table(arms, donors, donors))
            genes = root / "genes.tsv.gz"
            with gzip.open(genes, "wt") as handle:
                pd.DataFrame({"gene_id": GENE_IDS[:800], "gene_name": [f"G{i}" for i in range(800)]}).to_csv(handle, sep="\t", index=False)
            comparisons = root / "list.tsv"
            pd.DataFrame([comparison_row(comparison_id="oe"), comparison_row(comparison_id="ko", orientation="KO", category="ko"),
                          comparison_row(comparison_id="missing_arm", treated_arm="absent")]).to_csv(comparisons, sep="\t", index=False)
            finished = root / "finished"
            finished.mkdir()
            run_final(store, comparisons, genes, finished / "TEST", workers=2)
            (finished / "OLD").mkdir()
            shutil.copy(finished / "TEST" / "comparisons.tsv", finished / "OLD" / "comparisons.tsv")
            (finished / "OLD" / "SUPERSEDED").write_text("replaced by TEST\n")

            combined = run_relabel(finished, root / "relabelled", genes, workers=2).set_index("comparison_id")

            self.assertEqual(combined.status.to_dict(), {"oe": "fitted", "ko": "fitted", "missing_arm": "not_fitted"})
            self.assertFalse((root / "relabelled" / "OLD").exists())
            self.assertTrue((root / "relabelled" / "TEST" / "genes" / "missing_arm.tsv.gz").exists())
            down, up = pd.Series(GENE_IDS[:100]), pd.Series(GENE_IDS[100:200])
            for comparison_id, predicted, other in (("oe", down, up), ("ko", up, down)):
                table = pd.read_csv(root / "relabelled" / "TEST" / "genes" / f"{comparison_id}.tsv.gz", sep="\t")
                responsive = table.gene_id[table.state_d025 == "responsive"]
                self.assertGreater(responsive.isin(predicted).sum(), 50)
                self.assertEqual(responsive.isin(other).sum(), 0)
                self.assertEqual(int(combined.loc[comparison_id, "d05_responsive"]), int((table.state_d05 == "responsive").sum()))

            table_path = root / "relabelled" / "TEST" / "genes" / "oe.tsv.gz"
            before = table_path.read_bytes()
            run_relabel(finished, root / "relabelled", genes, workers=2)
            self.assertEqual(table_path.read_bytes(), before)

            tampered = pd.read_csv(table_path, sep="\t")
            tampered["logFC"] = tampered.logFC + 1
            tampered_path = root / "tampered.tsv.gz"
            tampered.to_csv(tampered_path, sep="\t", index=False)
            with self.assertRaisesRegex(ValueError, "logFC"):
                check_only_labels_changed(finished / "TEST" / "genes" / "oe.tsv.gz", tampered_path)

    def test_known_fold_changes_come_through_deseq2_and_the_labels(self):
        with tempfile.TemporaryDirectory() as temporary_directory, mock.patch.dict(
                os.environ, {"PATH": f"{R_ENV}:{os.environ['PATH']}", "OMP_NUM_THREADS": "1"}):
            from funmirbench.single_cell_final import run_final

            root = Path(temporary_directory)
            write_gtf(root)
            change = np.ones(N_GENES)
            change[:100] = 0.5  # targets: down two-fold in the treated arm
            counts, arms, donors = simulate_cells(change)
            # a fourth donor with control cells only: its batch is not shared, so it must not enter the fit
            counts = sparse.vstack([counts, counts[:60]]).tocsr()
            arms, donors = arms + ["ctrl"] * 60, donors + ["d4"] * 60
            store = write_store(root, counts, obs_table(arms, donors, donors))
            genes = root / "genes.tsv.gz"
            with gzip.open(genes, "wt") as handle:
                pd.DataFrame({"gene_id": GENE_IDS[:800], "gene_name": [f"G{i}" for i in range(800)]}).to_csv(handle, sep="\t", index=False)
            target_gene = GENE_IDS[int(np.asarray(counts[:, :100].sum(axis=0)).argmax())]
            rows = [comparison_row(comparison_id="oe_paired", target_gene_id=f"ENSG00000999999;{target_gene}"),
                    comparison_row(comparison_id="ko_paired", orientation="KO", category="ko"),
                    comparison_row(comparison_id="missing_arm", treated_arm="absent")]
            comparisons = root / "list.tsv"
            pd.DataFrame(rows).to_csv(comparisons, sep="\t", index=False)

            table = run_final(store, comparisons, genes, root / "out", workers=2).set_index("comparison_id")

            self.assertEqual(table.loc["oe_paired", "status"], "fitted", table.loc["oe_paired", "reason"])
            self.assertEqual(table.loc["oe_paired", "design"], "~ batch + condition")
            self.assertEqual(table.loc["oe_paired", "n_treated_units"], 3)
            self.assertEqual(table.loc["oe_paired", "n_control_units"], 3)
            self.assertEqual(table.loc["oe_paired", "n_control_cells"], 180)
            self.assertEqual(table.loc["oe_paired", "replicate_type"], "biological")
            self.assertEqual(table.loc["missing_arm", "status"], "not_fitted")

            final = pd.read_csv(root / "out" / "genes" / "oe_paired.tsv.gz", sep="\t")
            self.assertEqual(final.gene_id.tolist(), GENE_IDS[:800])
            target = final.gene_id.isin(GENE_IDS[:100])
            expressed = final.baseMean > 50
            self.assertAlmostEqual(final.logFC[target & expressed].median(), -1, delta=0.2)
            self.assertAlmostEqual(final.logFC[~target & expressed].median(), 0, delta=0.15)
            responsive = final.state_d025 == "responsive"
            self.assertGreater(responsive[target & expressed].mean(), 0.8)
            self.assertLess(responsive[~target & expressed].mean(), 0.05)

            self.assertEqual(table.loc["oe_paired", "target_used_gene_id"], target_gene)
            self.assertAlmostEqual(table.loc["oe_paired", "target_log2fc"], -1, delta=0.35)
            self.assertAlmostEqual(table.loc["oe_paired", "target_log2_ratio_of_sums"], -1, delta=0.35)
            self.assertAlmostEqual(table.loc["oe_paired", "target_log2fc_shrunk"], -1, delta=0.35)

            flipped = pd.read_csv(root / "out" / "genes" / "ko_paired.tsv.gz", sep="\t")
            self.assertTrue(np.allclose(flipped.oriented_logFC.dropna(), -final.oriented_logFC.dropna()))
            self.assertEqual(table.loc["oe_paired", "d05_responsive"], (final.state_d05 == "responsive").sum())


class CombineTests(unittest.TestCase):
    def test_combined_table_leaves_out_superseded_folders_and_rejects_duplicates(self):
        from funmirbench.single_cell_final import combine_comparisons

        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            for folder, comparison_id in (("a", "x"), ("b", "y"), ("old", "y")):
                (root / folder).mkdir()
                pd.DataFrame({"comparison_id": [comparison_id], "status": ["fitted"]}).to_csv(root / folder / "comparisons.tsv", sep="\t", index=False)
            with self.assertRaises(ValueError):
                combine_comparisons(root)
            (root / "old" / "SUPERSEDED").write_text("b\n")
            combine_comparisons(root)
            table = pd.read_csv(root / "comparisons.tsv", sep="\t")
            self.assertEqual(table.comparison_id.tolist(), ["x", "y"])
            self.assertEqual(table.final_dir.tolist(), ["a", "b"])


class CooccurrenceTests(unittest.TestCase):
    @staticmethod
    def _sparse_mirna_tables(root: Path) -> None:
        """12 cells: one miRNA read in half of them, whole libraries of 1,000 to 12,000 counts."""
        cells = [f"c{i}" for i in range(12)]
        pd.DataFrame({"miR": [1, 0] * 6}, index=pd.Index(cells, name="cell")).to_csv(root / "mirna.tsv", sep="\t")
        pd.DataFrame({"G": np.arange(1, 13) * 10, "H": np.arange(12, 0, -1) * 10}, index=pd.Index(cells, name="cell")).to_csv(root / "mrna.tsv", sep="\t")
        pd.DataFrame({"cell": cells, "group": "g"}).to_csv(root / "cells.tsv", sep="\t", index=False)
        pd.DataFrame({"cell": cells, "library_total": np.arange(1, 13) * 1000}).to_csv(root / "totals.tsv", sep="\t", index=False)

    def test_library_totals_keep_cells_without_mirna_and_set_the_denominator(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self._sparse_mirna_tables(root)
            summary = run_cooccurrence(root / "mirna.tsv", root / "mrna.tsv", root / "cells.tsv", "group", root / "out", library_totals_path=root / "totals.tsv")
            self.assertEqual(summary.n_cells.tolist(), [12])
            self.assertEqual(summary.n_cells_zero_mirna.tolist(), [6])
            self.assertEqual(summary.mirna_normalisation.tolist(), ["whole_library_total"])
            per_cell = pd.read_csv(root / "out/per_cell/g_mirna.tsv.gz", sep="\t", index_col=0)
            self.assertEqual(per_cell.library_size.tolist(), (np.arange(1, 13) * 1000).tolist())
            table = pd.read_csv(root / "out/correlations/g.tsv.gz", sep="\t").set_index("gene")
            share = np.array([1, 0] * 6) / (np.arange(1, 13) * 1000)
            gene = np.arange(1, 13) * 10 / (np.arange(1, 13) * 10 + np.arange(12, 0, -1) * 10)
            self.assertAlmostEqual(table.spearman_r["G"], stats.spearmanr(share, gene).statistic, places=4)

    def test_without_library_totals_the_mirna_share_rule_stays(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self._sparse_mirna_tables(root)
            summary = run_cooccurrence(root / "mirna.tsv", root / "mrna.tsv", root / "cells.tsv", "group", root / "out", min_cells=4)
            self.assertEqual(summary.n_cells.tolist(), [6])
            self.assertEqual(summary.n_cells_without_counts.tolist(), [6])
            self.assertEqual(summary.mirna_normalisation.tolist(), ["miRNA_table_total"])

    def test_library_totals_must_cover_every_cell_and_the_mirna_counts(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self._sparse_mirna_tables(root)
            totals = pd.read_csv(root / "totals.tsv", sep="\t")
            totals.iloc[1:].to_csv(root / "missing.tsv", sep="\t", index=False)
            with self.assertRaisesRegex(ValueError, "no library total"):
                run_cooccurrence(root / "mirna.tsv", root / "mrna.tsv", root / "cells.tsv", "group", root / "o1", library_totals_path=root / "missing.tsv")
            totals.assign(library_total=0.5).to_csv(root / "small.tsv", sep="\t", index=False)
            with self.assertRaisesRegex(ValueError, "exceed"):
                run_cooccurrence(root / "mirna.tsv", root / "mrna.tsv", root / "cells.tsv", "group", root / "o2", library_totals_path=root / "small.tsv")

    def test_allow_fractional_is_off_by_default(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self._sparse_mirna_tables(root)
            table = pd.read_csv(root / "mrna.tsv", sep="\t", index_col=0) + 0.5
            table.to_csv(root / "mrna.tsv", sep="\t")
            with self.assertRaisesRegex(ValueError, "fractional"):
                run_cooccurrence(root / "mirna.tsv", root / "mrna.tsv", root / "cells.tsv", "group", root / "o1", min_cells=4)
            summary = run_cooccurrence(root / "mirna.tsv", root / "mrna.tsv", root / "cells.tsv", "group", root / "o2", min_cells=4, allow_fractional=True)
            self.assertEqual(summary.status.tolist(), ["correlated"])
            table.iloc[0, 0] = -1
            table.to_csv(root / "mrna.tsv", sep="\t")
            with self.assertRaisesRegex(ValueError, "negative"):
                run_cooccurrence(root / "mirna.tsv", root / "mrna.tsv", root / "cells.tsv", "group", root / "o3", min_cells=4, allow_fractional=True)

    def test_correlations_within_groups_match_scipy_and_depth_does_not_correlate(self):
        rng = np.random.default_rng(3)
        n = 200
        groups = np.array(["line A"] * n + ["line B"] * n + ["tiny"] * 3)
        cells = [f"cell{i}" for i in range(len(groups))]
        activity = rng.lognormal(0, 1, len(groups))
        depth = rng.lognormal(0, 1, len(groups))  # shared by both libraries in each cell
        mirna = pd.DataFrame({"miR-1": rng.poisson(20 * activity * depth),
                              "miR-background": rng.poisson(200 * depth),
                              "miR-abundant": rng.poisson(20000 * depth),  # keeps the library from following miR-1
                              "miR-rare": (rng.random(len(groups)) < 0.02).astype(int)}, index=cells)
        anticorrelated = np.where(groups == "line A", rng.poisson(100 * depth / activity), rng.poisson(100 * depth))
        genes = pd.DataFrame({"G_down": anticorrelated, "G_other": rng.poisson(100 * depth),
                              "G_background": rng.poisson(50000 * depth)}, index=cells)

        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            mirna.rename_axis("cell").to_csv(root / "mirna.tsv", sep="\t")
            genes.rename_axis("cell").to_csv(root / "mrna.tsv", sep="\t")
            pd.DataFrame({"cell": cells, "cell_line": groups}).to_csv(root / "cells.tsv", sep="\t", index=False)

            summary = run_cooccurrence(root / "mirna.tsv", root / "mrna.tsv", root / "cells.tsv", "cell_line", root / "out")

            summary = summary.set_index("group")
            self.assertEqual(summary.loc["line A", "status"], "correlated")
            self.assertEqual(summary.loc["tiny", "status"], "not_correlated")
            self.assertIn("fewer than 10", summary.loc["tiny", "reason"])
            table = pd.read_csv(root / "out" / "correlations" / "line_A.tsv.gz", sep="\t")
            expressed = pd.read_csv(root / "out" / "expressed_mirnas.tsv", sep="\t")
            self.assertEqual(sorted(expressed.mirna[expressed.group == "line A"]), ["miR-1", "miR-abundant", "miR-background"])
            self.assertFalse((root / "out" / "correlations" / "tiny.tsv.gz").exists())

            pair = table[(table.mirna == "miR-1") & (table.gene == "G_down")].iloc[0]
            self.assertLess(pair.spearman_r, -0.3)
            self.assertLess(pair.ci95_high, 0)
            self.assertLess(pair.ci95_low, pair.spearman_r)
            mask = groups == "line A"
            mirna_a = (mirna["miR-1"] / mirna.sum(axis=1))[mask]
            gene_a = (genes["G_down"] / genes.sum(axis=1))[mask]
            self.assertAlmostEqual(pair.spearman_r, stats.spearmanr(mirna_a, gene_a).statistic, places=3)

            # same depth in both libraries, no relation after normalisation; raw counts would give r far above 0
            null = table[(table.mirna == "miR-background") & (table.gene == "G_other")].iloc[0]
            self.assertLess(abs(null.spearman_r), 0.2)
            self.assertLess(null.ci95_low, 0)
            self.assertGreater(null.ci95_high, 0)
            self.assertGreater(stats.spearmanr(mirna["miR-background"][mask], genes["G_other"][mask]).statistic, 0.5)

            per_cell = pd.read_csv(root / "out" / "per_cell" / "line_B_mirna.tsv.gz", sep="\t")
            self.assertEqual(len(per_cell), n)
            self.assertIn("library_size", per_cell)

    def test_h5ad_inputs_use_the_gene_id_and_reject_fractional_counts(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            cells = [f"c{i}" for i in range(20)]
            rng = np.random.default_rng(0)
            var = pd.DataFrame({"gene_id": ["ENSG00000000001", "", "ENSG00000000003"]}, index=["a", "b", "c"])
            ad.AnnData(sparse.csr_matrix(rng.poisson(5, (20, 3)).astype(np.int32)), obs=pd.DataFrame(index=cells), var=var).write_h5ad(root / "m.h5ad")
            pd.DataFrame(rng.poisson(5, (20, 2)), index=pd.Index(cells, name="cell"), columns=["miR-1", "miR-2"]).to_csv(root / "s.tsv", sep="\t")
            pd.DataFrame({"cell": cells, "group": "one"}).to_csv(root / "cells.tsv", sep="\t", index=False)

            run_cooccurrence(root / "s.tsv", root / "m.h5ad", root / "cells.tsv", "group", root / "out")

            table = pd.read_csv(root / "out" / "correlations" / "one.tsv.gz", sep="\t")
            self.assertEqual(sorted(table.gene.unique()), ["ENSG00000000001", "ENSG00000000003"])

            pd.DataFrame(rng.random((20, 2)), index=pd.Index(cells, name="cell"), columns=["miR-1", "miR-2"]).to_csv(root / "f.tsv", sep="\t")
            with self.assertRaisesRegex(ValueError, "raw integers"):
                run_cooccurrence(root / "f.tsv", root / "m.h5ad", root / "cells.tsv", "group", root / "out2")


if __name__ == "__main__":
    unittest.main()
