"""Size factors from the control samples: no sample's factor depends on the treated samples."""
import pathlib
import shutil
import subprocess
import tempfile
import unittest

import numpy as np
import pandas as pd

SRC = pathlib.Path(__file__).resolve().parents[1] / "pipelines" / "processing" / "src"
RSCRIPT = shutil.which("Rscript")
HAS_DESEQ2 = RSCRIPT is not None and subprocess.run(
    [RSCRIPT, "-e", "suppressMessages({library(DESeq2); library(apeglm)})"], capture_output=True).returncode == 0
GENES = [f"ENSG{i:011d}" for i in range(400)]


def toy_counts(columns, seed, depth=1.0, repressed=()):
    rng = np.random.default_rng(seed)
    base = np.exp(np.linspace(1, 8, len(GENES)))
    mean = np.outer(base, rng.uniform(0.5, 2.0, len(columns)) * depth)
    mean[list(repressed)] /= 2
    counts = pd.DataFrame(rng.poisson(mean), index=GENES, columns=columns)
    counts.iloc[0] = 0  # counted but zero in every sample
    return counts


def expected_factor(counts, controls, column):
    with np.errstate(divide="ignore"):
        log_reference = np.log(counts[controls]).mean(axis=1)
        used = np.isfinite(log_reference) & (counts[column] > 0)
        return np.exp(np.median(np.log(counts[column][used]) - log_reference[used]))


@unittest.skipIf(RSCRIPT is None, "Rscript not found")
class SizeFactorTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = pathlib.Path(tmp.name)

    def factors(self, counts, controls):
        path = self.root / "counts.tsv"
        counts.to_csv(path, sep="\t")
        code = (f'source("{SRC / "size_factors.R"}"); m <- as.matrix(read.delim("{path}", row.names = 1, check.names = FALSE)); '
                f'write.table(control_size_factors(m, c({", ".join(repr(c) for c in controls)})), stdout(), sep = "\\t", '
                'quote = FALSE, row.names = FALSE)')
        result = subprocess.run([RSCRIPT, "-e", code], capture_output=True, text=True)
        return result, (pd.read_csv(pd.io.common.StringIO(result.stdout), sep="\t").set_index("sample_id")
                        if result.returncode == 0 else None)

    def test_factors_use_only_the_controls_and_the_sample_itself(self):
        controls = ["c1", "c2"]
        counts = pd.concat([toy_counts(controls, seed=2), toy_counts(["t1", "t2"], seed=3, repressed=range(50, 90))], axis=1)
        _, together = self.factors(counts, controls)
        for column in counts:
            self.assertAlmostEqual(together.size_factor[column], expected_factor(counts, controls, column), places=10)
        _, alone = self.factors(counts[controls], controls)
        _, reordered = self.factors(counts[["t2", "c2", "t1", "c1"]], controls)
        _, scaled = self.factors(counts.assign(t1=counts.t1 * 11, t2=counts.t2 * 11), controls)
        self.assertEqual(alone.size_factor[controls].tolist(), together.size_factor[controls].tolist())
        self.assertTrue(reordered.size_factor.sort_index().equals(together.size_factor.sort_index()))
        self.assertEqual(scaled.size_factor[controls].tolist(), together.size_factor[controls].tolist())

    def test_a_sample_sharing_no_gene_with_the_controls_fails(self):
        counts = toy_counts(["c1", "c2", "t1"], seed=2)
        counts.loc[counts.c1 == 0, "t1"] = 5
        counts.loc[counts.c1 > 0, "t1"] = 0
        result, _ = self.factors(counts, ["c1", "c2"])
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("No gene is counted above 0 in every control sample and in: t1", result.stderr)

    @unittest.skipUnless(HAS_DESEQ2, "DESeq2 and apeglm not installed")
    def test_deseq2_keeps_the_factors_and_a_shared_control_is_unchanged(self):
        controls = toy_counts(["c1", "c2", "c3"], seed=4)
        runs = {}
        for name, treated in {"a": toy_counts(["t1", "t2", "t3"], seed=5, repressed=range(100, 160)),
                              "b": toy_counts(["u1", "u2", "u3"], seed=6, depth=5.0, repressed=range(200, 260))}.items():
            counts = pd.concat([controls, treated], axis=1).rename_axis("gene_id").reset_index()
            out = self.root / name
            out.mkdir()
            counts.to_csv(out / "counts.tsv", sep="\t", index=False)
            subprocess.run([RSCRIPT, str(SRC / "run_deseq2_counts.R"), "--counts", str(out / "counts.tsv"), "--gene-id-column", "gene_id",
                            "--control-columns", "c1,c2,c3", "--treated-columns", ",".join(treated.columns),
                            "--output", str(out / "de_results.tsv")], check=True, capture_output=True)
            runs[name] = out
            factors = pd.read_csv(out / "sample_size_factors.tsv", sep="\t").set_index("sample_id").size_factor
            for column in factors.index:
                self.assertAlmostEqual(factors[column], expected_factor(counts.set_index("gene_id"), ["c1", "c2", "c3"], column),
                                       places=10)
        normalized = {k: pd.read_csv(v / "normalized_counts.tsv", sep="\t", index_col=0) for k, v in runs.items()}
        pd.testing.assert_frame_equal(normalized["a"][["c1", "c2", "c3"]], normalized["b"][["c1", "c2", "c3"]])

        genes = pd.DataFrame({"gene_id": [*GENES, "ENSG99999999999"], "gene_name": "g"})
        genes.to_csv(self.root / "genes.tsv", sep="\t", index=False)
        older = self.root / "older.rds"
        subprocess.run([RSCRIPT, "-e", f'a <- readRDS("{runs["a"] / "complete_analysis.rds"}"); a$size_factor_method <- NULL; '
                        f'saveRDS(a, "{older}")'], check=True, capture_output=True)
        label = [RSCRIPT, str(SRC / "label_states.R")]
        refused = subprocess.run([*label, str(older), str(self.root / "genes.tsv"), str(self.root / "refused")],
                                 capture_output=True, text=True)
        self.assertNotEqual(refused.returncode, 0)
        self.assertIn("fitted before size factors from the control samples", refused.stderr)
        subprocess.run([*label, str(runs["a"] / "complete_analysis.rds"), str(self.root / "genes.tsv"), str(self.root / "table")],
                       check=True, capture_output=True)
        table = pd.read_csv(self.root / "table.tsv.gz", sep="\t").set_index("gene_id")
        self.assertEqual(table.loc[GENES[0], ["control_mean_normalized_count", "treated_nonzero_samples"]].tolist(), [0, 0])
        self.assertTrue(np.isnan(table.loc["ENSG99999999999", "control_mean_normalized_count"]))
        expected = normalized["a"][["c1", "c2", "c3"]].mean(axis=1)
        fitted = table.analysis_status != "not_fitted_all_zero"
        np.testing.assert_allclose(table.control_mean_normalized_count[fitted], expected[table.index[fitted]], rtol=1e-10)


if __name__ == "__main__":
    unittest.main()
