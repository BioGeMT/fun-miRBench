"""Size factors from the control samples: a control's normalised counts do not depend on the treated samples."""
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
    mean = np.outer(np.exp(np.linspace(1, 8, len(GENES))), rng.uniform(0.5, 2.0, len(columns)) * depth)
    mean[list(repressed)] /= 2
    return pd.DataFrame(rng.poisson(mean), index=GENES, columns=columns)


@unittest.skipUnless(HAS_DESEQ2, "DESeq2 and apeglm not installed")
class SizeFactorTests(unittest.TestCase):
    def test_a_shared_control_is_unchanged_and_old_fits_are_refused(self):
        root = pathlib.Path(self.enterContext(tempfile.TemporaryDirectory()))
        controls = toy_counts(["c1", "c2", "c3"], seed=4)
        normalized = {}
        for name, treated in {"a": toy_counts(["t1", "t2", "t3"], seed=5, repressed=range(100, 160)),
                              "b": toy_counts(["u1", "u2", "u3"], seed=6, depth=5.0, repressed=range(200, 260))}.items():
            out = root / name
            out.mkdir()
            pd.concat([controls, treated], axis=1).rename_axis("gene_id").to_csv(out / "counts.tsv", sep="\t")
            subprocess.run([RSCRIPT, str(SRC / "run_deseq2_counts.R"), "--counts", str(out / "counts.tsv"), "--gene-id-column",
                            "gene_id", "--control-columns", "c1,c2,c3", "--treated-columns", ",".join(treated.columns),
                            "--output", str(out / "de_results.tsv")], check=True, capture_output=True)
            normalized[name] = pd.read_csv(out / "normalized_counts.tsv", sep="\t", index_col=0)[["c1", "c2", "c3"]]
        pd.testing.assert_frame_equal(normalized["a"], normalized["b"], rtol=1e-12)

        pd.DataFrame({"gene_id": GENES, "gene_name": "g"}).to_csv(root / "genes.tsv", sep="\t", index=False)
        fit, older = root / "a" / "complete_analysis.rds", root / "older.rds"
        subprocess.run([RSCRIPT, "-e", f'a <- readRDS("{fit}"); a$size_factors <- NULL; saveRDS(a, "{older}")'],
                       check=True, capture_output=True)
        label = [RSCRIPT, str(SRC / "label_states.R")]
        refused = subprocess.run([*label, str(older), str(root / "genes.tsv"), str(root / "old")], capture_output=True, text=True)
        self.assertNotEqual(refused.returncode, 0)
        self.assertIn("fitted before size factors from the control samples", refused.stderr)
        subprocess.run([*label, str(fit), str(root / "genes.tsv"), str(root / "new")], check=True, capture_output=True)
        table = pd.read_csv(root / "new.tsv.gz", sep="\t", index_col="gene_id").dropna(subset="control_mean_normalized_count")
        for group, samples, mean in (("control", ["c1", "c2", "c3"], table.control_mean_normalized_count),
                                     ("treated", ["t1", "t2", "t3"], table.treated_mean_normalized_count)):
            per_sample = table[[f"{group}_normalized_count_{s}" for s in samples]]
            np.testing.assert_allclose(per_sample.mean(axis=1), mean, rtol=1e-12)
        np.testing.assert_allclose(table[[f"control_normalized_count_{s}" for s in ["c1", "c2", "c3"]]],
                                   normalized["a"].loc[table.index], rtol=1e-12)


if __name__ == "__main__":
    unittest.main()
