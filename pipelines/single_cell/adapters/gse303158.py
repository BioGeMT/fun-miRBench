"""GSE303158 hPSC knockout screen, primed and naive, deposited Seurat objects (no raw counts layer). Counts are reconstructed
as round(expm1(data) * nCount_RNA / 1e4). All cells with a guide are used; mixscape classes depend on response and are not.
Each batch is a batch and a unit. Usage: python gse303158.py <primed|naive>"""

import sys

import numpy as np
import pandas as pd
from scipy import io, sparse

from funmirbench.single_cell_store import raw_counts_from_lognorm
from screen_common import DATA, Targets, comparison_rows, plan_cells, write_screen

CONTROL = "CTRL"


def main(state: str) -> None:
    dataset = f"GSE303158_hPSC_{state}"
    root = DATA / "GSE303158" / ("x" if state == "primed" else "x_naive")
    targets = Targets()
    meta = pd.read_csv(root / "meta.tsv", sep="\t")
    features = (root / "features.tsv").read_text().splitlines()
    data = sparse.csr_matrix(io.mmread(root / "counts.mtx").T)
    batch = meta.orig_ident_batch.astype(str).to_numpy() if "orig_ident_batch" in meta else meta["batch"].astype(str).to_numpy()
    guide = meta.guide.astype(str)
    control = (guide == CONTROL).to_numpy()
    obs = pd.DataFrame({"sample": batch, "batch": batch, "replicate_unit": batch, "replicate_type": "lane_set", "arm": guide.str.split("_sg_").str[0],
                        "is_control": control, "control_type": np.where(control, "non_targeting", ""), "guide": guide, "n_guides": 1,
                        "passed_qc": True})
    kept, rows = plan_cells(obs, targets)
    counts = raw_counts_from_lognorm(data[rows], meta.nCount_RNA.to_numpy()[rows])
    caveats = ("Cas9 knockout, pluripotency library, MOI 0.3; replicate units are the two batches; counts reconstructed from the deposited "
               "log-normalised layer and nCount_RNA; all guide cells used, mixscape classes ignored")
    comparisons = comparison_rows(dataset, kept.arm[~kept.is_control], targets, control_arm=CONTROL, caveats=caveats)
    write_screen(dataset, counts, kept, pd.Series(features, dtype=str), comparisons=comparisons,
                 source=f"authors' Seurat object, GSE303158 {state}, counts reconstructed from data layer")


if __name__ == "__main__":
    main(sys.argv[1])
