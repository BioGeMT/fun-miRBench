"""Perturb-ME A375 (figshare 33264225): Cas9 knockout, Perturb-CITE-seq, cells sorted on surface HLA. The MOI1 object holds
sgRNA-singlet cells; X has raw RNA counts and normalised CITE columns, which are left out. Each sorted population (Low,
High, Control) is one batch and one replicate unit, as cells per guide per channel are too few for channels. Usage: python perturbme.py
"""

import h5py
import numpy as np
import pandas as pd

from screen_common import DATA, Targets, comparison_rows, h5_frame, h5_rows, plan_cells, write_screen

CONTROLS = ("NO_SITE", "ONE_NON-GENE_SITE")
CONTROL = "no_site_or_non_gene"
DATASET = "PerturbME_A375"


def main() -> None:
    targets = Targets()
    path = DATA / "perturbme_a375" / "PerturbME_A375_HLA_MOI1.h5ad"
    with h5py.File(path, "r") as handle:
        deposit = h5_frame(handle["obs"])
        var = h5_frame(handle["var"])
        target = deposit.sgRNA.astype(str).str.rsplit("_", n=1).str[0]
        control = target.isin(CONTROLS).to_numpy()
        condition = deposit.Condition.astype(str).to_numpy()
        obs = pd.DataFrame({"sample": deposit.Sample.astype(str).to_numpy(), "batch": condition, "replicate_unit": condition,
                            "replicate_type": "lane_set", "arm": np.where(control, CONTROL, target), "is_control": control,
                            "control_type": np.where(control, target, ""), "guide": deposit.sgRNA.astype(str).to_numpy(),
                            "n_guides": 1, "passed_qc": True})
        kept, rows = plan_cells(obs, targets)
        whole = h5_rows(handle, rows, integer=False)
    rna = (var.modality == "RNA").to_numpy()
    counts = whole[:, np.flatnonzero(rna)]
    if counts.data.size and not np.array_equal(counts.data, np.round(counts.data)):
        raise ValueError("RNA block is not integer")
    caveats = ("Cas9 knockout; cells were sorted on surface HLA (Low, High, Control populations), so cells are not a "
               "random sample; each sorted population is one batch and one unit; very few cells per gene")
    comparisons = comparison_rows(DATASET, kept.arm[~kept.is_control], targets, control_arm=CONTROL, caveats=caveats)
    write_screen(DATASET, counts.astype(np.int32), kept, pd.Series(var.index[rna], dtype=str), comparisons=comparisons,
                 source="authors' raw RNA counts, PerturbME_A375_HLA_MOI1.h5ad, figshare 33264225; cells as the authors filtered them")


if __name__ == "__main__":
    main()
