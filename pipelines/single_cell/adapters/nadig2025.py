"""Nadig et al. 2025 CRISPRi Perturb-seq (GSE264667), authors' raw single-cell files. Same layout as Replogle 2022: only
the cells and genes the authors kept; unit = lane set of eight GEM groups.

Usage: python nadig2025.py <jurkat|hepg2>
"""

import sys

import h5py
import numpy as np
import pandas as pd

from screen_common import DATA, Targets, comparison_rows, h5_frame, h5_rows, lane_sets, plan_cells, write_screen

CONTROL = "non-targeting"


def main(line: str) -> None:
    dataset = {"jurkat": "Nadig2025_Jurkat", "hepg2": "Nadig2025_HepG2"}[line]
    targets = Targets()
    name = f"GSE264667_{line}_raw_singlecell_01.h5ad"
    with h5py.File(DATA / "GSE264667" / name, "r") as handle:
        deposit = h5_frame(handle["obs"])
        var = h5_frame(handle["var"])
        units, rule = lane_sets(deposit.gem_group.astype(int).to_numpy())
        control = deposit.gene.to_numpy() == CONTROL
        obs = pd.DataFrame({
            "sample": "lane_" + deposit.gem_group.astype(str).to_numpy(), "batch": units, "replicate_unit": units,
            "replicate_type": "lane_set", "arm": deposit.gene.to_numpy(), "is_control": control,
            "control_type": np.where(control, "non_targeting", ""), "guide": deposit.gene_transcript.to_numpy(),
            "n_guides": 1, "passed_qc": True})
        kept, rows = plan_cells(obs, targets)
        counts = h5_rows(handle, rows)
    caveats = (f"CRISPRi knockdown; {rule}; authors' matrix holds only the genes they kept ({var.shape[0]}), other "
               "listed genes have no counts")
    comparisons = comparison_rows(dataset, kept.arm[~kept.is_control], targets, control_arm=CONTROL, caveats=caveats)
    write_screen(dataset, counts, kept, pd.Series(var.index, dtype=str), comparisons=comparisons,
                 symbols=pd.Series(var.gene_name.to_numpy(), dtype=str),
                 source=f"authors' raw counts, {name}, GEO GSE264667; cells as the authors filtered them")


if __name__ == "__main__":
    main(sys.argv[1])
