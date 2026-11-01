"""Replogle et al. 2022 CRISPRi Perturb-seq, authors' raw single-cell files (figshare+ 20029387). The files hold only
cells the authors kept and only the genes they kept (8,248 to 8,749), so every stored cell passes; unit = lane set.

Usage: python replogle2022.py <rpe1|k562_essential|k562_gw>
"""

import sys

import h5py
import numpy as np
import pandas as pd

from screen_common import DATA, Targets, comparison_rows, h5_frame, h5_rows, lane_sets, plan_cells, write_screen

FILES = {"rpe1": "rpe1_raw_singlecell_01.h5ad", "k562_essential": "K562_essential_raw_singlecell_01.h5ad",
         "k562_gw": "K562_gwps_raw_singlecell_01.h5ad"}
DATASETS = {"rpe1": "Replogle2022_RPE1", "k562_essential": "Replogle2022_K562_essential", "k562_gw": "Replogle2022_K562_GW"}
CONTROL = "non-targeting"


def main(screen: str) -> None:
    dataset = DATASETS[screen]
    targets = Targets()
    path = DATA / "replogle2022" / FILES[screen]
    with h5py.File(path, "r") as handle:
        deposit = h5_frame(handle["obs"])
        var = h5_frame(handle["var"])
        units, rule = lane_sets(deposit.gem_group.astype(int).to_numpy())
        control = deposit.gene.to_numpy() == CONTROL
        obs = pd.DataFrame({
            "sample": "lane_" + deposit.gem_group.astype(str).to_numpy(), "batch": units, "replicate_unit": units,
            "replicate_type": "lane_set", "arm": deposit.gene.to_numpy(), "is_control": control,
            "control_type": np.where(control, "non_targeting", ""), "guide": deposit.sgID_AB.to_numpy()
            if "sgID_AB" in deposit else deposit.gene_transcript.to_numpy(), "n_guides": 1, "passed_qc": True})
        kept, rows = plan_cells(obs, targets)
        counts = h5_rows(handle, rows)
    caveats = (f"CRISPRi knockdown; {rule}; authors' matrix holds only the genes they kept ({var.shape[0]}), other "
               "listed genes have no counts")
    comparisons = comparison_rows(dataset, kept.arm[~kept.is_control], targets, control_arm=CONTROL, caveats=caveats)
    write_screen(dataset, counts, kept, pd.Series(var.index, dtype=str), comparisons=comparisons,
                 symbols=pd.Series(var.gene_name.to_numpy(), dtype=str),
                 source=f"authors' raw counts, {FILES[screen]}, figshare+ 20029387; cells as the authors filtered them")


if __name__ == "__main__":
    main(sys.argv[1])
