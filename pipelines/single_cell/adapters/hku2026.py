"""Screens taken from the HKU uniform reprocessing (10.25442/hku.34002315.v1): X holds log1p(normalize_total(counts, 1e4)),
raw counts are round(expm1(X) * library_size / 1e4). Used where the authors' files exceed about 100 GB (X-Atlas/Orion,
Nourreddine) or hold no per-cell guide calls (Pan: only HKU has them). HKU drops perturbations below 20 cells, so absent targets may exist in the authors' files.

Usage: python hku2026.py <XAtlasOrion_HCT116|XAtlasOrion_HEK293T|Pan2026_H1|Nourreddine2026_KOLF>
"""

import sys

import h5py
import numpy as np
import pandas as pd

from funmirbench.single_cell_store import raw_counts_from_lognorm
from screen_common import DATA, Targets, comparison_rows, h5_frame, h5_rows, lane_sets, natural_order, plan_cells, write_screen

SCREENS = {
    "XAtlasOrion_HCT116": dict(prefix="Huang-HCT116", control="Non-Targeting", batch="sample", arm="gene_target",
                         deposit="X-Atlas/Orion HCT116, figshare+ 29190726"),
    "XAtlasOrion_HEK293T": dict(prefix="Huang-HEK293T", control="Non-Targeting", batch="sample", arm="gene_target",
                          deposit="X-Atlas/Orion HEK293T, figshare+ 29190726"),
    "Pan2026_H1": dict(prefix="Pan-GW-hESC", control="NTC", batch="batch", arm="perturbation",
                       deposit="Pan H1 hESC, GEO GSE295214"),
    "Nourreddine2026_KOLF": dict(prefix="Nourreddine-GW-ipsc", control="NTC", batch="batch", arm="perturbation",
                             deposit="Nourreddine KOLF2.1J iPSC, figshare+ 27261219"),
}


def main(dataset: str) -> None:
    screen = SCREENS[dataset]
    targets = Targets()
    path = DATA / "hku2026" / f"sc_{screen['prefix']}.h5ad"
    with h5py.File(path, "r") as handle:
        deposit = h5_frame(handle["obs"])
        var = h5_frame(handle["var"])
        names = deposit[screen["batch"]].astype(str).to_numpy()
        units, rule = lane_sets(natural_order(names))
        control = deposit[screen["arm"]].astype(str).to_numpy() == screen["control"]
        guide_column = next((column for column in ("guide_id", "gRNA") if column in deposit), None)
        guide = deposit[guide_column].astype(str).to_numpy() if guide_column else ""
        obs = pd.DataFrame({
            "sample": names, "batch": units, "replicate_unit": units, "replicate_type": "lane_set",
            "arm": deposit[screen["arm"]].astype(str).to_numpy(), "is_control": control,
            "control_type": np.where(control, "non_targeting", ""), "guide": guide, "n_guides": 1, "passed_qc": True})
        kept, rows = plan_cells(obs, targets)
        lognorm = h5_rows(handle, rows, integer=False)
    counts = raw_counts_from_lognorm(lognorm, deposit.library_size.to_numpy()[rows])
    caveats = (f"CRISPRi knockdown; {rule}; counts reconstructed from the HKU log-normalised file; HKU leaves out "
               "perturbations below 20 cells and genes it filtered")
    comparisons = comparison_rows(dataset, kept.arm[~kept.is_control], targets, control_arm=screen["control"],
                                  caveats=caveats)
    symbols = var.gene_name.astype(str).to_numpy() if "gene_name" in var else var.index.to_numpy()
    id_column = next((column for column in ("gene_ids", "gene_id", "ensembl_id") if column in var), None)
    features = var[id_column].astype(str).to_numpy() if id_column else var.index.to_numpy()
    write_screen(dataset, counts, kept, pd.Series(features, dtype=str), comparisons=comparisons,
                 symbols=pd.Series(symbols, dtype=str),
                 source=f"HKU 10.25442/hku.34002315.v1 sc_{screen['prefix']} reconstructed counts ({screen['deposit']}); "
                        "cells as HKU filtered them")


if __name__ == "__main__":
    main(sys.argv[1])
