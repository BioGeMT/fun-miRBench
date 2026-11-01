"""GSE181338 CROP-seq of 361 miRNA loci in K562 (Li et al. 2023). Cells and guide assignments rebuilt from the deposit in
single_cell_mirna/analysis_mirloci_cropseq (top guide with no other above 20% of its UMIs; rule fixed before any response).
Library 1 keeps the precursor (host_gene), library 2 deletes it (ko). Library 3 holds the non-targeting cells in its own
sample: control pseudo-replicates, treated replicates are the three screens, one batch. Usage: python gse181338.py <host|ko>"""

import pathlib
import sys

import numpy as np
import pandas as pd
from scipy import sparse

from funmirbench.single_cell_store import pseudo_replicate_units, read_10x_h5, write_cell_store
from screen_common import ADAPTERS, DATA, FINAL, GTF, MIN_CELLS_PER_UNIT

CONTROL = "non_targeting"
SCREENS = {"host": ("GSE181338_K562_host", "host only", "host_gene", "library 1: deletion downstream of the precursor, precursor kept, host transcript lost"),
           "ko": ("GSE181338_K562_ko", "miRNA and host", "ko", "library 2: precursor and host transcript deleted")}


def main(key: str) -> None:
    dataset, deletion, category, design = SCREENS[key]
    root = DATA / "GSE181338"
    meta = pd.read_csv(root / "processed" / "analysed_cells_obs.tsv.gz", sep="\t", dtype=str, keep_default_na=False)
    counts = sparse.load_npz(root / "processed" / "analysed_cells.npz").tocsr()
    assigned = meta.assigned == "True"
    treated = assigned & (meta.deletion_type == deletion) & (meta.locus != "") & (meta.locus != "unannotated")
    control = (meta.classification == "NTC").to_numpy()
    rows = np.flatnonzero(treated.to_numpy() | control)
    table = meta.iloc[rows]
    is_control = control[rows]
    replicate = table.replicate.to_numpy()
    obs = pd.DataFrame({"sample": table.lane_id.to_numpy(), "batch": "K562", "replicate_unit": "", "replicate_type": np.where(is_control, "technical", "biological"),
                        "arm": np.where(is_control, CONTROL, table.locus), "is_control": is_control,
                        "control_type": np.where(is_control, "non_targeting", ""), "guide": table.top_guide.to_numpy(),
                        "n_guides": 1, "passed_qc": True})
    units = np.array([f"screen_{value}" for value in replicate], dtype=object)
    units[is_control] = pseudo_replicate_units(table.lane_id.to_numpy()[is_control], obs.arm.to_numpy()[is_control],
                                               np.ones(is_control.sum(), dtype=bool))
    obs["replicate_unit"] = units
    features = read_10x_h5(next((root / "geo" / "h5").glob("*.h5")))[1]
    caveats = (f"{design}; paired-guide CRISPR deletion at MOI about 4, one assigned guide per cell, guides rebuilt from the guide "
               "FASTQ; non-targeting cells are one separate sample (library 3), so treatment is confounded with sample; treated "
               "units are the three screens, control units are three technical pseudo-replicates of that sample")
    loci = sorted(obs.arm[~obs.is_control].unique())
    comparisons = pd.DataFrame([{"comparison_id": f"{dataset}__{locus}_KO", "dataset": dataset, "category": category, "mirna_or_gene": locus,
                                 "orientation": "KO", "treated_arm": locus, "control_arm": CONTROL, "min_cells": MIN_CELLS_PER_UNIT,
                                 "caveats": caveats, "target_gene_id": "", "other_pathway_arms": ""} for locus in loci])
    store_path = FINAL / dataset / "store" / f"{dataset}.h5ad"
    store_path.parent.mkdir(parents=True, exist_ok=True)
    write_cell_store(store_path, counts[rows], obs.assign(dataset=dataset), features[0].astype(str), gtf=GTF, symbols=features[1].astype(str),
                     source="authors' 10x counts, GSE181338 raw h5, cells and guides as rebuilt in analysis_mirloci_cropseq")
    comparisons.to_csv(ADAPTERS.parent / "comparisons" / f"{dataset}.tsv", sep="\t", index=False)
    print(dataset, len(rows), len(loci), flush=True)


if __name__ == "__main__":
    main(sys.argv[1])
