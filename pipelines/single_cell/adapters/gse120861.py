"""GSE120861 Gasperini 2019 at-scale K562 enhancer CRISPRi, high MOI (scPerturb copy of the authors' counts, with per-cell guide
groups). Comparisons are the guide groups within 2 kb of a miRNA gene or non-coding host TSS (coordinates lifted to GRCh38 by the
Claude search). Controls are the other cells of the lane set. Usage: python gse120861.py"""

import re

import anndata as ad
import numpy as np
import pandas as pd

from screen_common import DATA, MIN_CELLS_PER_UNIT, Targets, lane_sets, natural_order, write_screen

DATASET = "GSE120861_K562"
OVERLAP = "/SCRATCH/dtzim01/fun-miRBind-work/single_cell_mirna/search_all_20261030/opus/workers/collections/scratch/gasperini_mirna_overlap.tsv"
CONTROL = "cells_without_these_guides"


def main() -> None:
    targets = Targets()
    overlap = pd.read_csv(OVERLAP, sep="\t", dtype=str, keep_default_na=False)
    overlap = overlap[overlap.screen == "at_scale"]
    groups = sorted(overlap.grna_group_or_site.unique())
    deposit = ad.read_h5ad(DATA / "GSE120861" / "gasperini_2019_atscale.h5ad")
    meta = deposit.obs
    pattern = re.compile("|".join(re.escape(group) for group in sorted(groups, key=len, reverse=True)))
    arms = [";".join(sorted(set(pattern.findall(value)))) or "no_listed_guide" for value in meta.gene.astype(str)]
    units, rule = lane_sets(natural_order(meta["sample"].astype(str).to_numpy()))
    obs = pd.DataFrame({"sample": meta["sample"].astype(str).to_numpy(), "batch": units, "replicate_unit": units, "replicate_type": "lane_set",
                        "arm": arms, "is_control": False, "control_type": "", "guide": "", "n_guides": meta.guide_count.astype(int).to_numpy(),
                        "passed_qc": True})
    rows = []
    for group in groups:
        names = overlap[overlap.grna_group_or_site == group]
        hosts = "; ".join(names.kind + " " + names.name)
        symbols = [name.split("[")[0] for name in names.name]
        known = [targets.table.loc[symbol] for symbol in symbols if symbol in targets.table.index]
        ids = ";".join(row.gene_ids for row in known)
        rows.append({"comparison_id": f"{DATASET}__{group}_KO", "dataset": DATASET, "category": "host_gene", "mirna_or_gene": ";".join(symbols),
                     "orientation": "KO", "treated_arm": group, "control_arm": CONTROL, "min_cells": MIN_CELLS_PER_UNIT,
                     "caveats": f"CRISPRi of a candidate enhancer guide group within 2 kb of {hosts}; not designed for miRNAs; about 28 guides "
                                f"per cell; {rule}; the mature miRNA was not measured", "target_gene_id": ids, "other_pathway_arms": ""})
    write_screen(DATASET, deposit.X.tocsr().astype(np.int32), obs, pd.Series(deposit.var.ensembl_id.astype(str).to_numpy()),
                 comparisons=pd.DataFrame(rows), symbols=pd.Series(deposit.var_names.astype(str).to_numpy()),
                 source="scPerturb GasperiniShendure2019_atscale via the pertpy mirror; counts as in GSE120861")


if __name__ == "__main__":
    main()
