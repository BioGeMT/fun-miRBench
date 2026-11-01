"""GSE249597 A549 CRISPRi mini-library, MOI about 0.3, cells hashed by oxygen condition. Authors' processed RNA with HTO counts and
guide UMI tables per library. Condition = hash with the most UMIs and at least half of the cell's hash UMIs; guide with
screen_common.call_guides. Comparisons within each condition; the two libraries are batches and units. Usage: python gse249597.py <RAW folder>"""

import pathlib
import sys

import numpy as np
import pandas as pd
from scipy import sparse

from funmirbench.single_cell_store import write_cell_store
from funmirbench.single_cell_pipeline import read_10x_mtx
from screen_common import ADAPTERS, FINAL, GTF, MIN_CELLS_PER_UNIT, Targets, call_guides

DATASET = "GSE249597_A549"
CONDITIONS = {"Normoxie": "normoxia", "Hypoxie_3H": "hypoxia3h", "Hypoxie_6H": "hypoxia6h", "Hypoxie_24h": "hypoxia24h"}
CONTROL = "Neg"
GUIDE_TARGET = {"LINC00152": "CYTOR"}


def read_library(root: pathlib.Path, library: str) -> tuple:
    folder = root / "data" / "truchi_data" / "DATA" / "CRISPRi_2"
    rna = folder / f"{library}_processed_RNA"
    import gzip
    import scipy.io as scipy_io
    whole = sparse.csr_matrix(scipy_io.mmread(gzip.open(rna / "matrix.mtx.gz"))).T.tocsr()
    features = pd.read_csv(rna / "features.tsv.gz", sep="\t", header=None, dtype=str)
    barcodes = pd.read_csv(rna / "barcodes.tsv.gz", header=None, dtype=str)[0]
    hto = (features[2] == "Antibody Capture").to_numpy()
    hto_counts = whole[:, np.flatnonzero(hto)].toarray()
    labels = features[0][hto].to_numpy()
    share = hto_counts.max(axis=1) / np.maximum(hto_counts.sum(axis=1), 1)
    condition = np.where(share >= 0.5, np.array([CONDITIONS[label] for label in labels])[hto_counts.argmax(axis=1)], "unassigned")
    guide_rows, guide_names = [], None
    lookup = pd.Series(np.arange(len(barcodes)), index=barcodes.str.replace("-1$", "", regex=True))
    total = None
    for lane in ("L001", "L002"):
        umi = folder / "gRNA_counts" / f"{library}_{lane}" / "umi_count"
        guide = sparse.csr_matrix(scipy_io.mmread(gzip.open(umi / "matrix.mtx.gz"))).T.tocsr()
        names = pd.read_csv(umi / "features.tsv.gz", sep="\t", header=None, dtype=str)[0].to_numpy()
        codes = pd.read_csv(umi / "barcodes.tsv.gz", header=None, dtype=str)[0].to_numpy()
        keep = names != "unmapped"
        guide = guide[:, np.flatnonzero(keep)]
        names = names[keep]
        where = lookup.reindex(codes).to_numpy()
        found = ~np.isnan(where)
        order = np.argsort(names)
        guide, names = guide[:, order], names[order]
        placed = sparse.csr_matrix((np.ones(found.sum(), dtype=np.int32), (where[found].astype(int), np.flatnonzero(found))),
                                   shape=(len(barcodes), guide.shape[0])) @ guide
        total = placed if total is None else total + placed
        guide_names = names
    return whole[:, np.flatnonzero(~hto)], features[~hto], barcodes, condition, total.tocsr(), guide_names


def main(root: pathlib.Path) -> None:
    targets = Targets()
    parts, tables, features = [], [], None
    for library in ("Lib1", "Lib2"):
        counts, genes, barcodes, condition, guide, names = read_library(root, library)
        features = genes
        target = np.array([name.split("-sg")[0] for name in names], dtype=object)
        target = np.array([GUIDE_TARGET.get(value, value) for value in target], dtype=object)
        arms, numbers, called = call_guides(guide, target)
        parts.append(counts)
        tables.append(pd.DataFrame({"sample": library, "batch": library, "replicate_unit": library, "replicate_type": "lane_set",
                                    "guide_arm": arms, "condition": condition, "guide": called, "n_guides": numbers}))
    table = pd.concat(tables, ignore_index=True)
    single = table.guide_arm.str.fullmatch(r"[^;]+").to_numpy() & (table.condition != "unassigned").to_numpy()
    arm = np.where(single, table.guide_arm + "|" + table.condition, "other")
    control = np.array([value.startswith(f"{CONTROL}|") for value in arm])
    obs = table.drop(columns=["guide_arm", "condition"]).assign(dataset=DATASET, arm=arm, is_control=control,
                                                                 control_type=np.where(control, "non_targeting", ""),
                                                                 passed_qc=single)
    listed = [symbol for symbol in pd.unique(table.guide_arm[single]) if symbol in targets.table.index]
    rows = []
    for symbol in listed:
        for condition in CONDITIONS.values():
            target_row = targets.table.loc[symbol]
            rows.append({"comparison_id": f"{DATASET}__{symbol}_{condition}_KO", "dataset": DATASET, "category": target_row.category,
                         "mirna_or_gene": symbol, "orientation": "KO", "treated_arm": f"{symbol}|{condition}",
                         "control_arm": f"{CONTROL}|{condition}", "min_cells": MIN_CELLS_PER_UNIT,
                         "caveats": f"CRISPRi knockdown in A549, {condition}; miRNAs inside: {target_row.mirnas_inside}; the CYTOR guide also silences MIR4435-2HG; "
                                    "units are the two libraries; hash and guide calls by a fixed rule on those counts; the mature miRNA was not measured",
                         "target_gene_id": target_row.gene_ids, "other_pathway_arms": ""})
    comparisons = pd.DataFrame(rows)
    store_path = FINAL / DATASET / "store" / f"{DATASET}.h5ad"
    store_path.parent.mkdir(parents=True, exist_ok=True)
    write_cell_store(store_path, sparse.vstack(parts).tocsr(), obs, features[0].astype(str), gtf=GTF, symbols=features[1].astype(str),
                     source="authors' processed RNA, HTO and guide counts, GEO GSE249597")
    comparisons.to_csv(ADAPTERS.parent / "comparisons" / f"{DATASET}.tsv", sep="\t", index=False)
    print(DATASET, obs.shape, comparisons.shape, flush=True)


if __name__ == "__main__":
    main(pathlib.Path(sys.argv[1]))
