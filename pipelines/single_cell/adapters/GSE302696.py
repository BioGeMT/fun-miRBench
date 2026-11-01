"""Engineered vesicles carrying PD-L1 and miR-27a-3p against no vesicles, human CD4+ T cells from humanised mice with colitis.
Each arm's two sublibraries share over 99% of their barcodes, so they are re-sequencing of one capture: counts are summed
per barcode and the cells split into three pseudo-replicates. QC: deposit's filter."""

import numpy as np
import pandas as pd
import pathlib
from scipy import sparse

from _sample_arms import DATA, Library, build_store, store_path
from funmirbench.single_cell_store import read_10x_mtx

RAW = DATA / "GSE302696" / "extracted"


def merged_capture(stems):
    parts = [read_10x_mtx(pathlib.Path(f"{stem}_matrix.mtx.gz"), pathlib.Path(f"{stem}_genes.tsv.gz"), pathlib.Path(f"{stem}_barcodes.tsv.gz")) for stem in stems]
    if len({part[0].shape[0] for part in parts}) != 1 or not all((part[1][0] == parts[0][1][0]).all() for part in parts):
        raise ValueError("sublibraries list different features")
    keys = sorted({barcode.split("-")[0] for part in parts for barcode in part[2]})
    index = pd.Index(keys)
    columns = [index.get_indexer([barcode.split("-")[0] for barcode in part[2]]) for part in parts]
    rows, cols, data = [], [], []
    for part, column in zip(parts, columns):
        coo = part[0].tocoo()
        rows.append(coo.row), cols.append(column[coo.col]), data.append(coo.data)
    summed = sparse.coo_matrix((np.concatenate(data), (np.concatenate(rows), np.concatenate(cols))), shape=(parts[0][0].shape[0], len(keys))).tocsc()
    summed.sum_duplicates()
    return summed, parts[0][1], np.array(keys, dtype=object)


libraries = []
for name, arm, control_type, stems in (("IBD", "no_vesicles", "no_vesicles", [RAW / "GSM9108990_TD_S1", RAW / "GSM9108990_TD_S2"]),
                                       ("IBD_EV", "PD-L1_miR-27a-3p_vesicles", "", [RAW / "GSM9108991_TDE_S1", RAW / "GSM9108991_TDE_merged"])):
    counts, features, barcodes = merged_capture(stems)
    libraries.append(Library(name, "series", arm, counts, features, barcodes, control_type))
build_store("GSE302696", libraries, qc={}, qc_source="deposit's filter: authors' barcode lists", source="authors' matrices, GEO GSE302696", out=store_path())
