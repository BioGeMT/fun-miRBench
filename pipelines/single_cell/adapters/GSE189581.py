"""HCMV TB40/E against mock infection, monocytes and monocyte-derived macrophages, 5 h and 144 h (the times with a mock arm);
MARS-seq wells of experiment 1, one plate per arm, three pseudo-replicates. The whole virus acts. QC: deposit's filter:
wells listed by the authors."""

import numpy as np
import pandas as pd
from scipy import sparse
from _sample_arms import DATA, Library, build_store, store_path

RAW = DATA / "GSE189581" / "raw"
cells = pd.read_csv(RAW / "GSE189581_exp1.cells.tsv.gz", sep="\t")
reads = pd.read_csv(RAW / "GSE189581_exp1.reads.tsv.gz", sep="\t")
features = pd.DataFrame({0: reads.Ensembl_ID.astype(str), 1: reads.gene_name.astype(str)})
libraries = []
for (cell_type, infection, hours), group in cells.groupby(["Cell_type", "Infection", "time (hr)"]):
    if hours not in (5, 144):
        continue
    wells = group.Well_ID.tolist()
    counts = sparse.csc_matrix(reads[wells].to_numpy())
    name = f"{cell_type}_{infection}_{hours}h"
    libraries.append(Library(name, f"{cell_type}_{hours}h", name, counts, features, np.array(wells, dtype=object), "mock" if infection == "Mock" else ""))
build_store("GSE189581", libraries, qc={}, qc_source="deposit's filter: wells listed by the authors", source="authors' MARS-seq tables, GEO GSE189581", out=store_path())
