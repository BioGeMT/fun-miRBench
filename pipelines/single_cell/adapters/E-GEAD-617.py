"""let-7 switch mRNA in imMKCL: cells sorted by switch output into let-7 high and let-7 low, one capture each, three
pseudo-replicates. Counts and cells from the authors' Seurat object, exported to matrix.mtx.gz, genes.txt, cells.txt with R:
cells in it passed their filter, stated in the record as at most 15% mitochondrial and at least 2,500 genes."""

import numpy as np
import pandas as pd
from scipy import io as scipy_io
from scipy import sparse
from _sample_arms import DATA, Library, build_store, store_path

RAW = DATA / "E-GEAD-617" / "extracted"
counts = sparse.csc_matrix(scipy_io.mmread(RAW / "matrix.mtx.gz"))
genes = pd.DataFrame({0: pd.read_csv(RAW / "genes.txt", header=None, dtype=str)[0]})
cells = pd.read_csv(RAW / "cells.txt", header=None, dtype=str)[0].to_numpy()
group = np.array(["_".join(cell.split("_")[:2]) for cell in cells])
libraries = []
for name, arm, control_type in (("Let7_Low", "let-7_low", "low_reporter_sort"), ("Let7_High", "let-7_high", "")):
    wanted = np.flatnonzero(group == name)
    libraries.append(Library(name, "transfection", arm, counts[:, wanted], genes, cells[wanted], control_type))
build_store("E-GEAD-617", libraries, qc={}, qc_source="authors' filter, record E-GEAD-617: over 15% mitochondrial or under 2,500 genes removed before the object", source="authors' Seurat object, DDBJ GEA E-GEAD-617", out=store_path())
