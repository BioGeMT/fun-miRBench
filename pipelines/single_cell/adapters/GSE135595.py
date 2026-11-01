"""miR-122 and miR-200a knockout, let-7i and miR-452 overexpression and both together (MiMAC) against empty vector in hiPSC
cardiomyocytes, glucose media. The series table is log-normalised (scale 10,000), so raw counts are rebuilt with the
metadata's nUMI. Single series table; cells listed in the authors' metadata passed their QC
(10.1038/s41467-019-12482-1: at most 40% mitochondrial, at least 200 genes, at least 2,000 UMIs). The fatty-acid-media arms
have no empty vector and are left out. One library per arm, three pseudo-replicates."""

import numpy as np
import pandas as pd
from _sample_arms import DATA, Library, build_store, read_dense_counts, store_path
from funmirbench.single_cell_store import raw_counts_from_lognorm

RAW = DATA / "GSE135595" / "raw"
ARMS = {"EV": ("EV", "empty_vector"), "122_200": ("miR-122_miR-200a_KO", ""), "let7_452": ("let-7i_miR-452_OE", ""), "MiMAC": ("MiMAC", "")}
meta = pd.read_csv(RAW / "GSE135595_cell_metadata.txt.gz", sep="\t")
counts, features, barcodes = read_dense_counts(RAW / "GSE135595_scRNAseq_counts.txt.gz")
# the two files write the arm prefix differently, so cells are matched by position and checked on the last 16 bases
if not np.array_equal(meta.cell_barcode.str[-16:].to_numpy(), np.array([name[-16:] for name in barcodes])):
    raise ValueError("metadata and count table list cells in a different order")
counts = raw_counts_from_lognorm(counts.T, meta.nUMI.to_numpy()).T.tocsc()
detected = np.diff(counts.indptr)
if not np.array_equal(detected, meta.nGene.to_numpy()):
    raise ValueError("rebuilt counts do not reproduce nGene of the metadata")
libraries = []
for orig, (arm, control_type) in ARMS.items():
    wanted = np.flatnonzero(meta["orig.ident"].to_numpy() == orig)
    libraries.append(Library(orig, "series", arm, counts[:, wanted], features, barcodes[wanted], control_type))
build_store("GSE135595", libraries, qc={}, qc_source="authors' metadata (10.1038/s41467-019-12482-1): only cells that passed their filter are listed",
            source="authors' log-normalised series table, raw counts rebuilt with nUMI, GEO GSE135595", out=store_path())
