"""miR-125b promoter EGFP reporter in MDA-MB-231: top and bottom 5% sorted, 12 and 9 libraries. Treated as single cells on the
GEO wording, not proven (GEO also says 3 ug RNA per sample). Arms split into three pseudo-replicates of whole libraries;
the libraries of one sort are not independent. QC: all deposited libraries."""

import pandas as pd
from _sample_arms import DATA, Library, build_store, read_dense_counts, store_path
import numpy as np

RAW = DATA / "GSE123358" / "raw"
counts, features, barcodes = read_dense_counts(RAW / "GSE123358_ReadCounts.txt.gz")
names = pd.Series(barcodes)
high = np.flatnonzero(names.str.fullmatch(r"H\d+").to_numpy())
low = np.flatnonzero(names.str.fullmatch(r"L\d+").to_numpy())
if len(high) + len(low) != len(barcodes):
    raise ValueError(f"columns not labelled high or low: {barcodes[:5]}")
libraries = [Library("high", "sort", "reporter_high", counts[:, high], features, barcodes[high]),
             Library("low", "sort", "reporter_low", counts[:, low], features, barcodes[low], "low_reporter_sort")]
build_store("GSE123358", libraries, qc={}, qc_source="all libraries in the deposit's table", source="authors' read count table, GEO GSE123358", out=store_path())
