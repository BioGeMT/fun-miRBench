"""miR-9/9*-124 conversion of adult fibroblasts, Coriell GM02171: miR only (GSE154359) and with ISL1 and LHX3 (GSE154356,
GSE154358), each day against day-0 fibroblasts; no same-day neutral arm, so time and treatment are confounded. One library
per day, three pseudo-replicates. GSE154356 is one merged table whose file lacks its last 3 entries (header 54,672,111, file 54,672,108), all in the last cell; read as
it stands. QC: authors, GEO processing: at least 200 genes; their mitochondrial cut-off has no number, so it is not applied."""

import numpy as np
import pandas as pd
import pathlib
from scipy import sparse

from _sample_arms import DATA, Library, build_store, store_path
from funmirbench.single_cell_store import read_10x_mtx

RAW = DATA
libraries = []
for path in sorted((RAW / "GSE154359" / "raw").glob("GSM*_barcodes.tsv.gz")):
    stem = str(path).removesuffix("_barcodes.tsv.gz")
    day = stem.rsplit("day", 1)[1]
    counts, features, barcodes = read_10x_mtx(pathlib.Path(f"{stem}_matrix.mtx.gz"), pathlib.Path(f"{stem}_features.tsv.gz"), pathlib.Path(f"{stem}_barcodes.tsv.gz"))
    libraries.append(Library(f"miR_only_d{day}", "all", f"miR_only_d{day}", counts, features, barcodes, "day0_untransduced" if day == "0" else ""))
merged = RAW / "GSE154356" / "raw"
table = pd.read_csv(merged / "GSE154356_IL3_D0toD20_matrix.mtx.gz", sep=" ", comment="%", header=None, dtype=np.int64)
header, entries = table.iloc[0].tolist(), table.iloc[1:]
print(f"GSE154356: header lists {header[2]} entries, file holds {len(entries)}", flush=True)
counts = sparse.coo_matrix((entries[2], (entries[0] - 1, entries[1] - 1)), shape=(header[0], header[1])).tocsc()
features = pd.read_csv(merged / "GSE154356_IL3_D0toD20_genes.tsv.gz", sep="\t", header=None, dtype=str)
barcodes = pd.read_csv(merged / "GSE154356_IL3_D0toD20_barcodes.tsv.gz", sep="\t", header=None, dtype=str)[0].to_numpy()
prefix = np.array([barcode.split("_")[0] for barcode in barcodes])
for name, day in (("Fibroblast", "0"), ("D5.IL3", "5"), ("D10.IL3", "10"), ("D15.IL3", "15"), ("D20.IL3", "20")):
    wanted = np.flatnonzero(prefix == name)
    libraries.append(Library(f"miR_I1L3_d{day}", "all", f"miR_I1L3_d{day}", counts[:, wanted], features, barcodes[wanted], "day0_untransduced" if day == "0" else ""))
for path in sorted((RAW / "GSE154358" / "raw").glob("GSM*_barcodes.tsv.gz")):
    stem = str(path).removesuffix("_barcodes.tsv.gz")
    day = stem.rsplit("day", 1)[1]
    counts, features, barcodes = read_10x_mtx(pathlib.Path(f"{stem}_matrix.mtx.gz"), pathlib.Path(f"{stem}_features.tsv.gz"), pathlib.Path(f"{stem}_barcodes.tsv.gz"))
    libraries.append(Library(f"miR_I1L3_rep2_d{day}", "all", f"miR_I1L3_rep2_d{day}", counts, features, barcodes))
build_store("GSE154178", libraries, qc={"min_genes": 200}, qc_source="authors, GEO data processing: at least 200 genes; mitochondrial cut-off not stated", source="authors' matrices, GEO GSE154356, GSE154358, GSE154359", out=store_path())
