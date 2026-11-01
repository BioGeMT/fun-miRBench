"""EBV infection of primary B cells from two donors against the same donor's uninfected cells at day 0; donor is the batch and
the real replicate; time and infection are confounded and no viral miRNA is isolated. QC: deposit's Cell Ranger filter."""

import re

from _sample_arms import DATA, Library, build_store, store_path
from funmirbench.single_cell_store import read_10x_mtx

RAW = DATA / "GSE189141" / "raw"
libraries = []
for path in sorted(RAW.glob("GSM*_GEX_EBV_barcodes.tsv.gz")):
    stem = str(path).removesuffix("barcodes.tsv.gz")
    donor, day = re.search(r"(TX\d+)_(\d+)_GEX", stem).groups()
    counts, features, barcodes = read_10x_mtx(f"{stem}matrix.mtx.gz", f"{stem}features.tsv.gz", f"{stem}barcodes.tsv.gz")
    sample = f"{donor}_d{day}"
    arm, control_type = ("uninfected_d0", "uninfected") if day == "0" else (f"EBV_d{day}", "")
    libraries.append(Library(sample, donor, arm, counts, features, barcodes, control_type, unit=sample))
build_store("GSE189141", libraries, qc={}, qc_source="deposit's Cell Ranger filter", source="authors' matrices with the EBV genome, GEO GSE189141", out=store_path())
