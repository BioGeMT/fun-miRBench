"""BK polyomavirus phases against control kidney allograft biopsies; each biopsy is a real replicate, arm = the record's bk phase.
Between patients, whole virus. Raw reads are withheld; processed matrices on GEO. QC: authors, GEO processing: at most 50% mitochondrial."""

import re

from _sample_arms import DATA, Library, build_store, store_path
from funmirbench.single_cell_store import read_10x_mtx

RAW = DATA / "GSE317012" / "extracted"
libraries = []
for line in open(DATA / "GSE317012" / "geo_record" / "summary.txt"):
    if not line.startswith("GSM"):
        continue
    fields = line.rstrip("\n").split("\t")
    gsm, title, phase = fields[0], fields[1], re.search(r"bk phase: ([^;]+)", fields[4]).group(1)
    stem = str(next(path for path in RAW.glob(f"{gsm}_*barcodes.tsv.gz") if "_raw_" not in path.name)).removesuffix("barcodes.tsv.gz")
    counts, features, barcodes = read_10x_mtx(f"{stem}matrix.mtx.gz", f"{stem}features.tsv.gz", f"{stem}barcodes.tsv.gz")
    arm = phase.replace(" ", "_")
    libraries.append(Library(title, "series", arm, counts, features, barcodes, "uninfected_biopsy" if phase == "Control" else "", unit=title))
build_store("GSE317012", libraries, qc={"max_mito": 0.50}, qc_source="authors, GEO GSE317012 data processing", source="authors' matrices, GEO GSE317012", out=store_path())
