"""KSHV against mock infection, oral epithelial 3D tissues, 3 and 6 days, three experiment sets. The set is the batch and the
real replicate; the whole virus changes the transcriptome, no viral miRNA is isolated. QC: deposit's filter."""

import pathlib
import re

from _sample_arms import DATA, Library, build_store, store_path
from funmirbench.single_cell_store import read_10x_mtx

RAW = DATA / "GSE184257" / "extracted"
TITLES = dict(line.split("\t")[:2] for line in open(DATA / "GSE184257" / "geo_record" / "summary.txt") if line.startswith("GSM"))
libraries = []
for gsm, title in TITLES.items():
    infection, day, batch = re.fullmatch(r"(Mock|KSHV) (\ddpi) Set (\d)", title).groups()
    stem = str(next(RAW.glob(f"{gsm}*barcodes*"))).split("barcodes")[0]
    genes = f"{stem}genes.tsv.gz" if pathlib.Path(f"{stem}genes.tsv.gz").exists() else f"{stem}features.tsv.gz"
    counts, features, barcodes = read_10x_mtx(f"{stem}matrix.mtx.gz", genes, f"{stem}barcodes.tsv.gz")
    sample = f"{infection}_{day}_set{batch}"
    libraries.append(Library(sample, f"set{batch}", f"{infection}_{day}", counts, features, barcodes, "mock" if infection == "Mock" else "", unit=sample))
build_store("GSE184257", libraries, qc={}, qc_source="deposit's filter", source="authors' matrices, GEO GSE184257", out=store_path())
