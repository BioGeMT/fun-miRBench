"""Exosomes of human embryonic stem cells (carrying miR-302b and everything else) against PBS, senescent LO2 cells; one library
per arm, three pseudo-replicates. QC: deposit's filter."""


from _sample_arms import DATA, Library, build_store, store_path
from funmirbench.single_cell_store import read_10x_mtx

RAW = DATA / "GSE247073" / "extracted"
libraries = []
for stem, arm, control_type in (("GSM7882928_Dox_pbs", "PBS", "saline"), ("GSM7882929_Dox_Exos", "hESC_exosomes", "")):
    counts, features, barcodes = read_10x_mtx(RAW / f"{stem}_matrix.mtx.gz", RAW / f"{stem}_features.tsv.gz", RAW / f"{stem}_barcodes.tsv.gz")
    libraries.append(Library(arm, "series", arm, counts, features, barcodes, control_type))
build_store("GSE247073", libraries, qc={}, qc_source="deposit's filter", source="authors' matrices, GEO GSE247073", out=store_path())
