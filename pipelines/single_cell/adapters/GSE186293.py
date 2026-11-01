"""Tetralogy of Fallot with 22q11.2 deletion (TOF-DG) against healthy lines and TOF without deletion (TOF-ND), hiPSC cardiac
progenitors on days 7, 8, 11 (day is the batch; one TOF-DG line) and cardiomyocytes on day 22. QC: authors
(10.1038/s42003-023-05344-6: at most 20 to 25% mitochondrial, I use 25%; at least 500 transcripts)."""

import pathlib

from _sample_arms import DATA, Library, build_store, store_path
from funmirbench.single_cell_store import read_10x_mtx

RAW = DATA / "GSE186293" / "extracted"
GROUP = {"Control": ("healthy", "healthy_line"), "TOF-DG": ("TOF-DG", ""), "TOF-ND": ("TOF-ND", "no_deletion_disease_line")}
libraries = []
for path in sorted(RAW.glob("GSM*_barcodes.tsv.gz")):
    stem = str(path).removesuffix("_barcodes.tsv.gz")
    line, stage = stem.rsplit("/", 1)[1].split("_", 1)[1].rsplit("_", 1)
    group, control_type = GROUP[line.rstrip("0123456789")]
    counts, features, barcodes = read_10x_mtx(pathlib.Path(f"{stem}_matrix.mtx.gz"), pathlib.Path(f"{stem}_features.tsv.gz"), pathlib.Path(f"{stem}_barcodes.tsv.gz"))
    arm = f"{group}_{'CM' if stage == 'CAS' else 'CP'}"
    libraries.append(Library(f"{line}_{stage}", stage, arm, counts, features, barcodes, control_type, unit=f"{line}_{stage}"))
build_store("GSE186293", libraries, qc={"min_counts": 500, "max_mito": 0.25}, qc_source="authors, 10.1038/s42003-023-05344-6", source="authors' matrices, GEO GSE186293", out=store_path())
