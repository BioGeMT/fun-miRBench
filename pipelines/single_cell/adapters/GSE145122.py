"""22q11.2 deletion patients against control donors, iPSC cortical spheroids; two donor lines per arm are real replicates.
QC: deposit's filter (paper full text not available)."""

import pathlib

from _sample_arms import DATA, Library, build_store, store_path
from funmirbench.single_cell_store import read_10x_mtx

RAW = DATA / "GSE145122" / "extracted"
libraries = []
for gsm, sample, arm, control_type in (("GSM4306931_Control_1", "Control_1", "control", "unaffected_donor"), ("GSM4306932_Control_2", "Control_2", "control", "unaffected_donor"),
                                       ("GSM4306933_Patient_1", "Patient_1", "22q11.2del", ""), ("GSM4306934_Patient_2", "Patient_2", "22q11.2del", "")):
    stem = RAW / gsm
    counts, features, barcodes = read_10x_mtx(pathlib.Path(f"{stem}_matrix.mtx.gz"), pathlib.Path(f"{stem}_features.tsv.gz"), pathlib.Path(f"{stem}_barcodes.tsv.gz"))
    libraries.append(Library(sample, "series", arm, counts, features, barcodes, control_type, unit=sample))
build_store("GSE145122", libraries, qc={}, qc_source="deposit's filter: authors' barcode lists", source="authors' matrices, GEO GSE145122", out=store_path())
