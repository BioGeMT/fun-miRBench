"""miR-449a plus miR-340-5p against vehicle (lipofectamine), RD18 rhabdomyosarcoma, 60 h, single nuclei; one capture per arm,
three pseudo-replicates. QC: deposit's Cell Ranger filter."""


from _sample_arms import DATA, Library, build_store, store_path
from funmirbench.single_cell_store import read_10x_mtx

RAW = DATA / "GSE282611" / "raw"
libraries = []
for stem, arm, control_type in (("GSM8647165_GC112500", "vehicle_60h", "vehicle"), ("GSM8647166_GC112501", "miR-449a_miR-340-5p_60h", "")):
    counts, features, barcodes = read_10x_mtx(RAW / f"{stem}_matrix.mtx.gz", RAW / f"{stem}_features.tsv.gz", RAW / f"{stem}_barcodes.tsv.gz")
    libraries.append(Library(stem, "series", arm, counts, features, barcodes, control_type))
build_store("GSE282611", libraries, qc={}, qc_source="deposit's Cell Ranger filtered matrix", source="authors' matrices, GEO GSE282611", out=store_path())
