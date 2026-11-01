"""miR-874-3p precursor mimic against negative-control miRNA, MCF-7, 36 h; one library per arm, three pseudo-replicates.
The siRNA arms of the series are other perturbations and are not stored. QC: authors (10.1038/s41598-022-23205-w)."""

from _sample_arms import DATA, Library, build_store, store_path
from funmirbench.single_cell_store import read_10x_mtx

RAW = DATA / "GSE214508" / "matrices"
libraries = []
for folder, arm, control_type in (("MCF7_control_mir1", "control_mirna", "control_mimic"), ("MCF7_mir874", "miR-874-3p_mimic", "")):
    counts, features, barcodes = read_10x_mtx(RAW / folder / "matrix.mtx.gz", RAW / folder / "features.tsv.gz", RAW / folder / "barcodes.tsv.gz")
    libraries.append(Library(folder, "series", arm, counts, features, barcodes, control_type))
build_store("GSE214508", libraries, qc={"max_mito": 0.10}, qc_source="authors, 10.1038/s41598-022-23205-w: more than 10% mitochondrial removed; deposit's Cell Ranger filter",
            source="authors' Cell Ranger matrices, GEO GSE214508", out=store_path())
