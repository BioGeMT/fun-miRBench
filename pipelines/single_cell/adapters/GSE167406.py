"""miR-200b-3p hairpin inhibitor against control inhibitor, adult dermal fibroblasts, by day; one library per arm and day,
three pseudo-replicates. Controls exist for days 1 and 7 only. QC: authors (10.1038/s41467-023-36665-z)."""

from _sample_arms import DATA, Library, build_store, store_path
from funmirbench.single_cell_store import read_10x_mtx

RAW = DATA / "GSE167406" / "matrices"
SAMPLES = [
    ("GSM5104820_Untreated_HADF", "GSM5104820_Unreated_HADF", "untreated", "d0", "untreated_d0", "untreated"),
    ("GSM5104821_MI-Treated_HADF_D1", None, "inhibitor_d1", "d1", "inhibitor_d1", ""),
    ("GSM5104822_MI-Treated_HADF_D3", None, "inhibitor_d3", "d3", "inhibitor_d3", ""),
    ("GSM5104823_MI-Treated_HADF_D5", None, "inhibitor_d5", "d5", "inhibitor_d5", ""),
    ("GSM5104824_MI-Treated_HADF_D7", None, "inhibitor_d7", "d7", "inhibitor_d7", ""),
    ("GSM6911548_CI-Treated_HADF_D1", None, "control_inhibitor_d1", "d1", "control_inhibitor_d1", "control_inhibitor"),
    ("GSM6911549_CI-Treated_HADF_D7", None, "control_inhibitor_d7", "d7", "control_inhibitor_d7", "control_inhibitor"),
]
libraries = []
for stem, barcode_stem, sample, batch, arm, control_type in SAMPLES:
    counts, features, barcodes = read_10x_mtx(RAW / f"{stem}_matrix.mtx.gz", RAW / f"{stem}_features.tsv.gz", RAW / f"{barcode_stem or stem}_barcodes.tsv.gz")
    libraries.append(Library(sample, batch, arm, counts, features, barcodes, control_type))
build_store("GSE167406", libraries, qc={"min_genes": 200, "max_genes": 8000, "min_counts": 500, "max_counts": 60000, "max_mito": 0.10},
            qc_source="authors, 10.1038/s41467-023-36665-z", source="authors' Cell Ranger matrices, GEO GSE167406", out=store_path())
