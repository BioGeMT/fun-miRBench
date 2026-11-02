"""miR-98-5p mimic against control mimic, H9 hESC during resetting to naive pluripotency. Days 2 and 3 share an arm name
and stand as two paired time points with the day as batch (GEO labels each rep1, so independent cultures are not established); days 8 and 14 are their own comparisons; day 5 has no control.
QC: authors' thresholds (preprint 10.1101/2025.10.14.682465)."""

from _sample_arms import DATA, Library, build_store, store_path
from funmirbench.single_cell_store import read_10x_mtx

RAW = DATA / "GSE280065" / "matrices"
SAMPLES = [
    # file stem, sample, batch, arm, control type, replicate unit
    ("GSM8586813_ctrl_d02", "ctrl_d02", "d02", "control_d02_d03", "control_mimic", "ctrl_d02"),
    ("GSM8586814_ctrl_d03", "ctrl_d03", "d03", "control_d02_d03", "control_mimic", "ctrl_d03"),
    ("GSM8586818_mimic_d02", "mimic_d02", "d02", "mimic_d02_d03", "", "mimic_d02"),
    ("GSM8586819_mimic_d03", "mimic_d03", "d03", "mimic_d02_d03", "", "mimic_d03"),
    ("GSM8586820_mimic_d05", "mimic_d05", "d05", "mimic_d05", "", None),
    ("GSM8586815_ctrl_d08", "ctrl_d08", "d08", "control_d08", "control_mimic", None),
    ("GSM8586821_mimic_d08", "mimic_d08", "d08", "mimic_d08", "", None),
    ("GSM8586816_ctrl_d14_1", "ctrl_d14_1", "d14", "control_d14", "control_mimic", "ctrl_d14_1"),
    ("GSM8586817_ctrl_d14_2", "ctrl_d14_2", "d14", "control_d14", "control_mimic", "ctrl_d14_2"),
    ("GSM8586822_mimic_d14_1", "mimic_d14_1", "d14", "mimic_d14", "", "mimic_d14_1"),
    ("GSM8586823_mimic_d14_2", "mimic_d14_2", "d14", "mimic_d14", "", "mimic_d14_2"),
]

libraries = []
for stem, sample, batch, arm, control_type, unit in SAMPLES:
    counts, features, barcodes = read_10x_mtx(RAW / f"{stem}_matrix.mtx.gz", RAW / f"{stem}_features.tsv.gz", RAW / f"{stem}_barcodes.tsv.gz")
    libraries.append(Library(sample, batch, arm, counts, features, barcodes, control_type, unit))
build_store("GSE280065", libraries, qc={"min_genes": 2000, "min_counts": 5000, "max_counts": 100000, "max_mito": 0.15},
            qc_source="authors, preprint 10.1101/2025.10.14.682465", source="authors' Cell Ranger matrices, GEO GSE280065",
            out=store_path())
