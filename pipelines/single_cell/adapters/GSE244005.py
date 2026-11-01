"""22q11.2 deletion patients against sex- and age-matched unaffected lines, cerebral organoids at 70 and 150 days; three lines
per arm are real replicates. QC: deposit's Cell Ranger filter; the paper's text gives no thresholds I could read."""

from _sample_arms import DATA, Library, build_store, store_path
from funmirbench.single_cell_store import read_10x_h5

RAW = DATA / "GSE244005" / "extracted"
SAMPLES = [("GSM7804653", "Q2", "WT", 70), ("GSM7804654", "Q1", "del", 70), ("GSM7804655", "Q5", "WT", 70), ("GSM7804656", "Q6", "del", 70), ("GSM7804657", "QR20", "WT", 70), ("GSM7804658", "QR27", "del", 70),
           ("GSM7804659", "Q2", "WT", 150), ("GSM7804660", "Q1", "del", 150), ("GSM7804661", "Q5", "WT", 150), ("GSM7804662", "Q6", "del", 150), ("GSM7804663", "QR20", "WT", 150), ("GSM7804664", "QR27", "del", 150)]
libraries = []
for gsm, line, genotype, day in SAMPLES:
    counts, features, barcodes = read_10x_h5(next(RAW.glob(f"{gsm}_*/**/filtered_feature_bc_matrix.h5")))
    arm = f"{'22q11.2del' if genotype == 'del' else 'control'}_DIV{day}"
    libraries.append(Library(f"{line}_DIV{day}", f"DIV{day}", arm, counts, features, barcodes,
                             "" if genotype == "del" else "unaffected_donor", unit=f"{line}_DIV{day}"))
build_store("GSE244005", libraries, qc={}, qc_source="deposit's Cell Ranger filtered matrix", source="authors' Cell Ranger outputs, GEO GSE244005", out=store_path())
