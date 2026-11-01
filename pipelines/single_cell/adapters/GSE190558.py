"""KSHV lytic reactivation (doxycycline 4 days) against latent cells, iSLK.219; one library per arm, three pseudo-replicates.
The latent library also holds iSLK.H cells with no per-cell call. Whole virus. QC: deposit's filter."""

from _sample_arms import DATA, Library, build_store, store_path
from funmirbench.single_cell_store import read_10x_mtx

RAW = DATA / "GSE190558" / "extracted"
libraries = []
for stem, arm, control_type in (("GSM5725695_iSLKH_iSLK219", "latent", "latent_infection"), ("GSM5725696_iSLK219_dox", "lytic_dox_4d", "")):
    counts, features, barcodes = read_10x_mtx(RAW / f"{stem}_matrix.mtx.gz", RAW / f"{stem}_features.tsv.gz", RAW / f"{stem}_barcodes.tsv.gz")
    libraries.append(Library(arm, "series", arm, counts, features, barcodes, control_type))
build_store("GSE190558", libraries, qc={}, qc_source="deposit's filter", source="authors' matrices, GEO GSE190558", out=store_path())
