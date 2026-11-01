"""CRISPRa reprogramming with the miR-302/367 promoter guide added (CRISPRa+ME) against the same background without it
(CRISPRa+E), LCL-derived, early stage and day 15. Same sequencing batch label; one library per arm, three
pseudo-replicates. The deposit is not cell-filtered. QC: authors, GEO processing: at least 2,000 genes, at most 50% mitochondrial."""


from _sample_arms import DATA, Library, build_store, store_path
from funmirbench.single_cell_store import read_10x_mtx

RAW = DATA / "GSE162528" / "raw"
libraries = []
for stem, arm, control_type in (("GSM4953914_CRISPRa_E_early", "CRISPRa_E_early", "same_reprogramming_background"), ("GSM4953917_CRISPRa_ME_early", "CRISPRa_ME_early", ""),
                                ("GSM4953916_CRISPRa_E_d15", "CRISPRa_E_d15", "same_reprogramming_background"), ("GSM4953918_CRISPRa_ME_d15", "CRISPRa_ME_d15", "")):
    counts, features, barcodes = read_10x_mtx(RAW / f"{stem}_matrix.mtx.gz", RAW / f"{stem}_features.tsv.gz", RAW / f"{stem}_barcodes.tsv.gz")
    libraries.append(Library(arm, "series", arm, counts, features, barcodes, control_type))
build_store("GSE162528", libraries, qc={"min_genes": 2000, "max_mito": 0.50}, qc_source="authors, GEO GSE162528 data processing", source="authors' matrices, GEO GSE162528", out=store_path())
