"""miR-18a knockout clone against wild-type A375. The deposit has one merged matrix per condition, so the two knockout and
eight wild-type libraries cannot be separated: three pseudo-replicates per arm. QC: authors (10.1002/advs.202515158)."""


from _sample_arms import DATA, Library, build_store, store_path
from funmirbench.single_cell_store import read_10x_mtx

RAW = DATA / "GSE318353" / "extracted"
libraries = []
for condition, arm, control_type in (("WT", "wild_type", "wild_type"), ("KO", "miR-18a_KO", "")):
    folder = RAW / condition / "filtered_feature_bc_matrix"
    counts, features, barcodes = read_10x_mtx(folder / "matrix.mtx.gz", folder / "features.tsv.gz", folder / "barcodes.tsv.gz")
    libraries.append(Library(condition, "series", arm, counts, features, barcodes, control_type))
build_store("GSE318353", libraries, qc={"min_genes": 200, "min_counts": 1000, "max_mito": 0.20}, qc_source="authors, 10.1002/advs.202515158",
            source="authors' per-condition Cell Ranger matrices, GEO GSE318353", out=store_path())
