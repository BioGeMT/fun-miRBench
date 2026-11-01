"""HSV-1 wild type and ICP0 deletion mutant against mock infection, neonatal fibroblasts, 5 h, Drop-seq; one library per arm,
three pseudo-replicates. The whole virus acts; no viral miRNA is isolated. QC: deposit's filter."""

from _sample_arms import DATA, Library, build_store, read_dense_counts, store_path

RAW = DATA / "GSE126042" / "raw"
libraries = []
for name, arm, control_type in (("mock", "mock", "mock"), ("wt", "HSV-1_wt", ""), ("d0", "HSV-1_dICP0", "")):
    counts, features, barcodes = read_dense_counts(RAW / f"GSE126042_{name}_merged3.dge.txt.gz")
    libraries.append(Library(name, "series", arm, counts, features, barcodes, control_type))
build_store("GSE126042", libraries, qc={}, qc_source="deposit's filter", source="authors' Drop-seq tables, GEO GSE126042", out=store_path())
