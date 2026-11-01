"""AGO2 knockout HeLa S3 clone against parental cells, Drop-seq, one library per arm, three pseudo-replicates. Exon counts.
QC: at least 500 genes, the authors' filter (10.15252/msb.20209946)."""

from _sample_arms import DATA, Library, build_store, read_dense_counts, store_path

RAW = DATA / "GSE142277" / "extracted"
libraries = []
for path, arm, control_type in ((next(RAW.glob("*exon*WT*")), "wild_type", "wild_type"), (next(RAW.glob("*exon*KO*")), "AGO2_KO", "")):
    counts, features, barcodes = read_dense_counts(path)
    libraries.append(Library(arm, "series", arm, counts, features, barcodes, control_type))
build_store("GSE142277", libraries, qc={"min_genes": 500}, qc_source="authors, 10.15252/msb.20209946", source="authors' Drop-seq exon matrices, GEO GSE142277", out=store_path())
