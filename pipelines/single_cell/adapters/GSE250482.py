"""Isogenic 22q11.2 deletion against wild-type MEL-1, forebrain organoids day 60, two independent differentiations per arm.
The 15q13.3 arms are another deletion and are left out. QC: deposit's filter (paper full text not available)."""

from _sample_arms import DATA, Library, build_store, read_dense_counts, store_path

RAW = DATA / "GSE250482" / "extracted"
libraries = []
for gsm, sample, arm, control_type in (("GSM7979673_WT_d60_1", "WT_d60_1", "control_d60", "isogenic_control"), ("GSM7979675_WT_d60_2", "WT_d60_2", "control_d60", "isogenic_control"),
                                       ("GSM7979674_22q11_del_d60_1", "del_d60_1", "22q11.2del_d60", ""), ("GSM7979676_22q11_del_d60_2", "del_d60_2", "22q11.2del_d60", "")):
    counts, features, barcodes = read_dense_counts(RAW / f"{gsm}.csv.gz", sep=",")
    libraries.append(Library(sample, "d60", arm, counts, features, barcodes, control_type, unit=sample))
build_store("GSE250482", libraries, qc={}, qc_source="deposit's filter: cells in the authors' tables", source="authors' count tables, GEO GSE250482", out=store_path())
