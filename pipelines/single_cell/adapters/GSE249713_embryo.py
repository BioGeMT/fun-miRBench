"""Co-seq of human preimplantation embryos: Small-seq and Smart-seq2 from the same cell, matched on the sample title."""
import pathlib, sys
sys.path.insert(0, str(pathlib.Path(__file__).parent))
import pandas as pd
from scipy import sparse
from cooccurrence_common import DATA, write_study

deposit = DATA / "GSE249713/deposit"
small = pd.read_csv(deposit / "GSE249708_smallseq.counts.tsv.gz", sep="\t")
all_small = small.drop(columns=["ID", "type"]).sum(axis=0)  # every small RNA biotype of the Small-seq library
small = small[small.type == "miRNA"].set_index("ID").drop(columns="type")
mrna = pd.read_csv(deposit / "GSE249493_coseq.smt2.scRNAseq.counts.tsv.gz", sep="\t", index_col=0)
common = [c for c in mrna.columns if c in set(small.columns)]
day = pd.Series(common).str.split("_").str[0]
cells = pd.DataFrame({"cell": common, "group": day, "embryonic_day": day, "embryo": pd.Series(common).str.split("_").str[1]})
print(small.shape, mrna.shape, len(common), day.value_counts().to_dict(), flush=True)
write_study("GSE249713_embryo", mirna=sparse.csr_matrix(small[common].to_numpy().T), mirna_features=small.index, mrna=sparse.csr_matrix(mrna[common].to_numpy().T), mrna_features=mrna.index, cells=cells, library_totals=all_small[common], library_total_source="sum of every row of the Small-seq table (miRNA, piRNA, tRNA, snRNA, snoRNA, rRNA and others)", info={
    "feature_level": "mature_arm", "same_cell": "co-seq: sample titles '<cell> (coseq-smallseq)' and '<cell> (coseq-mRNA)' (GSE249708, GSE249493)",
    "values": "deposited counts", "grouping": "embryonic day", "qc_source": "deposit (cells in the authors' matrices)", "source": "authors, GSE249713 subseries",
    "caveats": "Embryo and lineage vary within a day; no lineage labels in the deposit, so a correlation can reflect them. The abstract the Sol search found reports 189 paired cells; the public files pair fewer."})
