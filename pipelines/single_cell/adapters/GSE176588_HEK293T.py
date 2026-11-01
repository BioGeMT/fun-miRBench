"""VASA-seq, human arms only: HEK293T in the VASA-plate and VASA-drop mixes with mouse cells."""
import pathlib, sys
sys.path.insert(0, str(pathlib.Path(__file__).parent))
import pandas as pd
from scipy import sparse
from cooccurrence_common import DATA, write_study

MIN_HUMAN_SHARE = 0.9
root = DATA / "GSE176588/extracted"
parts, rows, library_totals = [], [], []
for assay, path in (("plate", root / "GSM5369495_vasaplate_HEK293T-mESC_split_total.UFICounts.tsv.gz"), ("drop", root / "GSM5369496_vasadrop_mm_hs_NextSeq-75bp_split_total.UFICounts.tsv.gz")):
    table = pd.read_csv(path, sep="\t", index_col=0)
    human = table.index.str.startswith("ENSG")
    totals, human_totals = table.sum(axis=0), table[human].sum(axis=0)
    share = human_totals / totals.where(totals > 0)
    chosen = share.index[share >= MIN_HUMAN_SHARE]
    print(assay, table.shape, "human cells", len(chosen), "of", table.shape[1], "total UFI quartiles of chosen", totals[chosen].quantile([.25, .5, .75]).tolist(), flush=True)
    part = table.loc[human, chosen]
    part.columns = [f"{assay}_{c}" for c in part.columns]
    parts.append(part)
    library_totals.append(totals[chosen].rename(lambda c: f"{assay}_{c}"))
    rows += [(f"{assay}_{c}", f"HEK293T_{assay}", "HEK293T", assay, float(share[c])) for c in chosen]
table = pd.concat(parts, axis=1).fillna(0)
cells = pd.DataFrame(rows, columns=["cell", "group", "cell_line", "assay", "human_share"])
n_rows_several_genes = int((table.index.str.count("_") != 2).sum())
table = table[table.index.str.count("_") == 2]  # one gene per row; the other rows hold reads shared between genes
is_mirna = table.index.str.endswith("_miRNA")
write_study("GSE176588_HEK293T", mirna=sparse.csr_matrix(table[is_mirna].to_numpy().T), mirna_features=table.index[is_mirna], mrna=sparse.csr_matrix(table[~is_mirna].to_numpy().T), mrna_features=pd.Series(table.index[~is_mirna]).str.split("_").str[0], cells=cells, library_totals=pd.concat(library_totals), library_total_source="sum of every row of the deposited table for the cell, including rows of reads shared between genes, tRNA rows and mouse rows", info={
    "feature_level": "precursor_locus", "same_cell": "one total RNA library per cell (GEO design)", "values": "deposited UFI counts, single-gene rows only",
    "grouping": "assay (plate, drop); all human cells are HEK293T", "qc_source": "no cell filter in the deposit: human cells called by at least 90% of counts on human genes (my rule, set before any correlation)", "source": "authors, GSE176588 RAW",
    "caveats": f"Only rows of one gene are used; {n_rows_several_genes} rows of reads shared between genes are left out. Every well or droplet with a human share of 90% or more is kept, including low-depth ones, since neither the deposit nor the record gives a filter. Random-primed total RNA, so reads at miRNA loci can be precursor or primary transcript."})
