"""scDASH MATQ-seq of HCT116: one total RNA library per cell, with and without rRNA depletion."""
import pathlib, sys
sys.path.insert(0, str(pathlib.Path(__file__).parent))
import pandas as pd
from scipy import sparse
from cooccurrence_common import DATA, is_mirna_symbol, write_study

table = pd.read_csv(DATA / "GSE158880/deposit/GSE158880_RawCounts_GeneFeatures.csv.gz", index_col=0)
single = [c for c in table.columns if "singlecell" in c]
table = table[single]
cells = pd.DataFrame({"cell": single, "group": ["HCT116_" + c.split("_")[0] for c in single], "cell_line": "HCT116"})
is_mirna = is_mirna_symbol(pd.Series(table.index))
write_study("GSE158880_scDASH", mirna=sparse.csr_matrix(table[is_mirna].to_numpy().T), mirna_features=table.index[is_mirna], mrna=sparse.csr_matrix(table[~is_mirna].to_numpy().T), mrna_features=table.index[~is_mirna], cells=cells, library_totals=table.sum(axis=0), library_total_source="sum of every row of the deposited table, rRNA included", info={
    "feature_level": "precursor_locus", "same_cell": "one total RNA library per cell (GEO design)", "values": "deposited raw counts, partly fractional",
    "grouping": "cell line and rRNA depletion", "qc_source": "deposit", "source": "authors, GSE158880 supplementary file",
    "caveats": "Four cells per group, below the step's 10, so groups are not correlated; per-cell tables only. Mini-bulk samples left out."})
