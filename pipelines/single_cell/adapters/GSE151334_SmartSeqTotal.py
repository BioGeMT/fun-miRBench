"""Smart-seq-total (fibroblasts, HEK293T, MCF7; a further HEK293T UMI table). One library per cell, gene-level counts."""
import pathlib, sys
sys.path.insert(0, str(pathlib.Path(__file__).parent))
import pandas as pd
from scipy import sparse
from cooccurrence_common import DATA, is_mirna_symbol, write_study

deposit = DATA / "GSE151334/deposit"
main = pd.read_csv(deposit / "GSE151334_counts.human.tsv.gz", sep="\t", index_col=0)
umi = pd.read_csv(deposit / "GSE151334_SST-UMI-raw-1.tsv.gz", sep="\t", index_col=0).drop(columns="biotype")
umi.columns = [f"HEK293T-UMI_{c}" for c in umi.columns]
table = pd.concat([main, umi], axis=1).fillna(0)
line = pd.Series(table.columns).str.extract(r"^([^_]+)_")[0]
cells = pd.DataFrame({"cell": table.columns, "group": line, "cell_line": line.str.replace("-UMI", "")})
is_mirna = is_mirna_symbol(pd.Series(table.index))
write_study("GSE151334_SmartSeqTotal", mirna=sparse.csr_matrix(table[is_mirna].to_numpy().T), mirna_features=table.index[is_mirna], mrna=sparse.csr_matrix(table[~is_mirna].to_numpy().T), mrna_features=table.index[~is_mirna], cells=cells, library_totals=table.sum(axis=0), library_total_source="sum of every row of the deposited table of the cell's library", info={
    "feature_level": "precursor_locus", "same_cell": "one total RNA library per cell (GEO design)",
    "values": "deposited counts", "grouping": "cell line; the UMI table's wells are their own group",
    "qc_source": "deposit (cells in the authors' matrices)", "source": "authors, GSE151334 supplementary files",
    "caveats": "miRNA features are gene symbols (MIR21, hsa-mir-3687-1), so they count reads at the precursor locus, not mature arms. Velut et al. 2025 report the miRNA depth too low for target analysis. The UMI wells are HEK293T per GEO sample titles; the search's reading as stem cells does not hold (those are mouse). The 15-well tn5 comparison table is left out as a protocol test."})
