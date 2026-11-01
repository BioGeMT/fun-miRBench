"""K562 half-cell small RNA-seq and SMART-seq v4 (Wang et al. 2019). Deposit: normalised log2 miRNA fractions, RPKM."""
import pathlib, sys
sys.path.insert(0, str(pathlib.Path(__file__).parent))
import numpy as np, pandas as pd
from scipy import sparse
from cooccurrence_common import DATA, write_study

deposit = DATA / "GSE114071/deposit"
small = pd.read_csv(deposit / "GSE114071_NW_scsmRNA_K562_norm_log2.gct.gz", sep="\t")
mrna = pd.read_csv(deposit / "GSE114071_NW_half_cell_K562_RNAseq_processed_data.gct.gz", sep="\t", skiprows=2)
small_cells = {int(c.rsplit("_", 1)[1]): c for c in small.columns[2:]}
mrna_cells = {c.rsplit("_", 1)[1]: c for c in mrna.columns[2:]}
paired = sorted(n for n in small_cells if str(n) in mrna_cells)
values = small[[small_cells[n] for n in paired]].to_numpy(float)
fractions = np.where(values <= values.min() + 1e-6, 0, 2 ** values)  # the authors' floor, log2(1e-4), means not detected
names = small.Name + "|" + small.Description
genes = mrna.drop(columns="Description").drop_duplicates("NAME").set_index("NAME")
genes = genes[[mrna_cells[str(n)] for n in paired]]
n_genes_without_value = int(genes.isna().any(axis=1).sum())
genes = genes.dropna()  # the deposit has no value for these genes in some cells
cells = pd.DataFrame({"cell": [f"K562_cell_{n:02d}" for n in paired], "group": "K562", "cell_line": "K562", "deposit_small": [small_cells[n] for n in paired], "deposit_mrna": [mrna_cells[str(n)] for n in paired]})
keep = fractions.sum(axis=1) > 0
write_study("GSE114071_K562", mirna=sparse.csr_matrix(fractions[keep].T), mirna_features=names[keep], mrna=sparse.csr_matrix(genes.to_numpy(float).T), mrna_features=genes.index, cells=cells, library_totals=None, library_total_source="none: the deposit holds normalised miRNA rows only, with no small RNA library total; the mRNA is a separate library", info={
    "feature_level": "mature_arm_per_precursor", "same_cell": "half of one cell's lysate to each library (GEO design)",
    "values": "miRNA: normalised fraction from the authors' log2 table; mRNA: RPKM. No integer counts in the deposit",
    "grouping": "cell line", "qc_source": "deposit (authors' matrices)", "source": "authors, GSE114071 supplementary files",
    "caveats": f"{n_genes_without_value} genes without a value in some paired cells are dropped. Cells paired by the number in the column names (HalfCell_N with hcRNAseq_N); 19 pairs, as the authors' QC. GEO lists 20 pairs, but HalfCell_23, HalfCell_24, hcRNAseq_21_1 and hcRNAseq_21_2 have no mapping in the deposit. Rows of one mature miRNA from several precursors are near-duplicates. Raw reads SRP144659 could give integers (tier 2)."})
