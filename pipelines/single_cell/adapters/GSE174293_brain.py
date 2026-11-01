"""MATQ-Drop total RNA-seq of nuclei and synapse-like compartments from frozen human brain (two donors, two regions)."""
import pathlib, re, sys
sys.path.insert(0, str(pathlib.Path(__file__).parent))
import pandas as pd
from scipy import sparse
from cooccurrence_common import DATA, mirna_ensembl_ids, write_study

parts, rows = [], []
for path in sorted((DATA / "GSE174293/deposit").glob("GSE174293_P*_transcript_count.tsv.gz")):
    sample = re.match(r"GSE174293_(P\d+_\w+?_(?:nucleus|synapse))_transcript", path.name).group(1)
    table = pd.read_csv(path, sep="\t", index_col=0)
    table.index = table.index.str.replace(r"\.\d+$", "", regex=True)
    table = table.groupby(level=0).sum()
    table.columns = [f"{sample}_{c}" for c in table.columns]
    parts.append(table)
    rows += [(c, sample, *sample.split("_")) for c in table.columns]
table = pd.concat(parts, axis=1).fillna(0)
cells = pd.DataFrame(rows, columns=["cell", "group", "donor", "region", "compartment"])
is_mirna = table.index.isin(mirna_ensembl_ids())
write_study("GSE174293_brain", mirna=sparse.csr_matrix(table[is_mirna].to_numpy().T), mirna_features=table.index[is_mirna], mrna=sparse.csr_matrix(table[~is_mirna].to_numpy().T), mrna_features=table.index[~is_mirna], cells=cells, library_totals=table.sum(axis=0), library_total_source="sum of every row of the deposited transcript_count table", info={
    "feature_level": "precursor_locus", "same_cell": "one total RNA library per nucleus or compartment (GEO design), not whole cells", "values": "deposited counts",
    "grouping": "donor, region and compartment (no cell-type labels in the deposit)", "qc_source": "deposit (barcodes in the authors' matrices)", "source": "authors, GSE174293 supplementary files",
    "caveats": "Nuclei hold nascent transcripts, not mature cytoplasmic miRNA. Groups mix cell types, so a correlation can reflect cell-type composition."})
