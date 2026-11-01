"""SMARTer total RNA-seq per cell, cell lines with and without treatment. Deposit: kallisto per transcript per well."""
import pathlib, re, sys
sys.path.insert(0, str(pathlib.Path(__file__).parent))
import numpy as np, pandas as pd
from scipy import sparse
from cooccurrence_common import DATA, SCRATCH, gtf_tables, mirna_ensembl_ids, write_study

rows = []
for block in (SCRATCH / "geo/GSE119984.gsm.soft").read_text().split("^SAMPLE")[1:]:
    traits = dict(x.split(": ", 1) for x in re.findall(r"!Sample_characteristics_ch1 = (.*)", block) if ": " in x)
    gsm = re.search(r"!Sample_geo_accession = (\S+)", block).group(1)
    rows.append((gsm, traits.get("cell line"), traits.get("treatment", ""), traits.get("mrna seq", ""), traits.get("library preparation kit", "")))
samples = pd.DataFrame(rows, columns=["gsm", "cell_line", "treatment", "library", "kit"])
samples = samples[samples.library != "polyA+"]  # poly(A) libraries cannot hold miRNA
_, transcripts = gtf_tables()
gene_of = transcripts.set_index("transcript_id").gene_id
files = {p.name.split("_")[0]: p for p in (DATA / "GSE119984/extracted").glob("*_abundance.tsv.gz")}
samples = samples[samples.gsm.isin(files)].reset_index(drop=True)
columns, unmapped, totals = [], [], {}
for gsm in samples.gsm:
    table = pd.read_csv(files[gsm], sep="\t", usecols=["target_id", "est_counts"])
    table["target_id"] = table.target_id.str.replace(r"\.\d+$", "", regex=True)
    known = table.target_id.isin(gene_of.index)
    totals[gsm] = float(table.est_counts.sum())
    unmapped.append(float(table.est_counts[~known].sum() / max(table.est_counts.sum(), 1)))
    columns.append(table[known].groupby(table.target_id[known].map(gene_of)).est_counts.sum().rename(gsm))
table = pd.concat(columns, axis=1).fillna(0)
cells = samples.rename(columns={"gsm": "cell"})
cells["group"] = cells.cell_line + "_" + cells.treatment.where(cells.treatment != "", cells.kit)
is_mirna = table.index.isin(mirna_ensembl_ids())
write_study("GSE119984_neuroblastoma", mirna=sparse.csr_matrix(table[is_mirna].to_numpy().T), mirna_features=table.index[is_mirna], mrna=sparse.csr_matrix(table[~is_mirna].to_numpy().T), mrna_features=table.index[~is_mirna], cells=cells, library_totals=pd.Series(totals), library_total_source="sum of est_counts over every transcript in the well's kallisto file, including those absent from Ensembl 115", info={
    "feature_level": "precursor_locus", "same_cell": "one total RNA library per well (GEO)", "values": "kallisto est_counts summed per Ensembl 115 gene, fractional",
    "grouping": "cell line and treatment (kit for A375 and Jurkat)", "qc_source": "deposit (all total RNA wells)", "source": "authors, GSE119984 RAW",
    "caveats": f"Transcriptome pseudoalignment against an annotation holding miRNA transcripts; random-primed total RNA, so not mature arms. Transcripts absent from Ensembl 115 are dropped (median {np.median(unmapped):.3f} of counts per well). Poly(A) wells left out. Vehicle and treated are separate groups."})
