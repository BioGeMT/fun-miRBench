"""Human preimplantation embryos, dexamethasone against none (E-MTAB-10098 small RNA, Small-seq; E-MTAB-10096 mRNA, Smart-seq2). Small RNA and mRNA come from
different cells of the same embryos (Genome Res 2022 methods; cell names do not match in the SDRFs), so this is population-level: no same-cell pairing and no miRNA perturbation.
No comparison. Writes the miRNA cell store from our read counts (the deposited mature.counts.tsv.gz is fractional), the mRNA cell store from the deposited
integer table (rpkmforgenes read counts, RefSeq symbols, hg38; symbols mapped through Ensembl 115), and embryo-level sums of both for embryos with both assays.
Usage: EMTAB10098_embryo.py prepare | build <final dir>"""
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import pandas as pd
from scipy import sparse
from funmirbench.single_cell_store import OBS_COLUMNS, write_cell_store
from tier2_plate import DATA, GTF

ROOT = DATA / "E-MTAB-10098"
ALIGNED = ROOT / "aligned"


def sdrf(name: str) -> pd.DataFrame:
    table = pd.read_csv(ROOT / "deposit" / f"{name}.sdrf.txt", sep="\t", dtype=str)
    return pd.DataFrame({"cell": table["Source Name"], "embryo": table["Characteristics[individual]"], "compound": table["Factor Value[compound]"],
                         "run": table["Comment[ENA_RUN]"], "well_quality": table.get("Characteristics[single cell well quality]", "OK")}).drop_duplicates("cell").set_index("cell")


def obs(table: pd.DataFrame, dataset: str) -> pd.DataFrame:
    no_compound = (table.compound == "none").to_numpy()
    return pd.DataFrame({"dataset": dataset, "sample": table.index, "batch": table.embryo, "replicate_unit": table.embryo, "replicate_type": "biological",
                         "arm": table.compound, "control_type": ["no_compound" if value else "" for value in no_compound], "guide": "",
                         "is_control": no_compound, "n_guides": 0, "passed_qc": (table.well_quality == "OK").to_numpy()}, index=table.index)[list(OBS_COLUMNS)]


def prepare() -> None:
    ALIGNED.mkdir(exist_ok=True)
    table = sdrf("E-MTAB-10098")
    (ALIGNED / "small.manifest").write_text("".join(f"{ROOT}/fastq/{row.run}.fastq.gz\t{cell}\n" for cell, row in table.iterrows()))


def build(final: pathlib.Path) -> None:
    (final / "store").mkdir(parents=True)
    (final / "tables").mkdir()
    small = sdrf("E-MTAB-10098")
    mirna = pd.read_csv(ALIGNED / "small/mirna_counts.tsv.gz", sep="\t", index_col=0).reindex(small.index)
    if mirna.isna().any(axis=None):
        raise ValueError("small RNA cells without counts")
    write_cell_store(final / "store/EMTAB10098_embryo_mirna.h5ad", sparse.csr_matrix(mirna.to_numpy("int64")), obs(small, "EMTAB10098_embryo"), pd.Series(mirna.columns), gtf=GTF,
                     gene_ids=pd.Series("", index=range(mirna.shape[1])), symbols=pd.Series(mirna.columns),
                     source="our count of E-MTAB-10098 reads: UMI 8 nt and CA removed, TGGAATTCTCGGGTGCCAAGG trimmed, inserts of 18-30 nt, one molecule per UMI and insert, bowtie 1 to miRBase 22.1 hairpins, reads on one arm only; cell QC: the SDRF's well quality")
    genes = pd.read_csv(ROOT / "processed/merge.rpkmforgenes_counts.tsv.gz", sep="\t", index_col=0)
    mrna_cells = sdrf("E-MTAB-10096").loc[genes.columns]
    write_cell_store(final / "store/EMTAB10098_embryo.h5ad", sparse.csr_matrix(genes.T.to_numpy("int64")), obs(mrna_cells, "EMTAB10098_embryo"), pd.Series(genes.index), gtf=GTF,
                     source="authors' integer read counts, E-MTAB-10096 merge.rpkmforgenes_counts.tsv.gz (hg38, RefSeq); cell QC: none beyond the deposit, all listed cells kept")
    shared = sorted(set(small.embryo) & set(mrna_cells.embryo))
    keep = small.embryo.isin(shared) & (small.well_quality == "OK")
    mirna.loc[keep].T.T.groupby(small.embryo[keep]).sum().T.to_csv(final / "tables/embryo_mirna_counts.tsv.gz", sep="\t")
    genes.T.loc[mrna_cells.embryo.isin(shared)].groupby(mrna_cells.embryo[mrna_cells.embryo.isin(shared)]).sum().T.to_csv(final / "tables/embryo_mrna_counts.tsv.gz", sep="\t")
    pd.DataFrame({"embryo": shared, "compound": [small.compound[small.embryo == embryo].iloc[0] for embryo in shared],
                  "small_rna_cells": [int((keep & (small.embryo == embryo)).sum()) for embryo in shared],
                  "mrna_cells": [int((mrna_cells.embryo == embryo).sum()) for embryo in shared]}).to_csv(final / "tables/embryos.tsv", sep="\t", index=False)
    print(f"{len(shared)} embryos with both assays; miRNA store {small.shape[0]} cells, mRNA store {mrna_cells.shape[0]} cells", flush=True)


if __name__ == "__main__":
    prepare() if sys.argv[1] == "prepare" else build(pathlib.Path(sys.argv[2]))
