"""Shared part of the co-occurrence adapters: split a deposit into the miRNA table, the mRNA table and cells.tsv, write
them under the study's final folder, and run funmirbench-single-cell-cooccurrence on them."""

from __future__ import annotations

import gzip
import json
import os
import pathlib
import re
import subprocess
import sys

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from funmirbench.single_cell_store import map_features_to_gene_ids

REPO = pathlib.Path(__file__).resolve().parents[3]
FINAL = pathlib.Path("/SCRATCH/dtzim01/fun-experiment-data/single-cell-final/20261031/cooccurrence_library_total")
DATA = pathlib.Path("/SCRATCH/dtzim01/single_cell_mirna_data")
GTF = pathlib.Path("/SCRATCH/dtzim01/fun-experiment-data/references/ensembl_v115/Homo_sapiens.GRCh38.115.gtf.gz")
SCRATCH = pathlib.Path("/SCRATCH/dtzim01/tmp/cooccurrence")
BENCHMARK_GENES = REPO / "metadata/benchmark_genes.tsv.gz"
# miRNA gene symbols, without host genes (MIR100HG), antisense (MIR4435-2-AS1) and other transcripts named after them.
MIRNA_SYMBOL = re.compile(r"^(MIR\d+[A-Z]*\d*(-\d+)?|MIRLET7[A-Z]*\d*(-\d+)?|hsa-mir-[0-9a-z-]+|hsa-let-[0-9a-z-]+)$")
NOT_MIRNA = re.compile(r"(HG|-AS\d*|-DT|-OT\d*|-IT\d*)$")


def gtf_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Genes (gene_id, symbol, biotype) and transcripts (transcript_id, gene_id) of Ensembl 115, cached in scratch."""
    genes_path, transcripts_path = SCRATCH / "ensembl115_genes.tsv", SCRATCH / "ensembl115_transcripts.tsv"
    if not genes_path.exists() or not transcripts_path.exists():
        genes, transcripts = [], []
        with gzip.open(GTF, "rt") as handle:
            for line in handle:
                fields = line.split("\t")
                if len(fields) < 9:
                    continue
                gene = re.search(r'gene_id "([^"]+)"', fields[8])
                if fields[2] == "gene":
                    name = re.search(r'gene_name "([^"]+)"', fields[8])
                    genes.append((gene.group(1), name.group(1) if name else "", re.search(r'gene_biotype "([^"]+)"', fields[8]).group(1)))
                elif fields[2] == "transcript":
                    transcripts.append((re.search(r'transcript_id "([^"]+)"', fields[8]).group(1), gene.group(1)))
        pd.DataFrame(genes, columns=["gene_id", "symbol", "biotype"]).to_csv(genes_path, sep="\t", index=False)
        pd.DataFrame(transcripts, columns=["transcript_id", "gene_id"]).to_csv(transcripts_path, sep="\t", index=False)
    return (pd.read_csv(genes_path, sep="\t", keep_default_na=False), pd.read_csv(transcripts_path, sep="\t"))


def benchmark_gene_ids() -> set[str]:
    return set(pd.read_csv(BENCHMARK_GENES, sep="\t", usecols=["gene_id"]).gene_id)


def is_mirna_symbol(names: pd.Series) -> np.ndarray:
    names = pd.Series(names, dtype=str)
    return (names.str.match(MIRNA_SYMBOL) & ~names.str.contains(NOT_MIRNA)).to_numpy()


def mirna_ensembl_ids() -> set[str]:
    genes, _ = gtf_tables()
    return set(genes.gene_id[genes.biotype == "miRNA"])


def write_study(dataset: str, *, mirna: sparse.spmatrix, mirna_features, mrna: sparse.spmatrix, mrna_features,
                cells: pd.DataFrame, library_totals: pd.Series | None, library_total_source: str, mirna_cells=None, mrna_cells=None,
                info: dict) -> None:
    """Write store/<dataset>_mirna.h5ad, store/<dataset>_mrna.h5ad, store/cells.tsv and study.json, then run the step
    with --out correlation. Cells are rows; the mRNA table keeps the features in the benchmark's 19,226-gene list.
    library_totals (index = cell): the sum of every deposited feature of the library holding the miRNAs, or None when
    the deposit has no such total; library_total_source says which, and why."""
    out = FINAL / dataset
    if out.exists():
        raise FileExistsError(f"{out} exists; a finished dataset is never overwritten")
    store = out / "store"
    store.mkdir(parents=True)
    mirna, mrna = sparse.csr_matrix(mirna), sparse.csr_matrix(mrna)
    mirna_features = pd.Index(pd.Series(mirna_features, dtype=str))
    mrna_features = pd.Series(mrna_features, dtype=str).reset_index(drop=True)
    if mirna_features.has_duplicates:
        raise ValueError(f"{dataset}: duplicated miRNA features {mirna_features[mirna_features.duplicated()][:5].tolist()}")
    for name, table in (("miRNA", mirna), ("mRNA", mrna)):
        if table.nnz and not (table.data >= 0).all():  # also catches NaN
            raise ValueError(f"{dataset}: negative or missing {name} values")
    gene_ids = map_features_to_gene_ids(mrna_features, GTF)
    keep = gene_ids.isin(benchmark_gene_ids()).to_numpy()
    if not keep.any():
        raise ValueError(f"{dataset}: no mRNA feature maps to the benchmark gene list")
    mirna_cells = pd.Index(cells.cell if mirna_cells is None else mirna_cells, dtype=str)
    mrna_cells = pd.Index(cells.cell if mrna_cells is None else mrna_cells, dtype=str)
    integer = {name: bool(np.array_equal(table.data, np.round(table.data))) for name, table in (("mirna", mirna), ("mrna", mrna))}
    miRNA_table = ad.AnnData(X=mirna.astype("int32" if integer["mirna"] else "float32"), obs=pd.DataFrame(index=mirna_cells),
                             var=pd.DataFrame({"feature_level": info["feature_level"]}, index=mirna_features))
    mRNA_table = ad.AnnData(X=mrna[:, np.flatnonzero(keep)].astype("int32" if integer["mrna"] else "float32"),
                            obs=pd.DataFrame(index=mrna_cells),
                            var=pd.DataFrame({"gene_id": gene_ids[keep].to_numpy()}, index=pd.Index(mrna_features[keep].to_numpy())))
    miRNA_table.var_names_make_unique(), mRNA_table.var_names_make_unique()
    miRNA_table.write_h5ad(store / f"{dataset}_mirna.h5ad", compression="gzip")
    mRNA_table.write_h5ad(store / f"{dataset}_mrna.h5ad", compression="gzip")
    cells.to_csv(store / "cells.tsv", sep="\t", index=False)
    allow_fractional = not (integer["mirna"] and integer["mrna"])
    if library_totals is not None:
        library_totals.rename("library_total").rename_axis("cell").reset_index().to_csv(store / "library_totals.tsv", sep="\t", index=False)
    (out / "study.json").write_text(json.dumps(
        {"dataset": dataset, **info, "mirna_values_integer": integer["mirna"], "mrna_values_integer": integer["mrna"],
         "allow_fractional": allow_fractional, "mirna_normalisation": "miRNA_table_total" if library_totals is None else "whole_library_total",
         "library_total_source": library_total_source,
         "cells_in_mirna_table": len(mirna_cells), "cells_in_mrna_table": len(mrna_cells), "cells_listed": len(cells),
         "mirna_features": int(mirna.shape[1]), "mrna_features_deposited": int(mrna.shape[1]),
         "mrna_features_kept": int(keep.sum())}, indent=2) + "\n")
    step = [str(out / "store" / f"{dataset}_mirna.h5ad"), str(out / "store" / f"{dataset}_mrna.h5ad"), str(store / "cells.tsv")]
    arguments = ["--mirna", step[0], "--mrna", step[1], "--cells", step[2], "--group-column", "group", "--out", str(out / "correlation")]
    if allow_fractional:
        arguments.append("--allow-fractional")
    if library_totals is not None:
        arguments += ["--library-totals", str(store / "library_totals.tsv")]
    command = [str(pathlib.Path(sys.executable).parent / "funmirbench-single-cell-cooccurrence"), *arguments]
    subprocess.run(command, check=True, cwd=SCRATCH, env={**os.environ, "OPENBLAS_NUM_THREADS": "4", "OMP_NUM_THREADS": "4"})
