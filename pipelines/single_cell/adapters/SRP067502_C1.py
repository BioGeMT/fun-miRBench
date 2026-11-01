"""Doxycycline-inducible miR-199a with GFP, HEK293 lines i199 and i199-KTN1 (Rzepiela et al., preprint 10.1101/035097), 597 Fluidigm C1 cells on seven chips, doses pooled on each chip.
Table S1 lists the doses of chips a to g; the deposit names chips by date, and neither source links the two, so no dose is stored per cell or chip.
Arms come from GFP as in the preprint's methods: control log2 GFP < 8, induced log2 GFP > 12, in TPM. Our TPM: reads per kb (union of exons, eGFP CDS 0.72 kb), GFP reads included in the total.
Cell QC from the same methods, adapted: the preprint removes cells whose median Pearson correlation of log2 expression with the other cells is below 0.64, 5-10% of cells. On our counts (log2(TPM + 1), genes of mean TPM >= 15, per line) the median cell sits at 0.54, so that value would remove every cell; we keep the stated fraction and remove the lowest QC_FRACTION of each line. The preprint's second rule (total gene expression below about 500,000 TPM) has no equivalent in our gene counts and is not applied.
GFP reads: at least 30 nt (the deposit's reads are adapter-trimmed), at most 2 mismatches.
Chips are replicate units within each line. Usage: SRP067502_C1.py prepare | build <store path>"""
import gzip
import pathlib
import re
import sys

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import numpy as np
import pandas as pd
from funmirbench.single_cell_store import read_10x_mtx
from tier2_plate import DATA, GTF, write_plate_store

ROOT = DATA / "SRP067502"
ALIGNED = ROOT / "aligned"
EGFP_KB = 0.72
LOW_LOG2, HIGH_LOG2, QC_FRACTION, MIN_MEAN_TPM = 8, 12, 0.075, 15


def report() -> pd.DataFrame:
    table = pd.read_csv(ROOT / "meta/SRP067502.tsv", sep="\t")
    table["chip"] = table.library_name.str.split("/").str[0]
    table["cell"] = table.library_name.str.replace("/", "__", regex=False)
    table["line"] = np.where(table.chip.str.contains("ktn1"), "i199-KTN1", "i199")
    return table.set_index("cell")


def prepare() -> None:
    ALIGNED.mkdir(exist_ok=True)
    table = report()
    (ALIGNED / "cells.manifest").write_text("".join(f"{ROOT}/fastq/{r.run_accession}.fastq.gz\t-\t{cell}\n" for cell, r in table.iterrows()))
    (ALIGNED / "markers.manifest").write_text("".join(f"{ROOT}/fastq/{r.run_accession}.fastq.gz\t{cell}\n" for cell, r in table.iterrows()))


def gene_kb() -> pd.Series:
    exons: dict[str, list[tuple[str, int, int]]] = {}
    with gzip.open(GTF, "rt") as handle:
        for line in handle:
            fields = line.split("\t")
            if len(fields) > 8 and fields[2] == "exon":
                exons.setdefault(re.search(r'gene_id "([^"]+)"', fields[8]).group(1), []).append((int(fields[3]), int(fields[4])))
    lengths = {}
    for gene, spans in exons.items():
        total, end = 0, 0
        for start, stop in sorted(spans):
            total += max(0, stop - max(start, end + 1) + 1)
            end = max(end, stop)
        lengths[gene] = total / 1000
    return pd.Series(lengths)


def build(out: pathlib.Path) -> None:
    table = report()
    solo = ALIGNED / "cells/Solo.out/Gene/raw"
    counts, features, barcodes = read_10x_mtx(solo / "matrix.mtx", solo / "features.tsv", solo / "barcodes.tsv")
    gfp = pd.read_csv(ALIGNED / "marker_min30/counts.tsv", sep="\t").set_index("cell").marker_reads.reindex(barcodes).to_numpy(float)
    kb = gene_kb().reindex(features[0]).to_numpy()
    if np.isnan(kb).any():
        raise ValueError("genes without exons in the GTF")
    matrix = counts.toarray().astype(float)  # genes x cells
    per_kb = matrix / kb[:, None]
    scale = per_kb.sum(axis=0) + gfp / EGFP_KB
    tpm = per_kb / scale * 1e6
    gfp_tpm = pd.Series(gfp / EGFP_KB / scale * 1e6, index=barcodes)
    median_r = pd.Series(np.nan, index=barcodes)
    for line in ("i199", "i199-KTN1"):
        cells = np.flatnonzero(table.line.reindex(barcodes).to_numpy() == line)
        keep = tpm[:, cells].mean(axis=1) >= MIN_MEAN_TPM
        r = np.corrcoef(np.log2(tpm[np.ix_(keep, cells)] + 1).T)
        np.fill_diagonal(r, np.nan)
        median_r.iloc[cells] = np.nanmedian(r, axis=0)
    log2_gfp = np.log2(gfp_tpm + 1)
    arm_of = np.where(log2_gfp < LOW_LOG2, "GFP_low", np.where(log2_gfp > HIGH_LOG2, "GFP_high", "GFP_middle"))
    info = table.reindex(barcodes)
    cells = pd.DataFrame({"sample": info.chip, "batch": info.chip, "replicate_unit": info.chip, "replicate_type": "biological",
                          "arm": [f"{arm}_{line}" for arm, line in zip(arm_of, info.line)],
                          "control_type": np.where(arm_of == "GFP_low", "low_GFP_same_chip", ""),
                          "passed_qc": median_r.to_numpy() >= median_r.groupby(info.line.to_numpy()).transform(lambda values: values.quantile(QC_FRACTION)).to_numpy()}, index=barcodes)
    print(pd.DataFrame({"chip": info.chip, "arm": arm_of, "passed": cells.passed_qc}).groupby(["chip", "arm"]).passed.agg(["size", "sum"]).to_string(), flush=True)
    print(f"fail the correlation rule: {int((~cells.passed_qc).sum())} of {len(cells)}; median correlation at the cut per line: {median_r[~cells.passed_qc].groupby(info.line[~cells.passed_qc.to_numpy()].to_numpy()).max().round(3).to_dict()}", flush=True)
    pd.DataFrame({"chip": info.chip, "line": info.line, "gfp_reads": gfp, "gfp_tpm": gfp_tpm, "median_r": median_r}).rename_axis("cell").to_csv(ALIGNED / "cell_table.tsv", sep="\t")
    write_plate_store("SRP067502_C1", solo, cells, out, "our STARsolo SmartSeq alignment against Ensembl 115, unique mappers, reads SRA SRP067502; cell QC: lowest 7.5% of the median Pearson correlation of log2 expression with the other cells removed per line (preprint: 5-10% removed, by a cut-off of 0.64 on its own values)")


if __name__ == "__main__":
    prepare() if sys.argv[1] == "prepare" else build(pathlib.Path(sys.argv[2]))
