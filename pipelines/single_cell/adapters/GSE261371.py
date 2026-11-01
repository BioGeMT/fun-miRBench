"""scComplete-seq: three hashed cell lines and hashed PBMCs. miRNA from the short fraction, mRNA from the long fraction."""
import pathlib, sys
sys.path.insert(0, str(pathlib.Path(__file__).parent))
import numpy as np, pandas as pd, scipy.io
from scipy import sparse
from cooccurrence_common import DATA, write_study

MIN_TOP_TAG_SHARE = 0.6
LINE_MARKERS = {"K562": ["HBG1", "HBG2", "GATA1", "HBZ"], "MCF7": ["ESR1", "TFF1", "KRT8", "KRT18"]}
PBMC_MARKERS = {"T": ["CD3D", "CD3E", "CD2"], "Monocyte": ["LYZ", "CD14", "S100A8"], "B": ["MS4A1", "CD79A"], "NK": ["NKG7", "GNLY", "KLRD1"]}


def read(stem: str):
    root = DATA / "GSE261371/extracted"
    features = pd.read_csv(root / f"{stem}_features.tsv.gz", sep="\t", header=None, names=["id", "symbol", "kind"])
    barcodes = pd.read_csv(root / f"{stem}_barcodes.tsv.gz", sep="\t", header=None)[0].astype(str)
    return scipy.io.mmread(root / f"{stem}_matrix.mtx.gz").T.tocsr(), features, barcodes


def demultiplex(stem: str) -> pd.Series:
    counts, features, barcodes = read(stem)
    tags = (~features.id.astype(str).str.startswith("unmapped")).to_numpy()
    dense = np.asarray(counts[:, np.flatnonzero(tags)].todense())
    total = dense.sum(axis=1)
    top = dense.argmax(axis=1)
    share = np.divide(dense.max(axis=1), total, out=np.zeros_like(total, dtype=float), where=total > 0)
    label = pd.Series(features.id[tags].astype(str).str.split("-").str[0].to_numpy()[top], index=barcodes.to_numpy())
    return label[share >= MIN_TOP_TAG_SHARE]


ARMS = {"lines": ("GSM8141547_CLtotalRNAshort", "GSM8141548_CLtotalRNAlong", "GSM8141549_CLhto"), "PBMC": ("GSM8141550_PBMCtotalRNAshort", "GSM8141551_PBMCtotalRNAlong", "GSM8141552_PBMChto")}
parts = []
for name, (short_stem, long_stem, hto_stem) in [(sys.argv[1], ARMS[sys.argv[1]])]:
    short, short_features, short_barcodes = read(short_stem)
    long, long_features, long_barcodes = read(long_stem)
    hash_tag = demultiplex(hto_stem)
    shared = pd.Index(short_barcodes).intersection(long_barcodes).intersection(hash_tag.index)
    print(name, "short", len(short_barcodes), "long", len(long_barcodes), "hashed", len(hash_tag), "shared", len(shared), flush=True)
    short = short[pd.Index(short_barcodes).get_indexer(shared)]
    long = long[pd.Index(long_barcodes).get_indexer(shared)]
    hash_tag = hash_tag.reindex(shared)
    total = np.asarray(long.sum(axis=1)).ravel()
    normalised = sparse.diags(1e4 / np.maximum(total, 1)) @ long

    def score(markers):
        columns = np.flatnonzero(long_features.symbol.isin(markers).to_numpy())
        return np.asarray(normalised[:, columns].mean(axis=1)).ravel()

    if name == "lines":
        scores = pd.DataFrame({line: score(markers) for line, markers in LINE_MARKERS.items()}).groupby(hash_tag.to_numpy()).mean()
        print(scores.round(2).to_string(), flush=True)
        line_of_tag = {scores[line].idxmax(): line for line in LINE_MARKERS}
        line_of_tag.update({tag: "HEK293T" for tag in scores.index if tag not in line_of_tag})
        print(line_of_tag, flush=True)
        groups = hash_tag.map(line_of_tag).to_numpy()
        annotation = groups
    else:
        scores = pd.DataFrame({kind: score(markers) for kind, markers in PBMC_MARKERS.items()})
        annotation = np.where(scores.max(axis=1) > 0, scores.idxmax(axis=1), "unassigned")
        groups = hash_tag.to_numpy() + "_" + annotation
    parts.append((name, shared, groups, annotation, hash_tag.to_numpy(), short, short_features, long, long_features))

name, shared, groups, annotation, hash_tag, short, short_features, long, long_features = parts[0]
mirna_rows = short_features.id.astype(str).str.startswith("MI0").to_numpy()
gene_rows = long_features.id.astype(str).str.startswith("ENSG").to_numpy()
cells = pd.DataFrame({"cell": shared, "group": groups, "hash_tag": hash_tag, "annotation": annotation})
write_study(f"GSE261371_{name}", mirna=short[:, np.flatnonzero(mirna_rows)], mirna_features=short_features.symbol[mirna_rows], mrna=long[:, np.flatnonzero(gene_rows)], mrna_features=long_features.id[gene_rows], cells=cells, library_totals=pd.Series(np.asarray(short.sum(axis=1)).ravel(), index=shared), library_total_source="sum of every feature of the short-fraction library, which holds the miRNA features", info={
    "feature_level": "precursor", "same_cell": "1x SPRI (short) and 0.6x SPRI (long) fractions of one droplet's cDNA, same barcodes (GEO)", "values": "deposited counts",
    "grouping": "cell line from hash tag" if name == "lines" else "hash tag and coarse marker annotation", "qc_source": "deposit's barcodes; hash tag kept when the top tag holds at least 60% of tag counts (my rule, set before any correlation)", "source": "authors, GSE261371 RAW",
    "caveats": "Features are miRBase precursors (MI accessions, hsa-mir-*) in the authors' custom annotation. " + ("The deposit does not map tags to cell lines: they are named from K562 and MCF7 marker genes. The scmRNA-seq arm has no small RNA and is left out." if name == "lines" else "The deposit does not map tags to stimuli, so groups keep the tag name (unstimulated, LPS or PMA-ionomycin, not resolved). Cell types are my coarse marker annotation, not the authors'.")})
