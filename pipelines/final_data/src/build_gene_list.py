"""Benchmark gene list: Ensembl 115 protein-coding genes with a MANE Select transcript.

The list depends only on the annotation, never on a comparison's counts, so every comparison is scored on the
same genes.

Usage:
  uv run python pipelines/final_data/src/build_gene_list.py \
      --gtf data/experiments/raw/refs/ensembl_v115/Homo_sapiens.GRCh38.115.gtf.gz \
      --out metadata/benchmark_genes.tsv.gz
"""
import argparse
import gzip
import pathlib
import re

import pandas as pd

ATTR = re.compile(r'(\w+) "([^"]*)"')


def read_gtf(path: pathlib.Path) -> tuple[pd.DataFrame, dict]:
    opener = gzip.open if path.suffix == ".gz" else open
    genes, mane = {}, {}
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[2] not in ("gene", "transcript"):
                continue
            a = dict(ATTR.findall(f[8]))
            if f[2] == "gene":
                genes[a["gene_id"]] = {"gene_id": a["gene_id"], "gene_version": a.get("gene_version", ""),
                                       "gene_name": a.get("gene_name", ""), "gene_biotype": a["gene_biotype"],
                                       "chromosome": f[0], "start": int(f[3]), "end": int(f[4]), "strand": f[6]}
            elif 'tag "MANE_Select"' in f[8]:
                if a["gene_id"] in mane:
                    raise ValueError(f"two MANE Select transcripts for {a['gene_id']}")
                mane[a["gene_id"]] = f"{a['transcript_id']}.{a.get('transcript_version', '')}".rstrip(".")
    return pd.DataFrame(genes.values()), mane


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--gtf", type=pathlib.Path, required=True)
    ap.add_argument("--out", type=pathlib.Path, required=True)
    args = ap.parse_args()
    genes, mane = read_gtf(args.gtf)
    pc = genes[genes.gene_biotype == "protein_coding"]
    out = pc[pc.gene_id.isin(mane)].copy()
    out["mane_select_transcript"] = out.gene_id.map(mane)
    out = out.drop(columns="gene_biotype").sort_values(["chromosome", "start", "gene_id"])
    args.out.parent.mkdir(parents=True, exist_ok=True)
    gz = {"method": "gzip", "mtime": 0} if args.out.suffix == ".gz" else None
    out.to_csv(args.out, sep="\t", index=False, compression=gz)
    print(f"{len(pc)} protein-coding genes; {len(out)} with a MANE Select transcript -> {args.out}")


if __name__ == "__main__":
    main()
