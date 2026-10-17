"""Download the FASTQ files of the v2 test comparisons and prepare split-mate runs.

Reads pipelines/experiments/configs/v2/fastq_manifest.tsv: every file is downloaded from ENA to its repository path
and checked against the manifest md5 and size. Rows marked prepare=split_mate_blocks are then split into _1/_2
files in prepared_dir (see split_mate_blocks.py). Existing files that already match their md5 are kept.

Usage:
  uv run python pipelines/experiments/download_v2_reads.py [--study GSE115646 ...] [--workers 4] [--dry-run]
"""
import argparse
import hashlib
import pathlib
from concurrent.futures import ThreadPoolExecutor

import pandas as pd
import requests

from split_mate_blocks import split

REPO = pathlib.Path(__file__).resolve().parents[2]
MANIFEST = REPO / "pipelines/experiments/configs/v2/fastq_manifest.tsv"


def md5(path: pathlib.Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 22), b""):
            h.update(block)
    return h.hexdigest()


def fetch(row) -> str:
    dest = REPO / row.path
    if dest.exists() and dest.stat().st_size == row.expected_bytes and md5(dest) == row.expected_md5:
        return f"ok      {row.path}"
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_name(dest.name + ".part")
    with requests.get(row.source_url, stream=True, timeout=120) as r:
        r.raise_for_status()
        with open(part, "wb") as fh:
            for chunk in r.iter_content(1 << 22):
                fh.write(chunk)
    if part.stat().st_size != row.expected_bytes or md5(part) != row.expected_md5:
        raise ValueError(f"size or md5 mismatch for {row.path}")
    part.replace(dest)
    return f"fetched {row.path}"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--study", nargs="*", help="only these GEO series")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    man = pd.read_csv(MANIFEST, sep="\t", keep_default_na=False)
    if args.study:
        man = man[man.study.isin(args.study)]
    print(f"{len(man)} files, {man.expected_bytes.sum() / 1e9:.1f} GB")
    if args.dry_run:
        return
    with ThreadPoolExecutor(args.workers) as pool:
        for line in pool.map(fetch, man.itertuples(index=False)):
            print(line, flush=True)
    for row in man[man.prepare == "split_mate_blocks"].itertuples(index=False):
        out_dir = REPO / row.prepared_dir
        if (out_dir / f"{row.run_accession}.audit.json").exists():
            continue
        report = split(REPO / row.path, row.run_accession, out_dir)
        print(f"split   {row.run_accession}: {report['outputs'][0]['records']} pairs", flush=True)


if __name__ == "__main__":
    main()
