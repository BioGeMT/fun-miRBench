#!/usr/bin/env python3
"""Build the MANE Select gene-level 3'UTR conservation table used by Figure 2F."""

from __future__ import annotations

import argparse
import gzip
import logging
import re
from pathlib import Path

import numpy as np
import pandas as pd
import pyBigWig

from funmirbench.logger import setup_logging


DEFAULT_GTF = Path("data/common_resources/ensembl/Homo_sapiens.GRCh38.115.gtf.gz")
DEFAULT_BIGWIG = Path("data/resources/conservation/hg38.phastCons100way.bw")
DEFAULT_OUT = Path("manuscript_assets/tables/figure2F_utr_conservation_mane_raw.tsv")
ATTR_RE = re.compile(r'([A-Za-z0-9_]+) "([^"]*)"')
UTR3_FEATURE_NAMES = {"three_prime_utr", "three_prime_UTR", "3UTR", "3utr"}
logger = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute mean phastCons100way score for each protein-coding MANE Select 3'UTR."
    )
    parser.add_argument("--gtf", type=Path, default=DEFAULT_GTF)
    parser.add_argument("--bigwig", type=Path, default=DEFAULT_BIGWIG)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--log-level", default="INFO")
    return parser.parse_args()


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def parse_attributes(value: str) -> dict[str, str]:
    return {match.group(1): match.group(2) for match in ATTR_RE.finditer(value)}


def strip_version(gene_id: str) -> str:
    return str(gene_id).split(".", 1)[0]


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not intervals:
        return []
    ordered = sorted(intervals)
    merged = [ordered[0]]
    for start, end in ordered[1:]:
        previous_start, previous_end = merged[-1]
        if start <= previous_end + 1:
            merged[-1] = (previous_start, max(previous_end, end))
        else:
            merged.append((start, end))
    return merged


def interval_length(intervals: list[tuple[int, int]]) -> int:
    return sum(end - start + 1 for start, end in intervals)


def normalize_chrom_for_bigwig(chrom: str, chrom_sizes: dict[str, int]) -> str | None:
    chrom = str(chrom).strip()
    candidates = []
    if chrom == "MT":
        candidates.extend(["chrM", "MT", "M"])
    else:
        candidates.append(chrom)
        candidates.append(chrom.removeprefix("chr") if chrom.startswith("chr") else f"chr{chrom}")
    for candidate in candidates:
        if candidate in chrom_sizes:
            return candidate
    return None


def parse_mane_select_utr3_intervals(gtf_path: Path) -> dict[str, dict[str, object]]:
    """Return annotated 3'UTRs for each protein-coding MANE Select transcript."""
    if not gtf_path.exists():
        raise FileNotFoundError(f"GTF does not exist: {gtf_path}")

    intervals_by_gene_tx: dict[str, dict[str, dict[str, list[tuple[int, int]]]]] = {}
    with open_text(gtf_path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] not in UTR3_FEATURE_NAMES:
                continue
            if 'tag "MANE_Select"' not in fields[8]:
                continue
            attrs = parse_attributes(fields[8])
            if attrs.get("transcript_biotype") != "protein_coding":
                continue
            gene_id = attrs.get("gene_id")
            transcript_id = attrs.get("transcript_id")
            if not gene_id or not transcript_id:
                continue
            key_gene = strip_version(gene_id)
            key_transcript = strip_version(transcript_id)
            intervals_by_gene_tx.setdefault(key_gene, {}).setdefault(key_transcript, {}).setdefault(fields[0], []).append(
                (int(fields[3]), int(fields[4]))
            )

    selected = {}
    for gene_id, transcript_intervals in intervals_by_gene_tx.items():
        for transcript_id, chrom_intervals in transcript_intervals.items():
            merged_by_chrom = {
                chrom: merge_intervals(intervals)
                for chrom, intervals in chrom_intervals.items()
                if intervals
            }
            length = sum(interval_length(intervals) for intervals in merged_by_chrom.values())
            if length > 0:
                selected[gene_id] = {
                    "transcript_id": transcript_id,
                    "utr3_length": int(length),
                    "intervals_by_chrom": merged_by_chrom,
                }
    if not selected:
        raise ValueError(f"No protein-coding MANE Select 3'UTR intervals parsed from {gtf_path}")
    return selected


def parse_longest_protein_coding_utr3_intervals(gtf_path: Path) -> dict[str, dict[str, object]]:
    """Return each gene's longest protein-coding 3'UTR (legacy helper)."""
    if not gtf_path.exists():
        raise FileNotFoundError(f"GTF does not exist: {gtf_path}")
    intervals_by_gene_tx: dict[str, dict[str, dict[str, list[tuple[int, int]]]]] = {}
    with open_text(gtf_path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] not in UTR3_FEATURE_NAMES:
                continue
            attrs = parse_attributes(fields[8])
            biotype = (
                attrs.get("transcript_biotype")
                or attrs.get("transcript_type")
                or attrs.get("gene_biotype")
                or attrs.get("gene_type")
            )
            if biotype and biotype != "protein_coding":
                continue
            gene_id = attrs.get("gene_id")
            transcript_id = attrs.get("transcript_id")
            if not gene_id or not transcript_id:
                continue
            intervals_by_gene_tx.setdefault(strip_version(gene_id), {}).setdefault(
                strip_version(transcript_id), {}
            ).setdefault(fields[0], []).append((int(fields[3]), int(fields[4])))
    selected = {}
    for gene_id, transcript_intervals in intervals_by_gene_tx.items():
        candidates = []
        for transcript_id, chrom_intervals in transcript_intervals.items():
            merged_by_chrom = {
                chrom: merge_intervals(intervals)
                for chrom, intervals in chrom_intervals.items()
                if intervals
            }
            length = sum(interval_length(intervals) for intervals in merged_by_chrom.values())
            if length > 0:
                candidates.append((transcript_id, length, merged_by_chrom))
        if candidates:
            transcript_id, length, merged_by_chrom = sorted(
                candidates, key=lambda item: (-item[1], item[0])
            )[0]
            selected[gene_id] = {
                "transcript_id": transcript_id,
                "utr3_length": int(length),
                "intervals_by_chrom": merged_by_chrom,
            }
    return selected


def mean_bigwig_signal(bw, chrom_intervals: dict[str, list[tuple[int, int]]]) -> tuple[float, int]:
    chrom_sizes = bw.chroms()
    chunks = []
    scored_bases = 0
    for chrom, intervals in chrom_intervals.items():
        bw_chrom = normalize_chrom_for_bigwig(chrom, chrom_sizes)
        if bw_chrom is None:
            continue
        chrom_size = chrom_sizes[bw_chrom]
        for start_1based, end_1based in intervals:
            start_0based = max(0, int(start_1based) - 1)
            end_0based = min(int(end_1based), chrom_size)
            if start_0based >= end_0based:
                continue
            values = bw.values(bw_chrom, start_0based, end_0based)
            arr = np.asarray([np.nan if value is None else value for value in values], dtype=float)
            valid = arr[np.isfinite(arr)]
            if valid.size:
                chunks.append(valid)
                scored_bases += int(valid.size)
    if not chunks:
        return np.nan, 0
    return float(np.mean(np.concatenate(chunks))), scored_bases


def build_table(gtf_path: Path, bigwig_path: Path) -> pd.DataFrame:
    if not bigwig_path.exists():
        raise FileNotFoundError(
            f"phastCons BigWig does not exist: {bigwig_path}. "
            "Download it with: uv run python scripts/download_phastcons100way.py"
        )

    selected_by_gene = parse_mane_select_utr3_intervals(gtf_path)
    logger.info("Parsed MANE Select 3'UTRs for %d genes", len(selected_by_gene))
    rows = []
    with pyBigWig.open(str(bigwig_path)) as bw:
        for index, (gene_id, selected) in enumerate(selected_by_gene.items(), start=1):
            mean_score, scored_bases = mean_bigwig_signal(bw, selected["intervals_by_chrom"])
            rows.append(
                {
                    "gene_id": gene_id,
                    "transcript_id": selected["transcript_id"],
                    "utr3_mean_conservation": mean_score,
                    "utr3_scored_bases": scored_bases,
                    "utr3_length": selected["utr3_length"],
                }
            )
            if index % 1000 == 0:
                logger.info("Scored %d/%d genes", index, len(selected_by_gene))
    table = pd.DataFrame(rows).sort_values("gene_id").reset_index(drop=True)
    return table


def main() -> int:
    args = parse_args()
    setup_logging(args.log_level)
    table = build_table(args.gtf, args.bigwig)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.out, sep="\t", index=False)
    logger.info("Wrote %d gene-level conservation rows: %s", len(table), args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
