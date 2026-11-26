"""Gene-level annotation features for predictor coverage groups."""
from __future__ import annotations

import gzip
import logging
import pathlib
import re
from dataclasses import dataclass

import numpy as np
import pandas as pd
import requests
from scipy.special import log_ndtr
from scipy.stats import mannwhitneyu

from funmirbench.gene_ids import strip_ensembl_version

DEFAULT_GTF = pathlib.Path("data/experiments/raw/refs/ensembl_v115/Homo_sapiens.GRCh38.115.gtf.gz")
DEFAULT_PHASTCONS = pathlib.Path("data/common_resources/conservation/hg38.phastCons100way.bw")
GTF_URL = "https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz"
PHASTCONS_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons100way/hg38.phastCons100way.bw"
UTR_FEATURES = {"three_prime_utr", "three_prime_UTR", "3UTR", "3utr"}
ATTRIBUTE_PATTERN = re.compile(r'([A-Za-z0-9_]+) "([^"]*)"')
FEATURE_COLUMNS = ["gene_id", "transcript_id", "utr3_length_nt", "utr3_mean_phastcons", "utr3_scored_bases"]
logger = logging.getLogger(__name__)


def _format_bytes(value: int) -> str:
    size = float(value)
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if size < 1024 or unit == "TiB":
            return f"{size:.1f} {unit}"
        size /= 1024
    return f"{size:.1f} TiB"


def format_pvalue(pvalue, log10_pvalue=None) -> str:
    """Format a p-value, retaining a log-domain value when the float underflows to zero."""
    if not pd.isna(pvalue) and pvalue > 0:
        return f"{pvalue:.2g}"
    if pd.isna(log10_pvalue):
        return "NA"
    exponent = int(np.floor(log10_pvalue))
    mantissa = 10 ** (log10_pvalue - exponent)
    return f"{mantissa:.1f}e{exponent}"


def _mann_whitney_log10_pvalue(x: np.ndarray, y: np.ndarray, statistic: float, pvalue: float) -> float:
    """Return the log10 asymptotic p-value, including values below double precision."""
    if pvalue > 0:
        return float(np.log10(pvalue))
    n1, n2 = len(x), len(y)
    _, ties = np.unique(np.concatenate([x, y]), return_counts=True)
    n = n1 + n2
    tie_term = np.sum(ties.astype(float) ** 3 - ties)
    standard_deviation = np.sqrt(n1 * n2 / 12 * ((n + 1) - tie_term / (n * (n - 1))))
    u = max(statistic, n1 * n2 - statistic)
    z = (u - n1 * n2 / 2 - 0.5) / standard_deviation
    return float((np.log(2) + log_ndtr(-z)) / np.log(10))


@dataclass(frozen=True)
class GeneSetFeatureSettings:
    enabled: bool
    gtf_path: pathlib.Path
    phastcons_path: pathlib.Path


def parse_settings(value, root: pathlib.Path) -> GeneSetFeatureSettings:
    if value is None:
        value = {}
    if not isinstance(value, dict):
        raise ValueError("reports.gene_set_features must be a mapping")
    enabled = value.get("enabled", False)
    if not isinstance(enabled, bool):
        raise ValueError("reports.gene_set_features.enabled must be a boolean")

    def resolve(path_value, default):
        path = pathlib.Path(str(path_value if path_value is not None else default)).expanduser()
        return path if path.is_absolute() else root / path

    return GeneSetFeatureSettings(
        enabled=enabled,
        gtf_path=resolve(value.get("gtf_path"), DEFAULT_GTF),
        phastcons_path=resolve(value.get("phastcons_path"), DEFAULT_PHASTCONS),
    )


def _download_resource(url: str, path: pathlib.Path, *, label: str, magic: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_name(f"{path.name}.part")
    try:
        logger.info("Downloading %s from %s", label, url)
        with requests.get(url, stream=True, timeout=(20, 120)) as response:
            response.raise_for_status()
            expected_size = int(response.headers.get("content-length", 0) or 0)
            if expected_size:
                logger.info("%s download size: %s", label, _format_bytes(expected_size))
            downloaded = 0
            next_progress = 10
            with partial.open("wb") as handle:
                for chunk in response.iter_content(chunk_size=8 * 1024 * 1024):
                    if not chunk:
                        continue
                    handle.write(chunk)
                    downloaded += len(chunk)
                    if expected_size:
                        progress = int(downloaded * 100 / expected_size)
                        if progress >= next_progress:
                            logger.info("%s download: %d%%", label, min(progress, 100))
                            next_progress = (progress // 10 + 1) * 10
            if expected_size and downloaded != expected_size:
                raise IOError(f"incomplete download ({downloaded} of {expected_size} bytes)")
        with partial.open("rb") as handle:
            if handle.read(len(magic)) != magic:
                raise IOError("downloaded file has an unexpected format")
        partial.replace(path)
        logger.info("Downloaded %s to %s", label, path)
    except Exception as exc:
        partial.unlink(missing_ok=True)
        raise RuntimeError(f"Could not download {label} from {url}: {exc}") from exc


def ensure_feature_resources(settings: GeneSetFeatureSettings) -> None:
    """Download the standard annotation files on demand and keep them in the configured paths."""
    if not settings.enabled:
        return
    try:
        import pyBigWig  # noqa: F401
    except ImportError as exc:
        raise RuntimeError("Install the optional feature dependency with `uv sync --extra gene-features`.") from exc
    if settings.gtf_path.is_file():
        logger.info("Ensembl 115 GTF already present: %s", settings.gtf_path)
    else:
        _download_resource(GTF_URL, settings.gtf_path, label="Ensembl 115 GTF", magic=b"\x1f\x8b")
    if settings.phastcons_path.is_file():
        logger.info("hg38 phastCons100way track already present: %s", settings.phastcons_path)
    else:
        _download_resource(PHASTCONS_URL, settings.phastcons_path,
                           label="hg38 phastCons100way BigWig", magic=b"\x26\xfc\x8f\x88")


def _merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
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


def _read_gtf(path: pathlib.Path, gene_ids: set[str]) -> dict[str, tuple[str, dict[str, list[tuple[int, int]]], int]]:
    opener = gzip.open if path.suffix == ".gz" else open
    by_transcript: dict[tuple[str, str], dict[str, list[tuple[int, int]]]] = {}
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] not in UTR_FEATURES:
                continue
            attrs = dict(ATTRIBUTE_PATTERN.findall(fields[8]))
            gene_id = strip_ensembl_version(attrs.get("gene_id", ""))
            transcript_id = strip_ensembl_version(attrs.get("transcript_id", ""))
            biotype = (attrs.get("transcript_biotype") or attrs.get("transcript_type")
                       or attrs.get("gene_biotype") or attrs.get("gene_type"))
            if (gene_id not in gene_ids or not transcript_id
                    or (biotype is not None and biotype != "protein_coding")):
                continue
            start, end = int(fields[3]), int(fields[4])
            by_transcript.setdefault((gene_id, transcript_id), {}).setdefault(fields[0], []).append((start, end))

    candidates: dict[str, list[tuple[str, dict[str, list[tuple[int, int]]], int]]] = {}
    for (gene_id, transcript_id), chrom_intervals in by_transcript.items():
        merged = {chrom: _merge_intervals(intervals) for chrom, intervals in chrom_intervals.items()}
        length = sum(end - start + 1 for intervals in merged.values() for start, end in intervals)
        if length:
            candidates.setdefault(gene_id, []).append((transcript_id, merged, length))
    return {gene_id: min(rows, key=lambda row: (-row[2], row[0])) for gene_id, rows in candidates.items()}


def _bigwig_chrom(chrom: str, chrom_sizes: dict[str, int]) -> str | None:
    candidates = [chrom]
    if chrom == "MT":
        candidates.extend(["chrM", "M"])
    elif chrom.startswith("chr"):
        candidates.append(chrom.removeprefix("chr"))
    else:
        candidates.append(f"chr{chrom}")
    return next((candidate for candidate in candidates if candidate in chrom_sizes), None)


def _mean_signal(bw, intervals_by_chrom: dict[str, list[tuple[int, int]]]) -> tuple[float, int]:
    chrom_sizes = bw.chroms()
    total, scored_bases = 0.0, 0
    for chrom, intervals in intervals_by_chrom.items():
        bw_chrom = _bigwig_chrom(chrom, chrom_sizes)
        if bw_chrom is None:
            continue
        for start, end in intervals:
            start0 = max(0, start - 1)
            end0 = min(end, chrom_sizes[bw_chrom])
            if start0 >= end0:
                continue
            values = np.asarray([np.nan if value is None else value
                                 for value in bw.values(bw_chrom, start0, end0)], dtype=float)
            valid = values[np.isfinite(values)]
            total += float(valid.sum())
            scored_bases += int(valid.size)
    return (total / scored_bases if scored_bases else np.nan), scored_bases


def build_utr3_feature_table(gtf_path: pathlib.Path, phastcons_path: pathlib.Path,
                             gene_ids: pd.Index) -> pd.DataFrame:
    try:
        import pyBigWig
    except ImportError as exc:
        raise RuntimeError("Install the optional feature dependency with `uv sync --extra gene-features`.") from exc

    requested = {strip_ensembl_version(gene_id) for gene_id in gene_ids}
    selected = _read_gtf(gtf_path, requested)
    bw = pyBigWig.open(str(phastcons_path))
    try:
        rows = []
        for gene_id in gene_ids:
            stable_id = strip_ensembl_version(gene_id)
            transcript_id, intervals, length = selected.get(stable_id, ("", {}, 0))
            mean_score, scored_bases = _mean_signal(bw, intervals) if intervals else (np.nan, 0)
            rows.append({"gene_id": stable_id, "transcript_id": transcript_id,
                         "utr3_length_nt": length if length else np.nan,
                         "utr3_mean_phastcons": mean_score, "utr3_scored_bases": scored_bases})
    finally:
        bw.close()
    return pd.DataFrame(rows, columns=FEATURE_COLUMNS)


def intersection_gene_set(predictions: dict[str, dict[str, pd.Series]], mirnas: set[str],
                          gene_universe: pd.Index) -> set[str]:
    universe = set(gene_universe.astype(str))
    covered_by_tool = []
    for by_mirna in predictions.values():
        covered = set()
        for mirna in mirnas:
            covered.update(str(gene_id) for gene_id in by_mirna.get(mirna, pd.Series(dtype=float)).index)
        covered_by_tool.append(covered & universe)
    return set.intersection(*covered_by_tool) if covered_by_tool else set()


def assign_coverage_groups(features: pd.DataFrame, gene_universe: pd.Index,
                           intersection: set[str]) -> pd.DataFrame:
    frame = pd.DataFrame({"gene_id": gene_universe.astype(str)}).merge(
        features.assign(gene_id=features.gene_id.astype(str)), on="gene_id", how="left", validate="one_to_one")
    frame.insert(1, "coverage_group", np.where(frame.gene_id.isin(intersection), "IGS", "non-IGS"))
    return frame


def summarize_features(frame: pd.DataFrame) -> pd.DataFrame:
    records = []
    labels = {"utr3_length_nt": "3′UTR length (nt)",
              "utr3_mean_phastcons": "Mean 3′UTR phastCons100way"}
    for column, label in labels.items():
        igs = pd.to_numeric(frame.loc[frame.coverage_group == "IGS", column], errors="coerce").dropna()
        non_igs = pd.to_numeric(frame.loc[frame.coverage_group == "non-IGS", column], errors="coerce").dropna()
        if len(igs) and len(non_igs):
            test = mannwhitneyu(igs, non_igs, alternative="two-sided", method="asymptotic")
            statistic, raw_pvalue = float(test.statistic), float(test.pvalue)
            log10_pvalue = _mann_whitney_log10_pvalue(igs.to_numpy(float), non_igs.to_numpy(float),
                                                       statistic, raw_pvalue)
            pvalue = raw_pvalue if raw_pvalue > 0 else np.nan
        else:
            statistic = pvalue = log10_pvalue = np.nan
        median_igs = float(igs.median()) if len(igs) else np.nan
        median_non_igs = float(non_igs.median()) if len(non_igs) else np.nan
        records.append({"feature": label, "feature_column": column,
                        "n_igs": int((frame.coverage_group == "IGS").sum()),
                        "n_non_igs": int((frame.coverage_group == "non-IGS").sum()),
                        "n_igs_with_value": len(igs), "n_non_igs_with_value": len(non_igs),
                        "median_igs": median_igs, "median_non_igs": median_non_igs,
                        "median_difference_igs_minus_non_igs": median_igs - median_non_igs,
                        "mann_whitney_u": statistic, "pvalue_two_sided": pvalue,
                        "log10_pvalue_two_sided": log10_pvalue,
                        "pvalue_two_sided_display": format_pvalue(pvalue, log10_pvalue)})
    return pd.DataFrame(records)
