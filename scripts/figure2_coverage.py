#!/usr/bin/env python3
"""Generate the coverage panels for manuscript Figure 2.

The script reads the ``joined.tsv`` files from one completed fun-miRBenvh run and
writes each panel as a separate publication-ready plot plus the exact summary
values used to draw it. Manuscript figures are written to
``manuscript_assets/figure2`` and manuscript tables are written to
``manuscript_assets/tables``.

Panel definitions
-----------------
A. Number of benchmark experiments retained for each predictor.
B. Predictor gene-set coverage across the full benchmark gene universe, plus
   the all-algorithm Intersection Gene Set (IGS).
C. Coverage of perturbation-consistent positive miRNA-gene pairs.
D. Coverage of background/non-positive miRNA-gene pairs.

Examples
--------
Generate all panels from a specific run::

    uv run python scripts/figure2_coverage.py \
        --run-dir results/20260709_122904 \
        --panel all

Outputs are written by default to::

    manuscript_assets/figure2/
    manuscript_assets/tables/
"""

from __future__ import annotations

import argparse
import gzip
import itertools
import logging
import math
import re
import shutil
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

from funmirbench.evaluate_common import CURVE_COLORS
from funmirbench.logger import setup_logging


logger = logging.getLogger(__name__)

DEFAULT_RESULTS_DIR = Path("results")
DEFAULT_METADATA_PATH = Path("metadata/predictions_info.tsv")
DEFAULT_MANUSCRIPT_OUTPUT_DIR = Path("manuscript_assets/figure2")
DEFAULT_MANUSCRIPT_TABLES_DIR = Path("manuscript_assets/tables")
DEFAULT_FORMATS = ("png", "svg")
PANEL_CHOICES = ("A", "B", "C", "D", "E", "F", "all")
DEFAULT_CONSERVATION_TABLE = Path("figure2F_utr_conservation_raw.tsv")
DEFAULT_GTF = Path("data/resources/ensembl/Homo_sapiens.GRCh38.115.gtf.gz")
ATTR_RE = re.compile(r'([A-Za-z0-9_]+) "([^"]*)"')
SVG_NS = "http://www.w3.org/2000/svg"
XLINK_NS = "http://www.w3.org/1999/xlink"
INTERSECTION_COLORS = {
    "all_algorithm_intersection": "#7A8798",
}
PREDICTOR_BAR_ALPHA = 0.35
IGS_BAR_ALPHA = 0.55
PANEL_LETTER_FONTSIZE = 18
PANEL_TITLE_FONTSIZE = 17
AXIS_LABEL_FONTSIZE = 14
TICK_LABEL_FONTSIZE = 12
BAR_LABEL_FONTSIZE = 16
ANNOTATION_FONTSIZE = 14
LEGEND_FONTSIZE = 13
GRID_ALPHA = 0.25

OE_ALIASES = {
    "OE",
    "OVEREXPRESSION",
    "OVER_EXPRESSION",
}
LOSS_ALIASES = {
    "KO",
    "KD",
    "KNOCKOUT",
    "KNOCK_OUT",
    "KNOCKDOWN",
    "KNOCK_DOWN",
}


@dataclass(frozen=True)
class BenchmarkDataset:
    """One benchmark dataset loaded from a completed run."""

    dataset_id: str
    frame: pd.DataFrame


@dataclass(frozen=True)
class Figure2Inputs:
    """Validated inputs shared by all Figure 2 panels."""

    run_dir: Path
    datasets: tuple[BenchmarkDataset, ...]
    tool_ids: tuple[str, ...]
    tool_labels: dict[str, str]
    tool_colors: dict[str, str]
    fdr_threshold: float | None
    effect_threshold: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate manuscript Figure 2 coverage panels from a fun-miRBenvh run."
    )
    parser.add_argument(
        "--run-dir",
        type=Path,
        default=None,
        help=(
            "Completed benchmark run directory. If omitted, the newest directory "
            "under results/ containing datasets/*/joined.tsv is used."
        ),
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=DEFAULT_RESULTS_DIR,
        help="Root used when auto-selecting the newest completed run.",
    )
    parser.add_argument(
        "--metadata",
        type=Path,
        default=DEFAULT_METADATA_PATH,
        help="Predictor metadata TSV used for display names and predictor order.",
    )
    parser.add_argument(
        "--panel",
        choices=PANEL_CHOICES,
        default="all",
        help="Generate one panel or all panels.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help=(
            "Optional diagnostic output directory for per-panel subdirectories. "
            "By default, only manuscript_assets outputs are written."
        ),
    )
    parser.add_argument(
        "--manuscript-out-dir",
        type=Path,
        default=DEFAULT_MANUSCRIPT_OUTPUT_DIR,
        help="Stable manuscript figure output directory. Default: manuscript_assets/figure2/.",
    )
    parser.add_argument(
        "--manuscript-tables-dir",
        type=Path,
        default=DEFAULT_MANUSCRIPT_TABLES_DIR,
        help="Stable manuscript table output directory. Default: manuscript_assets/tables/.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="Raster output resolution.",
    )
    parser.add_argument(
        "--gtf",
        type=Path,
        default=DEFAULT_GTF,
        help="GTF used to derive longest protein-coding transcript 3'UTR lengths.",
    )
    parser.add_argument(
        "--conservation-table",
        type=Path,
        default=None,
        help=(
            "Gene-level 3'UTR conservation table with gene_id and "
            "utr3_mean_conservation columns, used for panel F. Default: "
            "manuscript_assets/tables/figure2F_utr_conservation_raw.tsv."
        ),
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        help="Logging level. Default: INFO.",
    )
    return parser.parse_args()

def find_latest_completed_run(results_dir: Path) -> Path:
    if not results_dir.exists():
        raise FileNotFoundError(f"Results directory does not exist: {results_dir}")

    candidates = []
    for child in results_dir.iterdir():
        if child.is_dir() and list(child.glob("datasets/*/joined.tsv")):
            candidates.append(child)

    if not candidates:
        raise FileNotFoundError(
            f"No completed benchmark run with datasets/*/joined.tsv found under {results_dir}."
        )
    return max(candidates, key=lambda path: path.stat().st_mtime)


def load_run_thresholds(run_dir: Path) -> tuple[float | None, float]:
    config_path = run_dir / "benchmark_config.yaml"
    if not config_path.exists():
        raise FileNotFoundError(
            f"Completed run is missing benchmark config snapshot: {config_path}. "
            "Regenerate the benchmark run with the current pipeline so Figure 2 "
            "can reuse the run's ground-truth thresholds."
        )

    with config_path.open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle) or {}
    evaluation = config.get("evaluation") or {}
    raw_fdr_threshold = evaluation.get("fdr_threshold", 0.05)
    fdr_threshold = (
        None if raw_fdr_threshold is None else float(raw_fdr_threshold)
    )
    effect_threshold = float(evaluation.get("effect_threshold", 1.0))
    return fdr_threshold, effect_threshold


def default_conservation_table() -> Path:
    return DEFAULT_MANUSCRIPT_TABLES_DIR / DEFAULT_CONSERVATION_TABLE


def validate_panel_inputs(panels: tuple[str, ...], *, gtf_path: Path, conservation_table: Path) -> None:
    missing = []
    if "E" in panels and not gtf_path.exists():
        missing.append(f"Panel E GTF: {gtf_path}")
    if "F" in panels and not conservation_table.exists():
        missing.append(f"Panel F conservation table: {conservation_table}")
    if missing:
        details = "\n".join(f"- {item}" for item in missing)
        raise FileNotFoundError(
            "Missing required Figure 2 input files:\n"
            f"{details}\n"
            "For Panel F, generate the conservation table with:\n"
            "uv run python scripts/figure2_utr_conservation.py "
            f"--gtf {gtf_path} --out {conservation_table}"
        )


def load_tool_metadata(path: Path) -> tuple[tuple[str, ...], dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"Predictor metadata does not exist: {path}")

    metadata = pd.read_csv(path, sep="\t", dtype=str)
    required = {"tool_id", "official_name"}
    missing = required.difference(metadata.columns)
    if missing:
        raise ValueError(f"{path} is missing columns: {sorted(missing)}")

    metadata = metadata.dropna(subset=["tool_id"]).copy()
    metadata["tool_id"] = metadata["tool_id"].astype(str).str.strip()
    metadata["official_name"] = metadata["official_name"].fillna(metadata["tool_id"])
    metadata["official_name"] = metadata["official_name"].astype(str).str.strip()

    tool_ids = tuple(metadata["tool_id"].tolist())
    labels = dict(zip(metadata["tool_id"], metadata["official_name"]))
    return tool_ids, labels


def load_joined_datasets(run_dir: Path) -> tuple[BenchmarkDataset, ...]:
    joined_paths = sorted(run_dir.glob("datasets/*/joined.tsv"))
    if not joined_paths:
        raise FileNotFoundError(f"No datasets/*/joined.tsv files found in {run_dir}")

    datasets = []
    for path in joined_paths:
        frame = pd.read_csv(path, sep="\t", low_memory=False)
        dataset_id = path.parent.name
        if "dataset_id" in frame.columns and not frame.empty:
            dataset_id = str(frame["dataset_id"].iloc[0])
        datasets.append(BenchmarkDataset(dataset_id=dataset_id, frame=frame))
    return tuple(datasets)


def score_column(tool_id: str) -> str:
    return f"score_{tool_id}"


def available_tool_ids(
    datasets: Iterable[BenchmarkDataset], metadata_order: Iterable[str]
) -> tuple[str, ...]:
    datasets = tuple(datasets)
    common_score_columns: set[str] | None = None
    for dataset in datasets:
        current = {column for column in dataset.frame.columns if column.startswith("score_")}
        common_score_columns = (
            current if common_score_columns is None else common_score_columns & current
        )

    common_score_columns = common_score_columns or set()
    common_tool_ids = {column.removeprefix("score_") for column in common_score_columns}
    ordered = tuple(tool_id for tool_id in metadata_order if tool_id in common_tool_ids)
    extras = sorted(common_tool_ids.difference(ordered))
    resolved = ordered + tuple(extras)
    if not resolved:
        raise ValueError("No score_<tool_id> columns are shared by all joined datasets.")
    return resolved


def evaluation_tool_colors(tool_ids: Iterable[str]) -> dict[str, str]:
    """Reproduce the predictor colors used by the evaluation pipeline.

    ``evaluate_joined_dataframe`` discovers ``score_*`` columns in sorted order
    and assigns ``CURVE_COLORS`` sequentially. The mapping here deliberately
    follows that same ordering, independently of the display order in metadata.
    """

    evaluation_order = sorted(str(tool_id) for tool_id in tool_ids)
    return {
        tool_id: CURVE_COLORS[index % len(CURVE_COLORS)]
        for index, tool_id in enumerate(evaluation_order)
    }


def perturbation_series(frame: pd.DataFrame) -> pd.Series:
    if "perturbation" not in frame.columns:
        raise ValueError("Joined table is missing the perturbation column.")
    return (
        frame["perturbation"]
        .astype(str)
        .str.upper()
        .str.strip()
        .str.replace(r"[\s-]+", "_", regex=True)
    )


def annotate_ground_truth(
    frame: pd.DataFrame,
    *,
    fdr_threshold: float | None,
    effect_threshold: float,
) -> pd.DataFrame:
    required = {"gene_id", "logFC", "perturbation"}
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"Joined table is missing required columns: {sorted(missing)}")

    work = frame.copy()
    work["gene_id"] = work["gene_id"].astype(str)
    work["logFC"] = pd.to_numeric(work["logFC"], errors="coerce")

    perturbation = perturbation_series(work)
    oe_mask = perturbation.isin(OE_ALIASES)
    loss_mask = perturbation.isin(LOSS_ALIASES)
    unsupported = ~(oe_mask | loss_mask)
    if unsupported.any():
        values = sorted(perturbation.loc[unsupported].dropna().unique().tolist())
        raise ValueError(
            "Unsupported perturbation values after normalization: "
            f"{values}. Supported OE aliases: {sorted(OE_ALIASES)}; "
            f"loss-of-function aliases: {sorted(LOSS_ALIASES)}"
        )

    work["expected_effect"] = np.where(oe_mask, -work["logFC"], work["logFC"])
    usable = work["logFC"].notna()

    if fdr_threshold is not None:
        fdr_column = "benchmark_FDR" if "benchmark_FDR" in work.columns else "FDR"
        if fdr_column not in work.columns:
            raise ValueError(
                "An FDR threshold was requested, but neither benchmark_FDR nor FDR is available."
            )
        work[fdr_column] = pd.to_numeric(work[fdr_column], errors="coerce")
        usable &= work[fdr_column].between(0.0, 1.0, inclusive="both")
        significant = work[fdr_column] < float(fdr_threshold)
    else:
        significant = pd.Series(True, index=work.index)

    work = work.loc[usable].copy()
    significant = significant.loc[work.index]
    work["is_positive"] = (
        (work["expected_effect"] > float(effect_threshold)) & significant
    ).astype(bool)
    work["is_background"] = ~work["is_positive"]
    return work


def prepare_inputs(
    *,
    run_dir: Path,
    metadata_path: Path,
    fdr_threshold: float | None,
    effect_threshold: float,
) -> Figure2Inputs:
    datasets = load_joined_datasets(run_dir)
    metadata_order, labels = load_tool_metadata(metadata_path)
    tool_ids = available_tool_ids(datasets, metadata_order)

    annotated = tuple(
        BenchmarkDataset(
            dataset_id=dataset.dataset_id,
            frame=annotate_ground_truth(
                dataset.frame,
                fdr_threshold=fdr_threshold,
                effect_threshold=effect_threshold,
            ),
        )
        for dataset in datasets
    )
    resolved_labels = {
        tool_id: labels.get(tool_id, tool_id.replace("_", " "))
        for tool_id in tool_ids
    }
    return Figure2Inputs(
        run_dir=run_dir,
        datasets=annotated,
        tool_ids=tool_ids,
        tool_labels=resolved_labels,
        tool_colors=evaluation_tool_colors(tool_ids),
        fdr_threshold=fdr_threshold,
        effect_threshold=effect_threshold,
    )


def style_axis(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", alpha=GRID_ALPHA, linewidth=0.8)
    ax.set_axisbelow(True)
    ax.tick_params(axis="both", labelsize=TICK_LABEL_FONTSIZE)


def set_panel_title(ax: plt.Axes, letter: str, title: str) -> None:
    ax.text(
        -0.20,
        1.18,
        letter,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=PANEL_LETTER_FONTSIZE,
        fontweight="bold",
        clip_on=False,
    )
    ax.set_title(title, pad=18, fontsize=PANEL_TITLE_FONTSIZE, fontweight="bold")


def save_panel(
    fig: plt.Figure,
    *,
    out_dir: Path,
    stem: str,
    dpi: int,
) -> list[Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    for extension in DEFAULT_FORMATS:
        path = out_dir / f"{stem}.{extension}"
        fig.savefig(path, dpi=dpi, bbox_inches="tight", pad_inches=0.08)
        paths.append(path)
    plt.close(fig)
    return paths


def load_experiment_metric_table(run_dir: Path) -> pd.DataFrame | None:
    """Load the per-experiment coverage table used by the evaluation report.

    Missing values identify experiments that were not evaluable for a predictor,
    which is the definition used for panel A. This preserves the 31/29/27 counts
    reported in the manuscript draft instead of counting any raw score presence.
    """

    path = run_dir / "tables" / "per_experiment" / "coverage_per_experiment.tsv"
    if not path.exists():
        return None
    return pd.read_csv(path, sep="\t", low_memory=False)


def panel_a_experiment_coverage(inputs: Figure2Inputs) -> tuple[pd.DataFrame, plt.Figure]:
    metric_table = load_experiment_metric_table(inputs.run_dir)
    rows = []
    total_experiments = len(inputs.datasets)

    for tool_id in inputs.tool_ids:
        if metric_table is not None and tool_id in metric_table.columns:
            retained_mask = metric_table[tool_id].notna()
            retained = int(retained_mask.sum())
            excluded_datasets = (
                metric_table.loc[~retained_mask, "dataset_id"].astype(str).tolist()
                if "dataset_id" in metric_table.columns
                else []
            )
        else:
            retained = sum(
                bool(dataset.frame[score_column(tool_id)].notna().any())
                for dataset in inputs.datasets
            )
            excluded_datasets = []

        rows.append(
            {
                "tool_id": tool_id,
                "predictor": inputs.tool_labels[tool_id],
                "experiments_retained": retained,
                "experiments_total": total_experiments,
                "experiments_excluded": total_experiments - retained,
                "experiment_coverage": (
                    retained / total_experiments if total_experiments else np.nan
                ),
                "excluded_dataset_ids": ";".join(excluded_datasets),
            }
        )

    summary = pd.DataFrame(rows)
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    x = np.arange(len(summary))
    colors = [inputs.tool_colors[tool_id] for tool_id in summary["tool_id"]]
    bars = ax.bar(
        x,
        summary["experiments_retained"].to_numpy(),
        color=colors,
        edgecolor=colors,
        alpha=PREDICTOR_BAR_ALPHA,
        linewidth=1.0,
    )
    ax.axhline(total_experiments, linestyle="--", linewidth=1.2, color="black")
    ax.set_xticks(x, summary["predictor"], rotation=25, ha="right")
    ax.set_ylabel("Experiments retained", fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_ylim(0, max(total_experiments * 1.16, 1))
    set_panel_title(ax, "A", "Experiment coverage")
    for bar, retained in zip(bars, summary["experiments_retained"]):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height(),
            f"{int(retained)}/{total_experiments}",
            ha="center",
            va="bottom",
            fontsize=BAR_LABEL_FONTSIZE,
        )
    style_axis(ax)
    fig.tight_layout()
    return summary, fig


def panel_b_gene_set_coverage(inputs: Figure2Inputs) -> tuple[pd.DataFrame, plt.Figure]:
    all_pairs = pd.concat(
        [
            dataset.frame.assign(dataset_id=dataset.dataset_id)
            for dataset in inputs.datasets
        ],
        ignore_index=True,
    )
    fgs_genes = set(all_pairs["gene_id"].dropna().astype(str))
    scored_gene_sets = {
        tool_id: set(
            all_pairs.loc[all_pairs[score_column(tool_id)].notna(), "gene_id"]
            .dropna()
            .astype(str)
        )
        for tool_id in inputs.tool_ids
    }
    rows = []
    denominator = len(fgs_genes)
    for tool_id in inputs.tool_ids:
        count = len(scored_gene_sets[tool_id])
        rows.append(
            {
                "set_id": tool_id,
                "label": inputs.tool_labels[tool_id],
                "gene_count": count,
                "fgs_gene_count": denominator,
                "fgs_coverage": count / denominator if denominator else np.nan,
                "set_type": "predictor",
            }
        )

    available_sets = [scored_gene_sets[tool_id] for tool_id in inputs.tool_ids]
    igs_genes = set.intersection(*available_sets) if available_sets else set()
    rows.append(
        {
            "set_id": "all_algorithm_intersection",
            "label": "IGS",
            "gene_count": len(igs_genes),
            "fgs_gene_count": denominator,
            "fgs_coverage": len(igs_genes) / denominator if denominator else np.nan,
            "set_type": "intersection",
        }
    )
    summary = pd.DataFrame(rows)

    colors = [
        inputs.tool_colors[set_id]
        if set_id in inputs.tool_colors
        else INTERSECTION_COLORS[set_id]
        for set_id in summary["set_id"]
    ]
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    x = np.arange(len(summary))
    alphas = [
        PREDICTOR_BAR_ALPHA if set_type == "predictor" else IGS_BAR_ALPHA
        for set_type in summary["set_type"]
    ]
    bars = ax.bar(
        x,
        100.0 * summary["fgs_coverage"].to_numpy(),
        color=colors,
        edgecolor=colors,
        linewidth=1.0,
    )
    for bar, alpha in zip(bars, alphas):
        bar.set_alpha(alpha)
    ax.set_xticks(x, summary["label"], rotation=25, ha="right")
    ax.set_ylabel("Genes covered (%)", fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_ylim(0, 110)
    set_panel_title(ax, "B", "Gene-set coverage and IGS")
    ax.text(
        0.01,
        0.96,
        f"FGS: {denominator:,}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=ANNOTATION_FONTSIZE,
    )
    for bar, percentage in zip(bars, summary["fgs_coverage"]):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height(),
            f"{percentage:.1%}",
            ha="center",
            va="bottom",
            fontsize=BAR_LABEL_FONTSIZE,
        )
    style_axis(ax)
    fig.tight_layout()
    return summary, fig


def strip_version(gene_id: object) -> str:
    return str(gene_id).split(".", 1)[0]


def parse_attributes(attributes: str) -> dict[str, str]:
    return dict(ATTR_RE.findall(attributes))


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("rt")


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not intervals:
        return []
    ordered = sorted(intervals)
    merged = [ordered[0]]
    for start, end in ordered[1:]:
        prev_start, prev_end = merged[-1]
        if start <= prev_end + 1:
            merged[-1] = (prev_start, max(prev_end, end))
        else:
            merged.append((start, end))
    return merged


def interval_length(intervals: list[tuple[int, int]]) -> int:
    return sum(end - start + 1 for start, end in merge_intervals(intervals))


def load_max_utr_lengths(gtf_path: Path) -> pd.DataFrame:
    if not gtf_path.exists():
        raise FileNotFoundError(f"GTF does not exist: {gtf_path}")

    transcript_intervals: dict[tuple[str, str], list[tuple[int, int]]] = {}
    with open_text(gtf_path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "three_prime_utr":
                continue
            attrs = parse_attributes(fields[8])
            gene_id = attrs.get("gene_id")
            transcript_id = attrs.get("transcript_id")
            if not gene_id or not transcript_id:
                continue
            if attrs.get("transcript_biotype") != "protein_coding":
                continue
            key = (strip_version(gene_id), transcript_id)
            transcript_intervals.setdefault(key, []).append(
                (int(fields[3]), int(fields[4]))
            )

    gene_lengths: dict[str, int] = {}
    for (gene_id, _transcript_id), intervals in transcript_intervals.items():
        length = interval_length(intervals)
        gene_lengths[gene_id] = max(gene_lengths.get(gene_id, 0), length)

    return pd.DataFrame(
        {
            "gene_id": list(gene_lengths.keys()),
            "utr3_length": list(gene_lengths.values()),
        }
    )


def load_conservation_table(path: Path) -> pd.DataFrame:
    sep = "," if path.suffix.lower() == ".csv" else "\t"
    table = pd.read_csv(path, sep=sep)
    required = {"gene_id", "utr3_mean_conservation"}
    missing = required.difference(table.columns)
    if missing:
        raise ValueError(f"{path} is missing columns: {sorted(missing)}")
    out = table[["gene_id", "utr3_mean_conservation"]].copy()
    out["gene_id"] = out["gene_id"].map(strip_version)
    out["utr3_mean_conservation"] = pd.to_numeric(
        out["utr3_mean_conservation"],
        errors="coerce",
    )
    return out.dropna(subset=["gene_id", "utr3_mean_conservation"])


def annotate_sets(values: pd.DataFrame, gene_sets: SimpleNamespace) -> pd.DataFrame:
    work = values.copy()
    work["gene_id"] = work["gene_id"].map(strip_version)
    work = work[work["gene_id"].isin(gene_sets.fgs)].copy()
    work["gene_set"] = np.where(work["gene_id"].isin(gene_sets.igs), "IGS", "non-IGS")
    return work


def mann_whitney_pvalue(a: pd.Series, b: pd.Series) -> float:
    values = pd.concat(
        [
            pd.DataFrame({"value": a.astype(float), "group": "a"}),
            pd.DataFrame({"value": b.astype(float), "group": "b"}),
        ],
        ignore_index=True,
    ).dropna()
    n1 = int((values["group"] == "a").sum())
    n2 = int((values["group"] == "b").sum())
    if n1 == 0 or n2 == 0:
        return float("nan")

    values["rank"] = values["value"].rank(method="average")
    rank_sum_a = float(values.loc[values["group"] == "a", "rank"].sum())
    u1 = rank_sum_a - n1 * (n1 + 1) / 2.0
    mean_u = n1 * n2 / 2.0
    tie_counts = values["value"].value_counts()
    tie_term = float(((tie_counts**3) - tie_counts).sum())
    n = n1 + n2
    variance = n1 * n2 / 12.0 * ((n + 1) - tie_term / (n * (n - 1)))
    if variance <= 0:
        return float("nan")
    z = (u1 - mean_u) / math.sqrt(variance)
    return math.erfc(abs(z) / math.sqrt(2.0))


def metric_stats(frame: pd.DataFrame, *, metric: str, label: str) -> dict[str, object]:
    igs = frame.loc[frame["gene_set"] == "IGS", metric].dropna()
    non = frame.loc[frame["gene_set"] == "non-IGS", metric].dropna()
    pvalue = mann_whitney_pvalue(igs, non)
    return {
        "metric": label,
        "n_igs": int(igs.size),
        "n_non_igs": int(non.size),
        "mean_igs": float(igs.mean()) if not igs.empty else np.nan,
        "mean_non_igs": float(non.mean()) if not non.empty else np.nan,
        "mean_difference_igs_minus_non_igs": (
            float(igs.mean() - non.mean()) if not igs.empty and not non.empty else np.nan
        ),
        "median_igs": float(igs.median()) if not igs.empty else np.nan,
        "median_non_igs": float(non.median()) if not non.empty else np.nan,
        "mann_whitney_pvalue_two_sided_normal_approx": pvalue,
        "pvalue_label": "<1e-308" if pvalue == 0.0 else f"{pvalue:.3g}",
    }


def smoothed_density(values: pd.Series, *, bins: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    counts, edges = np.histogram(values.dropna().astype(float), bins=bins, density=True)
    centers = (edges[:-1] + edges[1:]) / 2.0
    kernel_x = np.linspace(-3.0, 3.0, 21)
    kernel = np.exp(-0.5 * kernel_x**2)
    kernel /= kernel.sum()
    density = np.convolve(counts, kernel, mode="same")
    return centers, density


def plot_metric(
    ax: plt.Axes,
    frame: pd.DataFrame,
    *,
    metric: str,
    title: str | None,
    xlabel: str,
    x_upper: float | None = None,
) -> None:
    non = frame.loc[frame["gene_set"] == "non-IGS", metric].dropna().astype(float)
    igs = frame.loc[frame["gene_set"] == "IGS", metric].dropna().astype(float)
    combined = pd.concat([non, igs], ignore_index=True)
    upper = float(x_upper) if x_upper is not None else float(combined.max())
    lower = max(0.0, float(combined.min())) if (combined >= 0).all() else float(combined.min())
    bins = np.linspace(lower, upper, 90)

    colors = {
        "non-IGS": "#8FB5D6",
        "IGS": "#D8766A",
    }
    igs_x, igs_y = smoothed_density(igs.clip(lower=lower, upper=upper), bins=bins)
    non_x, non_y = smoothed_density(non.clip(lower=lower, upper=upper), bins=bins)
    scale = 0.52 / max(float(igs_y.max()), float(non_y.max()), 1e-12)
    igs_y = igs_y * scale
    non_y = non_y * scale

    ax.fill_between(
        igs_x,
        0,
        igs_y,
        color=colors["IGS"],
        alpha=0.8,
        label=f"IGS genes (n={len(igs):,})",
    )
    ax.plot(igs_x, igs_y, color=colors["IGS"], linewidth=1.2)
    ax.fill_between(
        non_x,
        0,
        -non_y,
        color=colors["non-IGS"],
        alpha=0.8,
        label=f"non-IGS genes (n={len(non):,})",
    )
    ax.plot(non_x, -non_y, color=colors["non-IGS"], linewidth=1.2)
    ax.axhline(0, color="#888888", linewidth=0.8)
    ax.axvline(igs.median(), ymin=0.50, ymax=0.88, color="black", linewidth=1.2, linestyle="--", alpha=0.95)
    ax.axvline(non.median(), ymin=0.12, ymax=0.50, color="black", linewidth=1.2, linestyle="--", alpha=0.95)
    ax.set_xlabel(xlabel, fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_xlim(lower, upper)
    ax.set_ylim(-0.62, 0.62)
    ax.set_yticks([0.28, -0.28], ["IGS genes", "non-IGS genes"])
    ax.tick_params(axis="both", labelsize=TICK_LABEL_FONTSIZE)
    ax.set_ylabel("")
    if title:
        ax.set_title(title, pad=18, fontsize=PANEL_TITLE_FONTSIZE, fontweight="bold")
    ax.legend(frameon=False, fontsize=LEGEND_FONTSIZE, loc="upper right")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.grid(axis="x", alpha=GRID_ALPHA, linewidth=0.8)
    ax.set_axisbelow(True)


def figure2_gene_sets(inputs: Figure2Inputs) -> SimpleNamespace:
    all_pairs = pd.concat(
        [
            dataset.frame.assign(dataset_id=dataset.dataset_id)
            for dataset in inputs.datasets
        ],
        ignore_index=True,
    )
    fgs_genes = set(all_pairs["gene_id"].dropna().astype(str))
    scored_gene_sets = {
        tool_id: set(
            all_pairs.loc[all_pairs[score_column(tool_id)].notna(), "gene_id"]
            .dropna()
            .astype(str)
        )
        for tool_id in inputs.tool_ids
    }
    igs_genes = (
        set.intersection(*scored_gene_sets.values()) if scored_gene_sets else set()
    )
    return SimpleNamespace(fgs=fgs_genes, igs=igs_genes)


def panel_e_utr_length(
    inputs: Figure2Inputs,
    *,
    gtf_path: Path,
) -> tuple[pd.DataFrame, plt.Figure]:
    gene_sets = figure2_gene_sets(inputs)
    lengths = annotate_sets(load_max_utr_lengths(gtf_path), gene_sets)
    stats = pd.DataFrame(
        [metric_stats(lengths, metric="utr3_length", label="3'UTR length")]
    )
    stats.insert(1, "fgs_genes", len(gene_sets.fgs))
    stats.insert(2, "igs_genes", len(gene_sets.igs))
    stats.insert(3, "non_igs_genes", len(gene_sets.fgs - gene_sets.igs))

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    plot_metric(
        ax,
        lengths,
        metric="utr3_length",
        title=None,
        xlabel="3'UTR length (nt; displayed to 10,000 nt)",
        x_upper=10000.0,
    )
    set_panel_title(ax, "E", "3'UTR length")
    fig.tight_layout()
    return stats, fig


def panel_f_utr_conservation(
    inputs: Figure2Inputs,
    *,
    conservation_table: Path,
) -> tuple[pd.DataFrame, plt.Figure]:
    if not conservation_table.exists():
        raise FileNotFoundError(
            "Panel F requires a gene-level conservation table. "
            f"Missing: {conservation_table}"
        )
    gene_sets = figure2_gene_sets(inputs)
    conservation = annotate_sets(
        load_conservation_table(conservation_table),
        gene_sets,
    )
    stats = pd.DataFrame(
        [
            metric_stats(
                conservation,
                metric="utr3_mean_conservation",
                label="3'UTR mean conservation",
            )
        ]
    )
    stats.insert(1, "fgs_genes", len(gene_sets.fgs))
    stats.insert(2, "igs_genes", len(gene_sets.igs))
    stats.insert(3, "non_igs_genes", len(gene_sets.fgs - gene_sets.igs))

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    plot_metric(
        ax,
        conservation,
        metric="utr3_mean_conservation",
        title=None,
        xlabel="Mean 3'UTR conservation score",
        x_upper=1.0,
    )
    set_panel_title(ax, "F", "3'UTR conservation")
    fig.tight_layout()
    return stats, fig


def pooled_pair_coverage(
    inputs: Figure2Inputs,
    *,
    subset_column: str,
    subset_label: str,
) -> pd.DataFrame:
    rows = []
    for tool_id in inputs.tool_ids:
        denominator = 0
        numerator = 0
        for dataset in inputs.datasets:
            subset = dataset.frame[subset_column].astype(bool)
            denominator += int(subset.sum())
            numerator += int(
                dataset.frame.loc[subset, score_column(tool_id)].notna().sum()
            )
        rows.append(
            {
                "tool_id": tool_id,
                "predictor": inputs.tool_labels[tool_id],
                f"{subset_label}_pairs_scored": numerator,
                f"{subset_label}_pairs_total": denominator,
                f"{subset_label}_coverage": (
                    numerator / denominator if denominator else np.nan
                ),
            }
        )
    return pd.DataFrame(rows)


def pair_coverage_figure(
    summary: pd.DataFrame,
    *,
    inputs: Figure2Inputs,
    coverage_column: str,
    numerator_column: str,
    denominator_column: str,
    title: str,
    ylabel: str,
    total_label: str,
    panel: str,
) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    x = np.arange(len(summary))
    percentages = 100.0 * summary[coverage_column].to_numpy()
    colors = [inputs.tool_colors[tool_id] for tool_id in summary["tool_id"]]
    bars = ax.bar(
        x,
        percentages,
        color=colors,
        edgecolor=colors,
        alpha=PREDICTOR_BAR_ALPHA,
        linewidth=1.0,
    )
    ax.set_xticks(x, summary["predictor"], rotation=25, ha="right")
    ax.set_ylabel(ylabel, fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_ylim(0, 110)
    set_panel_title(ax, panel, title)
    denominator = int(summary[denominator_column].iloc[0]) if not summary.empty else 0
    ax.text(
        0.01,
        0.96,
        f"{total_label}: {denominator:,}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=ANNOTATION_FONTSIZE,
    )
    for bar, percentage in zip(bars, summary[coverage_column]):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height(),
            f"{percentage:.1%}",
            ha="center",
            va="bottom",
            fontsize=BAR_LABEL_FONTSIZE,
        )
    style_axis(ax)
    fig.tight_layout()
    return fig


def panel_c_positive_coverage(inputs: Figure2Inputs) -> tuple[pd.DataFrame, plt.Figure]:
    summary = pooled_pair_coverage(
        inputs,
        subset_column="is_positive",
        subset_label="positive",
    )
    fig = pair_coverage_figure(
        summary,
        inputs=inputs,
        coverage_column="positive_coverage",
        numerator_column="positive_pairs_scored",
        denominator_column="positive_pairs_total",
        title="Positive miRNA-gene pair coverage",
        ylabel="Positive miRNA-gene pairs covered (%)",
        total_label="Total positive",
        panel="C",
    )
    return summary, fig


def panel_d_background_coverage(inputs: Figure2Inputs) -> tuple[pd.DataFrame, plt.Figure]:
    summary = pooled_pair_coverage(
        inputs,
        subset_column="is_background",
        subset_label="background",
    )
    fig = pair_coverage_figure(
        summary,
        inputs=inputs,
        coverage_column="background_coverage",
        numerator_column="background_pairs_scored",
        denominator_column="background_pairs_total",
        title="Background miRNA-gene pair coverage",
        ylabel="Background miRNA-gene pairs covered (%)",
        total_label="Total background",
        panel="D",
    )
    return summary, fig


def predictor_mean_mirnas_per_scored_gene(inputs: Figure2Inputs) -> pd.DataFrame:
    all_pairs = pd.concat(
        [
            dataset.frame.assign(dataset_id=dataset.dataset_id)
            for dataset in inputs.datasets
        ],
        ignore_index=True,
    )
    rows = []
    for tool_id in inputs.tool_ids:
        scored = (
            all_pairs.loc[
                all_pairs[score_column(tool_id)].notna(),
                ["gene_id", "mirna"],
            ]
            .dropna()
            .copy()
        )
        total_pairs = int(scored.shape[0])
        per_gene = (
            scored.drop_duplicates(["gene_id", "mirna"])
            .groupby("gene_id")["mirna"]
            .nunique()
        )
        rows.append(
            {
                "tool_id": tool_id,
                "total_scored_pairs": total_pairs,
                "mean_mirnas_per_scored_gene": (
                    float(per_gene.mean()) if len(per_gene) else np.nan
                ),
            }
        )
    return pd.DataFrame(rows)


def write_supplementary_coverage_table(
    *,
    inputs: Figure2Inputs,
    manuscript_tables_dir: Path,
) -> Path:
    experiment, experiment_fig = panel_a_experiment_coverage(inputs)
    plt.close(experiment_fig)
    experiment = experiment.rename(
        columns={
            "experiments_retained": "panel_a_experiments_retained",
            "experiments_total": "panel_a_experiments_total",
            "experiments_excluded": "panel_a_experiments_excluded",
            "experiment_coverage": "panel_a_experiment_coverage",
        }
    )
    gene, gene_fig = panel_b_gene_set_coverage(inputs)
    plt.close(gene_fig)
    gene = gene.rename(
        columns={
            "set_id": "tool_id",
            "label": "predictor",
            "gene_count": "panel_b_gene_count",
            "fgs_gene_count": "panel_b_fgs_gene_count",
            "fgs_coverage": "panel_b_fgs_coverage",
        }
    )
    positive, positive_fig = panel_c_positive_coverage(inputs)
    plt.close(positive_fig)
    positive = positive.rename(
        columns={
            "positive_pairs_scored": "panel_c_positive_pairs_scored",
            "positive_pairs_total": "panel_c_positive_pairs_total",
            "positive_coverage": "panel_c_positive_pair_coverage",
        }
    )
    background, background_fig = panel_d_background_coverage(inputs)
    plt.close(background_fig)
    background = background.rename(
        columns={
            "background_pairs_scored": "panel_d_background_pairs_scored",
            "background_pairs_total": "panel_d_background_pairs_total",
            "background_coverage": "panel_d_background_pair_coverage",
        }
    )
    mirnas_per_gene = predictor_mean_mirnas_per_scored_gene(inputs)

    for frame in (experiment, positive, background):
        if "predictor" in frame.columns:
            frame.drop(columns=["predictor"], inplace=True)

    table = gene.merge(experiment, on="tool_id", how="left")
    table = table.merge(positive, on="tool_id", how="left")
    table = table.merge(background, on="tool_id", how="left")
    table = table.merge(mirnas_per_gene, on="tool_id", how="left")
    table.insert(
        1,
        "name",
        table["tool_id"].map(inputs.tool_labels).fillna(table["predictor"]),
    )
    table.loc[
        table["tool_id"] == "all_algorithm_intersection",
        "name",
    ] = "IGS (all algorithms)"

    columns = [
        "name",
        "set_type",
        "panel_a_experiments_retained",
        "panel_a_experiments_total",
        "panel_a_experiments_excluded",
        "panel_a_experiment_coverage",
        "panel_b_gene_count",
        "panel_b_fgs_gene_count",
        "panel_b_fgs_coverage",
        "panel_c_positive_pairs_scored",
        "panel_c_positive_pairs_total",
        "panel_c_positive_pair_coverage",
        "panel_d_background_pairs_scored",
        "panel_d_background_pairs_total",
        "panel_d_background_pair_coverage",
        "total_scored_pairs",
        "mean_mirnas_per_scored_gene",
    ]
    table = table[columns]
    path = manuscript_tables_dir / "figure2_supplementary_coverage_table.tsv"
    manuscript_tables_dir.mkdir(parents=True, exist_ok=True)
    table.to_csv(path, sep="\t", index=False)
    logger.info("Wrote supplementary Figure 2 table: %s", path)
    return path


def write_gene_set_overlap_table(
    *,
    inputs: Figure2Inputs,
    manuscript_tables_dir: Path,
) -> Path:
    all_pairs = pd.concat(
        [
            dataset.frame.assign(dataset_id=dataset.dataset_id)
            for dataset in inputs.datasets
        ],
        ignore_index=True,
    )
    scored_gene_sets = {
        tool_id: set(
            all_pairs.loc[all_pairs[score_column(tool_id)].notna(), "gene_id"]
            .dropna()
            .astype(str)
        )
        for tool_id in inputs.tool_ids
    }

    rows = []
    for combination_size in range(2, len(inputs.tool_ids) + 1):
        for combination in itertools.combinations(inputs.tool_ids, combination_size):
            gene_sets = [scored_gene_sets[tool_id] for tool_id in combination]
            intersection = set.intersection(*gene_sets)
            union = set.union(*gene_sets)
            predictor_sizes = {
                tool_id: len(scored_gene_sets[tool_id])
                for tool_id in combination
            }
            n_intersection = len(intersection)
            n_union = len(union)
            row = {
                "combination_size": combination_size,
                "predictors": ";".join(
                    inputs.tool_labels[tool_id] for tool_id in combination
                ),
                "n_predictors": ";".join(
                    str(predictor_sizes[tool_id]) for tool_id in combination
                ),
                "intersection": n_intersection,
                "union": n_union,
                "jaccard": n_intersection / n_union if n_union else np.nan,
                "overlap_coefficient": (
                    n_intersection / min(predictor_sizes.values())
                    if predictor_sizes and min(predictor_sizes.values())
                    else np.nan
                ),
            }
            if combination_size == 2:
                genes_a, genes_b = gene_sets
                row["only_first"] = len(genes_a - genes_b)
                row["only_second"] = len(genes_b - genes_a)
            else:
                row["only_first"] = np.nan
                row["only_second"] = np.nan
            for tool_id in inputs.tool_ids:
                if tool_id in predictor_sizes:
                    row[f"fraction_intersection_in_{tool_id}"] = (
                        n_intersection / predictor_sizes[tool_id]
                        if predictor_sizes[tool_id]
                        else np.nan
                    )
                else:
                    row[f"fraction_intersection_in_{tool_id}"] = np.nan
            rows.append(row)

    path = manuscript_tables_dir / "figure2_gene_set_overlap.tsv"
    manuscript_tables_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
    logger.info("Wrote Figure 2 gene-set overlap table: %s", path)
    return path


def write_gene_set_overlap_supplement_table(
    *,
    inputs: Figure2Inputs,
    manuscript_tables_dir: Path,
) -> Path:
    all_pairs = pd.concat(
        [
            dataset.frame.assign(dataset_id=dataset.dataset_id)
            for dataset in inputs.datasets
        ],
        ignore_index=True,
    )
    scored_gene_sets = {
        tool_id: set(
            all_pairs.loc[all_pairs[score_column(tool_id)].notna(), "gene_id"]
            .dropna()
            .astype(str)
        )
        for tool_id in inputs.tool_ids
    }
    mirna_means = predictor_mean_mirnas_per_scored_gene(inputs).set_index("tool_id")[
        "mean_mirnas_per_scored_gene"
    ]

    rows = []
    for combination_size in range(1, len(inputs.tool_ids) + 1):
        for combination in itertools.combinations(inputs.tool_ids, combination_size):
            gene_sets = [scored_gene_sets[tool_id] for tool_id in combination]
            intersection = set.intersection(*gene_sets)
            union = set.union(*gene_sets)
            predictor_sizes = {
                tool_id: len(scored_gene_sets[tool_id])
                for tool_id in combination
            }
            combination_means = {
                tool_id: float(mirna_means.loc[tool_id])
                for tool_id in combination
            }
            mean_values = list(combination_means.values())
            n_intersection = len(intersection)
            n_union = len(union)
            rows.append(
                {
                    "predictor_set": ";".join(
                        inputs.tool_labels[tool_id] for tool_id in combination
                    ),
                    "combination_size": combination_size,
                    "genes_per_predictor": ";".join(
                        f"{inputs.tool_labels[tool_id]}={predictor_sizes[tool_id]}"
                        for tool_id in combination
                    ),
                    "intersection_gene_count": n_intersection,
                    "union_gene_count": n_union,
                    "jaccard_index": n_intersection / n_union if n_union else np.nan,
                    "overlap_coefficient": (
                        n_intersection / min(predictor_sizes.values())
                        if predictor_sizes and min(predictor_sizes.values())
                        else np.nan
                    ),
                    "mean_mirnas_per_scored_gene_by_predictor": ";".join(
                        (
                            f"{inputs.tool_labels[tool_id]}="
                            f"{combination_means[tool_id]:.2f}"
                        )
                        for tool_id in combination
                    ),
                    "mean_mirnas_per_scored_gene_range": (
                        f"{min(mean_values):.2f}-{max(mean_values):.2f}"
                        if len(mean_values) > 1
                        else f"{mean_values[0]:.2f}"
                    ),
                }
            )

    path = manuscript_tables_dir / "figure2_supplement_gene_set_overlap_mirna_coverage.tsv"
    manuscript_tables_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
    logger.info("Wrote Figure 2 gene-set overlap supplement table: %s", path)
    return path


def write_panel(
    panel: str,
    *,
    inputs: Figure2Inputs,
    out_dir: Path | None,
    manuscript_figures_dir: Path,
    manuscript_tables_dir: Path,
    gtf_path: Path,
    conservation_table: Path,
    dpi: int,
) -> dict[str, Path]:
    builders = {
        "A": panel_a_experiment_coverage,
        "B": panel_b_gene_set_coverage,
        "C": panel_c_positive_coverage,
        "D": panel_d_background_coverage,
    }
    if panel == "E":
        summary, fig = panel_e_utr_length(inputs, gtf_path=gtf_path)
    elif panel == "F":
        summary, fig = panel_f_utr_conservation(
            inputs,
            conservation_table=conservation_table,
        )
    else:
        summary, fig = builders[panel](inputs)
    stem = {
        "A": "figure2A_experiment_coverage",
        "B": "figure2B_gene_set_coverage",
        "C": "figure2C_positive_coverage",
        "D": "figure2D_background_coverage",
        "E": "figure2E_utr_length",
        "F": "figure2F_utr_conservation",
    }[panel]
    manuscript_figures_dir.mkdir(parents=True, exist_ok=True)
    manuscript_tables_dir.mkdir(parents=True, exist_ok=True)
    manuscript_table = manuscript_tables_dir / f"{stem}.tsv"
    summary.to_csv(manuscript_table, sep="\t", index=False)
    manuscript_figure_paths = save_panel(
        fig,
        out_dir=manuscript_figures_dir,
        stem=stem,
        dpi=dpi,
    )
    manuscript_figures = {
        figure_path.suffix.lstrip("."): figure_path
        for figure_path in manuscript_figure_paths
    }
    if out_dir is not None:
        panel_dir = out_dir / f"panel_{panel}"
        panel_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(manuscript_table, panel_dir / f"{stem}.tsv")
        for figure_path in manuscript_figure_paths:
            shutil.copy2(figure_path, panel_dir / figure_path.name)
        logger.info("Wrote panel %s diagnostics: %s", panel, panel_dir / stem)
    logger.info("Wrote panel %s: %s", panel, manuscript_figures_dir / stem)
    return {
        "table": manuscript_table,
        **manuscript_figures,
    }


def trim_white_border(image: np.ndarray, *, threshold: float = 0.985, pad: int = 16) -> np.ndarray:
    rgb = image[..., :3]
    nonwhite = np.any(rgb < threshold, axis=2)
    if not np.any(nonwhite):
        return image
    rows = np.where(nonwhite.any(axis=1))[0]
    cols = np.where(nonwhite.any(axis=0))[0]
    top = max(int(rows[0]) - pad, 0)
    bottom = min(int(rows[-1]) + pad + 1, image.shape[0])
    left = max(int(cols[0]) - pad, 0)
    right = min(int(cols[-1]) + pad + 1, image.shape[1])
    return image[top:bottom, left:right]


def pad_to_shape(
    image: np.ndarray,
    shape: tuple[int, int],
    *,
    vertical: str = "center",
    horizontal: str = "center",
) -> np.ndarray:
    target_height, target_width = shape
    height, width = image.shape[:2]
    if height > target_height or width > target_width:
        raise ValueError("Target shape must be at least as large as the image.")

    channels = image.shape[2] if image.ndim == 3 else 1
    canvas_shape = (target_height, target_width, channels)
    canvas = np.ones(canvas_shape, dtype=image.dtype)
    if vertical == "top":
        top = 0
    elif vertical == "bottom":
        top = target_height - height
    else:
        top = (target_height - height) // 2
    if horizontal == "left":
        left = 0
    elif horizontal == "right":
        left = target_width - width
    else:
        left = (target_width - width) // 2
    canvas[top : top + height, left : left + width, ...] = image
    return canvas


def parse_svg_viewbox(path: Path) -> tuple[float, float, ET.Element]:
    root = ET.parse(path).getroot()
    viewbox = root.attrib.get("viewBox")
    if viewbox is None:
        raise ValueError(f"SVG is missing a viewBox: {path}")
    parts = [float(part) for part in viewbox.split()]
    if len(parts) != 4:
        raise ValueError(f"Unexpected SVG viewBox for {path}: {viewbox}")
    return parts[2], parts[3], root


def write_combined_png(panel_outputs: dict[str, dict[str, Path]], *, out_dir: Path, dpi: int) -> None:
    required = ("A", "B", "C", "D", "E", "F")
    missing = [panel for panel in required if panel not in panel_outputs or "png" not in panel_outputs[panel]]
    if missing:
        return

    panel_images = {
        panel: trim_white_border(plt.imread(panel_outputs[panel]["png"]))
        for panel in required
    }
    row_pairs = (("A", "B"), ("C", "D"), ("E", "F"))
    column_pairs = (("A", "C", "E"), ("B", "D", "F"))
    row_heights = [
        max(panel_images[panel].shape[0] for panel in row_pair)
        for row_pair in row_pairs
    ]
    column_widths = [
        max(panel_images[panel].shape[1] for panel in column_pair)
        for column_pair in column_pairs
    ]
    cells = {
        "A": pad_to_shape(panel_images["A"], (row_heights[0], column_widths[0]), vertical="top"),
        "B": pad_to_shape(panel_images["B"], (row_heights[0], column_widths[1]), vertical="top"),
        "C": pad_to_shape(panel_images["C"], (row_heights[1], column_widths[0]), vertical="top"),
        "D": pad_to_shape(panel_images["D"], (row_heights[1], column_widths[1]), vertical="top"),
        "E": pad_to_shape(panel_images["E"], (row_heights[2], column_widths[0]), vertical="top"),
        "F": pad_to_shape(panel_images["F"], (row_heights[2], column_widths[1]), vertical="top"),
    }

    gap_x = 90
    gap_y = 20
    channels = next(iter(cells.values())).shape[2]
    dtype = next(iter(cells.values())).dtype
    rows = []
    for left_panel, right_panel in row_pairs:
        rows.append(
            np.hstack(
                [
                    cells[left_panel],
                    np.ones(
                        (cells[left_panel].shape[0], gap_x, channels),
                        dtype=dtype,
                    ),
                    cells[right_panel],
                ]
            )
        )
    combined_parts = []
    for index, row in enumerate(rows):
        if index:
            combined_parts.append(
                np.ones((gap_y, rows[0].shape[1], channels), dtype=dtype)
            )
        combined_parts.append(row)
    combined = np.vstack(combined_parts)

    fig_width = 16.2
    fig_height = fig_width * combined.shape[0] / combined.shape[1]
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    ax.imshow(combined)
    ax.axis("off")
    fig.subplots_adjust(
        left=0,
        right=1,
        top=1,
        bottom=0,
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / "figure2_combined.png", dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def write_combined_svg(panel_outputs: dict[str, dict[str, Path]], *, out_dir: Path) -> None:
    required = ("A", "B", "C", "D", "E", "F")
    missing = [panel for panel in required if panel not in panel_outputs or "svg" not in panel_outputs[panel]]
    if missing:
        return

    ET.register_namespace("", SVG_NS)
    ET.register_namespace("xlink", XLINK_NS)
    row_pairs = (("A", "B"), ("C", "D"), ("E", "F"))
    panel_svgs = {
        panel: parse_svg_viewbox(panel_outputs[panel]["svg"])
        for panel in required
    }
    column_widths = [
        max(panel_svgs[panel][0] for panel in column_pair)
        for column_pair in (("A", "C", "E"), ("B", "D", "F"))
    ]
    row_heights = [
        max(panel_svgs[panel][1] for panel in row_pair)
        for row_pair in row_pairs
    ]
    gap_x = 18.0
    gap_y = 10.0
    combined_width = sum(column_widths) + gap_x
    combined_height = sum(row_heights) + gap_y * (len(row_heights) - 1)
    root = ET.Element(
        f"{{{SVG_NS}}}svg",
        {
            "width": f"{combined_width:g}pt",
            "height": f"{combined_height:g}pt",
            "viewBox": f"0 0 {combined_width:g} {combined_height:g}",
            "version": "1.1",
        },
    )
    ET.SubElement(
        root,
        f"{{{SVG_NS}}}rect",
        {
            "x": "0",
            "y": "0",
            "width": f"{combined_width:g}",
            "height": f"{combined_height:g}",
            "fill": "#ffffff",
        },
    )
    for row_index, (left_panel, right_panel) in enumerate(row_pairs):
        y = sum(row_heights[:row_index]) + gap_y * row_index
        for column_index, panel in enumerate((left_panel, right_panel)):
            panel_width, panel_height, panel_root = panel_svgs[panel]
            x = sum(column_widths[:column_index]) + gap_x * column_index
            nested = ET.SubElement(
                root,
                f"{{{SVG_NS}}}svg",
                {
                    "x": f"{x:g}",
                    "y": f"{y:g}",
                    "width": f"{panel_width:g}",
                    "height": f"{panel_height:g}",
                    "viewBox": panel_root.attrib["viewBox"],
                },
            )
            for child in list(panel_root):
                nested.append(child)

    out_dir.mkdir(parents=True, exist_ok=True)
    svg_path = out_dir / "figure2_combined.svg"
    ET.ElementTree(root).write(
        svg_path,
        encoding="utf-8",
        xml_declaration=True,
    )
    svg_text = "\n".join(
        line.rstrip()
        for line in svg_path.read_text(encoding="utf-8").splitlines()
    )
    svg_path.write_text(f"{svg_text}\n", encoding="utf-8")


def write_combined_figure(panel_outputs: dict[str, dict[str, Path]], *, out_dir: Path, dpi: int) -> None:
    write_combined_png(panel_outputs, out_dir=out_dir, dpi=dpi)
    write_combined_svg(panel_outputs, out_dir=out_dir)
    logger.info("Wrote combined Figure 2: %s", out_dir / "figure2_combined")


def main() -> int:
    args = parse_args()
    setup_logging(args.log_level)
    run_dir = (args.run_dir or find_latest_completed_run(args.results_dir)).resolve()
    metadata_path = args.metadata.resolve()
    out_dir = args.out_dir.resolve() if args.out_dir is not None else None
    manuscript_figures_dir = args.manuscript_out_dir.resolve()
    manuscript_tables_dir = args.manuscript_tables_dir.resolve()
    conservation_table = (
        args.conservation_table or default_conservation_table()
    ).resolve()
    fdr_threshold, effect_threshold = load_run_thresholds(run_dir)

    inputs = prepare_inputs(
        run_dir=run_dir,
        metadata_path=metadata_path,
        fdr_threshold=fdr_threshold,
        effect_threshold=effect_threshold,
    )

    panels = ("A", "B", "C", "D", "E", "F") if args.panel == "all" else (args.panel,)
    validate_panel_inputs(
        panels,
        gtf_path=args.gtf,
        conservation_table=conservation_table,
    )
    logger.info("Run directory: %s", run_dir)
    if out_dir is not None:
        logger.info("Diagnostic output directory: %s", out_dir)
    logger.info("Manuscript figure output directory: %s", manuscript_figures_dir)
    logger.info("Manuscript table output directory: %s", manuscript_tables_dir)
    logger.info("Conservation table: %s", conservation_table)
    logger.info("Datasets: %d", len(inputs.datasets))
    logger.info("Predictors: %s", ", ".join(inputs.tool_ids))
    gt_filter_text = (
        f"FDR<{inputs.fdr_threshold}, effect>{inputs.effect_threshold}"
        if inputs.fdr_threshold is not None
        else f"effect>{inputs.effect_threshold}"
    )
    logger.info("Ground-truth filters from run config: %s", gt_filter_text)
    panel_outputs = {}
    for panel in panels:
        panel_outputs[panel] = write_panel(
            panel,
            inputs=inputs,
            out_dir=out_dir,
            manuscript_figures_dir=manuscript_figures_dir,
            manuscript_tables_dir=manuscript_tables_dir,
            gtf_path=args.gtf,
            conservation_table=conservation_table,
            dpi=args.dpi,
        )
    if tuple(panels) == ("A", "B", "C", "D", "E", "F"):
        write_combined_figure(panel_outputs, out_dir=manuscript_figures_dir, dpi=args.dpi)
        write_supplementary_coverage_table(
            inputs=inputs,
            manuscript_tables_dir=manuscript_tables_dir,
        )
        write_gene_set_overlap_table(
            inputs=inputs,
            manuscript_tables_dir=manuscript_tables_dir,
        )
        write_gene_set_overlap_supplement_table(
            inputs=inputs,
            manuscript_tables_dir=manuscript_tables_dir,
        )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except FileNotFoundError as error:
        logger.error("%s", error)
        raise SystemExit(1) from None
