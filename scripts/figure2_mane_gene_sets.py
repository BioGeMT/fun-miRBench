#!/usr/bin/env python3
"""Generate Figure 2B/E/F gene-set panels from predictor scores over MANE genes.

Unlike the benchmark evaluation join, this reads standardized predictor tables
directly so every one of the 19,226 MANE Select genes is represented for every
tested experiment. The IGS is the intersection of genes scored by every
predictor in at least one benchmark experiment, matching the manuscript's
pooled Intersection Gene Set definition.
"""

from __future__ import annotations

import argparse
import json
import logging
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from figure2_coverage import (
    BenchmarkDataset,
    Figure2Inputs,
    evaluation_tool_colors,
    load_mane_select_gene_ids,
    load_tool_metadata,
    panel_b_gene_set_coverage,
    panel_e_utr_length,
    panel_f_utr_conservation,
    save_panel,
)
from funmirbench.logger import setup_logging


logger = logging.getLogger(__name__)
DEFAULT_EXPERIMENTS = Path("metadata/mirna_experiment_info.tsv")
DEFAULT_PREDICTORS = Path("metadata/predictions_info.tsv")
DEFAULT_GTF = Path("data/common_resources/ensembl/Homo_sapiens.GRCh38.115.gtf.gz")
DEFAULT_CONSERVATION = Path(
    "manuscript_assets/tables/figure2F_utr_conservation_mane_raw.tsv"
)
DEFAULT_OUT_ROOT = Path("results/mane_select_igs")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--experiments", type=Path, default=DEFAULT_EXPERIMENTS)
    parser.add_argument("--predictors", type=Path, default=DEFAULT_PREDICTORS)
    parser.add_argument("--gtf", type=Path, default=DEFAULT_GTF)
    parser.add_argument("--conservation-table", type=Path, default=DEFAULT_CONSERVATION)
    parser.add_argument("--out-dir", type=Path, default=None)
    parser.add_argument("--chunk-size", type=int, default=500_000)
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--log-level", default="INFO")
    return parser.parse_args()


def load_experiments(path: Path) -> pd.DataFrame:
    table = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    required = {"id", "mirna_name"}
    missing = required.difference(table.columns)
    if missing:
        raise ValueError(f"{path} is missing columns: {sorted(missing)}")
    table = table[["id", "mirna_name"]].copy()
    table["id"] = table["id"].str.strip()
    table["mirna_name"] = table["mirna_name"].str.strip()
    if table.empty or table[["id", "mirna_name"]].eq("").any().any():
        raise ValueError(f"{path} contains empty experiment IDs or miRNA names")
    if table["id"].duplicated().any():
        raise ValueError(f"{path} contains duplicate experiment IDs")
    return table


def collect_scored_genes(
    *,
    predictor_paths: dict[str, Path],
    mirna_names: set[str],
    mane_gene_ids: set[str],
    chunk_size: int,
) -> tuple[dict[str, dict[str, set[str]]], dict[str, dict[str, object]]]:
    scored: dict[str, dict[str, set[str]]] = {}
    file_audit: dict[str, dict[str, object]] = {}
    for tool_id, path in predictor_paths.items():
        if not path.exists():
            raise FileNotFoundError(f"Prediction file for {tool_id} does not exist: {path}")
        found = {name: set() for name in mirna_names}
        rows_seen = 0
        rows_for_benchmark_mirnas = 0
        logger.info("Scanning %s: %s", tool_id, path)
        for chunk in pd.read_csv(
            path,
            sep="\t",
            usecols=["Ensembl_ID", "miRNA_Name", "Score"],
            dtype={"Ensembl_ID": "string", "miRNA_Name": "string"},
            chunksize=chunk_size,
            low_memory=False,
        ):
            rows_seen += len(chunk)
            selected = chunk[chunk["miRNA_Name"].isin(mirna_names)].copy()
            if selected.empty:
                continue
            selected["Ensembl_ID"] = selected["Ensembl_ID"].str.split(".", n=1).str[0]
            selected = selected[
                selected["Ensembl_ID"].isin(mane_gene_ids)
                & selected["Score"].notna()
            ]
            rows_for_benchmark_mirnas += len(selected)
            for name, group in selected.groupby("miRNA_Name", sort=False):
                found[str(name)].update(group["Ensembl_ID"].dropna().astype(str))
        scored[tool_id] = found
        available = sorted(name for name, genes in found.items() if genes)
        missing = sorted(name for name, genes in found.items() if not genes)
        file_audit[tool_id] = {
            "path": str(path),
            "rows_scanned": rows_seen,
            "benchmark_score_rows_in_mane": rows_for_benchmark_mirnas,
            "miRNAs_with_mane_scores": available,
            "miRNAs_without_mane_scores": missing,
        }
        logger.info(
            "%s: %d MANE-scored miRNAs, %d selected score rows",
            tool_id,
            len(available),
            rows_for_benchmark_mirnas,
        )
    return scored, file_audit


def make_inputs(
    experiments: pd.DataFrame,
    *,
    tool_ids: tuple[str, ...],
    labels: dict[str, str],
    colors: dict[str, str],
    scored: dict[str, dict[str, set[str]]],
    mane_gene_ids: set[str],
) -> Figure2Inputs:
    ordered_genes = sorted(mane_gene_ids)
    datasets = []
    for experiment in experiments.itertuples(index=False):
        frame = pd.DataFrame({"gene_id": ordered_genes})
        for tool_id in tool_ids:
            genes = scored[tool_id][experiment.mirna_name]
            frame[f"score_{tool_id}"] = np.where(
                frame["gene_id"].isin(genes), 1.0, np.nan
            )
        datasets.append(BenchmarkDataset(dataset_id=experiment.id, frame=frame))
    return Figure2Inputs(
        run_dir=Path("."),
        datasets=tuple(datasets),
        tool_ids=tool_ids,
        tool_labels=labels,
        tool_colors=colors,
        fdr_threshold=None,
        effect_threshold=1.0,
    )


def main() -> int:
    args = parse_args()
    setup_logging(args.log_level)
    experiments_path = args.experiments.resolve()
    predictor_metadata = args.predictors.resolve()
    gtf_path = args.gtf.resolve()
    conservation_path = args.conservation_table.resolve()
    if not conservation_path.exists():
        raise FileNotFoundError(f"MANE conservation table not found: {conservation_path}")
    experiments = load_experiments(experiments_path)
    tool_ids, labels = load_tool_metadata(predictor_metadata)
    predictor_metadata_table = pd.read_csv(predictor_metadata, sep="\t", dtype=str).fillna("")
    paths = dict(
        zip(
            predictor_metadata_table["tool_id"].str.strip(),
            predictor_metadata_table["predictor_output_path"].str.strip().map(Path),
        )
    )
    missing_paths = set(tool_ids).difference(paths)
    if missing_paths:
        raise ValueError(f"Predictor metadata is missing paths for {sorted(missing_paths)}")
    predictor_paths = {tool_id: paths[tool_id] for tool_id in tool_ids}
    mane_gene_ids = load_mane_select_gene_ids(gtf_path)
    if len(mane_gene_ids) != 19_226:
        logger.warning(
            "Ensembl GTF contains %d protein-coding MANE Select genes (expected 19,226)",
            len(mane_gene_ids),
        )
    logger.info(
        "Experiments=%d, distinct miRNAs=%d, MANE genes=%d, predictors=%d",
        len(experiments),
        experiments["mirna_name"].nunique(),
        len(mane_gene_ids),
        len(tool_ids),
    )

    scored, file_audit = collect_scored_genes(
        predictor_paths=predictor_paths,
        mirna_names=set(experiments["mirna_name"]),
        mane_gene_ids=mane_gene_ids,
        chunk_size=args.chunk_size,
    )
    inputs = make_inputs(
        experiments,
        tool_ids=tool_ids,
        labels=labels,
        colors=evaluation_tool_colors(tool_ids),
        scored=scored,
        mane_gene_ids=mane_gene_ids,
    )

    out_dir = args.out_dir.resolve() if args.out_dir else (
        DEFAULT_OUT_ROOT
        / datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
    )
    figures_dir = out_dir / "plots"
    tables_dir = out_dir / "tables"
    figures_dir.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)

    outputs = {
        "B": panel_b_gene_set_coverage(inputs),
        "E": panel_e_utr_length(inputs, gtf_path=gtf_path),
        "F": panel_f_utr_conservation(
            inputs,
            gtf_path=gtf_path,
            conservation_table=conservation_path,
        ),
    }
    table_stems = {
        "B": "figure2B_gene_set_coverage_mane",
        "E": "figure2E_utr_length_mane",
        "F": "figure2F_utr_conservation_mane",
    }
    for panel, (summary, figure) in outputs.items():
        stem = table_stems[panel]
        summary.to_csv(tables_dir / f"{stem}.tsv", sep="\t", index=False)
        save_panel(figure, out_dir=figures_dir, stem=stem, dpi=args.dpi)

    membership_rows = []
    for experiment in experiments.itertuples(index=False):
        for tool_id in tool_ids:
            genes = scored[tool_id][experiment.mirna_name]
            membership_rows.append(
                {
                    "experiment_id": experiment.id,
                    "mirna_name": experiment.mirna_name,
                    "tool_id": tool_id,
                    "mane_gene_count": len(mane_gene_ids),
                    "scored_mane_genes": len(genes),
                    "score_coverage": len(genes) / len(mane_gene_ids),
                }
            )
    pd.DataFrame(membership_rows).to_csv(
        tables_dir / "per_experiment_score_coverage_mane.tsv", sep="\t", index=False
    )
    pooled_by_tool = {
        tool_id: set().union(*scored[tool_id].values())
        for tool_id in tool_ids
    }
    igs_genes = set.intersection(*pooled_by_tool.values()) if pooled_by_tool else set()
    pooled_membership = pd.DataFrame({"gene_id": sorted(mane_gene_ids)})
    for tool_id in tool_ids:
        pooled_membership[f"scored_{tool_id}"] = pooled_membership["gene_id"].isin(
            pooled_by_tool[tool_id]
        )
    pooled_membership["igs_gene"] = pooled_membership["gene_id"].isin(igs_genes)
    pooled_membership.to_csv(
        tables_dir / "gene_set_membership_mane.tsv", sep="\t", index=False
    )
    (out_dir / "input_audit.json").write_text(
        json.dumps(
            {
                "experiments": len(experiments),
                "distinct_mirnas": sorted(set(experiments["mirna_name"])),
                "mane_gene_count": len(mane_gene_ids),
                "predictors": file_audit,
                "definition": (
                    "A gene is scored for a predictor/miRNA when its standardized "
                    "score is non-null. Predictor gene sets pool over all listed "
                    "experiments; IGS is their intersection."
                ),
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    logger.info("Wrote MANE coverage and Figure 2B/E/F outputs to %s", out_dir)
    for panel in ("B", "E", "F"):
        summary = outputs[panel][0]
        logger.info("Panel %s summary:\n%s", panel, summary.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
