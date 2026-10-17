"""Join experiment DE tables with prediction tool scores."""

from __future__ import annotations

import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Collection

import pandas as pd

from funmirbench import DatasetMeta
from funmirbench.benchmark_config import clean_optional_string, resolve_predictor_output_path
from funmirbench.de_table import (
    DEFAULT_FDR_NA_POLICY,
    de_covariate_columns,
    find_gene_id_column,
    read_de_table,
)
from funmirbench.gene_ids import strip_ensembl_version


PREDICTOR_CHUNK_SIZE = 1_000_000


@dataclass(frozen=True)
class LoadedPredictor:
    tool_id: str
    path: Path
    score_direction: str
    rows_read: int
    rank_map: dict[float, float]
    scores: pd.DataFrame
    has_dataset_id: bool


def _emit_log(logger, message: str) -> None:
    if logger is not None:
        logger(message)


def _elapsed(start: float) -> float:
    return time.perf_counter() - start


def _global_rank_percentile_map(scores: pd.Series) -> dict[float, float]:
    values = pd.to_numeric(scores, errors="coerce").dropna().astype(float)
    if values.empty:
        return {}
    ordered = pd.Series(values.unique()).sort_values(kind="mergesort").reset_index(drop=True)
    if len(ordered) <= 1:
        return {float(value): 1.0 for value in ordered}
    denominator = float(len(ordered) - 1)
    return {float(value): float(index / denominator) for index, value in enumerate(ordered)}


def _normalize_scores(scores: pd.Series, *, score_direction: str, tool_id: str) -> pd.Series:
    normalized = pd.to_numeric(scores, errors="coerce").astype(float)
    if score_direction == "lower_is_stronger":
        return -normalized
    if score_direction != "higher_is_stronger":
        raise ValueError(
            f"Unsupported score_direction {score_direction!r} for tool {tool_id!r}."
        )
    return normalized


def load_experiment_table(
    meta: DatasetMeta,
    *,
    protein_coding_gene_ids: Collection[str] | None = None,
    fdr_na_policy: str = DEFAULT_FDR_NA_POLICY,
    logger=None,
) -> pd.DataFrame:
    start = time.perf_counter()
    de = read_de_table(meta.full_path, fdr_na_policy=fdr_na_policy)
    gene_src = find_gene_id_column(de)
    if gene_src == "__index__":
        de = de.copy()
        de.insert(0, "gene_id", de.index.astype(str))
    else:
        de = de.rename(columns={gene_src: "gene_id"})
    de["gene_id"] = de["gene_id"].map(strip_ensembl_version)
    missing = [col for col in ("logFC", "FDR") if col not in de.columns]
    if missing:
        raise ValueError(f"{meta.full_path} missing required columns: {missing}")
    if de["gene_id"].duplicated().any():
        raise ValueError(f"Duplicate gene_id values found in {meta.full_path}")
    keep = ["gene_id", "logFC", "FDR"]
    for optional in (
        "PValue",
        "plot_FDR",
        "benchmark_FDR",
    ):
        if optional in de.columns:
            keep.append(optional)
    keep.extend(column for column in de_covariate_columns(de.columns) if column not in keep)
    out = de[keep].copy()
    if protein_coding_gene_ids is not None:
        before = len(out)
        protein_coding_gene_ids = {
            strip_ensembl_version(gene_id)
            for gene_id in protein_coding_gene_ids
        }
        out = out.loc[
            out["gene_id"].isin(protein_coding_gene_ids)
        ].copy()
        _emit_log(
            logger,
            f"    Protein-coding filter | DE rows={before:,} -> {len(out):,}",
        )
    out.insert(0, "mirna", meta.miRNA)
    out.insert(0, "dataset_id", meta.id)
    out.insert(2, "perturbation", meta.perturbation)
    _emit_log(logger, f"    Join timing | DE table | rows={len(out):,} | { _elapsed(start):.2f}s")
    return out


def _read_relevant_tool_scores(
    path: Path,
    *,
    tool_id: str,
    score_direction: str,
    dataset_id: str,
    mirna: str,
    min_score: float | None,
) -> tuple[pd.DataFrame, dict[float, float], int]:
    header = pd.read_csv(path, sep="\t", nrows=0)
    columns = [str(c).strip() for c in header.columns]
    column_lookup = {str(c).strip(): c for c in header.columns}
    required = ("Ensembl_ID", "miRNA_Name", "Score")
    missing = [col for col in required if col not in columns]
    if missing:
        raise ValueError(f"{path} missing required columns: {missing}")

    has_dataset_id = "Dataset_ID" in columns
    usecols = [column_lookup[col] for col in required]
    if has_dataset_id:
        usecols.append(column_lookup["Dataset_ID"])

    chunks = []
    score_values = []
    rows_read = 0
    reader = pd.read_csv(path, sep="\t", usecols=usecols, chunksize=PREDICTOR_CHUNK_SIZE)
    for chunk in reader:
        chunk.columns = [str(c).strip() for c in chunk.columns]
        rows_read += len(chunk)
        normalized_score = _normalize_scores(
            chunk["Score"],
            score_direction=score_direction,
            tool_id=tool_id,
        )
        score_values.append(normalized_score.dropna())

        mask = chunk["miRNA_Name"].astype(str) == str(mirna)
        if has_dataset_id:
            mask &= chunk["Dataset_ID"].astype(str) == str(dataset_id)
        if min_score is not None:
            mask &= normalized_score >= float(min_score)
        if not bool(mask.any()):
            continue
        selected = chunk.loc[mask, ["Ensembl_ID"]].copy()
        selected["Score"] = normalized_score.loc[mask].astype(float)
        chunks.append(selected)

    all_scores = pd.concat(score_values, ignore_index=True) if score_values else pd.Series(dtype=float)
    rank_map = _global_rank_percentile_map(all_scores)
    relevant = pd.concat(chunks, ignore_index=True) if chunks else pd.DataFrame(columns=["Ensembl_ID", "Score"])
    return relevant, rank_map, rows_read


def load_predictor_scores(
    tool_id: str,
    tool_meta: dict,
    root: Path,
    *,
    logger=None,
) -> LoadedPredictor:
    start = time.perf_counter()
    path = resolve_predictor_output_path(tool_meta["predictor_output_path"], root)
    score_direction = clean_optional_string(
        tool_meta.get("score_direction"),
        default="higher_is_stronger",
    )

    header = pd.read_csv(path, sep="\t", nrows=0)
    columns = [str(c).strip() for c in header.columns]
    column_lookup = {str(c).strip(): c for c in header.columns}
    required = ("Ensembl_ID", "miRNA_Name", "Score")
    missing = [col for col in required if col not in columns]
    if missing:
        raise ValueError(f"{path} missing required columns: {missing}")

    has_dataset_id = "Dataset_ID" in columns
    usecols = [column_lookup[col] for col in required]
    if has_dataset_id:
        usecols.append(column_lookup["Dataset_ID"])

    chunks = []
    score_values = []
    rows_read = 0
    reader = pd.read_csv(path, sep="\t", usecols=usecols, chunksize=PREDICTOR_CHUNK_SIZE)
    for chunk in reader:
        chunk.columns = [str(c).strip() for c in chunk.columns]
        rows_read += len(chunk)
        normalized_score = _normalize_scores(
            chunk["Score"],
            score_direction=score_direction,
            tool_id=tool_id,
        )
        score_values.append(normalized_score.dropna())

        selected = pd.DataFrame(
            {
                "Ensembl_ID": chunk["Ensembl_ID"].map(strip_ensembl_version),
                "miRNA_Name": chunk["miRNA_Name"].astype(str),
                "Score": normalized_score.astype(float),
            }
        )
        if has_dataset_id:
            selected["Dataset_ID"] = chunk["Dataset_ID"].astype(str)
        chunks.append(selected)

    all_scores = pd.concat(score_values, ignore_index=True) if score_values else pd.Series(dtype=float)
    rank_map = _global_rank_percentile_map(all_scores)
    if chunks:
        scores = pd.concat(chunks, ignore_index=True)
    else:
        columns = ["Ensembl_ID", "miRNA_Name", "Score"]
        if has_dataset_id:
            columns.append("Dataset_ID")
        scores = pd.DataFrame(columns=columns)
    _emit_log(
        logger,
        (
            f"Loaded predictor | {tool_id} | file_rows={rows_read:,} | "
            f"unique_scores={len(rank_map):,} | {_elapsed(start):.2f}s"
        ),
    )
    return LoadedPredictor(
        tool_id=tool_id,
        path=path,
        score_direction=score_direction,
        rows_read=rows_read,
        rank_map=rank_map,
        scores=scores,
        has_dataset_id=has_dataset_id,
    )


def load_predictor_score_cache(
    tool_ids,
    predictions,
    root,
    *,
    max_workers: int = 1,
    logger=None,
) -> dict[str, LoadedPredictor]:
    tool_ids = list(tool_ids)
    if max_workers <= 1 or len(tool_ids) <= 1:
        return {
            tool_id: load_predictor_scores(tool_id, predictions[tool_id], root, logger=logger)
            for tool_id in tool_ids
        }

    loaded_by_tool = {}
    worker_count = min(int(max_workers), len(tool_ids))
    executor = ThreadPoolExecutor(max_workers=worker_count)
    should_wait = True
    future_to_tool = {
        executor.submit(
            load_predictor_scores,
            tool_id,
            predictions[tool_id],
            root,
            logger=logger,
        ): tool_id
        for tool_id in tool_ids
    }
    try:
        for future in as_completed(future_to_tool):
            tool_id = future_to_tool[future]
            try:
                loaded_by_tool[tool_id] = future.result()
            except Exception:
                for pending in future_to_tool:
                    if pending is not future:
                        pending.cancel()
                should_wait = False
                executor.shutdown(wait=False, cancel_futures=True)
                raise
    finally:
        if should_wait:
            executor.shutdown(wait=True)
    return {tool_id: loaded_by_tool[tool_id] for tool_id in tool_ids}


def _select_loaded_tool_scores(
    loaded: LoadedPredictor,
    *,
    dataset_id: str,
    mirna: str,
    col_name: str,
    rank_col_name: str,
    min_score: float | None,
) -> pd.DataFrame:
    scores = loaded.scores
    mask = scores["miRNA_Name"].astype(str) == str(mirna)
    if loaded.has_dataset_id:
        mask &= scores["Dataset_ID"].astype(str) == str(dataset_id)
    if min_score is not None:
        mask &= scores["Score"] >= float(min_score)

    if not bool(mask.any()):
        return pd.DataFrame(columns=["gene_id", col_name, rank_col_name])

    df = scores.loc[mask, ["Ensembl_ID", "Score"]].copy()
    df[rank_col_name] = df["Score"].map(loaded.rank_map)
    df = df.rename(columns={"Ensembl_ID": "gene_id"})
    if df["gene_id"].duplicated().any():
        keep_idx = df.groupby("gene_id")["Score"].idxmax()
        df = df.loc[keep_idx, ["gene_id", "Score", rank_col_name]].reset_index(drop=True)
    else:
        df = df[["gene_id", "Score", rank_col_name]].copy()
    return df.rename(columns={"Score": col_name})


def load_tool_scores(
    tool_id: str,
    tool_meta: dict,
    root: Path,
    dataset_id: str,
    mirna: str,
    col_name: str,
    rank_col_name: str,
    min_score: float | None = None,
    *,
    logger=None,
) -> tuple[pd.DataFrame, Path]:
    start = time.perf_counter()
    path = resolve_predictor_output_path(tool_meta["predictor_output_path"], root)
    score_direction = clean_optional_string(
        tool_meta.get("score_direction"),
        default="higher_is_stronger",
    )

    df, rank_map, rows_read = _read_relevant_tool_scores(
        path,
        tool_id=tool_id,
        score_direction=score_direction,
        dataset_id=dataset_id,
        mirna=mirna,
        min_score=min_score,
    )
    df[rank_col_name] = df["Score"].map(rank_map)
    df["gene_id"] = df["Ensembl_ID"].map(strip_ensembl_version)
    if df["gene_id"].duplicated().any():
        keep_idx = df.groupby("gene_id")["Score"].idxmax()
        df = df.loc[keep_idx, ["gene_id", "Score", rank_col_name]].reset_index(drop=True)
    else:
        df = df[["gene_id", "Score", rank_col_name]].copy()
    out = df.rename(columns={"Score": col_name})
    _emit_log(
        logger,
        (
            f"    Join timing | {tool_id} | file_rows={rows_read:,} | "
            f"matched_rows={len(out):,} | unique_scores={len(rank_map):,} | {_elapsed(start):.2f}s"
        ),
    )
    return out, path


def build_joined(
    meta,
    tool_ids,
    predictions,
    root,
    min_score: float | None = None,
    *,
    protein_coding_gene_ids: Collection[str] | None = None,
    predictor_cache: dict[str, LoadedPredictor] | None = None,
    fdr_na_policy: str = DEFAULT_FDR_NA_POLICY,
    logger=None,
):
    total_start = time.perf_counter()
    joined = load_experiment_table(
        meta,
        protein_coding_gene_ids=protein_coding_gene_ids,
        fdr_na_policy=fdr_na_policy,
        logger=logger,
    )
    paths = {}
    for tool_id in tool_ids:
        if tool_id not in predictions:
            raise ValueError(f"Unknown tool {tool_id!r}. Known: {sorted(predictions)}")
        tool_start = time.perf_counter()
        if predictor_cache is not None and tool_id in predictor_cache:
            loaded = predictor_cache[tool_id]
            scores = _select_loaded_tool_scores(
                loaded,
                dataset_id=meta.id,
                mirna=meta.miRNA,
                col_name=f"score_{tool_id}",
                rank_col_name=f"global_rank_{tool_id}",
                min_score=min_score,
            )
            predictor_output_path = loaded.path
            _emit_log(
                logger,
                (
                    f"    Join timing | {tool_id} | cached_file_rows={loaded.rows_read:,} | "
                    f"matched_rows={len(scores):,} | unique_scores={len(loaded.rank_map):,}"
                ),
            )
        else:
            scores, predictor_output_path = load_tool_scores(
                tool_id,
                predictions[tool_id],
                root,
                meta.id,
                meta.miRNA,
                f"score_{tool_id}",
                f"global_rank_{tool_id}",
                min_score=min_score,
                logger=logger,
            )
        joined = joined.merge(scores, on="gene_id", how="left")
        paths[tool_id] = str(predictor_output_path)
        scored = int(joined[f"score_{tool_id}"].notna().sum())
        _emit_log(
            logger,
            f"    Join timing | {tool_id} merge | scored_genes={scored:,} | {_elapsed(tool_start):.2f}s",
        )
    _emit_log(
        logger,
        f"    Join timing | {meta.id} total | rows={len(joined):,} | {_elapsed(total_start):.2f}s",
    )
    return joined, paths
