"""Common-prediction coverage summaries for benchmark reports."""

from __future__ import annotations

import itertools
import pathlib

import pandas as pd

from funmirbench.evaluate import REFERENCE_TOOL_IDS, SCORE_PREFIX


EXCLUDED_COMMON_SUMMARY_TOOL_IDS = set(REFERENCE_TOOL_IDS)


def _score_col(tool_id: str) -> str:
    return f"{SCORE_PREFIX}{tool_id}"


def _available_real_tools(joined, tool_ids, *, excluded_tool_ids=None):
    excluded_tool_ids = set(excluded_tool_ids or EXCLUDED_COMMON_SUMMARY_TOOL_IDS)
    available = []
    for tool_id in tool_ids:
        tool_id = str(tool_id)
        if tool_id in excluded_tool_ids:
            continue
        if _score_col(tool_id) in joined.columns:
            available.append(tool_id)
    return available


def _count_common(joined, tools):
    if not tools:
        return 0
    cols = [_score_col(tool_id) for tool_id in tools]
    missing = [col for col in cols if col not in joined.columns]
    if missing:
        return 0
    return int(joined[cols].notna().all(axis=1).sum())


def build_common_prediction_summary(
    joined,
    *,
    dataset_id,
    tool_ids,
    report_min_common_coverage=0.10,
    excluded_tool_ids=None,
):
    """Return dataset-level common prediction percentages.

    Percentages use the joined table row count as denominator.
    """
    total_rows = int(len(joined))
    real_tools = _available_real_tools(joined, tool_ids, excluded_tool_ids=excluded_tool_ids)
    rows = []

    for tool_id in real_tools:
        scored = int(joined[_score_col(tool_id)].notna().sum())
        rows.append(
            {
                "dataset_id": dataset_id,
                "summary_type": "single_predictor",
                "tools": tool_id,
                "tool_count": 1,
                "rows_total": total_rows,
                "rows_common": scored,
                "percent_common": float(scored / total_rows) if total_rows else float("nan"),
                "included_in_report_common_plots": bool(
                    total_rows and scored / total_rows >= float(report_min_common_coverage)
                ),
            }
        )

    eligible_tools = [
        row["tools"]
        for row in rows
        if row["summary_type"] == "single_predictor" and row["included_in_report_common_plots"]
    ]
    if eligible_tools:
        common = _count_common(joined, eligible_tools)
        rows.append(
            {
                "dataset_id": dataset_id,
                "summary_type": "report_common_set",
                "tools": ",".join(eligible_tools),
                "tool_count": len(eligible_tools),
                "rows_total": total_rows,
                "rows_common": common,
                "percent_common": float(common / total_rows) if total_rows else float("nan"),
                "included_in_report_common_plots": True,
            }
        )

    if len(real_tools) >= 2:
        common_all = _count_common(joined, real_tools)
        rows.append(
            {
                "dataset_id": dataset_id,
                "summary_type": "all_real_predictors_common_set",
                "tools": ",".join(real_tools),
                "tool_count": len(real_tools),
                "rows_total": total_rows,
                "rows_common": common_all,
                "percent_common": float(common_all / total_rows) if total_rows else float("nan"),
                "included_in_report_common_plots": False,
            }
        )

    for left, right in itertools.combinations(real_tools, 2):
        common = _count_common(joined, [left, right])
        rows.append(
            {
                "dataset_id": dataset_id,
                "summary_type": "pairwise_common_set",
                "tools": f"{left},{right}",
                "tool_count": 2,
                "rows_total": total_rows,
                "rows_common": common,
                "percent_common": float(common / total_rows) if total_rows else float("nan"),
                "included_in_report_common_plots": False,
            }
        )

    return pd.DataFrame(rows)


def write_common_prediction_summary(
    joined,
    reports_dir,
    *,
    dataset_id,
    tool_ids,
    report_min_common_coverage=0.10,
):
    reports_dir = pathlib.Path(reports_dir)
    reports_dir.mkdir(parents=True, exist_ok=True)
    summary = build_common_prediction_summary(
        joined,
        dataset_id=dataset_id,
        tool_ids=tool_ids,
        report_min_common_coverage=report_min_common_coverage,
    )
    path = reports_dir / f"{dataset_id}__common_prediction_summary.tsv"
    summary.to_csv(path, sep="\t", index=False)
    return path, summary


def write_combined_common_prediction_summary(dataset_summaries, out_tables_dir):
    out_tables_dir = pathlib.Path(out_tables_dir)
    out_tables_dir.mkdir(parents=True, exist_ok=True)
    path = out_tables_dir / "common_prediction_summary.tsv"
    if dataset_summaries:
        pd.concat(dataset_summaries, ignore_index=True).to_csv(path, sep="\t", index=False)
    else:
        pd.DataFrame().to_csv(path, sep="\t", index=False)
    return path
