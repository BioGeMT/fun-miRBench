"""Audit generated fun-miRBench report bundles for common review issues.

This script is intentionally read-only: it inspects a generated run directory and
reports missing artifacts or reporting inconsistencies that are easy to miss in a
manual PDF/plot review.

Usage:
    python scripts/audit_report_bundle.py results/<run_id>
"""

from __future__ import annotations

import argparse
import json
import pathlib
import re
import sys
from dataclasses import dataclass


PDF_MEDIA_BOX_PATTERN = re.compile(rb"/MediaBox\s*\[\s*0\s+0\s+([0-9.]+)\s+([0-9.]+)\s*\]")
EXPECTED_METRIC_TABLES = {
    "coverage",
    "positive_coverage",
    "aps",
    "pr_auc",
    "spearman",
    "auroc",
}
EXPECTED_COMBINED_PLOTS = {
    "cross_dataset_coverage_distribution",
    "cross_dataset_positive_coverage_distribution",
    "cross_dataset_aps_distribution",
    "cross_dataset_pr_auc_distribution",
    "cross_dataset_spearman_distribution",
    "cross_dataset_auroc_distribution",
    "positive_background_local_rank_distributions",
    "positive_background_local_rank_counts",
    "positive_background_global_rank_distributions",
    "positive_background_global_rank_counts",
    "positive_recovery_fraction_by_prediction_count",
    "predictor_combination_expanded_frontier",
}
DEFAULT_MIN_COVERAGE_FOR_HEADLINE = 0.10


@dataclass(frozen=True)
class AuditIssue:
    level: str
    message: str


def _read_json(path: pathlib.Path) -> dict:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError:
        raise SystemExit(f"Missing required file: {path}") from None
    except json.JSONDecodeError as exc:
        raise SystemExit(f"Invalid JSON in {path}: {exc}") from exc


def _pdf_media_boxes(path: pathlib.Path) -> list[tuple[float, float]]:
    if not path.is_file():
        return []
    return [
        (float(width), float(height))
        for width, height in PDF_MEDIA_BOX_PATTERN.findall(path.read_bytes())
    ]


def _media_boxes_share_page_size(boxes: list[tuple[float, float]]) -> bool:
    if not boxes:
        return True
    normalized = {
        tuple(round(value, 2) for value in sorted((width, height)))
        for width, height in boxes
    }
    return len(normalized) == 1


def _load_cross_dataset_rows(summary: dict) -> list[dict[str, str]]:
    table_path = summary.get("cross_dataset_outputs", {}).get("tables", {}).get(
        "cross_dataset_predictor_summary"
    )
    if not table_path:
        return []

    path = pathlib.Path(table_path)
    if not path.is_file():
        path = pathlib.Path(summary["out_dir"]) / table_path
    if not path.is_file():
        return []

    lines = path.read_text(encoding="utf-8").strip().splitlines()
    if not lines:
        return []
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"))) for line in lines[1:]]


def _float_or_none(value: object) -> float | None:
    try:
        return float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None


def audit_run(run_dir: pathlib.Path, *, min_headline_coverage: float) -> list[AuditIssue]:
    run_dir = run_dir.expanduser().resolve()
    summary_path = run_dir / "summary.json"
    summary = _read_json(summary_path)
    issues: list[AuditIssue] = []

    report_path = pathlib.Path(summary.get("report_pdf", run_dir / "REPORT.pdf"))
    if not report_path.is_absolute():
        report_path = run_dir / report_path
    if not report_path.is_file():
        issues.append(AuditIssue("ERROR", f"Missing run PDF report: {report_path}"))
    else:
        boxes = _pdf_media_boxes(report_path)
        if len(boxes) <= 1:
            issues.append(AuditIssue("WARN", f"Run PDF has only {len(boxes)} detected page(s)."))
        if boxes and not _media_boxes_share_page_size(boxes):
            issues.append(AuditIssue("WARN", "Run PDF pages have inconsistent MediaBox sizes."))

    metric_tables = set(summary.get("metric_tables", {}))
    missing_metric_tables = sorted(EXPECTED_METRIC_TABLES - metric_tables)
    if missing_metric_tables:
        issues.append(
            AuditIssue(
                "ERROR",
                "summary.json metric_tables is missing: " + ", ".join(missing_metric_tables),
            )
        )

    combined_plots = set(summary.get("cross_dataset_outputs", {}).get("plots", {}))
    missing_combined_plots = sorted(EXPECTED_COMBINED_PLOTS - combined_plots)
    if missing_combined_plots:
        issues.append(
            AuditIssue(
                "ERROR",
                "summary.json cross_dataset_outputs.plots is missing: "
                + ", ".join(missing_combined_plots),
            )
        )

    rows = _load_cross_dataset_rows(summary)
    sparse_headline_candidates = []
    for row in rows:
        tool_id = str(row.get("tool_id", ""))
        coverage = _float_or_none(row.get("coverage_mean"))
        aps = _float_or_none(row.get("aps_mean"))
        if coverage is not None and aps is not None and coverage < min_headline_coverage:
            sparse_headline_candidates.append((tool_id, coverage, aps))
    if sparse_headline_candidates:
        formatted = ", ".join(
            f"{tool_id} coverage={coverage:.1%}, APS={aps:.3f}"
            for tool_id, coverage, aps in sparse_headline_candidates
        )
        issues.append(
            AuditIssue(
                "WARN",
                "Sparse predictors should not be used for headline rankings: "
                + formatted,
            )
        )

    dataset_reports = sorted((run_dir / "datasets").glob("*/reports/*.pdf"))
    if not dataset_reports:
        issues.append(AuditIssue("WARN", "No per-dataset predictor report PDFs found."))
    for path in dataset_reports:
        boxes = _pdf_media_boxes(path)
        if boxes and not _media_boxes_share_page_size(boxes):
            issues.append(AuditIssue("WARN", f"Inconsistent page sizes in {path}"))

    return issues


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir", type=pathlib.Path)
    parser.add_argument(
        "--min-headline-coverage",
        type=float,
        default=DEFAULT_MIN_COVERAGE_FOR_HEADLINE,
        help="Minimum mean coverage required before a predictor can be used in headline ranking warnings.",
    )
    args = parser.parse_args(argv)

    issues = audit_run(args.run_dir, min_headline_coverage=args.min_headline_coverage)
    if not issues:
        print("No report bundle audit issues found.")
        return 0

    for issue in issues:
        print(f"[{issue.level}] {issue.message}")
    return 1 if any(issue.level == "ERROR" for issue in issues) else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
