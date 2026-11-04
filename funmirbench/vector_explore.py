"""Download a published snapshot, explore response vectors and write standalone reports."""
from __future__ import annotations

import argparse
import hashlib
from contextlib import contextmanager
import json
import logging
import pathlib
import shutil
import time

import pandas as pd
import yaml

from funmirbench.benchmark_config import filter_df, resolve_predictor_output_path
from funmirbench.logger import parse_log_level, setup_logging
from funmirbench.predictor_store import repo_root
from funmirbench.response_vectors import create_output_dir, read_response_table
from funmirbench.vector_engine import analyse
from funmirbench.vector_report import prepare_figures, read_results, write_reports
from funmirbench.vector_sources import atomic_json, digest, fingerprint, prepare_inputs

logger = logging.getLogger(__name__)


@contextmanager
def cache_lock(path):
    """Avoid concurrent writes to resumable downloads and compact caches on POSIX."""
    path.mkdir(parents=True, exist_ok=True)
    with (path / ".explore.lock").open("a") as handle:
        try:
            import fcntl
        except ImportError:
            logger.warning("cache locking is unavailable; do not run concurrent explorers with the same cache")
            yield
        else:
            logger.info("acquiring exploratory cache lock")
            fcntl.flock(handle, fcntl.LOCK_EX)
            try:
                yield
            finally:
                fcntl.flock(handle, fcntl.LOCK_UN)


def run_exploration(config_path, **kwargs):
    config_path = pathlib.Path(config_path).expanduser().resolve()
    config = yaml.safe_load(config_path.read_text())
    if not isinstance(config, dict):
        raise ValueError("exploration config must be a YAML mapping")
    cache = resolve_predictor_output_path(config.get("cache_dir", "data/vector_explore"), config_path.parent).resolve()
    with cache_lock(cache):
        return _run_exploration(config_path, **kwargs)


def verify_table(path, expected, checks):
    state = {**fingerprint(path), "sha256": expected}
    if checks.get(str(path)) != state:
        h = hashlib.sha256()
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(4 << 20), b""):
                h.update(block)
        if h.hexdigest() != expected:
            raise ValueError(f"published table checksum mismatch: {path}; remove the isolated record cache and retry")
        checks[str(path)] = state


def _run_exploration(config_path, *, record=None, predictors=None, published_only=False, refresh=False,
                    export_vectors=False, details=None, dpi=None, timeout=120):
    start = time.perf_counter()
    config_path = pathlib.Path(config_path).expanduser().resolve()
    config = yaml.safe_load(config_path.read_text())
    if not isinstance(config, dict):
        raise ValueError("exploration config must be a YAML mapping")
    root = config_path.parent
    resolve = lambda p: resolve_predictor_output_path(p, root).resolve()
    record = str(record or config.get("zenodo_record", "23839843"))
    if not record.isdigit():
        raise ValueError("zenodo_record must be a numeric record ID")
    options = config.get("exploration") or {}
    details = options.get("detail_experiments", 3) if details is None else details
    dpi = options.get("dpi", 130) if dpi is None else dpi
    if not isinstance(details, int) or details < 0 or not isinstance(dpi, int) or not 60 <= dpi <= 300:
        raise ValueError("details must be a nonnegative integer and dpi must be between 60 and 300")
    cache = resolve(config.get("cache_dir", "data/vector_explore"))
    tool_path = resolve(config["predictions_tsv"])
    tools = pd.read_csv(tool_path, sep="\t", keep_default_na=False)
    if predictors:
        missing = set(predictors) - set(tools.tool_id)
        if missing:
            raise ValueError(f"unknown predictors: {sorted(missing)}")
        tools = tools[tools.tool_id.isin(predictors)]
    elif config.get("predictors"):
        tools = filter_df(tools, config["predictors"])
    if tools.empty or tools.tool_id.duplicated().any() or tools.tool_id.eq("").any():
        raise ValueError("select at least one uniquely named predictor")
    if not tools.score_direction.isin(["higher_is_stronger", "lower_is_stronger"]).all():
        raise ValueError("invalid predictor score direction")
    mirna_path = resolve(config["mirnas_tsv"])
    mirnas = pd.read_csv(mirna_path, sep="\t", keep_default_na=False).set_index("mirna_name")
    if mirnas.index.has_duplicates or mirnas.index.hasnans or any(not str(m).strip() for m in mirnas.index):
        raise ValueError("miRNA names must be unique and nonempty")
    sequences = mirnas.mirna_sequence.str.strip().str.upper().str.replace("T", "U", regex=False)
    if not sequences.str.fullmatch("[ACGU]{8,}", na=False).all():
        raise ValueError("miRNA sequences must contain at least eight RNA/DNA bases")
    mirnas["mirna_sequence"] = sequences
    tables_dir, sources, provenance = prepare_inputs(record, cache, tools, repo=root,
                                                    published_only=published_only, refresh=refresh, timeout=timeout)
    comparisons_path, experiments_path = tables_dir / "comparisons.tsv", tables_dir / "experiments.tsv"
    comps = pd.read_csv(comparisons_path, sep="\t", keep_default_na=False)
    if config.get("experiments"):
        comps = filter_df(comps, config["experiments"])
    skipped = comps[(comps.perturbation != "OE") | (comps.action == "set_aside")]
    comps = comps[(comps.perturbation == "OE") & (comps.action != "set_aside")]
    if comps.empty or comps.experiment_id.duplicated().any():
        raise ValueError("select at least one unique OE comparison")
    if not set(comps.mirna_name).issubset(mirnas.index):
        raise ValueError("every comparison miRNA needs a sequence in mirnas_tsv")
    final = pd.read_csv(experiments_path, sep="\t", keep_default_na=False).set_index("experiment_id")
    if final.index.has_duplicates or not set(comps.experiment_id).issubset(final.index):
        raise ValueError("missing or duplicate experiment table entries")
    check_path = tables_dir.parent / "table_checks.json"
    checks = json.loads(check_path.read_text()) if check_path.is_file() else {}
    tables, files = {}, [comparisons_path, experiments_path, mirna_path, tool_path]
    for exp in comps.experiment_id:
        path = (tables_dir / final.loc[exp, "table"]).resolve()
        if not path.is_relative_to(tables_dir.resolve()):
            raise ValueError("table path points outside the downloaded snapshot")
        if "table_sha256" in final:
            verify_table(path, final.loc[exp, "table_sha256"], checks)
        tables[exp] = read_response_table(path)
        files.append(path)
    atomic_json(check_path, checks)
    genes = next(iter(tables.values())).index
    if any(set(t.index) != set(genes) for t in tables.values()):
        raise ValueError("all selected experiments must share the same gene universe")
    tables = {exp: table.reindex(genes) for exp, table in tables.items()}
    result_cache, key, hit = analyse(tables, comps, tools, mirnas, sources, cache, input_files=files,
                                     export_vectors=export_vectors or options.get("export_vectors", False))
    out = create_output_dir(resolve(config.get("out_dir", "results/vector_explore")))
    for path in result_cache.iterdir():
        if path.is_file():
            shutil.copy2(path, out / path.name)
    if (result_cache / "vectors").is_dir():
        shutil.copytree(result_cache / "vectors", out / "vectors")
    comps.to_csv(out / "comparisons.tsv", sep="\t", index=False)
    skipped.to_csv(out / "skipped_comparisons.tsv", sep="\t", index=False)
    tools.to_csv(out / "predictors.tsv", sep="\t", index=False)
    mirnas.to_csv(out / "mirnas.tsv", sep="\t")
    shutil.copy2(config_path, out / "exploration_config.yaml")
    provenance.update({"analysis_key": key, "analysis_cache_hit": hit, "plot_dpi": dpi,
                       "detail_experiments": details, "source_scope": "published snapshot; exploratory analysis"})
    atomic_json(out / "provenance.json", provenance)
    plot_start = time.perf_counter()
    specs = prepare_figures(read_results(out), result_cache, out / "plots", dpi=dpi, details=details)
    write_reports(out, specs, provenance)
    atomic_json(out / "timings.json", {"total_seconds": time.perf_counter() - start,
                                       "plots_and_reports_seconds": time.perf_counter() - plot_start,
                                       "analysis_cache_hit": hit})
    logger.info("wrote %s (%.1f seconds)", out, time.perf_counter() - start)
    return out


def replot(directory, *, dpi=130, details=3):
    directory = pathlib.Path(directory).expanduser().resolve()
    if not 60 <= dpi <= 300 or details < 0:
        raise ValueError("dpi must be 60-300 and details must be nonnegative")
    provenance = json.loads((directory / "provenance.json").read_text())
    provenance.update(plot_dpi=dpi, detail_experiments=details)
    # No network or source TSV parsing. The saved result tables are sufficient.
    inputs = [fingerprint(directory / f"{name}.tsv") for name in
              ("specificity", "candidate_alignment", "expression_bins", "comparisons", "predictors")]
    specs = prepare_figures(read_results(directory), directory / ".report_cache" / digest(inputs), directory / "plots",
                            dpi=dpi, details=details)
    write_reports(directory, specs, provenance)
    atomic_json(directory / "provenance.json", provenance)
    return directory


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=pathlib.Path, default=repo_root() / "benchmark_vectors_zenodo.yaml")
    parser.add_argument("--record", help="pinned benchmark-table record (default: 23839843)")
    parser.add_argument("--predictors", nargs="+", help="override the selected tool IDs")
    parser.add_argument("--published-only", action="store_true", help="use Zenodo predictions even when local files exist")
    parser.add_argument("--refresh-metadata", action="store_true")
    parser.add_argument("--export-vectors", action="store_true", help="also export per-gene vectors; slower and larger")
    parser.add_argument("--details", type=int, help="number of per-experiment decoy plots (default: first 3)")
    parser.add_argument("--dpi", type=int, help="plot raster resolution, 60-300 (default: 130)")
    parser.add_argument("--timeout", type=int, default=120)
    parser.add_argument("--replot", type=pathlib.Path, help="rebuild reports from a saved run, without network or predictor parsing")
    parser.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"])
    args = parser.parse_args()
    setup_logging(parse_log_level(args.log_level))
    if args.replot:
        out = replot(args.replot, dpi=args.dpi or 130, details=3 if args.details is None else args.details)
    else:
        out = run_exploration(args.config, record=args.record, predictors=args.predictors,
                              published_only=args.published_only, refresh=args.refresh_metadata,
                              export_vectors=args.export_vectors, details=args.details, dpi=args.dpi, timeout=args.timeout)
    print(out)


if __name__ == "__main__":
    main()
