"""Add one OE experiment, from raw reads or a count matrix, to the benchmark tables.

Both routes reach the same DESeq2 fit and the same labels as the published experiments. The table goes to
<benchmark data>/genes/, with its rows in experiments.tsv, comparisons.tsv and the manifest. A miRNA new to the benchmark
gets its miRBase sequence added to metadata/benchmark_mirnas.tsv.

Usage:
  uv run funmirbench-add --config my_experiment.yaml
"""
from __future__ import annotations

import argparse
import json
import logging
import os
import pathlib
import shutil
import sys

import pandas as pd

from funmirbench.experiments_pipeline import load_yaml, repo_root, run_ingestion_config, write_json
from funmirbench.logger import parse_log_level, setup_logging

sys.path.insert(0, str(repo_root() / "pipelines/processing/src"))  # the labelling step is shared with the full rebuild
from build_benchmark_data import SIZE_FACTORS, build_manifest, run_one, sha256  # noqa: E402

logger = logging.getLogger(__name__)

GENES = pathlib.Path("metadata/benchmark_genes.tsv.gz")
MIRNAS = pathlib.Path("metadata/benchmark_mirnas.tsv")
MATURE_FASTA = pathlib.Path("data/experiments/raw/refs/mirbase/mature.fa")
MIN_REPLICATES = 2


def check_config(config: dict) -> None:
    if config.get("experiment_type") != "Overexpression":
        raise ValueError("Only overexpression experiments can be added: set experiment_type: Overexpression.")
    if not (config.get("metadata") or {}).get("tested_cell_line"):
        raise ValueError("Config metadata.tested_cell_line is required.")
    comparison = config.get("comparison") or {}
    for side in ("control", "treated"):
        n = len(comparison.get(f"{side}_samples") or comparison.get(f"{side}_columns") or [])
        if n < MIN_REPLICATES:
            raise ValueError(f"DESeq2 needs at least {MIN_REPLICATES} {side} samples; the config has {n}.")


def mirbase_sequence(name: str, mature_fasta: pathlib.Path) -> str:
    lines = mature_fasta.read_text(encoding="utf-8").splitlines()
    for header, sequence in zip(lines, lines[1:]):
        if header.startswith(">") and header[1:].split()[0] == name:
            return sequence.strip()
    raise ValueError(f"{name} has no sequence in {mature_fasta}")


def replace_table(table: pd.DataFrame, path: pathlib.Path) -> None:
    partial = path.with_name(f"{path.name}.part")
    table.to_csv(partial, sep="\t", index=False)
    os.replace(partial, path)


def add_mirna(name: str, mirnas_path: pathlib.Path, mature_fasta: pathlib.Path) -> bool:
    """Adds the miRNA if it is new to the benchmark; returns whether it was added."""
    mirnas = pd.read_csv(mirnas_path, sep="\t", keep_default_na=False)
    if name in set(mirnas.mirna_name):
        return False
    sequence = mirbase_sequence(name, mature_fasta)
    row = {"mirna_name": name, "mirna_sequence": sequence, "seed_2_8": sequence[1:8], "seed_family": ""}
    replace_table(pd.concat([mirnas, pd.DataFrame([row])], ignore_index=True), mirnas_path)
    logger.info("Added %s to %s", name, mirnas_path)
    return True


def check_same_size_factors(benchmark_data: pathlib.Path) -> None:
    """A table is only added to a folder whose tables were normalised the same way."""
    manifest_path = benchmark_data / "manifest.json"
    if not manifest_path.is_file() and not (benchmark_data / "comparisons.tsv").is_file():
        return  # nothing added yet, or a retryable partial first addition
    manifest = json.loads(manifest_path.read_text(encoding="utf-8")) if manifest_path.is_file() else {}
    if manifest.get("size_factors") != SIZE_FACTORS:
        raise ValueError(f"The tables in {benchmark_data} were not built with size factors from the control samples. "
                         "Get the current tables (`uv run funmirbench-download-data`) or rebuild them "
                         "(pipelines/processing/README.md) first.")


def comparison_row(config: dict) -> dict:
    dataset_id = config["dataset_id"]
    return {"experiment_id": dataset_id, "study": config.get("gse") or dataset_id, "mirna_name": config["mirna_name"],
            "perturbation": "OE", "cell_line": config["metadata"]["tested_cell_line"], "split": "test",
            "design": "~ batch + condition" if config["comparison"].get("batch_by_sample") else "~ condition",
            "dependency_group": config.get("dependency_group") or f"group_{dataset_id}", "action": "keep", "flags": "",
            "source_mode": config["source"]["mode"]}


def with_row(path: pathlib.Path, row: dict) -> pd.DataFrame:
    """The table with this row; a row left by an interrupted earlier try is dropped."""
    table = pd.read_csv(path, sep="\t", keep_default_na=False) if path.is_file() else pd.DataFrame(columns=list(row))
    table = pd.concat([table[table.experiment_id != row["experiment_id"]], pd.DataFrame([row])], ignore_index=True)
    return table.fillna("")


def add_experiment(config_path: pathlib.Path, benchmark_data: pathlib.Path, repo: pathlib.Path | None = None) -> pathlib.Path:
    repo = (repo or repo_root()).resolve()
    config = load_yaml(config_path)
    check_config(config)
    dataset_id = config["dataset_id"]
    comparisons = benchmark_data / "comparisons.tsv"
    if comparisons.is_file() and dataset_id in set(pd.read_csv(comparisons, sep="\t").experiment_id):
        raise ValueError(f"{dataset_id} is already in {benchmark_data}; choose another dataset_id.")
    check_same_size_factors(benchmark_data)
    env_bin = pathlib.Path(os.environ.get("ENV_PREFIX", repo / "pipelines/processing/work/env")) / "bin"
    os.environ["PATH"] = f"{env_bin}{os.pathsep}{os.environ['PATH']}"
    rscript = shutil.which("Rscript")
    if rscript is None:
        raise RuntimeError("Rscript was not found; run `pipelines/processing/setup.sh` first.")

    run_dir = pathlib.Path(run_ingestion_config(config_path, repo)["run_dir"])
    (benchmark_data / "genes").mkdir(parents=True, exist_ok=True)
    (benchmark_data / "logs").mkdir(exist_ok=True)
    fit = pd.DataFrame([{"experiment_id": dataset_id, "fit_rds": str(run_dir / "complete_analysis.rds")}])
    row = run_one(next(fit.itertuples(index=False)), genes=repo / GENES, out=benchmark_data, rscript=rscript)
    if row["n_eligible"] == 0:
        (benchmark_data / row["table"]).unlink()
        raise ValueError("No gene of the benchmark gene list has a usable estimate. Counts need Ensembl gene IDs without "
                         "version suffixes, such as ENSG00000141510.")

    mirnas_before = (repo / MIRNAS).read_bytes()
    added_mirna = add_mirna(config["mirna_name"], repo / MIRNAS, repo / MATURE_FASTA)
    try:
        experiments = with_row(benchmark_data / "experiments.tsv", row)
        replace_table(experiments, benchmark_data / "experiments.tsv")
        # comparisons.tsv goes last: an experiment counts as added once it is listed there, so anything before can be retried
        partial = comparisons.with_name(f"{comparisons.name}.part")
        with_row(comparisons, comparison_row(config)).to_csv(partial, sep="\t", index=False)
        manifest_path = benchmark_data / "manifest.json"
        if manifest_path.is_file():
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            manifest["tables"][dataset_id] = row["table_sha256"]
            manifest["experiments"] = len(manifest["tables"])
            manifest["comparisons_sha256"] = sha256(partial)
        else:
            manifest = build_manifest(repo / GENES, partial, dict(zip(experiments.experiment_id, experiments.table_sha256)))
        manifest_partial = manifest_path.with_name(f"{manifest_path.name}.part")
        write_json(manifest, manifest_partial)
        os.replace(manifest_partial, manifest_path)
        os.replace(partial, comparisons)
    except BaseException:
        if added_mirna:  # without its experiment the new miRNA would only join the decoys of the others
            restore = (repo / MIRNAS).with_name(f"{MIRNAS.name}.part")
            restore.write_bytes(mirnas_before)
            os.replace(restore, repo / MIRNAS)
        raise
    logger.info("%s: %d responsive, %d no response, %d unresolved at d = 0.5", dataset_id, row["d05_responsive"],
                row["d05_no_response_above_margin"], row["d05_unresolved"])
    return benchmark_data / row["table"]


def main() -> None:
    parser = argparse.ArgumentParser(description="Add one OE experiment to the benchmark tables.")
    parser.add_argument("--config", type=pathlib.Path, required=True)
    parser.add_argument("--benchmark-data", type=pathlib.Path, default=pathlib.Path("data/benchmark_data"))
    parser.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"])
    args = parser.parse_args()
    setup_logging(parse_log_level(args.log_level))
    logger.info("Wrote %s", add_experiment(args.config.expanduser().resolve(), args.benchmark_data.expanduser().resolve()))


if __name__ == "__main__":
    main()
