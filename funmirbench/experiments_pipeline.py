"""Config-driven experiment ingestion pipeline for benchmark-ready DE tables."""

from __future__ import annotations

import argparse
import gzip
import json
import logging
import os
import pathlib
import re
import shutil
import subprocess
import sys
from datetime import datetime, timezone

import pandas as pd
import requests
import yaml


from funmirbench.logger import parse_log_level, setup_logging

GSE_URL_TEMPLATE = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse}"
DEFAULT_HSA_MATURE_MIRNAS_PATH = pathlib.Path(
    "data/experiments/raw/refs/mirbase/hsa_mature_mirnas.txt"
)
logger = logging.getLogger(__name__)


def repo_root() -> pathlib.Path:
    return pathlib.Path(__file__).resolve().parents[1]


def utc_now_stamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")


def normalize_space(value: str) -> str:
    return re.sub(r"\s+", " ", str(value).strip())


def gse_url(gse: str) -> str:
    return GSE_URL_TEMPLATE.format(gse=gse)


def ensure_clean_dir(path: pathlib.Path, *, force: bool) -> None:
    if path.exists():
        if not force:
            raise FileExistsError(f"{path} already exists; rerun with --force to replace it.")
        shutil.rmtree(path)
    path.mkdir(parents=True, exist_ok=True)


def write_json(obj, path: pathlib.Path) -> None:
    path.write_text(json.dumps(obj, indent=2) + "\n", encoding="utf-8")


def default_thread_count(*, cap: int, floor: int = 4) -> int:
    cpus = os.cpu_count() or 1
    return max(1, min(cap, max(floor, cpus // 2)))


def load_yaml(path: pathlib.Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def resolve_path(value: str, *, root: pathlib.Path, repo: pathlib.Path | None = None) -> pathlib.Path:
    path = pathlib.Path(value)
    if path.is_absolute():
        return path
    relative_to_config = (root / path).resolve()
    if repo is None:
        return relative_to_config
    relative_to_repo = (repo / path).resolve()
    if relative_to_config.exists() or not relative_to_repo.exists():
        return relative_to_config
    return relative_to_repo


def download_bytes(url: str) -> bytes:
    response = requests.get(url, timeout=120)
    response.raise_for_status()
    return response.content


def cached_download(url: str, dest: pathlib.Path, *, force: bool) -> pathlib.Path:
    if dest.exists() and not force:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_bytes(download_bytes(url))
    return dest


def read_table_auto(path: pathlib.Path) -> pd.DataFrame:
    return pd.read_csv(path, sep=None, engine="python", compression="infer")


def output_de_table_rel_path(dataset_id: str) -> pathlib.Path:
    return pathlib.Path("data") / "experiments" / "processed" / f"{dataset_id}.tsv"


def require_fields(config: dict, fields: list[str]) -> None:
    missing = [field for field in fields if not normalize_space(config.get(field, ""))]
    if missing:
        raise ValueError(f"Config is missing required top-level fields: {missing}")

def load_hsa_mature_mirna_names(repo: pathlib.Path) -> set[str]:
    path = repo / DEFAULT_HSA_MATURE_MIRNAS_PATH

    if not path.exists():
        raise ValueError(
            f"miRBase hsa mature miRNA list does not exist: {path}. "
            "Run `uv run funmirbench-experiments-download-examples --example mirbase-hsa-mature` first."
        )

    names = {
        normalize_space(line)
        for line in path.read_text(encoding="utf-8").splitlines()
        if normalize_space(line)
    }

    if not names:
        raise ValueError(f"miRBase hsa mature miRNA list is empty: {path}")

    return names


def validate_mirna_name(config: dict, *, repo: pathlib.Path) -> None:
    mirna_name = normalize_space(config.get("mirna_name", ""))

    if not mirna_name:
        raise ValueError("Config mirna_name is required.")

    mature_names = load_hsa_mature_mirna_names(repo)

    if mirna_name not in mature_names:
        query = mirna_name.lower()
        query_no_hsa = query.removeprefix("hsa-")

        suggestions = sorted(
            name for name in mature_names
            if query in name.lower()
            or query_no_hsa in name.lower()
            or name.lower() in query
        )[:10]

        message = (
            f"Config mirna_name {mirna_name!r} was not found in the miRBase "
            "Homo sapiens mature miRNA list. Use a canonical mature miRBase name "
            "with the 'hsa-' prefix and exact casing, for example 'hsa-miR-21-5p'."
        )

        if suggestions:
            message += f" Did you mean one of: {', '.join(suggestions)}?"

        raise ValueError(message)

def open_text_auto(path: pathlib.Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def load_local_source_table(
    source_cfg: dict,
    *,
    path_key: str,
    config_path: pathlib.Path,
    repo: pathlib.Path,
) -> tuple[pd.DataFrame, pathlib.Path]:
    path_value = normalize_space(source_cfg.get(path_key, ""))
    if not path_value:
        raise ValueError(f"Config source.{path_key} is required.")
    source_path = resolve_path(path_value, root=config_path.parent, repo=repo)
    if not source_path.exists():
        raise ValueError(f"Config source.{path_key} does not exist: {source_path}")
    return read_table_auto(source_path), source_path


def candidate_metadata_row(config: dict, *, de_table_rel_path: str) -> dict:
    gse = normalize_space(config.get("gse", ""))
    metadata_cfg = config.get("metadata", {})
    return {
        "id": config["dataset_id"],
        "mirna_name": config["mirna_name"],
        "mirna_sequence": metadata_cfg.get("mirna_sequence", ""),
        "article_pubmed_id": metadata_cfg.get("article_pubmed_id", ""),
        "organism": metadata_cfg.get("organism", ""),
        "tested_cell_line": metadata_cfg.get("tested_cell_line", ""),
        "treatment": metadata_cfg.get("treatment", ""),
        "tissue": metadata_cfg.get("tissue", ""),
        "method": metadata_cfg.get("method", "RNA-seq"),
        "experiment_type": config["experiment_type"],
        "gse_url": gse_url(gse) if gse else "",
        "de_table_path": de_table_rel_path,
    }


def write_candidate_metadata(config: dict, out_path: pathlib.Path, *, de_table_rel_path: str) -> None:
    pd.DataFrame([candidate_metadata_row(config, de_table_rel_path=de_table_rel_path)]).to_csv(
        out_path, sep="\t", index=False
    )


def validate_explicit_columns(df: pd.DataFrame, columns: list[str], *, key_name: str) -> list[str]:
    if not columns:
        raise ValueError(f"Config comparison.{key_name} is empty.")
    missing = [column for column in columns if column not in df.columns]
    if missing:
        raise ValueError(f"Config comparison.{key_name} references missing columns: {missing}")
    if len(set(columns)) != len(columns):
        raise ValueError(f"Config comparison.{key_name} contains duplicates: {columns}")
    return columns


def sample_entries_from_columns(columns: list[str], *, group_name: str) -> list[dict]:
    return [
        {
            "accession": column,
            "title": column,
            "group_label": group_name,
            "count_matrix_column": column,
        }
        for column in columns
    ]


def sample_column_mapping_rows(group_name: str, sample_entries: list[dict], columns: list[str]) -> list[dict]:
    rows = []
    for sample, column in zip(sample_entries, columns, strict=True):
        rows.append(
            {
                "group": group_name,
                "sample_accession": sample.get("accession", ""),
                "sample_title": sample.get("title", ""),
                "group_label": sample.get("group_label", ""),
                "resolved_count_column": column,
            }
        )
    return rows


def resolve_gene_id_column(source_cfg: dict, counts_df: pd.DataFrame) -> str:
    configured_gene_id = normalize_space(source_cfg.get("gene_id_column", ""))
    if not configured_gene_id:
        raise ValueError("Config source.gene_id_column is required for count-based modes.")
    if configured_gene_id not in counts_df.columns:
        raise ValueError(
            f"Configured source.gene_id_column {configured_gene_id!r} was not found in the count matrix."
        )
    return configured_gene_id


def require_local_binary(name: str) -> None:
    if shutil.which(name) is None:
        raise RuntimeError(
            f"Required executable {name!r} was not found on PATH. "
            "Install and activate the experiments environment from pipelines/experiments/environment.yml."
        )


def run_logged_command(
    command: list[str],
    *,
    cwd: pathlib.Path,
    stdout_path: pathlib.Path,
    stderr_path: pathlib.Path,
    error_label: str,
) -> None:
    env = os.environ.copy()
    env["LC_ALL"] = "C"
    env["LANG"] = "C"
    env["LANGUAGE"] = "C"
    env["LC_CTYPE"] = "C"
    completed = subprocess.run(
        command,
        cwd=str(cwd),
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    stdout_path.write_text(completed.stdout, encoding="utf-8")
    stderr_path.write_text(completed.stderr, encoding="utf-8")
    if completed.returncode != 0:
        raise RuntimeError(
            f"{error_label} failed with exit code {completed.returncode}. "
            f"See {stdout_path} and {stderr_path}."
        )


def run_deseq2_from_counts(
    *,
    repo: pathlib.Path,
    run_dir: pathlib.Path,
    counts_path: pathlib.Path,
    gene_id_column: str,
    control_cols: list[str],
    treated_cols: list[str],
    output_path: pathlib.Path,
) -> tuple[list[str], pathlib.Path, pathlib.Path]:
    require_local_binary("Rscript")
    logger.info("Running DESeq2 from count matrix...")
    r_script = repo / "pipelines" / "experiments" / "run_deseq2_counts.R"
    command = [
        "Rscript",
        str(r_script),
        "--counts",
        str(counts_path),
        "--gene-id-column",
        gene_id_column,
        "--control-columns",
        ",".join(control_cols),
        "--treated-columns",
        ",".join(treated_cols),
        "--output",
        str(output_path),
    ]
    stdout_path = run_dir / "deseq2.stdout.txt"
    stderr_path = run_dir / "deseq2.stderr.txt"
    run_logged_command(
        command,
        cwd=repo,
        stdout_path=stdout_path,
        stderr_path=stderr_path,
        error_label="DESeq2 run",
    )
    if not output_path.exists():
        raise ValueError(f"DESeq2 completed but did not create the expected DE table: {output_path}")
    return command, stdout_path, stderr_path


def run_de_from_counts(
    *,
    config: dict,
    config_path: pathlib.Path,
    repo: pathlib.Path,
    run_dir: pathlib.Path,
    counts_df: pd.DataFrame,
    counts_source: pathlib.Path,
    source_mode: str,
    control_cols: list[str],
    treated_cols: list[str],
    control_samples: list[dict],
    treated_samples: list[dict],
    gene_id_column_override: str | None = None,
    extra_manifest: dict | None = None,
) -> dict:
    overlap = sorted(set(control_cols) & set(treated_cols))
    if overlap:
        raise ValueError(f"Control and treated groups overlap: {overlap}")

    source_cfg = config.get("source", {})
    gene_id_col = gene_id_column_override or resolve_gene_id_column(source_cfg, counts_df)
    mapping_rows = sample_column_mapping_rows("control", control_samples, control_cols) + sample_column_mapping_rows(
        "treated", treated_samples, treated_cols
    )
    pd.DataFrame(mapping_rows).to_csv(run_dir / "sample_column_mapping.tsv", sep="\t", index=False)

    de_results_path = run_dir / "de_results.tsv"
    deseq2_command, stdout_path, stderr_path = run_deseq2_from_counts(
        repo=repo,
        run_dir=run_dir,
        counts_path=counts_source,
        gene_id_column=gene_id_col,
        control_cols=control_cols,
        treated_cols=treated_cols,
        output_path=de_results_path,
    )
    de_df = read_table_auto(de_results_path)
    required = {"gene_id", "logFC", "FDR"}
    if not required.issubset(de_df.columns):
        raise ValueError(
            f"DE output is missing required benchmark columns. Expected at least {sorted(required)}, "
            f"found {sorted(de_df.columns)}"
        )

    out_rel = output_de_table_rel_path(config["dataset_id"])
    out_path = repo / out_rel
    out_path.parent.mkdir(parents=True, exist_ok=True)
    de_df.to_csv(out_path, sep="\t", index=False)

    write_candidate_metadata(config, run_dir / "candidate_metadata.tsv", de_table_rel_path=str(out_rel))
    manifest = {
        "dataset_id": config["dataset_id"],
        "gse": normalize_space(config.get("gse", "")),
        "config_path": str(config_path),
        "source_mode": source_mode,
        "runtime": "local",
        "count_matrix_source": str(counts_source),
        "count_matrix_rows": int(counts_df.shape[0]),
        "count_matrix_columns": int(counts_df.shape[1]),
        "gene_id_column": gene_id_col,
        "control_columns": control_cols,
        "treated_columns": treated_cols,
        "sample_column_mapping": mapping_rows,
        "deseq2_command": deseq2_command,
        "deseq2_stdout": str(stdout_path),
        "deseq2_stderr": str(stderr_path),
        "de_gene_count": int(de_df.shape[0]),
        "output_de_table": str(out_path),
    }
    if extra_manifest:
        manifest.update(extra_manifest)
    write_json(manifest, run_dir / "run_manifest.json")
    pd.DataFrame(control_samples + treated_samples).to_csv(run_dir / "samples.tsv", sep="\t", index=False)
    return {"run_dir": str(run_dir), "de_table_path": str(out_path)}


def normalize_read_paths(
    value: str | list[str] | None,
    *,
    field_name: str,
    sample_id: str,
    root: pathlib.Path,
    repo: pathlib.Path,
) -> list[str]:
    if value is None:
        values = []
    elif isinstance(value, list):
        values = value
    elif isinstance(value, str):
        values = [value]
    else:
        raise ValueError(
            f"Sample {sample_id!r} has invalid {field_name}; "
            "expected a path or list of paths."
        )
    paths = [
        normalize_space(item)
        for item in values
        if item is not None and normalize_space(item)
    ]
    if not paths:
        raise ValueError(f"Sample {sample_id!r} must provide {field_name}.")

    resolved_paths = []
    for path_value in paths:
        resolved_path = resolve_path(path_value, root=root, repo=repo)
        if not resolved_path.exists():
            raise ValueError(
                f"{field_name} for sample {sample_id!r} does not exist: {resolved_path}"
            )
        resolved_paths.append(str(resolved_path))
    return resolved_paths


def normalize_sample_entry(sample: dict, *, group_name: str, root: pathlib.Path, repo: pathlib.Path) -> dict:
    sample_id = normalize_space(sample.get("sample_id") or sample.get("id") or sample.get("accession", ""))
    if not sample_id:
        raise ValueError(f"Each {group_name} sample must define sample_id.")

    reads_1 = normalize_read_paths(
        sample.get("reads_1", ""),
        field_name="reads_1",
        sample_id=sample_id,
        root=root,
        repo=repo,
    )
    reads_2_value = sample.get("reads_2", "")
    reads_2_values = reads_2_value if isinstance(reads_2_value, list) else [reads_2_value]
    has_reads_2 = any(
        item is not None and normalize_space(item) for item in reads_2_values
    )
    reads_2 = (
        normalize_read_paths(
            reads_2_value,
            field_name="reads_2",
            sample_id=sample_id,
            root=root,
            repo=repo,
        )
        if has_reads_2
        else []
    )
    if reads_2 and len(reads_1) != len(reads_2):
        raise ValueError(
            f"Sample {sample_id!r} must provide the same number of reads_1 and reads_2 files."
        )

    count_matrix_column = normalize_space(sample.get("count_matrix_column", "")) or sample_id
    return {
        "sample_id": sample_id,
        "group": group_name,
        "reads_1": reads_1,
        "reads_2": reads_2,
        "count_matrix_column": count_matrix_column,
        "accession": sample_id,
        "title": sample_id,
        "group_label": group_name,
    }


def open_reads_file(path: str):
    reads_path = pathlib.Path(path)
    if reads_path.suffix == ".gz":
        return gzip.open(reads_path, "rb")
    return reads_path.open("rb")


def combine_reads_files(paths: list[str], output_path: pathlib.Path) -> pathlib.Path:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with (
        output_path.open("wb") as output_handle,
        gzip.GzipFile(fileobj=output_handle, mode="wb", mtime=0) as compressed_handle,
    ):
        for path in paths:
            with open_reads_file(path) as input_handle:
                shutil.copyfileobj(input_handle, compressed_handle)
    return output_path


def materialize_sample_reads(sample: dict, *, run_dir: pathlib.Path) -> dict:
    reads_1_files = list(sample["reads_1"])
    reads_2_files = list(sample["reads_2"])
    materialized = dict(sample)
    materialized["source_reads_1"] = reads_1_files
    materialized["source_reads_2"] = reads_2_files

    if len(reads_1_files) == 1:
        materialized["reads_1"] = reads_1_files[0]
        materialized["reads_2"] = reads_2_files[0] if reads_2_files else ""
        return materialized

    output_dir = run_dir / "combined_reads" / sample["sample_id"]
    paired = bool(reads_2_files)
    reads_1_name = (
        f"{sample['sample_id']}_R1.fastq.gz"
        if paired
        else f"{sample['sample_id']}.fastq.gz"
    )
    materialized["reads_1"] = str(
        combine_reads_files(reads_1_files, output_dir / reads_1_name)
    )
    materialized["reads_2"] = (
        str(
            combine_reads_files(
                reads_2_files,
                output_dir / f"{sample['sample_id']}_R2.fastq.gz",
            )
        )
        if paired
        else ""
    )
    return materialized


def load_reads_samples(config: dict, *, config_path: pathlib.Path, repo: pathlib.Path) -> tuple[list[dict], list[dict]]:
    comparison_cfg = config.get("comparison", {})
    control_raw = comparison_cfg.get("control_samples", [])
    treated_raw = comparison_cfg.get("treated_samples", [])
    if not control_raw:
        raise ValueError("Config comparison.control_samples is empty.")
    if not treated_raw:
        raise ValueError("Config comparison.treated_samples is empty.")

    control_samples = [
        normalize_sample_entry(sample, group_name="control", root=config_path.parent, repo=repo)
        for sample in control_raw
    ]
    treated_samples = [
        normalize_sample_entry(sample, group_name="treated", root=config_path.parent, repo=repo)
        for sample in treated_raw
    ]
    sample_ids = [sample["sample_id"] for sample in control_samples + treated_samples]
    if len(set(sample_ids)) != len(sample_ids):
        raise ValueError(f"Reads config contains duplicate sample_id values: {sample_ids}")
    return control_samples, treated_samples


def write_reads_sample_sheet(run_dir: pathlib.Path, control_samples: list[dict], treated_samples: list[dict]) -> pathlib.Path:
    sample_sheet_path = run_dir / "reads_samples.tsv"
    pd.DataFrame(control_samples + treated_samples).to_csv(sample_sheet_path, sep="\t", index=False)
    return sample_sheet_path


def previous_run_dirs_for_dataset(
    *, repo: pathlib.Path, dataset_id: str, current_run_dir: pathlib.Path
) -> list[pathlib.Path]:
    runs_root = repo / "pipelines" / "experiments" / "runs"
    if not runs_root.exists():
        return []
    pattern = f"*_{dataset_id}"
    candidates = [path for path in runs_root.glob(pattern) if path.is_dir() and path.resolve() != current_run_dir.resolve()]
    return sorted(candidates, key=lambda path: path.stat().st_mtime, reverse=True)


def find_reusable_star_index(
    *, repo: pathlib.Path, dataset_id: str, current_run_dir: pathlib.Path
) -> pathlib.Path | None:
    for run_dir in previous_run_dirs_for_dataset(repo=repo, dataset_id=dataset_id, current_run_dir=current_run_dir):
        star_index = run_dir / "reference" / "star_index"
        if star_index.exists() and any(star_index.iterdir()):
            return star_index
    return None

def star_index_exists(path: pathlib.Path) -> bool:
    required = ["Genome", "SA", "SAindex", "genomeParameters.txt"]
    return path.exists() and path.is_dir() and all((path / name).exists() for name in required)


def shared_star_index_for_genome(genome_fasta: pathlib.Path) -> pathlib.Path:
    return genome_fasta.parent / "star_index"

def materialize_reference_file(path: pathlib.Path, *, dest_dir: pathlib.Path) -> pathlib.Path:
    if path.suffix != ".gz":
        return path
    dest_dir.mkdir(parents=True, exist_ok=True)
    out_path = dest_dir / path.stem
    if out_path.exists():
        return out_path
    with gzip.open(path, "rb") as src, out_path.open("wb") as dst:
        shutil.copyfileobj(src, dst)
    return out_path


def build_star_index(
    *,
    source_cfg: dict,
    repo: pathlib.Path,
    run_dir: pathlib.Path,
    genome_fasta: pathlib.Path,
    gtf_path: pathlib.Path,
    out_dir: pathlib.Path,
) -> tuple[pathlib.Path, list[str], pathlib.Path, pathlib.Path]:
    require_local_binary("STAR")
    threads = int(source_cfg.get("star_threads", default_thread_count(cap=32)))
    logger.info("Building STAR index with %s threads...", threads)
    out_dir.mkdir(parents=True, exist_ok=True)
    extra_args = [str(arg) for arg in source_cfg.get("star_index_extra_args", [])]
    command = [
        "STAR",
        "--runMode",
        "genomeGenerate",
        "--runThreadN",
        str(threads),
        "--genomeDir",
        str(out_dir),
        "--genomeFastaFiles",
        str(genome_fasta),
        "--sjdbGTFfile",
        str(gtf_path),
        *extra_args,
    ]
    stdout_path = run_dir / "star_index.stdout.txt"
    stderr_path = run_dir / "star_index.stderr.txt"
    run_logged_command(
        command,
        cwd=repo,
        stdout_path=stdout_path,
        stderr_path=stderr_path,
        error_label="STAR index build",
    )
    if not out_dir.exists() or not any(out_dir.iterdir()):
        raise ValueError(f"STAR index build completed but did not create files in {out_dir}")
    return out_dir, command, stdout_path, stderr_path


def prepare_reads_reference_assets(
    *,
    dataset_id: str,
    source_cfg: dict,
    config_path: pathlib.Path,
    repo: pathlib.Path,
    run_dir: pathlib.Path,
    force: bool,
) -> dict:
    genome_fasta_source = normalize_space(source_cfg.get("genome_fasta_path", ""))
    if not genome_fasta_source:
        raise ValueError("Config source.genome_fasta_path is required for reads mode.")
    gtf_source_value = normalize_space(source_cfg.get("gtf_path", ""))
    if not gtf_source_value:
        raise ValueError("Config source.gtf_path is required for reads mode.")
    genome_fasta_input = resolve_path(genome_fasta_source, root=config_path.parent, repo=repo)
    if not genome_fasta_input.exists():
        raise ValueError(f"Config source.genome_fasta_path does not exist: {genome_fasta_input}")
    gtf_input = resolve_path(gtf_source_value, root=config_path.parent, repo=repo)
    if not gtf_input.exists():
        raise ValueError(f"Config source.gtf_path does not exist: {gtf_input}")

    star_index = None
    genome_fasta_path = None
    generated_star_index = False
    reused_star_index = False
    command_paths: dict[str, str | list[str]] = {}


    materialized_gtf = materialize_reference_file(gtf_input, dest_dir=gtf_input.parent)
    genome_fasta_path = materialize_reference_file(genome_fasta_input, dest_dir=genome_fasta_input.parent)

    shared_star_index = shared_star_index_for_genome(genome_fasta_path)

    if star_index_exists(shared_star_index) and not force:
        logger.info("Reusing shared STAR index from %s...", shared_star_index)
        star_index = shared_star_index
        reused_star_index = True
    else:
        star_index, command, stdout_path, stderr_path = build_star_index(
            source_cfg=source_cfg,
            repo=repo,
            run_dir=run_dir,
            genome_fasta=genome_fasta_path,
            gtf_path=materialized_gtf,
            out_dir=shared_star_index,
        )
        generated_star_index = True
        command_paths.update(
            {
                "star_index_command": command,
                "star_index_stdout": str(stdout_path),
                "star_index_stderr": str(stderr_path),
            }
        )

    assert star_index is not None
    return {
        "star_index": star_index,
        "gtf_path": materialized_gtf,
        "genome_fasta": str(genome_fasta_path),
        "generated_star_index": generated_star_index,
        "reused_star_index": reused_star_index,
        **command_paths,
    }


def infer_library_layout(samples: list[dict]) -> str:
    paired_flags = [bool(sample["reads_2"]) for sample in samples]
    if all(paired_flags):
        return "paired"
    if not any(paired_flags):
        return "single"
    return "mixed"


def group_sample_ids_by_layout(
    samples: list[dict], sample_order: list[str]
) -> dict[str, list[str]]:
    samples_by_id = {sample["sample_id"]: sample for sample in samples}
    groups = {
        layout: [
            sample_id
            for sample_id in sample_order
            if bool(samples_by_id[sample_id]["reads_2"]) == (layout == "paired")
        ]
        for layout in ("single", "paired")
    }
    return {layout: sample_ids for layout, sample_ids in groups.items() if sample_ids}


def run_fastqc(
    *,
    stage_name: str,
    source_cfg: dict,
    repo: pathlib.Path,
    run_dir: pathlib.Path,
    sample: dict,
    reads_1: str,
    reads_2: str,
) -> tuple[dict, list[str], pathlib.Path, pathlib.Path]:
    require_local_binary("fastqc")
    threads = int(source_cfg.get("fastqc_threads", default_thread_count(cap=16)))
    logger.info(
        "Running FastQC (%s) for sample %s with %s threads...",
        stage_name,
        sample["sample_id"],
        threads,
    )
    out_dir = run_dir / "fastqc" / stage_name / sample["sample_id"]
    out_dir.mkdir(parents=True, exist_ok=True)
    inputs = [reads_1]
    if reads_2:
        inputs.append(reads_2)
    command = [
        "fastqc",
        "--threads",
        str(threads),
        "--outdir",
        str(out_dir),
        *inputs,
    ]
    stdout_path = run_dir / f"fastqc_{stage_name}_{sample['sample_id']}.stdout.txt"
    stderr_path = run_dir / f"fastqc_{stage_name}_{sample['sample_id']}.stderr.txt"
    run_logged_command(
        command,
        cwd=repo,
        stdout_path=stdout_path,
        stderr_path=stderr_path,
        error_label=f"FastQC {stage_name} for sample {sample['sample_id']}",
    )
    report_info = {
        "output_dir": str(out_dir),
        "html": sorted(str(path) for path in out_dir.glob("*_fastqc.html")),
        "zip": sorted(str(path) for path in out_dir.glob("*_fastqc.zip")),
    }
    return report_info, command, stdout_path, stderr_path


def run_fastp(
    *,
    source_cfg: dict,
    repo: pathlib.Path,
    run_dir: pathlib.Path,
    sample: dict,
) -> tuple[dict, list[str], pathlib.Path, pathlib.Path]:
    require_local_binary("fastp")
    threads = int(source_cfg.get("fastp_threads", default_thread_count(cap=16)))
    logger.info(
        "Running fastp for sample %s with %s threads...",
        sample["sample_id"],
        threads,
    )
    out_dir = run_dir / "trimmed" / sample["sample_id"]
    out_dir.mkdir(parents=True, exist_ok=True)
    paired = bool(sample["reads_2"])
    reads_1_out = out_dir / (f"{sample['sample_id']}_R1.fastq.gz" if paired else f"{sample['sample_id']}.fastq.gz")
    reads_2_out = out_dir / f"{sample['sample_id']}_R2.fastq.gz" if paired else None
    html_path = out_dir / "fastp.html"
    json_path = out_dir / "fastp.json"
    command = [
        "fastp",
        "--thread",
        str(threads),
        "-i",
        sample["reads_1"],
        "-o",
        str(reads_1_out),
        "--html",
        str(html_path),
        "--json",
        str(json_path),
        *[str(arg) for arg in source_cfg.get("fastp_extra_args", [])],
    ]
    if paired:
        assert reads_2_out is not None
        command.extend(["-I", sample["reads_2"], "-O", str(reads_2_out)])
    stdout_path = run_dir / f"fastp_{sample['sample_id']}.stdout.txt"
    stderr_path = run_dir / f"fastp_{sample['sample_id']}.stderr.txt"
    run_logged_command(
        command,
        cwd=repo,
        stdout_path=stdout_path,
        stderr_path=stderr_path,
        error_label=f"fastp trimming for sample {sample['sample_id']}",
    )
    if not reads_1_out.exists():
        raise ValueError(f"fastp completed but did not create {reads_1_out}")
    if paired and (reads_2_out is None or not reads_2_out.exists()):
        raise ValueError(f"fastp completed but did not create the expected mate file for {sample['sample_id']}")
    trimmed = dict(sample)
    trimmed["reads_1"] = str(reads_1_out)
    trimmed["reads_2"] = str(reads_2_out) if reads_2_out else ""
    return trimmed, command, stdout_path, stderr_path


def run_star_alignment(
    *,
    source_cfg: dict,
    repo: pathlib.Path,
    run_dir: pathlib.Path,
    sample: dict,
    star_index: pathlib.Path,
) -> tuple[pathlib.Path, list[str], pathlib.Path, pathlib.Path]:
    require_local_binary("STAR")
    threads = int(source_cfg.get("star_threads", default_thread_count(cap=32)))
    logger.info("Running STAR for sample %s with %s threads...", sample["sample_id"], threads)
    out_dir = run_dir / "star" / sample["sample_id"]
    out_dir.mkdir(parents=True, exist_ok=True)
    reads = [sample["reads_1"]]
    if sample["reads_2"]:
        reads.append(sample["reads_2"])
    command = [
        "STAR",
        "--runThreadN",
        str(threads),
        "--genomeDir",
        str(star_index),
        "--readFilesIn",
        *reads,
        "--outFileNamePrefix",
        f"{out_dir}/",
        "--outSAMtype",
        "BAM",
        "SortedByCoordinate",
        *[str(arg) for arg in source_cfg.get("star_extra_args", [])],
    ]
    if any(read.endswith(".gz") for read in reads):
        read_files_command = normalize_space(source_cfg.get("read_files_command", "")) or "zcat"
        command.extend(["--readFilesCommand", read_files_command])
    stdout_path = run_dir / f"star_{sample['sample_id']}.stdout.txt"
    stderr_path = run_dir / f"star_{sample['sample_id']}.stderr.txt"
    run_logged_command(
        command,
        cwd=repo,
        stdout_path=stdout_path,
        stderr_path=stderr_path,
        error_label=f"STAR alignment for sample {sample['sample_id']}",
    )
    bam_path = out_dir / "Aligned.sortedByCoord.out.bam"
    if not bam_path.exists():
        raise ValueError(f"STAR completed but did not create {bam_path}")
    return bam_path, command, stdout_path, stderr_path


def run_featurecounts(
    *,
    source_cfg: dict,
    repo: pathlib.Path,
    run_dir: pathlib.Path,
    gtf_path: pathlib.Path,
    bam_paths: dict[str, pathlib.Path],
    sample_order: list[str],
    paired_end: bool,
    output_label: str = "",
) -> tuple[pathlib.Path, list[str], pathlib.Path, pathlib.Path]:
    require_local_binary("featureCounts")
    threads = int(source_cfg.get("featurecounts_threads", default_thread_count(cap=16)))
    logger.info("Running featureCounts with %s threads...", threads)
    feature_type = normalize_space(source_cfg.get("featurecounts_feature_type", "")) or "exon"
    gene_attribute = normalize_space(source_cfg.get("featurecounts_gene_attribute", "")) or "gene_id"
    suffix = f"_{output_label}" if output_label else ""
    output_path = run_dir / f"featurecounts{suffix}_counts.tsv"
    command = [
        "featureCounts",
        "-T",
        str(threads),
        "-a",
        str(gtf_path),
        "-o",
        str(output_path),
        "-t",
        feature_type,
        "-g",
        gene_attribute,
    ]
    if paired_end:
        command.extend(["-p", "--countReadPairs"])
    command.extend(str(arg) for arg in source_cfg.get("featurecounts_extra_args", []))
    command.extend(str(bam_paths[sample_id]) for sample_id in sample_order)
    stdout_path = run_dir / f"featurecounts{suffix}.stdout.txt"
    stderr_path = run_dir / f"featurecounts{suffix}.stderr.txt"
    run_logged_command(
        command,
        cwd=repo,
        stdout_path=stdout_path,
        stderr_path=stderr_path,
        error_label="featureCounts run",
    )
    if not output_path.exists():
        raise ValueError(f"featureCounts completed but did not create {output_path}")
    return output_path, command, stdout_path, stderr_path


def build_featurecounts_matrix(
    *,
    featurecounts_path: pathlib.Path,
    bam_paths: dict[str, pathlib.Path],
    sample_order: list[str],
    out_path: pathlib.Path,
) -> pathlib.Path:
    df = pd.read_csv(featurecounts_path, sep="\t", comment="#")
    if "Geneid" not in df.columns:
        raise ValueError(f"featureCounts output is missing the Geneid column: {featurecounts_path}")

    matrix = pd.DataFrame({"gene_id": df["Geneid"].astype(str)})
    column_names = {str(column): column for column in df.columns}
    for sample_id in sample_order:
        bam_path = bam_paths[sample_id]
        candidates = [str(bam_path), bam_path.name]
        match = None
        for candidate in candidates:
            if candidate in column_names:
                match = column_names[candidate]
                break
        if match is None:
            suffix_matches = [column for column in df.columns if str(column).endswith(bam_path.name)]
            if len(suffix_matches) == 1:
                match = suffix_matches[0]
        if match is None:
            raise ValueError(
                f"Could not find the featureCounts column for sample {sample_id!r} and BAM {bam_path!s}."
            )
        matrix[sample_id] = pd.to_numeric(df[match], errors="raise")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    matrix.to_csv(out_path, sep="\t", index=False)
    return out_path


def combine_count_matrices(
    matrix_paths: list[pathlib.Path],
    sample_order: list[str],
    out_path: pathlib.Path,
) -> pathlib.Path:
    matrices = [pd.read_csv(path, sep="\t") for path in matrix_paths]
    if not matrices:
        raise ValueError("No count matrices were provided.")

    combined = matrices[0]
    for matrix in matrices[1:]:
        combined = combined.merge(matrix, on="gene_id", how="outer", validate="one_to_one")

    required_columns = ["gene_id", *sample_order]
    missing = [column for column in required_columns if column not in combined.columns]
    if missing:
        raise ValueError(f"Combined count matrix is missing sample columns: {missing}")
    if combined[required_columns].isna().any().any():
        raise ValueError("featureCounts runs produced different gene sets.")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    combined[required_columns].to_csv(out_path, sep="\t", index=False)
    return out_path


def run_count_matrix_mode(
    config: dict,
    *,
    config_path: pathlib.Path,
    repo: pathlib.Path,
    run_dir: pathlib.Path,
    force: bool,
) -> dict:
    logger.info("Mode: count_matrix (%s)", config["dataset_id"])
    source_cfg = config.get("source", {})
    logger.info("Loading count matrix...")
    counts_df, counts_source = load_local_source_table(
        source_cfg,
        path_key="count_matrix_path",
        config_path=config_path,
        repo=repo,
    )

    comparison_cfg = config.get("comparison", {})
    control_cols = validate_explicit_columns(
        counts_df,
        [normalize_space(value) for value in comparison_cfg.get("control_columns", []) if str(value).strip()],
        key_name="control_columns",
    )
    treated_cols = validate_explicit_columns(
        counts_df,
        [normalize_space(value) for value in comparison_cfg.get("treated_columns", []) if str(value).strip()],
        key_name="treated_columns",
    )
    control_samples = sample_entries_from_columns(control_cols, group_name="control")
    treated_samples = sample_entries_from_columns(treated_cols, group_name="treated")
    return run_de_from_counts(
        config=config,
        config_path=config_path,
        repo=repo,
        run_dir=run_dir,
        counts_df=counts_df,
        counts_source=counts_source,
        source_mode="count_matrix",
        control_cols=control_cols,
        treated_cols=treated_cols,
        control_samples=control_samples,
        treated_samples=treated_samples,
    )


def run_reads_mode(
    config: dict,
    *,
    config_path: pathlib.Path,
    repo: pathlib.Path,
    run_dir: pathlib.Path,
    force: bool,
) -> dict:
    logger.info("Mode: reads (%s)", config["dataset_id"])
    source_cfg = config.get("source", {})
    logger.info("Loading reads config...")
    control_samples, treated_samples = load_reads_samples(config, config_path=config_path, repo=repo)
    control_samples = [
        materialize_sample_reads(sample, run_dir=run_dir) for sample in control_samples
    ]
    treated_samples = [
        materialize_sample_reads(sample, run_dir=run_dir) for sample in treated_samples
    ]
    sample_sheet = write_reads_sample_sheet(run_dir, control_samples, treated_samples)
    sample_order = [sample["sample_id"] for sample in control_samples + treated_samples]
    library_layout = infer_library_layout(control_samples + treated_samples)

    logger.info("Preparing references...")
    reference_assets = prepare_reads_reference_assets(
        dataset_id=config["dataset_id"],
        source_cfg=source_cfg,
        config_path=config_path,
        repo=repo,
        run_dir=run_dir,
        force=force,
    )

    star_index = pathlib.Path(reference_assets["star_index"])
    gtf_path = pathlib.Path(reference_assets["gtf_path"])

    raw_fastqc_outputs: dict[str, dict] = {}
    raw_fastqc_commands: dict[str, list[str]] = {}
    raw_fastqc_stdout: dict[str, str] = {}
    raw_fastqc_stderr: dict[str, str] = {}
    trimmed_samples: list[dict] = []
    fastp_commands: dict[str, list[str]] = {}
    fastp_stdout: dict[str, str] = {}
    fastp_stderr: dict[str, str] = {}
    trimmed_fastqc_outputs: dict[str, dict] = {}
    trimmed_fastqc_commands: dict[str, list[str]] = {}
    trimmed_fastqc_stdout: dict[str, str] = {}
    trimmed_fastqc_stderr: dict[str, str] = {}

    for sample in control_samples + treated_samples:
        report_info, command, stdout_path, stderr_path = run_fastqc(
            stage_name="raw",
            source_cfg=source_cfg,
            repo=repo,
            run_dir=run_dir,
            sample=sample,
            reads_1=sample["reads_1"],
            reads_2=sample["reads_2"],
        )
        raw_fastqc_outputs[sample["sample_id"]] = report_info
        raw_fastqc_commands[sample["sample_id"]] = command
        raw_fastqc_stdout[sample["sample_id"]] = str(stdout_path)
        raw_fastqc_stderr[sample["sample_id"]] = str(stderr_path)

        current_sample, command, stdout_path, stderr_path = run_fastp(
            source_cfg=source_cfg,
            repo=repo,
            run_dir=run_dir,
            sample=sample,
        )
        fastp_commands[sample["sample_id"]] = command
        fastp_stdout[sample["sample_id"]] = str(stdout_path)
        fastp_stderr[sample["sample_id"]] = str(stderr_path)

        report_info, command, stdout_path, stderr_path = run_fastqc(
            stage_name="trimmed",
            source_cfg=source_cfg,
            repo=repo,
            run_dir=run_dir,
            sample=current_sample,
            reads_1=current_sample["reads_1"],
            reads_2=current_sample["reads_2"],
        )
        trimmed_fastqc_outputs[current_sample["sample_id"]] = report_info
        trimmed_fastqc_commands[current_sample["sample_id"]] = command
        trimmed_fastqc_stdout[current_sample["sample_id"]] = str(stdout_path)
        trimmed_fastqc_stderr[current_sample["sample_id"]] = str(stderr_path)

        trimmed_samples.append(current_sample)

    trimmed_samples_by_id = {sample["sample_id"]: sample for sample in trimmed_samples}
    control_samples = [trimmed_samples_by_id[sample["sample_id"]] for sample in control_samples]
    treated_samples = [trimmed_samples_by_id[sample["sample_id"]] for sample in treated_samples]

    bam_paths: dict[str, pathlib.Path] = {}
    star_commands: dict[str, list[str]] = {}
    star_stdout: dict[str, str] = {}
    star_stderr: dict[str, str] = {}
    for sample in control_samples + treated_samples:
        bam_path, command, stdout_path, stderr_path = run_star_alignment(
            source_cfg=source_cfg,
            repo=repo,
            run_dir=run_dir,
            sample=sample,
            star_index=star_index,
        )
        bam_paths[sample["sample_id"]] = bam_path
        star_commands[sample["sample_id"]] = command
        star_stdout[sample["sample_id"]] = str(stdout_path)
        star_stderr[sample["sample_id"]] = str(stderr_path)

    layout_groups = group_sample_ids_by_layout(
        control_samples + treated_samples, sample_order
    )

    featurecounts_runs = {}
    count_matrix_parts = []
    for layout, layout_sample_order in layout_groups.items():
        output_label = layout if library_layout == "mixed" else ""
        output, command, stdout_path, stderr_path = run_featurecounts(
            source_cfg=source_cfg,
            repo=repo,
            run_dir=run_dir,
            gtf_path=gtf_path,
            bam_paths=bam_paths,
            sample_order=layout_sample_order,
            paired_end=layout == "paired",
            output_label=output_label,
        )
        matrix_path = run_dir / (
            f"counts_matrix_{layout}.tsv" if library_layout == "mixed" else "counts_matrix.tsv"
        )
        build_featurecounts_matrix(
            featurecounts_path=output,
            bam_paths=bam_paths,
            sample_order=layout_sample_order,
            out_path=matrix_path,
        )
        count_matrix_parts.append(matrix_path)
        featurecounts_runs[layout] = {
            "command": command,
            "stdout": str(stdout_path),
            "stderr": str(stderr_path),
            "output": str(output),
        }

    counts_matrix_path = run_dir / "counts_matrix.tsv"
    if library_layout == "mixed":
        combine_count_matrices(count_matrix_parts, sample_order, counts_matrix_path)

    featurecounts_manifest = {"featurecounts_runs": featurecounts_runs}
    if len(featurecounts_runs) == 1:
        only_run = next(iter(featurecounts_runs.values()))
        featurecounts_manifest.update(
            {
                "featurecounts_command": only_run["command"],
                "featurecounts_stdout": only_run["stdout"],
                "featurecounts_stderr": only_run["stderr"],
                "featurecounts_output": only_run["output"],
            }
        )
    counts_df = read_table_auto(counts_matrix_path)
    control_cols = [sample["sample_id"] for sample in control_samples]
    treated_cols = [sample["sample_id"] for sample in treated_samples]
    trimmed_reads = {
        sample["sample_id"]: {"reads_1": sample["reads_1"], "reads_2": sample["reads_2"]}
        for sample in control_samples + treated_samples
    }
    return run_de_from_counts(
        config=config,
        config_path=config_path,
        repo=repo,
        run_dir=run_dir,
        counts_df=counts_df,
        counts_source=counts_matrix_path,
        source_mode="reads",
        control_cols=control_cols,
        treated_cols=treated_cols,
        control_samples=control_samples,
        treated_samples=treated_samples,
        gene_id_column_override="gene_id",
        extra_manifest={
            "reads_sample_sheet": str(sample_sheet),
            "library_layout": library_layout,
            "pipeline_stages": ["fastqc_raw", "fastp", "fastqc_trimmed", "star", "featurecounts", "deseq2"],
            "raw_fastqc_outputs": raw_fastqc_outputs,
            "raw_fastqc_commands": raw_fastqc_commands,
            "raw_fastqc_stdout": raw_fastqc_stdout,
            "raw_fastqc_stderr": raw_fastqc_stderr,
            "fastp_commands": fastp_commands,
            "fastp_stdout": fastp_stdout,
            "fastp_stderr": fastp_stderr,
            "trimmed_reads": trimmed_reads,
            "trimmed_fastqc_outputs": trimmed_fastqc_outputs,
            "trimmed_fastqc_commands": trimmed_fastqc_commands,
            "trimmed_fastqc_stdout": trimmed_fastqc_stdout,
            "trimmed_fastqc_stderr": trimmed_fastqc_stderr,
            "star_index": str(star_index),
            "generated_star_index": reference_assets["generated_star_index"],
            "reused_star_index": reference_assets["reused_star_index"],
            "genome_fasta": reference_assets["genome_fasta"],
            "gtf_path": str(gtf_path),
            "star_index_command": reference_assets.get("star_index_command", []),
            "star_index_stdout": reference_assets.get("star_index_stdout", ""),
            "star_index_stderr": reference_assets.get("star_index_stderr", ""),
            "star_commands": star_commands,
            "star_stdout": star_stdout,
            "star_stderr": star_stderr,
            "sample_bams": {sample_id: str(path) for sample_id, path in bam_paths.items()},
            **featurecounts_manifest,
        },
    )


def run_ingestion_config(config_path: pathlib.Path, repo: pathlib.Path | None = None, *, force: bool = False) -> dict:
    config_path = config_path.expanduser().resolve()
    repo = (repo or repo_root()).resolve()

    config = load_yaml(config_path)
    require_fields(config, ["dataset_id", "mirna_name", "experiment_type"])
    validate_mirna_name(config, repo=repo)

    dataset_id = config["dataset_id"]
    run_dir = repo / "pipelines" / "experiments" / "runs" / f"{utc_now_stamp()}_{dataset_id}"
    ensure_clean_dir(run_dir, force=force)
    shutil.copy2(config_path, run_dir / "config.yaml")
    logger.info("Run dir: %s", run_dir)

    source_cfg = config.get("source", {})
    mode = normalize_space(source_cfg.get("mode", ""))
    if not mode:
        raise ValueError("Config source.mode is required.")
    if mode == "count_matrix":
        return run_count_matrix_mode(config, config_path=config_path, repo=repo, run_dir=run_dir, force=force)
    if mode == "reads":
        return run_reads_mode(config, config_path=config_path, repo=repo, run_dir=run_dir, force=force)
    raise NotImplementedError(f"Unsupported source.mode {mode!r}.")


def parse_args() -> argparse.Namespace:
    argv = sys.argv[1:]
    if argv and argv[0] == "run":
        argv = argv[1:]

    parser = argparse.ArgumentParser(description="Run one experiment-ingestion config.")
    parser.add_argument("--config", type=pathlib.Path, required=True)
    parser.add_argument("--force", action="store_true")
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
    )
    return parser.parse_args(argv)


def main() -> None:
    args = parse_args()
    setup_logging(parse_log_level(args.log_level))
    result = run_ingestion_config(args.config, force=args.force)
    logger.info("Wrote DE table: %s", result["de_table_path"])
    logger.info("Run dir: %s", result["run_dir"])


if __name__ == "__main__":
    main()
