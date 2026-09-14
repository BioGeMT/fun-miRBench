"""
GEO Download Pipeline - Downloads FASTQ files from GEO/SRA and generates RNA-seq pipeline configs.

Reads experiment metadata from a TSV file (e.g. metadata/mirna_experiment_info.tsv),
downloads raw FASTQ files from GEO/SRA (GEO mode) or locates local files (local mode),
and auto-generates YAML configuration files for the RNA-seq pipeline (funmirbench-experiments).

Usage:
    python geo_download.py --tsv metadata/mirna_experiment_info.tsv

The TSV must contain the columns:
    id, mirna_name, experiment_type, gse_url, control_samples, condition_samples

Optional column:
    raw_data_dir  -- if set, local mode is used instead of GEO download.
                     control_samples/condition_samples are then interpreted as sample
                     base-names (without extension) rather than GSM accession IDs.
"""

import argparse
import csv
import json
import logging
import shutil
import subprocess
import time
import sys
import urllib.request
from collections import defaultdict
from urllib.parse import parse_qs, urlparse
from dataclasses import dataclass
from io import StringIO
from pathlib import Path

import yaml
from Bio import Entrez
import pandas as pd

# Set up logging
logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s - %(levelname)s - %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S")
logger = logging.getLogger(__name__)

# Repo root: two levels up from this file (pipelines/geo/geo_download.py)
REPO_ROOT = Path(__file__).resolve().parent.parent.parent

REQUIRED_TSV_COLUMNS = [
    "id",
    "mirna_name",
    "experiment_type",
    "gse_url",
    "control_samples",
    "condition_samples",
]

FASTQ_OUTPUT_DIR = REPO_ROOT / "data/experiments/raw"
CONFIG_OUTPUT_DIR = REPO_ROOT / "pipelines/experiments/configs"

# Default genome reference paths (Ensembl v115, downloaded by
# funmirbench-experiments-download-examples)
DEFAULT_GENOME_FASTA = str(
    REPO_ROOT / "data/experiments/raw/refs/ensembl_v115"
    / "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
)
DEFAULT_GTF = str(
    REPO_ROOT / "data/experiments/raw/refs/ensembl_v115"
    / "Homo_sapiens.GRCh38.115.gtf.gz"
)

REQUIRED_COUNT_MATRIX_COLUMNS = [
    "ENSEMBL",
    "name",
    "type",
    "chr",
    "start",
    "end",
    "str",
]


# ---------------------------------------------------------------------------
# Utility helpers
# ---------------------------------------------------------------------------

def _validate_tsv_header(fieldnames, tsv_path):
    if not fieldnames:
        raise ValueError(f"TSV has no header row: {tsv_path}")
    missing = [col for col in REQUIRED_TSV_COLUMNS if col not in fieldnames]
    if missing:
        raise ValueError(f"TSV missing required columns {missing}: {tsv_path}")


def extract_gse_accession(gse_url):
    """Extract a GSE accession from a GEO URL or accession-like string."""
    gse_url = str(gse_url).strip()
    if not gse_url:
        return ""

    parsed = urlparse(gse_url)
    accession = parse_qs(parsed.query).get("acc", [""])[0].strip()
    if accession:
        return accession

    tail = parsed.path.rstrip("/").split("/")[-1].strip()
    if tail.upper().startswith("GSE"):
        return tail

    if gse_url.upper().startswith("GSE"):
        return gse_url

    return ""


def parse_sample_list(value):
    """Parse a comma-separated sample list, handling empty/NaN values explicitly."""
    if pd.isna(value):
        return []
    return [sample.strip() for sample in str(value).split(",") if sample.strip()]


def log_step(current, total, message, prefix="STEP"):
    """Consistent numbered progress logging."""
    logger.info("[%s %d/%d] %s", prefix, current, total, message)


@dataclass
class SRRInfo:
    """Holds information about an SRA run."""
    srr: str
    gsm: str
    layout: str  # "SINGLE" or "PAIRED"

    def get_fastq_files(self):
        if self.layout == "PAIRED":
            return [f"{self.srr}_1.fastq.gz", f"{self.srr}_2.fastq.gz"]
        return [f"{self.srr}.fastq.gz"]


# ---------------------------------------------------------------------------
# GEO/SRA resolution
# ---------------------------------------------------------------------------

def _is_srr_accession(sample_id):
    """SRR, ERR, and DRR are all valid SRA run accession prefixes."""
    return sample_id.upper().startswith(("SRR", "ERR", "DRR"))


def resolve_srr_info(srr):
    """Fetch layout info for a known SRR/ERR/DRR accession directly."""
    logger.info(f"Fetching info for {srr}...")

    handle = Entrez.esearch(db="sra", term=srr)
    record = Entrez.read(handle)
    handle.close()

    sra_ids = record.get("IdList", [])
    if not sra_ids:
        raise ValueError(f"No SRA record found for {srr}")

    handle = Entrez.efetch(db="sra", id=sra_ids[0], rettype="runinfo", retmode="text")
    csv_text = handle.read()
    handle.close()

    if isinstance(csv_text, bytes):
        csv_text = csv_text.decode("utf-8")

    reader = csv.DictReader(StringIO(csv_text))
    for row in reader:
        if row.get("Run", "") == srr:
            layout = row.get("LibraryLayout", "SINGLE").upper()
            logger.info(f"  {srr} -> {layout}")
            return SRRInfo(srr=srr, gsm=srr, layout=layout)

    raise ValueError(f"Run {srr} not found in SRA runinfo")


def resolve_gsm_to_srr(gsm_id):
    logger.info(f"Resolving {gsm_id} to SRR...")

    handle = Entrez.esearch(db="gds", term=gsm_id)
    record = Entrez.read(handle)
    handle.close()

    gds_ids = record.get("IdList", [])
    if not gds_ids:
        raise ValueError(f"No GDS IDs found for {gsm_id}")

    handle = Entrez.elink(dbfrom="gds", db="sra", id=gds_ids)
    linkset = Entrez.read(handle)
    handle.close()

    if not linkset or not linkset[0].get("LinkSetDb"):
        raise ValueError(f"No SRA link found for {gsm_id}")

    sra_ids = [link["Id"] for link in linkset[0]["LinkSetDb"][0]["Link"]]

    results = []
    for sra_id in sra_ids:
        handle = Entrez.efetch(db="sra", id=sra_id, rettype="runinfo", retmode="text")
        csv_text = handle.read()
        handle.close()

        if isinstance(csv_text, bytes):
            csv_text = csv_text.decode("utf-8")

        reader = csv.DictReader(StringIO(csv_text))
        for row in reader:
            sample_name = row.get("SampleName", "")
            sample_id = row.get("Sample", "")

            if gsm_id in (sample_name, sample_id):
                srr = row.get("Run", "")
                layout = row.get("LibraryLayout", "SINGLE").upper()

                if srr:
                    results.append(SRRInfo(srr=srr, gsm=gsm_id, layout=layout))
                    logger.info(f"  {gsm_id} -> {srr} ({layout})")

    if not results:
        raise ValueError(f"No SRR runs found for {gsm_id}")

    return results


def resolve_gse_samples(sample_ids, max_retries=3, retry_delay=2.0):
    all_results = []
    failed_samples = []
    for sample_id in sample_ids:
        for attempt in range(1, max_retries + 1):
            try:
                if _is_srr_accession(sample_id):
                    all_results.append(resolve_srr_info(sample_id))
                else:
                    all_results.extend(resolve_gsm_to_srr(sample_id))
                break
            except ValueError as e:
                logger.error("Failed to resolve %s: %s", sample_id, e)
                failed_samples.append(sample_id)
                break
            except Exception as e:
                if attempt < max_retries:
                    logger.warning(
                        "Transient error resolving %s (attempt %d/%d): %s. Retrying in %.1fs.",
                        sample_id, attempt, max_retries, e, retry_delay,
                    )
                    time.sleep(retry_delay)
                else:
                    logger.error(
                        "Failed to resolve %s after %d attempts due to network/API error: %s",
                        sample_id, max_retries, e,
                    )
                    failed_samples.append(sample_id)
    return all_results, failed_samples


# ---------------------------------------------------------------------------
# Download helpers (NCBI / ENA)
# ---------------------------------------------------------------------------

def get_ena_fastq_urls(srr, layout):
    """
    Build ENA FTP URLs for FASTQ files.

    ENA URL pattern:
    ftp://ftp.sra.ebi.ac.uk/vol1/fastq/{first6}/{padding}/{srr}/{srr}_{1,2}.fastq.gz
    """
    base_url = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq"
    first6 = srr[:6]

    if len(srr) > 9:
        padding = srr[9:].zfill(3)
        path = f"{base_url}/{first6}/{padding}/{srr}"
    else:
        path = f"{base_url}/{first6}/{srr}"

    if layout == "PAIRED":
        return [f"{path}/{srr}_1.fastq.gz", f"{path}/{srr}_2.fastq.gz"]
    else:
        return [f"{path}/{srr}.fastq.gz"]


def download_from_ena(srr, output_dir, layout):
    """Download pre-compressed FASTQ files directly from ENA. Returns list of paths."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    urls = get_ena_fastq_urls(srr, layout)
    downloaded_files = []

    for url in urls:
        filename = url.split("/")[-1]
        output_path = output_dir / filename

        logger.info(f"  Downloading from ENA: {filename}...")
        try:
            urllib.request.urlretrieve(url, output_path)
            downloaded_files.append(output_path)
            logger.info(f"  Downloaded: {output_path.name} ({output_path.stat().st_size / 1024 / 1024:.1f} MB)")
        except Exception as e:
            for f in downloaded_files:
                if f.exists():
                    f.unlink()
            raise RuntimeError(f"ENA download failed for {url}: {e}")

    return downloaded_files


def download_from_ncbi(srr, output_dir, threads, layout):
    """Download and convert SRA files using prefetch + fasterq-dump. Returns compressed FASTQ paths."""
    output_dir = Path(output_dir)

    logger.info(f"  Prefetching {srr}...")
    cmd_prefetch = ["prefetch", srr, "--output-directory", str(output_dir)]
    result = subprocess.run(cmd_prefetch, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Prefetch failed for {srr}: {result.stderr}")

    logger.info(f"  Converting {srr} to FASTQ...")
    cmd_fasterq = [
        "fasterq-dump", srr,
        "--outdir", str(output_dir),
        "--temp", str(output_dir),
        "--threads", str(threads),
        "--split-files",
    ]
    result = subprocess.run(cmd_fasterq, capture_output=True, text=True, cwd=str(output_dir))
    if result.returncode != 0:
        raise RuntimeError(f"fasterq-dump failed for {srr}: {result.stderr}")

    sra_dir = output_dir / srr
    if sra_dir.is_dir():
        logger.info(f"  Cleaning up SRA directory {sra_dir}...")
        shutil.rmtree(sra_dir)

    for tmp_file in output_dir.glob("fasterq.tmp.*"):
        logger.info(f"  Removing temp file {tmp_file.name}...")
        if tmp_file.is_dir():
            shutil.rmtree(tmp_file)
        else:
            tmp_file.unlink()

    logger.info(f"  Compressing FASTQ files...")
    fastq_files = list(output_dir.glob(f"{srr}*.fastq"))
    compressed_files = []
    for fq in fastq_files:
        gz_path = Path(str(fq) + ".gz")
        if shutil.which("pigz"):
            subprocess.run(["pigz", "-p", str(threads), str(fq)], check=True)
        else:
            subprocess.run(["gzip", str(fq)], check=True)
        compressed_files.append(gz_path)

    if layout == "SINGLE":
        renamed_files = []
        for gz_file in compressed_files:
            if "_1.fastq.gz" in gz_file.name:
                new_name = gz_file.parent / gz_file.name.replace("_1.fastq.gz", ".fastq.gz")
                gz_file.rename(new_name)
                renamed_files.append(new_name)
                logger.info(f"  Renamed {gz_file.name} -> {new_name.name}")
            else:
                renamed_files.append(gz_file)
        compressed_files = renamed_files

    return compressed_files


def expected_fastq_filenames(srr, layout):
    """Return the exact compressed FASTQ filenames expected for this run."""
    if layout == "PAIRED":
        return [f"{srr}_1.fastq.gz", f"{srr}_2.fastq.gz"]
    return [f"{srr}.fastq.gz"]


def download_srr(srr, output_dir, threads, layout="PAIRED"):
    """Download FASTQ files for a single SRR accession (NCBI with ENA fallback)."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    expected_names = expected_fastq_filenames(srr, layout)
    expected_paths = [output_dir / name for name in expected_names]
    existing_gz = sorted(output_dir.glob(f"{srr}*.fastq.gz"))

    if all(path.exists() for path in expected_paths):
        logger.info(f"Skipping {srr} - already downloaded: {expected_names}")
        return expected_paths

    if existing_gz:
        logger.warning(
            "Partial FASTQ set for %s in %s; removing incomplete files before retry: %s",
            srr, output_dir, [path.name for path in existing_gz],
        )
        for path in existing_gz:
            path.unlink()

    logger.info(f"Downloading {srr} to {output_dir}")

    try:
        compressed_files = download_from_ncbi(srr, output_dir, threads, layout)
    except Exception as ncbi_error:
        logger.warning(f"  NCBI download failed: {ncbi_error}")
        logger.info(f"  Falling back to ENA...")
        try:
            compressed_files = download_from_ena(srr, output_dir, layout)
        except Exception as ena_error:
            raise RuntimeError(
                f"Both NCBI and ENA downloads failed for {srr}. "
                f"NCBI: {ncbi_error}, ENA: {ena_error}"
            )

    logger.info(f"  Done: {[f.name for f in compressed_files]}")
    return compressed_files


# ---------------------------------------------------------------------------
# TSV parsing
# ---------------------------------------------------------------------------

def parse_experiment_metadata(tsv_path):
    """
    Parse the experiments TSV and return a list of experiment dicts.

    Rows where control_samples or condition_samples are empty are skipped
    with a warning (they have not been filled in yet).
    """
    df = pd.read_csv(tsv_path, sep="\t", dtype=str, keep_default_na=False)
    _validate_tsv_header(list(df.columns), tsv_path)

    experiments = []
    for row_num, row in enumerate(df.to_dict(orient="records"), start=2):
        def get(col):
            return (row.get(col) or "").strip()

        control_samples = parse_sample_list(row.get("control_samples", ""))
        condition_samples = parse_sample_list(row.get("condition_samples", ""))

        # Skip rows not yet filled in — do this before any other validation
        # so that partially-filled placeholder rows don't cause errors.
        if not control_samples or not condition_samples:
            logger.warning(
                "Row %d: skipping — control_samples or condition_samples is empty.",
                row_num,
            )
            continue

        dataset_id = get("id")
        mirna = get("mirna_name")
        experiment_type = get("experiment_type")
        gse_url_val = get("gse_url")
        gse = extract_gse_accession(gse_url_val)
        raw_data_dir = get("raw_data_dir")
        count_matrix_path = get("count_matrix_path")
        gene_id_column = get("gene_id_column")

        # Hard errors for fields required on processable rows
        if not dataset_id:
            raise ValueError(f"Row {row_num} missing id in {tsv_path}")
        if not mirna:
            raise ValueError(f"Row {row_num} missing mirna_name in {tsv_path}")
        if not experiment_type:
            raise ValueError(f"Row {row_num} missing experiment_type in {tsv_path}")
        if not gse_url_val:
            raise ValueError(f"Row {row_num} missing gse_url in {tsv_path}")
        if not gse:
            raise ValueError(
                f"Row {row_num} has invalid gse_url with no GEO accession: {gse_url_val}"
            )

        # count_matrix mode requires gene_id_column
        if count_matrix_path and not gene_id_column:
            raise ValueError(
                f"Row {row_num} ({dataset_id}): count_matrix_path is set but gene_id_column is missing."
            )

        experiments.append({
            "id": dataset_id,
            "gse": gse,
            "mirna": mirna,
            "experiment_type": experiment_type,
            "control_samples": control_samples,
            "condition_samples": condition_samples,
            "raw_data_dir": raw_data_dir,
            "count_matrix_path": count_matrix_path,
            "gene_id_column": gene_id_column,
            # fields for the YAML metadata section
            "organism": get("organism"),
            "tested_cell_line": get("tested_cell_line"),
            "treatment": get("treatment"),
            "tissue": get("tissue"),
            "article_pubmed_id": get("article_pubmed_id"),
        })

    return experiments


# ---------------------------------------------------------------------------
# Local mode
# ---------------------------------------------------------------------------

def resolve_local_samples(sample_names, raw_data_dir):
    """
    Locate FASTQ files for each sample name inside raw_data_dir.

    Paired-end is auto-detected: looks for {name}_1.fastq.gz + {name}_2.fastq.gz.
    Single-end fallback: {name}.fastq.gz.

    Returns a list of dicts: {sample_id, reads_1, reads_2 (empty string if SE)}.
    """
    raw_data_dir = Path(raw_data_dir)
    if not raw_data_dir.is_dir():
        raise FileNotFoundError(f"raw_data_dir does not exist: {raw_data_dir}")

    entries = []
    for name in sample_names:
        reads_1_pe = raw_data_dir / f"{name}_1.fastq.gz"
        reads_2_pe = raw_data_dir / f"{name}_2.fastq.gz"
        reads_1_se = raw_data_dir / f"{name}.fastq.gz"

        if reads_1_pe.exists() and reads_2_pe.exists():
            entries.append({
                "sample_id": name,
                "reads_1": str(reads_1_pe),
                "reads_2": str(reads_2_pe),
            })
            logger.info("  %s -> paired-end (%s, %s)", name, reads_1_pe.name, reads_2_pe.name)
        elif reads_1_se.exists():
            entries.append({
                "sample_id": name,
                "reads_1": str(reads_1_se),
                "reads_2": "",
            })
            logger.info("  %s -> single-end (%s)", name, reads_1_se.name)
        else:
            raise FileNotFoundError(
                f"No FASTQ files found for sample {name!r} in {raw_data_dir}. "
                f"Expected '{name}.fastq.gz' (SE) or "
                f"'{name}_1.fastq.gz' + '{name}_2.fastq.gz' (PE)."
            )
    return entries


def process_local_experiment(experiment):
    """Local mode: resolve FASTQ files and return (control_entries, treated_entries)."""
    raw_data_dir = experiment["raw_data_dir"]
    logger.info("Local mode for %s — raw_data_dir: %s", experiment["id"], raw_data_dir)
    control_entries = resolve_local_samples(experiment["control_samples"], raw_data_dir)
    treated_entries = resolve_local_samples(experiment["condition_samples"], raw_data_dir)
    return control_entries, treated_entries


# ---------------------------------------------------------------------------
# GEO mode: build sample entries from SRR info
# ---------------------------------------------------------------------------

def build_geo_sample_entries(srr_infos, control_samples, condition_samples, output_dir):
    """
    Build YAML sample entry dicts from resolved SRR info.

    Uses each GSM accession as one biological sample. Multiple SRR runs for the
    same GSM are kept as multiple input files for that sample.

    Returns (control_entries, treated_entries).
    """
    gsm_to_srrs = defaultdict(list)
    for srr_info in srr_infos:
        gsm_to_srrs[srr_info.gsm].append(srr_info)

    control_entries = []
    treated_entries = []

    for gsm, srrs in gsm_to_srrs.items():
        sorted_srrs = sorted(srrs, key=lambda item: item.srr)
        layouts = {srr_info.layout for srr_info in sorted_srrs}
        if len(layouts) != 1:
            raise ValueError(f"GSM {gsm} has mixed SRR library layouts: {sorted(layouts)}")

        paired = next(iter(layouts)) == "PAIRED"
        reads_1 = [
            str(
                output_dir
                / (
                    f"{srr_info.srr}_1.fastq.gz"
                    if paired
                    else f"{srr_info.srr}.fastq.gz"
                )
            )
            for srr_info in sorted_srrs
        ]
        entry = {
            "sample_id": gsm,
            "reads_1": reads_1[0] if len(reads_1) == 1 else reads_1,
            "reads_2": "",
        }
        if paired:
            reads_2 = [
                str(output_dir / f"{srr_info.srr}_2.fastq.gz")
                for srr_info in sorted_srrs
            ]
            entry["reads_2"] = reads_2[0] if len(reads_2) == 1 else reads_2

        if gsm in control_samples:
            control_entries.append(entry)
        elif gsm in condition_samples:
            treated_entries.append(entry)
        else:
            raise ValueError(f"Resolved GSM {gsm} is not assigned to either comparison group.")

    return control_entries, treated_entries


def download_experiment(experiment, output_dir, threads):
    """
    GEO mode: resolve GSM accessions to SRR runs, download FASTQs,
    write manifest.json, and return (control_entries, treated_entries).
    """
    gse = experiment["gse"]
    control_samples = experiment["control_samples"]
    condition_samples = experiment["condition_samples"]

    gse_output_dir = Path(output_dir) / gse
    gse_output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Downloading experiment %s", gse)
    logger.info("  Control samples: %s", control_samples)
    logger.info("  Condition samples: %s", condition_samples)

    all_samples = control_samples + condition_samples
    step_total = 3

    log_step(1, step_total, f"{gse}: resolving GSM accessions to SRR runs", prefix="EXP")
    srr_infos, failed_resolutions = resolve_gse_samples(all_samples)
    logger.info("[%s] Resolved %d SRR run(s)", gse, len(srr_infos))
    if failed_resolutions:
        raise RuntimeError(
            f"{gse}: {len(failed_resolutions)} sample(s) could not be resolved to SRR: "
            f"{failed_resolutions}. Aborting to avoid an incomplete experiment."
        )

    log_step(2, step_total, f"{gse}: downloading SRR FASTQ files", prefix="EXP")
    failed_srrs = []
    for idx, srr_info in enumerate(srr_infos, start=1):
        try:
            logger.info(
                "[RUN %d/%d] %s -> %s (%s)",
                idx, len(srr_infos), srr_info.gsm, srr_info.srr, srr_info.layout,
            )
            download_srr(srr_info.srr, gse_output_dir, threads=threads, layout=srr_info.layout)
        except Exception as e:
            logger.error("Failed to download %s: %s", srr_info.srr, e)
            failed_srrs.append(srr_info.srr)

    if failed_srrs:
        raise RuntimeError(
            f"{gse}: {len(failed_srrs)} SRR download(s) failed: {failed_srrs}. "
            "Aborting manifest and config generation to avoid referencing missing files."
        )

    log_step(3, step_total, f"{gse}: writing manifest", prefix="EXP")
    control_entries, treated_entries = build_geo_sample_entries(
        srr_infos, control_samples, condition_samples, gse_output_dir
    )
    manifest = {
        "gse": gse,
        "samples": {
            **{e["sample_id"]: {"group": "control", "reads_1": e["reads_1"], "reads_2": e.get("reads_2", "")} for e in control_entries},
            **{e["sample_id"]: {"group": "condition", "reads_1": e["reads_1"], "reads_2": e.get("reads_2", "")} for e in treated_entries},
        },
    }
    manifest_path = gse_output_dir / "manifest.json"
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)
    logger.info("Manifest saved to %s", manifest_path)

    return control_entries, treated_entries


# ---------------------------------------------------------------------------
# Count matrix mode
# ---------------------------------------------------------------------------

def _read_count_matrix_header(path):
    """Read only the count matrix header without loading the full matrix."""
    try:
        return list(pd.read_csv(path, nrows=0).columns)
    except Exception as e:
        raise ValueError(f"Could not read count matrix header from {path}: {e}") from e


def _validate_count_matrix_header(path, experiment):
    """Validate required annotation, gene ID, and sample columns."""
    header = _read_count_matrix_header(path)
    header_set = set(header)

    required_columns = set(REQUIRED_COUNT_MATRIX_COLUMNS)
    required_columns.add(experiment["gene_id_column"])
    required_columns.update(experiment["control_samples"])
    required_columns.update(experiment["condition_samples"])

    missing = sorted(col for col in required_columns if col not in header_set)
    if missing:
        raise ValueError(
            f"Count matrix for {experiment['id']} is missing required column(s): "
            f"{missing}. Matrix path: {path}"
        )


def process_count_matrix_experiment(experiment):
    """Count matrix mode: verify the matrix file exists and has expected columns."""
    path = Path(experiment["count_matrix_path"])
    if not path.is_absolute():
        path = REPO_ROOT / path
    if not path.exists():
        raise FileNotFoundError(
            f"count_matrix_path does not exist: {path}"
        )

    _validate_count_matrix_header(path, experiment)
    logger.info("Count matrix found and validated: %s", path)


def generate_count_matrix_yaml_config(experiment, config_output_dir):
    """Write a count-matrix YAML config for funmirbench-experiments."""
    config_output_dir = Path(config_output_dir)
    config_output_dir.mkdir(parents=True, exist_ok=True)

    config = {
        "dataset_id": experiment["id"],
        "mirna_name": experiment["mirna"],
        "experiment_type": experiment["experiment_type"],
        "gse": experiment["gse"],
        "source": {
            "mode": "count_matrix",
            "count_matrix_path": experiment["count_matrix_path"],
            "gene_id_column": experiment["gene_id_column"],
        },
        "comparison": {
            "control_columns": experiment["control_samples"],
            "treated_columns": experiment["condition_samples"],
        },
        "metadata": {
            "organism": experiment.get("organism", ""),
            "tested_cell_line": experiment.get("tested_cell_line", ""),
            "treatment": experiment.get("treatment", ""),
            "tissue": experiment.get("tissue", ""),
            "method": "RNA-seq",
            "article_pubmed_id": experiment.get("article_pubmed_id", ""),
        },
    }

    config_path = config_output_dir / f"{experiment['id']}.yaml"
    yaml_str = yaml.dump(config, default_flow_style=False, allow_unicode=True, sort_keys=False)
    for key in ("source:", "comparison:", "metadata:"):
        yaml_str = yaml_str.replace(f"\n{key}", f"\n\n{key}")
    with open(config_path, "w", encoding="utf-8") as f:
        f.write(yaml_str)
    logger.info("YAML config written: %s", config_path)
    return config_path


# ---------------------------------------------------------------------------
# YAML config generation (reads mode)
# ---------------------------------------------------------------------------

def generate_yaml_config(experiment, control_entries, treated_entries, config_output_dir):
    """
    Write a YAML config for funmirbench-experiments (the RNA-seq pipeline).

    Output path: {config_output_dir}/{dataset_id}.yaml
    """
    config_output_dir = Path(config_output_dir)
    config_output_dir.mkdir(parents=True, exist_ok=True)

    def sample_block(entry):
        block = {"sample_id": entry["sample_id"], "reads_1": entry["reads_1"]}
        if entry.get("reads_2"):
            block["reads_2"] = entry["reads_2"]
        return block

    config = {
        "dataset_id": experiment["id"],
        "mirna_name": experiment["mirna"],
        "experiment_type": experiment["experiment_type"],
        "gse": experiment["gse"],
        "source": {
            "mode": "reads",
            "genome_fasta_path": DEFAULT_GENOME_FASTA,
            "gtf_path": DEFAULT_GTF,
            "fastqc_threads": 8,
            "fastp_threads": 8,
            "star_threads": 16,
            "featurecounts_threads": 8,
        },
        "comparison": {
            "control_samples": [sample_block(e) for e in control_entries],
            "treated_samples": [sample_block(e) for e in treated_entries],
        },
        "metadata": {
            "organism": experiment.get("organism", ""),
            "tested_cell_line": experiment.get("tested_cell_line", ""),
            "treatment": experiment.get("treatment", ""),
            "tissue": experiment.get("tissue", ""),
            "method": "RNA-seq",
            "article_pubmed_id": experiment.get("article_pubmed_id", ""),
        },
    }

    config_path = config_output_dir / f"{experiment['id']}.yaml"
    yaml_str = yaml.dump(config, default_flow_style=False, allow_unicode=True, sort_keys=False)
    for key in ("source:", "comparison:", "metadata:"):
        yaml_str = yaml_str.replace(f"\n{key}", f"\n\n{key}")
    with open(config_path, "w", encoding="utf-8") as f:
        f.write(yaml_str)
    logger.info("YAML config written: %s", config_path)
    return config_path


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Download FASTQ files from GEO/SRA (or locate local files) "
            "and generate RNA-seq pipeline YAML configs."
        )
    )
    parser.add_argument(
        "--tsv", required=True,
        help="Path to experiments TSV (e.g. metadata/mirna_experiment_info.tsv)",
    )
    parser.add_argument(
        "--threads", "-t", type=int, default=4,
        help="Threads for fasterq-dump (default: 4)",
    )
    parser.add_argument(
        "--entrez-email", required=True,
        help="Email address for NCBI Entrez API (required by NCBI)",
    )

    args = parser.parse_args()
    Entrez.email = args.entrez_email

    total_steps = 4
    log_step(1, total_steps, f"loading experiments metadata from {args.tsv}")
    try:
        experiments = parse_experiment_metadata(args.tsv)
    except (ValueError, FileNotFoundError) as e:
        logger.error("Failed to load TSV: %s", e)
        return 1
    log_step(2, total_steps, f"parsed {len(experiments)} processable experiment(s)")

    FASTQ_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    log_step(3, total_steps, f"output directory ready: {FASTQ_OUTPUT_DIR}")
    log_step(4, total_steps, "starting experiment processing")

    total_experiments = len(experiments)
    failed_experiments = []

    for exp_idx, exp in enumerate(experiments, start=1):
        exp_label = f"{exp['gse']} | {exp['mirna']} | {exp['experiment_type']}"
        try:
            logger.info(
                "[EXPERIMENT %d/%d] %s",
                exp_idx, total_experiments, exp_label,
            )
            if exp.get("count_matrix_path", "").strip():
                process_count_matrix_experiment(exp)
                generate_count_matrix_yaml_config(exp, CONFIG_OUTPUT_DIR)
            elif exp.get("raw_data_dir", "").strip():
                control_entries, treated_entries = process_local_experiment(exp)
                generate_yaml_config(exp, control_entries, treated_entries, CONFIG_OUTPUT_DIR)
            else:
                control_entries, treated_entries = download_experiment(
                    exp, FASTQ_OUTPUT_DIR, args.threads
                )
                generate_yaml_config(exp, control_entries, treated_entries, CONFIG_OUTPUT_DIR)
            logger.info("[EXPERIMENT %d/%d] %s — DONE", exp_idx, total_experiments, exp_label)
        except Exception as e:
            logger.error(
                "[EXPERIMENT %d/%d] %s — FAILED: %s",
                exp_idx, total_experiments, exp_label, e,
            )
            failed_experiments.append((exp_label, str(e)))

    sep = "=" * 60
    logger.info(sep)
    if failed_experiments:
        logger.error(
            "FINISHED WITH ERRORS — %d/%d experiment(s) failed:",
            len(failed_experiments), total_experiments,
        )
        for label, reason in failed_experiments:
            logger.error("  ✗ %s", label)
            logger.error("    Reason: %s", reason)
        logger.error(sep)
        return 1
    else:
        logger.info(
            "ALL DONE — %d/%d experiment(s) completed successfully.",
            total_experiments, total_experiments,
        )
        logger.info(sep)
        return 0


if __name__ == "__main__":
    sys.exit(main())
