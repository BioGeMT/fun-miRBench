import gzip
import logging
import math
import re
from pathlib import Path
from typing import Optional

import pandas as pd
import requests

from common import display_path

logger = logging.getLogger("miraw_utils")

MIRBASE_RELEASE = "22.1"
ENSEMBL_SOURCE_RELEASE = 112
ENSEMBL_TARGET_RELEASE = 115
ENSEMBL_GENE_ID_RE = re.compile(r"^ENSG\d{11}$")
MIMAT_ID_RE = re.compile(r"^MIMAT\d+$")


def repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def resolve_path_relative_to_root(path: Path) -> Path:
    try:
        return path.resolve().relative_to(repo_root())
    except ValueError:
        return path


def configure_logging(log_path: Path, log_level: str) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    level = getattr(logging, log_level.upper(), logging.INFO)
    logging.basicConfig(
        level=level,
        format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler(log_path, mode="w", encoding="utf-8"),
        ],
        force=True,
    )


def _log_row_count_change(label: str, before: int, after: int) -> None:
    logger.info("%s: %d -> %d rows", label, before, after)


def download_file(
    url: str,
    output_path: Path,
    params: Optional[dict[str, str]] = None,
    timeout: int = 120,
    resource_label: str = "file",
) -> Path:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    relative_path = resolve_path_relative_to_root(output_path)
    if output_path.exists():
        logger.info("Using cached %s: %s", resource_label, relative_path)
        return output_path

    logger.info("Downloading %s: %s", resource_label, relative_path)
    try:
        with requests.get(url, params=params, timeout=timeout, stream=True) as response:
            response.raise_for_status()
            with output_path.open("wb") as handle:
                wrote_bytes = 0
                for chunk in response.iter_content(chunk_size=1024 * 1024):
                    if not chunk:
                        continue
                    handle.write(chunk)
                    wrote_bytes += len(chunk)
    except requests.RequestException as exc:
        output_path.unlink(missing_ok=True)
        raise RuntimeError(f"Failed to download {resource_label} from {url}: {exc}") from exc
    except OSError:
        output_path.unlink(missing_ok=True)
        raise

    if wrote_bytes == 0:
        output_path.unlink(missing_ok=True)
        raise RuntimeError(f"Failed to download {resource_label}: empty response from {url}")

    logger.info("Saved %s: %s (%d bytes)", resource_label, relative_path, wrote_bytes)
    return output_path


def create_mirna_name_to_mimat_mapping(mature_fa_path: Path) -> dict[str, str]:
    if not mature_fa_path.exists():
        raise FileNotFoundError(f"Missing miRBase {MIRBASE_RELEASE} mature FASTA: {mature_fa_path}")

    mapping: dict[str, str] = {}
    human_headers = 0
    with mature_fa_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            parts = line[1:].strip().split()
            if len(parts) < 2:
                raise ValueError(f"{mature_fa_path}: invalid FASTA header: {line.strip()}")

            mirna_name, mimat_id = parts[0].strip(), parts[1].strip()
            if not mirna_name.startswith("hsa-"):
                continue
            human_headers += 1
            if not MIMAT_ID_RE.fullmatch(mimat_id):
                raise ValueError(
                    f"{mature_fa_path}: invalid miRBase {MIRBASE_RELEASE} mature accession: {mimat_id}"
                )
            previous = mapping.get(mirna_name)
            if previous is not None and previous != mimat_id:
                raise ValueError(
                    f"{mature_fa_path}: conflicting MIMAT accessions for {mirna_name}: "
                    f"{previous} versus {mimat_id}"
                )
            mapping[mirna_name] = mimat_id

    if not mapping:
        raise RuntimeError(
            f"No human mature miRNAs were parsed from miRBase {MIRBASE_RELEASE}: {mature_fa_path}"
        )
    logger.info(
        "Loaded %d unique human miRNA-name-to-MIMAT mappings from %d miRBase %s FASTA headers",
        len(mapping),
        human_headers,
        MIRBASE_RELEASE,
    )
    return mapping


def _strip_ensembl_version(value: object) -> str:
    return str(value).strip().split(".", 1)[0]


def _parse_gtf_attributes(raw_attrs: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for item in raw_attrs.strip().split(";"):
        item = item.strip()
        if not item or " " not in item:
            continue
        key, value = item.split(" ", 1)
        attrs[key] = value.strip().strip('"')
    return attrs


def create_ensembl_to_gene_name_mapping_from_gtf(gtf_gz_path: Path) -> dict[str, str]:
    if not gtf_gz_path.exists():
        raise FileNotFoundError(
            f"Missing Ensembl release {ENSEMBL_TARGET_RELEASE} GTF: {gtf_gz_path}"
        )

    mapping: dict[str, str] = {}
    gene_rows = 0
    invalid_ids = 0
    missing_names = 0
    conflicts = 0
    with gzip.open(gtf_gz_path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            gene_rows += 1
            attrs = _parse_gtf_attributes(fields[8])
            gene_id = _strip_ensembl_version(attrs.get("gene_id", ""))
            gene_name = attrs.get("gene_name", "").strip()
            if not ENSEMBL_GENE_ID_RE.fullmatch(gene_id):
                invalid_ids += 1
                continue
            if not gene_name:
                missing_names += 1
                continue
            previous = mapping.get(gene_id)
            if previous is not None and previous != gene_name:
                conflicts += 1
                continue
            mapping[gene_id] = gene_name

    if not mapping:
        raise RuntimeError(
            f"No valid human ENSG-to-gene-name mappings were parsed from Ensembl "
            f"release {ENSEMBL_TARGET_RELEASE}: {gtf_gz_path}"
        )
    logger.info(
        "Loaded %d Ensembl release %d ENSG-to-gene-name mappings from %d gene rows; "
        "invalid_ids=%d, missing_names=%d, conflicts=%d",
        len(mapping),
        ENSEMBL_TARGET_RELEASE,
        gene_rows,
        invalid_ids,
        missing_names,
        conflicts,
    )
    return mapping


def _extract_raw_gene_name_suffix(raw_gene_value: str) -> Optional[str]:
    if "__" not in raw_gene_value:
        return None
    suffix = raw_gene_value.split("__", 1)[1].strip()
    if not suffix or suffix.startswith("biotype="):
        return None
    return suffix


def _clean_loaded_miraw_predictions(
    df: pd.DataFrame,
    duplicate_subset: list[str],
) -> pd.DataFrame:
    before = len(df)
    df = df.dropna(subset=["Ensembl_ID", "miRNA_Name", "Score"]).copy()
    df["Ensembl_ID"] = df["Ensembl_ID"].map(_strip_ensembl_version)
    df["Raw_Gene_Name"] = df["Raw_Gene_Name"].fillna("").astype(str).str.strip()
    df["miRNA_Name"] = df["miRNA_Name"].astype(str).str.strip()
    df["Score"] = pd.to_numeric(df["Score"], errors="coerce")
    df = df.dropna(subset=["Score"])
    df = df[df["Score"].map(math.isfinite)].copy()
    df = df[
        df["Ensembl_ID"].map(lambda value: bool(ENSEMBL_GENE_ID_RE.fullmatch(value)))
        & df["miRNA_Name"].str.startswith("hsa-")
    ].copy()
    _log_row_count_change("Validate parsed miRAW rows", before, len(df))

    before = len(df)
    df = df.drop_duplicates(subset=duplicate_subset).copy()
    _log_row_count_change("Drop exact duplicate parsed miRAW rows", before, len(df))
    return df


def load_miraw_predictions(predictions_path: Path) -> pd.DataFrame:
    if not predictions_path.exists():
        raise FileNotFoundError(f"Missing prediction file: {predictions_path}")

    df = pd.read_csv(
        predictions_path,
        sep="\t",
        dtype={"Target_ENSG": "string", "GeneName": "string", "miRNA": "string"},
        keep_default_na=False,
    )
    required = {"Target_ENSG", "GeneName", "miRNA", "Prediction"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"{predictions_path}: missing required Figshare miRAW columns: {sorted(missing)}"
        )

    target_ensg = df["Target_ENSG"].astype(str).str.strip()
    gene_name_field = df["GeneName"].astype(str).str.strip()
    target_id = target_ensg.str.split("__", n=1).str[0]
    gene_id = gene_name_field.str.split("__", n=1).str[0]
    disagree = (target_id != "") & (gene_id != "") & (target_id != gene_id)
    if disagree.any():
        examples = df.loc[disagree, ["Target_ENSG", "GeneName"]].head(5).to_dict("records")
        raise ValueError(
            f"{predictions_path}: Target_ENSG and GeneName disagree for "
            f"{int(disagree.sum())} rows; examples={examples}"
        )

    ensembl_id = target_id.where(target_id != "", gene_id)
    out = pd.DataFrame(
        {
            "Raw_Gene_Field": gene_name_field,
            "Ensembl_ID": ensembl_id,
            "Raw_Gene_Name": gene_name_field.map(_extract_raw_gene_name_suffix),
            "miRNA_Name": df["miRNA"],
            "Score": df["Prediction"],
        }
    )
    logger.info(
        "Loaded %d miRAW prediction rows generated against Ensembl release %d and "
        "miRBase release %s from %s",
        len(out),
        ENSEMBL_SOURCE_RELEASE,
        MIRBASE_RELEASE,
        display_path(predictions_path),
    )
    return _clean_loaded_miraw_predictions(
        out,
        duplicate_subset=["Ensembl_ID", "miRNA_Name", "Score", "Raw_Gene_Field"],
    )


def collapse_to_best_score_per_pair(
    df: pd.DataFrame,
    ensembl_id_column: str,
    mirna_name_column: str,
    score_column: str,
) -> pd.DataFrame:
    before = len(df)
    out = (
        df.sort_values(
            [ensembl_id_column, mirna_name_column, score_column],
            ascending=[True, True, False],
        )
        .drop_duplicates(subset=[ensembl_id_column, mirna_name_column], keep="first")
        .copy()
    )
    _log_row_count_change(
        "Collapse site-level miRAW rows to best score per Ensembl_ID-miRNA_Name pair",
        before,
        len(out),
    )
    return out


def map_mirna_names_to_mimat(
    df: pd.DataFrame,
    mirna_name_to_id: dict[str, str],
    mirna_name_column: str,
    mimat_column: str,
) -> pd.DataFrame:
    out = df.copy()
    out[mimat_column] = out[mirna_name_column].map(mirna_name_to_id)
    before = len(out)
    missing_names = sorted(out.loc[out[mimat_column].isna(), mirna_name_column].unique())
    out = out.dropna(subset=[mimat_column]).copy()
    _log_row_count_change(
        f"Retain miRNAs present in miRBase release {MIRBASE_RELEASE}",
        before,
        len(out),
    )
    if missing_names:
        logger.info(
            "Dropped %d unique miRNA names absent from miRBase %s; first examples=%s",
            len(missing_names),
            MIRBASE_RELEASE,
            missing_names[:10],
        )
    if not out[mimat_column].map(lambda value: bool(MIMAT_ID_RE.fullmatch(str(value)))).all():
        raise ValueError(f"Invalid MIMAT accession produced from miRBase {MIRBASE_RELEASE}")
    return out


def map_ensembl_to_gene_name(
    df: pd.DataFrame,
    ensembl_to_gene_name_map: dict[str, str],
    ensembl_id_column: str,
    gene_name_column: str,
    raw_gene_name_column: str,
) -> pd.DataFrame:
    del raw_gene_name_column
    out = df.copy()
    out[gene_name_column] = out[ensembl_id_column].map(ensembl_to_gene_name_map)
    before = len(out)
    missing_ids = sorted(out.loc[out[gene_name_column].isna(), ensembl_id_column].unique())
    out = out.dropna(subset=[gene_name_column]).copy()
    _log_row_count_change(
        f"Retain genes present in Ensembl release {ENSEMBL_TARGET_RELEASE}",
        before,
        len(out),
    )
    if missing_ids:
        logger.info(
            "Dropped %d unique Ensembl release %d source IDs absent from release %d; "
            "first examples=%s",
            len(missing_ids),
            ENSEMBL_SOURCE_RELEASE,
            ENSEMBL_TARGET_RELEASE,
            missing_ids[:10],
        )
    out[gene_name_column] = out[gene_name_column].astype(str).str.strip()
    out = out[out[gene_name_column] != ""].copy()
    return out


def build_output_table(
    df: pd.DataFrame,
    final_columns: list[str],
    ensembl_id_column: str,
    mimat_column: str,
    score_column: str,
) -> pd.DataFrame:
    out = df.loc[:, final_columns].copy()
    before = len(out)
    out = out.dropna(subset=[ensembl_id_column, mimat_column, score_column]).copy()
    for column in final_columns:
        if column != score_column:
            out[column] = out[column].astype(str).str.strip()
    out[score_column] = pd.to_numeric(out[score_column], errors="coerce")
    out = out.dropna(subset=[score_column])
    out = out[
        out[ensembl_id_column].map(lambda value: bool(ENSEMBL_GENE_ID_RE.fullmatch(value)))
        & out[mimat_column].map(lambda value: bool(MIMAT_ID_RE.fullmatch(value)))
        & out["miRNA_Name"].str.startswith("hsa-")
        & (out["Gene_Name"] != "")
    ].copy()
    _log_row_count_change("Validate standardized miRAW identifiers and values", before, len(out))

    before = len(out)
    out = (
        out.sort_values(
            [ensembl_id_column, mimat_column, score_column, "miRNA_Name"],
            ascending=[True, True, False, True],
        )
        .drop_duplicates(subset=[ensembl_id_column, mimat_column], keep="first")
        .copy()
    )
    _log_row_count_change(
        "Finalize one highest-scoring row per Ensembl_ID-miRNA_ID pair",
        before,
        len(out),
    )

    if out.duplicated(subset=[ensembl_id_column, mimat_column]).any():
        raise ValueError("Final miRAW output contains duplicate Ensembl_ID-miRNA_ID pairs")
    if out.empty:
        raise RuntimeError("No miRAW rows remain after release-aware standardization")

    logger.info(
        "Final miRAW table uses Ensembl release %d gene IDs/names and miRBase release %s "
        "mature miRNA names/accessions",
        ENSEMBL_TARGET_RELEASE,
        MIRBASE_RELEASE,
    )
    return out.sort_values([ensembl_id_column, mimat_column]).reset_index(drop=True)
