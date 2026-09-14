"""Protein-coding gene-set helpers for benchmark evaluation."""

from __future__ import annotations

import gzip
import logging
import pathlib
import shutil
import tempfile
import time
import urllib.error
import urllib.request

from funmirbench.gene_ids import strip_ensembl_version


ENSEMBL_RELEASE = "115"
DEFAULT_ENSEMBL_GTF_URL = (
    "https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/"
    "Homo_sapiens.GRCh38.115.gtf.gz"
)
DEFAULT_GTF_REL_PATH = pathlib.Path(
    "data/common_resources/ensembl/Homo_sapiens.GRCh38.115.gtf.gz"
)
DEFAULT_CACHE_REL_PATH = pathlib.Path(
    "data/common_resources/ensembl/protein_coding_gene_ids.txt"
)
DOWNLOAD_RETRIES = 3
DOWNLOAD_TIMEOUT_SECONDS = 60

logger = logging.getLogger(__name__)


def _split_gtf_attributes(raw_attrs: str) -> list[str]:
    items = []
    current = []
    in_quotes = False
    escaped = False
    for char in raw_attrs.strip():
        if escaped:
            current.append(char)
            escaped = False
            continue
        if char == "\\":
            current.append(char)
            escaped = True
            continue
        if char == '"':
            in_quotes = not in_quotes
            current.append(char)
            continue
        if char == ";" and not in_quotes:
            items.append("".join(current))
            current = []
            continue
        current.append(char)
    if current:
        items.append("".join(current))
    return items


def _parse_gtf_attributes(raw_attrs: str) -> dict[str, str]:
    attrs = {}
    for item in _split_gtf_attributes(raw_attrs):
        item = item.strip()
        if not item or " " not in item:
            continue
        key, value = item.split(" ", 1)
        attrs[key] = value.strip().strip('"')
    return attrs


def _download_file(
    url: str,
    dest: pathlib.Path,
    *,
    retries: int = DOWNLOAD_RETRIES,
    timeout: int = DOWNLOAD_TIMEOUT_SECONDS,
) -> pathlib.Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    last_error = None
    for attempt in range(1, retries + 1):
        tmp_path = None
        try:
            with tempfile.NamedTemporaryFile(
                mode="wb",
                dir=dest.parent,
                prefix=f".{dest.name}.",
                suffix=".tmp",
                delete=False,
            ) as handle:
                tmp_path = pathlib.Path(handle.name)
                with urllib.request.urlopen(url, timeout=timeout) as response:
                    shutil.copyfileobj(response, handle)
            if tmp_path.stat().st_size == 0:
                raise RuntimeError(f"Downloaded empty file from {url}")
            tmp_path.replace(dest)
            tmp_path = None
            return dest
        except (OSError, RuntimeError, urllib.error.URLError, TimeoutError) as exc:
            last_error = exc
            if attempt >= retries:
                break
            logger.warning(
                "Download failed for %s (attempt %d/%d): %s",
                url,
                attempt,
                retries,
                exc,
            )
            time.sleep(min(2 ** (attempt - 1), 5))
        finally:
            if tmp_path is not None and tmp_path.exists():
                tmp_path.unlink()
    raise RuntimeError(
        f"Failed to download {url} to {dest} after {retries} attempts"
    ) from last_error


def ensure_ensembl_gtf(
    *,
    root: pathlib.Path,
    gtf_path: str | pathlib.Path | None = None,
    url: str = DEFAULT_ENSEMBL_GTF_URL,
) -> pathlib.Path:
    """Return a local Ensembl GTF path, downloading the default cache if needed."""
    path = pathlib.Path(gtf_path) if gtf_path is not None else DEFAULT_GTF_REL_PATH
    if not path.is_absolute():
        path = root / path
    if path.exists():
        return path
    logger.info("Downloading Ensembl GTF for protein-coding filter: %s", path)
    return _download_file(url, path)


def parse_protein_coding_gene_ids_from_gtf(gtf_gz_path: pathlib.Path) -> set[str]:
    """Parse protein-coding Ensembl gene IDs from a gzipped Ensembl GTF."""
    gene_ids: set[str] = set()
    with gzip.open(gtf_gz_path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            attrs = _parse_gtf_attributes(fields[8])
            biotype = attrs.get("gene_biotype") or attrs.get("gene_type")
            if biotype != "protein_coding":
                continue
            gene_id = attrs.get("gene_id")
            if gene_id:
                gene_ids.add(strip_ensembl_version(gene_id))
    if not gene_ids:
        raise ValueError(f"No protein-coding gene IDs parsed from {gtf_gz_path}")
    return gene_ids


def read_gene_id_cache(cache_path: pathlib.Path) -> set[str]:
    return {
        strip_ensembl_version(line)
        for line in cache_path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    }


def write_gene_id_cache(cache_path: pathlib.Path, gene_ids: set[str]) -> None:
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache_path.write_text("\n".join(sorted(gene_ids)) + "\n", encoding="utf-8")


def load_protein_coding_gene_ids(
    *,
    root: pathlib.Path,
    gtf_path: str | pathlib.Path | None = None,
    cache_path: str | pathlib.Path | None = None,
) -> set[str]:
    """Load cached protein-coding Ensembl IDs, building the cache from Ensembl GTF if needed."""
    root = root.resolve()
    cache = pathlib.Path(cache_path) if cache_path is not None else DEFAULT_CACHE_REL_PATH
    if not cache.is_absolute():
        cache = root / cache
    if cache.exists():
        gene_ids = read_gene_id_cache(cache)
        if gene_ids:
            logger.info("Loaded %d protein-coding gene IDs from %s", len(gene_ids), cache)
            return gene_ids

    gtf = ensure_ensembl_gtf(root=root, gtf_path=gtf_path)
    gene_ids = parse_protein_coding_gene_ids_from_gtf(gtf)
    write_gene_id_cache(cache, gene_ids)
    logger.info("Cached %d protein-coding gene IDs to %s", len(gene_ids), cache)
    return gene_ids
