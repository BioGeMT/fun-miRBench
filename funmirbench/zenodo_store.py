"""Shared helpers for accessing the fun-miRBench Zenodo record."""

from __future__ import annotations

import hashlib
import pathlib
from urllib.parse import quote

import requests


ZENODO_RECORD = "21671476"
ZENODO_API_RECORD_URL = f"https://zenodo.org/api/records/{ZENODO_RECORD}"


def compute_md5(path: str | pathlib.Path, *, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.md5()
    path = pathlib.Path(path)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_checksum(value: str) -> tuple[str, str]:
    algorithm, _, digest = str(value).partition(":")
    if not algorithm or not digest:
        raise ValueError(f"Unsupported checksum value: {value!r}")
    return algorithm.lower(), digest.lower()


def zenodo_download_url(filename: str) -> str:
    return f"{ZENODO_API_RECORD_URL}/files/{quote(str(filename), safe='')}/content"


def fetch_zenodo_file_registry(
    *,
    token: str | None = None,
    timeout: int = 120,
) -> dict[str, dict]:
    url = ZENODO_API_RECORD_URL
    if token:
        url = f"{url}?token={token}"

    response = requests.get(url, timeout=timeout)
    response.raise_for_status()
    payload = response.json()

    registry = {}
    for item in payload.get("files", []):
        key = str(item["key"])
        registry[key] = {
            "filename": key,
            "size": int(item.get("size", 0) or 0),
            "checksum": str(item.get("checksum", "")),
            "url": item.get("links", {}).get("content") or zenodo_download_url(key),
        }
    return registry
