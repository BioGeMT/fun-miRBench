"""Helpers for caching benchmark experiment DE tables from local files or Zenodo."""

from __future__ import annotations

import argparse
import hashlib
import pathlib
import tempfile
from urllib.parse import quote

import pandas as pd
import requests


ZENODO_RECORD = "21671476"
ZENODO_API_RECORD_URL = f"https://zenodo.org/api/records/{ZENODO_RECORD}"


def repo_root() -> pathlib.Path:
    return pathlib.Path(__file__).resolve().parents[1]


def experiments_cache_root_dir(*, repo: pathlib.Path | None = None) -> pathlib.Path:
    repo = (repo or repo_root()).resolve()
    return repo / "data" / "experiments" / "processed"


def experiments_processed_dir(*, repo: pathlib.Path | None = None) -> pathlib.Path:
    return experiments_cache_root_dir(repo=repo) / ZENODO_RECORD


def experiments_metadata_tsv(*, repo: pathlib.Path | None = None) -> pathlib.Path:
    return (repo or repo_root()).resolve() / "metadata" / "mirna_experiment_info.tsv"


def experiment_cache_relpath(filename: str | pathlib.Path) -> pathlib.Path:
    return pathlib.Path("data") / "experiments" / "processed" / ZENODO_RECORD / pathlib.Path(filename).name


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
    return (
        f"{ZENODO_API_RECORD_URL}/files/{quote(str(filename), safe='')}/content"
    )


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
            "url": item.get("links", {}).get("self") or zenodo_download_url(key),
        }
    return registry


def resolve_cached_experiment_path(
    de_table_path: str | pathlib.Path,
    *,
    repo: pathlib.Path | None = None,
) -> pathlib.Path:
    path = pathlib.Path(de_table_path)
    if path.is_absolute():
        return path
    return (repo or repo_root()).resolve() / path


def _verify_checksum_if_available(dest: pathlib.Path, filename: str, registry: dict[str, dict]) -> None:
    meta = registry.get(filename)
    if not meta:
        return
    checksum_value = str(meta.get("checksum", "") or "")
    if not checksum_value:
        return
    algorithm, expected_digest = parse_checksum(checksum_value)
    if algorithm != "md5":
        raise ValueError(f"Unsupported checksum algorithm: {algorithm}")
    if compute_md5(dest) != expected_digest:
        raise ValueError(
            f"Existing file {dest} failed checksum verification for {filename}."
        )


def ensure_zenodo_experiment_cached(
    de_table_path: str | pathlib.Path,
    *,
    repo: pathlib.Path | None = None,
    registry: dict[str, dict] | None = None,
    token: str | None = None,
    timeout: int = 120,
    force: bool = False,
) -> pathlib.Path:
    """Return a usable DE table path, using local files first and Zenodo as fallback.

    If ``de_table_path`` already exists locally, it is accepted as-is. When the
    local filename also exists in the Zenodo registry, its checksum is verified.
    If the local file is missing, the filename must be present in the Zenodo
    record and is downloaded to ``de_table_path``.
    """
    dest = resolve_cached_experiment_path(de_table_path, repo=repo)
    filename = dest.name

    if dest.exists() and not force:
        if registry is not None:
            _verify_checksum_if_available(dest, filename, registry)
        return dest

    registry = registry or fetch_zenodo_file_registry(token=token, timeout=timeout)
    if dest.exists() and not force:
        _verify_checksum_if_available(dest, filename, registry)
        return dest

    if filename not in registry:
        raise KeyError(
            f"{filename!r} is not present in Zenodo record {ZENODO_RECORD} and no local file exists at {dest}."
        )

    meta = registry[filename]
    checksum_value = str(meta.get("checksum", "") or "")

    dest.parent.mkdir(parents=True, exist_ok=True)
    response = requests.get(str(meta["url"]), stream=True, timeout=timeout)
    response.raise_for_status()
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
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    handle.write(chunk)

        if checksum_value:
            algorithm, expected_digest = parse_checksum(checksum_value)
            if algorithm != "md5":
                raise ValueError(f"Unsupported checksum algorithm: {algorithm}")
            actual_digest = compute_md5(tmp_path)
            if actual_digest != expected_digest:
                raise ValueError(
                    f"Checksum mismatch for {filename}: expected {expected_digest}, got {actual_digest}."
                )

        tmp_path.replace(dest)
        tmp_path = None
        return dest
    finally:
        if tmp_path is not None and tmp_path.exists():
            tmp_path.unlink()


def sync_all_zenodo_experiments(
    *,
    repo: pathlib.Path | None = None,
    token: str | None = None,
    timeout: int = 120,
    force: bool = False,
    ) -> list[pathlib.Path]:
    repo = (repo or repo_root()).resolve()
    df = pd.read_csv(experiments_metadata_tsv(repo=repo), sep="\t")
    de_table_paths = [str(value) for value in df["de_table_path"].dropna().tolist()]

    return sync_zenodo_experiments(
        de_table_paths,
        repo=repo,
        token=token,
        timeout=timeout,
        force=force,
    )


def sync_zenodo_experiments(
    de_table_paths: list[str | pathlib.Path],
    *,
    repo: pathlib.Path | None = None,
    registry: dict[str, dict] | None = None,
    token: str | None = None,
    timeout: int = 120,
    force: bool = False,
) -> list[pathlib.Path]:
    repo = (repo or repo_root()).resolve()
    registry_cache = registry

    saved = []
    seen = set()
    for de_table_path in de_table_paths:
        rel_path = pathlib.Path(de_table_path)
        key = str(rel_path)
        if key in seen:
            continue
        seen.add(key)

        dest = resolve_cached_experiment_path(rel_path, repo=repo)
        if dest.exists() and not force:
            print(f"local {rel_path.name}")
            saved.append(dest)
            continue

        if registry_cache is None:
            registry_cache = fetch_zenodo_file_registry(token=token, timeout=timeout)
        print(f"sync {rel_path.name}")
        saved.append(
            ensure_zenodo_experiment_cached(
                rel_path,
                repo=repo,
                registry=registry_cache,
                timeout=timeout,
                force=force,
            )
        )
    return saved



def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Sync benchmark experiment DE tables from local files or Zenodo."
    )
    parser.add_argument("--repo", type=pathlib.Path, default=None)
    parser.add_argument("--token", default=None)
    parser.add_argument("--timeout", type=int, default=120)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    saved = sync_all_zenodo_experiments(
        repo=args.repo,
        token=args.token,
        timeout=args.timeout,
        force=args.force,
    )
    print(
        "Prepared "
        f"{len(saved)} experiment tables under data/experiments/processed/"
    )
