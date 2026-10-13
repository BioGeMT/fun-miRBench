"""Helpers for caching benchmark experiment DE tables from local files or Zenodo."""

from __future__ import annotations

import argparse
import json
import logging
import pathlib
import tempfile

import pandas as pd
import requests

from funmirbench.zenodo_store import (
    ZENODO_RECORD,
    compute_md5,
    fetch_zenodo_file_registry,
    parse_checksum,
)


logger = logging.getLogger(__name__)


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
    try:
        remote_filename = _resolve_zenodo_filename(filename, registry)
    except KeyError:
        return
    meta = registry[remote_filename]
    checksum_value = str(meta.get("checksum", "") or "")
    if not checksum_value:
        return
    algorithm, expected_digest = parse_checksum(checksum_value)
    if algorithm != "md5":
        raise ValueError(f"Unsupported checksum algorithm: {algorithm}")

    if remote_filename != filename:
        _verify_cache_manifest(
            dest,
            remote_filename=remote_filename,
            remote_checksum=checksum_value,
        )
        return

    if compute_md5(dest) != expected_digest:
        raise ValueError(
            f"Existing file {dest} failed checksum verification for {filename}."
        )


def _resolve_zenodo_filename(filename: str, registry: dict[str, dict]) -> str:
    """Resolve a curated table filename against the Zenodo record.

    Metadata keeps stable local .tsv paths while the published artifact may
    store the same table either as .tsv or .tsv.gz.
    """
    candidates = [filename]
    if filename.endswith(".tsv"):
        candidates.append(f"{filename}.gz")
    elif filename.endswith(".tsv.gz"):
        candidates.append(filename[:-3])

    for candidate in candidates:
        if candidate in registry:
            return candidate
    raise KeyError(
        f"{filename!r} is not present in Zenodo record {ZENODO_RECORD} "
        f"(checked {candidates!r})."
    )


def _cache_manifest_path(dest: pathlib.Path) -> pathlib.Path:
    return dest.with_name(f"{dest.name}.funmirbench.json")


def _write_cache_manifest(
    dest: pathlib.Path,
    *,
    remote_filename: str,
    remote_checksum: str,
    remote_size: int,
) -> None:
    manifest = {
        "zenodo_record": ZENODO_RECORD,
        "local_filename": dest.name,
        "local_md5": compute_md5(dest),
        "remote_filename": remote_filename,
        "remote_checksum": remote_checksum,
        "remote_size": int(remote_size or 0),
    }
    _cache_manifest_path(dest).write_text(
        json.dumps(manifest, indent=2) + "\n",
        encoding="utf-8",
    )


def _verify_cache_manifest(
    dest: pathlib.Path,
    *,
    remote_filename: str,
    remote_checksum: str,
) -> None:
    manifest_path = _cache_manifest_path(dest)
    if not manifest_path.is_file():
        raise ValueError(
            f"Existing cached experiment {dest} cannot be checksum-verified "
            f"against compressed Zenodo artifact {remote_filename} without "
            f"sidecar manifest {manifest_path}. Rerun with --force to refresh "
            "this cache."
        )

    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    expected = {
        "zenodo_record": ZENODO_RECORD,
        "remote_filename": remote_filename,
        "remote_checksum": remote_checksum,
    }
    for key, value in expected.items():
        if manifest.get(key) != value:
            raise ValueError(
                f"Existing cache manifest {manifest_path} does not match Zenodo {key}: "
                f"expected {value!r}, got {manifest.get(key)!r}."
            )

    local_md5 = str(manifest.get("local_md5") or "")
    if not local_md5:
        raise ValueError(f"Existing cache manifest {manifest_path} is missing local_md5.")
    if compute_md5(dest) != local_md5:
        raise ValueError(
            f"Existing cached file {dest} failed sidecar checksum verification."
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

    remote_filename = _resolve_zenodo_filename(filename, registry)
    meta = registry[remote_filename]
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

        if remote_filename.endswith(".gz") and not filename.endswith(".gz"):
            import gzip
            import shutil

            with gzip.open(tmp_path, "rb") as source, tempfile.NamedTemporaryFile(
                mode="wb",
                dir=dest.parent,
                prefix=f".{dest.name}.",
                suffix=".tmp",
                delete=False,
            ) as handle:
                unpacked_path = pathlib.Path(handle.name)
                shutil.copyfileobj(source, handle)
            unpacked_path.replace(dest)
            tmp_path.unlink()
            tmp_path = None
        else:
            tmp_path.replace(dest)
            tmp_path = None
        _write_cache_manifest(
            dest,
            remote_filename=remote_filename,
            remote_checksum=checksum_value,
            remote_size=int(meta.get("size", 0) or 0),
        )
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
    published_root = experiments_processed_dir(repo=repo).resolve()
    for de_table_path in de_table_paths:
        rel_path = pathlib.Path(de_table_path)
        key = str(rel_path)
        if key in seen:
            continue
        seen.add(key)

        dest = resolve_cached_experiment_path(rel_path, repo=repo).resolve()
        if not dest.is_relative_to(published_root):
            logger.info("Skipping non-Zenodo experiment path: %s", dest)
            continue

        if registry_cache is None:
            registry_cache = fetch_zenodo_file_registry(token=token, timeout=timeout)

        if dest.exists() and not force:
            logger.info("Validating cached experiment: %s", dest)
        else:
            logger.info("Downloading experiment: %s", dest)
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
