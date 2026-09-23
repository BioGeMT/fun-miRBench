"""Cache published standardized predictor artifacts from the FuNmiRBench Zenodo record."""

from __future__ import annotations

import gzip
import pathlib
import shutil
import tempfile
import zipfile

import requests

from funmirbench.zenodo_store import (
    ZENODO_RECORD,
    compute_md5,
    fetch_zenodo_file_registry,
    parse_checksum,
)


def repo_root() -> pathlib.Path:
    return pathlib.Path(__file__).resolve().parents[1]


def _resolve_local_path(value: str | pathlib.Path, *, repo: pathlib.Path) -> pathlib.Path:
    path = pathlib.Path(value)
    if path.is_absolute():
        return path
    return repo / path


def _select_predictor_archive(registry: dict[str, dict]) -> dict:
    candidates = [
        meta
        for key, meta in registry.items()
        if key.lower().endswith(".zip")
        and "standardized" in key.lower()
        and ("prediction" in key.lower() or "predictor" in key.lower())
    ]
    if not candidates:
        zip_files = [
            meta for key, meta in registry.items() if key.lower().endswith(".zip")
        ]
        if len(zip_files) == 1:
            return zip_files[0]
        raise KeyError(
            f"Could not identify the standardized predictor archive in Zenodo record {ZENODO_RECORD}."
        )
    if len(candidates) > 1:
        names = sorted(str(meta["filename"]) for meta in candidates)
        raise ValueError(
            "Multiple standardized predictor archives found in Zenodo record "
            f"{ZENODO_RECORD}: {names}"
        )
    return candidates[0]


def _download_archive(meta: dict, *, directory: pathlib.Path, timeout: int) -> pathlib.Path:
    directory.mkdir(parents=True, exist_ok=True)
    response = requests.get(str(meta["url"]), stream=True, timeout=timeout)
    response.raise_for_status()

    with tempfile.NamedTemporaryFile(
        mode="wb",
        dir=directory,
        prefix=".funmirbench_predictors.",
        suffix=".zip",
        delete=False,
    ) as handle:
        archive_path = pathlib.Path(handle.name)
        for chunk in response.iter_content(chunk_size=1024 * 1024):
            if chunk:
                handle.write(chunk)

    checksum_value = str(meta.get("checksum", "") or "")
    if checksum_value:
        algorithm, expected_digest = parse_checksum(checksum_value)
        if algorithm != "md5":
            archive_path.unlink(missing_ok=True)
            raise ValueError(f"Unsupported checksum algorithm: {algorithm}")
        actual_digest = compute_md5(archive_path)
        if actual_digest != expected_digest:
            archive_path.unlink(missing_ok=True)
            raise ValueError(
                "Checksum mismatch for standardized predictor archive: "
                f"expected {expected_digest}, got {actual_digest}."
            )
    return archive_path


def _find_archive_member(
    archive: zipfile.ZipFile,
    filename: str,
) -> str:
    names = [name for name in archive.namelist() if not name.endswith("/")]
    target_names = {filename}
    if filename.endswith(".tsv"):
        target_names.add(f"{filename}.gz")

    matches = [
        name for name in names if pathlib.PurePosixPath(name).name in target_names
    ]
    if len(matches) != 1:
        raise KeyError(
            f"Expected exactly one archive member for {filename!r}, found {matches!r}."
        )
    return matches[0]


def _extract_member(
    archive: zipfile.ZipFile,
    member: str,
    destination: pathlib.Path,
) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="wb",
            dir=destination.parent,
            prefix=f".{destination.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            tmp_path = pathlib.Path(handle.name)
            with archive.open(member, "r") as source:
                if member.endswith(".gz") and not destination.name.endswith(".gz"):
                    with gzip.GzipFile(fileobj=source, mode="rb") as unpacked:
                        shutil.copyfileobj(unpacked, handle)
                else:
                    shutil.copyfileobj(source, handle)
        tmp_path.replace(destination)
        tmp_path = None
    finally:
        if tmp_path is not None and tmp_path.exists():
            tmp_path.unlink()


def sync_zenodo_predictors(
    predictions: dict[str, dict],
    *,
    repo: pathlib.Path | None = None,
    registry: dict[str, dict] | None = None,
    timeout: int = 120,
    force: bool = False,
) -> list[pathlib.Path]:
    """Ensure selected standardized predictor files are available locally."""
    repo = (repo or repo_root()).resolve()
    destinations = []
    for tool_id, meta in predictions.items():
        configured_path = meta.get("predictor_output_path")
        if not configured_path:
            raise ValueError(f"Predictor {tool_id!r} has no predictor_output_path.")
        destinations.append(
            (tool_id, _resolve_local_path(configured_path, repo=repo))
        )

    missing = [
        (tool_id, path)
        for tool_id, path in destinations
        if force or not path.is_file()
    ]
    if not missing:
        return [path for _, path in destinations]

    registry = registry or fetch_zenodo_file_registry(timeout=timeout)
    archive_meta = _select_predictor_archive(registry)
    archive_path = _download_archive(
        archive_meta,
        directory=repo / "data" / "predictions",
        timeout=timeout,
    )
    try:
        with zipfile.ZipFile(archive_path, "r") as archive:
            for tool_id, destination in missing:
                member = _find_archive_member(archive, destination.name)
                _extract_member(archive, member, destination)
                if not destination.is_file():
                    raise FileNotFoundError(
                        f"Failed to prepare predictor {tool_id!r} at {destination}."
                    )
    finally:
        archive_path.unlink(missing_ok=True)

    return [path for _, path in destinations]
