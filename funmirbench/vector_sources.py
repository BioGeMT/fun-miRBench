"""Isolated, pinned Zenodo inputs for exploratory vector analysis."""
from __future__ import annotations

import gzip
import hashlib
import json
import logging
import pathlib
import shutil
import tempfile
import time
import zipfile
from contextlib import contextmanager

import requests

from funmirbench.predictor_store import PUBLISHED_PREDICTOR_FILENAMES, _find_archive_member
from funmirbench.zenodo_store import compute_md5, fetch_zenodo_file_registry

logger = logging.getLogger(__name__)
PREDICTOR_RECORD = "21671476"
# File sizes/checksums verified against these published records on 2026-11-04.
# Public file routes avoid depending on the metadata API during an outage.
PINNED_FILES = {
    "23839843": {
        "benchmark_data.zip": (85741740, "42bd622a8760671ca80c477e7fa81b87"),
        "targetscan_family_expanded_standardized.tsv.gz":
            (120368347, "6110bb5d9b025793f84960c095114c7f"),
    },
    "21671476": {
        "fun-miRBench_standardized_predictions.zip":
            (2011298508, "9e619a17668ea780fa3f4de85a7f5476"),
    },
}


def atomic_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode="w", dir=path.parent, delete=False, encoding="utf-8") as handle:
        temporary = pathlib.Path(handle.name)
        json.dump(value, handle, indent=2)
    temporary.replace(path)


def fingerprint(path):
    path = pathlib.Path(path).resolve()
    stat = path.stat()
    return {"path": str(path), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns}


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def registry(record, cache, *, refresh=False, timeout=120):
    record = str(record)
    path = cache / "zenodo" / str(record) / "registry.json"
    if path.is_file() and not refresh:
        return json.loads(path.read_text())
    for attempt in range(3):
        try:
            result = fetch_zenodo_file_registry(record=record, timeout=timeout)
            if not result:
                raise ValueError(f"empty file registry for Zenodo {record}")
            break
        except (requests.RequestException, ValueError) as exc:
            # Invalid/private record IDs need correction, rather than a fallback.
            if isinstance(exc, requests.HTTPError) and exc.response is not None:
                status = exc.response.status_code
                if 400 <= status < 500 and status not in (408, 429):
                    raise
            if attempt < 2:
                delay = 2 ** (attempt + 1)
                logger.warning("Zenodo %s metadata unavailable (%s); retry %d/3 in %ss",
                               record, exc, attempt + 2, delay)
                time.sleep(delay)
                continue
            if record not in PINNED_FILES:
                raise
            logger.warning("Zenodo %s API unavailable; using bundled checksum-pinned file registry", record)
            result = {name: {"filename": name, "size": size, "checksum": "md5:" + checksum,
                             "url": f"https://zenodo.org/records/{record}/files/{name}?download=1"}
                      for name, (size, checksum) in PINNED_FILES[record].items()}
    atomic_json(path, result)
    return result


def download(meta, directory, *, timeout=120):
    """Cache a checksum-verified archive; retain partial transfers for HTTP Range retries."""
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / pathlib.PurePosixPath(meta["filename"]).name
    marker = path.with_name(path.name + ".verified.json")
    expected = meta["checksum"]
    if not expected.startswith("md5:"):
        raise ValueError("exploratory downloads require the Zenodo MD5 checksum")
    if path.is_file():
        state = {**fingerprint(path), "checksum": expected}
        if marker.is_file() and json.loads(marker.read_text()) == state:
            logger.info("cached %s", path.name)
            return path
        if path.stat().st_size == meta["size"] and compute_md5(path) == expected[4:]:
            atomic_json(marker, state)
            return path
    part = path.with_name(path.name + ".part")
    for attempt in range(3):
        try:
            if part.exists() and part.stat().st_size >= meta["size"]:
                if part.stat().st_size == meta["size"] and compute_md5(part) == expected[4:]:
                    part.replace(path)
                    atomic_json(marker, {**fingerprint(path), "checksum": expected})
                    return path
                part.unlink()
            offset = part.stat().st_size if part.exists() else 0
            logger.info("downloading %s (%.1f MB; resume at %.1f MB)", path.name,
                        meta["size"] / 1e6, offset / 1e6)
            headers = {"Range": f"bytes={offset}-"} if offset else {}
            with requests.get(meta["url"], stream=True, timeout=timeout, headers=headers) as response:
                response.raise_for_status()
                if response.status_code == 206:
                    if not response.headers.get("Content-Range", "").startswith(f"bytes {offset}-"):
                        raise ValueError("unexpected Content-Range")
                else:
                    offset = 0  # server ignored Range: restart, never append a full response
                written, last_percent = offset, -1
                with part.open("ab" if offset else "wb") as handle:
                    for block in response.iter_content(4 << 20):
                        handle.write(block)
                        written += len(block)
                        percent = int(written * 100 / meta["size"])
                        if percent // 10 > last_percent:
                            logger.info("%s: %d%%", path.name, min(100, percent))
                            last_percent = percent // 10
            if part.stat().st_size != meta["size"] or compute_md5(part) != expected[4:]:
                part.unlink(missing_ok=True)
                raise ValueError(f"size/checksum mismatch: {path.name}")
            part.replace(path)
            atomic_json(marker, {**fingerprint(path), "checksum": expected})
            return path
        except (requests.RequestException, ValueError):
            if attempt == 2:
                raise
            time.sleep(2 ** attempt)


def unpack_tables(archive, directory):
    marker = directory / "tables_ready.json"
    state = fingerprint(archive)
    destination = directory / "benchmark_data"
    if marker.is_file():
        saved = json.loads(marker.read_text())
        if saved.get("archive") == state:
            intact = all((directory / name).is_file() and fingerprint(directory / name) == meta
                         for name, meta in saved.get("files", {}).items())
            if intact and saved.get("files"):
                return destination
    with tempfile.TemporaryDirectory(dir=directory) as tmp:
        tmp = pathlib.Path(tmp)
        with zipfile.ZipFile(archive) as handle:
            for info in handle.infolist():
                path = pathlib.PurePosixPath(info.filename)
                if path.is_absolute() or ".." in path.parts or not path.parts or path.parts[0] != "benchmark_data":
                    raise ValueError(f"unexpected benchmark archive member: {info.filename}")
                if info.is_dir():
                    continue
                output = tmp.joinpath(*path.parts)
                output.parent.mkdir(parents=True, exist_ok=True)
                with handle.open(info) as source, output.open("wb") as target:
                    shutil.copyfileobj(source, target)
        for name in ("comparisons.tsv", "experiments.tsv"):
            if not (tmp / "benchmark_data" / name).is_file():
                raise ValueError(f"benchmark archive is missing {name}")
        if destination.exists():
            shutil.rmtree(destination)  # only our isolated record cache
        (tmp / "benchmark_data").replace(destination)
    atomic_json(marker, {"archive": state, "files": {
        str(p.relative_to(directory)): fingerprint(p) for p in destination.rglob("*") if p.is_file()}})
    return destination


@contextmanager
def prediction_stream(source):
    """Read a TSV or gzipped TSV directly from disk or a ZIP, without full extraction."""
    if source.get("member"):
        with zipfile.ZipFile(source["path"]) as archive, archive.open(source["member"]) as handle:
            if source["member"].endswith(".gz"):
                with gzip.GzipFile(fileobj=handle) as stream:
                    yield stream
            else:
                yield handle
    elif str(source["path"]).endswith(".gz"):
        with gzip.open(source["path"], "rb") as handle:
            yield handle
    else:
        with open(source["path"], "rb") as handle:
            yield handle


def prepare_inputs(record, cache, tools, *, repo, published_only=False, refresh=False, timeout=120):
    files = registry(record, cache, refresh=refresh, timeout=timeout)
    directory = cache / "zenodo" / str(record)
    if "benchmark_data.zip" not in files:
        raise ValueError(f"record {record} has no benchmark_data.zip")
    archive = download(files["benchmark_data.zip"], directory, timeout=timeout)
    tables = unpack_tables(archive, directory)
    sources = {}
    for tool in tools.itertuples(index=False):
        local = pathlib.Path(tool.predictor_output_path).expanduser()
        local = local if local.is_absolute() else repo / local
        if not published_only and local.is_file():
            logger.info("reusing local predictions for %s", tool.tool_id)
            sources[tool.tool_id] = {"path": str(local.resolve()), "origin": "local, not verified against Zenodo"}
            continue
        if tool.tool_id == "targetscan_family_expanded":
            name = "targetscan_family_expanded_standardized.tsv.gz"
            if name not in files:
                raise ValueError(f"record {record} has no {name}")
            path = download(files[name], directory, timeout=timeout)
            sources[tool.tool_id] = {"path": str(path), "record": str(record), "checksum": files[name]["checksum"]}
        elif tool.tool_id in PUBLISHED_PREDICTOR_FILENAMES:
            old = registry(PREDICTOR_RECORD, cache, refresh=refresh, timeout=timeout)
            name = "fun-miRBench_standardized_predictions.zip"
            path = download(old[name], cache / "zenodo" / PREDICTOR_RECORD, timeout=timeout)
            with zipfile.ZipFile(path) as handle:
                member = _find_archive_member(handle, PUBLISHED_PREDICTOR_FILENAMES[tool.tool_id])
            sources[tool.tool_id] = {"path": str(path), "member": member, "record": PREDICTOR_RECORD,
                                    "checksum": old[name]["checksum"]}
        else:
            raise ValueError(f"no local or published predictions for {tool.tool_id}")
    return tables, sources, {"benchmark_record": str(record), "benchmark_archive": files["benchmark_data.zip"],
                             "prediction_sources": sources}
