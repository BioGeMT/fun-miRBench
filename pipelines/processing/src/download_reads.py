"""Download the FASTQ files of the 21 OE test comparisons.

Reads pipelines/processing/configs/fastq_manifest.tsv: every file is downloaded from ENA to its repository path
and checked against the manifest md5 and size. Existing files that already match their md5 are kept.

Usage:
  uv run python pipelines/processing/src/download_reads.py [--study GSE115646 ...] [--workers 4] [--attempts 5] [--dry-run]
"""
import argparse
import hashlib
import pathlib
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd
import requests

REPO = pathlib.Path(__file__).resolve().parents[3]
MANIFEST = REPO / "pipelines/processing/configs/fastq_manifest.tsv"


def md5(path: pathlib.Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 22), b""):
            h.update(block)
    return h.hexdigest()


class DownloadError(RuntimeError):
    """A transfer still failed after retries, or returned a permanent HTTP error."""


def fetch(row, *, attempts: int = 5) -> str:
    if attempts < 1:
        raise ValueError("attempts must be at least 1")
    dest = REPO / row.path
    if dest.exists() and dest.stat().st_size == row.expected_bytes and md5(dest) == row.expected_md5:
        return f"ok      {row.path}"
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_name(dest.name + ".part")
    for attempt in range(1, attempts + 1):
        try:
            with requests.get(row.source_url, stream=True, timeout=120) as r:
                r.raise_for_status()
                with open(part, "wb") as fh:
                    for chunk in r.iter_content(1 << 22):
                        fh.write(chunk)
            if part.stat().st_size != row.expected_bytes or md5(part) != row.expected_md5:
                raise ValueError("size or md5 mismatch")
            part.replace(dest)
            return f"fetched {row.path}"
        except (requests.RequestException, ValueError) as exc:
            status = exc.response.status_code if isinstance(exc, requests.HTTPError) and exc.response is not None else None
            retryable = status is None or status in (408, 429) or status >= 500
            reason = (f"HTTP {status}" if status is not None else
                      "size or md5 mismatch" if isinstance(exc, ValueError) else type(exc).__name__)
            if not retryable or attempt == attempts:
                raise DownloadError(f"{row.path}: {reason} after {attempt} attempt(s)") from exc
            delay = min(30, 2 ** min(attempt, 5))
            print(f"retry   {row.path}: {reason}; attempt {attempt + 1}/{attempts} in {delay}s",
                  file=sys.stderr, flush=True)
        finally:
            # Only validated complete files are promoted. A retry starts this file again.
            part.unlink(missing_ok=True)
        time.sleep(delay)


def positive_int(value: str) -> int:
    number = int(value)
    if number < 1:
        raise argparse.ArgumentTypeError("must be at least 1")
    return number


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--study", nargs="*", help="only these GEO series")
    ap.add_argument("--workers", type=positive_int, default=4)
    ap.add_argument("--attempts", type=positive_int, default=5,
                    help="attempts per file, including the first (default: 5); failed transfers restart that file")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    man = pd.read_csv(MANIFEST, sep="\t", keep_default_na=False)
    if args.study:
        man = man[man.study.isin(args.study)]
    print(f"{len(man)} files, {man.expected_bytes.sum() / 1e9:.1f} GB")
    if args.dry_run:
        return
    failed = []
    with ThreadPoolExecutor(args.workers) as pool:
        jobs = [pool.submit(fetch, row, attempts=args.attempts) for row in man.itertuples(index=False)]
        for job in as_completed(jobs):
            try:
                print(job.result(), flush=True)
            except DownloadError as exc:
                failed.append(str(exc))
                print(f"FAILED  {exc}", file=sys.stderr, flush=True)
    if failed:
        raise SystemExit(f"{len(failed)} download(s) failed. Rerun the same command to retry; verified files are kept.")


if __name__ == "__main__":
    main()
