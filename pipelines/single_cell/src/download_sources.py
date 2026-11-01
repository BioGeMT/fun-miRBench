"""Download and unpack the source files of the selected single-cell builds (sources.tsv). A file whose size and MD5 match is
kept, a partial file resumes, and a finished unpack step is not repeated, so the script can be rerun after an interruption.

Unpack steps: untar:<dir>, untar_inner:<dir>:<name;name> (archives inside a folder already unpacked), unzip:<dir>,
seurat:<dir> (raw counts of a Seurat object, for the targets and controls of target_genes.tsv).
"""
from __future__ import annotations

import argparse
import hashlib
import pathlib
import subprocess
import sys
import tarfile

import pandas as pd

SINGLE_CELL = pathlib.Path(__file__).resolve().parents[1]
ADAPTERS = SINGLE_CELL / "adapters"
SEURAT_EXTRA_NAMES = ("NO-TARGET",)
RETRIES = 5


def md5sum(path: pathlib.Path) -> str:
    digest = hashlib.md5()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 22), b""):
            digest.update(block)
    return digest.hexdigest()


def fetch(row, path: pathlib.Path) -> str:
    """The file at path, checked. Returns its MD5."""
    path.parent.mkdir(parents=True, exist_ok=True)
    expected_bytes = int(row.bytes) if row.bytes else None
    if path.exists() and expected_bytes == path.stat().st_size:
        actual = md5sum(path)
        if not row.md5 or actual == row.md5:
            print(f"ok      {row.file}", flush=True)
            return actual
        path.unlink()
    part = path.with_name(path.name + ".part")
    for attempt in range(1, RETRIES + 1):
        done = subprocess.run(["curl", "-fsSL", "--retry", "5", "--retry-delay", "20", "-C", "-", "-o", str(part), row.url],
                              check=False)
        if done.returncode == 0 and (expected_bytes is None or part.stat().st_size == expected_bytes):
            break
        print(f"retry   {row.file}: attempt {attempt}, curl exit {done.returncode}", flush=True)
    else:
        raise RuntimeError(f"{row.file}: download failed after {RETRIES} attempts")
    actual = md5sum(part)
    if row.md5 and actual != row.md5:
        part.unlink()
        raise RuntimeError(f"{row.file}: MD5 {actual}, expected {row.md5}")
    part.replace(path)
    print(f"fetched {row.file}", flush=True)
    return actual


def safe_extract(archive: pathlib.Path, target: pathlib.Path) -> None:
    target.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive) as handle:
        handle.extractall(target, filter="data")


def seurat_counts(rds: pathlib.Path, target: pathlib.Path) -> None:
    table = pd.read_csv(ADAPTERS / "target_genes.tsv", sep="\t", dtype=str, keep_default_na=False)
    names = {*table.symbol, *SEURAT_EXTRA_NAMES}
    for aliases in table.aliases:
        names.update(name for name in aliases.split(";") if name)
    target.mkdir(parents=True, exist_ok=True)
    keep = target / "kept_values.txt"
    keep.write_text("\n".join(sorted(names)) + "\n")
    subprocess.run(["Rscript", str(ADAPTERS / "seurat_counts.R"), str(rds), "RNA", "gene", str(keep), str(target)], check=True)


def unpack(row, path: pathlib.Path, root: pathlib.Path) -> None:
    if not row.unpack:
        return
    marker = path.with_name(path.name + ".unpacked")
    if marker.exists():
        return
    for step in row.unpack.split(","):
        kind, directory, *rest = step.split(":")
        target = root / directory
        if kind == "untar":
            safe_extract(path, target)
        elif kind == "untar_inner":
            for name in rest[0].split(";"):
                safe_extract(target / name, target)
        elif kind == "unzip":
            subprocess.run(["unzip", "-o", "-q", str(path), "-d", str(target)], check=True)
        elif kind == "seurat":
            seurat_counts(path, target)
        else:
            raise ValueError(f"{row.file}: unknown unpack step {step}")
    marker.touch()
    print(f"unpacked {row.file}", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--data", type=pathlib.Path, required=True)
    parser.add_argument("--groups", default="main,machinery,supplement")
    parser.add_argument("--dataset", action="append", help="only these source datasets (a column of sources.tsv)")
    parser.add_argument("--checksums-out", type=pathlib.Path, help="write file, bytes and MD5 of what was fetched")
    args = parser.parse_args()
    table = pd.read_csv(SINGLE_CELL / "sources.tsv", sep="\t", dtype=str, keep_default_na=False)
    groups = set(args.groups.split(","))
    table = table[table.group.isin(groups)]
    if args.dataset:
        table = table[table.source_dataset.isin(args.dataset)]
    if table.empty:
        raise ValueError(f"no sources for groups {sorted(groups)} and datasets {args.dataset}")
    print(f"{len(table)} files, {sum(int(b) for b in table.bytes if b) / 1e9:.1f} GB listed with a size", flush=True)
    records, failed = [], []
    for count, row in enumerate(table.itertuples(index=False), start=1):
        path = args.data / row.file
        try:
            actual = fetch(row, path)
            unpack(row, path, args.data)
            records.append({"file": row.file, "bytes": path.stat().st_size, "md5": actual})
        except (RuntimeError, subprocess.CalledProcessError, tarfile.TarError, OSError) as error:
            failed.append(row.file)
            print(f"FAILED  {row.file}: {error}", flush=True)
        print(f"PROGRESS {count}/{len(table)} {row.file}", flush=True)
    if args.checksums_out:
        pd.DataFrame(records).to_csv(args.checksums_out, sep="\t", index=False)
    if failed:
        print(f"{len(failed)} files failed: {failed}", flush=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
