"""Download the FASTQ files of an ENA run report over HTTP, resuming partial files, and check each against ENA's md5.
Writes checksums.tsv (file, md5, bytes) in the output folder and exits non-zero if a file stays wrong after retries.
Usage: tier2_fetch.py <report.tsv> <out_dir> <jobs> [runs.txt]; runs.txt limits the download to the runs it lists."""

import concurrent.futures
import csv
import hashlib
import os
import pathlib
import subprocess
import sys
import time

RETRIES = int(os.environ.get("TIER2_FETCH_RETRIES", 8))  # GSA-Human drops connections often and needs far more
WAIT_SECONDS = 120  # ENA answers 403 for a while after a burst of requests


def md5_of(path: pathlib.Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        while chunk := handle.read(1 << 22):
            digest.update(chunk)
    return digest.hexdigest()


def fetch(url: str, md5: str, size: int, out: pathlib.Path) -> tuple[str, str]:
    """Returns (file name, error); a file already complete and correct is not fetched again."""
    for attempt in range(RETRIES):
        if attempt:
            time.sleep(WAIT_SECONDS)
        if not out.exists() or out.stat().st_size < size:
            subprocess.run(["curl", "-sS", "-f", "-L", "--retry", "5", "--retry-delay", "30", "--speed-time", "300",
                            "--speed-limit", "10000", "-C", "-", "-o", str(out), f"http://{url}"], check=False)
        if out.exists() and out.stat().st_size == size and md5_of(out) == md5:
            return out.name, ""
        if out.exists() and out.stat().st_size < 4096 and out.read_bytes().startswith(b"<!DOCTYPE"):
            out.unlink()
            return out.name, "the server lists a directory at this path, not a file"  # seen on ENA: SRR28574011_2, SRR3031892
        if out.exists() and out.stat().st_size >= size:
            out.unlink()  # complete but wrong: resuming cannot repair it
    return out.name, f"size or md5 wrong after {RETRIES} attempts"


def main() -> None:
    report, out_dir, jobs = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2]), int(sys.argv[3])
    wanted = set(pathlib.Path(sys.argv[4]).read_text().split()) if len(sys.argv) > 4 else None
    out_dir.mkdir(parents=True, exist_ok=True)
    tasks = []
    with report.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if wanted is None or row["run_accession"] in wanted:
                for url, md5, size in zip(row["fastq_ftp"].split(";"), row["fastq_md5"].split(";"), row["fastq_bytes"].split(";")):
                    tasks.append((url, md5, int(size), out_dir / url.rsplit("/", 1)[1]))
    if not tasks:
        sys.exit(f"nothing to download from {report}")
    print(f"{len(tasks)} files, {sum(task[2] for task in tasks) / 1e9:.1f} GB", flush=True)
    failed, done = [], 0
    with concurrent.futures.ThreadPoolExecutor(jobs) as pool:
        for name, error in pool.map(lambda task: fetch(*task), tasks):
            done += 1
            print(f"PROGRESS {done}/{len(tasks)} {name} {error or 'ok'}", flush=True)
            if error:
                failed.append(name)
    with (out_dir / "checksums.tsv").open("w") as handle:
        handle.write("file\tmd5\tbytes\n")
        for _, md5, size, out in tasks:
            if out.name not in failed:
                handle.write(f"{out.name}\t{md5}\t{size}\n")
    if failed:
        sys.exit(f"failed: {failed}")


if __name__ == "__main__":
    main()
