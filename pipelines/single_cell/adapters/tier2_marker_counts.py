"""Reads of a plate cell on a marker transgene (GFP, DsRed, tdTomato): bowtie 1, either strand, at most --mismatches, first hit only.
Usage: tier2_marker_counts.py <manifest: fastq, cell> <marker.fa> <out.tsv> <jobs> [--mismatches N] [--min-length N]. Writes cell, marker_reads, total_reads.
--min-length drops shorter reads first: adapter-trimmed deposits hold reads so short that they align anywhere with mismatches allowed."""
import argparse
import concurrent.futures
import pathlib
import subprocess

import pandas as pd

ENV = pathlib.Path("/SCRATCH/dtzim01/tools/envs/smallrna/bin")


def count(fastq: str, cell: str, index: pathlib.Path, mismatches: int, min_length: int) -> dict:
    keep = f"paste - - - - | awk -F'\\t' 'length($2) >= {min_length}' | tr '\\t' '\\n'"
    command = f"zcat '{fastq}' | {keep} | '{ENV}/bowtie' -p 1 -v {mismatches} -k 1 --suppress 1,2,3,4,5,6,7,8 '{index}' - 2>&1 >/dev/null"
    result = subprocess.run(["bash", "-c", f"set -o pipefail; {command}"], capture_output=True, text=True, check=True)
    fields = {line.split(":")[0].strip("# ").strip(): line.split(":")[1].split()[0] for line in result.stdout.splitlines() if ":" in line}
    return {"cell": cell, "marker_reads": int(fields["reads with at least one alignment"]), "total_reads": int(fields["reads processed"])}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest", type=pathlib.Path)
    parser.add_argument("marker", type=pathlib.Path)
    parser.add_argument("out", type=pathlib.Path)
    parser.add_argument("jobs", type=int)
    parser.add_argument("--mismatches", type=int, default=2)
    parser.add_argument("--min-length", type=int, default=0)
    args = parser.parse_args()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    index = args.out.parent / f"{args.marker.stem}_bowtie"
    subprocess.run([ENV / "bowtie-build", "-q", args.marker, index], check=True)
    files = pd.read_csv(args.manifest, sep="\t", names=["fastq", "cell"])
    with concurrent.futures.ThreadPoolExecutor(args.jobs) as pool:
        futures = [pool.submit(count, row.fastq, row.cell, index, args.mismatches, args.min_length) for row in files.itertuples()]
        rows = []
        for done, future in enumerate(futures, 1):
            rows.append(future.result())
            if done % 50 == 0:
                print(f"PROGRESS {done}/{len(futures)} markers", flush=True)
    pd.DataFrame(rows).to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    main()
