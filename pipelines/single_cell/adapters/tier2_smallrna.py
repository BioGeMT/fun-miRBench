"""Mature miRNA counts per cell from raw small RNA reads (miRBase 22.1). Same recipe as ko_hypotheses/B_abundance
quantify_smallrna.sh, with one change: a read whose best hits fall on several differently named arms is not counted,
so counts stay integers.

Each file: trim the 3' adapter with cutadapt (reads without it are dropped), remove the UMI and random bases, keep inserts
of 18-30 nt, collapse (with --umi-5p, one molecule per UMI and insert), align to the human hairpins with bowtie 1 (at most
1 mismatch, best stratum, all best hits, sense strand). A read counts for an arm when it starts within 2 nt of the arm's
start and ends at most 5 nt after the arm's end (miRDeep2's windows).

Usage: tier2_smallrna.py <manifest: fastq, cell> <out_dir> <jobs> --adapter SEQ|none [--umi-5p N] [--random-5p N] [--random-3p N]
Writes mirna_counts.tsv.gz (cells x arms) and library_stats.tsv, whose `inserts` column is the library total.
Prints one PROGRESS line per file and exits non-zero if any file fails."""

import argparse
import collections
import concurrent.futures
import pathlib
import subprocess
import tempfile

import pandas as pd

ENV = pathlib.Path("/SCRATCH/dtzim01/tools/envs/smallrna/bin")
REFERENCE = pathlib.Path("/SCRATCH/dtzim01/single_cell_mirna_data/tier2_reference/mirbase_22.1")
INDEX = REFERENCE / "bowtie/hsa_hairpin"
MIN_INSERT, MAX_INSERT = 18, 30
UP, DOWN = 2, 5


def build_index() -> None:
    if (INDEX.parent / "hsa_hairpin.1.ebwt").exists():
        return
    INDEX.parent.mkdir(exist_ok=True)
    keep, lines = False, []
    for line in (REFERENCE / "hairpin.fa").read_text().splitlines():
        if line.startswith(">"):
            keep = line.startswith(">hsa-")
            if keep:
                lines.append(line.split()[0])
        elif keep:
            lines.append(line.replace("U", "T"))
    fasta = INDEX.parent / "hsa_hairpin.fa"
    fasta.write_text("\n".join(lines) + "\n")
    subprocess.run([ENV / "bowtie-build", "-q", fasta, INDEX], check=True)


def read_arms() -> dict[str, list[tuple[int, int, str]]]:
    """Human hairpin -> [(start, end, mature name)], 1-based inclusive."""
    arms, hairpin, position = collections.defaultdict(list), None, None
    for line in (REFERENCE / "miRNA.dat").read_text().splitlines():
        if line.startswith("ID"):
            name = line.split()[1]
            hairpin = name if name.startswith("hsa-") else None
        elif hairpin and line.startswith("FT   miRNA"):
            position = tuple(int(value) for value in line.split()[2].split(".."))
        elif hairpin and position and "/product=" in line:
            arms[hairpin].append((*position, line.split('"')[1]))
            position = None
    return dict(arms)


def count_file(fastq: str, cell: str, arms: dict, adapter: str, umi: int, random_5p: int, random_3p: int) -> tuple[dict, collections.Counter]:
    trim = [] if adapter == "none" else ["-a", adapter, "--discard-untrimmed"]  # none: reads deposited already trimmed
    cutadapt = subprocess.Popen([ENV / "cutadapt", "-j", "1", *trim, "--quiet", fastq], stdout=subprocess.PIPE, text=True)
    reads_with_adapter, molecules, reads_per_insert = 0, set(), collections.Counter()
    for number, line in enumerate(cutadapt.stdout):
        if number % 4 != 1:
            continue
        reads_with_adapter += 1
        sequence = line.rstrip("\n")
        insert = sequence[umi + random_5p:len(sequence) - random_3p]
        if MIN_INSERT <= len(insert) <= MAX_INSERT and "N" not in sequence[:umi]:
            if umi:
                molecules.add((sequence[:umi], insert))
            else:
                reads_per_insert[insert] += 1
    if cutadapt.wait():
        raise RuntimeError(f"{cell}: cutadapt exit {cutadapt.returncode}")
    if umi:
        reads_per_insert = collections.Counter(insert for _, insert in molecules)
    inserts = sum(reads_per_insert.values())
    names, aligned = collections.defaultdict(set), set()
    if inserts:
        ids = {str(i): insert for i, insert in enumerate(reads_per_insert)}
        with tempfile.NamedTemporaryFile("w", suffix=".fa") as fasta:
            fasta.writelines(f">{i}\n{insert}\n" for i, insert in ids.items())
            fasta.flush()
            hits = subprocess.run([ENV / "bowtie", "-p", "1", "-f", "-v", "1", "-a", "--best", "--strata", "--norc", INDEX, fasta.name],
                                  check=True, capture_output=True, text=True).stdout
        for line in hits.splitlines():
            read, _, hairpin, offset, sequence = line.split("\t")[:5]
            aligned.add(read)
            start = int(offset) + 1
            end = start + len(sequence) - 1
            for arm_start, arm_end, name in arms.get(hairpin, []):
                if arm_start - UP <= start <= arm_start + UP and end <= arm_end + DOWN:
                    names[read].add(name)
    counts = collections.Counter()
    ambiguous = 0
    for read, hit in names.items():
        amount = reads_per_insert[ids[read]]
        if len(hit) == 1:
            counts[next(iter(hit))] += amount
        else:
            ambiguous += amount
    stats = {"cell": cell, "reads_with_adapter": reads_with_adapter, "inserts": inserts,
             "aligned_to_hairpins": sum(reads_per_insert[ids[read]] for read in aligned),
             "assigned_to_one_arm": sum(counts.values()), "not_counted_several_arms": ambiguous}
    return stats, counts


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest", type=pathlib.Path)
    parser.add_argument("out_dir", type=pathlib.Path)
    parser.add_argument("jobs", type=int)
    parser.add_argument("--adapter", required=True)
    parser.add_argument("--umi-5p", type=int, default=0)
    parser.add_argument("--random-5p", type=int, default=0)
    parser.add_argument("--random-3p", type=int, default=0)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    build_index()
    arms = read_arms()
    files = pd.read_csv(args.manifest, sep="\t", names=["fastq", "cell"])
    if files.cell.duplicated().any():
        raise ValueError(f"duplicated cells in {args.manifest}")
    stats, counts, failed = [], {}, []
    with concurrent.futures.ThreadPoolExecutor(args.jobs) as pool:
        futures = {pool.submit(count_file, row.fastq, row.cell, arms, args.adapter, args.umi_5p, args.random_5p, args.random_3p): row.cell
                   for row in files.itertuples()}
        for done, future in enumerate(concurrent.futures.as_completed(futures), 1):
            cell = futures[future]
            try:
                one, counted = future.result()
                stats.append(one)
                counts[cell] = counted
                print(f"PROGRESS {done}/{len(futures)} {cell} inserts={one['inserts']} assigned={one['assigned_to_one_arm']}", flush=True)
            except (RuntimeError, subprocess.CalledProcessError, OSError) as error:
                failed.append(cell)
                print(f"PROGRESS {done}/{len(futures)} {cell} FAILED {error}", flush=True)
    if failed:
        raise SystemExit(f"failed: {failed}")
    order = files.cell.tolist()
    arms_seen = sorted({name for counted in counts.values() for name in counted})
    table = pd.DataFrame([[counts[cell].get(name, 0) for name in arms_seen] for cell in order], index=order, columns=arms_seen, dtype="int64")
    table.index.name = "cell"
    table.to_csv(args.out_dir / "mirna_counts.tsv.gz", sep="\t")
    pd.DataFrame(stats).set_index("cell").loc[order].to_csv(args.out_dir / "library_stats.tsv", sep="\t")


if __name__ == "__main__":
    main()
