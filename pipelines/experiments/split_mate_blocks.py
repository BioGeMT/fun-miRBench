"""Split an ENA FASTQ that holds both mates as two positional blocks into _1 and _2 files.

Some ENA exports (for example GSE199332, BGI paired-end) put every mate-1 record first and every mate-2 record
second in one file, with SRA spot names that differ across mates. Equal block sizes alone do not prove the
pairing, so the split is accepted only if exact sequence overlaps between paired mates are far more frequent
than between shuffled pairs. The source file is never modified. A JSON audit is written next to the outputs.

Usage:
  uv run python pipelines/experiments/split_mate_blocks.py --input <run>.fastq.gz --run <SRR...> --out-dir <dir>
Needs pigz on PATH.
"""
import argparse
import hashlib
import json
import os
import pathlib
import subprocess

TRANSLATE = bytes.maketrans(b"ACGTNacgtn", b"TGCANtgcan")


def records(path: pathlib.Path):
    process = subprocess.Popen(["pigz", "-dc", "-p", "2", str(path)], stdout=subprocess.PIPE)
    try:
        stream = process.stdout
        while header := stream.readline():
            sequence, plus, quality = (stream.readline() for _ in range(3))
            if not header.startswith(b"@") or not plus.startswith(b"+"):
                raise ValueError(f"Malformed FASTQ in {path}")
            if not quality.endswith(b"\n") or len(sequence.rstrip()) != len(quality.rstrip()):
                raise ValueError(f"Truncated or unequal sequence/quality in {path}")
            yield header.rstrip(), sequence, quality
    finally:
        process.wait()
        process.stdout.close()


def overlaps(first: bytes, second: bytes) -> bool:
    """Exact >= 30-nt overlap between mate 1 and the reverse complement of mate 2."""
    second = second.translate(TRANSLATE)[::-1]
    position = second.find(first[-30:])
    if position < 0:
        return False
    length = position + 30
    return length <= len(first) and first[-length:] == second[:length]


def split(source: pathlib.Path, run: str, out_dir: pathlib.Path) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = {mate: out_dir / f"{run}_{mate}.fastq.gz" for mate in (1, 2)}
    counts, sampled, seen_blocks = {1: 0, 2: 0}, {1: [], 2: []}, []
    digests = {1: hashlib.sha256(), 2: hashlib.sha256()}
    handles = {mate: open(f"{p}.part", "wb") for mate, p in paths.items()}
    compressors = {mate: subprocess.Popen(["pigz", "-1", "-p", "2", "-c"], stdin=subprocess.PIPE, stdout=handle)
                   for mate, handle in handles.items()}
    buffers = {1: bytearray(), 2: bytearray()}
    try:
        for index, (header, sequence, quality) in enumerate(records(source), 1):
            tokens = header.split()
            if tokens[0] != f"@{run}.{index}".encode() or len(tokens) != 2:
                raise ValueError(f"Unexpected SRA record identity: {header!r}")
            if not tokens[1].endswith((b"/1", b"/2")):
                raise ValueError(f"Missing mate marker: {header!r}")
            mate = int(tokens[1][-1:])
            if not seen_blocks or seen_blocks[-1] != mate:
                seen_blocks.append(mate)
                if len(seen_blocks) > 2:
                    raise ValueError(f"Not two mate blocks in {run}")
            counts[mate] += 1
            digests[mate].update(sequence)
            digests[mate].update(quality)
            if counts[mate] <= 10000:
                sampled[mate].append(sequence.rstrip())
            buffers[mate].extend(f"@{run}.{counts[mate]}/{mate}\n".encode() + sequence + b"+\n" + quality)
            if len(buffers[mate]) >= 1 << 20:
                compressors[mate].stdin.write(buffers[mate])
                buffers[mate].clear()
        if counts[1] == 0 or counts[1] != counts[2] or set(seen_blocks) != {1, 2}:
            raise ValueError(f"Unequal or missing mate blocks in {run}: {counts}")
        n = min(len(sampled[1]), len(sampled[2]))
        matched = sum(overlaps(a, b) for a, b in zip(sampled[1], sampled[2]))
        shuffled = sum(overlaps(a, sampled[2][(i + 137) % n]) for i, a in enumerate(sampled[1]))
        if matched < 100 or matched <= 20 * (shuffled + 1):
            raise ValueError(f"Positional pairing lacks overlap support: {run}, {matched=}, {shuffled=}")
        for mate in (1, 2):
            compressors[mate].stdin.write(buffers[mate])
            compressors[mate].stdin.close()
            if compressors[mate].wait() != 0:
                raise ValueError(f"Compression failed for {run} mate {mate}")
            handles[mate].close()
        outputs = []
        for mate in (1, 2):
            part = pathlib.Path(f"{paths[mate]}.part")
            with open(part, "rb") as fh:
                sha = hashlib.file_digest(fh, "sha256").hexdigest()
            outputs.append({"mate": mate, "path": str(paths[mate]), "records": counts[mate],
                            "bytes": part.stat().st_size, "sha256": sha,
                            "sequence_quality_sha256": digests[mate].hexdigest()})
            os.replace(part, paths[mate])
        report = {"status": "validated", "run_accession": run, "source": str(source), "block_order": seen_blocks,
                  "sampled_pairs": n, "overlapping_pairs": matched, "shuffled_overlapping_pairs": shuffled,
                  "outputs": outputs}
        (out_dir / f"{run}.audit.json").write_text(json.dumps(report, indent=2) + "\n")
        return report
    finally:
        for process in compressors.values():
            if process.poll() is None:
                process.terminate()
            process.wait()
        for handle in handles.values():
            handle.close()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--input", type=pathlib.Path, required=True)
    ap.add_argument("--run", required=True)
    ap.add_argument("--out-dir", type=pathlib.Path, required=True)
    args = ap.parse_args()
    report = split(args.input, args.run, args.out_dir)
    print(f"{args.run}: {report['outputs'][0]['records']} pairs, {report['overlapping_pairs']} overlapping "
          f"vs {report['shuffled_overlapping_pairs']} shuffled")


if __name__ == "__main__":
    main()
