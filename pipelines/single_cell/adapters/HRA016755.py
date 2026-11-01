"""Perturb-Multiome (CRISPRi, three sgRNAs per vector, one 18 nt guide barcode per vector) in K562 and BJ, GSA-Human HRA016755; RNA part only, ATAC not downloaded. MIR155HG and MIR222HG against non-targeting.
Counts: our STARsolo alignment of the GEX runs. Guide calls: reprocessed from the guide-barcode runs with the barcode table of Supplementary Table 3 (Nature Aging 2026): read 1 holds cell and UMI, read 2 holds
the barcode at position 30 after TAGCAAACTGGGGCACAAGCTTAATTAAG, at most 1 mismatch. Most cell and UMI pairs carry several barcodes (55 to 84% per library, PCR chimeras), so each pair keeps the barcode with most reads, ties dropped, and it counts as a molecule with at least MIN_MOLECULE_READS reads (my rule; the authors' scripts are marked illustrative and do not show this step).
Single calls by the shape of the authors' rule (author_code/2_GBC_cellID.R): second barcode under SINGLE_BELOW of a cell's molecules, fourth under it. That script writes 1%, the next one calls its table oneGBClist_0.05, and 5% is what we use.
Non-targeting cells: a single call on nontarget-1 or nontarget-2, or a double call on both, as the authors' 4_K562_process_integration.R. No response-based filter (the authors' filter on perturbed cells uses the response, so it stays out).
Cell QC: 600 < RNA counts < 130,000, the RNA part of the thresholds in the authors' K562 script, on a commented-out line there, and used for BJ as well. Usage: HRA016755.py gbc <library> | build K562|BJ <store path>"""
import collections
import gzip
import itertools
import pathlib
import subprocess
import sys

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import numpy as np
import pandas as pd
from _sample_arms import Library, build_store
from funmirbench.single_cell_store import read_10x_mtx
from tier2_plate import DATA

ROOT = DATA / "HRA016755"
ALIGNED = ROOT / "aligned"
LIBRARIES = {"K562_rep1": ("K562", "HRR3070151", "fq"), "K562_rep2": ("K562", "HRR3070155", "fq"), "BJ": ("BJ", "HRR3070161", "fq")}
GBC_START, GBC_LENGTH = 29, 18
TARGETS = {"MIR155HG": "MIR155HG", "MIR222HG": "MIR222HG", "nontarget-1": "nontarget", "nontarget-2": "nontarget"}
SINGLE_BELOW, MIN_MOLECULE_READS, MIN_RNA, MAX_RNA = 0.05, 5, 600, 130000


def barcodes() -> dict[str, str]:
    table = pd.read_csv(ROOT / "meta/gbc_table.tsv", sep="\t")
    exact = dict(zip(table.gbc, table.target))
    neighbours: dict[str, set[str]] = collections.defaultdict(set)
    for gbc, target in exact.items():
        for position, base in itertools.product(range(GBC_LENGTH), "ACGTN"):
            neighbours[gbc[:position] + base + gbc[position + 1:]].add(target)
    return {seq: next(iter(names)) for seq, names in neighbours.items() if len(names) == 1}


def gbc(library: str) -> None:
    """Reads per cell, UMI and barcode, for the cells of the library's alignment."""
    run = LIBRARIES[library][1]
    cells = set(pd.read_csv(ALIGNED / library / "Solo.out/Gene/filtered/barcodes.tsv", header=None)[0])
    lookup = barcodes()
    r1 = subprocess.Popen(["pigz", "-dc", ROOT / f"fastq/{run}_f1.fq.gz"], stdout=subprocess.PIPE, text=True)
    r2 = subprocess.Popen(["pigz", "-dc", ROOT / f"fastq/{run}_r2.fq.gz"], stdout=subprocess.PIPE, text=True)
    counts, reads = collections.Counter(), 0
    for number, (line1, line2) in enumerate(zip(r1.stdout, r2.stdout)):
        if number % 4 != 1:
            continue
        reads += 1
        target = lookup.get(line2[GBC_START:GBC_START + GBC_LENGTH])
        if target and line1[:16] in cells:
            counts[(line1[:16], line1[16:28], target)] += 1
    if r1.wait() or r2.wait():
        raise RuntimeError(f"{library}: reading the guide-barcode run failed")
    with gzip.open(ALIGNED / library / "gbc_umi_reads.tsv.gz", "wt", compresslevel=1) as handle:
        handle.write("cell\tumi\ttarget\treads\n")
        handle.writelines(f"{cell}\t{umi}\t{target}\t{value}\n" for (cell, umi, target), value in counts.items())
    print(f"{library}: {reads} reads, {sum(counts.values())} with a known barcode in a cell, {len(counts)} cell, UMI and barcode sets", flush=True)


def molecules(library: str) -> pd.DataFrame:
    """Molecules per cell and barcode: per cell and UMI the barcode with most reads, ties dropped, at least MIN_MOLECULE_READS reads."""
    table = pd.read_csv(ALIGNED / library / "gbc_umi_reads.tsv.gz", sep="\t")
    tied = table.reads == table.groupby(["cell", "umi"]).reads.transform("max")
    top = table[tied]
    top = top[~top.duplicated(["cell", "umi"], keep=False)]
    top = top[top.reads >= MIN_MOLECULE_READS]
    return top.groupby(["cell", "target"]).size().rename("molecules").reset_index()


def calls(library: str) -> pd.DataFrame:
    """One row per cell with a barcode: call (single, double, triple or unclear), the two most common barcodes, n_guides."""
    table = molecules(library).sort_values(["cell", "molecules"], ascending=[True, False])
    rows = []
    for cell, group in table.groupby("cell", sort=False):
        share = group.molecules.to_numpy() / group.molecules.sum()
        second, third, fourth = (share[index] if len(share) > index else 0 for index in (1, 2, 3))
        call = "unclear" if fourth >= SINGLE_BELOW else "single" if second < SINGLE_BELOW else "double" if third < SINGLE_BELOW else "triple"
        rows.append((cell, call, group.target.iloc[0], group.target.iloc[1] if len(group) > 1 else "", {"single": 1, "double": 2, "triple": 3}.get(call, len(group)), int(group.molecules.sum())))
    return pd.DataFrame(rows, columns=["cell", "call", "top_barcode", "next_barcode", "n_guides", "molecules"]).set_index("cell")


def arm_of(row) -> str | None:
    if row.call == "single":
        return TARGETS.get(row.top_barcode)
    if row.call == "double" and {row.top_barcode, row.next_barcode} == {"nontarget-1", "nontarget-2"}:
        return "nontarget"
    return None


def build(line: str, out: pathlib.Path) -> None:
    libraries = []
    for library, (library_line, _, _) in LIBRARIES.items():
        if library_line != line:
            continue
        folder = ALIGNED / library / "Solo.out/Gene/filtered"
        counts, features, names = read_10x_mtx(folder / "matrix.mtx", folder / "features.tsv", folder / "barcodes.tsv")
        called = calls(library)
        called.to_csv(ALIGNED / library / "gbc_calls.tsv.gz", sep="\t")
        print(f"{library}: {len(names)} cells, {len(called)} with a barcode; calls {called.call.value_counts().to_dict()}; median molecules {called.molecules.median()}", flush=True)
        arm = [arm_of(called.loc[name]) if name in called.index else None for name in names]
        print(f"{library}: arms {pd.Series(arm).value_counts().to_dict()}", flush=True)
        control = ["non_targeting" if value == "nontarget" else "" for value in arm]
        unit = None if line == "BJ" else library
        libraries.append(Library(library, line, arm, counts, features, names, control, unit=unit, unit_type="biological"))
    build_store(f"HRA016755_{line}", libraries, qc={"min_counts": MIN_RNA, "max_counts": MAX_RNA}, qc_source="600 to 130,000 RNA counts, from a commented-out line of the authors' K562 script, also used for BJ; guide calls by the shape of the authors' rule at 5% on our molecule counts after resolving chimeras, no response-based filter",
                source="our STARsolo 2.7.11b alignment of GSA-Human HRA016755 GEX runs against Ensembl 115, unique mappers; guide barcodes from the GBC runs and Supplementary Table 3", out=out)


if __name__ == "__main__":
    gbc(sys.argv[2]) if sys.argv[1] == "gbc" else build(sys.argv[2], pathlib.Path(sys.argv[3]))
