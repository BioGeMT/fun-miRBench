"""Tier 2 adapters: small RNA counts stay integers and leave out reads on several arms, and guide calls follow the authors' rule."""

import gzip
import pathlib
import sys

import pandas as pd
import pytest

ADAPTERS = pathlib.Path(__file__).resolve().parents[1] / "pipelines" / "single_cell" / "adapters"
sys.path.insert(0, str(ADAPTERS))

import HRA016755  # noqa: E402
import tier2_smallrna  # noqa: E402

needs_reference = pytest.mark.skipif(not (tier2_smallrna.REFERENCE / "miRNA.dat").exists() or not (tier2_smallrna.ENV / "bowtie").exists(),
                                     reason="miRBase files and the bowtie environment exist only on biogemt-4")
ADAPTER = "TGGAATTCTCGGGTGCCAAGG"
LET7A = "TGAGGTAGTAGGTTGTATAGTT"  # hsa-let-7a-5p, three loci, one name
MIR21 = "TAGCTTATCAGACTGATGTTGA"  # hsa-miR-21-5p


def write_fastq(path: pathlib.Path, sequences: list[str]) -> None:
    with gzip.open(path, "wt") as handle:
        for number, sequence in enumerate(sequences):
            handle.write(f"@r{number}\n{sequence}\n+\n{'I' * len(sequence)}\n")


@needs_reference
def test_smallrna_counts_reads_on_one_arm_and_library_total(tmp_path):
    tier2_smallrna.build_index()
    fastq = tmp_path / "cell.fastq.gz"
    too_short = "ACGTACGTACGT"
    write_fastq(fastq, [LET7A + ADAPTER] * 3 + [MIR21 + ADAPTER] * 2 + [too_short + ADAPTER] + [LET7A])  # the last has no adapter
    stats, counts = tier2_smallrna.count_file(str(fastq), "cell", tier2_smallrna.read_arms(), ADAPTER, 0, 0, 0)
    assert counts == {"hsa-let-7a-5p": 3, "hsa-miR-21-5p": 2}
    assert stats["reads_with_adapter"] == 6
    assert stats["inserts"] == 5


@needs_reference
def test_smallrna_umis_collapse_to_molecules(tmp_path):
    tier2_smallrna.build_index()
    fastq = tmp_path / "cell.fastq.gz"
    write_fastq(fastq, ["AAAAAAAA" + "CA" + MIR21 + ADAPTER] * 2 + ["CCCCCCCC" + "CA" + MIR21 + ADAPTER])
    stats, counts = tier2_smallrna.count_file(str(fastq), "cell", tier2_smallrna.read_arms(), ADAPTER, 8, 2, 0)
    assert counts == {"hsa-miR-21-5p": 2}
    assert stats["inserts"] == 2


def test_guide_calls_resolve_chimeras_then_follow_the_authors_shares(tmp_path, monkeypatch):
    monkeypatch.setattr(HRA016755, "ALIGNED", tmp_path)
    (tmp_path / "library").mkdir()
    rows = [("single", f"U{n}", "MIR155HG", 20) for n in range(40)] + [("single", "U0", "nontarget-1", 2),  # a chimera on U0 loses to the barcode with most reads
            ("single", "low", "NEAT1", 4),  # under the read minimum
            ("single", "tie", "GAS5", 9), ("single", "tie", "XIST", 9)]  # a tie is dropped
    rows += [("double", f"A{n}", "nontarget-1", 9) for n in range(5)] + [("double", f"B{n}", "nontarget-2", 9) for n in range(5)]
    rows += [("double_other", f"A{n}", "MIR155HG", 9) for n in range(5)] + [("double_other", f"B{n}", "nontarget-1", 9) for n in range(5)]
    rows += [("triple", f"{target}{n}", target, 9) for target in ("MIR222HG", "NEAT1", "GAS5") for n in range(5)]
    rows += [("unclear", f"{target}{n}", target, 9) for target in ("MIR155HG", "GAS5", "NEAT1", "XIST") for n in range(5)]
    pd.DataFrame(rows, columns=["cell", "umi", "target", "reads"]).to_csv(tmp_path / "library" / "gbc_umi_reads.tsv.gz", sep="\t", index=False)
    calls = HRA016755.calls("library")
    assert calls.molecules["single"] == 40
    assert calls.call.to_dict() == {"single": "single", "double": "double", "double_other": "double", "triple": "triple", "unclear": "unclear"}
    arms = {cell: HRA016755.arm_of(row) for cell, row in calls.iterrows()}
    assert arms == {"single": "MIR155HG", "double": "nontarget", "double_other": None, "triple": None, "unclear": None}
