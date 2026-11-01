"""Human cardiac fibroblasts reprogrammed with MEF2C, GATA4, TBX5 and miR-133 (hMGT133, Zhou et al. 2019), 704 Fluidigm C1 wells on eight plates (GEO):
D0HCF uninfected day 0; D3M1, D3M2 hMGT133 day 3; D3UN DsRed + empty vector mixed 3:1 with uninfected, day 3; D3M3, D5M1, D7M1, D9M1
hMGT133 mixed 3:1 with tdTomato + empty vector. GEO holds only SCnorm tables, so counts come from our alignment. `prepare` writes the manifests;
`build` writes the store, using the authors' 652 cells (columns of the two GEO tables, saved from the xlsx) as QC.
Mixed plates: at least MARKER_POSITIVE reads on tdTomato is a control cell, below MARKER_NEGATIVE a treated cell, and cells between stay stored in their own arm, not compared. My rule, chosen from the marker reads alone (a gap in their distribution): it gives 14, 22, 18 and 19 control cells where GEO names 14, 21, 18 and 18.
D3UN: the authors' own table of 337 cells lists its 14 uninfected D3UN cells (GEO: 14 uninfected, 54 DsRed + empty vector, together the 68 of the other table)."""
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import pandas as pd
from tier2_plate import DATA, write_plate_store

ROOT = DATA / "GSE106888"
ALIGNED = ROOT / "aligned"
MARKER_POSITIVE, MARKER_NEGATIVE = 100, 10
UNINFECTED_TABLE = "GSE106888_SCnorm_337cells_d3579_d3UN.xlsx.norm_337cells_d3579_d3UN.columns.txt"
MIXED = {"D3M3": "d3", "D5M1": "d5", "D7M1": "d7", "D9M1": "d9"}
GEO_MARKER_POSITIVE = {"D3M3": 14, "D5M1": 21, "D7M1": 18, "D9M1": 18}  # GEO data processing: cells after QC


def report() -> pd.DataFrame:
    table = pd.read_csv(ROOT / "meta/SRP124967.tsv", sep="\t")
    table["cell"] = table.sample_title
    table["plate"] = table.cell.str.split("_").str[0]
    if table.cell.duplicated().any():
        raise ValueError("duplicated well names")
    return table


def prepare() -> None:
    ALIGNED.mkdir(exist_ok=True)
    table = report()
    fastq = ROOT / "fastq"
    (ALIGNED / "wells.manifest").write_text("".join(f"{fastq}/{r.run_accession}.fastq.gz\t-\t{r.cell}\n" for r in table.itertuples()))
    for plate in MIXED:
        rows = table[table.plate == plate]
        (ALIGNED / f"{plate}_tdtomato.manifest").write_text("".join(f"{fastq}/{r.run_accession}.fastq.gz\t{r.cell}\n" for r in rows.itertuples()))


def build(out: pathlib.Path) -> None:
    table = report().set_index("cell")
    kept = set()
    for path in (ROOT / "deposit").glob("*.columns.txt"):
        kept |= {line.strip() for line in path.read_text().splitlines() if line.strip() != "None"}
    if len(kept) != 652:
        raise ValueError(f"expected the 652 cells GEO names, found {len(kept)}")
    markers = pd.concat([pd.read_csv(ALIGNED / f"marker_{plate}/counts.tsv", sep="\t") for plate in MIXED]).set_index("cell")
    uninfected = {line.strip() for line in (ROOT / "deposit" / UNINFECTED_TABLE).read_text().splitlines() if line.startswith("D3UN")}
    if len(uninfected) != 14:
        raise ValueError(f"expected the 14 uninfected D3UN cells GEO names, found {len(uninfected)}")
    cells = pd.DataFrame({"sample": table.plate, "passed_qc": table.index.isin(kept)}, index=table.index)
    cells["marker_reads"] = markers.marker_reads.reindex(table.index)
    for plate, expected in GEO_MARKER_POSITIVE.items():
        print(f"{plate}: {int(((cells.marker_reads >= MARKER_POSITIVE) & cells.passed_qc & (cells['sample'] == plate)).sum())} tdTomato cells passing QC, GEO {expected}", flush=True)
    arm, control, batch = [], [], []
    for cell, row in zip(cells.index, cells.itertuples()):
        if row.sample == "D0HCF":
            arm.append("uninfected_d0"), control.append("uninfected"), batch.append("D0HCF")
        elif row.sample in ("D3M1", "D3M2"):
            arm.append("hMGT133_d3_pure_plates"), control.append(""), batch.append("d3_pure_plates")
        elif row.sample == "D3UN":
            if cell in uninfected:
                arm.append("uninfected_d3"), control.append("uninfected")
            elif row.passed_qc:
                arm.append("DsRed_EV_d3"), control.append("DsRed_empty_vector")
            else:
                arm.append("D3UN_not_kept_by_authors"), control.append("")
            batch.append("d3_pure_plates")
        else:
            day = MIXED[row.sample]
            if row.marker_reads >= MARKER_POSITIVE:
                arm.append(f"tdTomato_EV_{day}"), control.append("tdTomato_empty_vector")
            elif row.marker_reads < MARKER_NEGATIVE:
                arm.append(f"hMGT133_{day}"), control.append("")
            else:
                arm.append(f"tdTomato_between_{day}"), control.append("")
            batch.append(row.sample)
    cells["arm"], cells["control_type"], cells["batch"] = arm, control, batch
    write_plate_store("GSE106888_cardiac", ALIGNED / "wells/Solo.out/Gene/raw", cells, out,
                      "our STARsolo SmartSeq alignment against Ensembl 115, unique mappers, reads SRA SRP124967; cell QC: the authors' 652 cells (columns of the GEO SCnorm tables)")


if __name__ == "__main__":
    if sys.argv[1] == "prepare":
        prepare()
    else:
        build(pathlib.Path(sys.argv[2]))
