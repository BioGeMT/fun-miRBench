"""Utility functions for TargetScan standardization."""

from __future__ import annotations

import csv
import gzip as gz
import logging
import pathlib
import shutil
import urllib.request
import zipfile
from collections import Counter
from dataclasses import dataclass
from typing import Any, Dict, List, Tuple

from common import display_path

logger = logging.getLogger(__name__)

MIRBASE_RELEASE = "22.1"
MIRBASE_MATURE_URL = "https://mirbase.org/download_version_files/22.1/mature.fa"
PREDICTOR_NAME = "targetscan"
TARGETSCAN_SCORE_COL = "Cumulative weighted context++ score"


def repo_root() -> pathlib.Path:
    return pathlib.Path(__file__).resolve().parents[3]


def strip_version(value: str) -> str:
    s = value.strip()
    return s.split(".", 1)[0] if "." in s else s


def to_float_or_none(value: str):
    value = (value or "").strip()
    if value == "" or value.upper() == "NULL":
        return None
    return float(value)


def download_url(url: str, dest: pathlib.Path) -> None:
    tmp = dest.with_suffix(dest.suffix + ".part")
    try:
        with urllib.request.urlopen(url) as response, tmp.open("wb") as out:
            shutil.copyfileobj(response, out)
        tmp.replace(dest)
    finally:
        try:
            tmp.unlink()
        except FileNotFoundError:
            pass


def assert_fasta(path: pathlib.Path) -> None:
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            s = line.strip()
            if not s:
                continue
            if not s.startswith(">"):
                raise RuntimeError(f"Not FASTA (first line: {s[:120]})")
            return
    raise RuntimeError("FASTA check failed: file appears empty.")


def step1_download_targetscan_files(
    data_dir: pathlib.Path,
    resources_dir: pathlib.Path,
    *,
    force: bool = False,
) -> Dict[str, pathlib.Path]:
    data_dir = pathlib.Path(data_dir).resolve()
    data_dir.mkdir(parents=True, exist_ok=True)
    resources_dir = pathlib.Path(resources_dir).resolve()
    resources_dir.mkdir(parents=True, exist_ok=True)

    base = "https://www.targetscan.org/vert_80/vert_80_data_download"
    zip_destinations = {
        "Summary_Counts.all_predictions.txt.zip": data_dir,
        "miR_Family_Info.txt.zip": resources_dir,
        "Gene_info.txt.zip": resources_dir,
    }

    for fname, destination_dir in zip_destinations.items():
        url = f"{base}/{fname}"
        zip_path = destination_dir / fname

        if not zip_path.exists() or force:
            logger.info("Downloading %s", url)
            download_url(url, zip_path)
            logger.info("Downloaded %s", zip_path.name)
        else:
            logger.info("Skipping %s (already exists)", zip_path.name)

        with zipfile.ZipFile(zip_path, "r") as zf:
            for member in zf.namelist():
                extracted = destination_dir / member
                if extracted.exists() and not force:
                    continue
                zf.extract(member, destination_dir)
        logger.info("Unzipped %s", zip_path.name)

    expected = {
        "Summary_Counts.all_predictions.txt": data_dir / "Summary_Counts.all_predictions.txt",
        "miR_Family_Info.txt": resources_dir / "miR_Family_Info.txt",
        "Gene_info.txt": resources_dir / "Gene_info.txt",
    }

    logger.info("✔ TargetScan prediction data ready in %s", display_path(data_dir))
    logger.info("✔ TargetScan supporting resources ready in %s", display_path(resources_dir))
    return expected


def step2_build_representative_transcript_index(
    gene_info_path: pathlib.Path,
    *,
    species_id: str = "9606",
    report_top_n: int = 6,
) -> Dict[str, Any]:
    rep_txs_by_gene_id: Dict[str, List[Tuple[str, int, str]]] = {}
    gene_id_by_tx: Dict[str, str] = {}
    targetscan_gene_id_by_tx: Dict[str, str] = {}
    gene_symbol_by_gene_id: Dict[str, str] = {}

    tag_count_dist: Counter[int] = Counter()
    tx_count_by_gene: Counter[str] = Counter()
    rep_count_by_gene: Counter[str] = Counter()

    n_total = 0
    n_species = 0
    n_rep_rows = 0

    with pathlib.Path(gene_info_path).open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        for row in reader:
            n_total += 1
            if row["Species ID"].strip() != species_id:
                continue
            n_species += 1

            tx = row["Transcript ID"].strip()
            gene_id = row["Gene ID"].strip()
            gene_id_stable = strip_version(gene_id)
            gene_sym = row["Gene symbol"].strip()
            rep_flag = row["Representative transcript?"].strip()
            tags = int(row["3P-seq tags"].strip())

            gene_id_by_tx[tx] = gene_id
            targetscan_gene_id_by_tx[tx] = gene_id_stable
            gene_symbol_by_gene_id[gene_id_stable] = gene_sym
            tx_count_by_gene[gene_id_stable] += 1

            if rep_flag != "1":
                continue

            n_rep_rows += 1
            rep_count_by_gene[gene_id_stable] += 1
            tag_count_dist[tags] += 1

            rep_txs_by_gene_id.setdefault(gene_id_stable, []).append((tx, tags, gene_sym))

    genes_with_rep = len(rep_txs_by_gene_id)
    genes_with_multi_rep = sum(1 for count in rep_count_by_gene.values() if count > 1)

    logger.info(
        "Gene_Info stats: rows_total=%d | rows_species(%s)=%d | representative_rows=%d | genes_with_representative_tx=%d",
        n_total, species_id, n_species, n_rep_rows, genes_with_rep
    )

    top_bins = Counter(tx_count_by_gene.values()).most_common(report_top_n)
    logger.info(
        "Gene_Info transcripts-per-gene (top %d bins): %s",
        report_top_n,
        " | ".join([f"{k}tx:{v}" for k, v in top_bins]),
    )

    logger.info(
        "Representative-transcript QC: genes_with_multiple_representative_rows=%d",
        genes_with_multi_rep
    )
    if genes_with_multi_rep:
        logger.info("Genes with >1 representative transcript were found; all representative transcripts will be kept.")
        multi_rep_examples = 0
        for gene_id, reps in rep_txs_by_gene_id.items():
            if len(reps) <= 1:
                continue
            logger.info(
                "  multi-representative gene: %s (%s transcripts) -> %s",
                gene_id,
                len(reps),
                ", ".join(tx for tx, _tags, _sym in reps[:report_top_n]),
            )
            multi_rep_examples += 1
            if multi_rep_examples >= report_top_n:
                break

    top_tags = tag_count_dist.most_common(report_top_n)
    logger.info(
        "Representative transcript 3P-seq tag counts (top %d bins): %s",
        report_top_n,
        " | ".join([f"{k}:{v}" for k, v in top_tags]),
    )

    for gene_id, reps in list(rep_txs_by_gene_id.items())[:3]:
        tx, tags, sym = reps[0]
        logger.info("  sample representative transcript: %s -> %s (3P-seq tags=%d, symbol=%s)", gene_id, tx, tags, sym)

    return {
        "rep_txs_by_gene_id": rep_txs_by_gene_id,
        "gene_id_by_tx": gene_id_by_tx,
        "targetscan_gene_id_by_tx": targetscan_gene_id_by_tx,
        "gene_symbol_by_gene_id": gene_symbol_by_gene_id,
    }


def step3_download_ensembl115_gtf(
    destination: pathlib.Path,
    *,
    force: bool = False,
) -> pathlib.Path:
    url = "https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz"
    dest = pathlib.Path(destination).resolve()
    dest.parent.mkdir(parents=True, exist_ok=True)

    if dest.exists() and not force:
        logger.info("Skipping %s (already exists)", dest.name)
        return dest

    logger.info("Downloading %s", url)
    download_url(url, dest)
    logger.info("Downloaded %s", dest.name)
    return dest


def step4_build_and_cache_ensembl115_tables(
    ensembl_gtf_gz: pathlib.Path,
    *,
    cache_dir: pathlib.Path,
    force_rebuild: bool = False,
) -> Dict[str, Dict[str, str]]:
    cache_dir = pathlib.Path(cache_dir).resolve()
    cache_dir.mkdir(parents=True, exist_ok=True)

    tx2gene_path = cache_dir / "ensembl115_tx2gene.tsv.gz"
    gene2name_path = cache_dir / "ensembl115_gene2name.tsv.gz"

    if tx2gene_path.exists() and gene2name_path.exists() and not force_rebuild:
        tx_to_gene: Dict[str, str] = {}
        gene_to_name: Dict[str, str] = {}

        with gz.open(tx2gene_path, "rt", encoding="utf-8") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                tx_to_gene[row["transcript_id"]] = row["gene_id"]

        with gz.open(gene2name_path, "rt", encoding="utf-8") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                gene_to_name[row["gene_id"]] = row["gene_name"]

        logger.info(
            "Loaded cached Ensembl tables: tx_to_gene=%d | gene_to_name=%d",
            len(tx_to_gene), len(gene_to_name)
        )
        return {"tx_to_gene": tx_to_gene, "gene_to_name": gene_to_name}

    tx_to_gene: Dict[str, str] = {}
    gene_to_name: Dict[str, str] = {}

    n_total = 0
    n_transcript_rows = 0

    with gz.open(ensembl_gtf_gz, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            n_total += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "transcript":
                continue
            n_transcript_rows += 1

            attrs = {}
            for part in [item.strip() for item in fields[8].strip().split(";") if item.strip()]:
                if " " in part:
                    key, value = part.split(" ", 1)
                    attrs[key] = value.strip().strip('"')

            gene_id = attrs.get("gene_id")
            tx_id = attrs.get("transcript_id")
            gene_name = attrs.get("gene_name")

            if gene_id and tx_id:
                gene_id_s = strip_version(gene_id)
                tx_id_s = strip_version(tx_id)
                tx_to_gene[tx_id_s] = gene_id_s
                if gene_name:
                    gene_to_name.setdefault(gene_id_s, gene_name)

    logger.info(
        "Ensembl GTF parsed: lines=%d | transcript_features=%d | tx_mapped=%d | genes_named=%d",
        n_total, n_transcript_rows, len(tx_to_gene), len(gene_to_name)
    )

    with gz.open(tx2gene_path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["transcript_id", "gene_id"], delimiter="\t")
        writer.writeheader()
        for tx_id, gene_id in tx_to_gene.items():
            writer.writerow({"transcript_id": tx_id, "gene_id": gene_id})

    with gz.open(gene2name_path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["gene_id", "gene_name"], delimiter="\t")
        writer.writeheader()
        for gene_id, gene_name in gene_to_name.items():
            writer.writerow({"gene_id": gene_id, "gene_name": gene_name})

    logger.info("Wrote cached Ensembl tables:")
    logger.info("  - %s", display_path(tx2gene_path))
    logger.info("  - %s", display_path(gene2name_path))

    return {"tx_to_gene": tx_to_gene, "gene_to_name": gene_to_name}


def step5_qc_targetscan_vs_ensembl_transcripts(
    tx_index: Dict[str, Any],
    ensembl_tables: Dict[str, Dict[str, str]],
    *,
    report_n: int = 10,
) -> None:
    gene_id_by_tx = tx_index["gene_id_by_tx"]
    targetscan_gene_id_by_tx = tx_index["targetscan_gene_id_by_tx"]
    rep_txs_by_gene_id = tx_index["rep_txs_by_gene_id"]
    tx_to_gene_ensembl = ensembl_tables["tx_to_gene"]

    ts_all_tx = list(gene_id_by_tx.keys())
    ts_all_tx_stripped = [strip_version(value) for value in ts_all_tx]
    stripped_hits = sum(1 for value in ts_all_tx_stripped if value in tx_to_gene_ensembl)

    logger.info("TargetScan Gene_Info transcripts total: %d", len(ts_all_tx))
    logger.info("Ensembl v115 transcripts indexed: %d", len(tx_to_gene_ensembl))
    logger.info("Overlap (stable transcript ID): %d", stripped_hits)

    ts_rep_tx = [tx for reps in rep_txs_by_gene_id.values() for (tx, _tags, _sym) in reps]
    ts_rep_tx_stripped = [strip_version(value) for value in ts_rep_tx]
    rep_hits = sum(1 for value in ts_rep_tx_stripped if value in tx_to_gene_ensembl)

    logger.info("TargetScan representative transcripts: %d", len(ts_rep_tx))
    logger.info("Representative-transcript overlap with Ensembl (stable transcript ID): %d", rep_hits)

    rep_gene_matches = 0
    rep_gene_mismatches = 0
    mismatch_examples = []
    for tx_raw in ts_rep_tx:
        tx_stable = strip_version(tx_raw)
        ensembl_gene_id = tx_to_gene_ensembl.get(tx_stable)
        if ensembl_gene_id is None:
            continue
        targetscan_gene_id = targetscan_gene_id_by_tx.get(tx_raw)
        if targetscan_gene_id == ensembl_gene_id:
            rep_gene_matches += 1
        else:
            rep_gene_mismatches += 1
            if len(mismatch_examples) < report_n:
                mismatch_examples.append((tx_raw, targetscan_gene_id or "NA", ensembl_gene_id))

    logger.info(
        "Representative transcript gene-ID agreement vs Ensembl: match=%d | mismatch=%d",
        rep_gene_matches,
        rep_gene_mismatches,
    )
    if mismatch_examples:
        logger.info("Sample representative transcript gene-ID mismatches (tx -> TargetScan gene -> Ensembl gene):")
        for tx_raw, targetscan_gene_id, ensembl_gene_id in mismatch_examples:
            logger.info("  %s -> %s -> %s", tx_raw, targetscan_gene_id, ensembl_gene_id)

    misses = []
    for raw, stripped in zip(ts_rep_tx, ts_rep_tx_stripped):
        if stripped not in tx_to_gene_ensembl:
            misses.append((raw, stripped))
            if len(misses) >= report_n:
                break

    if misses:
        logger.info("Sample representative-transcript misses (raw -> stripped):")
        for raw, stripped in misses:
            logger.info("  %s -> %s", raw, stripped)

    non_enst = [value for value in ts_rep_tx_stripped if not value.startswith("ENST")]
    if non_enst:
        logger.info("Representative transcripts not ENST*: %d", len(non_enst))
        for value in non_enst[: min(report_n, len(non_enst))]:
            logger.info("  non-ENST example: %s", value)


def step6_download_mirbase_mature(
    destination: pathlib.Path,
    *,
    force: bool = False,
) -> pathlib.Path:
    dest = pathlib.Path(destination).resolve()
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and not force:
        logger.info("Skipping miRBase mature.fa (already exists): %s", display_path(dest))
        return dest

    logger.info("Downloading %s", MIRBASE_MATURE_URL)
    download_url(url=MIRBASE_MATURE_URL, dest=dest)
    assert_fasta(dest)
    logger.info("Downloaded miRBase mature.fa -> %s", display_path(dest))
    return dest


def parse_mirbase_mature(mature_fa: pathlib.Path) -> Dict[str, str]:
    logger.info("Parsing miRBase mature.fa")

    acc2name: Dict[str, str] = {}
    with pathlib.Path(mature_fa).open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if line.startswith(">"):
                parts = line[1:].strip().split()
                mir_name = parts[0]
                mir_acc = parts[1]
                acc2name[mir_acc] = mir_name

    logger.info("Parsed %d mature miRNAs from miRBase", len(acc2name))
    return acc2name


@dataclass(frozen=True)
class MirnaEntry:
    mirna_id: str
    mirna_name: str


def step_build_human_mirna_annotations(
    mir_family_info_path: pathlib.Path,
    *,
    mirbase_acc2name: Dict[str, str],
    species_id: str = "9606",
) -> Dict[str, MirnaEntry]:
    annotations: Dict[str, MirnaEntry] = {}

    n_total = 0
    n_species = 0
    n_missing_fields = 0
    n_acc_missing_in_mirbase = 0
    n_name_match = 0
    n_name_replaced = 0
    n_kept = 0
    n_duplicate_names = 0

    with pathlib.Path(mir_family_info_path).open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        for row in reader:
            n_total += 1
            if row["Species ID"].strip() != species_id:
                continue
            n_species += 1

            seed = row["Seed+m8"].strip()
            mir_name_ts = row["MiRBase ID"].strip()
            mir_acc = row["MiRBase Accession"].strip()

            if not seed or not mir_name_ts or not mir_acc:
                n_missing_fields += 1
                continue

            mir_name_v22 = mirbase_acc2name.get(mir_acc)
            if mir_name_v22 is None:
                n_acc_missing_in_mirbase += 1
                continue

            if mir_name_ts == mir_name_v22:
                n_name_match += 1
            else:
                n_name_replaced += 1

            if mir_name_ts in annotations:
                existing = annotations[mir_name_ts]
                if existing.mirna_id != mir_acc or existing.mirna_name != mir_name_v22:
                    raise ValueError(
                        f"Conflicting human mature miRNA annotation for {mir_name_ts!r}: "
                        f"{existing.mirna_id}/{existing.mirna_name} vs {mir_acc}/{mir_name_v22}"
                    )
                n_duplicate_names += 1
                continue

            annotations[mir_name_ts] = MirnaEntry(mirna_id=mir_acc, mirna_name=mir_name_v22)
            n_kept += 1

    logger.info(
        "miR_Family_Info stats: rows_total=%d | rows_species(%s)=%d | kept=%d | unique_human_mature_mirnas=%d",
        n_total, species_id, n_species, n_kept, len(annotations)
    )
    logger.info(
        "miRNA annotation QC vs miRBase v%s: name_match=%d | name_replaced=%d | duplicate_names=%d | dropped_missing_accession=%d | dropped_missing_fields=%d",
        MIRBASE_RELEASE, n_name_match, n_name_replaced, n_duplicate_names, n_acc_missing_in_mirbase, n_missing_fields
    )

    if annotations:
        key = next(iter(annotations.keys()))
        logger.info(
            "  sample annotation: %s -> %s",
            key, annotations[key].mirna_id
        )

    return annotations


def step_write_standardized_predictions(
    summary_counts_path: pathlib.Path,
    *,
    tx_index: Dict[str, Any],
    ensembl_tables: Dict[str, Dict[str, str]],
    mirna_annotations: Dict[str, MirnaEntry],
    output_path: pathlib.Path,
    species_id: str = "9606",
) -> None:
    rep_txs_by_gene_id = tx_index["rep_txs_by_gene_id"]
    targetscan_gene_id_by_tx = tx_index["targetscan_gene_id_by_tx"]
    tx_to_gene_ens = ensembl_tables["tx_to_gene"]
    gene_to_name_ens = ensembl_tables["gene_to_name"]

    representative_tx_set = {
        tx
        for reps in rep_txs_by_gene_id.values()
        for tx, _tags, _sym in reps
    }

    output_path = pathlib.Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    stats = Counter()
    mismatch_examples = []
    mismatch_example_keys = set()
    null_score_examples = []
    duplicate_examples = []
    unique_rows: Dict[Tuple[str, str], Dict[str, str | float]] = {}
    pre_mismatch_pair_counts: Counter[Tuple[str, str]] = Counter()
    post_mismatch_pair_counts: Counter[Tuple[str, str]] = Counter()

    with pathlib.Path(summary_counts_path).open("r", encoding="utf-8") as handle:
        summary_reader = csv.DictReader(handle, delimiter="\t")

        for row in summary_reader:
            if row["Species ID"].strip() != species_id:
                continue

            tx_raw = row["Transcript ID"].strip()
            if tx_raw not in representative_tx_set:
                continue

            tx_stable = strip_version(tx_raw)
            gene_id_ens = tx_to_gene_ens.get(tx_stable)
            if gene_id_ens is None:
                stats["rows_missing_ensembl_gene"] += 1
                continue

            representative_mirna = row["Representative miRNA"].strip()
            mirna_entry = mirna_annotations.get(representative_mirna)
            if mirna_entry is None:
                stats["rows_missing_mirna_annotation"] += 1
                continue

            pre_mismatch_pair_counts[(gene_id_ens, mirna_entry.mirna_id)] += 1

            targetscan_gene_id = targetscan_gene_id_by_tx.get(tx_raw)
            if targetscan_gene_id != gene_id_ens:
                stats["rows_gene_id_mismatch"] += 1
                sample_key = (tx_raw, targetscan_gene_id or "NA", gene_id_ens)
                if len(mismatch_examples) < 10 and sample_key not in mismatch_example_keys:
                    mismatch_example_keys.add(sample_key)
                    mismatch_examples.append((tx_raw, targetscan_gene_id or "NA", gene_id_ens))
                continue

            raw_score = to_float_or_none(row[TARGETSCAN_SCORE_COL])
            if raw_score is None:
                stats["rows_null_score"] += 1
                if len(null_score_examples) < 10:
                    null_score_examples.append((tx_raw, representative_mirna))
                continue

            stats["rows_after_filters"] += 1
            key = (gene_id_ens, mirna_entry.mirna_id)
            post_mismatch_pair_counts[key] += 1
            if key in unique_rows:
                stats["duplicate_gene_mirna_rows"] += 1
                if len(duplicate_examples) < 10:
                    duplicate_examples.append(
                        (
                            gene_id_ens,
                            mirna_entry.mirna_id,
                            representative_mirna,
                            float(unique_rows[key]["Score"]),
                            raw_score,
                            tx_raw,
                        )
                    )
                continue

            unique_rows[key] = {
                "Ensembl_ID": gene_id_ens,
                "Gene_Name": gene_to_name_ens.get(gene_id_ens) or row["Gene Symbol"].strip(),
                "miRNA_ID": mirna_entry.mirna_id,
                "miRNA_Name": mirna_entry.mirna_name,
                "Score": raw_score,
            }

    if stats["rows_null_score"]:
        details = ", ".join(f"{tx}/{mirna}" for tx, mirna in null_score_examples)
        raise ValueError(
            f"Unexpected NULL {TARGETSCAN_SCORE_COL!r} values after filtering: "
            f"{stats['rows_null_score']} rows. Sample rows: {details}"
        )

    if stats["rows_after_filters"] == 0:
        raise ValueError(f"No rows remain after filters for {PREDICTOR_NAME}.")

    duplicate_pairs_before_mismatch = sum(1 for count in pre_mismatch_pair_counts.values() if count > 1)
    duplicate_pairs_after_mismatch = sum(1 for count in post_mismatch_pair_counts.values() if count > 1)
    rep_gene_mismatch_count = sum(
        1 for reps in rep_txs_by_gene_id.values()
        for tx_raw, _tags, _sym in reps
        if (
            (targetscan_gene_id_by_tx.get(tx_raw) is not None)
            and (strip_version(tx_raw) in tx_to_gene_ens)
            and (targetscan_gene_id_by_tx.get(tx_raw) != tx_to_gene_ens[strip_version(tx_raw)])
        )
    )

    logger.info(
        "%s score column '%s': kept raw TargetScan score | rows_missing_mirna_annotation=%d | rows_missing_ensembl_gene=%d | rows_gene_id_mismatch=%d | NULL rows=%d",
        PREDICTOR_NAME,
        TARGETSCAN_SCORE_COL,
        stats["rows_missing_mirna_annotation"],
        stats["rows_missing_ensembl_gene"],
        stats["rows_gene_id_mismatch"],
        stats["rows_null_score"],
    )
    logger.info(
        "Filtering summary: %s NULL scores after filtering | %s representative transcripts whose TargetScan gene ID disagrees with Ensembl",
        f"{stats['rows_null_score']:,}",
        f"{rep_gene_mismatch_count:,}",
    )
    logger.info(
        "Filtering summary: %s summary rows dropped by the TargetScan/Ensembl gene mismatch filter | duplicate gene-miRNA pairs before mismatch filter=%s | surviving duplicate gene-miRNA pairs after mismatch filter=%s",
        f"{stats['rows_gene_id_mismatch']:,}",
        f"{duplicate_pairs_before_mismatch:,}",
        f"{duplicate_pairs_after_mismatch:,}",
    )
    if mismatch_examples:
        logger.info("Sample rows dropped for TargetScan/Ensembl gene-ID mismatch (tx -> TargetScan gene -> Ensembl gene):")
        for tx_raw, targetscan_gene_id, gene_id_ens in mismatch_examples:
            logger.info("  %s -> %s -> %s", tx_raw, targetscan_gene_id, gene_id_ens)

    if stats["duplicate_gene_mirna_rows"]:
        details = "; ".join(
            (
                f"{gene_id}/{mirna_id} ({mirna_name}) "
                f"existing_score={existing_score:.6g} new_score={new_score:.6g} "
                f"tx={tx_raw}"
            )
            for gene_id, mirna_id, mirna_name, existing_score, new_score, tx_raw in duplicate_examples
        )
        logger.error(
            "Unexpected duplicate gene-miRNA rows remain after all filters: duplicate_rows=%s | duplicate_pairs=%s",
            f"{stats['duplicate_gene_mirna_rows']:,}",
            f"{duplicate_pairs_after_mismatch:,}",
        )
        raise ValueError(
            "Unexpected duplicate gene-miRNA rows remain after all TargetScan filters: "
            f"{stats['duplicate_gene_mirna_rows']} duplicate rows across "
            f"{duplicate_pairs_after_mismatch} duplicate pairs. "
            f"Sample rows: {details}"
        )

    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["Ensembl_ID", "Gene_Name", "miRNA_ID", "miRNA_Name", "Score"],
            delimiter="\t",
        )
        writer.writeheader()
        for record in sorted(unique_rows.values(), key=lambda row: (str(row["Ensembl_ID"]), str(row["miRNA_ID"]))):
            writer.writerow(
                {
                    **record,
                    "Score": f"{float(record['Score']):.6g}",
                }
            )

    logger.info(
        "%s -> retained %d unique gene-miRNA rows from %d transcript-level rows",
        PREDICTOR_NAME,
        len(unique_rows),
        stats["rows_after_filters"],
    )
    logger.info("Output written to: %s", display_path(output_path))
    logger.info("Rows written: %d", len(unique_rows))


def compute_final_statistics(predictions_path: pathlib.Path) -> None:
    predictions_path = pathlib.Path(predictions_path)
    genes = set()
    mirs = set()
    n_rows = 0

    with predictions_path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            n_rows += 1
            genes.add(row["Ensembl_ID"])
            mirs.add(row["miRNA_ID"])

    logger.info("%s: %d rows | %d unique genes | %d unique miRNAs", PREDICTOR_NAME, n_rows, len(genes), len(mirs))
