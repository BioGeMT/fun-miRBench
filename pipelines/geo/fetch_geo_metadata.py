"""
Fetch GEO series metadata and append a review row to pipelines/geo/input_experiments.tsv.

Uses the GEO SOFT text API to retrieve series and sample-level metadata without
requiring SRA credentials. Optionally uses the Gemini Flash LLM for smarter field
extraction. Auto-fills what it can from GEO; prints a summary of what still needs
manual editing before running geo_download.py.

Usage:
    python pipelines/geo/fetch_geo_metadata.py --gse-url GSE93717
    python pipelines/geo/fetch_geo_metadata.py --gse-url GSE93717 --llm gemini
    python pipelines/geo/fetch_geo_metadata.py --gse-url https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE93717
"""

import argparse
import csv
import json
import os
import re
import sys
import time
from pathlib import Path
from urllib.parse import parse_qs, urlparse

import requests

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
INPUT_TSV = Path(__file__).resolve().parent / "input_experiments.tsv"
MIRBASE_CACHE_DIR = REPO_ROOT / "data" / "mirbase"

TSV_COLUMNS = [
    "id", "geo_accession", "mirna_name", "experiment_type",
    "tested_cell_line", "tissue", "organism", "method", "pubmed_id",
    "control_samples", "condition_samples", "raw_data_dir",
    "count_matrix_path", "gene_id_column", "treatment",
]


GEO_SOFT_URL = (
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
    "?acc={gse}&targ=all&form=text&view=brief"
)

MIRBASE_MATURE_FA_URL = "https://mirbase.org/download/mature.fa"
MIRBASE_CACHE_MAX_AGE_DAYS = 30  # re-download if cache is older than this many days


# Keywords used to classify samples as likely control or likely condition.
# Scoring: each match adds +1 to that side; highest score wins.
_CONTROL_KEYWORDS = [
    r"\bcontrol\b", r"\bctrl\b", r"\bmock\b", r"\bscramble[d]?\b",
    r"\bneg(?:ative)?\b", r"\bNC\b", r"\bempty.vector\b",
    r"\bmir.ctrl\b", r"\bnon.targeting\b", r"\bwild.?type\b", r"\bWT\b",
]
_CONDITION_KEYWORDS = [
    r"\bmimic\b", r"\boverexpression\b", r"\bOE\b",
    r"\binhibitor\b", r"\bantago(?:mir)?\b", r"\bknockout\b",
    r"\bKO\b", r"\bknockdown\b", r"\bKD\b", r"\bsuppression\b",
    r"\btransfect\b",
]


# ---------------------------------------------------------------------------
# GEO SOFT parser
# ---------------------------------------------------------------------------

def extract_gse_accession(gse_url: str) -> str:
    gse_url = str(gse_url).strip()
    parsed = urlparse(gse_url)
    accession = parse_qs(parsed.query).get("acc", [""])[0].strip()
    if accession:
        return accession
    tail = parsed.path.rstrip("/").split("/")[-1].strip()
    if tail.upper().startswith("GSE"):
        return tail
    if gse_url.upper().startswith("GSE"):
        return gse_url
    raise ValueError(f"Cannot extract GSE accession from: {gse_url!r}")


def fetch_soft(gse: str) -> str:
    url = GEO_SOFT_URL.format(gse=gse)
    try:
        response = requests.get(url, timeout=60)
        response.raise_for_status()
    except requests.RequestException as e:
        raise RuntimeError(f"Failed to fetch GEO SOFT for {gse}: {e}")
    return response.text


def parse_soft(soft_text: str) -> dict:
    """
    Parse GEO SOFT text into a dict with keys:
      series: {field: value_or_list}
      samples: {GSM_id: {field: value_or_list}}
    """
    result = {"series": {}, "samples": {}}
    current_block = None
    current_id = None

    for line in soft_text.splitlines():
        line = line.strip()
        if not line:
            continue

        # Block headers
        if line.startswith("^SERIES"):
            current_block = "series"
            current_id = None
            continue
        if line.startswith("^SAMPLE"):
            parts = line.split("=", 1)
            current_id = parts[1].strip() if len(parts) > 1 else None
            current_block = "sample"
            if current_id:
                result["samples"].setdefault(current_id, {})
            continue
        if line.startswith("^"):
            current_block = None
            current_id = None
            continue

        # Field lines
        if not line.startswith("!"):
            continue
        parts = line[1:].split("=", 1)
        if len(parts) != 2:
            continue
        key = parts[0].strip()
        value = parts[1].strip()

        if current_block == "series":
            existing = result["series"].get(key)
            if existing is None:
                result["series"][key] = value
            elif isinstance(existing, list):
                existing.append(value)
            else:
                result["series"][key] = [existing, value]
        elif current_block == "sample" and current_id:
            existing = result["samples"][current_id].get(key)
            if existing is None:
                result["samples"][current_id][key] = value
            elif isinstance(existing, list):
                existing.append(value)
            else:
                result["samples"][current_id][key] = [existing, value]

    return result


# ---------------------------------------------------------------------------
# Metadata extraction helpers
# ---------------------------------------------------------------------------

def _as_list(value) -> list:
    if value is None:
        return []
    if isinstance(value, list):
        return value
    return [value]


def _first(value, default="") -> str:
    lst = _as_list(value)
    return lst[0].strip() if lst else default


def extract_pubmed_id(series: dict) -> str:
    pmids = _as_list(series.get("Series_pubmed_id"))
    if pmids:
        return pmids[0].strip()
    return ""


_CELL_LINE_SUFFIXES = (" cells", " cell line", " cell")


def _clean_cell_line(name: str) -> str:
    """Strip generic trailing words from cell line names (e.g. 'HaCaT cells' → 'HaCaT')."""
    for suffix in _CELL_LINE_SUFFIXES:
        if name.lower().endswith(suffix):
            return name[:-len(suffix)].strip()
    return name


def extract_sample_info(sample_data: dict) -> dict:
    """Flatten a sample's SOFT fields into a clean dict."""
    chars = {}
    for val in _as_list(sample_data.get("Sample_characteristics_ch1")):
        if ":" in val:
            k, v = val.split(":", 1)
            chars[k.strip().lower()] = v.strip()

    # Cell line: try characteristics first, then source
    raw_cell_line = (
        chars.get("cell line")
        or chars.get("cell type")
        or _first(sample_data.get("Sample_source_name_ch1"))
    )
    cell_line = _clean_cell_line(raw_cell_line) if raw_cell_line else ""
    tissue = chars.get("tissue") or chars.get("tissue type") or ""
    organism = _first(sample_data.get("Sample_organism_ch1"))
    title = _first(sample_data.get("Sample_title"))

    return {
        "title": title,
        "organism": organism,
        "cell_line": cell_line,
        "tissue": tissue,
        "characteristics": chars,
    }


def _score_sample(title: str, chars: dict) -> tuple[int, int]:
    """Return (control_score, condition_score) for a sample."""
    text = (title + " " + " ".join(chars.values())).lower()
    ctrl = sum(1 for pat in _CONTROL_KEYWORDS if re.search(pat, text, re.IGNORECASE))
    cond = sum(1 for pat in _CONDITION_KEYWORDS if re.search(pat, text, re.IGNORECASE))
    return ctrl, cond


def classify_samples(samples: dict) -> dict:
    """
    Classify each GSM sample as 'control', 'condition', or 'uncertain'.

    Returns {gsm_id: {title, organism, cell_line, tissue, group, ctrl_score, cond_score}}.
    """
    classified = {}
    for gsm, data in samples.items():
        info = extract_sample_info(data)
        ctrl_score, cond_score = _score_sample(info["title"], info["characteristics"])

        if ctrl_score > cond_score:
            group = "control"
        elif cond_score > ctrl_score:
            group = "condition"
        else:
            group = "uncertain"

        classified[gsm] = {**info, "group": group, "ctrl_score": ctrl_score, "cond_score": cond_score}

    return classified


# ---------------------------------------------------------------------------
# miRBase validation
# ---------------------------------------------------------------------------

def _get_cached_mature_fa(
    cache_dir: Path | None = None,
    max_age_days: int = MIRBASE_CACHE_MAX_AGE_DAYS,
) -> Path:
    """
    Return path to a local copy of miRBase mature.fa, downloading/refreshing as needed.

    The file is re-downloaded if it does not exist or is older than max_age_days.
    If a refresh fails but a stale copy exists, the stale copy is used with a warning.
    """
    if cache_dir is None:
        cache_dir = MIRBASE_CACHE_DIR
    cache_path = Path(cache_dir) / "mature.fa"

    needs_download = True
    if cache_path.exists():
        age_days = (time.time() - cache_path.stat().st_mtime) / 86400
        if age_days < max_age_days:
            needs_download = False
        else:
            print(
                f"[miRBase] Cache is {age_days:.0f} days old (limit {max_age_days}). Refreshing...",
                file=sys.stderr,
            )

    if needs_download:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        print(f"[miRBase] Downloading mature.fa from {MIRBASE_MATURE_FA_URL} ...", file=sys.stderr)
        try:
            response = requests.get(MIRBASE_MATURE_FA_URL, timeout=120, stream=True)
            response.raise_for_status()
            tmp_path = cache_path.with_suffix(".fa.tmp")
            with open(tmp_path, "wb") as fh:
                for chunk in response.iter_content(chunk_size=65536):
                    fh.write(chunk)
            tmp_path.replace(cache_path)  # atomic rename
            print(f"[miRBase] Saved to {cache_path}", file=sys.stderr)
        except requests.RequestException as e:
            if cache_path.exists():
                print(
                    f"[miRBase] WARNING: refresh failed ({e}); using stale cache at {cache_path}",
                    file=sys.stderr,
                )
            else:
                raise RuntimeError(f"Failed to download miRBase mature.fa: {e}")

    return cache_path


def _parse_mature_fa(path: Path) -> dict:
    """Parse a miRBase mature.fa FASTA file into {mirna_name: sequence}."""
    sequences = {}
    current_name = None
    current_seq: list[str] = []

    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if current_name:
                    sequences[current_name] = "".join(current_seq)
                # Header format: >hsa-miR-21-5p MIMAT0000076 Homo sapiens miR-21-5p
                current_name = line[1:].split()[0]
                current_seq = []
            elif current_name:
                current_seq.append(line)
        if current_name:
            sequences[current_name] = "".join(current_seq)

    return sequences


def validate_mirna_in_mirbase(name: str) -> bool | None:
    """
    Return True if name is found in miRBase mature.fa, False if not found,
    or None if the check could not be performed (e.g. download failure).
    """
    try:
        cache_path = _get_cached_mature_fa()
        sequences = _parse_mature_fa(cache_path)
        return name in sequences
    except RuntimeError:
        return None


# ---------------------------------------------------------------------------
# LLM integration (Gemini Flash)
# ---------------------------------------------------------------------------

def call_gemini(
    title: str,
    summary: str,
    sample_titles: list[str],
    api_key: str,
) -> dict:
    """
    Use Gemini Flash to extract experiment metadata fields from GEO series text.

    Returns a dict with some/all of: mirna_name, experiment_type, treatment,
    tested_cell_line, tissue, organism. Values are strings or None.
    """
    # Uses the Gemini REST API directly — no extra package required.

    sample_list = "\n".join(f"  - {t}" for t in sample_titles[:20])

    prompt = f"""You are a bioinformatics expert. Given the following GEO (Gene Expression Omnibus) experiment metadata, extract the fields listed below. Return a JSON object with EXACTLY these keys:

{{
  "mirna_name": "exact miRBase mature name (e.g. hsa-miR-21-5p) or null",
  "experiment_type": "Overexpression, Knockout, Knockdown, or null",
  "treatment": "short description of the experimental treatment or null",
  "tested_cell_line": "cell line name without trailing 'cells' or 'cell line' (e.g. HaCaT, A549) or null",
  "tissue": "tissue or organ type (e.g. Lung, Skin, Ovarian) or null",
  "organism": "full organism name (e.g. Homo sapiens, Mus musculus) or null"
}}

Extraction rules:
- mirna_name: Use exact miRBase mature miRNA format with organism prefix and arm suffix (e.g. hsa-miR-21-5p). If the arm (-3p/-5p) is ambiguous, pick the more likely one based on context. Return null if no specific miRNA is mentioned.
- experiment_type: "Overexpression" for overexpression/mimic/gain-of-function, "Knockout" for knockout, and "Knockdown" for knockdown/inhibition/antagomir/loss-of-function. Return null if unclear.
- treatment: Short phrase describing what was done (e.g. "Overexpression of miR-21-5p in HaCaT cells").
- tested_cell_line: Clean name without "cells" or "cell line" suffix.
- tissue: The tissue/organ this cell line originates from or the tissue type studied.
- organism: Full Latin species name.

Series title: {title}

Series summary: {summary[:2000]}

Sample titles:
{sample_list}

Return ONLY valid JSON — no markdown code blocks, no explanations."""

    url = (
        "https://generativelanguage.googleapis.com/v1beta/models/"
        f"gemini-2.0-flash:generateContent?key={api_key}"
    )
    payload = {"contents": [{"parts": [{"text": prompt}]}]}

    try:
        resp = requests.post(url, json=payload, timeout=60)
        resp.raise_for_status()
        text = resp.json()["candidates"][0]["content"]["parts"][0]["text"].strip()
        # Strip markdown code fences if the model added them
        text = re.sub(r"^```(?:json)?\s*", "", text)
        text = re.sub(r"\s*```$", "", text)
        fields = json.loads(text)
        expected = {"mirna_name", "experiment_type", "treatment", "tested_cell_line", "tissue", "organism"}
        result = {}
        for k in expected:
            v = fields.get(k)
            result[k] = str(v).strip() if v and str(v).strip().lower() not in ("null", "none", "") else None
        return result
    except Exception as e:
        print(f"[LLM] Gemini call failed: {e}", file=sys.stderr)
        return {}


# ---------------------------------------------------------------------------
# Row builder
# ---------------------------------------------------------------------------

def build_row(
    gse: str,
    soft_parsed: dict,
    llm_fields: dict | None = None,
) -> tuple[dict, dict, dict, dict]:
    """
    Build a TSV row dict from parsed GEO SOFT data, optionally merging LLM fields.

    Returns (row, classified_samples, series_info, sources).
    sources maps each row field to one of:
      "rule-based" | "llm" | "both" | "default" | "manual" | "llm (was: <old>)"
    """
    series = soft_parsed["series"]
    classified = classify_samples(soft_parsed["samples"])

    pubmed_id = extract_pubmed_id(series)

    organisms = [s["organism"] for s in classified.values() if s["organism"]]
    cell_lines = [s["cell_line"] for s in classified.values() if s["cell_line"]]
    tissues = [s["tissue"] for s in classified.values() if s["tissue"]]

    def majority(lst):
        return max(set(lst), key=lst.count) if lst else "TO BE FILLED"

    control_gsms = [gsm for gsm, s in classified.items() if s["group"] == "control"]
    condition_gsms = [gsm for gsm, s in classified.items() if s["group"] == "condition"]

    row = {
        "id": "TO BE FILLED",
        "geo_accession": gse,
        "mirna_name": "TO BE FILLED",
        "experiment_type": "TO BE FILLED",
        "tested_cell_line": majority(cell_lines),
        "tissue": majority(tissues),
        "organism": majority(organisms),
        "method": "RNA-seq",
        "pubmed_id": pubmed_id,
        "control_samples": ",".join(control_gsms),
        "condition_samples": ",".join(condition_gsms),
        "raw_data_dir": "",
        "count_matrix_path": "",
        "gene_id_column": "",
        "treatment": "TO BE FILLED",
    }

    sources = {
        "id": "manual",
        "geo_accession": "rule-based",
        "mirna_name": "manual",
        "experiment_type": "manual",
        "tested_cell_line": "rule-based" if cell_lines else "manual",
        "tissue": "rule-based" if tissues else "manual",
        "organism": "rule-based" if organisms else "manual",
        "method": "default",
        "pubmed_id": "rule-based" if pubmed_id else "default",
        "control_samples": "rule-based",
        "condition_samples": "rule-based",
        "raw_data_dir": "default",
        "count_matrix_path": "default",
        "gene_id_column": "default",
        "treatment": "manual",
    }

    # Merge LLM fields where available
    if llm_fields:
        _llm_field_map = {
            "mirna_name": "mirna_name",
            "experiment_type": "experiment_type",
            "treatment": "treatment",
            "tested_cell_line": "tested_cell_line",
            "tissue": "tissue",
            "organism": "organism",
        }
        for llm_key, row_key in _llm_field_map.items():
            llm_val = llm_fields.get(llm_key)
            if not llm_val:
                continue
            current = row[row_key]
            if current == "TO BE FILLED":
                row[row_key] = llm_val
                sources[row_key] = "llm"
            elif current.lower() == llm_val.lower():
                sources[row_key] = "both"
            else:
                # LLM disagrees with rule-based; prefer LLM, note the original value
                sources[row_key] = f"llm (rule-based was: {current})"
                row[row_key] = llm_val

    series_info = {
        "title": _first(series.get("Series_title")),
        "summary": _first(series.get("Series_summary")),
    }

    return row, classified, series_info, sources


# ---------------------------------------------------------------------------
# TSV writer
# ---------------------------------------------------------------------------

def append_to_tsv(row: dict, tsv_path: Path) -> bool:
    """
    Append a review row using the current GEO input schema.

    Refuse to append if an existing TSV has a different header, preventing
    values from being silently written under stale column names.
    """
    if tsv_path.exists():
        with open(tsv_path, newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            existing_columns = reader.fieldnames or []
            if existing_columns != TSV_COLUMNS:
                raise ValueError(
                    f"{tsv_path} uses an incompatible header. "
                    f"Expected: {TSV_COLUMNS}; found: {existing_columns}"
                )
            for existing in reader:
                if existing.get("geo_accession") == row["geo_accession"]:
                    print(
                        f"WARNING: {row['geo_accession']} already exists in {tsv_path.name}. "
                        "Skipping. Remove the existing row first if you want to re-add it.",
                        file=sys.stderr,
                    )
                    return False

    write_header = not tsv_path.exists()
    with open(tsv_path, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=TSV_COLUMNS, delimiter="\t")
        if write_header:
            writer.writeheader()
        writer.writerow(row)
    return True


# ---------------------------------------------------------------------------
# Summary printer
# ---------------------------------------------------------------------------

_SOURCE_LABELS = {
    "rule-based": "[rule-based]",
    "both":        "[LLM+rule]",
    "default":     "[default]",
    "manual":      "",
    "llm":         "[LLM]",
}


def _source_label(src: str) -> str:
    return _SOURCE_LABELS.get(src, f"[{src}]")


def print_summary(
    gse: str,
    row: dict,
    classified: dict,
    series_info: dict,
    tsv_path: Path,
    sources: dict,
    mirbase_valid: bool | None = None,
    llm_used: bool = False,
):
    sep = "=" * 60
    warn = "!" * 60
    print(warn)
    print("  IMPORTANT — PLEASE READ BEFORE PROCEEDING")
    print()
    if llm_used:
        print("  This script uses rule-based heuristics AND the Gemini Flash")
        print("  LLM to auto-fill fields. LLM output can contain hallucinations.")
    else:
        print("  This script uses rule-based heuristics (keyword matching")
        print("  and majority voting), NOT an LLM or AI model.")
    print("  All auto-filled fields may contain errors.")
    print("  Go through each field carefully and verify manually")
    print("  before running geo_download.py.")
    print(warn)
    print()
    print(sep)
    print(f"  {gse} — {series_info['title']}")
    print(sep)
    print(f"\nRow appended to: {tsv_path}\n")

    # Fields with auto-filled values
    auto_fields = [
        f for f in TSV_COLUMNS
        if sources.get(f) not in ("manual", "default")
        and row.get(f) not in ("TO BE FILLED", "", None)
    ]
    if auto_fields:
        print("Auto-filled fields (verify each one):")
        for field in auto_fields:
            val = row[field]
            label = _source_label(sources.get(field, ""))
            extra = ""
            if field == "mirna_name" and mirbase_valid is not None:
                if mirbase_valid:
                    extra = "  ← miRBase: VALID ✓"
                else:
                    extra = "  ← miRBase: NOT FOUND — check name and arm (-3p/-5p)"
            print(f"  {field:<22} {val}  {label}{extra}")

    # Fields still needing manual input
    manual_fields = [f for f in TSV_COLUMNS if row.get(f) == "TO BE FILLED"]
    if manual_fields:
        print("\nFields still needing manual edit (open the TSV):")
        hints = {
            "id":              "e.g. {gse}_{experiment_type}_{mirna_name_safe}",
            "mirna_name":      "exact miRBase mature name (e.g. hsa-miR-21-5p, NOT hsa-mir-21-5p)\n"
                               "                           wrong arm (-3p/-5p) will break downstream analysis\n"
                               "                           verify at https://mirbase.org",
            "experiment_type": "Overexpression, Knockout, or Knockdown",
            "treatment":       "short description of the experiment",
        }
        for field in manual_fields:
            hint = hints.get(field, "")
            print(f"  {field:<22} → {hint}" if hint else f"  {field}")

    uncertain = [(gsm, s) for gsm, s in classified.items() if s["group"] == "uncertain"]
    if uncertain:
        print(f"\nWARNING: {len(uncertain)} sample(s) could not be auto-classified:")
        for gsm, s in uncertain:
            print(f"  {gsm}  ctrl_score={s['ctrl_score']}  cond_score={s['cond_score']}  \"{s['title']}\"")
        print("  → Please verify control_samples / condition_samples in the TSV.")

    print("\nSample classification:")
    for gsm, s in classified.items():
        flag = "  ← uncertain" if s["group"] == "uncertain" else ""
        print(f"  {gsm}  {s['group']:<10}  (ctrl={s['ctrl_score']}, cond={s['cond_score']})  \"{s['title']}\"{flag}")

    print(f"\nSeries summary:\n  {series_info['summary']}\n")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Fetch GEO series metadata and append a review row to input_experiments.tsv."
        )
    )
    parser.add_argument(
        "--gse-url", required=True,
        help="GEO series URL or accession (e.g. GSE93717 or full GEO URL)",
    )
    parser.add_argument(
        "--tsv", default=str(INPUT_TSV),
        help=f"Path to input_experiments.tsv (default: {INPUT_TSV})",
    )
    parser.add_argument(
        "--llm", choices=["gemini"], default=None,
        help="Use an LLM to help fill in fields (e.g. --llm gemini)",
    )
    parser.add_argument(
        "--gemini-key", default=None,
        help=(
            "Gemini API key. If omitted, reads GEMINI_API_KEY from the environment. "
            "Required when --llm gemini is set."
        ),
    )
    args = parser.parse_args()

    try:
        gse = extract_gse_accession(args.gse_url)
    except ValueError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1

    try:
        soft_text = fetch_soft(gse)
    except RuntimeError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1

    soft_parsed = parse_soft(soft_text)
    if not soft_parsed["series"]:
        print(f"ERROR: No series data found for {gse}. Check the accession.", file=sys.stderr)
        return 1

    # Optionally call LLM
    llm_fields = None
    llm_used = False
    if args.llm == "gemini":
        api_key = args.gemini_key or os.environ.get("GEMINI_API_KEY", "")
        if not api_key:
            print(
                "ERROR: Gemini API key required. Pass --gemini-key or set GEMINI_API_KEY.",
                file=sys.stderr,
            )
            return 1
        series = soft_parsed["series"]
        title = _first(series.get("Series_title"))
        summary = _first(series.get("Series_summary"))
        sample_titles = [
            _first(data.get("Sample_title"))
            for data in soft_parsed["samples"].values()
        ]
        print(f"[LLM] Calling Gemini Flash for {gse}...", file=sys.stderr)
        llm_fields = call_gemini(title, summary, sample_titles, api_key)
        llm_used = bool(llm_fields)
        if llm_used:
            print(f"[LLM] Received fields: {list(k for k, v in llm_fields.items() if v)}", file=sys.stderr)
        else:
            print("[LLM] No fields returned; falling back to rule-based only.", file=sys.stderr)

    row, classified, series_info, sources = build_row(gse, soft_parsed, llm_fields=llm_fields)

    # Validate mirna_name in miRBase if it was filled
    mirbase_valid = None
    if row.get("mirna_name") not in ("TO BE FILLED", "", None):
        print(f"[miRBase] Validating '{row['mirna_name']}'...", file=sys.stderr)
        mirbase_valid = validate_mirna_in_mirbase(row["mirna_name"])

    tsv_path = Path(args.tsv)
    appended = append_to_tsv(row, tsv_path)
    if appended:
        print_summary(
            gse, row, classified, series_info, tsv_path,
            sources=sources,
            mirbase_valid=mirbase_valid,
            llm_used=llm_used,
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
