# TargetScan

This directory contains the standardization pipeline for TargetScan v8 gene-level predictions.

## Files

- `pipeline.py`: CLI entrypoint for the pipeline.
- `utils.py`: shared helpers for logging, downloads, parsing, QC, mapping, and output construction.
- `data/predictions/targetscan/targetscan_pipeline.log`: example log from a completed run, relative to the repository root.
- `data/`: downloaded TargetScan prediction data.
- `resources/`: downloaded TargetScan metadata and derived Ensembl lookup tables.

Shared miRBase and Ensembl annotations are downloaded to:

```text
data/common_resources/mirbase/mature.fa
data/common_resources/ensembl/Homo_sapiens.GRCh38.115.gtf.gz
```

By default, `Summary_Counts.all_predictions.txt` and its ZIP archive are stored under `data/`.
`Gene_info.txt`, `miR_Family_Info.txt`, their ZIP archives, and derived Ensembl lookup tables are
stored under `resources/`.

## What The Pipeline Does

The pipeline:

1. Downloads or reuses the TargetScan v8 input files:
   - `Summary_Counts.all_predictions.txt`
   - `miR_Family_Info.txt`
   - `Gene_info.txt`
2. Builds a representative-transcript index from `Gene_info.txt` for human rows only (`Species ID = 9606`).
3. Keeps every transcript with `Representative transcript? == 1`.
4. Logs whether any gene has more than one representative transcript row.
5. Downloads or reuses the Ensembl v115 GTF for GRCh38.
6. Builds cached Ensembl lookup tables for:
   - transcript stable ID to gene stable ID
   - gene stable ID to gene name
7. Runs QC comparing TargetScan representative transcripts against Ensembl transcript and gene mappings.
8. Downloads or reuses miRBase `mature.fa` version 22.1.
9. Parses human family membership from `miR_Family_Info.txt`, indexed by both TargetScan
   `miR family` and `Seed+m8`.
10. Resolves every valid mature family member through its accession to the current mature name in
    miRBase v22.1. Family rows with missing fields or accessions absent from v22.1 are excluded and
    reported by reason.
11. Reads `Summary_Counts.all_predictions.txt`, resolves its `miRNA family` value to all valid
    human mature members, and expands each retained family-level `(gene, score)` prediction to one
    row per member. `Representative miRNA` is retained for provenance; it does not restrict which
    family members are written.
12. Keeps only human rows and drops rows where:
    - the transcript is not in the representative transcript set
    - the transcript does not map to an Ensembl gene
    - the TargetScan gene ID for that transcript disagrees with the Ensembl gene ID recovered from the transcript
    - the TargetScan family or seed has no valid human member resolved in miRBase v22.1
    - the score is null
13. Preserves the raw TargetScan `Cumulative weighted context++ score` as `Score`.
14. Labels every retained score in a companion provenance audit as either:
    - `official_representative_level_score` for TargetScan's named representative miRNA
    - `family_expanded_score` for another mature member receiving the family-level score
15. Handles duplicate `(Ensembl_ID, miRNA_ID)` pairs deterministically. Identical numeric scores are
    collapsed using a stable source ordering, with the official representative record preferred.
    Conflicting scores raise an error rather than being silently selected.
16. Writes the standardized output, score-provenance audit, and before/after miRNA gene-coverage
    audit.

## Output Schema

The output TSV contains:

- `Ensembl_ID`
- `Gene_Name`
- `miRNA_ID`
- `miRNA_Name`
- `Score`

`Score` is the raw numeric form of TargetScan `Cumulative weighted context++ score`.
For non-representative family members, this is a copied family-level score. It is not a
member-specific score and must not be interpreted or reported as one. The pipeline does not
reconstruct member-specific scores.

## Output Location

By default, the standardized file is written to:

```text
data/predictions/targetscan/targetscan_standardized.tsv
```

relative to the repository root.

Two report-facing audit tables are written alongside it:

- `targetscan_standardized_score_provenance.tsv` links every output pair and score to its
  TargetScan family, seed, named representative miRNA, source transcript, and provenance class.
- `targetscan_standardized_mirna_gene_coverage.tsv` reports gene counts per mature miRNA under the
  old representative-only interpretation and after family expansion, including the delta.

## Run

From the repository root:

```bash
uv run pipelines/standardized_predictors/targetscan/pipeline.py
```

## CLI Arguments

```bash
uv run pipelines/standardized_predictors/targetscan/pipeline.py \
  --common-resources-dir data/common_resources \
  --data-dir pipelines/standardized_predictors/targetscan/data \
  --resources-dir pipelines/standardized_predictors/targetscan/resources \
  --standardized-output-file data/predictions/targetscan/targetscan_standardized.tsv \
  --log-file data/predictions/targetscan/targetscan_pipeline.log \
  --log-level INFO
```

## Logging

Logging is written both to stdout and to the file passed via `--log-file`. By default, the log is written to `data/predictions/targetscan/targetscan_pipeline.log`.
QC logging reports the number of human TargetScan families and seeds, the members-per-family
distribution, unmapped members and reasons, row counts before and after expansion, identical and
conflicting duplicate pairs, and miRNA-level gene coverage before and after family expansion.
