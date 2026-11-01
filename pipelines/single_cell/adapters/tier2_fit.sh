#!/usr/bin/env bash
# Smoke fit of the first comparison in scratch, then all comparisons into the final folder, next to the store.
# Usage: tier2_fit.sh <dataset> <store.h5ad, moved to the final folder> <final folder, must not exist> <workers>
set -euo pipefail
W=/SCRATCH/dtzim01/funmirbench/workspaces/single-cell-processing
DATASET=$1; STORE=$2; FINAL=$3; WORKERS=$4
LIST=$W/pipelines/single_cell/comparisons/$DATASET.tsv
SMOKE=/SCRATCH/dtzim01/tmp/tier2/smoke/$DATASET
[ ! -e "$FINAL" ] || { echo "$FINAL exists"; exit 1; }
rm -rf "$SMOKE"; mkdir -p "$SMOKE"
head -2 "$LIST" > "$SMOKE/first.tsv"
cd "$W"
.venv/bin/funmirbench-single-cell-final --store "$STORE" --comparisons "$SMOKE/first.tsv" --genes metadata/benchmark_genes.tsv.gz --out "$SMOKE/out" --workers 1
mkdir -p "$FINAL/store"
mv "$STORE" "$FINAL/store/$DATASET.h5ad"
.venv/bin/funmirbench-single-cell-final --store "$FINAL/store/$DATASET.h5ad" --comparisons "$LIST" --genes metadata/benchmark_genes.tsv.gz --out "$FINAL" --workers "$WORKERS"
