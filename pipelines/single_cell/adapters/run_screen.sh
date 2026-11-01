#!/bin/bash
# One screen: adapter (unless the store exists), smoke fit of the first comparison in scratch, all fits, outcome row.
# Usage (as the command of run.sh): run_screen.sh <dataset> <workers> <run_dir> <source_files> <size_gb> [adapter args...]
set -euo pipefail
dataset=$1; workers=$2; run_dir=$3; source_files=$4; size_gb=$5; shift 5
repo=/SCRATCH/dtzim01/funmirbench/workspaces/single-cell-processing
final=/SCRATCH/dtzim01/fun-experiment-data/single-cell-final/20261031/$dataset
scratch=/SCRATCH/dtzim01/tmp/screens/smoke/$dataset
export PATH=/SCRATCH/dtzim01/fun-experiment-data/tools/envs/deseq2-shrinkage/bin:$PATH
cd "$repo"
store_dataset=${dataset%_pc_hosts}
store=/SCRATCH/dtzim01/fun-experiment-data/single-cell-final/20261031/$store_dataset/store/$store_dataset.h5ad
list=pipelines/single_cell/comparisons/$dataset.tsv
if [ ! -s "$store" ] || [ ! -s "$list" ]; then
  (cd pipelines/single_cell/adapters && uv run python "$@")
fi
rm -rf "$scratch"; mkdir -p "$(dirname "$scratch")"
head -2 "$list" > "$scratch.tsv"
uv run funmirbench-single-cell-final --store "$store" --comparisons "$scratch.tsv" --genes metadata/benchmark_genes.tsv.gz --out "$scratch" --workers 1 --log-level WARNING | grep -v '^PROGRESS'
rm -rf "$scratch" "$scratch.tsv"
uv run funmirbench-single-cell-final --store "$store" --comparisons "$list" --genes metadata/benchmark_genes.tsv.gz --out "$final" --workers "$workers" --units "$run_dir/units.tsv" --log-level WARNING
(cd pipelines/single_cell/adapters && uv run python record_screen.py "$dataset" "$source_files" "$size_gb")
