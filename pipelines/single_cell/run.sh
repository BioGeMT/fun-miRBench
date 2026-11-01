#!/usr/bin/env bash
# Build the cell stores and fit the selected comparisons, then write the combined table. Finished datasets are skipped.
# Arguments go to run_selected.py, for example --groups main. RUN_WORKERS is the number of fits at a time per dataset.
set -euo pipefail

repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
env_prefix="${ENV_PREFIX:-$repo/pipelines/single_cell/work/env}"
data="${SINGLE_CELL_DATA:-$repo/data/single_cell/raw}"
final="${SINGLE_CELL_FINAL:-$repo/data/single_cell/final}"

[[ -x "$env_prefix/bin/Rscript" ]] || { echo "ERROR: run setup.sh first" >&2; exit 1; }
export PATH="$env_prefix/bin:$PATH"
cd "$repo"
uv run --quiet python pipelines/single_cell/src/run_selected.py --data "$data" --final "$final" \
  --gtf "${SINGLE_CELL_GTF:-$data/refs/Homo_sapiens.GRCh38.115.gtf.gz}" --genes metadata/benchmark_genes.tsv.gz \
  --workers "${RUN_WORKERS:-8}" "$@"
