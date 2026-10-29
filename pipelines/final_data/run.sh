#!/usr/bin/env bash
# Run the experiments pipeline for each test comparison, then build the final per-gene tables.
# Arguments: the comparisons to run (default: all 30). Finished runs are skipped; the final tables are built once
# every comparison has a finished run.
set -euo pipefail

repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
env_prefix="${ENV_PREFIX:-$repo/pipelines/final_data/work/env}"
workers="${RUN_WORKERS:-2}"
final_dir="${FINAL_DIR:-$repo/results/final_data}"
work="$repo/pipelines/final_data/work"

[[ -x "$env_prefix/bin/uv" ]] || { echo "ERROR: run setup.sh first" >&2; exit 1; }
export PATH="$env_prefix/bin:$PATH" work
cd "$repo"
mkdir -p "$work/logs"

# The newest run of a comparison that wrote its fitted DESeq2 object.
latest_fit() {
  local dir fit=""
  for dir in pipelines/experiments/runs/????????_??????_"$1"; do
    if [[ -f "$dir/complete_analysis.rds" ]]; then fit="$dir/complete_analysis.rds"; fi
  done
  echo "$fit"
}

run_one() {
  if [[ -n "$(latest_fit "$1")" ]]; then echo "done    $1"; return; fi
  echo "running $1"
  uv run --quiet funmirbench-experiments --config "pipelines/final_data/configs/$1.yaml" >"$work/logs/$1.log" 2>&1 \
    || { echo "FAILED  $1 (log: $work/logs/$1.log)" >&2; return 1; }
  echo "ran     $1"
}
export -f latest_fit run_one

mapfile -t all_ids < <(awk -F'\t' 'NR > 1 {print $1}' metadata/comparisons.tsv)
ids=("$@")
if (( ${#ids[@]} == 0 )); then ids=("${all_ids[@]}"); fi
for id in "${ids[@]}"; do
  [[ -f "pipelines/final_data/configs/$id.yaml" ]] || { echo "ERROR: no config for $id" >&2; exit 1; }
done

echo "[1/3] Pipeline runs, $workers at a time..."
# The first run builds the shared STAR index when it is missing, so it runs alone.
run_one "${ids[0]}"
printf '%s\n' "${ids[@]:1}" | xargs -r -P "$workers" -I{} bash -c 'run_one "$1"' _ {}

echo "[2/3] Fit list..."
fits="$work/fits.tsv"
printf 'experiment_id\tperturbation\tfit_rds\n' >"$fits"
missing=0
while IFS=$'\t' read -r id perturbation; do
  fit="$(latest_fit "$id")"
  if [[ -z "$fit" ]]; then missing=$((missing + 1)); continue; fi
  printf '%s\t%s\t%s\n' "$id" "$perturbation" "$repo/$fit" >>"$fits"
done < <(awk -F'\t' 'NR == 1 {for (i = 1; i <= NF; i++) col[$i] = i; next}
                     {print $col["experiment_id"] "\t" $col["perturbation"]}' metadata/comparisons.tsv)
if (( missing > 0 )); then
  echo "$missing of ${#all_ids[@]} comparisons have no finished run yet; the final tables need all of them."
  exit 0
fi

echo "[3/3] Final tables in $final_dir..."
[[ ! -e "$final_dir" ]] || { echo "ERROR: $final_dir exists; move it or set FINAL_DIR" >&2; exit 1; }
uv run --quiet python pipelines/final_data/src/build_final_data.py --fits "$fits" \
  --genes metadata/benchmark_genes.tsv.gz --out "$final_dir" --rscript "$env_prefix/bin/Rscript"
