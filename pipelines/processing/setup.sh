#!/usr/bin/env bash
# Create the conda environment of pipelines/processing/environment.yml and the Python environment, and check the tools.
set -euo pipefail

repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
env_prefix="${ENV_PREFIX:-$repo/pipelines/processing/work/env}"
env_file="$repo/pipelines/processing/environment.yml"

command -v conda >/dev/null || { echo "ERROR: conda is required" >&2; exit 1; }

echo "[1/4] Conda environment in $env_prefix..."
if [[ ! -d "$env_prefix/conda-meta" ]]; then
  conda env create --quiet --prefix "$env_prefix" --file "$env_file"
fi
export PATH="$env_prefix/bin:$PATH"

echo "[2/4] Python environment..."
(cd "$repo" && uv sync --quiet)

echo "[3/4] Checking tools..."
for tool in fastqc fastp STAR featureCounts Rscript uv; do
  command -v "$tool" >/dev/null || {
    echo "ERROR: $tool is missing; run: conda env update --prefix $env_prefix --file $env_file" >&2; exit 1; }
done
Rscript -e 'stopifnot(requireNamespace("DESeq2", quietly = TRUE), requireNamespace("apeglm", quietly = TRUE))' \
  >/dev/null 2>&1 || {
  echo "ERROR: R needs DESeq2 and apeglm; run: conda env update --prefix $env_prefix --file $env_file" >&2; exit 1; }
echo "[4/4] miRBase names the pipeline checks configs against..."
mirbase="$repo/data/experiments/raw/refs/mirbase"
mkdir -p "$mirbase"
[[ -f "$mirbase/mature.fa" ]] || curl -fsSL -o "$mirbase/mature.fa" https://mirbase.org/download_version_files/22.1/mature.fa
awk '/^>hsa-/ {sub(/^>/, ""); print $1}' "$mirbase/mature.fa" | LC_ALL=C sort -u >"$mirbase/hsa_mature_mirnas.txt"
echo "Ready."
