#!/usr/bin/env bash
# Download the Ensembl 115 GTF and the source files of the selected experiments. Complete files are kept, so the script
# can be rerun after an interruption. Arguments go to download_sources.py, for example --groups main.
set -euo pipefail

repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
env_prefix="${ENV_PREFIX:-$repo/pipelines/single_cell/work/env}"
data="${SINGLE_CELL_DATA:-$repo/data/single_cell/raw}"
gtf="$data/refs/Homo_sapiens.GRCh38.115.gtf.gz"
ensembl="https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens"

[[ -x "$env_prefix/bin/Rscript" ]] || { echo "ERROR: run setup.sh first" >&2; exit 1; }
export PATH="$env_prefix/bin:$PATH"
mkdir -p "$(dirname "$gtf")"

echo "[1/2] GTF..."
# Ensembl CHECKSUMS hold the BSD sum (checksum and 1 kB blocks) of each .gz file.
if [[ ! -s "$gtf" ]]; then
  curl -fsSL --retry 5 -o "$gtf.part" "$ensembl/Homo_sapiens.GRCh38.115.gtf.gz"
  expected="$(curl -fsSL --retry 5 "$ensembl/CHECKSUMS" | awk '$3 == "Homo_sapiens.GRCh38.115.gtf.gz" {print $1 + 0, $2 + 0}')"
  actual="$(sum "$gtf.part" | awk '{print $1 + 0, $2 + 0}')"
  [[ -n "$expected" && "$expected" == "$actual" ]] || { echo "ERROR: checksum mismatch for the GTF" >&2; exit 1; }
  mv "$gtf.part" "$gtf"
fi

echo "[2/2] Source files..."
cd "$repo"
uv run --quiet python pipelines/single_cell/src/download_sources.py --data "$data" "$@"
