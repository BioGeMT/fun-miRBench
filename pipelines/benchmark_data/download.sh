#!/usr/bin/env bash
# Download the Ensembl 115 genome and GTF and the FASTQ files of the 21 OE test comparisons. Complete files are kept,
# so the script can be rerun after an interruption. Arguments go to download_reads.py (for example --study GSE124411).
set -euo pipefail

repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
env_prefix="${ENV_PREFIX:-$repo/pipelines/benchmark_data/work/env}"
refs="$repo/data/experiments/raw/refs/ensembl_v115"
ensembl="https://ftp.ensembl.org/pub/release-115"

[[ -x "$env_prefix/bin/uv" ]] || { echo "ERROR: run setup.sh first" >&2; exit 1; }
export PATH="$env_prefix/bin:$PATH"

# Ensembl CHECKSUMS hold the BSD sum (checksum and 1 kB blocks) of each .gz file.
fetch_reference() {
  local url="$1" source="$2" name="$3" expected actual
  if [[ -s "$refs/$name" ]]; then echo "ok      $name"; return; fi
  curl -fsSL --retry 5 -o "$refs/$name.gz.part" "$url/$source.gz"
  expected="$(curl -fsSL --retry 5 "$url/CHECKSUMS" | awk -v f="$source.gz" '$3 == f {print $1 + 0, $2 + 0}')"
  actual="$(sum "$refs/$name.gz.part" | awk '{print $1 + 0, $2 + 0}')"
  [[ -n "$expected" && "$expected" == "$actual" ]] || { echo "ERROR: checksum mismatch for $source.gz" >&2; exit 1; }
  gunzip -c "$refs/$name.gz.part" >"$refs/$name.part"
  mv "$refs/$name.part" "$refs/$name"
  rm "$refs/$name.gz.part"
  echo "fetched $name"
}

echo "[1/2] References and the miRBase names the pipeline checks configs against..."
mkdir -p "$refs"
# The reruns used the soft-masked genome under the name the configs point to.
fetch_reference "$ensembl/fasta/homo_sapiens/dna" Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
  Homo_sapiens.GRCh38.dna.primary_assembly.fa
fetch_reference "$ensembl/gtf/homo_sapiens" Homo_sapiens.GRCh38.115.gtf Homo_sapiens.GRCh38.115.gtf
cd "$repo"
uv run --quiet funmirbench-experiments-download-examples --example mirbase-hsa-mature --log-level WARNING

echo "[2/2] Reads..."
uv run --quiet python pipelines/benchmark_data/src/download_reads.py "$@"
