#!/usr/bin/env bash
# Plate data: one STARsolo SmartSeq process over a manifest (read 1, read 2 or -, cell), Ensembl 115 index, no BAM.
# Counts reads per cell in the exons of one gene: uniquely mapping reads only, no multi-mapper fractions, no UMI or
# position deduplication, so the counts are integers. Waits until at least 100 GB of memory are available.
# Usage: tier2_star_plate.sh <manifest.tsv> <out_dir> <threads> [Unstranded|Forward|Reverse]
set -euo pipefail
STAR=/SCRATCH/dtzim01/tools/envs/eisar/bin/STAR
GENOME=/SCRATCH/dtzim01/fun-experiment-data/references/ensembl_v115/star_index
MANIFEST=$1; OUT=$2; THREADS=$3; STRAND=${4:-Unstranded}
[ -s "$MANIFEST" ] || { echo "empty manifest $MANIFEST" >&2; exit 1; }
[ ! -e "$OUT/Solo.out" ] || { echo "$OUT already holds an alignment" >&2; exit 1; }
mkdir -p "$OUT"
until [ "$(free -g | awk '/^Mem:/{print $7}')" -ge 100 ]; do echo "waiting for memory $(date +%T)"; sleep 120; done
ulimit -n 65000 2>/dev/null || true
"$STAR" --runThreadN "$THREADS" --genomeDir "$GENOME" --readFilesManifest "$MANIFEST" --readFilesCommand zcat \
    --soloType SmartSeq --soloUMIdedup NoDedup --soloStrand "$STRAND" --soloFeatures Gene --soloMultiMappers Unique \
    --outSAMtype None --outFileNamePrefix "$OUT/" --outTmpDir "$OUT/_STARtmp"
[ -s "$OUT/Solo.out/Gene/raw/matrix.mtx" ]
