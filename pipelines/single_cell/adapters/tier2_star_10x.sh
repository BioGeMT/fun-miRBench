#!/usr/bin/env bash
# Droplet data: STARsolo against the Ensembl 115 index, the Cell Ranger-like settings of human_benchmark/source/align_srp150046.sh, no BAM.
# Counts UMIs in the exons of one gene, unique mappers only, so counts are integers. Waits until 100 GB of memory are available.
# Usage: tier2_star_10x.sh <out_dir> <threads> <whitelist> <umi_len> <read 1 files, comma-separated> <read 2 files, comma-separated>
set -euo pipefail
STAR=/SCRATCH/dtzim01/tools/envs/eisar/bin/STAR
GENOME=/SCRATCH/dtzim01/fun-experiment-data/references/ensembl_v115/star_index
OUT=$1; THREADS=$2; WHITELIST=$3; UMI_LEN=$4; R1=$5; R2=$6
[ ! -e "$OUT/Solo.out" ] || { echo "$OUT already holds an alignment" >&2; exit 1; }
mkdir -p "$OUT"
until [ "$(free -g | awk '/^Mem:/{print $7}')" -ge 100 ]; do echo "waiting for memory $(date +%T)"; sleep 120; done
"$STAR" --runThreadN "$THREADS" --genomeDir "$GENOME" --readFilesIn "$R2" "$R1" --readFilesCommand zcat \
    --soloType CB_UMI_Simple --soloCBwhitelist "$WHITELIST" --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 \
    --soloUMIlen "$UMI_LEN" --soloBarcodeReadLength 0 --soloFeatures Gene --soloMultiMappers Unique \
    --soloCellFilter EmptyDrops_CR --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR \
    --soloUMIdedup 1MM_CR --clipAdapterType CellRanger4 --outFilterScoreMin 30 --outSAMtype None \
    --outFileNamePrefix "$OUT/" --outTmpDir "$OUT/_STARtmp"
[ -s "$OUT/Solo.out/Gene/filtered/matrix.mtx" ]
