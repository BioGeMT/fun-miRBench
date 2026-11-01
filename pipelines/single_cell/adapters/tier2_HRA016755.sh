#!/usr/bin/env bash
# HRA016755: STARsolo on the GEX runs, guide barcodes from the GBC runs, two stores, smoke and fits, outcome rows. No ATAC runs are used.
set -euo pipefail
A=/SCRATCH/dtzim01/funmirbench/workspaces/single-cell-processing/pipelines/single_cell/adapters
PY=/SCRATCH/dtzim01/funmirbench/workspaces/single-cell-processing/.venv/bin/python
D=/SCRATCH/dtzim01/single_cell_mirna_data/HRA016755
WL=/SCRATCH/dtzim01/single_cell_mirna_data/tier2_reference/gex_737K-arc-v1.txt
T=/SCRATCH/dtzim01/tmp/tier2/stores; mkdir -p "$T"
export PATH=/SCRATCH/dtzim01/fun-experiment-data/tools/envs/deseq2-shrinkage/bin:$PATH OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
list(){ local out=""; for r in "$@"; do out="$out,$D/fastq/$r"; done; echo "${out#,}"; }
star(){ [ -s "$D/aligned/$1/Solo.out/Gene/filtered/matrix.mtx" ] || bash "$A/tier2_star_10x.sh" "$D/aligned/$1" 20 "$WL" 12 "$2" "$3"; }
echo "PROGRESS 0/5 K562 rep1"
star K562_rep1 "$(list HRR3070148_f1.fastq.gz)" "$(list HRR3070148_r2.fastq.gz)"
echo "PROGRESS 1/5 K562 rep2"
star K562_rep2 "$(list HRR3070152_f1.fastq.gz HRR3070153_f1.fastq.gz)" "$(list HRR3070152_r2.fastq.gz HRR3070153_r2.fastq.gz)"
echo "PROGRESS 2/5 BJ"
star BJ "$(list HRR3070156_f1.fastq.gz HRR3070157_f1.fastq.gz HRR3070158_f1.fastq.gz HRR3070159_f1.fastq.gz)" "$(list HRR3070156_r2.fastq.gz HRR3070157_r2.fastq.gz HRR3070158_r2.fastq.gz HRR3070159_r2.fastq.gz)"
echo "PROGRESS 3/5 guide barcodes"
for library in K562_rep1 K562_rep2 BJ; do [ -s "$D/aligned/$library/gbc_umi_reads.tsv.gz" ] || "$PY" -I "$A/HRA016755.py" gbc "$library"; done
echo "PROGRESS 4/5 stores and fits"
for line in K562 BJ; do
  ds=HRA016755_$line; F=/SCRATCH/dtzim01/fun-experiment-data/single-cell-final/20261031/${ds}_v2
  "$PY" -I "$A/HRA016755.py" build "$line" "$T/$ds.h5ad"
  bash "$A/tier2_fit.sh" "$ds" "$T/$ds.h5ad" "$F" 2
  if [ "$line" = K562 ]; then SOLO=(--solo "K562_rep1=$D/aligned/K562_rep1/Solo.out/Gene/filtered" --solo "K562_rep2=$D/aligned/K562_rep2/Solo.out/Gene/filtered"); else SOLO=(--solo "BJ=$D/aligned/BJ/Solo.out/Gene/filtered"); fi
  mkdir -p "$F/tables"; for l in K562_rep1 K562_rep2 BJ; do case $l in $line*) cp "$D/aligned/$l/gbc_calls.tsv.gz" "$F/tables/${l}_gbc_calls.tsv.gz";; esac; done
  "$PY" -I "$A/tier2_record.py" "$ds" host_gene "$D" done "guide calls reprocessed from the GBC runs with chimeras resolved, no response-based filter; ATAC runs not downloaded" --final-dir "$F" --store "$F/store/$ds.h5ad" "${SOLO[@]}"
done
echo "PROGRESS 5/5 done"
