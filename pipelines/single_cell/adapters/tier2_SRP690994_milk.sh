#!/usr/bin/env bash
# SRP690994: STARsolo per sample, one at a time; miRNA counts per sample (five runs deposited trimmed, two with the adapter); store; outcome row.
set -euo pipefail
A=/SCRATCH/dtzim01/funmirbench/workspaces/single-cell-processing/pipelines/single_cell/adapters
PY=/SCRATCH/dtzim01/funmirbench/workspaces/single-cell-processing/.venv/bin/python
D=/SCRATCH/dtzim01/single_cell_mirna_data/SRP690994
F=/SCRATCH/dtzim01/fun-experiment-data/single-cell-final/20261031/SRP690994_milk
WL=/SCRATCH/dtzim01/single_cell_mirna_data/tier2_reference/3M-february-2018.txt
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
"$PY" -I "$A/SRP690994_milk.py" prepare
n=0
while IFS=$'\t' read -r sample r1 r2; do
  echo "PROGRESS $n/9 $sample"
  [ -s "$D/aligned/$sample/Solo.out/Gene/filtered/matrix.mtx" ] || bash "$A/tier2_star_10x.sh" "$D/aligned/$sample" 20 "$WL" 12 "$r1" "$r2"
  n=$((n + 1))
done < "$D/aligned/samples.tsv"
echo "PROGRESS 7/9 miRNA"
"$PY" -I "$A/tier2_smallrna.py" "$D/aligned/mirna_trimmed.manifest" "$D/aligned/mirna_trimmed" 5 --adapter none
"$PY" -I "$A/tier2_smallrna.py" "$D/aligned/mirna_adapter.manifest" "$D/aligned/mirna_adapter" 2 --adapter AGATCGGAAGAGCACACGTCT
echo "PROGRESS 8/9 store"
[ ! -e "$F" ] || { echo "$F exists"; exit 1; }
"$PY" -I "$A/SRP690994_milk.py" build "$F"
SOLO=(); while IFS=$'\t' read -r sample r1 r2; do SOLO+=(--solo "$sample=$D/aligned/$sample/Solo.out/Gene/filtered"); done < "$D/aligned/samples.tsv"
"$PY" -I "$A/tier2_record.py" SRP690994_milk other "$D" done "store only, no comparison: no perturbation arm, gestational age only; miRNA is exosomal and per sample, so sample-level only" --final-dir "$F" --store "$F/store/SRP690994_milk.h5ad" "${SOLO[@]}"
echo "PROGRESS 9/9 done"
