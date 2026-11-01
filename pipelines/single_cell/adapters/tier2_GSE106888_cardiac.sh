#!/usr/bin/env bash
# GSE106888: wells through STAR, tdTomato reads on the mixed plates, store, one smoke fit, all fits, outcome row. Steps already done are kept.
set -euo pipefail
A=/SCRATCH/dtzim01/funmirbench/workspaces/single-cell-processing/pipelines/single_cell/adapters
W=/SCRATCH/dtzim01/funmirbench/workspaces/single-cell-processing
PY=$W/.venv/bin/python
D=/SCRATCH/dtzim01/single_cell_mirna_data/GSE106888
F=/SCRATCH/dtzim01/fun-experiment-data/single-cell-final/20261031/GSE106888_cardiac_v2
export PATH=/SCRATCH/dtzim01/fun-experiment-data/tools/envs/deseq2-shrinkage/bin:$PATH OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
"$PY" -I "$A/GSE106888_cardiac.py" prepare
echo "PROGRESS 0/4 STAR"
[ -s "$D/aligned/wells/Solo.out/Gene/raw/matrix.mtx" ] || bash "$A/tier2_star_plate.sh" "$D/aligned/wells.manifest" "$D/aligned/wells" 20
echo "PROGRESS 1/4 markers"
for p in D3M3 D5M1 D7M1 D9M1; do [ -s "$D/aligned/marker_$p/counts.tsv" ] || "$PY" -I "$A/tier2_marker_counts.py" "$D/aligned/${p}_tdtomato.manifest" "$D/marker/AY678269.fa" "$D/aligned/marker_$p/counts.tsv" 12 --mismatches 3; done
echo "PROGRESS 2/4 store"
T=/SCRATCH/dtzim01/tmp/tier2/stores; mkdir -p "$T"
"$PY" -I "$A/GSE106888_cardiac.py" build "$T/GSE106888_cardiac.h5ad"
echo "PROGRESS 3/4 fits"
bash "$A/tier2_fit.sh" GSE106888_cardiac "$T/GSE106888_cardiac.h5ad" "$F" 3
"$PY" -I "$A/tier2_record.py" GSE106888_cardiac other "$D" done "hMGT133 is miR-133 plus three factors, so a comparison is not a miR-133 effect alone" --final-dir "$F" --store "$F/store/GSE106888_cardiac.h5ad" --solo "=$D/aligned/wells/Solo.out/Gene/raw" 
echo "PROGRESS 4/4 done"
