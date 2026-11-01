#!/usr/bin/env bash
# SRP067502: one STAR process over the 597 cells, GFP reads, store, smoke and fits, outcome row. An alignment already done is kept.
set -euo pipefail
A=/SCRATCH/dtzim01/funmirbench/workspaces/single-cell-processing/pipelines/single_cell/adapters
PY=/SCRATCH/dtzim01/funmirbench/workspaces/single-cell-processing/.venv/bin/python
D=/SCRATCH/dtzim01/single_cell_mirna_data/SRP067502
F=/SCRATCH/dtzim01/fun-experiment-data/single-cell-final/20261031/SRP067502_C1_v2
T=/SCRATCH/dtzim01/tmp/tier2/stores; mkdir -p "$T"
export PATH=/SCRATCH/dtzim01/fun-experiment-data/tools/envs/deseq2-shrinkage/bin:$PATH OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
"$PY" -I "$A/SRP067502_C1.py" prepare
echo "PROGRESS 0/3 STAR"
[ -s "$D/aligned/cells/Solo.out/Gene/raw/matrix.mtx" ] || bash "$A/tier2_star_plate.sh" "$D/aligned/cells.manifest" "$D/aligned/cells" 20
echo "PROGRESS 1/3 GFP"
"$PY" -I "$A/tier2_marker_counts.py" "$D/aligned/markers.manifest" /SCRATCH/dtzim01/single_cell_mirna_data/SRP150046/reference/egfp.fa "$D/aligned/marker_min30/counts.tsv" 16 --mismatches 2 --min-length 30
echo "PROGRESS 2/3 store and fits"
"$PY" -I "$A/SRP067502_C1.py" build "$T/SRP067502_C1.h5ad"
bash "$A/tier2_fit.sh" SRP067502_C1 "$T/SRP067502_C1.h5ad" "$F" 2
"$PY" -I "$A/tier2_record.py" SRP067502_C1 raise "$D" done "arms from GFP per the preprint's methods, chips as replicate units; no dose per cell or chip: Table S1 names chips a to g, the deposit by date, and neither links them" --final-dir "$F" --store "$F/store/SRP067502_C1.h5ad" --solo "=$D/aligned/cells/Solo.out/Gene/raw"
echo "PROGRESS 3/3 done"
