#!/usr/bin/env bash
# PRJNA994505: one STAR process over the 643 single cells, two stores, smoke and fits, outcome rows. Calling agoTRIBE editing events is out of scope.
set -euo pipefail
A=/SCRATCH/dtzim01/funmirbench/workspaces/single-cell-processing/pipelines/single_cell/adapters
PY=/SCRATCH/dtzim01/funmirbench/workspaces/single-cell-processing/.venv/bin/python
D=/SCRATCH/dtzim01/single_cell_mirna_data/PRJNA994505
T=/SCRATCH/dtzim01/tmp/tier2/stores; mkdir -p "$T"
export PATH=/SCRATCH/dtzim01/fun-experiment-data/tools/envs/deseq2-shrinkage/bin:$PATH OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
echo "PROGRESS 0/3 STAR"
bash "$A/tier2_star_plate.sh" "$D/aligned/cells.manifest" "$D/aligned/cells" 20
for line in HEK293T K562; do
  ds=PRJNA994505_$line; F=/SCRATCH/dtzim01/fun-experiment-data/single-cell-final/20261031/$ds
  "$PY" -I "$A/PRJNA994505.py" build "$line" "$T/$ds.h5ad"
  bash "$A/tier2_fit.sh" "$ds" "$T/$ds.h5ad" "$F" 2
  "$PY" -I "$A/tier2_record.py" "$ds" other "$D" done "expression comparisons only, calling agoTRIBE editing events is out of scope; read counts, not UMIs, since ENA holds only read 2" --final-dir "$F" --store "$F/store/$ds.h5ad" --solo "=$D/aligned/cells/Solo.out/Gene/raw"
  echo "PROGRESS done $line"
done
