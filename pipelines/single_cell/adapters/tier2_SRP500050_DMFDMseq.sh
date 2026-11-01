#!/usr/bin/env bash
# SRP500050: align mRNA (16 threads), count miRNAs (8 files at a time), build the study. Reads in single_cell_mirna_data/SRP500050/fastq.
set -euo pipefail
A=/SCRATCH/dtzim01/funmirbench/workspaces/single-cell-processing/pipelines/single_cell/adapters
PY=/SCRATCH/dtzim01/funmirbench/workspaces/single-cell-processing/.venv/bin/python
D=/SCRATCH/dtzim01/single_cell_mirna_data/SRP500050/aligned
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4
"$PY" -I "$A/SRP500050_DMFDMseq.py" prepare
echo "PROGRESS 0/3 aligning mRNA"
bash "$A/tier2_star_plate.sh" "$D/mrna.manifest" "$D/mrna" 16
echo "PROGRESS 1/3 counting miRNAs"
"$PY" -I "$A/tier2_smallrna.py" "$D/mirna.manifest" "$D/mirna" 8 --adapter TGGAATTCTCGGGTGCCAAGG
echo "PROGRESS 2/3 building the study"
"$PY" -I "$A/SRP500050_DMFDMseq.py" build
echo "PROGRESS 3/3 done"
