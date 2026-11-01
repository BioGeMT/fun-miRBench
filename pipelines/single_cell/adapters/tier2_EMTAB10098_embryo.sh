#!/usr/bin/env bash
# E-MTAB-10098: miRNA counts from the small RNA reads, stores from them and from the deposited mRNA table, outcome row.
set -euo pipefail
A=/SCRATCH/dtzim01/funmirbench/workspaces/single-cell-processing/pipelines/single_cell/adapters
PY=/SCRATCH/dtzim01/funmirbench/workspaces/single-cell-processing/.venv/bin/python
D=/SCRATCH/dtzim01/single_cell_mirna_data/E-MTAB-10098
F=/SCRATCH/dtzim01/fun-experiment-data/single-cell-final/20261031/EMTAB10098_embryo
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
"$PY" -I "$A/EMTAB10098_embryo.py" prepare
echo "PROGRESS 0/2 small RNA"
"$PY" -I "$A/tier2_smallrna.py" "$D/aligned/small.manifest" "$D/aligned/small" 12 --adapter TGGAATTCTCGGGTGCCAAGG --umi-5p 8 --random-5p 2
echo "PROGRESS 1/2 stores"
[ ! -e "$F" ] || { echo "$F exists"; exit 1; }
"$PY" -I "$A/EMTAB10098_embryo.py" build "$F"
"$PY" -I "$A/tier2_record.py" EMTAB10098_embryo other "$D" done "stores only, no comparison: no miRNA perturbation (dexamethasone against none), small RNA and mRNA from different cells of the same embryos; embryo-level sums in tables/" --final-dir "$F" --store "$F/store/EMTAB10098_embryo_mirna.h5ad"
echo "PROGRESS 2/2 done"
