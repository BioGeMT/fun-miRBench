#!/usr/bin/env bash
# GSE154900: STARsolo per replicate (10x 3' v2; barcode 16, UMI 10, read 1 is 28 nt in rep1 and 101 nt in rep2), store, check, outcome row.
set -euo pipefail
A=/SCRATCH/dtzim01/funmirbench/workspaces/single-cell-processing/pipelines/single_cell/adapters
PY=/SCRATCH/dtzim01/funmirbench/workspaces/single-cell-processing/.venv/bin/python
D=/SCRATCH/dtzim01/single_cell_mirna_data/GSE154900
WL=/SCRATCH/dtzim01/single_cell_mirna_data/SRP150046/reference/737K-august-2016.txt
F=/SCRATCH/dtzim01/fun-experiment-data/single-cell-final/20261031/GSE154900_virus
echo "PROGRESS 0/3 rep1"
bash "$A/tier2_star_10x.sh" "$D/aligned/rep1" 20 "$WL" 10 "$D/fastq/SRR12287751_1.fastq.gz,$D/fastq/SRR12287752_1.fastq.gz" "$D/fastq/SRR12287751_2.fastq.gz,$D/fastq/SRR12287752_2.fastq.gz"
echo "PROGRESS 1/3 rep2"
bash "$A/tier2_star_10x.sh" "$D/aligned/rep2" 20 "$WL" 10 "$D/fastq/SRR12287753_1.fastq.gz,$D/fastq/SRR12287754_1.fastq.gz" "$D/fastq/SRR12287753_2.fastq.gz,$D/fastq/SRR12287754_2.fastq.gz"
echo "PROGRESS 2/3 store"
[ ! -e "$F" ] || { echo "$F exists"; exit 1; }
mkdir -p "$F/store"
"$PY" -I "$A/GSE154900_virus.py" --out "$F/store/GSE154900_virus.h5ad"
"$PY" -I "$A/tier2_record.py" GSE154900_virus other "$D" done "store only, no comparison: no sample-level arm; latent and lytic cells share one culture, and Ensembl 115 has no viral genes" --final-dir "$F" --store "$F/store/GSE154900_virus.h5ad" --solo "rep1=$D/aligned/rep1/Solo.out/Gene/raw" --solo "rep2=$D/aligned/rep2/Solo.out/Gene/raw"
echo "PROGRESS 3/3 done"
