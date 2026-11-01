#!/usr/bin/env bash
# One line per study: run its adapter, two at a time (each up to 4 threads). Writes units.tsv in the run directory.
set -uo pipefail
A=/SCRATCH/dtzim01/funmirbench/workspaces/single-cell-processing/pipelines/single_cell/adapters
PY=/SCRATCH/dtzim01/funmirbench/workspaces/single-cell-processing/.venv/bin/python
FINAL=/SCRATCH/dtzim01/fun-experiment-data/single-cell-final/20261031/cooccurrence_library_total
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 RUN=$PWD
one(){ d=$1; script=$2; shift 2; log=$RUN/$d.log
  echo "START $d $(date +%T)"
  if "$PY" -I "$A/$script.py" "$@" > "$log" 2>&1 && [ -s "$FINAL/$d/correlation/groups.tsv" ]; then printf '%s\tok\t\n' "$d" >> "$RUN/units.tsv"; echo "OK $d $(date +%T)"
  else printf '%s\tfailed\t%s\n' "$d" "$(tail -1 "$log" | tr '\t' ' ' | cut -c1-200)" >> "$RUN/units.tsv"; echo "FAILED $d $(date +%T)"; fi; }
export -f one; export A PY FINAL
printf 'unit\tstatus\tnote\n' > units.tsv
cat "$1" | xargs -P 2 -L 1 bash -c 'one "$@"' _
