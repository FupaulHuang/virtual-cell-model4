#!/usr/bin/env bash
# Invoke only after selecting the candidate when REFERENCE is the reserved Test.
set -euo pipefail
cd "$(dirname "$0")/.."
PREDICTION="${1:?prediction H5AD required}"
REFERENCE="${2:?reference H5AD required}"
OUTDIR="${3:?metric output directory required}"
mkdir -p .cache/numba .cache/matplotlib "$OUTDIR"
export NUMBA_CACHE_DIR="$PWD/.cache/numba"
export MPLCONFIGDIR="$PWD/.cache/matplotlib"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=2 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=2
LOCAL_INPUT="$OUTDIR/local_scoring_input.h5ad"
.venv/bin/python scripts/prepare_metric_input.py \
  "$PREDICTION" --reference "$REFERENCE" --output "$LOCAL_INPUT"
.venv/bin/cell-eval2 run \
  -ap "$LOCAL_INPUT" -ar "$REFERENCE" --write-degenes \
  --preset vcc2026 --pert-col target_gene --set device=cpu --set de.backend=scanpy \
  --set control_source=real --set num_threads=2 -o "$OUTDIR"
