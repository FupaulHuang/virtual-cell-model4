#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p .cache/numba .cache/matplotlib evaluations/h1_validation
export NUMBA_CACHE_DIR="$PWD/.cache/numba"
export MPLCONFIGDIR="$PWD/.cache/matplotlib"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=2 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=2
CELL_EVAL2=.venv/bin/cell-eval2
EVAL=evaluations/h1_validation
COMMON=(--preset vcc2026 --pert-col target_gene --set device=cpu --set de.backend=scanpy --set control_source=real --set num_threads=2)
"$CELL_EVAL2" baseline -ar "$EVAL/reference.h5ad" "${COMMON[@]}" \
  -o "$EVAL/baseline_cpu" --save-pred "$EVAL/baseline_cpu_prediction.h5ad"
"$CELL_EVAL2" prep-real-bundle --real "$EVAL/reference.h5ad" \
  --baseline "$EVAL/baseline_cpu_prediction.h5ad" "${COMMON[@]}" \
  --anchor-splits 5 -o "$EVAL/real_bundle_cpu"
