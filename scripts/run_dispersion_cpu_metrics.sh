#!/usr/bin/env bash
# Same frozen cohort, preset, backend and threads for all experimental arms.
set -euo pipefail
cd "$(dirname "$0")/.."
ARM="${1:?arm name required}"
PREDICTION="${2:?prediction H5AD required}"
OUTDIR="evaluations/dispersion_h1_training_400/${ARM}_cpu"
REFERENCE="evaluations/dispersion_h1_training_400/reference.h5ad"
mkdir -p "$OUTDIR" .cache/numba .cache/matplotlib
export NUMBA_CACHE_DIR="$PWD/.cache/numba"
export MPLCONFIGDIR="$PWD/.cache/matplotlib"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=4 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=4
LOCAL_INPUT="$OUTDIR/local_scoring_input.h5ad"
if [[ ! -f "$LOCAL_INPUT" ]]; then
  .venv/bin/python scripts/prepare_metric_input.py \
    "$PREDICTION" --reference "$REFERENCE" --output "$LOCAL_INPUT"
fi
EXTRA=()
if [[ "$ARM" != incumbent ]]; then
  REAL_DE="evaluations/dispersion_h1_training_400/incumbent_cpu/de_real.parquet"
  test -f "$REAL_DE"
  EXTRA=(--de-real "$REAL_DE")
fi
.venv/bin/cell-eval2 run \
  -ap "$LOCAL_INPUT" -ar "$REFERENCE" --write-degenes \
  --preset vcc2026 --pert-col target_gene --set device=cpu --set de.backend=scanpy \
  --set control_source=real --set num_threads=4 "${EXTRA[@]}" -o "$OUTDIR"
