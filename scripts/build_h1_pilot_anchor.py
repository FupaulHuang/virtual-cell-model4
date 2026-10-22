#!/usr/bin/env python3
"""One-split diagnostic calibration using official APIs; not the five-split result."""
import json
from pathlib import Path
import time
from cell_eval2.anchor import compute_replicate_anchor, build_meta, write_anchor, _derive_seeds
from cell_eval2.catalog import resolve_metrics
from cell_eval2.config import EvalConfig
from cell_eval2.io import load_anndata
from cell_eval2.run import metric_output_names

ROOT = Path(__file__).resolve().parents[1]
output = ROOT / 'evaluations/h1_validation/pilot_anchor1'
if output.exists():
    raise FileExistsError(output)
cfg = EvalConfig.from_yaml(str(ROOT / 'evaluations/h1_validation/pie_dual_cpu/run_params.yaml'))
assert cfg.device == 'cpu' and cfg.de.backend == 'scanpy' and cfg.control_source == 'real'
cfg.outdir = str(output)
real = load_anndata(ROOT / 'evaluations/h1_validation/reference.h5ad', backed=False)
started = time.monotonic()
splits, anchor = compute_replicate_anchor(real, config=cfg, base_seed=0, n_splits=1)
names, _ = resolve_metrics(cfg.metrics, version=cfg.version)
meta = build_meta(real_ad=real, cfg=cfg, names=list(names), base_seed=0, n_splits=1,
                  seeds=_derive_seeds(0, 1), metrics=metric_output_names(cfg))
meta['diagnostic_label'] = 'PILOT ONE SPLIT ONLY: provisional calibration; five-split job remains independent.'
write_anchor(output, splits, anchor, meta=meta)
(output / 'pilot_status.json').write_text(json.dumps({'status': 'complete', 'n_splits': 1,
    'purpose': 'Diagnostic provisional selection weights only; not five-split calibration',
    'elapsed_seconds': time.monotonic() - started}, indent=2) + '\n')
print(str(output), flush=True)
