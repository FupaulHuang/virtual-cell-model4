#!/usr/bin/env python3
"""Verify the structural experiment and score against the frozen H1 reference."""
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(key, '1')
import argparse
from dataclasses import asdict, replace
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess
import yaml
from cell_eval2 import EvalConfig
from build_prediction import ROOT, sha256
from run_bulk_cpu_metrics import verify_local_copy
from run_verified_dispersion_metrics import baseline_checks, frozen_files
from validate_prediction_receipt import validate
from summarize_dispersion_experiment import real_de_identity


def read(path): return json.loads(Path(path).read_text())


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('prediction', type=Path)
    args = ap.parse_args()
    prediction = args.prediction.resolve()
    frozen = read(ROOT/'reports/redesign_implementation_freeze.json')
    for name, digest in frozen['files'].items():
        if sha256(ROOT/name) != digest: raise ValueError('Frozen file changed: ' + name)
    frozen_files()
    baseline = baseline_checks()
    cache = read(ROOT/'reports/dispersion_real_de_cache.json')
    if any(cache.get(k) != v for k, v in baseline.items()): raise ValueError('Original real-DE/cache changed')
    receipt = read(prediction.with_suffix('.receipt.json'))
    spec = receipt['redesign']
    name = f'{spec["response"]}_{spec["emitter"]}_s{spec["seed"]}'
    if (spec['engineering_subset'] or spec['seed'] not in frozen['seeds']
            or spec['response'] not in ['current', 'adapter'] or spec['emitter'] not in ['current', 'preserving']
            or prediction != ROOT/'outputs/redesign_h1_400'/name/'prediction.h5ad'
            or spec['incumbent_sha256'] != baseline['prediction_sha256']
            or spec['profiles_sha256'] != frozen['files']['models/response_adapter_v1/profiles.npz']
            or spec['profile_receipt_sha256'] != frozen['files']['models/response_adapter_v1/profiles.json']):
        raise ValueError('Arm/seed/profiles differ from frozen experiment')
    for source in receipt['sources']:
        if 'h1' in source['path'].lower() or sha256(source['path']) != source['sha256']:
            raise ValueError('Forbidden or altered response source')
    for path, digest in spec['code_hashes'].items():
        if frozen['files'].get(path) != digest: raise ValueError('Candidate code differs: '+path)
    validation = validate(prediction, prediction.with_suffix('.receipt.json'), excluded_source_tokens=['H1'])
    null = read(ROOT/'reports/control_preserving_identity_gate.json')
    if not null['passed'] or null['emitter_sha256'] != sha256(ROOT/'scripts/control_preserving_emitter.py'):
        raise ValueError('Exact identity gate must pass before performance scoring')
    folder = ROOT/'evaluations/redesign_h1_400'/name
    if (folder/'terminal.json').exists() or (folder/'agg_results.csv').exists():
        raise FileExistsError('Existing scoring attempt; do not overwrite')
    folder.mkdir(parents=True, exist_ok=True)
    checked = verify_local_copy(prediction, folder)
    if checked['prediction_sha256'] != receipt['output_sha256']:
        raise ValueError('Generated/scored prediction hashes differ')
    checked.update(arm=name, validation=validation, real_de_sha256=baseline['real_de_sha256'],
                   freeze_sha256=sha256(ROOT/'reports/redesign_implementation_freeze.json'),
                   entrypoint_sha256=sha256(__file__), verified_at_utc=datetime.now(timezone.utc).isoformat())
    (folder/'preflight.json').write_text(json.dumps(checked, indent=2)+'\n')
    env = os.environ.copy()
    env.update(OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='4', MKL_NUM_THREADS='1', NUMBA_NUM_THREADS='4',
               NUMBA_CACHE_DIR=str(ROOT/'.cache/numba'), MPLCONFIGDIR=str(ROOT/'.cache/matplotlib'))
    command = [str(ROOT / '.venv/bin/cell-eval2'), 'run',
        '-ap', str(folder/'local_scoring_input.h5ad'),
        '-ar', str(ROOT/'evaluations/dispersion_h1_training_400/reference.h5ad'),
        '--de-real', str(ROOT/'evaluations/dispersion_h1_training_400/incumbent_cpu/de_real.parquet'),
        '--preset', 'vcc2026', '--pert-col', 'target_gene', '--write-degenes',
        '--set', 'device=cpu', '--set', 'de.backend=scanpy', '--set', 'control_source=real',
        '--set', 'num_threads=4', '-o', str(folder)]
    result = subprocess.run(command, cwd=ROOT, env=env)
    terminal = dict(exit_code=result.returncode, arm=name, completed_at_utc=datetime.now(timezone.utc).isoformat(),
                    entrypoint_sha256=sha256(__file__), preflight_sha256=sha256(folder/'preflight.json'))
    if result.returncode == 0:
        observed = yaml.safe_load((folder/'run_params.yaml').read_text()); observed.pop('outdir', None)
        cfg = replace(EvalConfig.from_preset('vcc2026'), pert_col='target_gene', device='cpu',
                      control_source='real', num_threads=4)
        cfg = replace(cfg, de=replace(cfg.de, backend='scanpy'))
        expected = asdict(cfg); expected.pop('outdir', None)
        meta = read(folder/'run_meta.json')
        baseline_folder = ROOT/'evaluations/dispersion_h1_training_400/incumbent_cpu'
        baseline_meta = read(baseline_folder/'run_meta.json')
        real_identity = real_de_identity(folder, baseline_folder)
        if (observed != expected or meta['resolved_device'] != 'cpu'
                or meta['cell_eval2_version'] != '0.16.0'
                or meta['resolved_de_backend'] not in ['scanpy', None]
                or meta.get('de_pred_fingerprint') is not None
                or meta.get('de_real_fingerprint') != real_identity['baseline_strict_fingerprint']
                or any(meta.get(k) != baseline_meta.get(k) for k in ['source_fingerprint', 'anchor_semantic_identity'])):
            terminal.update(exit_code=98, validation_error='Unexpected scorer configuration/backend/version')
        terminal['real_de_identity'] = real_identity
    (folder/'terminal.json').write_text(json.dumps(terminal, indent=2)+'\n')
    raise SystemExit(terminal['exit_code'])


if __name__ == '__main__': main()
