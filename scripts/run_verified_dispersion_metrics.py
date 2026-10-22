#!/usr/bin/env python3
"""Verify prediction/scoring-copy identities and frozen real-DE provenance before scoring."""
import argparse
from dataclasses import asdict, replace
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess
import sys
import anndata as ad
import numpy as np
import yaml
from cell_eval2 import EvalConfig
from build_prediction import sha256

ROOT = Path(__file__).resolve().parents[1]
EVAL = ROOT / 'evaluations/dispersion_h1_training_400'
CACHE = ROOT / 'reports/dispersion_real_de_cache.json'


def read(path):
    return json.loads(Path(path).read_text())


def frozen_files():
    record = read(ROOT / 'reports/dispersion_implementation_freeze.json')
    for relative, digest in record['files'].items():
        if sha256(ROOT / relative) != digest:
            raise ValueError(f'Frozen file changed: {relative}')


def check_scoring_copy(arm, prediction, create=False):
    local = EVAL / f'{arm}_cpu/local_scoring_input.h5ad'
    if not local.exists():
        if not create:
            raise FileNotFoundError(local)
        subprocess.run([sys.executable, str(ROOT / 'scripts/prepare_metric_input.py'),
                        str(prediction), '--reference', str(EVAL / 'reference.h5ad'),
                        '--output', str(local)], cwd=ROOT, check=True)
    receipt = read(local.with_suffix('.receipt.json'))
    reference_hash = sha256(EVAL / 'reference.h5ad')
    prediction_hash = sha256(prediction)
    if (Path(receipt['prediction']).resolve() != prediction.resolve()
            or receipt['prediction_sha256'] != prediction_hash
            or Path(receipt['reference']).resolve() != (EVAL / 'reference.h5ad').resolve()
            or receipt['reference_sha256'] != reference_hash
            or receipt['controls_added'] != 4000 or receipt['shape'] != [24000, 18080]):
        raise ValueError('Existing local scoring copy has stale or mismatched provenance')
    a = ad.read_h5ad(local, backed='r')
    p = ad.read_h5ad(prediction, backed='r')
    try:
        labels = a.obs.target_gene.astype(str).to_numpy()
        if (a.shape != (24000, 18080) or p.shape != (20000, 18080)
                or not np.array_equal(a.var_names, p.var_names)
                or not np.array_equal(labels[:20000], p.obs.target_gene.astype(str))
                or not np.all(labels[20000:] == 'non-targeting')):
            raise ValueError('Local scoring copy axes or labels differ from prediction plus controls')
    finally:
        a.file.close(); p.file.close()
    return {'prediction_sha256': prediction_hash, 'reference_sha256': reference_hash,
            'local_scoring_input_sha256': sha256(local),
            'local_scoring_receipt_sha256': sha256(local.with_suffix('.receipt.json'))}


def baseline_checks():
    terminal = read(ROOT / 'reports/dispersion_terminal_statuses.json')
    if terminal.get('incumbent', {}).get('exit_code') != 0:
        raise ValueError('Successful incumbent terminal status is required before real-DE reuse')
    folder = EVAL / 'incumbent_cpu'
    meta = read(folder / 'run_meta.json')
    cfg = yaml.safe_load((folder / 'run_params.yaml').read_text())
    expected = replace(EvalConfig.from_preset('vcc2026'), pert_col='target_gene',
                       device='cpu', control_source='real', num_threads=4)
    expected = replace(expected, de=replace(expected.de, backend='scanpy'))
    expected_cfg = asdict(expected)
    # The output location differs by arm; every biological/scoring field must match.
    expected_cfg.pop('outdir', None)
    observed_cfg = dict(cfg); observed_cfg.pop('outdir', None)
    if observed_cfg != expected_cfg:
        raise ValueError('Baseline scoring configuration differs from the frozen intended preset')
    if (meta['resolved_device'] != 'cpu' or meta['resolved_de_backend'] != 'scanpy'
            or meta['cell_eval2_version'] != '0.16.0'
            or Path(meta['source']).resolve() != (EVAL / 'reference.h5ad').resolve()
            or meta.get('de_real_fingerprint') is not None):
        raise ValueError('Baseline real-DE was not freshly computed using the intended reference/backend')
    for name in ['agg_results.csv', 'metric_aggregation.csv', 'results.csv', 'de_real.parquet']:
        if not (folder / name).is_file():
            raise FileNotFoundError(folder / name)
    manifest = read(EVAL / 'manifest.json')
    if sha256(EVAL / 'reference.h5ad') != manifest['reference_sha256']:
        raise ValueError('Reference differs from current audited manifest')
    if sha256(EVAL / 'context_H1.h5ad') != manifest['controls_sha256']:
        raise ValueError('Control file differs from audited manifest')
    source = check_scoring_copy('incumbent', ROOT / 'outputs/dispersion_h1_400_incumbent/prediction.h5ad')
    return dict(source, controls_sha256=manifest['controls_sha256'],
                baseline_run_meta_sha256=sha256(folder / 'run_meta.json'),
                baseline_run_params_sha256=sha256(folder / 'run_params.yaml'),
                real_de_sha256=sha256(folder / 'de_real.parquet'),
                source_fingerprint=meta['source_fingerprint'])


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('action', choices=['freeze-real', 'sd15', 'sd20'])
    args = p.parse_args()
    frozen_files()
    current = baseline_checks()
    if args.action == 'freeze-real':
        if CACHE.exists():
            raise FileExistsError(CACHE)
        current.update(frozen_at_utc=datetime.now(timezone.utc).isoformat(),
                       preflight_sha256=sha256(__file__),
                       purpose='Real-DE exported by the successful incumbent run on the exact same reference, controls, preset and backend')
        CACHE.write_text(json.dumps(current, indent=2) + '\n')
        print(json.dumps({'cache_frozen': str(CACHE), 'real_de_sha256': current['real_de_sha256']}), flush=True)
        return
    cache = read(CACHE)
    if any(cache.get(k) != v for k, v in current.items()) or cache['preflight_sha256'] != sha256(__file__):
        raise ValueError('Baseline/reference/config/real-DE/preflight identity changed after freeze')
    arm = args.action
    prediction = ROOT / f'outputs/dispersion_h1_400_{arm}/prediction.h5ad'
    validation = read(prediction.parent / 'validation.json')
    if not validation['valid'] or not validation['hash_verified'] or validation['excluded_source_tokens'] != ['H1']:
        raise ValueError('Candidate has no successful complete validation and H1-exclusion audit')
    candidate_receipt = read(prediction.with_suffix('.receipt.json'))
    derivation = candidate_receipt['derivation']
    if derivation['incumbent_sha256'] != cache['prediction_sha256']:
        raise ValueError('Candidate derives from a different incumbent')
    implementation = read(ROOT / 'reports/dispersion_implementation_freeze.json')
    expected_derivation = {
        'module_sha256': implementation['files']['scripts/constrained_dispersion.py'],
        'wrapper_sha256': implementation['files']['scripts/apply_dispersion_experiment.py'],
        'requested_residual_sd_multiplier': implementation['arms'][arm],
        'seed': implementation['seed'], 'cpm_rtol': implementation['cpm_rtol'],
        'cpm_atol': implementation['cpm_atol'], 'max_passes': 8, 'pair_window': 4,
    }
    if any(derivation.get(k) != v for k, v in expected_derivation.items()):
        raise ValueError('Candidate derivation differs from its frozen arm/module/parameters')
    identities = check_scoring_copy(arm, prediction, create=True)
    if candidate_receipt['output_sha256'] != identities['prediction_sha256']:
        raise ValueError('Candidate bytes changed after validation')
    identities.update(verified_utc=datetime.now(timezone.utc).isoformat(),
                      real_de_cache_sha256=sha256(CACHE), preflight_sha256=sha256(__file__),
                      real_de_sha256=cache['real_de_sha256'])
    (EVAL / f'{arm}_cpu/preflight.json').write_text(json.dumps(identities, indent=2) + '\n')
    print(json.dumps({'preflight_passed': arm, **identities}), flush=True)
    subprocess.run(['bash', str(ROOT / 'scripts/run_dispersion_cpu_metrics.sh'), arm,
                    str(prediction)], cwd=ROOT, check=True)


if __name__ == '__main__':
    main()
