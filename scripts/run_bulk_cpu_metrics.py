#!/usr/bin/env python3
"""Verify bulk-ablation identity and reuse the exact reference DE for CPU scoring."""
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(key, '1')
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess
import sys
import anndata as ad
import numpy as np
from scipy import sparse
from build_prediction import sha256
import run_verified_dispersion_metrics as original

ROOT = Path(__file__).resolve().parents[1]
REFERENCE = ROOT / 'evaluations/dispersion_h1_training_400/reference.h5ad'
BASELINE = ROOT / 'evaluations/dispersion_h1_training_400/incumbent_cpu'
EVAL = ROOT / 'evaluations/bulk_h1_training_400'
CACHE = ROOT / 'reports/bulk_real_de_cache.json'


def read(path):
    return json.loads(Path(path).read_text())


def check_frozen():
    freeze = read(ROOT / 'reports/bulk_implementation_freeze.json')
    for path, digest in freeze['files'].items():
        if sha256(ROOT / path) != digest:
            raise ValueError(f'Frozen file changed: {path}')
    return freeze


def old_baseline(freeze):
    old_cache = read(ROOT / 'reports/dispersion_real_de_cache.json')
    if sha256(Path(original.__file__)) != old_cache['preflight_sha256']:
        raise ValueError('Original verified-baseline helper changed')
    original.frozen_files()
    current = original.baseline_checks()
    if any(old_cache.get(k) != v for k, v in current.items()):
        raise ValueError('Original incumbent/reference/config/real-DE cache changed')
    identity_path = ROOT / 'reports/bulk_baseline_identity.json'
    identity = read(identity_path)
    if not identity['counts_bit_identical'] or not identity['axes_identical']:
        raise ValueError('New builder has not reproduced the incumbent exactly')
    if identity['old_incumbent_sha256'] != current['prediction_sha256']:
        raise ValueError('Identity audit used a different old incumbent')
    if identity['new_incumbent_sha256'] != sha256(identity['new_incumbent']):
        raise ValueError('New incumbent changed after full identity audit')
    new_receipt_hash = sha256(Path(identity['new_incumbent']).with_suffix('.receipt.json'))
    if new_receipt_hash != identity['new_receipt_sha256']:
        raise ValueError('New incumbent receipt changed after the identity audit')
    if identity['builder_sha256'] != freeze['files']['scripts/build_bulk_ablation.py']:
        raise ValueError('Identity audit used a different isolated builder')
    current.update(identity_audit_sha256=sha256(identity_path),
                   new_incumbent_sha256=identity['new_incumbent_sha256'],
                   new_incumbent_receipt_sha256=new_receipt_hash)
    return current


def verify_local_copy(prediction, folder):
    local = folder / 'local_scoring_input.h5ad'
    if not local.exists():
        subprocess.run([sys.executable, str(ROOT / 'scripts/prepare_metric_input.py'),
                        str(prediction), '--reference', str(REFERENCE), '--output', str(local)],
                       cwd=ROOT, check=True)
    receipt = read(local.with_suffix('.receipt.json'))
    prediction_hash, reference_hash = sha256(prediction), sha256(REFERENCE)
    if (Path(receipt['prediction']).resolve() != prediction.resolve()
            or receipt['prediction_sha256'] != prediction_hash
            or Path(receipt['reference']).resolve() != REFERENCE.resolve()
            or receipt['reference_sha256'] != reference_hash
            or receipt['shape'] != [24000, 18080] or receipt['controls_added'] != 4000):
        raise ValueError('Stale or mismatched local scoring-copy receipt')
    local_obj = ad.read_h5ad(local, backed='r')
    predicted = ad.read_h5ad(prediction, backed='r')
    real = ad.read_h5ad(REFERENCE, backed='r')
    try:
        labels = local_obj.obs.target_gene.astype(str).to_numpy()
        if (local_obj.shape != (24000, 18080) or predicted.shape != (20000, 18080)
                or not np.array_equal(local_obj.var_names, predicted.var_names)
                or not np.array_equal(local_obj.var_names, real.var_names)
                or not np.array_equal(labels[:20000], predicted.obs.target_gene.astype(str))
                or not np.all(labels[20000:] == 'non-targeting')):
            raise ValueError('Scoring-copy axes differ')
        for start in range(0, 20000, 1024):
            end = min(start + 1024, 20000)
            a = sparse.csr_matrix(local_obj.X[start:end])
            b = sparse.csr_matrix(predicted.X[start:end])
            if (a != b).nnz:
                raise ValueError('Scored treated counts differ from the validated candidate')
        rows = np.flatnonzero(real.obs.target_gene.astype(str).to_numpy() == 'non-targeting')
        if len(rows) != 4000:
            raise ValueError('Reference control count changed')
        for start in range(0, 4000, 1024):
            end = min(start + 1024, 4000)
            a = sparse.csr_matrix(local_obj.X[20000 + start:20000 + end])
            b = sparse.csr_matrix(real.X[rows[start:end]])
            if (a != b).nnz:
                raise ValueError('Scored control counts differ from the reference')
    finally:
        local_obj.file.close(); predicted.file.close(); real.file.close()
    return {'prediction_sha256': prediction_hash, 'reference_sha256': reference_hash,
            'local_scoring_input_sha256': sha256(local),
            'local_scoring_receipt_sha256': sha256(local.with_suffix('.receipt.json')),
            'scored_treated_and_control_counts_verified': True}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('arm', choices=['freeze-real', 'common05', 'bulk015', 'bulk060'])
    args = p.parse_args()
    freeze = check_frozen()
    current = old_baseline(freeze)
    if args.arm == 'freeze-real':
        if CACHE.exists():
            raise FileExistsError(CACHE)
        current.update(frozen_at_utc=datetime.now(timezone.utc).isoformat(),
                       scoring_entrypoint_sha256=sha256(__file__),
                       implementation_freeze_sha256=sha256(ROOT / 'reports/bulk_implementation_freeze.json'))
        CACHE.write_text(json.dumps(current, indent=2) + '\n')
        print(json.dumps({'cache_frozen': str(CACHE), 'real_de_sha256': current['real_de_sha256']}), flush=True)
        return
    cache = read(CACHE)
    if (any(cache.get(k) != v for k, v in current.items())
            or cache['scoring_entrypoint_sha256'] != sha256(__file__)
            or cache['implementation_freeze_sha256'] != sha256(ROOT / 'reports/bulk_implementation_freeze.json')):
        raise ValueError('Bulk baseline/cache/code identity changed')
    manifest = read(ROOT / 'configs/bulk_h1_400_experiment.json')
    arm = args.arm; arm_spec = manifest['arms'][arm]
    prediction = ROOT / arm_spec['output']
    validation = read(prediction.parent / 'validation.json')
    if not validation['valid'] or not validation['hash_verified'] or validation['excluded_source_tokens'] != ['H1']:
        raise ValueError('Complete validated prediction with H1 source exclusion is required')
    receipt = read(prediction.with_suffix('.receipt.json'))
    expected = {'arm': arm, 'raw_bulk_common_subtract': arm_spec['raw_bulk_common_subtract'],
                'amplitude_bulk': arm_spec['amplitude_bulk'],
                'base_config_sha256': sha256(ROOT / manifest['base_config']),
                'incumbent_prediction_sha256': cache['prediction_sha256']}
    if any(receipt['ablation'].get(k) != v for k, v in expected.items()):
        raise ValueError('Candidate arm/parameters/provenance differ from registration')
    if receipt['model_code_sha256']['bulk_builder'] != freeze['files']['scripts/build_bulk_ablation.py']:
        raise ValueError('Candidate used a different builder')
    old_receipt = read(ROOT / 'outputs/dispersion_h1_400_incumbent/prediction.receipt.json')
    old_sources = {(x['path'], x['sha256'], x['weight']) for x in old_receipt['sources']}
    new_sources = {(x['path'], x['sha256'], x['weight']) for x in receipt['sources']}
    if old_sources != new_sources:
        raise ValueError('Candidate source banks or weights changed')
    folder = EVAL / f'{arm}_cpu'
    if (folder / 'agg_results.csv').exists():
        raise FileExistsError('Successful arm results already exist; do not rerun them')
    folder.mkdir(parents=True, exist_ok=True)
    identities = verify_local_copy(prediction, folder)
    if identities['prediction_sha256'] != receipt['output_sha256']:
        raise ValueError('Candidate bytes changed after generation/validation')
    identities.update(verified_at_utc=datetime.now(timezone.utc).isoformat(), arm=arm,
                      real_de_cache_sha256=sha256(CACHE), real_de_sha256=cache['real_de_sha256'],
                      scoring_entrypoint_sha256=sha256(__file__))
    (folder / 'preflight.json').write_text(json.dumps(identities, indent=2) + '\n')
    print(json.dumps({'preflight_passed': arm, **identities}), flush=True)
    env = os.environ.copy()
    env.update(OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='4', MKL_NUM_THREADS='1',
               NUMBA_NUM_THREADS='4', NUMBA_CACHE_DIR=str(ROOT / '.cache/numba'),
               MPLCONFIGDIR=str(ROOT / '.cache/matplotlib'))
    cmd = [str(ROOT / '.venv/bin/cell-eval2'), 'run',
           '-ap', str(folder / 'local_scoring_input.h5ad'), '-ar', str(REFERENCE),
           '--de-real', str(BASELINE / 'de_real.parquet'), '--write-degenes',
           '--preset', 'vcc2026', '--pert-col', 'target_gene', '--set', 'device=cpu',
           '--set', 'de.backend=scanpy', '--set', 'control_source=real',
           '--set', 'num_threads=4', '-o', str(folder)]
    subprocess.run(cmd, cwd=ROOT, env=env, check=True)


if __name__ == '__main__':
    main()
