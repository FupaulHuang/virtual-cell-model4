#!/usr/bin/env python3
"""Freeze metadata-selected H1 Training development rows for a 400-cell experiment."""
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(key, '1')
import hashlib
import json
from pathlib import Path
import time
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from build_prediction import sha256

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'evaluations/dispersion_h1_training_400'
SEED = 20261021


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    freeze = OUT / 'cohort_freeze.json'
    if freeze.exists():
        raise FileExistsError('Development cohort already frozen; do not resample')
    source = ROOT / 'data/h1_2025/adata_Training.h5ad'
    a = ad.read_h5ad(source, backed='r')
    start = time.time()
    try:
        labels = a.obs.target_gene.astype(str).to_numpy()
        counts = pd.Series(labels).value_counts()
        genes = a.var_names.astype(str).to_numpy()
        if not a.obs_names.is_unique or not a.var_names.is_unique or len(genes) != 18080:
            raise ValueError('Unexpected or duplicate source axes')
        excluded = set()
        exclusion_files = []
        for name in ['h1_validation', 'h1_test', 'k562_confirmation']:
            path = ROOT / 'evaluations' / name / 'pert_counts.csv'
            excluded.update(pd.read_csv(path).target_gene.astype(str))
            exclusion_files.append({'path': str(path), 'sha256': sha256(path)})
        eligible = np.array(sorted(t for t, n in counts.items()
                                   if t != 'non-targeting' and n >= 400 and t not in excluded))
        if len(eligible) < 75:
            raise ValueError(f'Need 75 eligible disjoint targets; found {len(eligible)}')
        rng = np.random.default_rng(SEED)
        order = rng.permutation(eligible)
        development = sorted(order[:50].tolist())
        reserved = sorted(order[50:75].tolist())
        treated_rows = []
        for target in development:
            treated_rows.extend(rng.choice(np.flatnonzero(labels == target), 400,
                                           replace=False).tolist())
        controls = np.sort(rng.choice(np.flatnonzero(labels == 'non-targeting'),
                                     4000, replace=False))
        rows = np.sort(np.r_[treated_rows, controls].astype(np.int64))
        if set(reserved) & set(labels[rows]):
            raise ValueError('Reserved responses entered extraction')
        metadata = {
            'purpose': 'Development only: fixed 400-cell constrained dispersion comparison',
            'frozen_before_response_extraction_and_metrics': True,
            'source': str(source), 'source_bytes': source.stat().st_size,
            'source_mtime_ns': source.stat().st_mtime_ns,
            'seed': SEED, 'eligible_targets': eligible.tolist(),
            'development_targets': development, 'reserved_targets_not_extracted': reserved,
            'excluded_previous_target_files': exclusion_files,
            'selection': 'Uniform seeded selection among Training targets with >=400 genuine cells, excluding previous Validation and both confirmation target lists',
            'reference_cells_per_target': 400, 'prediction_cells_per_target': 400,
            'n_controls': 4000, 'n_genes': len(genes),
            'source_response_family_excluded': 'All H1 response banks',
            'limitation': 'New Model4 development cohort, not an untouched biological context; H1 Training was used in earlier source preparation. Control count differs from official contexts.',
            'metrics': 'Six paired raw official-definition metrics; no cross-panel anchors or inferred leaderboard score',
            'source_obs_sha256': hashlib.sha256('\n'.join(a.obs_names).encode()).hexdigest(),
            'source_gene_sha256': hashlib.sha256('\n'.join(genes).encode()).hexdigest(),
        }
        pd.DataFrame({'source_row': rows, 'source_cell': a.obs_names[rows],
                      'target_gene': labels[rows]}).to_csv(OUT / 'source_rows.csv', index=False)
        metadata['source_rows_sha256'] = sha256(OUT / 'source_rows.csv')
        freeze.write_text(json.dumps(metadata, indent=2) + '\n')
        print(json.dumps({'development_targets': len(development), 'reserved_targets': len(reserved),
                          'eligible_targets': len(eligible), 'rows_to_extract': len(rows)}), flush=True)
        blocks = []
        for offset in range(0, len(rows), 1024):
            x = sparse.csr_matrix(a.X[rows[offset:offset + 1024]])
            if not np.isfinite(x.data).all() or (x.data < 0).any() or (x.data != np.floor(x.data)).any():
                raise ValueError('Reference contains invalid raw counts')
            x.eliminate_zeros()
            blocks.append(x)
            if offset % 5120 == 0:
                print('Extracted', min(offset + 1024, len(rows)), '/', len(rows), flush=True)
        reference = ad.AnnData(sparse.vstack(blocks, format='csr'),
                              obs=a.obs.iloc[rows].copy(), var=a.var.copy())
    finally:
        a.file.close()
    reference.obs['context'] = 'H1'
    for column in reference.obs:
        if isinstance(reference.obs[column].dtype, pd.CategoricalDtype):
            reference.obs[column] = reference.obs[column].cat.remove_unused_categories()
    control = reference[reference.obs.target_gene.astype(str) == 'non-targeting'].copy()
    reference.write_h5ad(OUT / 'reference.h5ad', compression='gzip', compression_opts=1)
    control.write_h5ad(OUT / 'context_H1.h5ad', compression='lzf')
    pd.DataFrame({'gene_name': genes}).to_csv(OUT / 'gene_names.csv', index=False)
    pd.DataFrame({'target_gene': development, 'n_cells': 400}).to_csv(OUT / 'pert_counts.csv', index=False)
    metadata.update(reference_sha256=sha256(OUT / 'reference.h5ad'),
                    controls_sha256=sha256(OUT / 'context_H1.h5ad'),
                    shape=list(reference.shape), elapsed_seconds=time.time() - start,
                    cohort_freeze_sha256=sha256(freeze))
    (OUT / 'manifest.json').write_text(json.dumps(metadata, indent=2) + '\n')
    config = json.loads((ROOT / 'configs/h1_raw_dual_full.json').read_text())
    config.update(purpose='Frozen Model4 method applied to new 400-cell H1 Training development cohort; all H1 response sources excluded',
                  context=['H1=evaluations/dispersion_h1_training_400/context_H1.h5ad'],
                  genes='evaluations/dispersion_h1_training_400/gene_names.csv',
                  targets='evaluations/dispersion_h1_training_400/pert_counts.csv',
                  cells=400, workers=8, seed=SEED, name='dispersion_h1_400_incumbent',
                  output='outputs/dispersion_h1_400_incumbent/prediction.h5ad')
    path = ROOT / 'configs/dispersion_h1_400_incumbent.json'
    if path.exists():
        raise FileExistsError(path)
    path.write_text(json.dumps(config, indent=2) + '\n')
    print(json.dumps({'complete': True, 'shape': list(reference.shape),
                      'manifest': str(OUT / 'manifest.json'), 'elapsed_seconds': time.time() - start}), flush=True)


if __name__ == '__main__':
    main()
