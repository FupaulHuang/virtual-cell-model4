#!/usr/bin/env python3
"""Verify the new isolated builder reproduces every incumbent count before cache reuse."""
from datetime import datetime, timezone
import json
from pathlib import Path
import anndata as ad
import numpy as np
from scipy import sparse
from build_prediction import sha256

ROOT = Path(__file__).resolve().parents[1]


def main():
    old = ROOT / 'outputs/dispersion_h1_400_incumbent/prediction.h5ad'
    new = ROOT / 'outputs/bulk_h1_400_incumbent/prediction.h5ad'
    validation = json.loads((new.parent / 'validation.json').read_text())
    if not validation['valid'] or not validation['hash_verified']:
        raise ValueError('Complete new-incumbent validation is required first')
    a, b = ad.read_h5ad(old, backed='r'), ad.read_h5ad(new, backed='r')
    try:
        if a.shape != b.shape or a.shape != (20000, 18080):
            raise ValueError('Incumbent shape differs')
        if not np.array_equal(a.var_names, b.var_names):
            raise ValueError('Incumbent gene axis differs')
        for name in ['context', 'target_gene']:
            if not np.array_equal(a.obs[name].astype(str), b.obs[name].astype(str)):
                raise ValueError(f'Incumbent {name} axis differs')
        blocks = 0
        for start in range(0, a.n_obs, 512):
            left = sparse.csr_matrix(a.X[start:start + 512])
            right = sparse.csr_matrix(b.X[start:start + 512])
            if (left != right).nnz:
                raise ValueError(f'New builder differs from incumbent at row block {start}')
            blocks += 1
    finally:
        a.file.close(); b.file.close()
    receipt = json.loads(new.with_suffix('.receipt.json').read_text())
    record = {'verified_at_utc': datetime.now(timezone.utc).isoformat(),
              'counts_bit_identical': True, 'axes_identical': True,
              'all_count_blocks_checked': blocks, 'shape': [20000, 18080],
              'old_incumbent': str(old), 'new_incumbent': str(new),
              'old_incumbent_sha256': sha256(old), 'new_incumbent_sha256': sha256(new),
              'new_receipt_sha256': sha256(new.with_suffix('.receipt.json')),
              'builder_sha256': receipt['model_code_sha256']['bulk_builder'],
              'reuse': 'Existing incumbent metrics on the identical frozen reference may be reused; no repeat scoring required.'}
    output = ROOT / 'reports/bulk_baseline_identity.json'
    if output.exists():
        raise FileExistsError(output)
    output.write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps(record, indent=2), flush=True)


if __name__ == '__main__':
    main()
