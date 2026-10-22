#!/usr/bin/env python3
"""Prove complete current/current reproduction before reusing incumbent scores."""
import json
from datetime import datetime, timezone
import anndata as ad
import numpy as np
from scipy import sparse
from build_prediction import ROOT, sha256
from validate_prediction_receipt import validate


def main():
    old_path = ROOT/'outputs/dispersion_h1_400_incumbent/prediction.h5ad'
    new_path = ROOT/'outputs/redesign_h1_400/current_current_s20261021/prediction.h5ad'
    report_path = ROOT/'reports/redesign_baseline_identity.json'
    if report_path.exists(): raise FileExistsError(report_path)
    validate(new_path, new_path.with_suffix('.receipt.json'), excluded_source_tokens=['H1'])
    old, new = ad.read_h5ad(old_path, backed='r'), ad.read_h5ad(new_path, backed='r')
    try:
        if (old.shape != new.shape or not old.var_names.equals(new.var_names)
                or not old.obs_names.equals(new.obs_names)
                or not np.array_equal(old.obs[['context','target_gene']].astype(str), new.obs[['context','target_gene']].astype(str))):
            raise ValueError('Incumbent axes changed')
        for start in range(0, old.n_obs, 512):
            if (sparse.csr_matrix(old.X[start:start+512]) != sparse.csr_matrix(new.X[start:start+512])).nnz:
                raise ValueError('Incumbent counts changed')
    finally:
        old.file.close(); new.file.close()
    report = dict(counts_bit_identical=True, axes_identical=True,
        old_incumbent_sha256=sha256(old_path), new_incumbent_sha256=sha256(new_path),
        new_receipt_sha256=sha256(new_path.with_suffix('.receipt.json')),
        verifier_sha256=sha256(__file__), verified_at_utc=datetime.now(timezone.utc).isoformat())
    report_path.write_text(json.dumps(report, indent=2)+'\n')
    terminal_dir = ROOT/'evaluations/redesign_h1_400/current_current_s20261021'
    terminal_dir.mkdir(parents=True, exist_ok=True)
    terminal = terminal_dir/'terminal.json'
    if terminal.exists(): raise FileExistsError(terminal)
    terminal.write_text(json.dumps(dict(exit_code=0, arm='current_current_s20261021',
        reused=True, original_metrics=str(ROOT/'evaluations/dispersion_h1_training_400/incumbent_cpu'),
        identity_sha256=sha256(report_path), note='Original successful scorer result reused only after complete count/axis identity.'), indent=2)+'\n')
    print(json.dumps(report), flush=True)


if __name__ == '__main__': main()
