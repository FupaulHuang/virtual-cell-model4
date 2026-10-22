#!/usr/bin/env python3
"""Materialize frozen public H1 validation controls without touching H1 Test.

The 50-target Model3 reference has 200 real cells per target (MAX has 161)
and 2,000 real controls. Preserve that exact reference rather than introducing
another random evaluation draw. Predictions may contain 200 cells per target.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import anndata as ad
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_REFERENCE = ROOT / 'data/legacy_model3/evaluations/m3_e07_hybrid_validation/reference.h5ad'


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(8 * 1024**2), b''):
            h.update(chunk)
    return h.hexdigest()


def prepare(reference, output):
    reference, output = Path(reference).resolve(), Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    if (output / 'manifest.json').exists():
        raise FileExistsError('Frozen evaluation already prepared; choose a new output directory')
    original = ad.read_h5ad(reference, backed='r')
    try:
        labels = original.obs.target_gene.astype(str).to_numpy()
        genes = original.var_names.astype(str).to_numpy()
        targets = np.asarray(sorted(set(labels) - {'non-targeting'}))
        if len(targets) != 50 or len(genes) != 18080:
            raise ValueError('Expected the frozen H1 validation 50-target/18080-gene reference')
        rows = np.flatnonzero(labels == 'non-targeting')
        if len(rows) < 800:
            raise ValueError('At least 800 controls required for 200 cells with pooling=4')
        control = original[rows].to_memory()
        control.obs['context'] = 'H1'
        control_path = output / 'context_H1.h5ad'
        control.write_h5ad(control_path, compression='lzf')
        truth_counts = {t: int(np.sum(labels == t)) for t in targets}
    finally:
        original.file.close()
    # Copy into this workspace to keep downstream cache/sidecar writes isolated.
    destination_reference = output / 'reference.h5ad'
    shutil.copyfile(reference, destination_reference)
    pd.DataFrame({'gene_name': genes}).to_csv(output / 'gene_names.csv', index=False)
    pd.DataFrame({'target_gene': targets, 'n_cells': 200}).to_csv(output / 'pert_counts.csv', index=False)
    manifest = {
        'evaluation_context': 'H1', 'split': 'public H1 2025 Validation',
        'reference_source': str(reference), 'reference_sha256': digest(reference),
        'reference_copy': str(destination_reference),
        'control_path': str(control_path), 'control_sha256': digest(control_path),
        'n_control_cells': len(rows), 'genes': len(genes), 'targets': targets.tolist(),
        'reference_cells_per_target': truth_counts, 'prediction_cells_per_target': 200,
        'required_source_exclusion': ['all H1 response banks for leave-context-out evaluation'],
        'allowed_inference_information': 'H1 controls and perturbation identities only',
        'test_split_used': False,
        'historical_comparator': str(ROOT / 'data/legacy_model3/evaluations/m3_e07_hybrid_validation/prediction.h5ad'),
        'historical_comparator_limit': 'Trained on H1 Training, so it is a historical target-holdout comparator, not leave-context-out evidence.',
        'scoring': 'Recompute all candidates, mean baseline and five replicate anchors on one CPU scanpy backend. Do not compare new CPU metrics to saved GPU/gpudge scores.',
    }
    (output / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps({k: v for k, v in manifest.items() if k not in {'targets', 'reference_cells_per_target'}}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference', type=Path, default=DEFAULT_REFERENCE)
    parser.add_argument('--output', type=Path, default=ROOT / 'evaluations/h1_validation')
    args = parser.parse_args()
    prepare(args.reference, args.output)
