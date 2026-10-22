#!/usr/bin/env python3
"""Generate a neutral donor-based comparator on a prepared local context."""
import os
for key in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS']:
    os.environ.setdefault(key, '1')
import argparse
import json
from pathlib import Path
import time
import numpy as np
from build_prediction import Writer, control_template, axis_csv, atlas, sha256
from validate_prediction_receipt import validate


def build(directory, output, *, cells=200, pool=4, seed=20261020):
    directory, output = Path(directory).resolve(), Path(output).resolve()
    if output.exists():
        raise FileExistsError(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    genes_path, targets_path = directory / 'gene_names.csv', directory / 'pert_counts.csv'
    genes, targets = axis_csv(genes_path, 'gene_name'), axis_csv(targets_path, 'target_gene')
    control_path = directory / 'context_H1.h5ad'
    template, depths, mean, bulk = control_template(control_path, genes, cells, pool, seed)
    started, moments = time.monotonic(), []
    pending = output.with_suffix('.partial.h5ad')
    if pending.exists():
        raise FileExistsError(pending)
    writer = Writer(pending, genes)
    try:
        for index, target in enumerate(targets):
            counts = atlas.dual_moment_counts(template, mean, bulk, depths=depths,
                                              seed=atlas.seed_for('H1:' + target, seed))
            writer.append(counts, 'H1', target)
            moments.append({'context': 'H1', 'target': target,
                'nnz': int(np.count_nonzero(counts)), 'max_library': int(depths.max()),
                'mean_probability_l1': float(abs((counts / depths[:, None]).mean(0) - mean).sum()),
                'bulk_probability_l1': float(abs(counts.sum(0) / counts.sum() - bulk).sum())})
            if (index + 1) % 10 == 0:
                print(f'Neutral H1: {index + 1}/{len(targets)}, {time.monotonic() - started:.1f}s', flush=True)
        writer.close()
    except BaseException:
        writer.f.close()
        raise
    pending.rename(output)
    receipt = {'arguments': {'genes': str(genes_path), 'targets': str(targets_path),
                 'context': [f'H1={control_path}'], 'cells': cells, 'pool': pool,
                 'seed': seed, 'target_limit': None},
               'sources': [], 'targets': targets.tolist(), 'output_genes': len(genes),
               'moments': moments, 'cells': writer.n, 'stored_entries': writer.nnz,
               'bytes': output.stat().st_size, 'output_sha256': sha256(output),
               'elapsed_seconds': time.monotonic() - started,
               'note': 'Neutral control-moment comparator; no perturbation outcomes were used.'}
    receipt_path = output.with_suffix('.receipt.json')
    receipt_path.write_text(json.dumps(receipt, indent=2) + '\n')
    result = validate(output, receipt_path, excluded_source_tokens=['H1'])
    output.with_suffix('.validation.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory', type=Path, default=Path('evaluations/h1_validation'))
    parser.add_argument('--output', type=Path, default=Path('evaluations/h1_validation/neutral_prediction.h5ad'))
    parser.add_argument('--cells', type=int, default=200)
    parser.add_argument('--pool', type=int, default=4)
    parser.add_argument('--seed', type=int, default=20261020)
    args = parser.parse_args()
    build(args.directory, args.output, cells=args.cells, pool=args.pool, seed=args.seed)
