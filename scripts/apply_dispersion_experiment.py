#!/usr/bin/env python3
"""Apply an audited dispersion transform to a frozen incumbent, group by group."""
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(key, '1')
import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
import multiprocessing as mp
from pathlib import Path
import time
import anndata as ad
import numpy as np
import pandas as pd
from build_prediction import Writer, sha256

ROOT = Path(__file__).resolve().parents[1]


def apply_group(item):
    from constrained_dispersion import redistribute_dispersion
    context, target, original, multiplier, seed = item
    result = redistribute_dispersion(original, multiplier, seed=seed,
                                     cpm_rtol=1e-4, cpm_atol=1e-8,
                                     max_passes=8, pair_window=4)
    output = result.counts
    if not np.isfinite(output).all() or (output < 0).any() or (output != np.floor(output)).any():
        raise ValueError(f'Invalid integer counts for {target}')
    depths = original.sum(axis=1, dtype=np.int64)
    if not np.array_equal(output.sum(axis=1, dtype=np.int64), depths):
        raise ValueError(f'Cell depths changed for {target}')
    if not np.array_equal(original.sum(axis=0, dtype=np.int64), output.sum(axis=0, dtype=np.int64)):
        raise ValueError(f'Gene totals changed for {target}')
    before = (original.astype(np.float64) / depths[:, None]).mean(axis=0) * 1e6
    after = (output.astype(np.float64) / depths[:, None]).mean(axis=0) * 1e6
    error = np.abs(after - before)
    bound = np.maximum(1e-8, 1e-4 * before)
    if np.any(error > bound * (1 + 1e-8) + 1e-10):
        raise ValueError(f'Realized CPM preservation failed for {target}')
    audit = {'context': context, 'target': target, 'nnz': int(np.count_nonzero(output)),
             'max_library': int(depths.max()), 'rows_exact': True, 'gene_totals_exact': True,
             'max_mean_cpm_absolute_drift': float(error.max()),
             'max_mean_cpm_drift_fraction_of_bound': float((error / bound).max()),
             'mean_probability_l1_to_incumbent': float(error.sum() / 1e6),
             'bulk_probability_l1_to_incumbent': 0.0,
             'transform': result.diagnostics}
    return output, audit, result.gene_stats


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--incumbent', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--multiplier', type=float, required=True)
    p.add_argument('--seed', type=int, default=20261021)
    p.add_argument('--workers', type=int, default=4)
    args = p.parse_args()
    if args.multiplier not in (1.5, 2.0):
        raise ValueError('Only preregistered experimental multipliers are allowed')
    pending = args.output.with_suffix('.partial.h5ad')
    if args.output.exists() or pending.exists():
        raise FileExistsError('Choose a new output directory')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    baseline_receipt = args.incumbent.with_suffix('.receipt.json')
    receipt = json.loads(baseline_receipt.read_text())
    incumbent_hash = sha256(args.incumbent)
    if incumbent_hash != receipt['output_sha256']:
        raise ValueError('Incumbent hash does not match its frozen receipt')
    a = ad.read_h5ad(args.incumbent, backed='r')
    genes = a.var_names.astype(str).to_numpy()
    contexts = a.obs.context.astype(str).to_numpy()
    labels = a.obs.target_gene.astype(str).to_numpy()
    groups = sorted(set(zip(contexts, labels)))
    writer = Writer(pending, genes)
    start = time.time()
    audits = []
    gene_tables = []
    try:
        with ProcessPoolExecutor(max_workers=args.workers, mp_context=mp.get_context('spawn')) as pool:
            futures = {}
            launched = 0
            for i, (context, target) in enumerate(groups):
                while launched < min(len(groups), i + args.workers):
                    c, t = groups[launched]
                    rows = np.flatnonzero((contexts == c) & (labels == t))
                    if len(rows) != 400:
                        raise ValueError('Every experimental group must contain 400 cells')
                    x = a.X[rows].toarray().astype(np.int32)
                    seed = int.from_bytes(hashlib.sha256(f'{args.seed}:{c}:{t}'.encode()).digest()[:4], 'little')
                    futures[launched] = pool.submit(apply_group, (c, t, x, args.multiplier, seed))
                    launched += 1
                x, audit, stats = futures.pop(i).result()
                writer.append(x, context, target)
                audits.append(audit)
                frame = pd.DataFrame(stats)
                if len(frame) != len(genes):
                    raise ValueError('Per-gene audit has an incorrect axis')
                frame.insert(0, 'gene', genes)
                frame.insert(0, 'target', target)
                frame.insert(0, 'context', context)
                gene_tables.append(frame)
                print(json.dumps({'completed_groups': i + 1, 'total_groups': len(groups),
                                  'target': target, 'elapsed_seconds': time.time() - start,
                                  'diagnostics': audit['transform']}, default=str), flush=True)
        writer.close()
    except BaseException:
        writer.f.close()
        raise
    finally:
        a.file.close()
    pending.rename(args.output)
    pd.concat(gene_tables, ignore_index=True).to_parquet(args.output.parent / 'gene_profile_audit.parquet', index=False)
    receipt['arguments']['output'] = str(args.output.resolve())
    receipt['moments'] = audits
    receipt['derivation'] = {'incumbent': str(args.incumbent.resolve()), 'incumbent_sha256': incumbent_hash,
                             'incumbent_receipt_sha256': sha256(baseline_receipt),
                             'transform': 'constrained_dispersion.redistribute_dispersion',
                             'requested_residual_sd_multiplier': args.multiplier, 'seed': args.seed,
                             'cpm_rtol': 1e-4, 'cpm_atol': 1e-8, 'max_passes': 8, 'pair_window': 4,
                             'module_sha256': sha256(ROOT / 'scripts/constrained_dispersion.py'),
                             'wrapper_sha256': sha256(__file__),
                             'note': 'Exact realized integer cell/gene totals; bounded arithmetic-CPM drift. Achieved dispersion may fall short of requested multiplier.'}
    receipt.update(cells=writer.n, stored_entries=writer.nnz, elapsed_seconds=time.time() - start,
                   bytes=args.output.stat().st_size, output_sha256=sha256(args.output))
    args.output.with_suffix('.receipt.json').write_text(json.dumps(receipt, indent=2, default=str) + '\n')
    (args.output.parent / 'profile_audit.json').write_text(json.dumps(audits, indent=2, default=str) + '\n')
    print(json.dumps({'complete': True, 'output': str(args.output), 'groups': len(groups),
                      'elapsed_seconds': time.time() - start}), flush=True)


if __name__ == '__main__':
    main()
