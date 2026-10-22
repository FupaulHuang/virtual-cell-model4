#!/usr/bin/env python3
"""Cross-context atlas response transfer with calibrated raw-count emission.

Uses the separately pinned MIT AtlasShift numerical kernels. Sources and their
weights are explicit inputs. Evaluation must exclude the destination source.
"""
import os
for _key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_key, '1')
import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
import multiprocessing as mp
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'external/atlasshift'))
import model as atlas
from predict import control_template, axis_csv
import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(8 * 1024**2), b''):
            h.update(b)
    return h.hexdigest()


class Writer:
    def __init__(self, path, genes):
        self.f = h5py.File(path, 'w')
        self.f.attrs.update({'encoding-type': 'anndata', 'encoding-version': '0.1.0'})
        self.x = self.f.create_group('X')
        self.x.attrs.update({'encoding-type': 'csr_matrix', 'encoding-version': '0.1.0'})
        for k, dtype in [('data', 'int32'), ('indices', 'int32'), ('indptr', 'int64')]:
            self.x.create_dataset(k, shape=(1 if k == 'indptr' else 0,), maxshape=(None,),
                                  dtype=dtype, chunks=(65536,), compression='gzip', compression_opts=1,
                                  shuffle=True)
        self.x['indptr'][0] = 0
        self.n = 0
        self.nnz = 0
        self.obs = []
        self.genes = genes

    def append(self, x, context, target):
        if x.ndim != 2 or x.shape[1] != len(self.genes) or x.shape[0] == 0:
            raise ValueError('Count block has incorrect cell/gene axes')
        if not np.isfinite(x).all() or (x < 0).any() or (x != np.floor(x)).any():
            raise ValueError('Invalid raw counts')
        if np.max(x.sum(1)) > 1_000_000:
            raise ValueError('Cell library cap exceeded')
        x = sparse.csr_matrix(x, dtype=np.int32)
        if self.nnz + x.nnz > 4_750_000_000:
            raise ValueError('Submission sparse density cap exceeded')
        for k, a in [('data', x.data), ('indices', x.indices),
                     ('indptr', x.indptr[1:].astype(np.int64) + self.nnz)]:
            d = self.x[k]; n = len(d); d.resize((n + len(a),)); d[n:] = a
        self.n += x.shape[0]; self.nnz += x.nnz
        self.obs.extend([(context, target)] * x.shape[0])

    def close(self):
        self.x.attrs['shape'] = (self.n, len(self.genes))
        obs = pd.DataFrame(self.obs, columns=['context', 'target_gene'],
                           index=[f'm4_{i}' for i in range(self.n)])
        for k in obs:
            obs[k] = pd.Categorical(obs[k])
        ad.io.write_elem(self.f, 'obs', obs)
        ad.io.write_elem(self.f, 'var', pd.DataFrame(index=self.genes))
        for k in ['uns', 'layers', 'obsm', 'varm', 'obsp', 'varp']:
            ad.io.write_elem(self.f, k, {})
        self.f.close()


_TEMPLATE = None
_DEPTHS = None


def initialize_worker(template, depths):
    global _TEMPLATE, _DEPTHS
    _TEMPLATE, _DEPTHS = template, depths


def emit(item):
    target, context, p, b, seed, mode = item
    if mode == 'single':
        p = b
    x = atlas.dual_moment_counts(_TEMPLATE, p, b, depths=_DEPTHS, seed=seed)
    achieved_p = (x / _DEPTHS[:, None]).mean(0)
    achieved_b = x.sum(0) / x.sum()
    return x, {'context': context, 'target': target,
               'mean_probability_l1': float(np.abs(achieved_p - p).sum()),
               'bulk_probability_l1': float(np.abs(achieved_b - b).sum()),
               'max_library': int(x.sum(1).max()), 'nnz': int(np.count_nonzero(x))}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source', action='append', default=[], help='Raw source PATH=WEIGHT (repeat)')
    ap.add_argument('--de-source', action='append', default=[], help='DE-only source PATH=WEIGHT (repeat)')
    ap.add_argument('--cd4', type=Path)
    ap.add_argument('--cd4-weight', type=float, default=.5)
    ap.add_argument('--promoters', type=Path)
    ap.add_argument('--context', action='append', required=True, help='LABEL=CONTROL_H5AD')
    ap.add_argument('--genes', type=Path, default=ROOT/'data/gene_names.csv')
    ap.add_argument('--targets', type=Path, default=ROOT/'data/pert_counts.csv')
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--amplitude-cpm', type=float, default=.6)
    ap.add_argument('--amplitude-bulk', type=float, default=.3)
    ap.add_argument('--common-subtract', type=float, default=1.)
    ap.add_argument('--mode', choices=['dual', 'single'], default='dual')
    ap.add_argument('--prior-counts', type=float, default=100000.)
    ap.add_argument('--pool', type=int, default=4)
    ap.add_argument('--cells', type=int, default=400)
    ap.add_argument('--workers', type=int, default=4)
    ap.add_argument('--seed', type=int, default=20261020)
    ap.add_argument('--target-limit', type=int)
    ap.add_argument('--agreement', action='store_true', help='Require consistent measured source directions where at least two sources exist')
    ap.add_argument('--agreement-min-effect', type=float, default=.02)
    ap.add_argument('--agreement-fraction', type=float, default=.75)
    args = ap.parse_args()
    if not args.source and not args.de_source:
        ap.error('At least one measured source is required')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    pending = args.output.with_suffix('.partial.h5ad')
    if args.output.exists() or pending.exists():
        raise FileExistsError('Output exists; choose a new run path')
    targets = axis_csv(args.targets, 'target_gene')
    if args.target_limit:
        targets = targets[:args.target_limit]
    genes = axis_csv(args.genes, 'gene_name')
    sources, weights, provenance = [], [], []
    for spec in args.source:
        path, weight = spec.rsplit('=', 1)
        sources.append(atlas.load_xatlas(Path(path), prior_counts=args.prior_counts))
        weights.append(float(weight))
        provenance.append({'path': str(Path(path).resolve()), 'sha256': sha256(path), 'weight': float(weight)})
    effects = {space: atlas.fuse_source_centered(sources, weights, targets, genes,
                space=space, common_subtract=args.common_subtract)
               for space in ['log2fc', 'bulk_delta']}
    votes, vote_masks = [], []
    if args.agreement:
        for source, weight in zip(sources, weights):
            if weight <= 0:
                continue
            e, mask = atlas.aligned_effect(source, source.targets, genes, space='log2fc',
                                           common_subtract=args.common_subtract)
            ri = pd.Index(source.targets).get_indexer(targets)
            rr = np.flatnonzero(ri >= 0)
            aligned = np.zeros((len(targets),len(genes)),np.float32)
            measured = np.zeros(aligned.shape,bool)
            aligned[rr] = e[ri[rr]]; measured[rr] = mask[ri[rr]]
            votes.append(aligned); vote_masks.append(measured)
    for spec in args.de_source:
        path, weight = spec.rsplit('=', 1); weight = float(weight)
        if weight < 0:
            raise ValueError('Negative source weight')
        with np.load(path, allow_pickle=False) as d:
            ri = pd.Index(d['targets'].astype(str)).get_indexer(targets)
            ci = pd.Index(d['genes'].astype(str)).get_indexer(genes)
            rr, cc = np.flatnonzero(ri >= 0), np.flatnonzero(ci >= 0)
            available = d['measured'][np.ix_(ri[rr], ci[cc])].astype(bool)
            for space in effects:
                e = d['effects_'+space][np.ix_(ri[rr], ci[cc])].astype(float)
                e -= args.common_subtract * d['source_common_'+space][ci[cc]][None,:]
                # Missing gene-target pairs never contribute a false zero vote.
                e = np.where(available, e, 0.)
                old, cov = effects[space]
                numerator = old.astype(float) * cov
                numerator[np.ix_(rr, cc)] += weight * e
                cov = cov.copy(); cov[np.ix_(rr, cc)] += weight * available
                effects[space] = (np.divide(numerator, cov, out=np.zeros_like(numerator), where=cov>0), cov)
                if args.agreement and space == 'log2fc' and weight > 0:
                    aligned = np.zeros((len(targets),len(genes)),np.float32)
                    measured = np.zeros(aligned.shape,bool)
                    aligned[np.ix_(rr,cc)] = e; measured[np.ix_(rr,cc)] = available
                    votes.append(aligned); vote_masks.append(measured)
        provenance.append({'path': str(Path(path).resolve()), 'sha256': sha256(path), 'weight': weight,
          'kind': 'DE-only arithmetic CPM; bulk response approximated from arithmetic CPM, not measured pseudobulk'})
    cd4 = atlas.aligned_cd4(args.cd4, targets, genes, common_subtract=args.common_subtract,
                            center_scope='source') if args.cd4 else None
    if args.cd4:
        provenance.append({'path': str(args.cd4.resolve()), 'sha256': sha256(args.cd4),
                           'weight': args.cd4_weight, 'kind': 'CD4 publisher DE summaries'})
    gate = None
    agreement_diagnostics = None
    if args.agreement:
        from response_confidence import confidence_gate
        fused = effects['log2fc'][0]
        if cd4 and args.cd4_weight > 0:
            votes.append(cd4[0]); vote_masks.append(cd4[1])
            fused, _ = atlas.add_cd4_family(fused, effects['log2fc'][1], *cd4, args.cd4_weight,
                                            space='log2fc', control_probability=np.zeros(len(genes)))
        if not votes:
            raise ValueError('Agreement gating requires a positive-weight measured source')
        gate, agreement_diagnostics = confidence_gate(np.stack(votes), np.stack(vote_masks), fused,
                min_abs_effect=args.agreement_min_effect, min_agreement=args.agreement_fraction)
    receipt = {'arguments': {k: str(v) if isinstance(v, Path) else v for k,v in vars(args).items()},
               'sources': provenance, 'contexts': {}, 'moments': [],
               'output_genes': len(genes), 'targets': targets.tolist(),
               'note': 'Local emission audit is not an official challenge score.'}
    receipt['model_code_sha256'] = {'builder': sha256(__file__),
        'atlas_model': sha256(ROOT/'external/atlasshift/model.py'),
        'atlas_control_template': sha256(ROOT/'external/atlasshift/predict.py')}
    receipt['auxiliary_input_sha256'] = {'genes': sha256(args.genes), 'targets': sha256(args.targets)}
    if args.promoters:
        receipt['auxiliary_input_sha256']['promoters'] = sha256(args.promoters)
    if args.agreement:
        receipt['model_code_sha256']['response_confidence'] = sha256(ROOT/'scripts/response_confidence.py')
        receipt['agreement'] = agreement_diagnostics
    coverage = (effects['log2fc'][1] > 0).any(1)
    if cd4:
        coverage |= cd4[1].any(1)
    receipt['same_target_coverage'] = int(coverage.sum())
    receipt['uncovered_targets'] = targets[~coverage].tolist()
    print('Coverage', receipt['same_target_coverage'], '/', len(targets), flush=True)
    t0 = time.time(); writer = Writer(pending, genes)
    try:
        for ci, spec in enumerate(args.context):
            context, path = spec.split('=', 1)
            template, depths, mean, bulk = control_template(path, genes, args.cells, args.pool, args.seed + ci)
            desired = []
            promoter_changes = []
            for space, ctrl, amp in [('log2fc', mean, args.amplitude_cpm), ('bulk_delta', bulk, args.amplitude_bulk)]:
                effect, cov = effects[space]
                if cd4:
                    effect, cov = atlas.add_cd4_family(effect, cov, *cd4, args.cd4_weight,
                                                      space=space, control_probability=ctrl)
                if gate is not None:
                    effect = np.where(gate, effect, 0.)
                p = atlas.desired_mean(ctrl, effect, space=space, amplitude=amp, clip=3)
                if args.promoters:
                    p, changes = atlas.apply_promoter_prior(p, ctrl, targets, genes, args.promoters, .15)
                    promoter_changes.extend(changes)
                desired.append(p)
            receipt['contexts'][context] = {'path': str(Path(path).resolve()), 'sha256': sha256(path), 'pool': args.pool,
                'depth_min': int(depths.min()), 'depth_max': int(depths.max()),
                'promoter_changes': promoter_changes}
            items = [(str(t), context, desired[0][i], desired[1][i],
                      atlas.seed_for(context+':'+str(t), args.seed), args.mode)
                     for i,t in enumerate(targets)]
            with ProcessPoolExecutor(max_workers=args.workers, mp_context=mp.get_context('spawn'),
                 initializer=initialize_worker, initargs=(template, depths)) as pool:
                # Limit outstanding dense blocks to two per worker.
                pending_futures = {}
                launched = 0
                for i in range(len(items)):
                    while launched < min(len(items), i + 2*args.workers):
                        pending_futures[launched] = pool.submit(emit, items[launched]); launched += 1
                    x, audit = pending_futures.pop(i).result()
                    writer.append(x, context, str(targets[i]))
                    receipt['moments'].append(audit)
                    if (i+1) % 10 == 0 or i == len(items)-1:
                        print(context, i+1, '/', len(items), 'seconds', round(time.time()-t0),
                              'nnz', writer.nnz, flush=True)
        expected_cells = len(args.context) * len(targets) * args.cells
        if writer.n != expected_cells or len({s.split('=',1)[0] for s in args.context}) != len(args.context):
            raise ValueError('Incorrect total cells or duplicate context labels')
        writer.close()
    except BaseException:
        writer.f.close()
        raise
    pending.rename(args.output)
    receipt.update(cells=writer.n, stored_entries=writer.nnz, elapsed_seconds=time.time()-t0,
                   bytes=args.output.stat().st_size, output_sha256=sha256(args.output))
    args.output.with_suffix('.receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print('Saved', args.output, writer.n, 'cells', flush=True)


if __name__ == '__main__':
    main()
