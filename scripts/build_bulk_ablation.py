#!/usr/bin/env python3
"""Isolated four-arm bulk experiment using the unchanged frozen count emitter.

The base configuration and all inputs must match the incumbent receipt. Only
raw-source bulk centering or the final bulk amplitude changes. No dispersion
transform is applied. The read-only emitter profiler records the actual fitted
moments; it does not implement, replace, or relax the numerical projection.
"""
import os
for _key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_key, '1')
import argparse
from concurrent.futures import ProcessPoolExecutor
import copy
import hashlib
import json
import multiprocessing as mp
from pathlib import Path
import sys
import time

import anndata as ad
import numpy as np
import pandas as pd

from build_prediction import ROOT, Writer, atlas, axis_csv, control_template, sha256

ARMS = {
    'incumbent': {'raw_bulk_common_subtract': 1.0, 'bulk_amplitude': 0.3},
    'common05': {'raw_bulk_common_subtract': 0.5, 'bulk_amplitude': 0.3},
    'bulk015': {'raw_bulk_common_subtract': 1.0, 'bulk_amplitude': 0.15},
    'bulk060': {'raw_bulk_common_subtract': 1.0, 'bulk_amplitude': 0.6},
}


def array_sha256(value):
    a = np.ascontiguousarray(value)
    h = hashlib.sha256()
    h.update(str(a.dtype).encode())
    h.update(str(a.shape).encode())
    h.update(a.tobytes())
    return h.hexdigest()


def prepare_effects(sources, weights, targets, genes, arm):
    """Center on each source's full target population, never the output cohort."""
    settings = ARMS[arm]
    return {
        space: atlas.fuse_source_centered(
            sources, weights, targets, genes, space=space,
            common_subtract=1.0 if space == 'log2fc' else settings['raw_bulk_common_subtract'])
        for space in ('log2fc', 'bulk_delta')
    }


def make_profiles(effects, mean, bulk, cd4, cd4_weight, arm):
    """Return requested (arithmetic mean probability, bulk probability).

    ``cd4`` must already have been centered with common_subtract=1 and the
    source target population. It is shared unchanged by all four arms.
    """
    desired = []
    for space, ctrl, amplitude in (
        ('log2fc', mean, 0.6),
        ('bulk_delta', bulk, ARMS[arm]['bulk_amplitude']),
    ):
        effect, coverage = effects[space]
        if cd4 is not None:
            effect, coverage = atlas.add_cd4_family(
                effect, coverage, *cd4, cd4_weight,
                space=space, control_probability=ctrl)
        desired.append(atlas.desired_mean(ctrl, effect, space=space, amplitude=amplitude, clip=3))
    return tuple(desired)


def emit_audited(template, depths, probability, bulk_probability, seed):
    """Call the unchanged emitter and copy diagnostic vectors at its return.

    Profile hooks inspect only locals of this exact function's return frame.
    The previous hook is restored even if emission or audit fails. The copied
    arrays are small gene vectors; no live emitter array is mutated.
    """
    captured = {}
    previous = sys.getprofile()

    def observe(frame, event, result):
        if event == 'return' and frame.f_code is atlas.dual_moment_counts.__code__:
            local = frame.f_locals
            if result is None:
                for key in ('projection_error', 'cpm_error', 'bulk_error', 'step'):
                    if key in local:
                        captured[key] = float(local[key])
                return
            expected = local['expected']
            captured.update(
                bulk_probability_projected=local['bulk'].copy(),
                mean_probability_fitted=local['x'].mean(axis=0),
                bulk_probability_fitted=expected.sum(axis=0) / expected.sum(),
                expected_gene_count_sums=expected.sum(axis=0),
                projection_l1=float(local['projection_error']),
                fitting_mean_l1=float(local['cpm_error']),
                fitting_bulk_l1=float(local['bulk_error']),
                fitting_iterations=int(local['step']) + 1,
            )

    try:
        sys.setprofile(observe)
        counts = atlas.dual_moment_counts(template, probability, bulk_probability,
                                         depths=depths, seed=seed)
    except Exception as error:
        error.bulk_diagnostics = captured.copy()
        raise
    finally:
        sys.setprofile(previous)
    if not captured:
        raise RuntimeError('Frozen emitter diagnostics were not captured')
    realized_mean = (counts / np.asarray(depths)[:, None]).mean(axis=0)
    sums = counts.sum(axis=0, dtype=np.int64)
    realized_bulk = sums / sums.sum()
    gene_stats = {
        'mean_probability_requested': np.asarray(probability).copy(),
        'mean_probability_fitted': captured.pop('mean_probability_fitted'),
        'mean_probability_realized': realized_mean,
        'bulk_probability_requested': np.asarray(bulk_probability).copy(),
        'bulk_probability_projected': captured.pop('bulk_probability_projected'),
        'bulk_probability_fitted': captured.pop('bulk_probability_fitted'),
        'bulk_probability_realized': realized_bulk,
        'expected_gene_count_sums': captured.pop('expected_gene_count_sums'),
        'gene_count_sums_realized': sums,
        'mean_cpm_rounding_change': 1e6 * realized_mean,
    }
    gene_stats['mean_cpm_rounding_change'] -= 1e6 * gene_stats['mean_probability_fitted']
    captured.update(
        mean_probability_l1=float(np.abs(realized_mean - probability).sum()),
        bulk_probability_l1=float(np.abs(realized_bulk - bulk_probability).sum()),
        bulk_probability_l1_to_projected=float(np.abs(realized_bulk - gene_stats['bulk_probability_projected']).sum()),
        mean_cpm_rounding_max_abs=float(np.abs(gene_stats['mean_cpm_rounding_change']).max()),
        requested_mean_sha256=array_sha256(probability),
        requested_bulk_sha256=array_sha256(bulk_probability),
        projected_bulk_sha256=array_sha256(gene_stats['bulk_probability_projected']),
        counts_array_sha256=array_sha256(counts),
        max_library=int(counts.sum(axis=1).max()),
        nnz=int(np.count_nonzero(counts)),
    )
    return counts, captured, gene_stats


_TEMPLATE = None
_DEPTHS = None


def initialize_worker(template, depths):
    global _TEMPLATE, _DEPTHS
    _TEMPLATE, _DEPTHS = template, depths


def emit_worker(item):
    target, context, probability, bulk, seed = item
    try:
        x, audit, stats = emit_audited(_TEMPLATE, _DEPTHS, probability, bulk, seed)
    except Exception as error:
        return None, dict(context=context, target=target, seed=int(seed),
                          exception_type=type(error).__name__, exception=str(error),
                          emitter_diagnostics=getattr(error, 'bulk_diagnostics', {})), None
    audit.update(context=context, target=target, seed=int(seed))
    return x, audit, stats


def compare_incumbent(counts, original, depths, stats, *, require_identity=False):
    """Recompute realized margins and CPM from final integer counts."""
    for value in (counts, original):
        if not np.isfinite(value).all() or (value < 0).any() or (value != np.floor(value)).any():
            raise ValueError('Comparison requires finite nonnegative integer counts')
    if counts.shape != original.shape:
        raise ValueError('Incumbent and emitted axes differ')
    if not np.array_equal(counts.sum(axis=1, dtype=np.int64), depths):
        raise ValueError('Emission changed the frozen row depths')
    if not np.array_equal(original.sum(axis=1, dtype=np.int64), depths):
        raise ValueError('Incumbent and selected control depths differ')
    exact = np.array_equal(counts, original)
    if require_identity and not exact:
        raise ValueError('Incumbent arm is not bit-identical to frozen incumbent counts')
    before = (original / depths[:, None]).mean(axis=0)
    after = (counts / depths[:, None]).mean(axis=0)
    original_sums = original.sum(axis=0, dtype=np.int64)
    stats.update(mean_probability_incumbent=before,
                 bulk_probability_incumbent=original_sums / original_sums.sum(),
                 gene_count_sums_incumbent=original_sums,
                 mean_cpm_change_to_incumbent=1e6 * (after - before))
    absolute = np.abs(stats['mean_cpm_change_to_incumbent'])
    bound = np.maximum(1e-8, 1e-4 * before * 1e6)
    # This is a reported comparison, not a claim that the dispersion experiment's
    # CPM bound is imposed on a different mean-profile intervention.
    return dict(rows_exact=True, counts_bitidentical_to_incumbent=bool(exact),
                gene_totals_exact_to_incumbent=bool(np.array_equal(counts.sum(axis=0, dtype=np.int64), original_sums)),
                mean_probability_l1_to_incumbent=float(np.abs(after - before).sum()),
                bulk_probability_l1_to_incumbent=float(np.abs(stats['bulk_probability_realized'] - stats['bulk_probability_incumbent']).sum()),
                max_mean_cpm_absolute_drift=float(absolute.max()),
                max_mean_cpm_drift_fraction_of_dispersion_bound=float((absolute / bound).max()),
                genes_exceeding_dispersion_cpm_bound=int(np.count_nonzero(absolute > bound)),
                dispersion_cpm_bound_is_enforced=False)


def resolve(path):
    return (ROOT / path).resolve()


def verify_inputs(config, incumbent, receipt):
    """Require original model settings, dependencies, axes, and source bytes."""
    old = receipt['arguments']
    for key in ('cells', 'pool', 'prior_counts', 'common_subtract', 'seed',
                'amplitude_cpm', 'amplitude_bulk', 'mode', 'cd4_weight'):
        if config.get(key) != old.get(key):
            raise ValueError(f'Base configuration differs from incumbent: {key}')
    expected = dict(cells=400, pool=4, common_subtract=1.0, amplitude_cpm=0.6,
                    amplitude_bulk=0.3, mode='dual')
    if any(config.get(k) != v for k, v in expected.items()):
        raise ValueError('Only the preregistered 400-cell dual baseline is supported')
    if config.get('promoters') or config.get('agreement') or config.get('de_source'):
        raise ValueError('Promoters, agreement, and additional DE-only sources are out of scope')
    raw_specs = [f'{resolve(s["path"])}={float(s["weight"])}' for s in config['source']]
    old_specs = [f'{Path(p).resolve()}={float(w)}' for p, w in (x.rsplit('=', 1) for x in old['source'])]
    if raw_specs != old_specs or resolve(config['cd4']) != Path(old['cd4']).resolve():
        raise ValueError('Source identities, order, or weights differ from incumbent')
    for token in config.get('exclude_source_family', []):
        if any(token.lower() in s['path'].lower() for s in config['source']):
            raise ValueError(f'Forbidden source family: {token}')
    for s in receipt['sources']:
        if sha256(s['path']) != s['sha256']:
            raise ValueError(f'Source hash changed: {s["path"]}')
    for key, path in [('builder', ROOT/'scripts/build_prediction.py'),
                      ('atlas_model', ROOT/'external/atlasshift/model.py'),
                      ('atlas_control_template', ROOT/'external/atlasshift/predict.py')]:
        if sha256(path) != receipt['model_code_sha256'][key]:
            raise ValueError(f'Frozen dependency changed: {key}')
    for key in ('genes', 'targets'):
        if resolve(config[key]) != Path(old[key]).resolve() or sha256(resolve(config[key])) != receipt['auxiliary_input_sha256'][key]:
            raise ValueError(f'Frozen axis changed: {key}')
    contexts = [f'{c}={resolve(p)}' for c, p in (s.split('=', 1) for s in config['context'])]
    if contexts != old['context']:
        raise ValueError('Context identity/order changed')
    for context, source in receipt['contexts'].items():
        if sha256(source['path']) != source['sha256']:
            raise ValueError(f'Frozen control changed: {context}')
    if incumbent.resolve() != Path(old['output']).resolve() or sha256(incumbent) != receipt['output_sha256']:
        raise ValueError('Incumbent identity/hash changed')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--arm', choices=ARMS, required=True)
    parser.add_argument('--incumbent', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--workers', type=int)
    parser.add_argument('--target-limit', type=int, help='Engineering subset only; never a complete experiment arm')
    args = parser.parse_args()
    started = time.time()
    config = json.loads(args.config.read_text())
    workers = args.workers if args.workers is not None else config['workers']
    if workers < 1 or (args.target_limit is not None and args.target_limit < 1):
        raise ValueError('Workers and target limit must be positive')
    pending = args.output.with_suffix('.partial.h5ad')
    artifacts = [args.output, pending, args.output.with_suffix('.receipt.json'), args.output.parent/'failure.json',
                 args.output.parent/'profile_audit.json', args.output.parent/'gene_profile_audit.parquet']
    if any(p.exists() for p in artifacts):
        raise FileExistsError('Choose a new output directory; an artifact already exists')
    receipt_path = args.incumbent.with_suffix('.receipt.json')
    original_receipt = json.loads(receipt_path.read_text())
    verify_inputs(config, args.incumbent, original_receipt)
    targets = axis_csv(resolve(config['targets']), 'target_gene')
    if args.target_limit is not None:
        targets = targets[:args.target_limit]
    genes = axis_csv(resolve(config['genes']), 'gene_name')
    sources = [atlas.load_xatlas(resolve(s['path']), prior_counts=config['prior_counts']) for s in config['source']]
    weights = [float(s['weight']) for s in config['source']]
    effects = prepare_effects(sources, weights, targets, genes, args.arm)
    cd4 = atlas.aligned_cd4(resolve(config['cd4']), targets, genes, common_subtract=1.0, center_scope='source')
    incumbent = ad.read_h5ad(args.incumbent, backed='r')
    if not np.array_equal(incumbent.var_names.astype(str), genes):
        raise ValueError('Incumbent gene axis differs')
    labels = incumbent.obs.target_gene.astype(str).to_numpy()
    contexts = incumbent.obs.context.astype(str).to_numpy()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    writer = Writer(pending, genes)
    audits, tables = [], []
    context_profiles = {}
    failure_detail = {}
    try:
        for ci, spec in enumerate(config['context']):
            context, path = spec.split('=', 1)
            template, depths, mean, bulk = control_template(resolve(path), genes, config['cells'], config['pool'], config['seed'] + ci)
            p, b = make_profiles(effects, mean, bulk, cd4, config['cd4_weight'], args.arm)
            context_profiles[context] = dict(requested_mean_sha256=array_sha256(p), requested_bulk_sha256=array_sha256(b),
                template_sha256=array_sha256(template), depths_sha256=array_sha256(depths),
                cd4_centered_lfc_sha256=array_sha256(cd4[0]), cd4_available_sha256=array_sha256(cd4[1]))
            items = [(str(t), context, p[i], b[i], atlas.seed_for(context+':'+str(t), config['seed'])) for i, t in enumerate(targets)]
            with ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context('spawn'),
                                     initializer=initialize_worker, initargs=(template, depths)) as pool:
                futures, launched = {}, 0
                for i, target in enumerate(targets):
                    while launched < min(len(items), i + 2*workers):
                        futures[launched] = pool.submit(emit_worker, items[launched])
                        launched += 1
                    x, audit, stats = futures.pop(i).result()
                    failure_detail = dict(context=context, target=str(target))
                    if x is None:
                        failure_detail.update(audit)
                        raise ValueError(f'Frozen emitter failed for {context}/{target}: {audit["exception"]}')
                    rows = np.flatnonzero((contexts == context) & (labels == target))
                    if len(rows) != config['cells']:
                        raise ValueError(f'Incorrect incumbent group size: {context}/{target}')
                    original = incumbent.X[rows].toarray()
                    audit.update(compare_incumbent(x, original, depths, stats, require_identity=args.arm == 'incumbent'))
                    writer.append(x, context, str(target))
                    audits.append(audit)
                    table = pd.DataFrame(stats)
                    table.insert(0, 'gene', genes)
                    table.insert(0, 'target', str(target))
                    table.insert(0, 'context', context)
                    tables.append(table)
                    print(json.dumps(dict(arm=args.arm, completed_groups=len(audits), target=str(target),
                        elapsed_seconds=time.time()-started, projection_l1=audit['projection_l1'],
                        counts_identical=audit['counts_bitidentical_to_incumbent'])), flush=True)
        if writer.n != len(config['context'])*len(targets)*config['cells']:
            raise ValueError('Incorrect output row count')
        writer.close()
    except BaseException as error:
        writer.f.close()
        (args.output.parent/'failure.json').write_text(json.dumps(dict(
            status='failed_unpromoted', arm=args.arm, **ARMS[args.arm], completed_groups=len(audits),
            partial_output=str(pending), exception_type=type(error).__name__, exception=str(error),
            failed_group=failure_detail, unchanged_projection_l1_limit=0.03,
            unchanged_fit_l1_limit=0.001), indent=2)+'\n')
        raise
    finally:
        incumbent.file.close()
    pending.rename(args.output)
    pd.concat(tables, ignore_index=True).to_parquet(args.output.parent/'gene_profile_audit.parquet', index=False)
    (args.output.parent/'profile_audit.json').write_text(json.dumps(audits, indent=2)+'\n')
    receipt = copy.deepcopy(original_receipt)
    receipt['arguments'].update(output=str(args.output.resolve()), workers=workers,
        target_limit=args.target_limit, amplitude_bulk=ARMS[args.arm]['bulk_amplitude'],
        **ARMS[args.arm])
    receipt['targets'] = targets.tolist()
    coverage = (effects['log2fc'][1] > 0).any(axis=1) | cd4[1].any(axis=1)
    receipt['same_target_coverage'] = int(coverage.sum())
    receipt['uncovered_targets'] = targets[~coverage].tolist()
    receipt['moments'] = audits
    receipt['bulk_ablation'] = dict(arm=args.arm, **ARMS[args.arm], arithmetic_common_subtract=1.0,
        arithmetic_amplitude=0.6, cd4_common_subtract=1.0, center_scope='full source target populations',
        cd4_fixed_before_final_bulk_amplitude=True, dispersion_transform=None,
        status='engineering_subset' if args.target_limit is not None else 'unselected_development_candidate',
        base_config=str(args.config.resolve()), base_config_sha256=sha256(args.config),
        incumbent=str(args.incumbent.resolve()), incumbent_sha256=original_receipt['output_sha256'],
        incumbent_receipt_sha256=sha256(receipt_path), builder_sha256=sha256(__file__),
        projection_capture='Read-only return-frame profiler of unchanged pinned dual_moment_counts',
        unchanged_projection_l1_limit=0.03, unchanged_fit_l1_limit=0.001,
        context_profiles=context_profiles,
        gene_profile_audit_sha256=sha256(args.output.parent/'gene_profile_audit.parquet'),
        profile_audit_sha256=sha256(args.output.parent/'profile_audit.json'))
    receipt['ablation'] = dict(arm=args.arm, raw_bulk_common_subtract=ARMS[args.arm]['raw_bulk_common_subtract'],
        amplitude_bulk=ARMS[args.arm]['bulk_amplitude'], base_config_sha256=sha256(args.config),
        incumbent_prediction_sha256=original_receipt['output_sha256'])
    receipt['model_code_sha256']['bulk_builder'] = sha256(__file__)
    receipt.update(cells=writer.n, stored_entries=writer.nnz, elapsed_seconds=time.time()-started,
        bytes=args.output.stat().st_size, output_sha256=sha256(args.output),
        note='New isolated bulk ablation; unchanged dependencies recorded in model_code_sha256, new builder hash in bulk_ablation. No metric pass or selection is implied.')
    args.output.with_suffix('.receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps(dict(complete=True, arm=args.arm, output=str(args.output), groups=len(audits))), flush=True)


if __name__ == '__main__':
    main()
