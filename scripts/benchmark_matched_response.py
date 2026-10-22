#!/usr/bin/env python3
"""Preregistered matched raw-source response benchmark; no H1 or official DE labels.

Run `freeze` before `run`. Both phases refuse existing output artifacts. The
frozen v1 module supplies features and ridge fitting; it is never modified.
"""
from __future__ import annotations
import os
for _key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_key, '1')
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

import numpy as np
import pandas as pd
import response_transfer_adapter as adapter

ROOT = Path(__file__).resolve().parents[1]
PREFIX = ROOT/'reports/matched_response'
MANIFEST = PREFIX.with_name(PREFIX.name+'_manifest.json')
PAIRS = PREFIX.with_name(PREFIX.name+'_pairs.npz')
MODEL = ROOT/'models/response_adapter_v1/model.json'
PROTOCOL = ROOT/'docs/MATCHED_RESPONSE_BENCHMARK_PROTOCOL_20261021.md'
AUDIT = ROOT/'reports/matched_response_source_audit.json'
ADAPTER_SHA = '9eea25c6193418c2687b5f5fd8b5b04a62818c0e53cc54bf490658697276e497'
METHODS = ('zero', 'matched_fusion', 'adapter_v1')
THRESHOLDS = {'log2fc': .25, 'bulk_delta': .10}
POLICY = {'eval_pairs_per_context': 100000, 'train_pairs_per_task': 30000,
          'eval_seed': 20261023, 'train_seed': 20261021, 'bootstrap_seed': 20261023,
          'bootstrap_draws': 2000, 'minimum_cells': 20, 'input_prior_counts': 100000.,
          'label_prior_counts': 0., 'ridge': .03, 'clip': 3., 'thresholds': THRESHOLDS,
          'nonregression_tolerance': 1e-8, 'cd4': False, 'h1_responses': False,
          'label_kind': 'exploratory magnitude proxy, NOT differential-expression significance'}
PRIMITIVES = ('covered', 'predicted', 'reference', 'joint', 'correct', 'stable', 'stable_fp',
              'squared_error', 'stable_squared_error', 'responsive_squared_error')


def read(path):
    return json.loads(Path(path).read_text())


def banks_from_frozen_model():
    if adapter.sha256(ROOT/'scripts/response_transfer_adapter.py') != ADAPTER_SHA:
        raise ValueError('Immutable adapter implementation changed')
    model = read(MODEL)
    if model['implementation_sha256'] != ADAPTER_SHA:
        raise ValueError('Adapter model/module identity mismatch')
    banks = [adapter.load_raw_bank(p['path']) for p in model['source_banks']]
    if set(adapter.validate_families(banks, minimum=3)) != set(adapter.FAMILIES):
        raise ValueError('Exactly the three registered non-H1 raw families are required')
    for b, p in zip(banks, model['source_banks']):
        if b.digest != p['sha256']:
            raise ValueError('Frozen raw source changed')
    return sorted(banks, key=lambda b: adapter.FAMILIES.index(b.family))


def evaluation_pairs(destination, maximum, seed, minimum_cells=20):
    """Metadata-only sampling includes uncovered targets/genes; own gene excluded."""
    targets = destination.targets[destination.n_cells >= minimum_cells]
    if maximum < len(targets) or not len(targets):
        raise ValueError('Budget must represent every eligible target')
    rng = np.random.default_rng(seed); labels, genes, weights = [], [], []
    quotient, remainder = divmod(maximum, len(targets))
    for i, target in enumerate(targets):
        columns = np.flatnonzero(destination.measured & (destination.genes != target))
        n = min(len(columns), quotient+(i < remainder))
        if not n:
            raise ValueError('Eligible target has no measured non-own genes')
        chosen = np.sort(rng.choice(columns, size=n, replace=False))
        labels.extend([target]*n); genes.extend(destination.genes[chosen]); weights.extend([1/(len(targets)*n)]*n)
    return np.asarray(labels), np.asarray(genes), np.asarray(weights)


def put_pairs(arrays, key, pair):
    for suffix, value in zip(('targets', 'genes', 'weight'), pair):
        arrays[key+'_'+suffix] = value


def get_pairs(arrays, key):
    return tuple(arrays[key+'_'+s] for s in ('targets', 'genes', 'weight'))


def freeze():
    if MANIFEST.exists() or PAIRS.exists():
        raise FileExistsError('Benchmark pair/manifest freeze already exists')
    audit = read(AUDIT)
    if audit['cross_family_matched_DE_label_set_available']:
        raise ValueError('This preregistered implementation uses magnitude proxies only; review protocol before a different label estimand')
    banks = banks_from_frozen_model(); arrays = {}; folds = []
    for index, heldout in enumerate(banks):
        remaining = [b for b in banks if b.family != heldout.family]
        key = heldout.family+'_eval'
        pair = evaluation_pairs(heldout, POLICY['eval_pairs_per_context'], POLICY['eval_seed']+index)
        put_pairs(arrays, key, pair)
        tasks = []
        for destination in remaining:
            inputs = [b for b in remaining if b.family != destination.family]
            train_key = heldout.family+'_train_'+destination.family
            train_pair = adapter._sample_pairs(destination, inputs, POLICY['train_pairs_per_task'],
                                               POLICY['minimum_cells'], POLICY['train_seed'])
            put_pairs(arrays, train_key, train_pair)
            tasks.append({'destination': destination.family, 'inputs': [b.family for b in inputs],
                          'pairs_key': train_key, 'pairs': len(train_pair[0]), 'targets': len(set(train_pair[0]))})
        folds.append({'heldout': heldout.family, 'eval_key': key, 'eval_pairs': len(pair[0]),
                      'eval_targets': sorted(set(pair[0])), 'test_inputs': [b.family for b in remaining], 'training_tasks': tasks})
        arrays[heldout.family+'_native_genes'] = heldout.genes
        arrays[heldout.family+'_measured_genes'] = heldout.measured
    np.savez_compressed(PAIRS, **arrays)
    files = ['scripts/benchmark_matched_response.py', 'tests/test_matched_response_benchmark.py',
             'scripts/response_transfer_adapter.py', 'models/response_adapter_v1/model.json',
             str(PROTOCOL.relative_to(ROOT)), str(AUDIT.relative_to(ROOT))]
    manifest = {'frozen_at_utc': datetime.now(timezone.utc).isoformat(), 'policy': POLICY,
        'files': {p: adapter.sha256(ROOT/p) for p in files}, 'pairs_sha256': adapter.sha256(PAIRS),
        'sources': [b.provenance() for b in banks], 'folds': folds, 'outcomes_read_before_freeze': False,
        'history': 'Source contexts and v1 regularization were already examined previously; this is development, not untouched validation.',
        'label_audit': str(AUDIT.relative_to(ROOT)), 'matched_DE_labels_used': False,
        'fusion_definition': 'Specific feature constant coefficient1, shared/interactions0; equal-family centering, amplitude1, noCD4; not the production fusion recipe.'}
    MANIFEST.write_text(json.dumps(manifest, indent=2)+'\n')
    print(json.dumps({'manifest': str(MANIFEST), 'folds': [{k:f[k] for k in ('heldout','eval_pairs')} for f in folds]}), flush=True)


def task_from_pairs(destination, inputs, space, pair):
    if destination.family in {b.family for b in inputs}:
        raise ValueError('Destination response family cannot be its own input')
    targets, genes, weight = pair
    ri = np.asarray([destination._rows[str(t)] for t in targets]); ci = np.asarray([destination._cols[str(g)] for g in genes])
    if (np.any(targets == genes) or np.any(destination.n_cells[ri] < POLICY['minimum_cells'])
            or not destination.measured[ci].all() or not np.isclose(weight.sum(), 1)):
        raise ValueError('Frozen pair axis/eligibility/weight violation')
    x, confidence, count, source_mean = adapter.build_features(inputs, targets, genes,
        destination.control_mean_cpm[ri, ci], space, POLICY['input_prior_counts'], POLICY['minimum_cells'])
    y = destination.effects(space, 0, POLICY['minimum_cells'])[0][ri, ci].astype(float)
    return {'x': x, 'confidence_x': confidence, 'y': y, 'weight': weight, 'source_mean': source_mean,
            'confidence_y': ((source_mean*y)>0).astype(float), 'targets': targets, 'genes': genes,
            'coverage': count, 'destination': destination.family, 'input_families': [b.family for b in inputs]}


def predictions(task, head):
    # Shared basis/order is fixed by the immutable adapter implementation.
    if adapter.FEATURE_NAMES[0] != 'specific:constant':
        raise ValueError('Frozen specific/shared representation changed')
    result = {'zero': np.zeros(len(task['x'])), 'matched_fusion': task['x'][:, 0].copy(),
              'adapter_v1': adapter.apply_head(head, task['x'])}
    for name in result:
        result[name] = np.clip(result[name], -POLICY['clip'], POLICY['clip'])
        result[name][task['coverage'] == 0] = 0
    return result


def target_primitives(targets, truth, prediction, coverage, threshold):
    truth, prediction = np.asarray(truth), np.asarray(prediction)
    if not np.isfinite(truth).all() or not np.isfinite(prediction).all():
        raise ValueError('Finite measured effects and predictions required')
    target_axis, inverse = np.unique(targets, return_inverse=True)
    counts = np.bincount(inverse)
    reference = np.abs(truth) >= threshold; predicted = np.abs(prediction) >= threshold
    joint = reference & predicted; correct = joint & (np.sign(truth) == np.sign(prediction))
    squared = (prediction-truth)**2
    values = [coverage > 0, predicted, reference, joint, correct, ~reference, predicted & ~reference,
              squared, squared*~reference, squared*reference]
    matrix = np.column_stack([np.bincount(inverse, weights=np.asarray(v, dtype=float))/counts for v in values])
    return target_axis, matrix


def summarize_primitives(matrix):
    p = np.asarray(matrix).mean(axis=0)
    d = dict(zip(PRIMITIVES, p))
    def ratio(a,b): return float(d[a]/d[b]) if d[b] > 0 else None
    return {'source_coverage': float(d['covered']), 'predicted_call_fraction': float(d['predicted']),
        'reference_call_fraction': float(d['reference']), 'signed_precision': ratio('correct','predicted'),
        'signed_recovery': ratio('correct','reference'), 'joint_call_sign_accuracy': ratio('correct','joint'),
        'stable_false_call_fraction': ratio('stable_fp','stable'), 'unweighted_mse': float(d['squared_error']),
        'stable_mse': ratio('stable_squared_error','stable'), 'responsive_mse': ratio('responsive_squared_error','reference')}


def paired_intervals(base, candidate, draws=2000, seed=20261023):
    if base.shape != candidate.shape or len(base) < 1:
        raise ValueError('Paired target primitives required')
    rng = np.random.default_rng(seed); samples = {}
    for _ in range(draws):
        rows = rng.integers(0,len(base),len(base))
        a, b = summarize_primitives(base[rows]), summarize_primitives(candidate[rows])
        for name in a:
            value = b[name]-a[name] if a[name] is not None and b[name] is not None else np.nan
            samples.setdefault(name, []).append(value)
    result = {}
    for name, sample in samples.items():
        values = np.asarray(sample); finite = np.isfinite(values)
        result[name] = {'candidate_minus_fusion_95_ci': np.quantile(values,[.025,.975]).tolist() if finite.all() else None,
            'valid_draws': int(finite.sum()), 'undefined_draws': int((~finite).sum()), 'draws': draws,
            'note': 'Paired target-cluster bootstrap; no undefined draws silently removed.'}
    return result


def advancement(contexts):
    """Fusion-based prospective point gate; undefined precision cannot pass."""
    decisions = {}; under = False
    for space in adapter.SPACES:
        pairs = [(c[space]['matched_fusion'], c[space]['adapter_v1']) for c in contexts]
        complete = all(x['signed_precision'] is not None and x['signed_recovery'] is not None for p in pairs for x in p)
        under_count = sum(b['signed_recovery'] is not None and a['signed_recovery'] is not None
                          and b['signed_recovery'] < a['signed_recovery']-1e-8 for a,b in pairs)
        under |= under_count >= 2
        if complete:
            no_context_loss = all(b['signed_precision'] >= a['signed_precision']-1e-8 and
                                  b['signed_recovery'] >= a['signed_recovery']-1e-8 for a,b in pairs)
            gain = np.mean([b['signed_recovery']-a['signed_recovery'] for a,b in pairs])
            precision = np.mean([b['signed_precision']-a['signed_precision'] for a,b in pairs])
            mse = np.mean([b['unweighted_mse']-a['unweighted_mse'] for a,b in pairs])
            passed = no_context_loss and gain > 1e-8 and precision >= -1e-8 and mse <= 1e-8
        else:
            no_context_loss = False; gain = precision = mse = None; passed = False
        decisions[space] = {'passed': bool(passed), 'defined_precision_and_recovery': bool(complete),
            'no_context_precision_or_recovery_regression': bool(no_context_loss),
            'equal_context_recovery_gain': None if gain is None else float(gain),
            'equal_context_precision_gain': None if precision is None else float(precision),
            'equal_context_unweighted_mse_change': None if mse is None else float(mse),
            'contexts_with_under_recovery': int(under_count)}
    return {'heads': decisions, 'both_heads_pass': all(x['passed'] for x in decisions.values()),
            'under_recovery_established': bool(under), 'conditional_p2_design_allowed': bool(under),
            'h1_scoring_authorized_by_this_result': False}


def full_axis_diagnostic(destination, inputs, heads, targets):
    """Identical exact inverse on all measured destination genes, never sampled columns."""
    genes = destination.genes[destination.measured]; ci = np.flatnonzero(destination.measured)
    rows = []
    for start in range(0,len(targets),4):
        chunk = np.asarray(targets[start:start+4]); ri = np.asarray([destination._rows[t] for t in chunk])
        control_mean = destination.control_mean_cpm[np.ix_(ri,ci)]
        control_bulk = destination.control_probability[np.ix_(ri,ci)]
        if not np.allclose(control_bulk.sum(1),1,atol=1e-6):
            raise ValueError('Unmeasured coordinates unexpectedly carry bulk mass')
        native_mean = destination.mean_cpm[np.ix_(ri,ci)]
        native_mean = native_mean/native_mean.sum(1,keepdims=True)
        native_bulk = destination.counts[np.ix_(ri,ci)].astype(float)
        native_bulk /= native_bulk.sum(1,keepdims=True)
        effect = {method:{} for method in METHODS}
        for space in adapter.SPACES:
            t = np.repeat(chunk,len(genes)); g = np.tile(genes,len(chunk))
            x, _, coverage, _ = adapter.build_features(inputs,t,g,control_mean.ravel(),space,
                POLICY['input_prior_counts'],POLICY['minimum_cells'])
            pred = predictions({'x':x,'coverage':coverage}, heads[space])
            for method in METHODS:
                effect[method][space] = (pred[method].reshape(len(chunk),len(genes)), coverage.reshape(len(chunk),len(genes)))
        for method in METHODS:
            mean, bulk = adapter.effects_to_profiles(effect[method],control_mean,control_bulk,clip=POLICY['clip'])
            for i,target in enumerate(chunk):
                mask = genes != target
                rows.append({'context':destination.family,'target':str(target),'method':method,
                    'arithmetic_composition_log2_mse':float(np.mean((np.log2(1+1e6*mean[i,mask])-np.log2(1+1e6*native_mean[i,mask]))**2)),
                    'bulk_log1p_mse':float(np.mean((np.log1p(50000*bulk[i,mask])-np.log1p(50000*native_bulk[i,mask]))**2)),
                    'scored_genes':int(mask.sum())})
    return rows


def render(result):
    lines=['# Matched source-only response benchmark','',
        '**Labels are exploratory magnitude calls, not DE significance.** No H1 responses or official scoring were used. '
        'These source contexts were previously examined during v1 selection; this is a development benchmark.','',
        '| Context | Head | Predictor | Signed precision | Signed recovery | Call fraction | Unweighted MSE |',
        '|---|---|---|---:|---:|---:|---:|']
    def number(x): return 'undefined' if x is None else f'{x:.6f}'
    for context, heads in result['contexts'].items():
        for space, methods in heads.items():
            for method, values in methods.items():
                lines.append('| '+' | '.join([context,space,method,*[number(values[k]) for k in
                    ('signed_precision','signed_recovery','predicted_call_fraction','unweighted_mse')]])+' |')
    gate=result['advancement']
    lines += ['',f'Both-head advancement versus matched fusion: **{gate["both_heads_pass"]}**. '
              f'Under-recovery established: **{gate["under_recovery_established"]}**.', '',
              'Matched fusion fixes the common representation to specific coefficient1/shared0 with equal raw-family pooling, '
              'amplitude1 and no CD4. It is not the published production fusion. The held-out context is absent from every '
              'training input and label; the all-context fitted adapter is never used as a fold predictor.', '',
              'The JSON retains per-context stable/response errors, coverage, paired target-cluster intervals, undefined '
              'precision denominators and the fusion-based gate. Full measured-axis inverse diagnostics are separate '
              'normalized mean-composition/bulk errors, not reconstruction of per-cell CPM or official scores.', '',
              '[Source-label audit](matched_response_source_audit.md) rejected mismatched PIE q-values as raw-population labels. '
              '[Frozen protocol](../docs/MATCHED_RESPONSE_BENCHMARK_PROTOCOL_20261021.md) fixes thresholds, folds and criteria.']
    if gate['under_recovery_established']:
        lines += ['', 'One separately preregistered response-aware P2 objective may now be designed using clearly qualified '
                  'magnitude-proxy responsivity labels. No threshold sweep, new H1 scoring cycle, promotion or submission follows automatically.']
    return '\n'.join(lines)+'\n'


def run():
    result_path=PREFIX.with_name(PREFIX.name+'_results.json')
    if result_path.exists(): raise FileExistsError('Matched benchmark results already exist')
    manifest=read(MANIFEST)
    if manifest['policy'] != POLICY or adapter.sha256(PAIRS) != manifest['pairs_sha256']:
        raise ValueError('Frozen policy/pairs changed')
    for path,digest in manifest['files'].items():
        if adapter.sha256(ROOT/path) != digest: raise ValueError('Frozen benchmark artifact changed: '+path)
    banks=banks_from_frozen_model(); lookup={b.family:b for b in banks}
    with np.load(PAIRS,allow_pickle=False) as d: pairs={k:d[k] for k in d.files}
    contexts={}; intervals={}; fold_models={}; saved={}; flat=[]; inverse=[]
    for fold in manifest['folds']:
        heldout=lookup[fold['heldout']]; inputs=[lookup[f] for f in fold['test_inputs']]
        if heldout.family in {name for t in fold['training_tasks'] for name in [t['destination'],*t['inputs']]}:
            raise ValueError('Outer held-out source leaked into training')
        contexts[heldout.family]={}; intervals[heldout.family]={}; heads={}
        for space in adapter.SPACES:
            tasks=[task_from_pairs(lookup[t['destination']],[lookup[f] for f in t['inputs']],space,get_pairs(pairs,t['pairs_key']))
                   for t in fold['training_tasks']]
            if any(not (task['coverage']>0).all() for task in tasks): raise ValueError('Training pairs lack source support')
            head=adapter.fit_ridge(tasks,POLICY['ridge']); heads[space]=head
            task=task_from_pairs(heldout,inputs,space,get_pairs(pairs,fold['eval_key']))
            predictions_by_method=predictions(task,head); primitives={}; methods={}
            for method, value in predictions_by_method.items():
                axis, primitive=target_primitives(task['targets'],task['y'],value,task['coverage'],THRESHOLDS[space])
                primitives[method]=primitive; methods[method]=summarize_primitives(primitive)
                flat.append({'context':heldout.family,'head':space,'method':method,**methods[method]})
                saved[f'{heldout.family}_{space}_{method}']=value
            saved[f'{heldout.family}_{space}_observed']=task['y']
            saved[f'{heldout.family}_{space}_coverage']=task['coverage']
            saved[f'{heldout.family}_{space}_primitive_targets']=axis
            for method in METHODS: saved[f'{heldout.family}_{space}_{method}_target_primitives']=primitives[method]
            contexts[heldout.family][space]=methods
            intervals[heldout.family][space]=paired_intervals(primitives['matched_fusion'],primitives['adapter_v1'],
                POLICY['bootstrap_draws'],POLICY['bootstrap_seed'])
        fold_models[heldout.family]={'heads':heads,'training_tasks':fold['training_tasks']}
        inverse.extend(full_axis_diagnostic(heldout,inputs,heads,fold['eval_targets']))
        print(json.dumps({'completed_context':heldout.family}),flush=True)
    prediction_path=PREFIX.with_name(PREFIX.name+'_predictions.npz'); np.savez_compressed(prediction_path,**saved)
    model_path=PREFIX.with_name(PREFIX.name+'_fold_models.json'); model_path.write_text(json.dumps(fold_models,indent=2)+'\n')
    inverse_path=PREFIX.with_name(PREFIX.name+'_inverse_diagnostics.csv'); pd.DataFrame(inverse).to_csv(inverse_path,index=False)
    result={'manifest_sha256':adapter.sha256(MANIFEST),'analysis_sha256':adapter.sha256(__file__),'policy':POLICY,
        'contexts':contexts,'adapter_minus_fusion_intervals':intervals,'advancement':advancement(list(contexts.values())),
        'predictions_sha256':adapter.sha256(prediction_path),'fold_models_sha256':adapter.sha256(model_path),
        'inverse_diagnostics_sha256':adapter.sha256(inverse_path),'h1_responses_used':False,'official_DE_labels_used':False,
        'source_label_audit':str(AUDIT.relative_to(ROOT)),
        'inverse_diagnostic_mean_by_context_predictor':pd.DataFrame(inverse).groupby(['context','method'])[
            ['arithmetic_composition_log2_mse','bulk_log1p_mse']].mean().reset_index().to_dict(orient='records'),
        'limitations':['Previously examined source contexts; not untouched validation.',
            'Magnitude proxies do not establish significance, stable biology, or calibrated DE probability.',
            'Native panels, CPM denominators and intervention efficacy differ across source contexts.',
            'Outer-fold training has one input context per task while testing has two.',
            'Target-cluster intervals condition on measured noisy aggregate effects and shared controls.',
            'Common representation is matched; adapter coefficients may learn shared effects while fixed fusion subtracts the common component.',
            'No H1 scoring, reserve access, promotion or submission is authorized by this benchmark alone.']}
    result_path.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    result_path.with_suffix('.md').write_text(render(result)); pd.DataFrame(flat).to_csv(PREFIX.with_name(PREFIX.name+'_metrics.csv'),index=False)
    print(json.dumps(result['advancement'],indent=2),flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('action',choices=['freeze','run']); args=p.parse_args()
    freeze() if args.action=='freeze' else run()


if __name__=='__main__': main()
