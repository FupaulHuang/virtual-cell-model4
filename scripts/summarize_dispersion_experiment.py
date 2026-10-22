#!/usr/bin/env python3
"""Preregistered paired analysis of the 400-cell dispersion development arms.

No model fitting, scoring, reference resampling, anchor scaling or submission.
Run only after root records successful metric termination for all three arms.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
EVAL = ROOT / 'evaluations/dispersion_h1_training_400'
ARMS = ('incumbent', 'sd15', 'sd20')
METRICS = {
    'pds_cosine': ('PDS', 1),
    'expr_mse_unbiased_capped_norm': ('Expression MSE', -1),
    'de_wilcoxon_lfc_nmae': ('LFC NMAE', -1),
    'de_wilcoxon_direction_fidelity_yield_raw': ('Direction fidelity', 1),
    'de_wilcoxon_direction_reach_raw': ('Direction reach', 1),
    'de_wilcoxon_sig_jaccard': ('Significance Jaccard', 1),
}
PRIMARY = ('de_wilcoxon_direction_fidelity_yield_raw', 'de_wilcoxon_sig_jaccard')
PROTECTED = tuple(k for k in METRICS if k not in PRIMARY)
BOOTSTRAPS, SEED = 2000, 20261021
PDS_ATOL, NONREGRESSION_ATOL = 1e-10, 1e-8
CPM_RTOL, CPM_ATOL = 1e-4, 1e-8


def read_json(path):
    return json.loads(Path(path).read_text())


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024**2), b''):
            digest.update(block)
    return digest.hexdigest()


def paired_bootstrap(difference):
    difference = np.asarray(difference, dtype=float)
    if difference.ndim != 1 or not len(difference) or not np.isfinite(difference).all():
        raise ValueError('Paired bootstrap requires a nonempty finite difference vector')
    rng = np.random.default_rng(SEED)
    draws = rng.integers(0, len(difference), size=(BOOTSTRAPS, len(difference)))
    boot = difference[draws].mean(axis=1)
    return {'n_targets': len(difference), 'n_resamples': BOOTSTRAPS, 'seed': SEED,
            'mean_gain': float(difference.mean()),
            'percentile_95_ci': [float(x) for x in np.quantile(boot, [.025, .975])],
            'improved_targets': int((difference > 1e-12).sum()),
            'worsened_targets': int((difference < -1e-12).sum()),
            'tied_targets': int((abs(difference) <= 1e-12).sum())}


def metric_arm(arm):
    directory = EVAL / f'{arm}_cpu'
    meta = read_json(directory / 'run_meta.json')
    config = yaml.safe_load((directory / 'run_params.yaml').read_text())
    # Pinned build_run_meta reports no backend when *real* DE was supplied,
    # although predicted DE is computed with the explicitly configured engine.
    # The supplied real table is independently verified below before any verdict.
    backend_record_ok = (meta['resolved_de_backend'] == 'scanpy' or
                         (meta['resolved_de_backend'] is None and meta.get('de_real_fingerprint') is not None))
    if (meta['resolved_device'] != 'cpu' or not backend_record_ok or config['de']['backend'] != 'scanpy'
            or config['control_source'] != 'real' or config['metrics'] != 'vcc2026'
            or config['input_type'] != 'counts' or meta.get('de_pred_fingerprint') is not None):
        raise ValueError(f'{arm}: unexpected scoring semantics or supplied predicted DE')
    reference = Path(meta['source'])
    if not reference.is_absolute():
        reference = ROOT / reference
    if reference.resolve() != (EVAL / 'reference.h5ad').resolve():
        raise ValueError(f'{arm}: wrong reference path')
    aggregate = pd.read_csv(directory / 'agg_results.csv').set_index('statistic')
    aggregation = pd.read_csv(directory / 'metric_aggregation.csv').set_index('metric')
    raw = {metric: float(aggregate.loc['mean', metric]) for metric in METRICS}
    if not np.isfinite(list(raw.values())).all():
        raise ValueError(f'{arm}: all six finite raw metrics are required')
    cohort = {metric: None if bool(aggregation.loc[metric, 'derived'])
              else int(aggregation.loc[metric, 'n_used']) for metric in METRICS}
    auxiliary = {name: float(aggregate.loc['mean', name]) for name in
                 ('expr_distance_unbiased', 'expr_mse_unbiased', 'expr_mse_unbiased_capped')}
    return {'raw': raw, 'finite_cohorts': cohort, 'meta': meta, 'config': config,
            'auxiliary_means': auxiliary,
            'tidy': pd.read_csv(directory / 'results.csv'), 'directory': directory}


def real_de_identity(arm_directory, baseline_directory):
    from cell_eval2.cache import fingerprint_de_table
    from cell_eval2.de import load_de_table
    baseline_path = baseline_directory / 'de_real.parquet'
    path = arm_directory / 'de_real.parquet'
    baseline = pd.read_parquet(baseline_path).sort_values(['target', 'feature']).reset_index(drop=True)
    current = pd.read_parquet(path).sort_values(['target', 'feature']).reset_index(drop=True)
    pd.testing.assert_frame_equal(baseline, current, check_dtype=False, check_exact=True)
    fingerprint = fingerprint_de_table(load_de_table(baseline_path), strict=True)
    return {'identical_values': True, 'baseline_table_sha256': sha256(baseline_path),
            'exported_table_sha256': sha256(path), 'baseline_strict_fingerprint': fingerprint}


def reference_call_diagnostics(directory, targets):
    tables = {}
    for side in ('real', 'pred'):
        table = pd.read_parquet(directory / f'de_{side}.parquet')
        tables[side] = table[(table.p_adj < .05) & np.isfinite(table.log2_fold_change)
                              & (table.target != table.feature)][['target', 'feature']].drop_duplicates()
    real, pred = tables['real'], tables['pred']
    common = real.merge(pred, on=['target', 'feature'], how='inner', validate='one_to_one')
    counts = {}
    for name, table in [('real', real), ('pred', pred), ('tp', common)]:
        counts[name] = table.groupby('target').size().reindex(targets, fill_value=0)
    per_target = []
    for target in targets:
        r, p, tp = (int(counts[k].loc[target]) for k in ('real', 'pred', 'tp'))
        per_target.append({'target': target, 'real_calls': r, 'pred_calls': p,
                           'tp': tp, 'fp': p - tp, 'fn': r - tp})
    return {'p_adj_rule': 'strict p_adj < 0.05; finite LFC; own target gene excluded',
            'real_calls': len(real), 'pred_calls': len(pred), 'tp': len(common),
            'fp': len(pred) - len(common), 'fn': len(real) - len(common),
            'per_target': per_target,
            'interpretation': 'TP/FP/FN refer to the finite experimental DE-call set, not proof of biological truth. Totals are diagnostics, not a replacement for per-target official Jaccard.'}


def matrix_audit(arm, targets):
    base_path = ROOT / 'outputs/dispersion_h1_400_incumbent/prediction.h5ad'
    directory = ROOT / f'outputs/dispersion_h1_400_{arm}'
    candidate_path = directory / 'prediction.h5ad'
    receipt = read_json(directory / 'prediction.receipt.json')
    validation = read_json(directory / 'validation.json')
    if (not validation.get('valid') or not validation.get('hash_verified')
            or validation.get('excluded_source_tokens') != ['H1']
            or validation.get('shape') != [400 * len(targets), 18080]
            or validation.get('nnz', 4_750_000_001) > 4_750_000_000):
        raise ValueError(f'{arm}: complete format/hash/source-exclusion validation is required')
    if sha256(candidate_path) != receipt['output_sha256']:
        raise ValueError(f'{arm}: output hash mismatch')
    if receipt['derivation']['incumbent_sha256'] != sha256(base_path):
        raise ValueError(f'{arm}: wrong incumbent derivation')
    declared = read_json(directory / 'profile_audit.json')
    expected_pairs = {('H1', t) for t in targets}
    if {(x['context'], x['target']) for x in declared} != expected_pairs or len(declared) != len(expected_pairs):
        raise ValueError(f'{arm}: incomplete or duplicate profile audits')
    profile = pd.read_parquet(directory / 'gene_profile_audit.parquet')
    if profile.duplicated(['context', 'target', 'gene']).any():
        raise ValueError(f'{arm}: duplicate gene profile rows')
    before, after = ad.read_h5ad(base_path, backed='r'), ad.read_h5ad(candidate_path, backed='r')
    scoring = ad.read_h5ad(EVAL / f'{arm}_cpu/local_scoring_input.h5ad', backed='r')
    baseline_scoring = ad.read_h5ad(EVAL / 'incumbent_cpu/local_scoring_input.h5ad', backed='r')
    groups = []
    try:
        if (before.shape != after.shape or not before.obs_names.equals(after.obs_names)
                or not before.var_names.equals(after.var_names)
                or not np.array_equal(before.obs[['context', 'target_gene']].astype(str),
                                      after.obs[['context', 'target_gene']].astype(str))):
            raise ValueError(f'{arm}: cell, label or gene axes changed')
        if before.shape != (400 * len(targets), 18080):
            raise ValueError(f'{arm}: wrong complete development shape')
        labels = before.obs.target_gene.astype(str).to_numpy()
        if set(labels) != set(targets):
            raise ValueError(f'{arm}: unexpected target axis')
        for target in targets:
            rows = np.flatnonzero(labels == target)
            if len(rows) != 400:
                raise ValueError(f'{arm}: incorrect cells for {target}')
            a, b = before.X[rows].toarray(), after.X[rows].toarray()
            if (not np.array_equal(scoring.X[rows].toarray(), b)
                    or not np.array_equal(baseline_scoring.X[rows].toarray(), a)):
                raise ValueError(f'{arm}: scored treated counts differ from prediction for {target}')
            numeric = all(np.isfinite(x).all() and (x >= 0).all() and (x == np.floor(x)).all() for x in (a, b))
            if not numeric:
                raise ValueError(f'{arm}: invalid counts for {target}')
            a, b = a.astype(np.int64), b.astype(np.int64)
            depth = a.sum(1, dtype=np.int64)
            new_depth = b.sum(1, dtype=np.int64)
            row_exact = np.array_equal(depth, new_depth)
            column_exact = np.array_equal(a.sum(0, dtype=np.int64), b.sum(0, dtype=np.int64))
            if (depth <= 0).any() or (depth > 1_000_000).any():
                raise ValueError(f'{arm}: invalid incumbent cell depths')
            if (new_depth <= 0).any() or (new_depth > 1_000_000).any():
                raise ValueError(f'{arm}: invalid candidate cell depths')
            ca = a.astype(float) * (1e6 / depth)[:, None]
            cb = b.astype(float) * (1e6 / new_depth)[:, None]
            ma, mb = ca.mean(0), cb.mean(0)
            drift = abs(mb - ma)
            tolerance = np.maximum(CPM_ATOL, CPM_RTOL * ma)
            violations = int((drift > tolerance).sum())
            va, vb = ca.var(0, ddof=1), cb.var(0, ddof=1)
            p = profile[(profile.context == 'H1') & (profile.target == target)].set_index('gene').reindex(before.var_names)
            if len(p) != len(before.var_names) or p.mean_cpm_before.isna().any():
                raise ValueError(f'{arm}: incomplete per-gene evidence for {target}')
            # Independent final-count recomputation must agree with exported diagnostic values.
            for key, actual in [('mean_cpm_before', ma), ('mean_cpm_after', mb),
                                ('variance_cpm_before', va), ('variance_cpm_after', vb)]:
                np.testing.assert_allclose(p[key], actual, rtol=1e-11, atol=1e-8)
            base_trace = float(va.sum())
            achieved = float(np.sqrt(vb.sum() / base_trace)) if base_trace > 0 else None
            groups.append({'target': target, 'rows_exact': bool(row_exact),
                'gene_totals_exact': bool(column_exact), 'cpm_violations': violations,
                'max_cpm_drift_fraction_of_bound': float((drift / tolerance).max()),
                'max_abs_cpm_drift': float(drift.max()),
                'max_abs_lfc_drift': float(abs(np.log2((mb + 1e-9) / (ma + 1e-9))).max()),
                'achieved_aggregate_sd_multiplier': achieved,
                'nnz': int(np.count_nonzero(b)), 'genes_variance_decreased': int((vb < va - 1e-12).sum())})
    finally:
        before.file.close()
        after.file.close()
        scoring.file.close()
        baseline_scoring.file.close()
    passed = all(g['rows_exact'] and g['gene_totals_exact'] and g['cpm_violations'] == 0 for g in groups)
    achieved = [g['achieved_aggregate_sd_multiplier'] for g in groups if g['achieved_aggregate_sd_multiplier'] is not None]
    return {'passed': passed, 'independent_final_count_recomputation': True,
            'format_validation': validation,
            'output_sha256': receipt['output_sha256'], 'derivation': receipt['derivation'],
            'profile_audit_sha256': sha256(directory / 'profile_audit.json'),
            'gene_profile_audit_sha256': sha256(directory / 'gene_profile_audit.parquet'),
            'group_count': len(groups), 'groups': groups,
            'max_cpm_bound_fraction': max(g['max_cpm_drift_fraction_of_bound'] for g in groups),
            'max_abs_lfc_drift': max(g['max_abs_lfc_drift'] for g in groups),
            'achieved_sd_min_median_max': [float(x) for x in np.quantile(achieved, [0, .5, 1])] if achieved else None,
            'declared_target_reached_groups': sum(x['transform']['target_reached'] for x in declared),
            'duplicate_increase_groups': sum(x['transform']['duplicate_cells_increased'] for x in declared)}


def main():
    terminal_path = ROOT / 'reports/dispersion_terminal_statuses.json'
    terminal = read_json(terminal_path)
    statuses = terminal.get('arms', terminal)
    for arm in ARMS:
        if statuses[arm].get('exit_code') != 0:
            raise ValueError(f'{arm}: metric process is not authoritatively successful')
    freeze = read_json(ROOT / 'reports/dispersion_implementation_freeze.json')
    for name, expected in freeze['files'].items():
        if sha256(ROOT / name) != expected:
            raise ValueError(f'Implementation or preregistration changed since freeze: {name}')
    cohort = read_json(EVAL / 'cohort_freeze.json')
    targets = cohort['development_targets']
    if len(targets) != 50 or set(targets) & set(cohort['reserved_targets_not_extracted']):
        raise ValueError('Invalid development/reserve partition')
    manifest = read_json(EVAL / 'manifest.json')
    if sha256(EVAL / 'reference.h5ad') != manifest['reference_sha256']:
        raise ValueError('Reference file differs from the current repaired manifest')
    data = {arm: metric_arm(arm) for arm in ARMS}
    base = data['incumbent']
    cache_path = ROOT / 'reports/dispersion_real_de_cache.json'
    cache = read_json(cache_path)
    if cache['real_de_sha256'] != sha256(base['directory'] / 'de_real.parquet'):
        raise ValueError('Real-DE cache hash differs from the completed baseline export')
    comparisons = {}
    for arm in ('sd15', 'sd20'):
        candidate = data[arm]
        for key in ('source_fingerprint', 'anchor_semantic_identity', 'resolved_device'):
            if candidate['meta'][key] != base['meta'][key]:
                raise ValueError(f'{arm}: scoring/reference identity mismatch in {key}')
        ignored = {'outdir', 'cache_real', 'cache_pred'}
        if ({k: v for k, v in candidate['config'].items() if k not in ignored}
                != {k: v for k, v in base['config'].items() if k not in ignored}):
            raise ValueError(f'{arm}: metric configuration mismatch')
        de_identity = real_de_identity(candidate['directory'], base['directory'])
        preflight = read_json(candidate['directory'] / 'preflight.json')
        if (preflight['real_de_cache_sha256'] != sha256(cache_path)
                or preflight['real_de_sha256'] != cache['real_de_sha256']
                or preflight['reference_sha256'] != manifest['reference_sha256']
                or preflight['local_scoring_input_sha256'] != sha256(candidate['directory'] / 'local_scoring_input.h5ad')):
            raise ValueError(f'{arm}: preflight/cache/scoring-copy identity mismatch')
        supplied = candidate['meta'].get('de_real_fingerprint')
        if supplied != de_identity['baseline_strict_fingerprint']:
            raise ValueError(f'{arm}: supplied real DE fingerprint does not match frozen baseline table')
        metrics = {}
        for metric, (label, direction) in METRICS.items():
            gain = direction * (candidate['raw'][metric] - base['raw'][metric])
            item = {'label': label, 'better_direction': 'higher' if direction == 1 else 'lower',
                    'incumbent': base['raw'][metric], 'candidate': candidate['raw'][metric],
                    'gain_in_better_direction': gain, 'paired_bootstrap': None,
                    'finite_cohort_identical': True}
            if metric != 'expr_mse_unbiased_capped_norm':
                paired = []
                masks = []
                for d in (base, candidate):
                    part = d['tidy'].loc[d['tidy'].metric == metric, ['perturbation', 'value']]
                    if part.perturbation.duplicated().any() or not set(part.perturbation).issubset(targets):
                        raise ValueError(f'{arm}: invalid per-target rows for {metric}')
                    values = part.set_index('perturbation').value.reindex(targets).to_numpy(dtype=float)
                    masks.append(np.isfinite(values))
                    paired.append(values)
                counts = [int(mask.sum()) for mask in masks]
                if base['finite_cohorts'][metric] != counts[0] or candidate['finite_cohorts'][metric] != counts[1]:
                    raise ValueError(f'{arm}: reported/observed finite cohorts disagree')
                identical = np.array_equal(masks[0], masks[1]) and bool(masks[0].any())
                item['finite_cohort_identical'] = identical
                item['finite_target_ids'] = {'incumbent': np.asarray(targets)[masks[0]].tolist(),
                                             'candidate': np.asarray(targets)[masks[1]].tolist()}
                if identical:
                    delta = direction * (paired[1][masks[0]] - paired[0][masks[0]])
                    item['paired_bootstrap'] = paired_bootstrap(delta)
                    if abs(float(delta.mean()) - gain) > 1e-10:
                        raise ValueError(f'{arm}: raw aggregate is not the expected paired mean for {metric}')
                else:
                    item['uncertainty_note'] = 'Finite target cohorts differ; no intersection-only bootstrap or comparable-cohort gain claim is made. This fails the strict comparison.'
                if metric == 'pds_cosine':
                    item['max_abs_per_target_change'] = float(abs(delta).max()) if identical else None
            else:
                item['uncertainty_note'] = 'Global ratio-of-sums with whole-panel correction budget; no invalid per-target-normalized bootstrap is performed.'
            if metric in PROTECTED:
                item['point_nonregression'] = bool(gain >= -NONREGRESSION_ATOL)
                if item['paired_bootstrap'] is not None:
                    item['bootstrap_noninferiority_supported'] = item['paired_bootstrap']['percentile_95_ci'][0] >= -NONREGRESSION_ATOL
            else:
                item['positive_gain'] = bool(gain > 0)
                item['positive_bootstrap_interval'] = bool(item['paired_bootstrap'] is not None
                    and item['paired_bootstrap']['percentile_95_ci'][0] > 0)
            metrics[metric] = item
        audit = matrix_audit(arm, targets)
        reference_distance = base['auxiliary_means']['expr_distance_unbiased']
        if abs(candidate['auxiliary_means']['expr_distance_unbiased'] - reference_distance) > 1e-12:
            raise ValueError(f'{arm}: reference-only MSE denominator changed')
        for d in (base, candidate):
            if reference_distance <= 0 or abs(d['auxiliary_means']['expr_mse_unbiased_capped'] / reference_distance
                    - d['raw']['expr_mse_unbiased_capped_norm']) > 1e-10:
                raise ValueError(f'{arm}: cannot attribute MSE without matching finite component cohorts')
        bulks_fixed = all(g['gene_totals_exact'] for g in audit['groups'])
        mse_correction = {
            'identical_integer_gene_totals_verified': bulks_fixed,
            'reference_distance_mean': reference_distance,
            'inferred_change_in_mean_effective_prediction_correction':
                metrics['expr_mse_unbiased_capped_norm']['gain_in_better_direction'] * reference_distance if bulks_fixed else None,
            'interpretation': 'With verified identical whole-panel pseudobulks and fixed reference, MSE movement comes from the capped sampling-noise correction (apart from numerical effects), not improved predicted expression means or better mean-response biology.'
                if bulks_fixed else 'Raw gene totals changed, so a sampling-correction-only attribution is not justified; the structural invariant failed.'}
        cohorts_fixed = all(x['finite_cohort_identical'] for x in metrics.values())
        pds_fixed = (metrics['pds_cosine']['max_abs_per_target_change'] is not None
                     and abs(metrics['pds_cosine']['gain_in_better_direction']) <= PDS_ATOL
                     and metrics['pds_cosine']['max_abs_per_target_change'] <= PDS_ATOL)
        protected = all(metrics[m]['point_nonregression'] for m in PROTECTED)
        primary = all(metrics[m]['positive_gain'] for m in PRIMARY)
        primary_ci = all(metrics[m]['positive_bootstrap_interval'] for m in PRIMARY)
        strict = audit['passed'] and cohorts_fixed and pds_fixed and protected and primary
        comparisons[arm] = {'strict_point_verdict': 'PASS STRICT' if strict else 'FAIL STRICT',
            'structural_pass': audit['passed'], 'pds_invariance_pass': pds_fixed,
            'all_finite_target_cohorts_identical': cohorts_fixed,
            'protected_metrics_nonworse': protected, 'both_primary_point_gains_positive': primary,
            'both_primary_bootstrap_intervals_positive': primary_ci,
            'development_evidence_status': 'promising development evidence only' if strict and primary_ci
                else 'point dominance only; confidence inconclusive' if strict else 'strict criterion not met',
            'promoted': False, 'metrics': metrics, 'matrix_audit': audit,
            'mse_correction_interpretation': mse_correction,
            'real_de_identity': de_identity, 'preflight': preflight,
            'backend_verification': {'declared_predicted_de_backend': candidate['config']['de']['backend'],
                'recorded_resolved_de_backend': candidate['meta']['resolved_de_backend'],
                'predicted_de_was_supplied': candidate['meta'].get('de_pred_fingerprint') is not None,
                'cached_real_table_verified': True,
                'metadata_note': 'Pinned cell-eval2 baseline.py:1065-1067 records null backend when realDE is supplied; explicit scanpy config and fresh predictedDE are checked.'},
            'de_call_diagnostics': reference_call_diagnostics(candidate['directory'], targets)}
    result = {'protocol': 'docs/DISPERSION_EXPERIMENT_20261021.md',
        'protocol_sha256': freeze['files']['docs/DISPERSION_EXPERIMENT_20261021.md'],
        'implementation_freeze': freeze, 'cohort_freeze': cohort,
        'current_reference_manifest': manifest, 'terminal_statuses': terminal,
        'reference_category_repair': read_json(EVAL / 'reference_category_fix.json'),
        'parameters': {'bootstrap_resamples': BOOTSTRAPS, 'bootstrap_seed': SEED,
                       'pds_invariance_atol': PDS_ATOL, 'protected_nonregression_atol': NONREGRESSION_ATOL,
                       'cpm_rtol': CPM_RTOL, 'cpm_atol': CPM_ATOL},
        'incumbent': {'raw_metrics': base['raw'], 'finite_cohorts': base['finite_cohorts'],
                      'de_call_diagnostics': reference_call_diagnostics(base['directory'], targets)},
        'comparisons': comparisons, 'overall_raw_score_computed': False, 'anchors_used': False,
        'limitations': ['Development evidence in one previously used H1 context; not an untouched biological validation.',
            'The 25-target reserve remains unextracted; closed H1 Test/K562 confirmations are not reused.',
            'Requested aggregate residual SD is not each gene\'s multiplier; inspect achieved distributions and variance decreases.',
            'Strict numerical point nonregression is not statistical proof of absence of harm.',
            'Bootstrap intervals are paired target percentiles for the fixed finite cohorts, without familywise multiplicity adjustment; separate confirmation is required.',
            'Four thousand measured controls differ from the official context-control count; no leaderboard score is inferred.']}
    output = ROOT / 'reports/dispersion_experiment_results.json'
    output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    lines = ['# Constrained dispersion development results', '',
        'All arms use the frozen 50-target, 400-cell H1 Training development cohort. '
        'Values below are raw official-definition metrics; no heterogeneous overall score or cross-panel anchors are used.', '',
        '| Metric | Incumbent | SD1.5 | SD2.0 | Better |', '| --- | ---: | ---: | ---: | --- |']
    for metric, (label, direction) in METRICS.items():
        lines.append(f'| {label} | {base["raw"][metric]:.9f} | {data["sd15"]["raw"][metric]:.9f} | '
                     f'{data["sd20"]["raw"][metric]:.9f} | {"higher" if direction == 1 else "lower"} |')
    for arm, comparison in comparisons.items():
        lines += ['', f'## {arm}: {comparison["strict_point_verdict"]}', '',
            f'Independent matrix invariants: {comparison["structural_pass"]}; PDS invariant: {comparison["pds_invariance_pass"]}; '
            f'finite cohorts identical: {comparison["all_finite_target_cohorts_identical"]}; '
            f'all protected point estimates nonworse: {comparison["protected_metrics_nonworse"]}. '
            f'Both primary gains positive: {comparison["both_primary_point_gains_positive"]}; '
            f'both paired 95% intervals above zero: {comparison["both_primary_bootstrap_intervals_positive"]}.', '',
            f'Evidence status: {comparison["development_evidence_status"]}. No model promotion.', '',
            '| Metric | Gain in better direction | Paired 95% interval |', '| --- | ---: | --- |']
        for metric, item in comparison['metrics'].items():
            interval = item['paired_bootstrap']['percentile_95_ci'] if item['paired_bootstrap'] else None
            text = f'[{interval[0]:+.9f}, {interval[1]:+.9f}]' if interval else item['uncertainty_note']
            lines.append(f'| {item["label"]} | {item["gain_in_better_direction"]:+.9f} | {text} |')
        audit = comparison['matrix_audit']
        lines += ['', f'Achieved aggregate SD multiplier, min/median/max: {audit["achieved_sd_min_median_max"]}. '
            f'Target reached in {audit["declared_target_reached_groups"]}/50 groups; '
            f'duplicate-cell increase flags: {audit["duplicate_increase_groups"]}. '
            f'Maximum realized LFC drift: {audit["max_abs_lfc_drift"]:.9g}.']
        lines += ['', comparison['mse_correction_interpretation']['interpretation']]
    lines += ['', 'The JSON includes source/config identity checks, actual finite target IDs, paired win/loss/tie counts, '
              'strict p_adj<0.05 reference-relative TP/FP/FN, full independent count audits and terminal statuses.', '']
    lines += ['- ' + item for item in result['limitations']]
    output.with_suffix('.md').write_text('\n'.join(lines) + '\n')
    print(json.dumps({arm: {k: value[k] for k in ['strict_point_verdict', 'development_evidence_status']}
                      for arm, value in comparisons.items()}, indent=2), flush=True)


if __name__ == '__main__':
    main()
