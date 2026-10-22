#!/usr/bin/env python3
"""Read-only paired analysis of the preregistered four-arm bulk experiment.

Run after root records terminal results. No scoring, fitting, anchor scaling,
reserve extraction, model edits, or submission is performed by this script.
"""
from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
import yaml

from summarize_dispersion_experiment import (
    METRICS, paired_bootstrap, read_json, real_de_identity,
    reference_call_diagnostics, sha256,
)

ROOT = Path(__file__).resolve().parents[1]
REFERENCE_DIR = ROOT / 'evaluations/dispersion_h1_training_400'
ARMS = ('incumbent', 'common05', 'bulk015', 'bulk060')
MSE = 'expr_mse_unbiased_capped_norm'
PROTECTED = tuple(x for x in METRICS if x != MSE)
TOLERANCE = 1e-8
PROFILES = ('requested', 'projected', 'fitted', 'realized')


def resolve(path):
    path = Path(path)
    return path if path.is_absolute() else ROOT / path


def logbulk(counts_or_probability):
    """Normalize on the complete axis before any target-specific exclusion."""
    values = np.asarray(counts_or_probability, dtype=np.float64)
    if values.ndim != 1 or not np.isfinite(values).all() or (values < 0).any() or values.sum() <= 0:
        raise ValueError('A finite nonnegative nonempty whole-axis bulk is required')
    return np.log1p(50000.0 * values / values.sum())


def profile_error(prediction, reference, control, mask):
    p, r, c = prediction[mask], reference[mask], control[mask]
    response, residual = p - c, r - p
    norms = float(np.linalg.norm(response) * np.linalg.norm(residual))
    return {'uncorrected_mean_error': float(np.mean(residual ** 2)),
            'response_residual_dot_per_gene': float(np.mean(response * residual)),
            'response_residual_cosine': float(np.dot(response, residual) / norms) if norms > 0 else None,
            'response_squared_norm_per_gene': float(np.mean(response ** 2)),
            'n_genes': int(mask.sum())}


def metric_arm(directory):
    meta = read_json(directory / 'run_meta.json')
    config = yaml.safe_load((directory / 'run_params.yaml').read_text())
    backend_ok = meta['resolved_de_backend'] == 'scanpy' or (
        meta['resolved_de_backend'] is None and meta.get('de_real_fingerprint') is not None)
    if (meta['resolved_device'] != 'cpu' or not backend_ok or config['de']['backend'] != 'scanpy'
            or config['control_source'] != 'real' or config['metrics'] != 'vcc2026'
            or config['input_type'] != 'counts' or meta.get('de_pred_fingerprint') is not None
            or resolve(meta['source']).resolve() != (REFERENCE_DIR / 'reference.h5ad').resolve()):
        raise ValueError(f'{directory}: unexpected scoring semantics/reference')
    aggregate = pd.read_csv(directory / 'agg_results.csv').set_index('statistic')
    aggregation = pd.read_csv(directory / 'metric_aggregation.csv').set_index('metric')
    raw = {metric: float(aggregate.loc['mean', metric]) for metric in METRICS}
    if not np.isfinite(list(raw.values())).all():
        raise ValueError(f'{directory}: all six finite raw metrics required')
    return {'raw': raw, 'meta': meta, 'config': config, 'directory': directory,
            'aggregation': aggregation, 'aggregate': aggregate,
            'tidy': pd.read_csv(directory / 'results.csv')}


def values(data, metric, targets):
    part = data['tidy'].loc[data['tidy'].metric == metric, ['perturbation', 'value']]
    if part.perturbation.duplicated().any() or not set(part.perturbation).issubset(targets):
        raise ValueError(f'Invalid target rows for {metric}')
    result = part.set_index('perturbation').value.reindex(targets).to_numpy(dtype=float)
    finite = np.isfinite(result)
    if int(data['aggregation'].loc[metric, 'n_used']) != int(finite.sum()):
        raise ValueError(f'Reported and actual finite target counts differ: {metric}')
    if finite.any() and not np.isclose(result[finite].mean(),
            data['aggregate'].loc['mean', metric], rtol=0, atol=1e-10):
        raise ValueError(f'Per-target and aggregate mean disagree: {metric}')
    return result, finite


def compare_metrics(base, candidate, targets):
    result = {}
    for metric, (label, direction) in METRICS.items():
        item = {'label': label, 'incumbent': base['raw'][metric], 'candidate': candidate['raw'][metric],
                'gain_in_better_direction': direction * (candidate['raw'][metric] - base['raw'][metric]),
                'better_direction': 'higher' if direction == 1 else 'lower', 'paired_bootstrap': None}
        if metric == MSE:
            components = {}
            for name in ('expr_mse_unbiased_capped', 'expr_distance_unbiased'):
                a, fa = values(base, name, targets)
                b, fb = values(candidate, name, targets)
                components[name] = {'identical_finite_cohort': bool(np.array_equal(fa, fb)),
                    'incumbent_targets': np.asarray(targets)[fa].tolist(),
                    'candidate_targets': np.asarray(targets)[fb].tolist()}
                # The complete registered panel is needed for the shared budget.
                if not fa.all() or not fb.all():
                    raise ValueError('MSE component cohort is not the complete frozen panel')
                if name == 'expr_distance_unbiased':
                    np.testing.assert_allclose(a, b, rtol=0, atol=1e-12)
            for data in (base, candidate):
                num, _ = values(data, 'expr_mse_unbiased_capped', targets)
                den, _ = values(data, 'expr_distance_unbiased', targets)
                if den.sum() <= 0 or not np.isclose(num.sum() / den.sum(), data['raw'][MSE], rtol=0, atol=1e-10):
                    raise ValueError('Global MSE ratio does not match complete finite component cohorts')
            item.update(finite_cohort_identical=True, component_cohorts=components,
                uncertainty_note='Global ratio with whole-panel correction budget: no naive per-target MSE bootstrap.')
        else:
            a, fa = values(base, metric, targets)
            b, fb = values(candidate, metric, targets)
            identical = np.array_equal(fa, fb) and bool(fa.any())
            item.update(finite_cohort_identical=bool(identical),
                incumbent_finite_targets=np.asarray(targets)[fa].tolist(),
                candidate_finite_targets=np.asarray(targets)[fb].tolist())
            if identical:
                item['paired_bootstrap'] = paired_bootstrap(direction * (b[fa] - a[fa]))
            else:
                item['uncertainty_note'] = 'Different finite cohorts fail comparability; no intersection-only bootstrap.'
            item['point_nonregression'] = bool(item['gain_in_better_direction'] >= -TOLERANCE)
        result[metric] = item
    return result


def source_strata(config, targets, genes):
    result = {t: {'eligible_raw_sources': [], 'cd4_available': False} for t in targets}
    for source in config['source']:
        if 'h1' in source['path'].lower():
            raise ValueError('Forbidden H1 source')
        with np.load(resolve(source['path']), allow_pickle=False) as d:
            measured = set(d['genes'].astype(str)[d['measured_genes'].astype(bool)]) & set(genes)
            eligible = set(d['targets'].astype(str)[d['n_cells'] >= 20]) if measured else set()
        for target in targets:
            if target in eligible and source['weight'] > 0:
                result[target]['eligible_raw_sources'].append(source['path'])
    with np.load(resolve(config['cd4']), allow_pickle=False) as d:
        columns = np.isin(d['genes'].astype(str), genes)
        available = d['available'] & d['quality_pass'] & (d['n_cells'] >= 20)
        finite_gene = np.isfinite(d['log2fc'][:, :, columns]).any(axis=2)
        eligible = set(d['targets'].astype(str)[(available & finite_gene).any(axis=0)])
    for target in targets:
        result[target]['raw_source_count'] = len(result[target]['eligible_raw_sources'])
        result[target]['cd4_available'] = target in eligible and config['cd4_weight'] > 0
    return result


def reference_profiles(targets):
    real = ad.read_h5ad(REFERENCE_DIR / 'reference.h5ad', backed='r')
    try:
        labels = real.obs.target_gene.astype(str).to_numpy()
        if real.shape != (24000, 18080) or set(labels) != set(targets) | {'non-targeting'}:
            raise ValueError('Reference axes/group set changed')
        genes = real.var_names.astype(str).to_numpy()
        if len(set(genes)) != len(genes):
            raise ValueError('Duplicate reference genes')
        profiles = {}
        for target in [*targets, 'non-targeting']:
            rows = np.flatnonzero(labels == target)
            expected = 4000 if target == 'non-targeting' else 400
            if len(rows) != expected:
                raise ValueError('Reference group count changed')
            sums = np.asarray(real.X[rows].sum(axis=0)).ravel().astype(np.float64)
            profiles[target] = logbulk(sums)
    finally:
        real.file.close()
    return genes, profiles


def final_count_audit(arm, spec, targets, genes, reference, baseline_profiles):
    path = resolve(spec['output'])
    directory = path.parent
    receipt = read_json(path.with_suffix('.receipt.json'))
    validation = read_json(directory / 'validation.json')
    if (not validation.get('valid') or not validation.get('hash_verified')
            or validation.get('excluded_source_tokens') != ['H1']
            or validation.get('shape') != [20000, 18080]
            or receipt['targets'] != targets or receipt['arguments'].get('target_limit') is not None):
        raise ValueError(f'{arm}: complete validated development prediction required')
    if sha256(path) != receipt['output_sha256']:
        raise ValueError(f'{arm}: prediction bytes changed')
    audit = read_json(directory / 'profile_audit.json')
    if len(audit) != 50 or {(a['context'], a['target']) for a in audit} != {('H1', t) for t in targets}:
        raise ValueError(f'{arm}: incomplete/duplicate target diagnostics')
    by_target = {a['target']: a for a in audit}
    table = pd.read_parquet(directory / 'gene_profile_audit.parquet')
    if table.duplicated(['context', 'target', 'gene']).any() or len(table) != 50 * len(genes):
        raise ValueError('Invalid per-gene diagnostic axes')
    for name in ('gene_profile_audit', 'profile_audit'):
        suffix = '.parquet' if name == 'gene_profile_audit' else '.json'
        if sha256(directory / (name + suffix)) != receipt['bulk_ablation'][name + '_sha256']:
            raise ValueError('Profile artifact hash changed')
    predicted = ad.read_h5ad(path, backed='r')
    old = ad.read_h5ad(ROOT / 'outputs/dispersion_h1_400_incumbent/prediction.h5ad', backed='r')
    scored = ad.read_h5ad(resolve(spec['metrics']) / 'local_scoring_input.h5ad', backed='r')
    groups = []
    requested_means = {}
    try:
        if (predicted.shape != old.shape or not predicted.obs_names.equals(old.obs_names)
                or not predicted.var_names.equals(old.var_names)
                or not np.array_equal(predicted.var_names.astype(str), genes)
                or not np.array_equal(predicted.obs[['context', 'target_gene']].astype(str),
                                      old.obs[['context', 'target_gene']].astype(str))
                or scored.shape != (24000, len(genes)) or not scored.var_names.equals(old.var_names)
                or not np.array_equal(scored.obs.target_gene.astype(str).to_numpy()[:20000],
                                      predicted.obs.target_gene.astype(str))):
            raise ValueError(f'{arm}: cell, gene or scoring axes differ')
        labels = predicted.obs.target_gene.astype(str).to_numpy()
        if set(labels) != set(targets):
            raise ValueError('Prediction target panel changed')
        for target in targets:
            rows = np.flatnonzero(labels == target)
            if len(rows) != 400:
                raise ValueError('Incomplete target group')
            x = predicted.X[rows].toarray()
            original = old.X[rows].toarray()
            if not np.isfinite(x).all() or (x < 0).any() or (x != np.floor(x)).any():
                raise ValueError('Nonintegral, negative or nonfinite emitted counts')
            if (sparse.csr_matrix(scored.X[rows]) != sparse.csr_matrix(x)).nnz:
                raise ValueError('Scored treated counts differ from prediction')
            depths = x.sum(axis=1, dtype=np.int64)
            old_depths = original.sum(axis=1, dtype=np.int64)
            if not np.array_equal(depths, old_depths) or (depths < 1).any() or (depths > 1000000).any():
                raise ValueError('Frozen row depths changed or violate caps')
            identical = np.array_equal(x, original)
            if arm == 'incumbent' and not identical:
                raise ValueError('New full incumbent is not bit-identical to the existing incumbent')
            p = table[(table.context == 'H1') & (table.target == target)].set_index('gene').reindex(genes)
            if p.isna().any().any():
                raise ValueError('Incomplete/nonfinite per-gene diagnostic table')
            sums = x.sum(axis=0, dtype=np.int64)
            mean = (x / depths[:, None]).mean(axis=0)
            before = (original / depths[:, None]).mean(axis=0)
            bulk = sums / sums.sum()
            np.testing.assert_array_equal(p.gene_count_sums_realized.to_numpy(), sums)
            np.testing.assert_allclose(p.mean_probability_realized, mean, rtol=1e-12, atol=1e-14)
            np.testing.assert_allclose(p.bulk_probability_realized, bulk, rtol=1e-12, atol=1e-14)
            np.testing.assert_allclose(p.mean_cpm_change_to_incumbent, 1e6 * (mean - before), rtol=1e-10, atol=1e-8)
            requested_means[target] = p.mean_probability_requested.to_numpy()
            if baseline_profiles is not None:
                np.testing.assert_array_equal(requested_means[target], baseline_profiles[target])
            a = by_target[target]
            if a['projection_l1'] > .03 or max(a['fitting_mean_l1'], a['fitting_bulk_l1']) > .001:
                raise ValueError('Original emitter feasibility limit exceeded')
            mask = genes != target
            own = len(genes) - int(mask.sum())
            errors = {}
            for stage in PROFILES:
                profile = bulk if stage == 'realized' else p['bulk_probability_' + stage].to_numpy()
                if not np.isclose(profile.sum(), 1, rtol=0, atol=1e-8):
                    raise ValueError('Non-normalized diagnostic bulk probability')
                errors[stage] = profile_error(logbulk(profile), reference[target], reference['non-targeting'], mask)
            groups.append({'target': target, 'own_gene_excluded': bool(own), 'rows_exact': True,
                'counts_bitidentical_to_incumbent': bool(identical),
                'gene_totals_identical_to_incumbent': bool(np.array_equal(sums, original.sum(axis=0, dtype=np.int64))),
                'max_abs_realized_mean_cpm_drift': float(np.abs(1e6 * (mean - before)).max()),
                'mean_probability_l1_to_incumbent': float(np.abs(mean - before).sum()),
                'bulk_probability_l1_to_incumbent': float(np.abs(bulk - p.bulk_probability_incumbent).sum()),
                'projection_l1': float(a['projection_l1']),
                'fitting_mean_l1': float(a['fitting_mean_l1']), 'fitting_bulk_l1': float(a['fitting_bulk_l1']),
                'profiles': errors})
    finally:
        predicted.file.close(); old.file.close(); scored.file.close()
    profile_means = {stage: float(np.mean([g['profiles'][stage]['uncorrected_mean_error'] for g in groups]))
                     for stage in PROFILES}
    alignment = {}
    for stage in PROFILES:
        cosines = [g['profiles'][stage]['response_residual_cosine'] for g in groups
                   if g['profiles'][stage]['response_residual_cosine'] is not None]
        alignment[stage] = {
            'mean_response_residual_dot_per_gene': float(np.mean([
                g['profiles'][stage]['response_residual_dot_per_gene'] for g in groups])),
            'mean_finite_response_residual_cosine': float(np.mean(cosines)) if cosines else None,
            'finite_cosine_targets': len(cosines)}
    return {'passed': True, 'independently_recomputed_from_final_counts': True,
            'output_sha256': receipt['output_sha256'], 'groups': groups,
            'profile_mean_errors': profile_means, 'profile_alignment': alignment,
            'max_realized_mean_cpm_drift': max(g['max_abs_realized_mean_cpm_drift'] for g in groups),
            'maximum_projection_l1': max(g['projection_l1'] for g in groups),
            'context_profiles': receipt['bulk_ablation']['context_profiles'],
            'receipt_sha256': sha256(path.with_suffix('.receipt.json'))}, requested_means


def infeasible_arm(arm, status, spec):
    if status.get('exit_code') in (None, 0):
        raise ValueError('Infeasible arm must have an authoritative nonzero generation exit')
    path = resolve(status.get('failure_json', str(Path(spec['output']).parent / 'failure.json')))
    failure = read_json(path)
    group = failure.get('failed_group', {})
    error = group.get('exception', '')
    allowed = ('Requested bulk projection too large:', 'Moment fitting failed:',
               'Bulk support cannot satisfy the feasible moment bounds', 'Cannot bracket the bulk projection')
    if (failure.get('arm') != arm or not error.startswith(allowed)
            or failure.get('unchanged_projection_l1_limit') != .03
            or failure.get('unchanged_fit_l1_limit') != .001):
        raise ValueError('Failure is not an attributable unchanged-emitter feasibility rejection')
    if (resolve(spec['metrics']) / 'agg_results.csv').exists():
        raise ValueError('An infeasible/partial arm must not have a scored aggregate')
    return {'status': 'GENERATION_INFEASIBLE', 'strict_point_verdict': 'GENERATION_INFEASIBLE',
            'eligible_for_reserved_confirmation': False, 'failure': failure,
            'failure_sha256': sha256(path), 'terminal_status': status,
            'interpretation': 'Pinned emitter/settings failed; no complete score, replacement, or tolerance relaxation.'}


def render(result):
    base = result['incumbent']['raw_metrics']
    lines = ['# Bulk-response development results', '',
        'Fixed 50-target H1 Training panel, 400 real/predicted cells per target and 4,000 real controls. '
        'Values are raw official metrics. No overall raw score or cross-panel anchor scaling.', '',
        '| Metric | Incumbent | common05 | bulk015 | bulk060 | Better |',
        '| --- | ---: | ---: | ---: | ---: | --- |']
    for metric, (label, direction) in METRICS.items():
        vals = [f'{base[metric]:.9f}']
        for arm in ARMS[1:]:
            comp = result['comparisons'][arm]
            vals.append(f'{comp["metrics"][metric]["candidate"]:.9f}' if 'metrics' in comp else 'infeasible')
        lines.append('| ' + ' | '.join([label, *vals, 'higher' if direction == 1 else 'lower']) + ' |')
    lines += ['', '| Profile error (lower better) | Incumbent | common05 | bulk015 | bulk060 |',
              '| --- | ---: | ---: | ---: | ---: |']
    for stage in PROFILES:
        vals = [f'{result["incumbent"]["matrix_audit"]["profile_mean_errors"][stage]:.9f}']
        for arm in ARMS[1:]:
            comp = result['comparisons'][arm]
            vals.append(f'{comp["matrix_audit"]["profile_mean_errors"][stage]:.9f}' if 'matrix_audit' in comp else 'infeasible')
        lines.append('| ' + ' | '.join([stage, *vals]) + ' |')
    for arm, comp in result['comparisons'].items():
        lines += ['', f'## {arm}: {comp["strict_point_verdict"]}', '']
        if comp['status'] == 'GENERATION_INFEASIBLE':
            lines += [comp['interpretation'], '']
            continue
        ci = comp['mean_error_paired_bootstrap']
        lines += [f'Realized mean-error gain: {ci["mean_gain"]:+.9f}; paired 95% interval '
            f'[{ci["percentile_95_ci"][0]:+.9f}, {ci["percentile_95_ci"][1]:+.9f}]. '
            f'All five protected metrics nonworse: {comp["protected_metrics_nonworse"]}. '
            f'Identical finite cohorts: {comp["all_finite_target_cohorts_identical"]}. '
            f'Eligible for reserved confirmation: {comp["eligible_for_reserved_confirmation"]}.', '',
            '| Metric | Gain in better direction | Paired 95% interval |', '| --- | ---: | --- |']
        for metric, item in comp['metrics'].items():
            b = item['paired_bootstrap']
            interval = f'[{b["percentile_95_ci"][0]:+.9f}, {b["percentile_95_ci"][1]:+.9f}]' if b else item['uncertainty_note']
            lines.append(f'| {item["label"]} | {item["gain_in_better_direction"]:+.9f} | {interval} |')
        correction = comp['mse_decomposition']
        lines += ['', f'Maximum realized CPM drift: {comp["matrix_audit"]["max_realized_mean_cpm_drift"]:.6g}. '
            f'Mean uncorrected error change: {correction["delta_uncorrected_error"]:+.9g}; '
            f'capped numerator change: {correction["delta_capped_numerator"]:+.9g}; '
            f'effective prediction-correction change: {correction["delta_effective_prediction_correction"]:+.9g}. '
            'Both mean response and the global correction budget may move.']
        alignment = comp['matrix_audit']['profile_alignment']['realized']
        baseline_alignment = result['incumbent']['matrix_audit']['profile_alignment']['realized']
        lines += ['', f'Mean realized response/residual dot per gene: '
            f'{baseline_alignment["mean_response_residual_dot_per_gene"]:+.9g} to '
            f'{alignment["mean_response_residual_dot_per_gene"]:+.9g}. '
            'Positive alignment is a local directional diagnostic, not a fitted amplitude or an extra selection gate.']
    lines += ['', f'Deterministic reserve selection: {result["reserved_confirmation"]["selected_arm"] or "none"}. '
              + result['reserved_confirmation']['action'], '',
              'The JSON retains actual finite target IDs, source-coverage strata, response/residual alignment, '
              'strict DE calls/TP/FP/FN, terminal statuses and independent final-count/provenance checks.', '']
    lines += ['- ' + note for note in result['limitations']]
    return '\n'.join(lines) + '\n'


def main():
    manifest = read_json(ROOT / 'configs/bulk_h1_400_experiment.json')
    terminal = read_json(ROOT / 'reports/bulk_terminal_statuses.json')
    statuses = terminal.get('arms', terminal)
    freeze = read_json(ROOT / 'reports/bulk_implementation_freeze.json')
    for name, expected in freeze['files'].items():
        if sha256(resolve(name)) != expected:
            raise ValueError(f'Frozen artifact changed: {name}')
    cohort = read_json(REFERENCE_DIR / 'cohort_freeze.json')
    targets = cohort['development_targets']
    if len(targets) != 50 or set(targets) & set(cohort['reserved_targets_not_extracted']):
        raise ValueError('Development/reserve partition changed')
    repaired = read_json(REFERENCE_DIR / 'manifest.json')
    if sha256(REFERENCE_DIR / 'reference.h5ad') != repaired['reference_sha256']:
        raise ValueError('Reference bytes changed')
    identity = read_json(ROOT / 'reports/bulk_baseline_identity.json')
    if not identity['counts_bit_identical'] or not identity['axes_identical']:
        raise ValueError('Full incumbent reproduction is required for old metric reuse')
    if (identity['old_incumbent_sha256'] != sha256(resolve(manifest['existing_incumbent']))
            or identity['new_incumbent_sha256'] != sha256(resolve(manifest['arms']['incumbent']['output']))
            or identity['new_receipt_sha256'] != sha256(resolve(manifest['arms']['incumbent']['output']).with_suffix('.receipt.json'))
            or identity['builder_sha256'] != freeze['files']['scripts/build_bulk_ablation.py']):
        raise ValueError('Incumbent identity evidence changed')
    if statuses['incumbent'].get('exit_code') != 0:
        raise ValueError('Authoritative completed incumbent scoring status required')
    cache_path = ROOT / 'reports/bulk_real_de_cache.json'
    cache = read_json(cache_path)
    base = metric_arm(resolve(manifest['arms']['incumbent']['metrics']))
    for key, file in [('real_de_sha256', 'de_real.parquet'), ('baseline_run_meta_sha256', 'run_meta.json'),
                      ('baseline_run_params_sha256', 'run_params.yaml')]:
        if cache[key] != sha256(base['directory'] / file):
            raise ValueError('Frozen baseline/cache artifact changed')
    if (cache['identity_audit_sha256'] != sha256(ROOT / 'reports/bulk_baseline_identity.json')
            or cache['implementation_freeze_sha256'] != sha256(ROOT / 'reports/bulk_implementation_freeze.json')
            or cache['reference_sha256'] != repaired['reference_sha256']):
        raise ValueError('Bulk cache/freeze identity mismatch')
    genes, reference = reference_profiles(targets)
    config = read_json(resolve(manifest['base_config']))
    strata = source_strata(config, targets, genes)
    base_audit, base_requested = final_count_audit('incumbent', manifest['arms']['incumbent'],
                                                targets, genes, reference, None)
    comparisons = {}
    for arm in ARMS[1:]:
        status, spec = statuses[arm], manifest['arms'][arm]
        if status.get('status') == 'GENERATION_INFEASIBLE':
            comparisons[arm] = infeasible_arm(arm, status, spec)
            continue
        if status.get('exit_code') != 0:
            raise ValueError(f'{arm}: pending/failed metrics, not formally infeasible')
        data = metric_arm(resolve(spec['metrics']))
        for key in ('source_fingerprint', 'anchor_semantic_identity', 'resolved_device', 'cell_eval2_version'):
            if data['meta'][key] != base['meta'][key]:
                raise ValueError(f'{arm}: reference/scorer identity mismatch: {key}')
        ignored = {'outdir', 'cache_real', 'cache_pred'}
        if ({k: v for k, v in data['config'].items() if k not in ignored}
                != {k: v for k, v in base['config'].items() if k not in ignored}):
            raise ValueError(f'{arm}: scoring configuration mismatch')
        de = real_de_identity(data['directory'], base['directory'])
        preflight = read_json(data['directory'] / 'preflight.json')
        if (data['meta']['de_real_fingerprint'] != de['baseline_strict_fingerprint']
                or preflight['real_de_cache_sha256'] != sha256(cache_path)
                or preflight['real_de_sha256'] != cache['real_de_sha256']
                or preflight['reference_sha256'] != repaired['reference_sha256']
                or preflight['prediction_sha256'] != sha256(resolve(spec['output']))
                or preflight['local_scoring_input_sha256'] != sha256(data['directory'] / 'local_scoring_input.h5ad')
                or not preflight['scored_treated_and_control_counts_verified']):
            raise ValueError(f'{arm}: supplied real-DE/preflight/scoring-copy mismatch')
        audit, _ = final_count_audit(arm, spec, targets, genes, reference, base_requested)
        for context, declared in audit['context_profiles'].items():
            for key in ('requested_mean_sha256', 'template_sha256', 'depths_sha256',
                        'cd4_centered_lfc_sha256', 'cd4_available_sha256'):
                if declared[key] != base_audit['context_profiles'][context][key]:
                    raise ValueError(f'{arm}: isolation invariant changed: {key}')
        metrics = compare_metrics(base, data, targets)
        before = np.asarray([g['profiles']['realized']['uncorrected_mean_error'] for g in base_audit['groups']])
        after = np.asarray([g['profiles']['realized']['uncorrected_mean_error'] for g in audit['groups']])
        paired = paired_bootstrap(before - after)
        delta_u = float(after.mean() - before.mean())
        num_b = float(base['aggregate'].loc['mean', 'expr_mse_unbiased_capped'])
        num_a = float(data['aggregate'].loc['mean', 'expr_mse_unbiased_capped'])
        mse_decomposition = {'incumbent_uncorrected_error': float(before.mean()),
            'candidate_uncorrected_error': float(after.mean()), 'delta_uncorrected_error': delta_u,
            'delta_capped_numerator': num_a - num_b,
            'delta_effective_prediction_correction': delta_u - (num_a - num_b),
            'fixed_reference_denominator_mean': float(base['aggregate'].loc['mean', 'expr_distance_unbiased']),
            'interpretation': 'Candidate minus incumbent; reference correction is fixed. Both predicted mean error and the global spread/capped correction may change.'}
        cohorts_ok = all(x['finite_cohort_identical'] for x in metrics.values())
        protected = all(metrics[m]['point_nonregression'] for m in PROTECTED)
        primary = metrics[MSE]['gain_in_better_direction'] > 0 and paired['mean_gain'] > 0
        strict = audit['passed'] and cohorts_ok and protected and primary
        eligible = strict and paired['percentile_95_ci'][0] > 0
        stratum_results = {}
        for count in range(4):
            selected = np.asarray([strata[t]['raw_source_count'] == count for t in targets])
            stratum_results[str(count)] = {'n_targets': int(selected.sum()),
                'targets': np.asarray(targets)[selected].tolist(),
                'incumbent_mean_error': float(before[selected].mean()) if selected.any() else None,
                'candidate_mean_error': float(after[selected].mean()) if selected.any() else None,
                'mean_error_gain': float((before[selected] - after[selected]).mean()) if selected.any() else None}
        comparisons[arm] = {'status': 'SCORED', 'strict_point_verdict': 'PASS STRICT' if strict else 'FAIL STRICT',
            'protected_metrics_nonworse': bool(protected), 'both_primary_point_gains_positive': bool(primary),
            'all_finite_target_cohorts_identical': bool(cohorts_ok),
            'mean_error_paired_bootstrap': paired, 'eligible_for_reserved_confirmation': bool(eligible),
            'promoted': False, 'metrics': metrics, 'matrix_audit': audit, 'mse_decomposition': mse_decomposition,
            'source_coverage_strata': stratum_results, 'preflight': preflight,
            'real_de_identity': de, 'terminal_status': status,
            'de_call_diagnostics': reference_call_diagnostics(data['directory'], targets)}
    eligible = [a for a in ARMS[1:] if comparisons[a]['eligible_for_reserved_confirmation']]
    selected = min(eligible, key=lambda a: (comparisons[a]['metrics'][MSE]['candidate'],
        comparisons[a]['matrix_audit']['profile_mean_errors']['realized'], ARMS.index(a))) if eligible else None
    result = {'protocol': 'docs/BULK_ABLATION_PROTOCOL_20261021.md',
        'protocol_sha256': sha256(ROOT / 'docs/BULK_ABLATION_PROTOCOL_20261021.md'),
        'analysis_sha256': sha256(__file__), 'implementation_freeze': freeze,
        'analysis_helper_sha256': sha256(ROOT / 'scripts/summarize_dispersion_experiment.py'),
        'terminal_statuses': terminal, 'cohort_freeze': cohort, 'baseline_identity': identity,
        'reference_sha256': repaired['reference_sha256'], 'real_de_cache_sha256': sha256(cache_path),
        'source_coverage_definitions': {'eligible_raw': 'target n_cells>=20 and at least one measured output gene, positive frozen weight',
                                      'cd4': 'available and quality_pass and n_cells>=20 and a finite measured output gene in any condition'},
        'source_coverage_per_target': strata,
        'incumbent': {'raw_metrics': base['raw'], 'matrix_audit': base_audit,
                      'de_call_diagnostics': reference_call_diagnostics(base['directory'], targets)},
        'comparisons': comparisons, 'overall_raw_score_computed': False, 'anchors_used': False,
        'reserved_confirmation': {'selected_arm': selected, 'eligible_arms': eligible,
            'selection_rule': 'lowest raw normalized MSE, then lower realized uncorrected mean error, then registered arm order',
            'action': 'Freeze one selected recipe before the authorized single 25-target confirmation; no retuning or official submission.' if selected
                      else 'No candidate qualifies: retain incumbent, leave the reserve unextracted and stop.'},
        'limitations': ['Previously observed single-context Training development cohort; this is not untouched biological validation.',
            'Paired target intervals are conditional/descriptive with shared controls and correlated perturbations; no multiplicity adjustment.',
            'Point nonregression within 1e-8 is not statistical proof of no harm.',
            'Fixed requested arithmetic profiles do not ensure identical realized means or DE tests.',
            'Four thousand controls and the native 18,080-gene axis differ from the official challenge; no leaderboard gain is inferred.',
            'The global MSE correction is not bootstrapped as independent per-target normalized errors; mean residual and correction effects are separated.',
            'No model is promoted or submitted by this analysis. Closed H1 Test/K562 confirmations remain closed.']}
    path = ROOT / 'reports/bulk_ablation_results.json'
    path.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    path.with_suffix('.md').write_text(render(result))
    print(json.dumps({'comparisons': {a: {k: c[k] for k in ('strict_point_verdict', 'eligible_for_reserved_confirmation')}
                                      for a, c in comparisons.items()},
                      'reserved_confirmation': result['reserved_confirmation']}, indent=2), flush=True)


if __name__ == '__main__':
    main()
