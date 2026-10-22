#!/usr/bin/env python3
"""Analyze the frozen 2x2 response/emitter experiment after all eight jobs exit 0.

No scoring, training, second-context extraction, retuning, or submission. Each
recipe is compared with the current/current baseline at the SAME emission seed.
Both seeds must satisfy the unchanged prespecified point gate.
"""
from __future__ import annotations
import os
for _key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_key, '1')
from pathlib import Path
import json

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from build_prediction import atlas, control_template
from build_bulk_ablation import array_sha256
from summarize_bulk_ablation import metric_arm, compare_metrics, reference_profiles, profile_error, logbulk
from summarize_dispersion_experiment import (METRICS, paired_bootstrap, read_json, real_de_identity,
                                            reference_call_diagnostics, sha256)

ROOT = Path(__file__).resolve().parents[1]
REFERENCE = ROOT / 'evaluations/dispersion_h1_training_400'
OUTPUTS = ROOT / 'outputs/redesign_h1_400'
METRIC_ROOT = ROOT / 'evaluations/redesign_h1_400'
SEEDS = (20261021, 20261022)
RECIPES = ('current_current', 'current_preserving', 'adapter_current', 'adapter_preserving')
TOLERANCE = 1e-8
PRIMARY = ('de_wilcoxon_direction_fidelity_yield_raw', 'de_wilcoxon_sig_jaccard')
MSE = 'expr_mse_unbiased_capped_norm'


def arm_name(recipe, seed):
    return f'{recipe}_s{seed}'


def metric_directory(recipe, seed):
    if recipe == 'current_current' and seed == SEEDS[0]:
        return REFERENCE / 'incumbent_cpu'
    return METRIC_ROOT / arm_name(recipe, seed)


def require_all_terminals(metric_root=METRIC_ROOT, aggregate_path=None):
    """This must run before any metric or treated-response performance reads."""
    records = {}
    for recipe in RECIPES:
        for seed in SEEDS:
            name = arm_name(recipe, seed)
            path = metric_root / name / 'terminal.json'
            record = read_json(path)
            if record.get('exit_code') != 0 or record.get('arm') != name:
                raise ValueError(f'{name}: authoritative successful terminal required')
            records[name] = record
    if aggregate_path is not None and Path(aggregate_path).exists():
        aggregate = read_json(aggregate_path)
        rows = aggregate.get('arms', aggregate)
        for name in records:
            if rows.get(name, {}).get('exit_code') != 0:
                raise ValueError(f'{name}: aggregate terminal status is not successful')
    return records


def point_gate(comparison):
    if set(comparison) != set(METRICS):
        raise ValueError('Exactly all six registered raw metrics are required')
    comparable = all(item['finite_cohort_identical'] for item in comparison.values())
    finite = all(np.isfinite(item['gain_in_better_direction']) for item in comparison.values())
    nonregression = finite and all(item['gain_in_better_direction'] >= -TOLERANCE for item in comparison.values())
    improved = finite and any(comparison[m]['gain_in_better_direction'] > TOLERANCE for m in PRIMARY)
    return {'passed': bool(comparable and nonregression and improved),
            'all_six_nonregression': bool(nonregression), 'identical_finite_cohorts': bool(comparable),
            'fidelity_or_jaccard_gain_gt_tolerance': bool(improved), 'tolerance': TOLERANCE}


def recipe_gate(seed_gates):
    if set(seed_gates) != set(SEEDS):
        raise ValueError('Both preregistered emission seeds are required')
    return bool(all(seed_gates[seed]['passed'] for seed in SEEDS))


def verify_donor_metadata(audit, target, seed, controls, emitted_depths):
    donor_seed = atlas.seed_for('donor:H1:' + target, seed)
    emission_seed = atlas.seed_for('transport:H1:' + target, seed)
    rows = np.random.default_rng(donor_seed).choice(controls.shape[0], 400, replace=False)
    supplied = np.asarray(audit['donor_rows'])
    if (audit['donor_seed'] != donor_seed or audit['emission_seed'] != emission_seed
            or not np.array_equal(rows, supplied) or len(set(rows)) != 400
            or audit['donor_rows_sha256'] != array_sha256(rows)):
        raise ValueError('Preserving donor rows/seeds differ from deterministic frozen selection')
    selected = sparse.csr_matrix(controls)[rows].toarray().astype(np.int32)
    if audit['donor_count_sha256'] != array_sha256(selected):
        raise ValueError('Preserving donor counts differ from the frozen control population')
    if not np.array_equal(emitted_depths, selected.sum(axis=1, dtype=np.int64)):
        raise ValueError('Preserving emission changed measured donor library depths')
    return selected


def mean_errors(probability, reference, control, genes, target):
    transformed = logbulk(probability)
    return {'own_excluded': profile_error(transformed, reference, control, genes != target),
            'whole_axis': profile_error(transformed, reference, control, np.ones(len(genes), dtype=bool))}


def verify_metric_identity(data, baseline, name, terminal, freeze_hash, old_cache, prediction_hash):
    ignored = {'outdir', 'cache_real', 'cache_pred'}
    if ({k: v for k, v in data['config'].items() if k not in ignored}
            != {k: v for k, v in baseline['config'].items() if k not in ignored}):
        raise ValueError(f'{name}: scoring configuration differs from the baseline')
    for key in ('source_fingerprint', 'anchor_semantic_identity', 'resolved_device', 'cell_eval2_version'):
        if data['meta'][key] != baseline['meta'][key]:
            raise ValueError(f'{name}: scorer/reference identity mismatch: {key}')
    if name == arm_name('current_current', SEEDS[0]):
        for key, file in [('local_scoring_input_sha256', 'local_scoring_input.h5ad'),
                          ('baseline_run_meta_sha256', 'run_meta.json'),
                          ('baseline_run_params_sha256', 'run_params.yaml'), ('real_de_sha256', 'de_real.parquet')]:
            if old_cache[key] != sha256(data['directory'] / file):
                raise ValueError('Reused incumbent scoring artifact changed')
        return {'baseline_reused': True, 'terminal': terminal}
    preflight = read_json(data['directory'] / 'preflight.json')
    de = real_de_identity(data['directory'], baseline['directory'])
    if (terminal.get('preflight_sha256') != sha256(data['directory'] / 'preflight.json')
            or preflight['arm'] != name or preflight['freeze_sha256'] != freeze_hash
            or preflight['entrypoint_sha256'] != sha256(ROOT / 'scripts/run_redesign_metrics.py')
            or terminal['entrypoint_sha256'] != preflight['entrypoint_sha256']
            or preflight['prediction_sha256'] != prediction_hash
            or preflight['reference_sha256'] != old_cache['reference_sha256']
            or preflight['local_scoring_input_sha256'] != sha256(data['directory'] / 'local_scoring_input.h5ad')
            or preflight['local_scoring_receipt_sha256'] != sha256(data['directory'] / 'local_scoring_input.receipt.json')
            or preflight['real_de_sha256'] != old_cache['real_de_sha256']
            or not preflight['scored_treated_and_control_counts_verified']
            or data['meta']['de_real_fingerprint'] != de['baseline_strict_fingerprint']):
        raise ValueError(f'{name}: supplied DE, scored input or terminal provenance mismatch')
    return {'preflight': preflight, 'real_de_identity': de, 'terminal': terminal}


def audit_counts(recipe, seed, targets, genes, reference, controls, expected_depths, profiles,
                 frozen, baseline_path, metrics_directory):
    name = arm_name(recipe, seed); response, emitter = recipe.split('_')
    directory = OUTPUTS / name; path = directory / 'prediction.h5ad'
    receipt = read_json(path.with_suffix('.receipt.json')); spec = receipt['redesign']
    validation = read_json(directory / 'validation.json')
    if (not validation['valid'] or not validation['hash_verified'] or validation['shape'] != [20000, len(genes)]
            or validation['excluded_source_tokens'] != ['H1'] or receipt['targets'] != targets
            or spec['engineering_subset'] or spec['seed'] != seed or spec['response'] != response
            or spec['emitter'] != emitter or receipt['arguments'].get('target_limit') is not None
            or spec['profiles_sha256'] != frozen['files']['models/response_adapter_v1/profiles.npz']
            or spec['profile_receipt_sha256'] != frozen['files']['models/response_adapter_v1/profiles.json']):
        raise ValueError(f'{name}: incomplete or mismatched prediction/response/emitter')
    for code, digest in spec['code_hashes'].items():
        if frozen['files'].get(code) != digest:
            raise ValueError(f'{name}: unregistered generation code')
    for source in receipt['sources']:
        if 'h1' in source['path'].lower() or sha256(source['path']) != source['sha256']:
            raise ValueError(f'{name}: forbidden/altered source bank')
    actual_hash = sha256(path)
    if receipt['output_sha256'] != actual_hash:
        raise ValueError(f'{name}: prediction bytes changed')
    for file, key in [('profile_audit.json', 'profile_audit_sha256'),
                      ('gene_profile_audit.parquet', 'gene_profile_audit_sha256')]:
        if sha256(directory / file) != spec[key]:
            raise ValueError(f'{name}: generation diagnostic artifact changed')
    audits = read_json(directory / 'profile_audit.json')
    if len(audits) != 50 or {(a['context'], a['target']) for a in audits} != {('H1', t) for t in targets}:
        raise ValueError('Missing or duplicate target-level audit')
    audit_by_target = {a['target']: a for a in audits}
    table = pd.read_parquet(directory / 'gene_profile_audit.parquet')
    if len(table) != len(targets) * len(genes) or table.duplicated(['context', 'target', 'gene']).any():
        raise ValueError('Invalid diagnostic gene axes')
    predicted = ad.read_h5ad(path, backed='r'); before = ad.read_h5ad(baseline_path, backed='r')
    scored = ad.read_h5ad(metrics_directory / 'local_scoring_input.h5ad', backed='r')
    groups, drift_tables, total_nnz = [], [], 0
    try:
        labels = predicted.obs.target_gene.astype(str).to_numpy()
        if (predicted.shape != (20000, len(genes)) or set(labels) != set(targets)
                or not np.array_equal(predicted.var_names.astype(str), genes)
                or not predicted.obs_names.equals(before.obs_names) or not predicted.var_names.equals(before.var_names)
                or not np.array_equal(predicted.obs[['context', 'target_gene']].astype(str), before.obs[['context', 'target_gene']].astype(str))
                or scored.shape != (24000, len(genes)) or not scored.var_names.equals(predicted.var_names)
                or not np.array_equal(scored.obs.target_gene.astype(str).to_numpy()[:20000], labels)
                or not np.all(scored.obs.target_gene.astype(str).to_numpy()[20000:] == 'non-targeting')):
            raise ValueError(f'{name}: count or semantic axes differ')
        for start in range(0, 4000, 512):
            end = min(start + 512, 4000)
            if (sparse.csr_matrix(scored.X[20000+start:20000+end]) != controls[start:end]).nnz:
                raise ValueError('Scored controls differ from frozen measured controls')
        for i, target in enumerate(targets):
            rows = np.flatnonzero(labels == target)
            if len(rows) != 400:
                raise ValueError('Every target must have 400 cells')
            x = predicted.X[rows].toarray(); old = before.X[rows].toarray()
            if not np.isfinite(x).all() or (x < 0).any() or (x != np.floor(x)).any():
                raise ValueError('Invalid final integer counts')
            total_nnz += int(np.count_nonzero(x))
            depth = x.sum(axis=1, dtype=np.int64)
            if (depth < 1).any() or (depth > 1000000).any():
                raise ValueError('Invalid final library depths')
            if (sparse.csr_matrix(scored.X[rows]) != sparse.csr_matrix(x)).nnz:
                raise ValueError('Scored treated counts differ from final prediction')
            audit = audit_by_target[target]
            if audit['seed'] != seed or audit['counts_array_sha256'] != array_sha256(x):
                raise ValueError('Final count-array hash or seed differs from generation receipt')
            if emitter == 'preserving':
                selected = verify_donor_metadata(audit, target, seed, controls, depth)
                if (audit['edited_molecules'] != int(np.abs(x.astype(np.int64)-selected).sum()//2)
                        or audit['activated_donor_zero_entries'] != int(np.count_nonzero((selected == 0) & (x > 0)))):
                    raise ValueError('Donor edit diagnostics differ from final counts')
                fit_keys = ('fitted_mean_l1', 'fitted_bulk_l1')
            else:
                if not np.array_equal(depth, expected_depths[seed]):
                    raise ValueError('Current-emitter pooled template depths changed')
                fit_keys = ('fitting_mean_l1', 'fitting_bulk_l1')
            if audit['projection_l1'] > .03 or max(audit[k] for k in fit_keys) > .001:
                raise ValueError('Frozen generation feasibility limit exceeded')
            if recipe == 'current_current' and seed == SEEDS[0] and not np.array_equal(x, old):
                raise ValueError('Reused original incumbent differs from regenerated counts')
            p = table[(table.context == 'H1') & (table.target == target)].set_index('gene').reindex(genes)
            required = ['mean_probability_requested', 'bulk_probability_requested', 'mean_probability_sample_expected',
                        'bulk_probability_sample_expected', 'mean_probability_realized', 'bulk_probability_realized',
                        'gene_count_sums_realized']
            if not np.isfinite(p[required].to_numpy()).all():
                raise ValueError('Nonfinite/incomplete diagnostic profile')
            for kind in ('mean', 'bulk'):
                np.testing.assert_array_equal(p[kind+'_probability_requested'].to_numpy(), profiles[response+'_'+kind][i])
            mean = (x / depth[:, None]).mean(axis=0)
            sums = x.sum(axis=0, dtype=np.int64); bulk = sums / sums.sum()
            before_depth = old.sum(axis=1, dtype=np.int64)
            before_mean = (old / before_depth[:, None]).mean(axis=0)
            before_bulk = old.sum(axis=0, dtype=np.float64) / before_depth.sum()
            np.testing.assert_array_equal(p.gene_count_sums_realized.to_numpy(), sums)
            np.testing.assert_allclose(p.mean_probability_realized, mean, atol=1e-14, rtol=1e-12)
            np.testing.assert_allclose(p.bulk_probability_realized, bulk, atol=1e-14, rtol=1e-12)
            errors = {}
            for stage, profile in [('requested', p.bulk_probability_requested),
                                   ('sample_expected', p.bulk_probability_sample_expected), ('realized', bulk)]:
                errors[stage] = mean_errors(profile, reference[target], reference['non-targeting'], genes, target)
            if emitter == 'preserving':
                if not np.isfinite(p[['mean_probability_population_fitted', 'bulk_probability_population_fitted']].to_numpy()).all():
                    raise ValueError('Preserving population fitted moments must be finite')
                errors['population_fitted'] = mean_errors(p.bulk_probability_population_fitted,
                    reference[target], reference['non-targeting'], genes, target)
            drift = 1e6 * (mean-before_mean)
            groups.append({'target': target, 'seed': seed, 'profiles': errors,
                'own_gene_excluded': bool(np.any(genes == target)), 'max_realized_cpm_drift_to_same_seed_baseline': float(np.abs(drift).max()),
                'requested_mean_l1': float(np.abs(mean-profiles[response+'_mean'][i]).sum()),
                'requested_bulk_l1': float(np.abs(bulk-profiles[response+'_bulk'][i]).sum()),
                'min_depth': int(depth.min()), 'max_depth': int(depth.max()), 'mean_depth': float(depth.mean()),
                'mean_detected_genes': float((x > 0).sum(axis=1).mean()),
                'projection_l1': float(audit['projection_l1']), 'donor_rows_verified': emitter == 'preserving',
                'counts_bit_identical_to_same_seed_baseline': bool(np.array_equal(x, old))})
            drift_tables.append(pd.DataFrame({'target': target, 'gene': genes, 'own_target': genes == target,
                'mean_cpm_realized': 1e6*mean, 'mean_cpm_requested': 1e6*profiles[response+'_mean'][i],
                'mean_cpm_drift_to_request': 1e6*(mean-profiles[response+'_mean'][i]),
                'mean_cpm_drift_to_same_seed_baseline': drift,
                'bulk_log_drift_to_same_seed_baseline': logbulk(bulk)-logbulk(before_bulk),
                'bulk_probability_drift_to_request': bulk-profiles[response+'_bulk'][i]}))
    finally:
        predicted.file.close(); before.file.close(); scored.file.close()
    if total_nnz != receipt['stored_entries'] or total_nnz != validation['nnz'] or total_nnz > 4750000000:
        raise ValueError('Final matrix density differs from validated receipt/cap')
    drift_path = ROOT / 'reports/redesign_profile_drift' / (name+'.parquet')
    drift_path.parent.mkdir(parents=True, exist_ok=True)
    pd.concat(drift_tables, ignore_index=True).to_parquet(drift_path, index=False)
    coverage = {'same_target_coverage': receipt.get('same_target_coverage'),
                'uncovered_targets': receipt.get('uncovered_targets', []), 'source_families': [s['path'] for s in receipt['sources']]}
    if response == 'adapter':
        coverage['gene_target_coverage'] = float((profiles['adapter_coverage_log2fc'] > 0).mean())
    profile_errors = {stage: {mask: float(np.mean([g['profiles'][stage][mask]['uncorrected_mean_error'] for g in groups]))
                             for mask in ('own_excluded', 'whole_axis')} for stage in groups[0]['profiles']}
    return {'passed': True, 'prediction_sha256': actual_hash, 'receipt_sha256': sha256(path.with_suffix('.receipt.json')),
            'groups': groups, 'profile_mean_errors': profile_errors, 'coverage': coverage,
            'realized_profiles_recomputed_from_final_counts': True, 'requested_profiles_match_frozen_npz_exactly': True,
            'scored_treated_and_control_counts_verified': True, 'per_gene_drift': str(drift_path.relative_to(ROOT)),
            'per_gene_drift_sha256': sha256(drift_path),
            'max_cpm_drift_to_same_seed_baseline': max(g['max_realized_cpm_drift_to_same_seed_baseline'] for g in groups)}


def render(result):
    lines = ['# Response × emitter development experiment', '',
        'Eight frozen arms; 50 H1 Training targets, 400 real/predicted cells per target and 4,000 controls. '
        'H1 responses were excluded from adapter labels and inputs. These are raw metrics on a development cohort.', '',
        '| Arm | Seed | PDS ↑ | MSE ↓ | NMAE ↓ | Fidelity ↑ | Reach ↑ | Jaccard ↑ | Own-excluded mean error ↓ |',
        '|---|---:|---:|---:|---:|---:|---:|---:|---:|']
    for recipe in RECIPES:
        for seed in SEEDS:
            data = result['arms'][arm_name(recipe, seed)]
            numbers = [f'{data["raw_metrics"][m]:.9f}' for m in METRICS]
            u = data['matrix_audit']['profile_mean_errors']['realized']['own_excluded']
            lines.append('| '+' | '.join([recipe, str(seed), *numbers, f'{u:.9f}'])+' |')
    for recipe, comp in result['recipes'].items():
        lines.extend(['', f'## {recipe}: '+('PASS BOTH SEEDS' if comp['both_seeds_pass'] else 'FAIL STRICT'), ''])
        for seed, detail in comp['seeds'].items():
            lines.append(f'- Seed {seed}: {"PASS" if detail["gate"]["passed"] else "FAIL"}; '
                         f'all-six nonregression={detail["gate"]["all_six_nonregression"]}; '
                         f'fidelity/Jaccard improvement >1e-8={detail["gate"]["fidelity_or_jaccard_gain_gt_tolerance"]}.')
            for item in detail['metrics'].values():
                interval = item['paired_bootstrap']
                if interval:
                    low, high = interval['percentile_95_ci']
                    lines.append(f'  - {item["label"]}: gain {item["gain_in_better_direction"]:+.9f}, paired 95% [{low:+.9f}, {high:+.9f}].')
                else:
                    lines.append(f'  - {item["label"]}: gain {item["gain_in_better_direction"]:+.9f}; {item["uncertainty_note"]}')
            ci = detail['uncorrected_error_paired_bootstrap']
            lines.append(f'  - Own-excluded mean-error gain {ci["mean_gain"]:+.9f}; paired 95% {ci["percentile_95_ci"]}.')
    lines += ['', result['second_context_gate']['action'], '',
              'The point gate requires all six raw metrics nonworse within 1e-8, improvement greater than 1e-8 '
              'in fidelity or Jaccard, identical finite cohorts, and a pass at both seeds. Intervals are descriptive '
              'and cannot turn a failed point gate into a pass. No global normalized-MSE target bootstrap or cross-panel anchors are used.', '',
              'Adapter response includes learned transfer, exact pseudocount inversion and omission of CD4. '
              'The emitter changes donor construction and count transport. Requested profiles are held fixed across emitters; '
              'realized arithmetic/bulk drift is measured. The JSON retains donor/axis/hash checks, coverage, mean-error '
              'before and after own-target exclusion, correction decomposition, DE calls and per-gene drift artifacts. '
              'See [receipt metadata notes](redesign_receipt_metadata_notes.json) for inherited legacy source weights/common subtraction versus actual adapter semantics.', '']
    lines.extend('- '+s for s in result['limitations'])
    return '\n'.join(lines)+'\n'


def main():
    terminals = require_all_terminals(aggregate_path=ROOT / 'reports/redesign_terminal_statuses.json')
    # No metric/performance reads precede this all-terminal barrier.
    freeze_path = ROOT / 'reports/redesign_implementation_freeze.json'; frozen = read_json(freeze_path)
    if tuple(frozen['seeds']) != SEEDS:
        raise ValueError('Registered emission seeds changed')
    for name, digest in frozen['files'].items():
        if sha256(ROOT/name) != digest:
            raise ValueError('Frozen artifact changed: '+name)
    freeze_hash = sha256(freeze_path)
    old_cache = read_json(ROOT / 'reports/dispersion_real_de_cache.json')
    manifest = read_json(REFERENCE/'manifest.json'); cohort = read_json(REFERENCE/'cohort_freeze.json')
    targets = cohort['development_targets']
    if len(targets) != 50 or set(targets) & set(cohort['reserved_targets_not_extracted']):
        raise ValueError('Development/reserve partition changed')
    if (sha256(REFERENCE/'reference.h5ad') != manifest['reference_sha256']
            or manifest['reference_sha256'] != old_cache['reference_sha256']
            or sha256(REFERENCE/'context_H1.h5ad') != manifest['controls_sha256']
            or manifest['controls_sha256'] != old_cache['controls_sha256']):
        raise ValueError('Frozen reference/control bytes changed')
    identity = read_json(ROOT/'reports/redesign_baseline_identity.json')
    if not identity['counts_bit_identical'] or not identity['axes_identical']:
        raise ValueError('Full regenerated baseline identity required')
    old_prediction = ROOT/'outputs/dispersion_h1_400_incumbent/prediction.h5ad'
    new_baseline = OUTPUTS/arm_name('current_current', SEEDS[0])/'prediction.h5ad'
    if (identity['old_incumbent_sha256'] != sha256(old_prediction)
            or identity['new_incumbent_sha256'] != sha256(new_baseline)
            or identity['new_receipt_sha256'] != sha256(new_baseline.with_suffix('.receipt.json'))):
        raise ValueError('Full baseline identity provenance changed')
    profile_path = ROOT/'models/response_adapter_v1/profiles.npz'
    profile_receipt = read_json(profile_path.with_suffix('.json'))
    if profile_receipt['profiles_sha256'] != sha256(profile_path) or not profile_receipt['h1_responses_excluded']:
        raise ValueError('Requested profile identity/H1 exclusion changed')
    for name, digest in profile_receipt['input_hashes'].items():
        if sha256(ROOT/name) != digest:
            raise ValueError('Response-profile input changed: '+name)
    with np.load(profile_path, allow_pickle=False) as source:
        profiles = {k: source[k] for k in source.files}
    genes, reference = reference_profiles(targets)
    if not np.array_equal(profiles['genes'], genes) or profiles['targets'].tolist() != targets:
        raise ValueError('Frozen response profile axes differ from evaluation')
    context = ad.read_h5ad(REFERENCE/'context_H1.h5ad')
    if context.shape != (4000, len(genes)) or not np.array_equal(context.var_names, genes):
        raise ValueError('Incorrect measured control population')
    controls = sparse.csr_matrix(context.X)
    expected_depths = {}
    for seed in SEEDS:
        _, depth, _, _ = control_template(REFERENCE/'context_H1.h5ad', genes, 400, 4, seed)
        expected_depths[seed] = depth
    baseline = metric_arm(REFERENCE/'incumbent_cpu')
    arms, metric_data, table_rows = {}, {}, []
    for recipe in RECIPES:
        for seed in SEEDS:
            name = arm_name(recipe, seed); directory = metric_directory(recipe, seed)
            data = metric_arm(directory); metric_data[name] = data
            before_path = old_prediction if recipe == 'current_current' and seed == SEEDS[0] else OUTPUTS/arm_name('current_current', seed)/'prediction.h5ad'
            audit = audit_counts(recipe, seed, targets, genes, reference, controls, expected_depths,
                                 profiles, frozen, before_path, directory)
            provenance = verify_metric_identity(data, baseline, name, terminals[name], freeze_hash, old_cache, audit['prediction_sha256'])
            # Checks every per-target aggregate even for the two baselines.
            compare_metrics(data, data, targets)
            arms[name] = {'recipe': recipe, 'seed': seed, 'raw_metrics': data['raw'], 'matrix_audit': audit,
                          'provenance': provenance, 'de_call_diagnostics': reference_call_diagnostics(directory, targets)}
            for metric, (label, direction) in METRICS.items():
                aggregation = data['aggregation'].loc[metric]
                table_rows.append({'arm': name, 'response': recipe.split('_')[0], 'emitter': recipe.split('_')[1],
                    'seed': seed, 'metric': metric, 'label': label, 'value': data['raw'][metric],
                    'better_direction': direction, 'n_used': None if bool(aggregation['derived']) else int(aggregation['n_used']),
                    'uncorrected_mean_error_own_excluded': audit['profile_mean_errors']['realized']['own_excluded'],
                    'uncorrected_mean_error_whole_axis': audit['profile_mean_errors']['realized']['whole_axis']})
    recipes = {}
    for recipe in RECIPES[1:]:
        seed_results, gates = {}, {}
        for seed in SEEDS:
            name, base_name = arm_name(recipe, seed), arm_name('current_current', seed)
            data, base = metric_data[name], metric_data[base_name]
            metrics = compare_metrics(base, data, targets); gates[seed] = point_gate(metrics)
            before = np.asarray([g['profiles']['realized']['own_excluded']['uncorrected_mean_error'] for g in arms[base_name]['matrix_audit']['groups']])
            after = np.asarray([g['profiles']['realized']['own_excluded']['uncorrected_mean_error'] for g in arms[name]['matrix_audit']['groups']])
            delta_u = float(after.mean()-before.mean())
            delta_num = float(data['aggregate'].loc['mean','expr_mse_unbiased_capped']-base['aggregate'].loc['mean','expr_mse_unbiased_capped'])
            seed_results[str(seed)] = {'paired_baseline': base_name, 'gate': gates[seed], 'metrics': metrics,
                'uncorrected_error_paired_bootstrap': paired_bootstrap(before-after),
                'mse_decomposition': {'delta_uncorrected_error': delta_u, 'delta_capped_numerator': delta_num,
                    'delta_effective_prediction_correction': delta_u-delta_num,
                    'fixed_reference_denominator_mean': float(base['aggregate'].loc['mean','expr_distance_unbiased']),
                    'interpretation': 'Candidate minus same-seed baseline; mean response and capped global prediction correction can both change.'}}
        recipes[recipe] = {'both_seeds_pass': recipe_gate(gates), 'seeds': seed_results}
    eligible = [r for r in RECIPES[1:] if recipes[r]['both_seeds_pass']]
    result = {'implementation_freeze': frozen, 'analysis_sha256': sha256(__file__),
        'analysis_helpers_sha256': {p: sha256(ROOT/p) for p in ['scripts/summarize_bulk_ablation.py','scripts/summarize_dispersion_experiment.py']},
        'terminal_statuses': terminals, 'baseline_identity': identity, 'cohort_freeze': cohort,
        'profile_receipt': profile_receipt, 'arms': arms, 'recipes': recipes,
        'receipt_metadata_notes': {'path': 'reports/redesign_receipt_metadata_notes.json',
                                   'sha256': sha256(ROOT/'reports/redesign_receipt_metadata_notes.json')},
        'point_gate': 'Each seed: all six favorable gains >=-1e-8, fidelity or Jaccard gain >1e-8, identical finite cohorts. Both seeds must pass.',
        'second_context_gate': {'eligible_recipes': eligible, 'permitted': bool(eligible),
            'action': 'At least one recipe passes both development seeds; a recipe must be frozen before any authorized second-context test. No automatic promotion.' if eligible
                      else 'No recipe passes both seeds: retain incumbent and stop before second-context or reserved-response extraction.'},
        'overall_raw_score_computed': False, 'anchors_used': False, 'promotion_or_submission': False,
        'limitations': ['Previously observed H1 Training development cohort; not untouched destination validation.',
            'Paired target intervals are descriptive with shared controls/correlated perturbations and no multiplicity correction.',
            'Two emission seeds do not estimate all generation uncertainty; neither seed is selected after observing outcomes.',
            'The adapter uses only three non-H1 source contexts and LOSO-selected regularization; its confidence is not DE significance.',
            'Adapter versus current response also changes CD4 inclusion and exact arithmetic inverse; those are part of the registered response intervention.',
            'Adapter prediction receipts inherit legacy raw-source weights 2/1/1; actual adapter pooling is equal over valid raw families before learned feature coefficients, as fixed by adapter code/model/profiles.',
            'Point nonregression is an operational gate, not statistical proof of no harm.',
            'Global capped MSE has a shared correction budget and receives no naive target bootstrap.',
            'Native gene panels and this 18,080-gene reference differ from the challenge; no leaderboard gain is inferred.']}
    destination = ROOT/'reports/redesign_experiment_results.json'
    destination.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    destination.with_suffix('.md').write_text(render(result))
    pd.DataFrame(table_rows).to_csv(ROOT/'reports/redesign_experiment_metrics.csv', index=False)
    print(json.dumps({'both_seed_verdicts': {r:v['both_seeds_pass'] for r,v in recipes.items()},
                      'second_context_gate': result['second_context_gate']}, indent=2), flush=True)


if __name__ == '__main__':
    main()
