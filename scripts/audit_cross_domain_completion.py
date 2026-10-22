#!/usr/bin/env python3
"""Bounded completion/referee audit; never generates, fits, tests, or rescores.

Held-out recovery values are read only after the immutable confirmation report
exists and the recorded staged summaries have both completed successfully.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import hashlib
import json
import time

import numpy as np
import pandas as pd
import anndata as ad

ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT/'reports/cross_domain_completion_audit.json'
FAMILIES = ('composition', 'dirichlet_multinomial')
SEEDS = (20261026, 20261027)
CONTEXTS = ('H1', 'A', 'B', 'C')
HASHES = {}
CHECKS = {}


def path(value):
    value = Path(value)
    return value if value.is_absolute() else ROOT/value


def read(value):
    return json.loads(path(value).read_text())


def sha(value):
    file = path(value).resolve()
    if file not in HASHES:
        digest = hashlib.sha256()
        with file.open('rb') as handle:
            for chunk in iter(lambda: handle.read(8*1024**2), b''):
                digest.update(chunk)
        HASHES[file] = digest.hexdigest()
    return HASHES[file]


def require(name, condition):
    CHECKS[name] = bool(condition)
    if not condition:
        raise ValueError(name)


def same_hash(file, expected, label=None):
    require(label or f'hash:{file}', sha(file) == expected)


def check_freeze(file):
    record = read(file)
    for target, digest in record['files'].items():
        same_hash(target, digest, f'{file}:{target}')
    return record


def write_report(result, lines):
    REPORT.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    REPORT.with_suffix('.md').write_text('\n'.join(lines)+'\n')


def readiness(publish=True):
    required = ['reports/injection_recovery_de_execution.json',
                'reports/injection_recovery_summary_execution.json']
    for stage in ('development', 'confirmation'):
        required += [f'reports/injection_recovery_{stage}{suffix}' for suffix in ('.json', '.md', '.csv')]
    missing = [p for p in required if not path(p).is_file()]
    if not missing:
        execution = read(required[0])
        jobs = execution.get('terminal_jobs', [])
        if len(jobs) != 12 or any(j.get('exit_code') != 0 for j in jobs):
            missing.append('All 12 authoritative DE terminal exit codes must be zero')
        summaries = read(required[1])
        stages = summaries.get('stages', summaries)
        if isinstance(stages, list):
            stages = {x['stage']: x for x in stages}
        if any(stages.get(s, {}).get('exit_code') != 0 for s in ('development', 'confirmation')):
            missing.append('Both staged-summary terminal exit codes must be zero')
    if missing:
        result = {'checked_at_utc': datetime.now(timezone.utc).isoformat(), 'status': 'PENDING',
                  'missing_required_artifacts_or_receipts': missing,
                  'heldout_recovery_values_read': False,
                  'scope': 'Readiness only; final integrity/referee audit has not run.'}
        if publish:
            write_report(result, ['# Cross-domain completion audit: pending', '',
                'Final review waits for both immutable stage reports and authoritative successful terminal receipts. Held-out recovery values have not been read.', '',
                *['- '+x for x in missing]])
            print(json.dumps(result, indent=2), flush=True)
        return False
    return True


def balanced_audit():
    frozen = check_freeze('reports/balanced_sampling_freeze.json')
    result = read('reports/balanced_sampling_results.json')
    require('balanced embedded freeze unchanged', result['frozen'] == frozen)
    state = read('reports/cross_domain_experiment_state.json')
    for phase in ('balanced_feature_exit_codes', 'balanced_null_prepare_exit_codes',
                  'balanced_precision_exit_codes', 'balanced_null_exit_codes'):
        require(phase, set(state[phase]) == set(CONTEXTS) and all(v == 0 for v in state[phase].values()))
    require('balanced summary exit0', state['balanced_summary_exit_code'] == 0)
    numeric = {}
    for context in CONTEXTS:
        folder = ROOT/'evaluations/balanced_sampling_400'/context
        source = frozen['contexts'][context]
        same_hash(source['input']['path'], source['input']['sha256'])
        same_hash(source['precision_summary'], source['precision_summary_sha256'])
        old_summary = read(source['precision_summary'])
        old_folder = path(source['precision_summary']).parent
        for name, digest in old_summary['files'].items():
            same_hash(old_folder/name, digest)
        summary = read(folder/'precision_summary.json')
        require(f'{context} balanced precision exit0', read(folder/'precision_terminal.json')['exit_code'] == 0)
        for name, digest in summary['files'].items():
            same_hash(folder/name, digest)
        feature = read(folder/'features_receipt.json')
        same_hash(folder/'features.npy', feature['features_sha256'])
        same_hash(folder/'feature_state.json', feature['state_sha256'])
        old_draws, new_draws = pd.read_csv(old_folder/'draws.csv'), pd.read_csv(folder/'draws.csv')
        require(f'{context} balanced 200 fixed seeds', len(new_draws) == 200 and not new_draws.seed_id.duplicated().any()
                and set(new_draws.seed_id) == set(frozen['seed_ids']))
        saved = result['contexts'][context]
        for measure in ('arithmetic_log_mse', 'bulk_log_mse'):
            u = old_draws[old_draws.policy == 'uniform'].set_index('seed_id')[measure].reindex(frozen['seed_ids'])
            b = new_draws.set_index('seed_id')[measure].reindex(frozen['seed_ids'])
            require(f'{context} {measure} ratio matches draws', np.isclose(b.mean()/u.mean(), saved['measures'][measure]['uniform']['ratio'], rtol=0, atol=1e-12))
        a, b = (saved['measures'][m]['uniform'] for m in ('arithmetic_log_mse', 'bulk_log_mse'))
        precision = (b['ratio'] <= .95 and b['gain_95_ci'][0] > 0 and a['ratio'] <= 1.02 and a['ratio_95_ci'][1] <= 1.02) if context == 'H1' else all(v['ratio'] <= 1.05 and v['ratio_95_ci'][1] <= 1.05 for v in (a, b))
        require(f'{context} precision gate reconstructed', bool(precision) == saved['precision_gate_passed'])
        null = read(folder/'null/de/result.json')
        require(f'{context} null DE complete', null['complete'])
        same_hash(folder/'null/de/de.parquet', null['de_sha256'])
        receipt = read(folder/'null/prediction.receipt.json')
        same_hash(folder/'null/prediction.h5ad', receipt['output_sha256'])
        require(f'{context} all balanced null count identities', receipt['all_count_identity'])
        values = np.asarray([x['gate_nsig'] for x in null['counts']])
        require(f'{context} null alias completeness', len(values) == 10 and len({x['alias'] for x in null['counts']}) == 10)
        baseline = source['null_tables']['uniform']
        same_hash(baseline['path'], baseline['sha256']); same_hash(baseline['result'], baseline['result_sha256'])
        u = np.asarray([x['gate_nsig'] for x in read(baseline['result'])['counts']])
        um, up, cm, cp = np.median(u), np.quantile(u, .9), np.median(values), np.quantile(values, .9)
        null_pass = um <= 2 and up <= 5 and cm <= max(2, 2*um) and cp <= max(5, 2*up)
        require(f'{context} null gate reconstructed', bool(null_pass) == saved['null_gate_passed'])
        numeric[context] = {'arithmetic_error_ratio': a['ratio'], 'arithmetic_ratio_95_ci': a['ratio_95_ci'],
            'bulk_error_ratio': b['ratio'], 'bulk_ratio_95_ci': b['ratio_95_ci'],
            'precision_pass': bool(precision), 'null_pass': bool(null_pass),
            'balanced_null_max_calls': int(values.max()), 'balanced_null_median': float(cm), 'balanced_null_p90': float(cp)}
    combined = all(x['precision_pass'] and x['null_pass'] for x in numeric.values())
    require('balanced combined gate reconstructed', combined == result['combined_gate_passed'] == result['injection_balanced_arm_permitted'])
    return numeric, combined


def injection_artifacts(combined):
    truth_freeze = check_freeze('evaluations/injection_recovery_400/freeze.json')
    for p, h in [('input_manifest_path', 'input_manifest_sha256'), ('inference_path', 'inference_sha256'),
                 ('old_null_reference_path', 'old_null_reference_sha256')]:
        same_hash(truth_freeze['source_provenance'][p], truth_freeze['source_provenance'][h])
    truth = read('evaluations/injection_recovery_400/manifest.json')
    same_hash('evaluations/injection_recovery_400/freeze.json', truth['freeze_sha256'])
    same_hash(truth['population_parameters_path'], truth['population_parameters_sha256'])
    require('two independent reference families/two seeds', set(truth['families']) == set(FAMILIES) and truth['seeds'] == list(SEEDS))
    amendment = read('reports/injection_recovery_serialization_amendment.json')
    same_hash('evaluations/injection_recovery_400/freeze.json', amendment['amended_freeze_sha256'])
    require('serialization-only amendment preserved scientific parameters', amendment['all_four_generation_attempts_failed_before_any_count_output']
            and amendment['no_refit_or_reselection'] and amendment['no_seed_strength_kappa_or_threshold_change']
            and all(x['gene_axis_exactly_preserved'] and x['all_other_arrays_dtype_shape_bytes_identical'] for x in amendment['array_proofs']))
    truth_execution = read('reports/injection_recovery_truth_execution.json')
    require('all independent truth generation exit0/source audits', truth_execution['all_generation_exit_codes_zero']
            and truth_execution['independent_count_source_depth_audit_passed'] and len(truth_execution['runs']) == 4)
    candidate_freeze = check_freeze('outputs/injection_recovery_400/freeze.json')
    analysis_freeze = check_freeze('reports/injection_recovery_analysis_freeze.json')
    eligibility = read('reports/injection_recovery_eligibility.json')
    expected_methods = ['current', 'preserving_uniform'] + (['preserving_balanced'] if combined else [])
    require('exact prospective method eligibility', eligibility['planned_methods'] == candidate_freeze['methods'] == analysis_freeze['planned_methods'] == expected_methods)
    require('eligibility frozen before injection DE', eligibility['frozen_before_injection_DE'])
    same_hash(eligibility['balanced_gate_path'], eligibility['balanced_gate_sha256'])
    combined_manifest = read('outputs/injection_recovery_400/manifest.json')
    require('all planned candidate methods complete', combined_manifest['all_methods_complete'] and combined_manifest['methods'] == expected_methods)
    same_hash('outputs/injection_recovery_400/freeze.json', combined_manifest['candidate_freeze_sha256'])
    same_hash('outputs/injection_recovery_400/prepared.json', combined_manifest['prepared_manifest_sha256'])
    execution = read('reports/injection_prediction_execution_receipt.json')
    require('candidate preparation/generation/finalize terminal0', all(execution[k] == 0 for k in ('freeze_exit_code','prepare_exit_code','finalize_exit_code'))
            and len(execution['generation_processes']) == 4 and all(x['exit_code'] == 0 for x in execution['generation_processes']))
    require('eight exact candidate null identities', len(execution['null_identity_checks']) == 8 and all(x['exact_raw_identity'] for x in execution['null_identity_checks']))
    same_hash(execution['combined_manifest'], execution['combined_manifest_sha256'])
    candidate_records = {}; expected_jobs = set()
    for family, spec in truth['families'].items():
        same_hash(spec['profiles_path'], spec['profiles_sha256'])
        require(f'{family} eight distinct aliases', len(spec['aliases']) == len(set(spec['aliases'])) == 8)
        for seed in SEEDS:
            source = spec['seeds'][str(seed)]
            same_hash(source['receipt_path'], source['receipt_sha256'])
            for role in ('inference', 'comparison', 'reference'):
                same_hash(source[role+'_path'], source[role+'_sha256'])
            same_hash(source['empirical_profiles_path'], source['empirical_profiles_sha256'])
            same_hash(source['source_rows_path'], source['source_rows_sha256'])
            group = combined_manifest['groups'][f'{family}:{seed}']
            same_hash(group['generation_terminal_path'], group['generation_terminal_sha256'])
            terminal = read(group['generation_terminal_path'])
            require(f'{family}/{seed} all candidate aliases attempted', terminal['exit_code'] == 0 and terminal['all_aliases_attempted'])
            for method in ['reference', *expected_methods]:
                expected_jobs.add((family, seed, method))
            request_hashes = []
            for method, record in group['methods'].items():
                require(f'{family}/{seed}/{method} complete protocol axis', record['complete'] and record['DE_permitted']
                        and record['aliases'] == spec['aliases'] and record['successful_aliases'] == spec['aliases'])
                for key in ('prediction', 'receipt', 'gene_profiles', 'alias_audit', 'selected_rows', 'validation', 'requests'):
                    same_hash(record[key+'_path'], record[key+'_sha256'])
                require(f'{family}/{seed}/{method} original format validation', read(record['validation_path'])['valid'])
                request_hashes.append(record['requests_sha256'])
                candidate_records[(family, seed, method)] = record
            require(f'{family}/{seed} identical requests across emitters', len(set(request_hashes)) == 1)
    de_execution = read('reports/injection_recovery_de_execution.json')
    jobs = de_execution['terminal_jobs']
    require('exact complete DE execution roster', {(x['family'],x['seed'],x['method']) for x in jobs} == expected_jobs
            and len(jobs) == len(expected_jobs) and all(x['exit_code'] == 0 for x in jobs))
    for family, seed, method in sorted(expected_jobs):
        directory = ROOT/f'evaluations/injection_recovery_400/de/{family}/seed_{seed}/{method}'
        status, preflight = read(directory/'status.json'), read(directory/'preflight.json')
        require(f'{family}/{seed}/{method} DE status identity', status['status'] == 'complete'
                and (status['family'],status['seed'],status['method']) == (family,seed,method))
        same_hash(directory/'de.parquet', status['de_sha256'])
        same_hash('reports/injection_recovery_analysis_freeze.json', preflight['analysis_freeze_sha256'])
        same_hash('reports/injection_recovery_eligibility.json', preflight['eligibility_sha256'])
        require(f'{family}/{seed}/{method} pinned DE parameters', preflight['parameters'] == analysis_freeze['de_parameters'])
        source = truth['families'][family]['seeds'][str(seed)]
        expected_prediction = source['reference_sha256'] if method == 'reference' else candidate_records[(family,seed,method)]['prediction_sha256']
        require(f'{family}/{seed}/{method} exact DE inputs', preflight['prediction_sha256'] == expected_prediction
                and preflight['reference_sha256'] == source['reference_sha256'] and preflight['comparison_sha256'] == source['comparison_sha256'])
    return truth, candidate_records, expected_methods


def dm_sign_diagnostics(truth, candidate_records, candidates, tables, banks):
    """Explicitly post hoc sign split; no new threshold or decision gate."""
    rows = []
    family = 'dirichlet_multinomial'
    bank = banks[family]
    for seed in SEEDS:
        source = truth['families'][family]['seeds'][str(seed)]
        print(f'DM direct-count audit: materializing seed {seed}', flush=True)
        inference = ad.read_h5ad(source['inference_path'])
        reference = ad.read_h5ad(source['reference_path'])
        comparison = ad.read_h5ad(source['comparison_path'])
        objects = {}; selections = {}
        try:
            comparison_detection = np.asarray((comparison.X[:] > 0).mean(axis=0)).ravel()
            reference_labels = reference.obs.target_gene.astype(str).to_numpy()
            for row in candidates:
                if row['family'] != family or row['seed'] != seed or row['is_null']:
                    continue
                alias, method = row['alias'], row['method']
                i = list(bank['aliases']).index(alias)
                known = bank['truth_lfc_cpm_eps1e9'][i]
                key = (family, seed, method)
                ref_table = tables[(family,seed,'reference')]
                pred_table = tables[key]
                ref = ref_table[ref_table.target == alias].set_index('feature').reindex(bank['genes'])
                pred = pred_table[pred_table.target == alias].set_index('feature').reindex(bank['genes'])
                gate = pd.Index(ref_table.loc[ref_table.target == alias, 'feature']).get_indexer(bank['genes']) >= 0
                ref_call, pred_call = ref.p_adj.to_numpy() < .05, pred.p_adj.to_numpy() < .05
                ref_lfc, pred_lfc = ref.log2_fold_change.to_numpy(), pred.log2_fold_change.to_numpy()
                actual_counts = raw = reference_detection = candidate_detection = None
                if method == 'preserving_uniform':
                    if method not in objects:
                        manifest = candidate_records[key]
                        objects[method] = ad.read_h5ad(manifest['prediction_path'])
                        selections[method] = pd.read_csv(manifest['selected_rows_path'])
                    obj = objects[method]
                    selected = selections[method]
                    selected = selected[selected.alias == alias].sort_values('output_row')
                    require(f'{family}/{seed}/{alias} paired raw donor mapping', len(selected) == 400
                            and selected.output_row.tolist() == list(range(400))
                            and len(set(selected.inference_row)) == 400)
                    idx = selected.inference_row.to_numpy(); order = np.argsort(idx)
                    sorted_raw = inference.X[idx[order]].toarray(); raw = np.empty_like(sorted_raw); raw[order] = sorted_raw
                    labels = obj.obs.target_gene.astype(str).to_numpy()
                    actual_counts = obj.X[np.flatnonzero(labels == alias)].toarray()
                    require(f'{family}/{seed}/{alias} donor depth conservation', np.array_equal(raw.sum(1,dtype=np.int64),actual_counts.sum(1,dtype=np.int64)))
                    reference_detection = np.asarray((reference.X[np.flatnonzero(reference_labels == alias)] > 0).mean(axis=0)).ravel()
                    candidate_detection = (actual_counts > 0).mean(axis=0)
                for sign, mask in [('up', gate & (known >= .1)), ('down', gate & (known <= -.1))]:
                    ref_correct = mask & ref_call & np.isfinite(ref_lfc) & (ref_lfc*known > 0)
                    pred_correct = mask & pred_call & np.isfinite(pred_lfc) & (pred_lfc*known > 0)
                    count = int(mask.sum())
                    result = {k:row[k] for k in ('family','seed','alias','method','partition','regime')}
                    result.update(sign=sign, oracle_adjudicable_genes=count,
                        reference_correct_sign_calls=int(ref_correct.sum()), candidate_correct_sign_calls=int(pred_correct.sum()),
                        reference_signed_recovery=float(ref_correct.sum()/count) if count else None,
                        candidate_signed_recovery=float(pred_correct.sum()/count) if count else None)
                    if method == 'preserving_uniform':
                        zeros = raw[:,mask] == 0
                        activated = zeros & (actual_counts[:,mask] > 0)
                        lost = (raw[:,mask] > 0) & (actual_counts[:,mask] == 0)
                        group = 'oracle_mean_' + sign
                        saved = row['raw_donor_zero_edit_diagnostics']['groups'][group]
                        require(f'{family}/{seed}/{alias}/{sign} activation independently reproduced',
                            saved['donor_zero_entries'] == int(zeros.sum())
                            and saved['zero_to_positive_entries'] == int(activated.sum())
                            and saved['positive_to_zero_entries'] == int(lost.sum()))
                        detection = {
                            'reference_detection': float(reference_detection[mask].mean()) if count else None,
                            'candidate_detection': float(candidate_detection[mask].mean()) if count else None,
                            'comparison_detection': float(comparison_detection[mask].mean()) if count else None,
                            'reference_detection_change': float((reference_detection[mask]-comparison_detection[mask]).mean()) if count else None,
                            'candidate_detection_change': float((candidate_detection[mask]-comparison_detection[mask]).mean()) if count else None,
                            'candidate_reference_detection_mae': float(np.abs(candidate_detection[mask]-reference_detection[mask]).mean()) if count else None}
                        for name, value in detection.items():
                            recorded = row['detection'][group][name]
                            require(f'{family}/{seed}/{alias}/{sign}/{name} independently reproduced',
                                recorded is None if value is None else np.isclose(value,recorded,rtol=0,atol=1e-12))
                        result.update(donor_zero_entries=int(zeros.sum()), zero_to_positive_entries=int(activated.sum()),
                                      positive_to_zero_entries=int(lost.sum()), detection=detection)
                    rows.append(result)
            print(f'DM direct-count audit: seed {seed} completed', flush=True)
        finally:
            inference.file.close(); reference.file.close(); comparison.file.close()
            for obj in objects.values(): obj.file.close()
    return {'analysis_type':'Explicit descriptive post hoc per-sign explanation; not a new gate.',
            'sign_rule':'Existing |true LFC|>=0.1 within the realized comparison CPM>5 gate; q<0.05 unchanged.',
            'actual_preserving_raw_donor_edits_and_detection_independently_verified':True,
            'rows':rows,
            'causal_limit':'Different up/down recovery and zero observed donor-zero activation can support a distribution/support explanation in this oracle law. They do not isolate the support prior as a causal mechanism: no prior-strength or generator-component ablation was performed, and references are not per-cell counterfactual pairs.'}


def staged_results(truth, methods, candidate_records):
    # Caller only enters after readiness() verifies the immutable confirmation exists.
    development = read('reports/injection_recovery_development.json')
    confirmation = read('reports/injection_recovery_confirmation.json')
    execution = read('reports/injection_recovery_summary_execution.json')
    stages = execution.get('stages', execution)
    if isinstance(stages, list):
        stages = {x['stage']: x for x in stages}
    same_hash('reports/injection_recovery_analysis_freeze.json', execution['analysis_freeze_sha256'])
    for stage, report in [('development', development), ('confirmation', confirmation)]:
        require(f'{stage} summary process exit0', stages[stage]['exit_code'] == 0)
        expected_report_hash = stages[stage].get('report_sha256')
        if expected_report_hash is not None:
            same_hash(f'reports/injection_recovery_{stage}.json', expected_report_hash)
        require(f'{stage} stage identity', report['stage'] == stage and report['planned_methods'] == methods)
        same_hash('reports/injection_recovery_analysis_freeze.json', report['analysis_freeze_sha256'])
        same_hash('reports/injection_recovery_eligibility.json', report['eligibility_sha256'])
    same_hash('reports/injection_recovery_development.json', confirmation['development_report']['sha256'])
    require('development chronologically precedes confirmation', development['created_utc'] <= confirmation['created_utc'])
    require('development opens only template1/seed26', all(x['template'] == 1 and x['seed'] == 20261026 for x in development['reference_results']))
    require('confirmation excludes development partition', all(not (x['template'] == 1 and x['seed'] == 20261026) for x in confirmation['reference_results']))
    references = development['reference_results'] + confirmation['reference_results']
    candidates = development['candidate_results'] + confirmation['candidate_results']
    require('complete distinct staged reference roster', len(references) == 32 and len({(x['family'],x['seed'],x['alias']) for x in references}) == 32)
    require('complete distinct staged candidate roster', len(candidates) == 32*len(methods)
            and len({(x['family'],x['seed'],x['alias'],x['method']) for x in candidates}) == 32*len(methods))
    banks, tables = {}, {}
    for family, spec in truth['families'].items():
        with np.load(spec['profiles_path'], allow_pickle=False) as data:
            banks[family] = {k:data[k] for k in ('genes','aliases','truth_lfc_cpm_eps1e9','marginal_null_gene_mask','is_null')}
    def check_record(record, method):
        key = (record['family'],record['seed'],method)
        if key not in tables:
            tables[key] = pd.read_parquet(ROOT/f'evaluations/injection_recovery_400/de/{key[0]}/seed_{key[1]}/{key[2]}/de.parquet')
        bank = banks[key[0]]; i = list(bank['aliases']).index(record['alias'])
        table = tables[key]; sub = table[table.target == record['alias']].set_index('feature')
        require(f'{key}/{record["alias"]} unique DE rows', not sub.index.duplicated().any())
        aligned = sub.reindex(bank['genes']); gate = sub.index.get_indexer(bank['genes']) >= 0
        known = bank['truth_lfc_cpm_eps1e9'][i]; q = aligned.p_adj.to_numpy(); observed = aligned.log2_fold_change.to_numpy()
        called = gate & (q < .05); adjudicable = gate & (np.abs(known) >= .1)
        correct = called & adjudicable & np.isfinite(observed) & (known*observed > 0)
        marginal = gate & bank['marginal_null_gene_mask'][i]
        actual = {'significant_calls': int(called.sum()), 'oracle_adjudicable_gated': int(adjudicable.sum()),
                  'called_with_correct_oracle_sign': int(correct.sum()), 'marginal_null_calls': int((called&marginal).sum())}
        require(f'{key}/{record["alias"]} diagnostic values independently reproduced', all(record['metrics'][k] == v for k,v in actual.items()))
        strata = record['metrics']['effect_strata']
        require(f'{key}/{record["alias"]} nonoverlapping bins cover tested genes', sum(x['oracle_genes'] for x in strata) == int(gate.sum())
                and sum(x['significant_calls'] for x in strata) == int(called.sum()))
        recovery = None if not adjudicable.any() else float(correct.sum()/adjudicable.sum())
        require(f'{key}/{record["alias"]} signed recovery denominator', recovery == record['metrics']['signed_recovery'])
    for row in references:
        check_record(row, 'reference')
    for row in candidates:
        require('all candidates scored completely', row['status'] == 'complete')
        check_record(row, row['method'])
    nulls = [x for x in references if x['is_null']]
    require('all eight genuine reference NULLs assessed', len(nulls) == 8)
    validity = {}
    for family in FAMILIES:
        for seed in SEEDS:
            part = [x for x in nulls if x['family'] == family and x['seed'] == seed]
            passed = max(x['metrics']['significant_calls'] for x in part) <= 5
            saved = next(x for x in confirmation['reference_null_validity'] if x['family'] == family and x['seed'] == seed)
            require(f'{family}/{seed} complete NULL guard', saved['complete_family_seed_null_panel'] and saved['opened_null_aliases'] == 2
                    and saved['status'] == ('pass' if passed else 'inconclusive_reference_null'))
            validity[f'{family}:{seed}'] = {'passed': passed, 'maximum_calls': max(x['metrics']['significant_calls'] for x in part)}
    compact = []
    for row in candidates:
        if row['regime'] in ('NULL','STRONG'):
            compact.append({k:row[k] for k in ('family','seed','alias','method','partition','regime')} | {
                'reference_calls':row['reference_metrics']['significant_calls'],
                'candidate_calls':row['metrics']['significant_calls'],
                'reference_signed_recovery':row['reference_metrics']['signed_recovery'],
                'candidate_signed_recovery':row['metrics']['signed_recovery'],
                'marginal_null_calls':row['metrics']['marginal_null_calls']})
    return {'development_sha256': sha('reports/injection_recovery_development.json'),
            'confirmation_sha256': sha('reports/injection_recovery_confirmation.json'),
            'reference_rows': len(references), 'candidate_rows': len(candidates),
            'all_reported_calls_sign_matches_bins_and_recovery_verified': True,
            'reference_null_validity': validity, 'null_and_strong_rows': compact,
            'descriptive_posthoc_DM_sign_split': dm_sign_diagnostics(truth, candidate_records, candidates, tables, banks)}


def preservation():
    old = read('reports/completion_audit.json')
    model = read('models/m4_multisource_dual/model.json')
    same_hash('models/m4_multisource_dual/model.json', old['frozen_model_sha256'], 'published incumbent model unchanged')
    same_hash(model['configuration'], old['production_config_sha256'], 'published production configuration unchanged')
    bulk = read('reports/bulk_implementation_freeze.json')
    core = {}
    for name in ('scripts/build_prediction.py','external/atlasshift/model.py','external/atlasshift/predict.py'):
        same_hash(name, bulk['files'][name], f'published generator unchanged:{name}')
        core[name] = sha(name)
    redesign = read('reports/redesign_implementation_freeze.json')
    same_hash('scripts/control_preserving_emitter.py', redesign['files']['scripts/control_preserving_emitter.py'], 'preserving emitter unchanged')
    cohort = read('evaluations/dispersion_h1_training_400/cohort_freeze.json')
    cohort_path = 'evaluations/dispersion_h1_training_400/cohort_freeze.json'
    prior_dispersion = read('reports/dispersion_experiment_results.json')
    same_hash(cohort_path, prior_dispersion['current_reference_manifest']['cohort_freeze_sha256'],
              'original reserve partition unchanged')
    require('25 reserved targets remain metadata-only in partition', len(cohort['reserved_targets_not_extracted']) == 25)
    state = read('reports/cross_domain_experiment_state.json')
    require('recorded no reserve/submission action', state['reserved_responses_accessed'] is False
            and state['official_submission_made'] is False and state['submission_authorization']['used'] == 0)
    reserve_named = sorted(str(p.relative_to(ROOT)) for directory in ('evaluations','outputs')
                           for p in (ROOT/directory).iterdir() if 'reserve' in p.name.lower())
    declared_roots = ['evaluations/balanced_sampling_400','evaluations/injection_recovery_400',
                      'evaluations/injection_recovery_400_initial_string_encoding','outputs/injection_recovery_400']
    return {'published_model_sha256': sha('models/m4_multisource_dual/model.json'),
            'production_configuration_sha256': sha(model['configuration']), 'published_core_sha256': core,
            'preserving_emitter_sha256': sha('scripts/control_preserving_emitter.py'),
            'recorded_official_submission_made': False, 'recorded_reserved_responses_accessed': False,
            'declared_new_count_artifact_roots': declared_roots, 'reserve_named_top_level_artifacts_observed': reserve_named,
            'submission_inventory_current_names_only': sorted(p.name for p in (ROOT/'submissions').iterdir()),
            'claim_scope': 'No new submission/reserve action is recorded in this execution; incumbent/config/core hashes are unchanged and declared outputs are control-only/synthetic. No global pre-work filesystem inventory exists, so this is not a forensic proof of absence from timestamps or names. The published submission itself predates this study.'}


def main(structural_first=False):
    if not structural_first and not readiness():
        return
    print('Audit phase: balanced sampling integrity', flush=True)
    balanced, combined = balanced_audit()
    print(f'Audit phase: balanced complete; {len(HASHES)} files hashed', flush=True)
    print('Audit phase: injection artifact integrity', flush=True)
    truth, candidate_records, methods = injection_artifacts(combined)
    print(f'Audit phase: injection integrity complete; {len(HASHES)} files hashed', flush=True)
    print('Audit phase: preserved production artifacts', flush=True)
    preserved = preservation()
    if structural_first and not readiness(publish=False):
        checkpoint = {'checked_at_utc':datetime.now(timezone.utc).isoformat(),
            'status':'STRUCTURAL_PASS_WAITING_FOR_CONFIRMATION', 'checks_passed':sum(CHECKS.values()),
            'unique_files_hashed':len(HASHES), 'heldout_recovery_values_read':False,
            'pending':'Immutable confirmation report and authoritative successful staged-summary receipts.'}
        write_report(checkpoint, ['# Cross-domain completion audit: structural checks passed', '',
            'Frozen inputs, outputs, execution receipts and preserved production artifacts passed the structural audit. '
            'Final numeric/referee review is waiting for the immutable confirmation report and successful summary receipt. '
            'Held-out recovery values have not been read.'])
        print(json.dumps(checkpoint), flush=True)
        while not readiness(publish=False):
            time.sleep(15)
    print('Audit phase: staged numerical claims and direct counts', flush=True)
    staged = staged_results(truth, methods, candidate_records)
    print('Audit phase: staged numerical claims complete', flush=True)
    correction = read('reports/balanced_sampling_correction_audit.json')
    same_hash('scripts/audit_balanced_sampling_correction.py', correction['analysis_code_sha256'])
    require('correction diagnostic not official score or new gate', not correction['official_logbulk_or_global_capped_correction_computed'] and not correction['new_selection_gate'])
    result = {'checked_at_utc': datetime.now(timezone.utc).isoformat(), 'status': 'PASS',
        'meaning': 'Execution/provenance and bounded numerical-claim audit passed; this is not model promotion or proof of biological performance.',
        'analysis_code_sha256': sha(__file__), 'checks': CHECKS, 'checks_passed': sum(CHECKS.values()),
        'unique_files_hashed': len(HASHES), 'missing_required_artifacts_or_receipts': [],
        'heldout_recovery_values_read_only_after_confirmation_report_written': True,
        'balanced': balanced, 'balanced_combined_gate_passed': combined, 'exact_injection_methods': methods,
        'staged_recovery': staged, 'preservation_and_scope': preserved,
        'sampling_correction_diagnostic_sha256': sha('reports/balanced_sampling_correction_audit.json'),
        'tests_previously_passed_and_not_rerun': {'sampler':10,'truth_after_serialization_fix':14,'candidate':10,'staged_analysis':9,'total':43},
        'numerical_jobs_or_tests_rerun': False,
        'referee_limits': ['Balanced-sampling engineering gates and moment-oracle recovery are separate from biological model selection.',
            'Reference NULL failures make the relevant family/seed inconclusive; all rows remain reported.',
            'Known response moments do not specify complete target cell distributions; poor recovery is not automatically a numerical bug.',
            'No actual H1 perturbation outcome, reserve, official score, or submission is introduced by this audit.']}
    lines = ['# Cross-domain completion audit: PASS', '',
        f'{sum(CHECKS.values())} integrity/consistency checks passed across {len(HASHES)} hashed files. No numerical job or previously passed test was rerun. PASS refers to execution and provenance, not model advancement.', '',
        f'Balanced engineering gate: **{"PASS" if combined else "FAIL"}**. Exactly these injection methods were permitted: '+', '.join(methods)+'.', '',
        'Verified four truth generations, eight complete candidate arms, 12 successful DE jobs and immutable development-then-confirmation reports. Reproduced all32 reference and64 candidate call counts, known-sign matches, effect-bin totals and signed-recovery denominators from the frozen tables. Reference NULL guards were reconstructed per family/seed.', '',
        'Published incumbent, production configuration, original generator and preserving-emitter hashes are unchanged. The25-target reserve partition remains intact. '+preserved['claim_scope'], '',
        'The Unicode serialization-only amendment preserves original archived artifacts and all scientific arrays/parameters. The arithmetic sampling-correction diagnostic is explicitly not the official transformed-bulk/global-cap correction and changes no gate.', '',
        '[Detailed checks, numeric claims and provenance](cross_domain_completion_audit.json). The separate root final report owns the scientific synthesis; this audit does not mark a scaffold report complete.']
    write_report(result, lines)
    print(json.dumps({'status':'PASS','checks_passed':result['checks_passed'],'unique_files_hashed':len(HASHES),
                      'balanced_combined_gate_passed':combined,'methods':methods}, indent=2), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--structural-first', action='store_true', help='Hash finished inputs now, then wait without reading held-out values until confirmation is complete.')
    main(parser.parse_args().structural_first)
