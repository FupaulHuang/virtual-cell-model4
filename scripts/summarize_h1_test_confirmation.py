#!/usr/bin/env python3
"""Summarize the frozen H1 Test confirmation without fitting or anchor scaling."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
EVAL = ROOT / 'evaluations/h1_test'
METRICS = {
    'pds_cosine': ('PDS', 1),
    'expr_mse_unbiased_capped_norm': ('Expression MSE', -1),
    'de_wilcoxon_lfc_nmae': ('LFC NMAE', -1),
    'de_wilcoxon_direction_fidelity_yield_raw': ('Direction fidelity', 1),
    'de_wilcoxon_direction_reach_raw': ('Direction reach', 1),
    'de_wilcoxon_sig_jaccard': ('Significance Jaccard', 1),
}


def de_call_diagnostics(directory):
    result = {}
    for side in ['real', 'pred']:
        data = pd.read_parquet(directory / f'de_{side}.parquet')
        calls = data[(data.p_adj <= .05) & np.isfinite(data.log2_fold_change)
                     & (data.target != data.feature)]
        counts = calls.groupby('target').size()
        result[side] = {'significant_target_gene_pairs': len(calls),
                        'targets_with_calls': len(counts),
                        'mean_calls_per_target_over_all_100': float(len(calls) / 100),
                        'note': 'Descriptive p_adj<=0.05 calls, target gene excluded; not an alternative official metric.'}
    return result


def main():
    launch = json.loads((ROOT / 'reports/h1_test_confirmation_launch.json').read_text())
    terminals = json.loads((ROOT / 'reports/h1_test_terminal_statuses.json').read_text())
    if any(terminals[arm]['exit_code'] != 0 for arm in ['selected', 'neutral']):
        raise ValueError('Both metric jobs must have finished successfully')
    manifest = launch['reference_manifest']
    result = {'protocol': 'Frozen reserved H1 Test confirmation, one selected model versus neutral controls',
              'retuned_on_test': False, 'test_reference_scaled_score_computed': False,
              'validation_anchors_used': False, 'launch_provenance': launch,
              'frozen_model_snapshot': 'reports/h1_test_frozen_model.json',
              'terminal_statuses': terminals, 'arms': {}, 'comparisons': {}}
    per_target = {}
    reference_identity = None
    for arm in ['selected', 'neutral']:
        directory = EVAL / f'{arm}_cpu'
        meta = json.loads((directory / 'run_meta.json').read_text())
        cfg = yaml.safe_load((directory / 'run_params.yaml').read_text())
        if meta['resolved_device'] != 'cpu' or meta['resolved_de_backend'] != 'scanpy' or cfg['control_source'] != 'real':
            raise ValueError('Confirmation backend/control-source mismatch')
        if Path(meta['source']).resolve() != (EVAL / 'reference.h5ad').resolve():
            raise ValueError('Metric run used the wrong reference')
        if reference_identity is not None and reference_identity != meta['source_fingerprint']:
            raise ValueError('Metric runs used different reference data')
        reference_identity = meta['source_fingerprint']
        agg = pd.read_csv(directory / 'agg_results.csv').set_index('statistic')
        aggregation = pd.read_csv(directory / 'metric_aggregation.csv').set_index('metric')
        raw = {metric: float(agg.loc['mean', metric]) for metric in METRICS}
        if not all(np.isfinite(list(raw.values()))):
            raise ValueError('Missing/nonfinite scored raw metric')
        # Wide count includes NaN slots; n_used is the actual finite cohort.
        cohort = {metric: int(aggregation.loc[metric, 'n_used'])
                  if not bool(aggregation.loc[metric, 'derived']) else None for metric in METRICS}
        per_target[arm] = pd.read_csv(directory / 'results.csv')
        receipt = json.loads((ROOT / f'outputs/h1_test_{arm}/prediction.receipt.json').read_text())
        result['arms'][arm] = {'raw_metrics': raw, 'metric_target_cohorts': cohort,
            'source_target_coverage': receipt['same_target_coverage'],
            'uncovered_targets': receipt['uncovered_targets'],
            'mean_final_cpm_probability_l1': float(np.mean([x['mean_probability_l1'] for x in receipt['moments']])),
            'mean_final_bulk_probability_l1': float(np.mean([x['bulk_probability_l1'] for x in receipt['moments']])),
            'de_calls': de_call_diagnostics(directory), 'run_metadata': meta}
    for metric, (label, orientation) in METRICS.items():
        selected = result['arms']['selected']['raw_metrics'][metric]
        neutral = result['arms']['neutral']['raw_metrics'][metric]
        gain = orientation * (selected - neutral)
        entry = {'label': label, 'better_direction': 'higher' if orientation == 1 else 'lower',
                 'selected': selected, 'neutral': neutral, 'gain_in_better_direction': gain,
                 'outcome': 'improved' if gain > 1e-12 else 'worsened' if gain < -1e-12 else 'tied'}
        left = per_target['selected'].query('metric == @metric')[['perturbation', 'value']]
        right = per_target['neutral'].query('metric == @metric')[['perturbation', 'value']]
        pairs = left.merge(right, on='perturbation', suffixes=('_selected', '_neutral')).dropna()
        if len(pairs):
            difference = orientation * (pairs.value_selected.to_numpy() - pairs.value_neutral.to_numpy())
            difference = difference[np.isfinite(difference)]
            entry['paired_targets'] = {'n': len(difference), 'improved': int((difference > 1e-12).sum()),
                'worsened': int((difference < -1e-12).sum()), 'tied': int((abs(difference) <= 1e-12).sum())}
        else:
            entry['paired_targets'] = None
            entry['cohort_note'] = 'Expression normalized MSE is a global ratio of sums, not a mean of per-target normalized scores.'
        result['comparisons'][metric] = entry
    result['reference_identity'] = reference_identity
    result['improved_metric_count'] = sum(x['outcome'] == 'improved' for x in result['comparisons'].values())
    result['reference_cohort'] = {'targets': len(manifest['targets']), 'genes': manifest['n_genes'],
        'treated_cells': manifest['n_reference_treated'], 'controls': manifest['n_controls'],
        'min_cells_per_target': min(manifest['reference_cells_per_target'].values()),
        'max_cells_per_target': max(manifest['reference_cells_per_target'].values()),
        'targets_below_200_reference_cells': sum(x < 200 for x in manifest['reference_cells_per_target'].values())}
    result['limitations'] = [
        'This is one public H1 context and is not an estimate of hidden A/B/C leaderboard performance.',
        'H1 Validation was used for model selection; this Test confirms reserved target outcomes within H1, not a new cell context never used for selection.',
        'Both candidates were evaluated once after the selected configuration was frozen; no Test-driven parameter changes are made.',
        'Only raw metric values are compared on this Test reference; no Validation anchors or cross-panel scaled comparisons are used.',
        'Reference group sizes vary; predicted groups contain 200 cells each. NMAE may omit targets below its real-side significance gate.',
        'Usable cohorts come from metric_aggregation.csv n_used: direction reach uses 97 targets and NMAE uses 72; the wide aggregate count column includes the three undefined reach slots.',
        'Local confirmation emits 200 cells per target, while the challenge requires 400; DE detection power therefore differs from production.',
        'Coverage100/100 means some source evidence exists for every target, not every target-gene pair is measured in every source.',
        'Pooling and moment fitting alter variance; the official scorer retains its safeguards against suppressed sampling scatter.',
    ]
    real_calls = result['arms']['selected']['de_calls']['real']['significant_target_gene_pairs']
    predicted_calls = result['arms']['selected']['de_calls']['pred']['significant_target_gene_pairs']
    result['limitations'].append(f'The selected model calls {predicted_calls:,} significant target-gene pairs versus {real_calls:,} in the reference; improved metrics do not eliminate excess DE calls.')
    output = ROOT / 'reports/h1_test_confirmation.json'
    output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    lines = ['# Frozen H1 Test confirmation', '',
        f'The frozen selected model improves {result["improved_metric_count"]}/6 raw metrics versus the matched neutral comparator. '
        'Both scorer jobs completed successfully using CPU/scanpy and the same reserved Test reference. '
        'The model was frozen before these outcomes, and no retuning was performed.', '',
        '| Metric | Neutral | Selected | Gain in better direction | Direction |',
        '| --- | ---: | ---: | ---: | --- |']
    for item in result['comparisons'].values():
        lines.append(f'| {item["label"]} | {item["neutral"]:.6f} | {item["selected"]:.6f} | '
                     f'{item["gain_in_better_direction"]:+.6f} | {item["better_direction"]} |')
    lines += ['', f'The reference contains {manifest["n_reference_treated"]:,} treated cells across 100 targets, '
              f'{manifest["n_controls"]:,} controls and {manifest["n_genes"]:,} genes. '
              'Both predictions contain 200 cells per target. The selected model has source coverage for 100/100 targets; '
              'all H1 response sources are excluded.', '',
              'The JSON report records metric-specific target cohorts, paired target gains, DE call counts, '
              'output hashes, frozen-model provenance and terminal process statuses.', '']
    lines += ['- ' + text for text in result['limitations']]
    output.with_suffix('.md').write_text('\n'.join(lines) + '\n')
    print(json.dumps({'improved_metric_count': result['improved_metric_count'],
                      'comparisons': result['comparisons']}, indent=2), flush=True)


if __name__ == '__main__':
    main()
