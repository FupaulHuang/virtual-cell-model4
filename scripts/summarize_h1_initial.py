#!/usr/bin/env python3
"""Score completed initial H1 arms against one fresh official reference bundle."""
import json
from pathlib import Path
import subprocess
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
EVAL = ROOT / 'evaluations/h1_validation'
METRICS = {
    'pds_cosine': 'PDS', 'expr_mse_unbiased_capped_norm': 'MSE',
    'de_wilcoxon_lfc_nmae': 'NMAE',
    'de_wilcoxon_direction_fidelity_yield_raw': 'Fidelity',
    'de_wilcoxon_direction_reach_raw': 'Reach', 'de_wilcoxon_sig_jaccard': 'Jaccard',
}
ARMS = [
    ('neutral_pool4', 'neutral_cpu', EVAL / 'neutral_prediction.receipt.json'),
    ('pie_dual_pool4', 'pie_dual_cpu', ROOT / 'outputs/h1_pie_dual/prediction.receipt.json'),
    ('pie_single_pool4', 'pie_single_cpu', ROOT / 'outputs/h1_pie_single/prediction.receipt.json'),
    ('pie_dual_pool1', 'pie_dual_pool1_cpu', ROOT / 'outputs/h1_pie_dual_pool1/prediction.receipt.json'),
]


def main():
    bundle = EVAL / 'real_bundle_cpu'
    if not (bundle / 'manifest.json').exists():
        raise SystemExit('Fresh CPU reference anchors are not complete')
    report = {'context': 'public H1 Validation', 'source_exclusion': 'all H1 response banks',
              'backend': 'scanpy', 'device': 'cpu', 'control_source': 'real',
              'reference_manifest': json.loads((bundle / 'manifest.json').read_text()), 'arms': {}}
    for name, folder, receipt_path in ARMS:
        directory = EVAL / folder
        if not (directory / 'agg_results.csv').exists():
            raise SystemExit(f'Arm {name} has no completed aggregate metrics')
        score_path = directory / 'official_local_scores.csv'
        subprocess.run([str(ROOT / '.venv/bin/cell-eval2'),
                        'score', '--user-agg', str(directory / 'agg_results.csv'),
                        '--real-bundle', str(bundle), '-o', str(score_path)], check=True)
        scores = pd.read_csv(score_path).set_index('metric')
        components = {metric: float(scores.loc[metric, 'from_replicate']) for metric in METRICS}
        overall = float(scores.loc['avg_score', 'from_replicate'])
        if abs(overall - sum(components.values()) / 6) > 1e-12:
            raise ValueError('Official bundle avg_score disagrees with the six-component mean')
        aggregate = pd.read_csv(directory / 'agg_results.csv').set_index('statistic')
        receipt = json.loads(receipt_path.read_text())
        covered = set(receipt['targets']) - set(receipt.get('uncovered_targets', receipt['targets']))
        per_target = pd.read_csv(directory / 'results.csv')
        strata = {}
        for label, selected in [('covered', covered), ('uncovered', set(receipt['targets']) - covered)]:
            subset = per_target[per_target.perturbation.isin(selected)]
            strata[label] = {'targets': len(selected), 'raw_metric_means': {metric: {
                'mean': float(part.value.mean()), 'finite_targets': int(part.value.notna().sum())}
                for metric in METRICS if len(part := subset[subset.metric == metric])}}
        report['arms'][name] = {'overall': overall, 'six_component_mean_verified': True,
            'scaled': components,
            'raw': {metric: float(aggregate.loc['mean', metric]) for metric in METRICS},
            'source_coverage': len(covered), 'coverage_strata': strata,
            'receipt': str(receipt_path.relative_to(ROOT)),
            'mean_final_cpm_probability_l1': sum(m['mean_probability_l1'] for m in receipt['moments']) / len(receipt['moments']),
            'mean_final_bulk_probability_l1': sum(m['bulk_probability_l1'] for m in receipt['moments']) / len(receipt['moments'])}
    report['limits'] = [
        'Public H1 reference is not the hidden A/B/C panel; local scores do not estimate leaderboard scores.',
        'PIE sources cover only 22 of 50 H1 validation targets.',
        'Pool1 versus pool4 changes donor selection and depths as well as variance.',
        'Single-profile mode replaces the requested arithmetic mean with the bulk-derived profile, so its comparison also changes effective arithmetic-mean effect amplitudes.',
        'Initial runs did not request DE-table export; significant-set call counts are unavailable.',
        'No H1 Test result was used in this comparison.',
        'The expression scorer flagged suppressed across-perturbation pseudobulk scatter from moment-constrained counts; its bounded correction remains active.',
        'The initial evaluator processes completed and wrote all metric artifacts, but their shell wrappers subsequently exited 127 because the live wrapper file was edited. Final scoring validates persisted run metadata; metrics were not restarted.',
    ]
    output = ROOT / 'reports/h1_initial_comparison.json'
    output.write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    lines = ['# Initial H1 held-out-context comparison', '',
             'All arms exclude H1 response banks and use identical CPU/scanpy metric definitions, '
             'real controls and fresh five-split reference anchors. Scores are local to this public H1 bundle.', '',
             '| Arm | Covered | Overall | ' + ' | '.join(METRICS.values()) + ' |',
             '| --- | ---: | ---: | ' + ' | '.join(['---:'] * len(METRICS)) + ' |']
    for name, arm in report['arms'].items():
        lines.append(f'| {name} | {arm["source_coverage"]}/50 | {arm["overall"]:.6f} | ' +
                     ' | '.join(f'{arm["scaled"][m]:.6f}' for m in METRICS) + ' |')
    lines += ['', 'The table reports reference-scaled components. Full raw values, coverage-stratified '
              'raw summaries and post-rounding moment errors are in `h1_initial_comparison.json`.', '']
    lines += ['- ' + limit for limit in report['limits']]
    output.with_suffix('.md').write_text('\n'.join(lines) + '\n')
    print(json.dumps({k: v['overall'] for k, v in report['arms'].items()}, indent=2))


if __name__ == '__main__':
    main()
