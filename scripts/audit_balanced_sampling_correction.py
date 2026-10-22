#!/usr/bin/env python3
"""Post hoc linear-mean sampling diagnostic; not the official MSE correction."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CONTEXTS = ('H1', 'A', 'B', 'C')
N, SAMPLE, DRAWS = 4000, 400, 200
FPC = 1 - SAMPLE / N


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(8 * 1024**2), b''):
            digest.update(chunk)
    return digest.hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def main():
    frozen_path = ROOT / 'evaluations/balanced_sampling_400/freeze.json'
    frozen = read(frozen_path)
    output = {}
    for context in CONTEXTS:
        old_dir = ROOT / 'evaluations/sampling_precision_400' / context
        new_dir = ROOT / 'evaluations/balanced_sampling_400' / context
        old_summary_path, new_summary_path = old_dir/'summary.json', new_dir/'precision_summary.json'
        if sha(old_summary_path) != frozen['contexts'][context]['precision_summary_sha256']:
            raise ValueError('Frozen uniform/depth comparator summary changed')
        old_summary, new_summary = read(old_summary_path), read(new_summary_path)
        if (read(old_dir/'terminal.json')['exit_code'] != 0
                or read(new_dir/'precision_terminal.json')['exit_code'] != 0
                or old_summary['input'] != frozen['contexts'][context]['input']
                or new_summary['features']['input'] != old_summary['input']):
            raise ValueError('Incomplete or unmatched precision populations')
        for directory, summary in ((old_dir, old_summary), (new_dir, new_summary)):
            for filename in ('draws.csv', 'gene_sampling_audit.parquet'):
                if sha(directory/filename) != summary['files'][filename]:
                    raise ValueError('Sampling diagnostic artifact changed')
        prior = pd.read_parquet(old_dir/'gene_sampling_audit.parquet')
        balanced = pd.read_parquet(new_dir/'gene_sampling_audit.parquet')
        genes = prior.loc[prior.policy == 'uniform', 'gene'].to_numpy()
        if len(genes) != len(set(genes)) or not np.array_equal(balanced.gene, genes):
            raise ValueError('Duplicate or mismatched gene axes')
        tables = {p: prior.loc[prior.policy == p].reset_index(drop=True)
                  for p in ('uniform', 'depth_stratified')}
        for table in tables.values():
            if not np.array_equal(table.gene, genes):
                raise ValueError('Prior policy gene axes differ')
        population = tables['uniform'].arithmetic_population.to_numpy()
        np.testing.assert_array_equal(tables['depth_stratified'].arithmetic_population, population)
        gate = population * 1e6 > 5
        if int(gate.sum()) != old_summary['gated_genes']:
            raise ValueError('Population CPM gate differs from the variance-trace gate')
        uniform_exact = float(tables['uniform'].loc[gate, 'arithmetic_exact_variance'].sum())
        if not np.isfinite(uniform_exact) or uniform_exact <= 0:
            raise ValueError('A positive finite sampling-variance scale is required')
        prior_draws, new_draws = pd.read_csv(old_dir/'draws.csv'), pd.read_csv(new_dir/'draws.csv')
        policies = {}
        for policy in ('uniform', 'depth_stratified', 'balanced'):
            table = balanced if policy == 'balanced' else tables[policy]
            draw = new_draws if policy == 'balanced' else prior_draws[prior_draws.policy == policy]
            if (len(draw) != DRAWS or draw.seed_id.duplicated().any()
                    or set(draw.seed_id) != set(range(20261025, 20261025+DRAWS))):
                raise ValueError('Expected all 200 frozen unique seed draws')
            ratio = draw.variance_trace_ratio.to_numpy(dtype=float)
            mse = table.loc[gate, 'arithmetic_mc_mse'].to_numpy(dtype=float)
            bias = table.loc[gate, 'arithmetic_mc_bias'].to_numpy(dtype=float)
            if not np.isfinite(np.r_[ratio, mse, bias]).all() or (mse < 0).any():
                raise ValueError('Invalid sampling statistics')
            mc_mse = float(mse.sum())
            mc_bias_squared = float((bias*bias).sum())
            if mc_mse <= 0 or mc_bias_squared > mc_mse + 1e-15:
                raise ValueError('Inconsistent Monte Carlo mean decomposition')
            # vu = .9 * population sample variance /400. Each stored draw ratio
            # is its sample-variance trace divided by that same population trace.
            nominal_iid = float(ratio.mean() * uniform_exact / FPC)
            exact = None if policy == 'balanced' else float(table.loc[gate, 'arithmetic_exact_variance'].sum())
            policies[policy] = {
                'draws': DRAWS, 'mean_sample_to_population_variance_trace_ratio': float(ratio.mean()),
                'actual_mc_sum_arithmetic_mean_mse_CPM2': mc_mse * 1e12,
                'nominal_iid_sum_sample_variance_over_n_CPM2': nominal_iid * 1e12,
                'nominal_iid_to_actual_mc_mse_ratio': nominal_iid / mc_mse,
                'actual_mc_mse_to_uniform_exact_variance_ratio': mc_mse / uniform_exact,
                'exact_design_sum_arithmetic_mean_variance_CPM2': None if exact is None else exact * 1e12,
                'actual_mc_mse_to_exact_design_variance_ratio': None if exact is None else mc_mse / exact,
                'empirical_200_draw_mean_bias_squared_CPM2': mc_bias_squared * 1e12,
                'empirical_across_draw_variance_sum_ddof1_CPM2': max(0., mc_mse-mc_bias_squared) * DRAWS/(DRAWS-1) * 1e12,
                'bias_note': 'The finite-200-draw mean deviation is Monte Carlo noise, not proof of systematic design bias.'}
        output[context] = {
            'population_size': N, 'sample_size': SAMPLE, 'gated_genes': int(gate.sum()),
            'gene_axis_size': len(genes), 'gene_axis_sha256': hashlib.sha256('\0'.join(genes).encode()).hexdigest(),
            'gene_gate': 'Full supplied-axis population arithmetic CPM strictly greater than5; identical to stored variance-trace ratios.',
            'uniform_exact_sum_mean_variance_CPM2': uniform_exact * 1e12,
            'policies': policies,
            'inputs': {'population': old_summary['input'],
                'old_summary_sha256': sha(old_summary_path), 'balanced_summary_sha256': sha(new_summary_path),
                'old_gene_audit_sha256': sha(old_dir/'gene_sampling_audit.parquet'),
                'balanced_gene_audit_sha256': sha(new_dir/'gene_sampling_audit.parquet'),
                'old_draws_sha256': sha(old_dir/'draws.csv'), 'balanced_draws_sha256': sha(new_dir/'draws.csv')}}
    result = {'created_utc': datetime.now(timezone.utc).isoformat(), 'analysis_code_sha256': sha(__file__),
        'analysis_type': 'POST HOC descriptive arithmetic sampling-correction diagnostic; not a new gate.',
        'sampling_freeze_sha256': sha(frozen_path), 'finite_population_factor': FPC,
        'formulas': {'uniform_exact': 'vu_g = (1-400/4000) * population_sample_variance_g /400',
            'nominal_iid_trace': 'mean(draw_variance_trace_ratio) * sum_g(vu_g) /0.9',
            'actual_mc_error': 'sum_g mean_over_200_draws[(sample_mean_probability_g-population_probability_g)^2]',
            'units': 'CPM-squared values multiply probability-squared inputs by1e12. Gene sums are restricted to the same populationCPM>5 gate.'},
        'contexts': output,
        'official_logbulk_or_global_capped_correction_computed': False,
        'new_selection_gate': False, 'oracle_heldout_performance_read': False,
        'interpretation': 'Without-replacement and auxiliary-balanced draws are dependent designs. Individual-cell sample variance/n is an iid estimator, not the exact variance of these sample means. Uniform sampling already needs the0.9 finite-population correction; balancing changes additional covariance terms. A similar cell-variance trace does not imply the same group-mean sampling error.',
        'limitations': ['This is a linear arithmetic-CPM mean diagnostic, not the official log1p-bulk delete-one correction or its across-target cap.',
            'Monte Carlo ratios use200 frozen draws and have finite simulation error; no new confidence intervals or promotion threshold are inferred.',
            'Exact uniform/depth design variance is available; no closed-form balanced-design variance is asserted.',
            'No sampled counts, source files, old artifacts, candidate settings, or oracle outcomes were changed.']}
    path = ROOT/'reports/balanced_sampling_correction_audit.json'
    path.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    lines = ['# Post hoc arithmetic sampling-correction audit', '',
        'For each policy, compare the actual error over all200 frozen draws with the nominal independent-cell sample-variance/n term on population-CPM>5 genes. This is **not the official log-bulk/global-capped MSE correction** and changes no gate.', '',
        '| Context | Policy | Actual mean-error sum (CPM²) | Nominal iid term (CPM²) | Nominal / actual |',
        '| --- | --- | ---: | ---: | ---: |']
    for context, item in output.items():
        for policy, stats in item['policies'].items():
            lines.append(f'| {context} | {policy} | {stats["actual_mc_sum_arithmetic_mean_mse_CPM2"]:.3f} | '
                f'{stats["nominal_iid_sum_sample_variance_over_n_CPM2"]:.3f} | {stats["nominal_iid_to_actual_mc_mse_ratio"]:.4f} |')
    lines += ['', result['interpretation'], '',
        'The nominal term is recovered without rereading counts: the stored variance-trace ratio multiplies the known population variance. Specifically, `nominal = mean(trace_ratio) × sum(uniform_exact_variance) /0.9`. Actual error is the sum of per-gene Monte Carlo mean squared errors, with matching axes and gates verified.', '',
        'An iid correction can exceed actual balanced-sample mean error even when cell-level variances look unchanged. This does not demonstrate inflation of the official score: its transformed bulk jackknife and global cap are different calculations. The finite-draw mean bias in the JSON is simulation noise, not a finding of systematic selection bias.', '',
        '[Full formulas, values and input hashes](balanced_sampling_correction_audit.json). No held-out oracle performance was read.']
    path.with_suffix('.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({c: {p: s['nominal_iid_to_actual_mc_mse_ratio'] for p,s in x['policies'].items()}
                      for c,x in output.items()}, indent=2), flush=True)


if __name__ == '__main__':
    main()
