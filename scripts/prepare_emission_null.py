#!/usr/bin/env python3
"""Freeze and materialize a control-only diagnostic for emission-induced DE."""
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(key, '1')
from pathlib import Path
import json
from datetime import datetime, timezone
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.optimize import linprog
from build_prediction import ROOT, Writer, atlas, control_template, sha256

OUT = ROOT / 'evaluations/emission_null_400'
SEED = 20261023
GENERATOR_SEED = 20261021
ARMS = ('measured_null', 'donor_replay', 'pool4_no_fit', 'neutral_dual', 'response_probe')


def matched_marginal_counts(strata, original, capacities):
    """Minimum joint reallocation with exact batch and guide margins and capacities."""
    original, capacities = np.asarray(original, int), np.asarray(capacities, int)
    keys = [s.split('|', 1) for s in strata]
    batches = sorted({s[0] for s in keys}); guides = sorted({s[1] for s in keys})
    bi, gi = {x:i for i,x in enumerate(batches)}, {x:i for i,x in enumerate(guides)}
    a = np.zeros((len(batches) + len(guides), len(strata)))
    for j, (batch, guide) in enumerate(keys):
        a[bi[batch], j] = 1; a[len(batches) + gi[guide], j] = 1
    # Zero cost up to the desired joint count, unit cost above it. With fixed
    # total mass, minimizing excess also minimizes half the L1 deviation.
    bounds = [(0, int(min(o,c))) for o,c in zip(original,capacities)]
    bounds += [(0, int(max(c-o,0))) for o,c in zip(original,capacities)]
    fit = linprog(np.r_[np.zeros(len(strata)), np.ones(len(strata))],
                  A_eq=np.c_[a,a], b_eq=a @ original, bounds=bounds, method='highs')
    if not fit.success: raise ValueError('Exact separate batch/guide margins are infeasible: ' + fit.message)
    real = fit.x[:len(strata)] + fit.x[len(strata):]
    result = np.rint(real).astype(int)
    if (not np.allclose(real,result,atol=1e-6,rtol=0) or (result<0).any()
            or (result>capacities).any() or not np.array_equal(a@result,a@original)):
        raise ValueError('Integer marginal matching failed')
    return result


def split_quotas(counts, groups=10):
    counts = np.asarray(counts, dtype=int)
    q = np.tile(counts // groups, (groups, 1))
    room = np.full(groups, int(counts.sum() // groups - q[0].sum()))
    rng = np.random.default_rng(SEED)
    for j, remaining in enumerate(counts % groups):
        ties = rng.permutation(groups)
        chosen = ties[np.argsort(-room[ties], kind='stable')[:remaining]]
        q[chosen, j] += 1
        room[chosen] -= 1
    if (room != 0).any() or not np.array_equal(q.sum(0), counts) or not np.all(q.sum(1) == 400):
        raise ValueError('Could not balance fixed batch/guide strata into 400-cell groups')
    return q


def round_template(template, depths, seed):
    """Original emitter's rounding/repair, applied before its moment fitting."""
    expected = np.asarray(template, dtype=float) * depths[:, None]
    integer = np.floor(expected).astype(np.int32)
    fractions = expected - integer
    rng = np.random.default_rng(seed)
    for i in range(len(integer)):
        remaining = int(depths[i]) - int(integer[i].sum())
        if remaining:
            cumulative = np.cumsum(fractions[i])
            cumulative *= remaining / cumulative[-1]
            locations = np.searchsorted(cumulative, np.arange(remaining) + rng.random(), side='right')
            np.add.at(integer[i], locations, 1)
    result = atlas.repair_bulk_totals(integer, expected)
    if not np.array_equal(result.sum(1, dtype=np.int64), depths) or (result < 0).any():
        raise ValueError('Pooling-only emission lost integer depths')
    return result


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    if (OUT / 'cohort_freeze.json').exists():
        raise FileExistsError('Null cohort already frozen; do not resample')
    source_path = ROOT / 'data/h1_2025/adata_Training.h5ad'
    inference_path = ROOT / 'evaluations/dispersion_h1_training_400/context_H1.h5ad'
    source = ad.read_h5ad(source_path, backed='r')
    inference = ad.read_h5ad(inference_path)
    aliases = pd.read_csv(ROOT / 'evaluations/dispersion_h1_training_400/pert_counts.csv').target_gene.astype(str).tolist()[:10]
    if len(aliases) != 10 or inference.n_obs != 4000 or not np.array_equal(source.var_names, inference.var_names):
        raise ValueError('Unexpected inference controls or axes')
    if set(inference.obs.target_gene.astype(str)) != {'non-targeting'}:
        raise ValueError('Inference file includes perturbations')
    if not source.obs_names.is_unique or not inference.obs_names.is_unique:
        raise ValueError('Duplicate source IDs')
    source_labels = source.obs.target_gene.astype(str).to_numpy()
    strata_source = source.obs[['batch', 'guide_id']].astype(str).agg('|'.join, axis=1).to_numpy()
    strata_input = inference.obs[['batch', 'guide_id']].astype(str).agg('|'.join, axis=1)
    forbidden = np.isin(source.obs_names, inference.obs_names)
    available = pd.Series(strata_source[(source_labels == 'non-targeting') & ~forbidden]).value_counts()
    original_counts = strata_input.value_counts()
    strata = sorted(set(available.index) | set(original_counts.index))
    original_counts = original_counts.reindex(strata,fill_value=0).astype(int)
    capacities = (available.reindex(strata,fill_value=0).to_numpy() // 2).astype(int)
    chosen_counts = matched_marginal_counts(strata, original_counts.to_numpy(), capacities)
    counts = pd.Series(chosen_counts,index=strata)
    quotas = split_quotas(counts.to_numpy())
    rng = np.random.default_rng(SEED)
    comparison = []
    null_rows = [[] for _ in aliases]
    for j, (stratum, n) in enumerate(counts.items()):
        eligible = np.flatnonzero((source_labels == 'non-targeting') & (strata_source == stratum) & ~forbidden)
        if len(eligible) < 2 * n:
            raise ValueError(f'Insufficient unused controls in stratum {stratum}: {len(eligible)} < {2*n}')
        chosen = rng.choice(eligible, 2 * n, replace=False)
        comparison.extend(chosen[:n].tolist())
        at = n
        for i in range(len(aliases)):
            size = int(quotas[i, j]); null_rows[i].extend(chosen[at:at + size].tolist()); at += size
        if at != 2 * n:
            raise AssertionError('Stratum allocation lost cells')
    comparison = np.sort(np.asarray(comparison, dtype=np.int64))
    null_rows = [np.sort(np.asarray(x, dtype=np.int64)) for x in null_rows]
    all_rows = np.concatenate([comparison, *null_rows])
    if len(set(all_rows)) != 8000 or forbidden[all_rows].any() or not np.all(source_labels[all_rows] == 'non-targeting'):
        raise ValueError('Null split leaks inference cells or includes actual perturbations')
    metadata_rows = []
    for role, label, rows in [('comparison_control', 'non-targeting', comparison)] + [
        ('measured_null', alias, rows) for alias, rows in zip(aliases, null_rows)]:
        for row in rows:
            metadata_rows.append({'source_row': int(row), 'source_cell': source.obs_names[row],
                                  'role': role, 'technical_alias': label, 'source_target_gene': source_labels[row],
                                  'batch_guide_stratum': strata_source[row]})
    pd.DataFrame(metadata_rows).to_csv(OUT / 'source_rows.csv', index=False)
    freeze = {'created_at_utc': datetime.now(timezone.utc).isoformat(), 'seed': SEED,
              'generator_seed': GENERATOR_SEED, 'source': str(source_path),
              'source_bytes': source_path.stat().st_size, 'source_mtime_ns': source_path.stat().st_mtime_ns,
              'inference_controls': str(inference_path), 'inference_sha256': sha256(inference_path),
              'comparison_controls': 4000, 'null_groups': 10, 'cells_per_group': 400,
              'genes': inference.n_vars, 'technical_aliases': aliases,
              'source_rows_sha256': sha256(OUT / 'source_rows.csv'),
              'strata': counts.index.tolist(), 'stratum_counts_inference': original_counts.tolist(),
              'stratum_counts_comparison_and_pooled_null': counts.tolist(),
              'exact_separate_batch_and_guide_margins': True,
              'joint_reallocated_cells': int(np.abs(chosen_counts-original_counts.to_numpy()).sum()//2),
              'joint_total_variation_distance': float(np.abs(chosen_counts-original_counts.to_numpy()).sum()/8000),
              'stratum_counts_per_null_group': quotas.tolist(),
              'all_selected_real_cells_are_non_targeting': True, 'disjoint_from_inference': True,
              'reserved_response_data_accessed': False,
              'purpose': 'Control-only diagnostic; aliases are not gene knockdowns; never submit these artifacts.',
              'protocol_sha256': sha256(ROOT / 'docs/EMISSION_NULL_PROTOCOL_20261021.md')}
    (OUT / 'cohort_freeze.json').write_text(json.dumps(freeze, indent=2) + '\n')
    print('Metadata split frozen: 4000 comparison + 10x400 true null cells', flush=True)
    ordered = np.sort(all_rows)
    chunks = []
    for start in range(0, len(ordered), 1024):
        chunks.append(sparse.csr_matrix(source.X[ordered[start:start + 1024]]))
    raw_selected = sparse.vstack(chunks, format='csr')
    if not np.isfinite(raw_selected.data).all() or (raw_selected.data < 0).any() or (raw_selected.data != np.floor(raw_selected.data)).any():
        raise ValueError('Selected reference is not raw counts')
    positions = {int(row): i for i, row in enumerate(ordered)}
    selected_obs = source.obs.iloc[ordered].copy()
    selected_obs['source_target_gene'] = selected_obs.target_gene.astype(str)
    mapping = {int(row): alias for alias, rows in zip(aliases, null_rows) for row in rows}
    selected_obs['target_gene'] = pd.Categorical([mapping.get(int(r), 'non-targeting') for r in ordered])
    selected_obs['context'] = 'H1'
    selected_obs['diagnostic_role'] = ['measured_null' if int(r) in mapping else 'comparison_control' for r in ordered]
    for column in selected_obs:
        if isinstance(selected_obs[column].dtype, pd.CategoricalDtype):
            selected_obs[column] = selected_obs[column].cat.remove_unused_categories()
    reference = ad.AnnData(raw_selected, obs=selected_obs, var=source.var.copy())
    source.file.close()
    reference.uns['purpose'] = freeze['purpose']
    reference.write_h5ad(OUT / 'reference.h5ad', compression='gzip', compression_opts=1)
    pd.DataFrame({'gene_name': inference.var_names}).to_csv(OUT / 'gene_names.csv', index=False)
    pd.DataFrame({'target_gene': aliases, 'n_cells': 400}).to_csv(OUT / 'pert_counts.csv', index=False)
    template, depths, mean, bulk = control_template(inference_path, inference.var_names.to_numpy(), 400, 4, GENERATOR_SEED)
    raw_input = sparse.csr_matrix(inference.X, dtype=np.int32)
    library = np.asarray(raw_input.sum(1, dtype=np.int64)).ravel()
    selected = np.random.default_rng(GENERATOR_SEED).choice(4000, 1600, replace=False)
    selected = selected[np.argsort(library[selected], kind='stable')]
    pools = selected.reshape(400, 4)
    pd.DataFrame({'pool': np.repeat(np.arange(400), 4), 'within_pool': np.tile(np.arange(4), 400),
                  'inference_row': selected, 'source_cell': inference.obs_names[selected]}).to_csv(OUT / 'donor_pools.csv', index=False)
    response_path = ROOT / 'outputs/dispersion_h1_400_incumbent/prediction.h5ad'
    response_receipt = json.loads(response_path.with_suffix('.receipt.json').read_text())
    if sha256(response_path) != response_receipt['output_sha256']:
        raise ValueError('Frozen response probe bytes changed')
    response = ad.read_h5ad(response_path, backed='r')
    if not np.array_equal(response.var_names, inference.var_names):
        raise ValueError('Response probe gene axis differs')
    response_labels = response.obs.target_gene.astype(str).to_numpy()
    manifest = dict(freeze, reference_sha256=sha256(OUT / 'reference.h5ad'),
                    donor_pools_sha256=sha256(OUT / 'donor_pools.csv'), arms={})
    for arm in ARMS:
        folder = ROOT / f'outputs/emission_null_400_{arm}'
        folder.mkdir(parents=True, exist_ok=True)
        path = folder / 'prediction.h5ad'
        if path.exists(): raise FileExistsError(path)
        pending = path.with_suffix('.partial.h5ad')
        if pending.exists(): raise FileExistsError(pending)
        writer = Writer(pending, inference.var_names.to_numpy())
        moments = []; donors = []
        try:
            for i, alias in enumerate(aliases):
                seed = atlas.seed_for('H1:' + alias, GENERATOR_SEED)
                if arm == 'measured_null':
                    x = raw_selected[[positions[int(r)] for r in null_rows[i]]].toarray()
                elif arm == 'donor_replay':
                    chosen = pools[np.arange(400), np.random.default_rng(seed).integers(0, 4, 400)]
                    x = raw_input[chosen].toarray()
                    donors.extend({'alias': alias, 'row': j, 'inference_row': int(r),
                                   'source_cell': inference.obs_names[r]} for j, r in enumerate(chosen))
                elif arm == 'pool4_no_fit':
                    x = round_template(template, depths, seed)
                elif arm == 'neutral_dual':
                    x = atlas.dual_moment_counts(template, mean, bulk, depths=depths, seed=seed)
                else:
                    rows = np.flatnonzero(response_labels == alias)
                    if len(rows) != 400: raise ValueError('Response probe group size differs')
                    x = response.X[rows].toarray()
                writer.append(x, 'H1', alias)
                depth = x.sum(1, dtype=np.int64)
                moments.append({'context': 'H1', 'target': alias, 'nnz': int(np.count_nonzero(x)),
                                'max_library': int(depth.max()),
                                'arithmetic_probability_l1_to_inference': float(np.abs((x / depth[:, None]).mean(0) - mean).sum()),
                                'bulk_probability_l1_to_inference': float(np.abs(x.sum(0) / x.sum() - bulk).sum())})
            writer.close()
        except BaseException:
            writer.f.close(); raise
        pending.rename(path)
        if donors:
            pd.DataFrame(donors).to_csv(folder / 'donor_rows.csv', index=False)
        receipt = {'purpose': freeze['purpose'], 'arm': arm,
                   'arguments': {'genes': str(OUT / 'gene_names.csv'), 'targets': str(OUT / 'pert_counts.csv'),
                                 'context': [f'H1={inference_path}'], 'cells': 400, 'target_limit': None,
                                 'seed': GENERATOR_SEED},
                   'sources': response_receipt['sources'] if arm == 'response_probe' else [],
                   'targets': aliases, 'output_genes': inference.n_vars,
                   'cells': writer.n, 'stored_entries': writer.nnz, 'moments': moments,
                   'bytes': path.stat().st_size, 'output_sha256': sha256(path),
                   'inference_controls_sha256': freeze['inference_sha256'],
                   'cohort_freeze_sha256': sha256(OUT / 'cohort_freeze.json'),
                   'code_sha256': {'preparation': sha256(__file__), 'atlas': sha256(ROOT / 'external/atlasshift/model.py')},
                   'probe_source_sha256': sha256(response_path) if arm == 'response_probe' else None,
                   'no_response_effects': arm != 'response_probe'}
        path.with_suffix('.receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
        manifest['arms'][arm] = {'prediction': str(path), 'sha256': receipt['output_sha256']}
        print(json.dumps({'arm': arm, 'complete': True, 'shape': [writer.n, inference.n_vars]}), flush=True)
    response.file.close()
    (OUT / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print('All null-diagnostic arms complete', flush=True)


if __name__ == '__main__':
    main()
