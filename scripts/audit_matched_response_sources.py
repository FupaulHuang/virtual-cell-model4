#!/usr/bin/env python3
"""Audit non-H1 source-bank/DE compatibility without fitting or accessing H1 responses."""
from pathlib import Path
import hashlib
import json
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(8 * 1024**2), b''):
            h.update(chunk)
    return h.hexdigest()


def read(path):
    return json.loads((ROOT / path).read_text())


def load(path):
    with np.load(ROOT / path, allow_pickle=False) as data:
        return {k: data[k] for k in data.files}


def raw_record(family, path, receipt):
    data = load(path)
    if 'H1' in str(data['source']) or 'H1' in path:
        raise ValueError('H1 responses are forbidden')
    digest = sha(ROOT / path)
    assert digest == receipt['output_sha256']
    n = data['n_cells']
    return data, {'path': path, 'sha256': digest, 'source': str(data['source']),
        'stored_target_rows': len(n), 'targets_with_cells': int((n > 0).sum()),
        'targets_min20_cells': int((n >= 20).sum()), 'treated_cells': int(n.sum()),
        'gene_axis': len(data['genes']), 'measured_genes': int(data['measured_genes'].sum()),
        'available_arrays': sorted(data), 'contains_significance_values': any(k in data for k in ('p_adj', 'p_value', 'fdr')),
        'matched_control_policy': receipt.get('control_matching', receipt.get('matching')),
        'normalization': receipt['normalization'],
        'source_filter': receipt.get('source_filter', 'Fixed seeded subset: up to200 cells/target and50 controls/batch'),
        'actual_raw_counts_with_separate_arithmetic_CPM': True}


def compare_orion(family, raw, report):
    path = f'data/processed/PIE_{family}_DE_statistics.npz'
    de = load(path)
    release = read(f'reports/PIE_{family}_DE_statistics.json')
    assert np.array_equal(raw['genes'], de['genes'])
    target_index = {t: i for i, t in enumerate(raw['targets'].astype(str))}
    rows = np.asarray([target_index[t] for t in de['targets'].astype(str)])
    rn, dn = raw['n_cells'][rows], de['n_cells']
    measured = de['measured'].astype(bool) & raw['measured_genes'][None, :]
    nonself = de['targets'].astype(str)[:, None] != de['genes'].astype(str)[None, :]
    diagnostic = measured & nonself
    published_p = de['p_adj'][de['measured']]
    valid_p = np.isfinite(published_p) & (published_p >= 0) & (published_p <= 1)
    raw_mean = raw['target_mean_cpm'][rows].astype(float)
    raw_control = raw['matched_control_mean_cpm'][rows].astype(float)
    target_difference = np.abs(np.log2((de['target_mean_cpm'].astype(float) + 1) / (raw_mean + 1)))[diagnostic]
    control_difference = np.abs(np.log2((de['control_mean_cpm'].astype(float) + 1) / (raw_control + 1)))[diagnostic]
    raw_global = raw['global_control_mean_cpm'].astype(float)[None, :]
    global_difference = np.abs(np.log2((de['control_mean_cpm'].astype(float) + 1) / (raw_global + 1)))[diagnostic]
    raw_parquet = Path(release['input_path'])
    parquet = pq.ParquetFile(raw_parquet)
    example = parquet.read_row_group(0, columns=['target', 'feature', 'target_ncells', 'ref_ncells']).to_pandas()
    # Bounded metadata sample, not an assertion that every release row has the same reference count.
    sampled_ref_counts = sorted(int(x) for x in example.ref_ncells.unique())
    # The converted bank's own tested mask is authoritative; absent pairs are not negative labels.
    examples = [{'target': str(t), 'raw_n_cells': int(a), 'de_n_cells': int(b)}
                for t, a, b in zip(de['targets'], rn, dn)][:8]
    record = {'path': path, 'sha256': sha(ROOT / path), 'source': str(de['source']),
        'release_commit': release['release'], 'raw_release_commit': report['source_release'],
        'publisher_url': release['source_url'], 'source_parquet': str(raw_parquet),
        'source_parquet_sha256_recorded': release['input_sha256'],
        'source_parquet_schema': list(parquet.schema_arrow.names),
        'source_parquet_first_row_group_reference_ncells': sampled_ref_counts,
        'raw_total_control_cells': report['n_controls'], 'selected_DE_targets': len(dn),
        'shared_raw_targets': len(rows), 'same_saved_gene_axis': True,
        'DE_tested_pairs': int(de['measured'].sum()), 'raw_and_DE_measured_pairs': int(measured.sum()),
        'valid_p_adj_fraction_on_tested_pairs': float(valid_p.mean()),
        'missing_DE_pairs_preserved_by_mask': True,
        'same_target_cell_count': int((rn == dn).sum()),
        'different_target_cell_count': int((rn != dn).sum()),
        'DE_fewer_cells_than_raw': int((dn < rn).sum()), 'DE_more_cells_than_raw': int((dn > rn).sum()),
        'DE_to_raw_cell_count_ratio_min_median_max': np.quantile(dn / rn, [0, .5, 1]).tolist(),
        'target_count_examples': examples,
        'own_excluded_shared_pairs_for_mean_diagnostics': int(diagnostic.sum()),
        'median_abs_log2_ratio_of_CPM_plus1_target_means': float(np.median(target_difference)),
        'median_abs_log2_ratio_of_CPM_plus1_DEref_vs_raw_matched_controls': float(np.median(control_difference)),
        'median_abs_log2_ratio_of_CPM_plus1_DEref_vs_raw_global_controls': float(np.median(global_difference)),
        'recorded_DE_filter': release['efficacy_filter'],
        'filter_provenance_caution': 'Local DE-bank report attributes strong-KD filtering; vendored PIE Orion overlay enables it with strict group<0.30, cell<0.50, min30 defaults. Exact per-cell release membership/override command is not supplied here, so matching is not inferred from defaults alone.',
        'DE_control_policy': 'Publisher README: all non-targeting cells of each cell line; DE once per cell-line file. No target-specific batch-weighted control population is declared.',
        'raw_control_policy': report['control_matching'],
        'DE_assay_family_matches_raw': True,
        'same_perturbation_cell_population_verified': False,
        'same_target_matched_control_policy': False,
        'DE_bulk_is_arithmetic_approximation': bool(de['bulk_is_approximation']),
        'acceptable_as_confidence_labels_for_current_raw_effects': False,
        'verdict': 'REJECT AS MATCHED RAW-POPULATION DE LABELS; retain only as a separately defined efficacy-selected source view.'}
    return record


def main():
    reports = {'K562': read('data/K562_GWPS_CPM_full_statistics.receipt.json'),
               'HCT116': read('reports/HCT116_full_statistics.json'),
               'HEK293T': read('reports/HEK293T_full_statistics.json')}
    paths = {'K562': 'data/K562_GWPS_CPM_full_statistics.npz',
             'HCT116': 'data/processed/HCT116_full_statistics.npz',
             'HEK293T': 'data/processed/HEK293T_full_statistics.npz'}
    sources = {}
    for family in paths:
        raw, record = raw_record(family, paths[family], reports[family])
        if family == 'K562':
            record['DE_compatibility'] = {
                'matched_DE_table_available': False, 'acceptable_as_DE_confidence_labels': False,
                'observed_arithmetic_mean_variance_available': all(k in raw for k in ['target_mean_cpm_variance', 'matched_control_mean_cpm_variance']),
                'variance_limitation': 'Mean-variance estimates are not Wilcoxon/BH DE labels; no significance labels can be recovered from aggregate means alone.',
                'selection_receipt': 'data/K562_GWPS_CPM_full_statistics.receipt.json',
                'verdict': 'RAW RESPONSE MAGNITUDES AVAILABLE; MATCHED DE LABELS ABSENT.'}
        else:
            record['DE_compatibility'] = compare_orion(family, raw, reports[family])
        sources[family] = record
        print(family, record['DE_compatibility']['verdict'], flush=True)
    evidence = ['scripts/build_orion_raw_bank.py', 'scripts/build_orion_de_bank.py', 'scripts/gwps_subset_bank.py',
        'research/pie_orion_README.md', 'research/pie_orion_meta.json',
        'external/pie/src/pie/configs/process.yaml', 'external/pie/src/pie/configs/process/dataset/orion.yaml',
        'external/pie/src/pie/process/knockdown.py']
    result = {'created_utc': datetime.now(timezone.utc).isoformat(), 'audit_code_sha256': sha(__file__),
        'scope': 'Existing K562/HCT116/HEK293T source banks and publisher metadata only; no H1 responses, reserved responses, model fitting or score changes.',
        'sources': sources, 'evidence_sha256': {p: sha(ROOT / p) for p in evidence},
        'cross_family_matched_DE_label_set_available': False,
        'recommendation': 'Use predeclared native-effect magnitude labels for the matched source-only benchmark, explicitly not DE significance. Keep raw/DE views within one source family; do not impute untested DE pairs as negatives.',
        'requirements_for_future_true_DE_labels': ['Use exactly the raw-bank perturbation cell population and its documented matching/control policy, or rebuild raw and DE views jointly under a new frozen filter.',
            'Record per-cell/guide/batch inclusion, target and gene axes, normalization denominator, DE test/backend, multiple-testing universe and adjusted p-values.',
            'If a weighted matched-control estimand has no corresponding ordinary pooled-cell test, specify a batch-aware DE model rather than attaching global-control q-values.',
            'Do not infer significance from CPM means, approximate bulk deltas or available mean standard errors alone.'],
        'limitations': ['Raw source families have different native gene panels and arithmetic denominators; within-source response definitions are documented but cross-assay calibration is not established.',
            'The source target union contains target identities from public H1 lists, but this audit reads response measurements only from non-H1 source contexts.',
            'PIE control count metadata is sampled from one raw Parquet row group; exact matching is already rejected by population/control-policy discrepancies and is not asserted from this sample.',
            'Full original raw releases are not reread or rehashed; final bank hashes match existing acquisition/validation receipts, and each audited derived DE bank is hashed now.']}
    out = ROOT / 'reports/matched_response_source_audit.json'
    out.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    lines = ['# Matched source-response provenance audit', '',
        '**No coherent matched DE-label set is available for the three current raw banks.** Use the frozen raw effects for a response-only benchmark, with predeclared magnitude labels explicitly distinguished from DE significance.', '',
        '| Source | Raw targets with cells / eligible ≥20 | Raw gene axis | Candidate DE targets | Matched raw-population DE labels? |',
        '| --- | ---: | ---: | ---: | --- |']
    for family, source in sources.items():
        de = source['DE_compatibility']
        lines.append(f'| {family} | {source["targets_with_cells"]} / {source["targets_min20_cells"]} | {source["gene_axis"]} | {de.get("selected_DE_targets", "absent")} | No |')
    lines += ['', 'K562 contains genuine selected raw counts and separate arithmetic mean CPM/mean-variance estimates. It samples up to200 cells per target and50 controls per batch. There is no matching p-value table in this bank; means or standard errors must not be relabeled as Wilcoxon/BH confidence.', '']
    for family in ('HCT116', 'HEK293T'):
        d = sources[family]['DE_compatibility']
        lines += [f'**{family}:** {d["different_target_cell_count"]}/{d["shared_raw_targets"]} shared targets have different perturbation cell counts; median DE/raw count ratio {d["DE_to_raw_cell_count_ratio_min_median_max"][1]:.4f}. '
            f'On shared, tested, own-excluded pairs, median absolute log2(CPM+1) mean discrepancy is {d["median_abs_log2_ratio_of_CPM_plus1_target_means"]:.6f} for targets and {d["median_abs_log2_ratio_of_CPM_plus1_DEref_vs_raw_matched_controls"]:.6f} for matched controls.', '']
    lines += ['The raw Orion banks retain guide-passing cells without an efficacy filter and use target-weighted per-batch controls. Published PIE DE is a separately filtered view, with per-cell-line non-targeting controls. The local report attributes strong knockdown filtering; vendored defaults use strict group residual<0.30, cell residual<0.50 and at least30 cells, but exact release cell memberships and override commands are not recorded here. The observed cell-count and matched-control discrepancies independently prevent treating these q-values as labels for our raw-bank responses.', '',
        'PIE tables contain adjusted p-values and an explicit tested mask; missing records are not measured zeros or negative calls. Their arithmetic CP10k means can be converted to CPM, but the stored `effects_bulk_delta` is an arithmetic approximation and does not establish measured raw-count pseudobulk. Matching family names and target symbols does not resolve these population, control or estimand differences.', '',
        'For future true DE-confidence labels, build raw effects and DE jointly from identical cells, controls and normalization/test definitions. A batch-weighted control mean may require a batch-aware test, not a global-control q-value. Treat raw/DE views as one correlated source family.', '',
        'No H1 or reserved response measurements were accessed. This audit does not fit a model or alter frozen files. Full identities, comparison counts, settings and hashes are in [JSON](matched_response_source_audit.json).', '',
        'Evidence: [raw Orion builder](../scripts/build_orion_raw_bank.py), [PIE conversion](../scripts/build_orion_de_bank.py), [K562 builder](../scripts/gwps_subset_bank.py), [archived publisher README](../research/pie_orion_README.md), [PIE Orion configuration](../external/pie/src/pie/configs/process/dataset/orion.yaml).']
    out.with_suffix('.md').write_text('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()
