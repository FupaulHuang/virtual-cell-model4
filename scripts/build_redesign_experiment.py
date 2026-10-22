#!/usr/bin/env python3
"""Generate isolated response x emitter arms from frozen requested profiles."""
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(key, '1')
import argparse
from concurrent.futures import ProcessPoolExecutor
import copy
from datetime import datetime, timezone
import json
import multiprocessing as mp
from pathlib import Path
import time
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from build_prediction import ROOT, Writer, atlas, control_template, sha256
from build_bulk_ablation import array_sha256, emit_audited, verify_inputs
from control_preserving_emitter import emit_counts, fit_population, validate_counts
from validate_prediction_receipt import validate

STATE = {}


def initialize(raw, templates, seeds, emitter):
    STATE.update(raw=raw, templates=templates, seeds=seeds, emitter=emitter)


def worker(item):
    target, context, p, b = item
    raw, seeds = STATE['raw'], STATE['seeds']
    depth = raw.sum(1, dtype=np.int64)
    if STATE['emitter'] == 'preserving':
        fit, calibration = fit_population(raw, p, b)
        population_mean = fit.mean(0)
        population_bulk = depth @ fit / depth.sum()
    results = []
    for seed in seeds:
        if STATE['emitter'] == 'preserving':
            donor_seed = atlas.seed_for('donor:' + context + ':' + target, seed)
            emission_seed = atlas.seed_for('transport:' + context + ':' + target, seed)
            rows = np.random.default_rng(donor_seed).choice(len(raw), 400, replace=False)
            selected, expected = raw[rows], fit[rows]
            counts = emit_counts(selected, expected, emission_seed)
            expected_mean = expected.mean(0)
            expected_bulk = depth[rows] @ expected / depth[rows].sum()
            audit = dict(calibration, donor_seed=int(donor_seed), emission_seed=int(emission_seed),
                         donor_rows=rows.tolist(), donor_rows_sha256=array_sha256(rows),
                         donor_count_sha256=array_sha256(selected),
                         donor_depths_exact=bool(np.array_equal(counts.sum(1), depth[rows])),
                         edited_molecules=int(np.abs(counts.astype(np.int64)-selected).sum()//2),
                         activated_donor_zero_entries=int(np.count_nonzero((selected == 0) & (counts > 0))))
        else:
            template, selected_depths = STATE['templates'][seed]
            counts, audit, gene = emit_audited(template, selected_depths, p, b,
                                             atlas.seed_for(context + ':' + target, seed))
            expected_mean, expected_bulk = gene['mean_probability_fitted'], gene['bulk_probability_fitted']
            # Old emitter calibrates this 400-cell template, not the population.
            population_mean = population_bulk = np.full_like(p, np.nan)
        actual_depth = validate_counts(counts)
        actual_mean = (counts / actual_depth[:, None]).mean(0)
        actual_bulk = counts.sum(0, dtype=np.float64) / actual_depth.sum()
        audit.update(target=target, context=context, seed=int(seed),
                     nnz=int(np.count_nonzero(counts)), max_library=int(actual_depth.max()),
                     counts_array_sha256=array_sha256(counts),
                     mean_probability_l1=float(np.abs(actual_mean-p).sum()),
                     bulk_probability_l1=float(np.abs(actual_bulk-b).sum()),
                     sampled_expectation_mean_l1=float(np.abs(expected_mean-p).sum()),
                     sampled_expectation_bulk_l1=float(np.abs(expected_bulk-b).sum()))
        stats = dict(mean_probability_requested=p, bulk_probability_requested=b,
                     mean_probability_population_fitted=population_mean,
                     bulk_probability_population_fitted=population_bulk,
                     mean_probability_sample_expected=expected_mean,
                     bulk_probability_sample_expected=expected_bulk,
                     mean_probability_realized=actual_mean,
                     bulk_probability_realized=actual_bulk,
                     gene_count_sums_realized=counts.sum(0, dtype=np.int64))
        results.append((seed, counts, audit, stats))
    return results


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--profiles', type=Path, required=True)
    ap.add_argument('--response', choices=['current', 'adapter'], required=True)
    ap.add_argument('--emitter', choices=['current', 'preserving'], required=True)
    ap.add_argument('--output-root', type=Path, default=ROOT/'outputs/redesign_h1_400')
    ap.add_argument('--seeds', nargs='+', type=int, default=[20261021, 20261022])
    ap.add_argument('--workers', type=int, default=4)
    ap.add_argument('--target-limit', type=int)
    args = ap.parse_args()
    if args.workers < 1 or len(set(args.seeds)) != len(args.seeds):
        raise ValueError('Invalid workers/seeds')
    config_path = ROOT/'configs/dispersion_h1_400_incumbent.json'
    config = json.loads(config_path.read_text())
    incumbent = ROOT/'outputs/dispersion_h1_400_incumbent/prediction.h5ad'
    original = json.loads(incumbent.with_suffix('.receipt.json').read_text())
    verify_inputs(config, incumbent, original)
    profile_receipt = json.loads(args.profiles.with_suffix('.json').read_text())
    if profile_receipt['profiles_sha256'] != sha256(args.profiles):
        raise ValueError('Profile bytes changed')
    for name, digest in profile_receipt['input_hashes'].items():
        if sha256(ROOT/name) != digest:
            raise ValueError('Profile input changed: ' + name)
    values = np.load(args.profiles, allow_pickle=False)
    genes, targets = values['genes'].astype(str), values['targets'].astype(str)
    if args.target_limit is not None:
        if args.target_limit < 1: raise ValueError('Target limit must be positive')
        targets = targets[:args.target_limit]
    mean, bulk = values[args.response+'_mean'], values[args.response+'_bulk']
    context = 'H1'; control_path = ROOT/config['context'][0].split('=', 1)[1]
    control = ad.read_h5ad(control_path)
    if (not np.array_equal(control.var_names, genes) or set(control.obs.target_gene.astype(str)) != {'non-targeting'}
            or control.n_obs != 4000 or not np.array_equal(values['targets'], np.array(original['targets']))):
        raise ValueError('Incorrect control or target axes')
    raw = sparse.csr_matrix(control.X).toarray().astype(np.int32)
    validate_counts(raw)
    templates = {}
    if args.emitter == 'current':
        for seed in args.seeds:
            template, depths, _, _ = control_template(control_path, genes, 400, 4, seed)
            templates[seed] = (template, depths)
    writers, folders, audits, tables = {}, {}, {}, {}
    for seed in args.seeds:
        folder = args.output_root/f'{args.response}_{args.emitter}_s{seed}'
        if folder.exists(): raise FileExistsError(folder)
        folder.mkdir(parents=True)
        folders[seed] = folder
        writers[seed] = Writer(folder/'prediction.partial.h5ad', genes)
        audits[seed], tables[seed] = [], []
    start = time.time()
    try:
        with ProcessPoolExecutor(max_workers=args.workers, mp_context=mp.get_context('spawn'),
                                 initializer=initialize, initargs=(raw, templates, args.seeds, args.emitter)) as pool:
            futures, launched = {}, 0
            for i, target in enumerate(targets):
                while launched < min(len(targets), i + 2*args.workers):
                    futures[launched] = pool.submit(worker, (str(targets[launched]), context, mean[launched], bulk[launched]))
                    launched += 1
                for seed, counts, audit, stats in futures.pop(i).result():
                    writers[seed].append(counts, context, str(target))
                    audits[seed].append(audit)
                    table = pd.DataFrame(stats)
                    table.insert(0, 'gene', genes); table.insert(0, 'target', str(target)); table.insert(0, 'context', context)
                    tables[seed].append(table)
                print(json.dumps(dict(response=args.response, emitter=args.emitter, target=str(target),
                    completed=i+1, elapsed_seconds=time.time()-start)), flush=True)
        for writer in writers.values(): writer.close()
    except BaseException as error:
        for seed, writer in writers.items():
            if writer.f.id.valid: writer.f.close()
            (folders[seed]/'failure.json').write_text(json.dumps(dict(status='failed_unpromoted',
                completed_groups=len(audits[seed]), exception_type=type(error).__name__, exception=str(error)), indent=2)+'\n')
        raise
    for seed in args.seeds:
        folder, writer = folders[seed], writers[seed]
        path = folder/'prediction.h5ad'
        (folder/'prediction.partial.h5ad').rename(path)
        pd.concat(tables[seed], ignore_index=True).to_parquet(folder/'gene_profile_audit.parquet', index=False)
        (folder/'profile_audit.json').write_text(json.dumps(audits[seed], indent=2)+'\n')
        receipt = copy.deepcopy(original)
        receipt['arguments'].update(output=str(path.resolve()), seed=seed, target_limit=args.target_limit,
                                    redesign_response=args.response, redesign_emitter=args.emitter)
        if args.response == 'adapter':
            receipt['sources'] = [x for x in receipt['sources'] if 'cd4' not in x['path'].lower()]
            receipt['arguments'].update(cd4=None, cd4_weight=0., amplitude_cpm=1., amplitude_bulk=1.)
            coverage = (values['adapter_coverage_log2fc'][:len(targets)] > 0).any(1)
            receipt['same_target_coverage'] = int(coverage.sum())
            receipt['uncovered_targets'] = targets[~coverage].tolist()
        if args.emitter == 'preserving':
            receipt['arguments']['pool'] = 1
        receipt['targets'] = targets.tolist()
        receipt.update(moments=audits[seed], cells=writer.n, stored_entries=writer.nnz,
                       bytes=path.stat().st_size, output_sha256=sha256(path),
                       note='Unselected structural development candidate; no promotion or submission.')
        receipt['redesign'] = dict(response=args.response, emitter=args.emitter, seed=seed,
            engineering_subset=args.target_limit is not None, profile_path=str(args.profiles.resolve()),
            profiles_sha256=sha256(args.profiles), profile_receipt_sha256=sha256(args.profiles.with_suffix('.json')),
            base_config_sha256=sha256(config_path), incumbent_sha256=original['output_sha256'],
            population_controls=4000,
            donor_sampling=('400 without replacement per target' if args.emitter == 'preserving'
                            else '1600 without replacement per context/seed, depth-sorted and pooled by four; template shared across targets'),
            requested_profiles_identical_between_emitters=True,
            profile_audit_sha256=sha256(folder/'profile_audit.json'),
            gene_profile_audit_sha256=sha256(folder/'gene_profile_audit.parquet'),
            code_hashes={name:sha256(ROOT/name) for name in ['scripts/build_redesign_experiment.py',
                'scripts/control_preserving_emitter.py', 'scripts/build_bulk_ablation.py',
                'scripts/build_prediction.py', 'external/atlasshift/model.py', 'external/atlasshift/predict.py']},
            generated_at_utc=datetime.now(timezone.utc).isoformat())
        path.with_suffix('.receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
        checked = validate(path, path.with_suffix('.receipt.json'), excluded_source_tokens=['H1'])
        (folder/'validation.json').write_text(json.dumps(checked, indent=2)+'\n')
        print(json.dumps(dict(complete=True, output=str(path), validation=checked)), flush=True)


if __name__ == '__main__': main()
