#!/usr/bin/env python3
"""Pinned DE-only execution for independently generated moment-oracle references."""
from pathlib import Path
import os
ROOT=Path(__file__).resolve().parents[1]
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','POLARS_MAX_THREADS'):
    os.environ.setdefault(key,'1')
os.environ.setdefault('NUMBA_CACHE_DIR',str(ROOT/'.cache/numba'))
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'.cache/matplotlib'))
import argparse
from datetime import datetime,timezone
import importlib.metadata
import inspect
import json
import time
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from cell_eval2.de_compute import compute_de
from prepare_injection_recovery import sha,checked_counts

TRUTH=ROOT/'evaluations/injection_recovery_400'
ELIGIBILITY=ROOT/'reports/injection_recovery_eligibility.json'
FREEZE=ROOT/'reports/injection_recovery_analysis_freeze.json'
PARAMS=dict(backend='scanpy',input_type='counts',target_sum=1e6,mean_calc='arithmetic',epsilon=1e-9,
    clip_value=None,filter_gene_min_cpm_cell=5.,fdr_scope='per_pert',threads=1,groupby='target_gene',reference='non-targeting')


def read(path):return json.loads(Path(path).read_text())


def write(path,value):
    path=Path(path)
    if path.exists():raise FileExistsError(path)
    path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def validate_eligibility(record):
    methods=record['planned_methods']
    allowed=['current','preserving_uniform']+(['preserving_balanced'] if record['balanced_gate_passed'] else [])
    if methods!=allowed or not record['frozen_before_injection_DE']:
        raise ValueError('Invalid pre-outcome method eligibility')
    if sha(record['balanced_gate_path'])!=record['balanced_gate_sha256']:
        raise ValueError('Balanced gate artifact changed')
    return methods


def freeze(candidate_freeze):
    eligibility=read(ELIGIBILITY);methods=validate_eligibility(eligibility)
    if not candidate_freeze.is_file():raise FileNotFoundError(candidate_freeze)
    truth=read(TRUTH/'manifest.json')
    if (TRUTH/'de').exists() and any((TRUTH/'de').rglob('de.parquet')):
        raise ValueError('Cannot freeze after oracle DE outcomes exist')
    paths=[Path(__file__),ROOT/'scripts/summarize_injection_recovery.py',ROOT/'tests/test_injection_recovery_summary.py',
        ROOT/'docs/INJECTION_RECOVERY_PROTOCOL_20261021.md',ROOT/'scripts/prepare_injection_recovery.py',
        ELIGIBILITY,TRUTH/'manifest.json',candidate_freeze,Path(inspect.getsourcefile(compute_de)),
        Path(inspect.getsourcefile(compute_de)).with_name('de.py'),Path(inspect.getsourcefile(compute_de)).with_name('norm.py')]
    paths.append(ROOT/'docs/INJECTION_RECOVERY_ANALYSIS_PROTOCOL_20261021.md')
    paths.append(ROOT/'outputs/injection_recovery_400/manifest.json')
    if sha(truth['population_parameters_path'])!=truth['population_parameters_sha256']:
        raise ValueError('Oracle latent parameter file changed')
    paths.append(Path(truth['population_parameters_path']))
    for family,spec in truth['families'].items():
        if sha(spec['profiles_path'])!=spec['profiles_sha256']:raise ValueError('Oracle moment/effect profile file changed')
        paths.append(Path(spec['profiles_path']))
        for seed in spec['seeds']:
            for method in methods:
                manifest_path=ROOT/f'outputs/injection_recovery_400/{family}/{method}_s{seed}/manifest.json'
                paths.append(manifest_path);candidate=read(manifest_path)
                for key in ('receipt_path','gene_profiles_path','alias_audit_path','selected_rows_path','validation_path','requests_path'):
                    if candidate.get(key):
                        if sha(candidate[key])!=candidate[key.replace('_path','_sha256')]:raise ValueError('Candidate audit snapshot changed')
                        paths.append(Path(candidate[key]))
    record=dict(created_utc=datetime.now(timezone.utc).isoformat(),files={str(p):sha(p) for p in paths},
        planned_methods=methods,de_parameters=PARAMS,truth_manifest_sha256=sha(TRUTH/'manifest.json'),
        eligibility_sha256=sha(ELIGIBILITY),candidate_freeze_path=str(candidate_freeze),
        versions={n:importlib.metadata.version(n) for n in ('cell-eval2','scanpy','numpy','scipy','polars')},
        development='template1 AND seed20261026',confirmation=['template2/seed20261026','template1/seed20261027','template2/seed20261027'],
        sign_threshold=.1,effect_bins=[.1,.25,.5],numerical_zero_lfc_tolerance=1e-10,reference_null_max_significant_calls=5,
        no_oracle_outcomes_read_before_freeze=True,scope='Moment-oracle diagnostic only; no biological/official six-metric score or promotion')
    write(FREEZE,record);print(json.dumps({'analysis_frozen':True,'methods':methods}),flush=True)


def verify_freeze():
    f=read(FREEZE)
    for p,h in f['files'].items():
        if sha(p)!=h:raise ValueError('Frozen oracle analysis/input changed: '+p)
    validate_eligibility(read(ELIGIBILITY))
    return f


def assert_groups(obj,genes,aliases,context,controls):
    if not np.array_equal(obj.var_names.astype(str),genes) or not obj.obs_names.is_unique:
        raise ValueError('Incorrect gene/cell axes')
    labels=obj.obs.target_gene.astype(str).to_numpy()
    expected=set(aliases)|({'non-targeting'} if controls else set())
    if set(labels)!=expected or set(obj.obs.context.astype(str))!={context}:
        raise ValueError('Incorrect group/context labels')
    if any(np.count_nonzero(labels==a)!=400 for a in aliases):raise ValueError('Each alias must have400cells')
    if controls and np.count_nonzero(labels=='non-targeting')!=4000:raise ValueError('Expected4000comparison controls')
    checked_counts(obj.X)


def run(family,seed,method):
    f=verify_freeze()
    if method!='reference' and method not in f['planned_methods']:raise ValueError('Method not eligible before outcomes')
    truth=read(TRUTH/'manifest.json');spec=truth['families'][family];source=spec['seeds'][str(seed)]
    bank=np.load(spec['profiles_path'],allow_pickle=False);genes=bank['genes'].astype(str);aliases=bank['aliases'].astype(str)
    out=TRUTH/'de'/family/f'seed_{seed}'/method
    if out.exists():raise FileExistsError(out)
    for key in ('reference','comparison'):
        if sha(source[key+'_path'])!=source[key+'_sha256']:raise ValueError('Independent truth counts changed')
    candidate=None
    if method!='reference':
        mp=ROOT/f'outputs/injection_recovery_400/{family}/{method}_s{seed}/manifest.json';candidate=read(mp)
        if (candidate['family'],candidate['seed'],candidate['method'])!=(family,seed,method):raise ValueError('Candidate arm identity differs')
        if candidate['candidate_freeze_sha256']!=sha(f['candidate_freeze_path']) or candidate['planned_eligibility_sha256']!=f['eligibility_sha256']:
            raise ValueError('Candidate freeze or eligibility differs')
        for key in ('inference','comparison','reference'):
            if candidate[key+'_sha256']!=source[key+'_sha256']:raise ValueError('Candidate control/truth inputs differ')
        if not candidate['complete']:
            write(out/'status.json',dict(status='incomplete_candidate_not_scored',family=family,seed=seed,method=method,
                candidate_manifest=str(mp),candidate_manifest_sha256=sha(mp),failures=candidate.get('failures',[])))
            print(json.dumps({'family':family,'seed':seed,'method':method,'scored':False,'reason':'incomplete candidate'}),flush=True);return
        prediction=Path(candidate['prediction_path'])
        validation=read(candidate['validation_path'])
        if sha(candidate['validation_path'])!=candidate['validation_sha256'] or not validation['valid'] or not validation['hash_verified']:
            raise ValueError('Candidate lacks complete strict validation')
        if sha(prediction)!=candidate['prediction_sha256']:raise ValueError('Candidate prediction hash changed')
    reference=ad.read_h5ad(source['reference_path']);comparison=ad.read_h5ad(source['comparison_path'])
    assert_groups(reference,genes,aliases,spec['context'],True)
    if comparison.shape!=(4000,len(genes)) or not comparison.var_names.equals(reference.var_names) or set(comparison.obs.target_gene.astype(str))!={'non-targeting'}:
        raise ValueError('Incorrect comparison-control axis')
    control_rows=np.flatnonzero(reference.obs.target_gene.astype(str).to_numpy()=='non-targeting')
    if (sparse.csr_matrix(reference.X[control_rows])!=sparse.csr_matrix(comparison.X)).nnz:
        raise ValueError('Reference/comparison control count matrices differ')
    if method=='reference':obj=reference;prediction_sha=source['reference_sha256']
    else:
        pred=ad.read_h5ad(prediction);assert_groups(pred,genes,aliases,spec['context'],False)
        labels=pred.obs.target_gene.astype(str).tolist()+['non-targeting']*4000
        x,_=checked_counts(sparse.vstack([pred.X,comparison.X],format='csr'))
        obj=ad.AnnData(x,obs=pd.DataFrame({'target_gene':pd.Categorical(labels)},index=[f'run_{i}' for i in range(len(labels))]),var=pd.DataFrame(index=genes))
        prediction_sha=candidate['prediction_sha256']
    out.mkdir(parents=True);started=time.time()
    write(out/'preflight.json',dict(family=family,seed=seed,method=method,prediction_sha256=prediction_sha,
        reference_sha256=source['reference_sha256'],comparison_sha256=source['comparison_sha256'],
        truth_manifest_sha256=sha(TRUTH/'manifest.json'),analysis_freeze_sha256=sha(FREEZE),eligibility_sha256=sha(ELIGIBILITY),
        actual_counts_and_axes_verified=True,shape=list(obj.shape),parameters=PARAMS,
        candidate_manifest_sha256=None if candidate is None else sha(mp)))
    try:
        table=compute_de(obj,**PARAMS)
        if set(table['target'].to_list())!=set(aliases):raise ValueError('DE alias axis incomplete')
        table.write_parquet(out/'de.parquet')
        write(out/'status.json',dict(status='complete',family=family,seed=seed,method=method,de_sha256=sha(out/'de.parquet'),
            elapsed_seconds=time.time()-started,aliases=len(aliases),outcomes_not_printed=True))
    except Exception as error:
        write(out/'status.json',dict(status='failed',family=family,seed=seed,method=method,error_type=type(error).__name__,error=str(error)))
        raise
    print(json.dumps({'family':family,'seed':seed,'method':method,'DE_complete':True,'seconds':time.time()-started}),flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('stage',choices=['freeze','run'])
    p.add_argument('--candidate-freeze',type=Path);p.add_argument('--family',choices=['composition','dirichlet_multinomial'])
    p.add_argument('--seed',type=int,choices=[20261026,20261027]);p.add_argument('--method',choices=['reference','current','preserving_uniform','preserving_balanced'])
    a=p.parse_args()
    if a.stage=='freeze':
        if a.candidate_freeze is None:p.error('--candidate-freeze required')
        freeze(a.candidate_freeze.resolve())
    else:
        if a.family is None or a.seed is None or a.method is None:p.error('--family --seed --method required')
        run(a.family,a.seed,a.method)


if __name__=='__main__':main()
