#!/usr/bin/env python3
"""Freeze, run and summarize a two-policy control-only sampling experiment."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ.setdefault(key,'1')
from pathlib import Path
import argparse
from datetime import datetime, timezone
import hashlib
import json
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from build_prediction import ROOT, sha256
from stratified_donor_sampling import depth_blocks, select_donors, exact_mean_sampling_variance

OUT=ROOT/'evaluations/sampling_precision_400'
CONTEXTS=('H1','A','B','C')
POLICIES=('uniform','depth_stratified')
SEEDS=tuple(range(20261025,20261025+200))


def write(path,value):
    path=Path(path)
    if path.exists():raise FileExistsError(path)
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def seed_for(context,seed_id):
    return int.from_bytes(hashlib.sha256(f'sampling-precision-v1:{context}:{seed_id}'.encode()).digest()[:4],'little')


def freeze():
    if OUT.exists():raise FileExistsError(OUT)
    inputs={}
    for context in CONTEXTS:
        path=ROOT/'evaluations/dispersion_h1_training_400/context_H1.h5ad' if context=='H1' else ROOT/f'evaluations/control_preserving_null_abc/{context}/donor_population.h5ad'
        a=ad.read_h5ad(path,backed='r')
        try:
            if a.n_obs!=4000 or set(a.obs.target_gene.astype(str))!={'non-targeting'} or not a.obs_names.is_unique or not a.var_names.is_unique:
                raise ValueError('Unexpected donor population')
            inputs[context]={'path':str(path),'sha256':sha256(path),'shape':list(a.shape),
                'obs_axis_sha256':hashlib.sha256('\0'.join(a.obs_names).encode()).hexdigest(),
                'gene_axis_sha256':hashlib.sha256('\0'.join(a.var_names).encode()).hexdigest()}
        finally:a.file.close()
    files=['scripts/benchmark_sampling_precision.py','scripts/stratified_donor_sampling.py',
           'tests/test_stratified_donor_sampling.py','docs/SAMPLING_PRECISION_PROTOCOL_20261021.md']
    r={'created_at_utc':datetime.now(timezone.utc).isoformat(),'inputs':inputs,'seed_ids':list(SEEDS),
       'policies':list(POLICIES),'files':{p:sha256(ROOT/p) for p in files},
       'expression_outcomes_read_before_freeze':False,'perturbed_responses_accessed':False,
       'purpose':'Sampling precision on controls only; not a submission artifact'}
    write(OUT/'freeze.json',r);write(ROOT/'reports/sampling_precision_freeze.json',r)
    print(json.dumps({'frozen':True,'contexts':list(CONTEXTS),'draws_per_policy':len(SEEDS)}),flush=True)


def check_freeze():
    f=json.loads((OUT/'freeze.json').read_text())
    for path,digest in f['files'].items():
        if sha256(ROOT/path)!=digest:raise ValueError('Frozen code/protocol changed: '+path)
    return f


def run(context):
    f=check_freeze(); folder=OUT/context
    if folder.exists():raise FileExistsError(folder)
    inp=f['inputs'][context]
    if sha256(inp['path'])!=inp['sha256']:raise ValueError('Control bytes changed')
    a=ad.read_h5ad(inp['path']);x=sparse.csr_matrix(a.X).toarray()
    if not np.isfinite(x).all() or (x<0).any() or (x!=np.floor(x)).any():raise ValueError('Invalid counts')
    d=x.sum(1,dtype=np.int64)
    if (d<=0).any() or (d>1000000).any():raise ValueError('Invalid depth')
    q=x/d[:,None];p=q.mean(0);b=x.sum(0,dtype=np.float64)/d.sum()
    detection=(x>0).mean(0);variance=q.var(0,ddof=1);gate=p*1e6>5
    blocks=depth_blocks(d)
    vu,vs=exact_mean_sampling_variance(q,blocks)
    du,ds=exact_mean_sampling_variance(d[:,None],blocks)
    annotation='ntc_id' if 'ntc_id' in a.obs else 'guide_id'
    categories={}
    for column in [annotation,'batch']:
        if column in a.obs:
            labels,encoded=np.unique(a.obs[column].astype(str),return_inverse=True)
            categories[column]=(encoded,np.bincount(encoded,minlength=len(labels))/4000)
    all_rows=[];selected=[];gene_records=[]
    for policy in POLICIES:
        probabilities=[];bulks=[];inclusions=np.zeros(4000,int)
        for i,seed_id in enumerate(SEEDS):
            seed=seed_for(context,seed_id);rows=select_donors(d,seed,policy)
            if len(np.unique(rows))!=400:raise AssertionError('Repeated donor within draw')
            if policy=='depth_stratified' and not all(r in block for r,block in zip(rows,blocks)):
                raise AssertionError('Draw violates depth blocks')
            inclusions[rows]+=1;selected.append(rows)
            sub=q[rows];mean=sub.mean(0);counts=x[rows];depth=d[rows]
            bulk=counts.sum(0,dtype=np.float64)/depth.sum()
            detected=(counts>0).mean(0)
            sv=sub[:,gate].var(0,ddof=1).sum()
            r={'context':context,'policy':policy,'seed_id':seed_id,'rng_seed':seed,
               'arithmetic_log_mse':float(np.mean((np.log1p(1e6*mean)-np.log1p(1e6*p))**2)),
               'bulk_log_mse':float(np.mean((np.log1p(50000*bulk)-np.log1p(50000*b))**2)),
               'arithmetic_probability_l1':float(np.abs(mean-p).sum()),
               'bulk_probability_l1':float(np.abs(bulk-b).sum()),
               'depth_mean':float(depth.mean()),'depth_mean_error':float(depth.mean()-d.mean()),
               'depth_q05':float(np.quantile(depth,.05)),'depth_q50':float(np.quantile(depth,.5)),
               'depth_q95':float(np.quantile(depth,.95)),
               'detection_mae':float(np.mean(np.abs(detected[gate]-detection[gate]))),
               'variance_trace_ratio':float(sv/variance[gate].sum()),
               'selected_raw_count_sha256':hashlib.sha256(np.ascontiguousarray(counts).tobytes()).hexdigest()}
            for column,(encoded,proportion) in categories.items():
                empirical=np.bincount(encoded[rows],minlength=len(proportion))/400
                r[column+'_tv']=float(np.abs(empirical-proportion).sum()/2)
            all_rows.append(r);probabilities.append(mean);bulks.append(bulk)
            if (i+1)%50==0:print(json.dumps({'context':context,'policy':policy,'draws_completed':i+1}),flush=True)
        means=np.asarray(probabilities);bulk_arr=np.asarray(bulks)
        gene_records.append(pd.DataFrame({'gene':a.var_names,'policy':policy,
            'arithmetic_population':p,'arithmetic_mc_bias':means.mean(0)-p,
            'arithmetic_mc_mse':((means-p)**2).mean(0),
            'arithmetic_exact_variance':vu if policy=='uniform' else vs,
            'bulk_population':b,'bulk_mc_bias':bulk_arr.mean(0)-b,
            'bulk_mc_mse':((bulk_arr-b)**2).mean(0)}))
    folder.mkdir(parents=True)
    pd.DataFrame(all_rows).to_csv(folder/'draws.csv',index=False)
    np.savez_compressed(folder/'selected_rows.npz',rows=np.array(selected),policies=np.repeat(POLICIES,len(SEEDS)),seed_ids=np.tile(SEEDS,2),depth_blocks=blocks)
    pd.concat(gene_records,ignore_index=True).to_parquet(folder/'gene_sampling_audit.parquet',index=False)
    summary={'context':context,'input':inp,'draws':len(all_rows),'seed_ids':list(SEEDS),
        'exact_arithmetic_mean_variance_ratio':float(vs.sum()/vu.sum()),
        'exact_depth_mean_variance_ratio':float(ds.sum()/du.sum()) if du.sum()>0 else None,
        'gated_genes':int(gate.sum()),'depth_tied_cells':int(4000-len(np.unique(d))),
        'ties':'Stable frozen input-row order','population_gene_count':a.n_vars,
        'linear_estimators_unbiased_by_design':True,'bulk_ratio_exact_unbiasedness_claimed':False,
        'files':{n:sha256(folder/n) for n in ['draws.csv','selected_rows.npz','gene_sampling_audit.parquet']}}
    write(folder/'summary.json',summary)
    write(folder/'terminal.json',{'exit_code':0,'completed_at_utc':datetime.now(timezone.utc).isoformat(),'context':context})


def paired(u,s):
    u,s=np.asarray(u,float),np.asarray(s,float)
    rng=np.random.default_rng(20261025)
    ids=rng.integers(0,len(u),size=(5000,len(u)))
    um,sm=u[ids].mean(1),s[ids].mean(1)
    return {'uniform_mean':float(u.mean()),'stratified_mean':float(s.mean()),
        'ratio':float(s.mean()/u.mean()),'ratio_95_ci':np.quantile(sm/um,[.025,.975]).tolist(),
        'gain':float((u-s).mean()),'gain_95_ci':np.quantile(um-sm,[.025,.975]).tolist(),
        'uniform_p95':float(np.quantile(u,.95)),'stratified_p95':float(np.quantile(s,.95)),
        'stratified_better_draws':int((s<u).sum()),'draws':len(u)}


def summarize():
    f=check_freeze();contexts={};records=[]
    for context in CONTEXTS:
        folder=OUT/context
        if json.loads((folder/'terminal.json').read_text())['exit_code']!=0:raise ValueError('Incomplete context')
        summary=json.loads((folder/'summary.json').read_text())
        for name,digest in summary['files'].items():
            if sha256(folder/name)!=digest:raise ValueError('Sampling artifact changed')
        data=pd.read_csv(folder/'draws.csv');measures={}
        for key in ['arithmetic_log_mse','bulk_log_mse','arithmetic_probability_l1','bulk_probability_l1','detection_mae']:
            groups=[data[data.policy==policy].set_index('seed_id')[key].reindex(SEEDS).to_numpy() for policy in POLICIES]
            if not np.isfinite(groups).all():raise ValueError('Missing/nonfinite sampling result')
            measures[key]=paired(*groups)
        for policy in POLICIES:
            rows=data[data.policy==policy]
            if len(rows)!=200 or set(rows.seed_id)!=set(SEEDS):raise ValueError('Frozen seeds changed')
        a,b=measures['arithmetic_log_mse'],measures['bulk_log_mse']
        if context=='H1':
            passed=(b['ratio']<=.95 and b['gain_95_ci'][0]>0 and a['ratio']<=1.02 and a['ratio_95_ci'][1]<=1.02)
        else:passed=all(v['ratio']<=1.05 and v['ratio_95_ci'][1]<=1.05 for v in [a,b])
        contexts[context]={'precision_gate_passed':bool(passed),'measures':measures,'population_audit':summary}
        records.extend(data.to_dict('records'))
    result={'frozen_protocol':f,'contexts':contexts,'precision_gate_passed':all(v['precision_gate_passed'] for v in contexts.values()),
        'nonzero_experiment_requires_separate_null_pass':True,'perturbed_responses_used':False,
        'interpretation':'Monte Carlo uncertainty conditional on fixed controls; not biological prediction error or a score.'}
    write(ROOT/'reports/sampling_precision_results.json',result)
    pd.DataFrame(records).to_csv(ROOT/'reports/sampling_precision_draws.csv',index=False)
    lines=['# Control-only donor sampling precision','',
        '| Context | Arithmetic log-error ratio (95% interval) | Bulk log-error ratio (95% interval) | Exact arithmetic mean variance ratio | Precision gate |',
        '|---|---:|---:|---:|---|']
    for c,r in contexts.items():
        values=[]
        for key in ['arithmetic_log_mse','bulk_log_mse']:
            m=r['measures'][key];lo,hi=m['ratio_95_ci'];values.append(f'{m["ratio"]:.4f} ({lo:.4f}–{hi:.4f})')
        lines.append(f'| {c} | '+ ' | '.join(values)+f' | {r["population_audit"]["exact_arithmetic_mean_variance_ratio"]:.4f} | {"PASS" if r["precision_gate_passed"] else "FAIL"} |')
    lines+=['','Ratios are stratified/uniform; lower is better. There are200 predeclared draws per policy/context. Arithmetic and bulk log-errors use the complete gene axis. Exact variance applies only to the linear arithmetic mean, not log error or bulk ratios.',
        '',f'All-context precision gate: **{"PASS" if result["precision_gate_passed"] else "FAIL"}**. A separate exact-policy null pass is also required before a nonzero H1 experiment.',
        '',result['interpretation'],'','[Full results](sampling_precision_results.json); [all draws](sampling_precision_draws.csv).']
    (ROOT/'reports/sampling_precision_results.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({'precision_gate_passed':result['precision_gate_passed']}),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('action',choices=['freeze','run','summarize']);ap.add_argument('--context',choices=CONTEXTS)
    args=ap.parse_args()
    if args.action=='freeze':freeze()
    elif args.action=='run':
        if args.context is None:ap.error('--context required')
        run(args.context)
    else:summarize()
