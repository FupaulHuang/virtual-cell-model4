#!/usr/bin/env python3
"""Balanced raw-donor precision and null tests with verified frozen comparators."""
import os
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','POLARS_MAX_THREADS'):
    os.environ.setdefault(key,'1')
os.environ.setdefault('NUMBA_CACHE_DIR',str(ROOT/'.cache/numba'))
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'.cache/matplotlib'))
import argparse
from datetime import datetime,timezone
import hashlib
import inspect
import json
import time
import anndata as ad
import numpy as np
import pandas as pd
import polars as pl
from scipy import sparse
from cell_eval2.de_compute import compute_de
from build_prediction import Writer,sha256
from validate_prediction_receipt import validate
import benchmark_sampling_precision as old
import audit_stratified_sampling_null as null
import control_preserving_emitter as emitter
import balanced_donor_sampling as balanced

OUT=ROOT/'evaluations/balanced_sampling_400'
CONTEXTS=old.CONTEXTS


def read(path):return json.loads(Path(path).read_text())


def serial(value):
    if isinstance(value,np.ndarray):return value.tolist()
    if isinstance(value,np.generic):return value.item()
    if isinstance(value,dict):return {str(k):serial(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [serial(x) for x in value]
    return value


def write(path,value):
    path=Path(path)
    if path.exists():raise FileExistsError(path)
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(serial(value),indent=2,allow_nan=False)+'\n')


def freeze():
    if OUT.exists():raise FileExistsError(OUT)
    p0=old.check_freeze();nf=null.verify()
    files=['scripts/benchmark_balanced_sampling.py','scripts/balanced_donor_sampling.py',
        'tests/test_balanced_donor_sampling.py','docs/BALANCED_SAMPLING_PROTOCOL_20261021.md',
        'reports/balanced_donor_sampling_design.json',
        'scripts/benchmark_sampling_precision.py','scripts/stratified_donor_sampling.py',
        'scripts/audit_stratified_sampling_null.py','scripts/audit_control_preserving_null.py',
        'scripts/control_preserving_emitter.py','scripts/build_prediction.py','scripts/validate_prediction_receipt.py',
        'reports/sampling_precision_freeze.json','reports/sampling_precision_results.json',
        'reports/stratified_sampling_null_freeze.json','reports/stratified_sampling_null_results.json']
    records={}
    for c in CONTEXTS:
        folder=old.OUT/c;s=read(folder/'summary.json')
        if read(folder/'terminal.json')['exit_code']!=0:raise ValueError('Incomplete P0 comparator')
        for name,digest in s['files'].items():
            if sha256(folder/name)!=digest:raise ValueError('P0 comparator changed')
        if sha256(p0['inputs'][c]['path'])!=p0['inputs'][c]['sha256']:raise ValueError('Control changed')
        comp=null.OUT/c
        m=read(comp/'manifest.json')
        if p0['inputs'][c]['sha256']!=nf['sources'][c]['donor_sha256']:
            raise ValueError('Precision and null donor populations differ')
        if sha256(m['comparison'])!=m['comparison_sha256']:
            raise ValueError('Old null comparison changed')
        tables={}
        for policy in ['uniform','depth_stratified']:
            r=read(comp/(policy+'_de')/'result.json');pre=read(comp/(policy+'_de')/'preflight.json')
            if not r['complete'] or pre['parameters']!=null.PARAMETERS or sha256(comp/(policy+'_de')/'de.parquet')!=r['de_sha256']:
                raise ValueError('Null comparator/cache parameters differ')
            if (pre['prediction_sha256']!=m['policies'][policy]['prediction_sha256']
                    or pre['comparison_sha256']!=m['comparison_sha256']
                    or pre['helper_sha256']!=sha256(inspect.getsourcefile(compute_de))
                    or sha256(m['policies'][policy]['prediction'])!=pre['prediction_sha256']):
                raise ValueError('Old null preflight provenance differs')
            tables[policy]={'path':str(comp/(policy+'_de')/'de.parquet'),'sha256':r['de_sha256'],
                           'result':str(comp/(policy+'_de')/'result.json'),'result_sha256':sha256(comp/(policy+'_de')/'result.json'),
                           'preflight':str(comp/(policy+'_de')/'preflight.json'),
                           'preflight_sha256':sha256(comp/(policy+'_de')/'preflight.json')}
        records[c]={'input':p0['inputs'][c],'precision_summary':str(folder/'summary.json'),
            'precision_summary_sha256':sha256(folder/'summary.json'),
            'null_source':nf['sources'][c],'null_comparison':m['comparison'],
            'null_comparison_sha256':m['comparison_sha256'],'null_tables':tables,
            'null_manifest':str(comp/'manifest.json'),'null_manifest_sha256':sha256(comp/'manifest.json')}
    f={'created_at_utc':datetime.now(timezone.utc).isoformat(),'contexts':records,
       'files':{p:sha256(ROOT/p) for p in files},'seed_ids':list(old.SEEDS),
       'helper_sha256':sha256(inspect.getsourcefile(compute_de)),
       'new_candidate_outcomes_observed_before_freeze':False,'old_comparator_outcomes_previously_known':True,
       'feature_seed':20261026,'no_perturbed_response_access':True,'no_submission':True}
    write(OUT/'freeze.json',f);write(ROOT/'reports/balanced_sampling_freeze.json',f)
    print(json.dumps({'frozen':True,'contexts':list(CONTEXTS)}),flush=True)


def verify():
    f=read(OUT/'freeze.json')
    for p,h in f['files'].items():
        if sha256(ROOT/p)!=h:raise ValueError('Frozen code/input changed: '+p)
    if sha256(inspect.getsourcefile(compute_de))!=f['helper_sha256']:raise ValueError('DE helper changed')
    for c,r in f['contexts'].items():
        if (sha256(r['precision_summary'])!=r['precision_summary_sha256']
                or sha256(r['null_manifest'])!=r['null_manifest_sha256']
                or sha256(r['null_comparison'])!=r['null_comparison_sha256']):
            raise ValueError('Frozen comparator provenance changed: '+c)
        for spec in r['null_tables'].values():
            if sha256(spec['preflight'])!=spec['preflight_sha256']:
                raise ValueError('Old DE preflight changed')
    return f


def population(context):
    f=verify();spec=f['contexts'][context]['input']
    if sha256(spec['path'])!=spec['sha256']:raise ValueError('Frozen population changed')
    a=ad.read_h5ad(spec['path']);x=sparse.csr_matrix(a.X).toarray()
    emitter.validate_counts(x)
    if a.n_obs!=4000 or set(a.obs.target_gene.astype(str))!={'non-targeting'}:raise ValueError('Invalid population')
    return f,a,x


def prepare_features(context):
    f,a,x=population(context);folder=OUT/context
    if folder.exists():raise FileExistsError(folder)
    features,state=balanced.fit_features(x,seed=f['feature_seed'])
    if features.shape!=(4000,6) or not np.isfinite(features).all():raise ValueError('Expected fixed six balancing features')
    folder.mkdir();np.save(folder/'features.npy',features)
    write(folder/'feature_state.json',state)
    write(folder/'features_receipt.json',{'input':f['contexts'][context]['input'],
        'features_sha256':sha256(folder/'features.npy'),'state_sha256':sha256(folder/'feature_state.json'),
        'feature_seed':f['feature_seed'],'sampler_sha256':sha256(balanced.__file__),
        'created_before_candidate_draws':True})
    print(json.dumps({'features_prepared':context}),flush=True)


def load_features(context):
    folder=OUT/context;r=read(folder/'features_receipt.json')
    if (sha256(folder/'features.npy')!=r['features_sha256'] or sha256(folder/'feature_state.json')!=r['state_sha256']
            or sha256(balanced.__file__)!=r['sampler_sha256']):raise ValueError('Frozen features changed')
    return np.load(folder/'features.npy'),r


def draw_stats(x,q,d,rows,p,b,detection,variance,gate,categories):
    counts=x[rows];sub=q[rows];depth=d[rows];mean=sub.mean(0)
    bulk=counts.sum(0,dtype=np.float64)/depth.sum();detected=(counts>0).mean(0)
    r={'arithmetic_log_mse':float(np.mean((np.log1p(1e6*mean)-np.log1p(1e6*p))**2)),
       'bulk_log_mse':float(np.mean((np.log1p(50000*bulk)-np.log1p(50000*b))**2)),
       'arithmetic_probability_l1':float(np.abs(mean-p).sum()),'bulk_probability_l1':float(np.abs(bulk-b).sum()),
       'depth_mean':float(depth.mean()),'depth_mean_error':float(depth.mean()-d.mean()),
       'depth_q05':float(np.quantile(depth,.05)),'depth_q50':float(np.quantile(depth,.5)),
       'depth_q95':float(np.quantile(depth,.95)),
       'detection_mae':float(np.mean(np.abs(detected[gate]-detection[gate]))),
       'variance_trace_ratio':float(sub[:,gate].var(0,ddof=1).sum()/variance[gate].sum()),
       'selected_raw_count_sha256':hashlib.sha256(np.ascontiguousarray(counts).tobytes()).hexdigest()}
    for column,(encoded,proportion) in categories.items():
        empirical=np.bincount(encoded[rows],minlength=len(proportion))/400
        r[column+'_tv']=float(np.abs(empirical-proportion).sum()/2)
    return r,mean,bulk


def precision(context):
    f,a,x=population(context);folder=OUT/context
    if (folder/'draws.csv').exists():raise FileExistsError('Precision already attempted')
    features,receipt=load_features(context);sampler=balanced.prepare_sampler(features,cells=400)
    d=x.sum(1,dtype=np.int64);q=x/d[:,None];p=q.mean(0);b=x.sum(0,dtype=np.float64)/d.sum()
    detection=(x>0).mean(0);variance=q.var(0,ddof=1);gate=p*1e6>5;categories={}
    for column in [('ntc_id' if 'ntc_id' in a.obs else 'guide_id'),'batch']:
        if column in a.obs:
            labels,encoded=np.unique(a.obs[column].astype(str),return_inverse=True)
            categories[column]=(encoded,np.bincount(encoded,minlength=len(labels))/4000)
    prior_summary=read(f['contexts'][context]['precision_summary'])
    old_folder=Path(f['contexts'][context]['precision_summary']).parent
    for name,digest in prior_summary['files'].items():
        if sha256(old_folder/name)!=digest:raise ValueError('Old comparison artifact changed')
    comparison=pd.read_csv(old_folder/'draws.csv')
    first=old.SEEDS[0];rows=old.select_donors(d,old.seed_for(context,first),'uniform')
    fresh,_,_=draw_stats(x,q,d,rows,p,b,detection,variance,gate,categories)
    cached=comparison[(comparison.policy=='uniform')&(comparison.seed_id==first)].iloc[0]
    for key,value in fresh.items():
        if isinstance(value,str):
            if cached[key]!=value:raise ValueError('Fresh uniform count hash differs')
        elif not np.isclose(value,cached[key],atol=1e-12,rtol=1e-12):raise ValueError('Fresh uniform statistic differs: '+key)
    records=[];selections=[];audits=[];means=[];bulks=[];frequency=np.zeros(4000,int)
    for i,seed_id in enumerate(old.SEEDS):
        seed=old.seed_for(context,seed_id);rows,audit=sampler.select(seed)
        if len(rows)!=400 or len(set(rows))!=400:raise ValueError('Incorrect balanced sample size')
        record,mean,bulk=draw_stats(x,q,d,rows,p,b,detection,variance,gate,categories)
        balance=features[rows].mean(0)-features.mean(0)
        record.update(context=context,policy='balanced',seed_id=seed_id,rng_seed=seed,
                      max_feature_mean_residual=float(np.abs(balance).max()))
        records.append(record);selections.append(rows);audits.append(dict(seed_id=seed_id,**audit))
        frequency[rows]+=1;means.append(mean);bulks.append(bulk)
        if (i+1)%25==0:print(json.dumps({'context':context,'completed_draws':i+1}),flush=True)
    means=np.asarray(means);bulks=np.asarray(bulks)
    pd.DataFrame(records).to_csv(folder/'draws.csv',index=False)
    np.savez_compressed(folder/'selected_rows.npz',rows=np.asarray(selections),seed_ids=np.asarray(old.SEEDS))
    pd.DataFrame({'gene':a.var_names,'arithmetic_mc_bias':means.mean(0)-p,
        'arithmetic_mc_mse':((means-p)**2).mean(0),'bulk_mc_bias':bulks.mean(0)-b,
        'bulk_mc_mse':((bulks-b)**2).mean(0)}).to_parquet(folder/'gene_sampling_audit.parquet',index=False)
    pd.DataFrame({'source_cell':a.obs_names,'times_selected':frequency,'frequency':frequency/200}).to_csv(folder/'inclusion_frequencies.csv',index=False)
    write(folder/'sampler_audits.json',audits)
    write(folder/'precision_summary.json',{'context':context,'features':receipt,'draws':200,
        'uniform_recomputed_and_verified':True,'minimum_inclusion_frequency':float(frequency.min()/200),
        'maximum_inclusion_frequency':float(frequency.max()/200),'mean_inclusion_frequency':float(frequency.mean()/200),
        'inclusion_frequency_is_monte_carlo_not_exact_marginal':True,
        'files':{n:sha256(folder/n) for n in ['draws.csv','selected_rows.npz','gene_sampling_audit.parquet','inclusion_frequencies.csv','sampler_audits.json']}})
    write(folder/'precision_terminal.json',{'exit_code':0,'context':context})


def prepare_null(context):
    f,a,x=population(context);folder=OUT/context/'null';features,fr=load_features(context)
    if folder.exists():raise FileExistsError(folder)
    sampler=balanced.prepare_sampler(features,cells=400);source=f['contexts'][context]
    comp_path=source['null_comparison']
    if sha256(comp_path)!=source['null_comparison_sha256']:raise ValueError('Null comparison changed')
    comparison=ad.read_h5ad(comp_path)
    if not a.var_names.equals(comparison.var_names) or set(a.obs_names)&set(comparison.obs_names):raise ValueError('Control axes/disjointness failed')
    p,b=emitter.population_profiles(x);expected,fit=emitter.fit_population(x,p,b)
    if not fit['identity']:raise ValueError('Null population fit is not identity')
    aliases=source['null_source']['aliases'];folder.mkdir()
    pd.DataFrame({'gene_name':a.var_names}).to_csv(folder/'gene_names.csv',index=False)
    pd.DataFrame({'target_gene':aliases,'n_cells':400}).to_csv(folder/'pert_counts.csv',index=False)
    pending=folder/'prediction.partial.h5ad';writer=Writer(pending,a.var_names.to_numpy())
    audits=[];selections=[];dist=[];margins=[];cm=null.moments(comparison.X)
    try:
        for alias in aliases:
            seed=null.named_seed(context,alias);rows,sa=sampler.select(seed)
            counts=emitter.emit_counts(x[rows],expected[rows],(seed+1)%2**32)
            if not np.array_equal(counts,x[rows]):raise ValueError('Raw count identity failed')
            writer.append(counts,context,alias)
            audit={'context':context,'target':alias,'nnz':int(np.count_nonzero(counts)),
                   'max_library':int(counts.sum(1).max()),'every_count_identity':True,'sampler':sa}
            audits.append(audit);selections.extend({'alias':alias,'output_row':i,'donor_row':int(row),'source_cell':str(a.obs_names[row])} for i,row in enumerate(rows))
            _,stats=null.distribution(counts,cm);dist.append({'alias':alias,**stats})
            margins.append({'alias':alias,'margins':null.margins(a.obs,comparison.obs,rows)})
        writer.close()
    except BaseException:writer.f.close();raise
    path=folder/'prediction.h5ad';pending.rename(path)
    pd.DataFrame(selections).to_csv(folder/'selected_rows.csv',index=False)
    write(folder/'distribution.json',dist);write(folder/'margins.json',margins)
    receipt={'arguments':{'genes':str(folder/'gene_names.csv'),'targets':str(folder/'pert_counts.csv'),
        'context':[context+'='+str(comp_path)],'cells':400,'target_limit':None},'sources':[],
        'targets':aliases,'output_genes':a.n_vars,'cells':writer.n,'stored_entries':writer.nnz,
        'moments':audits,'bytes':path.stat().st_size,'output_sha256':sha256(path),'feature_receipt':fr,
        'freeze_sha256':sha256(OUT/'freeze.json'),'comparison_sha256':source['null_comparison_sha256'],
        'sampling_policy':'balanced','all_count_identity':True}
    write(path.with_suffix('.receipt.json'),receipt);write(folder/'validation.json',validate(path,path.with_suffix('.receipt.json')))
    print(json.dumps({'null_prepared':context,'groups':len(aliases),'all_identity':True}),flush=True)


def de(context):
    f=verify();folder=OUT/context/'null';out=folder/'de';r=read(folder/'prediction.receipt.json')
    if out.exists():raise FileExistsError(out)
    path=folder/'prediction.h5ad';spec=f['contexts'][context]
    if sha256(path)!=r['output_sha256'] or sha256(spec['null_comparison'])!=spec['null_comparison_sha256']:raise ValueError('Null counts changed')
    a=ad.read_h5ad(path);c=ad.read_h5ad(spec['null_comparison'])
    obj=ad.AnnData(sparse.vstack([a.X,c.X],format='csr'),
        obs=pd.DataFrame({'target_gene':pd.Categorical(a.obs.target_gene.astype(str).tolist()+['non-targeting']*4000)},index=[f'bnull_{i}' for i in range(8000)]),var=a.var.copy())
    out.mkdir();write(out/'preflight.json',{'parameters':null.PARAMETERS,'prediction_sha256':sha256(path),
        'comparison_sha256':sha256(spec['null_comparison']),'helper_sha256':f['helper_sha256']})
    table=compute_de(obj,**null.PARAMETERS);table.write_parquet(out/'de.parquet')
    write(out/'result.json',{'complete':True,'context':context,'counts':null.counts(table,r['targets']),
        'de_sha256':sha256(out/'de.parquet')})
    print(json.dumps({'null_DE_complete':context}),flush=True)


def compare(u,c):
    r=old.paired(u,c)
    for key in list(r):
        if key.startswith('stratified_'):r['candidate_'+key[len('stratified_'):]]=r.pop(key)
        elif key.startswith('uniform_'):r['baseline_'+key[len('uniform_'):]]=r.pop(key)
    return r


def summarize():
    f=verify();contexts={}
    for c in CONTEXTS:
        folder=OUT/c;s=read(folder/'precision_summary.json')
        if read(folder/'precision_terminal.json')['exit_code']!=0:raise ValueError('Incomplete precision run')
        for name,digest in s['files'].items():
            if sha256(folder/name)!=digest:raise ValueError('Balanced precision artifact changed')
        old_folder=Path(f['contexts'][c]['precision_summary']).parent
        old_summary=read(f['contexts'][c]['precision_summary'])
        for name,digest in old_summary['files'].items():
            if sha256(old_folder/name)!=digest:raise ValueError('Old precision artifact changed before analysis')
        before=pd.read_csv(old_folder/'draws.csv')
        after=pd.read_csv(folder/'draws.csv').set_index('seed_id')
        if len(after)!=200 or set(after.index)!=set(old.SEEDS):raise ValueError('Incorrect frozen draws')
        measures={}
        for key in ['arithmetic_log_mse','bulk_log_mse','arithmetic_probability_l1','bulk_probability_l1','detection_mae']:
            measures[key]={policy:compare(before[before.policy==policy].set_index('seed_id')[key].reindex(old.SEEDS).to_numpy(),after[key].reindex(old.SEEDS).to_numpy()) for policy in ['uniform','depth_stratified']}
        a=measures['arithmetic_log_mse']['uniform'];b=measures['bulk_log_mse']['uniform']
        pg=(b['ratio']<=.95 and b['gain_95_ci'][0]>0 and a['ratio']<=1.02 and a['ratio_95_ci'][1]<=1.02) if c=='H1' else all(v['ratio']<=1.05 and v['ratio_95_ci'][1]<=1.05 for v in [a,b])
        nd=read(folder/'null/de/result.json')
        if not nd['complete'] or sha256(folder/'null/de/de.parquet')!=nd['de_sha256']:raise ValueError('Incomplete null DE')
        urec=f['contexts'][c]['null_tables']['uniform']
        if sha256(urec['result'])!=urec['result_sha256'] or sha256(urec['path'])!=urec['sha256']:raise ValueError('Uniform null baseline changed')
        aliases=f['contexts'][c]['null_source']['aliases']
        if len(aliases)!=10 or len(set(aliases))!=10:raise ValueError('Expected ten distinct null aliases')
        uniform_result=read(urec['result'])
        for path,record in [(Path(urec['path']),uniform_result),(folder/'null/de/de.parquet',nd)]:
            table=pl.read_parquet(path)
            if (set(table['target'].unique().to_list())!=set(aliases)
                    or table.select(['target','feature']).is_duplicated().any()
                    or null.counts(table,aliases)!=record['counts']):
                raise ValueError('Null count table/alias completeness mismatch')
        u=np.array([v['gate_nsig'] for v in uniform_result['counts']]);v=np.array([v['gate_nsig'] for v in nd['counts']])
        um,up=float(np.median(u)),float(np.quantile(u,.9));cm,cp=float(np.median(v)),float(np.quantile(v,.9))
        identity=read(folder/'null/prediction.receipt.json')['all_count_identity']
        ng=identity and um<=2 and up<=5 and cm<=max(2,2*um) and cp<=max(5,2*up)
        contexts[c]={'precision_gate_passed':bool(pg),'null_gate_passed':bool(ng),'measures':measures,
            'null':{'uniform_median':um,'uniform_p90':up,'balanced_median':cm,'balanced_p90':cp,'counts':nd['counts']},
            'sampling_audit':s}
    passed=all(r['precision_gate_passed'] and r['null_gate_passed'] for r in contexts.values())
    result={'frozen':f,'contexts':contexts,'combined_gate_passed':passed,
        'injection_balanced_arm_permitted':passed,'biological_model_promoted':False,'official_submission':False,
        'interpretation':'Precision and null engineering only; nonlinear bulk ratio and downstream metrics need separate evaluation.'}
    write(ROOT/'reports/balanced_sampling_results.json',result)
    lines=['# Balanced raw-donor precision and null study','',
        '| Context | Arithmetic error ratio (95% interval) | Bulk error ratio (95% interval) | Precision | Null median / p90 | Null |',
        '|---|---:|---:|---|---:|---|']
    for c,r in contexts.items():
        nums=[]
        for key in ['arithmetic_log_mse','bulk_log_mse']:
            m=r['measures'][key]['uniform'];lo,hi=m['ratio_95_ci'];nums.append(f'{m["ratio"]:.4f} ({lo:.4f}–{hi:.4f})')
        lines.append(f'| {c} | '+' | '.join(nums)+f' | {"PASS" if r["precision_gate_passed"] else "FAIL"} | {r["null"]["balanced_median"]:g} / {r["null"]["balanced_p90"]:g} | {"PASS" if r["null_gate_passed"] else "FAIL"} |')
    lines+=['',f'Combined engineering gate: **{"PASS" if passed else "FAIL"}**. Ratios are balanced/uniform; lower is better. All200 frozen draws are used in each context.',
        '',result['interpretation'],'','[Full evidence](balanced_sampling_results.json). No biological model selection or submission follows from this diagnostic.']
    (ROOT/'reports/balanced_sampling_results.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({'combined_gate_passed':passed}),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('action',choices=['freeze','features','precision','prepare-null','de','summarize']);ap.add_argument('--context',choices=CONTEXTS)
    args=ap.parse_args()
    if args.action=='freeze':freeze()
    elif args.action=='summarize':summarize()
    elif args.context is None:ap.error('--context required')
    else:{'features':prepare_features,'precision':precision,'prepare-null':prepare_null,'de':de}[args.action](args.context)
