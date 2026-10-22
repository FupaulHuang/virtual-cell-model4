#!/usr/bin/env python3
"""Prospective control-only comparison of the actual uniform/stratified sampler."""
from pathlib import Path
import os
ROOT = Path(__file__).resolve().parents[1]
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','POLARS_MAX_THREADS'):
    os.environ.setdefault(key,'1')
os.environ.setdefault('NUMBA_CACHE_DIR',str(ROOT/'.cache/numba'))
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'.cache/matplotlib'))
import argparse
from datetime import datetime,timezone
import hashlib
import importlib.metadata
import inspect
import json
import time
import anndata as ad
import numpy as np
import pandas as pd
import polars as pl
import yaml
from scipy import sparse
from cell_eval2.de_compute import compute_de
import control_preserving_emitter as emitter
import stratified_donor_sampling as sampler
from audit_control_preserving_null import validate_counts,moments,distribution
from build_prediction import Writer,sha256
from validate_prediction_receipt import validate

OUT=ROOT/'evaluations/stratified_sampling_null'
PROTOCOL=ROOT/'docs/STRATIFIED_SAMPLING_NULL_PROTOCOL_20261021.md'
CONTEXTS=('H1','A','B','C')
POLICIES=('uniform','depth_stratified')
PARAMETERS=dict(backend='scanpy',input_type='counts',target_sum=1e6,mean_calc='arithmetic',
    epsilon=1e-9,clip_value=None,filter_gene_min_cpm_cell=5.0,fdr_scope='per_pert',threads=1,
    groupby='target_gene',reference='non-targeting')
GATE=dict(median_floor=2.,p90_floor=5.,relative_multiplier=2.,
    uniform_median_max=2.,uniform_p90_max=5.,quantile_method='linear',all_contexts_required=True,
    count='p_adj<0.05 after own-alias exclusion applied AFTER unchanged BH')


def read(path):return json.loads(Path(path).read_text())


def write(path,data):
    path=Path(path)
    if path.exists():raise FileExistsError(path)
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')


def named_seed(context,alias):
    return (int.from_bytes(hashlib.sha256(f'sampling-null-v1:{context}:{alias}'.encode()).digest()[:4],'little')+20261025)%2**32


def cached_source(context):
    """Verify old source/cache receipts; do not inspect cached call counts here."""
    if context!='H1':
        folder=ROOT/f'evaluations/control_preserving_null_abc/{context}'
        m=read(folder/'manifest.json');pre=read(folder/'measured_null_de/preflight.json');res=read(folder/'measured_null_de/result.json')
        if pre['parameters']!=PARAMETERS or pre['helper_sha256']!=sha256(inspect.getsourcefile(compute_de)):
            raise ValueError('ABC cached DE helper or parameters differ')
        donor=m['roles']['donor_population'];comp=m['roles']['comparison'];de=folder/'measured_null_de/de.parquet'
        if sha256(de)!=res['de_sha256']:raise ValueError('Cached measured-null DE changed')
        measured=m['arms']['measured_null']
        if sha256(measured['prediction'])!=measured['sha256'] or pre['input_hashes']!=[measured['sha256'],comp['sha256']]:
            raise ValueError('ABC cached DE inputs differ')
        aliases=read(Path(measured['prediction']).with_suffix('.receipt.json'))['targets']
        return dict(donor_path=donor['path'],donor_sha256=donor['sha256'],comparison_source=comp['path'],
            comparison_source_sha256=comp['sha256'],comparison_is_reference_subset=False,aliases=aliases,
            measured_de=str(de),measured_de_sha256=res['de_sha256'],cache_receipts={str(folder/'manifest.json'):sha256(folder/'manifest.json'),
            str(folder/'measured_null_de/preflight.json'):sha256(folder/'measured_null_de/preflight.json'),
            str(folder/'measured_null_de/result.json'):sha256(folder/'measured_null_de/result.json')})
    folder=ROOT/'evaluations/emission_null_400';m=read(folder/'manifest.json');cache=read(ROOT/'reports/emission_null_real_de_cache.json')
    de=folder/'measured_null_de/de_real.parquet';params=folder/'measured_null_de/run_params.yaml';meta=folder/'measured_null_de/run_meta.json'
    for path,key in [(de,'real_de_sha256'),(params,'run_params_sha256'),(meta,'run_meta_sha256')]:
        if sha256(path)!=cache[key]:raise ValueError('H1 cached measured DE changed')
    p=yaml.safe_load(params.read_text());v=read(meta)
    if (p['de']['backend']!='scanpy' or p['de']['mean_calc']!='arithmetic' or p['de']['epsilon']!=1e-9
        or p['de']['clip_value'] is not None or p['de']['fdr_scope']!='per_pert' or p['target_sum']!=1e6
        or p['filter']['filter_gene_min_cpm_cell']!=5 or p['input_type']!='counts'
        or p['pert_col']!='target_gene' or p['control']!='non-targeting' or p['control_source']!='real'
        or v['cell_eval2_version']!=importlib.metadata.version('cell-eval2') or v['resolved_de_backend']!='scanpy'):
        raise ValueError('H1 cache definitions differ')
    if cache['reference_sha256']!=m['reference_sha256'] or sha256(ROOT/'scripts/run_emission_null_de.py')!=cache['entrypoint_sha256']:
        raise ValueError('H1 cache provenance differs')
    return dict(donor_path=m['inference_controls'],donor_sha256=m['inference_sha256'],
        comparison_source=str(folder/'reference.h5ad'),comparison_source_sha256=m['reference_sha256'],
        comparison_is_reference_subset=True,aliases=m['technical_aliases'],measured_de=str(de),measured_de_sha256=cache['real_de_sha256'],
        cache_receipts={str(params):sha256(params),str(meta):sha256(meta),str(ROOT/'reports/emission_null_real_de_cache.json'):sha256(ROOT/'reports/emission_null_real_de_cache.json')},
        cache_provenance_limit='Older H1 cache records package/backend/semantic identity and hashed run parameters, not a retroactive helper source hash. It is reused descriptively; both new policy tests use the same newly pinned helper.',
        previous_joint_stratum_reallocation={k:m[k] for k in ['joint_reallocated_cells','joint_total_variation_distance','exact_separate_batch_and_guide_margins']})


def load_populations(record,backed=False):
    for p,h in [(record['donor_path'],record['donor_sha256']),(record['comparison_source'],record['comparison_source_sha256'])]:
        if sha256(p)!=h:raise ValueError('Frozen control input changed')
    donor=ad.read_h5ad(record['donor_path'],backed='r' if backed else None)
    comparison=ad.read_h5ad(record['comparison_source'],backed='r' if backed else None)
    rows=np.flatnonzero(comparison.obs.target_gene.astype(str).to_numpy()=='non-targeting')
    if donor.n_obs!=4000 or len(rows)!=4000 or set(donor.obs.target_gene.astype(str))!={'non-targeting'}:
        raise ValueError('Expected4000pure donor/comparison controls')
    if record['comparison_is_reference_subset'] and set(comparison.obs.source_target_gene.astype(str))!={'non-targeting'}:
        raise ValueError('Comparison reference contains real perturbations')
    if not donor.var_names.equals(comparison.var_names) or not donor.obs_names.is_unique or not comparison.obs_names.is_unique:
        raise ValueError('Invalid frozen control axes')
    if set(donor.obs_names)&set(comparison.obs_names[rows]):raise ValueError('Donor/comparison cell IDs overlap')
    return donor,comparison,rows


def margins(obs,comparison_obs,chosen):
    result={}
    for col in ('ntc_id','guide_id','batch'):
        if col not in obs or col not in comparison_obs:continue
        pop=obs[col].astype(str).value_counts();comp=comparison_obs[col].astype(str).value_counts();draw=obs.iloc[chosen][col].astype(str).value_counts()
        names=sorted(set(pop.index)|set(comp.index)|set(draw.index))
        a=pop.reindex(names,fill_value=0).to_numpy();b=comp.reindex(names,fill_value=0).to_numpy();c=draw.reindex(names,fill_value=0).to_numpy()
        result[col]=dict(donor_population_counts=dict(zip(names,a.tolist())),comparison_population_counts=dict(zip(names,b.tolist())),
            selected_counts=dict(zip(names,c.tolist())),population_margins_exact=bool(np.array_equal(a,b)),
            selected_vs_donor_total_variation=float(np.abs(c/400-a/4000).sum()/2),
            selected_vs_comparison_total_variation=float(np.abs(c/400-b/4000).sum()/2))
    return result


def freeze():
    if OUT.exists():raise FileExistsError('Sampling null directory exists; no resampling')
    sources={};all_rows=[];balance=[]
    for context in CONTEXTS:
        s=cached_source(context);donor,comp,rows=load_populations(s,backed=True)
        try:
            raw,depth=validate_counts(donor.X[:]);blocks=sampler.depth_blocks(depth,cells=400)
            counts=np.unique(depth,return_counts=True)[1]
            s.update(genes=donor.n_vars,comparison_rows=rows.tolist(),
                original_donor_row_order_sha256=hashlib.sha256('\0'.join(donor.obs_names).encode()).hexdigest(),
                depth_ties=dict(unique_depths=int(len(counts)),cells_in_tied_depths=int(counts[counts>1].sum()),
                    largest_tie=int(counts.max()),block_boundaries_with_equal_depth=int(np.sum(depth[blocks[:-1,-1]]==depth[blocks[1:,0]]))))
            for policy in POLICIES:
                for alias in s['aliases']:
                    seed=named_seed(context,alias);selected=sampler.select_donors(depth,seed,policy,cells=400)
                    if len(selected)!=400 or len(set(selected))!=400:raise ValueError('Sampler size/uniqueness failure')
                    if policy=='depth_stratified' and not np.all([selected[i] in blocks[i] for i in range(400)]):raise ValueError('Stratified block selection failed')
                    balance.append(dict(context=context,policy=policy,alias=alias,seed=seed,margins=margins(donor.obs,comp.obs.iloc[rows],selected)))
                    all_rows.extend(dict(context=context,policy=policy,alias=alias,seed=seed,output_row=i,donor_row=int(r),
                        source_cell=str(donor.obs_names[r]),depth=int(depth[r])) for i,r in enumerate(selected))
            sources[context]=s
        finally:donor.file.close();comp.file.close()
    OUT.mkdir(parents=True);pd.DataFrame(all_rows).to_csv(OUT/'sampled_rows.csv',index=False)
    write(OUT/'sampling_margins.json',balance)
    files=[Path(__file__),PROTOCOL,Path(sampler.__file__),Path(emitter.__file__),ROOT/'scripts/audit_control_preserving_null.py',
        ROOT/'scripts/build_prediction.py',ROOT/'scripts/validate_prediction_receipt.py',Path(inspect.getsourcefile(compute_de)),
        Path(inspect.getsourcefile(compute_de)).with_name('de.py'),Path(inspect.getsourcefile(compute_de)).with_name('norm.py')]
    record=dict(created_utc=datetime.now(timezone.utc).isoformat(),sources=sources,gate=GATE,de_parameters=PARAMETERS,
        named_seed_definition="(first4littleendian(SHA256('sampling-null-v1:'+context+':'+alias))+20261025)%2**32",
        files={str(p.resolve()):sha256(p) for p in files},sampled_rows_sha256=sha256(OUT/'sampled_rows.csv'),
        sampling_margins_sha256=sha256(OUT/'sampling_margins.json'),versions={n:importlib.metadata.version(n) for n in ('cell-eval2','scanpy','numpy','scipy','anndata')},
        new_DE_outcomes_observed_before_freeze=False,scope='Control-only pragmatic engineering gate; no performance score or reserve')
    write(OUT/'freeze.json',record);write(ROOT/'reports/stratified_sampling_null_freeze.json',record)
    print(json.dumps({'frozen':True,'selected_rows':len(all_rows),'groups':80}),flush=True)


def verify():
    f=read(OUT/'freeze.json')
    for p,h in f['files'].items():
        if sha256(p)!=h:raise ValueError('Frozen code/protocol changed: '+p)
    for name,key in [('sampled_rows.csv','sampled_rows_sha256'),('sampling_margins.json','sampling_margins_sha256')]:
        if sha256(OUT/name)!=f[key]:raise ValueError('Frozen selections changed')
    return f


def prepare(context):
    f=verify();s=f['sources'][context];folder=OUT/context
    if folder.exists():raise FileExistsError(folder)
    donor,comp,rows=load_populations(s);raw,depth=validate_counts(donor.X);dense=raw.toarray()
    comparison=comp[rows].copy()
    for col in comparison.obs:
        if isinstance(comparison.obs[col].dtype,pd.CategoricalDtype):comparison.obs[col]=comparison.obs[col].cat.remove_unused_categories()
    population=emitter.population_profiles(dense);expected,audit=emitter.fit_population(dense,*population)
    if not np.array_equal(expected,dense/depth[:,None]):raise ValueError('Zero-effect population calibration changed probabilities')
    folder.mkdir();comparison.write_h5ad(folder/'comparison.h5ad',compression='gzip',compression_opts=1)
    pd.DataFrame({'gene_name':donor.var_names}).to_csv(folder/'gene_names.csv',index=False)
    pd.DataFrame({'target_gene':s['aliases'],'n_cells':400}).to_csv(folder/'pert_counts.csv',index=False)
    all_rows=pd.read_csv(OUT/'sampled_rows.csv');all_rows=all_rows[all_rows.context==context]
    cm=moments(comparison.X);distribution_rows=[];policies={}
    for policy in POLICIES:
        directory=folder/policy;directory.mkdir();path=directory/'prediction.h5ad';pending=directory/'prediction.partial.h5ad'
        writer=Writer(pending,donor.var_names.to_numpy());audits=[]
        try:
            for alias in s['aliases']:
                selected=all_rows[(all_rows.policy==policy)&(all_rows.alias==alias)].sort_values('output_row').donor_row.to_numpy()
                repeat=sampler.select_donors(depth,named_seed(context,alias),policy,cells=400)
                if not np.array_equal(selected,repeat):raise ValueError('Sampler differs from frozen indices')
                x=emitter.emit_counts(dense[selected],expected[selected],seed=(named_seed(context,alias)+1)%2**32)
                if not np.array_equal(x,dense[selected]):raise ValueError('Zero-effect count identity failed')
                m,d=distribution(x,cm);d.update(context=context,policy=policy,alias=alias,
                    mean_detected_genes=float(np.mean(np.count_nonzero(x,axis=1))),
                    depth_quantiles=np.quantile(m['depth'],[0,.05,.5,.9,.95,1],method='linear').tolist())
                distribution_rows.append(d);writer.append(x,context,alias)
                audits.append(dict(context=context,target=alias,nnz=int(np.count_nonzero(x)),max_library=int(m['depth'].max()),
                    every_count_identity=True,selected_donors_unique=True,selected_donors=400))
            writer.close()
        except BaseException:writer.f.close();raise
        pending.rename(path)
        receipt=dict(arguments=dict(genes=str(folder/'gene_names.csv'),targets=str(folder/'pert_counts.csv'),
            context=[context+'='+str(folder/'comparison.h5ad')],cells=400,target_limit=None),sources=[],targets=s['aliases'],
            output_genes=donor.n_vars,cells=writer.n,stored_entries=writer.nnz,moments=audits,bytes=path.stat().st_size,
            output_sha256=sha256(path),policy=policy,freeze_sha256=sha256(OUT/'freeze.json'),
            every_count_identity=True,donor_comparison_disjoint=True,emitter_sha256=sha256(emitter.__file__),sampler_sha256=sha256(sampler.__file__))
        write(path.with_suffix('.receipt.json'),receipt);write(directory/'validation.json',validate(path,path.with_suffix('.receipt.json')))
        policies[policy]=dict(prediction=str(path),prediction_sha256=receipt['output_sha256'],receipt_sha256=sha256(path.with_suffix('.receipt.json')))
    write(folder/'distribution.json',distribution_rows)
    write(folder/'manifest.json',dict(context=context,comparison=str(folder/'comparison.h5ad'),comparison_sha256=sha256(folder/'comparison.h5ad'),
        policies=policies,all_identity=True,all_counts_validated=True,population_calibration_audit=audit,freeze_sha256=sha256(OUT/'freeze.json')))
    print(json.dumps({'prepared':context,'identity_groups':20,'validated':True}),flush=True)


def counts(de,aliases):
    rows=[]
    for alias in aliases:
        sub=de.filter(pl.col('target')==alias);off=sub.filter(pl.col('feature')!=alias)
        rows.append(dict(alias=alias,tested_genes=sub.height,tested_genes_after_own_exclusion=off.height,
            own_gene_rows_excluded=sub.height-off.height,raw_official_nsig=sub.filter(pl.col('p_adj')<.05).height,
            gate_nsig=off.filter(pl.col('p_adj')<.05).height))
    return rows


def run_de(context,policy):
    f=verify();folder=OUT/context;m=read(folder/'manifest.json');directory=folder/(policy+'_de')
    if directory.exists():raise FileExistsError(directory)
    p=m['policies'][policy]
    if sha256(p['prediction'])!=p['prediction_sha256'] or sha256(m['comparison'])!=m['comparison_sha256']:
        raise ValueError('Prepared null counts changed')
    prediction=ad.read_h5ad(p['prediction']);control=ad.read_h5ad(m['comparison'])
    if not prediction.var_names.equals(control.var_names):raise ValueError('Scoring gene axes differ')
    labels=prediction.obs.target_gene.astype(str).tolist()+['non-targeting']*4000
    raw,_=validate_counts(sparse.vstack([prediction.X,control.X],format='csr'))
    obj=ad.AnnData(raw,obs=pd.DataFrame({'target_gene':pd.Categorical(labels)},index=[f'null_{i}' for i in range(8000)]),var=pd.DataFrame(index=prediction.var_names))
    directory.mkdir();started=time.time()
    write(directory/'preflight.json',dict(freeze_sha256=sha256(OUT/'freeze.json'),prediction_sha256=p['prediction_sha256'],
        comparison_sha256=m['comparison_sha256'],shape=list(obj.shape),parameters=PARAMETERS,helper_sha256=sha256(inspect.getsourcefile(compute_de))))
    table=compute_de(obj,**PARAMETERS);table.write_parquet(directory/'de.parquet')
    write(directory/'result.json',dict(context=context,policy=policy,complete=True,counts=counts(table,f['sources'][context]['aliases']),
        seconds=time.time()-started,de_sha256=sha256(directory/'de.parquet')))
    print(json.dumps({'context':context,'policy':policy,'complete':True,'seconds':time.time()-started}),flush=True)


def summarize():
    f=verify();groups=[];decisions=[];cached={};identities={}
    for context in CONTEXTS:
        folder=OUT/context;m=read(folder/'manifest.json');s=f['sources'][context]
        distribution_map={(r['policy'],r['alias']):r for r in read(folder/'distribution.json')}
        vals={}
        for policy in POLICIES:
            r=read(folder/(policy+'_de')/'result.json')
            if not r['complete'] or sha256(folder/(policy+'_de')/'de.parquet')!=r['de_sha256']:raise ValueError('DE result changed')
            vals[policy]=np.asarray([x['gate_nsig'] for x in r['counts']])
            groups.extend(dict(context=context,policy=policy,**x,distribution=distribution_map[(policy,x['alias'])]) for x in r['counts'])
        if sha256(s['measured_de'])!=s['measured_de_sha256']:raise ValueError('Cached measured-null table changed')
        cached[context]=dict(path=s['measured_de'],sha256=s['measured_de_sha256'],counts=counts(pl.read_parquet(s['measured_de']),s['aliases']))
        u,st=vals['uniform'],vals['depth_stratified']
        um=float(np.median(u));up=float(np.quantile(u,.9,method='linear'));sm=float(np.median(st));sp=float(np.quantile(st,.9,method='linear'))
        identity=m['all_identity'] and m['all_counts_validated'];identities[context]=identity
        baseline=um<=2 and up<=5;median_limit=max(2,2*um);p90_limit=max(5,2*up)
        passed=identity and baseline and sm<=median_limit and sp<=p90_limit
        status='pass' if passed else ('inconclusive_uniform_baseline' if not baseline else 'fail')
        decisions.append(dict(context=context,status=status,passed=passed,identity=identity,uniform_baseline_passed=baseline,
            uniform_median=um,uniform_p90=up,stratified_median=sm,stratified_p90=sp,median_limit=median_limit,p90_limit=p90_limit))
    report=dict(created_utc=datetime.now(timezone.utc).isoformat(),scope='Pragmatic sampling-policy null triage only; no FDR guarantee or performance claim',
        freeze_sha256=sha256(OUT/'freeze.json'),sampler_sha256=sha256(sampler.__file__),emitter_sha256=sha256(emitter.__file__),
        gate=GATE,passed=all(d['passed'] for d in decisions),decisions=decisions,groups=groups,cached_measured_null=cached,
        sampling_margins=str(OUT/'sampling_margins.json'),sampling_margins_sha256=sha256(OUT/'sampling_margins.json'),
        no_response_or_reserve_access=True,alias_dependence='Shared donor/comparison populations make aliases correlated technical diagnostics')
    write(ROOT/'reports/stratified_sampling_null_results.json',report)
    flat=[{k:v for k,v in g.items() if k!='distribution'} for g in groups]
    csv=ROOT/'reports/stratified_sampling_null_counts.csv'
    if csv.exists():raise FileExistsError(csv)
    pd.DataFrame(flat).to_csv(csv,index=False)
    lines=['# Sampling-policy null audit','',
        'The actual uniform and global depth-stratified samplers each selected400of the same4000raw donors per technical group. No extra guide/batch quotas were imposed. Every count was preserved by the zero-effect emitter, and comparison controls were disjoint. Populations, selected indices, named seeds, code and thresholds were frozen before new DE outcomes.','',
        '| Context | Uniform median / p90 | Stratified median / p90 | Allowed median / p90 | Decision |',
        '| --- | ---: | ---: | ---: | --- |']
    for d in decisions:lines.append(f"| {d['context']} | {d['uniform_median']:g} / {d['uniform_p90']:g} | {d['stratified_median']:g} / {d['stratified_p90']:g} | {d['median_limit']:g} / {d['p90_limit']:g} | {d['status']} |")
    lines+=['', 'Overall engineering gate: **'+('PASS' if report['passed'] else 'NOT PASSED')+'**. Uniform median≤2/p90≤5 was separately required in every context; a failed uniform baseline is inconclusive. Quantiles use linear interpolation. These are pragmatic triage limits, not empirical FDR guarantees or promotion rules for perturbed-cell performance.','',
        'Counts above use q<0.05 after own-alias removal **after** the pinned helper’s reference-CPM>5 gate and per-alias BH. Raw official diagnostic counts are also retained. H1 cached measured-null results receive the same post-BH removal; A/B/C aliases match no genes.','',
        'The [JSON report](stratified_sampling_null_results.json) contains group distribution summaries, cached measured-null results and provenance; the [CSV](stratified_sampling_null_counts.csv) distinguishes raw and own-excluded call counts. Per-draw guide/batch counts and total-variation imbalances are in the frozen sampling-margins artifact. The earlier within-guide replay experiment is a different sampling policy. Shared populations make these ten aliases correlated technical diagnostics.']
    md=ROOT/'reports/stratified_sampling_null_results.md'
    if md.exists():raise FileExistsError(md)
    md.write_text('\n'.join(lines)+'\n')
    print(json.dumps({'complete':True,'passed':report['passed']}),flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('stage',choices=['freeze','prepare','de','summarize'])
    p.add_argument('--context',choices=CONTEXTS);p.add_argument('--policy',choices=POLICIES);a=p.parse_args()
    if a.stage in ('prepare','de') and a.context is None:p.error('--context required')
    if a.stage=='de' and a.policy is None:p.error('--policy required')
    if a.stage=='freeze':freeze()
    elif a.stage=='prepare':prepare(a.context)
    elif a.stage=='de':run_de(a.context,a.policy)
    else:summarize()


if __name__=='__main__':main()
