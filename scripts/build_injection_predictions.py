#!/usr/bin/env python3
"""Frozen candidate emission from independently generated oracle probability deltas.

No reference count generation, real perturbation response data, official scoring,
adaptive alias selection or limit relaxation. Incomplete methods are unscored.
"""
from __future__ import annotations
import os
for _key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ.setdefault(_key,'1')
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path
import sys
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from build_prediction import ROOT,Writer,atlas,control_template,sha256
from build_bulk_ablation import emit_audited,array_sha256
from control_preserving_emitter import fit_population,population_profiles,emit_counts,validate_counts
from validate_prediction_receipt import validate
import balanced_donor_sampling as balanced

TRUTH=ROOT/'evaluations/injection_recovery_400/manifest.json'
OUT=ROOT/'outputs/injection_recovery_400'
FREEZE=OUT/'freeze.json'
PREPARED=OUT/'prepared.json'
PROTOCOL=ROOT/'docs/INJECTION_PREDICTION_PROTOCOL_20261021.md'
BALANCED_RESULT=ROOT/'reports/balanced_sampling_results.json'
ELIGIBILITY=ROOT/'reports/injection_recovery_eligibility.json'
SEEDS=(20261026,20261027)
FAMILIES=('composition','dirichlet_multinomial')
BASE_METHODS=('current','preserving_uniform')
FEATURE_SEED=20261026
STAGES=('population_fitted_p','population_fitted_b','sample_expected_p','sample_expected_b',
        'bulk_projected_b','donor_baseline_p','donor_baseline_b','realized_p','realized_b',
        'realized_request_drift_p','realized_request_drift_b','realized_count_sums')


def read(path):return json.loads(Path(path).read_text())


def serial(value):
    if isinstance(value,np.ndarray):return value.tolist()
    if isinstance(value,np.generic):return value.item()
    if isinstance(value,dict):return {str(k):serial(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [serial(v) for v in value]
    return value


def write(path,value):
    path=Path(path)
    if path.exists():raise FileExistsError(path)
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(serial(value),indent=2,allow_nan=False)+'\n')


def seed_policy(family,alias,seed):
    return {'template_seed':int(atlas.seed_for('injection:template:'+family,seed)),
            'donor_seed':int(atlas.seed_for('injection:donor:'+family+':'+alias,seed)),
            'transport_seed':int(atlas.seed_for('injection:transport:'+family+':'+alias,seed))}


def checked_probability(value,name):
    x=np.asarray(value,dtype=np.float64)
    if x.ndim!=1 or not np.isfinite(x).all() or (x<0).any() or not np.isclose(x.sum(),1,atol=1e-8,rtol=0):
        raise ValueError('Invalid probability profile: '+name)
    return x


def anchor_requests(inference_p,inference_b,oracle_p,oracle_b,baseline_p,baseline_b,
                    delta_p,delta_b,is_null=False):
    """One shared request mapping. Exact zero deltas preserve input baseline bytes."""
    vectors={name:checked_probability(value,name) for name,value in
        [('inference_p',inference_p),('inference_b',inference_b),('oracle_p',oracle_p),('oracle_b',oracle_b),
         ('oracle_baseline_p',baseline_p),('oracle_baseline_b',baseline_b)]}
    shapes={x.shape for x in vectors.values()}
    if len(shapes)!=1:raise ValueError('Request profile gene axes differ')
    deltas={'p':np.asarray(delta_p,dtype=np.float64),'b':np.asarray(delta_b,dtype=np.float64)}
    for kind,delta in deltas.items():
        if delta.shape!=vectors['inference_p'].shape or not np.isfinite(delta).all():raise ValueError('Invalid declared oracle delta')
        expected=vectors['oracle_'+kind]-vectors['oracle_baseline_'+kind]
        if not np.array_equal(delta,expected):raise ValueError('Declared oracle delta differs from frozen oracle/baseline')
    zero=all(np.count_nonzero(delta)==0 for delta in deltas.values())
    if is_null and not zero:raise ValueError('NULL alias has a nonzero oracle delta')
    table={k:x.copy() for k,x in vectors.items()};audit={'zero_delta':bool(zero),'is_null':bool(is_null)}
    for kind,delta in deltas.items():
        inference=vectors['inference_'+kind];before=inference+delta
        removed=np.maximum(-before,0);after=np.maximum(before,0);mass=float(after.sum())
        if not np.isfinite(mass) or mass<=0:raise ValueError('Clipped request has no finite positive mass')
        request=inference.copy() if zero else after/mass
        actual_delta=request-inference;error=actual_delta-delta
        table.update({f'oracle_delta_{kind}':delta.copy(),f'anchored_before_clip_{kind}':before,
            f'clipped_negative_{kind}':removed,f'post_clip_{kind}':after,f'request_{kind}':request,
            f'request_delta_{kind}':actual_delta,f'request_delta_error_{kind}':error,
            f'request_minus_oracle_{kind}':request-vectors['oracle_'+kind]})
        audit[kind]={'clipped_negative_mass':float(removed.sum()),'clipped_genes':int(np.count_nonzero(removed)),
            'post_clip_mass':mass,'closure_adjustment_l1':float(np.abs(request-after).sum()),
            'request_delta_error_l1':float(np.abs(error).sum()),'request_minus_oracle_l1':float(np.abs(request-vectors['oracle_'+kind]).sum()),
            'probability_sign_audit_tolerance':1e-15,
            'delta_sign_flips':int(np.count_nonzero((delta*actual_delta<0)&(np.abs(delta)>1e-15)&(np.abs(actual_delta)>1e-15))),
            'introduced_changes_at_oracle_zero':int(np.count_nonzero((delta==0)&(np.abs(actual_delta)>1e-15))),
            'request_sha256':array_sha256(request),'exact_inference_baseline_copy':bool(np.array_equal(request,inference))}
    return table,audit


def eligible_methods(result):
    if not isinstance(result.get('combined_gate_passed'),bool) or not isinstance(result.get('injection_balanced_arm_permitted'),bool):
        raise ValueError('Final balanced precision/null eligibility is required')
    if result['combined_gate_passed']!=result['injection_balanced_arm_permitted']:
        raise ValueError('Balanced gate and injection eligibility disagree')
    return [*BASE_METHODS,*(['preserving_balanced'] if result['combined_gate_passed'] else [])]


def verified_eligibility():
    record=read(ELIGIBILITY);result=read(BALANCED_RESULT);methods=eligible_methods(result)
    if (Path(record['balanced_gate_path']).resolve()!=BALANCED_RESULT.resolve()
            or record['balanced_gate_sha256']!=sha256(BALANCED_RESULT)
            or record['balanced_gate_passed']!=result['combined_gate_passed']
            or record['planned_methods']!=methods or not record['frozen_before_injection_DE']):
        raise ValueError('Prospective injection-method eligibility differs from balanced gate')
    return methods,record


def verify_truth(manifest):
    if manifest['seeds']!=list(SEEDS) or set(manifest['families'])!=set(FAMILIES):raise ValueError('Unexpected truth families/seeds')
    for family,spec in manifest['families'].items():
        if len(spec['aliases'])!=8 or len(set(spec['aliases']))!=8:raise ValueError('Exactly eight unique aliases required')
        if sha256(spec['profiles_path'])!=spec['profiles_sha256']:raise ValueError('Oracle bank changed')
        for seed in SEEDS:
            record=spec['seeds'][str(seed)]
            if record['family']!=family or record['seed']!=seed or record['aliases']!=spec['aliases']:
                raise ValueError('Truth seed/alias identity differs')
            if not record['no_candidate_emitter_used'] or not record['no_actual_perturbed_response_used']:
                raise ValueError('Independent control-only/synthetic truth required')
            for role in ('inference','comparison','reference'):
                if sha256(record[role+'_path'])!=record[role+'_sha256']:raise ValueError('Truth role bytes changed')
            if sha256(record['receipt_path'])!=record['receipt_sha256']:raise ValueError('Truth receipt changed')


def load_inference(spec,seed,genes):
    source=spec['seeds'][str(seed)];obj=ad.read_h5ad(source['inference_path'])
    if (obj.shape!=(4000,len(genes)) or not np.array_equal(obj.var_names.astype(str),genes)
            or set(obj.obs.target_gene.astype(str))!={'non-targeting'} or set(obj.obs.context.astype(str))!={spec['context']}):
        raise ValueError('Inference controls/axes/context differ from frozen roles')
    raw=sparse.csr_matrix(obj.X).toarray()
    validate_counts(raw);raw=raw.astype(np.int32)
    if not obj.obs_names.is_unique:raise ValueError('Inference cell IDs must be unique')
    ids=np.asarray(obj.obs_names.astype(str),dtype=str)
    original=np.asarray(obj.obs['original_source_cell_id'].astype(str),dtype=str) if 'original_source_cell_id' in obj.obs else ids.copy()
    return raw,ids,original


def freeze():
    if OUT.exists():raise FileExistsError('Injection candidate output directory already exists')
    manifest=read(TRUTH);verify_truth(manifest);methods,eligibility=verified_eligibility()
    dependency_freezes={
        'reports/bulk_implementation_freeze.json':['scripts/build_bulk_ablation.py','scripts/build_prediction.py','external/atlasshift/model.py','external/atlasshift/predict.py'],
        'reports/redesign_implementation_freeze.json':['scripts/control_preserving_emitter.py'],
        'reports/balanced_sampling_freeze.json':['scripts/balanced_donor_sampling.py']}
    for receipt,names in dependency_freezes.items():
        approved=read(ROOT/receipt)
        for name in names:
            if sha256(ROOT/name)!=approved['files'][name]:raise ValueError('Previously frozen numerical dependency changed: '+name)
    files=['scripts/build_injection_predictions.py','tests/test_injection_prediction_requests.py',str(PROTOCOL.relative_to(ROOT)),
           'scripts/validate_prediction_receipt.py','scripts/control_preserving_emitter.py','scripts/balanced_donor_sampling.py',
           'scripts/build_bulk_ablation.py','scripts/build_prediction.py','external/atlasshift/model.py','external/atlasshift/predict.py',
           str(TRUTH.relative_to(ROOT)),str(BALANCED_RESULT.relative_to(ROOT)),str(ELIGIBILITY.relative_to(ROOT)),*dependency_freezes]
    for spec in manifest['families'].values():files.extend([spec['profiles_path'],spec['genes_path'],spec['targets_path']])
    record={'created_at_utc':datetime.now(timezone.utc).isoformat(),'methods':methods,'families':list(FAMILIES),'seeds':list(SEEDS),
        'files':{name:sha256(ROOT/name) for name in files},'truth_manifest_sha256':sha256(TRUTH),
        'balanced_eligibility':eligibility['balanced_gate_passed'],'balanced_result_sha256':sha256(BALANCED_RESULT),
        'planned_eligibility_path':str(ELIGIBILITY),'planned_eligibility_sha256':sha256(ELIGIBILITY),
        'requests':'One anchored array per family/seed/alias, identical across all methods; exactzero returns inferencebaselinecopy',
        'seed_policy':{'template':"atlas.seed_for('injection:template:'+family,seed)",
                       'donor':"atlas.seed_for('injection:donor:'+family+':'+alias,seed)",
                       'transport':"atlas.seed_for('injection:transport:'+family+':'+alias,seed)"},
        'feature_seed':FEATURE_SEED,'population':4000,'cells_per_alias':400,
        'candidate_outcomes_read_before_freeze':False,'no_official_six_metric_score':True,'no_submission':True}
    write(FREEZE,record);write(ROOT/'reports/injection_prediction_freeze.json',record)
    print(json.dumps({'frozen':True,'methods':methods}),flush=True)


def verify_freeze():
    record=read(FREEZE)
    for path,digest in record['files'].items():
        if sha256(ROOT/path)!=digest:raise ValueError('Frozen candidate artifact changed: '+path)
    if verified_eligibility()[0]!=record['methods']:raise ValueError('Balanced eligibility changed')
    manifest=read(TRUTH)
    return record,manifest


def prepare():
    frozen,truth=verify_freeze()
    if PREPARED.exists():raise FileExistsError(PREPARED)
    if any(OUT.glob('*/current_s*')):raise ValueError('All requests/features must be prepared before candidate generation')
    prepared={};feature_cache={}
    for family in FAMILIES:
        spec=truth['families'][family]
        with np.load(spec['profiles_path'],allow_pickle=False) as d:bank={k:d[k] for k in d.files}
        genes=bank['genes'].astype(str);aliases=bank['aliases'].astype(str)
        if aliases.tolist()!=spec['aliases']:raise ValueError('Oracle alias order differs')
        gene_csv=pd.read_csv(spec['genes_path'],dtype=str).gene_name.to_numpy()
        target_csv=pd.read_csv(spec['targets_path'])
        if (not np.array_equal(gene_csv,genes) or target_csv.target_gene.astype(str).tolist()!=aliases.tolist()
                or not np.all(target_csv.n_cells.to_numpy()==400)):
            raise ValueError('Frozen gene/target CSV axes or requested group sizes differ')
        for seed in SEEDS:
            raw,ids,original=load_inference(spec,seed,genes);own_p,own_b=population_profiles(raw)
            folder=OUT/'prepared'/family/f'seed_{seed}'
            if folder.exists():raise FileExistsError(folder)
            folder.mkdir(parents=True);arrays={'aliases':aliases,'genes':genes,'inference_p':own_p,'inference_b':own_b,
                                             'is_null':bank['is_null'],'inference_cell_ids':ids,'original_source_cell_ids':original}
            all_tables=[];audits=[]
            for i,alias in enumerate(aliases):
                table,audit=anchor_requests(own_p,own_b,bank['oracle_p'][i],bank['oracle_b'][i],bank['baseline_p'],bank['baseline_b'],
                    bank['delta_p'][i],bank['delta_b'][i],bool(bank['is_null'][i]))
                all_tables.append(table);audits.append({'alias':str(alias),'family':family,'seed':seed,**audit})
            for key in all_tables[0]:arrays[key]=np.stack([table[key] for table in all_tables])
            path=folder/'requests.npz';np.savez_compressed(path,**arrays);write(folder/'request_audit.json',audits)
            template_seed=seed_policy(family,str(aliases[0]),seed)['template_seed']
            template,depth,old_p,old_b=control_template(spec['seeds'][str(seed)]['inference_path'],genes,400,4,template_seed)
            if max(np.abs(old_p-own_p).max(),np.abs(old_b-own_b).max())>1e-12:raise ValueError('Current-template inference baseline differs')
            item={'family':family,'seed':seed,'requests_path':str(path),'requests_sha256':sha256(path),
                'request_audit_path':str(folder/'request_audit.json'),'request_audit_sha256':sha256(folder/'request_audit.json'),
                'inference_sha256':spec['seeds'][str(seed)]['inference_sha256'],
                'raw_counts_array_sha256':array_sha256(raw),'template_seed':template_seed,
                'template_sha256':array_sha256(template),'template_depths_sha256':array_sha256(depth),
                'current_baseline_arithmetic_max_abs_difference':float(np.abs(old_p-own_p).max()),
                'current_baseline_bulk_max_abs_difference':float(np.abs(old_b-own_b).max())}
            if 'preserving_balanced' in frozen['methods']:
                raw_hash=array_sha256(raw)
                if raw_hash in feature_cache:
                    features,state,reused=feature_cache[raw_hash]
                else:
                    features,state=balanced.fit_features(raw,seed=FEATURE_SEED);reused=None
                    feature_cache[raw_hash]=(features,state,f'{family}:{seed}')
                np.save(folder/'features.npy',features);write(folder/'feature_state.json',state)
                item['balancing_features']={'path':str(folder/'features.npy'),'sha256':sha256(folder/'features.npy'),
                    'state_path':str(folder/'feature_state.json'),'state_sha256':sha256(folder/'feature_state.json'),
                    'feature_matrix_sha256':array_sha256(features),'feature_seed':FEATURE_SEED,
                    'reused_identical_raw_population':reused,'prepared_before_candidate_outcomes':True}
            prepared[f'{family}:{seed}']=item
            print(json.dumps({'prepared':family,'seed':seed}),flush=True)
    record={'freeze_sha256':sha256(FREEZE),'groups':prepared,'all_prepared_before_any_candidate_generation':True}
    write(PREPARED,record)


def fit_population_audited(raw,p,b):
    captured={};previous=sys.getprofile()
    def observe(frame,event,result):
        if event=='return' and frame.f_code is fit_population.__code__:
            local=frame.f_locals
            if result is not None:
                captured['bulk_projected_b']=np.asarray(local.get('bulk',local['requested'])).copy()
    try:
        sys.setprofile(observe);fitted,audit=fit_population(raw,p,b)
    finally:sys.setprofile(previous)
    if 'bulk_projected_b' not in captured:raise RuntimeError('Actual population projection was not captured')
    depth=raw.sum(1,dtype=np.int64)
    stages={'population_fitted_p':fitted.mean(0),'population_fitted_b':depth@fitted/depth.sum(),**captured}
    return fitted,audit,stages


def failure_record(error):
    text=str(error)
    prefixes=('Requested bulk projection too large:','Moment fitting failed:','Bulk support cannot satisfy',
              'Cannot bracket the bulk projection','Bulk projection exceeds fixed limit:','Population fit failed:',
              'Cannot bracket bulk projection')
    return {'status':'FAILED','failure_kind':'FEASIBILITY' if text.startswith(prefixes) else 'ENGINEERING_ERROR',
            'exception_type':type(error).__name__,'exception':text,
            'diagnostics':getattr(error,'bulk_diagnostics',{}),'limits_relaxed':False}


def sample_rows(raw,method,seeds,sampler=None):
    if method=='preserving_uniform':
        rows=np.random.default_rng(seeds['donor_seed']).choice(len(raw),400,replace=False)
        audit={'policy':'uniform_without_replacement','inclusion_probability':400/len(raw)}
    elif method=='preserving_balanced' and sampler is not None:
        rows,audit=sampler.select(seeds['donor_seed'])
    else:raise ValueError('Invalid preserving sampling policy')
    if len(rows)!=400 or len(set(rows))!=400:raise ValueError('Incorrect donor cardinality')
    return rows,audit


def generate(family,seed):
    frozen,truth=verify_freeze();prepared=read(PREPARED)
    if prepared['freeze_sha256']!=sha256(FREEZE):raise ValueError('Prepared requests use a different freeze')
    spec=truth['families'][family];group=prepared['groups'][f'{family}:{seed}']
    for name in ('requests','request_audit'):
        if sha256(group[name+'_path'])!=group[name+'_sha256']:raise ValueError('Prepared request artifact changed')
    if sha256(spec['seeds'][str(seed)]['inference_path'])!=group['inference_sha256']:raise ValueError('Inference bytes changed')
    with np.load(group['requests_path'],allow_pickle=False) as d:requests={k:d[k] for k in d.files}
    aliases=requests['aliases'].astype(str);genes=requests['genes'].astype(str);context=spec['context']
    raw,ids,source_ids=load_inference(spec,seed,genes)
    if array_sha256(raw)!=group['raw_counts_array_sha256']:raise ValueError('Prepared raw inference population changed')
    if not np.array_equal(ids,requests['inference_cell_ids']):raise ValueError('Inference cell ordering changed')
    template,template_depth,_,_=control_template(spec['seeds'][str(seed)]['inference_path'],genes,400,4,group['template_seed'])
    if array_sha256(template)!=group['template_sha256'] or array_sha256(template_depth)!=group['template_depths_sha256']:
        raise ValueError('Frozen current template changed')
    all_depths=raw.sum(1,dtype=np.int64)
    template_rows=np.random.default_rng(group['template_seed']).choice(len(raw),1600,replace=False)
    template_rows=template_rows[np.argsort(all_depths[template_rows],kind='stable')].reshape(400,4)
    if not np.array_equal(np.rint(all_depths[template_rows].mean(1)).astype(np.int64),template_depth):
        raise ValueError('Current pooled-template donor provenance differs')
    sampler=None
    if 'preserving_balanced' in frozen['methods']:
        f=group['balancing_features']
        if sha256(f['path'])!=f['sha256'] or sha256(f['state_path'])!=f['state_sha256']:raise ValueError('Balanced feature freeze changed')
        features=np.load(f['path'],allow_pickle=False)
        if array_sha256(features)!=f['feature_matrix_sha256']:raise ValueError('Feature array identity changed')
        sampler=balanced.prepare_sampler(features,cells=400)
    folders={};writers={};audits={};tables={};donors={}
    for method in frozen['methods']:
        folder=OUT/family/f'{method}_s{seed}'
        if folder.exists():raise FileExistsError(folder)
        folder.mkdir(parents=True);folders[method]=folder;writers[method]=Writer(folder/'prediction.partial.h5ad',genes)
        audits[method]=[];tables[method]=[];donors[method]=[]
    request_audits=read(group['request_audit_path']);depths=raw.sum(1,dtype=np.int64)
    try:
        for i,alias in enumerate(aliases):
            seeds=seed_policy(family,str(alias),seed);p=requests['request_p'][i].copy();b=requests['request_b'][i].copy()
            p_hash,b_hash=array_sha256(p),array_sha256(b)
            if p_hash!=request_audits[i]['p']['request_sha256'] or b_hash!=request_audits[i]['b']['request_sha256']:
                raise ValueError('Alias request hashes differ from preparation')
            successes={};failures={}
            donors['current'].extend({'alias':str(alias),'output_row':j,'pool_member':k,'inference_row':int(r),
                'source_cell_id':str(ids[r]),'original_source_cell_id':str(source_ids[r]),'selection_role':'current_pooled_template'}
                for j,group_rows in enumerate(template_rows) for k,r in enumerate(group_rows))
            try:
                counts,calibration,stages=emit_audited(template,template_depth,p.copy(),b.copy(),seeds['transport_seed'])
                successes['current']=(counts,{**calibration,**seeds,'donor_seed_used':False,'template_seed_used':True,
                    'template_source_rows_sha256':array_sha256(template_rows),
                    'donor_baseline_semantics':'Pooled normalized template before fitting, not individual raw cells'},{
                    'sample_expected_p':stages['mean_probability_fitted'],'sample_expected_b':stages['bulk_probability_fitted'],
                    'bulk_projected_b':stages['bulk_probability_projected'],
                    'donor_baseline_p':template.mean(0),'donor_baseline_b':template_depth@template/template_depth.sum()})
            except Exception as error:failures['current']=failure_record(error)
            preserving=[m for m in frozen['methods'] if m!='current']
            try:
                fitted,calibration,population=fit_population_audited(raw,p.copy(),b.copy())
                if request_audits[i]['zero_delta'] and not calibration['identity']:raise ValueError('Zero-delta population fit is not identity')
            except Exception as error:
                for method in preserving:failures[method]=failure_record(error)
            else:
                for method in preserving:
                    try:
                        rows,selection=sample_rows(raw,method,seeds,sampler);selected=raw[rows];expected=fitted[rows]
                        counts=emit_counts(selected,expected,seeds['transport_seed'])
                        if request_audits[i]['zero_delta'] and not np.array_equal(counts,selected):raise ValueError('Zero-delta preserving counts differ from raw donors')
                        actual_depth=validate_counts(counts)
                        if not np.array_equal(actual_depth,depths[rows]):raise ValueError('Preserving donor depth changed')
                        donor_p,donor_b=population_profiles(selected)
                        stages={**population,'sample_expected_p':expected.mean(0),
                                'sample_expected_b':depths[rows]@expected/depths[rows].sum(),
                                'donor_baseline_p':donor_p,'donor_baseline_b':donor_b}
                        audit={**calibration,**seeds,'sampler':selection,'donor_rows':rows.tolist(),
                            'donor_seed_used':True,'template_seed_used':False,'donor_baseline_semantics':'Selected unmodified raw cells',
                            'donor_rows_sha256':array_sha256(rows),'donor_count_sha256':array_sha256(selected),
                            'zero_delta_raw_identity':bool(np.array_equal(counts,selected)) if request_audits[i]['zero_delta'] else None,
                            'donor_depths_exact':True}
                        donors[method].extend({'alias':str(alias),'output_row':j,'pool_member':0,'inference_row':int(r),
                            'source_cell_id':str(ids[r]),'original_source_cell_id':str(source_ids[r]),'selection_role':'preserving_raw_donor'} for j,r in enumerate(rows))
                        successes[method]=(counts,audit,stages)
                    except Exception as error:failures[method]=failure_record(error)
                del fitted
            for method in frozen['methods']:
                base={k:requests[k][i].copy() for k in requests if k not in
                    ('aliases','genes','is_null','inference_cell_ids','original_source_cell_ids') and requests[k].shape==(len(aliases),len(genes))}
                base.update({name:np.full(len(genes),np.nan) for name in STAGES})
                audit={'alias':str(alias),'target':str(alias),'context':context,'family':family,'seed':seed,'method':method,
                    'request_p_sha256':p_hash,'request_b_sha256':b_hash,'request_mapping':request_audits[i]}
                if method in successes:
                    counts,calibration,stages=successes[method];actual_depth=validate_counts(counts)
                    realized_p=(counts/actual_depth[:,None]).mean(0);sums=counts.sum(0,dtype=np.int64);realized_b=sums/sums.sum()
                    writers[method].append(counts,context,str(alias))
                    base.update(stages,realized_p=realized_p,realized_b=realized_b,realized_count_sums=sums,
                                realized_request_drift_p=realized_p-p,realized_request_drift_b=realized_b-b)
                    audit.update(calibration,status='COMPLETE',nnz=int(np.count_nonzero(counts)),max_library=int(actual_depth.max()),
                        counts_array_sha256=array_sha256(counts),realized_request_p_l1=float(np.abs(realized_p-p).sum()),
                        realized_request_b_l1=float(np.abs(realized_b-b).sum()),population_fit_applicable=method!='current')
                else:audit.update(failures[method])
                audit['requested_arrays_identical_across_methods']=True
                audits[method].append(audit);table=pd.DataFrame(base)
                for column,value in [('gene',genes),('alias',str(alias)),('status',audit['status'])]:table.insert(0,column,value)
                tables[method].append(table)
            print(json.dumps({'family':family,'seed':seed,'alias':str(alias),'status':{m:audits[m][-1]['status'] for m in frozen['methods']}}),flush=True)
        for writer in writers.values():writer.close()
    except BaseException:
        for writer in writers.values():
            if writer.f.id.valid:writer.f.close()
        raise
    arm_records={};engineering_errors=False
    for method in frozen['methods']:
        folder=folders[method];writer=writers[method];failed=[a for a in audits[method] if a['status']!='COMPLETE'];complete=not failed
        engineering_errors|=any(a.get('failure_kind')=='ENGINEERING_ERROR' for a in failed)
        path=folder/('prediction.h5ad' if complete else 'prediction.partial.h5ad')
        if complete:(folder/'prediction.partial.h5ad').rename(path)
        pd.concat(tables[method],ignore_index=True).to_parquet(folder/'gene_profiles.parquet',index=False)
        write(folder/'alias_audit.json',audits[method])
        pd.DataFrame(donors[method],columns=['alias','output_row','pool_member','inference_row','source_cell_id','original_source_cell_id','selection_role']).to_csv(folder/'selected_rows.csv',index=False)
        receipt={'arguments':{'genes':spec['genes_path'],'targets':spec['targets_path'],'context':[context+'='+spec['seeds'][str(seed)]['inference_path']],
            'cells':400,'target_limit':None},'sources':[],'targets':aliases.tolist(),'output_genes':len(genes),
            'cells':writer.n,'stored_entries':writer.nnz,'moments':[a for a in audits[method] if a['status']=='COMPLETE'],
            'bytes':path.stat().st_size,'output_sha256':sha256(path),'complete':complete,
            'note':'Independent synthetic/control-only diagnostic; never submit. Incomplete methods are unscored.',
            'truth_manifest_sha256':sha256(TRUTH),'candidate_freeze_sha256':sha256(FREEZE),'prepared_manifest_sha256':sha256(PREPARED),
            'requested_profiles_sha256':group['requests_sha256'],'family':family,'seed':seed,'method':method}
        write(folder/'prediction.receipt.json',receipt)
        checked=validate(path,folder/'prediction.receipt.json') if complete else {'valid':False,'complete':False,'reason':'Incomplete eight-alias protocol axis; DE/scoring prohibited'}
        write(folder/'validation.json',checked)
        record={'family':family,'seed':seed,'context':context,'method':method,'complete':complete,'aliases':aliases.tolist(),
            'candidate_freeze_path':str(FREEZE),'candidate_freeze_sha256':sha256(FREEZE),
            'planned_eligibility_path':str(ELIGIBILITY),'planned_eligibility_sha256':sha256(ELIGIBILITY),
            'successful_aliases':[a['alias'] for a in audits[method] if a['status']=='COMPLETE'],'failures':failed,
            'prediction_path':str(path) if complete else None,'partial_prediction_path':str(path) if not complete else None,
            'prediction_sha256':sha256(path),'receipt_path':str(folder/'prediction.receipt.json'),'receipt_sha256':sha256(folder/'prediction.receipt.json'),
            'gene_profiles_path':str(folder/'gene_profiles.parquet'),'gene_profiles_sha256':sha256(folder/'gene_profiles.parquet'),
            'alias_audit_path':str(folder/'alias_audit.json'),'alias_audit_sha256':sha256(folder/'alias_audit.json'),
            'selected_rows_path':str(folder/'selected_rows.csv'),'selected_rows_sha256':sha256(folder/'selected_rows.csv'),
            'validation_path':str(folder/'validation.json'),'validation_sha256':sha256(folder/'validation.json'),
            'inference_path':spec['seeds'][str(seed)]['inference_path'],'inference_sha256':group['inference_sha256'],
            'comparison_path':spec['seeds'][str(seed)]['comparison_path'],'comparison_sha256':spec['seeds'][str(seed)]['comparison_sha256'],
            'reference_path':spec['seeds'][str(seed)]['reference_path'],'reference_sha256':spec['seeds'][str(seed)]['reference_sha256'],
            'requests_path':group['requests_path'],'requests_sha256':group['requests_sha256'],'DE_permitted':complete,
            'full_population_fit_applicable':method!='current','no_selective_subset_performance_claim':True}
        write(folder/'manifest.json',record);arm_records[method]=record
    terminal={'family':family,'seed':seed,'exit_code':2 if engineering_errors else 0,'all_aliases_attempted':True,
              'all_methods_complete':all(r['complete'] for r in arm_records.values()),
              'methods':{m:{'manifest_path':str(folders[m]/'manifest.json'),'manifest_sha256':sha256(folders[m]/'manifest.json'),'complete':r['complete']} for m,r in arm_records.items()}}
    write(OUT/family/f'generation_s{seed}.json',terminal)
    print(json.dumps(terminal),flush=True)
    if engineering_errors:raise SystemExit(2)


def finalize():
    frozen,truth=verify_freeze();groups={};all_complete=True
    for family in FAMILIES:
        for seed in SEEDS:
            path=OUT/family/f'generation_s{seed}.json';terminal=read(path)
            if terminal['exit_code']!=0 or not terminal['all_aliases_attempted']:raise ValueError('Pending/failed generation engineering')
            if set(terminal['methods'])!=set(frozen['methods']):raise ValueError('Incomplete method attempt roster')
            records={}
            for method,metadata in terminal['methods'].items():
                if sha256(metadata['manifest_path'])!=metadata['manifest_sha256']:raise ValueError('Arm manifest changed')
                record=read(metadata['manifest_path']);records[method]=record;all_complete&=record['complete']
            groups[f'{family}:{seed}']={'generation_terminal_path':str(path),'generation_terminal_sha256':sha256(path),'methods':records}
    result={'candidate_freeze_sha256':sha256(FREEZE),'prepared_manifest_sha256':sha256(PREPARED),
            'planned_eligibility_path':str(ELIGIBILITY),'planned_eligibility_sha256':sha256(ELIGIBILITY),
            'truth_manifest_path':str(TRUTH),'truth_manifest_sha256':sha256(TRUTH),'methods':frozen['methods'],
            'balanced_eligibility':frozen['balanced_eligibility'],'groups':groups,'all_methods_complete':bool(all_complete),
            'incomplete_methods_unscored':True,'no_official_scores_or_submission':True}
    write(OUT/'manifest.json',result);write(ROOT/'reports/injection_prediction_manifest.json',result)
    print(json.dumps({'candidate_generation_finalized':True,'all_methods_complete':bool(all_complete),'methods':frozen['methods']}),flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('stage',choices=['freeze','prepare','generate','finalize'])
    p.add_argument('--family',choices=FAMILIES);p.add_argument('--seed',type=int,choices=SEEDS);args=p.parse_args()
    if args.stage=='generate':
        if args.family is None or args.seed is None:p.error('generate requires --family and --seed')
        generate(args.family,args.seed)
    else:{'freeze':freeze,'prepare':prepare,'finalize':finalize}[args.stage]()


if __name__=='__main__':main()
