#!/usr/bin/env python3
"""One preregistered response-balanced hurdle objective on the frozen P1 folds.

Magnitude-proxy probabilities only. No H1, DE significance, threshold search,
minimum call quota, alternate objective, scoring or submission.
"""
from __future__ import annotations
import os
for _key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ.setdefault(_key,'1')
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd

import response_transfer_adapter as adapter
import benchmark_matched_response as p1

ROOT=p1.ROOT
PREFIX=ROOT/'reports/matched_response_p2'
MANIFEST=PREFIX.with_name(PREFIX.name+'_manifest.json')
PROTOCOL=ROOT/'docs/MATCHED_RESPONSE_P2_PROTOCOL_20261021.md'
POLICY={'ridge':.03,'calibration_fraction':.2,'partition_seed':20261023,
        'probability_rule':'p > weighted recalibration prevalence + 1e-8',
        'decision_epsilon':1e-8,'sign_threshold':.5,'clip':3.,
        'response_class_mass':[.5,.5],'calibration_class_balanced':False,
        'calibration_transductive_common_features':True,'logistic_iterations':100,
        'diagnostic_control_cpm_threshold':5.,
        'probability_semantics':'magnitude proxy only; class-balanced decision costs, not DE/deployment probability'}


class DegenerateClasses(ValueError):
    pass


def split_targets(targets):
    targets=sorted(set(str(t) for t in targets))
    if len(targets)<2: raise DegenerateClasses('At least two training targets required')
    ordered=sorted(targets,key=lambda t:hashlib.sha256(f'{POLICY["partition_seed"]}:{t}'.encode()).hexdigest())
    n=max(1,min(len(targets)-1,int(round(len(targets)*POLICY['calibration_fraction']))))
    return sorted(set(targets)-set(ordered[:n])),sorted(ordered[:n])


def features(task):
    # The original basis has nine modifiers per pooled component.
    shared=adapter.FEATURE_NAMES.index('shared:constant')
    unsigned=np.column_stack([task['confidence_x'],np.log1p(np.abs(task['x'][:,0])),
                              np.log1p(np.abs(task['x'][:,shared]))])
    signed=np.column_stack([np.ones(len(task['x'])),task['x']])
    return unsigned,signed


def combine(tasks, selected_targets):
    selected=set(selected_targets); output={k:[] for k in ('unsigned','signed','y','weight','targets')}
    context_records=[]
    for task in tasks:
        mask=np.isin(task['targets'],list(selected))
        if not mask.any(): raise DegenerateClasses('A source context has no targets in this fitting/calibration subset')
        u,s=features(task); weight=task['weight'][mask].copy()
        original_mass=float(weight.sum());weight/=original_mass*len(tasks)
        for key,value in [('unsigned',u[mask]),('signed',s[mask]),('y',task['y'][mask]),
                          ('weight',weight),('targets',task['targets'][mask])]: output[key].append(value)
        context_records.append({'context':task['destination'],'targets':len(set(task['targets'][mask])),
                                'pairs':int(mask.sum()),'original_weight_mass':original_mass,
                                'renormalized_weight_mass':float(weight.sum())})
    combined={k:np.concatenate(v) for k,v in output.items()}
    combined['contexts']=context_records
    return combined


def sigmoid(x):
    return 1/(1+np.exp(-np.clip(x,-60,60)))


def scaler(x,weight):
    weight=np.asarray(weight,dtype=float);weight=weight/weight.sum()
    mean=(weight[:,None]*x).sum(0)
    scale=np.sqrt((weight[:,None]*(x-mean)**2).sum(0));scale=np.maximum(scale,1e-8)
    mean[0]=0;scale[0]=1
    return mean,scale


def fit_logistic(x,y,weight,*,scaler_weight=None):
    x=np.asarray(x,dtype=float);y=np.asarray(y,dtype=float);w=np.asarray(weight,dtype=float)
    if set(np.unique(y))!={0.,1.}: raise DegenerateClasses('Both logistic response/sign classes are required')
    if not np.isfinite(x).all() or not np.isfinite(w).all() or (w<=0).any(): raise ValueError('Invalid logistic inputs/weights')
    w=w/w.sum();mean,scale=scaler(x,w if scaler_weight is None else scaler_weight)
    z=(x-mean)/scale; penalty=np.full(z.shape[1],POLICY['ridge']);penalty[0]=0
    beta=np.zeros(z.shape[1]);prior=float(w@y);beta[0]=np.log(prior/(1-prior))
    def objective(b):
        eta=z@b
        return float(w@(np.logaddexp(0,eta)-y*eta)+.5*np.sum(penalty*b*b))
    converged=False
    for iteration in range(POLICY['logistic_iterations']):
        probability=sigmoid(z@beta)
        gradient=z.T@(w*(probability-y))+penalty*beta
        if np.max(np.abs(gradient))<1e-8:
            converged=True;break
        hessian=(z*(w*probability*(1-probability))[:,None]).T@z+np.diag(penalty)
        step=np.linalg.solve(hessian,gradient); old=objective(beta); rate=1.
        while rate>2**-30:
            candidate=beta-rate*step
            if objective(candidate)<=old+1e-12: break
            rate/=2
        if rate<=2**-30: raise ValueError('Fixed logistic solver failed to descend')
        beta=candidate
    if not converged:
        gradient=z.T@(w*(sigmoid(z@beta)-y))+penalty*beta
        converged=np.max(np.abs(gradient))<1e-7
    if not converged: raise ValueError('Fixed logistic solver did not converge')
    return {'coefficient':beta.tolist(),'mean':mean.tolist(),'scale':scale.tolist(),
            'iterations':iteration+1,'ridge':POLICY['ridge'],'intercept_penalized':False}


def score(model,x):
    return ((np.asarray(x)-np.asarray(model['mean']))/np.asarray(model['scale']))@np.asarray(model['coefficient'])


def fit_linear(x,y,weight):
    w=np.asarray(weight,dtype=float).copy();w/=w.sum();mean,scale=scaler(x,w);z=(x-mean)/scale
    penalty=np.full(z.shape[1],POLICY['ridge']);penalty[0]=0
    beta=np.linalg.solve((z*w[:,None]).T@z+np.diag(penalty),(z*w[:,None]).T@y)
    return {'coefficient':beta.tolist(),'mean':mean.tolist(),'scale':scale.tolist(),'ridge':POLICY['ridge']}


def fit_platt(logits,y,weight):
    w=np.asarray(weight,dtype=float);w=w/w.sum();prevalence=float(w@y)
    if not 0<prevalence<1: raise DegenerateClasses('Both proxy classes are required in recalibration')
    fitted=fit_logistic(np.column_stack([np.ones(len(logits)),logits]),y,w)
    slope=fitted['coefficient'][1]/fitted['scale'][1]
    intercept=fitted['coefficient'][0]-slope*fitted['mean'][1]
    if slope<=1e-12:
        slope=0.;intercept=float(np.log(prevalence/(1-prevalence)))
    return {'slope':float(slope),'intercept':float(intercept),'prevalence':prevalence,
            'constant':bool(slope==0),'ridge':POLICY['ridge'],'intercept_penalized':False,
            'class_balanced_weights':False,'standardized_slope_fit':fitted}


def calibrated_probability(calibrator,logits):
    if calibrator['constant']: return np.full(len(logits),calibrator['prevalence'])
    return sigmoid(calibrator['intercept']+calibrator['slope']*logits)


def reliability(y,probability,weight):
    w=np.asarray(weight,dtype=float);w=w/w.sum();p=np.asarray(probability);y=np.asarray(y)
    bins=np.minimum((p*10).astype(int),9);records=[]
    for i in range(10):
        selected=bins==i;mass=float(w[selected].sum())
        records.append({'bin':i,'lower':i/10,'upper':(i+1)/10,'pairs':int(selected.sum()),'weight_mass':mass,
            'mean_probability':float(w[selected]@p[selected]/mass) if mass else None,
            'observed_proxy_frequency':float(w[selected]@y[selected]/mass) if mass else None})
    return {'brier':float(w@((p-y)**2)),'weighted_prevalence':float(w@y),
            'mean_probability':float(w@p),'bins':records,'label':'magnitude proxy, not DE significance'}


def fit_hurdle(tasks,calibration_targets,space):
    all_targets=set(t for task in tasks for t in task['targets'])
    calibration_targets=set(calibration_targets)
    if not calibration_targets<=all_targets: raise ValueError('Unknown recalibration target')
    fitting_targets=sorted(all_targets-calibration_targets)
    fit=combine(tasks,fitting_targets);cal=combine(tasks,sorted(calibration_targets))
    if set(fit['targets'])&set(cal['targets']): raise ValueError('Direct fitting/recalibration labels overlap')
    threshold=p1.THRESHOLDS[space]
    responsive=np.abs(fit['y'])>=threshold;labels=responsive.astype(float)
    prevalence=float(fit['weight']@labels)
    if not 0<prevalence<1: raise DegenerateClasses('Both fitting magnitude-proxy classes are required')
    balanced=fit['weight']*np.where(responsive,.5/prevalence,.5/(1-prevalence))
    response=fit_logistic(fit['unsigned'],labels,balanced,scaler_weight=fit['weight'])
    conditional_weight=fit['weight'][responsive].copy();conditional_weight/=conditional_weight.sum()
    sign=fit_logistic(fit['signed'][responsive],(fit['y'][responsive]>0).astype(float),conditional_weight)
    magnitude=fit_linear(fit['unsigned'][responsive],np.log(np.abs(fit['y'][responsive])/threshold),conditional_weight)
    cal_labels=(np.abs(cal['y'])>=threshold).astype(float)
    logits=score(response,cal['unsigned']);calibrator=fit_platt(logits,cal_labels,cal['weight'])
    probability=calibrated_probability(calibrator,logits)
    return {'space':space,'threshold':threshold,'response':response,'sign':sign,'magnitude':magnitude,
        'calibrator':calibrator,'fitting_prevalence':prevalence,
        'fitting_target_ids':fitting_targets,'recalibration_target_ids':sorted(calibration_targets),
        'fitting_contexts':fit['contexts'],'recalibration_contexts':cal['contexts'],
        'balanced_class_weight_mass':[float(balanced[~responsive].sum()),float(balanced[responsive].sum())],
        'recalibration_subset_reliability':reliability(cal_labels,probability,cal['weight']),
        'reliability_caution':'Recalibration fit subset with transductive full-source common preprocessing; not independent unseen-target validation.'}


def decide(probability,prevalence,sign_probability,log_excess,coverage,threshold):
    calls=(np.asarray(probability)>prevalence+POLICY['decision_epsilon'])&(np.asarray(coverage)>0)
    magnitude=np.minimum(POLICY['clip'],threshold*np.exp(np.clip(log_excess,0,np.log(POLICY['clip']/threshold))))
    direction=np.where(np.asarray(sign_probability)>=POLICY['sign_threshold'],1.,-1.)
    return np.where(calls,direction*magnitude,0.),calls


def predict_hurdle(model,task):
    unsigned,signed=features(task)
    probability=calibrated_probability(model['calibrator'],score(model['response'],unsigned))
    probability[task['coverage']==0]=0
    sign_probability=sigmoid(score(model['sign'],signed))
    prediction,calls=decide(probability,model['calibrator']['prevalence'],sign_probability,
        score(model['magnitude'],unsigned),task['coverage'],model['threshold'])
    return prediction,probability,sign_probability,calls


def control_expression_strata(destination,task,fusion,prediction,threshold):
    """Descriptive strata retain original pair weights; never enter the gate."""
    ri=np.asarray([destination._rows[str(t)] for t in task['targets']])
    ci=np.asarray([destination._cols[str(g)] for g in task['genes']])
    expression=destination.control_mean_cpm[ri,ci]
    masks={'control_CPM_gt5':expression>5.,'control_CPM_le5':expression<=5.}
    result={}
    for name,mask in masks.items():
        record={'pairs':int(mask.sum()),'original_weight_mass':float(task['weight'][mask].sum()),
                'targets':sorted(set(task['targets'][mask])),'methods':{},'selection_or_advancement_gate':False}
        if mask.any():
            axis,inverse=np.unique(task['targets'][mask],return_inverse=True)
            mass=np.bincount(inverse,weights=task['weight'][mask])
            for method,values in [('matched_fusion',fusion),('response_hurdle',prediction)]:
                checked,primitive=p1.target_primitives(task['targets'][mask],task['y'][mask],values[mask],
                    task['coverage'][mask],threshold)
                np.testing.assert_array_equal(axis,checked)
                record['methods'][method]=p1.summarize_primitives(np.average(primitive,axis=0,weights=mass)[None,:])
        result[name]=record
    return result


def training_tasks(banks,fold,pairs,space):
    lookup={b.family:b for b in banks}; heldout=fold['heldout']
    if any(heldout in [t['destination'],*t['inputs']] for t in fold['training_tasks']):
        raise ValueError('Outer held-out context leaked into coefficient fitting or recalibration')
    return [p1.task_from_pairs(lookup[t['destination']],[lookup[f] for f in t['inputs']],space,p1.get_pairs(pairs,t['pairs_key']))
            for t in fold['training_tasks']]


def validate_p1_files(original):
    for path,digest in original['files'].items():
        if adapter.sha256(ROOT/path)!=digest: raise ValueError('Original P1 artifact changed: '+path)
    if adapter.sha256(p1.PAIRS)!=original['pairs_sha256']: raise ValueError('Original P1 pairs changed')


def freeze():
    if MANIFEST.exists(): raise FileExistsError(MANIFEST)
    p1_result=p1.read(ROOT/'reports/matched_response_results.json')
    if not p1_result['advancement']['under_recovery_established']: raise ValueError('Conditional P2 trigger did not fire')
    original=p1.read(p1.MANIFEST)
    validate_p1_files(original)
    with np.load(p1.PAIRS,allow_pickle=False) as d: pairs={k:d[k] for k in d.files}
    splits={}
    for fold in original['folds']:
        targets=set(t for task in fold['training_tasks'] for t in p1.get_pairs(pairs,task['pairs_key'])[0])
        fit,cal=split_targets(targets);splits[fold['heldout']]={'fitting_targets':fit,'recalibration_targets':cal}
    files=['scripts/benchmark_response_hurdle.py','tests/test_response_hurdle.py',str(PROTOCOL.relative_to(ROOT)),
           'scripts/benchmark_matched_response.py','scripts/response_transfer_adapter.py',
           'reports/matched_response_manifest.json','reports/matched_response_pairs.npz',
           'reports/matched_response_results.json','reports/matched_response_predictions.npz',
           'reports/matched_response_source_audit.json']
    record={'frozen_at_utc':datetime.now(timezone.utc).isoformat(),'policy':POLICY,'splits':splits,
        'files':{f:adapter.sha256(ROOT/f) for f in files},'p2_outcomes_read_before_freeze':False,
        'approved_scope':'One response-balanced hurdle procedure; no fallback, grid, new H1 cycle or submission.',
        'advancement':'Exact unchanged P1 signed precision/recovery/unweighted-error gate; stable false calls diagnostic only.'}
    MANIFEST.write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps({'frozen':str(MANIFEST),'targets':{f:{k:len(v) for k,v in s.items()} for f,s in splits.items()}}),flush=True)


def inverse_diagnostics(destination,inputs,models,targets):
    genes=destination.genes[destination.measured];ci=np.flatnonzero(destination.measured);rows=[]
    for start in range(0,len(targets),4):
        chunk=np.asarray(targets[start:start+4]);ri=np.asarray([destination._rows[t] for t in chunk])
        control_mean=destination.control_mean_cpm[np.ix_(ri,ci)];control_bulk=destination.control_probability[np.ix_(ri,ci)]
        effect={}
        for space in adapter.SPACES:
            t=np.repeat(chunk,len(genes));g=np.tile(genes,len(chunk))
            x,cx,n,_=adapter.build_features(inputs,t,g,control_mean.ravel(),space,p1.POLICY['input_prior_counts'],20)
            prediction,_,_,_=predict_hurdle(models[space],{'x':x,'confidence_x':cx,'coverage':n})
            effect[space]=(prediction.reshape(len(chunk),len(genes)),n.reshape(len(chunk),len(genes)))
        mean,bulk=adapter.effects_to_profiles(effect,control_mean,control_bulk,clip=3.)
        actual_mean=destination.mean_cpm[np.ix_(ri,ci)];actual_mean=actual_mean/actual_mean.sum(1,keepdims=True)
        actual_bulk=destination.counts[np.ix_(ri,ci)].astype(float);actual_bulk/=actual_bulk.sum(1,keepdims=True)
        for i,target in enumerate(chunk):
            mask=genes!=target
            rows.append({'context':destination.family,'target':str(target),
                'arithmetic_composition_log2_mse':float(np.mean((np.log2(1+1e6*mean[i,mask])-np.log2(1+1e6*actual_mean[i,mask]))**2)),
                'bulk_log1p_mse':float(np.mean((np.log1p(50000*bulk[i,mask])-np.log1p(50000*actual_bulk[i,mask]))**2))})
    return rows


def render(result):
    lines=['# One response-balanced hurdle objective: P2','',
        '**Magnitude proxies only; no DE-significance or deployment-probability claim.** '
        'The fixed prevalence operating point expresses class-balanced costs, not probability >0.5.','',
        '| Context | Head | Signed precision | Signed recovery | Call fraction | Unweighted MSE | Proxy Brier |',
        '|---|---|---:|---:|---:|---:|---:|']
    def fmt(x): return 'undefined' if x is None else f'{x:.6f}'
    for context,heads in result['contexts'].items():
        for space,data in heads.items():
            if data['status']!='COMPLETE':
                lines += [f'| {context} | {space} | failed: {data["reason"]} | | | | |'];continue
            metric=data['response_hurdle']
            lines.append('| '+' | '.join([context,space,*[fmt(metric[k]) for k in
                ('signed_precision','signed_recovery','predicted_call_fraction','unweighted_mse')],
                fmt(data['outer_context_reliability']['brier'])])+' |')
    lines += ['',f'Both-head advancement against matched fusion: **{result["advancement"]["both_heads_pass"]}**.', '',
        'The JSON records fitting/recalibration weights and target IDs, original prevalences, monotone Platt parameters, '
        'proxy reliability bins, conditional sign/magnitude models, fixed decisions, paired uncertainty and failures. '
        'Recalibration uses a target-label subset within training contexts with full-source common preprocessing; '
        'it is not independent unseen-target calibration. The outer held-out context is absent from every training/calibration input and label.', '',
        'The exact [P1 fusion gate](../docs/MATCHED_RESPONSE_BENCHMARK_PROTOCOL_20261021.md) is unchanged; '
        'stable false-call rate and matched-control CPM>5/<=5 strata remain diagnostic. [P2 protocol](../docs/MATCHED_RESPONSE_P2_PROTOCOL_20261021.md) '
        'was fixed before fitting. No method or threshold rescue is performed.']
    if not result['advancement']['both_heads_pass']:
        lines += ['', '**Stop the next H1 scoring cycle.** This one objective did not establish sufficient source-only response '
                  'evidence. Do not promote, extract reserved responses, or submit this candidate.']
    return '\n'.join(lines)+'\n'


def run():
    output=PREFIX.with_name(PREFIX.name+'_results.json')
    if output.exists(): raise FileExistsError(output)
    frozen=p1.read(MANIFEST)
    if frozen['policy']!=POLICY: raise ValueError('P2 policy changed')
    for path,digest in frozen['files'].items():
        if adapter.sha256(ROOT/path)!=digest: raise ValueError('P2 frozen artifact changed: '+path)
    original=p1.read(p1.MANIFEST);validate_p1_files(original)
    banks=p1.banks_from_frozen_model();lookup={b.family:b for b in banks}
    expected_sources={s['family']:s['sha256'] for s in original['sources']}
    if {b.family:b.digest for b in banks}!=expected_sources: raise ValueError('Original P1 source identities changed')
    with np.load(p1.PAIRS,allow_pickle=False) as d: pairs={k:d[k] for k in d.files}
    with np.load(ROOT/'reports/matched_response_predictions.npz',allow_pickle=False) as d: prior={k:d[k] for k in d.files}
    previous=p1.read(ROOT/'reports/matched_response_results.json')
    contexts={};gate_contexts=[];saved={};fold_models={};flat=[];inverse=[];complete=True
    for fold in original['folds']:
        family=fold['heldout'];destination=lookup[family];inputs=[lookup[f] for f in fold['test_inputs']]
        contexts[family]={};gate_heads={};models={}
        for space in adapter.SPACES:
            try:
                tasks=training_tasks(banks,fold,pairs,space)
                model=fit_hurdle(tasks,frozen['splits'][family]['recalibration_targets'],space)
            except DegenerateClasses as error:
                contexts[family][space]={'status':'DEGENERATE_CLASSES','reason':str(error),'no_fallback_or_artificial_calls':True}
                complete=False;continue
            models[space]=model
            task=p1.task_from_pairs(destination,inputs,space,p1.get_pairs(pairs,fold['eval_key']))
            np.testing.assert_array_equal(task['y'],prior[f'{family}_{space}_observed'])
            prediction,probability,sign_probability,calls=predict_hurdle(model,task)
            target_axis,primitive=p1.target_primitives(task['targets'],task['y'],prediction,task['coverage'],p1.THRESHOLDS[space])
            np.testing.assert_array_equal(target_axis,prior[f'{family}_{space}_primitive_targets'])
            metric=p1.summarize_primitives(primitive);base=previous['contexts'][family][space]['matched_fusion']
            rel=reliability((np.abs(task['y'])>=p1.THRESHOLDS[space]).astype(float),probability,task['weight'])
            intervals=p1.paired_intervals(prior[f'{family}_{space}_matched_fusion_target_primitives'],primitive,
                                          p1.POLICY['bootstrap_draws'],p1.POLICY['bootstrap_seed'])
            contexts[family][space]={'status':'COMPLETE','matched_fusion':base,'response_hurdle':metric,
                'outer_context_reliability':rel,'candidate_minus_fusion_intervals':intervals,
                'call_count':int(calls.sum()),'model':model,
                'control_expression_strata':control_expression_strata(destination,task,
                    prior[f'{family}_{space}_matched_fusion'],prediction,p1.THRESHOLDS[space])}
            gate_heads[space]={'matched_fusion':base,'adapter_v1':metric}
            flat.append({'context':family,'head':space,**metric,'proxy_brier':rel['brier']})
            for name,value in [('prediction',prediction),('probability',probability),('sign_probability',sign_probability),
                               ('calls',calls),('target_primitives',primitive)]: saved[f'{family}_{space}_{name}']=value
        if len(models)==2:
            inverse.extend(inverse_diagnostics(destination,inputs,models,fold['eval_targets']))
        gate_contexts.append(gate_heads);fold_models[family]=models
        print(json.dumps({'completed_context':family}),flush=True)
    gate=p1.advancement(gate_contexts) if complete else {'both_heads_pass':False,'reason':'At least one registered fold/head has degenerate classes; no fallback or rescue.'}
    gate.pop('conditional_p2_design_allowed',None)
    if 'under_recovery_established' in gate:
        gate['under_recovery_diagnostic_only']=gate.pop('under_recovery_established')
    gate['additional_response_objective_authorized']=False
    gate['one_prespecified_p2_objective_consumed']=True
    pred_path=PREFIX.with_name(PREFIX.name+'_predictions.npz');np.savez_compressed(pred_path,**saved)
    model_path=PREFIX.with_name(PREFIX.name+'_fold_models.json');model_path.write_text(json.dumps(fold_models,indent=2)+'\n')
    inverse_path=PREFIX.with_name(PREFIX.name+'_inverse_diagnostics.csv');pd.DataFrame(inverse).to_csv(inverse_path,index=False)
    result={'analysis_sha256':adapter.sha256(__file__),'manifest_sha256':adapter.sha256(MANIFEST),
        'policy':POLICY,'contexts':contexts,'advancement':gate,'all_fold_heads_completed':complete,
        'predictions_sha256':adapter.sha256(pred_path),'models_sha256':adapter.sha256(model_path),
        'inverse_diagnostics_sha256':adapter.sha256(inverse_path),'h1_responses_used':False,'DE_significance_labels_used':False,
        'candidate_label':'response_hurdle (gate helper field adapter_v1 is only a schema alias)',
        'next_action':'Source-only proxy gate passes; further registered upstream gates still required before any H1 work or promotion.' if gate['both_heads_pass']
                      else 'Stop the next H1 scoring cycle; no additional objective, threshold rescue, reserve extraction or submission.'}
    output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');output.with_suffix('.md').write_text(render(result))
    pd.DataFrame(flat).to_csv(PREFIX.with_name(PREFIX.name+'_metrics.csv'),index=False)
    print(json.dumps(gate,indent=2),flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=['freeze','run']);args=p.parse_args()
    freeze() if args.action=='freeze' else run()


if __name__=='__main__': main()
