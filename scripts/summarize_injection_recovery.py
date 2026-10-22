#!/usr/bin/env python3
"""Immutable development then confirmation summaries for the moment-oracle study."""
from pathlib import Path
import argparse
from datetime import datetime,timezone
import json
import anndata as ad
import numpy as np
import pandas as pd
import polars as pl
from run_injection_recovery_de import ROOT,TRUTH,FREEZE,ELIGIBILITY,read,write,verify_freeze
from prepare_injection_recovery import sha,population_moments

BIN_EDGES=(.1,.25,.5)


def stage_partition(seed,template):
    if seed not in (20261026,20261027) or template not in (1,2):raise ValueError('Unknown fixed seed/template')
    if seed==20261026 and template==1:return 'development'
    if seed==20261026:return 'heldout_template'
    if template==1:return 'heldout_seed'
    return 'heldout_template_and_seed'


def safe_ratio(a,b):return float(a/b) if b else None


def side_metrics(truth_lfc,observed_lfc,q,gate,marginal_null):
    truth_lfc=np.asarray(truth_lfc,float);observed_lfc=np.asarray(observed_lfc,float);q=np.asarray(q,float)
    gate=np.asarray(gate,bool);marginal_null=np.asarray(marginal_null,bool)
    if not all(a.shape==truth_lfc.shape for a in (observed_lfc,q,gate,marginal_null)):raise ValueError('Metric gene axes differ')
    significant=gate&(q<.05)
    adjudicable=gate&(np.abs(truth_lfc)>=.1)
    correct=significant&adjudicable&np.isfinite(observed_lfc)&(truth_lfc*observed_lfc>0)
    valid=gate&np.isfinite(observed_lfc)&np.isfinite(truth_lfc)
    null=gate&marginal_null
    result=dict(tested_genes=int(gate.sum()),significant_calls=int(significant.sum()),
        oracle_adjudicable_all_genes=int(np.count_nonzero(np.abs(truth_lfc)>=.1)),
        oracle_adjudicable_gated=int(adjudicable.sum()),oracle_adjudicable_outside_gate=int(np.count_nonzero((np.abs(truth_lfc)>=.1)&~gate)),
        called_adjudicable=int(np.count_nonzero(significant&adjudicable)),called_below_adjudication_threshold=int(np.count_nonzero(significant&~adjudicable)),
        called_with_correct_oracle_sign=int(correct.sum()),signed_recovery=safe_ratio(correct.sum(),adjudicable.sum()),
        direction_accuracy_among_adjudicable_calls=safe_ratio(correct.sum(),np.count_nonzero(significant&adjudicable)),
        marginal_null_gated_genes=int(null.sum()),marginal_null_calls=int(np.count_nonzero(significant&null)),
        marginal_null_call_fraction=safe_ratio(np.count_nonzero(significant&null),null.sum()),
        observed_lfc_mae_to_oracle=float(np.abs(observed_lfc[valid]-truth_lfc[valid]).mean()) if valid.any() else None,
        observed_lfc_rmse_to_oracle=float(np.sqrt(np.mean((observed_lfc[valid]-truth_lfc[valid])**2))) if valid.any() else None)
    bins=[('numerical_zero',0.,1e-10),('tiny',1e-10,.1),('small',.1,.25),('moderate',.25,.5),('large',.5,np.inf)]
    strata=[]
    for name,lo,hi in bins:
        absolute=np.abs(truth_lfc)
        mask=gate&((absolute<=hi) if name=='numerical_zero' else ((absolute>lo)&(absolute<hi)) if name=='tiny' else ((absolute>=lo)&(absolute<hi)))
        n=int(mask.sum());called=significant&mask
        signed=int(np.count_nonzero(called&np.isfinite(observed_lfc)&(truth_lfc*observed_lfc>0)))
        strata.append(dict(stratum=name,oracle_genes=n,significant_calls=int(called.sum()),
            correct_sign_calls=signed if lo>=.1 else None,
            call_recovery=safe_ratio(called.sum(),n),signed_recovery=safe_ratio(signed,n) if lo>=.1 else None))
    result['effect_strata']=strata
    return result


def call_overlap(reference_q,candidate_q,gate):
    a=np.asarray(gate,bool)&(np.asarray(reference_q)<.05);b=np.asarray(gate,bool)&(np.asarray(candidate_q)<.05)
    intersection=int(np.count_nonzero(a&b));union=int(np.count_nonzero(a|b))
    return dict(intersection=intersection,union=union,jaccard=safe_ratio(intersection,union),
        fraction_candidate_calls_in_reference=safe_ratio(intersection,b.sum()),fraction_reference_calls_recovered=safe_ratio(intersection,a.sum()))


def profile_error(p,b,oracle_p,oracle_b):
    p=np.asarray(p,float);b=np.asarray(b,float)
    if np.isnan(p).all() and np.isnan(b).all():return None
    if not np.isfinite(p).all() or not np.isfinite(b).all() or (p<0).any() or (b<0).any():raise ValueError('Invalid available profile stage')
    return dict(arithmetic_probability_l1=float(np.abs(p-oracle_p).sum()),
        arithmetic_cpm_rmse=float(1e6*np.sqrt(np.mean((p-oracle_p)**2))),
        bulk_probability_l1=float(np.abs(b-oracle_b).sum()),
        bulk_log1p50k_rmse=float(np.sqrt(np.mean((np.log1p(50000*b)-np.log1p(50000*oracle_b))**2))))


def detection_changes(reference,candidate,comparison,masks):
    result={}
    for name,mask in masks.items():
        mask=np.asarray(mask,bool);n=int(mask.sum())
        result[name]=dict(genes=n,reference_detection=float(reference[mask].mean()) if n else None,
            candidate_detection=float(candidate[mask].mean()) if n else None,comparison_detection=float(comparison[mask].mean()) if n else None,
            reference_detection_change=float((reference[mask]-comparison[mask]).mean()) if n else None,
            candidate_detection_change=float((candidate[mask]-comparison[mask]).mean()) if n else None,
            candidate_reference_detection_mae=float(np.abs(candidate[mask]-reference[mask]).mean()) if n else None)
    return result


def donor_zero_edits(donor,prediction,masks):
    if donor.shape!=prediction.shape:raise ValueError('Donor/output row axes differ')
    result={'meaning':'Generated edits only; independent reference cells provide no per-cell counterfactual activation truth','groups':{}}
    for name,mask in masks.items():
        a=donor[:,mask];b=prediction[:,mask];zeros=a==0;activated=zeros&(b>0)
        result['groups'][name]=dict(donor_zero_entries=int(zeros.sum()),zero_to_positive_entries=int(activated.sum()),
            activation_fraction_of_donor_zeros=safe_ratio(activated.sum(),zeros.sum()),
            genes_with_any_zero_to_positive=int(np.any(activated,axis=0).sum()),positive_to_zero_entries=int(np.count_nonzero((a>0)&(b==0))))
    return result


def read_de(family,seed,method,aliases):
    folder=TRUTH/'de'/family/f'seed_{seed}'/method;status=read(folder/'status.json')
    if (status['family'],status['seed'],status['method'])!=(family,seed,method):raise ValueError('DE status arm identity differs')
    if status['status']=='incomplete_candidate_not_scored':
        if method=='reference':raise ValueError('Independent reference cannot be an incomplete candidate')
        path=ROOT/f'outputs/injection_recovery_400/{family}/{method}_s{seed}/manifest.json';candidate=read(path)
        if candidate['complete'] or status['candidate_manifest_sha256']!=sha(path):raise ValueError('Invalid incomplete-candidate skip')
        return None,status
    if status['status']!='complete':raise ValueError('DE failed or has unknown status; do not convert failure into an unscored candidate')
    preflight=read(folder/'preflight.json');frozen=read(FREEZE)
    source=read(TRUTH/'manifest.json')['families'][family]['seeds'][str(seed)]
    prediction_sha=source['reference_sha256'] if method=='reference' else read(ROOT/f'outputs/injection_recovery_400/{family}/{method}_s{seed}/manifest.json')['prediction_sha256']
    for key,expected in [('analysis_freeze_sha256',sha(FREEZE)),('eligibility_sha256',sha(ELIGIBILITY)),
        ('truth_manifest_sha256',sha(TRUTH/'manifest.json')),('prediction_sha256',prediction_sha),
        ('reference_sha256',source['reference_sha256']),('comparison_sha256',source['comparison_sha256'])]:
        if preflight[key]!=expected:raise ValueError('DE input provenance differs: '+key)
    if preflight['parameters']!=frozen['de_parameters'] or (preflight['family'],preflight['seed'],preflight['method'])!=(family,seed,method):
        raise ValueError('DE preflight parameters/identity differ')
    if sha(folder/'de.parquet')!=status['de_sha256']:raise ValueError('DE output hash changed')
    table=pl.scan_parquet(folder/'de.parquet').filter(pl.col('target').is_in(list(aliases))).collect()
    if set(table['target'].to_list())!=set(aliases) or table.select(['target','feature']).is_duplicated().sum():
        raise ValueError('DE alias/feature table incomplete or duplicated')
    return table,status


def align_de(table,alias,genes):
    frame=table.filter(pl.col('target')==alias).to_pandas().set_index('feature')
    if frame.index.duplicated().any():raise ValueError('Duplicate DE features')
    aligned=frame.reindex(genes);gate=np.isin(genes,frame.index.to_numpy())
    return aligned['log2_fold_change'].to_numpy(float),aligned['p_adj'].to_numpy(float),gate


def profile_audit(candidate,alias,genes,bank,i,source,prediction):
    path=Path(candidate['gene_profiles_path'])
    if sha(path)!=candidate['gene_profiles_sha256']:raise ValueError('Candidate profile audit hash changed')
    table=pd.read_parquet(path,filters=[('alias','==',alias)]).set_index('gene')
    if table.index.duplicated().any() or set(table.index)!=set(genes):raise ValueError('Candidate profile gene axis differs')
    table=table.reindex(genes)
    def col(name):return table[name].to_numpy(float)
    np.testing.assert_array_equal(col('oracle_p'),bank['oracle_p'][i]);np.testing.assert_array_equal(col('oracle_b'),bank['oracle_b'][i])
    labels=prediction.obs.target_gene.astype(str).to_numpy();rows=np.flatnonzero(labels==alias)
    if len(rows)!=400:raise ValueError('Candidate actual count group is incomplete')
    actual=population_moments(prediction.X[rows])
    np.testing.assert_allclose(actual['p'],col('realized_p'),rtol=0,atol=1e-12)
    np.testing.assert_allclose(actual['b'],col('realized_b'),rtol=0,atol=1e-12)
    oracle_p=bank['oracle_p'][i];oracle_b=bank['oracle_b'][i]
    stages={}
    for name,prefix in [('requested','request'),('population_fitted','population_fitted'),('sample_expected','sample_expected')]:
        stages[name]=profile_error(col(prefix+'_p'),col(prefix+'_b'),oracle_p,oracle_b)
    if 'donor_baseline_p' in table:
        stages['donor_baseline']=profile_error(col('donor_baseline_p'),col('donor_baseline_b'),oracle_p,oracle_b)
    stages['realized']=profile_error(actual['p'],actual['b'],oracle_p,oracle_b)
    dp=col('request_p')-col('inference_p')-bank['delta_p'][i]
    db=col('request_b')-col('inference_b')-bank['delta_b'][i]
    diagnostics=dict(stages=stages,actual_realized_profiles_independently_verified=True,
        clipped_mass_p=float(-np.minimum(col('anchored_before_clip_p'),0).sum()),
        clipped_mass_b=float(-np.minimum(col('anchored_before_clip_b'),0).sum()),
        post_clip_mass_p=float(np.maximum(col('anchored_before_clip_p'),0).sum()),
        post_clip_mass_b=float(np.maximum(col('anchored_before_clip_b'),0).sum()),
        requested_delta_error_p_l1=float(np.abs(dp).sum()),requested_delta_error_b_l1=float(np.abs(db).sum()),
        requested_delta_error_p_max=float(np.abs(dp).max()),requested_delta_error_b_max=float(np.abs(db).max()),
        bulk_projection_l1_to_request=float(np.abs(col('bulk_projected_b')-col('request_b')).sum()),
        conditional_realization_noise_p_l1=float(np.abs(actual['p']-col('sample_expected_p')).sum()),
        conditional_realization_noise_b_l1=float(np.abs(actual['b']-col('sample_expected_b')).sum()),
        source_profile_audit=str(path),source_profile_audit_sha256=candidate['gene_profiles_sha256'])
    if np.isfinite(col('population_fitted_p')).all():
        diagnostics['donor_sampling_p_l1']=float(np.abs(col('sample_expected_p')-col('population_fitted_p')).sum())
        diagnostics['donor_sampling_b_l1']=float(np.abs(col('sample_expected_b')-col('population_fitted_b')).sum())
    else:diagnostics['donor_sampling_p_l1']=diagnostics['donor_sampling_b_l1']=None
    return diagnostics


def summarize(stage):
    frozen=verify_freeze();truth=read(TRUTH/'manifest.json')
    prefix=ROOT/f'reports/injection_recovery_{stage}';out=prefix.with_suffix('.json')
    if any(Path(str(prefix)+suffix).exists() for suffix in ('.json','.md','.csv')):raise FileExistsError('Immutable stage report exists')
    previous=None
    if stage=='confirmation':
        prior=ROOT/'reports/injection_recovery_development.json';previous=read(prior)
        if previous['analysis_freeze_sha256']!=sha(FREEZE) or previous['eligibility_sha256']!=sha(ELIGIBILITY):raise ValueError('Study changed between development and confirmation')
        previous={'path':str(prior),'sha256':sha(prior)}
    reference_rows=[];candidate_rows=[];statuses=[];baseline_offsets=[];full_null_checks=[];input_receipts=[]
    for family,spec in truth['families'].items():
        bank=np.load(spec['profiles_path'],allow_pickle=False);genes=bank['genes'].astype(str)
        for seed_text,source in spec['seeds'].items():
            seed=int(seed_text)
            indices=[i for i,t in enumerate(bank['template_id']) if (stage_partition(seed,int(t))=='development')==(stage=='development')]
            if not indices:continue
            aliases=bank['aliases'][indices].astype(str)
            for key in ('inference','comparison','reference'):
                if sha(source[key+'_path'])!=source[key+'_sha256']:raise ValueError('Truth/control count artifact changed')
            reference,status=read_de(family,seed,'reference',aliases)
            if reference is None:raise ValueError('Independent reference DE must complete before interpretation')
            empirical=np.load(source['empirical_profiles_path'],allow_pickle=False)
            if sha(source['empirical_profiles_path'])!=source['empirical_profiles_sha256']:raise ValueError('Empirical truth-role audit changed')
            baseline_offsets.append(dict(family=family,seed=seed,
                inference=profile_error(empirical['inference_p'],empirical['inference_b'],bank['baseline_p'],bank['baseline_b']),
                comparison=profile_error(empirical['comparison_p'],empirical['comparison_b'],bank['baseline_p'],bank['baseline_b']),
                note='Finite baseline offsets are nuisance, not injected biological effects'))
            candidates={};objects={};selections={}
            reference_cells=ad.read_h5ad(source['reference_path'],backed='r')
            comparison_cells=ad.read_h5ad(source['comparison_path'],backed='r')
            inference_cells=ad.read_h5ad(source['inference_path'],backed='r')
            comparison_detection=np.asarray((comparison_cells.X[:]>0).mean(axis=0)).ravel()
            reference_labels=reference_cells.obs.target_gene.astype(str).to_numpy()
            latent=np.load(truth['population_parameters_path'],allow_pickle=False)
            for method in frozen['planned_methods']:
                table,state=read_de(family,seed,method,aliases);statuses.append(dict(state))
                path=ROOT/f'outputs/injection_recovery_400/{family}/{method}_s{seed}/manifest.json';manifest=read(path)
                input_receipts.append(dict(path=str(path),sha256=sha(path)))
                candidates[method]=(table,manifest)
                if table is not None:
                    if not manifest['complete'] or sha(manifest['prediction_path'])!=manifest['prediction_sha256']:raise ValueError('Candidate completion/hash changed')
                    objects[method]=ad.read_h5ad(manifest['prediction_path'],backed='r')
                    if sha(manifest['selected_rows_path'])!=manifest['selected_rows_sha256']:raise ValueError('Selected donor provenance changed')
                    selections[method]=pd.read_csv(manifest['selected_rows_path'])
            try:
                for i in indices:
                    alias=str(bank['aliases'][i]);partition=stage_partition(seed,int(bank['template_id'][i]))
                    ref_lfc,ref_q,gate=align_de(reference,alias,genes)
                    known=bank['truth_lfc_cpm_eps1e9'][i];marginal_null=bank['marginal_null_gene_mask'][i]
                    reference_detection=np.asarray((reference_cells.X[np.flatnonzero(reference_labels==alias)]>0).mean(axis=0)).ravel()
                    program=np.zeros(len(genes),bool)
                    if family=='dirichlet_multinomial':
                        t=int(bank['template_id'][i]);program[np.r_[latent[f'program{t}_up'],latent[f'program{t}_down']]]=True
                    masks={'gated_all':gate,'oracle_mean_up':gate&(known>=.1),'oracle_mean_down':gate&(known<=-.1),'direct_program':gate&program}
                    rmetrics=side_metrics(known,ref_lfc,ref_q,gate,marginal_null)
                    ref_record=dict(family=family,seed=seed,alias=alias,template=int(bank['template_id'][i]),regime=str(bank['regime'][i]),
                        partition=partition,is_null=bool(bank['is_null'][i]),metrics=rmetrics,
                        reference_sample_profile_error=profile_error(empirical['reference_p'][i],empirical['reference_b'][i],bank['oracle_p'][i],bank['oracle_b'][i]))
                    reference_rows.append(ref_record)
                    if bank['is_null'][i]:full_null_checks.append(dict(family=family,seed=seed,alias=alias,significant_calls=rmetrics['significant_calls'],passed=rmetrics['significant_calls']<=5))
                    for method,(table,manifest) in candidates.items():
                        if table is None:
                            candidate_rows.append(dict(family=family,seed=seed,alias=alias,method=method,partition=partition,status='unscored_incomplete_candidate'))
                            continue
                        lfc,q,predgate=align_de(table,alias,genes)
                        if not np.array_equal(gate,predgate):raise ValueError('Reference and candidate control-expression gates differ')
                        metrics=side_metrics(known,lfc,q,gate,marginal_null)
                        actual_rows=np.flatnonzero(objects[method].obs.target_gene.astype(str).to_numpy()==alias)
                        actual_counts=objects[method].X[actual_rows]
                        candidate_detection=np.asarray((actual_counts>0).mean(axis=0)).ravel()
                        edits=None
                        selected=selections[method];selected=selected[selected.alias==alias].sort_values('output_row')
                        if set(selected.selection_role)=={'preserving_raw_donor'}:
                            if len(selected)!=400 or selected.output_row.tolist()!=list(range(400)):raise ValueError('Raw donor row mapping differs')
                            indices_raw=selected.inference_row.to_numpy();order=np.argsort(indices_raw)
                            sorted_raw=inference_cells.X[indices_raw[order]].toarray();raw=np.empty_like(sorted_raw);raw[order]=sorted_raw
                            edits=donor_zero_edits(raw,actual_counts.toarray(),masks)
                        candidate_rows.append(dict(family=family,seed=seed,alias=alias,method=method,template=int(bank['template_id'][i]),
                            regime=str(bank['regime'][i]),partition=partition,is_null=bool(bank['is_null'][i]),status='complete',
                            metrics=metrics,reference_metrics=rmetrics,call_overlap=call_overlap(ref_q,q,gate),
                            detection=detection_changes(reference_detection,candidate_detection,comparison_detection,masks),raw_donor_zero_edit_diagnostics=edits,
                            profile_audit=profile_audit(manifest,alias,genes,bank,i,source,objects[method])))
            finally:
                for obj in objects.values():obj.file.close()
                reference_cells.file.close();comparison_cells.file.close();inference_cells.file.close()
    null_checks=full_null_checks.copy()
    if stage=='confirmation':null_checks=read(ROOT/'reports/injection_recovery_development.json')['reference_null_checks_opened']+full_null_checks
    null_status=[]
    for family,seed in sorted({(r['family'],r['seed']) for r in null_checks}):
        group=[r for r in null_checks if r['family']==family and r['seed']==seed]
        complete=len(group)==2
        passed=all(r['passed'] for r in group)
        null_status.append(dict(family=family,seed=seed,opened_null_aliases=len(group),complete_family_seed_null_panel=complete,
            status='inconclusive_reference_null' if not passed else 'pass' if complete else 'provisional_opened_null_pass',
            max_significant_calls=max(r['significant_calls'] for r in group)))
    report=dict(created_utc=datetime.now(timezone.utc).isoformat(),stage=stage,study='moment-oracle diagnostic',
        analysis_freeze_sha256=sha(FREEZE),eligibility_sha256=sha(ELIGIBILITY),truth_manifest_sha256=sha(TRUTH/'manifest.json'),
        planned_methods=frozen['planned_methods'],development_report=previous,input_candidate_receipts=input_receipts,
        reference_results=reference_rows,candidate_results=candidate_rows,candidate_statuses=statuses,
        baseline_offsets=baseline_offsets,reference_null_checks_opened=full_null_checks,reference_null_validity=null_status,
        interpretation=['Candidates receive two population response profiles, not full target distributions or cell-state weights.',
        'Incomplete candidates are unscored; no successful-alias subset is selected.',
        'Sign adjudication uses absolute trueLFC>=0.1; tiny effects are not classified as stable.',
        'Composition nonnulls have no general marginal-null set; DM null genes require unchanged conditional law.',
        'Finite comparison/inference offsets, anchoring/clipping, fit, sampling and realization errors are distinct.',
        'A failed independent referenceNULL check makes that family/seed recovery inconclusive; no outcomes are discarded.',
        'No combined success score, biological-transfer conclusion, reserve access or official promotion is produced.'])
    write(out,report)
    flat=[]
    for r in candidate_rows:
        row={k:r[k] for k in ('family','seed','alias','method','partition','status')}
        if r['status']=='complete':
            row.update({k:v for k,v in r['metrics'].items() if k!='effect_strata'})
            row.update(reference_significant_calls=r['reference_metrics']['significant_calls'],reference_signed_recovery=r['reference_metrics']['signed_recovery'],
                call_jaccard=r['call_overlap']['jaccard'],realized_cpm_rmse=r['profile_audit']['stages']['realized']['arithmetic_cpm_rmse'],
                realized_bulk_log_rmse=r['profile_audit']['stages']['realized']['bulk_log1p50k_rmse'],
                requested_delta_error_p_l1=r['profile_audit']['requested_delta_error_p_l1'],requested_delta_error_b_l1=r['profile_audit']['requested_delta_error_b_l1'])
        flat.append(row)
    pd.DataFrame(flat).to_csv(prefix.with_suffix('.csv'),index=False)
    lines=[f'# Moment-oracle {stage} results','',
        'This stage evaluates known response moments against independently generated reference cells. It is not a full-distribution oracle, a biological-transfer test or an official score. All candidate methods and requests were frozen before DE outcomes; no adaptation is permitted between stages.','',
        '| Family | Seed | Alias | Method | Reference calls | Candidate calls | Reference signed recovery | Candidate signed recovery | NULL/marginal-null calls |',
        '| --- | ---: | --- | --- | ---: | ---: | ---: | ---: | ---: |']
    def fmt(x):return 'undefined' if x is None else f'{x:.3f}'
    for r in candidate_rows:
        if r['status']!='complete':
            lines.append(f"| {r['family']} | {r['seed']} | {r['alias']} | {r['method']} | — | unscored | — | — | — |");continue
        m=r['metrics'];ref=r['reference_metrics']
        lines.append(f"| {r['family']} | {r['seed']} | {r['alias']} | {r['method']} | {ref['significant_calls']} | {m['significant_calls']} | {fmt(ref['signed_recovery'])} | {fmt(m['signed_recovery'])} | {m['marginal_null_calls']} |")
    lines+=['','Signed recovery uses oracle |log2FC|≥0.1 within the realized comparison-CPM>5 gate. Below-threshold calls are reported without calling all of them false. NULL/marginal-null calls use the independent true-law mask, not merely a small aggregate mean. Undefined denominators remain undefined.','',
        'Reference NULL validity: '+ '; '.join(f"{x['family']}/{x['seed']}: {x['status']} (maximum{x['max_significant_calls']}calls)" for x in null_status)+'.','',
        'The detailed JSON/CSV retain effect strata, out-of-gate true signals, reference call power, overlap, finite-baseline nuisance, per-gene request/clipping provenance and population/sample/realized profile errors. Current-emitter population-fit statistics are unavailable rather than fabricated. Candidate profile values were checked against actual integer counts.','',
        'A weak result can reflect anchoring distortion, finite-control uncertainty, missing higher-order distribution information or generator-family limitations. It is not automatically an implementation bug. No scalar success metric or biological promotion decision is inferred.']
    if stage=='development':lines+=['','This immutable report opens only template1/seed20261026. Full family/seed reference-null validity remains pending confirmation; held-out templates/seeds are not used for a development decision.']
    else:lines+=['','Confirmation separately reports heldout-template, heldout-seed and joint-heldout rows. It does not repair a failed development or reference-null gate. The immutable development report hash is retained in JSON.']
    prefix.with_suffix('.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({'stage':stage,'written':str(out),'candidate_rows':len(candidate_rows)}),flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--stage',choices=['development','confirmation'],required=True)
    summarize(p.parse_args().stage)


if __name__=='__main__':main()
