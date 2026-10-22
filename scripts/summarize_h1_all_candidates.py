#!/usr/bin/env python3
"""Score assigned raw/HCT arms and compare frozen H1 Validation candidates.

Only the eight RAW_ARMS score files are written. Initial PIE/neutral directories
are read-only inputs owned by the initial-arm evaluator. No Test or K562
confirmation outcomes, prediction generation, or submission operations occur.
"""
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:
    os.environ[key]='1'
import argparse
from datetime import datetime,timezone
import hashlib,json,subprocess
from pathlib import Path
import numpy as np
import pandas as pd
import yaml
from cell_eval2.real_bundle import manifest_digest

ROOT=Path(__file__).resolve().parents[1]
EVAL=ROOT/'evaluations/h1_validation'
BUNDLE=EVAL/'real_bundle_cpu'
CLI=ROOT / '.venv/bin/cell-eval2'
METRICS={
    'pds_cosine':'PDS','expr_mse_unbiased_capped_norm':'MSE',
    'de_wilcoxon_lfc_nmae':'NMAE','de_wilcoxon_direction_fidelity_yield_raw':'Fidelity',
    'de_wilcoxon_direction_reach_raw':'Reach','de_wilcoxon_sig_jaccard':'Jaccard',
}
INITIAL_ARMS=['neutral_cpu','pie_dual_cpu','pie_single_cpu','pie_dual_pool1_cpu']
RAW_ARMS=['hct_raw_dual_cpu','raw_dual_cpu','raw_dual_half_cpu','raw_single_cpu',
          'raw_dual_agreement_cpu','raw_dual_hct2_cpu','raw_dual_cd4_cpu','raw_dual_full_cpu']
PEERS=['cell_eval2_version','config_digest','comparator','source_fingerprint','source_fingerprint_strict',
       'resolved_device','resolved_de_backend','input_type_real_effective','de_real_fingerprint']


def sha(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda:handle.read(8*1024**2),b''):digest.update(block)
    return digest.hexdigest()


def relative(path):
    try:return str(Path(path).resolve().relative_to(ROOT))
    except ValueError:return str(path)


def load_bundle():
    manifest=json.loads((BUNDLE/'manifest.json').read_text())
    if (manifest.get('n_splits')!=5 or manifest.get('resolved_device')!='cpu'
            or manifest.get('resolved_de_backend')!='scanpy' or manifest.get('profile')!='vcc2026'
            or not manifest.get('rule_digest') or manifest.get('rule_mismatches')
            or set(manifest.get('members',[]))!=set(METRICS)):
        raise ValueError('Require a fresh CPU/scanpy, competition-enrolled, five-split vcc2026 bundle')
    if len(manifest.get('derived_seeds',[]))!=5:
        raise ValueError('Bundle does not record all five replicate splits')
    return manifest


def validate_run(directory,manifest):
    for filename in ['agg_results.csv','results.csv','run_meta.json','run_params.yaml','local_scoring_input.receipt.json']:
        if not (directory/filename).exists():raise FileNotFoundError(directory/filename)
    meta=json.loads((directory/'run_meta.json').read_text())
    params=yaml.safe_load((directory/'run_params.yaml').read_text())
    mismatches={k:(manifest.get(k),meta.get(k)) for k in PEERS if manifest.get(k)!=meta.get(k)}
    if mismatches:raise ValueError(f'{directory.name} bundle identity mismatch: {mismatches}')
    if (params.get('device')!='cpu' or params.get('de',{}).get('backend')!='scanpy'
            or params.get('control_source')!='real' or params.get('metrics')!='vcc2026'
            or meta.get('de_pred_fingerprint') is not None):
        raise ValueError(f'{directory.name} is not an eligible measured-count CPU run')
    return meta,params


def read_scores(directory,manifest):
    frame=pd.read_csv(directory/'official_local_scores.csv')
    if frame.metric.duplicated().any():raise ValueError(f'Duplicate score rows: {directory}')
    indexed=frame.set_index('metric')
    values={m:float(indexed.loc[m,'from_replicate']) for m in METRICS}
    overall=float(indexed.loc['avg_score','from_replicate'])
    if not np.isfinite(list(values.values())+[overall]).all():raise ValueError('Nonfinite official score')
    np.testing.assert_allclose(overall,np.mean(list(values.values())),rtol=0,atol=1e-12)
    if pd.notna(indexed.loc['avg_score','from_baseline']):
        raise ValueError('Diagnostic baseline-relative avg_score is not a five-split competition score')
    if set(frame.anchor_source.dropna())!={'real_bundle'}:
        raise ValueError('Score was not produced through --real-bundle')
    if set(frame.real_bundle_digest.dropna())!={manifest_digest(manifest)}:
        raise ValueError('Score CSV uses a different reference bundle')
    return overall,values


def score_raw(manifest):
    # The full-source candidate is a required comparison; do not score a partial queue.
    for name in RAW_ARMS:validate_run(EVAL/name,manifest)
    bundle_hashes={p.name:sha(p) for p in BUNDLE.iterdir() if p.is_file()}
    for name in RAW_ARMS:
        directory=EVAL/name
        fingerprint={'aggregate_sha256':sha(directory/'agg_results.csv'),
                     'run_meta_sha256':sha(directory/'run_meta.json'),'bundle_artifact_sha256s':bundle_hashes}
        receipt_path=directory/'official_local_score_receipt.json'
        old=json.loads(receipt_path.read_text()) if receipt_path.exists() else {}
        scores_path=directory/'official_local_scores.csv'
        reusable=(scores_path.exists() and all(old.get(k)==v for k,v in fingerprint.items())
                  and old.get('score_csv_sha256')==sha(scores_path))
        if not reusable:
            command=[str(CLI),'score','--user-agg',str(directory/'agg_results.csv'),
                     '--real-bundle',str(BUNDLE),'-o',str(scores_path)]
            completed=subprocess.run(command,cwd=ROOT,capture_output=True,text=True)
            (directory/'official_local_score.log').write_text(completed.stdout+completed.stderr)
            if completed.returncode:raise RuntimeError(f'Official score failed: {name}: {completed.stderr[-1000:]}')
            overall,_=read_scores(directory,manifest)
            receipt_path.write_text(json.dumps({**fingerprint,'command':command,'exit_code':0,
                'scored_utc':datetime.now(timezone.utc).isoformat(),'score_csv_sha256':sha(scores_path),
                'n_anchor_splits':5,'six_component_mean_verified':True,'overall':overall},indent=2)+'\n')
        overall,_=read_scores(directory,manifest)
        print(name,overall,flush=True)


def collect_arm(name,manifest,panel,reference_hash):
    directory=EVAL/name;meta,params=validate_run(directory,manifest)
    overall,scaled=read_scores(directory,manifest)
    aggregate=pd.read_csv(directory/'agg_results.csv').set_index('statistic')
    raw={m:float(aggregate.loc['mean',m]) for m in METRICS}
    if not np.isfinite(list(raw.values())).all():raise ValueError('Nonfinite aggregate metric')
    link_path=directory/'local_scoring_input.receipt.json';link=json.loads(link_path.read_text())
    prediction=Path(link['prediction']);receipt_path=prediction.with_suffix('.receipt.json')
    receipt=json.loads(receipt_path.read_text());args=receipt.get('arguments',{})
    if link['reference_sha256']!=reference_hash or receipt['output_sha256']!=link['prediction_sha256']:
        raise ValueError(f'{name} reference/prediction receipt chain mismatch')
    targets=list(receipt.get('targets',panel))
    if set(targets)!=set(panel):raise ValueError(f'{name} target panel differs')
    sources=receipt.get('sources',[])
    if any('H1' in Path(s['path']).name.upper() for s in sources):
        raise ValueError(f'{name} contains a forbidden H1 response source')
    uncovered=receipt.get('uncovered_targets',targets if not sources else None)
    if uncovered is None:raise ValueError(f'{name} has no target coverage audit')
    coverage=len(targets)-len(uncovered)
    if 'same_target_coverage' in receipt and receipt['same_target_coverage']!=coverage:
        raise ValueError('Coverage audit count disagrees with uncovered-target list')
    per_target=pd.read_csv(directory/'results.csv')
    finite={m:int(np.isfinite(per_target.loc[per_target.metric==m,'value']).sum()) for m in METRICS}
    moments=receipt.get('moments',[])
    return dict(overall=overall,scaled=scaled,raw=raw,raw_finite_target_count=finite,
        coverage=dict(covered=coverage,total=len(panel),uncovered_targets=uncovered),
        source_configuration=dict(arguments=args,sources=sources,agreement=receipt.get('agreement')),
        artifact_paths=dict(run_directory=relative(directory),score_csv=relative(directory/'official_local_scores.csv'),
            local_input_receipt=relative(link_path),prediction_receipt=relative(receipt_path),prediction=relative(prediction)),
        artifact_sha256s=dict(aggregate=sha(directory/'agg_results.csv'),score_csv=sha(directory/'official_local_scores.csv'),
            local_input_receipt=sha(link_path),prediction_receipt=sha(receipt_path),prediction_recorded=link['prediction_sha256']),
        run_identity={k:meta.get(k) for k in PEERS},
        moment_errors={key:float(np.mean([m[key] for m in moments])) if moments else None
                       for key in ['mean_probability_l1','bulk_probability_l1']},
        avg_score_equals_six_component_mean=True)


def sources_label(arm):
    labels=[]
    names={'K562_GWPS_CPM_full_statistics.npz':'K562 raw','HCT116_full_statistics.npz':'HCT116 raw',
           'HEK293T_full_statistics.npz':'HEK293T raw','PIE_HCT116_DE_statistics.npz':'HCT116 PIE',
           'PIE_HEK293T_DE_statistics.npz':'HEK293T PIE','CD4_DE_statistics.npz':'CD4 DE'}
    for source in arm['source_configuration']['sources']:
        name=Path(source['path']).name
        labels.append(f"{names.get(name,name)} ×{source.get('weight',1):g}")
    return '; '.join(labels) if labels else 'Measured controls only'


def report(manifest):
    panel=pd.read_csv(EVAL/'pert_counts.csv').target_gene.astype(str).tolist()
    reference=json.loads((EVAL/'manifest.json').read_text());reference_hash=reference['reference_sha256']
    if len(panel)!=50:raise ValueError('Expected frozen 50-target H1 Validation panel')
    arms={name:collect_arm(name,manifest,panel,reference_hash) for name in INITIAL_ARMS+RAW_ARMS}
    ranking=sorted(arms,key=lambda n:arms[n]['overall'],reverse=True)
    winner,second=ranking[:2];gap=arms[winner]['overall']-arms[second]['overall']
    contrasts={}
    for name in RAW_ARMS:
        if name=='raw_dual_cpu':continue
        contrasts[name+'_vs_raw_dual_cpu']=dict(overall_delta=arms[name]['overall']-arms['raw_dual_cpu']['overall'],
            scaled_delta={m:arms[name]['scaled'][m]-arms['raw_dual_cpu']['scaled'][m] for m in METRICS})
    anchor=pd.read_parquet(BUNDLE/'anchor_agg.parquet').set_index('metric')
    baseline=pd.read_csv(BUNDLE/'baseline_agg.csv').set_index('statistic')
    result=dict(created_utc=datetime.now(timezone.utc).isoformat(),evaluation='Public H1 2025 Validation; 50 targets; held-out cell context',
        metric_backend='CPU scanpy',control_source='real',source_exclusion='All H1 response banks',
        anchor_type='Fresh five-split competition-enrolled real bundle',anchor_n_splits=5,
        reference_manifest=reference,bundle_manifest=manifest,
        bundle_artifact_sha256s={p.name:sha(p) for p in BUNDLE.iterdir() if p.is_file()},
        score_semantics='Official cell-eval2 score --real-bundle; avg_score from_replicate equals mean of exactly six official components; higher is better',
        scale_ends={m:dict(baseline=float(baseline.loc['mean',m]),replicate=float(anchor.loc[m,'replicate'])) for m in METRICS},
        arms=arms,ranking=ranking,numerical_leader=winner,runner_up=second,leader_margin=gap,
        score_gaps_from_leader={n:arms[winner]['overall']-arms[n]['overall'] for n in ranking},
        contrasts=contrasts,selection_made=False,selection_authority='Root agent/user',
        uncertainty=[
            'The numerical lead is a point-estimate ranking on 50 validation targets, not proof of superior generalization.',
            f'Doubling the HCT116 weight changes the Orion-only score by {contrasts["raw_dual_hct2_cpu_vs_raw_dual_cpu"]["overall_delta"]:+.6f}; adding CD4 changes it by {contrasts["raw_dual_cd4_cpu_vs_raw_dual_cpu"]["overall_delta"]:+.6f}. These small margins do not establish stable ordering beyond this panel.',
            'Small absolute score margins are uncertain: candidates are correlated and no independent candidate-generation replicates or performance confidence intervals were computed.',
            'Five replicate-anchor splits calibrate the score scale; they are not five independent candidate-performance trials.',
            'NMAE can use fewer finite targets than other components; effective counts are recorded separately. The normalized expression metric is aggregate-derived, not a per-target mean.',
            'Multiple candidates were compared on this Validation set; selection can overfit. No Test or K562 confirmation outcomes were read for this report.',
            'These are local scores on public H1 data and locally estimated anchors, not hidden 2026 A/B/C leaderboard estimates.',
        ],pilot_handling=dict(used_for_ranking=False,n_splits=1,path='reports/h1_initial_pilot1_comparison.json',
            note='Earlier one-split diagnostic values are provisional and separate; no pilot or old GPU anchors enter any score in this report.'))
    output=ROOT/'reports/h1_all_candidate_comparison.json';output.parent.mkdir(exist_ok=True)
    output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    lines=[f'The numerical leader on this frozen H1 Validation bundle is **{winner}** '
           f'({arms[winner]["overall"]:.6f}), ahead of {second} by **{gap:.6f}**. '
           'This report does not make the model-selection decision.', '',
           'All 12 arms use the same CPU/scanpy definitions, measured controls, and fresh five-split '
           'competition-enrolled reference bundle. Every overall score was checked against the mean '
           'of exactly six official `from_replicate` components. These are local validation scores.', '',
           '**Reference-scaled metrics (higher is better)**', '',
           '| Arm | Covered | Overall | '+' | '.join(METRICS.values())+' |',
           '| --- | ---: | ---: | '+' | '.join(['---:']*6)+' |']
    for name in ranking:
        a=arms[name];lines.append(f'| {name} | {a["coverage"]["covered"]}/50 | {a["overall"]:.6f} | '+
            ' | '.join(f'{a["scaled"][m]:.6f}' for m in METRICS)+' |')
    lines+=['','**Raw metrics**','',
            'PDS, fidelity, reach and Jaccard are higher-is-better; raw MSE and NMAE are lower-is-better.','',
            '| Arm | '+' | '.join(METRICS.values())+' |','| --- | '+' | '.join(['---:']*6)+' |']
    for name in ranking:lines.append('| '+name+' | '+' | '.join(f'{arms[name]["raw"][m]:.6f}' for m in METRICS)+' |')
    lines+=['','**Source and generation configuration**','',
            '| Arm | Sources and weights | Mode | CPM/bulk amplitudes | Pool | Agreement filter |',
            '| --- | --- | --- | --- | ---: | --- |']
    for name in ranking:
        a=arms[name];args=a['source_configuration']['arguments']
        lines.append(f'| {name} | {sources_label(a)} | {args.get("mode","neutral")} | '
            f'{args.get("amplitude_cpm","—")}/{args.get("amplitude_bulk","—")} | {args.get("pool","—")} | {args.get("agreement",False)} |')
    lines+=['','**Interpretation limits**','']+['- '+s for s in result['uncertainty']]
    lines+=['','Earlier one-split pilot values remain in `reports/h1_initial_pilot1_comparison.*` '
            'as explicitly provisional diagnostics. They were not used here. No old GPU anchors were used.', '',
            'The JSON records all raw/scaled values, finite-target counts, source arguments and hashes, '
            'coverage, paired variant contrasts, reference hashes, and the receipt chain from each '
            '`local_scoring_input.receipt.json` to its prediction and source configuration.']
    output.with_suffix('.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps(dict(numerical_leader=winner,overall=arms[winner]['overall'],runner_up=second,
        margin=gap,ranking={n:arms[n]['overall'] for n in ranking}),indent=2),flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report-only',action='store_true')
    parser.add_argument('--score-only',action='store_true')
    args=parser.parse_args()
    if args.report_only and args.score_only:parser.error('Choose at most one partial operation')
    manifest=load_bundle()
    if not args.report_only:score_raw(manifest)
    if not args.score_only:report(manifest)


if __name__=='__main__':main()
