#!/usr/bin/env python3
"""Compare six raw K562 confirmation metrics without anchors or retuning."""
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:
    os.environ[key]='1'
from datetime import datetime,timezone
import hashlib,json,math
from pathlib import Path
import numpy as np
import pandas as pd
import yaml

ROOT=Path(__file__).resolve().parents[1]
EVAL=ROOT/'evaluations/k562_confirmation'
METRICS={
    'pds_cosine':('PDS','higher'),
    'expr_mse_unbiased_capped_norm':('MSE','lower'),
    'de_wilcoxon_lfc_nmae':('NMAE','lower'),
    'de_wilcoxon_direction_fidelity_yield_raw':('Fidelity','higher'),
    'de_wilcoxon_direction_reach_raw':('Reach','higher'),
    'de_wilcoxon_sig_jaccard':('Jaccard','higher'),
}


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as handle:
        for b in iter(lambda:handle.read(8*1024**2),b''):h.update(b)
    return h.hexdigest()


def number(value):
    value=float(value)
    return value if math.isfinite(value) else None


def fmt(value):return 'undefined' if value is None else f'{value:.6f}'


def main():
    state_path=ROOT/'reports/k562_confirmation_job_state.json'
    state=json.loads(state_path.read_text())
    if state['status']!='complete' or any(j.get('exit_code')!=0 for j in state['jobs'].values()):
        raise RuntimeError('Both supervised metric jobs must have authoritatively exited zero')
    if sha(ROOT/'scripts/run_context_cpu_metrics.sh')!=state['wrapper_sha256']:
        raise ValueError('Scoring wrapper changed during confirmation')
    frozen=json.loads(Path(state['model_snapshot']).read_text())
    if sha(state['model_snapshot'])!=state['model_snapshot_sha256']:
        raise ValueError('Frozen model snapshot changed')
    manifest=json.loads((EVAL/'manifest.json').read_text())
    arms={};identities={}
    for arm in ['selected','neutral']:
        job=state['jobs'][arm];directory=Path(job['output'])
        aggregate=pd.read_csv(directory/'agg_results.csv').set_index('statistic')
        values={m:number(aggregate.loc['mean',m]) for m in METRICS}
        metadata=json.loads((directory/'run_meta.json').read_text())
        params=yaml.safe_load((directory/'run_params.yaml').read_text())
        if (metadata['resolved_device']!='cpu' or metadata['resolved_de_backend']!='scanpy'
                or metadata.get('de_pred_fingerprint') is not None or params.get('control_source')!='real'
                or params.get('metrics')!='vcc2026' or params.get('input_type')!='counts'):
            raise ValueError('Unexpected scorer protocol')
        link=json.loads((directory/'local_scoring_input.receipt.json').read_text())
        if link['reference_sha256']!=manifest['reference_sha256'] or link['prediction_sha256']!=job['prediction_sha256']:
            raise ValueError('Reference/prediction provenance chain differs')
        if sha(job['prediction'])!=job['prediction_sha256']:
            raise ValueError('Prediction changed after independent validation')
        receipt_path=Path(job['prediction']).with_suffix('.receipt.json')
        if sha(receipt_path)!=job['prediction_receipt_sha256']:
            raise ValueError('Prediction receipt changed after validation')
        receipt=json.loads(receipt_path.read_text())
        for source in receipt['sources']:
            if any(token in source['path'].lower() for token in ['k562','gwps']):
                raise ValueError('Forbidden held-out-context response source')
        if sha(directory/'agg_results.csv')!=job['aggregate_sha256']:
            raise ValueError('Aggregate changed after process completion')
        results=pd.read_csv(directory/'results.csv')
        finite={m:int(np.isfinite(results.loc[results.metric==m,'value']).sum()) for m in METRICS}
        finite_targets={m:sorted(results.loc[(results.metric==m)&np.isfinite(results.value),'perturbation'].astype(str).tolist()) for m in METRICS}
        identities[arm]={k:metadata.get(k) for k in ['cell_eval2_version','config_digest','source_fingerprint',
            'source_fingerprint_strict','resolved_device','resolved_de_backend','input_type_real_effective','de_real_fingerprint']}
        arms[arm]=dict(raw_metrics=values,finite_target_counts=finite,finite_target_sets=finite_targets,
            coverage=dict(covered=receipt['same_target_coverage'],total=50,uncovered_targets=receipt['uncovered_targets']),
            prediction_arguments=receipt['arguments'],sources=receipt['sources'],validation=job['validation'],
            terminal_status=dict(status=job['status'],exit_code=job['exit_code'],pid=job['pid'],
                started_utc=job['started_utc'],finished_utc=job['finished_utc']),
            provenance=dict(prediction=job['prediction'],prediction_sha256=job['prediction_sha256'],
                prediction_receipt=str(receipt_path),prediction_receipt_sha256=job['prediction_receipt_sha256'],
                local_scoring_input_receipt=str(directory/'local_scoring_input.receipt.json'),
                aggregate=str(directory/'agg_results.csv'),aggregate_sha256=job['aggregate_sha256'],
                run_meta_sha256=sha(directory/'run_meta.json')),scorer_identity=identities[arm])
    if identities['selected']!=identities['neutral']:
        raise ValueError('Two arms use different scorer/reference identities')
    comparison={};improved=[];worsened=[];ties=[];undefined=[]
    for metric,(label,direction) in METRICS.items():
        selected=arms['selected']['raw_metrics'][metric];neutral=arms['neutral']['raw_metrics'][metric]
        delta=None if selected is None or neutral is None else selected-neutral
        change=None if delta is None else delta*(1 if direction=='higher' else -1)
        status='undefined' if change is None else 'improved' if change>1e-12 else 'worsened' if change < -1e-12 else 'tied'
        {'improved':improved,'worsened':worsened,'tied':ties,'undefined':undefined}[status].append(metric)
        comparison[metric]=dict(label=label,direction=direction,selected=selected,neutral=neutral,
            selected_minus_neutral=delta,improvement_signed=change,status=status)
    limits=[
        'This is one public K562 context on 50 targets selected with a fixed seed among those with at least 100 cached cells; it does not represent every target or the hidden 2026 A/B/C contexts.',
        'The native panel has 8,246 unique measured genes, versus 18,080 in H1 and 18,533 in the challenge. HSPA14/TBCE duplicate source columns were summed; library normalization uses this native gene axis.',
        'Real groups contain 103–200 cells; predictions contain 200 cells per target. Exactly 2,000 measured controls are shared between inference and the reference.',
        f'NMAE is defined for only {arms["selected"]["finite_target_counts"]["de_wilcoxon_lfc_nmae"]}/50 targets and reach for {arms["selected"]["finite_target_counts"]["de_wilcoxon_direction_reach_raw"]}/50; their improvements apply to those eligible subsets. Finite target identities are recorded for both arms. Normalized expression MSE is aggregate-derived and has no per-target output rows.',
        'The reference retains 267 GEM batch labels, but this vcc2026 metric run compares against pooled real controls and does not adjust its DE model for batch effects.',
        'K562 raw responses were a public source for the earlier H1 validation comparison. All K562/GWPS response representations are excluded here; these metrics were unseen when parameters were frozen. This evaluates transfer across contexts using existing public data.',
        'The neutral comparator has zero response-source weights and preserves control moments. It is not the official generic-response baseline.',
        'No anchors were built or reused. There is no reference-scaled overall score, and heterogeneous raw metrics are not averaged.',
        'This confirmation is evaluated once after model freeze; no parameter retuning, model selection, or submission was performed by this workflow.',
    ]
    report=dict(created_utc=datetime.now(timezone.utc).isoformat(),context='K562',purpose='Reserved second-context confirmation of a previously frozen model',
        model_name=frozen['name'],model_frozen_at_utc=frozen['selected_at_utc'],
        frozen_model_snapshot=state['model_snapshot'],frozen_model_snapshot_sha256=state['model_snapshot_sha256'],
        confirmation_started_utc=state['started_utc'],confirmation_finished_utc=state['finished_utc'],
        supervisor_status=state['status'],job_state_path=str(state_path),job_state_sha256=sha(state_path),
        reference_manifest=manifest,reference_manifest_sha256=sha(EVAL/'manifest.json'),
        protocol=dict(metrics='six raw vcc2026 metrics',backend='scanpy',device='cpu',control_source='real',
            excluded_response_source_tokens=['K562','GWPS'],anchors_used=False,scaled_overall_score=None,
            retuned=False,model_selection_performed=False,submitted=False),
        arms=arms,comparison=comparison,improved_metrics=improved,worsened_metrics=worsened,tied_metrics=ties,undefined_metrics=undefined,
        finite_target_sets_match={m:arms['selected']['finite_target_sets'][m]==arms['neutral']['finite_target_sets'][m] for m in METRICS},
        cohort=dict(targets=manifest['targets'],eligible_target_count=manifest['eligible_target_count'],seed=manifest['seed'],
            reference_treated_cells=manifest['n_reference_treated'],controls=manifest['n_control_cells'],genes=manifest['n_genes'],
            real_cells_per_target=manifest['reference_cells_per_target'],predicted_cells_per_target=200,
            source_rows_path=manifest['source_rows_path'],source_rows_sha256=manifest['source_rows_sha256']),limits=limits)
    output=ROOT/'reports/k562_confirmation.json';output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    lines=[f'The frozen selected model improves **{len(improved)} of six raw metrics** over the zero-response control comparator on the reserved K562 cohort '
           f'({len(worsened)} worse, {len(ties)} tied, {len(undefined)} undefined). No anchor-scaled overall score is computed.', '',
           f'The model was frozen at `{frozen["selected_at_utc"]}`; confirmation started at `{state["started_utc"]}`. '
           'Both independent validations excluded K562 and GWPS response sources, and both metric processes exited 0.', '',
           '| Metric | Better | Selected | Neutral | Selected − neutral | Result |',
           '| --- | --- | ---: | ---: | ---: | --- |']
    for metric in METRICS:
        c=comparison[metric];lines.append(f'| {c["label"]} | {c["direction"]} | {fmt(c["selected"])} | {fmt(c["neutral"])} | {fmt(c["selected_minus_neutral"])} | {c["status"]} |')
    lines+=['',f'Selected source coverage is {arms["selected"]["coverage"]["covered"]}/50; neutral coverage is 0/50. '
            'The selected predictor uses H1 × 2, HCT116 × 1, HEK293T × 1 and CD4 DE × 0.5. '
            'The comparator keeps an HCT116 input at weight 0; it contributes no perturbation response.', '',
            'The frozen cohort contains 8,827 genuine treated cells, 2,000 measured controls, 50 targets and 8,246 unique genes. '
            'Target selection used only the fixed seed and the eligibility rule of at least 100 cells. Source row identities and GEM batches are preserved.', '',
            'Limits:', '']+['- '+s for s in limits]
    lines+=['','The JSON contains all six values, finite-target counts, coverage, frozen parameters, source/prediction/reference hashes, '
            'source-row provenance, and authoritative process exit codes. No H1 anchors or K562-based retuning were used.']
    output.with_suffix('.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps(dict(improved=len(improved),worsened=len(worsened),tied=len(ties),comparison=comparison),indent=2),flush=True)


if __name__=='__main__':main()
