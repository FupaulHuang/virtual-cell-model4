#!/usr/bin/env python3
"""Describe null-calibration failure without claiming biological model performance."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'): os.environ.setdefault(k,'1')
import json
from pathlib import Path
from datetime import datetime,timezone
import anndata as ad
import numpy as np
import pandas as pd
import yaml
from scipy import sparse
from build_prediction import ROOT,sha256

EVAL=ROOT/'evaluations/emission_null_400'
ARMS=('measured_null','donor_replay','pool4_no_fit','neutral_dual','response_probe')


def stats(x):
    x=sparse.csr_matrix(x,dtype=np.float64); depth=np.asarray(x.sum(1)).ravel()
    if (depth<=0).any():raise ValueError('Empty cell')
    y=x.multiply((1e6/depth)[:,None]).tocsr()
    mean=np.asarray(y.mean(0)).ravel(); variance=np.maximum(np.asarray(y.power(2).mean(0)).ravel()-mean**2,0)*len(depth)/(len(depth)-1)
    detect=np.asarray((x>0).mean(0)).ravel()
    return mean,variance,detect


def main():
    terminal=json.loads((ROOT/'reports/emission_null_terminal_statuses.json').read_text())
    if any(terminal.get(a,{}).get('exit_code')!=0 for a in ARMS):raise ValueError('All diagnostic jobs must complete')
    freeze=json.loads((ROOT/'reports/emission_null_freeze.json').read_text())
    for name,h in freeze['files'].items():
        if sha256(ROOT/name)!=h:raise ValueError('Frozen file changed')
    manifest=json.loads((EVAL/'manifest.json').read_text()); aliases=manifest['technical_aliases']
    ref=ad.read_h5ad(EVAL/'reference.h5ad',backed='r')
    ctrl=np.flatnonzero(ref.obs.target_gene.astype(str).to_numpy()=='non-targeting')
    mu,var,detect=stats(ref.X[ctrl]);genes=ref.var_names.to_numpy();ref.file.close()
    base_de=pd.read_parquet(EVAL/'measured_null_de/de_real.parquet').sort_values(['target','feature']).reset_index(drop=True)
    if sha256(EVAL/'reference.h5ad')!=manifest['reference_sha256']:raise ValueError('Reference changed')
    result={'completed_at_utc':datetime.now(timezone.utc).isoformat(),'scope':'Control-only emission diagnostic; no challenge score or candidate promotion',
            'cohort':manifest,'terminal_statuses':terminal,'arms':{},'overlap_with_response_probe':{},
            'limitations':['Ten aliases are correlated technical replicates with shared inference/comparison controls, not independent biological experiments.',
            'Real controls are independent of inference cells and comparison cells; batches/guides match marginally, with minimal recorded joint reallocation.',
            'Donor replay retains original donor depths; pooled arms use pooled depths. This compares full policies, not variance alone.',
            'Response-probe calls are not classified as false using null pseudo-outcomes; overlap with zero-response calls is only diagnostic.',
            'No real perturbation or reserved response outcomes were used to build this null reference.']}
    signatures={};expected_config=None;tested_sets=None
    for arm in ARMS:
        folder=EVAL/f'{arm}_de';cfg=yaml.safe_load((folder/'run_params.yaml').read_text());cfg.pop('outdir',None)
        if expected_config is None:expected_config=cfg
        elif cfg!=expected_config:raise ValueError('Diagnostic configuration mismatch')
        if cfg['de']['backend']!='scanpy' or cfg['control_source']!='real':raise ValueError('Unexpected DE settings')
        if cfg['metrics']!=['de_wilcoxon_nsig_counts_real','de_wilcoxon_nsig_counts_pred']:raise ValueError('Not diagnostic-only metrics')
        cur_real=pd.read_parquet(folder/'de_real.parquet').sort_values(['target','feature']).reset_index(drop=True)
        pd.testing.assert_frame_equal(cur_real,base_de,check_dtype=False,check_exact=True)
        table=pd.read_parquet(folder/'de_pred.parquet')
        if table.duplicated(['target','feature']).any():raise ValueError('Duplicate DE pairs')
        current_sets={t:set(table.loc[table.target==t,'feature'])-{t} for t in aliases}
        if tested_sets is None:tested_sets=current_sets
        elif current_sets!=tested_sets:raise ValueError('Tested gene universe differs across arms')
        sig=table[(table.p_adj<.05)&np.isfinite(table.log2_fold_change)&(table.target!=table.feature)]
        signature={t:set(sig.loc[sig.target==t,'feature']) for t in aliases};signatures[arm]=signature
        p=Path(manifest['arms'][arm]['prediction'])
        if sha256(p)!=manifest['arms'][arm]['sha256']:raise ValueError('Prediction bytes changed')
        a=ad.read_h5ad(p,backed='r');labels=a.obs.target_gene.astype(str).to_numpy()
        groups=[]
        for t in aliases:
            rows=np.flatnonzero(labels==t);m,v,d=stats(a.X[rows]);gate=(mu>5)&(genes!=t);common=gate&(var>0)&(v>0)
            groups.append({'alias':t,'significant_calls':len(signature[t]),'tested_genes':len(current_sets[t]),
                'call_fraction':len(signature[t])/max(len(current_sets[t]),1),
                'median_abs_significant_lfc':float(np.abs(sig.loc[sig.target==t,'log2_fold_change']).median()) if len(signature[t]) else None,
                'detection_mae_vs_comparison':float(np.abs(d[gate]-detect[gate]).mean()),
                'variance_trace_ratio_vs_comparison':float(v[gate].sum()/var[gate].sum()),
                'median_log2_variance_ratio':float(np.median(np.log2(v[common]/var[common]))),
                'mean_cpm_mae_vs_comparison':float(np.abs(m[gate]-mu[gate]).mean())})
        a.file.close();calls=[g['significant_calls'] for g in groups]
        result['arms'][arm]={'groups':groups,'median_calls':float(np.median(calls)),'mean_calls':float(np.mean(calls)),
            'min_calls':min(calls),'max_calls':max(calls),'total_calls':sum(calls),
            'median_call_fraction':float(np.median([g['call_fraction'] for g in groups])),
            'mean_detection_mae':float(np.mean([g['detection_mae_vs_comparison'] for g in groups])),
            'mean_variance_trace_ratio':float(np.mean([g['variance_trace_ratio_vs_comparison'] for g in groups]))}
    for arm in ('donor_replay','pool4_no_fit','neutral_dual'):
        intersections=[len(signatures[arm][t]&signatures['response_probe'][t]) for t in aliases]
        total=sum(len(signatures['response_probe'][t]) for t in aliases)
        result['overlap_with_response_probe'][arm]={'shared_calls':sum(intersections),'probe_calls':total,
            'fraction_of_probe_calls_shared':sum(intersections)/max(total,1),
            'not_biological_false_positive_attribution':True}
    measured=result['arms']['measured_null'];donor=result['arms']['donor_replay'];neutral=result['arms']['neutral_dual']
    paired=sum(a['significant_calls']>b['significant_calls'] for a,b in zip(neutral['groups'],donor['groups']))
    ratio=neutral['median_calls']/max(measured['median_calls'],donor['median_calls'],1)
    flag=ratio>=10 and neutral['median_call_fraction']>.01 and paired>=8
    result['registered_diagnostic']={'neutral_to_larger_empirical_null_median_ratio':ratio,'neutral_exceeds_donor_groups':paired,
        'flag':bool(flag),'criteria':'ratio>=10, neutralmedianfraction>1%, neutral>donor in>=8/10groups',
        'interpretation':'Prioritize a control-distribution-preserving zero-response emitter before further global response tuning.' if flag else 'No decisive emitter-excess flag; inspect empirical background before attributing failures to response transfer.'}
    out=ROOT/'reports/emission_null_results.json';out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    lines=['# Zero-response emission diagnostic','',result['registered_diagnostic']['interpretation'],'',
           '| Arm | Median DE calls | Range | Detection MAE | Variance-trace ratio |',
           '| --- | ---: | --- | ---: | ---: |']
    for arm,x in result['arms'].items():lines.append(f'| {arm} | {x["median_calls"]:.1f} | {x["min_calls"]}–{x["max_calls"]} | {x["mean_detection_mae"]:.5f} | {x["mean_variance_trace_ratio"]:.4f} |')
    lines+=['',f'Registered emitter-excess flag: **{flag}**; ratio versus the larger empirical-null median: **{ratio:.1f}×**; neutral exceeds donor replay in **{paired}/10** aliases.',
            '',f'Joint reallocation: {manifest["joint_reallocated_cells"]}/4000 cells; batch and guide margins remain exact.',
            '', '## Overlap with existing response calls','']
    for arm,x in result['overlap_with_response_probe'].items():lines.append(f'- {arm}: {x["shared_calls"]:,}/{x["probe_calls"]:,} response-probe calls shared ({100*x["fraction_of_probe_calls_shared"]:.1f}%). This is diagnostic overlap, not attribution of biological false positives.')
    lines+=['','## Limits','']+['- '+x for x in result['limitations']]
    out.with_suffix('.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({'flag':flag,'arms':{a:x['median_calls'] for a,x in result['arms'].items()}},indent=2))


if __name__=='__main__':main()
