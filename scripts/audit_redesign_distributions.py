#!/usr/bin/env python3
"""Posthoc descriptive distribution audit of the eight frozen redesign outputs.

No scorer is run and no promotion criterion is introduced. Only the extracted
H1 Training development reference and existing prediction artifacts are read.
"""
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(key, '1')
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
import argparse
import hashlib
import json
import multiprocessing as mp
from pathlib import Path
import time

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
REFERENCE = ROOT/'evaluations/dispersion_h1_training_400/reference.h5ad'
SEEDS = (20261021, 20261022)
QUANTILES = (0, .05, .25, .5, .75, .95, 1)
Q_NAMES = ('q00', 'q05', 'q25', 'q50', 'q75', 'q95', 'q100')


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(8*1024**2), b''):
            h.update(block)
    return h.hexdigest()


def summarize(values):
    a = np.asarray(values, dtype=float)
    return dict(n=len(a), mean=float(a.mean()), median=float(np.median(a)),
                min=float(a.min()), p05=float(np.quantile(a,.05)),
                p95=float(np.quantile(a,.95)), max=float(a.max()))


def cell_statistics(value):
    """Sparse sample moments, including implicit zeros; no dense cell×gene copy."""
    x = sparse.csr_matrix(value)
    if (not x.has_canonical_format or not np.isfinite(x.data).all()
            or (x.data <= 0).any() or (x.data != np.floor(x.data)).any()):
        raise ValueError('Counts must be canonical CSR with finite positive integer stored values')
    n, g = x.shape
    depth = np.asarray(x.sum(axis=1, dtype=np.int64)).ravel()
    if n < 2 or (depth <= 0).any() or (depth > 1_000_000).any():
        raise ValueError('Invalid group size or cell library depth')
    data = x.data.astype(np.float64)
    raw_sum = np.bincount(x.indices, weights=data, minlength=g)
    raw_square_sum = np.bincount(x.indices, weights=data*data, minlength=g)
    raw_variance = np.maximum(raw_square_sum-raw_sum*raw_sum/n, 0)/(n-1)
    data *= np.repeat(1e6/depth, np.diff(x.indptr))
    cpm_sum = np.bincount(x.indices, weights=data, minlength=g)
    cpm_square_sum = np.bincount(x.indices, weights=data*data, minlength=g)
    variance = np.maximum(cpm_square_sum-cpm_sum*cpm_sum/n, 0)/(n-1)
    detection = np.bincount(x.indices, minlength=g)/n
    return dict(mean_cpm=cpm_sum/n, variance_cpm=variance, raw_count_variance=raw_variance,
        detection=detection, mean_detected_genes_all=float(np.diff(x.indptr).mean()),
        depth_quantiles=np.quantile(depth, QUANTILES), depth_mean=float(depth.mean()),
        depth_sd=float(depth.std(ddof=1)), nnz=int(x.nnz), n_cells=n)


def small_formula_check():
    x = np.asarray([[0,2,3,8], [1,4,0,2], [3,0,5,2], [2,6,1,0]],dtype=np.int32)
    m = cell_statistics(sparse.csr_matrix(x))
    cpm = x/x.sum(axis=1)[:,None]*1e6
    np.testing.assert_allclose(m['mean_cpm'], cpm.mean(axis=0),rtol=1e-14,atol=1e-8)
    np.testing.assert_allclose(m['variance_cpm'], cpm.var(axis=0,ddof=1),rtol=1e-13)
    np.testing.assert_allclose(m['raw_count_variance'], x.var(axis=0,ddof=1),rtol=1e-14)
    np.testing.assert_array_equal(m['detection'], (x>0).mean(axis=0))
    assert m['mean_detected_genes_all'] == float((x>0).sum(axis=1).mean())


_GENES = _GATE = _REAL = _TARGETS = None


def initialize(genes, gate, real, targets):
    global _GENES, _GATE, _REAL, _TARGETS
    _GENES, _GATE, _REAL, _TARGETS = genes, gate, real, targets


def audit_arm(item):
    response, emitter, seed = item
    name=f'{response}_{emitter}_s{seed}'
    path=ROOT/'outputs/redesign_h1_400'/name/'prediction.h5ad'
    rp=path.with_suffix('.receipt.json');receipt=json.loads(rp.read_text())
    actual_sha=sha(path)
    if actual_sha!=receipt['output_sha256'] or path.stat().st_size!=receipt['bytes']:
        raise ValueError('Prediction hash or byte count differs from receipt: '+name)
    meta=receipt['redesign']
    if (meta['response'],meta['emitter'],meta['seed'])!=(response,emitter,seed) or meta['engineering_subset']:
        raise ValueError('Incorrect arm identity or incomplete experiment artifact')
    obj=ad.read_h5ad(path,backed='r');records=[];nnz=0
    try:
        if obj.shape!=(len(_TARGETS)*400,len(_GENES)) or not np.array_equal(obj.var_names.astype(str),_GENES):
            raise ValueError('Prediction axes differ: '+name)
        if not obj.obs_names.is_unique or set(obj.obs.context.astype(str))!={'H1'}:
            raise ValueError('Prediction cell IDs/context invalid')
        labels=obj.obs.target_gene.astype(str).to_numpy()
        if set(labels)!=set(_TARGETS) or set(receipt['targets'])!=set(_TARGETS):
            raise ValueError('Prediction target panel differs')
        for target in _TARGETS:
            rows=np.flatnonzero(labels==target)
            if len(rows)!=400:raise ValueError('Group does not contain400cells')
            pred=cell_statistics(obj.X[rows]);real=_REAL[target];nnz+=pred['nnz']
            gate=_GATE&(_GENES!=target)
            pv=pred['variance_cpm'][gate];rv=real['variance_cpm'][gate]
            trace_pred=float(pv.sum());trace_real=float(rv.sum())
            if trace_real<=0:raise ValueError('Reference CPM variance trace is zero')
            valid=(pv>0)&(rv>0)
            log_ratio=np.log2(pv[valid]/rv[valid])
            row=dict(arm=name,response=response,emitter=emitter,seed=seed,context='H1',target=target,
                n_pred=400,n_real=real['n_cells'],gated_genes=int(gate.sum()),
                own_gene_present=bool(np.any(_GENES==target)),
                pred_mean_detected_genes_all=pred['mean_detected_genes_all'],
                real_mean_detected_genes_all=real['mean_detected_genes_all'],
                pred_mean_detected_genes_gated=float(pred['detection'][gate].sum()),
                real_mean_detected_genes_gated=float(real['detection'][gate].sum()),
                gene_detection_mae=float(np.abs(pred['detection'][gate]-real['detection'][gate]).mean()),
                gene_detection_bias=float((pred['detection'][gate]-real['detection'][gate]).mean()),
                pred_cpm_variance_trace=trace_pred,real_cpm_variance_trace=trace_real,
                cpm_variance_trace_ratio=trace_pred/trace_real,
                median_gene_log2_cpm_variance_ratio=float(np.median(log_ratio)),
                mean_abs_gene_log2_cpm_variance_ratio=float(np.abs(log_ratio).mean()),
                positive_variance_pairs=int(valid.sum()),
                fraction_gene_variance_lower_than_real=float(np.mean(pv<rv)),
                pred_zero_real_positive_variance_genes=int(np.count_nonzero((pv==0)&(rv>0))),
                pred_raw_count_variance_trace=float(pred['raw_count_variance'][gate].sum()),
                real_raw_count_variance_trace=float(real['raw_count_variance'][gate].sum()),
                raw_count_variance_trace_ratio=float(pred['raw_count_variance'][gate].sum()/real['raw_count_variance'][gate].sum()),
                pred_depth_mean=pred['depth_mean'],real_depth_mean=real['depth_mean'],
                pred_depth_sd=pred['depth_sd'],real_depth_sd=real['depth_sd'])
            for prefix,stats in [('pred',pred),('real',real)]:
                row.update({prefix+'_depth_'+q:float(v) for q,v in zip(Q_NAMES,stats['depth_quantiles'])})
            row['depth_median_ratio']=row['pred_depth_q50']/row['real_depth_q50']
            records.append(row)
        if nnz!=receipt['stored_entries'] or receipt['cells']!=obj.n_obs or receipt['output_genes']!=obj.n_vars:
            raise ValueError('Actual count/axis totals differ from receipt')
    finally:obj.file.close()
    return dict(arm=name,response=response,emitter=emitter,seed=seed,
        prediction=str(path.relative_to(ROOT)),prediction_sha256=actual_sha,
        receipt=str(rp.relative_to(ROOT)),receipt_sha256=sha(rp),
        hash_verified=True,actual_axes_verified=True,all_counts_verified=True,
        groups=len(records),nnz=nnz,records=records)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--workers',type=int,default=4)
    args=parser.parse_args()
    if not 1<=args.workers<=4:raise ValueError('Bounded audit permits1–4workers')
    outputs=[ROOT/f'reports/redesign_distribution_audit.{suffix}' for suffix in ('json','csv','md')]
    if any(p.exists() for p in outputs):raise FileExistsError('Audit output exists; no overwrite')
    start=time.time();small_formula_check()
    frozen={}
    for f in ('redesign_implementation_freeze.json','redesign_analysis_freeze.json'):
        path=ROOT/'reports'/f;record=json.loads(path.read_text())
        checks={p:sha(ROOT/p)==h for p,h in record['files'].items()}
        if not all(checks.values()):raise ValueError('Frozen implementation/analysis dependency changed')
        frozen[f]={'sha256':sha(path),'all_listed_file_hashes_verified':True,'files':record['files']}
    fix=ROOT/'evaluations/dispersion_h1_training_400/reference_category_fix.json'
    reference_sha=sha(REFERENCE)
    if reference_sha!=json.loads(fix.read_text())['new_reference_sha256']:
        raise ValueError('Development reference differs from frozen repaired reference')
    ref=ad.read_h5ad(REFERENCE,backed='r');real={}
    try:
        genes=ref.var_names.astype(str).to_numpy()
        if not ref.var_names.is_unique or not ref.obs_names.is_unique or set(ref.obs.context.astype(str))!={'H1'}:
            raise ValueError('Invalid measured reference axes')
        labels=ref.obs.target_gene.astype(str).to_numpy();controls=np.flatnonzero(labels=='non-targeting')
        targets=sorted(set(labels)-{'non-targeting'})
        if len(controls)!=4000 or len(targets)!=50 or ref.shape!=(24000,18080):
            raise ValueError('Unexpected frozen development cohort')
        control_mean=np.zeros(len(genes))
        for left in range(0,len(controls),400):
            block=cell_statistics(ref.X[controls[left:left+400]])
            control_mean+=block['mean_cpm']*block['n_cells']
        control_mean/=len(controls);gate=control_mean>5
        for target in targets:
            rows=np.flatnonzero(labels==target)
            if len(rows)!=400:raise ValueError('Reference group must contain400cells')
            real[target]=cell_statistics(ref.X[rows])
    finally:ref.file.close()
    print(json.dumps({'reference_verified':True,'targets':len(targets),'control_gate_genes':int(gate.sum()),'seconds':time.time()-start}),flush=True)
    jobs=[(r,e,s) for r in ('current','adapter') for e in ('current','preserving') for s in SEEDS]
    results=[]
    with ProcessPoolExecutor(max_workers=args.workers,mp_context=mp.get_context('spawn'),
                             initializer=initialize,initargs=(genes,gate,real,targets)) as pool:
        futures=[pool.submit(audit_arm,job) for job in jobs]
        for future in as_completed(futures):
            result=future.result();results.append(result)
            print(json.dumps({'complete':result['arm'],'groups':result['groups'],'seconds':time.time()-start}),flush=True)
    results.sort(key=lambda x:x['arm'])
    frame=pd.DataFrame([row for result in results for row in result.pop('records')])
    if len(frame)!=400:raise ValueError('Incomplete posthoc audit')
    frame.to_csv(outputs[1],index=False)
    summary_keys=['cpm_variance_trace_ratio','raw_count_variance_trace_ratio','gene_detection_mae','gene_detection_bias',
        'pred_mean_detected_genes_all','real_mean_detected_genes_all','pred_mean_detected_genes_gated','real_mean_detected_genes_gated',
        'pred_depth_q05','real_depth_q05','pred_depth_q50','real_depth_q50','pred_depth_q95','real_depth_q95','depth_median_ratio',
        'median_gene_log2_cpm_variance_ratio','mean_abs_gene_log2_cpm_variance_ratio','fraction_gene_variance_lower_than_real']
    aggregates=[]
    for result in results:
        f=frame[frame.arm==result['arm']]
        aggregates.append(dict(arm=result['arm'],response=result['response'],emitter=result['emitter'],seed=result['seed'],
            targets=len(f),statistics={k:summarize(f[k]) for k in summary_keys}))
    pairs=[]
    for response in ('current','adapter'):
        for seed in SEEDS:
            a=frame[(frame.response==response)&(frame.emitter=='current')&(frame.seed==seed)].set_index('target').sort_index()
            b=frame[(frame.response==response)&(frame.emitter=='preserving')&(frame.seed==seed)].set_index('target').sort_index()
            pairs.append(dict(response=response,seed=seed,contrast='preserving minus current emitter; same response and seed',
                changes={k:summarize(b[k]-a[k]) for k in ('cpm_variance_trace_ratio','gene_detection_mae','gene_detection_bias','pred_mean_detected_genes_all','pred_depth_q50')}))
    report=dict(created_utc=datetime.now(timezone.utc).isoformat(),status='posthoc_descriptive_no_new_promotion_criteria',
        purpose='Describe actual nonzero400-cell prediction distributions across both frozen seeds; no rescoring, retuning or reserved data access.',
        script_sha256=sha(__file__),frozen_artifacts=frozen,reference=dict(path=str(REFERENCE.relative_to(ROOT)),sha256=reference_sha,
            reference_fix_receipt_sha256=sha(fix),shape=[24000,18080],controls=4000,targets=50,cells_per_target=400),
        methods=dict(counts='Actual canonical sparse integer counts validated per400-rowtarget; all8prediction hashes verified against receipts.',
            gate='Arithmetic mean CPM>5 in the4000realcontrols, then exclude each own target gene; same gate for every prediction.',
            gate_genes_before_own_target_exclusion=int(gate.sum()),variance='Sample variance ddof1 of per-cellCPM normalized on full18080-gene axis; variance trace is sum over gated genes. Units CPM squared.',
            raw_variance='Supplementary sample raw-count variance trace in squared-count units; differs from CPM variance and is depth-dependent.',
            detection='Positive-cell fractions per gene; MAE on the gated genes. Mean detectedgenes per cell also reported on full axis and gate.',
            depth='Raw integer library totals; pertarget quantiles0,.05,.25,.5,.75,.95,1 with linear interpolation.',
            aggregation='Each target weighted equally. Table depth quantiles are macro averages of pertarget quantiles, not pooled-cell quantiles.',
            mean_confound='Real/predicted means may differ; variance ratios are distribution descriptions, not pure dispersion causal effects.',
            self_check='Sparse moments compared with independent dense NumPy mean/samplevariance/detection on a small fixture before audit.'),
        predictions=results,aggregates=aggregates,paired_emitter_changes=pairs,csv=str(outputs[1].relative_to(ROOT)),csv_sha256=sha(outputs[1]),
        elapsed_seconds=time.time()-start,no_metric_rerun=True,no_reserved_response_access=True)
    outputs[0].write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    lines=['# Posthoc nonzero-response distribution audit','',
        'All eight frozen predictions were checked against their receipts and the actual 400-cell-per-target measured reference. All count values, target groups and ordered gene axes passed. This is descriptive analysis after the experiment outcomes; it adds no promotion criterion and does not change any previous failure or selection decision. No metric was rerun and no reserved responses were accessed.','',
        'The shared gene gate is real-control arithmetic mean CPM >5, excluding each target’s own gene. CPM variances use sample ddof1 and full-axis per-cell normalization. Variance-trace ratios compare the sum of predicted gene variances with the sum of measured gene variances on that gate. Means may differ, so these are not pure dispersion effects.','',
        '| Response | Emitter | Seed | CPM variance trace ratio, mean (median) | Gene detection MAE | Detected genes/cell, prediction / real | Depth q05 / q50 / q95, prediction |','| --- | --- | ---: | ---: | ---: | ---: | ---: |']
    for a in aggregates:
        s=a['statistics'];v=s['cpm_variance_trace_ratio']
        lines.append(f"| {a['response']} | {a['emitter']} | {a['seed']} | {v['mean']:.4f} ({v['median']:.4f}) | {s['gene_detection_mae']['mean']:.5f} | {s['pred_mean_detected_genes_all']['mean']:.1f} / {s['real_mean_detected_genes_all']['mean']:.1f} | {s['pred_depth_q05']['mean']:.0f} / {s['pred_depth_q50']['mean']:.0f} / {s['pred_depth_q95']['mean']:.0f} |")
    base=aggregates[0]['statistics']
    lines+=['',f"Measured depth quantiles averaged across targets are {base['real_depth_q05']['mean']:.0f} / {base['real_depth_q50']['mean']:.0f} / {base['real_depth_q95']['mean']:.0f}. Depth table values average per-target quantiles; they are not pooled-cell quantiles. Detection MAE is a fraction, so0.01 is one percentage point.",'',
        'The [CSV](redesign_distribution_audit.csv) contains all400arm/seed/target rows, full depth quantiles, gene detection, CPM variance traces and supplementary raw-count variance traces. Raw-count variance is explicitly in squared-count units and includes depth effects. The [JSON](redesign_distribution_audit.json) preserves artifact/receipt/freeze hashes, gate and normalization definitions, per-arm summaries and paired emitter differences.','',
        'This audit describes nonzero predictions, separately from the zero-effect identity/null audit. A more realistic distribution does not by itself establish improved DE directions, response specificity, expression error or the six-metric promotion gate.']
    outputs[2].write_text('\n'.join(lines)+'\n')
    print(json.dumps({'written':[str(p.relative_to(ROOT)) for p in outputs],'seconds':time.time()-start}),flush=True)


if __name__=='__main__':main()
