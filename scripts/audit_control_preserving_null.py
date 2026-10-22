#!/usr/bin/env python3
"""Frozen, control-only ABC null diagnostic; never a submission artifact.

Stages: freeze metadata/protocol/code; prepare counts; de per context/arm;
identity against the new emitter when available; summarize. Existing outputs
are never overwritten. No perturbed-cell reference or reserved target is read.
"""
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMBA_NUM_THREADS'):
    os.environ.setdefault(key, '1')
os.environ.setdefault('POLARS_MAX_THREADS', '2')
from pathlib import Path
import argparse
from datetime import datetime, timezone
import hashlib
import importlib
import importlib.metadata
import inspect
import json
import time

import anndata as ad
import numpy as np
import pandas as pd
import polars as pl
from scipy import sparse
from scipy.stats import ks_2samp
from cell_eval2.de_compute import compute_de
from build_prediction import ROOT, Writer, sha256, atlas
from prepare_emission_null import round_template
from validate_prediction_receipt import validate

OUT = ROOT / 'evaluations/control_preserving_null_abc'
REPORTS = ROOT / 'reports'
CONTEXTS = ('A', 'B', 'C')
ARMS = ('measured_null', 'donor_replay', 'pool4_no_fit')
ALIASES = [f'NULL_TECH_{i+1:02d}' for i in range(10)]
SEED = 20261024
PROTOCOL = {
    'purpose': 'Control-only technical null diagnostic. Aliases are not perturbations; never submit.',
    'contexts': list(CONTEXTS), 'seed': SEED, 'aliases': ALIASES,
    'cells_per_alias': 400, 'comparison_cells': 4000, 'donor_population_cells': 4000,
    'roles': 'Comparison4000, measured-null4000, donor-population4000: mutually disjoint within context.',
    'ntc_stratification': 'Fixed seeded 400-cell integer guide quota; identical per-alias quotas and comparison proportions. Each role has10times that quota.',
    'pooling': 'Sample4times each alias guide quota without replacement from donor population; depth-sort WITHIN guide; groups of4. Replay chooses1random raw cell per pool; pool4 averages normalized cells and uses original rounding/repair without moment fitting.',
    'difference_from_h1_null': 'Pools stay within ntc_id here, unlike prior H1 global depth pooling; this avoids a guide-mixing confound.',
    'reuse': 'Donor cells may recur across technical aliases; aliases and shared-comparison tests are correlated, not independent biological replicates.',
    'de': {'backend': 'scanpy', 'input_type': 'counts', 'target_sum': 1e6,
           'mean_calc': 'arithmetic', 'epsilon': 1e-9, 'clip_value': None,
           'filter_gene_min_cpm_cell': 5.0, 'fdr_scope': 'per_pert', 'threads': 1,
           'groupby': 'target_gene', 'reference': 'non-targeting'},
    'calls': 'Strict p_adj<0.05; reference-only arithmetic CPM>5 gate and recomputed per-alias BH are provided by the pinned compute_de helper. Aliases do not match genes, so own-target exclusion removes nothing.',
    'distribution': 'Full-row-depth CPM; sample variance ddof1; detection=positive fraction. Gene summaries on comparison-control meanCPM>5, positive-variance subset for log ratios. Depth KS is descriptive.',
    'identity_gate': 'fit_population at donor population own arithmetic mean/raw bulk, then emit the exact replay donor rows. Require expected probabilities and every emitted count to match raw donors exactly; no DE rerun if identity is proven.',
    'identity_status_at_freeze': 'Pending implementation/API availability; measured/replay/pool4 diagnostic is independently defined.',
    'compute_bound': 'Nine DE jobs, each10aliases versus4000 controls (90 contrasts); at most3 processes concurrently; no scorer or parameter sweep.',
    'no_hidden_or_reserved_response_access': True,
}


def seed_for(label):
    return (int.from_bytes(hashlib.sha256(label.encode()).digest()[:4], 'little') + SEED) % 2**32


def write_json(path, value):
    path = Path(path)
    if path.exists():
        raise FileExistsError(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def read_json(path):
    return json.loads(Path(path).read_text())


def validate_counts(x):
    x = sparse.csr_matrix(x)
    if not np.isfinite(x.data).all() or (x.data < 0).any() or (x.data != np.floor(x.data)).any():
        raise ValueError('Expected finite nonnegative integer counts')
    depth = np.asarray(x.sum(axis=1, dtype=np.int64)).ravel()
    if (depth <= 0).any() or (depth > 1_000_000).any():
        raise ValueError('Invalid control library depth')
    return x.astype(np.int32), depth


def freeze():
    if OUT.exists():
        raise FileExistsError('Null directory exists; do not resample or overwrite')
    source_rows, pool_rows, contexts = [], [], {}
    for context in CONTEXTS:
        path = ROOT/f'data/context_{context}.h5ad'
        obj = ad.read_h5ad(path, backed='r')
        try:
            if set(obj.obs.target_gene.astype(str)) != {'non-targeting'} or set(obj.obs.context.astype(str)) != {context}:
                raise ValueError('Source is not a pure named control context')
            if not obj.obs_names.is_unique or 'ntc_id' not in obj.obs:
                raise ValueError('Unique cell IDs and ntc_id required')
            if set(ALIASES) & set(obj.var_names):
                raise ValueError('Technical alias unexpectedly matches a gene')
            labels = obj.obs.ntc_id.astype(str).to_numpy()
            strata = np.asarray(sorted(set(labels)))
            if not 1 <= len(strata) <= 400:
                raise ValueError('Unsupported guide-stratum count')
            quota = np.full(len(strata), 400//len(strata), dtype=int)
            extra = np.random.default_rng(seed_for(context+':quota')).permutation(len(strata))[:400 % len(strata)]
            quota[extra] += 1
            for stratum, q in zip(strata, quota):
                eligible = np.flatnonzero(labels == stratum)
                if len(eligible) < 30*q:
                    raise ValueError(f'Insufficient controls for exact disjoint quotas: {context}/{stratum}')
                selected = np.random.default_rng(seed_for(context+':'+stratum)).choice(eligible, 30*q, replace=False)
                for role, chosen in [('comparison', selected[:10*q]), ('measured_null', selected[10*q:20*q]), ('donor_population', selected[20*q:])]:
                    for i, row in enumerate(chosen):
                        source_rows.append(dict(context=context, source_row=int(row), source_cell=str(obj.obs_names[row]),
                            ntc_id=str(stratum), role=role, alias=ALIASES[i//q] if role=='measured_null' else 'non-targeting'))
                donors = selected[20*q:]
                for alias in ALIASES:
                    chosen = np.random.default_rng(seed_for(context+':'+alias+':'+stratum)).choice(donors, 4*q, replace=False)
                    for row in chosen:
                        pool_rows.append(dict(context=context, alias=alias, ntc_id=str(stratum), source_row=int(row)))
            contexts[context] = dict(source=str(path), source_resolved=str(path.resolve()), source_sha256=sha256(path),
                shape=list(obj.shape), ntc_counts=pd.Series(labels).value_counts().sort_index().to_dict(),
                ntc_quota_per_alias=dict(zip(strata.tolist(), quota.tolist())),
                gene_axis_sha256=hashlib.sha256('\0'.join(obj.var_names).encode()).hexdigest())
        finally:
            obj.file.close()
    rows = pd.DataFrame(source_rows)
    for context in CONTEXTS:
        c = rows[rows.context == context]
        if c.source_row.duplicated().any() or len(c) != 12000:
            raise ValueError('Role allocation is not disjoint')
    OUT.mkdir(parents=True)
    rows.to_csv(OUT/'source_rows.csv', index=False)
    pd.DataFrame(pool_rows).to_csv(OUT/'pool_selections.csv', index=False)
    write_json(OUT/'protocol.json', PROTOCOL)
    files = [Path(__file__), ROOT/'scripts/prepare_emission_null.py', ROOT/'scripts/build_prediction.py',
             ROOT/'scripts/validate_prediction_receipt.py', ROOT/'external/atlasshift/model.py',
             Path(inspect.getsourcefile(compute_de)), Path(inspect.getsourcefile(compute_de)).with_name('de.py'),
             Path(inspect.getsourcefile(compute_de)).with_name('norm.py'),
             Path(inspect.getsourcefile(compute_de)).parent/'configs/vcc2026.yaml']
    record = dict(created_utc=datetime.now(timezone.utc).isoformat(), protocol=PROTOCOL,
        files={str(p.resolve()):sha256(p) for p in files}, contexts=contexts,
        source_rows_sha256=sha256(OUT/'source_rows.csv'), pool_selections_sha256=sha256(OUT/'pool_selections.csv'),
        protocol_sha256=sha256(OUT/'protocol.json'),
        versions={n:importlib.metadata.version(n) for n in ('cell-eval2','scanpy','numpy','scipy','anndata','polars')},
        expression_outcomes_read_before_freeze=False)
    write_json(OUT/'cohort_freeze.json', record)
    write_json(REPORTS/'control_preserving_null_freeze.json', record)
    protocol_path=REPORTS/'control_preserving_null_protocol.md'
    if protocol_path.exists():raise FileExistsError(protocol_path)
    protocol_path.write_text('# Frozen ABC control-only null protocol\n\n'+
        '\n\n'.join('**'+key.replace('_',' ')+'**: '+(json.dumps(value) if not isinstance(value,str) else value) for key,value in PROTOCOL.items())+
        '\n\nSource identities, all selected row identities, code/helper hashes and package versions are recorded in `control_preserving_null_freeze.json` and `evaluations/control_preserving_null_abc/`. No expression outcomes were inspected before the metadata freeze.\n')
    print(json.dumps({'frozen':True,'contexts':3,'selected_cells_per_context':12000,'source_rows_sha256':record['source_rows_sha256']}),flush=True)


def verify_freeze():
    f = read_json(OUT/'cohort_freeze.json')
    for p, digest in f['files'].items():
        if sha256(p) != digest:
            raise ValueError('Frozen code changed: '+p)
    for name, key in [('source_rows.csv','source_rows_sha256'),('pool_selections.csv','pool_selections_sha256'),('protocol.json','protocol_sha256')]:
        if sha256(OUT/name) != f[key]:
            raise ValueError('Frozen metadata changed: '+name)
    return f


def moments(x):
    x, depth = validate_counts(x)
    q = x.astype(np.float64).multiply((1e6/depth)[:,None]).tocsr()
    mean = np.asarray(q.mean(axis=0)).ravel()
    second = np.asarray(q.power(2).mean(axis=0)).ravel()
    var = np.maximum(second-mean*mean,0)*len(depth)/(len(depth)-1)
    detected = np.asarray((x>0).mean(axis=0)).ravel()
    return dict(mean=mean, variance=var, detection=detected, depth=depth)


def distribution(x, comparison):
    m = moments(x);gate=comparison['mean']>5
    valid=gate&(comparison['variance']>0)&(m['variance']>0)
    logs=np.log2(m['variance'][valid]/comparison['variance'][valid])
    report=dict(gated_genes=int(gate.sum()),variance_log_genes=int(valid.sum()),
        median_log2_variance_ratio=float(np.median(logs)),mean_abs_log2_variance_error=float(np.mean(np.abs(logs))),
        fraction_variance_lower=float(np.mean(logs<0)),
        detection_mae=float(np.mean(np.abs(m['detection'][gate]-comparison['detection'][gate]))),
        detection_bias=float(np.mean(m['detection'][gate]-comparison['detection'][gate])),
        arithmetic_cpm_mean_l1=float(np.abs(m['mean']-comparison['mean']).sum()),
        depth_mean=float(m['depth'].mean()),depth_median=float(np.median(m['depth'])),
        depth_sd=float(m['depth'].std(ddof=1)),depth_min=int(m['depth'].min()),depth_max=int(m['depth'].max()),
        depth_ks_to_comparison=float(ks_2samp(m['depth'],comparison['depth']).statistic))
    return m,report


def prepare(context):
    frozen=verify_freeze();folder=OUT/context
    if folder.exists():raise FileExistsError(folder)
    source_path=Path(frozen['contexts'][context]['source'])
    if sha256(source_path)!=frozen['contexts'][context]['source_sha256']:raise ValueError('Source changed')
    table=pd.read_csv(OUT/'source_rows.csv');table=table[table.context==context]
    selections=pd.read_csv(OUT/'pool_selections.csv');selections=selections[selections.context==context]
    source=ad.read_h5ad(source_path,backed='r');genes=source.var_names.astype(str).to_numpy()
    selected=np.sort(table.source_row.to_numpy());lookup={int(r):i for i,r in enumerate(selected)}
    blocks=[sparse.csr_matrix(source.X[selected[k:k+512]]) for k in range(0,len(selected),512)]
    raw,depth=validate_counts(sparse.vstack(blocks,format='csr'));del blocks
    def get(rows):return raw[[lookup[int(r)] for r in rows]]
    folder.mkdir()
    pd.DataFrame({'gene_name':genes}).to_csv(folder/'gene_names.csv',index=False)
    pd.DataFrame({'target_gene':ALIASES,'n_cells':400}).to_csv(folder/'pert_counts.csv',index=False)
    role_paths={}
    for role in ('comparison','donor_population'):
        rows=np.sort(table.loc[table.role==role,'source_row'].to_numpy())
        obs=source.obs.iloc[rows].copy()
        for col in obs:
            if isinstance(obs[col].dtype,pd.CategoricalDtype):obs[col]=obs[col].cat.remove_unused_categories()
        obj=ad.AnnData(get(rows),obs=obs,var=source.var.copy())
        path=folder/(role+'.h5ad');obj.write_h5ad(path,compression='gzip',compression_opts=1)
        role_paths[role]={'path':str(path),'sha256':sha256(path),'source_rows':rows.tolist()}
    comparison=moments(get(role_paths['comparison']['source_rows']))
    pool_table=[];replay_table=[];by_alias={}
    for alias in ALIASES:
        pooled=[];replay=[]
        for ntc in sorted(selections.ntc_id.unique()):
            rows=selections[(selections.alias==alias)&(selections.ntc_id==ntc)].source_row.to_numpy()
            dd=np.asarray(get(rows).sum(axis=1,dtype=np.int64)).ravel()
            rows=rows[np.argsort(dd,kind='stable')];pools=rows.reshape(-1,4)
            chosen=np.random.default_rng(seed_for(context+':'+alias+':replay:'+ntc)).integers(0,4,len(pools))
            for j,pool in enumerate(pools):
                group=len(pooled);pooled.append(pool);replay.append(pool[chosen[j]])
                pool_table.extend(dict(alias=alias,pool=group,within_pool=k,source_row=int(row),ntc_id=ntc) for k,row in enumerate(pool))
                replay_table.append(dict(alias=alias,row=group,source_row=int(pool[chosen[j]]),ntc_id=ntc))
        by_alias[alias]=(np.asarray(pooled),np.asarray(replay))
    pd.DataFrame(pool_table).to_csv(folder/'donor_pools.csv',index=False)
    pd.DataFrame(replay_table).to_csv(folder/'donor_replay_rows.csv',index=False)
    manifests={};gene_frames=[];distribution_reports=[]
    for arm in ARMS:
        out=folder/arm;out.mkdir();path=out/'prediction.h5ad';pending=out/'prediction.partial.h5ad'
        writer=Writer(pending,genes);audits=[]
        try:
            for alias in ALIASES:
                pools,replay=by_alias[alias]
                if arm=='measured_null':
                    rows=np.sort(table[(table.role=='measured_null')&(table.alias==alias)].source_row.to_numpy());x=get(rows).toarray()
                elif arm=='donor_replay':x=get(replay).toarray()
                else:
                    donors=get(pools.ravel()).toarray();d=donors.sum(axis=1,dtype=np.int64)
                    template=(donors/d[:,None]).reshape(400,4,len(genes)).mean(axis=1)
                    depths=np.rint(d.reshape(400,4).mean(axis=1)).astype(np.int64)
                    x=round_template(template,depths,seed_for(context+':'+alias+':pool_round'))
                if x.shape!=(400,len(genes)):raise ValueError('Incorrect group shape')
                m,report=distribution(x,comparison);report.update(context=context,arm=arm,alias=alias)
                distribution_reports.append(report)
                gene_frames.append(pd.DataFrame(dict(context=context,arm=arm,alias=alias,gene=genes,
                    mean_cpm=m['mean'],variance_cpm=m['variance'],detection=m['detection'],
                    comparison_mean_cpm=comparison['mean'],comparison_variance_cpm=comparison['variance'],comparison_detection=comparison['detection'])))
                writer.append(x,context,alias)
                audits.append(dict(context=context,target=alias,nnz=int(np.count_nonzero(x)),max_library=int(m['depth'].max())))
            writer.close()
        except BaseException:
            writer.f.close();raise
        pending.rename(path)
        receipt=dict(arm=arm,purpose=PROTOCOL['purpose'],arguments=dict(genes=str(folder/'gene_names.csv'),
            targets=str(folder/'pert_counts.csv'),context=[context+'='+role_paths['comparison']['path']],cells=400,target_limit=None),
            sources=[],targets=ALIASES,output_genes=len(genes),cells=writer.n,stored_entries=writer.nnz,moments=audits,
            bytes=path.stat().st_size,output_sha256=sha256(path),cohort_freeze_sha256=sha256(OUT/'cohort_freeze.json'),
            comparison_sha256=role_paths['comparison']['sha256'],identity_validation='pending' if arm=='donor_replay' else 'not_applicable')
        write_json(path.with_suffix('.receipt.json'),receipt)
        validation=validate(path,path.with_suffix('.receipt.json'))
        write_json(out/'validation.json',validation)
        manifests[arm]={'prediction':str(path),'sha256':receipt['output_sha256']}
        print(json.dumps({'context':context,'prepared':arm,'validated':True}),flush=True)
    source.file.close()
    pd.concat(gene_frames,ignore_index=True).to_parquet(folder/'gene_distributions.parquet',index=False)
    write_json(folder/'distribution_summary.json',distribution_reports)
    write_json(folder/'manifest.json',dict(context=context,roles=role_paths,arms=manifests,
        cohort_freeze_sha256=sha256(OUT/'cohort_freeze.json'),donor_pools_sha256=sha256(folder/'donor_pools.csv'),
        donor_replay_rows_sha256=sha256(folder/'donor_replay_rows.csv'),gene_distributions_sha256=sha256(folder/'gene_distributions.parquet'),
        distribution_summary_sha256=sha256(folder/'distribution_summary.json'),all_roles_disjoint=True,exact_ntc_quotas=True))


def run_de(context,arm):
    verify_freeze();folder=OUT/context;manifest=read_json(folder/'manifest.json')
    out=folder/(arm+'_de')
    if out.exists():raise FileExistsError(out)
    paths=[Path(manifest['arms'][arm]['prediction']),Path(manifest['roles']['comparison']['path'])]
    hashes=[manifest['arms'][arm]['sha256'],manifest['roles']['comparison']['sha256']]
    if any(sha256(p)!=h for p,h in zip(paths,hashes)):raise ValueError('Null counts changed')
    pred,control=[ad.read_h5ad(p) for p in paths]
    if not pred.var_names.equals(control.var_names):raise ValueError('Gene axes differ')
    counts=sparse.vstack([pred.X,control.X],format='csr');counts,_=validate_counts(counts)
    labels=pred.obs.target_gene.astype(str).tolist()+['non-targeting']*control.n_obs
    obj=ad.AnnData(counts,obs=pd.DataFrame({'target_gene':pd.Categorical(labels)}),var=pd.DataFrame(index=pred.var_names))
    out.mkdir();started=time.time()
    write_json(out/'preflight.json',dict(context=context,arm=arm,input_hashes=hashes,shape=list(obj.shape),
        cohort_freeze_sha256=sha256(OUT/'cohort_freeze.json'),helper_path=inspect.getsourcefile(compute_de),
        helper_sha256=sha256(inspect.getsourcefile(compute_de)),parameters=PROTOCOL['de']))
    de=compute_de(obj,**PROTOCOL['de'])
    de.write_parquet(out/'de.parquet')
    calls=[]
    for alias in ALIASES:
        sub=de.filter(pl.col('target')==alias)
        calls.append(dict(context=context,arm=arm,alias=alias,tested_genes=sub.height,
            significant_calls=sub.filter(pl.col('p_adj')<.05).height))
    write_json(out/'result.json',dict(context=context,arm=arm,complete=True,calls=calls,elapsed_seconds=time.time()-started,
        de_sha256=sha256(out/'de.parquet'),helper='Exact pinned cell_eval2.de_compute.compute_de',no_response_or_reserve_access=True))
    print(json.dumps({'context':context,'arm':arm,'DE_complete':True,'seconds':time.time()-started}),flush=True)


def identity(context):
    verify_freeze();folder=OUT/context;manifest=read_json(folder/'manifest.json');out=folder/'identity.json'
    if out.exists():raise FileExistsError(out)
    emitter=importlib.import_module('control_preserving_emitter')
    path=Path(manifest['roles']['donor_population']['path'])
    if sha256(path)!=manifest['roles']['donor_population']['sha256']:raise ValueError('Donor population changed')
    obj=ad.read_h5ad(path);raw,depth=validate_counts(obj.X);dense=raw.toarray()
    p=dense/depth[:,None];mean=p.mean(axis=0);bulk=dense.sum(axis=0)/dense.sum()
    expected,audit=emitter.fit_population(dense,mean,bulk)
    if not np.array_equal(expected,p):raise ValueError('Zero-effect population calibration is not exact identity')
    rows=manifest['roles']['donor_population']['source_rows'];lookup={r:i for i,r in enumerate(rows)}
    chosen=pd.read_csv(folder/'donor_replay_rows.csv')
    replay=ad.read_h5ad(manifest['arms']['donor_replay']['prediction']);checks=[]
    for alias in ALIASES:
        selected=np.asarray([lookup[int(r)] for r in chosen[chosen.alias==alias].sort_values('row').source_row])
        x=emitter.emit_counts(dense[selected],expected[selected],seed=seed_for(context+':'+alias+':identity'))
        reference=replay.X[replay.obs.target_gene.astype(str).to_numpy()==alias].toarray()
        if not np.array_equal(x,dense[selected]) or not np.array_equal(x,reference):raise ValueError('New zero-effect emitter differs from exact donor replay')
        checks.append(dict(alias=alias,rows=400,every_count_exact=True))
    write_json(out,dict(context=context,pass_exact_identity=True,population_rows=4000,genes=obj.n_vars,
        expected_population_probabilities_bitexact=True,checks=checks,emitter_path=emitter.__file__,emitter_sha256=sha256(emitter.__file__),
        replay_sha256=manifest['arms']['donor_replay']['sha256'],population_sha256=manifest['roles']['donor_population']['sha256'],
        population_audit=audit,DE_reused_only_because_all_counts_exact=True))
    print(json.dumps({'context':context,'identity_pass':True,'groups':10}),flush=True)


def identity_h1():
    """Prove identity on the existing H1 donor-replay artifact, without truth."""
    verify_freeze()
    out=REPORTS/'control_preserving_null_h1_identity.json'
    if out.exists():raise FileExistsError(out)
    path=ROOT/'outputs/emission_null_400_donor_replay/prediction.h5ad'
    receipt=read_json(path.with_suffix('.receipt.json'))
    if sha256(path)!=receipt['output_sha256']:raise ValueError('Existing H1 donor replay changed')
    obj=ad.read_h5ad(path);raw,depth=validate_counts(obj.X);counts=raw.toarray()
    emitter=importlib.import_module('control_preserving_emitter')
    p=counts/depth[:,None]
    expected,audit=emitter.fit_population(counts,p.mean(axis=0),counts.sum(axis=0)/counts.sum())
    if not np.array_equal(expected,p):raise ValueError('H1 identity population profile changed')
    actual=emitter.emit_counts(counts,expected,seed=seed_for('H1:existing_replay:identity'))
    if not np.array_equal(actual,counts):raise ValueError('H1 new-emitter zero effect differs from saved donor replay')
    write_json(out,dict(context='H1',pass_exact_identity=True,shape=list(obj.shape),every_count_exact=True,
        artifact=str(path),artifact_sha256=receipt['output_sha256'],population_audit=audit,
        emitter_sha256=sha256(emitter.__file__),emitter_path=emitter.__file__,
        meaning='The saved replay artifact itself serves as the4000-row population. This proves zero-effect identity and permits reuse of its already matched DE result; no perturbed truth was read.'))
    print(json.dumps({'context':'H1','identity_pass':True,'shape':list(obj.shape)}),flush=True)


def summarize():
    verify_freeze();groups=[];inputs={};identities={}
    for context in CONTEXTS:
        folder=OUT/context;manifest=read_json(folder/'manifest.json');dist=read_json(folder/'distribution_summary.json')
        dm={(x['arm'],x['alias']):x for x in dist}
        for arm in ARMS:
            p=folder/(arm+'_de')/'result.json';r=read_json(p)
            if not r['complete'] or sha256(p.with_name('de.parquet'))!=r['de_sha256']:raise ValueError('Incomplete or changed DE result')
            groups.extend(dict(dm[(arm,c['alias'])],**{k:v for k,v in c.items() if k not in ('context','arm','alias')}) for c in r['calls'])
        ip=folder/'identity.json';identities[context]=read_json(ip) if ip.exists() else {'status':'pending','pass_exact_identity':False}
        inputs[context]=dict(manifest_sha256=sha256(folder/'manifest.json'),source_sha256=read_json(OUT/'cohort_freeze.json')['contexts'][context]['source_sha256'])
    aggregate=[]
    for context in CONTEXTS:
        for arm in ARMS:
            g=[x for x in groups if x['context']==context and x['arm']==arm]
            a=dict(context=context,arm=arm,aliases=len(g))
            for key in ('significant_calls','median_log2_variance_ratio','mean_abs_log2_variance_error','detection_mae','detection_bias','depth_sd','depth_ks_to_comparison'):
                v=np.asarray([x[key] for x in g]);a[key]={'min':float(v.min()),'median':float(np.median(v)),'mean':float(v.mean()),'max':float(v.max())}
            aggregate.append(a)
    report=dict(created_utc=datetime.now(timezone.utc).isoformat(),purpose=PROTOCOL['purpose'],protocol=PROTOCOL,
        inputs=inputs,cohort_freeze_sha256=sha256(OUT/'cohort_freeze.json'),identity=identities,
        all_new_emitter_identity_checks_pass=all(x.get('pass_exact_identity',False) for x in identities.values()),
        aggregates=aggregate,groups=groups,limitations=['Technical aliases share donor populations/comparison controls; no independent-replicate inference.',
        'Within-guide depth pooling intentionally differs from prior H1 global depth pooling.',
        'Control-only results do not establish perturbed-response performance or a leaderboard improvement.'])
    write_json(REPORTS/'control_preserving_null_results.json',report)
    lines=['# Control-preserving null audit on supplied A/B/C controls','',
        'Each context uses disjoint 4,000-cell comparison, measured-null and donor populations, with exact ntc_id proportions. There are ten technical aliases of 400 cells per arm. All source cells are non-targeting controls; no hidden perturbed responses or reserved targets were accessed.','',
        'DE uses the pinned official scanpy helper: CPM normalization, arithmetic LFC, reference-control CPM>5 gene gate, tie-corrected Wilcoxon and per-alias BH; calls require q<0.05. Aliases are absent from the gene axis. These are diagnostic call counts, not a competition score.','',
        '| Context | Arm | Significant calls median (range) | Median log2 variance ratio | Detection MAE |','| --- | --- | ---: | ---: | ---: |']
    for a in aggregate:
        c=a['significant_calls'];lines.append(f"| {a['context']} | {a['arm']} | {c['median']:.0f} ({c['min']:.0f}–{c['max']:.0f}) | {a['median_log2_variance_ratio']['mean']:.4f} | {a['detection_mae']['mean']:.4f} |")
    lines+=['','Variance and detection compare each group with its disjoint comparison controls on the same comparison-mean-CPM>5 gate. The displayed variance statistic averages the per-alias median log2 variance ratio; detailed depth and gene-level records are retained.','',
        'New-emitter zero-effect identity: '+('PASS in all three contexts, including every count in all30replay groups. Reusing donor-replay DE is justified by exact count identity.' if report['all_new_emitter_identity_checks_pass'] else 'PENDING; donor replay is a measured raw-count comparator and is not yet certified as the new emitter.'),'',
        'Pools contain four cells from the same ntc_id and are depth-sorted within guide. This removes guide-mixing as an explanation but differs deliberately from the older H1 global-depth pooling experiment. Shared donors/comparison controls make aliases correlated; they are not independent biological replicates. Passing this null test alone does not establish a better perturbation model.','',
        'Full protocol, frozen hashes, group results and identity receipts: [JSON report](control_preserving_null_results.json).']
    path=REPORTS/'control_preserving_null_results.md'
    if path.exists():raise FileExistsError(path)
    path.write_text('\n'.join(lines)+'\n')
    print(json.dumps({'complete':True,'identity':report['all_new_emitter_identity_checks_pass']}),flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('stage',choices=['freeze','prepare','de','identity','identity-h1','summarize'])
    p.add_argument('--context',choices=CONTEXTS);p.add_argument('--arm',choices=ARMS)
    a=p.parse_args()
    if a.stage in ('prepare','de','identity') and a.context is None:p.error('--context required')
    if a.stage=='de' and a.arm is None:p.error('--arm required')
    if a.stage=='freeze':freeze()
    elif a.stage=='prepare':prepare(a.context)
    elif a.stage=='de':run_de(a.context,a.arm)
    elif a.stage=='identity':identity(a.context)
    elif a.stage=='identity-h1':identity_h1()
    else:summarize()


if __name__=='__main__':main()
