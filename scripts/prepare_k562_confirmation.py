#!/usr/bin/env python3
"""Freeze a second-context confirmation bundle from genuine cached GWPS cells.

Run only after the production K562 bank and receipt exist. Target selection
uses selected-cell counts only; no responses, fitted models, or outcomes enter
the selection. This command neither predicts nor scores.
"""
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:
    os.environ[key]='2'
from pathlib import Path
import hashlib,json,time
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'data/raw/gwps_subset'
SEED=20261020


def sha(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda:handle.read(8*1024**2),b''):digest.update(block)
    return digest.hexdigest()


def main():
    started=time.time(); out=ROOT/'evaluations/k562_confirmation'
    bank=ROOT/'data/K562_GWPS_CPM_full_statistics.npz'
    receipt=bank.with_suffix('.receipt.json')
    if not bank.exists() or not receipt.exists():
        raise FileNotFoundError('Complete and verify the production K562 bank first')
    if (out/'manifest.json').exists() or (out/'reference.h5ad').exists():
        raise FileExistsError('Frozen confirmation already exists; do not resample')
    provenance=json.loads(receipt.read_text())
    if sha(bank)!=provenance['output_sha256']:raise ValueError('Production bank hash mismatch')
    with np.load(bank,allow_pickle=False) as data:
        panel=data['targets'].astype(str); genes=data['genes'].astype(str); counts=data['n_cells']
    eligible=np.asarray(sorted(panel[counts>=100]))
    if len(eligible)<50:raise ValueError('Fewer than 50 eligible targets')
    rng=np.random.default_rng(SEED)
    targets=np.asarray(sorted(rng.choice(eligible,50,replace=False)))
    selection=np.load(BASE/'selection.npz',allow_pickle=False)
    if not np.array_equal(panel,selection['targets']):raise ValueError('Bank/selection target axes differ')
    tid=selection['target_ids']; rows=selection['rows']; batch=selection['batch_ids']
    chosen_ids=np.flatnonzero(np.isin(panel,targets))
    treated=np.flatnonzero(np.isin(tid,chosen_ids))
    controls=np.flatnonzero(tid<0)
    if len(controls)<2000:raise ValueError('Fewer than 2000 genuine cached control cells')
    chosen_controls=np.sort(rng.choice(controls,2000,replace=False))
    chosen=np.sort(np.concatenate([treated,chosen_controls]))
    wanted=np.zeros(len(rows),dtype=bool); wanted[chosen]=True
    source_genes=selection['genes'].astype(str)
    expected_genes=np.asarray(list(dict.fromkeys(source_genes)))
    if not np.array_equal(genes,expected_genes):raise ValueError('Unique source gene axis differs from bank')
    projection=sparse.csr_matrix((np.ones(len(source_genes),dtype=np.int64),
        (np.arange(len(source_genes)),pd.Index(genes).get_indexer(source_genes))),
        shape=(len(source_genes),len(genes)))
    plan=json.loads((BASE/'plan.json').read_text())
    blocks=[]; found=[]; source_count_total=0; used_ranges=0
    for i,p in enumerate(plan['ranges']):
        local=np.flatnonzero(wanted[p['first']:p['last']])+p['first']
        if not len(local):continue
        path=BASE/'rows'/f'{i:06d}.bin'
        if not path.exists() or path.stat().st_size!=p['end']-p['start']:
            raise ValueError(f'Missing cached range {i}')
        raw=np.fromfile(path,dtype='<f4').reshape(-1,len(source_genes))
        x=raw[rows[local]-rows[p['first']]]
        if not np.isfinite(x).all() or (x<0).any() or (x!=np.floor(x)).any():
            raise ValueError(f'Invalid reference raw counts in range {i}')
        source_count_total+=int(x.sum(dtype=np.float64))
        merged=sparse.csr_matrix(x,dtype=np.int64)@projection
        if merged.nnz and merged.data.max()>np.iinfo(np.int32).max:
            raise ValueError('A source count exceeds int32 storage')
        blocks.append(merged.astype(np.int32));found.extend(local.tolist());used_ranges+=1
        if used_ranges%1000==0:print('Read confirmation ranges',used_ranges,'cells',len(found),flush=True)
    np.testing.assert_equal(np.asarray(found),chosen)
    x=sparse.vstack(blocks,format='csr')
    if int(x.sum(dtype=np.int64))!=source_count_total:raise ValueError('Reference count total changed')
    source_targets=np.where(tid[chosen]>=0,panel[np.maximum(tid[chosen],0)],'non-targeting')
    source_batches=selection['batches'][batch[chosen]].astype(str)
    obs=pd.DataFrame(dict(target_gene=source_targets,context='K562',gem_group=source_batches,
        batch=source_batches,source_row=rows[chosen],source_selection_index=chosen),
        index=pd.Index([f'GWPS_row_{i}' for i in rows[chosen]],name=None))
    for column in ['target_gene','context','gem_group','batch']:obs[column]=pd.Categorical(obs[column])
    reference=ad.AnnData(x,obs=obs,var=pd.DataFrame(index=pd.Index(genes,name=None)))
    control_mask=source_targets=='non-targeting';control=reference[control_mask].copy()
    assert control.n_obs==2000 and reference.obs_names.is_unique and reference.var_names.is_unique
    actual=reference.obs.loc[~control_mask].groupby('target_gene',observed=True).size().to_dict()
    if set(actual)!=set(targets) or any(not 100<=n<=200 for n in actual.values()):
        raise ValueError('Reference target eligibility/counts changed')
    out.mkdir(parents=True,exist_ok=True)
    reference.write_h5ad(out/'reference.h5ad',compression='gzip',compression_opts=1)
    control.write_h5ad(out/'context_K562.h5ad',compression='gzip',compression_opts=1)
    pd.DataFrame({'gene_name':genes}).to_csv(out/'gene_names.csv',index=False)
    pd.DataFrame({'target_gene':targets,'n_cells':200}).to_csv(out/'pert_counts.csv',index=False)
    obs.to_csv(out/'source_rows.csv',index_label='reference_cell_id')
    manifest=dict(evaluation_context='K562',split='frozen public GWPS cell-context confirmation',
        seed=SEED,eligible_target_count=len(eligible),eligibility='At least 100 cells in the frozen raw-cell selection',
        target_selection='50 names drawn without replacement from the sorted eligible list, using fixed seed; no response/outcome filtering',
        control_selection='2000 measured non-targeting cells sampled uniformly without replacement from the frozen selected GWPS control pool',
        controls='The identical genuine control rows are both context inputs and reference controls; batch labels retained.',
        targets=targets.tolist(),reference_cells_per_target={str(k):int(v) for k,v in actual.items()},
        prediction_cells_per_target=200,n_reference_treated=int((~control_mask).sum()),n_control_cells=2000,
        n_genes=len(genes),n_reference_batches=len(set(source_batches)),n_control_batches=int(control.obs.gem_group.nunique()),
        reference_count_total=source_count_total,
        source_url=plan['source_url'],source_etag=plan['source_etag'],
        bank=str(bank),bank_sha256=provenance['output_sha256'],bank_receipt=str(receipt),
        selection_sha256=sha(BASE/'selection.npz'),plan_sha256=sha(BASE/'plan.json'),
        source_rows_path=str(out/'source_rows.csv'),source_rows_sha256=sha(out/'source_rows.csv'),
        reference_sha256=sha(out/'reference.h5ad'),control_sha256=sha(out/'context_K562.h5ad'),
        gene_names_sha256=sha(out/'gene_names.csv'),pert_counts_sha256=sha(out/'pert_counts.csv'),
        required_source_exclusion=['All K562/GWPS response banks, including raw, pseudobulk, and DE representations'],
        allowed_prediction_sources='H1, HCT116/HEK293T Orion, and optional CD4; use only K562 controls and perturbation identities for context conditioning',
        model_selection_rule='Freeze parameters using H1 Validation before predicting/scoring this bundle once; do not tune on K562 outcomes.',
        scored=False,predicted=False,all_counts_finite_nonnegative_integer=True,
        unique_gene_axis=True,elapsed_seconds=time.time()-started)
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({k:v for k,v in manifest.items() if k not in ['targets','reference_cells_per_target']},indent=2),flush=True)


if __name__=='__main__':main()
