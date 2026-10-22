#!/usr/bin/env python3
"""Prepare a seeded confirmation reference without scoring/tuning on outcomes."""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:
    os.environ.setdefault(k,'1')
from pathlib import Path
import json
import shutil
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from build_prediction import sha256

ROOT=Path(__file__).resolve().parents[1]


def main():
    out=ROOT/'evaluations/h1_test'
    out.mkdir(parents=True,exist_ok=True)
    if (out/'manifest.json').exists() or (out/'reference.h5ad').exists():
        raise FileExistsError('Test reference already exists; do not resample it')
    source=ROOT / 'data/h1_2025/adata_Test.h5ad'
    a=ad.read_h5ad(source,backed='r')
    try:
        labels=a.obs.target_gene.astype(str).to_numpy()
        targets=np.array(sorted(set(labels)-{'non-targeting'}))
        genes=a.var_names.astype(str).to_numpy()
        if len(targets)!=100 or len(genes)!=18080:
            raise ValueError('Unexpected public Test target/gene axes')
        val=pd.read_csv(ROOT/'evaluations/h1_validation/pert_counts.csv').target_gene.astype(str)
        if set(targets)&set(val):
            raise ValueError('Validation/Test target overlap')
        rng=np.random.default_rng(20261020)
        rows=[];counts={}
        for t in targets:
            ix=np.flatnonzero(labels==t);n=min(200,len(ix));counts[t]=n
            rows.extend(rng.choice(ix,n,replace=False).tolist())
        rows=np.array(sorted(rows),dtype=np.int64)
        blocks=[]
        for start in range(0,len(rows),1024):
            blocks.append(sparse.csr_matrix(a.X[rows[start:start+1024]]))
            print('Test reference rows',min(start+1024,len(rows)),'/',len(rows),flush=True)
        x=sparse.vstack(blocks,format='csr')
        if (x.data<0).any() or not np.isfinite(x.data).all() or (x.data!=np.floor(x.data)).any():
            raise ValueError('Reference is not raw integer counts')
        treated=ad.AnnData(x,obs=a.obs.iloc[rows].copy(),var=a.var.copy())
        treated.obs['context']='H1'
    finally:
        a.file.close()
    ctrl_source=ROOT/'evaluations/h1_validation/context_H1.h5ad'
    ctrl=ad.read_h5ad(ctrl_source)
    if not np.array_equal(ctrl.var_names.astype(str),genes):
        raise ValueError('Reference/control gene axes differ')
    ref=ad.concat([treated,ctrl],join='inner',merge='same')
    ref.write_h5ad(out/'reference.h5ad',compression='gzip',compression_opts=1)
    shutil.copyfile(ctrl_source,out/'context_H1.h5ad')
    pd.DataFrame({'gene_name':genes}).to_csv(out/'gene_names.csv',index=False)
    pd.DataFrame({'target_gene':targets,'n_cells':200}).to_csv(out/'pert_counts.csv',index=False)
    manifest={'split':'public H1 2025 Test','source':str(source),'seed':20261020,
       'reference_sha256':sha256(out/'reference.h5ad'),'control_sha256':sha256(ctrl_source),
       'n_controls':ctrl.n_obs,'n_reference_treated':treated.n_obs,'n_genes':len(genes),
       'targets':targets.tolist(),'reference_cells_per_target':counts,
       'prediction_cells_per_target':200,'scored':False,
       'source_exclusion':'All H1 response sources must be excluded during confirmation.',
       'selection_rule':'Freeze configuration using Validation before scoring Test; do not tune on Test.',
       'controls':'Same measured unperturbed H1 controls as Validation; these are inference inputs.'}
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({k:v for k,v in manifest.items() if k not in ['targets','reference_cells_per_target']},indent=2))


if __name__=='__main__':main()
