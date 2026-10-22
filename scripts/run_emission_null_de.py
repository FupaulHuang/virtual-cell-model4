#!/usr/bin/env python3
"""Run only official DE-count diagnostics with verified null data and real-DE reuse."""
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(key, '1')
from pathlib import Path
from datetime import datetime, timezone
import argparse
import json
import subprocess
import sys
import anndata as ad
import numpy as np
from scipy import sparse
from build_prediction import ROOT, sha256
from prepare_metric_input import prepare

EVAL = ROOT / 'evaluations/emission_null_400'
ARMS = ('measured_null', 'donor_replay', 'pool4_no_fit', 'neutral_dual', 'response_probe')
CACHE = ROOT / 'reports/emission_null_real_de_cache.json'


def read(path): return json.loads(Path(path).read_text())


def main():
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('arm', choices=ARMS)
    arm = parser.parse_args().arm
    for name, digest in read(ROOT / 'reports/emission_null_freeze.json')['files'].items():
        if sha256(ROOT / name) != digest: raise ValueError('Frozen code/protocol changed: ' + name)
    manifest = read(EVAL / 'manifest.json'); reference = EVAL / 'reference.h5ad'
    if sha256(reference) != manifest['reference_sha256']: raise ValueError('Null reference changed')
    prediction = Path(manifest['arms'][arm]['prediction'])
    validation = read(prediction.parent / 'validation.json')
    if not validation['valid'] or not validation['hash_verified']: raise ValueError('Full validation required')
    if sha256(prediction) != manifest['arms'][arm]['sha256']: raise ValueError('Prediction changed')
    folder = EVAL / f'{arm}_de'; folder.mkdir(parents=True,exist_ok=True)
    if (folder/'agg_results.csv').exists(): raise FileExistsError('Completed diagnostic already exists')
    local = folder / 'local_scoring_input.h5ad'
    if not local.exists(): prepare(prediction,reference,local)
    receipt = read(local.with_suffix('.receipt.json'))
    if (receipt['prediction_sha256'] != sha256(prediction)
            or receipt['reference_sha256'] != manifest['reference_sha256']
            or receipt['shape'] != [8000,18080] or receipt['controls_added'] != 4000):
        raise ValueError('Scoring copy provenance mismatch')
    a,p,r = [ad.read_h5ad(x,backed='r') for x in [local,prediction,reference]]
    try:
        if (a.shape!=(8000,18080) or p.shape!=(4000,18080)
                or not np.array_equal(a.var_names,p.var_names) or not np.array_equal(a.var_names,r.var_names)
                or not np.array_equal(a.obs.target_gene.astype(str).to_numpy()[:4000],p.obs.target_gene.astype(str))
                or not np.all(a.obs.target_gene.astype(str).to_numpy()[4000:]=='non-targeting')):
            raise ValueError('Scoring copy axes differ')
        for start in range(0,4000,512):
            end=min(start+512,4000)
            if (sparse.csr_matrix(a.X[start:end]) != sparse.csr_matrix(p.X[start:end])).nnz:
                raise ValueError('Scored prediction differs')
        rows=np.flatnonzero(r.obs.target_gene.astype(str).to_numpy()=='non-targeting')
        if len(rows)!=4000 or set(r.obs.source_target_gene.astype(str))!={'non-targeting'}:
            raise ValueError('Reference contains real perturbations or wrong controls')
        for start in range(0,4000,512):
            end=min(start+512,4000)
            if (sparse.csr_matrix(a.X[4000+start:4000+end]) != sparse.csr_matrix(r.X[rows[start:end]])).nnz:
                raise ValueError('Scored comparison controls differ')
    finally:
        a.file.close();p.file.close();r.file.close()
    extra=[]
    if arm!='measured_null':
        terminal=read(ROOT/'reports/emission_null_terminal_statuses.json')
        if terminal.get('measured_null',{}).get('exit_code')!=0: raise ValueError('Measured-null DE must complete first')
        cache=read(CACHE); baseline=EVAL/'measured_null_de'
        for key,name in [('real_de_sha256','de_real.parquet'),('run_params_sha256','run_params.yaml'),('run_meta_sha256','run_meta.json')]:
            if sha256(baseline/name)!=cache[key]: raise ValueError('Real-DE cache changed')
        if cache['reference_sha256']!=manifest['reference_sha256'] or cache['entrypoint_sha256']!=sha256(__file__):
            raise ValueError('Cache reference or diagnostic entrypoint differs')
        extra=['--de-real',str(baseline/'de_real.parquet')]
    preflight={'verified_at_utc':datetime.now(timezone.utc).isoformat(),'arm':arm,
               'prediction_sha256':sha256(prediction),'reference_sha256':manifest['reference_sha256'],
               'scored_input_sha256':sha256(local),'actual_scored_counts_verified':True,
               'entrypoint_sha256':sha256(__file__),'diagnostic_only':True}
    (folder/'preflight.json').write_text(json.dumps(preflight,indent=2)+'\n')
    env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='4',MKL_NUM_THREADS='1',NUMBA_NUM_THREADS='4',
                                   NUMBA_CACHE_DIR=str(ROOT/'.cache/numba'),MPLCONFIGDIR=str(ROOT/'.cache/matplotlib'))
    command=[str(ROOT / '.venv/bin/cell-eval2'),'run','-ap',str(local),'-ar',str(reference),
             '--preset','vcc2026','--set','metrics=[de_wilcoxon_nsig_counts_real,de_wilcoxon_nsig_counts_pred]',
             '--pert-col','target_gene','--set','device=cpu','--set','de.backend=scanpy','--set','control_source=real',
             '--set','num_threads=4','--write-degenes',*extra,'-o',str(folder)]
    subprocess.run(command,cwd=ROOT,env=env,check=True)
    if arm=='measured_null':
        if CACHE.exists(): raise FileExistsError(CACHE)
        CACHE.write_text(json.dumps({'reference_sha256':manifest['reference_sha256'],
            'real_de_sha256':sha256(folder/'de_real.parquet'),'run_params_sha256':sha256(folder/'run_params.yaml'),
            'run_meta_sha256':sha256(folder/'run_meta.json'),'entrypoint_sha256':sha256(__file__),
            'purpose':'Same null reference; diagnostic DE parameters identical for every arm'},indent=2)+'\n')


if __name__=='__main__': main()
