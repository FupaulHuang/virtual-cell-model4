#!/usr/bin/env python3
"""Run the frozen complete-source validation candidate after K562 verifies."""
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:
    os.environ.setdefault(key,'1')
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time
from datetime import datetime, timezone

ROOT=Path(__file__).resolve().parents[1]
STATE=ROOT/'reports/full_source_validation_job.json'


def save(state):
    state['updated']=datetime.now(timezone.utc).isoformat()
    tmp=STATE.with_suffix('.partial.json');tmp.write_text(json.dumps(state,indent=2)+'\n');tmp.replace(STATE)


def main():
    if STATE.exists():raise FileExistsError('Full-source validation already started; inspect its existing handle')
    state={'status':'waiting_for_verified_K562','pid':os.getpid()};save(state)
    bank=ROOT/'data/K562_GWPS_CPM_full_statistics.npz';receipt=bank.with_suffix('.receipt.json')
    started=time.monotonic()
    while True:
        try:
            r=json.loads(receipt.read_text())
            if r.get('all_counts_finite_nonnegative_integer') and r.get('output_sha256') and bank.exists():
                h=hashlib.sha256()
                with bank.open('rb') as f:
                    for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
                if h.hexdigest()!=r['output_sha256']:raise ValueError('K562 bank checksum mismatch')
                state['K562_sha256']=h.hexdigest();break
        except (FileNotFoundError,json.JSONDecodeError):pass
        if time.monotonic()-started>3600:raise TimeoutError('K562 verification deadline reached')
        time.sleep(3)
    state['status']='generating';save(state);print('Verified K562; generating full-source H1 validation',flush=True)
    subprocess.run([sys.executable,str(ROOT/'scripts/run_prediction_config.py'),str(ROOT/'configs/h1_raw_dual_full.json')],cwd=ROOT,check=True)
    pred=ROOT/'outputs/h1_raw_dual_full/prediction.h5ad'
    with pred.with_name('validation.json').open('w') as log:
        subprocess.run([sys.executable,str(ROOT/'scripts/validate_prediction_receipt.py'),str(pred),
                        '--exclude-source-token','H1'],cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,check=True)
    state['status']='scoring';save(state)
    out=ROOT/'evaluations/h1_validation/raw_dual_full_cpu'
    with out.with_suffix('.log').open('w') as log:
        completed=subprocess.run(['bash',str(ROOT/'scripts/run_context_cpu_metrics.sh'),str(pred),
                                 str(ROOT/'evaluations/h1_validation/reference.h5ad'),str(out)],
                                 cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
    state['exit_code']=completed.returncode;state['status']='complete' if completed.returncode==0 else 'failed';save(state)
    completed.check_returncode()


if __name__=='__main__':main()
