#!/usr/bin/env python3
"""Execute the prespecified raw-source validation queue after bank audits finish.

This process owns its child jobs and writes their exit codes. Do not edit its
scoring shell script while a child is running.
"""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:
    os.environ.setdefault(k,'1')
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
STATE=ROOT/'reports/raw_comparison_job_state.json'


def stamp():return datetime.now(timezone.utc).isoformat()


def write_state(state):
    state['updated']=stamp()
    tmp=STATE.with_suffix('.partial.json');tmp.write_text(json.dumps(state,indent=2)+'\n');tmp.replace(STATE)


def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
    return h.hexdigest()


def main():
    if STATE.exists():raise FileExistsError('Queue already started; inspect its existing process and outputs')
    state={'started':stamp(),'pid':os.getpid(),'status':'waiting_for_audited_banks','jobs':{}}
    write_state(state);start=time.monotonic()
    for line in ['HCT116','HEK293T']:
        report=ROOT/f'reports/{line}_full_statistics.json'
        while True:
            try:
                r=json.loads(report.read_text())
                if r.get('output_sha256') and r.get('audit'):
                    p=Path(r['output'])
                    if digest(p)!=r['output_sha256']:raise ValueError('Bank hash mismatch')
                    state.setdefault('verified_banks',{})[line]=r['output_sha256'];write_state(state)
                    print('Verified source bank',line,flush=True);break
            except (FileNotFoundError,json.JSONDecodeError):pass
            if time.monotonic()-start>3600:raise TimeoutError('Audited banks were not available within one hour')
            time.sleep(3)
    state['status']='running_comparisons';write_state(state)
    children=[]
    for name in ['raw_dual','raw_dual_half','raw_single']:
        cfg=ROOT/f'configs/h1_{name}.json';prediction=ROOT/f'outputs/h1_{name}/prediction.h5ad'
        out=ROOT/f'evaluations/h1_validation/{name}_cpu'
        state['jobs'][name]={'status':'generating','config':str(cfg),'prediction':str(prediction),'output':str(out)}
        write_state(state)
        with (ROOT/f'logs/h1_{name}.log').open('w') as log:
            subprocess.run([sys.executable,str(ROOT/'scripts/run_prediction_config.py'),str(cfg)],cwd=ROOT,
                           stdout=log,stderr=subprocess.STDOUT,check=True)
        with prediction.with_name('validation.json').open('w') as log:
            subprocess.run([sys.executable,str(ROOT/'scripts/validate_prediction_receipt.py'),str(prediction),
                            '--exclude-source-token','H1'],cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,check=True)
        logpath=out.with_suffix('.log');log=logpath.open('w')
        child=subprocess.Popen(['bash',str(ROOT/'scripts/run_context_cpu_metrics.sh'),str(prediction),
                               str(ROOT/'evaluations/h1_validation/reference.h5ad'),str(out)],cwd=ROOT,
                               stdout=log,stderr=subprocess.STDOUT)
        children.append((name,child,log))
        state['jobs'][name].update(status='scoring',pid=child.pid,log=str(logpath));write_state(state)
        print('Started metrics',name,'pid',child.pid,flush=True)
    failed=[]
    while children:
        for item in list(children):
            name,child,log=item;code=child.poll()
            if code is None:continue
            log.close();children.remove(item)
            state['jobs'][name].update(status='complete' if code==0 else 'failed',exit_code=code,finished=stamp())
            write_state(state);print('Metric job finished',name,code,flush=True)
            if code!=0:failed.append(name)
        if children:time.sleep(3)
    state['status']='failed' if failed else 'complete';write_state(state)
    if failed:raise RuntimeError(f'Metric jobs failed: {failed}')


if __name__=='__main__':main()
