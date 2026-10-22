#!/usr/bin/env python3
"""Validate and supervise exactly two frozen K562 confirmation metric runs."""
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:
    os.environ[key]='1'
from datetime import datetime,timezone
import hashlib,json,subprocess,time
from pathlib import Path
from validate_prediction_receipt import validate

ROOT=Path(__file__).resolve().parents[1]
EVAL=ROOT/'evaluations/k562_confirmation'
STATE=ROOT/'reports/k562_confirmation_job_state.json'


def stamp():return datetime.now(timezone.utc).isoformat()


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as handle:
        for b in iter(lambda:handle.read(8*1024**2),b''):h.update(b)
    return h.hexdigest()


def save(state):
    state['updated_utc']=stamp();p=STATE.with_suffix('.partial.json')
    p.write_text(json.dumps(state,indent=2)+'\n');p.replace(STATE)


def main():
    if STATE.exists():raise FileExistsError('Confirmation already started; inspect its existing process and outputs')
    model_path=ROOT/'models/m4_multisource_dual/model.json'
    frozen=json.loads(model_path.read_text());started=stamp()
    if frozen['selected_at_utc']>=started:raise ValueError('Model was not frozen before this confirmation')
    snapshot=ROOT/'reports/k562_confirmation_frozen_model.json'
    snapshot.write_text(json.dumps(frozen,indent=2)+'\n')
    cohort=json.loads((EVAL/'manifest.json').read_text())
    if sha(EVAL/'reference.h5ad')!=cohort['reference_sha256'] or sha(EVAL/'context_K562.h5ad')!=cohort['control_sha256']:
        raise ValueError('Frozen reference or controls changed')
    state=dict(started_utc=started,pid=os.getpid(),status='validating',model_snapshot=str(snapshot),
        model_snapshot_sha256=sha(snapshot),original_model_sha256=sha(model_path),
        model_frozen_at_utc=frozen['selected_at_utc'],reference_sha256=cohort['reference_sha256'],
        wrapper_sha256=sha(ROOT/'scripts/run_context_cpu_metrics.sh'),jobs={})
    save(state)
    allowed={str(Path(s['path']).resolve()):(s['sha256'],float(s['weight'])) for s in frozen['sources']
             if not any(t in str(s['path']).lower() for t in ['k562','gwps']) and float(s['weight'])>0}
    for arm in ['selected','neutral']:
        prediction=ROOT/f'outputs/k562_confirmation_{arm}/prediction.h5ad'
        rp=prediction.with_suffix('.receipt.json');r=json.loads(rp.read_text());args=r['arguments']
        result=validate(prediction,rp,excluded_source_tokens=['K562','GWPS'])
        for key in ['mode','amplitude_cpm','amplitude_bulk','pool','prior_counts','common_subtract','seed','agreement','promoters']:
            if args.get(key)!=frozen['parameters'].get(key):raise ValueError(f'{arm}: frozen parameter changed: {key}')
        if args['cells']!=200:raise ValueError('Confirmation requires 200 predicted cells per target')
        actual={str(Path(s['path']).resolve()):(s['sha256'],float(s['weight'])) for s in r['sources'] if float(s['weight'])>0}
        if arm=='selected' and actual!=allowed:raise ValueError('Selected confirmation source configuration differs from frozen source-exclusion protocol')
        if arm=='neutral' and actual:raise ValueError('Neutral comparator has nonzero response-source weights')
        for source in r['sources']:
            if sha(source['path'])!=source['sha256']:raise ValueError('Source bank hash changed')
        validation_path=prediction.with_name('confirmation_validation.json')
        validation_path.write_text(json.dumps(result,indent=2)+'\n')
        out=EVAL/f'{arm}_cpu'
        if out.exists() and any(out.iterdir()):raise FileExistsError(out)
        state['jobs'][arm]=dict(status='validated',prediction=str(prediction),prediction_sha256=r['output_sha256'],
            prediction_receipt_sha256=sha(rp),validation=result,validation_path=str(validation_path),
            output=str(out),same_target_coverage=r['same_target_coverage'])
        save(state)
    children=[]
    for arm in ['selected','neutral']:
        job=state['jobs'][arm];out=Path(job['output']);logpath=out.with_suffix('.log');log=logpath.open('w')
        command=['bash',str(ROOT/'scripts/run_context_cpu_metrics.sh'),job['prediction'],str(EVAL/'reference.h5ad'),str(out)]
        child=subprocess.Popen(command,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
        children.append((arm,child,log));job.update(status='running',pid=child.pid,command=command,log=str(logpath),started_utc=stamp())
        state['status']='running';save(state);print('Started',arm,'pid',child.pid,flush=True)
    failed=[]
    while children:
        for item in list(children):
            arm,child,log=item;code=child.poll()
            if code is None:continue
            log.close();children.remove(item)
            state['jobs'][arm].update(status='complete' if code==0 else 'failed',exit_code=code,finished_utc=stamp())
            if code==0:
                output=Path(state['jobs'][arm]['output'])
                for name in ['agg_results.csv','results.csv','run_meta.json','local_scoring_input.receipt.json']:
                    if not (output/name).exists():raise FileNotFoundError(output/name)
                state['jobs'][arm]['aggregate_sha256']=sha(output/'agg_results.csv')
            else:failed.append(arm)
            save(state);print('Finished',arm,'exit',code,flush=True)
        if children:time.sleep(3)
    state['status']='failed' if failed else 'complete';state['finished_utc']=stamp();save(state)
    if failed:raise RuntimeError(f'Metric runs failed: {failed}')


if __name__=='__main__':main()
