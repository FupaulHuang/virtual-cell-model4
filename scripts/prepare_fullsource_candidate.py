#!/usr/bin/env python3
"""Prepare and strictly package an explicitly UNSELECTED candidate. Never submit."""
from __future__ import annotations
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import time
from datetime import datetime, timezone

ROOT=Path(__file__).resolve().parents[1]
STATE=ROOT/'reports/fullsource_candidate_job.json'
CONFIG=ROOT/'configs/abc_fullsource_candidate.json'
PRED=ROOT/'outputs/m4_fullsource_candidate/prediction.h5ad'
PACKAGE=ROOT/'submissions/m4_fullsource_candidate.vcc'
CORE=['scripts/build_prediction.py','scripts/run_prediction_config.py',
      'scripts/validate_prediction_receipt.py','external/atlasshift/model.py',
      'external/atlasshift/predict.py']


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
    return h.hexdigest()


def save(state):
    state['updated']=datetime.now(timezone.utc).isoformat()
    tmp=STATE.with_suffix('.partial.json');tmp.write_text(json.dumps(state,indent=2)+'\n');tmp.replace(STATE)


def main():
    if STATE.exists():raise FileExistsError('Candidate preparation already started; inspect existing job')
    state={'status':'waiting_for_verified_K562','selection_status':'unselected_pending_holdout_confirmation',
           'selected_for_submission':False,'submitted':False,'pid':os.getpid(),
           'config':str(CONFIG),'config_sha256':digest(CONFIG),'prediction':str(PRED),'package':str(PACKAGE)}
    save(state)
    try:
        bank=ROOT/'data/K562_GWPS_CPM_full_statistics.npz';receipt=bank.with_suffix('.receipt.json');started=time.monotonic()
        while True:
            try:
                r=json.loads(receipt.read_text())
                if r.get('all_counts_finite_nonnegative_integer') and r.get('output_sha256') and bank.exists():
                    actual=digest(bank)
                    if actual!=r['output_sha256']:raise ValueError('K562 bank checksum mismatch')
                    state['K562_sha256']=actual;break
            except (FileNotFoundError,json.JSONDecodeError):pass
            if time.monotonic()-started>3600:raise TimeoutError('K562 verification deadline reached')
            time.sleep(3)
        state['frozen_core_sha256']={p:digest(ROOT/p) for p in CORE}
        state['status']='generating';save(state)
        print('Verified K562; generating provisional ABC candidate',flush=True)
        with (ROOT/'logs/m4_fullsource_candidate_generation.log').open('w') as log:
            subprocess.run([sys.executable,str(ROOT/'scripts/run_prediction_config.py'),str(CONFIG)],cwd=ROOT,
                           stdout=log,stderr=subprocess.STDOUT,check=True)
        if state['frozen_core_sha256']!={p:digest(ROOT/p) for p in CORE}:
            raise ValueError('Core code changed during generation; candidate must be reviewed')
        state['status']='validating';save(state)
        with PRED.with_name('validation.json').open('w') as log:
            subprocess.run([sys.executable,str(ROOT/'scripts/validate_prediction_receipt.py'),str(PRED)],cwd=ROOT,
                           stdout=log,stderr=subprocess.STDOUT,check=True)
        validation=json.loads(PRED.with_name('validation.json').read_text())
        if not validation.get('valid') or validation.get('shape')!=[360000,18533] or validation.get('groups')!=900:
            raise ValueError('Full official-panel validation did not pass')
        state['independent_validation']=validation
        tmpdir=ROOT/'scratch/vcc_fullsource_candidate';tmpdir.mkdir(parents=True,exist_ok=True)
        available=shutil.disk_usage(tmpdir).free
        if available<60*1024**3:raise RuntimeError('Packaging scratch has less than60GiB free')
        state['packaging_scratch']=str(tmpdir);state['packaging_scratch_free_bytes']=available
        PACKAGE.parent.mkdir(parents=True,exist_ok=True)
        if PACKAGE.exists():raise FileExistsError(PACKAGE)
        vcc=shutil.which('vcc')
        if not vcc:raise FileNotFoundError('Installed production vcc CLI missing')
        state['vcc_version']=subprocess.check_output([vcc,'--version'],text=True).strip()
        command=[vcc,'prep',str(PRED),'-g',str(ROOT/'data/gene_names.csv'),'--perts',str(ROOT/'data/pert_counts.csv'),
                 '-o',str(PACKAGE),'--contexts','A,B,C','--cells-per-pert','400','--encoding','32',
                 '--require-counts','--reject-controls','--verify-targets','--check-cell-counts',
                 '--expected-gene-dim','18533','--max-cell-dim','400000','--max-nnz','4750000000',
                 '--max-counts-per-cell','1000000','--json']
        state['prep_command']=command;state['status']='packaging';save(state)
        env=os.environ.copy();env['TMPDIR']=str(tmpdir)
        with (ROOT/'reports/m4_fullsource_candidate_prep.json').open('w') as out, (ROOT/'logs/m4_fullsource_candidate_prep.log').open('w') as err:
            subprocess.run(command,cwd=ROOT,env=env,stdout=out,stderr=err,check=True)
        with tarfile.open(PACKAGE,'r:') as archive:
            member=archive.getmember('pred.h5ad.zst')
            if not member.isfile() or member.size<=0:raise ValueError('Invalid package payload')
            state['compressed_payload_bytes']=member.size
        state['package_sha256']=digest(PACKAGE);state['package_bytes']=PACKAGE.stat().st_size
        PACKAGE.with_suffix('.vcc.sha256').write_text(f"{state['package_sha256']}  {PACKAGE.name}\n")
        state['prediction_sha256']=json.loads(PRED.with_suffix('.receipt.json').read_text())['output_sha256']
        state['status']='prepared_unselected_pending_holdout_confirmation';save(state)
        PACKAGE.with_suffix('.manifest.json').write_text(json.dumps(state,indent=2)+'\n')
        print(json.dumps(state,indent=2),flush=True)
    except BaseException as e:
        state['status']='failed';state['error']=repr(e);save(state);raise


if __name__=='__main__':main()
