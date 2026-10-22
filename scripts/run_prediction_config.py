#!/usr/bin/env python3
"""Translate an explicit model configuration to the audited prediction builder."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[1]


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('config',type=Path)
    ap.add_argument('--dry-run',action='store_true')
    args=ap.parse_args()
    config=args.config.resolve();c=json.loads(config.read_text())
    cmd=[sys.executable,str(ROOT/'scripts/build_prediction.py')]
    sources=c.get('source',[])+c.get('de_source',[])
    for token in c.get('exclude_source_family',[]):
        if any(token.lower() in x['path'].lower() for x in sources):
            raise ValueError(f'Forbidden response source family: {token}')
    for kind in ['source','de_source']:
        for x in c.get(kind,[]):
            p=ROOT/x['path']
            if not p.is_file():
                raise FileNotFoundError(f'Required bank is not complete: {p}')
            cmd += ['--'+kind.replace('_','-'),f'{p}={x["weight"]}']
    for x in c['context']:
        label,p=x.split('=',1);cmd += ['--context',f'{label}={ROOT/p}']
    for k in ['genes','targets','output','promoters','cd4']:
        if c.get(k) is not None:cmd += ['--'+k,str(ROOT/c[k])]
    for k in ['amplitude_cpm','amplitude_bulk','common_subtract','prior_counts',
              'mode','pool','cells','workers','seed','cd4_weight','agreement_min_effect','agreement_fraction']:
        if k in c:cmd += ['--'+k.replace('_','-'),str(c[k])]
    if c.get('agreement'):
        cmd += ['--agreement']
    print(json.dumps({'config':str(config),'command':cmd},indent=2),flush=True)
    if args.dry_run:return
    destination=ROOT/c['output']
    destination.parent.mkdir(parents=True,exist_ok=True)
    snapshot=destination.with_suffix('.config.json')
    if snapshot.exists():raise FileExistsError(snapshot)
    snapshot.write_text(json.dumps({'config':c,'config_sha256':hashlib.sha256(config.read_bytes()).hexdigest(),
                                  'source_config':str(config)},indent=2)+'\n')
    subprocess.run(cmd,cwd=ROOT,check=True)


if __name__=='__main__':main()
