#!/usr/bin/env python3
"""Emit one bounded batch using a freshly obtained public Figshare redirect.

First call direct curl to obtain headers for a one-byte public range. This
command reads that Location locally and emits only S3 range GETs, reducing
Figshare API load. It performs no network access and never prints the signed URL.
"""
import argparse
from datetime import datetime, timezone, timedelta
import json
from pathlib import Path
from urllib.parse import urlparse, parse_qs

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'data/raw/gwps_subset'
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--headers',type=Path,default=BASE/'fresh_headers.txt')
parser.add_argument('--start',type=int,default=0)
parser.add_argument('--limit',type=int,default=256)
parser.add_argument('--verify',action='store_true')
parser.add_argument('--output',type=Path,default=BASE/'signed_batch.curl')
args=parser.parse_args()
if args.verify:
    status=json.loads((BASE/'signed_batch_status.json').read_text())
    plan=json.loads((BASE/'plan.json').read_text())
    missing=[]
    for i in status['selected_indices']:
        p=plan['ranges'][i];path=BASE/'rows'/f'{i:06d}.bin'
        if not path.exists() or path.stat().st_size!=p['end']-p['start']:missing.append(i)
    status['incomplete_ranges']=len(missing)
    if missing:status['next_start']=min(missing)
    (BASE/'signed_batch_status.json').write_text(json.dumps(status))
    print(json.dumps({k:v for k,v in status.items() if k!='selected_indices'}))
    raise SystemExit(0)
url=None
for line in args.headers.read_text().splitlines():
    if line.lower().startswith('location:'):
        url=line.split(':',1)[1].strip()
if url is None:raise SystemExit('No public redirect in headers')
parsed=urlparse(url)
if parsed.scheme!='https' or not parsed.hostname.endswith('.amazonaws.com'):
    raise SystemExit('Unexpected public storage host')
query=parse_qs(parsed.query)
expires=datetime.strptime(query['X-Amz-Date'][0],'%Y%m%dT%H%M%SZ').replace(tzinfo=timezone.utc)+timedelta(seconds=int(query['X-Amz-Expires'][0]))
if (expires-datetime.now(timezone.utc)).total_seconds()<3:
    raise SystemExit('Public redirect expires too soon; obtain a fresh redirect')
plan=json.loads((BASE/'plan.json').read_text())
lines=[]; selected=[]; after=args.start
for i,p in enumerate(plan['ranges'][args.start:],start=args.start):
    after=i+1
    path=BASE/'rows'/f'{i:06d}.bin'
    if path.exists() and path.stat().st_size==p['end']-p['start']:continue
    lines += [f'url = "{url}"','location','fail','silent','show-error','connect-timeout = 5','max-time = 30',
        'header = "Accept-Encoding: identity"',
        f'header = "If-Match: \\"{plan["source_etag"][1:-1]}\\""',
        f'range = "{p["start"]}-{p["end"]-1}"',f'output = "{path}"',
        f'write-out = "%output{{>>{ROOT}/logs/gwps_subset_http.log}}%{{http_code}} %{{size_download}} %{{filename_effective}}\\n"','next']
    selected.append(i)
    if len(selected)>=args.limit:break
args.output.write_text('\n'.join(lines[:-1])+'\n' if lines else '')
status=dict(start=args.start,next_start=after,ranges=len(selected),bytes=sum(plan['ranges'][i]['end']-plan['ranges'][i]['start'] for i in selected),selected_indices=selected)
(BASE/'signed_batch_status.json').write_text(json.dumps(status))
print(json.dumps({k:v for k,v in status.items() if k!='selected_indices'}))
