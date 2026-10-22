#!/usr/bin/env python3
"""Run pinned AtlasShift aggregation on locally downloaded, verified public data."""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:
    os.environ.setdefault(k,'2')
from pathlib import Path
import argparse
import gzip
import hashlib
import json
import re
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'external/atlasshift'))
import prepare
import numpy as np
import pandas as pd


def union_targets():
    values = pd.read_csv(ROOT/'data/pert_counts.csv').target_gene.astype(str).tolist()
    for split in ['Training','Validation','Test']:
        values += pd.read_csv(ROOT / f'data/h1_2025/pert_counts_{split}.csv').target_gene.astype(str).tolist()
    return list(dict.fromkeys(t for t in values if t != 'non-targeting'))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source', required=True, choices=['k562','cd4','promoters'])
    args = ap.parse_args()
    data = ROOT/'data'; raw = data/'raw'
    verified = []
    def local_download(name, directory):
        spec = prepare.ASSETS[name]; path = directory/spec['filename']
        if path.stat().st_size != spec['bytes']:
            raise ValueError(f'Incomplete file {path}: {path.stat().st_size}/{spec["bytes"]}')
        h = hashlib.sha256()
        with path.open('rb') as f:
            for b in iter(lambda: f.read(8*1024**2), b''):
                h.update(b)
        if h.hexdigest() != spec['sha256']:
            raise ValueError(f'Checksum mismatch: {path}')
        verified.append({'path':str(path), **spec})
        print('Verified', name, h.hexdigest(), flush=True)
        return path
    prepare.download = local_download
    prepare.source_targets = lambda data,raw: union_targets()
    t0 = time.time()
    if args.source == 'k562':
        prepare.prepare_k562(data, raw)
        output = data/'K562_GWPS_CPM_full_statistics.npz'
    elif args.source == 'cd4':
        prepare.prepare_cd4(data, raw)
        output = data/'CD4_DE_statistics.npz'
    else:
        path = local_download('gencode47', raw); rows=[]
        with gzip.open(path,'rt') as f:
            for line in f:
                if line.startswith('#'):continue
                fields=line.rstrip('\n').split('\t')
                if fields[2]!='gene':continue
                attrs=dict(re.findall(r'(\w+) "([^"]*)"',fields[8]))
                if 'gene_name' in attrs:
                    rows.append(dict(gene=attrs['gene_name'],chromosome=fields[0],strand=fields[6],
                                     tss=int(fields[3] if fields[6]=='+' else fields[4])))
        genes=pd.read_csv(data/'gene_names.csv').gene_name.astype(str).tolist()
        hg=pd.read_csv(ROOT / 'data/h1_2025/gene_names.csv')
        genes=list(dict.fromkeys(genes+hg.iloc[:,0].astype(str).tolist()))
        output=data/'official_pairs.csv'
        prepare.pairs(pd.DataFrame(rows),union_targets(),genes).to_csv(output,index=False)
    receipt={'source':args.source,'inputs':verified,'output':str(output),
             'bytes':output.stat().st_size,'elapsed_seconds':time.time()-t0,
             'target_union_count':len(union_targets())}
    output.with_suffix('.receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt),flush=True)


if __name__=='__main__':main()
