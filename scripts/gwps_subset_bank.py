#!/usr/bin/env python3
"""Plan direct-curl GWPS ranges and aggregate the exact selected raw cells.

Python never performs network requests. Run emitted configurations with approved
direct ``curl -L --parallel --config ...`` commands. The existing full download
is read-only and optional; only complete selected ranges are aggregated.
"""
from pathlib import Path
import argparse, hashlib, json, os, time
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:
    os.environ[key] = '2'
import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'data/raw/gwps_subset'
URL = 'https://ndownloader.figshare.com/files/35775507'
ETAG = '"e85ededd1213eae9405a2f2422edef5d-252"'


def targets():
    values = pd.read_csv(ROOT/'data/pert_counts.csv').target_gene.astype(str).tolist()
    for split in ['Training','Validation','Test']:
        values += pd.read_csv(ROOT / f'data/h1_2025/pert_counts_{split}.csv').target_gene.astype(str).tolist()
    result = list(dict.fromkeys(t for t in values if t != 'non-targeting'))
    assert len(result) == 575
    return np.asarray(result)


def plan(args):
    meta = np.load(BASE/'metadata.npz')
    labels, batch = meta['labels'], meta['batches']
    panel = targets()
    lookup = {t:i for i,t in enumerate(panel)}
    # Group once; avoid 575 repeated scans of the two-million-row label axis.
    frame = pd.DataFrame({'target':labels, 'batch':batch})
    groups = frame.groupby('target',sort=False).indices
    rng = np.random.default_rng(args.seed)
    selected = []
    for target in panel:
        rows = groups.get(target, np.asarray([],dtype=np.int64))
        if len(rows): selected.extend(rng.choice(rows,min(args.cells_per_target,len(rows)),replace=False))
    control = frame[frame.target == 'non-targeting']
    for _, group in control.groupby('batch',sort=True):
        rows = group.index.to_numpy()
        selected.extend(rng.choice(rows,min(args.controls_per_batch,len(rows)),replace=False))
    rows = np.asarray(sorted(selected),dtype=np.int64)
    assert len(rows) == len(np.unique(rows))
    batch_names = np.asarray(sorted(set(batch[rows])))
    batch_ids = pd.Index(batch_names).get_indexer(batch[rows])
    target_ids = np.asarray([lookup.get(t,-1) for t in labels[rows]],dtype=np.int32)
    offset = int(meta['x_offset']); columns = len(meta['genes']); rowbytes = columns * 4
    plans = []
    for k,row in enumerate(rows):
        start = offset + int(row)*rowbytes; end = start+rowbytes
        if plans and start <= plans[-1]['end'] + rowbytes and end-plans[-1]['start'] <= 2*1024**2:
            plans[-1]['end'] = end; plans[-1]['last'] = k+1
        else:
            plans.append(dict(start=start,end=end,first=k,last=k+1))
    (BASE/'rows').mkdir(exist_ok=True)
    np.savez_compressed(BASE/'selection.npz',rows=rows,targets=panel,target_ids=target_ids,
        batch_ids=batch_ids,batches=batch_names,genes=meta['genes'],x_offset=offset,x_columns=columns)
    payload = dict(source_url=URL,source_etag=ETAG,seed=args.seed,cells_per_target=args.cells_per_target,
        controls_per_batch=args.controls_per_batch,selected_cells=len(rows),
        covered_targets=int(len(set(target_ids[target_ids>=0]))),
        selected_controls=int(sum(target_ids<0)),
        bytes=sum(p['end']-p['start'] for p in plans),ranges=plans)
    (BASE/'plan.json').write_text(json.dumps(payload))
    print(json.dumps({k:v for k,v in payload.items() if k!='ranges'},indent=2))
    print('ranges',len(plans),flush=True)


def config(args):
    payload=json.loads((BASE/'plan.json').read_text()); plans=payload['ranges']; lines=[];todo=[]
    request_tag=str(time.time_ns())
    prefix=ROOT/'data/raw/K562_gwps_raw_singlecell_01.h5ad'
    prefix_size=prefix.stat().st_size if prefix.exists() else 0
    reused=0
    end=len(plans) if args.limit is None else min(len(plans),args.start+args.limit)
    for i in range(args.start,end):
        p=plans[i]; path=BASE/'rows'/f'{i:06d}.bin'
        if path.exists() and path.stat().st_size==p['end']-p['start']: continue
        if p['end']<=prefix_size:
            with prefix.open('rb') as handle:
                handle.seek(p['start']); value=handle.read(p['end']-p['start'])
            if len(value)!=p['end']-p['start']:raise ValueError('Full-download prefix changed')
            path.write_bytes(value);reused+=1
            continue
        lines += [f'url = "{URL}?subset={request_tag}_{i}"','location','fail','silent','show-error',
            'retry = 2','connect-timeout = 15','max-time = 60',
            'header = "Accept-Encoding: identity"',
            f'header = "If-Match: \\{ETAG[0]}{ETAG[1:-1]}\\{ETAG[-1]}"',
            f'range = "{p["start"]}-{p["end"]-1}"',f'output = "{path}"',
            f'write-out = "%output{{>>{ROOT}/logs/gwps_subset_http.log}}%{{http_code}} %{{size_download}} %{{filename_effective}}\\n"','next']
        todo.append(i)
    Path(args.output).write_text('\n'.join(lines[:-1])+'\n' if lines else '')
    print(json.dumps(dict(config=args.output,ranges=len(todo),reused_prefix_ranges=reused,bytes=sum(plans[i]['end']-plans[i]['start'] for i in todo))))


def build(args):
    started=time.time(); plan=json.loads((BASE/'plan.json').read_text()); sel=np.load(BASE/'selection.npz')
    rows=sel['rows']; target_ids=sel['target_ids']; batch_ids=sel['batch_ids']; panel=sel['targets']
    source_genes=sel['genes']; genes=np.asarray(list(dict.fromkeys(source_genes)))
    projection=sparse.csr_matrix((np.ones(len(source_genes),dtype=np.int64),
        (np.arange(len(source_genes)),pd.Index(genes).get_indexer(source_genes))),shape=(len(source_genes),len(genes)))
    nt=len(panel); nb=len(sel['batches']); ng=len(genes); shape=(nt+nb,ng)
    sums=np.zeros(shape,dtype=np.int64); cpm=np.zeros(shape); squares=np.zeros(shape)
    n=np.zeros(nt+nb,dtype=np.int64); bn=np.zeros((nt,nb),dtype=np.int64); bl=np.zeros((nt,nb),dtype=np.int64)
    matrix_total=0; digests=[]
    for i,p in enumerate(plan['ranges']):
        path=BASE/'rows'/f'{i:06d}.bin'
        if args.wait_for_download:
            wait_started=time.time()
            while not path.exists() or path.stat().st_size!=p['end']-p['start']:
                if time.time()-wait_started>600:
                    raise TimeoutError(f'No complete download for range {i} after ten minutes')
                time.sleep(2)
        if not path.exists() or path.stat().st_size!=p['end']-p['start']:
            raise ValueError(f'Missing or truncated range {i}: {path}')
        raw=path.read_bytes(); digests.append(hashlib.sha256(raw).digest())
        k=np.arange(p['first'],p['last']); relative=int(rows[k[0]])
        block=np.frombuffer(raw,dtype='<f4').reshape(-1,len(source_genes))[rows[k]-relative]
        if not np.isfinite(block).all() or (block<0).any() or (block!=np.floor(block)).any():
            raise ValueError(f'Invalid integer raw counts in range {i}')
        x=sparse.csr_matrix(block,dtype=np.int64)@projection
        depth=np.asarray(x.sum(axis=1)).ravel()
        if (depth<=0).any():raise ValueError('Zero-depth cell')
        matrix_total+=int(depth.sum())
        group=np.where(target_ids[k]>=0,target_ids[k],nt+batch_ids[k])
        gg,inv=np.unique(group,return_inverse=True)
        member=sparse.csr_matrix((np.ones(len(k),dtype=np.int64),(inv,np.arange(len(k)))),shape=(len(gg),len(k)))
        sums[gg]+=(member@x).toarray()
        q=x.multiply((1e6/depth)[:,None]).tocsr(); cpm[gg]+=(member@q).toarray(); q.data**=2
        squares[gg]+=(member@q).toarray(); n+=np.bincount(group,minlength=nt+nb)
        pert=target_ids[k]>=0
        np.add.at(bn,(target_ids[k][pert],batch_ids[k][pert]),1)
        np.add.at(bl,(target_ids[k][pert],batch_ids[k][pert]),depth[pert])
        if i%5000==0: print('aggregated',i,'/',len(plan['ranges']),flush=True)
    assert int(sums.sum())==matrix_total
    expected_groups=np.where(target_ids>=0,target_ids,nt+batch_ids)
    np.testing.assert_equal(n,np.bincount(expected_groups,minlength=nt+nb))
    assert int(n.sum())==len(rows)
    assert (n[nt:]>0).all()
    np.testing.assert_equal(bn.sum(axis=1),n[:nt])
    np.testing.assert_allclose(cpm[n>0].sum(axis=1)/n[n>0],1e6,rtol=1e-12)
    means=np.divide(cpm,n[:,None],out=np.zeros_like(cpm),where=n[:,None]>0)
    second=np.divide(squares,n[:,None],out=np.zeros_like(squares),where=n[:,None]>0)
    mean_var=np.maximum(second-means**2,0)/np.maximum(n[:,None]-1,1)
    ctrl=sums[nt:]; cp=ctrl/ctrl.sum(axis=1,keepdims=True)
    lw=np.divide(bl,bl.sum(axis=1,keepdims=True),out=np.zeros_like(bl,dtype=float),where=bl.sum(axis=1,keepdims=True)>0)
    cw=np.divide(bn,n[:nt,None],out=np.zeros_like(bn,dtype=float),where=n[:nt,None]>0)
    gp=ctrl.sum(axis=0)/ctrl.sum(); gm=cpm[nt:].sum(axis=0)/n[nt:].sum()
    mp=lw@cp; mm=cw@means[nt:]; mp[n[:nt]==0]=gp;mm[n[:nt]==0]=gm
    np.testing.assert_allclose(mp.sum(axis=1),1,rtol=1e-12)
    np.testing.assert_allclose(mm.sum(axis=1),1e6,rtol=1e-12)
    cv=(cw*cw)@mean_var[nt:]
    output=Path(args.output);output.parent.mkdir(parents=True,exist_ok=True)
    arrays=dict(source=np.asarray('K562_GWPS_raw_subset'),targets=panel,genes=genes,n_cells=n[:nt],
        target_count_sums=sums[:nt],target_mean_cpm=means[:nt].astype('f4'),
        matched_control_probability=mp.astype('f4'),matched_control_mean_cpm=mm.astype('f4'),
        global_control_probability=gp.astype('f4'),global_control_mean_cpm=gm.astype('f4'),
        measured_genes=np.ones(ng,dtype=bool),target_mean_cpm_variance=mean_var[:nt].astype('f4'),
        matched_control_mean_cpm_variance=cv.astype('f4'),delta_cpm_variance=(mean_var[:nt]+cv).astype('f4'),
        arithmetic_cpm_observed=np.asarray(True),batches=sel['batches'],target_batch_n=bn,target_batch_library=bl,
        control_batch_n=n[nt:],control_batch_count_sums=ctrl)
    for name,value in arrays.items():
        if value.dtype.kind in 'fiu' and (not np.isfinite(value).all() or (value<0).any()):
            raise ValueError(f'Invalid bank array: {name}')
    pending=output.with_suffix('.partial.npz');np.savez_compressed(pending,**arrays);pending.replace(output)
    receipt={k:v for k,v in plan.items() if k!='ranges'}
    receipt.update(output=str(output),output_bytes=output.stat().st_size,genes=ng,observed_raw_count_total=matrix_total,
        sampled_cells=len(rows),coverage_at_20_cells=int(sum(n[:nt]>=20)),elapsed_seconds=time.time()-started,
        all_counts_finite_nonnegative_integer=True,range_sha256_digest=hashlib.sha256(b''.join(digests)).hexdigest(),
        source_integrity='Public byte ranges identified by matching ETag and exact length; whole-source SHA256 requires full download.',
        normalization='Arithmetic mean of per-cell CPM over all 8248 measured genes; duplicate symbols sum before aggregation.',
        matching='Target library weights for batch control probabilities; target cell weights for batch control arithmetic CPM.',
        metadata_sha256=hashlib.sha256((BASE/'metadata.npz').read_bytes()).hexdigest(),
        selection_sha256=hashlib.sha256((BASE/'selection.npz').read_bytes()).hexdigest(),
        plan_sha256=hashlib.sha256((BASE/'plan.json').read_bytes()).hexdigest(),
        builder_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        output_sha256=hashlib.sha256(output.read_bytes()).hexdigest())
    output.with_suffix('.receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2),flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__); sub=parser.add_subparsers(dest='action',required=True)
    p=sub.add_parser('plan');p.add_argument('--cells-per-target',type=int,default=200);p.add_argument('--controls-per-batch',type=int,default=50);p.add_argument('--seed',type=int,default=20261020)
    p=sub.add_parser('config');p.add_argument('--start',type=int,default=0);p.add_argument('--limit',type=int);p.add_argument('--output',default=str(BASE/'counts.curl'))
    p=sub.add_parser('build');p.add_argument('--output',default=str(ROOT/'data/K562_GWPS_CPM_full_statistics.npz'));p.add_argument('--wait-for-download',action='store_true',help='Overlap aggregation with the independently live downloader; wait up to ten minutes for each missing range')
    args=parser.parse_args();globals()[args.action](args)


if __name__=='__main__':main()
