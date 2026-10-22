#!/usr/bin/env python3
"""Verify public S3 range parts, assemble separately, and build a masked CD4 bank.

No network calls occur here. The original concurrent full-file download is never
opened or modified. Published quality gates follow the MIT AtlasShift reader.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import re
import time

import h5py
import numpy as np
from anndata.io import read_elem
from build_orion_de_bank import ROOT, axes, sha256


def assemble(plan):
    directory=ROOT/'data/raw/cd4_chunks'
    destination=ROOT/'data/raw/GWCD4i.DE_stats.ranged.h5ad'
    if destination.exists():
        if destination.stat().st_size!=plan['source']['bytes'] or sha256(destination)!=plan['source']['sha256']:
            raise ValueError('Existing separate assembled file has incorrect checksum')
        return destination
    pending=destination.with_suffix('.partial.h5ad')
    digest=hashlib.sha256()
    receipts=[]
    with pending.open('wb') as writer:
        for part in plan['parts']:
            stem=f"{part['index']:06d}"
            data=directory/(stem+'.bin'); header=directory/(stem+'.headers')
            if data.stat().st_size!=part['bytes']:
                raise ValueError(f'Incomplete range part {stem}')
            response=header.read_text()
            statuses=re.findall(r'HTTP/\S+\s+(\d{3})',response)
            ranges=re.findall(r'Content-Range:\s*bytes\s+(\d+)-(\d+)/(\d+)',response,re.I)
            etags=re.findall(r'ETag:\s*"([^"]+)"',response,re.I)
            if not statuses or statuses[-1]!='206' or not ranges or not etags:
                raise ValueError(f'Range response is not complete HTTP 206: {stem}')
            if tuple(map(int,ranges[-1]))!=(part['start'],part['end']-1,plan['source']['bytes']):
                raise ValueError(f'Incorrect Content-Range: {stem}')
            if etags[-1]!=plan['etag']:
                raise ValueError(f'Source ETag changed: {stem}')
            chunk_digest=hashlib.sha256()
            with data.open('rb') as reader:
                for block in iter(lambda:reader.read(8*1024**2),b''):
                    writer.write(block);digest.update(block);chunk_digest.update(block)
            receipts.append({**part,'sha256':chunk_digest.hexdigest(),'http_status':206,'etag':etags[-1]})
    if pending.stat().st_size!=plan['source']['bytes'] or digest.hexdigest()!=plan['source']['sha256']:
        raise ValueError('Whole assembled source SHA256 mismatch')
    pending.replace(destination)
    (ROOT/'reports/cd4_range_acquisition.json').write_text(json.dumps(dict(
        source_url=plan['source']['url'],source_version=plan['source_version'],etag=plan['etag'],
        source_sha256=digest.hexdigest(),bytes=destination.stat().st_size,
        assembled_path=str(destination),range_count=len(receipts),parts=receipts),indent=2)+'\n')
    print('Assembled and full-SHA256 verified',destination,flush=True)
    return destination


def build(path,plan):
    targets,_,official,h1=axes()
    with h5py.File(path,'r') as f:
        obs=read_elem(f['obs']);var=read_elem(f['var'])
        genes=var.gene_name.astype(str).to_numpy(dtype=str) if 'gene_name' in var else var.index.astype(str).to_numpy(dtype=str)
        if len(set(genes))!=len(genes):
            raise ValueError('Duplicate source symbols require explicit resolution')
        selected=np.flatnonzero(obs.target_contrast_gene_name.astype(str).isin(targets))
        chosen=obs.iloc[selected].copy()
        for key in ['single_guide_estimate','ontarget_significant','distal_offtarget_flag','low_target_gex']:
            if chosen[key].isna().any() or not set(chosen[key].unique()) <= {True,False,0,1}:
                raise ValueError(f'Unexpected nonboolean quality flag: {key}')
        layers={}
        for name in ['log_fc','adj_p_value','lfcSE']:
            ds=f['layers'][name]
            if not isinstance(ds,h5py.Dataset) or ds.ndim!=2:
                raise ValueError('Expected dense publisher DE layer')
            # Individual contiguous row reads avoid slow HDF5 fancy indexing.
            layers[name]=np.stack([np.asarray(ds[int(row),:],dtype=np.float32) for row in selected])
            print('Read',name,layers[name].shape,flush=True)
    conditions=sorted(chosen.culture_condition.astype(str).unique())
    shape=(len(conditions),len(targets),len(genes))
    effects=np.zeros(shape,np.float32);q=np.ones(shape,np.float32);se=np.zeros(shape,np.float32)
    available=np.zeros(shape[:2],bool);quality=np.zeros(shape[:2],bool)
    n=np.zeros(shape[:2],np.int64);rows_per=np.zeros(shape[:2],np.int32)
    culture=chosen.culture_condition.astype(str).to_numpy()
    names=chosen.target_contrast_gene_name.astype(str).to_numpy()
    for ci,condition in enumerate(conditions):
        for ti,target in enumerate(targets):
            rows=np.flatnonzero((culture==condition)&(names==target))
            if not len(rows):continue
            metadata=chosen.iloc[rows];values=layers['log_fc'][rows]
            if not np.isfinite(values).all():raise ValueError('Nonfinite source log2 fold changes')
            effects[ci,ti]=values.mean(0);q[ci,ti]=layers['adj_p_value'][rows].max(0);se[ci,ti]=layers['lfcSE'][rows].mean(0)
            available[ci,ti]=True;rows_per[ci,ti]=len(rows);n[ci,ti]=int(metadata.n_cells_target.min())
            quality[ci,ti]=bool(((metadata.n_guides>=2)&~metadata.single_guide_estimate.astype(bool)
                &metadata.ontarget_significant.astype(bool)&~metadata.distal_offtarget_flag.astype(bool)
                &~metadata.low_target_gex.astype(bool)).all())
    destination=ROOT/'data/CD4_DE_statistics.npz';pending=destination.with_suffix('.partial.npz')
    if destination.exists():raise FileExistsError(destination)
    np.savez_compressed(pending,targets=targets,genes=genes,conditions=np.asarray(conditions),log2fc=effects,
        adjusted_p=q,lfcSE=se,available=available,quality_pass=quality,n_cells=n,rows_per_result=rows_per)
    with np.load(pending,allow_pickle=False) as z:
        assert np.array_equal(z['targets'],targets) and np.array_equal(z['genes'],genes)
        assert z['log2fc'].shape==shape and np.isfinite(z['log2fc']).all()
        assert ((z['n_cells']>0)==z['available']).all()
    pending.replace(destination)
    usable=(available&quality&(n>=20)).any(0)
    report=dict(source_url=plan['source']['url'],source_sha256=plan['source']['sha256'],source_etag=plan['etag'],
        source_version=plan['source_version'],source_path=str(path),output=str(destination),output_sha256=sha256(destination),
        targets=len(targets),genes=len(genes),conditions=conditions,selected_source_rows=len(selected),
        official_targets_available=int(available.any(0)[np.isin(targets,official)].sum()),
        h1_targets_available=int(available.any(0)[np.isin(targets,h1)].sum()),
        official_coverage_quality_min20=int(usable[np.isin(targets,official)].sum()),
        h1_coverage_quality_min20=int(usable[np.isin(targets,h1)].sum()),
        quality_pass_pairs=int((available&quality&(n>=20)).sum()),
        quality_gates='>=2 guides, not single-guide estimate, significant on-target, no distal-off-target flag, adequate target expression; minimum20cells applied by model',
        aggregation='Within condition/target, arithmetic mean publisher log2FC and SE, maximum adjustedP, minimum cellcount; all replicates must pass quality gates',
        center_scope='575-target union; effects uncentered',
        audit='Full source SHA256 matched; every range HTTP206/Content-Range/ETag/length verified; finite effects, exact axes, and availability/cell counts checked')
    p=ROOT/'reports/CD4_DE_statistics.json';tmp=p.with_suffix('.partial.json');tmp.write_text(json.dumps(report,indent=2)+'\n');tmp.replace(p)
    print(json.dumps(report,indent=2),flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--wait-seconds',type=int,default=900)
    args=p.parse_args();plan=json.loads((ROOT/'research/cd4_chunks_plan.json').read_text());started=time.monotonic();last=0
    while True:
        incomplete=[r for r in plan['parts'] if not (ROOT/f'data/raw/cd4_chunks/{r["index"]:06d}.bin').exists()
                    or (ROOT/f'data/raw/cd4_chunks/{r["index"]:06d}.bin').stat().st_size!=r['bytes']]
        if not incomplete:break
        if time.monotonic()-started>args.wait_seconds:raise TimeoutError('Optional CD4 acquisition exceeded bounded wait')
        if time.monotonic()-last>60:print('Waiting for range parts:',len(plan['parts'])-len(incomplete),'/',len(plan['parts']),flush=True);last=time.monotonic()
        time.sleep(3)
    path=assemble(plan);build(path,plan)


if __name__=='__main__':main()
