#!/usr/bin/env python3
"""Audit a completed Orion raw source bank and its per-file checksum receipts."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from build_orion_de_bank import ROOT, axes, sha256


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--line',required=True,choices=['HCT116','HEK293T'])
    args=p.parse_args()
    path=ROOT/f'data/processed/{args.line}_full_statistics.npz'
    report_path=ROOT/f'reports/{args.line}_full_statistics.json'
    report=json.loads(report_path.read_text())
    requested, genes, official, h1=axes()
    expected=[r for r in json.loads((ROOT/'research/orion_raw_listing.json').read_text())
              if Path(r['path']).name.startswith(args.line+'_Batch')]
    assert len(report['files'])==len(expected)
    for entry in expected:
        receipt=ROOT/'data/processed/orion_batches'/args.line/(Path(entry['path']).stem+'.json')
        r=json.loads(receipt.read_text())
        assert r['source_sha256']==entry['lfs']['oid'] and r['bytes']==entry['size']
    with np.load(path,allow_pickle=False) as z:
        targets=z['targets'].astype(str); n=z['n_cells']
        assert len(set(targets))==len(targets) and set(targets)<=set(requested)
        assert np.array_equal(z['genes'].astype(str),genes)
        assert (n>0).all() and np.equal(n,np.floor(n)).all()
        for key in ['target_count_sums','target_mean_cpm','matched_control_probability','matched_control_mean_cpm']:
            assert z[key].shape==(len(targets),len(genes)),key
            assert np.isfinite(z[key]).all() and (z[key]>=0).all(),key
        counts=z['target_count_sums']
        assert np.equal(counts,np.floor(counts)).all()
        assert (counts.sum(1)>0).all()
        assert np.allclose(z['matched_control_probability'].sum(1),1,atol=2e-6)
        assert np.isclose(z['global_control_probability'].sum(),1,atol=2e-6)
        assert (z['target_mean_cpm'].sum(1)<=1e6+1).all()
        assert (z['matched_control_mean_cpm'].sum(1)<=1e6+1).all()
        assert z['measured_genes'].shape==(len(genes),)
        assert (counts[:,~z['measured_genes']]==0).all()
        eligible=set(targets[n>=20]); present=set(targets)
        report.update(official_coverage_min20=len(set(official)&eligible),
            h1_coverage_min20=len((set(h1)-{'non-targeting'})&eligible),
            official_targets_missing=sorted(set(official)-present),
            official_targets_below20=sorted(set(official)-eligible),
            target_cells_min=int(n.min()),target_cells_max=int(n.max()),
            output_sha256=sha256(path),verified_source_files=len(expected),
            source_release=json.loads((ROOT/'research/orion_raw_release.json').read_text())['sha'],
            center_scope='575-target official/H1 union; raw banks are uncentered')
        meta=pd.read_parquet(ROOT/'data/raw/orion_gene_metadata.parquet')
        duplicates=meta[meta.gene_name.isin(set(genes))&meta.gene_name.duplicated(keep=False)]
        report['duplicate_symbol_tokens_summed']=duplicates.to_dict(orient='records')
    report['audit']='Passed finite/nonnegative/integer counts, positive cell groups, exact gene axis, source/matched probability sums, CPM mass bounds, unmeasured-gene zeros, and every source checksum receipt'
    pending=report_path.with_suffix('.audited.tmp')
    pending.write_text(json.dumps(report,indent=2)+'\n')
    pending.replace(report_path)
    print(json.dumps({k:v for k,v in report.items() if k!='files'},indent=2))


if __name__=='__main__':
    main()
