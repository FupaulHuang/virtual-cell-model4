#!/usr/bin/env python3
"""Report measured target coverage, source efficacy and basal-context similarity.

These diagnostics do not use hidden challenge responses and are not leaderboard
scores. Efficacy is measured on the targeted gene; this gene is excluded from
official prediction scoring.
"""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]


def summarize(path, official):
    with np.load(path, allow_pickle=False) as z:
        targets = z['targets'].astype(str)
        genes = z['genes'].astype(str)
        n = z['n_cells']
        if n.ndim != 1:
            raise ValueError(f'Expected per-target n_cells: {path}')
        ci = pd.Index(genes).get_indexer(targets)
        raw = 'target_count_sums' in z
        mean = z['target_mean_cpm']
        ctrl = z['matched_control_mean_cpm'] if raw else z['control_mean_cpm']
        if ctrl.ndim == 1:
            ctrl = np.broadcast_to(ctrl, mean.shape)
        if mean.shape != (len(targets), len(genes)) or ctrl.shape != mean.shape:
            raise ValueError('Response mean axes differ')
        mask = np.broadcast_to(z['measured_genes'], mean.shape) if raw else z['measured']
        if mask.shape != mean.shape:
            raise ValueError('Measured-pair mask axes differ')
        available = (n >= 20) & mask.any(1)
        rr = np.flatnonzero((ci >= 0) & available)
        rr = rr[mask[rr, ci[rr]]]
        ratio = np.full(len(targets), np.nan)
        ratio[rr] = (mean[rr, ci[rr]] + 1.) / (ctrl[rr, ci[rr]] + 1.)
        on_panel = np.isin(targets, official)
        rows = pd.DataFrame({'source': Path(path).stem, 'target': targets, 'n_cells': n,
            'available_min20': available, 'in_challenge_panel': on_panel,
            'own_gene_cpm_ratio_pseudocount1': ratio})
        covered = set(targets[on_panel & available])
        finite_ratio = ratio[on_panel & np.isfinite(ratio)]
        report = {'path': str(Path(path).resolve()), 'kind': 'raw' if raw else 'DE-only',
            'source_targets':len(targets), 'source_genes':len(genes),
            'challenge_covered_min20':len(covered), 'challenge_targets_missing':sorted(set(official)-covered),
            'source_cells_min':int(n.min()), 'source_cells_median':float(np.median(n)),
            'challenge_targets_with_own_gene_measurement':len(finite_ratio),
            'challenge_own_gene_ratio_median':float(np.median(finite_ratio)) if len(finite_ratio) else None,
            'challenge_fraction_own_gene_decreased':float(np.mean(finite_ratio<1)) if len(finite_ratio) else None,
            'challenge_fraction_own_gene_ratio_below_half':float(np.mean(finite_ratio<.5)) if len(finite_ratio) else None,
            'caveat':'Pseudocount regularizes low-expression own-gene ratios; efficacy is a source diagnostic, not a performance score.'}
        if raw:
            gp=z['global_control_probability'].astype(float)
            similarities={}
            for c in ['A','B','C']:
                p=ROOT/f'data/legacy_model3/features/m3_e04_context_features/context_{c}.npz'
                # Trusted local Model3 cache stores gene labels as object dtype.
                with np.load(p,allow_pickle=True) as d:
                    ix=pd.Index(genes).get_indexer(d['genes'].astype(str));good=ix>=0
                    src=gp[ix[good]];dst=d['mean_counts'][good].astype(float)
                    src/=src.sum();dst/=dst.sum()
                    similarities[c]={'shared_genes':int(good.sum()),
                        'spearman_basal':float(spearmanr(src,dst).statistic),
                        'pearson_log1p_cp50k':float(np.corrcoef(np.log1p(50000*src),np.log1p(50000*dst))[0,1])}
            report['basal_similarity_unsupervised']=similarities
    return report, rows, covered


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('banks',nargs='+',type=Path)
    ap.add_argument('--output',type=Path,required=True)
    a=ap.parse_args()
    official=pd.read_csv(ROOT/'data/pert_counts.csv').target_gene.astype(str).tolist()
    reports=[];tables=[];union=set()
    for path in a.banks:
        r,t,c=summarize(path,official);reports.append(r);tables.append(t);union |= c
    out={'sources':reports,'union_coverage_min20':len(union),'panel_size':len(official),
         'uncovered':sorted(set(official)-union)}
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(out,indent=2)+'\n')
    pd.concat(tables,ignore_index=True).to_csv(a.output.with_suffix('.targets.csv'),index=False)
    print(json.dumps(out,indent=2))


if __name__=='__main__':main()
