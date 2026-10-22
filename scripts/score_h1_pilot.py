#!/usr/bin/env python3
"""Provisional six-component summary against the explicitly one-split pilot.

Official score APIs compute each component; their unweighted diagnostic mean
here is NOT a five-split result or a production challenge score.
"""
import argparse
import json
from pathlib import Path
import subprocess
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
EVAL = ROOT / 'evaluations/h1_validation'
METRICS = ['pds_cosine', 'expr_mse_unbiased_capped_norm', 'de_wilcoxon_lfc_nmae',
           'de_wilcoxon_direction_fidelity_yield_raw', 'de_wilcoxon_direction_reach_raw',
           'de_wilcoxon_sig_jaccard']


def score(directory):
    directory = Path(directory).resolve()
    anchor = EVAL / 'pilot_anchor1'
    meta = json.loads((anchor / 'anchor_meta.json').read_text())
    if meta['n_splits'] != 1:
        raise ValueError('Pilot utility requires exactly one split, explicitly labeled')
    destination = directory / 'pilot1_diagnostic_scores.csv'
    subprocess.run([str(ROOT / '.venv/bin/cell-eval2'), 'score',
        '--user-agg', str(directory / 'agg_results.csv'),
        '--baseline-agg', str(EVAL / 'baseline_cpu/baseline_agg.csv'),
        '--anchor', str(anchor), '-o', str(destination)], check=True)
    scores = pd.read_csv(destination).set_index('metric')
    values = scores.loc[METRICS, 'from_replicate'].to_numpy(dtype=float)
    if len(values) != 6 or not np.isfinite(values).all():
        raise ValueError('Pilot must supply all six finite diagnostic components')
    result = {'label': 'PROVISIONAL ONE-SPLIT DIAGNOSTIC ONLY', 'n_splits': 1,
              'run_directory': str(directory), 'anchor_directory': str(anchor),
              'provisional_six_component_mean': float(values.mean()),
              'scaled_components': dict(zip(METRICS, map(float, values))),
              'rule': 'Arithmetic mean of exactly six official from_replicate component calculations; not the CLI avg_score, which remains baseline-relative on this diagnostic path.',
              'pending_crosscheck': 'Recompute against real_bundle_cpu when its independent five-split job completes.'}
    destination.with_suffix('.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2), flush=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run_directory', type=Path)
    score(parser.parse_args().run_directory)
