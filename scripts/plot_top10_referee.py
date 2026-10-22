#!/usr/bin/env python3
"""Plot observed score gaps; this is arithmetic attribution, not model causality."""
import os
os.environ.setdefault('MPLCONFIGDIR', str(__import__('pathlib').Path(__file__).resolve().parents[1] / '.cache/matplotlib'))
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
record = json.loads((ROOT / 'reports/top10_referee_20261021.json').read_text())
keys = ['score_mse', 'score_reach', 'score_pds', 'score_nmae', 'score_fid', 'score_jac']
labels = ['Expression MSE', 'Direction reach', 'PDS', 'LFC NMAE', 'Direction fidelity', 'DE overlap']
incumbent = np.array([record['model4']['scaled'][k] for k in keys])
top = np.array([[row['scaled'][k] for k in keys] for row in record['top10_unique_teams']])
average = top.mean(axis=0)
contribution = (average - incumbent) / 6
share = contribution / contribution.sum() * 100
y = np.arange(len(keys))
plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
fig, axes = plt.subplots(1, 2, figsize=(12, 5.2), gridspec_kw={'width_ratios': [1.45, 1]})
ax = axes[0]
ax.barh(y - .16, incumbent, height=.3, color='#D87530', label='Model4')
ax.barh(y + .16, average, height=.3, color='#267CA6', label='Top-10 mean')
for i in range(len(keys)):
    ax.plot([top[:, i].min(), top[:, i].max()], [i + .16, i + .16], color='#193C50', lw=1.4)
ax.set_yticks(y, labels)
ax.invert_yaxis()
ax.axvline(0, color='#555555', linewidth=.7)
ax.set_xlim(-.06, .95)
ax.set_xlabel('Reference-scaled component (higher is better)')
ax.set_title('Where current top teams outperform Model4', loc='left', fontweight='bold')
ax.legend(frameon=False, loc='lower right', title='Lines: top-10 min–max', title_fontsize=8)
ax = axes[1]
ax.barh(y, share, color=['#267CA6', '#267CA6', '#8FA5AF', '#8FA5AF', '#8FA5AF', '#8FA5AF'])
ax.set_yticks(y, labels)
ax.invert_yaxis()
ax.set_xlim(0, 50)
ax.set_xlabel('Share of overall score gap (%)')
ax.set_title('MSE + reach account for 62.3% of the gap', loc='left', fontweight='bold')
for i, value in enumerate(share):
    ax.text(value + .7, i, f'{value:.1f}%', va='center')
fig.text(.01, .015,
         'Official snapshot: 21 Oct 2026, 10:18 HKT • Same panel and anchors • Score decomposition does not identify team methods',
         fontsize=8, color='#555555')
fig.tight_layout(rect=(0, .05, 1, 1))
for suffix in ['png', 'svg']:
    fig.savefig(ROOT / f'reports/top10_gap_20261021.{suffix}', dpi=180, bbox_inches='tight')
print('Saved reports/top10_gap_20261021.png and .svg')
