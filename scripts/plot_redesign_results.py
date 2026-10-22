#!/usr/bin/env python3
"""Render the preregistered two-seed six-metric comparison as exportable figures."""
from pathlib import Path
import os
ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR', str(ROOT/'.cache/matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main():
    path = ROOT/'reports/redesign_experiment_metrics.csv'
    data = pd.read_csv(path)
    if not {'arm','response','emitter','seed','metric','value'}.issubset(data.columns):
        raise ValueError('Expected complete tidy metric table')
    metrics = [
        ('pds_cosine', 'PDS ↑'),
        ('expr_mse_unbiased_capped_norm', 'Expression MSE ↓'),
        ('de_wilcoxon_lfc_nmae', 'LFC NMAE ↓'),
        ('de_wilcoxon_direction_fidelity_yield_raw', 'Direction fidelity ↑'),
        ('de_wilcoxon_direction_reach_raw', 'Direction reach ↑'),
        ('de_wilcoxon_sig_jaccard', 'Significance overlap ↑')]
    recipes = [('current','current'), ('current','preserving'), ('adapter','current'), ('adapter','preserving')]
    labels = ['Current response\nCurrent emitter', 'Current response\nNew emitter',
              'Adapter\nCurrent emitter', 'Adapter\nNew emitter']
    colors = ['#657284', '#167d9a', '#a66b19', '#6551b8']
    seeds = [20261021,20261022]
    fig, axes = plt.subplots(2,3,figsize=(14,8), layout='constrained')
    for ax,(metric,title) in zip(axes.ravel(),metrics):
        for j,((response,emitter),color) in enumerate(zip(recipes,colors)):
            selected = data[(data.metric==metric)&(data.response==response)&(data.emitter==emitter)]
            if len(selected)!=2 or set(selected.seed)!=set(seeds):
                raise ValueError('Two complete seeds required for every recipe/metric')
            values = selected.set_index('seed').value.reindex(seeds).to_numpy()
            ax.plot([j-.07,j+.07], values, color=color, linewidth=1, alpha=.6)
            ax.scatter(j-.07,values[0],color=color,marker='o',s=48,zorder=3)
            ax.scatter(j+.07,values[1],color=color,marker='D',s=42,zorder=3)
            ax.plot([j-.22,j+.22], [values.mean()]*2, color=color,linewidth=2)
        ax.set_title(title,loc='left',weight='bold',fontsize=12)
        ax.set_xticks(range(4),labels,fontsize=8)
        ax.grid(axis='y',alpha=.2);ax.set_axisbelow(True)
        ax.spines[['top','right']].set_visible(False)
        ax.tick_params(axis='y',labelsize=9)
        ax.margins(y=.18)
    fig.suptitle('H1 development: response model × cell emitter',fontsize=18,weight='bold')
    fig.supxlabel('50 targets × 400 cells; raw local metrics. Circle/diamond = two fixed seeds; short line = their mean.\n'
                  'Arrows mark the better direction. Panels use separate scales; no composite score is implied.',fontsize=10)
    for extension in ['png','svg']:
        fig.savefig(ROOT/f'reports/redesign_metrics_20261021.{extension}',dpi=180,bbox_inches='tight')
    plt.close(fig)


if __name__=='__main__':main()
