#!/usr/bin/env python3
"""Export the completed control-only sampling comparison, without changing gates."""
from pathlib import Path
import os,json
ROOT=Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'.cache/matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def main():
    depth=json.loads((ROOT/'reports/sampling_precision_results.json').read_text())
    cube=json.loads((ROOT/'reports/balanced_sampling_results.json').read_text())
    contexts=['H1','A','B','C'];fig,axes=plt.subplots(1,2,figsize=(11,4.7),layout='constrained')
    for ax,key,title in zip(axes,['arithmetic_log_mse','bulk_log_mse'],['Arithmetic profile error','Bulk profile error']):
        for offset,label,color,results in [(-.10,'Depth only','#b57623',depth),(.10,'Depth + four expression PCs','#187b9a',cube)]:
            values=[];lower=[];upper=[]
            for context in contexts:
                r=results['contexts'][context]['measures'][key]
                if results is cube:r=r['uniform']
                values.append(r['ratio']);lower.append(r['ratio_95_ci'][0]);upper.append(r['ratio_95_ci'][1])
            values=np.asarray(values)
            ax.errorbar(np.arange(4)+offset,values,yerr=[values-lower,np.asarray(upper)-values],
                        fmt='o',capsize=4,color=color,label=label,markersize=7,linewidth=1.7)
        ax.axhline(1,color='#697386',linewidth=1,linestyle='--')
        ax.set_xticks(range(4),contexts);ax.set_title(title,loc='left',weight='bold')
        ax.set_ylabel('MSE / uniform-sampling MSE  ↓');ax.set_ylim(.91,1.025)
        ax.grid(axis='y',alpha=.18);ax.spines[['top','right']].set_visible(False)
    axes[1].plot([-.34,.34],[.95,.95],color='#aa3232',linewidth=2)
    axes[1].annotate('H1 requirement: ≤0.95',xy=(.04,.95),xytext=(.7,.925),fontsize=9,
                     arrowprops={'arrowstyle':'-','color':'#aa3232'},color='#aa3232')
    handles,labels=axes[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='outside upper center',ncol=2,frameon=False)
    fig.supxlabel('400 actual cells per draw; 200 fixed draws per context.\n95% paired-seed intervals describe Monte Carlo uncertainty, not biological generalization.',fontsize=10)
    for suffix in ['png','svg']:fig.savefig(ROOT/f'reports/balanced_sampling_20261021.{suffix}',dpi=180,bbox_inches='tight')
    plt.close(fig)


if __name__=='__main__':main()
