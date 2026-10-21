from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
P=Path(__file__).resolve().parent
d=json.loads((P/'results/resource_assessment.json').read_text())
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.titlesize':12,'savefig.dpi':180})
colors=['#3c668a','#218b83','#c67c28','#a74c73'];masses=[.5,1,2,4]
fig,axes=plt.subplots(2,2,figsize=(11,8),layout='constrained')
for ax,label,title in zip(axes[0],['same_active_volume','30mm_available_depth'],['Same active volume','Port fills a 30 mm depth']):
    for m,color in zip(masses,colors):
        r=[x for x in d['halo_corrected'] if x['mass_meV']==m and x['geometry']==label and not x['uniform_profile']]
        ax.plot([x['pickups'] for x in r],[x['R_halo'] for x in r],'o-',color=color,label=f'{m:g} meV')
    ax.axhline(1,color='#555555',lw=1);ax.set(xscale='log',yscale='log',ylim=(1e-4,1500),xticks=[1,4,16,64],xticklabels=['1','4','16','64'],title=title,xlabel='Pickup columns feeding one ideal coherent receiver',ylabel='Scan rate ratio: radiative / port')
    ax.grid(axis='y',alpha=.15);ax.text(.98,.97,'Above 1: radiative advantage',transform=ax.transAxes,va='top',ha='right',fontsize=9)
axes[0,0].legend(ncol=2,loc='lower left',frameon=False)
ax=axes[1,0]
for lab,ls,col,label in [('same_active_volume','-','#218b83','Same active volume'),('30mm_available_depth','--','#a74c73','Port fills 30 mm')]:
    rr=[x for x in d['loss_thresholds'] if x['geometry']==lab]
    ax.plot([x['mass_meV'] for x in rr],[x['port_loss_for_parity_dB'] for x in rr],'o'+ls,color=col,label=label)
ax.set(xscale='log',xticks=masses,xticklabels=['0.5','1','2','4'],xlabel='Axion mass [meV]',ylabel='Port-chain power loss needed for parity [dB]',title='Measured loss is the useful criterion',ylim=(0,21));ax.legend(frameon=False,fontsize=9);ax.grid(axis='y',alpha=.15)
ax.text(.04,.22,'80% additional optical transmission\nStrongest tested port, including continuous pickup',transform=ax.transAxes,fontsize=9)
ax=axes[1,1]
rr=[x for x in d['loss_thresholds'] if x['geometry']=='same_active_volume']
ax.plot(masses,[x['guide_length_for_parity_cm'] for x in rr],'o-',color='#218b83',label='Nominal copper loss')
ax.plot(masses,[2*x['guide_length_for_parity_cm'] for x in rr],'o--',color='#888888',label='Half the surface resistance')
ax.set(xscale='log',yscale='log',xticks=masses,xticklabels=['0.5','1','2','4'],xlabel='Axion mass [meV]',ylabel='Common output-guide length at parity [cm]',title='Illustration: single-mode copper transport',ylim=(2,400));ax.grid(axis='y',alpha=.15);ax.legend(frameon=False,fontsize=9)
ax.text(.04,.05,'Guide aspect ratio 2:1; frequency / cutoff = 1.7\nNot a lower bound on all possible readout losses',transform=ax.transAxes,fontsize=8.5)
for ax in axes[1]:ax.xaxis.set_minor_formatter(NullFormatter())
fig.suptitle('Where radiative extraction can help',fontsize=17,fontweight='bold')
fig.savefig(P/'figures/comparison_overview.pdf');fig.savefig(P/'figures/comparison_overview.png');plt.close(fig)

# More compact figure for the short note: loss threshold and transport length.
fig,axes=plt.subplots(1,2,figsize=(10.5,3.5),layout='constrained')
for ax,lab,title in zip(axes,['same_active_volume','30mm_available_depth'],['Same active volume','Port fills a 30 mm depth']):
    for m,color in zip(masses,colors):
        r=[x for x in d['halo_corrected'] if x['mass_meV']==m and x['geometry']==lab and not x['uniform_profile']]
        ax.plot([x['pickups'] for x in r],[x['R_halo'] for x in r],'o-',color=color,label=f'{m:g} meV')
    ax.axhline(1,color='#555555',lw=1);ax.set(xscale='log',yscale='log',xticks=[1,4,16,64],xticklabels=['1','4','16','64'],xlabel='Pickup columns',ylabel='Radiative / port scan rate',title=title,ylim=(1e-4,1500));ax.grid(axis='y',alpha=.15)
axes[0].legend(ncol=2,frameon=False,fontsize=9,loc='lower left')
fig.savefig(P/'figures/pickup_comparison.pdf');plt.close(fig)
