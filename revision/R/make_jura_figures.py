"""Render the manuscript's Jura figures from the canonical saved results."""
from pathlib import Path
import json
import os
import sys
os.environ.setdefault('MPLCONFIGDIR','/tmp/spleakage-matplotlib')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[2]
raw=pd.read_csv(ROOT/'revision/results/application-jura.csv')
selection=pd.read_csv(ROOT/'revision/results/summaries/application-jura-selection-repetitions.csv')
geometry=json.loads((ROOT/'revision/results/summaries/application-jura-geometry.json').read_text())
out=ROOT/'revision/figures'
schemes=['random','block','buffered','knndm']; names=['Random','Blocks','Buffered','kNNDM']
colors=['#0072b2','#d55e00','#cc79a7','#009e73']
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                     'pdf.fonttype':42,'axes.titleweight':'bold'})

def distribution(ax,values,labels):
    obj=ax.boxplot(values,labels=labels,patch_artist=True,showfliers=True,
        flierprops={'markersize':2,'alpha':.6},medianprops={'color':'black'})
    for p,c in zip(obj['boxes'],colors): p.set_facecolor(c);p.set_alpha(.5)
    ax.scatter(np.arange(1,len(values)+1),[np.mean(v) for v in values],marker='D',s=18,c='black',zorder=4)
    ax.tick_params(axis='x',labelsize=9)

if '--supplement' not in sys.argv:
    g=next(g for g in geometry if g['target']=='published')
    f=raw[(raw.target=='published')&(raw.model=='ordinary')]
    fig=plt.figure(figsize=(7.5,6.8),layout='constrained')
    rows=fig.add_gridspec(2,1,height_ratios=[1.25,1])
    top=rows[0].subgridspec(1,2,width_ratios=[1,1.3],wspace=.12)
    distance_axes=top[0,1].subgridspec(3,1,hspace=.08)
    bottom=rows[1].subgridspec(1,3,wspace=.08)
    xy=np.asarray(g['xy']); q=np.asarray(g['target_xy'])
    a=fig.add_subplot(top[0,0])
    a.scatter(xy[:,0],xy[:,1],s=11,color='#4a718c',alpha=.7,label='Training (259)')
    a.scatter(q[:,0],q[:,1],s=18,color='#c17725',marker='x',label='Target (100)')
    a.set(xlabel='Easting (km)',ylabel='Northing (km)',title='A  Training and target sites')
    a.set_aspect('equal',adjustable='datalim');a.legend(fontsize=8.5,loc='upper left')
    xmax=max(max(v) for v in g['cv_distances'].values())*1.02
    assert np.array_equal(g['cv_distances']['block'],g['cv_distances']['buffered'])
    target_x=np.sort(g['target_distances'])
    target_y=np.arange(1,len(target_x)+1)/len(target_x)
    for i,(scheme,label,color) in enumerate(zip(
            ['random','block','knndm'],['Random','Blocks = buffered','kNNDM'],
            [colors[0],colors[1],colors[3]])):
        a=fig.add_subplot(distance_axes[i,0])
        x=np.sort(g['cv_distances'][scheme]);y=np.arange(1,len(x)+1)/len(x)
        a.step(np.r_[0,x,xmax],np.r_[0,y,1],where='post',color=color,linewidth=1.8,
               label=label+' (solid)')
        a.step(np.r_[0,target_x,xmax],np.r_[0,target_y,1],where='post',
               color='black',linestyle='--',linewidth=1.6,label='Target (dashed)')
        a.set(xlim=(0,xmax),ylim=(0,1.06),yticks=[0,.5,1],xticks=[0,.5,1,1.5])
        a.legend(loc='lower right',fontsize=8,handlelength=2.5,framealpha=.95,
                 edgecolor='white',borderpad=.2,labelspacing=.25)
        a.tick_params(labelsize=8.5,labelbottom=i==2)
        if i==0: a.set_title('B  Nearest-training distances',fontsize=10)
        if i==1: a.set_ylabel('Cumulative proportion')
        if i==2: a.set_xlabel('Nearest training distance (km)')
    short_names=['Random','Blocks','Buffered','kNNDM']
    for j,(column,title,ylabel) in enumerate([
            ('W','C  Distance discrepancy',r'$W$ (km)'),
            ('J2','D  Absolute contrast',r'$|J_2|$'),
            ('cv_mse','E  Prediction error','MSE (log scale)')]):
        a=fig.add_subplot(bottom[0,j])
        values=[f.loc[f.scheme==scheme,column] for scheme in schemes]
        if column=='J2': values=[abs(v) for v in values]
        distribution(a,values,short_names)
        a.tick_params(axis='x',labelsize=8.5,rotation=40)
        for tick in a.get_xticklabels(): tick.set_ha('right')
        a.set_ylabel(ylabel)
        a.set_title(title,fontsize=10)
        if column=='cv_mse':
            a.axhline(f.target_mse.iloc[0],color='black',linestyle='--',linewidth=1.4)
            a.text(.98,f.target_mse.iloc[0]+.01,'Target MSE',
                   transform=a.get_yaxis_transform(),ha='right',fontsize=8.5)
        else:
            a.set_ylim(bottom=0)
    fig.savefig(out/'F4-application.pdf',bbox_inches='tight')
    fig.savefig(out/'F4-application.png',dpi=170,bbox_inches='tight')
    plt.close(fig)

if '--main' not in sys.argv:
    rules=['W','J2','fixed_random','fixed_knndm'];labels=['Select W',r'Select $|J_2|$','Always\nrandom','Always\nkNNDM']
    fig,axes=plt.subplots(3,2,figsize=(7.2,9),layout='constrained')
    for i,target in enumerate(['published','uniform','north']):
        for j,model in enumerate(['ordinary','trend']):
            d=selection[(selection.target==target)&(selection.model==model)]
            a=axes[i,j]
            distribution(a,[d.loc[d.rule==rule,'absolute_error'] for rule in rules],labels)
            a.set_ylabel('Absolute CV discrepancy')
            a.set_title({'published':'Original split','uniform':'Random targets','north':'Northern targets'}[target]+
                        ' | '+('Constant mean' if model=='ordinary' else 'Linear trend'),fontsize=10)
            a.set_ylim(bottom=0)
    fig.suptitle('Jura: target and mean-model sensitivities\nAbsolute difference between observed target MSE and CV score',fontsize=12)
    fig.savefig(out/'S5-jura-targets.pdf',bbox_inches='tight')
    fig.savefig(out/'S5-jura-targets.png',dpi=150,bbox_inches='tight')
    plt.close(fig)
print('Generated requested Jura figures from saved results.')
