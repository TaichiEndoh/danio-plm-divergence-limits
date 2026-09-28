#!/usr/bin/env python3
"""Build revised scientific figures from saved outputs; no model inference."""
from pathlib import Path
import json,hashlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from tier_stratification import read_plddt, AHR2_POCKET
from ahr_pocket_axis import COMPARISONS, POCKET
R=Path(__file__).resolve().parents[1]; O=R/'reports/figures/final'; O.mkdir(exist_ok=True)
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'white','savefig.facecolor':'white','font.family':'DejaVu Sans'})
C=['#0072B2','#009E73','#E69F00','#D55E00','#777777']; used={};out={};metrics={}
def track(p):
 f=R/p;used[p]=hashlib.sha256(f.read_bytes()).hexdigest();return f
def csv(p):return pd.read_csv(track(p))
def save(fig,name):
 fig.savefig(O/name,dpi=300,bbox_inches='tight');plt.close(fig);out[name]=hashlib.sha256((O/name).read_bytes()).hexdigest()
def label(ax,s):ax.set_title(s,loc='left',fontweight='bold',fontsize=11)
m=csv('reports/paper_revision/pair_table.csv'); primary=m[m.max_length<=2046]
r=json.loads(track('reports/paper_revision/ranking.json').read_text())
track('analysis/tier_stratification.py');track('analysis/ahr_pocket_axis.py')
# Study design
fig,ax=plt.subplots(figsize=(9,4.6));ax.set(xlim=(0,10),ylim=(0,6));ax.axis('off')
boxes=[(2,5,'Recorded ortholog pairs\n50,921 distinct pairs'),(7.5,5,'Selected overlapping subsets\n3,428 pairs; 5,614 proteins'),(2,3,'Whole-protein ranks\n4 ESM-2 checkpoints vs ESM-3\nPrimary: 2,662 untruncated pairs'),(7.5,3,'Two residue case studies\nAHR2 and Kcnj13\nSequences, features, predicted structures'),(4.8,.9,'Interpretation checks\nSequence uncertainty · coordinates · shared proteins\nModel agreement is not functional validation')]
for x,y,t in boxes:ax.text(x,y,t,ha='center',va='center',bbox=dict(boxstyle='round,pad=0.7',fc='#edf4f8',ec='#7292a3'),fontsize=10)
for xy,xytext in [((5.5,5),(4,5)),((2,3.7),(2,4.4)),((7.5,3.7),(7.5,4.4)),((3.8,1.6),(2.3,2.3)),((6,1.6),(7.2,2.3))]:ax.annotate('',xy=xy,xytext=xytext,arrowprops=dict(arrowstyle='->',color='#555555'))
save(fig,'Figure1_study_design.png')
# Correlation
fig,axs=plt.subplots(1,2,figsize=(10,4.4),layout='constrained')
axs[0].scatter(primary['8m'].rank(pct=True),primary.esm3.rank(pct=True),s=5,alpha=.3,c=C[0],rasterized=True);axs[0].set(xlabel='ESM-2 8M rank percentile',ylabel='ESM-3 rank percentile',xlim=(0,1.02),ylim=(0,1.02));label(axs[0],'a  Primary sample (n = 2,662)');axs[0].text(.03,.95,f"ρ = {r['sets']['primary_untruncated']['rho']['8m']:.3f}\n95% CI {r['sets']['primary_untruncated']['ci_8m_component'][0]:.3f}–{r['sets']['primary_untruncated']['ci_8m_component'][1]:.3f}",transform=axs[0].transAxes,va='top',bbox=dict(fc='white',ec='none',alpha=.85))
for name,col,lab in [('primary_untruncated',C[0],'Primary (2,662)'),('full',C[3],'Full union (3,428)')]:axs[1].plot(range(4),list(r['sets'][name]['rho'].values()),'o-',color=col,label=lab)
axs[1].set(xticks=range(4),xticklabels=['8M','35M','150M','650M'],ylabel='Spearman ρ with ESM-3',xlabel='ESM-2 checkpoint',ylim=(.8,1));axs[1].legend(frameon=False,loc='lower left');label(axs[1],'b  Observed model concordance');save(fig,'Figure2_scale_robustness.png')
# Top-k
fig,axs=plt.subplots(1,2,figsize=(10,4.5),layout='constrained')
for ax,name,title in zip(axs,['primary_untruncated','full'],['a  Primary sample','b  Full selected union']):
 s=r['sets'][name];x=np.arange(3)
 for j,(c,lab) in enumerate(zip(['8m','35m','150m','650m','seqdiv'],['8M','35M','150M','650M','Sequence identity'])):ax.bar(x+(j-2)*.15,[t['overlap'][c] for t in s['topk']],width=.145,color=C[j],label=lab)
 for i,t in enumerate(s['topk']):ax.plot([i-.4,i+.4],[100*t['k']/s['n']]*2,'k--',lw=1)
 ax.set(xticks=x,xticklabels=[f"{t['fraction']:.0%}\nk = {t['k']}" for t in s['topk']],ylabel='Overlap with ESM-3 top-k (%)',ylim=(0,100));label(ax,title)
axs[1].legend(ncol=2,frameon=False,fontsize=8,loc='upper left');save(fig,'Figure3_prioritization.png')
# Profiles
D={p:csv(f'reports/{p}_residue_divergence_fixed.csv') for p in ['ahr2','kcnj13']}; E={p:csv(f'reports/{p}_residue_divergence_esm3.csv') for p in D};X=[96,106,581,620,872,981]
fig,axs=plt.subplots(2,2,figsize=(10,6),layout='constrained')
for j,p in enumerate(D):
 for i,model in enumerate(['ESM-2 8M (complete windows)','ESM-3 (full context)']):
  ax=axs[i,j];d=D[p] if i==0 else E[p];col='divergence' if i==0 else 'divergence_esm3';ax.plot(d.residue,d[col],lw=.9,color=C[i]);ax.set(xlabel='Reference residue',ylabel='Local distance');label(ax,f"{'abcd'[i*2+j]}  {p.upper() if p=='ahr2' else 'Kcnj13'} · {model}")
  if p=='ahr2':
   for xx in X:ax.axvline(xx,color=C[3],alpha=.4,lw=.7)
  else:
   for lo,hi in [(11,29),(166,186)]:ax.axvspan(lo,hi,color=C[2],alpha=.2)
save(fig,'Figure4_hotspots_and_axes.png')
# Confidence tiers from original PDB, not colored derivatives
fig,axs=plt.subplots(2,2,figsize=(10,7),layout='constrained')
for j,p in enumerate(D):
 pl=read_plddt(track(f'reports/{p}_esmfold.pdb'));d=D[p].copy();d['plddt']=d.residue.map(pl);assert d.plddt.notna().all()
 rho=float(spearmanr(d.divergence,d.plddt).statistic);ax=axs[0,j];ax.scatter(d.plddt,d.divergence,s=7,alpha=.5,color=C[j]);ax.set(xlabel='Reference pLDDT',ylabel='Local ESM-2 distance');label(ax,f"{'ab'[j]}  {p.upper() if p=='ahr2' else 'Kcnj13'} (ρ = {rho:.3f})")
 pocket=set(AHR2_POCKET) if p=='ahr2' else set();d['tier']=['Pocket' if r.residue in pocket else '≥70' if r.plddt>=70 else '50–<70' if r.plddt>=50 else '<50' for r in d.itertuples()];order=[t for t in ['Pocket','≥70','50–<70','<50'] if t in set(d.tier)];means=[d.loc[d.tier==t,'divergence'].mean() for t in order];ns=[sum(d.tier==t) for t in order];ax=axs[1,j];ax.bar(range(len(order)),means,color=C[:len(order)]);ax.set(xticks=range(len(order)),xticklabels=[f'{t}\nn = {n}' for t,n in zip(order,ns)],ylabel='Mean local ESM-2 distance',xlabel='Annotation / reference pLDDT tier',ylim=(0,max(means)*1.25));label(ax,f"{'cd'[j]}  Tier summaries")
 for i,v in enumerate(means):ax.text(i,v,f'{v:.2g}',ha='center',va='bottom',fontsize=9)
 metrics[p+'_tiers']={t:{'n':int(n),'mean':float(v)} for t,n,v in zip(order,ns,means)}
save(fig,'Figure5_structural_grounding.png')
# Features: labels include mapped coordinates and distinct helices; small values stay numeric
f=csv('reports/feature_divergence.csv');fig,axs=plt.subplots(1,2,figsize=(11,4.7),layout='constrained')
labels={'AHR2':['bHLH 25–86','PAS-A 119–231','PAS fold 308–385','C-terminal 386–1027','Unannotated'], 'Kcnj13':['TM domain 42–171','TM1 62–85','TM2 142–166','C-terminal 178–315','Unannotated']}
for j,p in enumerate(labels):
 d=f[f.protein==p];ax=axs[j];y=np.arange(len(d));ax.barh(y,d['mean'],color=C[j]);ax.set(yticks=y,yticklabels=labels[p],xlabel='Mean local ESM-2 distance',xlim=(0,d['mean'].max()*1.35));ax.invert_yaxis();label(ax,f"{'ab'[j]}  {p}")
 for yy,v in zip(y,d['mean']):ax.text(v+d['mean'].max()*.025,yy,f'{v:.2g}',va='center',fontsize=9)
metrics['features']=f.fillna('').to_dict('records');save(fig,'Figure6_feature_annotation.png')
# Structure corrected table
fig,axs=plt.subplots(1,2,figsize=(10,4.4),layout='constrained')
for j,p in enumerate(D):
 d=csv(f'reports/figures/{p}_plm_vs_structure.csv');core=(d.plddt_a>70)&(d.plddt_b>70);ax=axs[j];ax.scatter(d.loc[~core,'ca_dev'],d.loc[~core,'divergence'],s=8,alpha=.4,c='#999999',label='Outside confident core');ax.scatter(d.loc[core,'ca_dev'],d.loc[core,'divergence'],s=9,alpha=.6,c=C[j],label='Both pLDDT >70');rho=float(spearmanr(d.ca_dev,d.divergence).statistic);rc=float(spearmanr(d.loc[core,'ca_dev'],d.loc[core,'divergence']).statistic);ax.set(xlabel='Cα deviation after core fit (Å)',ylabel='Local ESM-2 distance');label(ax,f"{'ab'[j]}  {p.upper() if p=='ahr2' else 'Kcnj13'} (n = {len(d)})");ax.text(.98,.95,f'Overall ρ = {rho:.3f}\nCore ρ = {rc:.3f}\nCore n = {sum(core)}',transform=ax.transAxes,ha='right',va='top');ax.legend(frameon=False,fontsize=8,loc='upper left');metrics[p+'_structure']={'n':len(d),'core_n':int(sum(core)),'rho':rho,'core_rho':rc}
save(fig,'FigureS1_rmsd_control.png')
# Archived summaries only, absence displayed as missing
v=csv('reports/variant_rediscovery_benchmark.csv');v['label']=v.protein.str.replace('_HUMAN','',regex=False)+' '+v.variant;models=['esm2_t6_8M_UR50D','esm2_t12_35M_UR50D','esm2_t30_150M_UR50D','esm2_t33_650M_UR50D'];vl=list(dict.fromkeys(v.label));fig,axs=plt.subplots(1,2,figsize=(10,4.5),layout='constrained');cmap=plt.cm.viridis_r.copy();cmap.set_bad('#dddddd')
for ax,score,t in zip(axs,['wt-marginal','masked-marginal'],['a  Wild-type marginal','b  Masked marginal']):
 arr=v[v.scoring==score].pivot(index='label',columns='model',values='percentile_most_impactful').reindex(index=vl,columns=models).to_numpy();im=ax.imshow(arr,cmap=cmap,vmin=0,vmax=100,aspect='auto');ax.set(xticks=range(4),xticklabels=['8M','35M','150M','650M'],yticks=range(len(vl)),yticklabels=vl);label(ax,t)
 for row in range(arr.shape[0]):
  for col in range(4):
   x=arr[row,col];ax.text(col,row,'—' if np.isnan(x) else f'{x:.1f}',ha='center',va='center',fontsize=9,color='white' if x>55 else 'black')
fig.colorbar(im,ax=axs,shrink=.8,label='Archived percentile (lower = more negative relative score)');save(fig,'FigureS2_variant_rediscovery.png')
# Baseline distribution
b=csv('data/reference/esm2_8m_pair_distances_dedup.csv');col='plm_cosine_distance';assert len(b)==50921
fig,ax=plt.subplots(figsize=(8,4));ax.hist(b[col],bins=100,color=C[0]);ax.set(xlabel='Archived ESM-2 8M cosine distance',ylabel='Pair count',yscale='log');label(ax,'Recorded baseline distribution (50,921 unique pairs)');save(fig,'FigureS3_baseline_distribution.png')
# X sensitivity profile, pocket ratio inset deliberately omitted to avoid unneeded panels
fig,ax=plt.subplots(figsize=(10,3.8));d=D['ahr2'];ax.plot(d.residue,d.divergence,color=C[0],lw=1);domains=[(25,86,'bHLH'),(119,231,'PAS-A'),(308,385,'PAS fold')]
for a,b,n in domains:ax.axvspan(a,b,color=C[1],alpha=.12);ax.text((a+b)/2,d.divergence.max()*1.04,n,ha='center',fontsize=8)
for i,x in enumerate(X):ax.axvspan(x-10,x+10,color=C[3],alpha=.18,label='Centers within 10 residues of X' if i==0 else None)
ax.set(xlabel='D. rerio reference residue',ylabel='Local ESM-2 distance',ylim=(min(d.divergence),max(d.divergence)*1.15));ax.legend(frameon=False,fontsize=9,loc='upper right');save(fig,'FigureS4_ahr2_residue_profile.png')
manifest={'generator':'analysis/revision_figures.py','generator_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'inputs_sha256':used,'outputs_sha256':out,'metrics':metrics,'scope':'Saved-output reconstruction; no new model inference or biological validation.'};(R/'reports/paper_revision/figure_manifest.json').write_text(json.dumps(manifest,indent=2));print('Generated',len(out),'figures')
