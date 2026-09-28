#!/usr/bin/env python3
"""Deterministic descriptive ranking reanalysis; no inference or network access.

Primary: the previously designated untruncated set. Both species' shared proteins
are retained together through connected-component resampling. Results are
conditional on the selected union, not a population performance estimate.
"""
from pathlib import Path
import hashlib,json,sys
import numpy as np
import pandas as pd
from Bio import SeqIO
from scipy.stats import spearmanr
R=Path(__file__).resolve().parents[1]; O=R/'reports/paper_revision'; O.mkdir(exist_ok=True)
used={}
def read(p):
 f=R/p;used[p]=hashlib.sha256(f.read_bytes()).hexdigest();return pd.read_csv(f)
def rho(x,y):return float(spearmanr(x,y).statistic)
def run():
 keys=['zebrafish_id','aescallii_id']; m=read('reports/esm3_subset_pairs.csv')
 for tag in ['8m','35m','150m','650m','esm3']:
  file='reports/esm3_subset_'+('esm3_1p4b' if tag=='esm3' else 'esm2_'+tag)+'.csv'
  d=read(file);col='esm3_cosine_distance' if tag=='esm3' else 'esm2_cosine_distance'
  m=m.merge(d[keys+[col]].rename(columns={col:tag}),on=keys,validate='one_to_one')
 identity=read('reports/subset_identity.csv').rename(columns={'z':'zebrafish_id','a':'aescallii_id'})
 m=m.merge(identity[keys+['identity']],on=keys,validate='one_to_one');m['seqdiv']=100-m.identity
 seq={}
 for p in sorted((R/'data/reference/sequences/subset').glob('*.fasta')):
  used[str(p.relative_to(R))]=hashlib.sha256(p.read_bytes()).hexdigest()
  for r in SeqIO.parse(p,'fasta'):
   key=r.id.split('.')[0];s=str(r.seq);assert key not in seq or seq[key]==s;seq[key]=s
 m['max_length']=[max(len(seq[z]),len(seq[a])) for z,a in m[keys].itertuples(index=False,name=None)]
 m['nonstandard']=[any(set(seq[x])-set('ACDEFGHIKLMNPQRSTVWY') for x in (z,a)) for z,a in m[keys].itertuples(index=False,name=None)]
 assert len(m)==3428 and not m.duplicated(keys).any() and np.isfinite(m[['8m','35m','150m','650m','esm3']]).all().all()
 # Preserve the input row order for historical nlargest tie handling; report sensitivity separately.
 m['input_order']=np.arange(len(m));m.to_csv(O/'pair_table.csv',index=False)
 sets={'primary_untruncated':m[m.max_length<=2046].copy(),'full':m.copy(),'clean':m[~m.nonstandard].copy(),'clean_untruncated':m[(m.max_length<=2046)&~m.nonstandard].copy()}
 def components(d):
  parent={}
  def find(x):
   parent.setdefault(x,x)
   if parent[x]!=x:parent[x]=find(parent[x])
   return parent[x]
  for z,a in d[keys].itertuples(index=False,name=None):parent[find('z:'+z)]=find('a:'+a)
  labels=[find('z:'+z) for z in d.zebrafish_id]
  return list(pd.Series(np.arange(len(d))).groupby(labels,sort=True).apply(list))
 def overlaps(d,k,tie='input'):
  x=d if tie=='input' else d.sort_values(keys,ascending=(tie=='lexical'))
  target=set(x.nlargest(k,'esm3').index)
  return {c:100*len(target&set(x.nlargest(k,c).index))/k for c in ['8m','35m','150m','650m','seqdiv']}
 result={'definition':{'primary':'max recorded sequence length <= 2046; selected before this revision','bootstrap':'connected components of the bipartite protein graph within each analysis set; components resampled with replacement retaining all rows','replicates':2000,'seed':0,'ci':'percentile 2.5% and 97.5%; selected-sample uncertainty','topk':'k=round(N*f), pandas nlargest keep=first in saved input order; boundary sensitivity also reported','inferential_superiority_test':False},'sets':{}}
 for name,d in sets.items():
  gr=components(d);x=d['8m'].to_numpy();y=d.esm3.to_numpy();rng=np.random.default_rng(0)
  boot=[]
  for _ in range(2000):
   ix=np.concatenate([gr[i] for i in rng.integers(len(gr),size=len(gr))]);boot.append(rho(x[ix],y[ix]))
  stats={'n':len(d),'rho':{c:rho(d[c],d.esm3) for c in ['8m','35m','150m','650m']},'components':len(gr),'max_component_pairs':max(map(len,gr)),'ci_8m_component':np.quantile(boot,[.025,.975]).tolist(),'topk':[]}
  for f in [.01,.05,.1]:
   k=round(len(d)*f);row={'fraction':f,'k':k,'overlap':overlaps(d,k),'sensitivity':[]}
   for k2 in sorted({int(len(d)*f),k,int(np.ceil(len(d)*f))}):
    for tie in ['input','lexical','reverse_lexical']:row['sensitivity'].append({'k':k2,'tie':tie,'overlap':overlaps(d,k2,tie)})
   stats['topk'].append(row)
  result['sets'][name]=stats
 result['selection_subsets']={k:{'n':int(m[k].sum()),'rho_8m':rho(m.loc[m[k],'8m'],m.loc[m[k],'esm3'])} for k in ['in_random','in_stratified','in_top','in_target']}
 result['inputs_sha256']=used;result['script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest();result['environment']={'python':sys.version.split()[0],'numpy':np.__version__,'pandas':pd.__version__}
 (O/'ranking.json').write_text(json.dumps(result,indent=2))
 for name,s in result['sets'].items():print(name,s['n'],s['rho']['8m'],s['ci_8m_component'])
if __name__=='__main__':run()
