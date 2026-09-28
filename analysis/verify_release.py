#!/usr/bin/env python3
"""Verify frozen saved outputs. No network or model inference. Nonzero on mismatch."""
from pathlib import Path
import argparse, hashlib, json, sys
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
R=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def verify():
 checks=[]
 def need(ok,name):
  if not ok:raise AssertionError(name)
  checks.append(name)
 manifest=R/'FILE_MANIFEST.json'
 if manifest.exists():
  for rel,h in json.loads(manifest.read_text())['sha256'].items():need((R/rel).is_file() and sha(R/rel)==h,'file '+rel)
 else:raise AssertionError('FILE_MANIFEST.json missing')
 for file in ['ranking.json','figure_manifest.json']:
  x=json.loads((R/'reports/paper_revision'/file).read_text())
  for rel,h in x['inputs_sha256'].items():need(sha(R/rel)==h,'source '+rel)
  if file=='ranking.json':need(sha(R/'analysis/revision_statistics.py')==x['script_sha256'],'ranking code')
  else:
   need(sha(R/x['generator'])==x['generator_sha256'],'figure code')
   for rel,h in x['outputs_sha256'].items():need(sha(R/'reports/figures/final'/rel)==h,'figure '+rel)
 m=pd.read_csv(R/'reports/paper_revision/pair_table.csv'); keys=['zebrafish_id','aescallii_id']
 need(len(m)==3428 and not m.duplicated(keys).any(),'3428 unique pair keys')
 need(np.isfinite(m[['8m','35m','150m','650m','esm3']].to_numpy()).all(),'all pair distances finite')
 sets={'primary_untruncated':m[m.max_length<=2046],'full':m,'clean':m[~m.nonstandard],'clean_untruncated':m[(m.max_length<=2046)&~m.nonstandard]}
 ranking=json.loads((R/'reports/paper_revision/ranking.json').read_text())
 for name,d in sets.items():
  s=ranking['sets'][name];need(len(d)==s['n'],name+' count')
  for col in ['8m','35m','150m','650m']:need(abs(float(spearmanr(d[col],d.esm3).statistic)-s['rho'][col])<1e-12,name+' rho '+col)
  for row in s['topk']:
   k=row['k']; target=set(sorted(d.index,key=lambda i:-d.loc[i,'esm3'])[:k])
   for col in ['8m','35m','150m','650m','seqdiv']:
    selected=set(sorted(d.index,key=lambda i:-d.loc[i,col])[:k]);need(abs(100*len(target&selected)/k-row['overlap'][col])<1e-12,name+f' top{k} '+col)
 vector_errors={}
 for col,archive in [('8m','esm2_8m'),('35m','esm2_35m'),('150m','esm2_150m'),('650m','esm2_650m'),('esm3','esm3_1p4b')]:
  with np.load(R/'data/embeddings'/f'{archive}.npz',allow_pickle=False) as z:
   vectors={key:z[key].astype(np.float64) for key in z.files}; errors=[]
   for a,b,expected in m[keys+[col]].itertuples(index=False,name=None):
    x,y=vectors[a],vectors[b]
    need(np.isfinite(x).all() and np.isfinite(y).all(),'finite vectors '+col+' '+a+' '+b)
    denom=float(np.sqrt(np.sum(x*x))*np.sqrt(np.sum(y*y)))
    need(denom>0 and np.isfinite(denom),'nonzero vector norm '+col+' '+a+' '+b)
    calc=1-float(np.sum(x*y))/denom;errors.append(abs(calc-expected))
   vector_errors[col]={'pairs':len(errors),'max_abs_error':max(errors),'n_vectors':len(vectors)}
   need(len(errors)==3428 and np.isfinite(errors).all() and max(errors)<2e-6,'saved vector distance '+col)
 for name in ['ahr2','kcnj13']:
  fixed=pd.read_csv(R/f'reports/{name}_residue_divergence_fixed.csv').set_index('residue');d=pd.read_csv(R/f'reports/figures/{name}_plm_vs_structure.csv')
  need(d.ri.is_unique and np.allclose(d.divergence,fixed.loc[d.ri,'divergence'],rtol=0,atol=1e-12),name+' structure aligned to corrected profile')
 need(len(list((R/'reports/figures/final').glob('*.png')))==10,'ten final figures')
 need(len(list((R/'tables').glob('Table*.csv')))==6,'six manuscript tables')
 result={'status':'PASS','checks':len(checks),'scope':'saved-output numerical consistency, not biological validation or historical inference provenance','vector_distance_checks':vector_errors,'environment':{'python':sys.version.split()[0],'numpy':np.__version__,'pandas':pd.__version__}}
 print(json.dumps(result,indent=2));return result
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path);a=ap.parse_args();r=verify()
 if a.out:a.out.write_text(json.dumps(r,indent=2))
