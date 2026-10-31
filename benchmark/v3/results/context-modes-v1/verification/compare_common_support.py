from pathlib import Path
import json,sys,hashlib
import numpy as np
repo=Path('/home/leechuck/Public/software/gspa');sys.path.insert(0,str(repo/'scripts'))
import train_context_model as T
root=Path.home()/'.cache/gspa-v3-data/cohort-context-modes-v2'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
out={};bindings={}
for kind in ('logistic','gbdt'):
 arms={};supported={}
 for mode in ('UNIVERSAL','PATHWAY','ORGANISM'):
  table=root/mode/'table/training_table.jsonl';model=root/mode/'models'/f'model_{kind}.json'
  contract=T.load_contract(next((root/mode/'features').glob('*/feature_contract.json')))
  rows={r['example_id']:r for r in T.read_table(table,contract) if r['split']=='test'}
  artifact=json.loads(model.read_text());bindings[str(table)]=sha(table);bindings[str(model)]=sha(model)
  supported[mode]={k for k,r in rows.items() if T.unsupported_input_reason(artifact,r['x'],r['sequence_sources'],cohort_context_mode=r.get('cohort_context_mode','UNIVERSAL')) is None}
  arms[mode]=(rows,artifact,model,contract,table)
 common=set.intersection(*supported.values());metrics={};preds={};labels={}
 for mode,(rows,artifact,model,contract,table) in arms.items():
  selected=[rows[k] for k in sorted(common)]
  metrics[mode]=T.evaluate_frozen(model,contract,selected,sha(model),sha(table))
  preds[mode]={k:T.predict_proba(artifact,rows[k]['x'],rows[k]['sequence_sources'],cohort_context_mode=rows[k].get('cohort_context_mode','UNIVERSAL')) for k in sorted(common)}
  for k in common:
   identity=(rows[k]['genome_id'],rows[k]['label'])
   if k in labels:assert labels[k]==identity
   labels[k]=identity
 genomes=sorted({v[0] for v in labels.values()});rng=np.random.default_rng(20261027)
 draws=rng.integers(0,len(genomes),size=(10000,len(genomes)));deltas={}
 for mode in ('PATHWAY','ORGANISM'):
  per_genome={}
  for g in genomes:
   ks=[k for k in common if labels[k][0]==g]
   per_genome[g]=float(np.mean([(preds[mode][k]-labels[k][1])**2-(preds['UNIVERSAL'][k]-labels[k][1])**2 for k in ks]))
  values=np.asarray([per_genome[g] for g in genomes]);boot=values[draws].mean(axis=1)
  deltas[mode]=dict(per_genome_brier_difference=per_genome,mean_genome_difference=float(values.mean()),
                    paired_genome_bootstrap_95=list(np.quantile(boot,[.025,.975])))
 out[kind]=dict(common_test_rows=len(common),individual_supported_rows={m:len(v) for m,v in supported.items()},
                test_genomes=genomes,metrics=metrics,brier_difference_vs_universal=deltas)
assert all(sha(Path(p))==h for p,h in bindings.items())
result=dict(status='verified',inputs=bindings,results=out,
 scope='Descriptive comparison on the intersection of train/calibration-supported test profiles. Four previously examined test genomes with known reference overlap; not independent accuracy or a biological significance claim.',
 bootstrap=dict(unit='genome',paired=True,replicates=10000,seed=20261027,metric='mean per-genome Brier difference; negative favors comparator'))
Path(__file__).with_name('common-support-comparison.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
print(json.dumps({k:dict(common=v['common_test_rows'],deltas=v['brier_difference_vs_universal']) for k,v in out.items()}))
