from pathlib import Path
import json,sys,hashlib,itertools
from collections import Counter
repo=Path('/home/leechuck/Public/software/gspa');sys.path.insert(0,str(repo/'benchmark/v3'))
import migrate_cohort_feature_contract as M
root=Path.home()/'.cache/gspa-v3-data'
a=root/'cohort-signature-before-v1';b=root/'cohort-signature-availability-v1'
M.verified(a);M.verified(b)
inputs={};sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
for base in [a,b]:
 for name in ['export.json','inputs.json','inference/execution_receipt.json']:
  p=base/name
  if p.exists():inputs[str(p)]=sha(p)
summary={};allchanges=Counter();total_changed=0
for p in sorted((a/'features').glob('*/candidate_features.jsonl')):
 gid=p.parent.name;q=b/'features'/gid/p.name
 inputs[str(p)]=sha(p);inputs[str(q)]=sha(q)
 n=changed=local_lost=effective_lost=0;columns=Counter();representatives=[]
 with p.open() as f,q.open() as g:
  for left,right in itertools.zip_longest(f,g):
   assert left is not None and right is not None
   x,y=json.loads(left),json.loads(right);n+=1
   for key in ['candidate_key','genome_id','protein_id','function_type','function_id','role','sequence_sources','notes']:
    assert x[key]==y[key],(gid,key)
   old,new=x['features'],y['features'];assert old.keys()==new.keys()
   differences=[k for k in old if old[k]!=new[k]]
   assert all(k.startswith('cohort.w.') for k in differences)
   changed+=bool(differences);columns.update(differences)
   local_lost+=old['cohort.w.local_available']==1 and new['cohort.w.local_available']==0
   effective_lost+=old['cohort.w.effective_n']>0 and new['cohort.w.effective_n']==0
   for k in x['available_contexts']:
    if k!='local_reaction':assert x['available_contexts'][k]==y['available_contexts'][k]
 for suffix in ['.integrated.tsv']:
  old=a/'inference'/f'{gid}{suffix}';new=b/'inference'/f'{gid}{suffix}'
  inputs[str(old)]=sha(old);inputs[str(new)]=sha(new);assert sha(old)==sha(new)
 oldreport=json.loads((a/'inference'/f'{gid}.integrated.tsv.cohort.json').read_text())
 newreport=json.loads((b/'inference'/f'{gid}.integrated.tsv.cohort.json').read_text())
 oldamb={r['function_key']:r for r in oldreport['reaction_ambiguity']}
 for r in newreport['reaction_ambiguity']:
  k=r['function_key']
  if k in oldamb and oldamb[k]['representative_reaction']!=r['representative_reaction']:
   representatives.append(dict(function_key=k,before=oldamb[k],after=r))
 summary[gid]=dict(candidates=n,changed_candidates=changed,local_context_withheld=local_lost,
    effective_peer_evidence_withheld=effective_lost,changed_columns=dict(columns),
    changed_representatives_in_bounded_report=representatives,predictions_byte_identical=True)
 allchanges.update(columns);total_changed+=changed
assert inputs=={p:sha(Path(p)) for p in inputs}
result=dict(status='verified',genomes=len(summary),candidates=sum(r['candidates'] for r in summary.values()),
            changed_candidates=total_changed,changed_columns=dict(allchanges),summaries=summary,inputs=inputs,
            representative_limit='Production ambiguity examples are bounded; this is not an exhaustive representative-change count.')
Path(__file__).with_name('real-comparison.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ['summaries','inputs']},indent=2))
