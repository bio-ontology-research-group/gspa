from pathlib import Path
from collections import Counter
import json,hashlib
cache=Path.home()/'.cache/gspa-v3-data';before=cache/'pathway-gap-before-v1';after=cache/'pathway-gap-after-v1'
paths=[p/'integrated/features/candidate_features.jsonl' for p in [before,after]]
rows=[{r['candidate_key']:r for r in map(json.loads,p.read_text().splitlines())} for p in paths]
a,b=rows
assert a.keys()==b.keys()
changes=Counter()
for k in a:
 assert a[k]['contract']==b[k]['contract']
 for f,v in a[k]['features'].items():
  if v!=b[k]['features'][f]:
   assert f in ['ctx.gap.available','ctx.gap.function_listed','ctx.gap.gapseq_guessed_listed'],(k,f)
   changes[f]+=1
for relative in ['integrated/raw.jsonl','integrated/protein_integrated.tsv']:
 x,y=before/relative,after/relative
 assert x.read_bytes()==y.read_bytes()
 paths.extend([x,y])
r=json.loads((after/'verification.json').read_text())
assert len(r['runs'])==10 and r['observed_gap_records']==50
result=dict(status='verified',candidates=len(a),changed_columns=dict(changes),raw_claims_byte_identical=True,
            predictions_byte_identical=True,gap_feature_counts=r['gap_feature_counts'],
            inputs={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
(Path.home()/'tmp/gspa-pathway-gap-bridge-v1/comparison.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
