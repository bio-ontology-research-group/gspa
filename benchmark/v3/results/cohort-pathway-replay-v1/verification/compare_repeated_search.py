from pathlib import Path
from collections import Counter
import hashlib,json,sys
repo=Path('/home/leechuck/Public/software/gspa')
sys.path.insert(0,str(repo/'benchmark/v3'))
from replay_cohort_pathways import rows,expected_candidates,eligible_groups,canonical,jsonl
base=Path.home()/'.cache/gspa-v3-data';gid='GCA_000012525.1'
a=base/'gapsmith-observed-v2';b=base/'cohort-pathway-search-v1/observations'/gid
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
inputs={}
def bind(p): inputs[str(p)]=sha(p);return p
observations=[json.loads(bind(p/'observation.json').read_text()) for p in [a,b]]
references=[json.loads(bind(p/'reference-hashes.json').read_text()) for p in [a,b]]
rxns=[rows(bind(next((p/'outputs').glob('*-Reactions.tbl')))) for p in [a,b]]
pwys=[rows(bind(next((p/'outputs').glob('*-Pathways.tbl')))) for p in [a,b]]
cs=[Counter(canonical(r)) for r in rxns]
candidates=[expected_candidates(r)[0] for r in rxns]
oldgap=bind(base/'rlc-mapping-gap-bridge-v1/direct.tsv.gaps.jsonl')
newgap=bind(base/'cohort-pathway-replay-v1'/gid/'automatic/predictions.tsv.gaps.jsonl')
gaps=[jsonl(p) for p in [oldgap,newgap]]
# Only the old route's basename-derived genome identifier differs.
normalized=[canonical([{k:v for k,v in row.items() if k!='genome_id'} for row in records]) for records in gaps]
result=dict(schema='gspa.repeated-search-diagnostic/1',scope='Observed repeat difference; upstream cause is not established.',
    observations=[{'job_id':v['job_id'],'command':v['command'],'inputs':v['inputs']} for v in observations],
    references_identical=references[0]==references[1],old_rows=len(rxns[0]),new_rows=len(rxns[1]),
    added=[json.loads(s) for s in (cs[1]-cs[0]).elements()],removed=[json.loads(s) for s in (cs[0]-cs[1]).elements()],
    old_candidates=len(candidates[0]),new_candidates=len(candidates[1]),
    candidates_added=sorted(set(candidates[1])-set(candidates[0])),candidates_removed=sorted(set(candidates[0])-set(candidates[1])),
    shared_candidates_with_changed_evidence=sorted(k for k in candidates[0].keys()&candidates[1].keys() if canonical(candidates[0][k])!=canonical(candidates[1][k])),
    pathway_rows_identical=canonical(pwys[0])==canonical(pwys[1]),
    eligible_groups_identical=eligible_groups(rxns[0],pwys[0])==eligible_groups(rxns[1],pwys[1]),
    production_gaps_equal_except_genome_id=normalized[0]==normalized[1],gap_records=[len(x) for x in gaps],
    inputs=inputs)
assert inputs=={p:sha(Path(p)) for p in inputs}
out=Path(__file__).with_name('repeated-search-difference.json');out.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ['inputs','observations','added','removed']},indent=2))
