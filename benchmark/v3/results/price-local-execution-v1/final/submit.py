"""Register resource waves for every cell, with dependencies on validation."""
import hashlib,json,subprocess
from pathlib import Path
root=Path('/data/hohndor/gspa-v3/price-production-v3')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
plan=json.loads((root/'collection/plan.json').read_text())
pilot=json.loads((root/'pilot-verification.json').read_text())
assert pilot['status']=='verified' and pilot['job']=='10292'
groups=dict(pilot=[0],small=list(range(1,66)),medium=list(range(66,101)),large=list(range(101,131)),largest=[131])
assert sorted(i for cells in groups.values() for i in cells)==list(range(132))
assert sum(c['claims'] for c in plan['cells'])==1526919081
record=dict(schema='gspa.price-resource-waves/1',groups=groups,
 claims_by_group={k:sum(plan['cells'][i]['claims'] for i in v) for k,v in groups.items()},
 registration_sha256=sha(root/'collection/registration.json'),
 files={str(root/name):sha(root/name) for name in ('small.sbatch','medium.sbatch','large.sbatch','submit.py')},
 dependencies=dict(pilot=10292,compact_dense_parity=10279),submissions=[],
 interpretation='All132cells andbotharms remainrequired. Resourcewaves do not select bylabels. No ensemble/cohort or accuracycompletion follows.')
with (root/'schedule-started.json').open('x') as stream:json.dump(record,stream,indent=2)
def submit(kind,args):
 command=['sbatch','--parsable']+args
 result=subprocess.run(command,capture_output=True,text=True)
 entry=dict(kind=kind,command=command,returncode=result.returncode,stdout=result.stdout,stderr=result.stderr)
 record['submissions'].append(entry)
 with (root/'schedule-progress.json').open('w') as stream:json.dump(record,stream,indent=2)
 result.check_returncode()
 job=result.stdout.strip().split(';')[0]
 assert job.isdigit(),job
 return job
largest=submit('largest',['--dependency=afterok:10292:10279',str(root/'large.sbatch'),'131'])
small=submit('small',['--dependency=afterok:10292:10279','--array=1-65%2',str(root/'small.sbatch')])
medium=submit('medium',['--dependency=afterok:10292:10279','--array=66-100%2',str(root/'medium.sbatch')])
large=submit('large',['--dependency=afterok:'+largest,'--array=101-130%1',str(root/'large.sbatch')])
record['jobs']=dict(largest=largest,small=small,medium=medium,large=large)
with (root/'schedule.json').open('x') as stream:json.dump(record,stream,indent=2)
print(json.dumps(record['jobs']))
