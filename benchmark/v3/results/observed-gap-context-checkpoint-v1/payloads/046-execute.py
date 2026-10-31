import json,subprocess,os
from datetime import datetime,timezone
from pathlib import Path
root=Path(__file__).parent;repo=Path('/home/leechuck/Public/software/gspa');cache=Path.home()/'.cache/gspa-v3-data'
python='/home/leechuck/tmp/gspa-v3-testenv/bin/python';out=cache/'observed-gap-context-modes-v1'
jar=cache/'engines/b23a7f6d94eb06c375e63ad344cfbf4c15d7644e33c1848766a36c4287cfd090.jar'
stages=[dict(name=action,command=[python,str(repo/'benchmark/v3/run_observed_gap_context_modes.py'),action,'--out',str(out)],status='queued') for action in ('export','train','evaluate')]
for transfer in ('homology','family_excluded'):
 stages.append(dict(name='parity-'+transfer,command=[python,str(repo/'benchmark/v3/check_context_mode_models.py'),'--exports',str(out/transfer),'--jar',str(jar),'--out',str(cache/('observed-gap-context-modes-v1-parity-'+transfer))],status='queued'))
def save():(root/'workflow.json').write_text(json.dumps(stages,indent=2)+'\n')
save()
for stage in stages:
 stage.update(status='running',started_at=datetime.now(timezone.utc).isoformat());save()
 with (root/(stage['name']+'.log')).open('x') as log:
  r=subprocess.run(stage['command'],cwd=repo,stdout=log,stderr=subprocess.STDOUT)
 stage.update(status='completed' if r.returncode==0 else 'failed',returncode=r.returncode,finished_at=datetime.now(timezone.utc).isoformat());save()
 print(json.dumps(stage),flush=True)
 if r.returncode:raise SystemExit(r.returncode)
