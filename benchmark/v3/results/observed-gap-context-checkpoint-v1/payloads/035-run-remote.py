"""Resume registered context experiment on a Slurm compute node."""
import json,os,subprocess,sys
from datetime import datetime,timezone
from pathlib import Path
repo=Path('/home/leechuck/Public/software/gspa')
cache=Path('/home/leechuck/.cache/gspa-v3-data')
training=cache/'observed-gap-context-modes-v1'
replay=cache/'observed-gap-context-model-replay-v1'
analysis=cache/'observed-gap-context-model-analysis-v1'
logroot=Path('/home/leechuck/tmp/gspa-context-remote-execution');logroot.mkdir(exist_ok=False)
python='/home/leechuck/tmp/gspa-v3-testenv/bin/python'
jar=cache/'engines/b23a7f6d94eb06c375e63ad344cfbf4c15d7644e33c1848766a36c4287cfd090.jar'
commands=[(a,[python,str(repo/'benchmark/v3/run_observed_gap_context_modes.py'),a,'--out',str(training)]) for a in ('export','train','evaluate')]
for transfer in ('homology','family_excluded'):
 commands.append(('parity-'+transfer,[python,str(repo/'benchmark/v3/check_context_mode_models.py'),'--exports',str(training/transfer),'--jar',str(jar),'--out',str(cache/('observed-gap-context-modes-v1-parity-'+transfer))]))
commands += [('production-replay',[python,str(repo/'benchmark/v3/check_observed_gap_context_models.py'),'--training',str(training),'--out',str(replay)]),
 ('analysis-registration',[python,str(repo/'benchmark/v3/analyze_observed_gap_context_models.py'),'register','--training',str(training),'--replay',str(replay),'--out',str(analysis)]),
 ('analysis',[python,str(repo/'benchmark/v3/analyze_observed_gap_context_models.py'),'run','--out',str(analysis)])]
stages=[dict(name=n,command=c,status='queued') for n,c in commands]
def save():(logroot/'workflow.json').write_text(json.dumps(dict(slurm_job_id=os.environ.get('SLURM_JOB_ID'),hostname=os.uname().nodename,stages=stages),indent=2)+'\n')
save()
for stage in stages:
 stage.update(status='running',started_at=datetime.now(timezone.utc).isoformat());save()
 with (logroot/(stage['name']+'.log')).open('x') as log:
  result=subprocess.run(stage['command'],cwd=repo,stdout=log,stderr=subprocess.STDOUT)
 stage.update(status='completed' if result.returncode==0 else 'failed',returncode=result.returncode,finished_at=datetime.now(timezone.utc).isoformat());save()
 print(json.dumps(stage),flush=True)
 if result.returncode:raise SystemExit(result.returncode)
