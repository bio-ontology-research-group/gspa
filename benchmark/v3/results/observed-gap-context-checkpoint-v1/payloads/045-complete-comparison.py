"""Continue only after the registered five-stage context pipeline completes."""
import hashlib,json,os,subprocess,sys,time
from datetime import datetime,timezone
from pathlib import Path
repo=Path('/home/leechuck/Public/software/gspa'); scratch=Path(__file__).parent
sys.path.insert(0,str(repo/'benchmark/v3'))
import check_observed_gap_context_models as C
import analyze_observed_gap_context_models as A
R=C.R; F=R.F
cache=Path.home()/'.cache/gspa-v3-data'
training=cache/'observed-gap-context-modes-v1'
replay=cache/'observed-gap-context-model-replay-v1'
analysis=cache/'observed-gap-context-model-analysis-v1'
python='/home/leechuck/tmp/gspa-v3-testenv/bin/python'
pid=2926626
proc=Path(f'/proc/{pid}/stat')
start_identity=proc.read_text().split(') ',1)[1].split()[19]
inputs=F.hashes(R.sources()+[Path(__file__),scratch/'execute.py'])
stages=[dict(name='production-replay',command=[python,str(repo/'benchmark/v3/check_observed_gap_context_models.py'),'--training',str(training),'--out',str(replay)],status='queued'),
 dict(name='analysis-registration',command=[python,str(repo/'benchmark/v3/analyze_observed_gap_context_models.py'),'register','--training',str(training),'--replay',str(replay),'--out',str(analysis)],status='queued'),
 dict(name='analysis',command=[python,str(repo/'benchmark/v3/analyze_observed_gap_context_models.py'),'run','--out',str(analysis)],status='queued')]
state=dict(status='waiting',dependency_pid=pid,dependency_start_identity=start_identity,source_inputs=inputs,stages=stages)
state_path=scratch/'comparison-workflow.json'
def save():
 temp=state_path.with_suffix('.pending');temp.write_text(json.dumps(state,indent=2)+'\n');temp.replace(state_path)
def alive():
 try:return proc.read_text().split(') ',1)[1].split()[19]==start_identity
 except FileNotFoundError:return False
save(); deadline=time.monotonic()+6*3600
try:
 while True:
  try:previous=json.loads((scratch/'workflow.json').read_text())
  except json.JSONDecodeError:
   if not alive():raise ValueError('dependency exited with malformed status')
   time.sleep(1);continue
  if [s['name'] for s in previous]!=['export','train','evaluate','parity-homology','parity-family_excluded']:raise ValueError('dependency stage set changed')
  if any(s['status']=='failed' for s in previous):raise ValueError('dependency stage failed')
  if all(s['status']=='completed' and s.get('returncode')==0 for s in previous):break
  if not alive():raise ValueError('dependency exited before completing all stages')
  if time.monotonic()>deadline:raise TimeoutError('dependency still incomplete after six hours; no replay launched')
  time.sleep(15)
 F.check_hashes(inputs);R.verify(training)
 for transfer in R.G.MODES:
  root=cache/('observed-gap-context-modes-v1-parity-'+transfer)
  receipt=json.loads((root/'verification.json').read_text())
  if receipt.get('status')!='verified' or len(receipt['results'])!=6:raise ValueError('incomplete table parity')
  F.check_hashes(receipt['inputs']);F.check_hashes({str(root/p):v for p,v in receipt['outputs'].items()})
 state.update(status='running',dependency_completed_at=datetime.now(timezone.utc).isoformat());save()
 for stage in stages:
  F.check_hashes(inputs)
  stage.update(status='running',started_at=datetime.now(timezone.utc).isoformat());save()
  with (scratch/(stage['name']+'.log')).open('x') as log:
   result=subprocess.run(stage['command'],cwd=repo,stdout=log,stderr=subprocess.STDOUT)
  stage.update(status='completed' if result.returncode==0 else 'failed',returncode=result.returncode,finished_at=datetime.now(timezone.utc).isoformat());save()
  if result.returncode:raise RuntimeError('comparison stage failed: '+stage['name'])
 F.check_hashes(inputs);state['status']='completed';save()
except BaseException as exc:
 state.update(status='failed',error=repr(exc));save();raise
