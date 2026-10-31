"""Label-free largest-cell capacity diagnostic; never a complete collection cell."""
import json, os, subprocess, sys, threading, time, traceback
from pathlib import Path
import run_price_matrix_capacity as C
root=Path('/data/hohndor/gspa-v3/price-production-v1')
plan=json.loads((root/'collection/plan.json').read_text())
reg=json.loads((root/'collection/registration.json').read_text())
assert C.V.file_hashes(reg['inputs'])==reg['inputs']
assert C.V.sha256_file(str(root/'collection/plan.json'))==reg['plan_sha256']
cell=max(plan['cells'],key=lambda c:c['claims'])
assert cell['predictor']=='clean' and cell['claims']==45290880
out=root/'dense'/('job-'+os.environ['SLURM_JOB_ID']);out.mkdir(parents=True,exist_ok=False)
profile=out/'profile.jfr'
os.environ['JAVA_TOOL_OPTIONS']=f'-XX:StartFlightRecording=filename={profile},settings=profile,maxsize=256m,dumponexit=true'
(out/'started.json').write_text(json.dumps(dict(cell=cell,java_tool_options=os.environ['JAVA_TOOL_OPTIONS'],
 heap='168g',allocation_mib=196608,timeout_seconds=36000,arms=['stage0'],
 inputs=C.V.file_hashes([root/'collection/registration.json',Path(__file__),root/'dense.sbatch']),
 scope='Largest dense single-arm capacity diagnostic. Does not complete either the cell or the collection.'),indent=2)+'\n')
stop=threading.Event()
def monitor():
    cg={}
    for line in Path('/proc/self/cgroup').read_text().splitlines():
        _,controllers,path=line.split(':',2)
        if 'memory' in controllers.split(','):
            folder=Path('/sys/fs/cgroup/memory')/path.lstrip('/')
            cg={name:folder/name for name in ('memory.usage_in_bytes','memory.max_usage_in_bytes','memory.limit_in_bytes','memory.failcnt','memory.oom_control')}
    with (out/'memory.jsonl').open('w',buffering=1) as log:
        while not stop.is_set():
            mem={s.split(':')[0]:s.split(':')[1].strip() for s in Path('/proc/meminfo').read_text().splitlines() if s.startswith(('MemAvailable:','MemFree:','SwapFree:'))}
            values={k:p.read_text().strip() for k,p in cg.items() if p.is_file()}
            log.write(json.dumps(dict(time=time.time(),node_memory=mem,cgroup_memory=values))+'\n')
            stop.wait(15)
thread=threading.Thread(target=monitor,daemon=True);thread.start()
status='failed'
try:
    C.run(Path(cell['part']),Path(plan['layout']),Path(plan['jar']),Path(plan['reference']),out/'production',
          scratch=Path(os.environ['TMPDIR']),java_heap='168g',arms=('stage0',),timeout=36000)
    status='verified'
except BaseException:
    (out/'failure.log').write_text(traceback.format_exc())
    raise
finally:
    stop.set();thread.join(timeout=20)
    (out/'finished.json').write_text(json.dumps(dict(status=status,time=time.time(),profile_exists=profile.exists(),
        scope='Capacity only; no labels, no function accuracy and no complete two-arm cell.'),indent=2)+'\n')
