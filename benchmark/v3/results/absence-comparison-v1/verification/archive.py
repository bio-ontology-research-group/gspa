from pathlib import Path
import gzip, hashlib, json, shutil

repo = Path('/home/leechuck/Public/software/gspa')
cache = Path.home()/'.cache/gspa-v3-data'
scratch = Path(__file__).parent

def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as stream:
        for block in iter(lambda:stream.read(1<<20),b''):h.update(block)
    return h.hexdigest()

def archive(name, source, files, include_tests):
    out=repo/'benchmark/v3/results'/name
    out.mkdir(exist_ok=False)
    receipt=json.loads((source/'verification.json').read_text())
    for p,h in receipt['inputs'].items():
        assert digest(Path(p)) == h, p
    for p,h in receipt['outputs'].items():
        assert digest(source/p) == h, p
    index={}
    def add(p,key):
        size=p.stat().st_size;h=digest(p)
        stored=key+'.gz' if size>512*1024 else key
        dest=out/stored;dest.parent.mkdir(parents=True,exist_ok=True)
        if stored!=key:
            with p.open('rb') as src,dest.open('wb') as dst,gzip.GzipFile(filename='',mode='wb',compresslevel=1,mtime=0,fileobj=dst) as z:
                shutil.copyfileobj(src,z,1<<20)
        else:shutil.copyfile(p,dest)
        index[key]=dict(original=str(p),bytes=size,sha256=h,stored=stored)
    for p in sorted(source.rglob('*')):
        if p.is_file():add(p,'cli/'+str(p.relative_to(source)))
    for name in files:add(repo/name,'source/'+name)
    if include_tests:
        for p in sorted(scratch.rglob('*')):
            if p.is_file() and not p.name.startswith('archive.log') and p.name not in ['full-pathway-replay.log','full-pathway-summary.json']:
                add(p,'verification/'+str(p.relative_to(scratch)))
    else:
        for p in ['full-pathway-replay.log','full-pathway-summary.json','archive.py']:
            add(scratch/p,'verification/'+p)
    (out/'storage-index.json').write_text(json.dumps(index,indent=2,sort_keys=True)+'\n')
    for key,record in index.items():
        h=hashlib.sha256();opener=gzip.open if key!=record['stored'] else open
        with opener(out/record['stored'],'rb') as f:
            for block in iter(lambda:f.read(1<<20),b''):h.update(block)
        assert h.hexdigest()==record['sha256'],key
    report=dict(status='verified',archive=name,payloads=len(index),original_bytes=sum(r['bytes'] for r in index.values()),
                stored_payload_bytes=sum((out/r['stored']).stat().st_size for r in index.values()))
    (out/'archive-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True)

archive('absence-comparison-v1',cache/'absence-mechanism-v1',[
    'scripts/absence_experiments.py','scripts/metabolic_experiments.py','scripts/tests/test_absence_experiments.py',
    'benchmark/v3/check_absence_cobra.py','benchmark/v3/run_absence_mechanism.py','benchmark/v3/tests/test_absence_cobra.py',
    'gspa-core/src/main/groovy/gspa/integration/contextbridge/AbsenceComparisonBridge.groovy',
    'gspa-core/src/main/groovy/gspa/integration/contextbridge/ExperimentStage.groovy',
    'gspa-core/src/test/groovy/gspa/integration/workflow/ExperimentWorkflowSpec.groovy',
    'gspa-core/src/test/resources/experiments/make_fixtures.py','docs/benchmarks/v3-absence-comparison.md',
    'tasks/todo.md','tasks/v3-acceptance.md'],True)
archive('cohort-pathway-replay-v2',cache/'cohort-pathway-replay-v2',[
    'benchmark/v3/replay_cohort_pathways.py','benchmark/v3/audit_cohort_pathways.py','benchmark/v3/observe_gapsmith.py',
    'benchmark/v3/tests/test_replay_cohort_pathways.py','docs/benchmarks/v3-cohort-pathway-replay.md'],False)
