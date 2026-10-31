"""Check stored sources, whole pilot outputs, publication and full scheduling scope.

Large dense artifacts remain remote. Their stored receipts are checked here, not
replaced with a claim that the remote bytes were independently re-read offline.
"""
import hashlib,json,re,tarfile
from pathlib import Path
root=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
manifest=json.loads((root/'manifest.json').read_text())
files={str(p.relative_to(root)) for p in root.rglob('*') if p.is_file() and p.name!='manifest.json'}
assert files==set(manifest)
for name,digest in manifest.items():assert sha(root/name)==digest,name
sources={}
for filename,hashfile,prefix in [('code.tar.gz','code-hashes.json','code/'),('test-source.tar.gz','test-source-hashes.json','test-source/')]:
    expected=json.loads((root/'source'/hashfile).read_text())
    with tarfile.open(root/'source'/filename) as archive:
        assert {m.name for m in archive if m.isfile()}=={prefix+n for n in expected}
        for name,digest in expected.items():
            assert hashlib.sha256(archive.extractfile(prefix+name).read()).hexdigest()==digest,name
    sources[filename]=len(expected)
final=root/'final'; collection=final/'collection'
reg=json.loads((collection/'registration.json').read_text());plan=json.loads((collection/'plan.json').read_text())
assert reg['plan_sha256']==sha(collection/'plan.json')
assert len(plan['cells'])==132
assert len({(c['predictor'],c['assembly']) for c in plan['cells']})==132
assert len({c['predictor'] for c in plan['cells']})==6 and len({c['assembly'] for c in plan['cells']})==22
assert sum(c['claims'] for c in plan['cells'])==plan['expected_claims']==1526919081
assert plan['arms']==['stage0','stage1']
code=json.loads((root/'source/code-hashes.json').read_text())
for path,digest in reg['inputs'].items():
    prefix='/data/hohndor/gspa-v3/price-production-v3/code/'
    if path.startswith(prefix):assert code[path[len(prefix):]]==digest
assert (final/'logs/prepare-10292.exit-code.txt').read_text().strip()=='0'
log=(final/'logs/tests-10292.log').read_text()
assert re.search(r'Ran 64 tests in [0-9.]+s\s+OK\n',log)
assert (root/'failed-package/prepare-10287.exit-code.txt').read_text().strip()=='1'
pilot=json.loads((final/'pilot-verification.json').read_text())
attempt=collection/'attempts/0/job-10292';publication=json.loads((attempt/'publication.json').read_text())
assert publication['status']=='verified' and publication['execution']['returncode']==0
assert publication['execution']['registration_sha256']==sha(collection/'registration.json')
assert pilot['publication_sha256']==sha(attempt/'publication.json')
assert pilot['status']=='verified'
actual={str(p.relative_to(attempt)) for p in attempt.rglob('*') if p.is_file()}
assert actual==set(publication['artifacts'])|{'publication.json'}
for relative,digest in publication['artifacts'].items():assert sha(attempt/relative)==digest,relative
for arm in plan['arms']:
    after=attempt/'production'/arm;before=root/'before-pilot'/arm
    name='GCF_000236665.1.integrated.tsv'
    assert sha(after/name)==sha(before/name)
    ledgers=[json.loads((p/(name+'.stages.json')).read_text()) for p in (before,after)]
    for ledger in ledgers:
        for stage in ledger['stages']:stage.pop('elapsed_seconds',None)
    assert ledgers[0]==ledgers[1]
    assert sum(1 for _ in (after/name).open())-1==3481
    receipt=json.loads((after/'verified_run.json').read_text())
    assert pilot['comparisons'][arm]['receipt_sha256']==sha(after/'verified_run.json')
    for path,digest in receipt['artifacts'].items():assert sha(after/Path(path).name)==digest
schedule=json.loads((final/'schedule.json').read_text())
assert schedule['registration_sha256']==sha(collection/'registration.json')
assert sorted(i for group in schedule['groups'].values() for i in group)==list(range(132))
assert sum(schedule['claims_by_group'].values())==1526919081
for name,indices in schedule['groups'].items():
    assert sum(plan['cells'][i]['claims'] for i in indices)==schedule['claims_by_group'][name]
for path,digest in schedule['files'].items():assert sha(final/Path(path).name)==digest
assert schedule['dependencies']==dict(pilot=10292,compact_dense_parity=10279)
assert all(x['returncode']==0 for x in schedule['submissions'])
large=root/'largest-baseline-10230'
assert json.loads((large/'finished.json').read_text())['status']=='verified'
capacity=json.loads((large/'production/verification.json').read_text())
assert capacity['status']=='verified' and set(capacity['runs'])=={'stage0'}
assert capacity['runs']['stage0']['prediction_rows']==capacity['runs']['stage0']['input_candidates']==45290880
assert capacity['runs']['stage0']['candidate_universe_preserved'] is True
assert capacity['runs']['stage0']['receipt_sha256']==sha(large/'production/stage0/verified_run.json')
print(json.dumps(dict(status='verified',files=len(manifest),sources=sources,remote_tests=64,
 pilot_rows_per_arm=3481,registered_cells=132,scheduled_cells=132,claims=1526919081,
 largest_baseline_rows=45290880,scope='Storage correctness and scheduling checkpoint; full runs and biological accuracy remain incomplete.'),indent=2))
