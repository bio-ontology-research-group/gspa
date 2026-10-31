"""Verify portable compact-audit evidence; large production artifacts remain remote."""
import hashlib,json,tarfile,xml.etree.ElementTree as ET
from pathlib import Path
root=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
manifest=json.loads((root/'manifest.json').read_text())
actual={str(p.relative_to(root)) for p in root.rglob('*') if p.is_file() and p.name!='manifest.json'}
assert actual==set(manifest)
for name,digest in manifest.items():assert sha(root/name)==digest,name
sources=json.loads((root/'source/source-files-v2.json').read_text())
with tarfile.open(root/'source/source-v2.tar.gz') as archive:
    names={m.name for m in archive if m.isfile()}
    assert names=={'home/leechuck/Public/software/gspa/'+p for p in sources}
    for name,digest in sources.items():
        f=archive.extractfile('home/leechuck/Public/software/gspa/'+name)
        assert hashlib.sha256(f.read()).hexdigest()==digest,name
job=root/'job-10275'
assert (job/'exit-code.txt').read_text().strip()=='0'
assert (root/'failed-10273/exit-code.txt').read_text().strip()=='1'
counts={}
for phase in ('focused-tests','final-tests'):
    totals=dict(tests=0,failures=0,errors=0,skipped=0)
    files=list((job/phase).rglob('TEST-*.xml'));assert files
    for p in files:
        node=ET.parse(p).getroot()
        for k in totals:totals[k]+=int(node.get(k,0))
    assert totals['failures']==totals['errors']==0
    counts[phase]=totals
assert counts['focused-tests']['tests']==11
before=json.loads((job/'before-memory.json').read_text());after=json.loads((job/'after-memory.json').read_text())
assert before['count']==after['count']==1000000
assert before['json_sha256']==after['json_sha256']
assert 0<after['delta_heap_bytes']<before['delta_heap_bytes']
parity=json.loads((job/'parity/parity.json').read_text())
assert parity['status']=='verified' and set(parity['reports'])=={'stage0','stage1'}
for arm,report in parity['reports'].items():
    assert report['tsv_byte_identical'] is True
    assert report['complete_ledger_equal_except_stage_elapsed_seconds'] is True
    assert report['rows']==3481
    folder=job/'parity'/arm
    old=root/'before-pilot'/arm
    assert sha(old/'verified_run.json')==report['old_receipt_sha256']
    name='GCF_000236665.1.integrated.tsv'
    assert sha(folder/name)==sha(old/name)
    ledgers=[json.loads((d/(name+'.stages.json')).read_text()) for d in (old,folder)]
    for ledger in ledgers:
        for stage in ledger['stages']:stage.pop('elapsed_seconds',None)
    assert ledgers[0]==ledgers[1]
    receipt=json.loads((folder/'verified_run.json').read_text())
    assert sha(folder/'verified_run.json')==report['new_receipt_sha256']
    for path,digest in receipt['artifacts'].items():
        assert sha(folder/Path(path).name)==digest,path
print(json.dumps(dict(status='verified',files=len(manifest),source_files=len(sources),tests=counts,
 memory_delta_bytes=dict(before=before['delta_heap_bytes'],after=after['delta_heap_bytes']),
 full_serialized_million_record_hash_equal=True,pilot_rows_per_arm=3481,
 scope='Portable source/tests/small-output evidence. Dense replay and complete biological benchmark remain separate.'),indent=2))
