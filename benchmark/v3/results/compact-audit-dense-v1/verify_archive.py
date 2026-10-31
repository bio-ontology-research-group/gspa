"""Verify retained dense-parity receipt chains, not absent large remote artifacts."""
import hashlib
import json
from pathlib import Path

root=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
read=lambda p:json.loads(p.read_text())
manifest=read(root/'archive-hashes.json')
assert set(manifest)=={str(p.relative_to(root)) for p in root.rglob('*') if p.is_file() and p.name!='archive-hashes.json'}
for relative,digest in manifest.items():assert sha(root/relative)==digest,relative
state=(root/'scheduler.txt').read_text()
assert 'JobState=COMPLETED' in state and 'ExitCode=0:0' in state
comparison=read(root/'job-10279/comparison.json')
assert comparison['status']=='verified'
old=root/'baseline-dpz-10235/production';new=root/'job-10279/production'
reports={'old':read(old/'verification.json'),'new':read(new/'verification.json')}
resources={}
for arm,record in comparison['comparisons'].items():
    resources[arm]={}
    for label,directory in [('old',old),('new',new)]:
        path=directory/arm/'verified_run.json';receipt=read(path)
        assert sha(path)==record[label+'_receipt_sha256']
        assert reports[label]['runs'][arm]['receipt_sha256']==sha(path)
        assert receipt['returncode']==0
        path=receipt['spec']['outputs']['GCF_000236665.1']
        assert receipt['artifacts'][path]==record['tsv_sha256']
        report=reports[label]['runs'][arm]
        assert report['candidate_universe_preserved'] is True
        assert report['prediction_rows']==report['input_candidates']==9790188
        resources[arm][label]=int(report['resources']['Maximum resident set size (kbytes)'])
    resources[arm]['rss_reduction_percent']=100*(1-resources[arm]['new']/resources[arm]['old'])
assert set(resources)=={'stage0','stage1'}
print(json.dumps(dict(status='verified',files=len(manifest),candidates_per_arm=9790188,resources=resources,
    scope='Portable receipt chains only; full candidate/TSV/normalized-ledger equality was checked by completed remote job10279.'),indent=2))
