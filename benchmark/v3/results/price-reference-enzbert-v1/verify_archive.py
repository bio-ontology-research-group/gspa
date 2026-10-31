"""Verify the portable archive and independently recompute its reference metrics."""
import hashlib
import json
from pathlib import Path
import sys
sys.dont_write_bytecode=True
from check_retained import check

root=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
manifest=json.loads((root/'archive-hashes.json').read_text())
actual={str(p.relative_to(root)) for p in root.rglob('*') if p.is_file() and p.name!='archive-hashes.json'}
assert actual==set(manifest), 'archive inventory differs'
for relative,digest in manifest.items():assert sha(root/relative)==digest, relative
for line in (root/'source.sha256').read_text().splitlines():
    digest,relative=line.split('  ',1)
    assert sha(root/'source'/relative)==digest, relative
assert (root/'job-10313/exit-code.txt').read_text().strip()=='0'
tests=(root/'job-10313/tests.log').read_text()
assert 'Ran 15 tests' in tests and tests.rstrip().endswith('OK')
result=check(root/'job-10313/evaluation',root/'ec-audit.json',root/'source/benchmark/manifests/suite6_price149_manifest.tsv')
evaluation=root/'job-10313/evaluation'
registration=json.loads((evaluation/'registration.json').read_text())
report=json.loads((evaluation/'report.json').read_text())
assert report['predictors']==['enzbert'] and report['cells']==22
assert len(registration['inputs'])==506
assert sum(c['claims'] for c in registration['cells'])==11553500
changed={'raw_export_to_stage0':0,'stage0_to_stage1':0}
for relative in report['retained_files']:
    data=json.loads((evaluation/relative).read_text())
    for before,after in (('raw_export','stage0'),('stage0','stage1')):
        for protein,left in data['scores'][before].items():
            right=data['scores'][after][protein]
            top=lambda scores:{k for k,v in scores.items() if v==max(scores.values())}
            assert len(top(left))==1 and len(top(right))==1
            changed[before+'_to_'+after]+=top(left)!=top(right)
assert all(n==0 for n in changed.values())
for mode in ('exact','current_unambiguous'):
    summary=report['summaries'][mode]['enzbert']
    for arm in ('raw_export','stage0','stage1'):
        assert summary['arms'][arm]['top1_match']==53
        assert summary['arms'][arm]['candidate_reference_matches']==96
    for pair,n in changed.items():
        assert summary['paired'][pair]['changed_top1']==n
        assert summary['paired'][pair]['top1_gained_reference_match']==0
        assert summary['paired'][pair]['top1_lost_reference_match']==0
result.update(whole_proteome_claims=11553500,bound_remote_files=506,changed_top1=changed)
result.update(archive_files=len(manifest),remote_tests=15)
print(json.dumps(result,indent=2))
