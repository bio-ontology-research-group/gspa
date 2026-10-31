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
assert (root/'job-10306/exit-code.txt').read_text().strip()=='0'
tests=(root/'job-10306/tests.log').read_text()
assert 'Ran 15 tests' in tests and tests.rstrip().endswith('OK')
result=check(root/'job-10306/evaluation',root/'ec-audit.json',root/'source/benchmark/manifests/suite6_price149_manifest.tsv')
result.update(archive_files=len(manifest),remote_tests=15)
print(json.dumps(result,indent=2))
