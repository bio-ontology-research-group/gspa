"""Verify the supplemental supervision binding before every resumed stage."""
import hashlib,json,sys
from pathlib import Path
p=Path('/home/leechuck/tmp/context-input-supplement.json')
r=json.loads(p.read_text())
for path,expected in {**r['inputs'],**r['reference_manifests']}.items():
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for block in iter(lambda:f.read(1048576),b''):h.update(block)
 if h.hexdigest()!=expected:raise ValueError('supplemental input changed: '+path)
print(json.dumps({'status':'verified','inputs':len(r['inputs']),'references':len(r['reference_manifests'])}))
