import hashlib,json
from pathlib import Path
root=Path('/data/hohndor/gspa-v3/price-production-v3')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
r=root/'collection/registration.json'
assert sha(r)=='454f466cb635b70a9bbb65cbb70d1022065568d1931ac944be95e27a754b02dd'
reg=json.loads(r.read_text())
assert sha(root/'collection/plan.json')==reg['plan_sha256']=='2407147c51dd0fe25cea0ea8193feefcac75d8cdbe46f0009d9001474339c4a4'
files=list((root/'code').glob('*.py'))
expected={Path(name) for name in reg['inputs'] if Path(name).parent==root/'code' and Path(name).suffix=='.py'}
assert expected and set(files)==expected
for p in files:assert sha(p)==reg['inputs'][str(p)],str(p)
assert not (root/'collection/collection-verification.json').exists()
print('Original registration, plan and all',len(files),'registered Python source files verified.')
