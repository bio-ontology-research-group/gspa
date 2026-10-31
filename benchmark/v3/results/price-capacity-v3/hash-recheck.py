import hashlib,json,datetime
from pathlib import Path

def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024*1024), b''): h.update(chunk)
    return h.hexdigest()
root=Path('/data/hohndor/gspa-v3/price-capacity-v2/dpz-GCF_000236665.1-stage0')
receipt=json.loads((root/'stage0/verified_run.json').read_text())
expected={**receipt['spec']['inputs'],**receipt['artifacts']}
for name in ['verification.json','stage0/verified_run.json','stage0/command.json','stage0/resources.txt']:
    expected[str(root/name)]=digest(root/name)
expected[str(Path(__file__).resolve())]=digest(__file__)
verification=json.loads((root/'verification.json').read_text())
if verification['status']!='verified': raise ValueError('capacity measurement did not pass')
for name,sha in expected.items():
    if digest(name)!=sha: raise ValueError('changed recorded file: '+name)
out=Path('/data/hohndor/gspa-v3/price-capacity-v2/dpz-postrun-hashes.json')
if out.exists(): raise SystemExit('recheck already exists')
out.write_text(json.dumps({'status':'verified','captured_at':datetime.datetime.now(datetime.timezone.utc).isoformat(),
  'hashes':expected,'scope':'Post-measurement hash recheck of all recorded inputs and output TSV/ledger; candidate-universe comparison was performed by the completed capacity job.'},indent=2)+'\n')
print('Verified',len(expected),'recorded hashes after completed measurement')
