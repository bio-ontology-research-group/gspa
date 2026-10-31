import json, os, sys, tempfile
from pathlib import Path
sys.path.insert(0, '/data/hohndor/gspa-v3/price-capacity-v2/code')
import verified_batch as V
import verify_candidate_universe as U
part=Path('/data/hohndor/gspa-v3/price-full-matrix-v2/parts/enzbert/GCF_000236665.1')
run=Path('/data/hohndor/gspa-v3/price-capacity-v1/enzbert-GCF_000236665.1')
out=Path('/data/hohndor/gspa-v3/price-capacity-v2/enzbert-independent-recheck.json')
if out.exists(): raise SystemExit('recheck record already exists')
publication=json.loads((part/'job-receipt.json').read_text())
expected={str(part/name):sha for name,sha in publication['artifacts'].items()}
for arm in ['stage0','stage1']:
    receipt=json.loads((run/arm/'verified_run.json').read_text())
    for name in ['verified_run.json','resources.txt','command.json']:
        path=run/arm/name
        expected[str(path)]=V.sha256_file(str(path))
    for path,sha in {**receipt['spec']['inputs'],**receipt['artifacts']}.items():
        if path in expected and expected[path]!=sha: raise ValueError('inconsistent receipts')
        expected[path]=sha
for path in [part/'job-receipt.json',run/'verification.json',Path(__file__),Path(U.__file__),U.sort_executable()]:
    expected[str(path)]=V.sha256_file(str(path))
if V.file_hashes(expected)!=expected: raise ValueError('recorded input/output changed')
counts={}
with tempfile.TemporaryDirectory(dir=os.environ['TMPDIR']) as temp:
    baseline,count=U.expected_candidates(part/'claims.jsonl',Path(temp))
    for arm in ['stage0','stage1']:
        counts[arm]=U.compare_candidates(baseline,run/arm/'GCF_000236665.1.integrated.tsv',Path(temp))
if V.file_hashes(expected)!=expected: raise ValueError('input/output changed during verification')
report={'status':'verified','input_candidates':count,'arm_counts':counts,'hashes':expected,
        'method':'Every candidate key externally sorted and compared exactly; no sampling; all receipt inputs/outputs rehashed.'}
out.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
print(json.dumps({'status':report['status'],'input_candidates':count,'arm_counts':counts,'hashed_files':len(expected)}))
