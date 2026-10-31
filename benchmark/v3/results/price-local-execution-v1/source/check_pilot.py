import json,sys
from pathlib import Path
import run_price_production_collection as P
root=Path('/data/hohndor/gspa-v3/price-production-v3')
job=sys.argv[1];directory=root/'collection/attempts/0'/('job-'+job)
old=Path('/data/hohndor/gspa-v3/price-production-v1/collection/attempts/0/job-10228')
publication=json.loads((directory/'publication.json').read_text())
assert publication['status']=='verified' and publication['execution']['returncode']==0
assert publication['published_path']==str(directory)
assert P.C.V.sha256_file(str(root/'collection/registration.json'))==publication['execution']['registration_sha256']
for relative,digest in publication['artifacts'].items():
 assert P.C.V.sha256_file(str(directory/relative))==digest,relative
assert not Path(publication['local_path']).exists()
comparisons={}
for arm in ('stage0','stage1'):
 a=directory/'production'/arm;b=old/'production'/arm
 receipt=json.loads((a/P.C.V.RECEIPT).read_text())
 assert P.C.V.output_hashes(receipt['spec'])==receipt['artifacts']
 assert P.C.V.file_hashes(receipt['spec']['inputs'])==receipt['spec']['inputs']
 name='GCF_000236665.1.integrated.tsv'
 assert P.C.V.sha256_file(str(a/name))==P.C.V.sha256_file(str(b/name))
 ledgers=[json.loads((d/(name+'.stages.json')).read_text()) for d in (a,b)]
 for ledger in ledgers:
  for stage in ledger['stages']:stage.pop('elapsed_seconds',None)
 assert ledgers[0]==ledgers[1]
 comparisons[arm]=dict(tsv_byte_identical=True,ledger_equal_except_durations=True,
  receipt_sha256=P.C.V.sha256_file(str(a/P.C.V.RECEIPT)))
P.save(root/'pilot-verification.json',dict(status='verified',job=job,comparisons=comparisons,
 publication_sha256=P.C.V.sha256_file(str(directory/'publication.json')),
 interpretation='Full small production pair after private execution and checked publication. No labels or biological accuracy.'))
