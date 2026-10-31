import json, resource, sys, time
from pathlib import Path
import verified_batch as V
import v3_runner as R
import stage_ledger as L
root=Path(__file__).resolve().parent
out=root/('job-'+sys.argv[1]);out.mkdir(exist_ok=False)
inputs=[p for p in (root/'source').rglob('*') if p.is_file() and '__pycache__' not in p.parts]
inputs += [Path(__file__),root/'replay.sbatch',root/'source.tar.gz',*L.streaming_parser_sources()]
start_hashes=V.file_hashes(inputs)
old=Path('/data/hohndor/gspa-v3/price-index-capacity-v1/dpz-GCF_000236665.1-stage1/stage1/verified_run.json')
receipt=json.loads(old.read_text());gid='GCF_000236665.1';prediction=receipt['spec']['outputs'][gid]
start=time.monotonic();batch=V.output_hashes(receipt['spec']);batch_seconds=time.monotonic()-start
assert batch==receipt['artifacts'],'verified-batch output hashes differ'
start=time.monotonic();generic=R.stage_artifacts(gid,prediction,{'gspa':{'shared_args':[]}}, {'args':[]});generic_seconds=time.monotonic()-start
assert {p:r['sha256'] for p,r in generic.items()}=={p:h for p,h in batch.items() if p!=prediction},'generic sidecar hashes differ'
assert R.sha256_file(prediction)==receipt['artifacts'][prediction]
assert V.file_hashes(inputs)==start_hashes,'replay implementation changed'
result=dict(status='verified',old_receipt=str(old),old_receipt_sha256=V.sha256_file(str(old)),
 input_hashes=start_hashes,parser=L.parser_identity(),batch_artifacts=batch,generic_artifacts=generic,
 elapsed_seconds={'verified_batch':batch_seconds,'generic_runner':generic_seconds},
 peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
 interpretation='Both new validation entry points agree with the existing real dense production receipt. No inference, labels, model fitting or score changes.')
(out/'result.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ('input_hashes','batch_artifacts','generic_artifacts')},sort_keys=True))
