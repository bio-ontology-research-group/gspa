import json, os, resource, time, sys
from pathlib import Path
import verified_batch as V
import ijson
out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=False)
old=Path('/data/hohndor/gspa-v3/price-index-capacity-v1/dpz-GCF_000236665.1-stage1/stage1/verified_run.json')
r=json.loads(old.read_text())
started=time.monotonic()
actual=V.output_hashes(r['spec'])
assert actual==r['artifacts'], 'retained output receipts changed'
summary=dict(status='verified',receipt=str(old),receipt_sha256=V.sha256_file(str(old)),
 artifacts=actual,elapsed_seconds=time.monotonic()-started,
 max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
 parser_version=ijson.__version__,parser_backend=ijson.backend,
 code_hashes=V.file_hashes([Path(V.__file__),Path(__file__),*V.streaming_parser_sources()]),
 interpretation='Existing production TSV/ledger hashes preserved. No inference or accuracy measurement.')
(out/'result.json').write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n')
print(json.dumps(summary,sort_keys=True))
