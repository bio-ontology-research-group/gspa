"""Complete DPZ replay with compact audit; no candidate sampling or labels."""
import hashlib,json,re,sys
from pathlib import Path
root=Path('/data/hohndor/gspa-v3/compact-audit-v1')
source,jar,out,scratch=map(Path,sys.argv[1:])
sys.path.insert(0,str(source/'benchmark/v3'))
import run_price_matrix_capacity as C
old=Path('/data/hohndor/gspa-v3/price-production-v1/collection/attempts/66/job-10235/production')
reference=Path('/data/hohndor/gspa-v3/price-production-v1/reference')
part=Path('/data/hohndor/gspa-v3/price-full-matrix-v2/parts/dpz/GCF_000236665.1')
layout=Path('/data/hohndor/gspa-v3/price-full-matrix-v2/layout')
old_receipts={arm:json.loads((old/arm/C.V.RECEIPT).read_text()) for arm in ('stage0','stage1')}
# The wrapper and environment are frozen explicitly, in addition to run_batch's command receipts.
inputs=C.V.file_hashes([Path(__file__),root/'dense-parity.sbatch',root/'source-v2.tar.gz',jar])
out.mkdir(parents=True,exist_ok=False)
C.save(out/'registration.json',dict(inputs=inputs,arms=['stage0','stage1'],
    java_tool_options=__import__('os').environ.get('JAVA_TOOL_OPTIONS'),
    comparisons='Byte-identical TSV and every ledger line except exactly two stage elapsed_seconds lines. Same candidate universe; no labels.'))
for arm,receipt in old_receipts.items():
    if C.V.output_hashes(receipt['spec'])!=receipt['artifacts']: raise ValueError('old output changed')
    if C.V.file_hashes(receipt['spec']['inputs'])!=receipt['spec']['inputs']: raise ValueError('old input changed')
C.run(part,layout,jar,reference,out/'production',scratch=scratch,java_heap='64g',arms=('stage0','stage1'),timeout=10800)
def ledger_hash(path):
    h=hashlib.sha256(); skipped=0
    with path.open('rb') as stream:
        for line in stream:
            if re.fullmatch(rb'    "elapsed_seconds" : [0-9.eE+\-]+,?\r?\n',line):
                skipped+=1
            else:h.update(line)
    if skipped!=2:raise ValueError(f'Unexpected timed stage count: {skipped}')
    return h.hexdigest()
comparisons={}
for arm,receipt in old_receipts.items():
    name='GCF_000236665.1.integrated.tsv'
    before=old/arm/name;after=out/'production'/arm/name
    if C.V.sha256_file(str(before))!=C.V.sha256_file(str(after)):raise ValueError('prediction changed')
    lh=ledger_hash(Path(str(before)+'.stages.json'));rh=ledger_hash(Path(str(after)+'.stages.json'))
    if lh!=rh:raise ValueError('ledger changed beyond elapsed time')
    comparisons[arm]=dict(tsv_sha256=C.V.sha256_file(str(after)),normalized_ledger_sha256=rh,
        old_receipt_sha256=C.V.sha256_file(str(old/arm/C.V.RECEIPT)),
        new_receipt_sha256=C.V.sha256_file(str(out/'production'/arm/C.V.RECEIPT)))
if C.V.file_hashes(inputs)!=inputs:raise ValueError('comparison inputs changed')
C.save(out/'comparison.json',dict(status='verified',comparisons=comparisons,
    interpretation='Complete paired implementation parity. Not an accuracy result or a completed 132-cell collection.'))
