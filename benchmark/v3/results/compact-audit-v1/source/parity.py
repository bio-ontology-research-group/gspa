"""Full small production replay against retained verified pilot; no labels."""
import csv, hashlib, json, subprocess, sys
from pathlib import Path
repo, jar, out, scratch = map(Path,sys.argv[1:])
sys.path.insert(0,str(repo/'benchmark/v3'))
import verified_batch as V
root=Path('/data/hohndor/gspa-v3/price-production-v1')
oldroot=root/'collection/attempts/0/job-10228/production'
reports={}
for arm in ('stage0','stage1'):
    old=oldroot/arm
    receipt=json.loads((old/V.RECEIPT).read_text())
    assert V.file_hashes(receipt['spec']['inputs'])==receipt['spec']['inputs']
    assert V.output_hashes(receipt['spec'])==receipt['artifacts']
    directory=out/arm; directory.mkdir(parents=True,exist_ok=False)
    with (old/'manifest.tsv').open() as stream:
        reader=csv.DictReader(stream,delimiter='\t');fields=reader.fieldnames;rows=list(reader)
    assert len(rows)==1
    oldoutput=Path(rows[0]['out_path']); newoutput=directory/oldoutput.name
    rows[0]['out_path']=str(newoutput); manifest=directory/'manifest.tsv'
    with manifest.open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=fields,delimiter='\t');writer.writeheader();writer.writerows(rows)
    cmd=json.loads((old/'command.json').read_text())
    for flag,value in (('-o',directory/'resources.txt'),('-jar',jar),('--batch-manifest',manifest),('--batch-stats-out',directory/'batch_stats.tsv')):
        cmd[cmd.index(flag)+1]=str(value)
    cmd.insert(cmd.index('java')+1,'-Djava.io.tmpdir='+str(scratch))
    (directory/'command.json').write_text(json.dumps(cmd,indent=2)+'\n')
    fresh=V.run_batch(cmd,manifest,directory,sources=[Path(__file__)],timeout=1800)
    assert V.sha256_file(str(newoutput))==V.sha256_file(str(oldoutput))
    left=json.loads(Path(str(oldoutput)+'.stages.json').read_text())
    right=json.loads(Path(str(newoutput)+'.stages.json').read_text())
    for ledger in (left,right):
        for stage in ledger['stages']:stage.pop('elapsed_seconds',None)
    assert left==right
    reports[arm]=dict(tsv_byte_identical=True,complete_ledger_equal_except_stage_elapsed_seconds=True,
        old_receipt_sha256=V.sha256_file(str(old/V.RECEIPT)),
        new_receipt_sha256=V.sha256_file(str(directory/V.RECEIPT)),
        rows=sum(1 for _ in newoutput.open())-1)
(out/'parity.json').write_text(json.dumps(dict(status='verified',reports=reports),indent=2)+'\n')
