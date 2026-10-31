"""Verify retained checkpoint bytes and pilot receipts; bulk collection stays remote."""
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parent
REMOTE='/data/hohndor/gspa-v3/price-production-v1/'
def sha(path):
    digest=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):digest.update(block)
    return digest.hexdigest()
def read(path):return json.loads(path.read_text())

def verify():
    expected=read(ROOT/'archive-hashes.json')
    actual={str(p.relative_to(ROOT)):sha(p) for p in ROOT.rglob('*') if p.is_file() and p.name!='archive-hashes.json'}
    assert actual==expected,'archive file inventory/hash mismatch'
    plan=read(ROOT/'collection/plan.json');reg=read(ROOT/'collection/registration.json')
    assert sha(ROOT/'collection/plan.json')==reg['plan_sha256']
    cells=plan['cells'];keys={(c['predictor'],c['assembly']) for c in cells}
    assert len(cells)==len(keys)==132 and len({g for _,g in keys})==22
    assert len({p for p,_ in keys})==6
    assert keys=={(p,g) for p,_ in keys for _,g in keys}
    assert plan['arms']==['stage0','stage1']
    assert sum(c['claims'] for c in cells)==plan['expected_claims']==1526919081
    sparse=[c for c in cells if c['predictor'] in ('topec','mapred','enzbert')]
    assert len(sparse)==66 and sum(c['claims'] for c in sparse)==13682870
    for remote,digest in reg['inputs'].items():
        if remote.startswith(REMOTE+'code/'):
            assert sha(ROOT/'source'/remote.removeprefix(REMOTE))==digest,remote
    replay=read(ROOT/'replay-10226/result.json')
    assert replay['status']=='verified' and replay['parser_backend']=='yajl2_c'
    assert replay['max_rss_kib']==23424
    assert replay['artifacts']==read(ROOT/'replay-10225/result.json')['artifacts']
    for remote,digest in replay['code_hashes'].items():
        if remote.startswith(REMOTE+'code/') or remote==REMOTE+'replay_ledger.py':
            assert sha(ROOT/'source'/remote.removeprefix(REMOTE))==digest
    pilot=ROOT/'pilot-10228';start=read(pilot/'started.json');done=read(pilot/'completed.json')
    assert start['cell']==done['cell']==cells[0]
    assert start['registration_sha256']==sha(ROOT/'collection/registration.json')
    assert done['status']=='verified' and not (pilot/'failed.json').exists()
    result=read(pilot/'production/verification.json')
    assert sha(pilot/'production/verification.json')==done['verification_sha256']
    assert result['status']=='verified' and set(result['runs'])==set(plan['arms'])
    for arm,report in result['runs'].items():
        assert report['prediction_rows']==report['input_candidates']==cells[0]['claims']==3481
        assert report['candidate_universe_preserved'] is True
        directory=pilot/'production'/arm;receipt=read(directory/'verified_run.json')
        assert sha(directory/'verified_run.json')==report['receipt_sha256']
        assert receipt['returncode']==0
        for remote,digest in receipt['artifacts'].items():
            local=directory/Path(remote).name
            assert sha(local)==digest,remote
        command=receipt['spec']['command']
        assert (arm=='stage0')==('--enable-priors' in command)
    print(json.dumps(dict(status='verified',archive_files=len(expected),registered_cells=len(cells),
        completed_cells_retained=1,full_collection_complete=False,
        interpretation='Retained bytes and pilot receipts checked; remote bulk inputs and ongoing outputs not rehashed here.'),sort_keys=True))
if __name__=='__main__':verify()
