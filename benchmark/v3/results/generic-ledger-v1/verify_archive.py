"""Portable checks of retained source, tests and real-output audit evidence."""
import hashlib
import json
from pathlib import Path
import tarfile
ROOT=Path(__file__).resolve().parent
PREFIX='/data/hohndor/gspa-v3/generic-ledger-v2/'
def sha_bytes(data):return hashlib.sha256(data).hexdigest()
def sha(path):
    with path.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()
def read(path):return json.loads(path.read_text())
def verify():
    expected=read(ROOT/'archive-hashes.json')
    actual={str(p.relative_to(ROOT)):sha(p) for p in ROOT.rglob('*') if p.is_file() and p.name!='archive-hashes.json'}
    assert actual==expected,'archive inventory or hash changed'
    result=read(ROOT/'remote/job-10244/result.json'); prior=read(ROOT/'prior-replay.json')
    assert result['status']==prior['status']=='verified'
    assert result['batch_artifacts']==prior['artifacts']
    assert result['old_receipt_sha256']==prior['receipt_sha256']
    assert {p:r['sha256'] for p,r in result['generic_artifacts'].items()}=={
        p:h for p,h in result['batch_artifacts'].items() if p.endswith('.stages.json')}
    assert result['peak_rss_kib']==24400
    assert result['parser']=={'name':'ijson','version':'3.4.0.post0','backend':'yajl2_c'}
    count=0
    with tarfile.open(ROOT/'remote/source.tar.gz','r:gz') as archive:
        members={m.name:m for m in archive if m.isfile()}
        for path,digest in result['input_hashes'].items():
            if path.startswith(PREFIX+'source/'):
                name=path.removeprefix(PREFIX+'source/')
                with archive.extractfile(members[name]) as stream:
                    assert sha_bytes(stream.read())==digest,name
                count+=1
            elif path.startswith(PREFIX):
                assert sha(ROOT/'remote'/path.removeprefix(PREFIX))==digest,path
        assert count==len(members)==133
    for path in (ROOT/'logs/final-tests.log',ROOT/'remote/logs/10244.tests.log'):
        text=path.read_text();assert 'Ran 106 tests' in text and 'OK (skipped=1)' in text
    assert 'FAILED (failures=2)' in (ROOT/'logs/before.log').read_text()
    assert 'ModuleNotFoundError' in (ROOT/'failed-10241/10241.tests.log').read_text()
    print(json.dumps(dict(status='verified',archive_files=len(expected),bound_source_files=count,
        test_cases=106,passed=105,skipped=1,peak_rss_kib=result['peak_rss_kib'],
        scope='Retained evidence only; no bulk-input rehash or biological accuracy claim.'),sort_keys=True))
if __name__=='__main__':verify()
