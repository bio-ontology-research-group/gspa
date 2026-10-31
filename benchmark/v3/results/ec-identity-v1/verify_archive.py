"""Reconstruct every entry and disposition without importing ec_identity."""
from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path
import re
ROOT=Path(__file__).resolve().parent
EC=r'\d+\.\d+\.\d+\.(?:n)?\d+'
def sha_bytes(data):return hashlib.sha256(data).hexdigest()
def read(path):return json.loads(path.read_text())
def verify():
    expected=read(ROOT/'archive-hashes.json')
    actual={str(p.relative_to(ROOT)):sha_bytes(p.read_bytes()) for p in ROOT.rglob('*')
            if p.is_file() and p.name!='archive-hashes.json'}
    assert actual==expected,'archive inventory/hash mismatch'
    source=gzip.decompress((ROOT/'inputs/enzyme.dat.gz').read_bytes())
    acquisition=read(ROOT/'inputs/source.json');audit=read(ROOT/'audit.json')
    assert sha_bytes((ROOT/'source/ec_identity.py').read_bytes())==audit['generator_sha256']
    assert len(source)==acquisition['bytes'] and sha_bytes(source)==acquisition['sha256']
    assert sha_bytes(source) in audit['inputs'].values()
    assert sha_bytes((ROOT/'inputs/legacy-projection.json').read_bytes()) in audit['inputs'].values()
    descriptions={}
    for block in source.decode().split('//\n'):
        ids=re.findall(r'^ID   (.+)$',block,re.M)
        if not ids:continue
        assert len(ids)==1 and ids[0] not in descriptions
        descriptions[ids[0]]=' '.join(re.findall(r'^DE   (.+)$',block,re.M))
    entries={}
    for ec,de in descriptions.items():
        if de.startswith('Transferred entry:'):
            targets=re.findall(EC,de);status='transferred'
        elif de=='Deleted entry.':targets=[];status='deleted'
        else:targets=[];status='active'
        entries[ec]=dict(description=de,status=status,targets=targets)
    assert entries==audit['entries'] and len(entries)==8472
    def resolve(ec):
        path=[]
        while True:
            assert ec not in path,'cycle'
            path.append(ec)
            if ec not in entries:return dict(status='unknown',canonical_id=None,path=path)
            r=entries[ec]
            if r['status']=='active':
                return dict(status='resolved',canonical_id=ec,path=path,preliminary=bool(re.search(r'\.n\d+$',ec)))
            if r['status']=='deleted':return dict(status='deleted',canonical_id=None,path=path)
            if len(r['targets'])>1:return dict(status='ambiguous',canonical_id=None,path=path,alternatives=r['targets'])
            assert len(r['targets'])==1
            ec=r['targets'][0]
    resolutions={ec:resolve(ec) for ec in entries}
    assert resolutions==audit['resolutions']
    assert Counter(r['status'] for r in resolutions.values())=={'resolved':8186,'ambiguous':102,'deleted':184}
    aliases=read(ROOT/'inputs/legacy-projection.json')['aliases']
    rows=audit['legacy_audit']['rows']
    assert len(rows)==len(aliases)==1309 and len({r['source_id'] for r in rows})==1309
    for row in rows:
        old=row['source_id'];target=aliases[old];a=resolve(old);b=resolve(target)
        assert row['legacy_target']==target and row['resolution']==a and row['legacy_target_resolution']==b
        if a['status']!='resolved':status=a['status']
        elif b['status']!='resolved':status='legacy_target_unresolved'
        elif a['canonical_id']!=b['canonical_id']:status='different_active_function'
        else:status='unique_active_equivalence'
        assert row['disposition']==status
    counts=dict(Counter(r['disposition'] for r in rows))
    assert counts==audit['legacy_audit']['counts']=={'unique_active_equivalence':1202,'ambiguous':101,'deleted':6}
    assert entries['1.1.1.5']['targets']==['1.1.1.303','1.1.1.304']
    assert '(R)-acetoin' in descriptions['1.1.1.303']
    assert '(S)-acetoin' in descriptions['1.1.1.304']
    assert resolve('1.1.1.5')['canonical_id'] is None
    assert resolve('1.6.5.3')['canonical_id']=='7.1.1.2'
    assert resolve('1.5.98.3')['canonical_id']=='1.5.98.3'
    assert resolve('1.6.2.1')['status']=='deleted'
    print(json.dumps(dict(status='verified',archive_files=len(expected),source_entries=len(entries),
                         legacy_replacements=len(rows),dispositions=counts,labels_read=False),sort_keys=True))
if __name__=='__main__':verify()
