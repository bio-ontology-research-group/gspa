"""Bind completed prerequisites before any full-comparison inference."""
import argparse
from datetime import datetime, timezone
import hashlib,json
from pathlib import Path

PROBE=Path('/data/hohndor/gspa-v3/multisource-comparison-probe-v4/jobs/10403/results/verification.json')
CLEAN=Path('/data/hohndor/gspa-v3/clean-cohort-v3/jobs/10368/audit.json')
JAR=Path('/data/hohndor/gspa-v3/ec-prefix-identity-v3/jobs/10402/gspa-3.0.0.jar')

def digest(p):
    with p.open('rb') as f:
        h=hashlib.sha256()
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
        return h.hexdigest()

def seal(root,probe=PROBE,clean=CLEAN):
    static=root/'static.sha256'
    records={}
    for line in static.read_text().splitlines():
        h,p=line.split('  ',1)
        if p in records or digest(Path(p))!=h: raise ValueError('changed or duplicate static input: '+p)
        records[p]=h
    if records.get(str(JAR))!='37aeee2b565951dd3a4e049a23aacd2e8155f29bd8b6216fed7006e4273965e6':
        raise ValueError('wrong tested JAR')
    exits=[probe.parent.parent/'exit-code.txt',clean.parent/'exit-code.txt']
    for p in exits:
        if p.read_text().strip()!='0': raise ValueError('prerequisite job did not exit successfully')
    p=json.loads(probe.read_text());c=json.loads(clean.read_text())
    expected=dict(shared_source_candidates=511,verified_both_source_candidates=511,
                  combined_candidates=29512,raw_claim_files_unchanged=True)
    if p.get('schema')!='gspa.source-context-probe/2' or p.get('status')!='verified' or p.get('ec_identity')!=expected:
        raise ValueError('corrected complete EC identity probe required')
    if p['inputs'].get(str(JAR))!=records[str(JAR)]: raise ValueError('probe JAR differs')
    if c.get('schema')!='gspa.clean-cohort-audit/2' or (c.get('genomes'),c.get('proteins'),c.get('native_proteins'),c.get('unsupported_proteins'))!=(18,49463,49462,1):
        raise ValueError('complete original CLEAN population required')
    for item in [static,probe,clean,*exits]: records[str(item)]=digest(item)
    # Exclusive creation: an interrupted/previous sealing attempt is retained.
    with (root/'sealed-prerequisites.json').open('x') as f:
        json.dump(dict(sealed_at=datetime.now(timezone.utc).isoformat(),inputs=records,
            scope='Frozen inputs before fitting; full runner still verifies all content and populations.'),f,indent=2)
        f.write('\n')
    records[str(root/'sealed-prerequisites.json')]=digest(root/'sealed-prerequisites.json')
    with (root/'packages.sha256').open('x') as f:
        f.write(''.join(h+'  '+p+'\n' for p,h in sorted(records.items())))
    return records

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('root',type=Path)
    seal(parser.parse_args().root)
