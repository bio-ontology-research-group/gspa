#!/usr/bin/env python3
"""Check retained byte inventory, test counts and controlled route score semantics."""
import csv
import hashlib
import json
import math
from pathlib import Path
import xml.etree.ElementTree as ET

ROOT=Path(__file__).resolve().parent

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def read(path):return json.loads(path.read_text())
def jsonl(path):return [json.loads(line) for line in path.read_text().splitlines()]


def check_typed(folder,tag,input_fasta,kind):
    receipt=read(folder/(tag+'.clean.receipt.json'))
    tsv=folder/(tag+'.clean.tsv');native=folder/(tag+'.clean.native.csv')
    assert receipt['input_fasta_sha256']==sha(input_fasta)
    assert receipt['output_sha256']==sha(tsv) and receipt['native_output_sha256']==sha(native)
    assert receipt['min_score']==0 and receipt['filtered_calls']==0
    ids=receipt['native_to_input_ids'];calls={};seen=set()
    with native.open() as stream:
        for row in csv.reader(stream):
            assert row[0] not in seen;seen.add(row[0])
            for token in row[1:]:
                term,value=token.split('/')
                key=(ids[row[0]],term)
                assert key not in calls;calls[key]=(float(value),value,row[0])
    assert seen==set(ids) and len(ids)==2
    emitted=set()
    for line in tsv.read_text().splitlines()[1:]:
        pid,term,score,atype,metadata=line.split('\t');metadata=json.loads(metadata)
        key=(pid,term);value,text,native_id=calls[key]
        assert key not in emitted;emitted.add(key)
        assert atype=='EC' and metadata['score_kind']==kind
        assert metadata['native_score']==value and metadata['native_score_text']==text
        assert metadata['native_protein_id']==native_id
        assert metadata['score_available']==(kind=='gmm_confidence')
        assert metadata['native_score_direction']==('higher' if kind=='gmm_confidence' else 'lower')
        assert float(score)==(value if kind=='gmm_confidence' else 0)
    assert emitted==set(calls) and len(calls)==receipt['emitted_calls']==4


def main():
    listed={}
    for line in (ROOT/'SHA256SUMS').read_text().splitlines():
        digest,name=line.split('  ',1);listed[name]=digest
        assert sha(ROOT/name)==digest,name
    actual={str(p.relative_to(ROOT)) for p in ROOT.rglob('*') if p.is_file() and p.name!='SHA256SUMS' and '__pycache__' not in p.parts}
    assert set(listed)==actual,'archive inventory changed'
    report=read(ROOT/'report.json')
    for gate in report['jvm_gates']:
        folder=ROOT/'jobs'/str(gate['job'])/'final-tests'/gate['module']
        suites=[ET.parse(p).getroot() for p in folder.rglob('TEST-*.xml')]
        counts={key:sum(int(r.attrib.get(key,0)) for r in suites) for key in ['tests','failures','errors','skipped']}
        assert counts==gate['counts'] and counts['failures']==counts['errors']==0
    route=ROOT/'jobs'/str(report['route_job'])/'route'
    assert read(route/'report.json')['status']=='controlled_route_parity_verified'
    commands=read(route/'commands.json')
    assert len(commands)==11
    for kind in ('gmm_confidence','distance'):
        capture=route/('capture-'+kind)
        check_typed(capture,'query',capture/'input.faa',kind)
        assert (route/('annotate-'+kind+'.jsonl')).read_bytes()==(route/('import-'+kind+'.jsonl')).read_bytes()
        nf=route/('nf-output-'+kind)/'q'
        check_typed(nf/'clean','q',nf/'gene_calling/q.faa',kind)
        assert (nf/'q_claims.jsonl').read_bytes()==(nf/'clean/q_claims_clean.jsonl').read_bytes()
        for path in (route/('import-'+kind+'.jsonl'),nf/'q_claims.jsonl'):
            rows=jsonl(path)
            assert len(rows)==4 and all(r['metadata']['score_available']==(kind=='gmm_confidence') for r in rows)
            assert {r['metadata']['native_score'] for r in rows}==({1e-8,.8} if kind=='gmm_confidence' else {.01,8.0627})
        if kind=='distance':
            assert len((nf/'q_integrated.tsv').read_text().splitlines())==1
            assert not (nf/'ensemble/q.ensemble.tsv').exists()
    for fault in ('missing','unknown','malformed'):
        status=read(route/('annotation-'+fault)/'q_predictor_status.json')
        assert status['summary']['FAILED']==1 and status['summary']['SUCCESS_ZERO']==0
        assert (route/(fault+'.jsonl')).read_bytes()==b''
    assert next(c for c in commands if c['name']=='nextflow-missing')['exit_code']!=0
    dbcan=ROOT/'jobs'/str(report['route_job'])/'dbcan-route'
    assert read(dbcan/'report.json')['status']=='controlled_route_parity_verified'
    assert len(read(dbcan/'commands.json'))==10
    assert (dbcan/'nextflow-output/q/q_claims.jsonl').read_bytes()==(dbcan/'import-5.jsonl').read_bytes()
    if report.get('native_audit'):
        from audit_native import audit
        expected=read(ROOT/'native-audit.json')
        native_job=ROOT/'jobs'/str(report['native_job'])
        assert audit(native_job)==expected
        imported=ROOT/'jobs'/str(report['native_import_job'])
        assert (imported/'exit-code.txt').read_text().strip()=='0'
        parity=read(imported/'import-parity.json')
        for dataset in ['new','price']:
            predictions={}
            for line in (native_job/'predictions'/dataset/(dataset+'.clean.tsv')).read_text().splitlines()[1:]:
                pid,term,score,atype,metadata=line.split('\t')
                predictions[pid,term]=(float(score),json.loads(metadata))
            rows=jsonl(imported/(dataset+'.claims.jsonl'));seen=set()
            assert len(rows)==len(predictions)==parity['datasets'][dataset]['claims']
            for row in rows:
                key=(row['protein_id'],row['function_id']);assert key not in seen;seen.add(key)
                score,metadata=predictions[key]
                assert row['raw_score']==score and row['metadata']==metadata
                assert row['source']=='clean' and row['function_type']=='EC'
                assert math.isfinite(row['calibrated_prob']) and 0<=row['calibrated_prob']<=1
            assert seen==set(predictions)
            assert [min(r['calibrated_prob'] for r in rows),max(r['calibrated_prob'] for r in rows)]==parity['datasets'][dataset]['derived_probability_range']
    print(json.dumps({'status':'verified','files':len(listed),'controlled_invocations':21,
        'native_audit':bool(report.get('native_audit'))},sort_keys=True))

if __name__=='__main__':main()
