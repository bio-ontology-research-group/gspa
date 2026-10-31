"""Independent source-feature membership audit; does not import the adapter/checker.

Uses native region SeqFeature locations and qualifiers, not records.areas.
This is a format/identity check, not independently measured biological truth.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re


def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1048576),b''):h.update(block)
    return h.hexdigest()


def check(native, cli):
    data=json.loads(native.read_text())
    expected={};regions=[];protein_ids=set()
    location=lambda text:re.fullmatch(r'\[(\d+):(\d+)\]\(([+-])\)',text)
    for record in data['records']:
        cid=record.get('original_id') or record['id']
        local=[]
        for feature in record['features']:
            if feature['type']!='region':continue
            m=location(feature['location']);assert m
            start,end=int(m[1]),int(m[2]);q=feature['qualifiers']
            number=int(q['region_number'][0]);rid=f'{cid}:region{number}'
            products=q['product'];local.append((start,end,rid,products));regions.append((cid,start+1,end))
        for feature in record['features']:
            if feature['type']!='CDS':continue
            m=location(feature['location']);assert m
            start,end=int(m[1]),int(m[2]);pid='GSPA_'+feature['qualifiers']['locus_tag'][0]
            assert pid not in protein_ids;protein_ids.add(pid)
            for left,right,rid,products in local:
                if max(left,start)<min(right,end):
                    expected[pid,rid]=dict(products=products,overlap='contained' if left<=start and end<=right else 'partial')
    assert len(data['records'])==145 and len(protein_ids)==4915
    assert len(regions)==3 and len(expected)==70
    for route in ('annotate','annotate_integrated'):
        rows=[json.loads(line) for line in (cli/route/'raw.jsonl').read_text().splitlines()]
        observed={}
        for row in rows:
            meta=row['metadata'];key=(row['protein_id'],meta['region_id'])
            assert key not in observed
            assert row['function_type']=='SECONDARY_METABOLITE' and meta['score_available'] is False
            assert row.get('calibrated_prob') is None and meta['source_sha256']==sha(native)
            observed[key]=dict(products=meta['products'],overlap=meta['cdsOverlap'])
        assert observed==expected,route
    before=[json.loads(line) for line in (cli/'before/raw.jsonl').read_text().splitlines()]
    assert not before
    report=json.loads((cli/'verification.json').read_text())
    assert report['status']=='verified' and report['input_proteins']==4915 and report['raw_memberships']==70
    return dict(status='verified',native_sha256=sha(native),cli_report_sha256=sha(cli/'verification.json'),
        source_records=145,input_proteins=4915,source_feature_regions=3,expected_memberships=70,
        checked_routes=['annotate','annotate_integrated'],
        scope='Independent native SeqFeature versus output membership/identity check; not experimental function or BGC accuracy.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--native',type=Path,required=True);p.add_argument('--cli',type=Path,required=True)
    a=p.parse_args();print(json.dumps(check(a.native,a.cli),indent=2))
