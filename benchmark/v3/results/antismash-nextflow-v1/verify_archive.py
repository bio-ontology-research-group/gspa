#!/usr/bin/env python3
"""Verify archived hashes, source bindings, test gates and route observations (no compute)."""
import hashlib
import gzip
import re
import json
from pathlib import Path
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parent


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def main():
    entries = (ROOT / 'archive.sha256').read_text().splitlines()
    for line in entries:
        expected, relative = line.split('  ', 1)
        assert sha(ROOT / relative) == expected, relative
    manifest = json.loads((ROOT / 'source-files-v2.json').read_text())
    sources = [p for p in (ROOT / 'source').rglob('*') if p.is_file()]
    for source in sources:
        relative = str(source.relative_to(ROOT/'source'))
        assert sha(source) == manifest[relative], relative
    with gzip.open(ROOT/'GCA_003598085.1.json.gz', 'rb') as stream:
        native_bytes = stream.read()
    native_hash = hashlib.sha256(native_bytes).hexdigest()
    assert native_hash == '423d21fab26bfcfacb0dafd4a21f10a7172c659ee0dc3bd36c2b4cfcbc2ec9a6'
    data = json.loads(native_bytes)
    expected = set()
    for record in data['records']:
        cid = record.get('original_id') or record['id']
        for pi, cds in enumerate((f for f in record['features'] if f['type']=='CDS'), 1):
            match = re.fullmatch(r'\[(\d+):(\d+)\]\(([+-])\)', cds['location'])
            assert match
            start,end = int(match[1])+1,int(match[2])
            for region in (f for f in record['features'] if f['type']=='region'):
                rm = re.fullmatch(r'\[(\d+):(\d+)\]\(([+-])\)', region['location'])
                assert rm
                rs,re_ = int(rm[1])+1,int(rm[2])
                if start <= re_ and end >= rs:
                    expected.add((f'{cid}_{pi}',cid+':region'+region['qualifiers']['region_number'][0],
                                  rs,re_,'contained' if start>=rs and end<=re_ else 'partial'))
    assert len(expected) == 70
    overlays = json.loads((ROOT/'route-overlay-hashes.json').read_text())
    for name,digest in overlays.items():
        assert sha(ROOT/Path(name).name) == digest
    successes = []
    for job in (ROOT/'jobs').iterdir():
        code = (job/'exit-code.txt').read_text().strip()
        if code != '0':
            continue
        suites = [ET.parse(p).getroot() for p in (ROOT/'jobs/10481/final-tests').rglob('TEST-*.xml')]
        assert (job/'parent-gate.log').is_file()
        assert sum(int(s.get('tests', 0)) for s in suites) == 1521
        assert sum(int(s.get('skipped', 0)) for s in suites) == 6
        assert suites and not any(int(s.get('failures', 0)) or int(s.get('errors', 0)) for s in suites)
        report = json.loads((job/'cli/report.json').read_text())
        assert report['status'] == 'observed_replay_route_parity_verified'
        assert report['memberships'] == {'q':70,'zero':0}
        assert report['protein_counts'] == {'q':4915,'zero':6}
        assert report['functional_predictions'] == 0
        assert report['jar_sha256'] == '069263f7b876222815928a637d02106b9c3b31a756edd1cf692432acc2e52303'
        assert report['checker_sha256'] == overlays['benchmark/v3/check_antismash_nextflow.py']
        statuses = {r['name']:r['exit_code'] for r in report['commands']}
        assert set(statuses) == {'before-import','import','annotate','nextflow','integrate','nextflow-missing'}
        assert statuses['before-import'] != 0 and statuses['nextflow-missing'] != 0
        assert all(statuses[name] == 0 for name in ('import','annotate','nextflow','integrate'))
        cli = job/'cli'
        raw = (cli/'import.jsonl').read_bytes()
        rows = [json.loads(line) for line in raw.splitlines()]
        assert len(rows) == 70
        assert {(r['protein_id'],r['metadata']['region_id'],r['metadata']['native_start']+1,
                 r['metadata']['native_end'],r['metadata']['cdsOverlap']) for r in rows} == expected
        assert all(r['metadata']['score_available'] is False and r['calibrated_prob'] is None
                   and r['metadata']['source_sha256'] == native_hash for r in rows)
        assert raw == (cli/'annotation.jsonl').read_bytes()
        assert raw == (cli/'nextflow-output/q/q_claims.jsonl').read_bytes()
        assert (cli/'direct.tsv').read_bytes() == (cli/'nextflow-output/q/q_integrated.tsv').read_bytes()
        assert len((cli/'direct.tsv').read_text().splitlines()) == 1
        assert not (cli/'nextflow-output/zero/zero_claims.jsonl').read_bytes()
        for sample in ('q','zero'):
            receipt = json.loads((cli/f'nextflow-output/{sample}/antismash/{sample}_claims_antismash.jsonl.import.json').read_text())
            assert receipt['coverage'] == 'complete'
            assert receipt['region_count'] == (3 if sample == 'q' else 0)
            assert receipt['output_sha256'] == sha(cli/f'nextflow-output/{sample}/{sample}_claims.jsonl')
        successes.append(dict(job=job.name, cases=sum(int(s.get('tests',0)) for s in suites),
                              skips=sum(int(s.get('skipped',0)) for s in suites)))
    assert successes, 'no successful full validation archived'
    print(json.dumps(dict(files=len(entries), scoped_sources=len(sources), successes=successes),indent=2))


if __name__ == '__main__':
    main()
