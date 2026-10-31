#!/usr/bin/env python3
"""Verify retained dbCAN route evidence without executing models or tools."""
import csv
import hashlib
import json
from pathlib import Path
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parent


def totals(folder):
    suites = [ET.parse(p).getroot() for p in folder.rglob('TEST-*.xml')]
    assert suites, folder
    return {key: sum(int(s.get(key, 0)) for s in suites)
            for key in ('tests', 'failures', 'errors', 'skipped')}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify():
    summary = json.loads((ROOT / 'validation.json').read_text())
    for row in summary['test_results']:
        assert totals(ROOT / row['path']) == row['totals'], row
    replay = ROOT / summary['accepted_replay']
    report = json.loads((replay / 'report.json').read_text())
    assert report['status'] == 'controlled_route_parity_verified'
    for version in ('4', '5'):
        native = replay / ('annotate-' + version + '.jsonl')
        imported = replay / ('import-' + version + '.jsonl')
        assert native.read_bytes() == imported.read_bytes()
        rows = [json.loads(line) for line in native.read_text().splitlines()]
        assert len(rows) == 4
        assert all(r['metadata']['score_available'] is False and r['calibrated_prob'] is None for r in rows)
        assert all(isinstance(json.loads(r['metadata']['rawFieldsJson']), dict) for r in rows)
        assert all(isinstance(json.loads(r['metadata']['supportingMethodsJson']), list) for r in rows)
        assert 'Loaded 4 claims' in (replay / ('integrate-' + version + '.log')).read_text()
        assert len((replay / ('integrated-' + version + '.tsv')).read_text().splitlines()) == 1
        receipt = json.loads((replay / ('import-' + version + '.jsonl.import.json')).read_text())
        assert receipt['overview_sha256'] == sha(replay / ('native-' + version + '.tsv'))
        assert receipt['proteins_sha256'] == sha(replay / 'q.faa')
        assert receipt['output_sha256'] == sha(imported)
        assert receipt['annotation_count'] == 4 and receipt['protein_count'] == 3
    assert 'ERROR: no claims produced' in (replay / 'empty-builtin-default.log').read_text()
    assert (replay / 'empty-builtin.jsonl').read_bytes() == b''
    commands = json.loads((replay / 'commands.json').read_text())
    assert len(commands) == report['commands'] == 10
    assert next(c for c in commands if c['name'] == 'empty-builtin-default')['exit_code'] != 0
    assert next(c for c in commands if c['name'] == 'empty-builtin-allowed')['exit_code'] == 0
    merged = replay / 'nextflow-output/q/q_claims.jsonl'
    assert merged.read_bytes() == (replay / 'import-5.jsonl').read_bytes()
    assert report['nextflow_merged_claims_sha256'] == sha(merged)
    with (replay / 'trace.tsv').open() as stream:
        trace = list(csv.DictReader(stream, delimiter='\t'))
    for process in ('DBCAN', 'DBCAN_CLAIMS', 'MERGE_CLAIMS', 'INTEGRATE'):
        matches = [r for r in trace if r['name'].split('(')[0].strip() == process]
        assert len(matches) == 1 and matches[0]['status'] == 'COMPLETED', (process, matches)
    with (replay / 'missing-trace.tsv').open() as stream:
        trace = list(csv.DictReader(stream, delimiter='\t'))
    failures = [r for r in trace if r['status'] == 'FAILED']
    assert failures and all(r['name'].split('(')[0].strip() == 'DBCAN' for r in failures), failures
    assert not any(r['name'].split('(')[0].strip() in ('DBCAN_CLAIMS', 'MERGE_CLAIMS', 'INTEGRATE')
                   and r['status'] == 'COMPLETED' for r in trace)
    print('Retained test counts, import receipts, raw route parity and dbCAN-specific failure verified.')


if __name__ == '__main__':
    verify()
