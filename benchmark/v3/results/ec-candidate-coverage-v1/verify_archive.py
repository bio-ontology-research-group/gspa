#!/usr/bin/env python3
"""Portable arithmetic/provenance check; no models or production evaluator imported."""
import csv
import gzip
import hashlib
import io
import json
from collections import Counter, defaultdict
from pathlib import Path
import tarfile

ROOT = Path(__file__).resolve().parent


def read(name):
    return json.loads((ROOT/name).read_text())


def digest(data):
    return hashlib.sha256(data).hexdigest()


def main():
    hashes = dict(line.split('  ', 1)[::-1] for line in (ROOT/'SHA256SUMS').read_text().splitlines())
    files = {str(p.relative_to(ROOT)) for p in ROOT.rglob('*') if p.is_file() and p.name != 'SHA256SUMS'}
    assert files == hashes.keys()
    for name, expected in hashes.items(): assert digest((ROOT/name).read_bytes()) == expected, name
    trace = read('trace-report.json'); assert trace['status'] == 'verified'
    inputs = {}
    for source, expected in trace['inputs'].items():
        if source.endswith('.claims.jsonl'): continue
        path = ROOT/'trace-inputs'/Path(source).name
        data = path.read_bytes() if path.exists() else gzip.decompress(Path(str(path)+'.gz').read_bytes())
        assert digest(data) == expected, source
        inputs[path.name] = data.decode()
    labels = list(csv.DictReader(io.StringIO(inputs['labels.tsv']), delimiter='\t'))
    truth = {(r['genome_id'], r['protein_id'], r['function_id']) for r in labels
        if r['evidence'] == 'ECO:0000269' and r['label'] == '1' and r['function_id'].startswith('EC:')}
    assert len(truth) == trace['reference_pairs'] == 335
    splits = {r['genome_id']: r['split'] for r in csv.DictReader(io.StringIO(inputs['splits.tsv']), delimiter='\t')}
    seqs = {(r['genome_id'], r['genome_protein_id']): r['protein_id'] for r in csv.DictReader(io.StringIO(inputs['proteins.tsv']), delimiter='\t')}
    refs = defaultdict(set)
    for r in csv.DictReader(io.StringIO(inputs['reference_labels.tsv']), delimiter='\t'): refs[r['function_id']].add(r['accession'])
    hits = defaultdict(list)
    for line in inputs['hits.tsv'].splitlines():
        q, s, identity, coverage, evalue, bitscore = line.split('\t')
        hits[q].append(dict(subject=s, identity=float(identity), coverage=float(coverage), evalue=float(evalue), bitscore=float(bitscore)))
    protocol = read('trace-inputs/context-protocol.json'); claims = {}; seen = set()
    for line in (ROOT/'trace-inputs/reference-locus-claims.jsonl').read_text().splitlines():
        r = json.loads(line); c = r['claim']; key = r['mode'], r['genome_id'], c['protein_id'], 'EC:'+c['function_id'].removeprefix('EC:')
        assert key not in claims; claims[key] = c['raw_score']
        source = protocol['claims'][r['mode']][r['genome_id']]
        assert source['path'] == r['source_path'] and source['sha256'] == trace['inputs'][r['source_path']]
    assert sum(v['records'] for files in protocol['claims'].values() for v in files.values()) == trace['claims_scanned'] == 756568
    counts = defaultdict(Counter)
    for line in (ROOT/'trace-rows.jsonl').read_text().splitlines():
        r = json.loads(line); key = r['mode'], r['genome_id'], r['protein_id'], r['function_id']
        assert key not in seen; seen.add(key)
        assert key[1:] in truth and seqs[key[1:3]] == r['sequence_id'] and splits[key[1]] == r['split']
        assert sorted(refs[key[3]]) == r['reference_accessions']
        all_hits = hits[r['sequence_id']]; support = [h for h in all_hits if h['subject'] in refs[key[3]]]
        eligible = [h for h in support if key[0] == 'homology' or not(h['identity'] >= 30 and h['coverage'] >= 80)]
        assert len(all_hits) == r['query_hit_count'] and support == r['matching_hits'] and eligible == r['eligible_hits']
        score = max([round(h['identity']/100 * min(h['coverage'],100)/100, 6) for h in eligible]+[0.0])
        assert (key in claims) == (score > 0)
        assert r['claim_score'] == claims.get(key) == (score if score > 0 else None)
        status = 'emitted' if score > 0 else 'no_reference_annotation' if not refs[key[3]] else 'no_matching_search_hit' if not support else 'family_excluded' if not eligible else 'zero_score'
        assert status == r['status']; counts[key[0], r['split']][status] += 1
    assert len(seen) == trace['rows'] == 670
    for (mode, split), count in counts.items(): assert dict(count) == trace['summary'][mode][split]
    assert {k[1:] for k in seen if k[0] == 'homology'} == truth
    assert {k[1:] for k in seen if k[0] == 'family_excluded'} == truth

    job = ROOT/'job-10224'; assert (job/'exit-code.txt').read_text().strip() == '0'
    old = read('job-10224/inputs/results.json')['results']; new = read('job-10224/results/results.json')['results']
    audit = read('job-10224/inputs/audit-report.json')['regimes']
    assert read('job-10224/results/thresholds.json') == read('job-10224/inputs/thresholds.json')
    assert read('job-10224/results/protocol.json')['schema'] == 'gspa.observed-gap-context-comparison/2'
    total = 0
    for transfer, result in new.items():
        for policy, summary in result['policies'].items():
            panels = summary.pop('full_ec_reference')
            assert summary == old[transfer]['policies'][policy]
            for arm, panel in panels.items():
                expected = audit[transfer]['policies'][policy]['full_ec_reference'][arm]
                assert all(panel[k] == v for k, v in expected.items())
                assert panel['reference_pairs_by_genome'] == {'GCA_000007065.1':21, 'GCA_000008765.1':13, 'GCA_000009985.1':0, 'GCA_000147055.1':5}
                for r in panel['records']:
                    emitted = {ec.removeprefix('EC:') for (m,g,p,ec) in claims if (m,g,p) == (transfer,r['genome_id'],r['protein_id'])}
                    assert r['candidates'] == len(emitted)
                    assert set(r['candidate_reference_matches']) == emitted & set(r['reference_ecs'])
                total += 1
    assert total == 36
    with tarfile.open(ROOT/'source-v1.tar.gz') as tar:
        for name in ['analyze_observed_gap_context_models.py','ec_reference_panel.py']:
            frozen = tar.extractfile('home/leechuck/Public/software/gspa/benchmark/v3/'+name).read()
            assert frozen == (ROOT/'source'/name).read_bytes(), 'production code changed after real replay'
    print(json.dumps(dict(status='verified', files=len(files), traced_assertions=335, trace_rows=670,
        complete_ec_panels=36, unchanged_thresholds=18, unchanged_primary_summaries=36, unchanged_paired_statistics=48)))


if __name__ == '__main__':
    main()
