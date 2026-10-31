#!/usr/bin/env python3
"""Verify retained evidence and independently recompute report arithmetic.

No model execution. Raw production TSV reconstruction and the 87 GB file audit
ran in job 10223; this portable check verifies that receipt, not remote bytes.
"""
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def read(name):
    return json.loads((ROOT / name).read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def same(a, b):
    if a is None or b is None:
        assert a is b, (a, b)
    else:
        assert math.isclose(a, b, rel_tol=1e-9, abs_tol=1e-11), (a, b)


def metrics(rows, arm, cut):
    cells = Counter(); squared = []; bins = defaultdict(list)
    for r in rows:
        p, y = r['scores'][arm], r['label']
        assert type(y) is int and y in (0, 1)
        assert type(p) in (int, float) and math.isfinite(p) and 0 <= p <= 1
        cells[(y, p >= cut)] += 1
        squared.append((p-y)**2); bins[min(9, int(10*p))].append(p-y)
    tp, fp, fn, tn = (cells[k] for k in [(1, True), (0, True), (1, False), (0, False)])
    return dict(n=len(rows), positives=tp+fn, tp=tp, fp=fp, fn=fn, tn=tn,
        precision=tp/(tp+fp) if tp+fp else None, recall=tp/(tp+fn) if tp+fn else None,
        f1=2*tp/(2*tp+fp+fn) if 2*tp+fp+fn else None,
        brier=math.fsum(squared)/len(rows),
        ece_10bin=math.fsum(abs(math.fsum(v)) for v in bins.values())/len(rows), threshold=cut)


def main():
    sums = dict(line.split('  ', 1)[::-1] for line in (ROOT/'SHA256SUMS').read_text().splitlines())
    actual = {str(p.relative_to(ROOT)) for p in ROOT.rglob('*') if p.is_file() and p.name != 'SHA256SUMS'}
    assert actual == sums.keys(), 'incomplete archive manifest'
    for name, digest in sums.items():
        assert sha(ROOT/name) == digest, name
    report = read('job-10223/audit/report.json')
    ledger = [json.loads(l) for l in (ROOT/'job-10223/audit/file-checks.jsonl').read_text().splitlines()]
    checks = {r['path']: r for r in ledger}
    assert len(checks) == len(ledger) == report['source_files'] == 6981
    assert all(r['agrees'] and r['actual'] == r['expected'] for r in ledger)
    assert sum(r['bytes'] for r in ledger) == report['source_bytes'] == 87052664741
    assert report['status'] == 'verified' and report['consumed_files_rechecked'] == 337
    assert report['root_verification_sha256'] == sha(ROOT/'analysis/verification.json')
    bindings = read('retained-source-bindings.json')
    for name, b in bindings.items():
        assert sha(ROOT/name) == b['sha256'] == checks[b['original_path']]['expected']
    parent = read('references/parent-registration.json')['inputs']
    for name, digest in read('original-job/supplement.json')['inputs'].items():
        assert parent[name] == digest == checks[name]['expected']
    for job in ('original-job', 'job-10223'):
        assert (ROOT/job/'exit-code.txt').read_text().strip() == '0'
    stages = read('logs/workflow.json')['stages']
    assert len(stages) == 7 and all(s['status'] == 'completed' and s['returncode'] == 0 for s in stages)
    assert sha(ROOT/'job-10223/audit-v1.py') == sha(ROOT/'audit-source/audit_context_comparison.py')

    rows = [json.loads(l) for l in (ROOT/'analysis/scored-supervision.jsonl').read_text().splitlines()]
    assert len(rows) == 6706
    old = read('analysis/results.json')['results']; thresholds = read('analysis/thresholds.json')
    with (ROOT/'references/splits.tsv').open() as f:
        test_genomes = {r['genome_id'] for r in csv.DictReader(f, delimiter='\t') if r['split'] == 'test'}
    truth = defaultdict(set)
    with (ROOT/'references/labels.tsv').open() as f:
        for r in csv.DictReader(f, delimiter='\t'):
            if r['genome_id'] in test_genomes and r['evidence'] == 'ECO:0000269' and r['label'] == '1' and r['function_id'].startswith('EC:'):
                truth[r['genome_id'], r['protein_id']].add(r['function_id'][3:])
    assert len(test_genomes) == 4 and len(truth) == 37 and sum(map(len, truth.values())) == 39
    metrics_checked = panels_checked = 0
    for regime, result in report['regimes'].items():
        subset = [r for r in rows if r['transfer'] == regime]
        assert len(subset) == result['supervision_rows']
        cal = [r for r in subset if r['split'] == 'calibration']
        test = [r for r in subset if r['split'] == 'test']
        assert len(test) == result['test_rows']
        for arm, t in thresholds[regime].items():
            best, cut = max((metrics(cal, arm, c)['f1'], c) for c in {r['scores'][arm] for r in cal})
            same(cut, t['value']); same(best, t['calibration_f1'])
        for policy, data in result['policies'].items():
            for arm, original in old[regime]['policies'][policy]['arms'].items():
                cut = thresholds[regime][arm]['value'] if policy == 'calibration_selected' else .5
                for key, value in metrics(test, arm, cut).items():
                    same(value, original[key])
                metrics_checked += 1
                for kind, by_arm in data['supervised_by_type'].items():
                    typed = [r for r in test if r['example_id'].split('|')[2] == kind]
                    for key, value in metrics(typed, arm, cut).items(): same(value, by_arm[arm][key])
                panel = data['full_ec_reference'][arm]; counts = Counter(); seen = set()
                for r in panel['records']:
                    key = r['genome_id'], r['protein_id']; assert key not in seen; seen.add(key)
                    refs, calls = set(r['reference_ecs']), set(r['calls'])
                    assert refs == truth[key] and refs & calls == set(r['matched_calls'])
                    assert calls - refs == set(r['unconfirmed_calls']) and refs - calls == set(r['missed_reference'])
                    assert set(r['candidate_reference_matches']) <= refs
                    assert len(r['ranked_top5']) == min(5, r['candidates'])
                    counts['reference_matches'] += len(refs & calls)
                    counts['unconfirmed_calls'] += len(calls - refs)
                    counts['missed_reference'] += len(refs - calls)
                    counts['emitted_reference_pairs'] += len(r['candidate_reference_matches'])
                    counts['proteins_with_candidates'] += r['candidates'] > 0
                    counts['top1_tied'] += len(r['top_tie']) > 1
                    for k in (1, 3, 5):
                        hit = bool(refs & set(r['ranked_top5'][:k]))
                        assert hit == r[f'top{k}_reference_match']
                        counts[f'top{k}_matches'] += hit
                assert seen == truth.keys()
                for key, value in counts.items(): same(value, panel[key])
                tp, un = counts['reference_matches'], counts['unconfirmed_calls']
                for key, value in dict(proteins=37, reference_pairs=39, threshold=cut,
                    reference_precision=tp/(tp+un) if tp+un else None,
                    reference_recall=tp/39, reference_micro_f1=2*tp/(tp+un+39)).items(): same(value, panel[key])
                panels_checked += 1
    assert metrics_checked == panels_checked == 36
    audits = list((ROOT/'metadata/observed-gap-context-model-replay-v1').glob('*/*/*/audit.json'))
    assert len(audits) == 12
    assert sum(len(json.loads(p.read_text())) for p in audits) == 216
    assert sum(v['candidates'] for p in audits for v in json.loads(p.read_text()).values()) == 4539408
    parity = list((ROOT/'metadata').glob('*parity*/verification.json'))
    assert len(parity) == 2
    assert sum(v['rows'] for p in parity for v in json.loads(p.read_text())['results'].values()) == 40236
    print(json.dumps(dict(status='verified', files=len(sums), retained_original_bindings=len(bindings),
        original_summaries=metrics_checked, full_ec_panels=panels_checked, reference_proteins=37,
        reference_assertions=39, remote_hash_receipt_files=6981), sort_keys=True))


if __name__ == '__main__':
    main()
