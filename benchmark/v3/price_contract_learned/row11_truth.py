#!/usr/bin/env python3
"""Plan step G2: truth intersection of the frozen public case bundle (evaluate_stage3_agent.py prepare output).

Truth sources: Price-148 complete ECs (`<assembly>.ec.tsv`) and the Swiss-Prot experimental Price ECs (exact
sequence). A case target that is an EC is matched exactly (positive / unknown; the existing wrapper's negative rule
is reported separately as a sensitivity count, never used). A GO target is "positive" when ec2go maps one of the
candidate protein's labelled ECs to that exact GO term; no GO negative is derived. Labels from Swiss-Prot entries that
are themselves reference entries of the DIAMOND predictor (leakage report) are counted separately as circular.
The public bundle's SHA-256 is recorded before any label file is opened. Writes per-split truth files in the
`evaluate_stage3_agent.py` format (status, source) for the methods, with circular labels excluded.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

COMPLETE = re.compile(r'^(?:EC:)?(\d+\.\d+\.\d+\.\d+)$')


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def ec2go(path):
    m = defaultdict(set)
    for line in Path(path).read_text().splitlines():
        if line.startswith('EC:') and ';' in line:
            ec = line.split()[0].removeprefix('EC:')
            m[ec].add(line.rsplit(';', 1)[1].strip())
    return m


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--public', type=Path, required=True)
    ap.add_argument('--price-labels', type=Path, required=True)
    ap.add_argument('--swissprot-ecs', type=Path, required=True)
    ap.add_argument('--ec2go', type=Path, required=True)
    ap.add_argument('--leakage', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    a = ap.parse_args(argv)
    a.out.mkdir(parents=True, exist_ok=False)
    public_sha = sha(a.public)
    (a.out / 'public.sha256').write_text(f'{public_sha}  {a.public}\n')
    cases = json.loads(a.public.read_text())['cases']
    labels = defaultdict(set)          # (genome, protein) -> {(ec, source)}
    for p in sorted(a.price_labels.glob('*.ec.tsv')):
        g = p.name.split('.ec.tsv')[0]
        for r in csv.DictReader(p.open(), delimiter='\t'):
            m = COMPLETE.match(r['ec'])
            if m:
                labels[(g, r['protein_id'])].add((m.group(1), 'price'))
    circ = {(c['genome_id'], c['protein_id'], c['function_id']) for c in json.loads(a.leakage.read_text())['swissprot_price_labels_in_reference']}
    for r in csv.DictReader(a.swissprot_ecs.open(), delimiter='\t'):
        m = COMPLETE.match(r['function_id'])
        if m:
            src = 'swissprot_circular' if (r['genome_id'], r['protein_id'], m.group(1)) in circ else 'swissprot'
            labels[(r['genome_id'], r['protein_id'])].add((m.group(1), src))
    mapping = ec2go(a.ec2go)
    counts = Counter()
    by_split = defaultdict(Counter)
    truth_rows = defaultdict(list)
    detail = []
    target_kinds = Counter()
    for case in cases:
        an = case['anomaly']
        g, split, target = case['genome_id'], case['split'], an['target_function_id']
        kind = 'GO' if target.startswith('GO:') else 'EC'
        target_kinds[(an.get('anomaly_type') or an.get('type'), kind)] += 1
        for c in an['candidates']:
            lab = labels.get((g, c['protein_id']), set())
            status, source, neg_rule = 'unknown', '', False
            if lab:
                if kind == 'EC':
                    t = target.removeprefix('EC:')
                    hit = [s for e, s in lab if e == t]
                    neg_rule = not hit and COMPLETE.match(target) is not None
                else:
                    hit = [s for e, s in lab if target in mapping.get(e, ())]
                if hit:
                    status = 'positive'
                    source = 'swissprot_circular' if set(hit) == {'swissprot_circular'} else sorted(set(hit) - {'swissprot_circular'})[0]
            key = status if status == 'unknown' else f'positive/{source}'
            counts[key] += 1
            by_split[split][key] += 1
            if neg_rule:
                counts['negative_under_existing_rule_not_used'] += 1
            if lab:
                counts['candidate_pairs_on_labelled_proteins'] += 1
                by_split[split]['candidate_pairs_on_labelled_proteins'] += 1
                detail.append(dict(case=an['anomaly_id'], split=split, genome=g, protein=c['protein_id'], target=target,
                                   labels=sorted(lab), status=status, source=source))
            if status == 'positive' and source != 'swissprot_circular':
                truth_rows[split].append(dict(genome_id=g, protein_id=c['protein_id'], function_id=target, status='positive',
                                              source=source))
    for split in ('calibration', 'test'):
        (a.out / f'truth_{split}.json').write_text(json.dumps(truth_rows[split], indent=1) + '\n')
    rep = dict(schema='gspa.price-contract-row11-truth/1', public_sha256=public_sha, cases=len(cases),
               cases_by_split=dict(Counter(c['split'] for c in cases)),
               candidate_pairs=sum(len(c['anomaly']['candidates']) for c in cases),
               case_types={f'{t}/{k}': n for (t, k), n in sorted(target_kinds.items(), key=str)},
               pairs=dict(counts), pairs_by_split={k: dict(v) for k, v in by_split.items()},
               label_sources=dict(price_tables=len(list(a.price_labels.glob('*.ec.tsv'))), swissprot=sha(a.swissprot_ecs), ec2go=sha(a.ec2go)),
               labelled_pairs_detail=detail,
               methods_applicable=all(any(r for r in truth_rows[s]) for s in ('calibration', 'test')))
    (a.out / 'truth_report.json').write_text(json.dumps(rep, indent=1, sort_keys=True) + '\n')
    print(json.dumps({k: rep[k] for k in ('cases', 'cases_by_split', 'candidate_pairs', 'case_types', 'pairs', 'pairs_by_split', 'methods_applicable')}, indent=1))


if __name__ == '__main__':
    main()
