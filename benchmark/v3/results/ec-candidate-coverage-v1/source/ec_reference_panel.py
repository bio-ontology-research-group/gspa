"""Full positive-reference EC evaluation, including functions never proposed.

This module never selects thresholds or generates predictions. Nonmatching calls
are unconfirmed, not experimentally established negatives. Complex membership
annotations do not establish isolated-subunit catalytic activity.
"""
import csv
import math
import re
from collections import Counter, defaultdict

POLICY = {
    'schema': 'gspa.full-ec-reference-panel/1',
    'evidence': 'ECO:0000269',
    'population': 'Every positive exact EC assertion on every reference protein in the registered test genomes, including absent candidates.',
    'thresholds': 'Use the frozen mixed GO/EC calibration cutoffs and 0.5 without EC-specific retuning.',
    'ranking': 'Descending probability then lexical EC; report top-score ties and top-1/3/5 over all emitted candidates.',
    'interpretation': 'Positive-reference support only; nonmatches are unconfirmed, and missing reference coverage is not biological absence.',
}


def probability(value):
    if type(value) not in (int, float) or not math.isfinite(value) or not 0 <= value <= 1:
        raise ValueError('invalid probability')
    return value


def load_references(labels, splits, genomes):
    """Bind population to explicit genome splits, never to emitted candidates."""
    assignments = {}
    with splits.open() as f:
        for row in csv.DictReader(f, delimiter='\t'):
            gid, split = row['genome_id'], row['split']
            if not gid or gid in assignments or split not in ('train', 'calibration', 'test'):
                raise ValueError('invalid or duplicate genome split')
            assignments[gid] = split
    if set(assignments) != set(genomes):
        raise ValueError('reference split/genome population differs')
    test = {g for g, s in assignments.items() if s == 'test'}
    if not test:
        raise ValueError('no registered test genomes')
    refs = defaultdict(set); labels_seen = defaultdict(set)
    with labels.open() as f:
        for row in csv.DictReader(f, delimiter='\t'):
            if row['genome_id'] not in test or row['evidence'] != POLICY['evidence'] or not row['function_id'].startswith('EC:'):
                continue
            term = row['function_id'][3:]; pid = row['protein_id']
            if not pid or not re.fullmatch(r'[0-9]+\.[0-9]+\.[0-9]+\.[0-9]+', term):
                raise ValueError('invalid exact EC reference identity')
            label = row.get('label', '1')
            if label not in ('0', '1'):
                raise ValueError('invalid reference label')
            key = row['genome_id'], pid
            labels_seen[key, term].add(label)
            if label == '1':
                refs[key].add(term)
    if any(len(v) > 1 for v in labels_seen.values()):
        raise ValueError('contradictory experimental EC reference')
    return dict(refs), test


def evaluate(references, by_genome, test_genomes, cut):
    """Production keys are (protein_id, annotation_type, unprefixed EC id)."""
    probability(cut)
    if not set(test_genomes) <= by_genome.keys() or any(g not in test_genomes for g, _ in references):
        raise ValueError('missing test-genome predictions or foreign reference')
    records = []
    for (gid, protein), truth in sorted(references.items()):
        candidates = {ec: probability(score) for (pid, kind, ec), score in by_genome[gid].items()
                      if pid == protein and kind == 'EC'}
        ranked = sorted(candidates, key=lambda ec: (-candidates[ec], ec))
        calls = {ec for ec, p in candidates.items() if p >= cut}
        top = [ec for ec in ranked if candidates[ec] == candidates[ranked[0]]] if ranked else []
        records.append(dict(genome_id=gid, protein_id=protein, reference_ecs=sorted(truth),
            candidates=len(ranked), candidate_reference_matches=sorted(truth & candidates.keys()),
            calls=sorted(calls), matched_calls=sorted(calls & truth), unconfirmed_calls=sorted(calls - truth),
            missed_reference=sorted(truth - calls), top_tie=top, ranked_top5=ranked[:5],
            top1_reference_match=bool(truth & set(ranked[:1])), top3_reference_match=bool(truth & set(ranked[:3])),
            top5_reference_match=bool(truth & set(ranked[:5]))))
    tp = sum(len(r['matched_calls']) for r in records); unconfirmed = sum(len(r['unconfirmed_calls']) for r in records)
    positive = sum(len(t) for t in references.values())
    by_genome_count = Counter(g for (g, _), terms in references.items() for _ in terms)
    return dict(policy=POLICY, test_genomes=sorted(test_genomes),
        reference_pairs_by_genome={g: by_genome_count[g] for g in sorted(test_genomes)},
        proteins=len(records), reference_pairs=positive,
        emitted_reference_pairs=sum(len(r['candidate_reference_matches']) for r in records),
        proteins_with_candidates=sum(r['candidates'] > 0 for r in records), reference_matches=tp,
        unconfirmed_calls=unconfirmed, missed_reference=positive - tp,
        reference_precision=tp / (tp + unconfirmed) if tp + unconfirmed else None,
        reference_recall=tp / positive if positive else None,
        reference_micro_f1=2 * tp / (tp + unconfirmed + positive) if tp + unconfirmed + positive else None,
        top1_matches=sum(r['top1_reference_match'] for r in records), top3_matches=sum(r['top3_reference_match'] for r in records),
        top5_matches=sum(r['top5_reference_match'] for r in records), top1_tied=sum(len(r['top_tie']) > 1 for r in records),
        threshold=cut, records=records)
