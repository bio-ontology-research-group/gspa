"""Price positive-reference metrics; no fitting, inference or threshold selection.

Unlabelled functions are unconfirmed, not experimentally negative. The exact-ID
analysis is primary. A second analysis merges only unambiguous ENZYME transfers;
unresolved IDs remain in both the candidate and reference populations.
"""
from collections import Counter, defaultdict
import math
import random

import ec_identity as E

POLICY = {
    'primary': 'exact',
    'sensitivity': 'current_unambiguous',
    'projection': 'Single-destination transfers only; retain unresolved original IDs; deduplicate references and take maximum candidate score.',
    'ranking': 'All emitted four-level EC IDs, including preliminary IDs. Primary top-1 gives the fraction of top-tied IDs that match the reference; report worst/best tie bounds. Lexical ranks/top-k/MRR are explicitly supplementary.',
    'population': 'All declared reference proteins, including proteins with no candidates. Other proteins have no accuracy labels.',
    'calls': 'Integrated scores >= 0.5 are a fixed diagnostic, not a calibrated operating point. No thresholded raw-score claims.',
    'interpretation': 'Reference agreement, not enzyme detection or experimentally established correctness of every nonmatching call.',
}


def ec(value):
    if not isinstance(value, str):
        raise ValueError('EC must be a string')
    result = value.removeprefix('EC:')
    if not E.ID.fullmatch(result):
        raise ValueError(f'invalid four-level EC: {value!r}')
    return result


def score(value):
    if type(value) not in (int, float) or not math.isfinite(value) or not 0 <= value <= 1:
        raise ValueError('invalid score')
    return value


def project(identifier, entries, mode):
    identifier = ec(identifier)
    result = E.resolve(entries, identifier)
    if mode not in ('exact', 'current_unambiguous'):
        raise ValueError('unknown EC identity policy')
    return (result['canonical_id'] if mode == 'current_unambiguous' and
            result['status'] == 'resolved' else identifier), result['status']


def protein(reference, candidates, entries, mode, *, integrated):
    """One complete reference protein; missing candidates never remove its truth."""
    if not reference:
        raise ValueError('empty protein reference')
    reference = {ec(x) for x in reference}
    truth = {project(x, entries, mode)[0] for x in reference}
    values, origins = {}, defaultdict(list)
    for original, value in candidates.items():
        original = ec(original)
        target, _ = project(original, entries, mode)
        values[target] = max(values.get(target, -1), score(value))
        origins[target].append(original)
    ranked = sorted(values, key=lambda x: (-values[x], x))
    ties = [x for x in ranked if values[x] == values[ranked[0]]] if ranked else []
    ranks = {x: i+1 for i, x in enumerate(ranked)}
    matching_ranks = [ranks[x] for x in truth if x in ranks]
    if matching_ranks:
        first_score = max(values[x] for x in truth if x in values)
        higher = sum(v > first_score for v in values.values())
        block = {x for x, v in values.items() if v == first_score}
        rr_low, rr_high = 1/(higher+len(block-truth)+1), 1/(higher+1)
    else:
        rr_low = rr_high = 0.0
    result = dict(reference_ecs=sorted(truth), original_reference_ecs=sorted(reference),
        reference_resolution=dict(Counter(E.resolve(entries, x)['status'] for x in reference)),
        candidate_resolution=dict(Counter(E.resolve(entries, x)['status'] for x in candidates)),
        candidates=len(values), original_candidates=len(candidates),
        merged_candidates={x: sorted(v) for x, v in origins.items() if len(v) > 1},
        candidate_reference_matches=sorted(truth & values.keys()),
        absent_reference_candidates=sorted(truth - values.keys()),
        ranked_top5=[dict(ec=x, score=values[x]) for x in ranked[:5]], top_tie=ties,
        top1_match=bool(truth & set(ranked[:1])), top3_match=bool(truth & set(ranked[:3])),
        top5_match=bool(truth & set(ranked[:5])),
        top_tie_any_match=bool(truth & set(ties)),
        top_tie_all_match=bool(ties) and set(ties) <= truth,
        top1_credit=len(truth & set(ties))/len(ties) if ties else 0.0,
        reciprocal_rank_low=rr_low, reciprocal_rank_high=rr_high,
        reciprocal_rank=1/min(matching_ranks) if matching_ranks else 0.0)
    if integrated:
        calls = {x for x, value in values.items() if value >= 0.5}
        result['fixed_0_5'] = dict(matched=sorted(calls & truth),
            unconfirmed=sorted(calls - truth), missed=sorted(truth - calls))
    return result


def summary(records):
    if not records:
        raise ValueError('empty reference population')
    n = len(records)
    result = dict(proteins=n, reference_pairs=sum(len(r['reference_ecs']) for r in records),
        original_reference_pairs=sum(len(r['original_reference_ecs']) for r in records),
        candidate_reference_matches=sum(len(r['candidate_reference_matches']) for r in records),
        proteins_without_candidates=sum(r['candidates'] == 0 for r in records),
        tied_top_proteins=sum(len(r['top_tie']) > 1 for r in records),
        lexical_mean_reciprocal_rank=sum(r['reciprocal_rank'] for r in records)/n,
        mean_reciprocal_rank_low=sum(r['reciprocal_rank_low'] for r in records)/n,
        mean_reciprocal_rank_high=sum(r['reciprocal_rank_high'] for r in records)/n)
    for key in ('top1_match', 'top3_match', 'top5_match', 'top_tie_any_match', 'top_tie_all_match'):
        result[key] = sum(r[key] for r in records)
    result['lexical_top1_reference_accuracy'] = result['top1_match']/n
    result['top1_reference_accuracy'] = sum(r['top1_credit'] for r in records)/n
    if any('fixed_0_5' in r for r in records):
        if not all('fixed_0_5' in r for r in records):
            raise ValueError('cannot mix raw and integrated threshold records')
        counts = {k: sum(len(r['fixed_0_5'][k]) for r in records)
                  for k in ('matched', 'unconfirmed', 'missed')}
        tp, unconfirmed, missed = (counts[k] for k in ('matched', 'unconfirmed', 'missed'))
        result['fixed_0_5'] = dict(**counts,
            reference_precision=tp/(tp+unconfirmed) if tp+unconfirmed else None,
            reference_recall=tp/(tp+missed),
            reference_f1=2*tp/(2*tp+unconfirmed+missed))
    return result


def paired(before, after):
    """Reference outcomes, retaining ties and every missing-candidate protein."""
    if set(before) != set(after):
        raise ValueError('different paired protein populations')
    if any(before[k]['reference_ecs'] != after[k]['reference_ecs'] for k in before):
        raise ValueError('different paired reference functions')
    return dict(proteins=len(before),
        top1_credit_increase=sum(max(0,after[k]['top1_credit']-before[k]['top1_credit']) for k in before),
        top1_credit_decrease=sum(max(0,before[k]['top1_credit']-after[k]['top1_credit']) for k in before),
        top1_gained_reference_match=sum(not before[k]['top1_match'] and after[k]['top1_match'] for k in before),
        top1_lost_reference_match=sum(before[k]['top1_match'] and not after[k]['top1_match'] for k in before),
        changed_top1=sum([x['ec'] for x in before[k]['ranked_top5'][:1]] !=
                         [x['ec'] for x in after[k]['ranked_top5'][:1]] for k in before))


def bootstrap_top1(by_genome, *, baseline=None, repetitions=2000, seed=13):
    """Paired genome resampling; retain all genomes, even with zero references."""
    if not by_genome or repetitions < 1 or (baseline is not None and set(baseline) != set(by_genome)):
        raise ValueError('invalid bootstrap population')
    if baseline is not None and any(len(baseline[g]) != len(rows) for g, rows in by_genome.items()):
        raise ValueError('unmatched bootstrap protein populations')
    groups = sorted(by_genome)
    counts = {g: (sum(by_genome[g]) - (sum(baseline[g]) if baseline is not None else 0),
                  len(by_genome[g])) for g in groups}
    total = sum(n for _, n in counts.values())
    if total == 0:
        raise ValueError('bootstrap has no reference proteins')
    rng = random.Random(seed); values = []; empty = 0
    for _ in range(repetitions):
        sample = [counts[rng.choice(groups)] for _ in groups]
        denominator = sum(n for _, n in sample)
        if denominator:
            values.append(sum(n for n, _ in sample)/denominator)
        else:
            empty += 1
    values.sort()
    return dict(estimate=sum(n for n, _ in counts.values())/total, n_genomes=len(groups),
        repetitions=repetitions, seed=seed, empty_resamples=empty,
        ci_low=values[int(0.025*(len(values)-1))] if values else None,
        ci_high=values[int(0.975*(len(values)-1))] if values else None)
