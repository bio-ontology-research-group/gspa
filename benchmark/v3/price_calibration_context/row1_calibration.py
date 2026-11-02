#!/usr/bin/env python3
"""Row 1 of the acceptance bar on the Price benchmark: genus-blocked calibration/test split,
calibration-only call thresholds (the WP5 rule: >= 10 calls at >= 90 % support), calibration
metrics, corrections/harms and genome-bootstrap intervals, from the retained per-protein score
files of the existing release-candidate collections. No inference is rerun.

Frozen plan: tasks/v3-row1-row5-price-calibration-context-2026-11-01.md (Part A).

Unconfirmed calls are never negatives. Where a statistic needs a binary outcome (threshold
selection, top-call Brier/ECE) an unconfirmed call scores 0 and the statistic is labelled a
reference-support statistic, exactly as unknown decisions count against precision in the WP5 rule.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import random
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
V3 = HERE.parent
sys.path.insert(0, str(V3))
import ec_identity as E  # noqa: E402
import price_reference_metrics as M  # noqa: E402

SCHEMA = 'gspa.price-row1-calibration/1'
PRECISION, MIN_DECISIONS = 0.90, 10           # the WP5 rule, unchanged
SELECTION_RULE = 'max_supported_then_min_unsupported_then_highest_threshold; unconfirmed calls count against support'
SPLIT_SALT = 'gspa-price-row1-2026-11-01|'
MERGED_BLOCKS = {'ncbi_genus_2901164': 'ncbi_genus_286'}   # Stutzerimonas -> Pseudomonas block
MIN_ARM_PROTEINS, MIN_ARM_GENOMES = 40, 5
BINS = 10
BOOT_N, BOOT_SEED = 2000, 13                  # the Price evaluator's own bootstrap settings
PREDICTORS = ('clean', 'dpz', 'enzbert', 'graphec', 'mapred', 'topec')
ARMS = ('raw_export', 'stage0', 'stage1_taxon', 'stage2b')
CONTRASTS = (('raw_export', 'stage0'), ('stage0', 'stage1_taxon'), ('stage0', 'stage2b'))


def sha(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def save(path, obj):
    path = Path(path)
    with path.open('x') as f:
        json.dump(obj, f, indent=2, sort_keys=True, allow_nan=False)
        f.write('\n')
    return sha(path)


def load(path):
    return json.loads(Path(path).read_text())


# ----------------------------------------------------------------------------- data

def read_retained(directory, arm_name, bound):
    """Return {predictor: {genome: {pid: (reference set, {ec: score})}}} for one arm of one evaluation directory."""
    out = defaultdict(dict)
    files = sorted(Path(directory).glob('*.json'))
    if not files:
        raise ValueError(f'no retained files in {directory}')
    for path in files:
        rec = load(path)
        bound[str(path)] = sha(path)
        cell = rec['cell']
        if arm_name not in rec['scores']:
            raise ValueError(f'{path}: arm {arm_name} missing (has {sorted(rec["scores"])})')
        genome = cell['assembly']
        proteins = {}
        for pid, refs in rec['reference'].items():
            scores = rec['scores'][arm_name].get(pid)
            if scores is None:
                raise ValueError(f'{path}: no score vector for reference protein {pid}')
            proteins[pid] = ({M.ec(x) for x in refs}, {M.ec(k): M.score(v) for k, v in scores.items()})
        if genome in out[cell['predictor']]:
            raise ValueError(f'duplicate genome {genome} for {cell["predictor"]}')
        out[cell['predictor']][genome] = proteins
    return out


def gather(args, bound):
    """All arms for all predictors: data[arm][predictor][genome][pid] = (refs, scores)."""
    data = {}
    data['raw_export'] = read_retained(args.rc1 / 'retained', 'raw_export', bound)
    data['stage0'] = read_retained(args.rc1 / 'retained', 'stage0', bound)
    stage1 = read_retained(args.rc1 / 'retained', 'stage1', bound)
    if stage1 != data['stage0']:
        raise ValueError('stage1 differs from stage0 in the rc1 retained files; the plan assumes identity')
    taxon = defaultdict(dict)
    for p in PREDICTORS:
        taxon.update(read_retained(args.taxon / f'report-{p}' / 'retained', 'stage1_taxon', bound))
    data['stage1_taxon'] = taxon
    s2b = defaultdict(dict)
    for p in ('clean', 'dpz', 'graphec'):
        s2b.update(read_retained(args.stage2b_a / f'report-{p}' / 'retained', 'stage2b', bound))
    for p in ('enzbert', 'mapred', 'topec'):
        s2b.update(read_retained(args.stage2b_b / f'report-{p}' / 'retained', 'stage2b', bound))
    data['stage2b'] = s2b
    # one population everywhere
    reference = None
    for arm, by_pred in data.items():
        if set(by_pred) != set(PREDICTORS):
            raise ValueError(f'{arm}: predictors {sorted(by_pred)}')
        for p, by_genome in by_pred.items():
            pop = {(g, pid): frozenset(v[0]) for g, prot in by_genome.items() for pid, v in prot.items()}
            if reference is None:
                reference = pop
            elif pop != reference:
                raise ValueError(f'{arm}/{p}: reference population differs')
    if len(reference) != 148 or sum(len(v) for v in reference.values()) != 151 or len({g for g, _ in reference}) != 22:
        raise ValueError('reference population is not 148 proteins / 151 pairs / 22 genomes')
    return data, reference


# ----------------------------------------------------------------------------- split

def genus_split(metadata, reference, order_rule='hash'):
    """order_rule 'hash': the pre-registered A2 rule. 'size_then_hash': the dated deviation (blocks walked from the
    largest reference-protein count down, ties by the same hash), applied only after the A2 rule failed its minimum."""
    genomes = sorted({g for g, _ in reference})
    counts = Counter(g for g, _ in reference)
    blocks = defaultdict(list)
    for g in genomes:
        clade = metadata['genomes'][g]['clade']
        blocks[MERGED_BLOCKS.get(clade, clade)].append(g)
    h = lambda b: hashlib.sha256((SPLIT_SALT + b).encode()).hexdigest()
    if order_rule == 'hash':
        order = sorted(blocks, key=h)
    elif order_rule == 'size_then_hash':
        order = sorted(blocks, key=lambda b: (-sum(counts[g] for g in blocks[b]), h(b)))
    else:
        raise ValueError(order_rule)
    arms = {'calibration': [], 'test': []}
    totals = {'calibration': 0, 'test': 0}
    for b in order:
        arm = 'calibration' if totals['calibration'] <= totals['test'] else 'test'
        arms[arm].append(b)
        totals[arm] += sum(counts[g] for g in blocks[b])
    assignment = {g: arm for arm, bs in arms.items() for b in bs for g in blocks[b]}
    summary = {arm: dict(blocks=bs, genomes=sorted(g for b in bs for g in blocks[b]),
                         reference_proteins=totals[arm]) for arm, bs in arms.items()}
    ok = all(s['reference_proteins'] >= MIN_ARM_PROTEINS and len(s['genomes']) >= MIN_ARM_GENOMES for s in summary.values())
    return dict(schema=SCHEMA + ':split', rule=dict(block='NCBI genus (cohort metadata clade), Stutzerimonas merged into Pseudomonas',
        order=('sha256(' + repr(SPLIT_SALT) + ' + block) ascending' if order_rule == 'hash' else
               'reference-protein count descending, ties by sha256(' + repr(SPLIT_SALT) + ' + block)'),
        order_rule=order_rule, assignment='greedy to the arm with fewer reference proteins, calibration on ties',
        minimum=dict(reference_proteins=MIN_ARM_PROTEINS, genomes=MIN_ARM_GENOMES)),
        block_order=order, blocks={b: sorted(gs) for b, gs in blocks.items()},
        reference_proteins_per_genome=dict(sorted(counts.items())), arms=summary, assignment=assignment,
        requirement_met=ok)


# ----------------------------------------------------------------------------- calls and thresholds

def pairs(proteins, entries):
    """All candidate pairs of reference proteins: list of (pid, ec, score, supported) and the truth set."""
    rows, truth = [], set()
    for pid, (refs, scores) in proteins.items():
        t = {M.project(x, entries, 'exact')[0] for x in refs}
        truth |= {(pid, x) for x in t}
        merged = {}
        for ec, s in scores.items():
            tgt, _ = M.project(ec, entries, 'exact')
            merged[tgt] = max(merged.get(tgt, -1.0), s)
        for ec, s in merged.items():
            rows.append((pid, ec, s, ec in t))
    return rows, truth


def choose_threshold(rows):
    """WP5 rule on the calibration calls: >= MIN_DECISIONS calls, support >= PRECISION,
    max (supported, -unsupported, threshold)."""
    ordered = sorted(rows, key=lambda r: -r[2])
    best = None
    most_supported = None   # diagnostic: the admissible-size threshold with the highest support share
    first = None            # diagnostic: the first threshold reaching MIN_DECISIONS calls
    sup = unsup = 0
    i = 0
    n = len(ordered)
    while i < n:
        t = ordered[i][2]
        j = i
        while j < n and ordered[j][2] == t:
            sup += ordered[j][3]
            unsup += not ordered[j][3]
            j += 1
        calls = sup + unsup
        if calls >= MIN_DECISIONS:
            rec = dict(threshold=t, calls=calls, supported=sup, unsupported=unsup, support=sup / calls)
            if first is None:
                first = rec
            if most_supported is None or (rec['support'], t) > (most_supported['support'], most_supported['threshold']):
                most_supported = rec
            if sup / calls >= PRECISION:
                key = (sup, -unsup, t)
                if best is None or key > best[0]:
                    best = (key, rec)
        i = j
    diagnostic = dict(calibration_pairs=n, supported_pairs=sum(r[3] for r in rows),
                      first_threshold_with_min_calls=first, highest_support_with_min_calls=most_supported)
    if best is None:
        return dict(threshold=None, status='no_supported_threshold_abstain', calibration_calls={}, diagnostic=diagnostic)
    return dict(threshold=best[1]['threshold'], status='selected_on_calibration', calibration_calls=best[1], diagnostic=diagnostic)


def call_counts(rows, truth, t):
    """Per-genome-agnostic counts at threshold t."""
    called = {(pid, ec) for pid, ec, s, _ in rows if s >= t}
    supported = called & truth
    return dict(calls=len(called), supported=len(supported), unconfirmed=len(called - truth),
                reference_pairs=len(truth), supported_reference_pairs=len(supported),
                missed=len(truth - called))


# ----------------------------------------------------------------------------- bootstrap

def genome_bootstrap(vectors, stat, n=BOOT_N, seed=BOOT_SEED):
    """Paired genome resampling; `vectors` maps genome -> list of additive counts; stat(summed) -> value or None."""
    groups = sorted(vectors)
    rng = random.Random(seed)
    values = []
    empty = 0
    for _ in range(n):
        total = [0.0] * len(next(iter(vectors.values())))
        for _ in groups:
            v = vectors[rng.choice(groups)]
            for k, x in enumerate(v):
                total[k] += x
        value = stat(total)
        if value is None:
            empty += 1
        else:
            values.append(value)
    values.sort()
    point = stat([sum(vectors[g][k] for g in groups) for k in range(len(total))])
    return dict(estimate=point, n_genomes=len(groups), repetitions=n, seed=seed, empty_resamples=empty,
                ci_low=values[int(0.025 * (len(values) - 1))] if values else None,
                ci_high=values[int(0.975 * (len(values) - 1))] if values else None)


def ratio(i, j):
    return lambda v: (v[i] / v[j]) if v[j] else None


# ----------------------------------------------------------------------------- calibration metrics

def bin_of(s):
    return min(int(s * BINS), BINS - 1)


def calibration_metrics(by_genome, entries):
    """by_genome: {genome: proteins}. Returns positive-pair Brier, top-call Brier/ECE, reliability tables."""
    pos_sq, pos_n = 0.0, 0
    top = []   # (score, credit)
    pair_bins = [[0, 0, 0.0] for _ in range(BINS)]   # pairs, reference pairs, score sum
    for proteins in by_genome.values():
        rows, truth = pairs(proteins, entries)
        best = {}
        for pid, ec, s, supported in rows:
            b = bin_of(s)
            pair_bins[b][0] += 1
            pair_bins[b][1] += supported
            pair_bins[b][2] += s
        for pid, (refs, scores) in proteins.items():
            rec = M.protein(refs, scores, entries, 'exact', integrated=True)
            for ec in rec['reference_ecs']:
                s = max([v for k, v in scores.items() if M.project(k, entries, 'exact')[0] == ec], default=0.0)
                pos_sq += (1.0 - s) ** 2
                pos_n += 1
            if rec['ranked_top5']:
                top.append((rec['ranked_top5'][0]['score'], rec['top1_credit']))
            else:
                top.append((0.0, 0.0))
    top_bins = [[0, 0.0, 0.0] for _ in range(BINS)]   # n, score sum, credit sum
    for s, y in top:
        b = bin_of(s)
        top_bins[b][0] += 1
        top_bins[b][1] += s
        top_bins[b][2] += y
    n_top = len(top)
    ece = sum(abs(b[1] / b[0] - b[2] / b[0]) * b[0] / n_top for b in top_bins if b[0]) if n_top else None
    brier_top = sum((s - y) ** 2 for s, y in top) / n_top if n_top else None
    return dict(
        convention='positive-pair Brier is label-safe; top-call Brier/ECE treat an unconfirmed top call as 0 (reference support, not error rate)',
        positive_pair_brier=pos_sq / pos_n if pos_n else None, reference_pairs=pos_n,
        top_call_brier=brier_top, top_call_ece=ece, bins=BINS, top_calls=n_top,
        top_call_reliability=[dict(bin=f'[{i/BINS:.1f},{(i+1)/BINS:.1f}{")" if i < BINS-1 else "]"}', n=b[0],
            mean_score=b[1] / b[0] if b[0] else None, supported_credit=b[2],
            unconfirmed_credit=b[0] - b[2]) for i, b in enumerate(top_bins)],
        pair_reliability=[dict(bin=f'[{i/BINS:.1f},{(i+1)/BINS:.1f}{")" if i < BINS-1 else "]"}', pairs=b[0],
            reference_pairs=b[1], mean_score=b[2] / b[0] if b[0] else None,
            reference_share=b[1] / b[0] if b[0] else None) for i, b in enumerate(pair_bins)])


# ----------------------------------------------------------------------------- per-arm evaluation

def arm_vectors(by_genome, entries, t):
    """Per-genome additive vectors [calls, supported, unconfirmed, reference pairs, supported reference pairs, top1 credit, proteins]."""
    vectors, records = {}, {}
    for g, proteins in by_genome.items():
        rows, truth = pairs(proteins, entries)
        c = call_counts(rows, truth, t) if t is not None else dict(calls=0, supported=0, unconfirmed=0,
                                                                    reference_pairs=len(truth), supported_reference_pairs=0)
        recs = {pid: M.protein(refs, scores, entries, 'exact', integrated=True) for pid, (refs, scores) in proteins.items()}
        records[g] = recs
        vectors[g] = [c['calls'], c['supported'], c['unconfirmed'], c['reference_pairs'], c['supported_reference_pairs'],
                      sum(r['top1_credit'] for r in recs.values()), len(recs)]
    return vectors, records


def evaluate_arm(by_genome, entries, policy):
    t = policy['threshold']
    vectors, records = arm_vectors(by_genome, entries, t)
    total = [sum(v[k] for v in vectors.values()) for k in range(7)]
    result = dict(threshold=t, status=policy['status'], genomes=len(vectors), proteins=total[6],
                  top1=dict(credit=total[5], share=genome_bootstrap(vectors, ratio(5, 6))),
                  top1_match=sum(r['top1_match'] for recs in records.values() for r in recs.values()),
                  top3_match=sum(r['top3_match'] for recs in records.values() for r in recs.values()),
                  fixed_0_5=M.summary([r for recs in records.values() for r in recs.values()])['fixed_0_5'])
    if t is not None:
        result['calibrated_calls'] = dict(calls=total[0], supported=total[1], unconfirmed=total[2],
            reference_pairs=total[3], supported_reference_pairs=total[4],
            support_share=genome_bootstrap(vectors, ratio(1, 0)),
            reference_recall=genome_bootstrap(vectors, ratio(4, 3)))
    return result, records, vectors


def contrast(before, after, entries, t_before, t_after):
    """Reference pairs called only after (corrections) / only before (harms); top-1 gained/lost; genome bootstrap."""
    vectors = {}
    flat_b, flat_a = {}, {}
    for g in sorted(before):
        rb, truth = pairs(before[g], entries)
        ra, truth_a = pairs(after[g], entries)
        if truth != truth_a:
            raise ValueError('contrast truth differs')
        def called(rows, t):
            return {(p, e) for p, e, s, _ in rows if t is not None and s >= t} & truth
        cb, ca = called(rb, t_before), called(ra, t_after)
        cb5 = {(p, e) for p, e, s, _ in rb if s >= 0.5} & truth
        ca5 = {(p, e) for p, e, s, _ in ra if s >= 0.5} & truth
        recs_b = {pid: M.protein(r, s, entries, 'exact', integrated=True) for pid, (r, s) in before[g].items()}
        recs_a = {pid: M.protein(r, s, entries, 'exact', integrated=True) for pid, (r, s) in after[g].items()}
        for pid in recs_b:
            flat_b[(g, pid)] = recs_b[pid]
            flat_a[(g, pid)] = recs_a[pid]
        gained = sum(not recs_b[p]['top1_match'] and recs_a[p]['top1_match'] for p in recs_b)
        lost = sum(recs_b[p]['top1_match'] and not recs_a[p]['top1_match'] for p in recs_b)
        vectors[g] = [len(ca - cb), len(cb - ca), len(ca5 - cb5), len(cb5 - ca5), gained, lost]
    tot = [sum(v[k] for v in vectors.values()) for k in range(6)]
    diff = lambda i, j: (lambda v: v[i] - v[j])
    return dict(thresholds=dict(before=t_before, after=t_after,
                                rule='each arm at its own calibrated threshold; None means the arm abstained and no calibrated calls exist'),
        calibrated=dict(corrected=tot[0], harmed=tot[1], net=genome_bootstrap(vectors, diff(0, 1))),
        fixed_0_5=dict(corrected=tot[2], harmed=tot[3], net=genome_bootstrap(vectors, diff(2, 3))),
        top1=dict(gained=tot[4], lost=tot[5], net=genome_bootstrap(vectors, diff(4, 5))),
        paired=M.paired(flat_b, flat_a))


# ----------------------------------------------------------------------------- integrity (C6)

def published_top1(results_root, arm, predictor):
    folder = {'raw_export': 'price-reference-rc1', 'stage0': 'price-reference-rc1',
              'stage1_taxon': 'price-taxon-rc2', 'stage2b': 'price-stage2b-rc1'}[arm]
    rep = load(results_root / folder / predictor / 'report.json')
    return rep['summaries']['exact'][predictor]['arms'][arm]['top1_match']


# ----------------------------------------------------------------------------- main

def run(args):
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=False)
    ledger = {}
    bound = {}
    audit = load(args.ec_audit)
    entries = audit['entries']
    E.validate_transfer_graph(entries)
    bound[str(args.ec_audit)] = sha(args.ec_audit)
    bound[str(args.cohort_metadata)] = sha(args.cohort_metadata)
    for path in (Path(__file__), Path(M.__file__), Path(E.__file__)):
        bound[str(path)] = sha(path)
    data, reference = gather(args, bound)

    # C1: split, written and hashed before any threshold is chosen
    split = genus_split(load(args.cohort_metadata), reference, 'hash')
    if not split['requirement_met']:
        # The pre-registered rule failed its own minimum: keep its record and apply the dated deviation rule.
        ledger['split_a2_rule_failed.json'] = save(out / 'split_a2_rule_failed.json', split)
        split = genus_split(load(args.cohort_metadata), reference, 'size_then_hash')
        split['deviation'] = 'A2 hash-order rule gave fewer than 40 calibration proteins (see split_a2_rule_failed.json); size-then-hash order applied, plan deviation dated 2026-11-02'
    ledger['split.json'] = save(out / 'split.json', split)
    if not split['requirement_met']:
        raise SystemExit('split requirement not met under both rules; see split.json')
    cal = {g for g, a in split['assignment'].items() if a == 'calibration'}
    test = {g for g, a in split['assignment'].items() if a == 'test'}

    def subset(by_genome, genomes):
        return {g: p for g, p in by_genome.items() if g in genomes}

    # C2: thresholds on calibration genomes only, for the primary and the swapped split
    policies = {}
    for name, cal_set in (('primary', cal), ('swapped', test)):
        pol = {}
        for arm in ARMS:
            for p in PREDICTORS:
                rows = []
                for g, proteins in subset(data[arm][p], cal_set).items():
                    rows += pairs(proteins, entries)[0]
                pol[f'{arm}/{p}'] = choose_threshold(rows)
        policies[name] = pol
    policy = dict(schema=SCHEMA + ':policy', target_support=PRECISION, min_calls=MIN_DECISIONS,
                  selection_rule=SELECTION_RULE, split_sha256=ledger['split.json'],
                  calibration_genomes=dict(primary=sorted(cal), swapped=sorted(test)), arms=policies)
    ledger['policy.json'] = save(out / 'policy.json', policy)

    # C3-C5: test metrics with the frozen thresholds
    report = dict(schema=SCHEMA, plan='tasks/v3-row1-row5-price-calibration-context-2026-11-01.md',
                  split_sha256=ledger['split.json'], policy_sha256=ledger['policy.json'],
                  convention='unconfirmed calls are not negatives; see each table',
                  analyses={}, integrity={}, inputs=bound)
    for name, test_set, cal_set in (('primary', test, cal), ('swapped', cal, test)):
        analysis = dict(test_genomes=sorted(test_set), calibration_genomes=sorted(cal_set), arms={}, contrasts={},
                        calibration_set_metrics={})
        recs = {}
        for arm in ARMS:
            analysis['arms'][arm] = {}
            analysis['calibration_set_metrics'][arm] = {}
            recs[arm] = {}
            for p in PREDICTORS:
                pol = policies[name][f'{arm}/{p}']
                res, records, _ = evaluate_arm(subset(data[arm][p], test_set), entries, pol)
                res['calibration_metrics'] = calibration_metrics(subset(data[arm][p], test_set), entries)
                analysis['arms'][arm][p] = res
                recs[arm][p] = records
                cal_res, _, _ = evaluate_arm(subset(data[arm][p], cal_set), entries, pol)
                analysis['calibration_set_metrics'][arm][p] = dict(top1_match=cal_res['top1_match'],
                    calibrated_calls=cal_res.get('calibrated_calls'), fixed_0_5=cal_res['fixed_0_5'])
        for before, after in CONTRASTS:
            analysis['contrasts'][f'{before}_to_{after}'] = {
                p: contrast(subset(data[before][p], test_set), subset(data[after][p], test_set), entries,
                            policies[name][f'{before}/{p}']['threshold'], policies[name][f'{after}/{p}']['threshold'])
                for p in PREDICTORS}
        report['analyses'][name] = analysis

    # C6: recomputed top-1 over all 22 genomes equals the published summaries
    for arm in ARMS:
        for p in PREDICTORS:
            recomputed = sum(M.protein(r, s, entries, 'exact', integrated=arm != 'raw_export')['top1_match']
                             for proteins in data[arm][p].values() for r, s in proteins.values())
            published = published_top1(args.results_root, arm, p)
            report['integrity'][f'{arm}/{p}'] = dict(recomputed_top1=recomputed, published_top1=published,
                                                     equal=recomputed == published)
    report['integrity_all_equal'] = all(v['equal'] for v in report['integrity'].values())
    ledger['report.json'] = save(out / 'report.json', report)
    write_tables(report, out)
    save(out / 'ledger.json', dict(schema=SCHEMA + ':ledger', files=ledger))
    print(json.dumps(dict(split=ledger['split.json'], policy=ledger['policy.json'], report=ledger['report.json'],
                          integrity_all_equal=report['integrity_all_equal']), indent=1))


def fmt(x, nd=3):
    return '' if x is None else (f'{x:.{nd}f}' if isinstance(x, float) else str(x))


def ci(b, nd=3):
    return f"{fmt(b['estimate'], nd)} [{fmt(b['ci_low'], nd)}, {fmt(b['ci_high'], nd)}]"


def write_tables(report, out):
    for name, analysis in report['analyses'].items():
        with (out / f'thresholds_{name}.tsv').open('x') as f:
            w = csv.writer(f, delimiter='\t')
            w.writerow(['arm', 'predictor', 'status', 'threshold', 'cal_calls', 'cal_supported', 'cal_unsupported',
                        'test_calls', 'test_supported', 'test_unconfirmed', 'test_support_share_ci', 'test_reference_recall_ci',
                        'test_top1_match', 'test_top1_share_ci', 'test_top3_match', 'fixed05_matched', 'fixed05_unconfirmed'])
            for arm, by_p in analysis['arms'].items():
                for p, r in by_p.items():
                    cc = r.get('calibrated_calls')
                    calc = analysis['calibration_set_metrics'][arm][p]['calibrated_calls'] or {}
                    w.writerow([arm, p, r['status'], fmt(r['threshold'], 6), calc.get('calls', ''), calc.get('supported', ''),
                                calc.get('unconfirmed', ''),
                                cc['calls'] if cc else '', cc['supported'] if cc else '', cc['unconfirmed'] if cc else '',
                                ci(cc['support_share']) if cc else '', ci(cc['reference_recall']) if cc else '',
                                r['top1_match'], ci(r['top1']['share']), r['top3_match'],
                                r['fixed_0_5']['matched'], r['fixed_0_5']['unconfirmed']])
        with (out / f'calibration_{name}.tsv').open('x') as f:
            w = csv.writer(f, delimiter='\t')
            w.writerow(['arm', 'predictor', 'reference_pairs', 'positive_pair_brier', 'top_calls', 'top_call_brier', 'top_call_ece_10bins'])
            for arm, by_p in analysis['arms'].items():
                for p, r in by_p.items():
                    c = r['calibration_metrics']
                    w.writerow([arm, p, c['reference_pairs'], fmt(c['positive_pair_brier'], 4), c['top_calls'],
                                fmt(c['top_call_brier'], 4), fmt(c['top_call_ece'], 4)])
        with (out / f'reliability_{name}.tsv').open('x') as f:
            w = csv.writer(f, delimiter='\t')
            w.writerow(['arm', 'predictor', 'table', 'bin', 'n', 'mean_score', 'supported_or_reference', 'unconfirmed_or_other'])
            for arm, by_p in analysis['arms'].items():
                for p, r in by_p.items():
                    c = r['calibration_metrics']
                    for b in c['top_call_reliability']:
                        w.writerow([arm, p, 'top_call', b['bin'], b['n'], fmt(b['mean_score']), fmt(b['supported_credit']), fmt(b['unconfirmed_credit'])])
                    for b in c['pair_reliability']:
                        w.writerow([arm, p, 'all_pairs', b['bin'], b['pairs'], fmt(b['mean_score']), b['reference_pairs'], b['pairs'] - b['reference_pairs']])
        with (out / f'contrasts_{name}.tsv').open('x') as f:
            w = csv.writer(f, delimiter='\t')
            w.writerow(['contrast', 'predictor', 't_before', 't_after', 'cal_corrected', 'cal_harmed', 'cal_net_ci',
                        'fixed05_corrected', 'fixed05_harmed', 'fixed05_net_ci', 'top1_gained', 'top1_lost', 'top1_net_ci', 'changed_top1'])
            for cname, by_p in analysis['contrasts'].items():
                for p, c in by_p.items():
                    w.writerow([cname, p, fmt(c['thresholds']['before'], 6), fmt(c['thresholds']['after'], 6),
                                c['calibrated']['corrected'], c['calibrated']['harmed'], ci(c['calibrated']['net'], 2),
                                c['fixed_0_5']['corrected'], c['fixed_0_5']['harmed'], ci(c['fixed_0_5']['net'], 2),
                                c['top1']['gained'], c['top1']['lost'], ci(c['top1']['net'], 2), c['paired']['changed_top1']])


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--rc1', type=Path, required=True, help='rc1 evaluation directory (contains retained/)')
    p.add_argument('--taxon', type=Path, required=True, help='taxon-rc2 evaluation directory (report-<p>/retained)')
    p.add_argument('--stage2b-a', type=Path, required=True, help='stage2b evaluation dir with report-{clean,dpz,graphec}')
    p.add_argument('--stage2b-b', type=Path, required=True, help='stage2b evaluation dir with report-{enzbert,mapred,topec}')
    p.add_argument('--cohort-metadata', type=Path, required=True)
    p.add_argument('--ec-audit', type=Path, required=True)
    p.add_argument('--results-root', type=Path, required=True, help='benchmark/v3/results (published report.json files)')
    p.add_argument('--out', type=Path, required=True)
    run(p.parse_args(argv))


if __name__ == '__main__':
    main()
