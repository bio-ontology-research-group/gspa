#!/usr/bin/env python3
"""WP5 step 5: join hidden truth and measured phenotypes, calibrate, evaluate, report.

Phases (design §6–§8 and supervisor rulings §12–§13)::

  join       --bundle BUNDLE --runs RUNS --out OUT
  calibrate  --bundle BUNDLE --runs RUNS --out OUT
  evaluate   --bundle BUNDLE --runs RUNS --out OUT --policy-sha256 SHA256
  report     --out OUT

Only this script opens private truth and measured labels. ``join`` reads calibration files only.
``evaluate`` verifies the exact frozen policy file, protocol, run files and calibration join before
opening test truth or test labels. Unconfirmed phenotypes are never treated as negatives.
"""
from __future__ import annotations

import argparse
import json
import math
import random
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import prepare_wp5_cases as P5  # noqa: E402
import run_wp5_recovery as R  # noqa: E402

SCHEMA = 'gspa.wp5.evaluation/1'
PRECISION, MIN_DECISIONS = 0.90, 10
BOOTSTRAP = {'n': 10000, 'seed': 20261029}
FLUX_TOL = 1e-7
GROWTH_TOL = 1e-6
PRIMARY_STRATA = ('catalytic', 'complex_subunit', 'transport_only')
CALLS = ('correct_call', 'wrong_call')
REFERENCE_ARMS = ('sequence_only', 'exhaustive')
KEIO_NOTE = ('Keio (calibration) measured pairs were chosen by an iML1515 design that used the target GPR; '
             'the menu offered the full substrate grid plus glucose.')
PROTOCOL_NOTES = [
    'Nutrient assays use oxygen uptake 20 mmol gDW^-1 h^-1; printed is primary and additional_inorganic is sensitivity.',
    'For supplemented nutrient conditions, the wild type is the corresponding unsupplemented wild type.',
    'Sulfonate sulfur uptake ceiling is 1.0; tauA and ssuA share a target-free background menu.',
    'No oxygen sweep, product objective, reaction-level condition knockout, learned predictor or Tier B case is admitted.',
]


# --------------------------------------------------------------------------- measured classes and cells

def measured_class(label):
    """Return growth/no_growth/partial/unknown/unmeasured/wt_unverified (§6)."""
    if label is None:
        return 'unmeasured'
    if label.get('wt_status') not in ('growth', 'stated_by_source'):
        return 'wt_unverified'
    status = label.get('mutant_status')
    if status in ('growth', 'no_growth'):
        return status if label.get('stable') is True else 'unknown'
    return 'partial' if status == 'partial' else 'unknown'


def prediction(cell):
    contrast = cell.get('contrast')
    if contrast == 'both_grow':
        return 'growth'
    if contrast in ('knockout_no_growth', 'knockout_infeasible_with_maintenance'):
        return 'no_growth'
    return None


def score_of(h):
    values = list((h.get('evidence') or {}).values())
    return max(values, default=0.0)


def fva_class(interval):
    if not interval or interval.get('minimum') is None or interval.get('maximum') is None:
        return 'unassessed'
    lo, hi = interval['minimum'], interval['maximum']
    if abs(lo) <= FLUX_TOL and abs(hi) <= FLUX_TOL:
        return 'no_flux_at_floor'
    if lo > FLUX_TOL or hi < -FLUX_TOL:
        return 'required'
    return 'optional'


def truth_match(h, truth):
    return bool(h.get('assignments')) and set(h['assignments']) == set(truth['reactions'])


def call(h, truth):
    correct = truth_match(h, truth) and truth['stratum'] != 'transport_only'
    overlap = bool(set(h['assignments']) & set(truth['reactions']))
    return {'decision': 'correct_call' if correct else 'wrong_call', 'called': h['id'],
            'called_reactions': sorted(h['assignments']), 'partial_overlap': overlap and not correct,
            'transport_catalytic_call': truth['stratum'] == 'transport_only'}


def _cell_predictions(public, run):
    """Validate the visited rectangle and return condition -> hypothesis -> binary prediction."""
    hyps = {h['id'] for h in public['hypotheses']}
    conditions = {c['id'] for c in public['conditions']}
    visited = run.get('visited_conditions')
    cells = run.get('cells')
    if not isinstance(visited, list) or not isinstance(cells, list) or len(set(visited)) != len(visited):
        raise P5.WP5Error('invalid visited-condition list')
    if not set(visited) <= conditions:
        raise P5.WP5Error('run visited a condition outside the public menu')
    expected = {(cid, hid) for cid in visited for hid in hyps}
    seen = set()
    preds = defaultdict(dict)
    for cell in cells:
        pair = (cell.get('condition_id'), cell.get('hypothesis_id'))
        if pair not in expected or pair in seen:
            raise P5.WP5Error('incomplete or duplicated hypothesis-by-condition rectangle')
        seen.add(pair)
        preds[pair[0]][pair[1]] = prediction(cell)
    if seen != expected:
        raise P5.WP5Error('run lacks a complete hypothesis-by-condition rectangle')
    return preds


def decide(public, run, truth, measured):
    """Apply the frozen §6 case rule; ``measured`` maps condition ID to measured class."""
    hyps = {h['id']: h for h in public['hypotheses']}
    null_rows = [h for h in public['hypotheses'] if not h.get('assignments')]
    cands = [h for h in public['hypotheses'] if h.get('assignments')]
    base = {'truth_in_candidates': any(truth_match(h, truth) for h in cands), 'n_hypotheses': len(hyps)}
    if len(null_rows) != 1:
        raise P5.WP5Error('exactly one null hypothesis is required')
    null = null_rows[0]['id']
    if run is None or run.get('status') != 'completed':
        return dict(base, decision='inconclusive', reason=(run or {}).get('status', 'missing_run'),
                    pairs=[], discriminating_visited=False, visited=[])
    try:
        preds = _cell_predictions(public, run)
    except P5.WP5Error as error:
        return dict(base, decision='inconclusive', reason=str(error), pairs=[],
                    discriminating_visited=False, visited=[])
    cells = {(c['condition_id'], c['hypothesis_id']): c for c in run['cells']}
    if any(cells[key].get('contrast') == 'inconclusive' for key in cells):
        return dict(base, decision='inconclusive', reason='solver_or_budget_inconclusive', pairs=[],
                    discriminating_visited=False, visited=list(run['visited_conditions']))

    contradicted, support = Counter(), defaultdict(list)
    pairs, discriminating = [], False
    for cid in run['visited_conditions']:
        p = preds[cid]
        values = {v for v in p.values() if v is not None}
        differs = len(values) > 1
        discriminating |= differs
        m = measured.get(cid, 'unmeasured')
        for hid, value in p.items():
            # Partial, unknown, unmeasured and WT-unverified labels neither confirm nor contradict.
            if value is None or m not in ('growth', 'no_growth'):
                continue
            if value != m:
                contradicted[hid] += 1
            elif differs:
                support[hid].append(cid)
        pairs.append({'condition_id': cid, 'measured': m, 'predictions': dict(p)})
    out = dict(base, pairs=pairs, discriminating_visited=discriminating,
               contradicted=sorted(contradicted), visited=list(run['visited_conditions']))
    consistent = [h for h in cands if not contradicted[h['id']]]
    if all(contradicted[h] for h in hyps):
        return dict(out, decision='no_consistent_hypothesis')
    if consistent:
        top = max(score_of(h) for h in consistent)
        best = [h for h in consistent if score_of(h) == top]
        if len(best) > 1:
            return dict(out, decision='abstain_tied', score=top)
        chosen = best[0]
        if not support[chosen['id']]:
            return dict(out, decision='abstain_unsupported', score=top)
        return dict(out, **call(chosen, truth), score=top, supporting=support[chosen['id']])
    if not contradicted[null] and support[null]:
        negative = truth['stratum'] == 'transport_only' or not base['truth_in_candidates']
        return dict(out, decision='null_call', null_outcome='correct_negative' if negative else 'false_null',
                    called=null, supporting=support[null])
    return dict(out, decision='abstain_unsupported')


def sequence_only(public, truth):
    cands = [h for h in public['hypotheses'] if h.get('assignments')]
    base = {'truth_in_candidates': any(truth_match(h, truth) for h in cands),
            'n_hypotheses': len(public['hypotheses']), 'pairs': [],
            'discriminating_visited': False, 'visited': []}
    if not cands:
        return dict(base, decision='no_admissible_experiment')
    top = max(score_of(h) for h in cands)
    best = [h for h in cands if score_of(h) == top]
    if len(best) > 1:
        return dict(base, decision='abstain_tied', score=top)
    return dict(base, **call(best[0], truth), score=top)


def capacity(public, run, truth, decision):
    """C2 FVA classes for called and truth-matching assignments in growing-control cells."""
    out = {'called': Counter(), 'truth': Counter()}
    if run is None or run.get('status') != 'completed':
        return {k: dict(v) for k, v in out.items()}
    truth_ids = {h['id'] for h in public['hypotheses'] if truth_match(h, truth)}
    by_id = {h['id']: h for h in public['hypotheses']}
    for cell in run.get('cells', []):
        wt = (cell.get('phenotypes') or {}).get('wild_type') or {}
        if not (wt.get('status') == 'optimal' and (wt.get('growth') or 0) > GROWTH_TOL):
            continue
        for role, ids in (('called', {decision.get('called')}), ('truth', truth_ids)):
            hid = cell.get('hypothesis_id')
            if hid in ids and hid in by_id:
                for rid in sorted(by_id[hid]['assignments']):
                    out[role][fva_class((wt.get('flux_intervals') or {}).get(rid))] += 1
    return {k: dict(v) for k, v in out.items()}


# --------------------------------------------------------------------------- public/run validation and joining

def load_inputs(bundle, runs):
    bundle, runs = Path(bundle), Path(runs)
    public = bundle / 'public'
    _, _, _, records = R.load_public(public)
    protocol = P5.read_json(runs / 'protocol.json')
    if protocol.get('schema') != R.SCHEMA + ':protocol':
        raise P5.WP5Error('invalid recovery protocol schema')
    if protocol['public_manifest_sha256'] != P5.sha_file(public / 'manifest.json'):
        raise P5.WP5Error('runs were made on different public inputs')
    audit = P5.read_json(public / 'masking_audit.json')
    if audit.get('status') != 'pass' or protocol['masking_audit_sha256'] != P5.sha_file(public / 'masking_audit.json'):
        raise P5.WP5Error('masking audit failed or changed: the run is invalid')
    cases = {r['problem']: r for r in P5.read_jsonl(public / 'cases.jsonl')}
    if len(cases) != sum(1 for _ in P5.read_jsonl(public / 'cases.jsonl')):
        raise P5.WP5Error('duplicate public problem')
    if set(cases) != set(records) or set(cases) != set(protocol['problems']):
        raise P5.WP5Error('protocol problems differ from the prepared cases')
    return protocol, cases, records


def _invalid(run, status):
    return dict(run, status=status)


def run_index(runs, protocol):
    """Return planned runs and exact hashes; provenance/budget/audit failures remain inconclusive."""
    runs = Path(runs)
    protocol_path = runs / 'protocol.json'
    protocol_sha = P5.sha_file(protocol_path)
    audited, audited_protocol = {}, None
    tpath = runs / 'transcript_audit.json'
    if tpath.is_file():
        audit = P5.read_json(tpath)
        if audit.get('status') == 'pass':
            audited = audit.get('transcripts_sha256') or {}
            audited_protocol = audit.get('protocol_sha256')
    index, files = {}, {}
    for stage in ('exhaustive', 'deterministic', 'qwen'):
        for name, arm, k, seed in R.plan_runs(protocol, stage):
            path = runs / 'runs' / R.run_name(name, arm, k, seed)
            run = P5.read_json(path) if path.is_file() else None
            if run is None:
                index[(name, arm, k, seed)] = None
                continue
            files[str(path.relative_to(runs))] = P5.sha_file(path)
            entry = protocol['problems'][name]
            expected_cap = entry['n_conditions'] * entry['cost_per_condition'] if arm == 'exhaustive' \
                else entry['caps'][str(k)]['max_solver_calls']
            expected_exp = entry['n_conditions'] if arm == 'exhaustive' else entry['caps'][str(k)]['max_experiments']
            if (run.get('protocol_sha256') != protocol_sha
                    or run.get('public_problem_sha256') != entry['public_problem_sha256']
                    or run.get('problem') != name or run.get('arm') != arm
                    or run.get('k') != k or run.get('seed') != seed):
                run = _invalid(run, 'invalid_provenance')
            elif run.get('max_solver_calls') != expected_cap or run.get('max_experiments') != expected_exp:
                run = _invalid(run, 'invalid_budget')
            else:
                calls = run.get('solver_calls')
                used = run.get('solver_calls_used')
                if (not isinstance(calls, list) or used != len(calls) or not 0 <= used <= expected_cap
                        or run.get('solver_calls_unused') != expected_cap - used):
                    run = _invalid(run, 'invalid_solver_ledger')
            if arm == 'qwen':
                transcript_name = R.transcript_name(name, arm, k, seed)
                transcript = runs / transcript_name
                if (run.get('transcript') != transcript_name or not transcript.is_file()
                        or audited_protocol != protocol_sha):
                    run = _invalid(run, 'transcript_unaudited')
                else:
                    sha = P5.sha_file(transcript)
                    files[transcript_name] = sha
                    if audited.get(transcript_name) != sha or run.get('transcript_sha256') != sha:
                        run = _invalid(run, 'transcript_unaudited')
            index[(name, arm, k, seed)] = run
    # Fail closed on runs that did not come from the production selection loop (ruling §13.5) and on arm
    # cells that differ from the exhaustive cell for the same (case, hypothesis, condition) (design §5).
    for key, run in list(index.items()):
        if run is None or run.get('status') != 'completed':
            continue
        if run.get('engine_entry') != R.ENGINE_ENTRY:
            index[key] = _invalid(run, 'invalid_engine')
            continue
        if key[1] == 'exhaustive':
            continue
        exhaustive = index.get((key[0], 'exhaustive', 0, 0))
        if exhaustive is None or exhaustive.get('status') != 'completed' or \
                exhaustive.get('engine_entry') != R.ENGINE_ENTRY or R.check_against_exhaustive(run, exhaustive):
            index[key] = _invalid(run, 'cell_mismatch')
    return index, files


def _unique(rows, key, description):
    out = {}
    for row in rows:
        value = key(row)
        if value in out:
            raise P5.WP5Error('duplicate ' + description)
        out[value] = row
    return out


def measured_map(record, labels):
    variant = record['problem'].split('-', 1)[1]
    rows = labels.get((record['case'], variant), [])
    by_index = _unique(rows, lambda r: r['menu_index'], 'measured label')
    # Label rows carry integer menu indices; the menu map is keyed by their string form.
    if not {str(i) for i in by_index} <= set(record.get('menu_map', {})):
        raise P5.WP5Error('label refers to a condition outside the frozen menu')
    return {cid: measured_class(by_index.get(int(i))) for i, cid in record.get('menu_map', {}).items()}


def arm_keys(protocol):
    keys = [('sequence_only', 0, 0), ('exhaustive', 0, 0)]
    for k in protocol['budgets']:
        for arm in protocol['arms']:
            for seed in (protocol['random_seeds'] if arm == 'random' else [0]):
                keys.append((arm, k, seed))
    return keys


def arm_label(row):
    return row['arm'] if row['arm'] in REFERENCE_ARMS else f'{row["arm"]}:k{row["k"]}'


def join_split(bundle, runs, split, policy_verified=False):
    """Join one private split. Test joins require the evaluate phase's verified-policy token."""
    if split not in ('calibration', 'test'):
        raise P5.WP5Error('unknown split')
    if split == 'test' and not policy_verified:
        raise P5.WP5Error('test labels are joined only by evaluate after policy verification')
    bundle = Path(bundle)
    protocol, cases, records = load_inputs(bundle, runs)
    truth_rows = P5.read_jsonl(bundle / f'private/truth.{split}.jsonl')
    truth = _unique(truth_rows, lambda t: t['case'], 'truth case')
    expected_cases = {r['case'] for r in cases.values() if r['split'] == split}
    if set(truth) != expected_cases:
        raise P5.WP5Error('truth file does not cover exactly its split')
    labels = defaultdict(list)
    for row in P5.read_jsonl(bundle / f'private/labels.{split}.jsonl'):
        labels[(row['case'], row['variant'])].append(row)
    if {case for case, _ in labels} - expected_cases:
        raise P5.WP5Error('label file contains another split')
    index, files = run_index(runs, protocol)
    rows = []
    for name, case in sorted(cases.items()):
        if case['split'] != split:
            continue
        record, t = records[name], truth[case['case']]
        common = {'problem': name, 'case': case['case'], 'source': case['source'], 'species': case['species'],
                  'variant': case['variant'], 'primary_variant': case['primary_variant'], 'split': split,
                  'stratum': t['stratum'], 'eligible': t['stratum'] in PRIMARY_STRATA}
        if record['status'] != 'constructed':
            for arm, k, seed in arm_keys(protocol):
                rows.append(dict(common, arm=arm, k=k, seed=seed, decision='no_admissible_experiment',
                                 truth_in_candidates=False, n_hypotheses=0, pairs=[],
                                 discriminating_visited=False, visited=[], capacity={'called': {}, 'truth': {}}))
            continue
        public = record['public_problem']
        measured = measured_map(record, labels)
        exhaustive = index.get((name, 'exhaustive', 0, 0))
        for arm, k, seed in arm_keys(protocol):
            if arm == 'sequence_only':
                decision, run = sequence_only(public, t), None
            else:
                run = exhaustive if arm == 'exhaustive' else index.get((name, arm, k, seed))
                decision = decide(public, run, t, measured)
            decision['capacity'] = capacity(public, run, t, decision)
            if run is not None:
                decision['cost'] = {
                    'solver_calls_used': run.get('solver_calls_used'), 'max_solver_calls': run.get('max_solver_calls'),
                    'solver_calls_unused': run.get('solver_calls_unused'), 'wall_seconds': run.get('wall_seconds'),
                    'agent_usage': run.get('agent_usage')}
            rows.append(dict(common, arm=arm, k=k, seed=seed, **decision))
    return rows, files


def join(bundle, runs, out):
    out = Path(out)
    if (out / 'joined.calibration.jsonl').exists() or (out / 'join.receipt.json').exists():
        raise P5.WP5Error('calibration join exists; outputs are immutable')
    rows, files = join_split(bundle, runs, 'calibration')
    P5.write_jsonl(out / 'joined.calibration.jsonl', rows)
    receipt = {'schema': SCHEMA + ':join', 'split': 'calibration', 'run_files': files,
               'joined_sha256': P5.sha_file(out / 'joined.calibration.jsonl'),
               'protocol_sha256': P5.sha_file(Path(runs) / 'protocol.json')}
    P5.write_json(out / 'join.receipt.json', receipt)
    return rows


# --------------------------------------------------------------------------- calibration and frozen policy

def gate(row, threshold):
    if row['decision'] in CALLS and (threshold is None or row.get('score', 0) < threshold):
        return dict(row, decision='abstain_policy', gated_from=row['decision'])
    return row


def weighted_counts(rows, threshold):
    """Decision counts; random seeds have per-case weights summing to one."""
    seeds = Counter(r['problem'] for r in rows)
    counts = Counter()
    for row in rows:
        counts[gate(row, threshold)['decision']] += 1 / seeds[row['problem']]
    return counts


def calibrate_rows(rows):
    by_arm = defaultdict(list)
    for row in rows:
        if row['eligible'] and row['primary_variant']:
            by_arm[arm_label(row)].append(row)
    policy = {}
    for arm, items in sorted(by_arm.items()):
        scores = sorted({r['score'] for r in items if r['decision'] in CALLS})
        choices = []
        for threshold in scores:
            counts = weighted_counts(items, threshold)
            # Unknown active decisions, if a future adapter emits them, reduce demonstrated precision.
            n = counts['correct_call'] + counts['wrong_call'] + counts['unknown_call']
            precision = counts['correct_call'] / n if n else 0.0
            if n >= MIN_DECISIONS and precision >= PRECISION:
                choices.append((round(counts['correct_call'], 9), -round(counts['wrong_call'], 9), threshold,
                                {k: round(v, 6) for k, v in counts.items()}))
        best = max(choices, default=None)
        policy[arm] = {'threshold': best[2] if best else None,
                       'calibration_counts': best[3] if best else {},
                       'status': 'selected_on_calibration' if best else 'no_supported_threshold_abstain'}
    return policy


def calibrate(bundle, runs, out):
    out = Path(out)
    receipt = P5.read_json(out / 'join.receipt.json')
    if receipt.get('schema') != SCHEMA + ':join' or receipt.get('split') != 'calibration':
        raise P5.WP5Error('invalid calibration join receipt')
    if P5.sha_file(out / 'joined.calibration.jsonl') != receipt['joined_sha256']:
        raise P5.WP5Error('calibration join changed')
    if (out / 'policy.json').exists():
        raise P5.WP5Error('a policy is already frozen')
    protocol = P5.read_json(Path(runs) / 'protocol.json')
    if receipt['protocol_sha256'] != P5.sha_file(Path(runs) / 'protocol.json'):
        raise P5.WP5Error('protocol changed after calibration join')
    _, current_files = run_index(runs, protocol)
    if current_files != receipt['run_files']:
        raise P5.WP5Error('run files changed after calibration join')
    rows = P5.read_jsonl(out / 'joined.calibration.jsonl')
    policy = {
        'schema': SCHEMA + ':policy', 'target_precision': PRECISION, 'min_decisions': MIN_DECISIONS,
        'selection_rule': 'max_correct_then_min_harm_then_highest_threshold; unknowns count against precision',
        'score': 'max evidence channel of the called hypothesis', 'arms': calibrate_rows(rows),
        'protocol_sha256': receipt['protocol_sha256'], 'run_files': receipt['run_files'],
        'calibration_join_sha256': receipt['joined_sha256'],
        'truth_calibration_sha256': P5.sha_file(Path(bundle) / 'private/truth.calibration.jsonl'),
        'labels_calibration_sha256': P5.sha_file(Path(bundle) / 'private/labels.calibration.jsonl'),
        'note': KEIO_NOTE, 'limitation': 'empirical calibration precision, not a population guarantee'}
    P5.write_json(out / 'policy.json', policy)
    return P5.sha_file(out / 'policy.json')


def verify_policy(bundle, runs, out, policy_sha):
    out = Path(out)
    if not isinstance(policy_sha, str) or not re.fullmatch(r'[0-9a-f]{64}', policy_sha):
        raise P5.WP5Error('a concrete lowercase SHA-256 policy hash is required')
    path = out / 'policy.json'
    if not path.is_file() or P5.sha_file(path) != policy_sha:
        raise P5.WP5Error('policy hash mismatch: refusing to evaluate')
    policy = P5.read_json(path)
    if (policy.get('schema') != SCHEMA + ':policy'
            or policy['protocol_sha256'] != P5.sha_file(Path(runs) / 'protocol.json')):
        raise P5.WP5Error('policy was frozen for another protocol')
    protocol = P5.read_json(Path(runs) / 'protocol.json')
    _, files = run_index(runs, protocol)
    if files != policy['run_files']:
        raise P5.WP5Error('run files changed after the policy was frozen')
    if P5.sha_file(out / 'joined.calibration.jsonl') != policy['calibration_join_sha256']:
        raise P5.WP5Error('calibration join changed after the policy was frozen')
    for kind in ('truth', 'labels'):
        if P5.sha_file(Path(bundle) / f'private/{kind}.calibration.jsonl') != policy[f'{kind}_calibration_sha256']:
            raise P5.WP5Error('calibration truth or labels changed after the policy was frozen')
    return policy


# --------------------------------------------------------------------------- metrics and statistics

def wilson(x, n, z=1.959963984540054):
    if not n:
        return [None, None]
    p = x / n
    centre = (p + z * z / (2 * n)) / (1 + z * z / n)
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return [max(0.0, centre - half), min(1.0, centre + half)]


def mcnemar(a, b):
    shared = sorted(set(a) & set(b))
    only_a = sum(1 for case in shared if a[case] and not b[case])
    only_b = sum(1 for case in shared if b[case] and not a[case])
    n = only_a + only_b
    if n == 0:
        return {'only_a': 0, 'only_b': 0, 'p': 1.0, 'n_pairs': len(shared)}
    tail = sum(math.comb(n, i) for i in range(min(only_a, only_b) + 1)) / 2 ** n
    return {'only_a': only_a, 'only_b': only_b, 'p': min(1.0, 2 * tail), 'n_pairs': len(shared)}


def bootstrap_difference(a, b, n=BOOTSTRAP['n'], seed=BOOTSTRAP['seed']):
    cases = sorted(set(a) & set(b))
    if not cases:
        return [None, None]
    rng = random.Random(seed)
    diffs = sorted(sum(a[c] - b[c] for c in (rng.choice(cases) for _ in cases)) / len(cases)
                   for _ in range(n))
    return [diffs[int(0.025 * n)], diffs[min(n - 1, int(0.975 * n))]]


def per_case(rows, threshold):
    groups = defaultdict(list)
    for row in rows:
        groups[row['case']].append(gate(row, threshold))
    return groups


def summarize(items, threshold):
    groups = per_case(items, threshold)
    n = len(groups)
    counts = Counter()
    for group in groups.values():
        for row in group:
            key = row['decision'] if row['decision'] != 'null_call' else row['null_outcome']
            counts[key] += 1 / len(group)
    recall = sum(1 for group in groups.values() if group[0]['truth_in_candidates'])
    with_truth = [group for group in groups.values() if group[0]['truth_in_candidates']]
    conditional = sum(sum(row['decision'] == 'correct_call' for row in group) / len(group)
                      for group in with_truth)
    abstain = sum(value for key, value in counts.items()
                  if key.startswith('abstain') or key == 'no_consistent_hypothesis')
    transport = sum(sum(bool(row.get('transport_catalytic_call')) and row['decision'] == 'wrong_call'
                        for row in group) / len(group) for group in groups.values())
    partial = sum(sum(bool(row.get('partial_overlap')) and row['decision'] == 'wrong_call'
                      for row in group) / len(group) for group in groups.values())
    return {
        'denominator': n, 'counts': {k: round(v, 6) for k, v in sorted(counts.items())},
        'recovery': counts['correct_call'] / n if n else None,
        'recovery_wilson95': wilson(counts['correct_call'], n),
        'harm': counts['wrong_call'] / n if n else None,
        'abstention': abstain / n if n else None,
        'false_null': counts['false_null'] / n if n else None,
        'correct_negative': counts['correct_negative'] / n if n else None,
        'transport_catalytic_calls': round(transport, 6), 'partial_overlap_calls': round(partial, 6),
        'candidate_recall': recall / n if n else None, 'candidate_recall_wilson95': wilson(recall, n),
        'conditional_recovery': conditional / len(with_truth) if with_truth else None,
        'conditional_denominator': len(with_truth)}


def _agreement(items, threshold, pair_rows):
    out = Counter()
    for row in items:
        gated = gate(row, threshold)
        hid = gated.get('called') if gated['decision'] in ('correct_call', 'wrong_call', 'null_call') else None
        for pair in pair_rows(row):
            measured = pair['measured']
            if measured == 'partial':
                out['partial_pairs'] += 1
                continue
            if measured not in ('growth', 'no_growth'):
                out['not_evaluable_pairs'] += 1
                continue
            out['binary_pairs'] += 1
            out['constant_growth_agree'] += measured == 'growth'
            out['constant_no_growth_agree'] += measured == 'no_growth'
            if hid is None:
                continue
            value = pair['predictions'].get(hid)
            if value is None:
                out['not_evaluable_called_pairs'] += 1
                continue
            out['called_pairs'] += 1
            agrees = value == measured
            out['called_agree'] += agrees
            out['agree' if agrees else 'contradict'] += 1
    result = dict(out)
    result['called_agreement'] = out['called_agree'] / out['called_pairs'] if out['called_pairs'] else None
    result['called_coverage'] = out['called_pairs'] / out['binary_pairs'] if out['binary_pairs'] else None
    result['constant_growth_agreement'] = out['constant_growth_agree'] / out['binary_pairs'] if out['binary_pairs'] else None
    result['constant_no_growth_agreement'] = out['constant_no_growth_agree'] / out['binary_pairs'] if out['binary_pairs'] else None
    return result


def phenotype_agreement(items, threshold, exhaustive_items=None):
    """Visited-pair agreement plus the common exhaustive measured set; partials are separate."""
    result = _agreement(items, threshold, lambda row: row.get('pairs', []))
    references = {r['problem']: r for r in (exhaustive_items or []) if r['arm'] == 'exhaustive'}

    def exhaustive_pairs(row):
        reference = references.get(row['problem'])
        return reference.get('pairs', []) if reference and reference['decision'] != 'inconclusive' else []

    result['common_exhaustive'] = _agreement(items, threshold, exhaustive_pairs)
    result['unmasked_original'] = {
        'status': 'not_available',
        'reason': 'The §8 public-only runner has no unmasked-original reference stage; no private model is introduced here.'}
    return result


def capacity_summary(items):
    eligible = [r for r in items if r.get('n_hypotheses', 0) >= 2 and r['decision'] != 'no_admissible_experiment']
    fva = {'called': Counter(), 'truth': Counter()}
    for row in items:
        for role in fva:
            fva[role].update((row.get('capacity') or {}).get(role, {}))
    seeds = Counter(r['case'] for r in eligible)
    hits = sum(r.get('discriminating_visited', False) / seeds[r['case']] for r in eligible)
    cases = len(seeds)
    return {'c1_cases': cases, 'c1_discriminating': round(hits, 6),
            'c1_fraction': hits / cases if cases else None,
            'c2_fva_classes': {k: dict(v) for k, v in fva.items()}}


def cost_summary(items):
    total = Counter()
    for row in items:
        cost = row.get('cost') or {}
        for key in ('solver_calls_used', 'max_solver_calls', 'solver_calls_unused', 'wall_seconds'):
            total[key] += cost.get(key) or 0
        usage = cost.get('agent_usage') or {}
        for key in ('requests', 'failures', 'seconds'):
            total['qwen_' + key] += usage.get(key, 0)
        for key, value in (usage.get('transport') or {}).items():
            total['qwen_' + key] += value
        total['runs'] += 1
        total['failed_runs'] += row['decision'] == 'inconclusive'
    return dict(total)


def precision_recovery_curve(items):
    """Descriptive test curve; never used to select the frozen policy threshold."""
    scores = sorted({r['score'] for r in items if r['decision'] in CALLS})
    curve = []
    for threshold in scores:
        gated = [gate(r, threshold) for r in items]
        summary = summarize(gated, threshold)
        correct = summary['counts'].get('correct_call', 0)
        wrong = summary['counts'].get('wrong_call', 0)
        decisions = correct + wrong
        curve.append({'threshold': threshold, 'decisions': round(decisions, 6),
                      'precision': correct / decisions if decisions else None,
                      'recovery': summary['recovery'], 'harm': summary['harm']})
    return curve


def outcomes(items, threshold):
    return {case: sum(row['decision'] == 'correct_call' for row in group) / len(group)
            for case, group in per_case([r for r in items if r['eligible']], threshold).items()}


def paired_tests(by_arm, policy):
    primary = {arm: items for (arm, primary), items in by_arm.items() if primary}
    scored = {arm: outcomes(items, policy['arms'].get(arm, {}).get('threshold'))
              for arm, items in primary.items()}
    pairs = []
    for arm in sorted(scored):
        k = arm.split(':k')[-1] if ':k' in arm else None
        comparators = ['sequence_only'] + ([f'minimal_medium:k{k}'] if k else [])
        if arm.startswith('qwen:') and k:
            comparators += [f'{other}:k{k}' for other in ('fixed_panel', 'random', 'evidence')]
        for other in comparators:
            if other == arm or other not in scored:
                continue
            a, b = scored[arm], scored[other]
            row = {'arm': arm, 'comparator': other,
                   'bootstrap95_recovery_difference': bootstrap_difference(a, b)}
            if not arm.startswith('random') and not other.startswith('random'):
                row['mcnemar_exact'] = mcnemar({c: v == 1 for c, v in a.items()},
                                               {c: v == 1 for c, v in b.items()})
            pairs.append(row)
    return pairs


def evaluate(bundle, runs, out, policy_sha):
    out = Path(out)
    policy = verify_policy(bundle, runs, out, policy_sha)
    if (out / 'evaluation.json').exists() or (out / 'joined.test.jsonl').exists():
        raise P5.WP5Error('evaluation exists; outputs are immutable')
    # This is the first and only test-label join, after every policy/provenance check above.
    rows, _ = join_split(bundle, runs, 'test', policy_verified=True)
    P5.write_jsonl(out / 'joined.test.jsonl', rows)
    by_arm = defaultdict(list)
    for row in rows:
        by_arm[(arm_label(row), row['primary_variant'])].append(row)
    tables = {}
    for (arm, primary), items in sorted(by_arm.items()):
        threshold = policy['arms'].get(arm, {}).get('threshold')
        status = policy['arms'].get(arm, {}).get('status', 'no_policy_abstain')
        eligible = [r for r in items if r['eligible']]
        strata = {}
        keys = {('stratum', r['stratum']) for r in eligible} | {('species', r['species']) for r in eligible}
        for key in sorted(keys):
            strata[f'{key[0]}:{key[1]}'] = summarize([r for r in eligible if r[key[0]] == key[1]], threshold)
        exhaustive_items = [r for r in rows if r['primary_variant'] == primary]
        tables[f'{arm}|{"primary" if primary else "sensitivity"}'] = {
            'arm': arm, 'variant_scope': 'primary' if primary else 'sensitivity',
            'policy_status': status, 'threshold': threshold, **summarize(eligible, threshold),
            'reported_only': dict(Counter(r['stratum'] for r in items if not r['eligible'])),
            'precision_recovery_curve': precision_recovery_curve(eligible),
            'phenotype_agreement': phenotype_agreement(eligible, threshold, exhaustive_items),
            'capacity': capacity_summary(eligible), 'cost': cost_summary(items), 'strata': strata}
    result = {
        'schema': SCHEMA, 'policy_sha256': policy_sha, 'tables': tables, 'tests': paired_tests(by_arm, policy),
        'bootstrap': BOOTSTRAP, 'calibration_policy': policy['arms'],
        'test_truth_sha256': P5.sha_file(Path(bundle) / 'private/truth.test.jsonl'),
        'test_labels_sha256': P5.sha_file(Path(bundle) / 'private/labels.test.jsonl'),
        'notes': [KEIO_NOTE] + PROTOCOL_NOTES + [
            'Small denominators: every proportion carries its Wilson interval and denominator.',
            'Random arm: per-case mean over seeds; it enters bootstrap differences, not exact McNemar.',
            'Candidate-generation failures and transport controls remain in the all-eligible denominator.',
            'The unmasked-original phenotype comparator has no §8 runner stage and is marked unavailable, never invented.',
        ]}
    P5.write_json(out / 'evaluation.json', result)
    return result


# --------------------------------------------------------------------------- Markdown report

def fmt(value, digits=3):
    if value is None:
        return '–'
    if isinstance(value, float):
        return f'{value:.{digits}f}'
    return str(value)


def report(out):
    out = Path(out)
    result = P5.read_json(out / 'evaluation.json')
    if result.get('schema') != SCHEMA:
        raise P5.WP5Error('invalid evaluation schema')
    lines = ['# WP5 masked recovery: evaluation tables', '',
             'Function recovery, phenotype agreement and model-conditional capacity are separate. '
             'Simulation agreement is not biological confirmation.', '',
             '## Function recovery at fixed precision (test split)', '',
             'Denominators include every eligible test case, including candidate-generation failures. '
             'Random seeds are averaged per case.', '',
             '| Arm | Scope | Policy | Threshold | N | Recovery | 95% CI | Harm | Abstention | False null | Correct negative | Candidate recall | Conditional recovery |',
             '|---|---|---|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|']
    for _, table in sorted(result['tables'].items()):
        ci = table['recovery_wilson95']
        lines.append(
            f'| {table["arm"]} | {table["variant_scope"]} | {table["policy_status"]} | {fmt(table["threshold"])} | '
            f'{table["denominator"]} | {fmt(table["recovery"])} | {fmt(ci[0])}–{fmt(ci[1])} | '
            f'{fmt(table["harm"])} | {fmt(table["abstention"])} | {fmt(table["false_null"])} | '
            f'{fmt(table["correct_negative"])} | {fmt(table["candidate_recall"])} | {fmt(table["conditional_recovery"])} |')

    lines += ['', '## Descriptive test precision–recovery curves (primary scope)', '',
              'These curves do not select or replace the calibration-frozen operating point.', '',
              '| Arm | Threshold | Decisions | Precision | Recovery | Harm |', '|---|---:|---:|---:|---:|---:|']
    for _, table in sorted(result['tables'].items()):
        if table['variant_scope'] != 'primary':
            continue
        for point in table['precision_recovery_curve']:
            lines.append(f'| {table["arm"]} | {fmt(point["threshold"])} | {fmt(point["decisions"])} | '
                         f'{fmt(point["precision"])} | {fmt(point["recovery"])} | {fmt(point["harm"])} |')

    lines += ['', '## Paired comparisons (primary scope)', '',
              '| Arm | Comparator | McNemar only-arm | only-comparator | p | Bootstrap 95% recovery difference |',
              '|---|---|---:|---:|---:|---|']
    for row in result['tests']:
        m = row.get('mcnemar_exact', {})
        interval = row['bootstrap95_recovery_difference']
        lines.append(f'| {row["arm"]} | {row["comparator"]} | {fmt(m.get("only_a"))} | {fmt(m.get("only_b"))} | '
                     f'{fmt(m.get("p"))} | {fmt(interval[0])} to {fmt(interval[1])} |')

    lines += ['', '## Phenotype agreement', '',
              'Called-agreement uses binary-evaluable pairs with a called prediction; coverage shows the '
              'called-pair fraction of all binary pairs. Partial growth is excluded and reported separately.', '',
              '| Arm | Scope | Set | Binary pairs | Called pairs | Called agree | Called agreement | Coverage | Constant growth | Constant no-growth | Partial pairs |',
              '|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|']
    for _, table in sorted(result['tables'].items()):
        agreement = table['phenotype_agreement']
        for label, values in (('visited', agreement), ('common exhaustive', agreement['common_exhaustive'])):
            lines.append(
                f'| {table["arm"]} | {table["variant_scope"]} | {label} | {values.get("binary_pairs", 0)} | '
                f'{values.get("called_pairs", 0)} | {values.get("called_agree", 0)} | '
                f'{fmt(values.get("called_agreement"))} | {fmt(values.get("called_coverage"))} | '
                f'{fmt(values.get("constant_growth_agreement"))} | {fmt(values.get("constant_no_growth_agreement"))} | '
                f'{values.get("partial_pairs", 0)} |')
    lines += ['', 'The unmasked-original comparator is unavailable: the public-only §8 runner has no such stage. '
              'No private original model is silently introduced into evaluation.']

    lines += ['', '## Model-conditional capacity and cost', '',
              '| Arm | Scope | C1 cases | C1 discriminating | C1 fraction | LPs used | LP cap | Unused LPs | Runs | Failed runs | Qwen requests | Qwen failures |',
              '|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
    for _, table in sorted(result['tables'].items()):
        capacity, cost = table['capacity'], table['cost']
        lines.append(
            f'| {table["arm"]} | {table["variant_scope"]} | {capacity["c1_cases"]} | '
            f'{fmt(capacity["c1_discriminating"])} | {fmt(capacity["c1_fraction"])} | '
            f'{fmt(cost.get("solver_calls_used"))} | {fmt(cost.get("max_solver_calls"))} | '
            f'{fmt(cost.get("solver_calls_unused"))} | {fmt(cost.get("runs"))} | {fmt(cost.get("failed_runs"))} | '
            f'{fmt(cost.get("qwen_requests"))} | {fmt(cost.get("qwen_failures"))} |')

    lines += ['', '## Per-stratum and per-species tables (primary scope)', '',
              '| Arm | Group | N | Recovery | Harm | False null | Candidate recall | Conditional recovery |',
              '|---|---|---:|---:|---:|---:|---:|---:|']
    for _, table in sorted(result['tables'].items()):
        if table['variant_scope'] != 'primary':
            continue
        for group, values in sorted(table['strata'].items()):
            lines.append(f'| {table["arm"]} | {group} | {values["denominator"]} | {fmt(values["recovery"])} | '
                         f'{fmt(values["harm"])} | {fmt(values["false_null"])} | '
                         f'{fmt(values["candidate_recall"])} | {fmt(values["conditional_recovery"])} |')

    lines += ['', '## Frozen protocol and limitations', '']
    lines += ['- ' + note for note in result['notes']]
    lines += ['', '## Calibration policy', '', '> ' + KEIO_NOTE, '',
              '| Arm | Status | Threshold | Calibration counts |', '|---|---|---:|---|']
    for arm, policy in sorted(result['calibration_policy'].items()):
        lines.append(f'| {arm} | {policy["status"]} | {fmt(policy["threshold"])} | '
                     f'`{json.dumps(policy["calibration_counts"], sort_keys=True)}` |')
    text = '\n'.join(lines) + '\n'
    (out / 'report.md').write_text(text, encoding='utf-8')
    P5.write_json(out / 'report.json', {'schema': SCHEMA + ':report', 'evaluation_sha256': P5.sha_file(out / 'evaluation.json'),
                                       'report_sha256': P5.sha_file(out / 'report.md')})
    return text


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest='phase', required=True)
    for phase in ('join', 'calibrate', 'evaluate'):
        command = sub.add_parser(phase)
        command.add_argument('--bundle', type=Path, required=True)
        command.add_argument('--runs', type=Path, required=True)
        command.add_argument('--out', type=Path, required=True)
        if phase == 'evaluate':
            command.add_argument('--policy-sha256', required=True)
    command = sub.add_parser('report')
    command.add_argument('--out', type=Path, required=True)
    args = parser.parse_args(argv)
    if args.phase == 'join':
        rows = join(args.bundle, args.runs, args.out)
        print(json.dumps({'calibration_rows': len(rows)}))
    elif args.phase == 'calibrate':
        print(calibrate(args.bundle, args.runs, args.out))
    elif args.phase == 'evaluate':
        result = evaluate(args.bundle, args.runs, args.out, args.policy_sha256)
        print(json.dumps({'tables': len(result['tables']), 'policy_sha256': args.policy_sha256}))
    else:
        report(args.out)
        print(str(args.out / 'report.md'))
    return 0


if __name__ == '__main__':
    sys.exit(main())
