#!/usr/bin/env python3
"""Truth-separated comparison of production-selected Stage-3 anomalies.

prepare reads only content-verified production runs. suggest reads only the
public bundle. calibrate reads calibration truth; evaluate reads test truth only
after requiring a frozen policy. Unknown labels are never converted to negatives.
The invalid v1 generator is available in git history, not as an executable mode.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
from pathlib import Path
import sys
import time
import concurrent.futures
import threading
from collections import Counter, defaultdict

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'benchmark/agent'))
import run_gspa_agent as agent
import v3_runner as runner

MODES = ('abstain_only', 'rank_by_raw_seq', 'rank_by_stage1_posterior',
         'rank_by_context_support', 'learned_logistic', 'learned_gbdt')
FORBIDDEN = {'truth', 'label', 'labels', 'ground_truth', 'expected_action',
             'expected_protein_id', 'expected_function_id', 'is_true_target',
             'valid_ancestors_in_truth', 'y_true', 'is_true'}


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
                                     allow_nan=False).encode()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def write(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    # Outputs are immutable: a test run cannot silently replace its frozen policy.
    with path.open('x') as fh:
        json.dump(data, fh, indent=2, sort_keys=True, allow_nan=False)
        fh.write('\n')


def reject_truth(value):
    if isinstance(value, dict):
        if FORBIDDEN & set(value):
            raise ValueError('truth fields in public evidence')
        for child in value.values():
            reject_truth(child)
    elif isinstance(value, list):
        for child in value:
            reject_truth(child)


def probability(value):
    if isinstance(value, bool) or not isinstance(value, (float, int)) or not math.isfinite(value) or not 0 <= value <= 1:
        raise ValueError(f'invalid probability: {value!r}')
    return float(value)


def verified_run(path):
    old = read(path)
    identity = old.get('inference_identity')
    if not identity:
        raise ValueError('production run lacks input identity')
    current = runner.load_run_record(str(Path(path).parent), identity['genomes'],
                                    {'gspa': identity['launch']}, identity['arm'])
    if len(current['genomes']) != len(identity['genomes']):
        raise ValueError(f'unverified production run: {current["failed"]}')
    return old


def prepare(source):
    """No truth-path argument and no reading of truth files in this phase."""
    if source.get('schema') != 'gspa.agent-evidence-sources/2':
        raise ValueError('expected production evidence sources v2')
    runs = {name: verified_run(path) for name, path in source['runs'].items()}
    if set(runs) != {'stage1', 'learned_logistic', 'learned_gbdt'}:
        raise ValueError('require baseline plus actual logistic and GBDT production runs')
    base = runs['stage1']
    for name, run in runs.items():
        if run['inference_identity']['genomes'] != base['inference_identity']['genomes']:
            raise ValueError(f'{name}: unmatched claims or genome context')
        if run['inference_identity']['launch'] != base['inference_identity']['launch']:
            raise ValueError(f'{name}: unmatched production configuration')
        if name != 'stage1':
            args = run['inference_identity']['arm']['args']
            if len(args) != 2 or args[0] != '--learned-model':
                raise ValueError('model arm must differ only by its frozen learned artifact')
            model_path = str(Path(args[1]).resolve())
            expected_files = {k: v for k, v in run['inference_identity']['files'].items() if k != model_path}
            if expected_files != base['inference_identity']['files']:
                raise ValueError('model arm reference or executable differs from baseline')
    cases, coverage = [], []
    for gid in sorted(base['genomes']):
        split = source['splits'][gid]['split']
        if split not in ('calibration', 'test'):
            raise ValueError('only calibration and test cases are admitted')
        output = base['genomes'][gid]['path']
        anomaly_file = output + '.anomalies.jsonl'
        if anomaly_file not in base['artifacts'].get(gid, {}):
            raise ValueError('anomaly export not covered by production receipt')
        predictions = {mode: {**runner.parse_integrated(r['genomes'][gid]['path'], 'go'),
                               **runner.parse_integrated(r['genomes'][gid]['path'], 'ec')}
                       for mode, r in runs.items()}
        records = [json.loads(line) for line in Path(anomaly_file).read_text().splitlines() if line.strip()]
        coverage.append({'genome_id': gid, 'split': split, 'n_anomalies': len(records),
                         'extractor_limit': 150, 'possibly_capped': len(records) >= 150})
        for ordinal, record in enumerate(records):
            reject_truth(record)
            a = copy.deepcopy(record)
            # Ignore source identifiers; the new ID has no exposed protein/function fields.
            a.pop('anomaly_id')
            aid = 'A_' + digest([gid, ordinal, a])[:24]
            a['anomaly_id'] = aid
            a['genome_id'] = gid
            target = a['target_function_id'].removeprefix('EC:')
            ids = [c['protein_id'] for c in a['candidates']]
            if len(set(ids)) != len(ids):
                raise ValueError('duplicate candidate protein')
            for c in a['candidates']:
                pid = c['protein_id']
                c['function_id'] = a['target_function_id']
                # A suggestion score is not the baseline's actual called probability.
                c['baseline_probability'] = predictions['stage1'].get((pid, target), 0.0)
                c['model_probabilities'] = {
                    mode: predictions[mode].get((pid, target))
                    for mode in ('learned_logistic', 'learned_gbdt')}
                c['raw_sequence_score'] = max((float(e['raw_score']) for e in c.get('supporting_claims', [])
                                                if e.get('raw_score') is not None), default=None)
                for value in c['model_probabilities'].values():
                    if value is not None:
                        probability(value)
                probability(c['baseline_probability'])
            cases.append({'split': split, 'genome_id': gid, 'anomaly': a})
    return {'schema': 'gspa.agent-public/2', 'cases': cases, 'coverage': coverage,
            'source_sha256': digest(source), 'source_receipts': {
                name: runner.hash_record(path) for name, path in source['runs'].items()},
            'shared_production_wall_seconds': sum(r['runtime_sec'] for r in runs.values()),
            'limitations': source.get('limitations', [])}


def abstain(a, **extra):
    return {'anomaly_id': a['anomaly_id'], 'action': 'ABSTAIN', 'protein_id': '',
            'function_id': a['target_function_id'], 'confidence': 0.0, **extra}


def rank(mode, a):
    if mode not in MODES:
        raise ValueError('unknown comparator')
    if mode == 'abstain_only':
        return abstain(a)
    def score(c):
        if mode == 'rank_by_raw_seq':
            return c.get('raw_sequence_score')
        if mode == 'rank_by_stage1_posterior':
            return c['baseline_probability']
        if mode == 'rank_by_context_support':
            # Frozen single-feature context comparator; no truth-derived weights.
            return float(bool(c.get('in_pathway_operon')))
        return c['model_probabilities'].get(mode)
    eligible = [(score(c), c['protein_id']) for c in a['candidates'] if score(c) is not None]
    if not eligible:
        return abstain(a, reason='no measured score')
    if any(not math.isfinite(s) for s, _ in eligible):
        raise ValueError('nonfinite candidate score')
    score_value, pid = sorted(eligible, key=lambda p: (-p[0], p[1]))[0]
    return {'anomaly_id': a['anomaly_id'], 'action': 'PROMOTE', 'protein_id': pid,
            'function_id': a['target_function_id'], 'confidence': score_value}


def shared_guardrail(a, proposal, guardrail):
    """Same final constraints for every arm; rejects rather than inventing another action."""
    if proposal.get('anomaly_id') != a['anomaly_id']:
        raise ValueError('resolution has wrong case ID')
    r = copy.deepcopy(proposal)
    action = r.get('action')
    if action == 'ABSTAIN':
        return r
    target, fid = a['target_function_id'], r.get('function_id')
    reason = None
    if action not in ('PROMOTE', 'DEMOTE', 'RELAX'):
        reason = 'invalid action'
    elif r.get('protein_id') not in {c['protein_id'] for c in a['candidates']}:
        reason = 'outside candidate set'
    elif action != 'RELAX' and fid != target:
        reason = 'changed target function'
    elif action == 'RELAX' and (not target.startswith('GO:') or fid == target or fid not in guardrail.go_ancestors(target)):
        reason = 'not a strict verified ancestor'
    elif action in ('PROMOTE', 'RELAX') and fid.startswith('GO:') and not guardrail.check_go_term(fid, a.get('organism_taxon', ''))['valid']:
        reason = 'symbolic constraint'
    if reason:
        return abstain(a, guardrail_rejected=reason, proposal=r)
    confidence = r.get('confidence')
    if isinstance(confidence, bool) or not isinstance(confidence, (int, float)) or not math.isfinite(confidence):
        return abstain(a, error='invalid decision score', proposal=r)
    return r


def public_cases(public):
    reject_truth(public)
    if public.get('schema') != 'gspa.agent-public/2':
        raise ValueError('expected public bundle v2')
    cases = public['cases']
    ids = [c['anomaly']['anomaly_id'] for c in cases]
    if len(set(ids)) != len(ids):
        raise ValueError('duplicate case ID')
    return cases


def suggest(public, mode, guardrail, solve=None):
    cases = public_cases(public)
    start = time.monotonic()
    rows = []
    for case in cases:
        a = copy.deepcopy(case['anomaly'])
        proposal = solve(a) if solve else rank(mode, a)
        rows.append(shared_guardrail(a, proposal, guardrail))
    return {'schema': 'gspa.agent-proposals/2', 'public_sha256': digest(public), 'mode': mode,
            'guardrail_policy': 'shared-reject-v2', 'resolutions': rows,
            'decision_wall_seconds': time.monotonic() - start,
            'shared_production_wall_seconds': public['shared_production_wall_seconds']}


def suggest_agent(public, guardrail, model='qwen3.8-27b', api_base=None, max_turns=2, workers=4, seed=42):
    """Bounded production agent, same public evidence and final constraints as all rankers.

    Counts HTTP attempts including retries, tokens actually reported by the
    endpoint, and errors. No fabricated zero cost when usage is absent.
    """
    cases = public_cases(public)
    if max_turns < 2 or workers < 1:
        raise ValueError('at least two turns needed for the full-evidence agent prompt')
    base, resolved = agent.resolve_endpoint_and_model(model, api_base)
    key = agent.load_api_key()
    original_call, original_open = agent.call_chat_completion, agent.urllib.request.urlopen
    stats, lock = Counter(), threading.Lock()
    def counted_open(*args, **kwargs):
        with lock:
            stats['http_attempts'] += 1
        return original_open(*args, **kwargs)
    def counted_call(**kwargs):
        try:
            response = original_call(**kwargs)
        except Exception:
            with lock:
                stats['failed_calls'] += 1
            raise
        with lock:
            stats['completed_calls'] += 1
            usage = response.get('usage')
            if not usage or any(usage.get(k) is None for k in ('prompt_tokens', 'completion_tokens')):
                stats['calls_without_token_usage'] += 1
            else:
                for k in ('prompt_tokens', 'completion_tokens'):
                    stats[k] += usage[k]
        return response
    def solve(case):
        a = copy.deepcopy(case['anomaly'])
        raw = agent.solve_single_anomaly(a, guardrail, set(), base, key, resolved, model,
                                        False, max_turns, seed, 0.0)
        return shared_guardrail(a, raw, guardrail)
    start = time.monotonic()
    try:
        agent.call_chat_completion, agent.urllib.request.urlopen = counted_call, counted_open
        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
            rows = list(pool.map(solve, cases))
    finally:
        agent.call_chat_completion, agent.urllib.request.urlopen = original_call, original_open
    return {'schema': 'gspa.agent-proposals/2', 'public_sha256': digest(public), 'mode': 'qwen_agent',
            'guardrail_policy': 'shared-reject-v2', 'resolutions': rows,
            'decision_wall_seconds': time.monotonic() - start,
            'shared_production_wall_seconds': public['shared_production_wall_seconds'],
            'endpoint': base, 'model': resolved, 'max_turns': max_turns, 'workers': workers,
            'seed': seed, 'temperature': 0.0, 'api_usage': dict(stats),
            'billing_cost': None,
            'usage_limitation': 'Tokens for failed calls or responses without usage are unknown; no billing rate assumed.'}


def checked_predictions(public, proposals):
    cases = public_cases(public)
    if proposals.get('public_sha256') != digest(public) or proposals.get('guardrail_policy') != 'shared-reject-v2':
        raise ValueError('proposals do not match public evidence/constraint policy')
    rows = proposals['resolutions']
    by_id = {r['anomaly_id']: r for r in rows}
    if len(by_id) != len(rows) or set(by_id) != {c['anomaly']['anomaly_id'] for c in cases}:
        raise ValueError('missing, duplicate or extra resolutions; failures must be explicit')
    return by_id


def truth_map(records, cases, split):
    allowed = {c['genome_id'] for c in cases if c['split'] == split}
    truth = {}
    for row in records:
        if row['genome_id'] not in allowed:
            raise ValueError(f'{split} truth file contains another split or unknown genome')
        if row['status'] not in ('positive', 'negative', 'unknown') or not row.get('source'):
            raise ValueError('truth requires explicit status and traceable source')
        key = row['genome_id'], row['protein_id'], row['function_id']
        if key in truth and truth[key] != row['status']:
            raise ValueError('conflicting truth')
        truth[key] = row['status']
    return truth


def outcome(case, resolution, truth):
    """Classify actual baseline changes. Absence of a label stays unknown."""
    if resolution.get('error') or resolution.get('guardrail_status') == 'ERROR':
        return 'error'
    if resolution['action'] == 'ABSTAIN':
        return 'abstained'
    a = case['anomaly']
    c = next(c for c in a['candidates'] if c['protein_id'] == resolution['protein_id'])
    called = c['baseline_probability'] >= .5
    action = resolution['action']
    if (action == 'PROMOTE' and called) or (action in ('DEMOTE', 'RELAX') and not called):
        return 'unchanged'
    target = truth.get((case['genome_id'], c['protein_id'], a['target_function_id']), 'unknown')
    if action == 'RELAX':
        broad = truth.get((case['genome_id'], c['protein_id'], resolution['function_id']), 'unknown')
        if target == 'positive' or broad == 'negative':
            return 'harmful'
        return 'correct' if target == 'negative' and broad == 'positive' else 'unknown'
    if target == 'unknown':
        return 'unknown'
    correct = (target == 'positive') == (action == 'PROMOTE')
    return 'correct' if correct else 'harmful'


def gated(a, proposal, threshold):
    if proposal['action'] != 'ABSTAIN' and (threshold is None or proposal['confidence'] < threshold):
        return abstain(a, reason='calibration policy abstention')
    return proposal


def calibrate(public, proposals, labels, precision=.9, min_decisions=10):
    if not 0 < precision <= 1 or min_decisions < 1:
        raise ValueError('invalid calibration requirements')
    cases = public_cases(public)
    rows = checked_predictions(public, proposals)
    truth = truth_map(labels, cases, 'calibration')
    calibration = [c for c in cases if c['split'] == 'calibration']
    scores = sorted({rows[c['anomaly']['anomaly_id']]['confidence'] for c in calibration
                     if rows[c['anomaly']['anomaly_id']]['action'] != 'ABSTAIN'})
    choices = []
    for threshold in scores:
        counts = Counter(outcome(c, gated(c['anomaly'], rows[c['anomaly']['anomaly_id']], threshold), truth) for c in calibration)
        n = sum(counts[k] for k in ('correct', 'harmful', 'unknown'))
        p = counts['correct'] / n if n else 0
        # Unknown active decisions count against demonstrable precision.
        if n >= min_decisions and p >= precision:
            choices.append((counts['correct'], -counts['harmful'], threshold, dict(counts)))
    best = max(choices, default=None)
    return {'schema': 'gspa.agent-policy/2', 'mode': proposals['mode'],
            'public_sha256': digest(public), 'calibration_truth_sha256': digest(labels),
            'calibration_proposals_sha256': digest([rows[c['anomaly']['anomaly_id']] for c in calibration]),
            'target_precision': precision, 'min_decisions': min_decisions,
            'threshold': best[2] if best else None,
            'calibration_counts': best[3] if best else {},
            'status': 'selected_on_calibration' if best else 'no_supported_threshold_abstain',
            'selection_rule': 'max_correct_then_min_harm_then_highest_threshold; unknowns count against precision',
            'limitation': 'empirical calibration precision, not a population guarantee'}


def evaluate(public, proposals, policy, labels):
    if policy.get('schema') != 'gspa.agent-policy/2' or policy['public_sha256'] != digest(public) or policy['mode'] != proposals['mode']:
        raise ValueError('missing or mismatched frozen calibration policy')
    cases, rows = public_cases(public), checked_predictions(public, proposals)
    cal_rows = [rows[c['anomaly']['anomaly_id']] for c in cases if c['split'] == 'calibration']
    if policy['calibration_proposals_sha256'] != digest(cal_rows):
        raise ValueError('calibration predictions changed after policy freeze')
    truth = truth_map(labels, cases, 'test')
    details, counts, by_genome = [], Counter(), defaultdict(Counter)
    for case in cases:
        if case['split'] != 'test':
            continue
        a = case['anomaly']
        raw = rows[a['anomaly_id']]
        resolution = gated(a, raw, policy['threshold'])
        status = outcome(case, resolution, truth)
        # Do not conceal a failed inference call behind a confidence gate.
        if raw.get('error') or raw.get('guardrail_status') == 'ERROR':
            status = 'error'
        counts[status] += 1
        by_genome[case['genome_id']][status] += 1
        details.append({'anomaly_id': a['anomaly_id'], 'genome_id': case['genome_id'],
                        'outcome': status, 'resolution': resolution})
    active = sum(counts[k] for k in ('correct', 'harmful', 'unknown'))
    known = counts['correct'] + counts['harmful']
    return {'schema': 'gspa.agent-evaluation/2', 'mode': proposals['mode'],
            'public_sha256': digest(public), 'policy_sha256': digest(policy),
            'test_truth_sha256': digest(labels), 'counts': dict(counts), 'n_cases': len(details),
            'known_decision_precision': counts['correct'] / known if known else None,
            'demonstrated_precision_lower_bound': counts['correct'] / active if active else None,
            'net_known_gain': counts['correct'] - counts['harmful'],
            'by_genome': {g: dict(c) for g, c in by_genome.items()}, 'details': details,
            'decision_wall_seconds': proposals['decision_wall_seconds'],
            'shared_production_wall_seconds': public['shared_production_wall_seconds'],
            'api_usage': proposals.get('api_usage'),
            'cost_note': 'All arms receive both learned scores; shared production includes both model runs. Decision time is additional.',
            'limitations': public.get('limitations', []) + ['Unknown truth is not a negative.',
                'Repeated case decisions are clustered within genome; counts are case-level.']}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest='phase', required=True)
    prep = sub.add_parser('prepare')
    prep.add_argument('--sources', required=True, type=Path)
    for phase in ('suggest', 'agent', 'calibrate', 'evaluate'):
        cmd = sub.add_parser(phase)
        cmd.add_argument('--public', required=True, type=Path)
        if phase in ('suggest', 'agent'):
            cmd.add_argument('--go-obo', required=True, type=Path)
            if phase == 'suggest':
                cmd.add_argument('--mode', required=True, choices=MODES)
            else:
                cmd.add_argument('--model', default='qwen3.8-27b')
                cmd.add_argument('--api-base')
                cmd.add_argument('--max-turns', type=int, default=2)
                cmd.add_argument('--workers', type=int, default=4)
        else:
            cmd.add_argument('--proposals', required=True, type=Path)
            cmd.add_argument('--truth', required=True, type=Path)
        if phase == 'calibrate':
            cmd.add_argument('--precision', type=float, default=.9)
            cmd.add_argument('--min-decisions', type=int, default=10)
        if phase == 'evaluate':
            cmd.add_argument('--policy', required=True, type=Path)
    for cmd in sub.choices.values():
        cmd.add_argument('--out', required=True, type=Path)
    args = p.parse_args()
    code_identity = {str(path.relative_to(ROOT)): runner.sha256_file(str(path))
                     for path in (Path(__file__).resolve(), Path(agent.__file__).resolve())}
    if args.phase == 'prepare':
        result = prepare(read(args.sources))
    elif args.phase in ('suggest', 'agent'):
        guard = agent.SymbolicGuardrail(ROOT / 'gspa-core/src/main/resources', args.go_obo)
        if not guard.resources_loaded or not guard.go_parents:
            raise ValueError('complete symbolic resources and GO hierarchy are required')
        result = (suggest(read(args.public), args.mode, guard) if args.phase == 'suggest' else
                  suggest_agent(read(args.public), guard, args.model, args.api_base, args.max_turns, args.workers))
        resources = ROOT / 'gspa-core/src/main/resources'
        result['decision_code_sha256'] = code_identity
        result['constraint_inputs'] = [runner.hash_record(str(args.go_obo))] + [
            runner.hash_record(str(path)) for directory in ('taxon-constraints', 'anatomical-constraints')
            for path in sorted((resources / directory).glob('*.tsv'))]
    elif args.phase == 'calibrate':
        result = calibrate(read(args.public), read(args.proposals), read(args.truth), args.precision, args.min_decisions)
    else:
        result = evaluate(read(args.public), read(args.proposals), read(args.policy), read(args.truth))
    if code_identity != {str(path.relative_to(ROOT)): runner.sha256_file(str(path))
                         for path in (Path(__file__).resolve(), Path(agent.__file__).resolve())}:
        raise RuntimeError('decision code changed during execution')
    write(args.out, result)


if __name__ == '__main__':
    main()
