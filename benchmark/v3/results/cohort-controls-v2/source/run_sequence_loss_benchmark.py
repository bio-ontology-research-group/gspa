#!/usr/bin/env python3
"""Matched intact / coordinate-loss / nucleotide-loss production cohort comparison.

Calibration policies are frozen before test scoring. Positive-only reference
annotations support recovery/loss counts; unmatched calls remain unverified,
not experimentally false. This development cohort is not a new blind test.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import prepare_sequence_loss as P
import run_cohort_perturbations as R
import verified_batch as V
import v3_metrics as M

ARMS = ('intact', 'coordinate_loss', 'sequence_loss')


def read_splits(path):
    with path.open() as fh:
        rows = list(csv.DictReader(fh, delimiter='\t'))
    if len({r['genome_id'] for r in rows}) != len(rows):
        raise ValueError('duplicate split genome')
    return {r['genome_id']: r['split'] for r in rows}


def run(inputs, work, model):
    preparation = P.verify(inputs)
    splits = read_splits(R.RES / 'splits.tsv')
    if set(preparation['genomes']) != set(splits):
        raise ValueError('prepared and split genomes differ')
    R.PERT_DIR = work
    common = dict(genomes=sorted(splits), learned_model_path=model, rlc_suggester=False,
                  sources=(Path(__file__), Path(P.__file__), inputs / 'preparation.json'))
    for arm in ARMS:
        R.run_perturbation_batch(
            arm, layout_dir=R.CACHE / 'layout' if arm == 'intact' else inputs / 'layout',
            cohort_meta_path=R.CACHE / 'cohort_inputs/cohort_metadata.json' if arm == 'intact'
            else inputs / f'metadata_{arm}.json',
            cohort_og_path=inputs / 'orthogroups.tsv' if arm == 'sequence_loss'
            else R.CACHE / 'cohort_inputs/cohort_orthogroups.tsv',
            claims_override={g: inputs / 'claims' / f'{g}.claims.jsonl' for g in splits}
            if arm == 'sequence_loss' else None, **common)
    P.verify(inputs)


def verified_predictions(work, genomes, *, arms=ARMS, extra_genomes=None):
    predictions, hashes = {}, {}
    for arm in arms:
        out = work / arm / 'out'
        receipt = json.loads((out / V.RECEIPT).read_text())
        spec = receipt['spec']
        if V.file_hashes(spec['inputs']) != spec['inputs'] or not V.reusable(out, spec):
            raise ValueError(f'{arm}: execution receipt is stale or invalid')
        expected = set(genomes) | set((extra_genomes or {}).get(arm, []))
        if set(spec['outputs']) != expected:
            raise ValueError(f'{arm}: output genome coverage differs')
        hashes[arm] = P.digest(out / V.RECEIPT)
        predictions[arm] = read_predictions(out, genomes)
    return predictions, hashes


def read_predictions(out, genomes):
    result = {}
    for gid in genomes:
        result[gid] = {'GO': {}, 'EC': {}}
        with (out / f'{gid}.integrated.tsv').open() as fh:
            for row in csv.DictReader(fh, delimiter='\t'):
                kind = row['type']
                if kind not in ('GO', 'EC'):
                    continue
                value = row['function_id']
                value = M.normalize_ec(value) if kind == 'EC' else M.normalize_go(value)
                score = float(row['posterior_prob'])
                if not math.isfinite(score) or not 0 <= score <= 1:
                    raise ValueError(f'{gid}: invalid score')
                key = (row['protein_id'], value)
                if key in result[gid][kind]:
                    raise ValueError(f'{gid}: duplicate candidate {key}')
                result[gid][kind][key] = score
    return result


def read_truth(labels, genomes, go):
    truth = {g: {'GO': set(), 'EC': set()} for g in genomes}
    with labels.open() as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            if row['genome_id'] not in truth or row['tier'] != 'experimental':
                continue
            fid = row['function_id']
            kind = 'GO' if fid.startswith('GO:') else 'EC'
            fid = M.normalize_ec(fid) if kind == 'EC' else M.normalize_go(fid)
            truth[row['genome_id']][kind].add((row['protein_id'], fid))
    for g in truth:
        truth[g]['GO'] = go.propagate_truth(truth[g]['GO'])
    return truth


def calls(prediction, threshold):
    return set() if threshold is None else {k for k, s in prediction.items() if s >= threshold}


def select_policy(predictions, truth, calibration_genomes, *, minimum=10, precision=0.9):
    policies = {}
    for kind in ('GO', 'EC'):
        choices = []
        for tau in M.GRID:
            if tau <= 0:
                continue
            supported = total = 0
            for g in calibration_genomes:
                labelled = {p for p, _ in truth[g][kind]}
                active = {k for k in calls(predictions[g][kind], tau) if k[0] in labelled}
                total += len(active)
                supported += len(active & truth[g][kind])
            if total >= minimum and supported / total >= precision:
                choices.append((supported, supported / total, tau, total))
        if choices:
            n, fraction, tau, total = max(choices)
            policies[kind] = {'threshold': tau, 'supported': n, 'called': total,
                              'supported_fraction': fraction}
        else:
            policies[kind] = {'threshold': None, 'reason': 'insufficient calibration support; abstain'}
    return {'minimum_calls': minimum, 'required_reference_support': precision,
            'selection': 'maximize supported calls, then supported fraction, then threshold; intact calibration only',
            'kinds': policies}


def propagated(predictions, go):
    for by_genome in predictions.values():
        for by_kind in by_genome.values():
            by_kind['GO'] = go.propagate_scores(by_kind['GO'])
    return predictions


def calibrate(inputs, work, policy_path, labels, splits_path, go_path):
    if policy_path.exists():
        raise FileExistsError('policy is frozen; refusing overwrite')
    P.verify(inputs)
    splits = read_splits(splits_path)
    predictions, receipts = verified_predictions(work, splits)
    go = M.load_obo(str(go_path))
    propagated(predictions, go)
    calibration = [g for g, s in splits.items() if s == 'calibration']
    if not calibration:
        raise ValueError('no calibration genomes')
    truth = read_truth(labels, calibration, go)
    policy = select_policy(predictions['intact'], truth, calibration)
    policy.update(schema='gspa.sequence-loss-policy/1', calibration_genomes=calibration,
                  receipts=receipts, inputs={str(p.resolve()): P.digest(p) for p in
                  (inputs / 'preparation.json', labels, splits_path, go_path, Path(__file__))})
    P.write_json(policy_path, policy)
    return policy


def evaluate(inputs, work, policy_path, labels, splits_path, go_path, n_boot=2000):
    preparation = P.verify(inputs)
    splits = read_splits(splits_path)
    predictions, receipts = verified_predictions(work, splits)
    policy = json.loads(policy_path.read_text())
    if receipts != policy['receipts'] or any(P.digest(p) != h for p, h in policy['inputs'].items()):
        raise ValueError('inputs changed after calibration')
    if policy['inputs'] != {str(p.resolve()): P.digest(p) for p in
                           (inputs / 'preparation.json', labels, splits_path, go_path, Path(__file__))}:
        raise ValueError('evaluation inputs differ from frozen calibration provenance')
    test = [g for g, s in splits.items() if s == 'test']
    if set(test) & set(policy['calibration_genomes']):
        raise ValueError('test/calibration overlap')
    go = M.load_obo(str(go_path))
    truth = read_truth(labels, test, go)
    propagated(predictions, go)
    # A removed protein must never reappear, even if an identical sequence survives elsewhere.
    for g in splits:
        retained = set(preparation['genomes'][g]['retained_proteins'])
        for kind in ('GO', 'EC'):
            if any(p not in retained for p, _ in predictions['sequence_loss'][g][kind]):
                raise ValueError(f'{g}: removed protein reappeared in production output')
    return summarize(predictions, truth, preparation, policy, test, n_boot=n_boot) | {
        'schema': 'gspa.sequence-loss-evaluation/1', 'receipts': receipts,
        'policy_sha256': P.digest(policy_path), 'preparation_sha256': P.digest(inputs / 'preparation.json'),
        'limits': ['previously examined development cohort; not a new blind test',
                   'unmatched calls are unverified, not experimentally false; reference coverage is incomplete',
                   'GO ancestors are correlated reference pairs, not independent observations',
                   'no reaction gaps supplied: cohort stage runs; reaction-context stage is not requested',
                   'whole-CDS deletion only; boundary CDS excluded; unplaced proteins retained',
                   'known sequence loss is an intervention label, not inferred absence attribution']}


def summarize(predictions, truth, preparation, policy, test, *, n_boot=2000,
              arms=ARMS, scopes=('all_original', 'retained', 'removed')):
    result = {'test_genomes': test, 'scopes': {}}
    for scope in scopes:
        report = result['scopes'][scope] = {}
        for kind in ('GO', 'EC'):
            tau = policy['kinds'][kind]['threshold']
            by_arm = report[kind] = {}
            for arm in arms:
                vectors, top1, changes = {}, {}, {}
                for g in test:
                    keep = set(preparation['genomes'][g]['retained_proteins'])
                    t = {k for k in truth[g][kind] if scope == 'all_original' or
                         ((k[0] in keep) == (scope == 'retained'))}
                    if not t:
                        continue
                    labelled = {p for p, _ in t}
                    active = calls(predictions[arm][g][kind], tau)
                    baseline = calls(predictions['intact'][g][kind], tau)
                    known, unknown = len(active & t), len({k for k in active - t if k[0] in labelled})
                    vectors[g] = [known, unknown, len(t)]
                    changes[g] = [len((active - baseline) & t), len((baseline - active) & t)]
                    if kind == 'EC':
                        a = M.ec_top1_outcomes(predictions[arm][g][kind], t)
                        b = M.ec_top1_outcomes(predictions['intact'][g][kind], t)
                        top1[g] = M.top1_vector(a) + M.top1_vector(b)
                pooled = M.sum_vectors(vectors.values()) if vectors else [0, 0, 0]
                change = M.sum_vectors(changes.values()) if changes else [0, 0]
                recall = lambda v: v[0] / v[2] if v and v[2] else float('nan')
                by_arm[arm] = {'threshold': tau, 'supported_calls': pooled[0], 'unverified_calls': pooled[1],
                              'reference_pairs': pooled[2], 'new_supported_calls': change[0],
                              'lost_supported_calls': change[1],
                              'reference_recall': M.genome_bootstrap(vectors, recall, n_boot=n_boot),
                              'net_supported_calls': M.genome_bootstrap(changes, lambda v: v[0]-v[1] if v else float('nan'), n_boot=n_boot)}
                if top1:
                    rate = lambda v: v[2] / v[0] if v[0] else float('nan')
                    by_arm[arm]['ec_top1_reference_agreement'] = M.genome_bootstrap(top1, rate, n_boot=n_boot)
                    by_arm[arm]['ec_top1_delta_vs_intact'] = M.genome_bootstrap(top1, M.paired_stat(rate, 3), n_boot=n_boot)
                    pooled_top = M.sum_vectors(top1.values())
                    by_arm[arm]['ec_top1_counts'] = dict(zip(('labelled', 'called', 'reference_match'), pooled_top[:3]))
    return result


def clean_json(value):
    if isinstance(value, dict):
        return {k: clean_json(v) for k, v in value.items()}
    if isinstance(value, list):
        return [clean_json(v) for v in value]
    return None if isinstance(value, float) and not math.isfinite(value) else value


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('action', choices=('run', 'calibrate', 'evaluate'))
    ap.add_argument('--inputs', type=Path, required=True)
    ap.add_argument('--work', type=Path, required=True)
    ap.add_argument('--model', type=Path, default=R.RES / 'learned_models/model_logistic.json')
    ap.add_argument('--policy', type=Path)
    ap.add_argument('--out', type=Path)
    args = ap.parse_args()
    if args.action == 'run':
        run(args.inputs, args.work, args.model)
    else:
        if args.policy is None:
            ap.error('--policy is required')
        params = (args.inputs, args.work, args.policy, R.RES / 'labels.tsv', R.RES / 'splits.tsv', R.CACHE / 'go-basic.obo')
        if args.action == 'calibrate':
            calibrate(*params)
        else:
            if args.out is None or args.out.exists():
                ap.error('--out must name a new report')
            P.write_json(args.out, clean_json(evaluate(*params)))
