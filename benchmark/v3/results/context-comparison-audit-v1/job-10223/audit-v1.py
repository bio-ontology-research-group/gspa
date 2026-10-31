#!/usr/bin/env python3
"""Independent audit of the completed context-mode comparison; never fits models.

Reconstructs all scored supervision from production TSVs, checks every registered
file, reproduces frozen thresholds and statistics, and evaluates all experimental
EC-reference loci at those unchanged thresholds. The EC expansion is post-hoc,
descriptive, and positive-reference-only; unconfirmed calls are not proven errors.
Run the complete audit on a compute node, not on a login node or laptop.
"""
import argparse
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
import csv
import hashlib
import json
import math
from pathlib import Path
import random

ARMS = tuple(f'{mode}/{kind}' for mode in ('UNIVERSAL', 'PATHWAY', 'ORGANISM')
             for kind in ('rule', 'logistic', 'gbdt'))
PREFIX = Path('/home/leechuck/.cache/gspa-v3-data')
LABELS = '/home/leechuck/.cache/gspa-v3-data/reviewed-thrs-locus-v1/labels.tsv'
SPLITS = '/home/leechuck/Public/software/gspa/benchmark/v3/results/training-cohort/splits.tsv'
PROTOCOL = {
    'schema': 'gspa.context-comparison-independent-audit/1',
    'thresholds': 'Recompute original max-F1 cuts from calibration rows only, ties prefer higher cut. Never retune on EC subset or test.',
    'population': 'All 6706 registered supervision rows. Additional EC panel contains every ECO:0000269 positive EC assertion in the four test genomes, including unemitted functions/proteins.',
    'ec_policy': 'Exact EC identifiers; no ontology remapping. Unconfirmed predictions are reference nonmatches, not experimentally proven negatives.',
    'ranking': 'Descending production probability, then lexical EC. Report all ties at the top, and hits at 1,3,5 over the full emitted candidate list.',
    'bootstrap': {'seed': 20261027, 'replicates': 10000, 'unit': 'paired genome'},
    'interpretation': 'Post-hoc descriptive audit of a previously examined cohort; no untouched-test, biological-calibration or causal context-improvement claim.',
}


def read(path):
    return json.loads(Path(path).read_text(), parse_constant=lambda x: (_ for _ in ()).throw(ValueError(x)))


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def probability(value):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or not 0 <= value <= 1:
        raise ValueError('invalid probability')
    return value


def equivalent(a, b):
    if isinstance(a, dict):
        return isinstance(b, dict) and a.keys() == b.keys() and all(equivalent(a[k], b[k]) for k in a)
    if isinstance(a, list):
        return isinstance(b, list) and len(a) == len(b) and all(equivalent(x, y) for x, y in zip(a, b))
    if isinstance(a, (int, float)) and not isinstance(a, bool):
        return isinstance(b, (int, float)) and math.isclose(a, b, rel_tol=1e-9, abs_tol=1e-11)
    return a == b


def metrics(rows, arm, cut):
    counts = Counter(); bins = defaultdict(list); squared = []
    for r in rows:
        p = probability(r['scores'][arm]); y = r['label']; yes = p >= cut
        if type(y) is not int or y not in (0, 1):
            raise ValueError('invalid label')
        counts[('t' if y else 'f') + 'p' if yes else ('f' if y else 't') + 'n'] += 1
        bins[min(9, int(p * 10))].append(p - y); squared.append((p - y) ** 2)
    tp, fp, fn, tn = (counts[k] for k in ('tp', 'fp', 'fn', 'tn'))
    return dict(n=len(rows), positives=tp + fn, tp=tp, fp=fp, fn=fn, tn=tn,
        precision=tp / (tp + fp) if tp + fp else None, recall=tp / (tp + fn) if tp + fn else None,
        f1=2 * tp / (2 * tp + fp + fn) if 2 * tp + fp + fn else None,
        brier=math.fsum(squared) / len(rows) if rows else None,
        ece_10bin=math.fsum(abs(math.fsum(v)) for v in bins.values()) / len(rows) if rows else None, threshold=cut)


def threshold(rows, arm):
    if not rows or any(r['split'] != 'calibration' for r in rows):
        raise ValueError('calibration rows required')
    if {r['label'] for r in rows} != {0, 1}:
        raise ValueError('both calibration classes required')
    candidates = sorted({probability(r['scores'][arm]) for r in rows}, reverse=True)
    f1, cut = max((metrics(rows, arm, cut)['f1'], cut) for cut in candidates)
    return dict(rule='max_f1', value=cut, selected_on='calibration', advisory_only=True, calibration_f1=f1)


def bootstrap(values):
    values = [values[k] for k in sorted(values)]
    n = len(values); point = sum(values) / n if n else None
    if n < 2:
        return dict(estimate=point, ci_low=None, ci_high=None, n_genomes=n, n_boot=0)
    rng = random.Random(PROTOCOL['bootstrap']['seed'])
    draws = sorted(sum(values[rng.randrange(n)] for _ in range(n)) / n
                   for _ in range(PROTOCOL['bootstrap']['replicates']))
    return dict(estimate=point, ci_low=draws[math.floor(.025 * (len(draws) - 1))],
                ci_high=draws[math.ceil(.975 * (len(draws) - 1))], n_genomes=n, n_boot=len(draws))


def paired(rows, arm, base, cuts):
    per = {}
    for gid in sorted({r['genome_id'] for r in rows}):
        subset = [r for r in rows if r['genome_id'] == gid]
        changes = [(r['scores'][arm] >= cuts[arm], r['scores'][base] >= cuts[base], r['label']) for r in subset]
        per[gid] = dict(rows=len(subset), brier_delta=metrics(subset, arm, cuts[arm])['brier'] - metrics(subset, base, cuts[base])['brier'],
            reference_corrected=sum(a == y and b != y for a, b, y in changes),
            reference_harmed=sum(a != y and b == y for a, b, y in changes), changed_calls=sum(a != b for a, b, y in changes))
    corrected = sum(v['reference_corrected'] for v in per.values()); harmed = sum(v['reference_harmed'] for v in per.values())
    return dict(reference_corrected=corrected, reference_harmed=harmed, net_reference_corrections=corrected - harmed,
        changed_calls=sum(v['changed_calls'] for v in per.values()), per_genome=per,
        mean_genome_brier_delta=bootstrap({g: v['brier_delta'] for g, v in per.items()}))


def ec_panel(references, scores, cut):
    """All reference loci remain in the denominator, including empty predictions."""
    records = []
    for (gid, protein), truth in sorted(references.items()):
        candidates = scores.get((gid, protein), {})
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
    return dict(proteins=len(records), reference_pairs=positive, emitted_reference_pairs=sum(len(r['candidate_reference_matches']) for r in records),
        proteins_with_candidates=sum(r['candidates'] > 0 for r in records), reference_matches=tp, unconfirmed_calls=unconfirmed,
        missed_reference=positive - tp, reference_precision=tp / (tp + unconfirmed) if tp + unconfirmed else None,
        reference_recall=tp / positive if positive else None, reference_micro_f1=2 * tp / (tp + unconfirmed + positive) if tp + unconfirmed + positive else None,
        top1_matches=sum(r['top1_reference_match'] for r in records), top3_matches=sum(r['top3_reference_match'] for r in records),
        top5_matches=sum(r['top5_reference_match'] for r in records), top1_tied=sum(len(r['top_tie']) > 1 for r in records),
        threshold=cut, records=records)


class Inputs:
    def __init__(self, package, base, final):
        self.package, self.base, self.final = map(Path, (package, base, final)); self.used = set()

    def locate(self, path):
        original = Path(path)
        if not original.is_absolute() or '..' in original.parts:
            raise ValueError('absolute registered path required')
        if original.is_relative_to(PREFIX):
            rel = original.relative_to(PREFIX)
            for root in (self.final, self.base):
                p = root / rel
                if p.is_file():
                    return p
        p = self.package / str(original).lstrip('/')
        if not p.is_file():
            raise FileNotFoundError(path)
        return p

    def consume(self, path):
        self.used.add(str(path)); return self.locate(path)


def audit(args):
    out = args.out; out.mkdir(parents=True, exist_ok=False); save(out / 'protocol.json', PROTOCOL)
    source = Inputs(args.package, args.base, args.final)
    analysis = str(PREFIX / 'observed-gap-context-model-analysis-v1')
    verification = read(source.consume(analysis + '/verification.json'))
    if verification.get('status') != 'verified':
        raise ValueError('analysis not verified')
    bindings = dict(verification['inputs'])
    for name, digest in verification['outputs'].items():
        if name in bindings and bindings[name] != digest:
            raise ValueError('conflicting registered hashes')
        bindings[name] = digest
    # Include the seven supplemental supervision files and their original parent
    # registration, rather than treating an unhashed registry as prior evidence.
    parent = str(PREFIX / 'observed-gap-models-v1/registration.json')
    parent_record = read(source.consume(parent))
    supplement = read(args.supplement)
    for name, digest in supplement['inputs'].items():
        if parent_record['inputs'].get(name) != digest:
            raise ValueError('supplement differs from original parent registration: ' + name)
        if name in bindings and bindings[name] != digest:
            raise ValueError('supplement hash conflict')
        bindings[name] = digest
    def check(item):
        path, expected = item; actual_path = source.locate(path); actual = sha(actual_path)
        return dict(path=path, physical_path=str(actual_path), bytes=actual_path.stat().st_size,
                    expected=expected, actual=actual, agrees=expected == actual)
    checked = []
    with (out / 'file-checks.jsonl').open('x') as log, ThreadPoolExecutor(max_workers=args.workers) as pool:
        for n, row in enumerate(pool.map(check, sorted(bindings.items())), 1):
            log.write(json.dumps(row) + '\n'); log.flush(); checked.append(row)
            if not row['agrees']:
                raise ValueError('changed registered file: ' + row['path'])
            if n % 100 == 0:
                print(json.dumps({'hashed_files': n, 'total': len(bindings)}), flush=True)
    protocol = read(source.consume(analysis + '/protocol.json')); original = read(source.consume(analysis + '/results.json'))['results']
    thresholds = read(source.consume(analysis + '/thresholds.json'))
    scored = [json.loads(line) for line in source.consume(analysis + '/scored-supervision.jsonl').read_text().splitlines()]
    with source.consume(SPLITS).open() as f:
        test_genomes = {r['genome_id'] for r in csv.DictReader(f, delimiter='\t') if r['split'] == 'test'}
    references = defaultdict(set)
    with source.consume(LABELS).open() as f:
        for row in csv.DictReader(f, delimiter='\t'):
            if row['genome_id'] in test_genomes and row['function_id'].startswith('EC:') and row['evidence'] == 'ECO:0000269' and row['label'] == '1':
                references[(row['genome_id'], row['protein_id'])].add(row['function_id'][3:])
    report = dict(protocol=PROTOCOL, source_files=len(checked), source_bytes=sum(r['bytes'] for r in checked),
                  original_registry_bound=True, root_verification_sha256=sha(source.locate(analysis + '/verification.json')), regimes={})
    for transfer, plan in protocol['plans'].items():
        table = [json.loads(line) for line in source.consume(plan['table']).read_text().splitlines()]
        for mode in ('PATHWAY', 'ORGANISM'):
            other_path = plan['table'].replace('/UNIVERSAL/table/', '/' + mode + '/table/')
            other = [json.loads(line) for line in source.consume(other_path).read_text().splitlines()]
            strip = lambda r: {k: v for k, v in r.items() if k not in ('features', 'cohort_context_mode')}
            if [strip(r) for r in table] != [strip(r) for r in other]:
                raise ValueError('context modes changed supervision metadata')
        rows = [r for r in scored if r['transfer'] == transfer]
        if len(table) != len(rows) or len({r['example_id'] for r in rows}) != len(rows):
            raise ValueError('supervision population differs')
        by_id = {r['example_id']: r for r in rows}; wanted = defaultdict(dict)
        for t in table:
            r = by_id[t['example_id']]
            if any(r[k] != t[k] for k in ('genome_id', 'split', 'label', 'label_basis')):
                raise ValueError('supervision metadata differs')
            r['function_type'] = t['function_type']
            term = t['function_id'].removeprefix('EC:') if t['function_type'] == 'EC' else t['function_id']
            wanted[t['genome_id']][(t['protein_id'], t['function_type'], term)] = r
        all_ec = {}; predictions_checked = 0
        for arm, paths in plan['paths'].items():
            if arm not in ARMS:
                raise ValueError('unexpected arm')
            all_ec[arm] = defaultdict(dict)
            for gid, path in paths.items():
                seen = set(); matched = set()
                with source.consume(path).open() as f:
                    for prediction in csv.DictReader(f, delimiter='\t'):
                        kind = prediction['type']; term = prediction['function_id'].removeprefix('EC:') if kind == 'EC' else prediction['function_id']
                        key = (prediction['protein_id'], kind, term)
                        if key in seen:
                            raise ValueError('duplicate prediction')
                        seen.add(key); p = probability(float(prediction['posterior_prob'])); predictions_checked += 1
                        if key in wanted[gid]:
                            if p != wanted[gid][key]['scores'][arm]:
                                raise ValueError('scored supervision differs from production TSV')
                            matched.add(key)
                        if kind == 'EC' and (gid, key[0]) in references:
                            all_ec[arm][(gid, key[0])][term] = p
                if matched != wanted[gid].keys():
                    raise ValueError('missing scored supervision')
        cal = [r for r in rows if r['split'] == 'calibration']; test = [r for r in rows if r['split'] == 'test']
        if {r['genome_id'] for r in test} != test_genomes or set(plan['paths']) != set(ARMS):
            raise ValueError('incomplete test genomes or arms')
        for arm in ARMS:
            if not equivalent(threshold(cal, arm), thresholds[transfer][arm]):
                raise ValueError('calibration threshold mismatch: ' + arm)
        result = dict(supervision_rows=len(rows), production_rows_checked=predictions_checked, test_rows=len(test),
                      test_label_basis=dict(Counter(r['label_basis'] for r in test)), policies={})
        if set(original[transfer]['policies']) != {'calibration_selected', 'fixed_0.5'}:
            raise ValueError('unexpected threshold policies')
        for policy, original_policy in original[transfer]['policies'].items():
            cuts = {a: thresholds[transfer][a]['value'] if policy == 'calibration_selected' else .5 for a in ARMS}
            for arm in ARMS:
                if not equivalent(metrics(test, arm, cuts[arm]), original_policy['arms'][arm]):
                    raise ValueError('reported metrics differ')
            for pair, old in original_policy['comparisons'].items():
                arm, base = pair.split('_vs_')
                if not equivalent(paired(test, arm, base, cuts), old):
                    raise ValueError('reported paired statistics differ: ' + pair)
            result['policies'][policy] = dict(
                supervised_by_type={kind: {a: metrics([r for r in test if r['function_type'] == kind], a, cuts[a]) for a in ARMS} for kind in ('EC', 'GO')},
                full_ec_reference={a: ec_panel(references, all_ec[a], cuts[a]) for a in ARMS})
        report['regimes'][transfer] = result
    if source.used - bindings.keys() != {analysis + '/verification.json'}:
        raise ValueError('unregistered consumed inputs')
    for path in sorted(source.used):
        if path in bindings and sha(source.locate(path)) != bindings[path]:
            raise ValueError('consumed input changed during audit')
    report['consumed_files_rechecked'] = len(source.used & bindings.keys())
    report['status'] = 'verified'
    save(out / 'report.json', report)
    print(json.dumps({k: v for k, v in report.items() if k != 'regimes'}), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('package', 'base', 'final', 'supplement', 'out'):
        p.add_argument('--' + name, type=Path, required=True)
    p.add_argument('--workers', type=int, default=4)
    audit(p.parse_args())
