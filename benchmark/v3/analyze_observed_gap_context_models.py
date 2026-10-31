"""Matched final-score comparison of observed-gap context modes.

All frozen supervision rows, including unsupported learned rule fallback, remain
included. Calibration selects thresholds before test metrics. This reused cohort
measures reference support, not biological accuracy; unconfirmed is not a proven
negative. Sparse gap-linked supervision and reference overlap prevent claims of
missing-enzyme recovery or unseen-family generalization. No model is promoted.
"""
import argparse
from collections import Counter
import csv
import json
from pathlib import Path

import analyze_observed_gap_models as A
import run_observed_gap_context_modes as R
import ec_reference_panel as EC
F = R.F
KINDS = ('rule',) + tuple(F.KINDS)
ARMS = tuple(f'{mode}/{kind}' for mode in R.MODES for kind in KINDS)
POLICY = dict(schema='gspa.observed-gap-context-comparison/2', scope=__doc__,
    arms=list(ARMS), thresholds=A.POLICY['thresholds'],
    population=A.POLICY['population'], bootstrap=A.POLICY['bootstrap'],
    uncertainty=A.POLICY['uncertainty'], full_ec_reference=EC.POLICY,
    comparisons='Each restricted context against UNIVERSAL with the same scorer; each learned scorer against rules within the same context.')


def attach_scores(rows, by_arm):
    if set(by_arm) != set(ARMS):
        raise ValueError('missing or unexpected comparison arm')
    result = []; seen = set()
    for row in rows:
        if row['split'] not in F.T.SPLITS:
            raise ValueError('unknown split')
        if type(row['label']) is not int or row['label'] not in (0, 1):
            raise ValueError('invalid label')
        gid = row['genome_id']; key = A.function_key(row)
        if (gid, key) in seen:
            raise ValueError('duplicate supervised function identity')
        seen.add((gid, key)); scores = {}
        for arm in ARMS:
            if gid not in by_arm[arm] or key not in by_arm[arm][gid]:
                raise ValueError('missing supervised prediction')
            scores[arm] = A.probability(by_arm[arm][gid][key])
        result.append(dict(example_id=row['example_id'], genome_id=gid,
            split=row['split'], label=row['label'], label_basis=row['label_basis'],
            scores=scores))
    return result


def select_thresholds(rows):
    if not rows or any(r['split'] != 'calibration' for r in rows):
        raise ValueError('threshold selection requires calibration rows only')
    if {r['label'] for r in rows} != {0, 1}:
        raise ValueError('both calibration classes required')
    return {arm: F.T.choose_threshold(F.T.np.asarray([r['scores'][arm] for r in rows]),
        F.T.np.asarray([r['label'] for r in rows]), 'max_f1') for arm in ARMS}


def evaluate_rows(rows, thresholds, n_boot=10000):
    if not rows or any(r['split'] != 'test' for r in rows):
        raise ValueError('test rows only')
    pairs = [(f'{mode}/{kind}', f'UNIVERSAL/{kind}')
        for mode in R.MODES if mode != 'UNIVERSAL' for kind in KINDS]
    pairs += [(f'{mode}/{kind}', f'{mode}/rule')
        for mode in R.MODES for kind in F.KINDS]
    result = dict(rows=len(rows), genomes=sorted({r['genome_id'] for r in rows}),
        label_basis=dict(Counter(r['label_basis'] for r in rows)), policies={})
    for name, cuts in [('calibration_selected', {a: thresholds[a]['value'] for a in ARMS}),
                       ('fixed_0.5', {a: .5 for a in ARMS})]:
        result['policies'][name] = dict(
            arms={a: A.metrics(rows, a, cuts[a]) for a in ARMS},
            comparisons={f'{a}_vs_{b}': A.paired(rows, a, b, cuts, n_boot) for a, b in pairs})
    return result


def register(training, replay, out):
    protocol = R.verify(training)
    production = json.loads((replay / 'verification.json').read_text())
    if production.get('status') != 'verified' or production.get('schema') != 'gspa.observed-gap-context-model-verification/1':
        raise ValueError('complete verified context production replay required')
    F.check_hashes(production['inputs']); F.check_hashes(production['outputs'])
    frozen = json.loads((training / 'training.json').read_text())
    exported = json.loads((training / 'export.json').read_text())
    if frozen.get('status') != 'verified' or exported.get('status') != 'verified':
        raise ValueError('complete frozen models and exports required')
    F.check_hashes(frozen['outputs']); F.check_hashes(exported['outputs'])
    inputs = F.hashes(R.sources() + [replay / 'verification.json', training / 'training.json', training / 'export.json'])
    inputs.update(production['inputs']); inputs.update(production['outputs'])
    inputs.update(frozen['outputs']); inputs.update(exported['outputs'])
    # R.verify checks the supervision files against the original fit registration.
    # Bind the actual reference population directly in this analysis as well.
    supervision = {t: protocol['recipes'][t]['supervision'] for t in R.G.MODES}
    reference_inputs = {t: {k: str(Path(s[k]).resolve()) for k in ('labels', 'splits')}
                        for t, s in supervision.items()}
    inputs.update(R.supervision_inputs(protocol['recipes'],
        json.loads((training / 'registration.json').read_text())['inputs']))
    plans = {}
    for transfer in R.G.MODES:
        table = training / transfer / 'UNIVERSAL/table/training_table.jsonl'
        for mode in R.MODES:
            R.compare_supervision(table, training / transfer / mode / 'table/training_table.jsonl', mode)
        paths = {}
        for mode in R.MODES:
            with Path(protocol['plans'][transfer][mode]['manifest']).open() as f:
                paths[f'{mode}/rule'] = {r['genome_id']: r['out_path'] for r in csv.DictReader(f, delimiter='\t')}
            for kind in F.KINDS:
                paths[f'{mode}/{kind}'] = {g: str(replay / transfer / mode / kind / 'inference' / (g + '.integrated.tsv')) for g in protocol['genomes']}
        if any(set(p) != set(protocol['genomes']) for p in paths.values()):
            raise ValueError('incomplete genome coverage')
        # Every consumed file must be among the verified artifacts, not merely
        # present beside them in a similarly named directory.
        consumed = [table] + [Path(p) for arm in paths.values() for p in arm.values()]
        if any(str(p.resolve()) not in inputs for p in consumed):
            raise ValueError('unverified analysis input')
        plans[transfer] = dict(table=str(table), paths=paths,
            reference_inputs=reference_inputs[transfer], genomes=protocol['genomes'])
    out.mkdir(parents=True, exist_ok=False)
    A.save(out / 'protocol.json', dict(POLICY, plans=plans))
    A.save(out / 'registration.json', dict(status='registered',
        registered_at=F.datetime.now(F.timezone.utc).isoformat(),
        protocol_sha256=F.T.sha256_file(out / 'protocol.json'), inputs=inputs))


def run(out):
    reg = json.loads((out / 'registration.json').read_text()); F.check_hashes(reg['inputs'])
    if F.T.sha256_file(out / 'protocol.json') != reg['protocol_sha256']:
        raise ValueError('comparison protocol changed')
    protocol = json.loads((out / 'protocol.json').read_text())
    if any(protocol[k] != v for k, v in POLICY.items()):
        raise ValueError('analysis policy changed')
    if (out / 'thresholds.json').exists():
        raise FileExistsError('refusing analysis overwrite')
    prepared = {}; thresholds = {}; panels = {}
    for transfer, plan in protocol['plans'].items():
        with Path(plan['table']).open() as f:
            rows = [json.loads(s) for s in f]
        scores = {a: {g: A.read_scores(Path(p)) for g, p in paths.items()} for a, paths in plan['paths'].items()}
        prepared[transfer] = attach_scores(rows, scores)
        thresholds[transfer] = select_thresholds([r for r in prepared[transfer] if r['split'] == 'calibration'])
        refs = plan['reference_inputs']
        if any(str(Path(p).resolve()) not in reg['inputs'] for p in refs.values()):
            raise ValueError('unregistered reference population input')
        references, test_genomes = EC.load_references(Path(refs['labels']), Path(refs['splits']), plan['genomes'])
        if {r['genome_id'] for r in prepared[transfer] if r['split'] == 'test'} != test_genomes:
            raise ValueError('supervision/reference test genomes differ')
        # Keep only the reference loci's EC rows, without filtering their scores.
        panels[transfer] = (references, test_genomes, {a: {g: {k: v for k, v in values.items()
            if k[1] == 'EC' and (g, k[0]) in references} for g, values in by_genome.items() if g in test_genomes}
            for a, by_genome in scores.items()})
    A.save(out / 'thresholds.json', thresholds)
    threshold_hash = F.T.sha256_file(out / 'thresholds.json')
    results = {t: evaluate_rows([r for r in rows if r['split'] == 'test'], thresholds[t]) for t, rows in prepared.items()}
    for transfer, (references, test_genomes, scores) in panels.items():
        for policy, summary in results[transfer]['policies'].items():
            cuts = {a: thresholds[transfer][a]['value'] if policy == 'calibration_selected' else .5 for a in ARMS}
            summary['full_ec_reference'] = {a: EC.evaluate(references, scores[a], test_genomes, cuts[a]) for a in ARMS}
    A.save(out / 'results.json', dict(scope=__doc__, results=results))
    with (out / 'scored-supervision.jsonl').open('x') as f:
        for transfer, rows in prepared.items():
            for row in rows:
                f.write(json.dumps(dict(transfer=transfer, **row), allow_nan=False) + '\n')
    F.check_hashes(reg['inputs'])
    if F.T.sha256_file(out / 'thresholds.json') != threshold_hash or F.T.sha256_file(out / 'protocol.json') != reg['protocol_sha256']:
        raise ValueError('analysis policy changed during evaluation')
    A.save(out / 'verification.json', dict(status='verified', inputs=reg['inputs'], outputs=F.hashes([p for p in out.iterdir() if p.is_file()])))
    return results


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['register', 'run'])
    for key in ('training', 'replay', 'out'):
        parser.add_argument('--' + key, type=Path, required=key == 'out')
    args = parser.parse_args()
    if args.action == 'register':
        if not args.training or not args.replay:
            parser.error('training and replay required')
        register(args.training.resolve(), args.replay.resolve(), args.out.resolve())
    else:
        run(args.out.resolve())
