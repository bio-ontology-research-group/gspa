#!/usr/bin/env python3
"""Controlled false enzyme assignments, and their effect on the measured coculture predictions.

Part B of the frozen plan ``tasks/v3-row8-row9-absence-crossfeeding-2026-11-01.md``.

``inject``  rewrites gene-protein-reaction rules of a model: it draws ``round(rate * n)`` of the
            reactions that have a GPR and replaces each drawn reaction's whole rule with one single
            gene drawn from the model's other gene products. The draw is a SHA-256 stream keyed by
            (model SHA-256, rate, seed), so it is reproducible from the receipt alone. Rate 0.00 is
            the control: it changes nothing.
``run``     reruns the frozen coculture plans through the production community engine with the
            perturbed rules. Only the members' knockout sets change, which is the effect under test.
``score``   compares the perturbed predictions with the baseline predictions and with the measured
            endpoints, using the frozen endpoint definitions.

The GPR rewrite is deliberately crude and explicit: it represents an enzyme assignment given to the
wrong gene. It is not a model of any particular annotation pipeline's error.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
import platform
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'scripts'))
sys.path.insert(0, str(ROOT / 'benchmark' / 'v3'))

import metabolic_experiments as ME  # noqa: E402
import community_metabolism as C  # noqa: E402
import run_nutrient_community as RC  # noqa: E402

SCHEMA = 'gspa.false_assignment_sensitivity/1'
POLICIES = RC.PROTOCOL['policies']
SOLVE_ORDER = RC.PROTOCOL['solve_order']


def sha(path) -> str:
    digest = hashlib.sha256()
    with open(path, 'rb') as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b''):
            digest.update(chunk)
    return digest.hexdigest()


def digest(value) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, allow_nan=False).encode()).hexdigest()


def rows(path):
    with open(path) as stream:
        for line in stream:
            if line.strip():
                yield json.loads(line)


def write_json(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')


def code_hashes() -> dict:
    return {str(p.relative_to(ROOT)): sha(p) for p in
            (Path(__file__), ROOT / 'scripts/metabolic_experiments.py',
             ROOT / 'scripts/community_metabolism.py', ROOT / 'benchmark/v3/run_nutrient_community.py')}


# --------------------------------------------------------------------------- injection

def _stream(key: str):
    """Deterministic uniform [0,1) stream from a SHA-256 counter keyed by the run identity."""
    counter = 0
    while True:
        block = hashlib.sha256(f'{key}|{counter}'.encode()).digest()
        for offset in range(0, 32, 8):
            yield int.from_bytes(block[offset:offset + 8], 'big') / 2 ** 64
        counter += 1


def _shuffle(items: list, draws) -> list:
    items = list(items)
    for i in range(len(items) - 1, 0, -1):
        j = int(next(draws) * (i + 1))
        items[i], items[j] = items[j], items[i]
    return items


def inject(model, rate: float, seed: int, mode: str = 'uniform', focal_genes: list | None = None) -> tuple[dict, dict]:
    """Return (perturbed GPR mapping, receipt). The model object is not modified.

    ``uniform`` (the frozen scheme) draws from every reaction that has a GPR. ``targeted``
    (deviation D1 of the plan, declared before it ran) draws only from the reactions whose rule
    mentions one of ``focal_genes`` - the genes the benchmark actually deletes - because a uniform
    draw almost never touches them and would measure a different question.
    """
    if not 0.0 <= rate <= 1.0:
        raise ValueError('rate must be in [0, 1]')
    if mode not in ('uniform', 'targeted'):
        raise ValueError('mode must be uniform or targeted')
    gpr = {rid: [list(clause) for clause in dnf] for rid, dnf in model.gpr.items() if dnf}
    genes = sorted({g for dnf in gpr.values() for clause in dnf for g in clause})
    if mode == 'targeted':
        focal = set(focal_genes or [])
        if not focal:
            raise ValueError('targeted mode requires focal genes')
        unknown = sorted(focal - set(genes))
        eligible = sorted(rid for rid, dnf in gpr.items()
                          if any(g in focal for clause in dnf for g in clause))
    else:
        focal, unknown, eligible = set(), [], sorted(gpr)
    count = int(round(rate * len(eligible)))
    draws = _stream(f'{model.sha256}|{mode}|{rate:.6f}|{seed}|{",".join(sorted(focal))}')
    order = _shuffle(eligible, draws)
    chosen = order[:count]
    perturbed = {rid: [list(c) for c in dnf] for rid, dnf in gpr.items()}
    changes = []
    for rid in sorted(chosen):
        own = {g for clause in gpr[rid] for g in clause}
        pool = [g for g in genes if g not in own and g not in focal]
        if not pool:
            changes.append({'reaction': rid, 'skipped': 'no other gene product available'})
            continue
        replacement = pool[int(next(draws) * len(pool))]
        perturbed[rid] = [[replacement]]
        changes.append({'reaction': rid, 'old_rule': gpr[rid], 'new_gene': replacement})
    receipt = {'schema': SCHEMA + ':injection', 'model_sha256': model.sha256, 'rate': rate, 'seed': seed,
               'mode': mode, 'focal_genes': sorted(focal), 'focal_genes_absent_from_rules': unknown,
               'eligible_reactions': len(eligible), 'reactions_with_gpr': len(gpr),
               'gene_products_in_rules': len(genes),
               'reactions_perturbed': sum(1 for c in changes if 'new_gene' in c),
               'draw': 'SHA-256 counter stream keyed by (model sha256, mode, rate, seed, focal genes); '
                       'Fisher-Yates over the eligible reactions, then one replacement gene per reaction; '
                       'replacement genes exclude the focal genes so the wrong assignment is not self-cancelling',
               'rule_semantics': "the drawn reaction's whole rule is replaced by a single wrong gene",
               'changes': changes}
    return perturbed, receipt


def apply_gpr(model, perturbed: dict):
    """A copy of the model whose GPR mapping is the perturbed one (stoichiometry untouched)."""
    clone = copy.copy(model)
    merged = dict(model.gpr)
    merged.update(perturbed)
    clone.gpr = merged
    return clone


# --------------------------------------------------------------------------- run

def member_knockouts(model, genes: list) -> list:
    return sorted(model.blocked_by_genes(set(genes)))


def simulate(plan: dict, context: dict, model, model_path: str, baseline_model) -> dict:
    """The frozen plan through the production community engine, with perturbed knockout sets."""
    members, knockout_changes = [], []
    for member in plan['members']:
        blocked = member_knockouts(model, member['genes'])
        frozen = sorted(member['blocked_reactions'])
        if blocked != frozen:
            knockout_changes.append({'member': member['id'], 'frozen': frozen, 'perturbed': blocked,
                                     'added': sorted(set(blocked) - set(frozen)),
                                     'removed': sorted(set(frozen) - set(blocked))})
        members.append({**member, 'blocked_reactions': blocked})
    perturbed_plan = {**plan, 'members': members}
    manifest = RC.make_manifest(perturbed_plan, context, model_path)
    cfg = {**C.DEFAULT_CONFIG, **manifest['config']}
    orgs = [C.build_organism(spec['id'], baseline_model, spec, context['pool_metabolites'], cfg)
            for spec in manifest['organisms']]
    analysis = C.Analysis(manifest, orgs, cfg)
    for index, organism in enumerate(orgs):
        analysis.growth(f'monoculture_{index}', 'monoculture', [organism.id], analysis.lp([organism.id]))
    analysis.growth('permissive', 'community', analysis.order, analysis.lp(analysis.order))
    native_orgs = []
    for spec in manifest['organisms']:
        native = copy.deepcopy(spec)
        for rid, port in native['ports'].items():
            if plan['supply_exchanges'].get(rid, 0.) <= 0:
                port.pop('open_uptake')
        native_orgs.append(C.build_organism(spec['id'], baseline_model, native, context['pool_metabolites'], cfg))
    native_lp = C.build_lp(native_orgs, analysis.ratio, analysis.supply, analysis.removal)
    analysis.growth('native_unsupplied', 'community', analysis.order, native_lp)
    closures = [(o.id, rid, 'uptake') for o in orgs for rid in o.ports
                if plan['supply_exchanges'].get(rid, 0.) <= 0]
    analysis.growth('no_unsupplied_uptake', 'community', analysis.order, analysis.lp(analysis.order, closures))
    for record in analysis.records:
        record.pop('_x', None)
        record.pop('_lp', None)
    return {'case_id': plan['case_id'], 'arm': plan['arm'], 'plan_sha256': digest(plan),
            'manifest_sha256': digest(manifest),
            'predictions': {r['id']: RC.prediction(r) for r in analysis.records},
            'knockout_changes': knockout_changes,
            'solver_calls': len(analysis.solver.log)}


def focal_genes_for(protocol: dict, species: str) -> list:
    """The canonical deletion genes of this species in the frozen registration."""
    return sorted({g for gene in protocol['models'][species].get('genes', []) for g in [gene['canonical']]})


def run(plan_dir: Path, out: Path, models: dict, rate: float, seed: int,
        shard_index: int = 0, shard_count: int = 1, mode: str = 'uniform') -> dict:
    if shard_count < 1 or not 0 <= shard_index < shard_count:
        raise ValueError('invalid shard index/count')
    protocol = json.loads((plan_dir / 'protocol.json').read_text())
    if protocol['plans_sha256'] != sha(plan_dir / 'plans.jsonl'):
        raise ValueError('frozen plans changed')
    contexts = json.loads((plan_dir / 'contexts.json').read_text())
    loaded, perturbed, injection = {}, {}, {}
    for species, path in models.items():
        model = ME.load_sbml_model(Path(path))
        if model.sha256 != protocol['models'][species]['sha256']:
            raise ValueError(f'{species}: model differs from the frozen registration')
        focal = focal_genes_for(protocol, species) if mode == 'targeted' else None
        mapping, receipt = inject(model, rate, seed, mode, focal)
        loaded[species] = model
        perturbed[species] = apply_gpr(model, mapping)
        injection[species] = receipt
    out.mkdir(parents=True, exist_ok=False)
    started, counts, changed_cases = time.time(), Counter(), 0
    with (out / 'predictions.jsonl').open('w') as stream:
        for index, plan in enumerate(rows(plan_dir / 'plans.jsonl')):
            if index % shard_count != shard_index:
                continue
            species = plan['species']
            try:
                row = simulate(plan, contexts[species], perturbed[species], models[species], loaded[species])
            except Exception as exc:
                row = {'case_id': plan['case_id'], 'arm': plan['arm'], 'plan_sha256': digest(plan),
                       'predictions': {k: None for k in SOLVE_ORDER},
                       'failure': f'{type(exc).__name__}: {exc}'}
            counts[str(row['predictions'].get('permissive'))] += 1
            changed_cases += bool(row.get('knockout_changes'))
            stream.write(json.dumps(row, sort_keys=True, allow_nan=False) + '\n')
            if (index + 1) % 500 == 0:
                stream.flush()
                print(json.dumps({'rate': rate, 'seed': seed, 'shard': shard_index,
                                  'seen': index + 1, 'elapsed_s': round(time.time() - started, 1)}), flush=True)
    receipt = {'schema': SCHEMA + ':run', 'rate': rate, 'seed': seed, 'mode': mode,
               'protocol_sha256': sha(plan_dir / 'protocol.json'),
               'plans_sha256': protocol['plans_sha256'],
               'contexts_sha256': sha(plan_dir / 'contexts.json'),
               'predictions_sha256': sha(out / 'predictions.jsonl'),
               'injection': injection, 'permissive_counts': dict(counts),
               'cases_with_changed_knockouts': changed_cases,
               'implementation': code_hashes(), 'python': sys.version, 'platform': platform.platform(),
               'shard_index': shard_index, 'shard_count': shard_count,
               'elapsed_s': time.time() - started, 'reads_measurements': False}
    write_json(out / 'run.json', receipt)
    return receipt


def merge(plan_dir: Path, shards: list, out: Path) -> dict:
    protocol = json.loads((plan_dir / 'protocol.json').read_text())
    receipts, seen = {}, {}
    rate = seed = mode = None
    for directory in shards:
        receipt = json.loads((directory / 'run.json').read_text())
        if (receipt['plans_sha256'] != protocol['plans_sha256']
                or receipt['implementation'] != code_hashes()
                or receipt['predictions_sha256'] != sha(directory / 'predictions.jsonl')):
            raise ValueError(f'{directory}: shard provenance mismatch')
        if receipt['shard_count'] != len(shards) or receipt['shard_index'] in receipts:
            raise ValueError(f'{directory}: duplicate or incomplete shard set')
        if mode is None:
            rate, seed, mode = receipt['rate'], receipt['seed'], receipt['mode']
        elif (rate, seed, mode) != (receipt['rate'], receipt['seed'], receipt['mode']):
            raise ValueError('shards mix different perturbations')
        receipts[receipt['shard_index']] = {'dir': str(directory), 'run_sha256': sha(directory / 'run.json')}
        for row in rows(directory / 'predictions.jsonl'):
            key = (row['case_id'], row['arm'])
            if key in seen:
                raise ValueError(f'duplicate prediction for {key}')
            seen[key] = row
    if set(receipts) != set(range(len(shards))):
        raise ValueError('incomplete shard set')
    order = [(p['case_id'], p['arm']) for p in rows(plan_dir / 'plans.jsonl')]
    if set(order) != set(seen):
        raise ValueError('shards do not cover every plan exactly once')
    out.mkdir(parents=True, exist_ok=False)
    counts, changed = Counter(), 0
    with (out / 'predictions.jsonl').open('w') as stream:
        for key in order:
            counts[str(seen[key]['predictions'].get('permissive'))] += 1
            changed += bool(seen[key].get('knockout_changes'))
            stream.write(json.dumps(seen[key], sort_keys=True, allow_nan=False) + '\n')
    receipt = {'schema': SCHEMA + ':run', 'rate': rate, 'seed': seed, 'mode': mode,
               'protocol_sha256': sha(plan_dir / 'protocol.json'), 'plans_sha256': protocol['plans_sha256'],
               'predictions_sha256': sha(out / 'predictions.jsonl'), 'permissive_counts': dict(counts),
               'cases_with_changed_knockouts': changed, 'merged_shards': receipts,
               'injection': {k: v for k, v in json.loads((shards[0] / 'run.json').read_text())['injection'].items()},
               'implementation': code_hashes(), 'plans': len(order), 'reads_measurements': False}
    write_json(out / 'run.json', receipt)
    return receipt


# --------------------------------------------------------------------------- score

def wilson(successes: int, total: int) -> list:
    if total == 0:
        return [None, None]
    z, p, n = 1.959963984540054, successes / total, total
    centre = (p + z * z / (2 * n)) / (1 + z * z / n)
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return [max(0.0, centre - half), min(1.0, centre + half)]


def accuracy(pairs):
    scorable = [(o, p) for o, p in pairs if o is not None and p is not None]
    correct = sum(1 for o, p in scorable if o == p)
    return {'n_total': len(pairs), 'n_scorable': len(scorable), 'n_correct': correct,
            'accuracy': (correct / len(scorable)) if scorable else None,
            'accuracy_ci': wilson(correct, len(scorable))}


def endpoints(observations: Path):
    obs = list(rows(observations))
    by_id = {o['case_id']: o for o in obs}
    mono = {(o['species'], o['environment'], o['strains'][0]): o for o in obs if o['kind'] == 'monoculture'}
    return by_id, mono


def score(plan_dir: Path, baseline: Path, arms: list, observations: Path, out: Path) -> dict:
    protocol = json.loads((plan_dir / 'protocol.json').read_text())
    plans = {(p['case_id'], p['arm']): p for p in rows(plan_dir / 'plans.jsonl')}
    by_id, mono = endpoints(observations)
    base = {(r['case_id'], r['arm']): r['predictions'] for r in rows(baseline)}
    if set(base) != set(plans):
        raise ValueError('baseline predictions do not cover the frozen plans')

    summary, per_arm, control = [], [], None
    for directory in arms:
        receipt = json.loads((directory / 'run.json').read_text())
        if receipt['predictions_sha256'] != sha(directory / 'predictions.jsonl'):
            raise ValueError(f'{directory}: predictions changed since the run')
        rate, seed, mode = receipt['rate'], receipt['seed'], receipt.get('mode', 'uniform')
        strata = defaultdict(list)
        flips = Counter()
        identical = 0
        total = 0
        changed_members = 0
        for row in rows(directory / 'predictions.jsonl'):
            key = (row['case_id'], row['arm'])
            plan = plans[key]
            total += 1
            changed_members += bool(row.get('knockout_changes'))
            if row['predictions'] == base[key]:
                identical += 1
            else:
                for policy in SOLVE_ORDER:
                    if row['predictions'].get(policy) != base[key].get(policy):
                        flips[policy] += 1
            observed = by_id.get(plan['case_id'], {}).get('mean_od600')
            observed_bulk = None if observed is None else observed > 0.02
            member_obs = [mono.get((plan['species'], plan['environment'], m['id'].split(':')[-1]), {}).get('mean_od600')
                          for m in plan['members']]
            member_bulk = [None if v is None else v > 0.02 for v in member_obs]
            observed_rescue = RC.rescue(observed_bulk, *member_bulk)
            scope = 'primary' if plan['primary_environment'] else 'no_nitrogen_control'
            group = (plan['species'], plan['arm'], scope)
            strata[group + ('permissive', 'bulk_growth')].append((observed_bulk, row['predictions']['permissive']))
            strata[group + ('permissive', 'paired_rescue_proxy')].append(
                (observed_rescue, RC.rescue(row['predictions']['permissive'],
                                            row['predictions']['monoculture_0'],
                                            row['predictions']['monoculture_1'])))
        entry = {'rate': rate, 'seed': seed, 'mode': mode, 'directory': str(directory),
                 'plans': total, 'identical_to_baseline': identical,
                 'plans_with_changed_knockouts': changed_members,
                 'prediction_flips_by_policy': dict(flips),
                 'reactions_perturbed': {s: r['reactions_perturbed'] for s, r in receipt['injection'].items()},
                 'strata': [{'species': k[0], 'arm': k[1], 'scope': k[2], 'method': k[3],
                             'endpoint': k[4], **accuracy(v)} for k, v in sorted(strata.items())]}
        summary.append(entry)
        if rate == 0.0 and mode == 'uniform':
            control = {'plans': total, 'identical_to_baseline': identical,
                       'reproduces_baseline_exactly': identical == total}
    report = {'schema': SCHEMA + ':score',
              'interpretation': 'Sensitivity of the frozen coculture endpoint predictions to controlled wrong '
                                'gene-reaction assignments. A prediction that does not change is not evidence '
                                'that the assignment is harmless biologically; it means the endpoint under this '
                                'design does not depend on it.',
              'inputs': {'protocol': sha(plan_dir / 'protocol.json'), 'plans': protocol['plans_sha256'],
                         'baseline_predictions': sha(baseline), 'observations': sha(observations)},
              'implementation': code_hashes(), 'rate_zero_control': control, 'arms': summary}
    write_json(out, report)
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    runner = sub.add_parser('run')
    for key in ('plan-dir', 'out', 'ec-model', 'bs-model'):
        runner.add_argument('--' + key, type=Path, required=True)
    runner.add_argument('--rate', type=float, required=True)
    runner.add_argument('--seed', type=int, required=True)
    runner.add_argument('--shard-index', type=int, default=0)
    runner.add_argument('--shard-count', type=int, default=1)
    runner.add_argument('--mode', choices=('uniform', 'targeted'), default='uniform')
    merger = sub.add_parser('merge')
    merger.add_argument('--plan-dir', type=Path, required=True)
    merger.add_argument('--shards', type=Path, nargs='+', required=True)
    merger.add_argument('--out', type=Path, required=True)
    scorer = sub.add_parser('score')
    scorer.add_argument('--plan-dir', type=Path, required=True)
    scorer.add_argument('--baseline', type=Path, required=True)
    scorer.add_argument('--arms', type=Path, nargs='+', required=True)
    scorer.add_argument('--observations', type=Path, required=True)
    scorer.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.command == 'run':
        receipt = run(args.plan_dir, args.out, {'Ec': str(args.ec_model), 'Bs': str(args.bs_model)},
                      args.rate, args.seed, args.shard_index, args.shard_count, args.mode)
        print(json.dumps({k: receipt[k] for k in ('rate', 'seed', 'mode', 'permissive_counts',
                                                  'cases_with_changed_knockouts', 'elapsed_s')}))
    elif args.command == 'merge':
        receipt = merge(args.plan_dir, list(args.shards), args.out)
        print(json.dumps({k: receipt[k] for k in ('rate', 'seed', 'mode', 'plans', 'permissive_counts',
                                                  'cases_with_changed_knockouts')}))
    else:
        report = score(args.plan_dir, args.baseline, list(args.arms), args.observations, args.out)
        print(json.dumps({'arms': len(report['arms']), 'rate_zero_control': report['rate_zero_control']}))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
