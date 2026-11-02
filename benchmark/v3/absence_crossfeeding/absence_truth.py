#!/usr/bin/env python3
"""Score the four absence explanations against known truth on the measured nutrient auxotrophs.

Three separate stages, as in the frozen plan
``tasks/v3-row8-row9-absence-crossfeeding-2026-11-01.md`` (Part A):

``build``  reads the frozen nutrient plan and the two curated models and writes one absence
           comparison per (deletion strain, environment, medium arm), plus wild-type controls.
           It reads no measurement.
``run``    executes each comparison through the production engine
           (``scripts/absence_experiments.run_absence_experiments``). It reads no measurement.
``score``  joins the measured endpoint OD600 and scores, per explanation, whether the known
           explanation is in the feasible set. Growth agreement is not re-scored here.

Truth (frozen before any run): for an unsupplemented deletion assay the focal reaction really is
absent, so the true explanation is ``pathway_loss``; for the matched supplemented assay the measured
growth is really due to the supplied amino acid, so the true explanation is ``external_supply``.
Wild-type assays have no true absence explanation and are negative controls.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
import multiprocessing as mp
import os
import platform
import sys
import time
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'scripts'))
sys.path.insert(0, str(ROOT / 'benchmark' / 'v3'))

import absence_experiments as AE  # noqa: E402
import metabolic_experiments as ME  # noqa: E402

SCHEMA = 'gspa.absence_truth/1'
HYPOTHESES = AE.HYPOTHESES
GROWTH_OD = 0.02  # the nutrient benchmark's own monoculture endpoint threshold
TRUTH = {'unsupplemented': 'pathway_loss', 'supplemented': 'external_supply'}


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
            (Path(__file__), ROOT / 'scripts/absence_experiments.py', ROOT / 'scripts/metabolic_experiments.py')}


# --------------------------------------------------------------------------- build

def medium_difference(unsupp: dict, supp: dict) -> tuple[str, float]:
    """The one exchange the matched supplemented assay adds, and its frozen rate."""
    a, b = unsupp['manifest']['baseline']['medium'], supp['manifest']['baseline']['medium']
    added = {k: v for k, v in b.items() if k not in a}
    raised = {k: v for k, v in b.items() if k in a and v > a[k]}
    if len(added) == 1 and not raised:
        return next(iter(added.items()))
    if not added and len(raised) == 1:
        return next(iter(raised.items()))
    raise ValueError(f'supplement is not a single exchange: added={sorted(added)} raised={sorted(raised)}')


def absence_manifest(model, base_plan: dict, model_path: str, focal: str,
                     supplement: str, rate: float) -> dict:
    """One four-scenario comparison: baseline = the unsupplemented medium, supply = the supplement."""
    manifest = copy.deepcopy(base_plan['manifest'])
    manifest.pop('context', None)  # the absence comparison has no candidate menu
    manifest['model'] = model_path
    exchanges = sorted(set(manifest['baseline']['exchanges']) | {supplement})
    manifest['baseline']['exchanges'] = exchanges  # declared, not supplied: uptake stays closed
    ports = sorted(AE.transport_reactions(model, supplement, set(
        manifest['baseline'].get('extracellular_compartments', ['e']))))
    manifest['absence_comparison'] = {
        'reaction': focal,
        'candidate_function_key': None,
        'assembly': {'status': 'COMPLETE', 'completeness_fraction': 1.0,
                     'source': 'curated published model of a laboratory strain; not an assembly assessment of this experiment'},
        'function_search': {'status': 'UNKNOWN',
                            'source': 'no whole-proteome functional search is part of this phenotype benchmark'},
        'external_supply': {'exchange': supplement, 'uptake_limit': float(rate),
                            'transport_reactions': ports,
                            'source': 'the measured supplement of the matched assay, at the frozen protocol rate'},
        'conditions': [{'id': 'unsupplemented', 'mods': []},
                       {'id': 'supplemented', 'mods': [{'op': 'open', 'exchange': supplement, 'rate': float(rate)}]}],
    }
    manifest['config'] = {**manifest.get('config', {}), 'strategy': 'fixed', 'panel': None,
                          'max_experiments': 2, 'max_solver_calls': 64}
    return manifest


def build(plan_dir: Path, models: dict, out: Path) -> dict:
    protocol = json.loads((plan_dir / 'protocol.json').read_text())
    plans = list(rows(plan_dir / 'plans.jsonl'))
    loaded = {species: ME.load_sbml_model(Path(path)) for species, path in models.items()}
    for species, model in loaded.items():
        if model.sha256 != protocol['models'][species]['sha256']:
            raise ValueError(f'{species}: model differs from the frozen registration')

    by_key = {}
    for plan in plans:
        by_key[(plan['species'], plan['strain_id'], plan['environment'], plan['arm'], plan['kind'])] = plan

    cases, skipped = [], []
    focal_by_strain: dict[tuple, list] = {}
    for (species, strain, environment, arm, kind), plan in sorted(by_key.items()):
        if kind != 'monoculture':
            continue
        supp_plan = by_key.get((species, strain, environment, arm, 'supplemented_monoculture'))
        if supp_plan is None:
            skipped.append({'reason': 'no matched supplemented assay', 'species': species,
                            'strain_id': strain, 'environment': environment, 'arm': arm})
            continue
        model = loaded[species]
        blocked = sorted(model.blocked_by_genes(set(plan['genes'])))
        focal_by_strain[(species, strain)] = blocked
        if not blocked:
            skipped.append({'reason': 'no blocked reaction (isozyme bypass): not constructible',
                            'species': species, 'strain_id': strain, 'environment': environment,
                            'arm': arm, 'genes': plan['genes']})
            continue
        scope = 'primary' if len(blocked) == 1 else 'sensitivity'
        supplement, rate = medium_difference(plan, supp_plan)
        for focal in blocked:
            manifest = absence_manifest(model, plan, models[species], focal, supplement, rate)
            cases.append({
                'case_id': f'ABS:{species}:{strain}:{environment}:{arm}:{focal}',
                'kind': 'auxotroph', 'scope': scope, 'species': species, 'strain_id': strain,
                'environment': environment, 'arm': arm, 'genes': plan['genes'],
                'focal_reaction': focal, 'blocked_reactions': blocked,
                'supplement_exchange': supplement, 'supplement_rate': float(rate),
                'assays': {'unsupplemented': plan['case_id'], 'supplemented': supp_plan['case_id']},
                'manifest': manifest, 'manifest_sha256': digest(manifest)})

    # Wild-type negative controls: the focal reaction of each primary auxotroph of the same species,
    # in the wild type where that reaction is really present. Truth: no absence explanation.
    for (species, strain, environment, arm, kind), plan in sorted(by_key.items()):
        if kind != 'wild_type':
            continue
        model = loaded[species]
        for (sp, dstrain), blocked in sorted(focal_by_strain.items()):
            if sp != species or len(blocked) != 1:
                continue
            donor = by_key.get((species, dstrain, environment, arm, 'monoculture'))
            donor_supp = by_key.get((species, dstrain, environment, arm, 'supplemented_monoculture'))
            if donor is None or donor_supp is None:
                continue
            supplement, rate = medium_difference(donor, donor_supp)
            manifest = absence_manifest(model, plan, models[species], blocked[0], supplement, rate)
            cases.append({
                'case_id': f'ABS-WT:{species}:{strain}:{environment}:{arm}:{blocked[0]}',
                'kind': 'wild_type_control', 'scope': 'control', 'species': species,
                'strain_id': strain, 'environment': environment, 'arm': arm, 'genes': [],
                'focal_reaction': blocked[0], 'blocked_reactions': blocked,
                'donor_strain_id': dstrain, 'supplement_exchange': supplement,
                'supplement_rate': float(rate),
                'assays': {'unsupplemented': plan['case_id'], 'supplemented': None},
                'manifest': manifest, 'manifest_sha256': digest(manifest)})

    out.mkdir(parents=True, exist_ok=False)
    with (out / 'cases.jsonl').open('w') as stream:
        for case in cases:
            stream.write(json.dumps(case, sort_keys=True, allow_nan=False) + '\n')
    counts = defaultdict(int)
    for case in cases:
        counts[f"{case['kind']}:{case['scope']}"] += 1
    receipt = {
        'schema': SCHEMA + ':build', 'truth': TRUTH,
        'scope_rules': 'one blocked reaction = primary; several = sensitivity, one comparison each; none = not constructible',
        'inputs': {'plans': sha(plan_dir / 'plans.jsonl'), 'protocol': sha(plan_dir / 'protocol.json'),
                   **{f'model_{s}': sha(p) for s, p in models.items()}},
        'model_sha256': {s: m.sha256 for s, m in loaded.items()},
        'implementation': code_hashes(), 'cases': len(cases), 'counts': dict(counts),
        'skipped': skipped, 'cases_sha256': sha(out / 'cases.jsonl'),
        'reads_measurements': False}
    write_json(out / 'build.json', receipt)
    return receipt


# --------------------------------------------------------------------------- run

_WORKER: dict = {}


def _init(models: dict) -> None:
    _WORKER.update({'models': {s: ME.load_sbml_model(Path(p)) for s, p in models.items()},
                    'paths': models})


def _execute(case: dict) -> dict:
    model = _WORKER['models'][case['species']]
    manifest = case['manifest']
    started = time.time()
    try:
        result = AE.run_absence_experiments(manifest, model)
        cells = {}
        for cell in result['cells']:
            entry = {'status': cell['status'], 'reason': cell.get('reason')}
            for perturbation, record in cell['phenotypes'].items():
                entry[perturbation] = {'status': record['status'], 'class': record['class'],
                                       'growth': record['growth']}
            cells[f"{cell['condition_id']}|{cell['hypothesis']}"] = entry
        return {'case_id': case['case_id'], 'status': 'completed',
                'manifest_sha256': case['manifest_sha256'],
                'input_model_sha256': result['input_model_sha256'],
                'problem_sha256': result['problem_sha256'],
                'gene_attribution': result['gene_attribution'],
                'selected_explanation': result['selected_explanation'],
                'solver_calls_used': len(result['solver_calls']),
                'visited_conditions': result['visited_conditions'],
                'unvisited_conditions': result['unvisited_conditions'],
                'comparisons': result['comparisons'], 'cells': cells,
                'wall_seconds': time.time() - started}
    except Exception as exc:  # recorded, never silently dropped
        return {'case_id': case['case_id'], 'status': 'failed',
                'manifest_sha256': case['manifest_sha256'],
                'failure': f'{type(exc).__name__}: {exc}', 'wall_seconds': time.time() - started}


def run(cases_dir: Path, out: Path, models: dict, workers: int, keep_manifests: Path | None = None,
        shard_index: int = 0, shard_count: int = 1) -> dict:
    """Execute the comparisons of one shard. Shards are independent processes; `merge` joins them."""
    if shard_count < 1 or not 0 <= shard_index < shard_count:
        raise ValueError('invalid shard index/count')
    build_receipt = json.loads((cases_dir / 'build.json').read_text())
    if build_receipt['cases_sha256'] != sha(cases_dir / 'cases.jsonl'):
        raise ValueError('cases changed since build')
    if build_receipt['implementation'] != code_hashes():
        raise ValueError('implementation changed since build')
    cases = [c for i, c in enumerate(rows(cases_dir / 'cases.jsonl')) if i % shard_count == shard_index]
    out.mkdir(parents=True, exist_ok=False)
    started = time.time()
    counts = defaultdict(int)
    with (out / 'results.jsonl').open('w') as stream:
        if workers > 1:
            with mp.get_context('fork').Pool(workers, initializer=_init, initargs=(models,)) as pool:
                for row in pool.imap(_execute, cases, chunksize=8):
                    counts[row['status']] += 1
                    stream.write(json.dumps(row, sort_keys=True, allow_nan=False) + '\n')
        else:
            _init(models)
            for index, case in enumerate(cases):
                row = _execute(case)
                counts[row['status']] += 1
                stream.write(json.dumps(row, sort_keys=True, allow_nan=False) + '\n')
                if (index + 1) % 100 == 0:
                    stream.flush()
                    print(json.dumps({'shard': shard_index, 'completed': index + 1,
                                      'of': len(cases), 'elapsed_s': round(time.time() - started, 1)}), flush=True)
    if keep_manifests is not None:
        keep_manifests.mkdir(parents=True, exist_ok=True)
        for case in cases:
            write_json(keep_manifests / f"{digest(case['manifest'])[:16]}.json", case['manifest'])
    receipt = {'schema': SCHEMA + ':run', 'build_sha256': sha(cases_dir / 'build.json'),
               'cases_sha256': build_receipt['cases_sha256'], 'counts': dict(counts),
               'results_sha256': sha(out / 'results.jsonl'), 'implementation': code_hashes(),
               'python': sys.version, 'platform': platform.platform(), 'workers': workers,
               'shard_index': shard_index, 'shard_count': shard_count,
               'elapsed_s': time.time() - started, 'reads_measurements': False}
    write_json(out / 'run.json', receipt)
    return receipt


def merge(cases_dir: Path, shards: list, out: Path) -> dict:
    """Join shard outputs into one results file, checking shard provenance and completeness."""
    build_receipt = json.loads((cases_dir / 'build.json').read_text())
    seen, receipts = {}, {}
    for directory in shards:
        receipt = json.loads((directory / 'run.json').read_text())
        if (receipt['cases_sha256'] != build_receipt['cases_sha256']
                or receipt['implementation'] != code_hashes()
                or receipt['results_sha256'] != sha(directory / 'results.jsonl')):
            raise ValueError(f'{directory}: shard provenance mismatch')
        if receipt['shard_count'] != len(shards) or receipt['shard_index'] in receipts:
            raise ValueError(f'{directory}: duplicate or incomplete shard set')
        receipts[receipt['shard_index']] = {'dir': str(directory), 'run_sha256': sha(directory / 'run.json')}
        for row in rows(directory / 'results.jsonl'):
            if row['case_id'] in seen:
                raise ValueError(f"duplicate result for {row['case_id']}")
            seen[row['case_id']] = row
    if set(receipts) != set(range(len(shards))):
        raise ValueError('incomplete shard set')
    order = [c['case_id'] for c in rows(cases_dir / 'cases.jsonl')]
    if set(order) != set(seen):
        raise ValueError('shards do not cover every case exactly once')
    out.mkdir(parents=True, exist_ok=False)
    counts = defaultdict(int)
    with (out / 'results.jsonl').open('w') as stream:
        for case_id in order:
            counts[seen[case_id]['status']] += 1
            stream.write(json.dumps(seen[case_id], sort_keys=True, allow_nan=False) + '\n')
    receipt = {'schema': SCHEMA + ':run', 'build_sha256': sha(cases_dir / 'build.json'),
               'cases_sha256': build_receipt['cases_sha256'], 'counts': dict(counts),
               'results_sha256': sha(out / 'results.jsonl'), 'implementation': code_hashes(),
               'python': sys.version, 'platform': platform.platform(),
               'merged_shards': receipts, 'reads_measurements': False}
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


def observed_class(observation: dict | None) -> str | None:
    if observation is None:
        return None
    value = observation.get('mean_od600')
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError('malformed observed OD')
    return 'growth' if value > GROWTH_OD else 'no_growth'


def comparable(cls: str | None) -> str | None:
    if cls in ('growth', 'no_growth'):
        return cls
    if cls == 'infeasible_under_constraints':
        return 'no_growth'
    return None


def score(cases_dir: Path, run_dir: Path, observations: Path, out: Path, n_boot: int = 2000) -> dict:
    cases = {c['case_id']: c for c in rows(cases_dir / 'cases.jsonl')}
    results = {r['case_id']: r for r in rows(run_dir / 'results.jsonl')}
    if set(cases) != set(results):
        raise ValueError('cases and results differ')
    obs = {o['case_id']: o for o in rows(observations)}

    scored, unavailable, noncomparable, failures = [], defaultdict(int), defaultdict(int), []
    for case_id, case in sorted(cases.items()):
        result = results[case_id]
        if result['status'] != 'completed':
            failures.append({'case_id': case_id, 'failure': result.get('failure')})
            continue
        for condition, assay_id in case['assays'].items():
            if assay_id is None:
                continue
            truth = None if case['kind'] == 'wild_type_control' else TRUTH[condition]
            measured = observed_class(obs.get(assay_id))
            if measured is None:
                noncomparable['no_measurement'] += 1
                continue
            predicted, feasible = {}, []
            for hypothesis in HYPOTHESES:
                cell = result['cells'].get(f'{condition}|{hypothesis}')
                if cell is None or cell['status'] == 'unavailable':
                    predicted[hypothesis] = 'unavailable'
                    unavailable[hypothesis] += 1
                    continue
                cls = comparable(cell.get('wild_type', {}).get('class'))
                predicted[hypothesis] = cls or 'not_comparable'
                if cls is None:
                    noncomparable[hypothesis] += 1
                elif cls == measured:
                    feasible.append(hypothesis)
            scored.append({'case_id': case_id, 'condition': condition, 'assay_case_id': assay_id,
                           'kind': case['kind'], 'scope': case['scope'], 'species': case['species'],
                           'strain_id': case['strain_id'], 'environment': case['environment'],
                           'arm': case['arm'], 'focal_reaction': case['focal_reaction'],
                           'truth': truth, 'measured_class': measured, 'predicted': predicted,
                           'feasible': feasible, 'feasible_size': len(feasible),
                           'truth_in_feasible': None if truth is None else truth in feasible,
                           'unique': len(feasible) == 1,
                           'uniquely_truth': None if truth is None else feasible == [truth]})

    def summarise(subset, label):
        out_rows = []
        for hypothesis in HYPOTHESES:
            with_truth = [s for s in subset if s['truth'] == hypothesis]
            in_feasible = [s for s in with_truth if s['truth_in_feasible']]
            selected = [s for s in subset if hypothesis in s['feasible']]
            unique_truth = [s for s in with_truth if s['uniquely_truth']]
            out_rows.append({
                'group': label, 'explanation': hypothesis,
                'cases_with_truth': len(with_truth), 'truth_in_feasible': len(in_feasible),
                'recall': (len(in_feasible) / len(with_truth)) if with_truth else None,
                'recall_ci': wilson(len(in_feasible), len(with_truth)),
                'cases_feasible_any_truth': len(selected),
                'precision': (len(in_feasible) / len(selected)) if selected else None,
                'precision_ci': wilson(len(in_feasible), len(selected)),
                'unique_identification': (len(unique_truth) / len(with_truth)) if with_truth else None,
                'unique_identification_ci': wilson(len(unique_truth), len(with_truth)),
                'available_cells': sum(1 for s in subset if s['predicted'].get(hypothesis) not in ('unavailable', None)),
                'unavailable_cells': sum(1 for s in subset if s['predicted'].get(hypothesis) == 'unavailable'),
            })
        return out_rows

    primary = [s for s in scored if s['scope'] == 'primary']
    summary = []
    summary += summarise(primary, 'primary:all_arms')
    for arm in sorted({s['arm'] for s in primary}):
        summary += summarise([s for s in primary if s['arm'] == arm], f'primary:{arm}')
    for species in sorted({s['species'] for s in primary}):
        summary += summarise([s for s in primary if s['species'] == species], f'primary:{species}')
    summary += summarise([s for s in scored if s['scope'] == 'sensitivity'], 'sensitivity:all_arms')

    controls = [s for s in scored if s['kind'] == 'wild_type_control']
    control_summary = {
        'cases': len(controls),
        'any_feasible_explanation': sum(1 for s in controls if s['feasible']),
        'uniquely_selected_explanation': sum(1 for s in controls if s['unique']),
        'unique_rate': (sum(1 for s in controls if s['unique']) / len(controls)) if controls else None,
        'unique_rate_ci': wilson(sum(1 for s in controls if s['unique']), len(controls)),
        'by_explanation': {h: sum(1 for s in controls if h in s['feasible']) for h in HYPOTHESES},
    }

    # Descriptive clustered bootstrap over (species, strain, environment) for the primary recalls.
    clusters = defaultdict(list)
    for s in primary:
        clusters[(s['species'], s['strain_id'], s['environment'])].append(s)
    keys = sorted(clusters)
    boot = {}
    if keys:
        import random
        rng = random.Random(42)
        for hypothesis in HYPOTHESES:
            draws = []
            for _ in range(n_boot):
                num = den = 0
                for _ in keys:
                    for s in clusters[keys[rng.randrange(len(keys))]]:
                        if s['truth'] == hypothesis:
                            den += 1
                            num += bool(s['truth_in_feasible'])
                if den:
                    draws.append(num / den)
            draws.sort()
            boot[hypothesis] = ({'low': draws[int(.025 * len(draws))], 'high': draws[int(.975 * len(draws)) - 1],
                                 'n_boot': len(draws), 'n_clusters': len(keys)} if draws else None)

    report = {
        'schema': SCHEMA + ':score', 'truth': TRUTH,
        'interpretation': 'Whether the known explanation is in the feasible set of the four-scenario comparison, '
                          'scored separately from growth agreement. Feasibility is conditional on the model, '
                          'medium and bounds; it is not biological proof of an explanation.',
        'inputs': {'cases': sha(cases_dir / 'cases.jsonl'), 'build': sha(cases_dir / 'build.json'),
                   'results': sha(run_dir / 'results.jsonl'), 'run': sha(run_dir / 'run.json'),
                   'observations': sha(observations)},
        'implementation': code_hashes(),
        'scored_assay_cells': len(scored), 'failures': failures,
        'unavailable_by_explanation': dict(unavailable), 'non_comparable': dict(noncomparable),
        'summary': summary, 'wild_type_controls': control_summary,
        'clustered_bootstrap_recall_primary': boot, 'cases': scored,
    }
    write_json(out, report)
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    builder = sub.add_parser('build')
    builder.add_argument('--plan-dir', type=Path, required=True)
    builder.add_argument('--ec-model', type=Path, required=True)
    builder.add_argument('--bs-model', type=Path, required=True)
    builder.add_argument('--out', type=Path, required=True)
    runner = sub.add_parser('run')
    runner.add_argument('--cases', type=Path, required=True)
    runner.add_argument('--out', type=Path, required=True)
    runner.add_argument('--ec-model', type=Path, required=True)
    runner.add_argument('--bs-model', type=Path, required=True)
    runner.add_argument('--workers', type=int, default=1)
    runner.add_argument('--keep-manifests', type=Path)
    runner.add_argument('--shard-index', type=int, default=0)
    runner.add_argument('--shard-count', type=int, default=1)
    merger = sub.add_parser('merge')
    merger.add_argument('--cases', type=Path, required=True)
    merger.add_argument('--shards', type=Path, nargs='+', required=True)
    merger.add_argument('--out', type=Path, required=True)
    scorer = sub.add_parser('score')
    scorer.add_argument('--cases', type=Path, required=True)
    scorer.add_argument('--run', type=Path, required=True)
    scorer.add_argument('--observations', type=Path, required=True)
    scorer.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.command == 'build':
        receipt = build(args.plan_dir, {'Ec': str(args.ec_model), 'Bs': str(args.bs_model)}, args.out)
        print(json.dumps({k: receipt[k] for k in ('cases', 'counts')}))
    elif args.command == 'run':
        receipt = run(args.cases, args.out, {'Ec': str(args.ec_model), 'Bs': str(args.bs_model)},
                      args.workers, args.keep_manifests, args.shard_index, args.shard_count)
        print(json.dumps({k: receipt[k] for k in ('counts', 'elapsed_s')}))
    elif args.command == 'merge':
        receipt = merge(args.cases, list(args.shards), args.out)
        print(json.dumps({k: receipt[k] for k in ('counts',)}))
    else:
        report = score(args.cases, args.run, args.observations, args.out)
        print(json.dumps({'scored_assay_cells': report['scored_assay_cells'],
                          'failures': len(report['failures'])}))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
