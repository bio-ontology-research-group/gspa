#!/usr/bin/env python3
"""WP5 independent COBRApy/GLPK replay on the public, already-masked models.

Replays every available arm-visited cell through ``check_masked_assignment_cobra.py``. Each cell gets
its own model copy; wild-type and target-knockout endpoints cannot contaminate another cell. Optionally
adds all exhaustive cells or a seeded fraction of exhaustive cells not visited by any arm.

  python benchmark/v3/check_wp5_cobra.py --public BUNDLE/public --runs RUNS --out ORACLE \
      [--include-exhaustive | --rest-sample 0.1] [--shard i/n] [--arm-stages deterministic]

No original SBML, mask-gene identity, private truth or measured label is accepted. Agreement checks only
execution semantics, never biological correctness.
"""
from __future__ import annotations

import argparse
import json
import random
import sys
from pathlib import Path

import cobra

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import prepare_wp5_cases as P5  # noqa: E402
import run_wp5_recovery as R  # noqa: E402
import check_masked_assignment_cobra as O  # noqa: E402

SCHEMA = 'gspa.wp5.cobra_oracle/1'
REST_SEED = 20261029


def _load_model(path):
    model = cobra.io.read_sbml_model(str(path), f_replace={})
    model.solver = 'glpk'
    return model


def replay_cell(model, public, run, cell):
    """Replay one cell from one isolated copy of the masked COBRA model."""
    hypotheses = {h['id']: h for h in public['hypotheses']}
    conditions = {c['id']: c for c in public['conditions']}
    hypothesis = hypotheses[cell['hypothesis_id']]
    condition = conditions[cell['condition_id']]
    cell_model = model.copy()  # mandatory cell isolation (§8.6)
    endpoints, comparisons = {}, []
    for kind in ('wild_type', 'target_knockout'):
        oracle = O.simulate(cell_model, public, hypothesis, condition, run['config'], kind == 'target_knockout')
        endpoints[kind] = oracle
        observed = {k: v for k, v in cell['phenotypes'][kind].items() if k != 'knocked_out_reactions'}
        comparisons.append({'hypothesis_id': hypothesis['id'], 'condition_id': condition['id'],
                            'perturbation': kind,
                            'agrees': set(observed) == set(oracle) and O.agrees(observed, oracle),
                            'oracle': oracle})
    expected = O.demand_assessment(endpoints, run['config'])
    cell_agrees = (cell.get('contrast') == expected['contrast']
                   and cell.get('status') == expected['growth_status'])
    cell_comparison = {'hypothesis_id': hypothesis['id'], 'condition_id': condition['id'],
                       'agrees': cell_agrees,
                       'oracle': {'contrast': expected['contrast'], 'status': expected['growth_status']}}
    return comparisons, cell_comparison


def compare(public, runs, include_exhaustive=False, rest_sample=0.0, shard=None, loader=None,
            arm_stages=('deterministic', 'qwen')):
    """arm_stages names the arm stages whose runs must exist (Shaheen has no Qwen endpoint)."""
    if not set(arm_stages) <= {'deterministic', 'qwen'} or not arm_stages:
        raise P5.WP5Error('arm stages must be a non-empty subset of deterministic, qwen')
    if not 0.0 <= rest_sample <= 1.0:
        raise P5.WP5Error('rest sample must be in [0,1]')
    if include_exhaustive and rest_sample:
        raise P5.WP5Error('choose full exhaustive replay or a rest sample, not both')
    public, problems, protocol = R.load_protocol(public, runs)
    runs = Path(runs)
    protocol_sha = P5.sha_file(runs / 'protocol.json')
    loader = loader or _load_model
    selected = []
    for stage in arm_stages:
        selected.extend(R.plan_runs(protocol, stage))
    if include_exhaustive or rest_sample:
        selected.extend(R.plan_runs(protocol, 'exhaustive'))
    selected = R.shard_of(selected, shard)
    loaded = {}
    for name, arm, k, seed in selected:
        path = runs / 'runs' / R.run_name(name, arm, k, seed)
        if not path.is_file():
            raise P5.WP5Error('missing run for COBRA replay: ' + str(path))
        run = P5.read_json(path)
        entry = protocol['problems'][name]
        if (run.get('protocol_sha256') != protocol_sha
                or run.get('public_problem_sha256') != entry['public_problem_sha256']
                or run.get('status') == 'cell_mismatch'):
            raise P5.WP5Error('invalid run provenance for COBRA replay: ' + str(path))
        record = problems[name]
        model_path = public / record['masked_model']
        if P5.sha_file(model_path) != entry['masked_model_sha256']:
            raise P5.WP5Error('masked model changed before COBRA replay')
        if name not in loaded:
            loaded[name] = (loader(model_path), record['public_problem'], P5.sha_file(model_path))
        cells = run.get('cells')
        if not isinstance(cells, list):
            raise P5.WP5Error('run lacks a cell list: ' + str(path))
        loaded.setdefault('runs', [])
        loaded['runs'].append((name, arm, k, seed, run, cells))

    planned = loaded.pop('runs')
    # For a rest sample, remove exhaustive cells visited by any selected arm, then seed-shuffle the rest.
    if rest_sample:
        visited = {(cell['hypothesis_id'], cell['condition_id'])
                   for name, arm, k, seed, run, cells in planned if arm != 'exhaustive' for cell in cells}
        sampled = []
        for name, arm, k, seed, run, cells in planned:
            if arm != 'exhaustive':
                sampled.append((name, arm, k, seed, run, cells))
                continue
            rest = [c for c in cells if (c['hypothesis_id'], c['condition_id']) not in visited]
            random.Random(f'{REST_SEED}:{name}').shuffle(rest)
            take = int(len(rest) * rest_sample)
            sampled.append((name, arm, k, seed, run, rest[:take]))
        planned = sampled

    comparisons, cell_comparisons, run_rows = [], [], []
    model_hashes = {}
    for name, arm, k, seed, run, cells in planned:
        model, public_problem, model_sha = loaded[name]
        model_hashes[name] = model_sha
        before = len(comparisons)
        for cell in cells:
            rows, cell_row = replay_cell(model, public_problem, run, cell)
            comparisons.extend(rows)
            cell_comparisons.append(cell_row)
        run_rows.append({'run': R.run_name(name, arm, k, seed), 'cells': len(cells),
                         'endpoint_comparisons': len(comparisons) - before})
    report = {
        'schema': SCHEMA, 'cobra_version': cobra.__version__, 'solver': 'glpk',
        'protocol_sha256': protocol_sha, 'masked_model_sha256': model_hashes,
        'include_exhaustive': include_exhaustive, 'rest_sample': rest_sample, 'rest_seed': REST_SEED,
        'arm_stages': list(arm_stages),
        'runs': run_rows, 'comparisons': comparisons, 'n_comparisons': len(comparisons),
        'n_agree': sum(r['agrees'] for r in comparisons),
        'cell_comparisons': cell_comparisons, 'n_cell_comparisons': len(cell_comparisons),
        'n_cell_agree': sum(r['agrees'] for r in cell_comparisons),
        'interpretation': 'Independent masked-model execution agreement, not function or phenotype accuracy.'}
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--public', type=Path, required=True)
    parser.add_argument('--runs', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--include-exhaustive', action='store_true')
    parser.add_argument('--rest-sample', type=float, default=0.0)
    parser.add_argument('--shard')
    parser.add_argument('--arm-stages', nargs='+', default=['deterministic', 'qwen'],
                        help='arm stages whose runs are replayed and must exist (default: deterministic qwen)')
    args = parser.parse_args(argv)
    report = compare(args.public, args.runs, args.include_exhaustive, args.rest_sample, args.shard,
                     arm_stages=tuple(args.arm_stages))
    P5.write_json(args.out, report)
    print(json.dumps({'cells': report['n_cell_comparisons'], 'agree': report['n_cell_agree']}))
    return 0 if (report['n_cell_comparisons'] and report['n_cell_agree'] == report['n_cell_comparisons']
                 and report['n_agree'] == report['n_comparisons']) else 1


if __name__ == '__main__':
    sys.exit(main())
