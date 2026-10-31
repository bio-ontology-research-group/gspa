#!/usr/bin/env python3
"""Merge sharded ``check_wp5_cobra.py --shard i/n`` reports into the report an unsharded run writes.

  python benchmark/v3/merge_wp5_cobra.py --runs RUNS --shard-dir DIR --shards N --out ORACLE.json \
      [--arm-stages deterministic]

Shard i must be ``DIR/shard-<i>-of-<N>.json``. The merge recomputes the run plan from ``RUNS/protocol.json``,
checks that every shard holds exactly the runs ``shard_of`` assigns to it (same order), that all shards agree on
every run-independent field, and then interleaves runs, endpoint comparisons and cell comparisons back into plan
order. The merged report is therefore identical to an unsharded replay. A receipt with the shard hashes is
written next to it (``<out>.merge.json``). Only rest-sample 0 is supported (the sampled rest is defined on the
full plan, not per shard).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import prepare_wp5_cases as P5  # noqa: E402
import run_wp5_recovery as R  # noqa: E402

SCHEMA = 'gspa.wp5.cobra_oracle/1'
RECEIPT_SCHEMA = 'gspa.wp5.cobra_oracle_merge/1'
SHARED = ('schema', 'cobra_version', 'solver', 'protocol_sha256', 'include_exhaustive', 'rest_sample',
          'rest_seed', 'arm_stages', 'interpretation')


class MergeError(Exception):
    pass


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def merge(protocol, protocol_sha, shard_reports, arm_stages=('deterministic',)):
    n = len(shard_reports)
    if n < 1:
        raise MergeError('no shards')
    plan = []
    for stage in arm_stages:
        plan.extend(R.plan_runs(protocol, stage))
    first = shard_reports[0]
    if first.get('schema') != SCHEMA:
        raise MergeError('unexpected schema: %r' % first.get('schema'))
    if first.get('protocol_sha256') != protocol_sha:
        raise MergeError('shard protocol hash differs from runs/protocol.json')
    if list(first.get('arm_stages', [])) != list(arm_stages):
        raise MergeError('shard arm stages differ from the requested ones')
    if first.get('rest_sample') or first.get('include_exhaustive'):
        raise MergeError('only arm-visited replays (rest sample 0, no exhaustive) can be merged')
    models = {}
    chunks = {}
    for i, rep in enumerate(shard_reports):
        for key in SHARED:
            if rep.get(key) != first.get(key):
                raise MergeError('shard %d differs from shard 0 in %s' % (i, key))
        expected = [R.run_name(*r) for r in R.shard_of(plan, '%d/%d' % (i, n))]
        got = [row['run'] for row in rep['runs']]
        if got != expected:
            raise MergeError('shard %d holds runs %d, expected %d (or a different order)' % (i, len(got), len(expected)))
        for name, sha in rep['masked_model_sha256'].items():
            if models.setdefault(name, sha) != sha:
                raise MergeError('masked model hash differs between shards: ' + name)
        e = c = 0
        for row in rep['runs']:
            ne, nc = row['endpoint_comparisons'], row['cells']
            if ne != 2 * nc:
                raise MergeError('run %s: %d endpoint comparisons for %d cells' % (row['run'], ne, nc))
            chunks[row['run']] = (row, rep['comparisons'][e:e + ne], rep['cell_comparisons'][c:c + nc])
            e += ne
            c += nc
        if e != len(rep['comparisons']) or c != len(rep['cell_comparisons']):
            raise MergeError('shard %d comparison counts do not add up' % i)
        if rep['n_comparisons'] != e or rep['n_cell_comparisons'] != c:
            raise MergeError('shard %d recorded counts are inconsistent' % i)
    runs, comparisons, cell_comparisons = [], [], []
    for r in plan:
        row, comp, cells = chunks[R.run_name(*r)]
        runs.append(row)
        comparisons.extend(comp)
        cell_comparisons.extend(cells)
    report = {key: first[key] for key in ('schema', 'cobra_version', 'solver', 'protocol_sha256')}
    report['masked_model_sha256'] = _plan_order(models, plan)
    report.update({key: first[key] for key in ('include_exhaustive', 'rest_sample', 'rest_seed', 'arm_stages')})
    report.update({'runs': runs, 'comparisons': comparisons, 'n_comparisons': len(comparisons),
                   'n_agree': sum(r['agrees'] for r in comparisons),
                   'cell_comparisons': cell_comparisons, 'n_cell_comparisons': len(cell_comparisons),
                   'n_cell_agree': sum(r['agrees'] for r in cell_comparisons),
                   'interpretation': first['interpretation']})
    return report


def _plan_order(models, plan):
    """Key order of an unsharded run: first appearance in plan order."""
    out = {}
    for r in plan:
        if r[0] in models and r[0] not in out:
            out[r[0]] = models[r[0]]
    return out


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--runs', type=Path, required=True)
    p.add_argument('--shard-dir', type=Path, required=True)
    p.add_argument('--shards', type=int, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--arm-stages', nargs='+', default=['deterministic'])
    a = p.parse_args(argv)
    protocol_path = a.runs / 'protocol.json'
    protocol = json.loads(protocol_path.read_text())
    paths = [a.shard_dir / ('shard-%d-of-%d.json' % (i, a.shards)) for i in range(a.shards)]
    missing = [str(x) for x in paths if not x.is_file()]
    if missing:
        raise MergeError('missing shard reports: ' + ', '.join(missing[:5]))
    report = merge(protocol, _sha(protocol_path), [json.loads(x.read_text()) for x in paths], tuple(a.arm_stages))
    a.out.parent.mkdir(parents=True, exist_ok=True)
    tmp = a.out.with_suffix(a.out.suffix + '.partial')
    P5.write_json(tmp, report)  # same writer as check_wp5_cobra.py: byte-identical to an unsharded replay
    tmp.replace(a.out)
    receipt = {'schema': RECEIPT_SCHEMA, 'merged': str(a.out), 'merged_sha256': _sha(a.out),
               'protocol_sha256': report['protocol_sha256'], 'shards': a.shards,
               'shard_sha256': {x.name: _sha(x) for x in paths}, 'merger_sha256': _sha(__file__)}
    Path(str(a.out) + '.merge.json').write_text(json.dumps(receipt, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'cells': report['n_cell_comparisons'], 'agree': report['n_cell_agree'],
                      'endpoints': report['n_comparisons'], 'endpoint_agree': report['n_agree']}))
    ok = report['n_cell_comparisons'] and report['n_cell_agree'] == report['n_cell_comparisons'] \
        and report['n_agree'] == report['n_comparisons']
    return 0 if ok else 1


if __name__ == '__main__':
    try:
        sys.exit(main())
    except MergeError as exc:
        print('merge error: %s' % exc, file=sys.stderr)
        sys.exit(2)
