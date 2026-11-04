#!/usr/bin/env python3
"""One Price genome with the contract-bearing DIAMOND claims through the production engine in the arms of plan
`tasks/v3-price-contract-learned-2026-11-02.md` (step I).

The commands are those of the row-5 runner (`price_calibration_context/run_context_cell.py`, imported, not copied):

  stage0        --enable-priors none
  stage1        default priors
  learned       default priors + --learned-model <logistic>
  learned_gbdt  default priors + --learned-model <gbdt>
  gaps          row-5 `gaps` arm (pathway observations, reaction-local-context suggester, reaction graph/aliases)
  gaps_learned  gaps + --learned-model <logistic> (no cohort)
  universal / pathway / organism   the row-5 mode arms unchanged (gaps + cohort inputs + mode + logistic model)

Per arm: command, input hashes, exit code, time, output hash, rows, the EC scores of the reference proteins
(maximum over duplicate rows, duplicates counted), the stage ledger, and from the learned report its counts,
unsupported reasons by candidate type, and every EC candidate of a reference protein with status, probability and
reason. Nothing here reads a label beyond the reference protein identifiers.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
V3 = HERE.parent
sys.path.insert(0, str(V3))
from price_calibration_context import run_context_cell as R  # noqa: E402

ARMS = ('stage0', 'stage1', 'learned', 'learned_gbdt', 'gaps', 'gaps_learned', 'universal', 'pathway', 'organism',
        'obs_learned', 'rlc_learned', 'cohort_universal', 'cohort_pathway', 'cohort_organism')
LEARNED_ARMS = {'learned', 'learned_gbdt', 'gaps_learned', 'universal', 'pathway', 'organism',
                'obs_learned', 'rlc_learned', 'cohort_universal', 'cohort_pathway', 'cohort_organism'}
# Added after the pilot (plan deviation 2026-11-02, no label read): with pathway observations and the reaction-local-
# context suggester the model abstains on every candidate (unseen availability profile). Diagnostic arms separate the
# two inputs; the cohort_* arms are the three context modes without either, so that the learned stage can score.
DERIVED = {'obs_learned': ('gaps_learned', ('--rlc-suggester',), ('--suggestions-out',)),
           'rlc_learned': ('gaps_learned', (), ('--pathway-observations',)),
           'cohort_universal': ('universal', ('--rlc-suggester',), ('--suggestions-out', '--pathway-observations')),
           'cohort_pathway': ('pathway', ('--rlc-suggester',), ('--suggestions-out', '--pathway-observations')),
           'cohort_organism': ('organism', ('--rlc-suggester',), ('--suggestions-out', '--pathway-observations'))}


def strip(cmd, flags, valued):
    out, skip = [], False
    for x in cmd:
        if skip:
            skip = False
            continue
        if x in flags:
            continue
        if x in valued:
            skip = True
            continue
        out.append(x)
    return out
SCHEMA = 'gspa.price-contract-learned-cell/1'


def command(args, cell, arm, directory, scratch):
    if arm in DERIVED:
        parent, flags, valued = DERIVED[arm]
        cmd, output, bound = command(args, cell, parent, directory, scratch)
        return strip(cmd, flags, valued), output, bound
    base = {'stage1': 'stage0', 'learned': 'stage0', 'learned_gbdt': 'stage0', 'gaps_learned': 'gaps'}.get(arm, arm)
    margs = argparse.Namespace(**vars(args))
    cmd, output, bound = R.command(margs, cell, base, directory, args.inputs, scratch)
    if base == 'stage0' and arm != 'stage0':
        i = cmd.index('--enable-priors')
        del cmd[i:i + 2]
    if arm in ('learned', 'gaps_learned'):
        cmd += ['--learned-model', str(args.model)]
        bound[str(args.model)] = None
    if arm == 'learned_gbdt':
        cmd += ['--learned-model', str(args.model_gbdt)]
        bound[str(args.model_gbdt)] = None
    return cmd, output, bound


def retain(output, reference):
    h = hashlib.sha256()
    rows, dup = 0, 0
    types = Counter()
    scores = {pid: {} for pid in reference}
    with Path(output).open('rb') as f:
        header = None
        for raw in f:
            h.update(raw)
            line = raw.decode('utf-8').rstrip('\n')
            if header is None:
                header = line.split('\t')
                idx = {k: header.index(k) for k in ('protein_id', 'type', 'function_id', 'posterior_prob')}
                continue
            rows += 1
            x = line.split('\t')
            types[x[idx['type']]] += 1
            pid = x[idx['protein_id']]
            if pid in scores and x[idx['type']] == 'EC':
                ec = x[idx['function_id']].removeprefix('EC:')
                s = float(x[idx['posterior_prob']])
                if ec in scores[pid]:
                    dup += 1
                    s = max(s, scores[pid][ec])
                scores[pid][ec] = s
    return dict(rows=rows, rows_by_type=dict(types), duplicate_reference_ec_rows=dup, sha256=h.hexdigest()), scores


def learned_summary(path, reference):
    p = Path(path)
    if not p.is_file():
        return None
    d = json.loads(p.read_text())
    by_type = Counter()
    reasons = Counter()
    ref = {}
    for c in d.get('candidates', []):
        pid, kind, fid = c['function_key'].split('|', 2)
        by_type[(kind, c['score_status'])] += 1
        if c['score_status'] != 'scored':
            reasons[(kind, c.get('unsupported_reason'))] += 1
        if pid in reference and kind == 'EC':
            ref.setdefault(pid, []).append(dict(ec=fid.removeprefix('EC:'), role=c.get('role'), status=c['score_status'],
                                                probability=c.get('probability'), reason=c.get('unsupported_reason')))
    return dict(model=dict(id=d['model']['id'], sha256=d['model']['sha256'], kind=d['model']['kind']), counts=d.get('counts'),
                locus_layer_unavailable_reason=d.get('locus_layer_unavailable_reason'),
                status_by_type={f'{k}/{s}': n for (k, s), n in sorted(by_type.items())},
                unsupported_reasons={f'{k}/{r}': n for (k, r), n in sorted(reasons.items(), key=str)},
                contract_unavailable=sum(n for (k, r), n in reasons.items() if r == 'candidate_predictor_contract_unavailable'),
                reference_ec_candidates=ref, sha256=R.sha(p))


def run(args):
    plan = json.loads(Path(args.plan).read_text())
    cell = next(c for c in plan['cells'] if c['index'] == args.index)
    out = Path(args.out) / str(args.index)
    out.mkdir(parents=True, exist_ok=False)
    scratch = Path(args.scratch)
    scratch.mkdir(parents=True, exist_ok=True)
    reference = R.read_manifest_reference(args.manifest, cell['assembly'])
    receipt = dict(schema=SCHEMA, plan='tasks/v3-price-contract-learned-2026-11-02.md', cell=cell, arms={},
                   started=time.strftime('%Y-%m-%dT%H:%M:%S%z'), host=os.uname().nodename, heap=args.heap)
    retained = dict(cell=cell, reference=reference, scores={}, learned={})
    outputs = {}
    for arm in args.arms:
        directory = out / arm
        directory.mkdir()
        (scratch / arm).mkdir(exist_ok=True)
        cmd, output, bound = command(args, cell, arm, directory, scratch / arm)
        for p in bound:
            if bound[p] is None:
                bound[p] = R.sha(p)
        t0 = time.time()
        with (directory / 'engine.log').open('w') as log:
            rc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, timeout=args.timeout).returncode
        record = dict(command=cmd, inputs=bound, returncode=rc, seconds=round(time.time() - t0, 1))
        if rc == 0 and output.is_file():
            info, scores = retain(output, reference)
            retained['scores'][arm] = scores
            outputs[arm] = output
            record.update(output=str(output), output_sha256=info['sha256'], prediction_rows=info['rows'],
                          rows_by_type=info['rows_by_type'], duplicate_reference_ec_rows=info['duplicate_reference_ec_rows'],
                          published_claims=cell['claims'], rows_equal_claims=info['rows'] == cell['claims'],
                          reference_proteins_with_ec_candidates=sum(bool(v) for v in scores.values()))
            record['stages'] = R.ledger_scalars(str(output) + '.stages.json')
            gaps = Path(str(output) + '.gaps.jsonl')
            if gaps.is_file():
                record['gap_records'] = sum(1 for line in gaps.read_text().splitlines() if line.strip())
            ls = learned_summary(str(output) + '.learned.json', reference)
            if ls is not None:
                retained['learned'][arm] = ls.pop('reference_ec_candidates')
                record['learned'] = ls
            gap_context = R.small_json(str(output) + '.gap-context.json')
            if isinstance(gap_context, dict):
                record['gap_context'] = {k: v for k, v in gap_context.items() if not isinstance(v, (list, dict))}
            p = directory / 'resources.txt'
            if p.is_file():
                record['resources.txt'] = p.read_text()[-3000:]
        receipt['arms'][arm] = record
    receipt['differences'] = {}
    names = [a for a in args.arms if a in outputs]
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            try:
                receipt['differences'][f'{a}_vs_{b}'] = R.count_differences(outputs[a], outputs[b])
            except Exception as e:  # noqa: BLE001 - recorded
                receipt['differences'][f'{a}_vs_{b}'] = dict(error=repr(e))
    with (out / 'retained.json').open('x') as f:
        json.dump(retained, f, sort_keys=True)
        f.write('\n')
    receipt['retained_sha256'] = R.sha(out / 'retained.json')
    receipt['finished'] = time.strftime('%Y-%m-%dT%H:%M:%S%z')

    def learned_ok(arm, r):
        if arm not in LEARNED_ARMS:
            return True
        st = next((s for s in r.get('stages', []) if s.get('name') == 'learned_integration'), {})
        return st.get('status') == 'executed' and (r.get('learned') or {}).get('contract_unavailable') == 0
    receipt['learned_executed'] = {a: learned_ok(a, r) for a, r in receipt['arms'].items() if a in LEARNED_ARMS}
    receipt['status'] = 'complete' if all(r.get('returncode') == 0 and 'output' in r for r in receipt['arms'].values()) else 'failed'
    with (out / 'receipt.json').open('x') as f:
        json.dump(receipt, f, indent=1, sort_keys=True)
        f.write('\n')
    print(json.dumps(dict(index=args.index, status=receipt['status'], learned_executed=receipt['learned_executed'],
                          arms={a: (r.get('returncode'), r.get('seconds')) for a, r in receipt['arms'].items()})))
    return 0 if receipt['status'] == 'complete' else 1


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan', type=Path, required=True)
    p.add_argument('--index', type=int, required=True)
    p.add_argument('--arms', nargs='+', default=list(ARMS), choices=ARMS)
    p.add_argument('--layout', type=Path, required=True)
    p.add_argument('--reference', type=Path, required=True)
    p.add_argument('--jar', type=Path, required=True)
    p.add_argument('--inputs', type=Path, required=True)
    p.add_argument('--orthogroups', type=Path, required=True)
    p.add_argument('--pathway-groups', type=Path, required=True)
    p.add_argument('--reaction-graph', type=Path, required=True)
    p.add_argument('--reaction-aliases', type=Path, required=True)
    p.add_argument('--model', type=Path, required=True)
    p.add_argument('--model-gbdt', type=Path, required=True)
    p.add_argument('--layout-metadata', type=Path, default=None)
    p.add_argument('--manifest', type=Path, required=True)
    p.add_argument('--heap', default='16g')
    p.add_argument('--timeout', type=int, default=7200)
    p.add_argument('--scratch', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--feature-export', action='store_true')
    return run(p.parse_args(argv))


if __name__ == '__main__':
    sys.exit(main())
