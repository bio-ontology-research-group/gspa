#!/usr/bin/env python3
"""Score the row-5 context cells (plan Part B5) with the unchanged Price metrics.

Reads `cells/<index>/{receipt,retained}.json` written by run_context_cell.py, checks every arm's
candidate-universe preservation and the stage0 replica against rc1's retained scores, and reports per
predictor and arm: top-1/3/5, tie-aware top-1 share with genome bootstrap, fixed-0.5 calls, the
stage0 -> gaps contrast (corrections/harms at 0.5 and at the row-1 calibrated test thresholds when a
policy is given), the mode contrasts (differing scores over the whole output and on the 148 loci) and
the learned-stage abstention counts. Unconfirmed calls are never negatives.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
V3 = HERE.parent
sys.path.insert(0, str(V3))
import ec_identity as E  # noqa: E402
import price_reference_metrics as M  # noqa: E402
from price_calibration_context.row1_calibration import contrast, evaluate_arm, sha, save, load  # noqa: E402

SCHEMA = 'gspa.price-context-evaluation/1'
ARMS = ('stage0', 'gaps', 'universal', 'pathway', 'organism')


def collect(cells_dir):
    data = defaultdict(lambda: defaultdict(dict))   # arm -> predictor -> genome -> proteins
    receipts = {}
    for d in sorted(Path(cells_dir).iterdir(), key=lambda p: int(p.name) if p.name.isdigit() else -1):
        rec = d / 'receipt.json'
        if not rec.is_file():
            continue
        receipt = load(rec)
        retained = load(d / 'retained.json')
        if sha(d / 'retained.json') != receipt['retained_sha256']:
            raise ValueError(f'{d}: retained file changed')
        cell = receipt['cell']
        receipts[cell['index']] = receipt
        for arm, scores in retained['scores'].items():
            proteins = {pid: ({M.ec(x) for x in refs}, {M.ec(k): M.score(v) for k, v in scores[pid].items()})
                        for pid, refs in retained['reference'].items()}
            data[arm][cell['predictor']][cell['assembly']] = proteins
    return data, receipts


def rc1_scores(rc1_retained):
    out = defaultdict(dict)
    for path in sorted(Path(rc1_retained).glob('*.json')):
        rec = load(path)
        cell = rec['cell']
        out[cell['predictor']][cell['assembly']] = {pid: {M.ec(k): M.score(v) for k, v in rec['scores']['stage0'][pid].items()}
                                                    for pid in rec['reference']}
    return out


def run(args):
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=False)
    audit = load(args.ec_audit)
    entries = audit['entries']
    E.validate_transfer_graph(entries)
    data, receipts = collect(args.cells)
    policy = load(args.row1_policy) if args.row1_policy else None
    split = load(args.row1_split) if args.row1_split else None
    test_genomes = {g for g, a in split['assignment'].items() if a == 'test'} if split else None
    rc1 = rc1_scores(args.rc1_retained) if args.rc1_retained else None

    report = dict(schema=SCHEMA, cells=len(receipts), arms_present=sorted(data), predictors={}, integrity={},
                  completeness={}, resources={}, learned={}, differences={}, gaps={})
    # Gap derivation depends on the pathway observations only: the gap records written by different predictors' cells
    # of the same genome must be identical files.
    by_genome = defaultdict(set)
    for idx, r in receipts.items():
        g = r['cell']['assembly']
        f = Path(args.cells) / str(idx) / f'{g}.gaps.jsonl'
        if f.is_file():
            by_genome[g].add(sha(f))
            report['gaps'].setdefault(g, {})['records'] = sum(1 for l in f.read_text().splitlines() if l.strip())
    for g, hashes in by_genome.items():
        report['gaps'][g]['cells'] = sum(1 for r in receipts.values() if r['cell']['assembly'] == g)
        report['gaps'][g]['identical_across_predictors'] = len(hashes) == 1
    report['gaps_identical_everywhere'] = all(v.get('identical_across_predictors', True) for v in report['gaps'].values())
    # completeness and resources
    by_pred_arm = Counter()
    for idx, r in receipts.items():
        for arm, a in r['arms'].items():
            ok = a.get('returncode') == 0 and a.get('candidate_universe_preserved') is True
            by_pred_arm[(r['cell']['predictor'], arm, 'complete' if ok else 'failed')] += 1
            report['resources'].setdefault(r['cell']['predictor'], {}).setdefault(arm, []).append(a.get('seconds'))
            # The learned stage's counts come from the stage ledger (the .learned.json of dense cells exceeds the
            # size the cell runner keeps inline).
            for st in a.get('stages', []):
                if st.get('name') == 'learned_integration' and st.get('requested'):
                    acc = report['learned'].setdefault(r['cell']['predictor'], {}).setdefault(arm, dict(cells=0, input_candidates=0, candidates_scored=0, candidates_unsupported=0, rescored=0, reasons=Counter()))
                    acc['cells'] += 1
                    for k in ('input_candidates', 'candidates_scored', 'candidates_unsupported', 'rescored'):
                        acc[k] += int(st.get(k) or 0)
                    acc['reasons'][st.get('reason') or st.get('status')] += 1
        for k, v in r.get('differences', {}).items():
            report['differences'].setdefault(r['cell']['predictor'], {}).setdefault(k, []).append(v)
    report['completeness'] = {f'{p}/{a}/{s}': n for (p, a, s), n in sorted(by_pred_arm.items())}
    for p in report['learned']:
        for a in report['learned'][p]:
            report['learned'][p][a]['reasons'] = dict(report['learned'][p][a]['reasons'])
    report['resources'] = {p: {a: dict(cells=len(v), max_seconds=max(x for x in v if x is not None), total_seconds=sum(x for x in v if x is not None))
                               for a, v in arms.items()} for p, arms in report['resources'].items()}
    report['differences'] = {p: {k: dict(cells=len(v), rows=sum(x.get('rows', 0) for x in v), differing=sum(x.get('differing', 0) for x in v),
                                         errors=[x['error'] for x in v if 'error' in x]) for k, v in d.items()}
                             for p, d in report['differences'].items()}
    # stage0 replica identity against rc1
    if rc1 is not None:
        for p, by_g in data.get('stage0', {}).items():
            equal = all(by_g[g][pid][1] == rc1[p][g][pid] for g in by_g for pid in by_g[g])
            report['integrity'][f'stage0/{p}'] = dict(genomes=len(by_g), equal_to_rc1=equal)
    # per predictor and arm metrics; locus-level differences between arms
    predictors = sorted({p for arm in data.values() for p in arm})
    for p in predictors:
        res = dict(arms={}, contrasts={}, locus_differences={})
        for arm in ARMS:
            if p not in data.get(arm, {}):
                continue
            by_g = data[arm][p]
            pol = dict(threshold=None, status='not_calibrated')
            metrics, _, _ = evaluate_arm(by_g, entries, pol)
            metrics['genomes'] = len(by_g)
            res['arms'][arm] = metrics
            if policy and test_genomes:
                key = f'stage0/{p}'   # row-1 thresholds were chosen for stage0; applied to every row-5 arm as pre-registered
                tp = policy['arms']['primary'][key]
                sub = {g: v for g, v in by_g.items() if g in test_genomes}
                if sub:
                    m_test, _, _ = evaluate_arm(sub, entries, tp)
                    res['arms'][arm]['test_genomes_at_row1_threshold'] = dict(threshold=tp['threshold'], status=tp['status'],
                        calibrated_calls=m_test.get('calibrated_calls'), top1_match=m_test['top1_match'], proteins=m_test['proteins'])
        for before, after in (('stage0', 'gaps'), ('gaps', 'universal'), ('gaps', 'pathway'), ('gaps', 'organism'),
                              ('universal', 'pathway'), ('universal', 'organism'), ('pathway', 'organism')):
            if p in data.get(before, {}) and p in data.get(after, {}):
                common = sorted(set(data[before][p]) & set(data[after][p]))
                b = {g: data[before][p][g] for g in common}
                a = {g: data[after][p][g] for g in common}
                t = policy['arms']['primary'][f'stage0/{p}']['threshold'] if policy else None
                if test_genomes:
                    bt = {g: v for g, v in b.items() if g in test_genomes}
                    at = {g: v for g, v in a.items() if g in test_genomes}
                else:
                    bt, at = b, a
                res['contrasts'][f'{before}_to_{after}'] = dict(all_genomes=contrast(b, a, entries, None, None),
                    test_genomes=contrast(bt, at, entries, t, t) if bt else None)
                diff = sum(b[g][pid][1] != a[g][pid][1] for g in common for pid in b[g])
                res['locus_differences'][f'{before}_vs_{after}'] = dict(reference_proteins=sum(len(b[g]) for g in common),
                                                                       proteins_with_any_score_change=diff)
        report['predictors'][p] = res
    save(out / 'report.json', report)
    with (out / 'summary.tsv').open('x') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['predictor', 'arm', 'genomes', 'proteins', 'top1_match', 'top1_share_ci', 'top3_match', 'fixed05_matched', 'fixed05_unconfirmed',
                    'test_threshold', 'test_calls', 'test_supported', 'test_unconfirmed'])
        for p, res in report['predictors'].items():
            for arm, m in res['arms'].items():
                t = m.get('test_genomes_at_row1_threshold') or {}
                cc = t.get('calibrated_calls') or {}
                b = m['top1']['share']
                w.writerow([p, arm, m['genomes'], m['proteins'], m['top1_match'],
                            f"{b['estimate']:.3f} [{b['ci_low']:.3f}, {b['ci_high']:.3f}]", m['top3_match'],
                            m['fixed_0_5']['matched'], m['fixed_0_5']['unconfirmed'], t.get('threshold', ''),
                            cc.get('calls', ''), cc.get('supported', ''), cc.get('unconfirmed', '')])
    with (out / 'contrasts.tsv').open('x') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['predictor', 'contrast', 'genomes_scope', 'fixed05_corrected', 'fixed05_harmed', 'fixed05_net_ci', 'top1_gained', 'top1_lost',
                    'changed_top1', 'cal_threshold', 'cal_corrected', 'cal_harmed', 'loci_with_score_change', 'whole_output_rows', 'whole_output_differing'])
        for p, res in report['predictors'].items():
            for name, c in res['contrasts'].items():
                ld = res['locus_differences'][name.replace('_to_', '_vs_')]
                wd = report['differences'].get(p, {}).get(name.replace('_to_', '_vs_'), {})
                for scope in ('all_genomes', 'test_genomes'):
                    cc = c[scope]
                    if cc is None:
                        continue
                    n = cc['fixed_0_5']['net']
                    w.writerow([p, name, scope, cc['fixed_0_5']['corrected'], cc['fixed_0_5']['harmed'],
                                f"{n['estimate']:.2f} [{n['ci_low']:.2f}, {n['ci_high']:.2f}]", cc['top1']['gained'], cc['top1']['lost'],
                                cc['paired']['changed_top1'], cc['thresholds']['after'] if scope == 'test_genomes' else '',
                                cc['calibrated']['corrected'] if scope == 'test_genomes' else '', cc['calibrated']['harmed'] if scope == 'test_genomes' else '',
                                ld['proteins_with_any_score_change'] if scope == 'all_genomes' else '',
                                wd.get('rows', '') if scope == 'all_genomes' else '', wd.get('differing', '') if scope == 'all_genomes' else ''])
    print(json.dumps(dict(cells=report['cells'], completeness=report['completeness'], integrity=report['integrity'],
                          differences=report['differences'], gaps_identical_everywhere=report['gaps_identical_everywhere']), indent=1))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cells', type=Path, required=True)
    p.add_argument('--manifest', type=Path, required=True)
    p.add_argument('--ec-audit', type=Path, required=True)
    p.add_argument('--rc1-retained', type=Path, default=None)
    p.add_argument('--row1-policy', type=Path, default=None)
    p.add_argument('--row1-split', type=Path, default=None)
    p.add_argument('--out', type=Path, required=True)
    run(p.parse_args(argv))


if __name__ == '__main__':
    main()
