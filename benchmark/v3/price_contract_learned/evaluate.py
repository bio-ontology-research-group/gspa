#!/usr/bin/env python3
"""Plan steps E (rows 1 and 4) and F (row 5) for the contract-bearing DIAMOND claims on the 22 Price genomes.

Reads `cells/<index>/{receipt,retained}.json` (run_cell.py) and the claims files (arm `raw`), and scores with the
unchanged functions of the Price evaluator and the row-1 analysis (`price_reference_metrics`,
`row1_calibration.pairs/choose_threshold/evaluate_arm/calibration_metrics/contrast`). Complete four-level EC rows
are the scored population (as in every Price analysis); GO rows and partial ECs are counted, not scored.

Order (plan E1): thresholds are chosen on the calibration genomes of the frozen row-1 split and `policy.json` is
written and hashed before any test-genome metric is computed. Unconfirmed calls are never negatives; statistics
that need a binary outcome treat them as unsupported and are labelled reference-support statistics.
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
from price_calibration_context import row1_calibration as R1  # noqa: E402

SCHEMA = 'gspa.price-contract-learned-evaluation/1'
ROW1_ARMS = ('raw', 'stage0', 'stage1', 'learned', 'learned_gbdt')
ROW5_ARMS = ('gaps', 'gaps_learned', 'universal', 'pathway', 'organism', 'obs_learned', 'rlc_learned', 'cohort_universal', 'cohort_pathway', 'cohort_organism')
LEARNED = ('learned', 'learned_gbdt', 'gaps_learned', 'universal', 'pathway', 'organism', 'obs_learned', 'rlc_learned', 'cohort_universal', 'cohort_pathway', 'cohort_organism')
CONTRASTS = (('raw', 'stage1'), ('stage0', 'stage1'), ('stage1', 'learned'), ('stage1', 'learned_gbdt'), ('learned', 'learned_gbdt'))
MODE_PAIRS = (('stage1', 'gaps'), ('learned', 'gaps_learned'), ('gaps_learned', 'universal'), ('gaps_learned', 'pathway'),
              ('gaps_learned', 'organism'), ('universal', 'pathway'), ('universal', 'organism'), ('pathway', 'organism'),
              ('learned', 'obs_learned'), ('learned', 'rlc_learned'), ('learned', 'cohort_universal'), ('learned', 'cohort_pathway'),
              ('learned', 'cohort_organism'), ('cohort_universal', 'cohort_pathway'), ('cohort_universal', 'cohort_organism'),
              ('cohort_pathway', 'cohort_organism'))


def complete(ec):
    return E.ID.fullmatch(ec.removeprefix('EC:')) is not None


def read_raw(claims, reference):
    out = {pid: {} for pid in reference}
    partial = 0
    with Path(claims).open() as f:
        for line in f:
            r = json.loads(line)
            if r['protein_id'] in out and not r['function_id'].startswith('GO:'):
                if not complete(r['function_id']):
                    partial += 1
                    continue
                k = M.ec(r['function_id'])
                out[r['protein_id']][k] = max(out[r['protein_id']].get(k, 0.0), M.score(float(r['raw_score'])))
    return out, partial


def collect(cells, claims_dir):
    data = defaultdict(dict)           # arm -> genome -> pid -> (refs, scores)
    learned = defaultdict(dict)        # arm -> genome -> pid -> list of candidate records
    receipts, dropped = {}, Counter()
    for d in sorted(Path(cells).iterdir(), key=lambda p: int(p.name) if p.name.isdigit() else -1):
        if not (d / 'receipt.json').is_file():
            continue
        rec = R1.load(d / 'receipt.json')
        if R1.sha(d / 'retained.json') != rec['retained_sha256']:
            raise ValueError(f'{d}: retained file changed')
        ret = R1.load(d / 'retained.json')
        g = rec['cell']['assembly']
        receipts[g] = rec
        refs = {pid: {M.ec(x) for x in r} for pid, r in ret['reference'].items()}
        for arm, scores in ret['scores'].items():
            prot = {}
            for pid in refs:
                s = {}
                for k, v in scores[pid].items():
                    if not complete(k):
                        dropped[arm] += 1
                        continue
                    s[M.ec(k)] = M.score(v)
                prot[pid] = (refs[pid], s)
            data[arm][g] = prot
        for arm, by_pid in ret.get('learned', {}).items():
            learned[arm][g] = by_pid
        raw, partial = read_raw(Path(claims_dir) / f'{g}.claims.jsonl', refs)
        data['raw'][g] = {pid: (refs[pid], raw[pid]) for pid in refs}
        dropped['raw'] += partial
    return data, learned, receipts, dropped


def subset(by_genome, genomes):
    return {g: p for g, p in by_genome.items() if g in genomes}


def learned_scored_only(by_genome, learned_by_genome):
    """Restrict each reference protein's candidates to those the learned model scored (probability not null)."""
    out = {}
    for g, prot in by_genome.items():
        out[g] = {}
        for pid, (refs, scores) in prot.items():
            ok = {M.ec(c['ec']) for c in learned_by_genome.get(g, {}).get(pid, []) if c['status'] == 'scored' and complete(c['ec'])}
            out[g][pid] = (refs, {k: v for k, v in scores.items() if k in ok})
    return out


def at_threshold(by_genome, entries, t):
    pol = dict(threshold=t, status='fixed')
    m, _, _ = R1.evaluate_arm(by_genome, entries, pol)
    return m.get('calibrated_calls')


def learned_counts(receipts, learned, arm, reference_population):
    acc = dict(cells=0, executed=0, candidates=0, scored=0, unsupported=0, rescored=0, contract_unavailable=0,
               status_by_type=Counter(), unsupported_reasons=Counter(), reference_ec=Counter(), reference_ec_reasons=Counter())
    for g, rec in receipts.items():
        a = rec['arms'].get(arm)
        if not a or 'learned' not in a:
            continue
        acc['cells'] += 1
        st = next((s for s in a.get('stages', []) if s.get('name') == 'learned_integration'), {})
        acc['executed'] += st.get('status') == 'executed'
        c = a['learned']['counts'] or {}
        acc['candidates'] += c.get('candidates', 0)
        acc['rescored'] += c.get('rescored', 0)
        acc['unsupported'] += c.get('unsupported', 0)
        acc['scored'] += c.get('candidates', 0) - c.get('unsupported', 0)
        acc['contract_unavailable'] += a['learned']['contract_unavailable']
        acc['status_by_type'].update(a['learned']['status_by_type'])
        acc['unsupported_reasons'].update(a['learned']['unsupported_reasons'])
        for pid, cands in learned.get(arm, {}).get(g, {}).items():
            if (g, pid) not in reference_population:
                continue
            for x in cands:
                acc['reference_ec'][x['status']] += 1
                if x['status'] != 'scored':
                    acc['reference_ec_reasons'][x['reason']] += 1
    for k in ('status_by_type', 'unsupported_reasons', 'reference_ec', 'reference_ec_reasons'):
        acc[k] = dict(acc[k])
    return acc


def run(args):
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=False)
    entries = R1.load(args.ec_audit)['entries']
    E.validate_transfer_graph(entries)
    data, learned, receipts, dropped = collect(args.cells, args.claims)
    split = R1.load(args.split)
    display = {arm: R1.load(args.model_gbdt if arm == 'learned_gbdt' else args.model)['display_threshold']['value'] for arm in LEARNED}
    genomes = sorted(receipts)
    if len(genomes) != 22:
        raise SystemExit(f'expected 22 genomes, found {len(genomes)}')
    pop = {(g, pid) for g in genomes for pid in data['stage0'][g]}
    if len(pop) != 148:
        raise SystemExit(f'reference population {len(pop)} != 148')
    leakage = R1.load(args.leakage)
    excluded = set(leakage['sensitivity_exclusion'])
    cal = {g for g, a in split['assignment'].items() if a == 'calibration'}
    test = {g for g, a in split['assignment'].items() if a == 'test'}
    arms = [a for a in ROW1_ARMS + ROW5_ARMS if a in data and len(data[a]) == 22]

    # ---- E1: thresholds on calibration genomes only (both splits), written and hashed first
    policies = {}
    for name, cal_set in (('primary', cal), ('swapped', test)):
        pol = {}
        for arm in arms:
            rows = []
            for g, prot in subset(data[arm], cal_set).items():
                rows += R1.pairs(prot, entries)[0]
            pol[arm] = R1.choose_threshold(rows)
        policies[name] = pol
    policy = dict(schema=SCHEMA + ':policy', target_support=R1.PRECISION, min_calls=R1.MIN_DECISIONS,
                  selection_rule=R1.SELECTION_RULE, split_sha256=R1.sha(args.split),
                  calibration_genomes=dict(primary=sorted(cal), swapped=sorted(test)), arms=policies,
                  model_display_threshold=dict(values=display, source='model artifact, chosen on its training calibration split (advisory)'))
    policy_sha = R1.save(out / 'policy.json', policy)

    # ---- E2-E6 and F
    populations = {'all': None, 'without_GCA_000006965.1': {'GCA_000006965.1'}}
    if excluded:
        populations['without_leaked_proteins'] = excluded
    report = dict(schema=SCHEMA, plan='tasks/v3-price-contract-learned-2026-11-02.md', policy_sha256=policy_sha,
                  split_sha256=policy['split_sha256'], arms=arms, dropped_non_complete_ec_rows=dict(dropped),
                  convention='unconfirmed calls are not negatives; support share counts them as unsupported (reference support)',
                  analyses={}, all_genomes={}, learned={}, differences={}, rows={})

    def drop(by_genome, popname):
        rule = populations[popname]
        if rule is None:
            return by_genome
        if popname.startswith('without_GCA'):
            return {g: p for g, p in by_genome.items() if g not in rule}
        return {g: {pid: v for pid, v in p.items() if f'{g}__{pid}' not in rule} for g, p in by_genome.items()}

    for popname in populations:
        for name, test_set, cal_set in (('primary', test, cal), ('swapped', cal, test)):
            key = f'{name}/{popname}'
            an = dict(test_genomes=sorted(test_set), arms={}, contrasts={})
            for arm in arms:
                by = drop(subset(data[arm], test_set), popname)
                pol = policies[name][arm]
                res, _, _ = R1.evaluate_arm(by, entries, pol)
                res['calibration_metrics'] = R1.calibration_metrics(by, entries)
                res['at_0_5'] = at_threshold(by, entries, 0.5)
                if arm in LEARNED:
                    display_t = display[arm]
                    res['at_model_display_threshold'] = dict(threshold=display_t, **at_threshold(by, entries, display_t))
                    ls = learned_scored_only(by, learned.get(arm, {}))
                    res['learned_scored_only'] = dict(calibration_metrics=R1.calibration_metrics(ls, entries),
                                                      at_0_5=at_threshold(ls, entries, 0.5),
                                                      at_model_display_threshold=at_threshold(ls, entries, display_t))
                an['arms'][arm] = res
            for before, after in CONTRASTS:
                if before in arms and after in arms:
                    an['contrasts'][f'{before}_to_{after}'] = R1.contrast(
                        drop(subset(data[before], test_set), popname), drop(subset(data[after], test_set), popname), entries,
                        policies[name][before]['threshold'], policies[name][after]['threshold'])
            report['analyses'][key] = an
    # all 22 genomes (row 5 F1 and descriptive row 1)
    for arm in arms:
        pol = dict(threshold=None, status='not_calibrated')
        m, _, _ = R1.evaluate_arm(data[arm], entries, pol)
        m['at_0_5'] = at_threshold(data[arm], entries, 0.5)
        report['all_genomes'][arm] = m
    for before, after in MODE_PAIRS:
        if before in arms and after in arms:
            b, a = data[before], data[after]
            report['all_genomes'][f'contrast/{before}_to_{after}'] = R1.contrast(b, a, entries, None, None)
            report['all_genomes'][f'loci_changed/{before}_vs_{after}'] = sum(b[g][pid][1] != a[g][pid][1] for g in b for pid in b[g])
    # learned counts (E5) and whole-output differences (F2)
    for arm in LEARNED:
        if arm in arms:
            report['learned'][arm] = learned_counts(receipts, learned, arm, {(g, pid) for g, pid in pop})
    diffs = defaultdict(lambda: dict(cells=0, rows=0, differing=0, errors=[]))
    for g, rec in receipts.items():
        for k, v in rec.get('differences', {}).items():
            d = diffs[k]
            d['cells'] += 1
            d['rows'] += v.get('rows', 0)
            d['differing'] += v.get('differing', 0)
            if 'error' in v:
                d['errors'].append(v['error'])
        for arm, a in rec['arms'].items():
            r = report['rows'].setdefault(arm, dict(rows=0, claims=0, rows_equal_claims=0, seconds=0.0, by_type=Counter(), taxon=Counter()))
            r['rows'] += a.get('prediction_rows', 0)
            r['claims'] += a.get('published_claims', 0)
            r['rows_equal_claims'] += bool(a.get('rows_equal_claims'))
            r['seconds'] += a.get('seconds', 0)
            r['by_type'].update(a.get('rows_by_type', {}))
            st = next((s for s in a.get('stages', []) if s.get('name') == 'symbolic_enforcement'), {})
            r['taxon'][str(st.get('taxon_consistency_state'))] += 1
    report['rows'] = {k: {kk: (dict(vv) if isinstance(vv, Counter) else vv) for kk, vv in v.items()} for k, v in report['rows'].items()}
    report['differences'] = dict(diffs)
    # F3: power gate restated with the unchanged gap-linked loci; gap records compared with the row-5 run
    gap = R1.load(args.gap_report)
    row5 = R1.load(args.row5_report)
    rec_counts = {g: rec['arms'].get('gaps', {}).get('gap_records') for g, rec in receipts.items()}
    loci = {}
    for key, v in gap['loci'].items():
        g, pid = key.split('/', 1)
        entry = dict(reference_ecs=v['reference_ecs'], context_available=v['context_available'])
        for arm in ('raw', 'stage1', 'learned', 'gaps_learned', 'universal', 'cohort_universal'):
            if arm in data:
                refs, scores = data[arm][g][pid]
                r = M.protein(refs, scores, entries, 'exact', integrated=True)
                entry[arm] = dict(candidates=r['candidates'], reference_proposed=bool(r['candidate_reference_matches']),
                                  top1=r['top1_match'], reference_scores={e: scores.get(e) for e in refs})
        entry['learned_status'] = [c for c in learned.get('learned', {}).get(g, {}).get(pid, []) if c['ec'] in v['reference_ecs']]
        loci[key] = entry
    report['power'] = dict(gap_linked_loci=gap['gap_linked_loci'], gate=gap['power_gate'], powered=gap['gap_linked_loci'] >= gap['power_gate'],
                           loci_with_reference_ec_proposed=sum(1 for e in loci.values() if e.get('raw', {}).get('reference_proposed')),
                           loci=loci, gap_records_per_genome=rec_counts,
                           gap_records_equal_row5=all(rec_counts[g] == row5['gaps'][g]['records'] for g in rec_counts))
    report['exit_codes'] = {g: {a: r.get('returncode') for a, r in rec['arms'].items()} for g, rec in receipts.items()}
    report['learned_executed'] = {g: rec.get('learned_executed') for g, rec in receipts.items()}
    R1.save(out / 'report.json', report)
    write_tables(report, out)
    print(json.dumps(dict(policy=policy_sha, arms=arms, dropped=dict(dropped)), indent=1))


def f(x, nd=3):
    return '' if x is None else (f'{x:.{nd}f}' if isinstance(x, float) else str(x))


def ci(b, nd=3):
    return f"{f(b['estimate'], nd)} [{f(b['ci_low'], nd)}, {f(b['ci_high'], nd)}]" if b else ''


def write_tables(report, out):
    with (out / 'test_metrics.tsv').open('x') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['analysis', 'arm', 'proteins', 'threshold_status', 'threshold', 'cal_calls', 'cal_supported', 'cal_unconfirmed',
                    'cal_support_share_ci', 'cal_reference_recall_ci', 'top1_match', 'top1_share_ci', 'top3_match',
                    'calls_0_5', 'supported_0_5', 'unconfirmed_0_5', 'calls_display', 'supported_display', 'unconfirmed_display',
                    'positive_pair_brier', 'top_call_brier', 'top_call_ece',
                    'ls_calls_0_5', 'ls_supported_0_5', 'ls_calls_display', 'ls_supported_display', 'ls_positive_pair_brier', 'ls_top_call_ece'])
        for key, an in report['analyses'].items():
            for arm, r in an['arms'].items():
                cc = r.get('calibrated_calls') or {}
                c5 = r['at_0_5'] or {}
                cd = r.get('at_model_display_threshold') or {}
                cm = r['calibration_metrics']
                ls = r.get('learned_scored_only') or {}
                lc = ls.get('calibration_metrics') or {}
                w.writerow([key, arm, r['proteins'], r['status'], f(r['threshold'], 6), cc.get('calls', ''), cc.get('supported', ''),
                            cc.get('unconfirmed', ''), ci(cc.get('support_share')), ci(cc.get('reference_recall')),
                            r['top1_match'], ci(r['top1']['share']), r['top3_match'],
                            c5.get('calls'), c5.get('supported'), c5.get('unconfirmed'),
                            cd.get('calls', ''), cd.get('supported', ''), cd.get('unconfirmed', ''),
                            f(cm['positive_pair_brier'], 4), f(cm['top_call_brier'], 4), f(cm['top_call_ece'], 4),
                            (ls.get('at_0_5') or {}).get('calls', ''), (ls.get('at_0_5') or {}).get('supported', ''),
                            (ls.get('at_model_display_threshold') or {}).get('calls', ''),
                            (ls.get('at_model_display_threshold') or {}).get('supported', ''),
                            f(lc.get('positive_pair_brier'), 4), f(lc.get('top_call_ece'), 4)])
    with (out / 'contrasts.tsv').open('x') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['analysis', 'contrast', 't_before', 't_after', 'cal_corrected', 'cal_harmed', 'fixed05_corrected', 'fixed05_harmed',
                    'fixed05_net_ci', 'top1_gained', 'top1_lost', 'top1_net_ci', 'changed_top1'])
        for key, an in report['analyses'].items():
            for name, c in an['contrasts'].items():
                w.writerow([key, name, f(c['thresholds']['before'], 6), f(c['thresholds']['after'], 6), c['calibrated']['corrected'],
                            c['calibrated']['harmed'], c['fixed_0_5']['corrected'], c['fixed_0_5']['harmed'], ci(c['fixed_0_5']['net'], 2),
                            c['top1']['gained'], c['top1']['lost'], ci(c['top1']['net'], 2), c['paired']['changed_top1']])
    with (out / 'reliability.tsv').open('x') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['analysis', 'arm', 'population', 'table', 'bin', 'n', 'mean_score', 'supported_or_reference', 'unconfirmed_or_other'])
        for key, an in report['analyses'].items():
            for arm, r in an['arms'].items():
                for popn, cm in (('all_candidates', r['calibration_metrics']), ('learned_scored_only', (r.get('learned_scored_only') or {}).get('calibration_metrics'))):
                    if not cm:
                        continue
                    for b in cm['top_call_reliability']:
                        w.writerow([key, arm, popn, 'top_call', b['bin'], b['n'], f(b['mean_score']), f(b['supported_credit']), f(b['unconfirmed_credit'])])
                    for b in cm['pair_reliability']:
                        w.writerow([key, arm, popn, 'all_pairs', b['bin'], b['pairs'], f(b['mean_score']), b['reference_pairs'], b['pairs'] - b['reference_pairs']])
    with (out / 'all_genomes.tsv').open('x') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['arm', 'proteins', 'top1_match', 'top1_share_ci', 'top3_match', 'matched_0_5', 'unconfirmed_0_5'])
        for arm, m in report['all_genomes'].items():
            if '/' in arm:
                continue
            w.writerow([arm, m['proteins'], m['top1_match'], ci(m['top1']['share']), m['top3_match'],
                        m['fixed_0_5']['matched'], m['fixed_0_5']['unconfirmed']])


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cells', type=Path, required=True)
    p.add_argument('--claims', type=Path, required=True)
    p.add_argument('--split', type=Path, required=True)
    p.add_argument('--ec-audit', type=Path, required=True)
    p.add_argument('--model', type=Path, required=True)
    p.add_argument('--model-gbdt', type=Path, required=True)
    p.add_argument('--leakage', type=Path, required=True)
    p.add_argument('--gap-report', type=Path, required=True, help='price-calibration-context-v1/row5/gap_report.json')
    p.add_argument('--row5-report', type=Path, required=True, help='price-calibration-context-v1/row5/evaluation/report.json')
    p.add_argument('--out', type=Path, required=True)
    run(p.parse_args(argv))


if __name__ == '__main__':
    main()
