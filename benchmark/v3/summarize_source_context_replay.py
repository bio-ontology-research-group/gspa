"""Tabulate the verified source x context production replay and its EC analysis.

Reads analysis/results.json (gspa.source-context-analysis/1) and
replay/verification.json (gspa.source-context-production/1) from a results
directory and writes

* source_context_replay_results.tsv - one row per (operating point, metric
  family, arm or comparison) with every numeric field present, nested
  dictionaries flattened with dots;
* source_context_replay_integrity.tsv - one row per (source, context, kind)
  with replay actions, unsupported reasons and the largest Python/JVM
  probability error;
* source_context_replay_summary.md - compact tables for both.

Stdlib only, deterministic (no timestamps, fixed ordering), and strict: any
missing, extra or inconsistent structure raises SummaryError.
"""
import argparse
import csv
import json
from pathlib import Path

SOURCES = ('diamond_homology', 'diamond_family_excluded', 'clean', 'clean_diamond_homology',
           'clean_diamond_family_excluded')
CONTEXTS = ('SEQONLY', 'UNIVERSAL', 'PATHWAY', 'ORGANISM')
KINDS = ('rule', 'logistic', 'gbdt')
REPLAY_KINDS = ('logistic', 'gbdt')
NATIVE = 'clean_native'
ARM_OPERATING_POINTS = ('calibration_selected', 'ec_calibration_selected', 'fixed_0.5')
NATIVE_OPERATING_POINT = 'native_emitted_calls'
OP_ORDER = ARM_OPERATING_POINTS + (NATIVE_OPERATING_POINT,)
ARM_SECTIONS = {'comparisons', 'comparisons_to_native_emitted_calls',
                'conditional_candidate_metrics', 'full_ec_reference'}
EC_KEYS = {'emitted_reference_pairs', 'missed_reference', 'policy', 'proteins', 'proteins_with_candidates',
           'records', 'reference_matches', 'reference_micro_f1', 'reference_pairs',
           'reference_pairs_by_genome', 'reference_precision', 'reference_recall', 'test_genomes',
           'threshold', 'top1_matches', 'top1_tied', 'top3_matches', 'top5_matches', 'unconfirmed_calls'}
COND_KEYS = {'brier_denominator', 'brier_on_available_scores', 'label_basis', 'missed_reference',
             'missing_predictions', 'positive', 'reference_f1', 'reference_matched_calls', 'rows',
             'threshold', 'uncalled_unconfirmed_or_negative', 'unconfirmed_or_negative_calls'}
CMP_KEYS = {'per_genome', 'proteins', 'reference_calls_gained', 'reference_calls_lost', 'top1_gained',
            'top1_lost', 'unconfirmed_call_change'}
CMP_GENOME_KEYS = CMP_KEYS - {'per_genome'}
REPLAY_CELL_KEYS = {'actions', 'added', 'candidates', 'functions', 'max_probability_error', 'reasons', 'removed'}
REPLAY_TIE_KEYS = {'display_tie_candidates', 'display_ties'}  # present together on cells replayed after the tie fix
FAMILIES = ('full_ec_reference', 'conditional_candidate_metrics', 'comparison', 'comparison_to_native_emitted')
RESULTS_TSV = 'source_context_replay_results.tsv'
INTEGRITY_TSV = 'source_context_replay_integrity.tsv'
MD_NAME = 'source_context_replay_summary.md'
LEAD = ['operating_point', 'family', 'name', 'source', 'context', 'kind', 'baseline', 'comparison_type']


class SummaryError(ValueError):
    pass


def require(cond, msg):
    if not cond:
        raise SummaryError(msg)


def is_num(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool)


def keys_equal(obj, expected, where):
    require(isinstance(obj, dict), f'{where}: expected an object')
    require(set(obj) == set(expected), f'{where}: keys differ: missing {sorted(set(expected) - set(obj))}, '
            f'unexpected {sorted(set(obj) - set(expected))}')


def arm_names():
    return [f'{s}/{c}/{k}' for s in SOURCES for c in CONTEXTS for k in KINDS]


def split_arm(name):
    parts = name.split('/')
    require(len(parts) == 3 and parts[0] in SOURCES and parts[1] in CONTEXTS and parts[2] in KINDS,
            f'unrecognised arm {name!r}')
    return parts


def classify(key):
    """Return (subject, baseline, comparison_type) for a 'A_vs_B' comparison key."""
    require(key.count('_vs_') == 1, f'unrecognised comparison key {key!r}')
    a, b = key.split('_vs_')
    sa = split_arm(a)
    if b == NATIVE:
        return a, b, 'vs_clean_native'
    sb = split_arm(b)
    if sa[0] == sb[0] and sa[1] == sb[1] and sa[2] != 'rule' and sb[2] == 'rule':
        return a, b, 'learned_vs_rule'
    if sa[0] == sb[0] and sa[2] == sb[2] and sa[1] != 'SEQONLY' and sb[1] == 'SEQONLY':
        return a, b, 'context_vs_seqonly'
    if sa[0] != sb[0] and sa[1] == sb[1] and sa[2] == sb[2]:
        return a, b, 'source_vs_source'
    raise SummaryError(f'comparison {key!r} matches no known comparison type')


def flatten(obj, prefix, out, where):
    """Flatten numeric leaves of nested dicts; short string lists are joined."""
    for k in sorted(obj):
        v = obj[k]
        name = f'{prefix}.{k}' if prefix else k
        if isinstance(v, dict):
            flatten(v, name, out, where)
        elif is_num(v):
            out[name] = v
        elif v is None:
            out[name] = None  # undefined metric (e.g. precision with no calls); empty TSV cell
        elif isinstance(v, list) and all(isinstance(x, str) for x in v):
            out[name] = ';'.join(v)
        elif isinstance(v, str):
            continue  # descriptive text (policy); not a metric
        else:
            raise SummaryError(f'{where}: field {name!r} has unsupported type {type(v).__name__}')


def ec_fields(entry, where):
    keys_equal(entry, EC_KEYS, where)
    require(isinstance(entry['records'], list) and len(entry['records']) == entry['proteins'],
            f'{where}: records count differs from proteins')
    require(entry['proteins_with_candidates'] <= entry['proteins'], f'{where}: more proteins with candidates than proteins')
    require(sum(entry['reference_pairs_by_genome'].values()) == entry['reference_pairs'],
            f'{where}: per-genome reference pairs do not sum to reference_pairs')
    require(entry['reference_matches'] + entry['missed_reference'] == entry['reference_pairs'],
            f'{where}: matches + missed != reference pairs')
    body = {k: v for k, v in entry.items() if k not in ('records', 'policy')}
    out = {}
    flatten(body, '', out, where)
    return out


def cond_fields(entry, where):
    keys_equal(entry, COND_KEYS, where)
    out = {}
    flatten(entry, '', out, where)
    return out


def cmp_fields(entry, where):
    keys_equal(entry, CMP_KEYS, where)
    for g, v in entry['per_genome'].items():
        keys_equal(v, CMP_GENOME_KEYS, f'{where}/{g}')
    for k in CMP_GENOME_KEYS - {'proteins'}:
        require(sum(v[k] for v in entry['per_genome'].values()) == entry[k], f'{where}: per-genome {k} does not sum')
    require(sum(v['proteins'] for v in entry['per_genome'].values()) == entry['proteins'],
            f'{where}: per-genome proteins do not sum')
    out = {}
    flatten(entry, '', out, where)
    return out


def result_rows(results):
    """Yield dict rows in fixed order; validate the whole results object."""
    require(isinstance(results, dict), 'analysis has no results object')
    require(set(results) == set(OP_ORDER), f'operating points differ: {sorted(results)}')
    arms = arm_names() + [NATIVE]
    for op in ARM_OPERATING_POINTS:
        sec = results[op]
        keys_equal(sec, ARM_SECTIONS, op)
        ref = sec['full_ec_reference']
        cond = sec['conditional_candidate_metrics']
        tonat = sec['comparisons_to_native_emitted_calls']
        for name, sect in (('full_ec_reference', ref), ('conditional_candidate_metrics', cond),
                           ('comparisons_to_native_emitted_calls', tonat)):
            require(isinstance(sect, dict), f'{op}/{name}: not an object')
            keys_equal(sect, arms, f'{op}/{name}')
        for a in arm_names():
            s, c, k = split_arm(a)
            yield dict(operating_point=op, family='full_ec_reference', name=a, source=s, context=c, kind=k,
                       **{'_': ec_fields(ref[a], f'{op}/full_ec_reference/{a}')})
        yield dict(operating_point=op, family='full_ec_reference', name=NATIVE, source='clean_native',
                   context='', kind='native', **{'_': ec_fields(ref[NATIVE], f'{op}/ref/{NATIVE}')})
        for a in arm_names():
            s, c, k = split_arm(a)
            yield dict(operating_point=op, family='conditional_candidate_metrics', name=a, source=s, context=c,
                       kind=k, **{'_': cond_fields(cond[a], f'{op}/conditional/{a}')})
        yield dict(operating_point=op, family='conditional_candidate_metrics', name=NATIVE,
                   source='clean_native', context='', kind='native',
                   **{'_': cond_fields(cond[NATIVE], f'{op}/cond/{NATIVE}')})
        for a in arm_names():
            s, c, k = split_arm(a)
            yield dict(operating_point=op, family='comparison_to_native_emitted', name=a, source=s, context=c,
                       kind=k, baseline=NATIVE_OPERATING_POINT, comparison_type='vs_native_emitted_calls',
                       **{'_': cmp_fields(tonat[a], f'{op}/to_native/{a}')})
        yield dict(operating_point=op, family='comparison_to_native_emitted', name=NATIVE, source='clean_native',
                   context='', kind='native', baseline=NATIVE_OPERATING_POINT,
                   comparison_type='vs_native_emitted_calls', **{'_': cmp_fields(tonat[NATIVE], f'{op}/to_native/{NATIVE}')})
        for key in sorted(sec['comparisons']):
            subject, baseline, ctype = classify(key)
            s, c, k = split_arm(subject)
            yield dict(operating_point=op, family='comparison', name=key, source=s, context=c, kind=k,
                       baseline=baseline, comparison_type=ctype,
                       **{'_': cmp_fields(sec['comparisons'][key], f'{op}/{key}')})
    nat = results[NATIVE_OPERATING_POINT]
    keys_equal(nat, {'conditional_candidate_metrics', 'full_ec_reference'}, NATIVE_OPERATING_POINT)
    yield dict(operating_point=NATIVE_OPERATING_POINT, family='full_ec_reference', name=NATIVE,
               source='clean_native', context='', kind='native',
               **{'_': ec_fields(nat['full_ec_reference'], 'native/full_ec_reference')})
    yield dict(operating_point=NATIVE_OPERATING_POINT, family='conditional_candidate_metrics', name=NATIVE,
               source='clean_native', context='', kind='native',
               **{'_': cond_fields(nat['conditional_candidate_metrics'], 'native/conditional')})


def build_result_table(results):
    rows = []
    for r in result_rows(results):
        flat = r.pop('_')
        for k in LEAD:
            r.setdefault(k, '')
        for k, v in flat.items():
            require(k not in LEAD, f'metric field {k!r} collides with a label column')
            r[k] = v
        rows.append(r)
    seen = set()
    for r in rows:
        key = (r['operating_point'], r['family'], r['name'])
        require(key not in seen, f'duplicate row {key}')
        seen.add(key)
    metric_cols = sorted({k for r in rows for k in r if k not in LEAD})
    return LEAD + metric_cols, rows


def fmt_cell(v):
    if v == '' or v is None:
        return ''
    return repr(v) if isinstance(v, float) else str(v)


def write_tsv(path, columns, rows):
    with open(path, 'w', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh, delimiter='\t', lineterminator='\n')
        w.writerow(columns)
        for r in rows:
            w.writerow([fmt_cell(r.get(c, '')) for c in columns])


# ---- replay integrity ------------------------------------------------------------------

def integrity_rows(verification):
    require(isinstance(verification, dict), 'verification is not an object')
    require(verification.get('status') == 'verified', f'replay status is {verification.get("status")!r}, not verified')
    require(verification.get('schema') == 'gspa.source-context-production/1', 'unexpected replay schema')
    summ = verification.get('summaries')
    require(isinstance(summ, dict) and set(summ) == set(SOURCES), 'replay summaries differ from the five sources')
    rows = []
    genomes_by_source = {}
    for s in SOURCES:
        require(set(summ[s]) == set(CONTEXTS), f'{s}: replay contexts differ')
        for c in CONTEXTS:
            require(set(summ[s][c]) == set(REPLAY_KINDS), f'{s}/{c}: replay kinds differ')
            for k in REPLAY_KINDS:
                cells = summ[s][c][k]
                gset = frozenset(cells)
                require(gset, f'{s}/{c}/{k}: no genomes')
                genomes_by_source.setdefault(s, gset)
                require(genomes_by_source[s] == gset, f'{s}/{c}/{k}: genome set differs within source')
                tot = dict(genomes=0, candidates=0, functions=0, added=0, removed=0, rescored=0,
                           unchanged_unsupported=0, scored_only=0, display_ties=0, cells_with_tie_record=0,
                           max_probability_error=0.0)
                reasons = {}
                for g in sorted(cells):
                    cell = cells[g]
                    where = f'{s}/{c}/{k}/{g}'
                    has_ties = REPLAY_TIE_KEYS <= set(cell)
                    keys_equal(cell, REPLAY_CELL_KEYS | (REPLAY_TIE_KEYS if has_ties else set()), where)
                    a = cell['actions']
                    require(isinstance(a, dict) and {'rescored'} <= set(a)
                            <= {'rescored', 'unchanged_unsupported', 'scored_only'},
                            f'{where}: unexpected actions {sorted(a)}')
                    scored_only = a.get('scored_only', 0)
                    require(sum(a.values()) == cell['candidates'], f'{where}: actions do not sum to candidates')
                    if has_ties:
                        require(isinstance(cell['display_tie_candidates'], list)
                                and cell['display_ties'] == len(cell['display_tie_candidates']),
                                f'{where}: display_ties disagrees with display_tie_candidates')
                        tot['display_ties'] += cell['display_ties']
                        tot['cells_with_tie_record'] += 1
                    require(sum(cell['reasons'].values()) == cell['candidates'],
                            f'{where}: reasons do not sum to candidates')
                    require(cell['reasons'].get('scored') == a['rescored'] + scored_only, f'{where}: scored != rescored + scored_only')
                    require(is_num(cell['max_probability_error']) and cell['max_probability_error'] >= 0,
                            f'{where}: bad max_probability_error')
                    tot['genomes'] += 1
                    for f in ('candidates', 'functions', 'added', 'removed'):
                        tot[f] += cell[f]
                    tot['rescored'] += a['rescored']
                    tot['scored_only'] += scored_only
                    tot['unchanged_unsupported'] += a.get('unchanged_unsupported', 0)
                    tot['max_probability_error'] = max(tot['max_probability_error'], cell['max_probability_error'])
                    for r, n in cell['reasons'].items():
                        if r != 'scored':
                            reasons[r] = reasons.get(r, 0) + n
                rows.append(dict(source=s, context=c, kind=k, **tot,
                                 reasons=';'.join(f'{r}={reasons[r]}' for r in sorted(reasons)), _reasons=reasons))
    return rows


def source_totals(rows):
    out = {}
    for s in SOURCES:
        sel = [r for r in rows if r['source'] == s]
        t = dict(genomes=sel[0]['genomes'], cells=len(sel), candidates=0, rescored=0, scored_only=0, unchanged_unsupported=0,
                 display_ties=0, cells_with_tie_record=0, added=0, removed=0, max_probability_error=0.0, reasons={})
        for r in sel:
            require(r['genomes'] == t['genomes'], f'{s}: genome count differs between cells')
            for f in ('candidates', 'rescored', 'scored_only', 'unchanged_unsupported', 'added', 'removed',
                      'display_ties', 'cells_with_tie_record'):
                t[f] += r[f]
            t['max_probability_error'] = max(t['max_probability_error'], r['max_probability_error'])
            for k, n in r['_reasons'].items():
                t['reasons'][k] = t['reasons'].get(k, 0) + n
        out[s] = t
    return out


# ---- markdown --------------------------------------------------------------------------

def f3(v):
    if v is None:
        return 'n/a'
    return f'{v:.3f}' if isinstance(v, float) else str(v)


def table(header, body):
    lines = ['| ' + ' | '.join(header) + ' |', '|' + '|'.join(['---'] * len(header)) + '|']
    lines += ['| ' + ' | '.join(str(x) for x in r) + ' |' for r in body]
    return lines


def render(results, rows_by, integrity, totals):
    md = ['# Source x context production replay: summary', '',
          'Generated by `benchmark/v3/summarize_source_context_replay.py` from `analysis/results.json` and '
          '`replay/verification.json`. The replay verifies execution, input preservation and Python/JVM score '
          'agreement; it does not measure biological accuracy. The reference is the full experimental EC '
          'population with exact identifiers; unmatched calls are unconfirmed, never negatives.', '',
          '## Replay integrity per source', '',
          'Totals over 4 contexts x 2 learned kinds (logistic, gbdt) x genomes; candidates are replayed '
          'candidate functions.', '']
    body = []
    for s in SOURCES:
        t = totals[s]
        body.append([f'`{s}`', t['genomes'], t['cells'], t['candidates'], t['rescored'], t['scored_only'],
                     t['unchanged_unsupported'], t['added'], t['removed'], t['display_ties'],
                     f'{t["max_probability_error"]:.2e}'])
    md += table(['source', 'genomes', 'cells', 'candidates', 'rescored', 'scored only', 'unchanged unsupported',
                 'added', 'removed', 'display ties', 'max abs Python-JVM error'], body)
    md += ['', '### Unsupported reasons (candidates left unchanged)', '']
    reasons = sorted({r for t in totals.values() for r in t['reasons']})
    md += table(['source'] + [f'`{r}`' for r in reasons],
                [[f'`{s}`'] + [totals[s]['reasons'].get(r, 0) for r in reasons] for s in SOURCES])
    md += ['']
    for op in OP_ORDER:
        md += [f'## Operating point `{op}`', '']
        if op == NATIVE_OPERATING_POINT:
            nat = results[op]
            md += ['Native CLEAN emitted calls, unchanged confidence; every emitted call is kept.', '']
            md += table(['arm', 'proteins', 'reference pairs', 'matched', 'precision', 'recall', 'micro-F1',
                         'top-1', 'top-3', 'top-5', 'unconfirmed calls'], [ec_row(NATIVE, nat['full_ec_reference'])])
            md += ['', 'Conditional candidate table:', '']
            c = nat['conditional_candidate_metrics']
            md += table(['rows', 'positive', 'closed-world unconfirmed', 'F1', 'matched calls', 'missed reference',
                         'Brier'], [[c['rows'], c['positive'], c['label_basis']['closed_world_unconfirmed'],
                                     f3(c['reference_f1']), c['reference_matched_calls'], c['missed_reference'],
                                     f3(c['brier_on_available_scores'])]])
            md += ['']
            continue
        sec = results[op]
        ref = sec['full_ec_reference']
        md += ['Full experimental EC reference, per arm. Matched/pairs counts exact reference EC assertions '
               'recovered; top-k counts proteins whose ranked list hits the reference; denominators are '
               'in the columns.', '']
        body = []
        for a in arm_names():
            e = ref[a]
            body.append(ec_row(a, e))
        body.append(ec_row(NATIVE, ref[NATIVE]))
        md += table(['arm', 'proteins', 'reference pairs', 'matched', 'precision', 'recall', 'micro-F1',
                     'top-1', 'top-3', 'top-5', 'unconfirmed calls'], body)
        for ctype, title in (('learned_vs_rule', 'Learned versus rule (same source and context)'),
                             ('context_vs_seqonly', 'Context versus SEQONLY (same source and kind)'),
                             ('vs_clean_native', 'Arm versus native CLEAN at the same operating point'),
                             ('source_vs_source', 'Source versus source (same context and kind)')):
            md += ['', f'### {title}', '',
                   'Paired over the proteins shown; gained/lost are exact reference EC calls, top-1 counts '
                   'proteins, unconfirmed change is the net change in unmatched calls.', '']
            body = []
            for r in rows_by[op]:
                if r['family'] == 'comparison' and r['comparison_type'] == ctype:
                    body.append([f'`{r["name"].replace("_vs_", "` vs `")}`', r['proteins'], r['reference_calls_gained'],
                                 r['reference_calls_lost'], r['top1_gained'], r['top1_lost'],
                                 r['unconfirmed_call_change']])
            md += table(['comparison', 'proteins', 'reference calls gained', 'lost', 'top-1 gained', 'top-1 lost',
                         'unconfirmed change'], body)
        md += ['', '### Arm versus native CLEAN emitted calls (no threshold)', '',
               'Baseline is every call native CLEAN emits at its unchanged confidence, '
               'including zero confidence.', '']
        body = []
        for r in rows_by[op]:
            if r['family'] == 'comparison_to_native_emitted':
                body.append([f'`{r["name"]}`', r['proteins'], r['reference_calls_gained'], r['reference_calls_lost'],
                             r['top1_gained'], r['top1_lost'], r['unconfirmed_call_change']])
        md += table(['arm', 'proteins', 'reference calls gained', 'lost', 'top-1 gained', 'top-1 lost',
                     'unconfirmed change'], body)
        md += ['']
    return '\n'.join(md).rstrip('\n') + '\n'


def ec_row(name, e):
    p = e['proteins']
    return [f'`{name}`', p, e['reference_pairs'], e['reference_matches'], f3(e['reference_precision']),
            f3(e['reference_recall']), f3(e['reference_micro_f1']), f'{e["top1_matches"]}/{p}',
            f'{e["top3_matches"]}/{p}', f'{e["top5_matches"]}/{p}', e['unconfirmed_calls']]


def run(results_dir, out_dir):
    results_dir, out_dir = Path(results_dir), Path(out_dir)
    analysis = json.loads((results_dir / 'analysis' / 'results.json').read_text(encoding='utf-8'))
    verification = json.loads((results_dir / 'replay' / 'verification.json').read_text(encoding='utf-8'))
    keys_equal(analysis, {'policy', 'results'}, 'analysis')
    require(analysis['policy'].get('schema') == 'gspa.source-context-analysis/1', 'unexpected analysis schema')
    columns, rows = build_result_table(analysis['results'])
    irows = integrity_rows(verification)
    totals = source_totals(irows)
    rows_by = {}
    for r in rows:
        rows_by.setdefault(r['operating_point'], []).append(r)
    md = render(analysis['results'], rows_by, irows, totals)
    out_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(out_dir / RESULTS_TSV, columns, rows)
    icols = ['source', 'context', 'kind', 'genomes', 'candidates', 'functions', 'added', 'removed', 'rescored',
             'unchanged_unsupported', 'scored_only', 'display_ties', 'cells_with_tie_record',
             'max_probability_error', 'reasons']
    write_tsv(out_dir / INTEGRITY_TSV, icols, irows)
    (out_dir / MD_NAME).write_text(md, encoding='utf-8')
    return len(rows), len(irows)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--results-dir', required=True, help='directory holding analysis/ and replay/')
    ap.add_argument('--out-dir', required=True)
    a = ap.parse_args(argv)
    try:
        n, m = run(a.results_dir, a.out_dir)
    except SummaryError as e:
        raise SystemExit(f'error: {e}')
    print(f'wrote {n} result rows and {m} integrity rows to {a.out_dir}')


if __name__ == '__main__':
    main()
