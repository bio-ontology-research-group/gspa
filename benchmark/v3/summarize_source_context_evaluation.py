"""Tabulate the source-arm x context evaluation written by
run_source_arm_context_comparison.py evaluate.

Writes one TSV row per (source, context, kind) with every metric and row count,
and a Markdown summary whose metric cells carry supported/total test rows.
Metrics are on supported rows only (the JSON's metric_scope); counts are shown
so that cells with tiny supported counts are not read as comparable. Stdlib only,
deterministic output; fails on any missing or unexpected cell.
"""
import argparse
import csv
import json
from pathlib import Path

SOURCES = ('diamond_homology', 'diamond_family_excluded', 'clean', 'clean_diamond_homology',
           'clean_diamond_family_excluded')
CONTEXTS = ('SEQONLY', 'UNIVERSAL', 'PATHWAY', 'ORGANISM')
KINDS = ('logistic', 'gbdt')
NO_ROWS = 'no_supported_test_rows'
TSV_NAME = 'source_context_evaluation.tsv'
MD_NAME = 'source_context_evaluation.md'
COLUMNS = ['source', 'context', 'kind', 'status', 'total_rows', 'supported_rows', 'unsupported_rows',
           'unsupported_reasons', 'locus_groups_with_unsupported_members', 'n', 'positive_rate',
           'n_positive_supported', 'auroc', 'logloss', 'brier', 'ece_10bin', 'display_threshold',
           'n_predicted_positive', 'precision', 'recall', 'recall_all_test_positives', 'locus_trained',
           'locus_reason', 'probability_target', 'evaluated_on', 'metric_scope', 'model_sha256',
           'source_file_sha256']


class SummaryError(ValueError):
    pass


def require(condition, message):
    if not condition:
        raise SummaryError(message)


def cells(evaluation):
    """Yield (source, context, kind, result) in fixed order; fail on gaps or extras."""
    results = evaluation.get('results')
    require(isinstance(results, dict), 'evaluation has no results object')
    require(set(results) == set(SOURCES), f'sources differ: missing {sorted(set(SOURCES) - set(results))}, '
            f'unexpected {sorted(set(results) - set(SOURCES))}')
    for s in SOURCES:
        require(set(results[s]) == set(CONTEXTS), f'{s}: contexts differ: missing {sorted(set(CONTEXTS) - set(results[s]))}, '
                f'unexpected {sorted(set(results[s]) - set(CONTEXTS))}')
        for c in CONTEXTS:
            require(set(results[s][c]) == set(KINDS), f'{s}/{c}: kinds differ: missing {sorted(set(KINDS) - set(results[s][c]))}, '
                    f'unexpected {sorted(set(results[s][c]) - set(KINDS))}')
            for k in KINDS:
                yield s, c, k, results[s][c][k]


def row(s, c, k, r):
    where = f'{s}/{c}/{k}'
    cov = r.get('input_coverage')
    require(isinstance(cov, dict), f'{where}: missing input_coverage')
    for key in ('total_rows', 'supported_rows', 'unsupported_reasons'):
        require(key in cov, f'{where}: input_coverage lacks {key}')
    total, supported, reasons = cov['total_rows'], cov['supported_rows'], cov['unsupported_reasons']
    require(total - supported == sum(reasons.values()), f'{where}: unsupported reasons do not sum to total - supported')
    empty = r.get('status') == NO_ROWS
    require(empty == (supported == 0), f'{where}: status {r.get("status")!r} inconsistent with {supported} supported rows')
    for key in ('auroc', 'logloss', 'brier', 'ece_10bin', 'locus', 'evaluated_on'):
        require(key in r, f'{where}: missing {key}')
    thr = r.get('at_display_threshold')
    if not empty:
        require(isinstance(thr, dict), f'{where}: missing at_display_threshold')
        for key in ('threshold', 'precision', 'recall', 'n_predicted_positive'):
            require(key in thr, f'{where}: at_display_threshold lacks {key}')
        require(r.get('n') == supported, f'{where}: n {r.get("n")} differs from supported_rows {supported}')
    thr = thr or {}
    rate = r.get('positive_rate')
    locus = r['locus']
    return dict(source=s, context=c, kind=k, status=r.get('status', 'evaluated'), total_rows=total,
                supported_rows=supported, unsupported_rows=total - supported,
                unsupported_reasons=';'.join(f'{name}={reasons[name]}' for name in sorted(reasons)),
                locus_groups_with_unsupported_members=cov.get('locus_groups_with_unsupported_members'),
                n=r.get('n'), positive_rate=rate,
                n_positive_supported=None if rate is None else round(rate * r['n']),
                auroc=r['auroc'], logloss=r['logloss'], brier=r['brier'], ece_10bin=r['ece_10bin'],
                display_threshold=thr.get('threshold'), n_predicted_positive=thr.get('n_predicted_positive'),
                precision=thr.get('precision'), recall=thr.get('recall'),
                recall_all_test_positives=thr.get('recall_all_test_positives'),
                locus_trained=locus.get('trained', 'n_groups' in locus), locus_reason=locus.get('reason'),
                probability_target=r.get('probability_target'), evaluated_on=r['evaluated_on'],
                metric_scope=r.get('metric_scope'), model_sha256=r.get('model_sha256'),
                source_file_sha256=r.get('source_file_sha256'), _reasons=reasons)


def tsv_value(v):
    return '' if v is None else repr(v) if isinstance(v, float) else str(v)


def num(v):
    return 'n/a' if v is None else f'{v:.3f}'


def counts(x):
    return f'{x["supported_rows"]}/{x["total_rows"]}'


def grid(rows, kind, fmt):
    by = {(x['source'], x['context']): x for x in rows if x['kind'] == kind}
    out = ['| source | ' + ' | '.join(CONTEXTS) + ' |', '|---|' + '---|' * len(CONTEXTS)]
    out += ['| ' + s + ' | ' + ' | '.join(fmt(by[s, c]) for c in CONTEXTS) + ' |' for s in SOURCES]
    return out


def markdown(evaluation, rows, input_name):
    scopes = sorted({x['metric_scope'] for x in rows if x['metric_scope']})
    tables = [('AUROC', lambda x: f'{num(x["auroc"])} ({counts(x)})'),
              ('Log-loss', lambda x: f'{num(x["logloss"])} ({counts(x)})'),
              ('Precision / recall at the display threshold',
               lambda x: f'{num(x["precision"])} / {num(x["recall"])} '
                         f'({x["n_predicted_positive"] if x["n_predicted_positive"] is not None else "n/a"} pred.+; {counts(x)})')]
    md = [f'# Source x context evaluation: {input_name}', '',
          'Cells: metric (supported/total test rows). Precision/recall cells also give the number of supported rows '
          'predicted positive at the frozen display threshold. Metric scope as recorded: '
          + '; '.join(f'"{t}"' for t in scopes) + '.', '']
    for title, fmt in tables:
        for kind in KINDS:
            md += [f'## {title}: {kind}', ''] + grid(rows, kind, fmt) + ['']
    reasons = sorted({name for x in rows for name in x['_reasons']})
    md += ['## Unsupported test rows by reason', '',
           '| source | context | kind | total | supported | ' + ' | '.join(f'`{n}`' for n in reasons) + ' |',
           '|---|---|---|---|---|' + '---|' * len(reasons)]
    md += [f'| {x["source"]} | {x["context"]} | {x["kind"]} | {x["total_rows"]} | {x["supported_rows"]} | '
           + ' | '.join(str(x['_reasons'].get(n, 0)) for n in reasons) + ' |' for x in rows]
    pending = evaluation.get('pending_stages') or {}
    md += ['', '## Limitations', '', 'Quoted from the evaluation JSON.', '',
           f'- scope: "{evaluation.get("scope", "(absent)")}"',
           f'- limitations: "{evaluation.get("limitations", "(absent)")}"',
           '- pending_stages: ' + (', '.join(f'`{k}`: {pending[k]}' for k in sorted(pending)) or '(none)'), '']
    return '\n'.join(md)


def summarize(evaluation_path, out_dir):
    evaluation_path, out_dir = Path(evaluation_path), Path(out_dir)
    evaluation = json.loads(evaluation_path.read_text())
    rows = [row(*cell) for cell in cells(evaluation)]
    out_dir.mkdir(parents=True, exist_ok=True)
    with open(out_dir / TSV_NAME, 'w', newline='') as h:
        w = csv.writer(h, delimiter='\t', lineterminator='\n')
        w.writerow(COLUMNS)
        w.writerows([tsv_value(x[col]) for col in COLUMNS] for x in rows)
    (out_dir / MD_NAME).write_text(markdown(evaluation, rows, evaluation_path.name))
    return rows


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--evaluation', type=Path, required=True)
    ap.add_argument('--out-dir', type=Path, required=True)
    a = ap.parse_args()
    summarize(a.evaluation, a.out_dir)
