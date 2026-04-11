#!/usr/bin/env python3
"""Compile final benchmark tables from bench9 results."""

import argparse
import json
import os


def fmt_f(v, ci=None):
    if ci and (ci[0] or ci[1]):
        return f'{v:.3f} [{ci[0]:.3f}, {ci[1]:.3f}]'
    return f'{v:.3f}'


def load_json(path):
    try:
        with open(path) as fh:
            return json.load(fh)
    except Exception:
        return None


def print_table(rows, title):
    print(f'### {title}\n')
    if not rows or len(rows) < 2:
        print('(no data)\n')
        return
    widths = [max(len(str(r[i])) for r in rows) for i in range(len(rows[0]))]
    for i, r in enumerate(rows):
        print('  '.join(str(v).ljust(widths[j]) for j, v in enumerate(r)))
        if i == 0:
            print('-' * (sum(widths) + 2 * (len(widths) - 1)))
    print()


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--bench-dir', required=True)
    p.add_argument('--gaef-dir', required=True)
    p.add_argument('--genomes', nargs='+', required=True)
    p.add_argument('--strip-dir', default=None)
    args = p.parse_args()

    print('# GSPA 9-Genome Benchmark Results\n')
    print(f'generated: $(date)')
    print(f'bench: {args.bench_dir}')
    print(f'gaef:  {args.gaef_dir}')
    print(f'genomes: {", ".join(args.genomes)}\n')

    # ===== F-max vs experimental-only and all-GOA truth =====
    for truth_name in ('exp', 'all'):
        rows = [['genome', 'method', 'F-max (95% CI)', 'MF', 'BP', 'CC', 'IC-recall', 'coverage']]
        for tag in args.genomes:
            d = load_json(os.path.join(args.bench_dir, f'{tag}_fmax.json'))
            if not d:
                continue
            block = d.get('by_truth', {}).get(truth_name)
            if not block:
                continue
            for r in block.get('results', []):
                rows.append([
                    tag, r['method'],
                    fmt_f(r['fmax_overall'], r.get('fmax_ci')),
                    fmt_f(r['fmax_MF']),
                    fmt_f(r['fmax_BP']),
                    fmt_f(r['fmax_CC']),
                    fmt_f(r.get('ic_recall', 0)),
                    fmt_f(r.get('coverage', 0)),
                ])
        label = ('Experimental-only truth (EXP/IDA/IMP/...)' if truth_name == 'exp'
                 else 'Full-GOA truth (all evidence incl. IEA)')
        print_table(rows, f'F-max — {label}')

    # ===== GAEF (Completeness / Coherence / Consistency) =====
    gaef_rows = [['genome', 'method', 'completeness', 'proc_coh', 'path_coh', 'cplx_coh', 'consist', 'IC', 'composite']]
    for tag in args.genomes:
        for method in ('gspa', 'pgap'):
            path = os.path.join(args.gaef_dir, f'{tag}_{method}_quality.json')
            d = load_json(path)
            if not d:
                continue
            comp = d.get('completeness', {}) or {}
            coh = d.get('coherence', {}) or {}
            cons = d.get('consistency', {}) or {}
            ic = d.get('information_content', {}) or {}
            summary = d.get('summary', {}) or {}

            def g(x, default='-'):
                return (f'{x:.3f}' if isinstance(x, (int, float)) and x >= 0 else default)

            gaef_rows.append([
                tag, method.upper(),
                g(comp.get('score') if isinstance(comp, dict) else comp),
                g(coh.get('process_coherence')),
                g(coh.get('pathway_coherence')),
                g(coh.get('complex_coherence')),
                str(cons.get('consistent', '-')),
                g(ic.get('mean_ic')),
                g(summary.get('composite_score') or d.get('composite_score')),
            ])
    print_table(gaef_rows, 'GAEF metrics (Completeness / Coherence / Consistency)')

    # ===== Strip test =====
    if args.strip_dir:
        strip_rows = [['genome', 'n_stripped', 'singleton_rate', 'any_suggestion_rate', 'survived_rate', 'pgap_floor']]
        for tag in args.genomes:
            d = load_json(os.path.join(args.strip_dir, f'{tag}_strip_report.json'))
            if not d:
                continue
            s = d.get('summary', {})
            strip_rows.append([
                tag, s.get('n_stripped', '-'),
                fmt_f(s.get('singleton_recovery_rate', 0)),
                fmt_f(s.get('any_suggestion_recovery_rate', 0)),
                fmt_f(s.get('integrated_survived_recovery_rate', 0)),
                fmt_f(s.get('pgap_baseline_recovery_rate', 0)),
            ])
        print_table(strip_rows, 'Dark-matter strip test (Phase 8 suggester)')


if __name__ == '__main__':
    main()
