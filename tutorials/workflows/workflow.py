#!/usr/bin/env python3
"""The actual GSPA commands behind the biological lessons (inside the container)."""
import argparse
import csv
from pathlib import Path
import shlex
import subprocess

import results

CASES = ('bacteria', 'archaea', 'eukaryote', 'mag', 'metagenome', 'virus', 'prokaryote')


def commands(case):
    """Keep biological assumptions visible and reusable by the command preview."""
    base = f'/work/{case}'
    source = 'bacteria' if case == 'prokaryote' else case
    supplied = case == 'eukaryote'
    annotate = ['gspa', 'annotate']
    annotate += ['-k', 'virus', '--mag'] if case == 'virus' else ['--preset', case]
    if case == 'mag':
        annotate += ['-k', 'archaea']
    annotate += ['--proteins' if supplied else '--genome',
                 f'/work/data/{source}.{"faa" if supplied else "fna"}',
                 '-c', f'/work/data/{"virus" if case == "virus" else "cpu"}.json',
                 '--db', '/work/data/reference.dmnd', '-o', f'{base}/annot', '-t', '2']
    integrate = ['gspa', 'integrate']
    integrate += (['--lite', '--enable-priors', '', '--no-enforce-stage2']
                  if case == 'virus' else ['--preset', case])
    if case == 'mag':
        integrate += ['--taxon', 'archaea']
    if case in ('bacteria', 'archaea', 'mag', 'prokaryote'):
        integrate += ['--operons', f'{base}/annot/operons_for_integrate.tsv']
    integrate += ['--claims', f'{base}/claims.jsonl', '--out', f'{base}/gspa_out/integrated.tsv',
                  '--provenance', f'{base}/gspa_out/provenance.json',
                  '--enforcement-actions-out', f'{base}/gspa_out/actions.tsv']
    visualize = ['gspa', 'visualize', '-w', base, '-g', case,
                 '--go-obo', '/work/data/browser-labels.obo', '--ec2go', '/work/data/no-ec-labels.txt']
    if not supplied:
        visualize += ['--fasta', f'/work/data/{source}.fna']
    return [('annotate', annotate), ('integrate', integrate), ('visualize', visualize)]


def logged(command, log):
    with log.open('w') as handle:
        handle.write('$ ' + shlex.join(command) + '\n')
        handle.flush()
        result = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT)
    if result.returncode:
        print('\n'.join(log.read_text(errors='replace').splitlines()[-60:]), flush=True)
        raise SystemExit(f'Step failed. Full log: {log}')


def run(case):
    base = Path('/work') / case
    (base / 'gspa_out').mkdir(parents=True, exist_ok=True)
    # Do not leave a previous successful browser/summary looking like this run's output.
    for path in (base / 'summary.json', base / 'gspa_out' / f'{case}_browser.html',
                 base / 'gspa_out/actions.tsv'):
        path.unlink(missing_ok=True)
    for step, command in commands(case):
        print({'annotate': '1/3 Finding sequence evidence',
               'integrate': '2/3 Applying the biological assumptions',
               'visualize': '3/3 Preparing the results browser'}[step], flush=True)
        if step == 'visualize':
            results.browser(case)
        logged(command, base / f'{step}.log')
        if step == 'annotate':
            results.claims(case)
        if step == 'integrate':
            results.summary(case)
    browser = base / 'gspa_out' / f'{case}_browser.html'
    if not browser.is_file() or browser.stat().st_size == 0:
        raise SystemExit('The results browser was not generated; inspect visualize.log.')
    results.inspect(case)


def neural(device):
    import neural as helper
    helper.device(device)
    base = Path('/work/neural')
    if device == 'cpu':
        helper.prepare()
        logged(['python3', '/src/benchmark/neural/build_esm2_centroids.py',
                '--swissprot-fasta', str(base / 'reference.faa'),
                '--swissprot-tsv', str(base / 'labels.tsv'),
                '--model', 'esm2_t12_35M_UR50D', '--batch-size', '4', '--min-class-size', '2',
                '--out', str(base / 'centroids.npz')], base / 'build-centroids.log')
    logged(['gspa', 'annotate', '--preset', 'metagenome', '--proteins', str(base / 'query.faa'),
            '-c', str(base / 'config.json'), '-o', str(base / device), '-t', '2'], base / f'{device}.log')
    path = next((base / device).glob('*_annotations.tsv'))
    with path.open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    if len(rows) != 48 or {r['source'] for r in rows} != {'esm2-centroid'}:
        raise SystemExit(f'Unexpected neural results; inspect {path}')
    print(f'{device.upper()}: 48 ESM-2 candidate annotations for 16 proteins.')
    if device == 'cuda':
        helper.compare()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('case', choices=CASES + ('neural-cpu', 'neural-cuda'))
    args = parser.parse_args()
    if args.case.startswith('neural-'):
        neural(args.case.removeprefix('neural-'))
    else:
        run(args.case)
