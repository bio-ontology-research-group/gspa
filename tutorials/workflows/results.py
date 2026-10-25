#!/usr/bin/env python3
"""Convert real annotation output and summarize it; no invented predictor scores."""
import csv
import html
import json
from pathlib import Path
import shutil
import sys
from collections import Counter
from data import fasta

ROOT = Path('/work') if Path('/work').is_dir() else Path('.tutorial-work')


def annotations(case):
    paths = list((ROOT / case / 'annot').glob('*_annotations.tsv'))
    assert len(paths) == 1, paths
    return list(csv.DictReader(paths[0].open(), delimiter='\t'))


def claims(case):
    rows = annotations(case)
    out = ROOT / case / 'claims.jsonl'
    count = 0
    seen = set()
    with out.open('w') as handle:
        for row in rows:
            if row['type'] not in ('GO', 'EC') or row['source'] != 'diamond':
                continue
            key = row['protein_id'], row['type'], row['value'], row['source']
            if key in seen:
                continue
            seen.add(key)
            # DIAMOND currently exposes identity/100, not a fitted calibration.
            score = min(.99, max(.01, float(row['score'])))
            handle.write(json.dumps({'protein_id': row['protein_id'], 'function_type': row['type'],
                'function_id': row['value'], 'source': row['source'],
                'raw_score': float(row['score']), 'calibrated_prob': score}) + '\n')
            count += 1
    assert count > 0, 'No functional claims: inspect the predictor log/database'
    print(f'{count:,} unique DIAMOND GO/EC claims from {len({r[0] for r in seen}):,} proteins')
    print('Teaching score: DIAMOND identity fraction, clipped to [0.01, 0.99]; not a fitted probability.')


def features(case):
    gff = next((ROOT / case / 'annot').glob('*.gff3'))
    result = []
    for line in gff.read_text().splitlines():
        if line == '##FASTA': break
        if line.startswith('#'): continue
        f = line.split('\t')
        if len(f) != 9 or f[2] != 'CDS': continue
        attrs = dict(p.split('=',1) for p in f[8].split(';') if '=' in p)
        result.append((attrs.get('locus_tag', attrs.get('ID')), f))
    return result


def summary(case):
    rows = list(csv.DictReader((ROOT / case / 'gspa_out/integrated.tsv').open(), delimiter='\t'))
    assert rows, 'No integrated output'
    raw = annotations(case)
    fts = features(case)
    n = len(fts) or len(list(fasta(ROOT / 'data/eukaryote.faa')))
    annotated = len({r['protein_id'] for r in raw if r['type'] in ('GO','EC') and r['source'] == 'diamond'})
    counts = Counter(r['type'] for r in rows)
    priors = Counter(p.split(':')[0] for r in rows for p in r['priors_fired'].split(';') if p)
    stat = {'case': case, 'proteins': n, 'with_go_ec': annotated,
            'go_rows': counts['GO'], 'ec_rows': counts['EC'],
            'accepted_rows': sum(float(r['posterior_prob']) >= .5 for r in rows),
            'priors': dict(sorted(priors.items()))}
    (ROOT / case / 'summary.json').write_text(json.dumps(stat, indent=2) + '\n')
    print(f'Proteins: {n:,}; with GO/EC search evidence: {annotated:,}')
    print(f'Integrated hypotheses: GO={counts["GO"]:,}; EC={counts["EC"]:,}')
    print(f'Hypotheses at posterior >= 0.5: {stat["accepted_rows"]:,}')
    print('These counts describe this teaching reference; they are not accuracy or completeness estimates.')
    chart(case, stat)


def chart(case, stat):
    # Standalone vector figure drawn from the actual output counts.
    total = max(stat['proteins'],1)
    covered = stat['with_go_ec'] / total
    title = {'bacteria':'Bacterial genome','archaea':'Archaeal genome','eukaryote':'Eukaryotic proteome',
             'mag':'Incomplete archaeal bin','metagenome':'Mixed community','virus':'Lambda phage',
             'prokaryote':'Prokaryote, domain unspecified'}.get(case, case)
    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" width="900" height="320" viewBox="0 0 900 320">
<rect width="900" height="320" rx="18" fill="#101f33"/>
<g font-family="sans-serif" fill="#e8f0f8"><text x="36" y="42" font-size="13" fill="#76d9cd">GSPA / EXECUTED WORKFLOW</text>
<text x="36" y="82" font-size="29">{html.escape(title)}</text>
<text x="36" y="126" font-size="16">{stat['proteins']:,} proteins • {stat['with_go_ec']:,} with GO/EC search evidence</text>
<rect x="36" y="148" width="828" height="24" rx="6" fill="#34465d"/>
<rect x="36" y="148" width="{828*covered:.2f}" height="24" rx="6" fill="#57c7bb"/>
<text x="36" y="214" font-size="26">{stat['go_rows']:,} GO</text><text x="300" y="214" font-size="26">{stat['ec_rows']:,} EC</text>
<text x="575" y="214" font-size="26">{stat['accepted_rows']:,} ≥ 0.5</text>
<text x="36" y="242" font-size="13" fill="#aebed0">integrated hypotheses</text><text x="575" y="242" font-size="13" fill="#aebed0">posterior threshold</text>
<text x="36" y="288" font-size="13" fill="#aebed0">Teaching reference includes close/self matches. This is an execution demonstration, not an accuracy test.</text></g></svg>'''
    (ROOT / case / 'summary.svg').write_text(svg)


def browser(case):
    base = ROOT / case
    folder = base / 'prokka_out'
    folder.mkdir(exist_ok=True)
    fts = features(case)
    with (folder / f'{case}.tsv').open('w') as out:
        out.write('locus_tag\tftype\tlength_bp\tgene\tEC_number\tCOG\tproduct\n')
        if fts:
            for pid, f in fts:
                out.write(f'{pid}\tCDS\t{int(f[4])-int(f[3])+1}\t\t\t\tProdigal CDS\n')
        else:
            for pid, seq in fasta(ROOT / 'data/eukaryote.faa'):
                out.write(f'{pid}\tCDS\t{len(seq)*3}\t\t\t\tSupplied protein; genomic coordinates unavailable\n')
    with (folder / f'{case}.gff').open('w') as out:
        out.write('##gff-version 3\n')
        for pid,f in fts:
            f[8] += ';locus_tag=' + pid
            out.write('\t'.join(f) + '\n')
    operons = base / 'annot/operons.tsv'
    if operons.exists(): shutil.copyfile(operons, base / 'gspa_out/operons.tsv')
    # Explicit empty lookups avoid machine-specific defaults. GO IDs remain visible.
    (ROOT / 'data/no-labels.obo').write_text('format-version: 1.2\n')
    (ROOT / 'data/no-ec-labels.txt').write_text('! No EC labels installed in this compact tutorial\n')
    print('Browser metadata prepared from real CDS coordinates / supplied protein identifiers.')


if __name__ == '__main__':
    {'claims': claims, 'summary': summary, 'browser': browser}[sys.argv[1]](sys.argv[2])
