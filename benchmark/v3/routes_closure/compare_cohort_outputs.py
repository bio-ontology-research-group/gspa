#!/usr/bin/env python3
"""Compare the outputs of several cohort routes genome by genome, file by file (plan 2026-11-01, check A6).

Routes are given as ``name=<root>=<pattern>`` where the pattern names the integrated table of a genome relative to
the root with ``{g}`` for the genome id (for example ``direct=/job/direct={g}/integrated.tsv``,
``nextflow=/job/published/cohort={g}/{g}_integrated.tsv``, ``web=/job/web/jobs/<id>/gspa_out/cohort={g}/integrated.tsv``).
For each genome every artifact next to the table is compared: the table bytes; the stage ledger (names, statuses,
reasons, notes and change lists; ``elapsed_seconds`` and the route-specific feature-export ``directory`` excluded);
the learned report (whole document); the cohort report; provenance; suggestions; enforcement actions; and the
feature-export files (``<features-root>/{g}/*``, given as ``name=<dir>`` with ``--features``). Absent-on-all-routes
is recorded, not counted as agreement. Nothing is summarised away; the result is one JSON record.
"""
import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

SUFFIXES = {'ledger': '.stages.json', 'learned': '.learned.json', 'cohort': '.cohort.json', 'provenance': '.provenance.json',
            'suggestions': '.suggestions.tsv', 'actions': '.enforcement.tsv'}


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b''):
            h.update(chunk)
    return h.hexdigest()


def ledger_view(path):
    d = json.loads(Path(path).read_text())
    stages = []
    for s in d['stages']:
        s = {k: v for k, v in s.items() if k not in ('elapsed_seconds', 'directory')}   # directory: the route's own feature-export path
        for block in (s.get('notes'), s.get('summary')):
            if isinstance(block, dict):
                block.pop('directory', None)
        stages.append(s)
    return {'top': {k: v for k, v in d.items() if k != 'stages'}, 'stages': stages}


def learned_facts(path):
    rep = json.loads(Path(path).read_text())
    reasons = Counter(c.get('unsupported_reason') for c in rep['candidates'] if c.get('unsupported_reason'))
    scored = [c for c in rep['candidates'] if c.get('probability') is not None]
    return {'model_sha256': rep['model']['sha256'], 'counts': rep['counts'], 'scored_candidates': len(scored),
            'unsupported_reasons': dict(reasons)}


def stage_facts(path):
    d = json.loads(Path(path).read_text())
    return {'scoring_mode': d.get('scoring_mode'), 'stages': {s['name']: s.get('status') for s in d['stages']},
            'cohort_reason': next((s.get('reason') for s in d['stages'] if s['name'] == 'cohort_context'), None),
            'taxon_state': next((s.get('notes', {}).get('taxon_consistency_state') or s.get('taxon_consistency_state') for s in d['stages'] if s['name'] == 'reference_preparation'), None)}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--genomes', nargs='+', required=True)
    ap.add_argument('--route', action='append', required=True, help='name=<root>=<pattern>')
    ap.add_argument('--features', action='append', default=[], help='name=<features root containing {g}/>')
    ap.add_argument('--out', required=True)
    args = ap.parse_args()
    routes = {}
    for spec in args.route:
        name, root, pattern = spec.split('=', 2)
        routes[name] = (Path(root), pattern)
    features = dict(f.split('=', 1) for f in args.features)
    names = list(routes)
    result = {'routes': {n: str(r[0]) for n, r in routes.items()}, 'genomes': {}, 'identical': True}
    for g in args.genomes:
        rec = {'artifacts': {}, 'identical': True}
        tables = {n: root / pattern.format(g=g) for n, (root, pattern) in routes.items()}
        for art, suffix in [('table', '')] + list(SUFFIXES.items()):
            paths = {n: Path(str(t) + suffix) for n, t in tables.items()}
            present = {n: p.exists() for n, p in paths.items()}
            entry = {'present': present}
            if all(present.values()):
                if art == 'ledger':
                    views = {n: ledger_view(p) for n, p in paths.items()}
                    entry['identical'] = all(views[n] == views[names[0]] for n in names)
                    entry['bytes_identical'] = len({sha(p) for p in paths.values()}) == 1
                    if not entry['identical']:
                        v0 = views[names[0]]
                        entry['differing_stages'] = {n: [a.get('name') for a, b in zip(views[n]['stages'], v0['stages']) if a != b] for n in names[1:]}
                        entry['top_identical'] = all(views[n]['top'] == v0['top'] for n in names)
                    entry['facts'] = stage_facts(paths[names[0]])
                elif art == 'learned':
                    docs = {n: json.loads(p.read_text()) for n, p in paths.items()}
                    entry['identical'] = all(docs[n] == docs[names[0]] for n in names)
                    entry['bytes_identical'] = len({sha(p) for p in paths.values()}) == 1
                    if not entry['identical']:
                        entry['differing_keys'] = {n: sorted(k for k in set(docs[n]) | set(docs[names[0]]) if docs[n].get(k) != docs[names[0]].get(k)) for n in names[1:]}
                    entry['facts'] = learned_facts(paths[names[0]])
                else:
                    hashes = {n: sha(p) for n, p in paths.items()}
                    entry['identical'] = len(set(hashes.values())) == 1
                    entry['sha256'] = hashes
                    if art == 'table':
                        with paths[names[0]].open() as fh:
                            entry['rows'] = sum(1 for _ in fh) - 1
            elif any(present.values()):
                entry['identical'] = False
                entry['note'] = 'written by some routes only'
            else:
                entry['identical'] = None
                entry['note'] = 'written by no route'
            if entry['identical'] is False:
                rec['identical'] = False
            rec['artifacts'][art] = entry
        if features:
            sets = {}
            for n, root in features.items():
                d = Path(root) / g
                sets[n] = {p.relative_to(d).as_posix(): sha(p) for p in d.rglob('*') if p.is_file()} if d.is_dir() else None
            if all(s is not None for s in sets.values()):
                same = all(sets[n] == sets[list(sets)[0]] for n in sets)
                rec['artifacts']['features'] = {'identical': same, 'files': len(sets[list(sets)[0]]),
                                                'differing': sorted(k for k in set().union(*sets.values()) if len({s.get(k) for s in sets.values()}) > 1)}
                if not same:
                    rec['identical'] = False
            else:
                rec['artifacts']['features'] = {'identical': False if any(s is not None for s in sets.values()) else None, 'present': {n: s is not None for n, s in sets.items()}}
                if any(s is not None for s in sets.values()):
                    rec['identical'] = False
        result['genomes'][g] = rec
        if not rec['identical']:
            result['identical'] = False
    Path(args.out).write_text(json.dumps(result, indent=2) + '\n')
    for g, rec in result['genomes'].items():
        print(g, 'identical' if rec['identical'] else 'DIFFERENT', {a: e.get('identical') for a, e in rec['artifacts'].items()})
    print('ALL IDENTICAL' if result['identical'] else 'DIFFERENCES FOUND')
    return 0 if result['identical'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
