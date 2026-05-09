#!/usr/bin/env python3
"""Convert sidecar TSVs (protein_id, term, score, annotation_type) to JSONL
claims compatible with gspa integrate.

Multiple input files can be passed; each gets a --source label assigned
to its claims.
"""
import argparse, json, sys

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--input', action='append', nargs=2, metavar=('SOURCE', 'TSV'),
                    required=True, help='Pair: --input <source_label> <tsv> ; repeatable')
    ap.add_argument('--output', required=True)
    ap.add_argument('--go-aspect-map', help='go.obo path; if given, look up aspect for each GO term')
    a = ap.parse_args()

    aspect = {}
    if a.go_aspect_map:
        cur = None; ns = None
        with open(a.go_aspect_map) as f:
            for line in f:
                line = line.rstrip()
                if line.startswith('[Term]'): cur = ns = None
                elif line.startswith('id: GO:'): cur = line.split(' ', 1)[1].strip()
                elif line.startswith('namespace: '):
                    n = line.split(' ', 1)[1].strip()
                    if cur:
                        if n == 'molecular_function': aspect[cur] = 'MF'
                        elif n == 'biological_process': aspect[cur] = 'BP'
                        elif n == 'cellular_component': aspect[cur] = 'CC'

    n_total = 0
    seen = set()
    with open(a.output, 'w') as out_fh:
        for source, tsv in a.input:
            try:
                fh = open(tsv)
            except FileNotFoundError:
                print(f'WARN: skipping {tsv} (missing)', file=sys.stderr)
                continue
            n = 0
            with fh:
                header = fh.readline().rstrip('\n').split('\t')
                # Find columns
                try:
                    pi = header.index('protein_id')
                except ValueError:
                    pi = 0
                try:
                    ti = header.index('term')
                except ValueError:
                    ti = 1
                try:
                    si = header.index('score')
                except ValueError:
                    si = 2
                try:
                    ai = header.index('annotation_type')
                except ValueError:
                    ai = 3
                for line in fh:
                    fs = line.rstrip('\n').split('\t')
                    if len(fs) <= max(pi, ti, si, ai):
                        continue
                    pid, term, score_s, ann_type = fs[pi], fs[ti], fs[si], fs[ai]
                    if not term:
                        continue
                    try:
                        score = float(score_s)
                    except ValueError:
                        continue
                    # Normalize annotation type
                    if ann_type in ('GO',) or term.startswith('GO:'):
                        ftype = 'GO'
                        fid = term if term.startswith('GO:') else 'GO:' + term
                    elif ann_type in ('EC',) or term.startswith('EC:'):
                        ftype = 'EC'
                        fid = term if term.startswith('EC:') else 'EC:' + term.replace('EC:','')
                    else:
                        ftype = ann_type or 'GO'
                        fid = term
                    asp = aspect.get(fid, '') if ftype == 'GO' else ''
                    key = (pid, fid, source)
                    if key in seen: continue
                    seen.add(key)
                    out_fh.write(json.dumps({
                        'protein_id': pid,
                        'function_type': ftype,
                        'function_id': fid,
                        'go_aspect': asp,
                        'source': source,
                        'raw_score': round(score, 6),
                        'metadata': {},
                    }) + '\n')
                    n += 1
            print(f'  {source}: {n:,} claims from {tsv}', file=sys.stderr)
            n_total += n
    print(f'wrote {n_total:,} claims → {a.output}', file=sys.stderr)

if __name__ == '__main__':
    main()
