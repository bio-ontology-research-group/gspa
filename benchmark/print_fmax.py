#!/usr/bin/env python3
import json
import sys

for path in sys.argv[1:]:
    try:
        d = json.load(open(path))
    except Exception as e:
        print(f'{path}: ERROR {e}')
        continue
    tag = d.get('tag', path)
    print(f'=== {tag} ===')
    for t, b in d.get('by_truth', {}).items():
        n = b.get('truth_proteins', 0)
        print(f'  truth={t} n={n}')
        for r in b.get('results', []):
            method = r['method']
            ci = r.get('fmax_ci', [0, 0])
            print(f'    {method:15s} Fmax={r["fmax_overall"]:.3f} [{ci[0]:.3f},{ci[1]:.3f}]  '
                  f'MF={r["fmax_MF"]:.3f}  BP={r["fmax_BP"]:.3f}  CC={r["fmax_CC"]:.3f}  '
                  f'icR={r["ic_recall"]:.3f}  cov={r["coverage"]:.3f}')
