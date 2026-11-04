"""Supervisor check: every number quoted in the report's tier tables comes from the shipped tables."""
import csv, json

base = '/home/leechuck/tmp/gspa-v3-offgrowth2'
checks = [
    ('A', 'N1', 'precision', 'test', '24', '29'),
    ('A', 'N1', 'gene_recovery', 'masked_offgrowth', '97', '800'),
    ('A', 'S1', 'precision', 'test', '71', '106'),
    ('A', 'S0', 'precision', 'test', '6', '6'),
    ('A', 'Q', 'precision', 'test', '26', '29'),
    ('B', 'N1', 'gene_recovery', 'masked_offgrowth', '3', '60'),
    ('C', 'N1', 'precision', 'test', '16', '21'),
    ('C', 'N1', 'gene_recovery', 'masked_offgrowth', '42', '280'),
    ('C', 'S1', 'precision', 'test', '58', '72'),
]
ok = True
for tier, arm, metric, label, num, den in checks:
    rows = list(csv.DictReader(open(f'{base}/score-tier{tier}-clean/out/arms.tsv'), delimiter='\t'))
    hit = [r for r in rows if r['group'] == 'pooled' and r['arm'] == arm and r['metric'] == metric and r['label_set'] == label]
    good = hit and hit[0]['numerator'] == num and hit[0]['denominator'] == den
    ok &= bool(good)
    print(('OK  ' if good else 'BAD '), tier, arm, metric, (hit[0]['numerator'], hit[0]['denominator']) if hit else None, 'expected', (num, den))

# pooled B+C recomputed from decisions and cases
sup = harm = gene = rxn = n = 0
for tier in ('B', 'C'):
    for d in csv.DictReader(open(f'{base}/score-tier{tier}-clean/out/decisions.tsv'), delimiter='\t'):
        if d['arm'] == 'N1' and d['label_set'] == 'test':
            sup += d['classification'] == 'supported'
            harm += d['classification'] == 'harmful'
    for c in csv.DictReader(open(f'{base}/score-tier{tier}-clean/out/cases.tsv'), delimiter='\t'):
        if c['arm'] == 'N1' and c['stratum'] == 'off_growth':
            n += 1
            gene += int(c['gene_recovered'])
            rxn += int(c['reaction_recovered'])
print(f'pooled B+C: supported {sup} harmful {harm} precision {sup/(sup+harm):.3f} | gene {gene}/{n} reaction {rxn}/{n}')
print('expected in report: 17 / 5 / 0.773 | 45/340 and 89/340')
s = json.load(open(f'{base}/score-tierA-clean/out/summary.json'))
print('tier A own-CLEAN additional functions:', len(s['additional_functions']),
      sum(1 for x in s['additional_functions'] if x['classification'] == 'supported'), 'supported (report says 17 / 13)')
print('ALL ARM CHECKS OK' if ok else 'MISMATCH FOUND')
