#!/usr/bin/env python3
"""Draw a compact overview from executed results, without inferring accuracy."""
import html
import json
from pathlib import Path

root=Path('.tutorial-work')
cases=[('bacteria','Bacterial genome'),('archaea','Archaeal genome'),('eukaryote','Yeast proteome'),
       ('mag','Constructed archaeal bin'),('metagenome','Constructed mixed community'),
       ('virus','Lambda phage'),('prokaryote','Prokaryote, domain unspecified')]
rows=[(label,json.loads((root/case/'summary.json').read_text())) for case,label in cases]
svg=['<svg xmlns="http://www.w3.org/2000/svg" width="1000" height="590" viewBox="0 0 1000 590">',
     '<rect width="1000" height="590" rx="18" fill="#101f33"/>',
     '<g font-family="sans-serif" fill="#e8f0f8">',
     '<text x="36" y="42" font-size="13" fill="#76d9cd">GSPA / SEVEN EXECUTED BIOLOGICAL WORKFLOWS</text>',
     '<text x="36" y="85" font-size="30">From sequences to inspectable annotations</text>',
     '<text x="36" y="119" font-size="15" fill="#aebed0">Bars show the fraction of proteins with GO/EC search evidence in the teaching reference.</text>']
for i,(label,stat) in enumerate(rows):
    y=166+i*49
    fraction=stat['with_go_ec']/stat['proteins']
    svg += [f'<text x="36" y="{y+15}" font-size="16">{html.escape(label)}</text>',
            f'<rect x="340" y="{y}" width="370" height="20" rx="5" fill="#34465d"/>',
            f'<rect x="340" y="{y}" width="{370*fraction:.2f}" height="20" rx="5" fill="#57c7bb"/>',
            f'<text x="740" y="{y+15}" font-size="16">{stat["with_go_ec"]:,} / {stat["proteins"]:,} proteins</text>']
svg += ['<text x="36" y="551" font-size="14" fill="#aebed0">Real sequences and searches; close/self matches included. These are execution results, not accuracy estimates.</text>',
        '</g></svg>']
(root/'overview.svg').write_text('\n'.join(svg)+'\n')
print('Seven biological workflows have nonempty, inspectable results.')
print('Overview chart generated from their recorded protein counts.')
