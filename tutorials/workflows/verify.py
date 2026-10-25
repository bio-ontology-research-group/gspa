#!/usr/bin/env python3
"""Rerun every CPU lesson and check the evidence used in its worked interpretation."""
import csv
import json
from pathlib import Path
import subprocess

from run import CASES, ROOT, WORK


def rows(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def main():
    observed = {}
    expected = {
        'bacteria': (4319, 3876), 'archaea': (2102, 675), 'eukaryote': (6733, 6051),
        'mag': (1307, 384), 'metagenome': (6749, 4594), 'virus': (62, 45),
        'prokaryote': (4319, 3876),
    }
    for case in CASES:
        subprocess.run([str(ROOT / 'tutorial'), 'run', case], cwd=ROOT, check=True)
        stat = json.loads((WORK / case / 'summary.json').read_text())
        assert (stat['proteins'], stat['with_go_ec']) == expected[case], (case, stat)
        observed[case] = rows(WORK / case / 'gspa_out/integrated.tsv')
        assert (WORK / case / 'gspa_out/provenance.json').stat().st_size
        browser = WORK / case / 'gspa_out' / f'{case}_browser.html'
        assert '<html' in browser.read_text().lower(), browser
    def candidates(case, term):
        return {r['protein_id'] for r in observed[case] if r['function_id'] == term}
    threonine = candidates('bacteria', 'GO:0009088')
    assert threonine == {'bacteria_00002', 'bacteria_00003', 'bacteria_00004',
                         'bacteria_03361', 'bacteria_03861'}
    assert candidates('prokaryote', 'GO:0009088') == threonine
    assert not candidates('archaea', 'GO:0005634')
    actions = rows(WORK / 'archaea/gspa_out/actions.tsv')
    assert len({r['protein_id'] for r in actions if r['term'] == 'GO:0005634' and r['action'] == 'remove'}) == 17
    assert len(candidates('eukaryote', 'GO:0005634')) == 1939
    assert len(candidates('eukaryote', 'GO:0005739')) == 1153
    assert candidates('mag', 'GO:0009088') == {'mag_00678', 'mag_00680'}
    assert candidates('virus', 'GO:0019073') == {'virus_00001', 'virus_00002', 'virus_00006'}
    assert candidates('metagenome', 'GO:0019073') == {'metagenome_06685', 'metagenome_06686', 'metagenome_06690'}
    # Verify the central neighbourhood and community attribution claims against CDS coordinates.
    import results
    results.ROOT = WORK
    coords = dict(results.features('bacteria'))
    assert [(coords[p][3], coords[p][4]) for p in sorted(threonine)[:3]] == [
        ('337', '2799'), ('2801', '3733'), ('3734', '5020')]
    mixed = dict(results.features('metagenome'))
    from collections import Counter
    assert Counter(mixed[p][0].rsplit('_', 1)[0] for p in candidates('metagenome', 'GO:0009088')) == {'bacteria': 5, 'archaea': 2}
    assert all(mixed[p][0].startswith('virus_') for p in candidates('metagenome', 'GO:0019073'))
    assert sum(r['type'] == 'GO' for r in observed['bacteria']) == 21226
    assert sum(r['type'] == 'GO' for r in observed['prokaryote']) == 21211
    print('PASS: seven CPU lessons, browsers, and the biological observations used in the text.')


if __name__ == '__main__':
    main()
