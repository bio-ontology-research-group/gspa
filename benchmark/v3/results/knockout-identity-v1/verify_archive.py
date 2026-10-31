#!/usr/bin/env python3
"""Offline byte, source, table and numerical checks; no model or solver execution."""
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys
import tarfile
import xml.etree.ElementTree as ET

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent
expected = {}
for line in (ROOT / 'SHA256SUMS').read_text().splitlines():
    digest, name = line.split('  ', 1)
    assert name not in expected
    expected[name] = digest
actual = {str(p.relative_to(ROOT)) for p in ROOT.rglob('*') if p.is_file() and p.name != 'SHA256SUMS'}
assert actual == set(expected)
for name, digest in expected.items():
    assert hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == digest, name

def read(name):
    return json.loads((ROOT / name).read_text())

job = ROOT / 'remote/jobs/10218'
assert (job / 'exit-code.txt').read_text().strip() == '0'
for name, count in [('tests.log', 144), ('benchmark-tests.log', 118)]:
    text = (job / name).read_text()
    assert f'Ran {count} tests in' in text and text.rstrip().endswith('OK')
before = (job / 'before-tests.log').read_text()
assert 'Ran 14 tests in' in before and 'FAILED (failures=15, errors=5)' in before
package = ROOT / 'source/source-v4.tar.gz'
assert hashlib.sha256(package.read_bytes()).hexdigest() in (job / 'packages.sha256').read_text()
with tarfile.open(package) as archive:
    for name in ['scripts/metabolic_experiments.py', 'scripts/community_metabolism.py',
                 'scripts/tests/test_knockout_identity.py', 'scripts/tests/test_community_knockout_identity.py',
                 'benchmark/v3/prepare_germicidin_observations.py', 'benchmark/v3/run_germicidin_diagnostic.py']:
        assert archive.extractfile('home/leechuck/Public/software/gspa/' + name).read() == (ROOT / 'source/final' / name).read_bytes()
    model_bytes = archive.extractfile('home/leechuck/tmp/germicidin-input/Sco-GEM.xml').read()
    medium_bytes = archive.extractfile('home/leechuck/tmp/germicidin-input/minimal_medium.csv').read()
base = 'remote/jobs/10218/diagnostic/'
protocol = read(base + 'protocol.json')
assert hashlib.sha256(model_bytes).hexdigest() == protocol['model_sha256']
assert hashlib.sha256(medium_bytes).hexdigest() == protocol['medium_sha256']
for key, name in [('before_source_sha256', 'before/metabolic_experiments.py'),
                  ('after_source_sha256', 'source/final/scripts/metabolic_experiments.py'),
                  ('diagnostic_source_sha256', 'source/final/benchmark/v3/run_germicidin_diagnostic.py')]:
    assert hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == protocol[key]
rows = read(base + 'results.json')
keys = {(r['condition'], r['genotype'], r['requirement']) for r in rows}
assert len(rows) == len(keys) == 18
assert keys == {(c, g, q) for c in protocol['conditions'] for g in protocol['genotypes'] for q in protocol['requirements']}

def agree(a, b):
    if isinstance(a, dict):
        return isinstance(b, dict) and set(a) == set(b) and all(agree(a[k], b[k]) for k in a)
    if isinstance(a, (int, float)):
        return isinstance(b, (int, float)) and math.isfinite(a) and math.isfinite(b) and math.isclose(a, b, rel_tol=1e-6, abs_tol=1e-7)
    return a == b

by_key = {tuple(r[k] for k in ('condition', 'genotype', 'requirement')): r for r in rows}
for row in rows:
    assert agree(row['production'], row['oracle'])
    assert row['production'] == row['before'] and row['knocked_out_reactions'] == row['before_knocked_out_reactions']
    record = row['production']
    if row['genotype'] == 'gcs_deleted':
        assert row['knocked_out_reactions'] == ['R_GERMA', 'R_GERMB', 'R_GERMC', 'R_GERMD']
    if row['condition'] == 'carbon_closed':
        assert record == {'growth_status': 'infeasible', 'growth': None, 'capacities': {}}
    elif row['condition'] == 'phosphate_closed' and row['requirement'] == 'growth_fraction':
        assert record['growth_status'] == 'optimal' and abs(record['growth']) < 1e-6 and not record['capacities']
    else:
        assert record['growth_status'] == 'optimal' and set(record['capacities']) == set(protocol['product_targets'])
        for v in record['capacities'].values():
            assert v['statuses'] == ['optimal', 'optimal'] and abs(v['minimum']) < 1e-7
            assert (abs(v['maximum']) < 1e-7) if row['genotype'] == 'gcs_deleted' else (v['maximum'] > 1e-7)
    if row['genotype'] == 'native_repeat':
        assert record == by_key[row['condition'], 'native', row['requirement']]['production']
assert sum(len(r['solver_calls']) for r in rows) == 90

# Check the label-only fixture independently from the archived SBML bytes.
ns = {'s': 'http://www.sbml.org/sbml/level3/version1/core', 'f': 'http://www.sbml.org/sbml/level3/version1/fbc/version2'}
f = '{' + ns['f'] + '}'
def structure(data):
    root = ET.fromstring(data)
    params = {p.attrib['id']: float(p.attrib['value']) for p in root.findall('.//s:parameter', ns)}
    genes = {g.attrib[f + 'id']: g.attrib.get(f + 'label') for g in root.findall('.//f:geneProduct', ns)}
    mets = {m.attrib['id']: (m.attrib['compartment'], m.attrib.get(f + 'chemicalFormula')) for m in root.findall('.//s:species', ns)}
    def gpr(node):
        if node is None: return None
        return (node.tag, tuple(sorted(node.attrib.items())), tuple(gpr(c) for c in node))
    reactions = {}
    for r in root.findall('.//s:reaction', ns):
        sides = [tuple((s.attrib['species'], float(s.attrib.get('stoichiometry', 1))) for s in r.findall(path, ns))
                 for path in ('s:listOfReactants/s:speciesReference', 's:listOfProducts/s:speciesReference')]
        reactions[r.attrib['id']] = (sides, params[r.attrib[f + 'lowerFluxBound']], params[r.attrib[f + 'upperFluxBound']], gpr(r.find('f:geneProductAssociation', ns)))
    return genes, mets, reactions
label_bytes = (ROOT / base / 'Sco-GEM-labelled.xml').read_bytes()
old, new = structure(model_bytes), structure(label_bytes)
assert old[1:] == new[1:]
assert new[0] == {**old[0], 'G_SCO7221': 'SCO7221'}
aliases = read(base + 'label-regressions.json')
assert len(aliases) == 2
for row in aliases:
    assert hashlib.sha256(label_bytes).hexdigest() == row['derived_label_only_model_sha256']
    assert row['before_knocked_out_reactions'] == []
    assert row['after'] == by_key['published_model_minimal', 'gcs_deleted', row['requirement']]['production']
    assert all(v['maximum'] > 1e-7 for v in row['before']['capacities'].values())
assert sum(len(r['solver_calls']) + len(r['before_solver_calls']) for r in rows) + sum(
    len(r['before_solver_calls']) + len(r['after_solver_calls']) for r in aliases) == 216

spec = importlib.util.spec_from_file_location('primary_observations', ROOT / 'source/final/benchmark/v3/prepare_germicidin_observations.py')
module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
observations = read('observations.json')
assert module.prepare(ROOT / 'primary/article.xml') == observations
assert len(observations['records']) == 24
assert [sum(r['status'] == status for r in observations['records']) for status in
        ('detected', 'not_detected', 'not_determined')] == [11, 5, 8]
assert all(r['value'] is None for r in observations['records'] if r['status'] != 'detected')
print(f'Verified {len(expected)} files; 262 tests; 18 native pairs and oracle scenarios; two label regressions; all 24 primary cells.')
