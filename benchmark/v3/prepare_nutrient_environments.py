#!/usr/bin/env python3
"""Audit raw 2026 nutrient-environment measurements without executing author code.

Genotypes and assay conditions are separate from measured outcomes. No model,
simulation result, filtered author table, or predicted auxotrophy selects cases.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path
import statistics
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

import numpy as np
import openpyxl

COMMIT = 'ee3a2079a4a7e47c090a7ad9945b053c72723490'
REPOSITORY = 'amichaibk/nutrient-dependent-auxotrophy_and_CF'
DOI = '10.64898/2026.05.20.726445'
PAPER = 'https://www.biorxiv.org/content/10.64898/2026.05.20.726445v3'
PAPER_XML = 'https://www.biorxiv.org/content/early/2026/08/06/2026.05.20.726445.source.xml'
AUX = {'M': 'methionine', 'F': 'phenylalanine', 'R': 'arginine',
       'C': 'cysteine', 'S': 'serine', 'L': 'leucine'}
SUPPLEMENT = {'M': 'Met', 'F': 'Phe', 'R': 'Arg', 'C': 'Cys', 'S': 'Ser', 'L': 'Leu'}
GENES = {'Ec': dict(zip(AUX, ['metA', 'pheA', 'argA', 'cysE', 'serA', 'leuB'])),
         'Bs': dict(zip(AUX, ['metA', 'pheA', 'argH', 'cysE', 'serA', 'leuB']))}
CARBON = {'Gluc': 'D-glucose', 'Gly': 'glycerol', 'Fum': 'sodium fumarate', 'NoC': None}
NITROGEN = {'Glut': 'L-glutamic acid', 'AmSu': 'ammonium sulfate', 'Urea': 'urea',
            'Asp': 'L-aspartic acid', 'Asn': 'L-asparagine', 'Gln': 'L-glutamine',
            'GNAc': 'N-acetylglucosamine', 'Glcn': 'D-glucosamine hydrochloride',
            'Ala': 'L-alanine', 'Pro': 'L-proline', 'NoN': None}
BASE_MEDIUM = {'potassium phosphate': '5 mM', 'MOPS': '100 mM, pH 7.05', 'MgCl2': '2 mM',
               'MnCl2': '50 uM', 'ZnCl2': '1 uM', 'thiamine': '2 uM', 'L-tryptophan': '50 ug/ml',
               'CaCl2': '700 uM', 'FeCl3': '50 uM'}
FILES = ['README.md', 'Consolidated Analysis and Figures.ipynb', 'Data/Sample Names.xlsx',
         'Data/AA_synthesis_pathways.csv', 'Generated Tables/big_norm_df.csv',
         'Generated Tables/Supplemented_AA_data.csv', 'Generated Tables/filtered_df.csv'] + [
             f'Data/Crossfeeding Experiments/{species}_data_Ex{rep}.xlsx'
             for species in GENES for rep in (1, 2, 3)]


class DataError(ValueError):
    pass


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def fetch(url):
    with urllib.request.urlopen(urllib.request.Request(url, headers={'User-Agent': 'gspa-validation'}), timeout=45) as response:
        return response.read()


def acquire(root):
    root.mkdir(parents=True, exist_ok=True)
    commit_data = fetch(f'https://api.github.com/repos/{REPOSITORY}/git/commits/{COMMIT}')
    commit = json.loads(commit_data)
    if commit.get('sha') != COMMIT:
        raise DataError('unexpected commit metadata')
    url = f'https://api.github.com/repos/{REPOSITORY}/git/trees/{commit["tree"]["sha"]}?recursive=1'
    payload = fetch(url)
    tree = json.loads(payload)
    if tree.get('truncated') or tree.get('sha') != commit['tree']['sha']:
        raise DataError('incomplete or changed pinned source tree')
    entries = {e['path']: e for e in tree['tree'] if e['type'] == 'blob'}
    receipt = {'schema': 'gspa.nutrient_acquisition/1', 'commit': COMMIT, 'repository': REPOSITORY,
               'tree_url': url, 'tree_sha256': hashlib.sha256(payload).hexdigest(),
               'commit_metadata_sha256': hashlib.sha256(commit_data).hexdigest(), 'files': []}
    (root / 'git-commit.json').write_bytes(commit_data)
    (root / 'tree.json').write_bytes(payload)
    for name in FILES:
        entry = entries[name]
        url = f'https://raw.githubusercontent.com/{REPOSITORY}/{COMMIT}/' + urllib.parse.quote(name)
        path = root / 'raw' / name
        data = path.read_bytes() if path.is_file() else fetch(url)
        blob = hashlib.sha1(b'blob ' + str(len(data)).encode() + b'\0' + data).hexdigest()
        if blob != entry['sha'] or len(data) != entry['size']:
            raise DataError('pinned Git blob mismatch: ' + name)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
        receipt['files'].append({'path': name, 'url': url, 'git_blob': blob,
                                 'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()})
    paper_path = root / 'paper-v3.xml'
    paper_data = paper_path.read_bytes() if paper_path.exists() else fetch(PAPER_XML)
    paper = ET.fromstring(paper_data)
    identities = [''.join(e.itertext()) for e in paper.findall('.//article-id')]
    if DOI not in identities or 'biorxiv;2026.05.20.726445v3' not in identities:
        raise DataError('paper identity mismatch')
    paper_path.write_bytes(paper_data)
    receipt['paper'] = {'url': PAPER_XML, 'version': 3, 'bytes': len(paper_data),
                        'sha256': hashlib.sha256(paper_data).hexdigest(), 'path': 'paper-v3.xml'}
    (root / 'acquisition.json').write_text(json.dumps(receipt, indent=2) + '\n')
    return receipt


def verify_inputs(root):
    receipt = json.loads((root / 'acquisition.json').read_text())
    if receipt['commit'] != COMMIT:
        raise DataError('unexpected source version')
    for filename, key in [('tree.json', 'tree_sha256'), ('git-commit.json', 'commit_metadata_sha256')]:
        if sha(root / filename) != receipt[key]:
            raise DataError('acquisition metadata changed: ' + filename)
    tree = json.loads((root / 'tree.json').read_text())
    commit = json.loads((root / 'git-commit.json').read_text())
    if commit['sha'] != COMMIT or tree['sha'] != commit['tree']['sha'] or tree.get('truncated'):
        raise DataError('acquisition metadata identity mismatch')
    entries = {e['path']: e for e in tree['tree'] if e['type'] == 'blob'}
    by_path = {r['path']: r for r in receipt['files']}
    if len(by_path) != len(receipt['files']) or set(by_path) != set(FILES):
        raise DataError('missing, duplicated or extra source files')
    for name, record in by_path.items():
        p = root / 'raw' / name
        if p.stat().st_size != record['bytes'] or sha(p) != record['sha256']:
            raise DataError('source content changed: ' + name)
        data = p.read_bytes()
        blob = hashlib.sha1(b'blob ' + str(len(data)).encode() + b'\0' + data).hexdigest()
        if blob != entries[name]['sha'] or blob != record['git_blob']:
            raise DataError('source differs from pinned Git tree: ' + name)
    if receipt['paper']['version'] != 3 or sha(root / receipt['paper']['path']) != receipt['paper']['sha256']:
        raise DataError('paper provenance changed')
    return receipt


def plate_grid(sheet):
    """Locate the actual labelled 8x12 plate, not a fixed export header length."""
    rows = list(sheet.iter_rows(values_only=True))
    matches = []
    for i, row in enumerate(rows):
        for j, value in enumerate(row):
            if value != 'A' or i == 0 or i + 7 >= len(rows):
                continue
            if list(rows[i - 1][j + 1:j + 13]) != list(range(1, 13)):
                continue
            if any(len(rows[i + k]) <= j + 12 or rows[i + k][j] != letter for k, letter in enumerate('ABCDEFGH')):
                continue
            matches.append({f'{letter}{col}': (rows[i + k][j + col], i + k + 1, j + col + 1)
                            for k, letter in enumerate('ABCDEFGH') for col in range(1, 13)})
    if len(matches) != 1:
        raise DataError(f'{sheet.title}: expected exactly one labelled plate, found {len(matches)}')
    return matches[0]


def sample_map(path):
    workbook = openpyxl.load_workbook(path, data_only=True, read_only=True)
    try:
        if len(workbook.sheetnames) != 1:
            raise DataError('ambiguous sample map')
        cells = plate_grid(workbook.active)
    finally:
        workbook.close()
    result = {}
    for well, (value, _, _) in cells.items():
        if not isinstance(value, str) or len(value.split('+')) != 2:
            raise DataError('invalid sample identity: ' + well)
        nitrogen, carbon = [x.strip() for x in value.split('+')]
        if nitrogen not in NITROGEN or carbon not in CARBON:
            raise DataError('unknown medium ingredient')
        result[well] = {'nitrogen': nitrogen, 'carbon': carbon, 'environment': nitrogen + ' + ' + carbon}
    expected = Counter({'Glut': 8, 'AmSu': 8, 'Urea': 8, 'Asp': 8, 'Asn': 8, 'Gln': 8,
                        'GNAc': 8, 'Glcn': 8, 'Ala': 8, 'Pro': 8, 'NoN': 16})
    if Counter(r['nitrogen'] for r in result.values()) != expected:
        raise DataError('sample-map condition coverage differs from the registered 40 environments plus controls')
    if Counter((r['nitrogen'], r['carbon']) for r in result.values()) != Counter({
            (n, c): 4 if n == 'NoN' else 2 for n in NITROGEN for c in CARBON}):
        raise DataError('sample-map environment replication differs from the protocol')
    if any(result[w]['nitrogen'] != 'NoN' or result[w]['carbon'] != 'NoC' for w in ('G11', 'G12', 'H11', 'H12')):
        raise DataError('registered blank wells are not no-carbon/no-nitrogen controls')
    return result


def pair_kind(pair):
    if pair == 'Blank':
        return 'blank', []
    if pair == 'WT':
        return 'wild_type', []
    parts = pair.split('+')
    if len(parts) != 2 or parts[0] not in AUX:
        raise DataError('unknown strain pair: ' + pair)
    first, second = parts
    if second == first:
        return 'monoculture', [first]
    if second == SUPPLEMENT[first]:
        return 'supplemented_monoculture', [first]
    if second in AUX:
        return 'coculture', parts
    raise DataError('unknown supplementation: ' + pair)


def number(value):
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, (float, int)) or not math.isfinite(value):
        raise DataError('plate measurement must be finite numeric or explicitly missing')
    return float(value)


def expected_sheets():
    expected = {'WT', 'Blank'} | {a + '+' + a for a in AUX} | {a + '+' + SUPPLEMENT[a] for a in AUX}
    expected |= {a + '+' + b for i, a in enumerate(AUX) for b in list(AUX)[i + 1:]}
    return expected


def read_measurements(root, mapping):
    records, coverage = [], []
    expected = expected_sheets()
    for species in GENES:
        for experiment in (1, 2, 3):
            relative = f'Data/Crossfeeding Experiments/{species}_data_Ex{experiment}.xlsx'
            path = root / 'raw' / relative
            workbook = openpyxl.load_workbook(path, data_only=True, read_only=True)
            try:
                present = set(workbook.sheetnames)
                # Author sheets use S+L before R pairs; pair order itself is still explicit.
                if present - expected:
                    raise DataError('unexpected assay sheet: ' + str(sorted(present - expected)))
                coverage.append({'species': species, 'experiment': experiment,
                                 'sheets_present': sorted(present), 'sheets_missing': sorted(expected - present)})
                for pair in workbook.sheetnames:
                    kind, strains = pair_kind(pair)
                    cells = plate_grid(workbook[pair])
                    values = {well: number(cell[0]) for well, cell in cells.items()}
                    blanks = [values[w] for w in ('G11', 'G12', 'H11', 'H12')]
                    blank = statistics.mean(blanks) if all(v is not None for v in blanks) else None
                    for well, value in values.items():
                        corrected = float(np.round(value - blank, 3)) if value is not None and blank is not None else None
                        records.append({'species': species, 'experiment': experiment, 'pair': 'WT+WT' if pair == 'WT' else pair,
                                        'kind': kind, 'strains': strains, **mapping[well], 'well': well,
                                        'raw_od600': value, 'plate_blank_od600': blank, 'corrected_od600_unclipped': corrected,
                                        'corrected_od600': max(.005, corrected) if corrected is not None else None,
                                        'status': 'measured' if corrected is not None else 'unavailable',
                                        'source_file': relative, 'source_sheet': pair,
                                        'source_cell': openpyxl.utils.get_column_letter(cells[well][2]) + str(cells[well][1])})
            finally:
                workbook.close()
    return records, coverage


def sem(values):
    return statistics.stdev(values) / math.sqrt(len(values)) if len(values) > 1 else None


def summarize(records):
    groups = defaultdict(list)
    for row in records:
        if row['kind'] != 'blank':
            groups[row['species'], row['environment'], row['pair']].append(row)
    result = []
    # Full assay catalogue is independent of observed growth and missing sheets.
    catalogue = {}
    for species in GENES:
        for nitrogen in NITROGEN:
            for carbon in CARBON:
                environment = nitrogen + ' + ' + carbon
                for sheet in sorted(expected_sheets() - {'Blank'}):
                    kind, strains = pair_kind(sheet)
                    pair = 'WT+WT' if sheet == 'WT' else sheet
                    catalogue[species, environment, pair] = {'kind': kind, 'strains': strains, 'nitrogen': nitrogen, 'carbon': carbon}
    for (species, environment, pair), first in sorted(catalogue.items()):
        group = groups[species, environment, pair]
        observed = [r for r in group if r['status'] == 'measured']
        values = [r['corrected_od600'] for r in observed]
        reps = defaultdict(list)
        for r in observed:
            reps[r['experiment']].append(r['corrected_od600'])
        rep_means = [statistics.mean(v) for v in reps.values()]
        mean = statistics.mean(values) if values else None
        result.append({'case_id': ':'.join(['NUTRIENT2026', species, environment, pair]),
                       'species': species, 'environment': environment, 'pair': pair,
                       'kind': first['kind'], 'strains': first['strains'],
                       'nitrogen': first['nitrogen'], 'carbon': first['carbon'],
                       'primary_environment': first['nitrogen'] != 'NoN',
                       'n_wells_measured': len(values), 'n_wells_missing': (12 if first['nitrogen'] == 'NoN' else 6) - len(values),
                       'experiments_present': sorted(reps), 'experiments_missing': sorted(set((1, 2, 3)) - set(reps)),
                       'mean_od600': mean, 'pooled_well_sem': sem(values),
                       'experiment_mean_od600': statistics.mean(rep_means) if rep_means else None,
                       'experiment_mean_sem': sem(rep_means),
                       'above_author_monoculture_threshold': mean > .02 if mean is not None else None,
                       'replicate_means': {str(k): statistics.mean(v) for k, v in sorted(reps.items())},
                       'endpoint_growth_is_not_gene_function_truth': True})
    return result


def case_inputs(observations):
    """Only assay identities enter inputs. Outcomes, coverage and OD never do."""
    result = []
    for row in observations:
        supplement = row['strains'][0] if row['kind'] == 'supplemented_monoculture' else None
        background = 'EcNR1 derivative of MG1655' if row['species'] == 'Ec' else 'Bacillus subtilis 168'
        organisms = [{'strain_id': row['species'] + ':' + strain, 'species': row['species'],
                      'background': background, 'deleted_genes': [GENES[row['species']][strain]]}
                     for strain in row['strains']]
        if not organisms:
            organisms = [{'strain_id': row['species'] + ':WT', 'species': row['species'],
                          'background': background, 'deleted_genes': []}]
        result.append({'case_id': row['case_id'], 'species': row['species'], 'kind': row['kind'],
                       'environment': row['environment'], 'pair': row['pair'],
                       'organisms': organisms,
                       'medium': {'base': BASE_MEDIUM, 'carbon_source': CARBON[row['carbon']],
                                  'nitrogen_source': NITROGEN[row['nitrogen']], 'source_concentration': '0.5% w/v each when present',
                                  'supplement': AUX[supplement] if supplement else None,
                                  'supplement_concentration_uM': (1000 if supplement == 'S' else 200) if supplement else None,
                                  'uptake_rates': None, 'oxygen_uptake': None},
                       'assay': {'temperature_c': 30, 'duration_hours': 96, 'volume_ul': 200, 'initial_od600': .005,
                                 'coculture_od_ratio': [1, 1] if row['kind'] == 'coculture' else None},
                       'primary_environment': row['nitrogen'] != 'NoN', 'protocol_source': PAPER})
    return result


def audit_author_exports(root, observations):
    observed = {(r['species'], r['environment'], r['pair']): r for r in observations}
    reports = {}
    for name in ('big_norm_df.csv', 'Supplemented_AA_data.csv'):
        with (root / 'raw/Generated Tables' / name).open() as stream:
            rows = list(csv.DictReader(stream))
        differences, seen = [], set()
        for row in rows:
            key = row['Species'], ' + '.join(p.strip() for p in row['Environment'].split('+')), row['Pair']
            if key in seen or key not in observed:
                raise DataError('duplicate or unmapped author summary: ' + str(key))
            seen.add(key)
            ours = observed[key]
            for column, own in [('Mean', 'mean_od600'), ('SEM', 'pooled_well_sem')]:
                value = float(row[column])
                if ours[own] is None or not math.isclose(value, ours[own], rel_tol=1e-9, abs_tol=1e-12):
                    differences.append({'key': list(key), 'column': column, 'author': value, 'raw_recomputed': ours[own]})
        reports[name] = {'rows': len(rows), 'differences': differences,
                         'mismatches_by_column': dict(Counter(d['column'] for d in differences))}
    return reports


def prepare(root, out):
    receipt = verify_inputs(root)
    out.mkdir(parents=True, exist_ok=False)
    mapping = sample_map(root / 'raw/Data/Sample Names.xlsx')
    measurements, coverage = read_measurements(root, mapping)
    observations = summarize(measurements)
    inputs = case_inputs(observations)
    audit = audit_author_exports(root, observations)
    for name, rows in [('measurements.jsonl', measurements), ('observations.jsonl', observations), ('inputs.jsonl', inputs)]:
        with (out / name).open('w') as stream:
            for row in rows:
                stream.write(json.dumps(row, sort_keys=True) + '\n')
    report = {'schema': 'gspa.nutrient_environment_data/1', 'source_commit': COMMIT, 'paper': PAPER,
              'source_hashes': {r['path']: r['sha256'] for r in receipt['files']}, 'paper_source': receipt['paper'],
              'implementation_sha256': sha(Path(__file__)), 'coverage': coverage,
              'measurements': len(measurements), 'observations': len(observations),
              'primary_observations': sum(o['primary_environment'] for o in observations),
              'observations_by_kind': dict(Counter(o['kind'] for o in observations)),
              'author_export_audit': audit,
              'interpretation': 'Measured 96-hour endpoint biomass, not maximum flux or independent gene attribution. No OD value selects a simulation input.',
              'unresolved_model_inputs': ['Concentrations do not determine flux bounds.', 'Oxygen transfer is unspecified.',
                                          'Biomass, maintenance, model strain differences and transport require explicit assumptions.'],
              'outputs': {p.name: sha(p) for p in sorted(out.iterdir()) if p.is_file()}}
    (out / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    verify_inputs(root)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('acquire', 'prepare'))
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    if args.action == 'acquire':
        acquire(args.root)
    elif args.out is None:
        parser.error('prepare requires --out')
    else:
        result = prepare(args.root, args.out)
        print(json.dumps({k: result[k] for k in ('measurements', 'observations', 'primary_observations', 'observations_by_kind')}))


if __name__ == '__main__':
    main()
