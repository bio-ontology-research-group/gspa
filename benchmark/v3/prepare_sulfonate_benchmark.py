#!/usr/bin/env python3
"""Extract the complete Eichhorn 2000 Table 3, retaining missing/ordinal outcomes.

Inputs are retained, line-numbered publisher web extractions, not invented tables.
This is a developmental phenotype dataset; it is not masked gene-function truth.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import re

DOI = '10.1128/JB.182.10.2687-2695.2000'
GENOTYPES = ['None', 'tauA', 'tauBC', 'ssuA', 'ssuCB', 'ssuABC', 'ssuA tauA',
             'ssuCB tauBC', 'ssuA tauBC', 'ssuCB tauA', 'ssuAC tauB', 'ssuAB tauC', 'ssuB tauAC']
# Table 3's complete row catalogue, plus the tested compound omitted from its table.
SUBSTRATES = ['Sulfate', 'Taurine', 'PIPES', 'MOPS', 'HEPES', '4-Phenyl-1-butanesulfonate',
              '2-(4-Pyridyl)ethanesulfonate', 'Isethionate', '1,3-Dioxo-2-isoindolineethanesulfonate',
              '3-Aminopropanesulfonate', 'Ethanesulfonate', 'Propanesulfonate', 'Butanesulfonate',
              'Pentanesulfonate', 'Hexanesulfonate', 'Octanesulfonate', 'Decanesulfonate',
              'N-Phenyltaurine', 'Sulfoacetate', 'Dodecanesulfonate']
CATEGORIES = {'−': {'ordinal': 0, 'relative_yield_interval': '<5%'},
              '+': {'ordinal': 1, 'relative_yield_interval': '5–33%'},
              '++': {'ordinal': 2, 'relative_yield_interval': '33–66%'},
              '+++': {'ordinal': 3, 'relative_yield_interval': '>66%'}}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def genes(label):
    if label == 'None':
        return []
    result = []
    for token in label.split():
        match = re.fullmatch(r'(ssu|tau)([ABC]+)', token)
        if not match:
            raise ValueError('unrecognized Table 3 genotype: ' + label)
        result.extend(match[1] + suffix for suffix in match[2])
    if len(result) != len(set(result)):
        raise ValueError('duplicated gene in genotype: ' + label)
    return sorted(result)


def source_lines(extraction):
    text = json.loads(Path(extraction).read_text())
    if not isinstance(text, str) or 'jb.182.10.2687-2695.2000' not in text.lower():
        raise ValueError('wrong publisher extraction')
    return [(int(m[1]), m[2]) for line in text.splitlines()
            if (m := re.fullmatch(r'L(\d+): ?(.*)', line))]


def parse_table(lines):
    selected, started, headers = {}, False, None
    for number, text in lines:
        if text.startswith('Table 3. Growth phenotypes'):
            if started:
                raise ValueError('duplicated table')
            started = True
        if not started:
            continue
        if text.startswith('Net growth relative'):
            if '<5%' not in text or 'three independent experiments' not in text:
                raise ValueError('unrecognized category definition')
            break
        cells = [s.strip() for s in text.split('|')]
        if len(cells) == 13 and cells[0] == 'None':
            if headers is not None or cells != GENOTYPES:
                raise ValueError('unexpected or duplicated genotype header')
            headers = cells
        elif cells[0] in SUBSTRATES:
            if headers is None or len(cells) != 14:
                raise ValueError('table row width/header mismatch')
            if cells[0] in selected or any(c not in CATEGORIES for c in cells[1:]):
                raise ValueError('duplicate row or invalid growth category')
            selected[cells[0]] = {'line': number, 'values': cells[1:]}
    else:
        raise ValueError('missing table or category footnote')
    if headers != GENOTYPES or set(selected) != set(SUBSTRATES[:-1]):
        raise ValueError('incomplete or changed Table 3 row catalogue')
    if not any('dodecanesulfonate' in text and 'Growth was tested' in text for _, text in lines):
        raise ValueError('missing evidence for tested compound omitted from Table 3')
    return selected


def extract(table_path, methods_path):
    lines = source_lines(table_path)
    table = parse_table(lines)
    methods = '\n'.join(text for _, text in source_lines(methods_path))
    for required in ('MC4100 trp-1', '4 μg of tryptophan per ml', '250 μM', '37°C',
                     '150 μl of culture', 'every 5 min', 'aerobic growth conditions'):
        if required not in methods:
            raise ValueError('missing study method: ' + required)
    inputs, observations = [], []
    for si, substrate in enumerate(SUBSTRATES):
        for gi, genotype in enumerate(GENOTYPES):
            case_id = f'SULFONATE2000:S{si:02d}:G{gi:02d}'
            inputs.append({'case_id': case_id, 'organism': 'Escherichia coli',
                           'strain_background': 'EC1250 (MC4100 trp-1)',
                           'genotype_label': genotype,
                           'reported_deleted_genes': genes(genotype),
                           'deleted_genes': None if genotype == 'ssuAB tauC' else genes(genotype),
                           'genotype_status': 'table_text_conflict' if genotype == 'ssuAB tauC' else 'table_reported',
                           'sulfur_source': substrate, 'sulfur_source_uM': 250,
                           'medium': 'sulfur-free M63', 'tryptophan_ug_ml': 4,
                           'temperature_C': 37, 'culture_volume_ul': 150,
                           'oxygenation': 'aerobic; 30-second shaking before each five-minute measurement'})
            row = table.get(substrate)
            symbol = row['values'][gi] if row else None
            observations.append({'case_id': case_id,
                                 'status': 'observed_ordinal' if row else 'not_reported_in_table',
                                 'genotype_status': 'table_text_conflict' if genotype == 'ssuAB tauC' else 'table_reported',
                                 'observation_role': 'wild_type_normalization_reference' if genotype == 'None' else
                                                     'sulfate_control' if substrate == 'Sulfate' else 'deletion_response',
                                 'symbol': symbol, **(CATEGORIES[symbol] if symbol else {'ordinal': None, 'relative_yield_interval': None}),
                                 'source': {'doi': DOI, 'table': 3, 'row': substrate,
                                            'column': genotype, 'rendered_line': row['line'] if row else None},
                                 'normalization': 'net growth relative to same-substrate wild type; ssuD tauD double mutant is background',
                                 'minimum_independent_experiments': 3 if row else None,
                                 'raw_replicates_available': False,
                                 'interpretation': 'conditional growth phenotype; not enzyme-function or universal negative truth'})
    report = {'schema': 'gspa.sulfonate_table3/1', 'doi': DOI,
              'source_kind': 'retained rendered publisher extraction; original HTML bytes unavailable',
              'sources': {str(p): sha(p) for p in (table_path, methods_path)},
              'catalogue_cases': len(inputs), 'reported_cells': sum(o['symbol'] is not None for o in observations),
              'missing_cells': sum(o['symbol'] is None for o in observations),
              'uncertain_genotype_cases': sum(i['genotype_status'] != 'table_reported' for i in inputs),
              'categories': dict(Counter(o['symbol'] or 'unreported' for o in observations)),
              'scope': 'developmental known-genotype phenotype benchmark; no independent function-recovery labels',
              'limits': ['Table header ssuAB tauC contradicts the text saying this construct was not obtained; its genotype is unresolved.',
                         'Dodecanesulfonate is named among tested compounds but absent from Table 3.',
                         'Qualitative categories are normalized yields, not absolute growth rates or replicate measurements.',
                         'The ssuD/tauD controls are described but are not extra measured Table 3 columns.',
                         'MC4100-derived EC1250 is not the MG1655 background of iML1515.',
                         'The study predates model curation; independence from model/reference knowledge is unproven.',
                         'No gene/function negatives, null-locus labels, catalytic roles or cross-feeding labels are inferred.']}
    return inputs, observations, report


def write_rows(path, rows):
    path.write_text(''.join(json.dumps(r, sort_keys=True, ensure_ascii=False) + '\n' for r in rows))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--table-extraction', type=Path, required=True)
    p.add_argument('--methods-extraction', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    rows, outcomes, report = extract(a.table_extraction, a.methods_extraction)
    a.out.mkdir(parents=True, exist_ok=False)
    write_rows(a.out / 'inputs.jsonl', rows)
    write_rows(a.out / 'observations.jsonl', outcomes)
    report['outputs'] = {name: sha(a.out / name) for name in ('inputs.jsonl', 'observations.jsonl')}
    (a.out / 'report.json').write_text(json.dumps(report, indent=2, ensure_ascii=False) + '\n')
    print(json.dumps({k: report[k] for k in ('catalogue_cases', 'reported_cells', 'missing_cells', 'categories')}))


if __name__ == '__main__':
    main()
