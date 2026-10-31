#!/usr/bin/env python3
"""Retain every Table 1 cell from Thanapipatsiri et al. 2015, without model labels."""
import argparse
from decimal import Decimal
import hashlib
import json
from pathlib import Path
import re
import xml.etree.ElementTree as ET


def prepare(path):
    raw = path.read_bytes()
    root = ET.fromstring(raw)
    ids = {x.attrib.get('pub-id-type'): x.text for x in root.findall('.//article-meta/article-id')}
    if ids.get('doi') != '10.1186/s12934-015-0335-0' or ids.get('pmid') != '26376792':
        raise ValueError('unexpected primary article')
    table = root.find(".//table-wrap[@id='Tab1']")
    if table is None:
        raise ValueError('missing production table')
    rows = [[''.join(c.itertext()).strip() for c in row] for row in table.findall('.//tr')]
    if rows[1] != ['Germicidin A', 'Germicidin B7', 'Total germicidins', 'Flaviolin'] or len(rows) != 8:
        raise ValueError('unexpected table columns or population')
    expected_heading = 'Relative production (mAU × 103)/mg dry weight'
    if ' '.join(rows[0][1].split()) != expected_heading:
        raise ValueError('unexpected production units')
    paragraphs = [''.join(p.itertext()) for p in root.findall('.//body//p')]
    checks = ['Cultures were grown for 5', 'Each experiment was carried out in duplicate.',
              '50\u00a0ml of medium in 250\u00a0ml flasks', '30\u00a0°C', '250\u00a0rpm']
    if any(not any(phrase in p for p in paragraphs) for phrase in checks):
        raise ValueError('curated culture metadata needs source reconciliation')
    metadata_sources = [{'body_paragraph_index': i, 'sha256': hashlib.sha256(p.encode()).hexdigest()}
                        for i, p in enumerate(paragraphs) if any(phrase in p for phrase in checks)]
    records = []
    for i, row in enumerate(rows[2:]):
        if len(row) != 5:
            raise ValueError('unexpected row width')
        for product, value in zip(rows[1], row[1:]):
            record = {'strain': row[0], 'product': product, 'raw_value': value,
                      'source_table': 'Tab1', 'source_data_row': i + 1,
                      'units': '(mAU × 10^3)/mg dry weight', 'replicates_reported': 2}
            if value == 'ND':
                record.update(status='not_determined', value=None, reported_uncertainty=None)
            elif value == 'None':
                record.update(status='not_detected', value=None, reported_uncertainty=None,
                              detection_limit=None)
            else:
                match = re.fullmatch(r'([0-9.]+)\s*±\s*([0-9.]+)', value)
                if not match:
                    raise ValueError('unexpected measured value: ' + value)
                mean, uncertainty = (Decimal(v) for v in match.groups())
                if not mean.is_finite() or not uncertainty.is_finite() or mean <= 0 or uncertainty < 0:
                    raise ValueError('invalid measured value')
                record.update(status='detected', value=float(mean), reported_uncertainty=float(uncertainty),
                              uncertainty_definition='not inferred; preserve reported ± value')
            records.append(record)
    if len({(r['strain'], r['product']) for r in records}) != 24:
        raise ValueError('duplicated or missing table cells')
    genotype_table = root.find(".//table-wrap[@id='Tab2']")
    if genotype_table is None:
        raise ValueError('missing genotype table')
    genotypes = []
    for row in genotype_table.findall('.//tr'):
        cells = [''.join(c.itertext()).strip() for c in row]
        if cells and any(s in cells[0] for s in ('M1152', 'M1314', 'M1316', 'M1317')):
            if len(cells) != 3:
                raise ValueError('unexpected genotype row')
            genotypes.append({'strain': cells[0], 'description': cells[1], 'source': cells[2]})
    if len(genotypes) != 4:
        raise ValueError('unexpected genotype population')
    return {'schema': 'gspa.germicidin_primary_observations/1', 'doi': ids['doi'],
            'source_sha256': hashlib.sha256(raw).hexdigest(), 'records': records, 'genotypes': genotypes,
            'metadata_provenance': {'method': 'manual curation from full primary article; checked phrases and source hashes',
                                    'paragraphs': metadata_sources, 'raw_units_header': rows[0][1]},
            'culture': {'medium': 'GYM (glucose-yeast extract-malt extract)', 'duration_days': 5,
                        'temperature_celsius': 30, 'shaking_rpm': 250, 'culture_ml': 50, 'flask_ml': 250,
                        'matrix': 'culture supernatant'},
            'interpretation': 'Production observations, not enzyme negatives or steady-state fluxes. '
                'GYM composition, secretion transport and congener mappings are not represented by the model diagnostic.'}


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--article', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    args = p.parse_args()
    if args.out.exists():
        p.error('output already exists')
    args.out.write_text(json.dumps(prepare(args.article), indent=2, ensure_ascii=False) + '\n')
