#!/usr/bin/env python3
"""Audit the complete HAD workbook population without accepting function labels."""
import argparse
from collections import Counter, defaultdict
from decimal import Decimal
import hashlib
import json
from pathlib import Path

import openpyxl


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--intake', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError(args.out)
    raw = args.intake / 'mirror/data/raw/phosphatase'
    metadata_path = raw / 'pnas.1423570112.sd01.xlsx'
    assays_path = raw / 'pnas.1423570112.sd02.xlsx'
    inventory_path = args.intake / 'raw-observations.jsonl'
    metadata = openpyxl.load_workbook(metadata_path, read_only=True, data_only=True)
    assays = openpyxl.load_workbook(assays_path, data_only=False)
    m2, m3 = defaultdict(list), defaultdict(list)
    for number, row in enumerate(metadata['Supplemental 2'].values, 1):
        if number > 1 and row[0] is not None:
            m2[str(row[0]).strip()].append(dict(row=number, accession=row[1],
                organism=row[2], screened=row[7]))
    for number, row in enumerate(metadata['Supplemental 3'].values, 1):
        if number > 1 and isinstance(row[6], (int, float)):
            m3[str(row[6]).strip()].append(dict(row=number,
                shared_accessions=row[3], gi=row[4], screened=row[9]))
    assay_ids = {s.title.strip() for s in assays if s.title.strip().isdigit()}
    lookup = {str(row[0]).strip(): row[1] for row in
              list(assays['EFI-UniProt lookup '].values)[1:] if row[0] is not None}
    yes2 = {k for k, rows in m2.items() if any(r['screened'] == 'Y' for r in rows)}
    yes3 = {k for k, rows in m3.items() if any(r['screened'] == 'Y' for r in rows)}
    conflicts = []
    for efi in sorted(set(m2) | set(m3)):
        flags2 = {r['screened'] for r in m2[efi]} & {'Y', 'N'}
        flags3 = {r['screened'] for r in m3[efi]} & {'Y', 'N'}
        if flags2 and flags3 and flags2 != flags3:
            conflicts.append(dict(efi_id=efi, accession=lookup.get(efi),
                has_assay=efi in assay_ids, supplemental_2=m2[efi], supplemental_3=m3[efi]))
    # Compare every saved observation to its originating workbook cell, including blanks.
    saved = {}
    with inventory_path.open() as handle:
        for line in handle:
            record = json.loads(line)
            key = (record['efi_id'], record['cell'])
            assert key not in saved, key
            saved[key] = record
    samples, extras, formulas, matched = {}, [], [], 0
    totals = Counter()
    for sheet in assays:
        efi = sheet.title.strip()
        if efi not in assay_ids:
            continue
        values, counts = [], Counter()
        for row in sheet['C7:N20']:
            for cell in row:
                record = saved.pop((efi, cell.coordinate))
                assert record['reported_absorbance'] == cell.value
                assert record['substrate_text'] == cell.comment.text.strip()
                assert record['uniprot'] == lookup[efi]
                matched += 1
                if record['blank']:
                    counts['blank_cells'] += 1
                    continue
                # These are descriptive signal strata, not accepted training labels.
                value = Decimal(str(cell.value))
                values.append(value)
                counts['substrate_cells'] += 1
                counts['below_0.17'] += value < Decimal('0.17')
                counts['within_0.17_to_0.23_inclusive'] += Decimal('0.17') <= value <= Decimal('0.23')
                counts['above_0.23'] += value > Decimal('0.23')
                counts['above_nominal_0.2'] += value > Decimal('0.2')
                counts['above_linear_response_0.7'] += value > Decimal('0.7')
                counts['exact_zero'] += value == 0
        assert len(values) == 167
        assert counts['blank_cells'] == 1
        for row in sheet:
            for cell in row:
                if cell.data_type == 'f':
                    formulas.append(dict(sheet=sheet.title, cell=cell.coordinate, formula=cell.value))
                if cell.row >= 7 and not (3 <= cell.column <= 14 and cell.row <= 20):
                    if isinstance(cell.value, (int, float)) and not isinstance(cell.value, bool):
                        extras.append(dict(sheet=sheet.title, cell=cell.coordinate, value=cell.value))
        samples[efi] = dict(accession=lookup[efi], supplemental_2=m2[efi],
            supplemental_3=m3[efi], signal_counts=dict(counts),
            metadata_says_screened=efi in yes2 | yes3,
            header_identity={c: sheet[c].value for c in ('C2', 'D2', 'E2', 'C4')})
        totals.update(counts)
    assert not saved
    assert matched == len(assay_ids) * 168
    average_formulas = {c.coordinate: c.value for row in assays['Averages'] for c in row
                        if c.data_type == 'f'}
    report = dict(status='quarantined_no_labels_accepted',
        source_doi='10.1073/pnas.1423570112',
        source_mirror_commit='627556e265e2a52d39753e684c05b550d53a9be4',
        publisher_byte_equality_verified=False,
        paper_reported_proteins=217, assay_sheets=len(assay_ids),
        supplemental_2_screened=len(yes2), supplemental_3_screened=len(yes3),
        screened_union=len(yes2 | yes3),
        screened_union_with_assays=len((yes2 | yes3) & assay_ids),
        screened_union_without_assays=sorted((yes2 | yes3) - assay_ids),
        assays_without_screened_metadata=sorted(assay_ids - (yes2 | yes3)),
        intertable_screening_conflicts=conflicts,
        duplicate_efi_rows={name: {k: v for k, v in table.items() if len(v) > 1}
                            for name, table in [('supplemental_2', m2), ('supplemental_3', m3)]},
        verified_original_inventory_cells=matched, signal_counts=dict(totals),
        assay_formulas=formulas, numeric_cells_outside_assay_block_from_row_7=extras,
        averages_sheet_formulas=average_formulas,
        samples=samples,
        interpretation=[
            '217 sheets have a Y in at least one metadata table; this is a metadata intersection, not proof of the paper population.',
            'Both tables mark EFI 502329 N despite its measured assay sheet; keep the sheet quarantined, not silently excluded.',
            'EFI 900190 is Y only in Supplemental 3 but has no assay sheet; do not synthesize an observation.',
            'Supplemental 3 is also a propagated annotation table; only screening flags and identities were audited, not imported functional truth.',
            'No individual replicate readings or replicate-specific errors are identifiable in the retained assay workbook.',
            'Averages contains across-enzyme formulas and is not a second experimental replicate.',
            'The 0.17/0.23 band is descriptive sensitivity analysis using the paper-wide variability; it is not a per-pair confidence interval.',
            'Low assay signal does not establish NOT-GO, absence of an EC activity, or a no-valid-gene locus.',
            'High endpoint signal does not establish physiological substrate preference or an EC assignment.',
        ])
    report['inputs_sha256'] = {str(p.relative_to(args.intake)): digest(p)
        for p in (metadata_path, assays_path, inventory_path, args.intake / 'article.html',
                  args.intake / 'metadata.json')}
    report['audit_script_sha256'] = digest(Path(__file__))
    report['openpyxl_version'] = openpyxl.__version__
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open('x') as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write('\n')
    print(json.dumps({k: report[k] for k in ('status', 'assay_sheets',
        'screened_union', 'screened_union_with_assays', 'screened_union_without_assays',
        'assays_without_screened_metadata', 'verified_original_inventory_cells', 'signal_counts')}, indent=2))


if __name__ == '__main__':
    main()
