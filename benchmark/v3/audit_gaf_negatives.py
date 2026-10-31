#!/usr/bin/env python3
"""Stream a pinned GAF snapshot and preserve explicit negatives for source review.

This inventories assertions, not false-function labels. Relations, annotation
extensions, forms and references remain attached. An accession join to an older
sequence-verified mapping is a candidate join, not a proof of snapshot identity.
"""
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
from pathlib import Path
import re
import subprocess
import time

import evidence_codes as EV


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(4 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def negative_record(raw, line_number):
    fields = raw.decode('utf-8').rstrip('\r\n').split('\t')
    if len(fields) != 17:
        raise ValueError(f'negative GAF row {line_number}: expected 17 fields, got {len(fields)}')
    qualifiers = fields[3].split('|')
    if 'NOT' not in qualifiers:
        return None
    if len(qualifiers) != 2 or qualifiers[0] != 'NOT' or not qualifiers[1]:
        raise ValueError(f'negative GAF row {line_number}: expected NOT|relation for GAF 2.2')
    if not re.fullmatch(r'GO:\d{7}', fields[4]) or fields[8] not in ('F', 'P', 'C'):
        raise ValueError(f'negative GAF row {line_number}: invalid GO identifier/aspect')
    taxa = fields[12].split('|')
    if not 1 <= len(taxa) <= 2 or any(not re.fullmatch(r'taxon:\d+', t) for t in taxa):
        raise ValueError(f'negative GAF row {line_number}: invalid taxa')
    if not fields[5] or not fields[6] or not fields[1]:
        raise ValueError(f'negative GAF row {line_number}: missing identity, evidence or reference')
    return {'line': line_number, 'db': fields[0], 'accession': fields[1], 'symbol': fields[2],
            'relation': qualifiers[1], 'go_id': fields[4], 'references': fields[5].split('|'),
            'evidence': fields[6], 'evidence_category': EV.category(fields[6]), 'with_from': fields[7],
            'aspect': fields[8], 'object_name': fields[9], 'object_type': fields[11],
            'taxid': taxa[0][6:], 'interacting_taxid': taxa[1][6:] if len(taxa) == 2 else None,
            'date': fields[13], 'assigned_by': fields[14], 'extension': fields[15], 'form_id': fields[16],
            'raw_sha256': hashlib.sha256(raw).hexdigest(), 'gaf_row': fields}


def match_index(path):
    result = defaultdict(list)
    if path is not None:
        with Path(path).open() as fh:
            for row in csv.DictReader(fh, delimiter='\t'):
                result[row['swissprot_accession']].append(row)
    return result


def candidate_joins(record, matches):
    if record['db'] != 'UniProtKB':
        return []
    out = []
    for match in matches.get(record['accession'], []):
        reasons = []
        if match['status'] != 'matched':
            reasons.append('prior_sequence_match_' + match['status'])
        if match['taxid'] != record['taxid']:
            reasons.append('taxon_differs')
        if record['object_type'] != 'protein':
            reasons.append('not_a_protein')
        if record['form_id'] and record['form_id'] != 'UniProtKB:' + record['accession']:
            reasons.append('specific_gene_product_form')
        if record['extension']:
            reasons.append('annotation_extension_requires_context')
        if record['interacting_taxid']:
            reasons.append('interaction_context')
        if record['evidence'] not in EV.EXPERIMENTAL_GO:
            reasons.append('not_experimental')
        out.append({'genome_id': match['genome_id'], 'protein_id': match['protein_id'],
                    'prior_identity_basis': match['basis'], 'exclusion_reasons': reasons,
                    'status': 'excluded' if reasons else 'needs_snapshot_and_reference_review'})
    return out


def audit(source, out, matches_path=None):
    out.mkdir(parents=True, exist_ok=False)
    stamp = source.stat()
    code_sha = sha(__file__)
    evidence_code_sha = sha(EV.__file__)
    print('Hashing compressed source before scan', flush=True)
    source_sha = sha(source)
    matches_sha = sha(matches_path) if matches_path else None
    matches = match_index(matches_path)
    counts = {k: Counter() for k in ('evidence', 'evidence_category', 'relation', 'aspect', 'taxid', 'db', 'object_type')}
    headers, negatives, experimental, joined = [], 0, 0, 0
    started, last = time.monotonic(), time.monotonic()
    line_number, uncompressed_bytes = 0, 0
    proc = subprocess.Popen(['gzip', '-dc', '--', str(source)], stdout=subprocess.PIPE)
    try:
        with (out / 'negatives.jsonl').open('w') as all_out, (out / 'cohort-candidates.jsonl').open('w') as cohort_out:
            for line_number, raw in enumerate(proc.stdout, 1):
                uncompressed_bytes += len(raw)
                if raw.startswith(b'!'):
                    headers.append(raw.decode('utf-8').rstrip('\r\n'))
                    continue
                # The delimiter-qualified prefilter avoids splitting billions of positive rows.
                if b'\tNOT|' in raw or b'\tNOT\t' in raw:
                    record = negative_record(raw, line_number)
                    if record is not None:
                        negatives += 1
                        experimental += record['evidence'] in EV.EXPERIMENTAL_GO
                        for key, count in counts.items():
                            count[record[key]] += 1
                        record['candidate_joins'] = candidate_joins(record, matches)
                        all_out.write(json.dumps(record, sort_keys=True) + '\n')
                        if record['candidate_joins']:
                            joined += 1
                            cohort_out.write(json.dumps(record, sort_keys=True) + '\n')
                if line_number % 1000000 == 0 and time.monotonic() - last >= 30:
                    last = time.monotonic()
                    print(json.dumps({'lines': line_number, 'decompressed_bytes': uncompressed_bytes,
                                      'negative_records': negatives, 'experimental': experimental,
                                      'cohort_records': joined, 'elapsed_seconds': round(last - started)}), flush=True)
        if proc.wait() != 0:
            raise ValueError('gzip failed or source CRC is invalid; scan cannot be certified')
    finally:
        proc.stdout.close()
        if proc.poll() is None:
            proc.terminate()
            proc.wait()
    if '!gaf-version: 2.2' not in headers:
        raise ValueError('audit requires declared GAF 2.2')
    print('Verifying compressed source after complete scan', flush=True)
    if source.stat().st_size != stamp.st_size or source.stat().st_mtime_ns != stamp.st_mtime_ns or sha(source) != source_sha:
        raise ValueError('source changed during scan')
    if sha(__file__) != code_sha or sha(EV.__file__) != evidence_code_sha or (matches_path and sha(matches_path) != matches_sha):
        raise ValueError('audit code or identity mapping changed during scan')
    report = {'schema': 'gspa.gaf-negative-inventory/1', 'status': 'verified', 'headers': headers,
              'source': {'path': str(source), 'bytes': stamp.st_size, 'sha256': source_sha},
              'code_sha256': code_sha, 'evidence_codes_sha256': evidence_code_sha,
              'match_table': str(matches_path) if matches_path else None,
              'match_table_sha256': matches_sha, 'lines': line_number, 'uncompressed_bytes': uncompressed_bytes,
              'negative_records': negatives, 'experimental_negative_records': experimental,
              'cohort_candidate_records': joined, 'counts': {k: dict(sorted(v.items())) for k, v in counts.items()},
              'output_sha256': {name: sha(out / name) for name in ('negatives.jsonl', 'cohort-candidates.jsonl')},
              'seconds': time.monotonic() - started,
              'interpretation': 'Inventory, not accepted training truth. Candidate joins require snapshot/reference review. Relations, forms and extensions are not flattened.'}
    (out / 'inventory.json').write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    return report


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source', required=True, type=Path)
    ap.add_argument('--out', required=True, type=Path)
    ap.add_argument('--matches', type=Path)
    args = ap.parse_args()
    result = audit(args.source.resolve(), args.out.resolve(), args.matches.resolve() if args.matches else None)
    print(json.dumps({k: result[k] for k in ('status', 'lines', 'negative_records', 'experimental_negative_records', 'cohort_candidate_records')}), flush=True)
