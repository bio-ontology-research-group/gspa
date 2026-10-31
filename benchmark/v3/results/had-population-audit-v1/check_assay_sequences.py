#!/usr/bin/env python3
"""Check all sheet GI/accession identities against independently archived sequences."""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

import openpyxl


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--intake', type=Path, required=True)
    p.add_argument('--cache', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--fetch', action='store_true')
    args = p.parse_args()
    if args.out.exists():
        raise FileExistsError(args.out)
    args.cache.mkdir(parents=True, exist_ok=True)
    raw = args.intake / 'mirror/data/raw/phosphatase'
    workbook = raw / 'pnas.1423570112.sd02.xlsx'
    mirror = raw / 'phosphatase_sequence_file.txt'
    with mirror.open() as handle:
        sequence_rows = list(csv.DictReader(handle, delimiter='\t'))
    sequences = {r['Entry']: r['Sequence'] for r in sequence_rows}
    assert len(sequence_rows) == len(sequences) == 218
    # Recheck the previously retrieved primary sequences, not just their summary.
    primary_sequences, primary_paths = {}, []
    for path in sorted((args.intake / 'sequence-verification-v1').glob('batch-*.json')):
        primary_paths.append(path)
        for r in json.loads(path.read_text())['results']:
            accession = r['primaryAccession']
            if r['entryType'] == 'Inactive':
                fasta = args.intake / 'sequence-resolution-v1' / (r['extraAttributes']['uniParcId'] + '.fasta')
                primary_paths.append(fasta)
                primary_sequences[accession] = ''.join(fasta.read_text().splitlines()[1:])
            else:
                primary_sequences[accession] = r['sequence']['value']
    assert all(primary_sequences[a] == seq for a, seq in sequences.items())
    w = openpyxl.load_workbook(workbook, read_only=True, data_only=True)
    lookup = {str(row[0]).strip(): row[1] for row in
              list(w['EFI-UniProt lookup '].values)[1:] if row[0] is not None}
    sheets = [(s.title.strip(), s['D2'].value) for s in w if s.title.strip().isdigit()]
    gis = sorted({str(identity) for _, identity in sheets if isinstance(identity, int)})
    fetched, receipts, failures = {}, [], []
    for start in range(0, len(gis), 32):
        batch = gis[start:start + 32]
        path = args.cache / ('batch-%02d.xml' % (start // 32))
        receipt_path = path.with_suffix('.receipt.json')
        url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?' + urllib.parse.urlencode(
            dict(db='protein', id=','.join(batch), rettype='gp', retmode='xml'))
        if not path.exists() and args.fetch:
            receipt = dict(url=url, requested_gis=batch,
                           retrieved_utc=datetime.now(timezone.utc).isoformat())
            try:
                with urllib.request.urlopen(url, timeout=40) as response:
                    data = response.read()
                    receipt.update(http_status=response.status, final_url=response.url,
                                   content_type=response.headers.get('Content-Type'))
                with path.open('xb') as handle:
                    handle.write(data)
                receipt.update(sha256=sha(path), bytes=len(data))
            except Exception as error:
                receipt['error'] = str(error)
            with receipt_path.open('x') as handle:
                json.dump(receipt, handle, indent=2)
                handle.write('\n')
            time.sleep(0.4)
        if not path.exists():
            failures.append(dict(batch=start // 32, reason='no_response_file'))
            continue
        receipt = json.loads(receipt_path.read_text())
        assert receipt['url'] == url and receipt['requested_gis'] == batch
        assert receipt['sha256'] == sha(path)
        receipts.append(dict(file=path.name, sha256=sha(path), receipt=receipt))
        try:
            root = ET.fromstring(path.read_bytes())
            for record in root.findall('GBSeq'):
                for identifier in record.findall('GBSeq_other-seqids/GBSeqid'):
                    if identifier.text.startswith('gi|'):
                        gi = identifier.text[3:]
                        if gi not in batch:
                            continue
                        assert gi not in fetched, gi
                        fetched[gi] = dict(accession_version=record.findtext('GBSeq_accession-version'),
                            organism=record.findtext('GBSeq_organism'),
                            sequence=(record.findtext('GBSeq_sequence') or '').upper(), file=path.name)
        except ET.ParseError as error:
            failures.append(dict(batch=start // 32, reason=str(error)))
        print('batch', start // 32, 'cumulative GI records', len(fetched), flush=True)
    results = []
    for efi, identity in sheets:
        accession = lookup[efi]
        row = dict(efi_id=efi, assay_identity=identity, lookup_accession=accession,
                   lookup_sequence_length=len(sequences[accession]))
        if isinstance(identity, int):
            record = fetched.get(str(identity))
            if record is None:
                row['status'] = 'unresolved_gi'
            else:
                seq = record['sequence']
                row.update({k: v for k, v in record.items() if k != 'sequence'})
                row.update(assay_sequence_length=len(seq),
                    sequence_sha256=hashlib.sha256(seq.encode()).hexdigest(),
                    status='exact_sequence_match' if seq == sequences[accession] else 'sequence_mismatch')
        else:
            row['status'] = 'direct_accession_match' if identity == accession else 'accession_mismatch'
        results.append(row)
    from collections import Counter
    report = dict(status='identity_audit_only_no_function_labels',
        counts=dict(Counter(r['status'] for r in results)), results=results, receipts=receipts,
        failures=failures, source_doi='10.1073/pnas.1423570112',
        inputs_sha256={str(path.relative_to(args.intake)): sha(path)
                       for path in [workbook, mirror] + primary_paths},
        script_sha256=sha(Path(__file__)), limitations=[
            'Identity matches do not verify actual recombinant constructs, tags, mutations or sample handling.',
            'The workbook remains a third-party mirror without publisher byte-equality verification.',
            'Protein/organism names are retained from NCBI records, not accepted as genome-coordinate mappings.',
            'No database functional descriptions are imported as experimental labels.'])
    with args.out.open('x') as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write('\n')
    print(json.dumps(report['counts'], sort_keys=True))


if __name__ == '__main__':
    main()
