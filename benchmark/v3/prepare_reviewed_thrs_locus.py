#!/usr/bin/env python3
"""Add the manually reviewed PMID:15240874 ThrRS competition to a cohort copy.

This is one training locus, not an independent benchmark or sufficient null-model
training. The two source assertions below were reviewed against the full paper:
ThrRS-cat aminoacylates tRNAThr; ThrRS-ed lacks threonylation activity. All other
GAF assertions remain unaccepted. Original candidates and splits are unchanged.
"""
import argparse
import csv
from datetime import datetime
import io
import hashlib
import json
from pathlib import Path
import sys

import audit_gaf_negatives as A
import match_swissprot_labels as M
import prepare_training_data as P

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
import train_context_model as T

GENOME = 'GCA_000007005.1'
TERM = 'GO:0004829'
MEMBERS = {'Q97VW8': ('AAK42622.1', 1), 'Q980D1': ('AAK40712.1', 0)}
REFERENCE = 'PMID:15240874'
SOURCE = ('PMID:15240874, Materials and Methods, Results (Elements in M. jannaschii ThrRS and '
          'S. solfataricus ThrRS-cat That Are Responsible for the Specificity of Threonylation), '
          'Table 1 and Fig. 3A inset; cloned SSO2486/SSO0384 from ATCC 35092D. '
          'Complete only for the two exported GO:0004829/UNSPECIFIED candidates, '
          'not all genomic proteins or every assay condition.')


def verify_entry(entry_path, accession, protein, proteins, snapshot_date):
    entries = list(M.iter_entries(entry_path))
    if len(entries) != 1:
        raise ValueError('expected one UniProt record')
    e = entries[0]
    if e['accessions'][0] != accession or e['taxid'] != '273057':
        raise ValueError('UniProt identity/taxon differs')
    if protein not in e['xrefs'] or proteins.get(protein) != e['sequence'] or not e['sequence']:
        raise ValueError('full sequence and versioned accession must both match')
    dates = [l for l in entry_path.read_text().splitlines() if l.startswith('DT   ') and 'sequence version' in l]
    if len(dates) != 1 or datetime.strptime(dates[0][5:16], '%d-%b-%Y').date() > snapshot_date:
        raise ValueError('sequence version postdates the annotation snapshot')
    return {'accession': accession, 'protein_id': protein, 'taxid': e['taxid'],
            'sequence_sha256': hashlib.sha256(e['sequence'].encode()).hexdigest(),
            'length': len(e['sequence']), 'sequence_version': dates[0],
            'basis': 'versioned_accession_and_full_sequence; last sequence update predates GAF'}


def build(original_manifest, provenance_path, inventory_path, references, out):
    out.mkdir(parents=True, exist_ok=False)
    old = json.loads(original_manifest.read_text())
    inputs = old['inputs']
    inventory = json.loads(inventory_path.read_text())
    if inventory.get('schema') != 'gspa.gaf-negative-inventory/1' or inventory.get('status') != 'verified':
        raise ValueError('complete GAF inventory required')
    candidates = inventory_path.parent / 'cohort-candidates.jsonl'
    if A.sha(candidates) != inventory['output_sha256']['cohort-candidates.jsonl']:
        raise ValueError('GAF cohort records changed')
    found = [r for r in map(json.loads, candidates.read_text().splitlines())
             if r['accession'] == 'Q980D1' and r['go_id'] == TERM and r['relation'] == 'enables'
             and r['evidence'] == 'IDA' and REFERENCE in r['references']]
    if len(found) != 1 or not any(j['genome_id'] == GENOME and j['protein_id'] == MEMBERS['Q980D1'][0]
                                and j['status'] == 'needs_snapshot_and_reference_review' for j in found[0]['candidate_joins']):
        raise ValueError('expected unqualified exact GAF assertion absent')
    assembly = json.loads(provenance_path.read_text())['organisms']['273057']['assembly']
    if assembly['used_accession'] != GENOME:
        raise ValueError('assembly identity differs')
    fasta = Path(assembly['directory']) / 'protein.faa'
    if A.sha(fasta) != assembly['member_sha256']['protein.faa']:
        raise ValueError('acquired genome proteins changed')
    proteins = M.PL.fasta_sequences(fasta)
    snapshot = [s for s in inventory['headers'] if s.startswith('!date-generated: ')][0][17:27]
    snapshot_date = datetime.strptime(snapshot, '%Y-%m-%d').date()
    sequence_checks = [verify_entry(references / (a + '.txt'), a, p, proteins, snapshot_date)
                       for a, (p, _) in MEMBERS.items()]
    # Full-text source is an archived, manually read retrieval, not a search snippet.
    paper = references / 'PMC478561.html'
    if A.sha(paper) != '194c62be69bf7225fec93286c8e26cd2e9e9218714354e061569590f87e21842':
        raise ValueError('reviewed full-text artifact differs; source needs review again')
    files = [original_manifest, provenance_path, inventory_path, candidates, fasta, paper,
             Path(__file__), Path(A.__file__), Path(M.__file__), Path(M.PL.__file__), Path(P.__file__), Path(T.__file__),
             Path(M.EV.__file__), references / 'Q97VW8.txt', references / 'Q980D1.txt', Path(old['table'])]
    if A.sha(old['table']) != old['table_sha256']:
        raise ValueError('original training table changed')
    for key in ('contract', 'labels', 'splits', 'label_policy', 'registry', 'go_obo', 'example_clusters'):
        if inputs.get(key):
            p = Path(inputs[key])
            if inputs.get(key + '_sha256') and A.sha(p) != inputs[key + '_sha256']:
                raise ValueError('original input hash differs: ' + key)
            files.append(p)
    files.extend(p for p in Path(inputs['feature_dir']).rglob('*') if p.is_file())
    identity = {str(p.resolve()): A.sha(p) for p in files}
    with Path(inputs['labels']).open() as fh:
        labels = list(csv.DictReader(fh, delimiter='\t'))
    if any(r.get('label', '1') != '1' or r.get('qualifier', '') for r in labels):
        raise ValueError('this augmentation requires the original positive-only labels')
    positive = [r for r in labels if r['genome_id'] == GENOME and r['protein_id'] == MEMBERS['Q97VW8'][0] and r['function_id'] == TERM]
    if len(positive) != 1 or positive[0]['evidence'] != 'IDA':
        raise ValueError('expected original positive assertion absent')
    for r in labels:
        r['label'] = '1'
        r['reference'] = REFERENCE if r is positive[0] else ''
    labels.append(dict(positive[0], protein_id=MEMBERS['Q980D1'][0], swissprot_accession='Q980D1', label='0'))
    stream = io.StringIO()
    writer = csv.DictWriter(stream, fieldnames=list(labels[0]), delimiter='\t', lineterminator='\n')
    writer.writeheader(); writer.writerows(labels)
    (out / 'labels.tsv').write_text(stream.getvalue())
    policy = json.loads(Path(inputs['label_policy']).read_text())
    policy['datasets']['swissprot_label_rich_v1']['GO']['explicit_negative_truth'] = SOURCE
    policy['locus_truth_complete'] = SOURCE
    (out / 'label_policy.json').write_text(json.dumps(policy, indent=2) + '\n')
    locus = out / 'locus_groups.tsv'
    locus.write_text('genome_id\tlocus_group_id\tcandidate_key\ttruth_complete\n' + ''.join(
        f'{GENOME}\tPMID15240874-threonylation\t{p}|GO|{TERM}|UNSPECIFIED\ttrue\n' for p, _ in MEMBERS.values()))
    result = P.build_table(Path(inputs['feature_dir']), Path(inputs['contract']), out / 'labels.tsv', Path(inputs['splits']),
        out / 'label_policy.json', out / 'table', registry_path=Path(inputs['registry']), go_obo=inputs['go_obo'],
        example_clusters=inputs['example_clusters'], locus_groups=locus, disjoint=tuple(old['split_audit']['disjoint_keys']))
    before = [json.loads(s) for s in Path(old['table']).read_text().splitlines()]
    after = [json.loads(s) for s in Path(result['table']).read_text().splitlines()]
    keys = ('example_id', 'features', 'label', 'split', 'group')
    if [{k:r[k] for k in keys} for r in before] != [{k:r[k] for k in keys} for r in after]:
        raise ValueError('candidate features, labels, identities or splits changed')
    contract = T.load_contract(Path(inputs['contract']))
    rows = T.read_table(Path(result['table']), contract)
    T.validate_locus_truth(rows)
    members = [r for r in rows if r.get('locus_group_id')]
    if result['locus_audit']['eligible_groups'] != 1 or len(members) != 2 or {r['split'] for r in members} != {'train'}:
        raise ValueError('reviewed pair did not survive complete-universe/family-mask checks')
    for p, digest in identity.items():
        if A.sha(p) != digest:
            raise ValueError('source changed: ' + p)
    report = {'schema': 'gspa.reviewed-thrs-locus/1', 'status': 'verified', 'source': SOURCE,
              'sequence_checks': sequence_checks, 'identity_sha256': identity, 'rows': len(rows),
              'original_features_labels_splits_unchanged': True, 'locus_audit': result['locus_audit'],
              'explicit_negative_rows': sum(r.get('label_basis') == 'explicit_negative_truth' for r in rows),
              'complete_groups_by_split': {'train': 1, 'calibration': 0, 'test': 0},
              'scope': 'One real training competition; no all-negative group, independent evaluation, null fit, or production model change.',
              'output_sha256': {str(p.relative_to(out)): A.sha(p) for p in out.rglob('*') if p.is_file()}}
    (out / 'verification.json').write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    return report


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('original-manifest', 'provenance', 'inventory', 'references', 'out'):
        p.add_argument('--' + name, type=Path, required=True)
    a = p.parse_args()
    r = build(a.original_manifest.resolve(), a.provenance.resolve(), a.inventory.resolve(), a.references.resolve(), a.out.resolve())
    print(json.dumps({k:r[k] for k in ('status', 'rows', 'explicit_negative_rows', 'complete_groups_by_split')}))
