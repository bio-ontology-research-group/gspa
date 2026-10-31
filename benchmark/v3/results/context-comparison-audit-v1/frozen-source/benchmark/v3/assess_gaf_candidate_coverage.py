#!/usr/bin/env python3
"""Intersect an audited negative inventory with existing production candidates.

This is a source-review queue, never a training-label exporter. Matching GO IDs
are exact: no obsolete/alternate-ID normalization or negative propagation.
"""
import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path

import audit_gaf_negatives as A


def assess(inventory_path, feature_dir, out):
    inventory = json.loads(inventory_path.read_text())
    if inventory.get('schema') != 'gspa.gaf-negative-inventory/1' or inventory.get('status') != 'verified':
        raise ValueError('a verified complete inventory is required')
    candidate_path = inventory_path.parent / 'cohort-candidates.jsonl'
    if A.sha(candidate_path) != inventory['output_sha256']['cohort-candidates.jsonl']:
        raise ValueError('candidate inventory hash mismatch')
    identity = {str(p.resolve()): A.sha(p) for p in (inventory_path, candidate_path, Path(__file__), Path(A.__file__))}
    records = [json.loads(s) for s in candidate_path.read_text().splitlines() if s]
    requested = defaultdict(set)
    for r in records:
        for j in r['candidate_joins']:
            if j['status'] == 'needs_snapshot_and_reference_review':
                requested[j['genome_id']].add((j['protein_id'], r['go_id']))
    matches, groups = defaultdict(list), defaultdict(set)
    exports = sorted(feature_dir.glob('*/export_manifest.json'))
    if not exports:
        raise ValueError('no production feature exports found')
    total_rows = 0
    genomes = set()
    for manifest_path in exports:
        manifest = json.loads(manifest_path.read_text())
        if manifest.get('schema') != 'gspa.context-feature-export/1':
            raise ValueError('unsupported feature export: ' + str(manifest_path))
        rows_path = manifest_path.parent / manifest['rows_file']
        digest = A.sha(rows_path)
        if digest != manifest['rows_sha256']:
            raise ValueError('production feature hash mismatch: ' + str(rows_path))
        identity[str(manifest_path.resolve())] = A.sha(manifest_path)
        identity[str(rows_path.resolve())] = digest
        seen, n = set(), 0
        with rows_path.open() as fh:
            for line in fh:
                r = json.loads(line)
                n += 1
                genome = r['genome_id']
                if genome != manifest_path.parent.name:
                    raise ValueError('feature genome differs from export directory')
                genomes.add(genome)
                key = (genome, r['candidate_key'])
                if key in seen:
                    raise ValueError('duplicate production candidate')
                seen.add(key)
                if r['function_type'] != 'GO':
                    continue
                group = (genome, r['function_id'], r['role'])
                groups[group].add(r['protein_id'])
                if (r['protein_id'], r['function_id']) in requested[genome]:
                    matches[(genome, r['protein_id'], r['function_id'])].append({
                        'candidate_key': r['candidate_key'], 'role': r['role']})
        if n != manifest['rows']:
            raise ValueError('production feature count mismatch')
        total_rows += n
    queue = []
    status_counts, relation_counts = Counter(), Counter()
    for r in records:
        for j in r['candidate_joins']:
            key = (j['genome_id'], j['protein_id'], r['go_id'])
            found = matches.get(key, []) if j['status'] == 'needs_snapshot_and_reference_review' else []
            status = ('excluded_by_inventory' if j['status'] == 'excluded' else
                      'genome_not_exported' if j['genome_id'] not in genomes else
                      'exact_candidate_found_needs_review' if found else 'no_exact_candidate')
            status_counts[status] += 1
            if found:
                relation_counts[r['relation']] += 1
            queue.append({'gaf_line': r['line'], 'accession': r['accession'], 'go_id': r['go_id'],
                          'relation': r['relation'], 'evidence': r['evidence'], 'references': r['references'],
                          'genome_id': j['genome_id'], 'protein_id': j['protein_id'],
                          'status': status, 'exclusion_reasons': j['exclusion_reasons'],
                          'candidates': [dict(c, competing_proteins=len(groups[(key[0], key[2], c['role'])])) for c in found],
                          'accepted_training_label': False})
    for path, digest in identity.items():
        if A.sha(path) != digest:
            raise ValueError('coverage input changed: ' + path)
    out.mkdir(parents=True, exist_ok=False)
    queue_path = out / 'source-review-queue.jsonl'
    queue_path.write_text(''.join(json.dumps(r, sort_keys=True) + '\n' for r in queue))
    report = {'schema': 'gspa.gaf-negative-candidate-coverage/1', 'status': 'verified',
              'genomes': len(genomes), 'production_candidate_rows': total_rows,
              'cohort_gaf_records': len(records), 'candidate_join_counts': dict(status_counts),
              'matched_relation_counts': dict(relation_counts), 'accepted_training_labels': 0,
              'identity_sha256': identity, 'queue_sha256': A.sha(queue_path),
              'scope': 'Exact GO identifier intersection only; no normalization/propagation. Accession joins, relation/role compatibility and references still need review. No locus completeness or truth is certified.'}
    (out / 'coverage.json').write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    return report


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--inventory', type=Path, required=True)
    p.add_argument('--features', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    args = p.parse_args()
    result = assess(args.inventory.resolve(), args.features.resolve(), args.out.resolve())
    print(json.dumps({k: result[k] for k in ('status', 'production_candidate_rows', 'candidate_join_counts')}))
