"""Audit every verified Price matrix against exact assembly proteins and loci.

The input verification manifest must be trusted: its pickles execute Python on
loading. SHA-256 is checked again before deserialization. No functional labels are
read. Counts of retained candidates are planning diagnostics, not predictions.
"""
import argparse
import json
import re
from pathlib import Path

import numpy as np

try:
    from .prepare_price import load_matrix, norm_ec, FLOAT_RTOL
    from .prepare_price_layout import digest, fasta_sequences, read_coordinates
except ImportError:
    from prepare_price import load_matrix, norm_ec, FLOAT_RTOL
    from prepare_price_layout import digest, fasta_sequences, read_coordinates


def check_matrix(index, columns, values, protein_ids, locus_ids):
    if values.shape != (len(index), len(columns)) or len(set(index)) != len(index):
        raise ValueError('matrix shape or duplicate protein IDs')
    ecs = [norm_ec(c) for c in columns]
    if len(set(ecs)) != len(ecs):
        raise ValueError('duplicate normalized EC columns')
    if any(not re.fullmatch(r'\d+(?:\.(?:\d+|-)){0,3}', ec) for ec in ecs):
        raise ValueError('invalid EC identifier')
    complete = np.array([bool(re.fullmatch(r'\d+\.\d+\.\d+\.\d+', ec)) for ec in ecs])
    # Versionless IDs have an explicit, unique accession match, not sequence verification.
    versions = {}
    for pid in protein_ids:
        versions.setdefault(re.sub(r'\.\d+$', '', pid), []).append(pid)
    mapped, versionless, absent, ambiguous = {}, {}, [], []
    for pid in index:
        if pid in protein_ids:
            mapped[pid] = pid
        elif '.' not in pid and pid in versions:
            if len(versions[pid]) == 1:
                mapped[pid] = versionless[pid] = versions[pid][0]
            else:
                ambiguous.append(pid)
        else:
            absent.append(pid)
    canonical = list(mapped.values())
    if len(set(canonical)) != len(canonical):
        raise ValueError('multiple matrix rows resolve to the same acquired protein')
    total, threshold, top100, nonzero_proteins = 0, 0, 0, 0
    for row in values:
        row = np.asarray(row, dtype=float)
        if not np.isfinite(row).all() or (row < 0).any() or (row > 1 + FLOAT_RTOL).any():
            raise ValueError('nonfinite or out-of-range score')
        nz = (row != 0) & complete
        n = int(nz.sum())
        total += n
        nonzero_proteins += bool(n)
        above = int(((row >= .01) & nz).sum())
        threshold += above
        top100 += max(above, min(100, n))
    return {'shape': list(values.shape), 'complete_ec_columns': int(complete.sum()),
            'partial_ec_columns': int((~complete).sum()), 'nonzero_complete_ec_scores': total,
            'proteins_with_nonzero_complete_ec': nonzero_proteins,
            'nonzero_at_least_0_01': threshold, 'threshold_0_01_plus_top100_count': top100,
            'exact_protein_ids': len(mapped) - len(versionless), 'versionless_unique_matches': versionless,
            'absent_protein_ids': absent, 'ambiguous_versionless_ids': ambiguous,
            'matched_proteins_with_single_locus': sum(p in locus_ids for p in mapped.values()),
            'sequence_identity_to_original_predictor_input': 'unverified: original predictor FASTA not supplied'}


def audit(verification, genomes, out):
    if Path(out).exists():
        raise ValueError('output already exists')
    records = json.loads(Path(verification).read_text())['files']
    result = {'schema': 'gspa.price_matrix_audit/1', 'verification_sha256': digest(verification), 'files': []}
    for item in records:
        rec = {k: item[k] for k in ('predictor', 'assembly', 'destination')}
        try:
            if item.get('verified') is not True:
                raise ValueError('source transfer was not verified')
            if digest(item['destination']) != item['destination_sha256']:
                raise ValueError('matrix hash differs from verified transfer')
            rec['sha256'] = item['destination_sha256']
            folder = Path(genomes) / item['assembly']
            fasta, gff = folder / 'protein.faa', folder / 'genomic.gff'
            proteins = set(fasta_sequences(fasta))
            coords, _, _ = read_coordinates(gff, item['assembly'], proteins)
            index, columns, values = load_matrix(item['destination'], trust_pickle=True)
            rec.update(check_matrix(index, columns, values, proteins, {r[0] for r in coords}))
            del values
            rec['input_hashes'] = {p.name: digest(p) for p in (fasta, gff)}
            rec['status'] = 'unresolved_ids' if rec['absent_protein_ids'] or rec['ambiguous_versionless_ids'] else 'validated_id_join'
        except Exception as exc:
            rec.update(status='failed', error=f'{type(exc).__name__}: {exc}')
        result['files'].append(rec)
        print(json.dumps({'predictor': rec['predictor'], 'assembly': rec['assembly'], 'status': rec['status'],
                          'nonzero': rec.get('nonzero_complete_ec_scores'), 'error': rec.get('error')}), flush=True)
    result['all_ids_joined'] = all(r['status'] == 'validated_id_join' for r in result['files']) and bool(result['files'])
    Path(out).write_text(json.dumps(result, indent=2) + '\n')
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--verification', required=True)
    p.add_argument('--genomes', required=True)
    p.add_argument('--out', required=True)
    a = p.parse_args()
    result = audit(a.verification, a.genomes, a.out)
    raise SystemExit(0 if result['all_ids_joined'] else 2)
