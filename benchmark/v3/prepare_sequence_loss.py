#!/usr/bin/env python3
"""Label-blind nucleotide-block deletion controls for the acquired training cohort.

Keep complete original CDS translations only. Boundary-crossing CDSs are excluded;
this is controlled loss of sequence, not a simulation of gene calling on drafts.
Unplaced proteins are retained unchanged, outside the spatial perturbation.
"""
from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

from prepare_price_layout import fasta_sequences


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')


def write_fasta(path, records):
    with path.open('w') as fh:
        for key, seq in sorted(records.items()):
            fh.write(f'>{key}\n{seq}\n')


def prepare(cache: Path, out: Path, *, seed=42, block_bp=10000, drop_fraction=0.5):
    if isinstance(block_bp, bool) or not isinstance(block_bp, int) or block_bp < 1:
        raise ValueError('block_bp must be a positive integer')
    if not math.isfinite(drop_fraction) or not 0 < drop_fraction < 1:
        raise ValueError('drop_fraction must be strictly between zero and one')
    cache, out = cache.resolve(), out.resolve()
    # Immutable output: a second preparation needs a new directory, even with the same seed.
    out.mkdir(parents=True, exist_ok=False)
    for name in ('layout', 'claims', 'proteins', 'assemblies'):
        (out / name).mkdir()
    sources = [Path(__file__), cache / 'cohort_inputs/cohort_metadata.json',
               cache / 'cohort_inputs/cohort_orthogroups.tsv']
    original = json.loads(sources[1].read_text())
    missing, coordinates = copy.deepcopy(original), copy.deepcopy(original)
    report = {'schema': 'gspa.sequence-loss-inputs/1', 'seed': seed, 'block_bp': block_bp,
              'drop_fraction_blocks': drop_fraction, 'genomes': {},
              'selection': 'SHA256(seed, genome, contig, block start), smallest hashes removed; no labels read',
              'boundary_policy': 'exclude CDS crossing a block boundary, even between two retained blocks',
              'unplaced_policy': 'retain unchanged outside the spatial perturbation',
              'scope': 'whole-CDS sequence loss; not gene-calling or partial-protein predictor simulation'}
    retained = {}
    for gid, meta in sorted(original['genomes'].items()):
        dna_path = cache / 'genomes' / gid / 'genomic.fna'
        prot_path = cache / 'genomes' / gid / 'protein.faa'
        layout_path = cache / 'layout' / f'{gid}.layout.tsv'
        claim_path = cache / 'baseline/claims/homology' / f'{gid}.claims.jsonl'
        sources.extend([dna_path, prot_path, layout_path, claim_path])
        dna, proteins = fasta_sequences(dna_path), fasta_sequences(prot_path)
        with layout_path.open() as fh:
            reader = csv.DictReader(fh, delimiter='\t')
            header, rows = reader.fieldnames, list(reader)
        placed = {r['protein_id'] for r in rows}
        unplaced = set(meta.get('unplacedProteins', []))
        if len(placed) != len(rows) or placed & unplaced or placed | unplaced != set(proteins):
            raise ValueError(f'{gid}: protein FASTA must match unique layout plus declared unplaced proteins')
        if set(dna) != set(meta['contigs']):
            raise ValueError(f'{gid}: assembly contigs do not match metadata')
        blocks = []
        for contig, sequence in sorted(dna.items()):
            if len(sequence) != meta['contigs'][contig]['length']:
                raise ValueError(f'{gid}/{contig}: assembly length mismatch')
            for start in range(0, len(sequence), block_bp):
                end = min(start + block_bp, len(sequence))
                rank = hashlib.sha256(f'{seed}\t{gid}\t{contig}\t{start}'.encode()).hexdigest()
                blocks.append((rank, contig, start, end))
        n_drop = int(len(blocks) * drop_fraction)
        if not 0 < n_drop < len(blocks):
            raise ValueError(f'{gid}: block size/fraction does not yield both retained and removed sequence')
        deleted = {(c, s) for _, c, s, _ in sorted(blocks)[:n_drop]}
        fragments, lookup, contigs = {}, {}, {}
        for _, contig, start, end in blocks:
            if (contig, start) in deleted:
                continue
            cid = f'{contig}__block_{start + 1}_{end}'
            if cid in fragments:
                raise ValueError('fragment identifier collision')
            fragments[cid] = dna[contig][start:end]
            lookup[(contig, start)] = cid
            contigs[cid] = {'length': end - start, 'closed': False}
        keep, lost, new_rows = set(unplaced), {}, []
        for r in rows:
            pid, contig = r['protein_id'], r['contig']
            start, end = int(r['start']) - 1, int(r['end'])
            if contig not in dna or not 0 <= start < end <= len(dna[contig]):
                raise ValueError(f'{gid}/{pid}: invalid CDS coordinates')
            block_start = start // block_bp * block_bp
            if (end - 1) // block_bp != start // block_bp:
                lost[pid] = 'boundary_crossing_CDS'
            elif (contig, block_start) in deleted:
                lost[pid] = 'deleted_nucleotide_block'
            else:
                keep.add(pid)
                new_rows.append({**r, 'contig': lookup[(contig, block_start)],
                                 'start': start - block_start + 1, 'end': end - block_start})
        retained[gid] = keep
        with (out / 'layout' / layout_path.name).open('w') as fh:
            writer = csv.DictWriter(fh, fieldnames=header, delimiter='\t', lineterminator='\n')
            writer.writeheader()
            writer.writerows(sorted(new_rows, key=lambda r: (r['contig'], int(r['start']), r['protein_id'])))
        write_fasta(out / 'proteins' / f'{gid}.faa', {p: proteins[p] for p in keep})
        write_fasta(out / 'assemblies' / f'{gid}.fna', fragments)
        n_input = n_kept = 0
        with claim_path.open() as src, (out / 'claims' / claim_path.name).open('w') as dst:
            for line in src:
                if not line.strip():
                    continue
                c = json.loads(line)
                if c['protein_id'] not in proteins:
                    raise ValueError(f'{gid}: claim for undeclared protein {c["protein_id"]}')
                n_input += 1
                if c['protein_id'] in keep:
                    dst.write(line if line.endswith('\n') else line + '\n')
                    n_kept += 1
        dna_fraction = sum(map(len, fragments.values())) / sum(map(len, dna.values()))
        missing['genomes'][gid].update(assemblyStatus='DRAFT', contigs=contigs,
                                      completenessFraction=dna_fraction, unplacedProteins=sorted(unplaced))
        # Complete protein sequences, same surviving positions: isolate loss of spatial context.
        coordinates['genomes'][gid].update(assemblyStatus='DRAFT', contigs=contigs,
                                          completenessFraction=1.0,
                                          unplacedProteins=sorted(unplaced | set(lost)))
        report['genomes'][gid] = {
            'n_blocks': len(blocks), 'n_blocks_removed': n_drop,
            'retained_dna_fraction': dna_fraction, 'retained_protein_fraction': len(keep) / len(proteins),
            'n_original_proteins': len(proteins), 'n_retained_proteins': len(keep),
            'n_unplaced_outside_perturbation': len(unplaced),
            'claims_original': n_input, 'claims_retained': n_kept,
            'removed_proteins': lost, 'retained_proteins': sorted(keep)}
    og_src = cache / 'cohort_inputs/cohort_orthogroups.tsv'
    counts = defaultdict(lambda: {'original': 0, 'retained': 0})
    with og_src.open() as src, (out / 'orthogroups.tsv').open('w') as dst:
        reader = csv.DictReader(src, delimiter='\t')
        writer = csv.DictWriter(dst, fieldnames=reader.fieldnames, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for r in reader:
            gid, pid = r['genome_id'], r['protein_id']
            if gid not in retained or pid not in retained[gid] | set(report['genomes'][gid]['removed_proteins']):
                raise ValueError(f'orthogroup member not in original proteome: {gid}/{pid}')
            counts[gid]['original'] += 1
            if pid in retained[gid]:
                writer.writerow(r)
                counts[gid]['retained'] += 1
    for gid, value in counts.items():
        report['genomes'][gid]['orthogroup_memberships'] = value
    write_json(out / 'metadata_sequence_loss.json', missing)
    write_json(out / 'metadata_coordinate_loss.json', coordinates)
    report['inputs'] = {str(p): digest(p) for p in sources}
    report['outputs'] = {str(p.relative_to(out)): digest(p) for p in sorted(out.rglob('*')) if p.is_file()}
    write_json(out / 'preparation.json', report)
    return report


def verify(out):
    report = json.loads((out / 'preparation.json').read_text())
    for filename, expected in report['inputs'].items():
        if digest(filename) != expected:
            raise ValueError(f'changed sequence-loss source: {filename}')
    for filename, expected in report['outputs'].items():
        if digest(out / filename) != expected:
            raise ValueError(f'changed sequence-loss output: {filename}')
    return report


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--cache', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--seed', type=int, default=42)
    ap.add_argument('--block-bp', type=int, default=10000)
    args = ap.parse_args()
    prepare(args.cache, args.out, seed=args.seed, block_bp=args.block_bp)
