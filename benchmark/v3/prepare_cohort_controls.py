"""Freeze label-blind neighborhood-shuffle and exact-strain-duplicate controls."""
import copy
import csv
import hashlib
import json
from pathlib import Path
import random

import prepare_sequence_loss as P
from run_sequence_loss_benchmark import read_splits


def prepare(cache, out, splits_path, source_genome, *, seed=42, n_dups=5):
    cache, out = cache.resolve(), out.resolve()
    if type(seed) is not int or type(n_dups) is not int or n_dups < 1:
        raise ValueError('seed must be an integer and n_dups a positive integer')
    splits = read_splits(splits_path)
    meta_path = cache / 'cohort_inputs/cohort_metadata.json'
    og_path = cache / 'cohort_inputs/cohort_orthogroups.tsv'
    meta = json.loads(meta_path.read_text())
    if set(splits) != set(meta['genomes']):
        raise ValueError('split and metadata genome coverage differs')
    if splits.get(source_genome) != 'train':
        raise ValueError('duplicate source must be a training genome')
    source_meta = meta['genomes'][source_genome]
    if not source_meta.get('clade') or not source_meta.get('strainCluster'):
        raise ValueError('duplicate source needs known clade and strain cluster')
    duplicates = [f'{source_genome}_dup{i}' for i in range(1, n_dups + 1)]
    if set(duplicates) & set(splits) or any('_dup' in g for g in splits):
        raise ValueError('duplicate genome identity collision')
    sources = [Path(__file__).resolve(), Path(P.__file__).resolve(),
               Path(P.fasta_sequences.__code__.co_filename).resolve(),
               Path(read_splits.__code__.co_filename).resolve(), splits_path.resolve(), meta_path, og_path]
    for gid in sorted(splits):
        sources += [cache / 'layout' / f'{gid}.layout.tsv',
                    cache / 'genomes' / gid / 'protein.faa',
                    cache / 'baseline/claims/homology' / f'{gid}.claims.jsonl']
    identity = {str(p): P.digest(p) for p in sources}
    out.mkdir(parents=True, exist_ok=False)
    (out / 'layout_shuffled').mkdir()
    genomes, original_ids = {}, {}
    for gid, genome in sorted(meta['genomes'].items()):
        proteins = P.fasta_sequences(cache / 'genomes' / gid / 'protein.faa')
        with (cache / 'layout' / f'{gid}.layout.tsv').open() as fh:
            reader = csv.DictReader(fh, delimiter='\t')
            fields, rows = reader.fieldnames, list(reader)
        placed = {r['protein_id'] for r in rows}
        unplaced = set(genome.get('unplacedProteins', []))
        if len(placed) != len(rows) or placed & unplaced or placed | unplaced != set(proteins):
            raise ValueError('FASTA must equal unique placed plus declared unplaced proteins')
        original_ids[gid] = set(proteins)
        with (cache / 'baseline/claims/homology' / f'{gid}.claims.jsonl').open() as fh:
            for line in fh:
                if json.loads(line)['protein_id'] not in proteins:
                    raise ValueError('claim names an undeclared protein')
        by_contig = {}
        for row in rows:
            contig = row['contig']
            if contig not in genome['contigs'] or not 0 < int(row['start']) <= int(row['end']) <= genome['contigs'][contig]['length']:
                raise ValueError('invalid layout position')
            by_contig.setdefault(contig, []).append(row)
        shuffled, moved = [], 0
        for contig, slots in sorted(by_contig.items()):
            slots = sorted(slots, key=lambda r: (int(r['start']), int(r['end']), r['protein_id']))
            ids = [r['protein_id'] for r in slots]
            local_seed = hashlib.sha256(f'{seed}\t{gid}\t{contig}'.encode()).hexdigest()
            random.Random(local_seed).shuffle(ids)
            for pid, row in zip(ids, slots):
                shuffled.append({**row, 'protein_id': pid})
                moved += pid != row['protein_id']
        with (out / 'layout_shuffled' / f'{gid}.layout.tsv').open('w') as fh:
            writer = csv.DictWriter(fh, fieldnames=fields, delimiter='\t', lineterminator='\n')
            writer.writeheader(); writer.writerows(shuffled)
        genomes[gid] = {'retained_proteins': sorted(proteins), 'placed': len(rows),
                        'unplaced': len(unplaced), 'moved_by_shuffle': moved}
    with og_path.open() as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        fields, rows = reader.fieldnames, list(reader)
    for row in rows:
        if row['genome_id'] not in original_ids or row['protein_id'] not in original_ids[row['genome_id']]:
            raise ValueError('orthogroup names undeclared genome/protein')
    duplicate_meta = copy.deepcopy(meta)
    original_rows = [r for r in rows if r['genome_id'] == source_genome]
    for gid in duplicates:
        duplicate_meta['genomes'][gid] = copy.deepcopy(source_meta)
        rows.extend({**row, 'genome_id': gid} for row in original_rows)
    with (out / 'orthogroups_duplicated.tsv').open('w') as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)
    P.write_json(out / 'metadata_duplicated.json', duplicate_meta)
    if any(P.digest(path) != sha for path, sha in identity.items()):
        raise ValueError('control source changed during preparation')
    report = {'schema': 'gspa.cohort-controls-inputs/1', 'cache': str(cache), 'seed': seed,
              'source_genome': source_genome, 'duplicate_genomes': duplicates, 'genomes': genomes,
              'shuffle': 'permute placed protein IDs within each contig; retain positions/strand slots and all unplaced proteins',
              'duplicate': 'reuse exact source claims/layout; copy metadata and memberships within the same clade/strain cluster',
              'selection': 'explicit training genome; no function labels read', 'inputs': identity,
              'outputs': {str(p.relative_to(out)): P.digest(p) for p in out.rglob('*') if p.is_file()}}
    P.write_json(out / 'preparation.json', report)
    return report


def verify(out):
    report = json.loads((out / 'preparation.json').read_text())
    if report.get('schema') != 'gspa.cohort-controls-inputs/1':
        raise ValueError('invalid cohort control schema')
    for path, expected in report['inputs'].items():
        if P.digest(path) != expected:
            raise ValueError('changed cohort control source: ' + path)
    for path, expected in report['outputs'].items():
        if P.digest(out / path) != expected:
            raise ValueError('changed cohort control output: ' + path)
    return report
