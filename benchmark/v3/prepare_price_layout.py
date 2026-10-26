"""Build label-free Price genome layouts from exact-version NCBI packages.

Only protein identifiers, coordinates, strand and contig bounds are extracted.
Functional GFF annotations never enter the output. A protein with several loci,
split CDS coordinates or an unknown strand has no invented single locus; it stays
in the sequence input and is recorded as unavailable for local context.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path
from urllib.parse import unquote


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for block in iter(lambda: fh.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def fasta_sequences(path):
    sequences, pid, parts = {}, None, []
    def finish():
        if pid is None:
            return
        if pid in sequences:
            raise ValueError(f'duplicate FASTA protein ID: {pid}')
        sequence = ''.join(parts).upper().removesuffix('*')
        if not sequence:
            raise ValueError(f'empty sequence: {pid}')
        sequences[pid] = sequence
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                finish()
                pid, parts = line[1:].split()[0], []
            elif pid is None:
                raise ValueError('sequence before FASTA header')
            else:
                parts.append(line)
    finish()
    if not sequences:
        raise ValueError('empty FASTA')
    return sequences


def attributes(text):
    return {k: unquote(v) for k, v in (x.split('=', 1) for x in text.split(';') if '=' in x)}


def read_coordinates(path, expected_assembly, protein_ids):
    loci, contigs, declared = defaultdict(set), {}, None
    with open(path) as fh:
        for line_no, line in enumerate(fh, 1):
            if line.startswith('##FASTA'):
                break
            if line.startswith('#!genome-build-accession '):
                declared = line.strip().split(':')[-1]
            elif line.startswith('##sequence-region '):
                _, contig, lo, hi = line.split()
                if int(lo) != 1 or int(hi) < 1:
                    raise ValueError(f'unsupported sequence region: {line.strip()}')
                contigs[contig] = {'length': int(hi)}
            elif line.startswith('#') or not line.strip():
                continue
            else:
                f = line.rstrip('\n').split('\t')
                if len(f) != 9:
                    raise ValueError(f'{path}:{line_no}: expected 9 GFF fields')
                if f[2] != 'CDS':
                    continue
                a = attributes(f[8])
                pid = a.get('protein_id')
                if pid not in protein_ids:
                    continue
                start, end = int(f[3]), int(f[4])
                if start < 1 or end < start:
                    raise ValueError(f'{path}:{line_no}: invalid coordinates')
                loci[pid].add((f[0], start, end, f[6]))
    if declared != expected_assembly:
        raise ValueError(f'GFF assembly {declared!r} != requested {expected_assembly!r}')
    rows, unavailable = [], {}
    for pid in sorted(protein_ids):
        found = loci[pid]
        if len(found) != 1:
            unavailable[pid] = 'missing_CDS' if not found else 'multiple_or_split_CDS_loci'
            continue
        contig, start, end, strand = next(iter(found))
        if contig not in contigs or end > contigs[contig]['length']:
            unavailable[pid] = 'coordinates_outside_declared_contig'
        elif strand not in ('+', '-'):
            unavailable[pid] = 'unknown_strand'
        else:
            rows.append((pid, contig, start, end, strand))
    return rows, unavailable, contigs


def prepare(genomes_dir, outdir):
    genomes_dir, outdir = Path(genomes_dir), Path(outdir)
    outdir.mkdir(parents=True, exist_ok=False)
    records, metadata = {}, {'genomes': {}}
    with (outdir / 'all_proteins.faa').open('w') as all_fasta:
        for folder in sorted(genomes_dir.glob('GC[AF]_*')):
            if not folder.is_dir():
                continue
            asm = folder.name
            fasta, gff = folder / 'protein.faa', folder / 'genomic.gff'
            seqs = fasta_sequences(fasta)
            rows, missing, contigs = read_coordinates(gff, asm, set(seqs))
            report_path = folder / 'assembly_data_report.json'
            report = json.loads(report_path.read_text())
            if report.get('accession') != asm:
                raise ValueError(f'assembly report accession does not match {asm}')
            complete = report.get('assemblyInfo', {}).get('assemblyLevel') == 'Complete Genome'
            # Complete replicons may be linear. Closure denotes assembled sequence,
            # not circularity; no family-search or phylogenetic independence is inferred.
            if complete:
                for c in contigs.values():
                    c['closed'] = True
            metadata['genomes'][asm] = {'assemblyStatus': 'COMPLETE' if complete else 'DRAFT', 'contigs': contigs}
            layout = outdir / f'{asm}.layout.tsv'
            with layout.open('w') as fh:
                writer = csv.writer(fh, delimiter='\t', lineterminator='\n')
                writer.writerow(['protein_id', 'contig', 'start', 'end', 'strand'])
                writer.writerows(rows)
            seqfile = outdir / f'{asm}.sequence_hashes.tsv'
            with seqfile.open('w') as fh:
                fh.write('protein_id\tlength\tsequence_sha256\n')
                for pid, seq in sorted(seqs.items()):
                    fh.write(f'{pid}\t{len(seq)}\t{hashlib.sha256(seq.encode()).hexdigest()}\n')
                    all_fasta.write(f'>{asm}__{pid}\n{seq}\n')
            records[asm] = {'proteins': len(seqs), 'single_locus_proteins': len(rows),
                            'unavailable_loci': missing, 'taxid': report.get('organism', {}).get('taxId'),
                            'input_hashes': {p.name: digest(p) for p in (fasta, gff, report_path)},
                            'layout_sha256': digest(layout), 'sequence_hashes_sha256': digest(seqfile)}
    if not records:
        raise ValueError('no GCA/GCF assembly directories')
    (outdir / 'cohort_metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')
    result = {'schema': 'gspa.price_layout/1', 'assemblies': records,
              'all_proteins_sha256': digest(outdir / 'all_proteins.faa'),
              'interpretation': 'Coordinates only. No functional labels, orthology, clade independence or family-search completeness inferred. Sequence hashes describe acquired NCBI proteins; equality to original predictor FASTA remains unverified.'}
    (outdir / 'report.json').write_text(json.dumps(result, indent=2) + '\n')
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--genomes', required=True)
    p.add_argument('--out', required=True)
    a = p.parse_args()
    result = prepare(a.genomes, a.out)
    print(json.dumps({'assemblies': len(result['assemblies']), 'proteins': sum(r['proteins'] for r in result['assemblies'].values()),
                      'single_locus_proteins': sum(r['single_locus_proteins'] for r in result['assemblies'].values())}))
