"""Find exact full-length protein alignments in the registered assembly DNA.

This is an annotation-drift diagnostic, not a gene/function predictor. Inputs are
the matrix ID audit, explicitly versioned NCBI protein sequences and exact-version
assembly DNA. No partial match or alternate protein version supplies a coordinate.
Uniqueness here means one distinct exact alignment returned by TBLASTN, not an
exhaustive proof that no other genomic location could encode the protein.
"""
import argparse
import json
import subprocess
from pathlib import Path

try:
    from .prepare_price_layout import digest, fasta_sequences
except ImportError:
    from prepare_price_layout import digest, fasta_sequences


def exact_hits(path, sequences):
    hits = {pid: set() for pid in sequences}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) != 12 or f[0] not in sequences:
                raise ValueError('unexpected TBLASTN row or query ID')
            pid, qlen = f[0], len(sequences[f[0]])
            if int(f[4]) != qlen:
                raise ValueError('alignment query length differs from acquired sequence')
            start, end = sorted([int(f[7]), int(f[8])])
            frame = int(f[11])
            if (int(f[5]) == 1 and int(f[6]) == qlen and int(f[9]) == qlen and
                    int(f[3]) == qlen and int(f[10]) == 0 and end - start + 1 == 3 * qlen and
                    start > 0 and frame in (-3, -2, -1, 1, 2, 3)):
                hits[pid].add((f[1], start, end, '+' if frame > 0 else '-'))
    return {pid: [dict(zip(('contig', 'start', 'end', 'strand'), h)) for h in sorted(found)] for pid, found in hits.items()}


def recover(audit, proteins, dna_dir, out, executable='tblastn'):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=False)
    seqs = fasta_sequences(proteins)
    by = {}
    for r in json.loads(Path(audit).read_text())['files']:
        by.setdefault(r['assembly'], set()).update(pid for pid in r.get('absent_protein_ids', []) if '.' in pid)
    records = []
    version = subprocess.run([executable, '-version'], capture_output=True, text=True, check=True).stdout.strip()
    for asm, ids in sorted(by.items()):
        if not ids:
            continue
        if ids - seqs.keys():
            raise ValueError(f'{asm}: requested exact-version protein sequences missing: {sorted(ids - seqs.keys())}')
        folder = out / asm
        folder.mkdir()
        query, alignments = folder / 'queries.faa', folder / 'tblastn.tsv'
        query.write_text(''.join(f'>{pid}\n{seqs[pid]}\n' for pid in sorted(ids)))
        subject = Path(dna_dir) / asm / 'genomic.fna'
        cmd = [executable, '-query', str(query), '-subject', str(subject), '-out', str(alignments),
               '-outfmt', '6 qseqid sseqid pident length qlen qstart qend sstart send nident gaps sframe',
               '-evalue', '1e-10', '-seg', 'no', '-db_gencode', '11']
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        rec = {'assembly': asm, 'command': cmd, 'returncode': proc.returncode, 'stderr': proc.stderr,
               'input_sha256': {str(p): digest(p) for p in (query, subject)}}
        if proc.returncode == 0:
            rec.update(exact_full_length_hits=exact_hits(alignments, {pid: seqs[pid] for pid in ids}),
                       alignments_sha256=digest(alignments))
        records.append(rec)
        (folder / 'report.json').write_text(json.dumps(rec, indent=2) + '\n')
    result = {'schema': 'gspa.price_coordinate_recovery/1', 'tblastn_version': version,
              'input_sha256': {'audit': digest(audit), 'protein_fasta': digest(proteins), 'script': digest(__file__)},
              'genetic_code': 11, 'interpretation': __doc__, 'records': records}
    (out / 'report.json').write_text(json.dumps(result, indent=2) + '\n')
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('audit', 'proteins', 'dna-dir', 'out'):
        p.add_argument('--' + name, required=True)
    a = p.parse_args()
    r = recover(a.audit, a.proteins, a.dna_dir, a.out)
    print(json.dumps({'assemblies': len(r['records']), 'unique_exact_alignments': sum(len(h) == 1 for rec in r['records'] for h in rec.get('exact_full_length_hits', {}).values())}))
    raise SystemExit(any(rec['returncode'] for rec in r['records']))
