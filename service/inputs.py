"""Validate annotated input and retain existing CDS identities and products."""
import re
from urllib.parse import unquote
from service import jobs


def gff3(folder):
    path = folder / 'sample.gff3'
    lines = path.read_text().splitlines()
    if not lines or not lines[0].startswith('##gff-version 3'):
        raise ValueError('Use an uncompressed GFF3 beginning with ##gff-version 3.')
    if '##FASTA' in lines:
        boundary = lines.index('##FASTA')
        if (folder / 'sample.fna').exists():
            raise ValueError('Supply either embedded FASTA or a companion genome FASTA, not both.')
        (folder / 'sample.fna').write_text('\n'.join(lines[boundary+1:])+'\n')
        lines = lines[:boundary]
    fasta = folder / 'sample.fna'
    if not fasta.exists():
        raise ValueError('This GFF3 has no embedded sequence. Add its matching genome FASTA in the companion field.')
    count, length = jobs.validate(fasta, 'genome')
    lengths = {pid: len(seq) for pid, seq in jobs.fasta(fasta)}
    cds = 0
    for line in lines:
        if not line or line.startswith('#'): continue
        f = line.split('\t')
        if len(f) != 9:
            raise ValueError('GFF3 feature rows must have nine tab-separated columns.')
        try:
            start, end = int(f[3]), int(f[4])
            if f[0] not in lengths or start < 1 or end < start or end > lengths[f[0]]: raise ValueError()
        except ValueError:
            raise ValueError('GFF3 coordinates must match the identifiers and lengths in the genome FASTA.')
        if f[2] == 'CDS':
            cds += 1
            attrs = dict(part.split('=', 1) for part in f[8].split(';') if '=' in part)
            if not any(attrs.get(k) for k in ('protein_id', 'ID', 'locus_tag')) or f[6] not in ('+', '-') or f[7] not in ('0','1','2'):
                raise ValueError('Each CDS needs a stable ID, strand, and phase (0, 1, or 2).')
            key = next((unquote(attrs[k]) for k in ('protein_id','ID','locus_tag') if attrs.get(k)), '')
            if not re.fullmatch(r'[A-Za-z0-9_.:|+-]{1,120}', key):
                raise ValueError('CDS protein identifiers must use letters, numbers, or _ . : | + - .')
            if attrs.get('transl_table', '11') not in ('1', '11'):
                raise ValueError('Currently the annotated-input workflow supports genetic codes 1 and 11.')
    if not cds or cds > jobs.MAX_SEQUENCES:
        raise ValueError('GFF3 must contain 1–20,000 CDS features.')
    path.write_text('\n'.join(lines)+'\n')
    return count, length
