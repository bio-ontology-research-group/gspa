#!/usr/bin/env python3
"""One streaming pass over uniprot_sprot.dat.gz to extract, per accession,
the fields needed for the CAFA6 literature (BM25 text-kNN) and Net-KNN
(STRING PPI) components:

  accession \t taxon \t string_id \t name_text \t char_text

  * name_text  = DE (protein names) + GN (gene names) -- *identification*
                 fields, the temporally-safe text we allow for QUERY proteins
                 (a novel protein's name is usually assigned at entry creation
                 from homology/structure, not from the post-t0 functional study).
  * char_text  = CC -!- FUNCTION free text + KW keywords -- *characterization*
                 text; used only for the training CORPUS (whose labels are the
                 pre-t0 train_terms), never for query proteins, to avoid leaking
                 post-t0 function into the no-knowledge evaluation.

Every accession on the AC line (primary + secondary) is indexed to the entry,
because CAFA6 ground-truth IDs can be secondary accessions.
"""
from __future__ import annotations

import gzip
import re
import sys

DAT = sys.argv[1] if len(sys.argv) > 1 else \
    '/home/leechuck/Public/software/cafa6/uniprot_sprot.dat.gz'
OUT = sys.argv[2] if len(sys.argv) > 2 else 'text_string_index.tsv'

_TAX = re.compile(r'NCBI_TaxID=(\d+)')
_DE_NAME = re.compile(r'(?:RecName|SubName|AltName):\s*Full=([^;{]+)')
_GN_NAME = re.compile(r'Name=([^;{]+)')
_STRING = re.compile(r'^DR\s+STRING;\s+(\S+);')


def clean(s: str) -> str:
    return s.replace('\t', ' ').replace('\n', ' ').strip()


def flush(out, acs, taxon, string_id, names, chars):
    if not acs:
        return
    name_text = clean(' '.join(names))
    char_text = clean(' '.join(chars))
    for ac in acs:
        out.write(f'{ac}\t{taxon}\t{string_id}\t{name_text}\t{char_text}\n')


def main():
    n = 0
    with gzip.open(DAT, 'rt', errors='replace') as fh, open(OUT, 'w') as out:
        acs: list[str] = []
        taxon = ''
        string_id = ''
        names: list[str] = []
        chars: list[str] = []
        in_function = False
        for line in fh:
            tag = line[:2]
            if tag == 'AC':
                for a in line[5:].split(';'):
                    a = a.strip()
                    if a:
                        acs.append(a)
            elif tag == 'OX':
                m = _TAX.search(line)
                if m:
                    taxon = m.group(1)
            elif tag == 'DE':
                for m in _DE_NAME.finditer(line):
                    names.append(m.group(1).strip())
            elif tag == 'GN':
                for m in _GN_NAME.finditer(line):
                    names.append(m.group(1).strip())
            elif tag == 'KW':
                for kw in line[5:].split(';'):
                    kw = kw.strip().rstrip('.')
                    if kw:
                        chars.append(kw)
            elif tag == 'CC':
                body = line[5:].strip()
                if body.startswith('-!- FUNCTION:'):
                    in_function = True
                    chars.append(body[len('-!- FUNCTION:'):].strip())
                elif body.startswith('-!-'):
                    in_function = False
                elif in_function and not body.startswith('---') \
                        and 'Copyrighted by the UniProt' not in body:
                    chars.append(body)
            elif tag == 'DR':
                m = _STRING.match(line)
                if m and not string_id:
                    string_id = m.group(1)
            elif line.startswith('//'):
                flush(out, acs, taxon, string_id, names, chars)
                n += 1
                if n % 100000 == 0:
                    print(f'  {n:,} entries', file=sys.stderr, flush=True)
                acs, taxon, string_id, names, chars = [], '', '', [], []
                in_function = False
    print(f'done: {n:,} entries -> {OUT}', file=sys.stderr)


if __name__ == '__main__':
    main()
