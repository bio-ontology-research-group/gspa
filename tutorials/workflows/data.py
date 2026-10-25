#!/usr/bin/env python3
"""Small public teaching data, with checksummed downloads and explicit derivations."""
import csv
import hashlib
import io
import json
from pathlib import Path
import re
import sys
import time
import urllib.parse
import urllib.request

ROOT = Path('/work') if Path('/work').is_dir() else Path('.tutorial-work')
DATA = ROOT / 'data'
DATA.mkdir(parents=True, exist_ok=True)
LOCK_PATH = Path(__file__).with_name('data-lock.json')
LOCK = json.loads(LOCK_PATH.read_text()) if LOCK_PATH.exists() else {}
RECEIPTS = {}


def fetch(name, url):
    dest = DATA / name
    expected = LOCK.get(name, {}).get('sha256')
    if not dest.exists():
        for attempt in range(3):
            try:
                with urllib.request.urlopen(url, timeout=180) as response:
                    body = response.read()
                assert body and not body.lstrip().lower().startswith(b'<!doctype html'), url
                dest.with_suffix(dest.suffix + '.part').write_bytes(body)
                dest.with_suffix(dest.suffix + '.part').replace(dest)
                break
            except Exception:
                if attempt == 2:
                    raise
                time.sleep(2)
    digest = hashlib.sha256(dest.read_bytes()).hexdigest()
    if expected and expected != digest:
        raise SystemExit(f'{name}: source changed. Expected {expected}, got {digest}. Preserve the recorded snapshot or deliberately refresh the tutorial lock and outputs.')
    RECEIPTS[name] = {'url': url, 'sha256': digest, 'bytes': dest.stat().st_size}
    return dest


def fasta(path):
    name, chunks = None, []
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):
            if name is not None:
                yield name, ''.join(chunks)
            name, chunks = line[1:].split()[0], []
        else:
            chunks.append(line.strip())
    if name is not None:
        yield name, ''.join(chunks)


def write_fasta(path, entries):
    with Path(path).open('w') as out:
        for pid, sequence in entries:
            out.write(f'>{pid}\n')
            for i in range(0, len(sequence), 70):
                out.write(sequence[i:i+70] + '\n')


def main():
    for name, accession in [('bacteria', 'U00096.3'), ('archaea', 'AE009950.1'), ('virus', 'J02459.1')]:
        path = fetch(f'{accession}.fasta', f'https://www.ebi.ac.uk/ena/browser/api/fasta/{accession}')
        entries = list(fasta(path))
        assert len(entries) == 1 and len(entries[0][1]) > 40000
        write_fasta(DATA / f'{name}.fna', [(accession, entries[0][1])])
    # Organism names are resolved by UniProt's taxonomy search, including descendants.
    query = 'reviewed:true AND (organism_id:83333 OR organism_id:186497 OR organism_id:559292 OR taxonomy_name:"Lambdavirus lambda")'
    params = urllib.parse.urlencode({'query': query, 'format': 'tsv',
        'fields': 'accession,sequence,go_id,ec,protein_name,organism_id', 'sort': 'accession asc'})
    table = fetch('uniprot-reviewed.tsv', 'https://rest.uniprot.org/uniprotkb/stream?' + params)
    rows = list(csv.DictReader(table.open(), delimiter='\t'))
    assert len(rows) > 5000, 'Expected whole reviewed E. coli/yeast proteomes plus archaeal and phage references'
    references, yeast = [], []
    for row in rows:
        terms = re.findall(r'GO:\d{7}', row['Gene Ontology IDs'])
        ecs = ['EC:' + e.strip() for e in row['EC number'].split(';') if e.strip()]
        references.append((row['Entry'] + ' ' + ' '.join(terms + ecs), row['Sequence']))
        if row['Organism (ID)'] == '559292':
            yeast.append((row['Entry'], row['Sequence']))
    write_fasta(DATA / 'reference.faa', references)
    write_fasta(DATA / 'eukaryote.faa', yeast)
    # Deterministic teaching bin: retain 60% of one archaeal genome in 20 kb chunks.
    archaeal = next(fasta(DATA / 'archaea.fna'))[1]
    pieces = [(f'archaeal_bin_{i//20000+1:03d}', archaeal[i:i+20000])
              for i in range(0, len(archaeal), 20000) if (i//20000) % 5 < 3]
    write_fasta(DATA / 'mag.fna', pieces)
    # Constructed unbinned community, keeping sequences separate (never concatenate genomes).
    mixed = []
    for source in ['bacteria', 'archaea', 'virus']:
        for pid, seq in fasta(DATA / f'{source}.fna'):
            mixed.extend((f'{source}_{i//20000+1:03d}', seq[i:i+20000]) for i in range(0, len(seq), 20000))
    write_fasta(DATA / 'metagenome.fna', mixed)
    config = {'predictors': {'domains': {'enabled': False}, 'neural': {'basePredictor': 'none'},
                            'similarity': {'maxTargetSeqs': 1}},
              'integration': {'enabled': False}}
    (DATA / 'cpu.json').write_text(json.dumps(config, indent=2) + '\n')
    virus = json.loads(json.dumps(config))
    virus['predictors']['operons'] = {'enabled': False}
    virus['predictors']['pathway'] = {'enabled': False}
    virus['predictors']['localization'] = {'enabled': False}
    virus['quality'] = {'enabled': False, 'completeness': {'profile': 'skip'},
                        'consistency': {'enforce': False}, 'coherence': {'enforce': False}}
    (DATA / 'virus.json').write_text(json.dumps(virus, indent=2) + '\n')
    (DATA / 'downloads.json').write_text(json.dumps(RECEIPTS, indent=2) + '\n')
    for name, suffix in [('bacteria','fna'),('archaea','fna'),('eukaryote','faa'),('mag','fna'),('metagenome','fna'),('virus','fna')]:
        seqs = list(fasta(DATA / f'{name}.{suffix}'))
        print(f'{name}: {len(seqs):,} sequences; {sum(len(s) for _,s in seqs):,} {"amino acids" if suffix == "faa" else "bp"}')
    print(f'DIAMOND teaching reference: {len(references):,} reviewed proteins with GO/EC headers')
    print('MAG and community inputs are constructed teaching datasets, not environmental samples.')
    print('Self/close matches are present: these runs demonstrate operation, not prediction accuracy.')

if __name__ == '__main__':
    main()
