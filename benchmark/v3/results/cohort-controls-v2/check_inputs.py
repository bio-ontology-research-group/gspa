"""Independent input audit; does not import the preparation implementation."""
import csv
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def check(folder):
    preparation = json.loads((folder / 'preparation.json').read_text())
    cache = Path(preparation['cache'])
    for path, expected in preparation['inputs'].items():
        assert sha(path) == expected, path
    for path, expected in preparation['outputs'].items():
        assert sha(folder / path) == expected, path
    original = json.loads((cache / 'cohort_inputs/cohort_metadata.json').read_text())['genomes']
    duplicated = json.loads((folder / 'metadata_duplicated.json').read_text())['genomes']
    source = preparation['source_genome']
    copies = preparation['duplicate_genomes']
    assert set(duplicated) == set(original) | set(copies)
    assert all(duplicated[g] == original[g] for g in original)
    assert all(duplicated[g] == original[source] for g in copies)
    total = moved = nproteins = unplaced = 0
    for gid, metadata in original.items():
        a = rows(cache / 'layout' / f'{gid}.layout.tsv')
        b = rows(folder / 'layout_shuffled' / f'{gid}.layout.tsv')
        slots = lambda data: Counter(tuple(sorted((k, v) for k, v in row.items() if k != 'protein_id')) for row in data)
        membership = lambda data: Counter((r['contig'], r['protein_id']) for r in data)
        assert slots(a) == slots(b), gid
        assert membership(a) == membership(b), gid
        fasta_ids = [line[1:].split()[0] for line in (cache / 'genomes' / gid / 'protein.faa').read_text().splitlines() if line.startswith('>')]
        declared = [r['protein_id'] for r in a] + metadata.get('unplacedProteins', [])
        assert Counter(fasta_ids) == Counter(declared)
        assert len(set(fasta_ids)) == len(fasta_ids)
        assert set(fasta_ids) == set(preparation['genomes'][gid]['retained_proteins'])
        a_by_id = {r['protein_id']: r for r in a}
        moved += sum(r != a_by_id[r['protein_id']] for r in b)
        total += len(a)
        nproteins += len(fasta_ids)
        unplaced += len(metadata.get('unplacedProteins', []))
    a = rows(cache / 'cohort_inputs/cohort_orthogroups.tsv')
    b = rows(folder / 'orthogroups_duplicated.tsv')
    canonical = lambda data: Counter(tuple(sorted(r.items())) for r in data)
    assert canonical(a) == canonical([r for r in b if r['genome_id'] in original])
    source_rows = [r for r in a if r['genome_id'] == source]
    for copy in copies:
        assert canonical(source_rows) == canonical([{**r, 'genome_id': source} for r in b if r['genome_id'] == copy])
    return dict(status='verified', genomes=len(original), original_proteins=nproteins,
                unchanged_contig_slots_and_memberships=total, unplaced=unplaced,
                proteins_moved=moved, duplicates=len(copies), source_gene_memberships=len(source_rows),
                inputs_rehashed=len(preparation['inputs']), outputs_rehashed=len(preparation['outputs']))


if __name__ == '__main__':
    print(json.dumps(check(Path(sys.argv[1])), indent=2))
