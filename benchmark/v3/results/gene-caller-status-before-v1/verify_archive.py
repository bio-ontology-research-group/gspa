"""Independent, standard-library verification of retained pre-repair evidence."""
from pathlib import Path
import csv
import gzip
import hashlib
import io
import json

ROOT = Path(__file__).resolve().parent

def sha(data):
    return hashlib.sha256(data).hexdigest()

def verify():
    for line in (ROOT / 'archive.sha256').read_text().splitlines():
        digest, name = line.split('  ', 1)
        assert sha((ROOT / name).read_bytes()) == digest, name
    bindings = json.loads((ROOT / 'retained-bindings.json').read_text())
    def data(name):
        b = bindings[name]
        value = (ROOT / b['retained']).read_bytes()
        if b['retained'].endswith('.gz') and not name.endswith('.gz'):
            value = gzip.decompress(value)
        assert sha(value) == b['sha256'], name
        return value
    def text(name):
        return data(name).decode()
    def obj(name):
        return json.loads(text(name))
    for name in bindings:
        data(name)
    for name, digest in json.loads((ROOT / 'source-bindings.json').read_text()).items():
        assert sha((ROOT / 'source' / name).read_bytes()) == digest, name
    for job in ('10507', '10509', '10510', '10511'):
        assert text(f'jobs/{job}/exit-code.txt').strip() == '0'
    def records(name):
        output = []
        for block in text(name).split('>'):
            if not block:
                continue
            lines = block.splitlines()
            output.append((lines[0], ''.join(lines[1:])))
        return output
    def native(name):
        output = {}
        for header, sequence in records(name):
            parts = header.split(' # ')
            key = (parts[0].rsplit('_', 1)[0], parts[1], parts[2], '+' if parts[3] == '1' else '-')
            assert key not in output
            output[key] = sequence
        return output
    def tsv(name):
        return list(csv.DictReader(io.StringIO(text(name)), delimiter='\t'))
    def key(row):
        return tuple(row[n] for n in ('contig_id', 'start', 'end', 'strand'))
    outcomes = obj('jobs/10507/observed/observations.json')
    byname = {r['name']: r for r in outcomes['cases']}
    assert len(byname) == 7
    assert byname['short-meta']['exit_code'] == 0
    assert byname['short-single']['exit_code'] == 1
    assert byname['all-n-single']['exit_code'] == 1
    assert byname['all-n-meta']['exit_code'] == 0
    assert byname['empty-meta']['exit_code'] == 0
    for rec in byname.values():
        for ext in ('gff', 'faa'):
            content = data(f"jobs/10507/observed/{rec['name']}.{ext}")
            assert len(content) == rec[ext]['bytes'] and sha(content) == rec[ext]['sha256']
        rows = [line.split('\t') for line in text(f"jobs/10507/observed/{rec['name']}.gff").splitlines() if line and not line.startswith('#')]
        proteins = records(f"jobs/10507/observed/{rec['name']}.faa")
        assert len(rows) == rec['gff']['cds_rows'] == len(proteins) == rec['faa']['protein_headers']
    assert not records('jobs/10507/observed/short-meta.faa')
    assert text('jobs/10507/observed/short-meta.gff').startswith('##gff-version')
    nproteins = records('jobs/10507/observed/all-n-meta.faa')
    assert len(nproteins) == 1 and set(nproteins[0][1]) == {'X'}
    counts = {}
    summary = obj('jobs/10509/probe/summary.json')
    for mode, expected in (('single', 4359), ('meta', 4291)):
        reference = native(f'jobs/10507/observed/genome-{mode}.faa')
        actual = tsv(f'jobs/10509/probe/native-{mode}/genes.tsv')
        assert len(reference) == len(actual) == expected
        assert len({key(r) for r in actual}) == expected
        assert all(reference[key(r)] == r['sequence'] for r in actual)
        counts[mode] = expected
    for mode in ('single', 'meta'):
        result = obj(f'jobs/10509/probe/short-{mode}/result.json')
        assert result['protein_count'] == 0
        assert result['statuses'][0]['outcome'] == 'FAILED'
        assert result['statuses'][0]['exit_code'] is None
    reference = native('jobs/10507/observed/genome-single.faa')
    for name, expected in (('reordered', 4286), ('missing-protein', 190)):
        actual = tsv(f'jobs/10509/probe/{name}/genes.tsv')
        wrong = sum(reference.get(key(r)) != r['sequence'] for r in actual)
        result = obj(f'jobs/10509/probe/{name}/result.json')
        assert wrong == expected == summary[name]['sequence_coordinate_mismatches']
        assert result['statuses'][0]['outcome'] == 'SUCCESS'
        assert result['feature_count'] == 4359
        assert len(actual) == result['protein_count'] == (4358 if name == 'missing-protein' else 4359)
        counts[name + '_mismatches'] = wrong
    strict = tsv('jobs/10509/probe/reordered/strict-genes.tsv')
    assert len(strict) == len(reference) and all(reference[key(r)] == r['sequence'] for r in strict)
    source_sequences = records('jobs/10507/observed/inputs/genome.fna')
    numeric_sequences = records('jobs/10510/observed/inputs/numeric.fna')
    assert [h for h, _ in numeric_sequences] == ['2', '1']
    assert [s for _, s in numeric_sequences] == [s for _, s in source_sequences[:2]]
    rejected = obj('jobs/10511/probe/numeric-reader.json')
    for mode, expected in (('single', 294), ('meta', 291)):
        rows = [l.split('\t') for l in text(f'jobs/10510/observed/numeric-{mode}.gff').splitlines() if l and not l.startswith('#')]
        faa = records(f'jobs/10510/observed/numeric-{mode}.faa')
        assert len(rows) == len(faa) == expected
        aliases = {dict(item.split('=', 1) for item in r[8].split(';') if '=' in item)['ID']: r for r in rows}
        collisions = 0
        for header, sequence in faa:
            parts = header.split(' # ')
            numeric_id = dict(item.split('=', 1) for item in parts[4].split(';') if '=' in item)['ID']
            identifiers = {parts[0], numeric_id} & aliases.keys()
            candidates = [aliases[i] for i in identifiers]
            coord = (parts[0].rsplit('_', 1)[0], parts[1], parts[2], '+' if parts[3] == '1' else '-')
            matching = [r for r in candidates if (r[0], r[3], r[4], r[6]) == coord]
            assert len(matching) == 1
            collisions += len(candidates) > 1
        assert collisions > 0 and not rejected[mode]['accepted']
        assert rejected[mode]['error'] == 'missing or ambiguous CDS identity for 2_1'
        counts['numeric_' + mode + '_alias_collisions'] = collisions
    return dict(verified=True, retained_bindings=len(bindings), independently_checked=counts,
                limitation='Observed contract defects, no repair or biological accuracy claim.')

if __name__ == '__main__':
    print(json.dumps(verify(), indent=2))
