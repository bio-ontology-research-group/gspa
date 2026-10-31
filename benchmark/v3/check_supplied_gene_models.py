#!/usr/bin/env python3
"""Check supplied native gene models through the actual CLI; not an accuracy benchmark."""
import argparse
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import time


def sha(path):
    with Path(path).open('rb') as stream:
        digest = hashlib.sha256()
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def fasta(path):
    opener = gzip.open if path.suffix == '.gz' else open
    records = {}
    with opener(path, 'rt') as stream:
        for line in stream:
            if line.startswith('>'):
                header = line[1:].strip()
                pid = header.split()[0]
                assert pid not in records
                records[pid] = [header, '']
            elif line.strip():
                records[pid][1] += line.strip()
    return records


def native_reference(gff, faa):
    records = fasta(faa)
    ref = {}
    for pid, (header, sequence) in records.items():
        fields = header.split(' # ')
        ref[pid] = dict(contig=pid.rsplit('_', 1)[0], start=int(fields[1]), end=int(fields[2]),
                        strand='+' if fields[3] == '1' else '-', sequence=sequence)
    cds = [line.split('\t') for line in gff.read_text().splitlines() if line and not line.startswith('#')]
    coords = [(r[0], int(r[3]), int(r[4]), r[6]) for r in cds if r[2] == 'CDS']
    assert len(coords) == len(set(coords)) == len(ref)
    assert set(coords) == {(r['contig'], r['start'], r['end'], r['strand']) for r in ref.values()}
    return ref


def expected_groups(reference, excluded=()):
    contigs = {}
    for pid, r in reference.items():
        if pid not in excluded:
            contigs.setdefault(r['contig'], []).append((pid, r))
    groups = set()
    for contig, proteins in contigs.items():
        proteins.sort(key=lambda item: item[1]['start'])
        current = []
        for pid, r in proteins:
            if current:
                previous = reference[current[-1]]
                if previous['strand'] != r['strand'] or r['start'] - previous['end'] - 1 > 300:
                    if len(current) >= 2:
                        groups.add(tuple(current))
                    current = []
            current.append(pid)
        if len(current) >= 2:
            groups.add(tuple(current))
    return groups


def run(args):
    root = args.root.resolve()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    (out / 'tmp').mkdir()
    bindings = {}
    runs, results = [], {}

    def bind(path):
        bindings[str(path)] = sha(path)

    for p in (args.jar.resolve(), args.java.resolve(), Path(__file__).resolve()):
        bind(p)
    cases = []
    for caller in ('pyrodigal', 'prodigal'):
        for mode in ('single', 'meta'):
            folder = root / caller
            cases.append(dict(name=caller + '-' + mode, gff=folder / ('genome-' + mode + '.gff'),
                              faa=folder / ('genome-' + mode + '.faa'), genome=folder / 'inputs/genome.fna'))
    for mode in ('single', 'meta'):
        folder = root / 'numeric'
        cases.append(dict(name='numeric-' + mode, gff=folder / ('numeric-' + mode + '.gff'),
                          faa=folder / ('numeric-' + mode + '.faa'), genome=folder / 'inputs/numeric.fna'))
    prod = cases[2]
    gz = out / 'prodigal-single.faa.gz'
    with gzip.open(gz, 'wb') as stream:
        stream.write(prod['faa'].read_bytes())
    records = fasta(prod['faa'])
    reordered = out / 'reordered.faa'
    reordered.write_text(''.join('>' + header + '\n' + seq + '\n' for header, seq in reversed(list(records.values()))))
    changed = out / 'conflicting-coordinate.faa'
    altered = list(records)
    first = altered[0]
    header, seq = records[first]
    fields = header.split(' # ')
    fields[1] = str(int(fields[1]) + 1)
    changed.write_text('>' + ' # '.join(fields) + '\n' + seq + '\n' + ''.join(
        '>' + h + '\n' + s + '\n' for h, s in list(records.values())[1:]))
    cases.extend([
        dict(prod, name='prodigal-gzip', faa=gz),
        dict(prod, name='prodigal-without-genome', genome=None),
        dict(prod, name='prodigal-reordered', faa=reordered),
        dict(prod, name='prodigal-header-conflict', faa=changed, excluded=[first]),
    ])
    bind(prod['faa'])
    for case in cases:
        name, gff, faa, genome = (case[k] for k in ('name', 'gff', 'faa', 'genome'))
        for p in (gff, faa, genome):
            if p:
                bind(p)
        excluded = case.get('excluded', [])
        reference = native_reference(gff, prod['faa'] if excluded else faa)
        config = out / (name + '.config.json')
        config.write_text(json.dumps(dict(input=dict(kingdom='bacteria', type='genome'),
            predictors={k: dict(enabled=k == 'operons') for k in ('similarity', 'domains', 'pathway', 'operons', 'localization')},
            integration=dict(enabled=False), quality=dict(enabled=False))) + '\n')
        bind(config)
        command = [str(args.java), '-XX:ActiveProcessorCount=3', '-jar', str(args.jar), 'annotate',
            '--gff3', str(gff), '--proteins', str(faa), '--config', str(config), '--output', str(out / name),
            '--threads', '1', '--no-metrics']
        if genome:
            command += ['--genome', str(genome)]
        started = time.monotonic()
        timed_out = False
        with (out / (name + '.log')).open('w') as log:
            proc = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT, start_new_session=True,
                env=dict(os.environ, TMPDIR=str(out / 'tmp'), JAVA_TOOL_OPTIONS='-Djava.io.tmpdir=' + str(out / 'tmp')))
            try:
                proc.wait(timeout=300)
            except subprocess.TimeoutExpired:
                timed_out = True
                try:
                    os.killpg(proc.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                proc.wait()
        runs.append(dict(name=name, command=command, exit_code=proc.returncode,
                         timed_out=timed_out, elapsed_seconds=time.monotonic() - started))
        (out / 'runs.json').write_text(json.dumps(runs, indent=2) + '\n')
        assert not timed_out and proc.returncode == 0, name
        folder = out / name
        with (folder / (gff.stem + '_proteins.tsv')).open() as stream:
            rows = list(csv.DictReader(stream, delimiter='\t'))
        assert len(rows) == len(reference) and {r['protein_id'] for r in rows} == set(reference), name
        emitted = fasta(folder / (gff.stem + '.faa'))
        assert {pid: seq for pid, (_, seq) in emitted.items()} == {pid: r['sequence'] for pid, r in reference.items()}, name
        for r in rows:
            ref = reference[r['protein_id']]
            if r['protein_id'] in excluded:
                assert int(r['start']) == int(r['end']) == 0, r
            else:
                assert (r['contig'], int(r['start']), int(r['end']), r['strand']) == (
                    ref['contig'], ref['start'], ref['end'], ref['strand']), r
        expected = expected_groups(reference, excluded)
        lines = (folder / 'operons_for_integrate.tsv').read_text().splitlines()
        observed = [tuple(line.split('\t')) for line in lines]
        assert len(observed) == len(set(observed)) and set(observed) == expected, name
        with (folder / 'protein_to_operon.tsv').open() as stream:
            members = list(csv.DictReader(stream, delimiter='\t'))
        assert len(members) == sum(map(len, expected))
        assert {r['protein_id'] for r in members} == {p for group in expected for p in group}
        # Scoped baseline parity: Pyrodigal already joined correctly before this repair.
        parity = []
        if name.startswith('pyrodigal-'):
            old = root / 'before' / name
            for suffix in ('.faa', '.gff3', '_proteins.tsv', '_annotations.tsv'):
                assert (folder / (gff.stem + suffix)).read_bytes() == (old / (gff.stem + suffix)).read_bytes(), (name, suffix)
                bind(old / (gff.stem + suffix))
                parity.append(suffix)
            for filename in ('operons.tsv', 'operons_for_integrate.tsv', 'protein_to_operon.tsv'):
                assert (folder / filename).read_bytes() == (old / filename).read_bytes(), (name, filename)
                bind(old / filename)
                parity.append(filename)
        results[name] = dict(proteins=len(rows), placed=len(rows)-len(excluded), unplaced=len(excluded),
                            operons=len(expected), members=sum(map(len, expected)), byte_identical_baseline=parity)
        (out / 'summary.json').write_text(json.dumps(results, indent=2) + '\n')
    assert all(sha(Path(p)) == digest for p, digest in bindings.items())
    (out / 'bindings.json').write_text(json.dumps(bindings, indent=2) + '\n')
    print(json.dumps(results, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('root', 'out', 'jar', 'java'):
        parser.add_argument('--' + key, type=Path, required=True)
    run(parser.parse_args())
