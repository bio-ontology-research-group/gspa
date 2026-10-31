"""Exact candidate identity comparison with bounded-memory external sorting.

All keys are compared, including duplicate detection; hashes or samples are not
substitutes for equality. Scratch belongs on the caller's local job disk.
"""
import csv
import itertools
import json
import os
from pathlib import Path
import shutil
import subprocess


def sort_executable():
    executable = shutil.which('sort')
    if executable is None:
        raise RuntimeError('GNU sort is required for bounded-memory verification')
    return Path(executable).resolve()


def _key(protein, kind, function):
    if not all(isinstance(x, str) and x for x in (protein, kind, function)):
        raise ValueError('candidate identity must contain three nonempty strings')
    # The production TSV adds this prefix for ECs. Other identifiers are exact.
    if kind == 'EC':
        function = function.removeprefix('EC:')
    return json.dumps([protein, kind, function], ensure_ascii=True, separators=(',', ':')) + '\n'


def _sort_keys(keys, destination):
    raw = destination.with_suffix('.unsorted')
    with raw.open('w', encoding='ascii') as stream:
        for key in keys:
            stream.write(key)
    subprocess.run([str(sort_executable()), '--buffer-size=64M', '--parallel=1',
                    '--temporary-directory', str(destination.parent),
                    '--output', str(destination), str(raw)],
                   env={**os.environ, 'LC_ALL': 'C'}, check=True)
    raw.unlink()


def _unique_lines(path):
    previous = None
    with path.open(encoding='ascii') as stream:
        for line in stream:
            if line == previous:
                raise ValueError(f'duplicate candidate in {path.name}: {line.strip()}')
            previous = line
            yield line


def expected_candidates(claims, scratch):
    """Write the sorted, unique complete input universe; return its path/count."""
    def keys():
        with Path(claims).open() as stream:
            for line in stream:
                row = json.loads(line)
                yield _key(row['protein_id'], row['function_type'], row['function_id'])
    result = Path(scratch) / 'expected.keys'
    _sort_keys(keys(), result)
    return result, sum(1 for _ in _unique_lines(result))


def compare_candidates(expected, prediction, scratch):
    """Accept arbitrary output order, but reject every missing/extra/duplicate key."""
    def keys():
        with Path(prediction).open() as stream:
            reader = csv.DictReader(stream, delimiter='\t')
            required = {'protein_id', 'type', 'function_id'}
            if not required <= set(reader.fieldnames or []):
                raise ValueError('prediction header lacks candidate identity')
            for row in reader:
                yield _key(row['protein_id'], row['type'], row['function_id'])
    actual = Path(scratch) / 'actual.keys'
    _sort_keys(keys(), actual)
    count = 0
    try:
        for left, right in itertools.zip_longest(_unique_lines(expected), _unique_lines(actual)):
            if left != right:
                raise ValueError('production output candidate universe differs from full matrix input')
            count += 1
        return count
    finally:
        actual.unlink(missing_ok=True)
