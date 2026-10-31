"""Rehash the production receipts and inspect all required stage records."""
import json
import sys
from pathlib import Path
from check_inputs import sha


def check(work):
    result = {}
    for arm in ('intact', 'shuffled_neighborhoods', 'duplicated_strain_5x'):
        path = work / arm / 'out/verified_run.json'
        receipt = json.loads(path.read_text())
        assert receipt['returncode'] == 0
        spec = receipt['spec']
        for source, expected in spec['inputs'].items():
            assert sha(source) == expected, source
        for source, expected in receipt['artifacts'].items():
            assert sha(source) == expected, source
        coverage = {}
        for gid, output in spec['outputs'].items():
            ledger = json.loads(Path(output + '.stages.json').read_text())
            assert ledger['genome_id'] == gid and ledger['schema_version'] == 2
            states = {s['name']: s['status'] for s in ledger['stages']}
            assert len(states) == len(ledger['stages'])
            assert all(states[s] == 'executed' for s in spec['required_stages'])
            assert not set(states.values()) & {'failed', 'not_run'}
            coverage[gid] = states
        result[arm] = dict(receipt_sha256=sha(path), input_files=len(spec['inputs']),
                           output_files=len(receipt['artifacts']), genomes=len(coverage), stages=coverage)
    return result


if __name__ == '__main__':
    print(json.dumps(check(Path(sys.argv[1])), indent=2, sort_keys=True))
