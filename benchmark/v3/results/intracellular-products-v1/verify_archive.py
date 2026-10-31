#!/usr/bin/env python3
"""Offline integrity and retained-result checks. No models, LPs or endpoints run."""
import hashlib
import json
import math
from pathlib import Path
import tarfile
import xml.etree.ElementTree as ET

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]


def read(path):
    return json.loads(path.read_text(), parse_constant=lambda s: (_ for _ in ()).throw(ValueError(s)))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def same(a, b):
    if isinstance(a, dict):
        return isinstance(b, dict) and a.keys() == b.keys() and all(same(a[k], b[k]) for k in a)
    if isinstance(a, list):
        return isinstance(b, list) and len(a) == len(b) and all(same(x, y) for x, y in zip(a, b))
    if isinstance(a, (int, float)) and not isinstance(a, bool):
        return isinstance(b, (int, float)) and math.isclose(a, b, rel_tol=1e-6, abs_tol=1e-7)
    return a == b


def masked(folder, stem='masked', expected_cells=None):
    result = read(folder / (stem + '-result.json'))
    oracle = read(folder / (stem + '-oracle.json'))
    cells = {(c['condition_id'], c['hypothesis_id']): c for c in result['cells']}
    if expected_cells is not None:
        assert len(cells) == len(result['cells']) == expected_cells
    assert result['promotion'] == 'none' and result['calibrated_probability'] is None
    assert len(result['solver_calls']) <= result['max_solver_calls']
    expected = {(c, h['id']) for c in result['visited_conditions'] for h in result['public_problem']['hypotheses']}
    assert set(cells) == expected
    assert len(oracle['comparisons']) == len(cells) * 2 == oracle['n_comparisons'] == oracle['n_agree']
    for row in oracle['comparisons']:
        actual = cells[(row['condition_id'], row['hypothesis_id'])]['phenotypes'][row['perturbation']]
        actual = {k: v for k, v in actual.items() if k != 'knocked_out_reactions'}
        assert row['agrees'] and same(actual, row['oracle'])
    assert oracle['n_cell_comparisons'] == oracle['n_cell_agree'] == len(oracle['cell_comparisons'])
    for row in oracle['cell_comparisons']:
        actual = cells[(row['condition_id'], row['hypothesis_id'])]
        assert row['agrees'] and all(actual[k] == v for k, v in row['oracle'].items())
    return result


def main():
    inventory = {}
    for line in (HERE / 'SHA256SUMS').read_text().splitlines():
        sha, name = line.split('  ', 1)
        assert name not in inventory and digest(HERE / name) == sha, name
        inventory[name] = sha
    assert set(inventory) == {str(p.relative_to(HERE)) for p in HERE.rglob('*')
                              if p.is_file() and p.name != 'SHA256SUMS' and '__pycache__' not in p.parts}
    record = read(HERE / 'record.json')
    for source in record['source_inputs']:
        assert digest(ROOT / source['path']) == source['sha256']
    counts = {}
    with tarfile.open(HERE / 'jvm-test-results.tar.gz') as archive:
        for member in archive:
            if not member.name.endswith('.xml'):
                continue
            module = member.name.split('/')[0]
            totals = counts.setdefault(module, dict(tests=0, failures=0, errors=0, skipped=0))
            root = ET.fromstring(archive.extractfile(member).read())
            for key in totals:
                totals[key] += int(root.get(key, 0))
    assert counts == record['jvm_counts']
    assert sum(c['tests'] for c in counts.values()) == 1436
    assert sum(c['skipped'] for c in counts.values()) == 6
    assert all(c['failures'] == c['errors'] == 0 for c in counts.values())
    for name, n in [('tests.log', 159), ('benchmark-tests.log', 121)]:
        log = (HERE / 'final-gates' / name).read_text()
        assert f'Ran {n} tests' in log and log.rstrip().endswith('OK')
    assert (HERE / 'final-gates/exit-code.txt').read_text().strip() == '0'
    parity = read(HERE / 'parity/summary.json')
    assert len(parity['non_agent_pairs']) == 6
    for row in parity['non_agent_pairs']:
        prefix = f"{row['strategy']}-{row['seed']}"
        assert read(HERE / 'parity' / (prefix + '-before.json')) == read(HERE / 'parity' / (prefix + '-after.json'))
    folder = HERE / 'diagnostic/sco-gem'
    standard = read(folder / 'standard-result.json')
    model_sha = next(x['sha256'] for x in record['source_inputs'] if x['path'].endswith('Sco-GEM.xml'))
    assert standard['hashes']['model_sha256'] == model_sha
    experiments = {r['experiment_id']: r for r in standard['experiments']}
    comparisons = read(folder / 'standard-oracle.json')
    assert len(comparisons) == 12
    for row in comparisons:
        product = experiments[row['experiment']]['measurements']['product']
        actual = {'minimum': product['min_synthesis_at_required_state'], 'maximum': product['max_synthesis'],
                  'statuses': [product['lp_status'][s] for s in ('min', 'max')]}
        assert row['agrees'] and same(actual, row['observed']) and same(actual, row['expected'])
        assert product['target_kind'] == 'intracellular_synthesis_capacity' and 'max_secretion' not in product
        if '@knockout:' in row['experiment']:
            assert actual['maximum'] <= 1e-7
        else:
            assert actual['maximum'] > 0 and abs(actual['minimum']) <= 1e-7
    result = masked(folder, expected_cells=36)
    assert result['input_model_sha256'] == model_sha
    for cell in result['cells']:
        if int(cell['condition_id'][1:]) % 3 == 2:
            assert cell['product_status'] == 'inconclusive' and cell['status'] == 'inconclusive'
        elif cell['hypothesis_id'] == 'H1':
            assert cell['product_contrast'] == 'synthesis_capacity_lost'
    masked(folder, 'masked-qwen')
    for kind in ('standard', 'masked'):
        result = read(folder / (kind + '-qwen-result.json'))
        requests = read(folder / (kind + '-qwen-requests.json'))
        transcript = read(folder / (kind + '-qwen-transcript.json'))
        assert len(requests) == len(transcript) and 1 <= len(requests) <= 3 and all(r['ok'] for r in requests)
        choices = result['selection_log' if kind == 'standard' else 'selection']
        assert len(choices) == len(transcript)
        for choice, turn in zip(choices, transcript):
            observation = json.loads(turn['messages'][1]['content'])
            if kind == 'standard':
                assert observation['schema'] == 'gspa.metabolic_planner_observation/3'
                encoded = json.dumps(turn['messages'], sort_keys=True, allow_nan=False).encode()
                assert choice['valid'] and choice['source'] == 'qwen' and choice['prompt_sha256'] == hashlib.sha256(encoded).hexdigest()
            else:
                assert choice['status'] == 'agent_selected' and observation['product_requirement'] == 'maintenance_only'
    cli = HERE / 'diagnostic/cli'
    # Production CLI directories retain the actual generated manifests/results.
    for kind in ('standard', 'masked'):
        result = read(cli / (kind + '.tsv.experiments/result.json'))
        assert result['promotion'] == 'none'
        ledger = read(cli / (kind + '.tsv.stages.json'))
        stage = next(s for s in ledger['stages'] if s['name'] == 'metabolic_experiments')
        assert stage['status'] == 'executed'
    oracle = read(cli / 'masked-oracle.json')
    assert oracle['n_agree'] == oracle['n_comparisons'] == 24
    assert oracle['n_cell_agree'] == oracle['n_cell_comparisons'] == 6
    print(json.dumps({'files_verified': len(inventory), 'jvm_tests': 1436, 'python_tests': 280,
                      'standard_capacity_comparisons': 12, 'masked_cells': 36, 'status': 'verified'}))


if __name__ == '__main__':
    main()
