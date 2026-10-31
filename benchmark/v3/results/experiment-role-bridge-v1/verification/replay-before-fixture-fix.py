#!/usr/bin/env python3
"""Before/after CLI proof for raw candidate roles in masked and absence experiments.

Uses synthetic claims and a synthetic network; this does not measure role accuracy.
"""
import argparse
import copy
import json
import os
from pathlib import Path
import subprocess
import sys

from check_absence_cobra import compare as absence_compare, sha
from check_masked_assignment_cobra import compare as masked_compare

ROOT = Path(__file__).resolve().parents[2]


def save(path, data):
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + '\n')


def run(old_jar, new_jar, out):
    out.mkdir(parents=True, exist_ok=False)
    scratch = out / 'scratch'
    scratch.mkdir()
    env = dict(os.environ, TMPDIR=str(scratch), JAVA_TOOL_OPTIONS='-Djava.io.tmpdir=' + str(scratch))
    fixture = ROOT / 'gspa-core/src/test/resources/experiments/make_fixtures.py'
    subprocess.run([sys.executable, '-I', str(fixture), str(out / 'models')], check=True, capture_output=True, env=env)
    template = json.loads((out / 'models/xref/xref_template.json').read_text())
    template['model'] = str(out / 'models/xref/xref.xml')
    template['assignment_recovery'] = {'mask_gene': 'pD', 'candidate_protein': 'query'}
    support = out / 'support.json'
    save(support, [{'id': 'complex', 'function_key': 'query|EC|EC:6.2.1.1', 'type': 'domain',
                   'source': 'synthetic independent partner evidence', 'score': .6, 'required_genes': ['pB']}])
    cases = [
        ('default', [], False, 'masked', 'executed'),
        ('accessory', ['ACCESSORY'], True, 'masked', 'no_admissible_experiment'),
        ('regulator', ['REGULATOR'], True, 'masked', 'no_admissible_experiment'),
        ('transporter', ['TRANSPORTER'], True, 'masked', 'no_admissible_experiment'),
        ('subunit', ['COMPLEX_SUBUNIT'], False, 'masked', 'no_admissible_experiment'),
        ('catalytic_subunit', ['CATALYTIC', 'COMPLEX_SUBUNIT'], False, 'masked', 'no_admissible_experiment'),
        ('supported_subunit', ['COMPLEX_SUBUNIT'], True, 'masked', 'executed'),
        ('absence_default', [], False, 'absence', 'executed'),
        ('absence_catalytic_subunit', ['CATALYTIC', 'COMPLEX_SUBUNIT'], False, 'absence', 'failed'),
    ]
    commands, results = [], []
    for name, roles, partners, kind, expected in cases:
        folder = out / name
        folder.mkdir()
        claims = folder / 'claims.jsonl'
        claims.write_text(''.join(json.dumps({'protein_id': 'query', 'function_type': 'EC',
            'function_id': 'EC:6.2.1.1', 'source': 'synthetic_predictor', 'evidence_type': 'SEQUENCE_SIMILARITY',
            'raw_score': .001, 'metadata': {'candidate_role': role} if role else {}}) + '\n' for role in roles or [None]))
        manifest = copy.deepcopy(template)
        if kind == 'absence':
            manifest.pop('assignment_recovery')
            manifest['absence_comparison'] = {'reaction': 'ACS', 'candidate_function_key': 'query|EC|EC:6.2.1.1',
                'assembly': {'status': 'UNKNOWN', 'source': 'synthetic test'},
                'function_search': {'status': 'UNKNOWN', 'source': 'synthetic test'}, 'external_supply': None,
                'conditions': [{'id': 'control', 'mods': []}]}
        source = folder / 'template.json'
        save(source, manifest)
        observations = {}
        for label, jar in [('before', old_jar), ('after', new_jar)]:
            output = folder / (label + '.tsv')
            command = ['java', '-jar', str(jar), 'integrate', '--lite', '--enable-priors', 'none',
                '--claims', str(claims), '--out', str(output), '--experiment-manifest', str(source),
                '--experiment-mode', 'evidence_guided' if kind == 'masked' else 'manifest',
                '--experiment-strategy', 'fixed', '--experiment-python-bin', sys.executable,
                '--experiment-sidecar', str(ROOT / 'scripts/metabolic_experiments.py'),
                '--experiment-max-experiments', '2', '--experiment-max-solver-calls', '80']
            if partners:
                command += ['--experiment-evidence', str(support)]
            proc = subprocess.run(command, cwd=ROOT, env=env, capture_output=True, text=True, timeout=120)
            (folder / (label + '.log')).write_text(proc.stdout + proc.stderr)
            commands.append(dict(case=name, version=label, argv=command, returncode=proc.returncode))
            save(out / 'commands.json', commands)
            ledger = json.loads(Path(str(output) + '.stages.json').read_text())
            stage = next(x for x in ledger['stages'] if x['name'] == 'metabolic_experiments')
            observed = 'failed' if proc.returncode else stage.get('outcome') if stage.get('outcome') == 'no_admissible_experiment' else 'executed'
            assert observed == ('executed' if label == 'before' else expected), (name, label, observed)
            artifact = Path(str(output) + '.experiments')
            comparisons = 0
            if observed == 'executed':
                result = json.loads((artifact / 'result.json').read_text())
                generated = json.loads((artifact / 'manifest.json').read_text())
                if kind == 'masked':
                    oracle = masked_compare(Path(generated['model']), generated, result)
                    assert oracle['n_agree'] == oracle['n_comparisons']
                    if label == 'after' and roles:
                        assert all(h['assignments'].get('ACS') in (None, [['MASKED_TARGET', 'gB']]) for h in result['public_problem']['hypotheses'])
                        assert stage['context']['hypotheses'][0]['candidates'][0]['roles'] == roles
                        assert stage['context']['hypotheses'][0]['candidates'][0]['role_source'] == 'raw_claims'
                else:
                    oracle = absence_compare(artifact / 'manifest.json', artifact / 'result.json')
                    assert oracle['agrees']
                comparisons = len(oracle['comparisons'])
                save(folder / (label + '.oracle.json'), oracle)
            else:
                assert not (artifact / 'result.json').exists()
                if observed == 'failed':
                    assert stage['status'] == 'failed'
                else:
                    assert stage['solver_calls_used'] == 0
            observations[label] = dict(outcome=observed, oracle_comparisons=comparisons)
        if expected != 'failed':
            assert (folder / 'before.tsv').read_bytes() == (folder / 'after.tsv').read_bytes()
        if expected == 'executed':
            old_result = json.loads((folder / 'before.tsv.experiments/result.json').read_text())
            new_result = json.loads((folder / 'after.tsv.experiments/result.json').read_text())
            assert old_result['cells'] == new_result['cells']
            if kind == 'masked':
                assert old_result['public_problem'] == new_result['public_problem']
        results.append(dict(case=name, observations=observations, annotation_bytes_unchanged=expected != 'failed'))
        print(json.dumps(results[-1]), flush=True)
    report = dict(status='verified', scope='synthetic raw-role contract, not measured biological accuracy', cases=results,
        inputs={str(p): sha(p) for p in [old_jar, new_jar, Path(__file__).resolve(), fixture,
            ROOT / 'benchmark/v3/check_absence_cobra.py', ROOT / 'benchmark/v3/check_masked_assignment_cobra.py']},
        outputs={str(p.relative_to(out)): sha(p) for p in sorted(out.rglob('*')) if p.is_file()})
    save(out / 'verification.json', report)
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--old-jar', required=True, type=Path)
    parser.add_argument('--new-jar', required=True, type=Path)
    parser.add_argument('--out', required=True, type=Path)
    args = parser.parse_args()
    run(args.old_jar.resolve(), args.new_jar.resolve(), args.out.resolve())
