#!/usr/bin/env python3
"""Production CLI/independent solver check for evidence-guided masked assignments.

Uses a declared synthetic network and synthetic predictor claims. This is a
mechanism test, never a biological function-recovery or planner-accuracy result.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import libsbml
import check_masked_assignment_cobra as oracle

ROOT = Path(__file__).resolve().parents[2]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')


def run(out, jar):
    out.mkdir(parents=True, exist_ok=False)
    scratch = out / 'scratch'
    scratch.mkdir()
    env = dict(os.environ, TMPDIR=str(scratch), JAVA_TOOL_OPTIONS=f'-Djava.io.tmpdir={scratch}')
    engine = out / 'engine.jar'
    shutil.copyfile(jar, engine)
    sources = [Path(__file__).resolve(), ROOT / 'benchmark/v3/check_masked_assignment_cobra.py',
               ROOT / 'scripts/metabolic_experiments.py', ROOT / 'scripts/masked_assignment_experiments.py',
               ROOT / 'gspa-core/src/main/groovy/gspa/integration/contextbridge/MaskedAssignmentBuilder.groovy',
               ROOT / 'gspa-core/src/main/groovy/gspa/integration/contextbridge/ExperimentStage.groovy',
               ROOT / 'gspa-core/src/main/groovy/gspa/integration/workflow/IntegrationWorkflow.groovy',
               ROOT / 'gspa-core/src/test/resources/experiments/make_fixtures.py', engine]
    identity = {str(p): sha(p) for p in sources}
    fixtures = out / 'fixtures'
    subprocess.run([sys.executable, str(sources[-2]), str(fixtures)], check=True, env=env, capture_output=True)
    model = fixtures / 'xref/xref.xml'
    template = json.loads((fixtures / 'xref/xref_template.json').read_text())
    template['model'] = str(model)
    template['assignment_recovery'] = {'mask_gene': 'pD', 'candidate_protein': 'query'}
    manifest = out / 'template.json'
    save(manifest, template)
    claims = out / 'claims.jsonl'
    records = [
        {'protein_id': 'query', 'function_type': 'EC', 'function_id': 'EC:6.2.1.1', 'source': 'fixture_diamond',
         'evidence_type': 'SEQUENCE_SIMILARITY', 'raw_score': .00001},
        {'protein_id': 'query', 'function_type': 'EC', 'function_id': 'EC:4.4.4.4', 'source': 'fixture_domain',
         'evidence_type': 'SEQUENCE_DOMAIN', 'raw_score': .3}]
    claims.write_text(''.join(json.dumps(r, sort_keys=True) + '\n' for r in records))
    # Change only the private target's original GPR, keeping every candidate input fixed.
    document = libsbml.readSBMLFromFile(str(model))
    for reaction in document.getModel().getListOfReactions():
        fbc = reaction.getPlugin('fbc')
        if fbc.isSetGeneProductAssociation() and 'gD' in fbc.getGeneProductAssociation().toXMLNode().toXMLString():
            fbc.unsetGeneProductAssociation()
    document.getModel().getReaction('VITS').getPlugin('fbc').createGeneProductAssociation().createGeneProductRef().setGeneProduct('gD')
    changed_model = out / 'changed-hidden-gpr.xml'
    libsbml.writeSBMLToFile(document, str(changed_model))
    changed_manifest = out / 'changed-template.json'
    save(changed_manifest, {**template, 'model': str(changed_model)})
    for p in (model, manifest, claims, changed_model, changed_manifest):
        identity[str(p)] = sha(p)
    runs = {}
    for name, strategy, max_experiments, source in [
            ('single_medium', 'fixed', 1, manifest), ('fixed', 'fixed', 2, manifest),
            ('random', 'random', 2, manifest), ('evidence', 'evidence', 2, manifest),
            ('changed_hidden_assignment', 'fixed', 2, changed_manifest)]:
        folder = out / name
        folder.mkdir()
        output = folder / 'predictions.tsv'
        command = ['java', '-jar', str(engine), 'integrate', '--lite', '--enable-priors', 'none',
                   '--claims', str(claims), '--out', str(output), '--experiment-manifest', str(source),
                   '--experiment-mode', 'evidence_guided', '--experiment-strategy', strategy,
                   '--experiment-python-bin', sys.executable,
                   '--experiment-sidecar', str(ROOT / 'scripts/metabolic_experiments.py'),
                   '--experiment-max-experiments', str(max_experiments), '--experiment-max-solver-calls', '20',
                   '--experiment-seed', '42']
        save(folder / 'command.json', command)
        with (folder / 'execution.log').open('w') as log:
            subprocess.run(command, check=True, env=env, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
        artifact = Path(str(output) + '.experiments')
        result = json.loads((artifact / 'result.json').read_text())
        constructed = json.loads((artifact / 'manifest.json').read_text())
        report = oracle.compare(constructed['model'], constructed, result)
        save(folder / 'independent-oracle.json', report)
        assert report['n_agree'] == report['n_comparisons']
        ledger = json.loads(Path(str(output) + '.stages.json').read_text())
        stage = next(s for s in ledger['stages'] if s['name'] == 'metabolic_experiments')
        assert stage['status'] == 'executed' and stage['context_source'] == 'evidence_guided_masked_assignments'
        assert stage['promotion'] == 'none' and stage['candidates_with_conditional_features'] == 0
        assert len(result['solver_calls']) == (10 if max_experiments == 1 else 20)
        assert set(h['reaction'] for h in stage['context']['hypotheses']) == {'ACS', 'SEC'}
        runs[name] = {'result': result, 'summary': {
            'solver_calls': len(result['solver_calls']), 'cells': len(result['cells']),
            'oracle_agreement': report['n_agree'], 'oracle_comparisons': report['n_comparisons'],
            'public_problem_sha256': result['public_problem_sha256'], 'visited_conditions': result['visited_conditions'],
            'result_sha256': sha(artifact / 'result.json'), 'stage_ledger_sha256': sha(Path(str(output) + '.stages.json'))}}
    base = runs['fixed']['result']
    for name in runs:
        assert runs[name]['result']['public_problem'] == base['public_problem'], name
    assert runs['changed_hidden_assignment']['result']['cells'] == base['cells']
    for path, expected in identity.items():
        assert sha(path) == expected, path
    save(out / 'input-sha256.json', identity)
    summary = {'schema': 'gspa.guided-masked-assignment-smoke/1', 'status': 'verified',
               'data': 'synthetic network and predictor claims; no biological accuracy claim',
               'hidden_assignment_invariance': True, 'same_public_problem_all_arms': True,
               'runs': {k: v['summary'] for k, v in runs.items()}}
    save(out / 'summary.json', summary)
    return summary


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out', required=True, type=Path)
    ap.add_argument('--jar', type=Path, default=ROOT / 'gspa-cli/build/libs/gspa-3.0.0.jar')
    args = ap.parse_args()
    print(json.dumps(run(args.out.resolve(), args.jar.resolve()), indent=2))
