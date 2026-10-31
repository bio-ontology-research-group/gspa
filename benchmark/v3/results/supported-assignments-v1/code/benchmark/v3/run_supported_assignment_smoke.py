"""Exercise supported complex/multifunction construction through CLI and independent GLPK.

All sequences, evidence and network changes are synthetic mechanism fixtures.
No measured biological recovery or strategy-superiority claim is made.
"""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import libsbml
from run_guided_assignment_smoke import ROOT, save, sha
import check_masked_assignment_cobra as oracle


def run(out, jar):
    out.mkdir(parents=True, exist_ok=False)
    scratch = out / 'scratch'
    scratch.mkdir()
    env = dict(os.environ, TMPDIR=str(scratch), JAVA_TOOL_OPTIONS=f'-Djava.io.tmpdir={scratch}')
    engine = out / 'engine.jar'
    shutil.copyfile(jar, engine)
    generator = ROOT / 'gspa-core/src/test/resources/experiments/make_fixtures.py'
    sources = [Path(__file__).resolve(), Path(__file__).with_name('run_guided_assignment_smoke.py'),
        ROOT / 'benchmark/v3/check_masked_assignment_cobra.py', generator, engine,
        ROOT / 'scripts/metabolic_experiments.py', ROOT / 'scripts/masked_assignment_experiments.py',
        ROOT / 'gspa-core/src/main/groovy/gspa/integration/contextbridge/MaskedAssignmentBuilder.groovy',
        ROOT / 'gspa-core/src/main/groovy/gspa/integration/contextbridge/SupportedAssignmentHypotheses.groovy',
        ROOT / 'gspa-core/src/main/groovy/gspa/integration/contextbridge/ExperimentStage.groovy']
    identity = {str(p): sha(p) for p in sources}
    fixtures = out / 'fixtures'
    subprocess.run([sys.executable, str(generator), str(fixtures)], check=True, env=env, capture_output=True)
    model = out / 'supported.xml'
    document = libsbml.readSBMLFromFile(str(fixtures / 'xref/xref.xml'))
    network = document.getModel()
    for gid in ('gPartner', 'gIsozyme'):
        gene = network.getPlugin('fbc').createGeneProduct()
        gene.setId(gid)
        gene.setLabel(gid)
    fbc = network.getReaction('ACS').getPlugin('fbc')
    fbc.unsetGeneProductAssociation()
    either = fbc.createGeneProductAssociation().createOr()
    for gene in ('gD', 'gIsozyme'):
        either.createGeneProductRef().setGeneProduct(gene)
    libsbml.writeSBMLToFile(document, str(model))
    template = json.loads((fixtures / 'xref/xref_template.json').read_text())
    template['model'] = str(model)
    template['assignment_recovery'] = {'mask_gene': 'pD', 'candidate_protein': 'query', 'conditions': [
        {'mods': [{'op': 'carbon_swap', 'exchange': 'EX_ac_e'},
                  {'op': 'knockout', 'reactions': [], 'genes': ['gIsozyme']}], 'evidence_channels': ['domain']},
        {'mods': [{'op': 'carbon_swap', 'exchange': 'EX_ac_e'},
                  {'op': 'knockout', 'reactions': [], 'genes': ['gIsozyme']},
                  {'op': 'knockout', 'reactions': [], 'genes': ['gPartner']}], 'evidence_channels': ['domain']}]}
    manifest = out / 'template.json'
    save(manifest, template)
    evidence = out / 'evidence.json'
    save(evidence, [{'id': 'complex', 'function_key': 'query|EC|EC:6.2.1.1', 'type': 'domain',
                    'source': 'synthetic independent complex support', 'score': .6,
                    'required_genes': ['gPartner'], 'assignment_group': 'synthetic_architecture'},
                   {'id': 'second_activity', 'function_key': 'query|EC|EC:1.9.3.1', 'type': 'domain',
                    'source': 'synthetic independent multifunction support', 'score': .6,
                    'required_genes': [], 'assignment_group': 'synthetic_architecture'}])
    claims = out / 'claims.jsonl'
    claims.write_text(''.join(json.dumps({'protein_id': 'query', 'function_type': 'EC', 'function_id': 'EC:' + ec,
        'source': 'fixture_predictor', 'evidence_type': 'SEQUENCE_SIMILARITY', 'raw_score': .1}) + '\n'
        for ec in ('6.2.1.1', '1.9.3.1')))
    # Change only the target-containing GPR clauses; retain the independent isozyme.
    fbc.unsetGeneProductAssociation()
    fbc.createGeneProductAssociation().createGeneProductRef().setGeneProduct('gIsozyme')
    network.getReaction('RESP').getPlugin('fbc').unsetGeneProductAssociation()
    network.getReaction('VITS').getPlugin('fbc').createGeneProductAssociation().createGeneProductRef().setGeneProduct('gD')
    changed_model = out / 'changed-hidden-gpr.xml'
    libsbml.writeSBMLToFile(document, str(changed_model))
    changed_manifest = out / 'changed-template.json'
    save(changed_manifest, {**template, 'model': str(changed_model)})
    for p in (model, manifest, evidence, claims, changed_model, changed_manifest):
        identity[str(p)] = sha(p)
    runs = {}
    for name, strategy, limit, source in [
            ('single_medium', 'fixed', 1, manifest), ('fixed', 'fixed', 20, manifest),
            ('random', 'random', 20, manifest), ('evidence', 'evidence', 20, manifest),
            ('changed_hidden_assignment', 'fixed', 20, changed_manifest)]:
        folder = out / name
        folder.mkdir()
        output = folder / 'predictions.tsv'
        command = ['java', '-jar', str(engine), 'integrate', '--lite', '--enable-priors', 'none',
            '--claims', str(claims), '--out', str(output), '--experiment-manifest', str(source),
            '--experiment-evidence', str(evidence), '--experiment-mode', 'evidence_guided',
            '--experiment-strategy', strategy, '--experiment-python-bin', sys.executable,
            '--experiment-sidecar', str(ROOT / 'scripts/metabolic_experiments.py'),
            '--experiment-max-experiments', str(limit), '--experiment-max-solver-calls', '256', '--experiment-seed', '42']
        save(folder / 'command.json', command)
        with (folder / 'execution.log').open('w') as log:
            subprocess.run(command, check=True, env=env, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT, timeout=120)
        artifact = Path(str(output) + '.experiments')
        result = json.loads((artifact / 'result.json').read_text())
        generated = json.loads((artifact / 'manifest.json').read_text())
        check = oracle.compare(Path(generated['model']), generated, result)
        save(folder / 'independent_check.json', check)
        assert check['n_agree'] == check['n_comparisons'], check
        assert result['promotion'] == 'none' and result['calibrated_probability'] is None
        assert len(result['solver_calls']) <= 256
        hypotheses = result['public_problem']['hypotheses']
        assert len(hypotheses) == 4, hypotheses
        assert any(h['assignments'] == {'ACS': [['MASKED_TARGET', 'gPartner']], 'RESP': [['MASKED_TARGET']]}
                   for h in hypotheses)
        assert not any(h['assignments'].get('ACS') == [['MASKED_TARGET']] for h in hypotheses)
        runs[name] = {'result': result, 'summary': {'comparisons': check['n_comparisons'],
            'agree': check['n_agree'], 'conditions': len(result['visited_conditions']),
            'solver_calls': len(result['solver_calls']), 'result_sha256': sha(artifact / 'result.json')}}
    for entry in runs.values():
        assert entry['result']['public_problem'] == runs['fixed']['result']['public_problem']
    assert runs['fixed']['result']['cells'] == runs['changed_hidden_assignment']['result']['cells']
    for path, expected in identity.items():
        assert sha(path) == expected, path
    save(out / 'input-sha256.json', identity)
    summary = {'schema': 'gspa.supported-assignment-smoke/1', 'status': 'verified',
        'scope': 'Synthetic complex/multifunction and combined-KO mechanism; no biological recovery claim.',
        'hidden_assignment_invariance': True, 'same_public_problem_all_arms': True,
        'runs': {k: v['summary'] for k, v in runs.items()}}
    save(out / 'summary.json', summary)
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--jar', type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(run(args.out.resolve(), args.jar.resolve()), indent=2))
