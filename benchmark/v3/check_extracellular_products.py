#!/usr/bin/env python3
"""Real-model secretion route check; not measured production or function recovery."""
import argparse
import copy
import hashlib
import json
from pathlib import Path
import sys
import subprocess

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / 'scripts'), str(Path(__file__).parent)]
import metabolic_experiments as ME
import masked_assignment_experiments as MA
import check_masked_assignment_cobra as CO
import run_germicidin_diagnostic as GD
import cobra

PRODUCTS = ['R_EX_act_e', 'R_EX_gACT_e']
GENE = 'G_SCO5089'


def save(path, value):
    with path.open('x') as f:
        f.write(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')


def check(model_path, medium, out):
    out.mkdir(exist_ok=False)
    model = ME.load_sbml_model(model_path)
    base = GD.manifest(model, medium, 'published_model_minimal')
    base['config'].update(strategy='fixed', product_requirement='maintenance_only',
        max_experiments=40, max_solver_calls=500, max_focus_reactions=0, max_baseline_focus_reactions=0)
    base['context'] = {'hypotheses': [{'id': 'H', 'gene': GENE}], 'evidence': [
        {'id': f'E{i}_{j}', 'type': 'bgc', 'hypothesis_id': 'H', 'product_exchange': product,
         'condition_hint': 'supplement_removal', 'exchange': nutrient}
        for i, product in enumerate(PRODUCTS) for j, nutrient in enumerate(['R_EX_pi_e', 'R_EX_glc__D_e'])]}
    masked = copy.deepcopy(base)
    masked['context'] = {}
    masked['assignment_recovery'] = {'mask_gene': GENE,
        'hypotheses': [{'id': 'H0', 'assignments': {}},
                       {'id': 'H1', 'assignments': {'R_ACTS1': [[MA.TARGET]]}},
                       {'id': 'H2', 'assignments': {'R_GERMA': [[MA.TARGET]]}}],
        'conditions': [{'id': 'C' + str(i * 3 + j), 'mods': mods, 'product_exchange': product,
                        'evidence_channels': ['bgc']}
            for i, product in enumerate(PRODUCTS) for j, mods in enumerate([
                [], [{'op': 'close', 'exchange': 'R_EX_pi_e'}],
                [{'op': 'close', 'exchange': 'R_EX_glc__D_e'}]])]}
    save(out / 'registration.json', {'model_sha256': hashlib.sha256(model_path.read_bytes()).hexdigest(),
        'medium_sha256': hashlib.sha256(medium.read_bytes()).hexdigest(), 'standard': base, 'masked': masked,
        'scope': 'Pinned Sco-GEM mechanism check with controlled hypotheses; no assay outcomes, independent function truth, or discovery score.'})
    original = cobra.io.read_sbml_model(model_path, f_replace={})
    original.solver = 'glpk'
    result = ME.run_experiments(base, model)
    save(out / 'standard-result.json', result)
    items = {r['id']: r for r in result['menu']}
    comparisons = []
    for row in result['experiments']:
        if not row['experiment_id'].startswith('product:'):
            continue
        item = items[row['experiment_id']]
        expected = CO.simulate(original, base, {'assignments': {}}, item, result['config'], False)
        measured = row['measurements'].get('product', {})
        interval = expected.get('product_interval')
        if interval and interval['statuses'] == ['optimal', 'optimal']:
            match = (CO.agrees(row['measurements']['max_growth'], expected['growth']) and
                     CO.agrees(measured.get('max_secretion'), max(0, interval['maximum'])) and
                     CO.agrees(measured.get('min_secretion_at_required_state'), max(0, interval['minimum'])))
        else:
            match = (expected['status'] == 'infeasible' and row.get('outcome') == 'infeasible_with_maintenance'
                     and 'max_secretion' not in measured)
        comparisons.append({'experiment': row['experiment_id'], 'agrees': match, 'oracle': expected,
                            'production': row})
    save(out / 'standard-oracle.json', comparisons)
    assert len(comparisons) == 8 and all(r['agrees'] for r in comparisons), 'standard secretion comparison failed'
    result = ME.run_experiments(masked, model)
    save(out / 'masked-result.json', result)
    report = CO.compare(model_path, masked, result)
    save(out / 'masked-oracle.json', report)
    assert len(result['cells']) == 18
    assert report['n_agree'] == report['n_comparisons'] == 36
    assert report['n_cell_agree'] == report['n_cell_comparisons'] == 18
    summary = {'status': 'verified', 'standard_product_comparisons': len(comparisons),
               'masked_phenotype_comparisons': report['n_comparisons'], 'masked_cell_comparisons': report['n_cell_comparisons'],
               'scope': 'Independent numerical agreement only; no measured-production accuracy or recovered-function claim.'}
    save(out / 'verification.json', summary)
    return summary


def cli_check(out, jar, before_jar):
    """Compare actual evidence-guided CLI routes; fixture chemistry is synthetic."""
    import libsbml
    sys.path.insert(0, str(ROOT / 'scripts/tests'))
    from test_metabolic_experiments import manifest, write_sbml
    out.mkdir(exist_ok=False)
    model_path = out / 'model.xml'
    write_sbml(model_path)
    doc = libsbml.readSBMLFromFile(str(model_path))
    for rid, ec in [('SEC', '4.4.4.4'), ('ACS', '6.2.1.1')]:
        reaction = doc.getModel().getReaction(rid)
        reaction.setMetaId('meta_' + rid)
        cv = libsbml.CVTerm(libsbml.BIOLOGICAL_QUALIFIER)
        cv.setBiologicalQualifierType(libsbml.BQB_IS)
        cv.addResource('https://identifiers.org/ec-code/' + ec)
        reaction.addCVTerm(cv)
    libsbml.writeSBMLToFile(doc, str(model_path))
    claims = out / 'claims.jsonl'
    claims.write_text(''.join(json.dumps({'protein_id': 'locus_g6', 'function_type': 'EC',
        'function_id': 'EC:' + ec, 'source': 'diamond', 'evidence_type': 'SEQUENCE_SIMILARITY',
        'raw_score': 0.6, 'calibrated_prob': 0.6,
        'metadata': {'product_exchange': 'EX_sec_e', 'substrate_exchange': 'EX_ac_e', 'substrate_role': 'carbon'}})
        + '\n' for ec in ['4.4.4.4', '6.2.1.1']))
    reports = []
    for version, artifact in [('before', before_jar), ('after', jar)]:
        for route in ['standard', 'masked']:
            template = manifest({}, strategy='fixed', max_focus_reactions=0,
                max_baseline_focus_reactions=0, product_requirement='maintenance_only')
            template['model'] = str(model_path)
            if route == 'masked':
                template['assignment_recovery'] = {'mask_gene': 'g6', 'candidate_protein': 'locus_g6'}
            stem = version + '-' + route
            path = out / (stem + '-template.json')
            save(path, template)
            target = out / (stem + '.tsv')
            command = ['java', '-jar', str(artifact), 'integrate', '--claims', str(claims), '--out', str(target),
                '--lite', '--enable-priors', 'none', '--experiment-manifest', str(path),
                '--experiment-mode', 'evidence_guided', '--experiment-strategy', 'fixed',
                '--experiment-max-solver-calls', '200', '--experiment-max-experiments', '30',
                '--experiment-python-bin', sys.executable,
                '--experiment-sidecar', str(ROOT / 'scripts/metabolic_experiments.py')]
            save(out / (stem + '-command.json'), command)
            with (out / (stem + '.log')).open('x') as log:
                process = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
            assert process.returncode == 0, (stem, process.returncode)
            folder = Path(str(target) + '.experiments')
            result = json.loads((folder / 'result.json').read_text())
            generated = json.loads((folder / 'manifest.json').read_text())
            assert result['promotion'] == 'none'
            if route == 'masked':
                conditions = result['public_problem']['conditions']
                combined = [c for c in conditions if c.get('product_exchange') and c['mods']]
                assert bool(combined) == (version == 'after'), (stem, combined)
                report = CO.compare(model_path, generated, result)
                save(out / (stem + '-oracle.json'), report)
                assert report['n_comparisons'] > 0 and report['n_agree'] == report['n_comparisons']
                assert report['n_cell_comparisons'] > 0 and report['n_cell_agree'] == report['n_cell_comparisons']
                reports.append({'run': stem, 'combined_conditions': len(combined),
                                'oracle_phenotypes': report['n_comparisons'], 'oracle_product_cells': report['n_cell_comparisons']})
            else:
                combined = [r for r in result['experiments'] if r.get('mods') and r['measurements'].get('product')]
                assert combined
                reports.append({'run': stem, 'combined_experiments': len(combined)})
    save(out / 'verification.json', {'status': 'verified', 'runs': reports,
        'scope': 'Before/after Java builder with the same repaired Python engine. Synthetic CLI integration, not function accuracy.'})


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', required=True, type=Path)
    p.add_argument('--model', type=Path)
    p.add_argument('--medium', type=Path)
    p.add_argument('--jar', type=Path)
    p.add_argument('--before-jar', type=Path)
    args = p.parse_args()
    if args.jar or args.before_jar:
        if not args.jar or not args.before_jar:
            p.error('--jar and --before-jar must be supplied together')
        cli_check(args.out, args.jar, args.before_jar)
    else:
        if not args.model or not args.medium:
            p.error('--model and --medium are required for the real-model check')
        print(json.dumps(check(args.model, args.medium, args.out), sort_keys=True))
