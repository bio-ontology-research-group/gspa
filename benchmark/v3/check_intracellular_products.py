#!/usr/bin/env python3
"""Intracellular-target workflow diagnostic, not a function accuracy benchmark.

Run on a compute node. The synthetic CLI fixture checks evidence-guided routing;
unmodified Sco-GEM checks native and masked conditional capacity. No assay labels
are read. Supplied candidate assignments are deliberately controlled hypotheses.
"""
import argparse
import copy
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / 'scripts'), str(ROOT / 'scripts/tests'), str(Path(__file__).parent)]
import metabolic_experiments as ME
import masked_assignment_experiments as MA
import check_masked_assignment_cobra as CO
import run_germicidin_diagnostic as GD
from run_planner_benchmark import UsageChat
from test_metabolic_experiments import REACTIONS, manifest, write_sbml
import cobra
import libsbml


def save(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')


def cli_fixture(out, jar):
    """No transport for the product; exact EC annotations are synthetic fixtures."""
    out.mkdir()
    reactions = copy.deepcopy(REACTIONS)
    del reactions['T_sec']
    reactions['DM_sec'] = ({'sec_c': 1}, {}, 0, 1000, None)
    model = out / 'toy.xml'
    write_sbml(model, reactions=reactions)
    doc = libsbml.readSBMLFromFile(str(model))
    for rid, ec in [('SEC', '4.4.4.4'), ('ACS', '6.2.1.1')]:
        r = doc.getModel().getReaction(rid); r.setMetaId('meta_' + rid)
        cv = libsbml.CVTerm(libsbml.BIOLOGICAL_QUALIFIER)
        cv.setBiologicalQualifierType(libsbml.BQB_IS)
        cv.addResource('https://identifiers.org/ec-code/' + ec); r.addCVTerm(cv)
    libsbml.writeSBMLToFile(doc, str(model))
    claims = out / 'claims.jsonl'
    claims.write_text(''.join(json.dumps({'protein_id': 'locus_g6', 'function_type': 'EC',
        'function_id': 'EC:' + ec, 'source': 'diamond', 'evidence_type': 'SEQUENCE_SIMILARITY',
        'raw_score': .6, 'calibrated_prob': .6,
        'metadata': {'product_demand': 'DM_sec', 'substrate_exchange': 'EX_ac_e', 'substrate_role': 'carbon'}}) + '\n'
        for ec in ('4.4.4.4', '6.2.1.1')))
    summary = []
    for kind in ('standard', 'masked'):
        template = manifest({}, product_requirement='maintenance_only', max_focus_reactions=0,
                            max_baseline_focus_reactions=0)
        template['model'] = str(model)
        template['baseline']['product_demands'] = ['DM_sec']
        if kind == 'masked':
            template['assignment_recovery'] = {'mask_gene': 'g6', 'candidate_protein': 'locus_g6'}
        path = out / (kind + '-template.json'); save(path, template)
        result_path = out / (kind + '.tsv')
        command = ['java', '-jar', str(jar), 'integrate', '--claims', str(claims), '--out', str(result_path),
            '--lite', '--enable-priors', 'none', '--experiment-manifest', str(path),
            '--experiment-mode', 'evidence_guided', '--experiment-strategy', 'fixed',
            '--experiment-max-solver-calls', '200', '--experiment-max-experiments', '30',
            '--experiment-python-bin', sys.executable, '--experiment-sidecar', str(ROOT / 'scripts/metabolic_experiments.py')]
        save(out / (kind + '-command.json'), command)
        with (out / (kind + '.log')).open('w') as log:
            process = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
        assert process.returncode == 0, (kind, process.returncode)
        folder = Path(str(result_path) + '.experiments')
        result = json.loads((folder / 'result.json').read_text())
        generated = json.loads((folder / 'manifest.json').read_text())
        assert result['promotion'] == 'none'
        if kind == 'standard':
            products = [r for r in result['experiments'] if r['experiment_id'].startswith('synthesis:')]
            assert products and any(r['measurements'].get('product', {}).get('max_synthesis', 0) > 0 for r in products)
            assert all('max_secretion' not in r['measurements'].get('product', {}) for r in products)
            summary.append({'kind': kind, 'synthesis_experiments': len(products), 'solver_calls': result['solver']['used']})
        else:
            report = CO.compare(model, generated, result); save(out / 'masked-oracle.json', report)
            assert report['n_comparisons'] > 0 and report['n_agree'] == report['n_comparisons']
            assert report['n_cell_agree'] == report['n_cell_comparisons']
            products = [c for c in result['cells'] if 'product_contrast' in c]
            assert {c['product_contrast'] for c in products} >= {'synthesis_capacity_lost', 'synthesis_capacity_retained'}
            assert any(c.get('product_demand') == 'DM_sec' and c['mods'] for c in result['public_problem']['conditions'])
            summary.append({'kind': kind, 'synthesis_cells': len(products), 'oracle_comparisons': report['n_comparisons']})
    return summary


def real_model(out, model_path, medium, live):
    out.mkdir()
    model = ME.load_sbml_model(model_path)
    original = cobra.io.read_sbml_model(model_path, f_replace={}); original.solver = 'glpk'
    base = GD.manifest(model, medium, 'published_model_minimal')
    base['baseline']['product_demands'] = GD.PRODUCTS
    base['config'].update(product_requirement='maintenance_only', max_baseline_focus_reactions=0,
                          max_experiments=40, max_solver_calls=500, strategy='fixed')
    base['context'] = {'hypotheses': [{'id': 'H', 'gene': GD.GENE}], 'evidence': [
        {'id': 'E' + str(i), 'type': 'bgc', 'hypothesis_id': 'H', 'product_demand': demand,
         'condition_hint': 'supplement_removal', 'exchange': 'R_EX_pi_e'} for i, demand in enumerate(GD.PRODUCTS)]}
    save(out / 'standard-manifest.json', base)
    result = ME.run_experiments(base, model); save(out / 'standard-result.json', result)
    independent = {}
    for condition in ('published_model_minimal', 'phosphate_closed'):
        for genotype in ('native', 'gcs_deleted'):
            m = GD.manifest(model, medium, condition)
            independent[(condition, genotype)] = GD.oracle(original, m, genotype, 'maintenance_only')
    comparisons = []
    for row in result['experiments']:
        if not row['experiment_id'].startswith('synthesis:'):
            continue
        condition = 'phosphate_closed' if any(m['op'] == 'close' for m in row['mods']) else 'published_model_minimal'
        genotype = 'gcs_deleted' if any(m['op'] == 'knockout' for m in row['mods']) else 'native'
        product = row['measurements']['product']; rid = product['demand']
        expected = independent[(condition, genotype)]['capacities'][rid]
        observed = {'minimum': product['min_synthesis_at_required_state'], 'maximum': product['max_synthesis'],
                    'statuses': [product.get('lp_status', {}).get(sense, 'optimal') for sense in ('min', 'max')]}
        comparisons.append({'experiment': row['experiment_id'], 'agrees': GD.agree(expected, observed),
                            'expected': expected, 'observed': observed})
    save(out / 'standard-oracle.json', comparisons)
    assert len(comparisons) == 12 and all(r['agrees'] for r in comparisons)
    # Known positive and decoy hypotheses are explicit controls, never inferred
    # from hidden target assignments by a planner or scored as discovery.
    masked = copy.deepcopy(base)
    masked['context'] = {}
    masked['assignment_recovery'] = {'mask_gene': GD.GENE,
        'hypotheses': [{'id': 'H0', 'assignments': {}},
            {'id': 'H1', 'assignments': {r: [[MA.TARGET]] for r in ('R_GERMA', 'R_GERMB', 'R_GERMC', 'R_GERMD')}},
            {'id': 'H2', 'assignments': {'R_GERMA': [[MA.TARGET]]}}],
        'conditions': [{'id': 'C' + str(i * 3 + j), 'mods': mods, 'product_demand': demand, 'evidence_channels': ['bgc']}
            for i, demand in enumerate(GD.PRODUCTS)
            for j, mods in enumerate(([], [{'op': 'close', 'exchange': 'R_EX_pi_e'}],
                                      [{'op': 'close', 'exchange': 'R_EX_glc__D_e'}]))]}
    save(out / 'masked-manifest.json', masked)
    result = ME.run_experiments(masked, model); save(out / 'masked-result.json', result)
    report = CO.compare(model_path, masked, result); save(out / 'masked-oracle.json', report)
    assert len(result['cells']) == 36 and report['n_agree'] == report['n_comparisons'] == 72
    assert report['n_cell_agree'] == report['n_cell_comparisons'] == 36
    for cell in result['cells']:
        cid = int(cell['condition_id'][1:])
        if cid % 3 == 2:
            assert cell['status'] == 'inconclusive' and cell['product_contrast'] == 'inconclusive'
        elif cell['hypothesis_id'] == 'H1':
            assert cell['status'] == 'completed' and cell['product_contrast'] == 'synthesis_capacity_lost'
    summary = {'standard_product_comparisons': len(comparisons), 'masked_cells': len(result['cells']),
               'masked_oracle_comparisons': report['n_comparisons'], 'live': []}
    if live:
        for kind, full in [('standard', base), ('masked', masked)]:
            m = copy.deepcopy(full)
            m['config'].update(strategy='qwen', max_experiments=3, max_solver_calls=120)
            raw_chat = UsageChat(ME.resolve_config(m))
            class Chat:
                usage = None
                def __call__(self, messages):
                    reply = raw_chat(messages); self.usage = raw_chat.usage(); return reply
            result = ME.run_experiments(m, model, chat=Chat())
            save(out / (kind + '-qwen-manifest.json'), m)
            save(out / (kind + '-qwen-result.json'), result)
            save(out / (kind + '-qwen-transcript.json'), raw_chat.transcript)
            save(out / (kind + '-qwen-requests.json'), raw_chat.log)
            if kind == 'masked':
                oracle = CO.compare(model_path, m, result); save(out / (kind + '-qwen-oracle.json'), oracle)
                assert oracle['n_agree'] == oracle['n_comparisons']
                assert oracle['n_cell_agree'] == oracle['n_cell_comparisons']
                valid = all(r['status'] == 'agent_selected' for r in result['selection'])
            else:
                valid = all(r.get('valid') and r['source'] == 'qwen' for r in result['selection_log'])
                for entry, transcript in zip(result['selection_log'], raw_chat.transcript):
                    assert entry['prompt_sha256'] == ME.canonical_hash(transcript['messages'])
            summary['live'].append({'kind': kind, 'valid': valid, 'usage': raw_chat.usage()})
            assert 1 <= raw_chat.usage()['ok'] <= 3 and raw_chat.usage()['failed'] == 0 and valid
    return summary


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('out', 'model', 'medium'):
        p.add_argument('--' + name, type=Path, required=True)
    p.add_argument('--jar', type=Path)
    p.add_argument('--live-qwen', action='store_true')
    p.add_argument('--key-file', type=Path)
    args = p.parse_args(); args.out.mkdir(parents=True, exist_ok=False)
    if args.key_file:
        os.environ['BORG_LLM_API_KEY'] = args.key_file.read_text().strip()
    summary = {'scope': 'mechanism and integration diagnostics; no biological accuracy or agent benefit claim'}
    if args.jar:
        summary['cli'] = cli_fixture(args.out / 'cli', args.jar.resolve())
    summary['sco_gem'] = real_model(args.out / 'sco-gem', args.model, args.medium, args.live_qwen)
    save(args.out / 'summary.json', summary); print(json.dumps(summary))


if __name__ == '__main__':
    main()
