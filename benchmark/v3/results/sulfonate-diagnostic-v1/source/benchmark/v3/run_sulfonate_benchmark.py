#!/usr/bin/env python3
"""Complete known-genotype sulfonate diagnostic, with unsupported cases retained.

No training, gene assignment recovery, or original-study reproduction is claimed.
Preparation/run never read observations. Evaluation joins them only afterward.
"""
import argparse
from collections import Counter
import json
from pathlib import Path
import sys
import time

import run_nutrient_phenotypes as N
import prepare_metabolic_benchmarks as P
ME = N.ME
ROOT = Path(__file__).resolve().parents[2]
SULFUR = {'Sulfate': 'so4', 'Taurine': 'taur', 'Isethionate': 'isetac',
          'Ethanesulfonate': 'ethso3', 'Butanesulfonate': 'butso3', 'Sulfoacetate': 'sulfac',
          'MOPS': 'mops', 'HEPES': 'hepes', 'PIPES': 'pipes'}
PROTOCOL = {'schema': 'gspa.sulfonate_diagnostic/1', 'sulfur_uptake_ceilings': [0.1, 1.0, 10.0],
            'oxygen_uptake': 20.0, 'growth_threshold': 1e-6,
            'base_rates': {'glc__D': 10., 'nh4': 10., 'trp__L': .1, 'thm': .01,
                           **{x: 1000. for x in ('ca2', 'cl', 'cobalt2', 'cu2', 'fe2', 'fe3', 'h', 'h2o',
                                                 'k', 'mg2', 'mn2', 'mobd', 'na1', 'ni2', 'pi', 'zn2')}},
            'assumptions': ['M63 carbon source is glucose according to the cited 1996 methods; 2000 specifies Trp supplementation.',
                            'Flux ceilings are declared modeling assumptions, not conversion of measured concentrations.',
                            'Trace ion availability is a model diagnostic assumption; the complete sulfur-free trace recipe is unresolved.',
                            'Native maintenance and all internal reactions/GPRs remain unchanged.',
                            'In-frame experimental deletions are approximated by loss of the indicated gene product.',
                            'MG1655 model is not strain-matched EC1250; background trp-1 lesion is not mapped.',
                            'Wild-type column defines normalization; sulfate is a control, not a recovered function.',
                            'iML1515 may encode this literature; this is not an independent or masked recovery benchmark.']}


def implementation():
    paths = [Path(__file__), Path(N.__file__), Path(P.__file__), Path(ME.__file__),
             Path(__file__).with_name('check_nutrient_cobra.py')]
    return {str(p.relative_to(ROOT)): N.sha(p) for p in paths}


def plan_case(case, model, index, sulfur_rate):
    reasons = []
    if case['genotype_status'] != 'table_reported' or case['deleted_genes'] is None:
        reasons.append('genotype_uncertain')
    token = SULFUR.get(case['sulfur_source'])
    exchange = 'R_EX_' + token + '_e' if token else None
    if exchange not in model.rxn_index:
        reasons.append('sulfur_exchange_unrepresented')
    genes, mappings = [], []
    for symbol in case['deleted_genes'] or []:
        locus, method = P.map_gene(symbol, index['symbol_index'])
        canonical = model.gene_alias.get(locus)
        mappings.append({'symbol': symbol, 'locus': locus, 'canonical': canonical, 'method': method})
        if not canonical:
            reasons.append('unmapped_gene:' + symbol)
        else:
            genes.append(canonical)
    base = {'case_id': case['case_id'], 'arm': f'sulfur_uptake={sulfur_rate:g}',
            'sulfur_source': case['sulfur_source'], 'genotype_label': case['genotype_label'],
            'status': 'unsupported' if reasons else 'ready', 'reasons': reasons,
            'sulfur_exchange': exchange, 'genes': genes, 'gene_mappings': mappings}
    if reasons:
        return base
    medium = {'R_EX_' + k + '_e': v for k, v in PROTOCOL['base_rates'].items()}
    missing = set(medium) - set(model.rxn_ids)
    if missing:
        raise ValueError('base exchange missing: ' + repr(sorted(missing)))
    medium[exchange] = sulfur_rate
    medium['R_EX_o2_e'] = PROTOCOL['oxygen_uptake']
    if len(model.objective_terms) != 1 or model.objective_sense != 'maximize':
        raise ValueError('requires one maximized biomass objective')
    manifest = {'model': model.path, 'model_sha256_expected': model.sha256,
                'objectives': {'biomass': model.objective_terms[0][0],
                               'maintenance': {'reaction': 'R_ATPM', 'lower_bound': float(model.lb[model.rxn_index['R_ATPM']])}},
                'baseline': {'exchanges': sorted(medium), 'medium': medium},
                'context': {'hypotheses': [{'id': 'observed_genotype', 'genes': genes}] if genes else []},
                'config': {'strategy': 'fixed', 'max_experiments': 1, 'max_solver_calls': 2,
                           'growth_threshold': PROTOCOL['growth_threshold'],
                           'max_focus_reactions': 0, 'max_baseline_focus_reactions': 0}}
    ME.validate_manifest(manifest)
    ME.validate_boundary(model, manifest)
    base['manifest'] = manifest
    return base


def prepare(inputs, model_path, out):
    cases = N.read_rows(inputs)
    if len({r['case_id'] for r in cases}) != len(cases):
        raise ValueError('duplicate case identity')
    model, index = ME.load_sbml_model(model_path), P.index_sbml(model_path)
    plans = [plan_case(case, model, index, rate) for rate in PROTOCOL['sulfur_uptake_ceilings'] for case in cases]
    out.mkdir(parents=True, exist_ok=False)
    N.write_json(out / 'protocol.json', {'protocol': PROTOCOL, 'implementation': implementation(),
                  'inputs_sha256': N.sha(inputs), 'model_sha256': model.sha256,
                  'cases': len(cases), 'plans': len(plans),
                  'statuses': dict(Counter(p['status'] for p in plans)),
                  'plans_sha256': N.digest(plans)})
    (out / 'plans.jsonl').write_text(''.join(json.dumps(p, sort_keys=True) + '\n' for p in plans))
    print(json.dumps({'plans': len(plans), 'statuses': dict(Counter(p['status'] for p in plans))}), flush=True)


def verify(plan, model_path=None):
    frozen = json.loads((plan / 'protocol.json').read_text())
    plans = N.read_rows(plan / 'plans.jsonl')
    if (frozen['protocol'] != PROTOCOL or frozen['implementation'] != implementation()
            or N.digest(plans) != frozen['plans_sha256']):
        raise ValueError('frozen code, protocol or plans changed')
    if model_path and N.sha(model_path) != frozen['model_sha256']:
        raise ValueError('model changed')
    return frozen, plans


def run(plan, model_path, out):
    frozen, plans = verify(plan, model_path)
    model = ME.load_sbml_model(model_path)
    out.mkdir(parents=True, exist_ok=False)
    start = time.time()
    with (out / 'predictions.jsonl').open('w') as stream:
        for i, entry in enumerate(plans):
            if entry['status'] == 'unsupported':
                row = {'case_id': entry['case_id'], 'arm': entry['arm'], 'plan_sha256': N.digest(entry),
                       'prediction': None, 'status': 'unsupported', 'reasons': entry['reasons']}
            else:
                try:
                    row = N.simulate(entry, model)
                    row['status'] = 'executed'
                except Exception as exc:
                    row = {'case_id': entry['case_id'], 'arm': entry['arm'], 'plan_sha256': N.digest(entry),
                           'prediction': None, 'status': 'failed', 'error': f'{type(exc).__name__}: {exc}'}
            stream.write(json.dumps(row, sort_keys=True, allow_nan=False) + '\n')
            if (i+1) % 100 == 0:
                stream.flush();print(json.dumps({'completed': i+1, 'elapsed_s': time.time()-start}), flush=True)
    N.write_json(out / 'run.json', {'model_sha256': model.sha256, 'protocol_sha256': N.sha(plan / 'protocol.json'),
                                   'predictions_sha256': N.sha(out / 'predictions.jsonl'), 'elapsed_s': time.time()-start})


def verify_run(plan, run):
    frozen, plans = verify(plan)
    receipt = json.loads((run / 'run.json').read_text())
    if (receipt['model_sha256'] != frozen['model_sha256']
            or receipt['protocol_sha256'] != N.sha(plan / 'protocol.json')
            or receipt['predictions_sha256'] != N.sha(run / 'predictions.jsonl')):
        raise ValueError('run receipt mismatch')
    predictions = N.read_rows(run / 'predictions.jsonl')
    by_id = {(p['case_id'], p['arm']): p for p in plans}
    if len(predictions) != len(by_id) or {(p['case_id'], p['arm']) for p in predictions} != set(by_id):
        raise ValueError('missing/duplicated predictions')
    for p in predictions:
        if p['plan_sha256'] != N.digest(by_id[p['case_id'], p['arm']]):
            raise ValueError('prediction-plan mismatch')
    return frozen, by_id, predictions


def oracle(plan, run, model_path, out):
    import cobra
    import check_nutrient_cobra as O
    frozen, plans, rows = verify_run(plan, run)
    if N.sha(model_path) != frozen['model_sha256']:
        raise ValueError('oracle model mismatch')
    model = cobra.io.read_sbml_model(model_path, f_replace={});model.solver = 'glpk'
    results = []
    for row in rows:
        entry = plans[row['case_id'], row['arm']]
        checks = []
        for record in row.get('production', {}).get('experiments', []):
            genes = [] if record['kind'] == 'baseline' else entry['genes']
            other = O.simulate(model, entry['manifest'], genes)
            checks.append({'kind': record['kind'], 'oracle': other, 'flux_agrees': O.agreement(record, other),
                           'gpr_agrees': record['measurements'].get('knocked_out_reactions') == other['blocked_reactions']})
        expected = (2 if entry['genes'] else 1) if entry['status'] == 'ready' else 0
        results.append({'case_id': row['case_id'], 'arm': row['arm'], 'checks': checks,
                        'status': entry['status'],
                        'agrees': (len(checks) == expected and all(c['flux_agrees'] and c['gpr_agrees'] for c in checks))
                                  if entry['status'] == 'ready' else None})
    report = {'schema': 'gspa.sulfonate_oracle/1', 'run_sha256': N.sha(run / 'run.json'),
              'model_sha256': N.sha(model_path), 'solver': 'COBRApy/GLPK; independently parsed SBML/GPRs',
              'comparison_policy': '1e-6 absolute/relative; feasibility 1e-9; standard basis before each LP',
              'comparisons': sum(len(r['checks']) for r in results),
              'unsupported': sum(r['status'] == 'unsupported' for r in results),
              'disagreements': sum(r['agrees'] is False for r in results),
              'results': results}
    N.write_json(out, report)
    return report


def evaluate(plan, run, observations, out):
    _, plans, rows = verify_run(plan, run)
    obs = N.read_rows(observations);truth = {r['case_id']: r for r in obs}
    if len(truth) != len(obs) or set(truth) != {r['case_id'] for r in rows}:
        raise ValueError('observation catalogue mismatch')
    cross = Counter();records = []
    for row in rows:
        outcome = truth[row['case_id']]
        cross[(row['arm'], outcome['observation_role'], outcome['symbol'] or 'unreported', str(row['prediction']))] += 1
        records.append({'case_id': row['case_id'], 'arm': row['arm'], 'observed_symbol': outcome['symbol'],
                        'role': outcome['observation_role'], 'prediction': row['prediction'],
                        'prediction_status': row['status'], 'genotype_status': outcome['genotype_status']})
    N.write_json(out, {'schema': 'gspa.sulfonate_phenotype_diagnostic/1',
                      'observations_sha256': N.sha(observations), 'run_sha256': N.sha(run / 'run.json'),
                      'cases_by_arm': dict(Counter(r['arm'] for r in records)),
                      'cross_tabulation': [{'arm': k[0], 'role': k[1], 'observed_symbol': k[2],
                                            'predicted_growth': k[3], 'count': v} for k, v in sorted(cross.items())],
                      'records': records,
                      'limits': ['All controls, intermediate categories, uncertain genotypes and unrepresented substrates retained.',
                                 'FBA capacity is not normalized endpoint yield; no invented quantitative ordinal accuracy.',
                                 'Known model GPRs retained: no function recovery or independent biological acceptance.']})


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest='command', required=True)
    for command, arguments in [('prepare', ['inputs', 'model', 'out']), ('run', ['plan', 'model', 'out']),
                                ('oracle', ['plan', 'run', 'model', 'out']), ('evaluate', ['plan', 'run', 'observations', 'out'])]:
        s = sub.add_parser(command)
        for name in arguments:s.add_argument('--' + name, type=Path, required=True)
    a = p.parse_args()
    if a.out.exists():p.error('output exists')
    if a.command == 'prepare':prepare(a.inputs, a.model, a.out)
    elif a.command == 'run':run(a.plan, a.model, a.out)
    elif a.command == 'oracle':
        report = oracle(a.plan, a.run, a.model, a.out)
        print(json.dumps({k: report[k] for k in ('comparisons', 'disagreements')}))
        if report['disagreements']:sys.exit(1)
    else:evaluate(a.plan, a.run, a.observations, a.out)


if __name__ == '__main__':main()
