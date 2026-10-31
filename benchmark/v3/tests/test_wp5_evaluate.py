"""WP5 step 5 (design §6, §7, §8.5; ruling §13.7): decision truth table, denominators, transport harm,
partial-growth exclusion, policy gating and the calibrate-before-test order."""
import copy
import json
import os
import shutil
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import wp5_fixtures as F  # noqa: E402
import prepare_wp5_cases as P5  # noqa: E402
import run_wp5_recovery as R  # noqa: E402
import evaluate_wp5_recovery as E  # noqa: E402

HYPS = [{'id': 'H0', 'assignments': {}},
        {'id': 'H1', 'assignments': {'R_A': [['MASKED_TARGET']]}, 'evidence': {'sequence': 0.9}},
        {'id': 'H2', 'assignments': {'R_B': [['MASKED_TARGET']]}, 'evidence': {'sequence': 0.6}}]
CONDS = [{'id': 'C0', 'mods': [], 'evidence_channels': []},
         {'id': 'C1', 'mods': [{'op': 'close', 'exchange': 'R_EX_a_e'}], 'evidence_channels': []},
         {'id': 'C2', 'mods': [{'op': 'close', 'exchange': 'R_EX_b_e'}], 'evidence_channels': []}]
KO = {('H1', 'C1'): 'knockout_no_growth', ('H2', 'C2'): 'knockout_infeasible_with_maintenance'}


def public_problem(hyps=HYPS):
    return {'baseline': {'exchanges': ['R_EX_a_e', 'R_EX_b_e'], 'medium': {'R_EX_a_e': 10.0}},
            'objectives': {'biomass': 'R_BIO', 'maintenance': 'none'}, 'model': 'masked:' + '0' * 64,
            'hypotheses': copy.deepcopy(hyps), 'conditions': copy.deepcopy(CONDS)}


def cell(h, c, contrast=None, assignments=None):
    contrast = contrast or KO.get((h, c), 'both_grow')
    rxns = assignments if assignments is not None else next(x['assignments'] for x in HYPS if x['id'] == h)
    return {'hypothesis_id': h, 'condition_id': c, 'contrast': contrast,
            'status': 'inconclusive' if contrast in ('inconclusive', 'no_growing_control') else 'completed',
            'phenotypes': {'wild_type': {'status': 'optimal', 'growth': 0.8, 'knocked_out_reactions': [],
                                         'flux_intervals': {r: {'minimum': 0.4, 'maximum': 1.0} for r in rxns}},
                           'target_knockout': {'status': 'optimal', 'growth': 0.0, 'knocked_out_reactions': [],
                                               'flux_intervals': {}}},
            'solver_calls': []}


def run_of(visited, status='completed', contrasts=None):
    contrasts = contrasts or {}
    return {'status': status, 'visited_conditions': list(visited),
            'cells': [cell(h['id'], c, contrasts.get((h['id'], c))) for c in visited for h in HYPS]}


def label(status, stable=True, wt='growth'):
    return {'mutant_status': status, 'stable': stable, 'wt_status': wt}


TRUTH_A = {'reactions': ['R_A'], 'stratum': 'catalytic'}
TRUTH_B = {'reactions': ['R_B'], 'stratum': 'complex_subunit'}


class DecisionTruthTable(unittest.TestCase):
    def decide(self, measured, truth=TRUTH_A, visited=('C0', 'C1', 'C2'), run=None, hyps=HYPS):
        return E.decide(public_problem(hyps), run or run_of(visited), truth, measured)

    def test_correct_call(self):
        d = self.decide({'C1': 'no_growth', 'C2': 'growth'})
        # C1 and C2 both discriminate and agree with H1.
        self.assertEqual((d['decision'], d['called'], d['supporting']), ('correct_call', 'H1', ['C1', 'C2']))
        self.assertEqual(d['score'], 0.9)
        d = self.decide({'C1': 'no_growth'}, visited=('C0', 'C1'))
        self.assertEqual((d['decision'], d['supporting']), ('correct_call', ['C1']))
        self.assertTrue(d['truth_in_candidates'] and d['discriminating_visited'])

    def test_wrong_call_and_partial_overlap(self):
        d = self.decide({'C1': 'no_growth', 'C2': 'growth'}, truth={'reactions': ['R_C'], 'stratum': 'catalytic'})
        self.assertEqual(d['decision'], 'wrong_call')
        self.assertFalse(d['partial_overlap'] or d['truth_in_candidates'])
        d = self.decide({'C1': 'no_growth', 'C2': 'growth'}, truth={'reactions': ['R_A', 'R_Z'], 'stratum': 'catalytic'})
        self.assertEqual(d['decision'], 'wrong_call')
        self.assertTrue(d['partial_overlap'])      # reported separately, never correct

    def test_any_call_on_a_transport_only_target_is_harmful(self):
        for reactions in (['R_T'], ['R_A']):
            d = self.decide({'C1': 'no_growth'}, truth={'reactions': reactions, 'stratum': 'transport_only'})
            self.assertEqual(d['decision'], 'wrong_call')
            self.assertTrue(d['transport_catalytic_call'])

    def test_null_call_correct_negative_and_false_null(self):
        growth_everywhere = {'C1': 'growth', 'C2': 'growth'}
        d = self.decide(growth_everywhere, truth={'reactions': ['R_T'], 'stratum': 'transport_only'})
        self.assertEqual((d['decision'], d['null_outcome'], d['called']), ('null_call', 'correct_negative', 'H0'))
        d = self.decide(growth_everywhere, truth={'reactions': ['R_C'], 'stratum': 'catalytic'})
        self.assertEqual(d['null_outcome'], 'correct_negative')         # truth has no candidate
        d = self.decide(growth_everywhere, truth=TRUTH_A)
        self.assertEqual(d['null_outcome'], 'false_null')

    def test_abstain_tied(self):
        tied = copy.deepcopy(HYPS)
        tied[2]['evidence'] = {'sequence': 0.9}
        d = self.decide({'C1': 'growth', 'C2': 'growth'}, hyps=[tied[0], tied[1], tied[2],
                        {'id': 'H3', 'assignments': {'R_C': [['MASKED_TARGET']]}, 'evidence': {'sequence': 0.9}}],
                        run={'status': 'completed', 'visited_conditions': ['C0'],
                             'cells': [cell(h, 'C0', assignments={}) for h in ('H0', 'H1', 'H2', 'H3')]})
        self.assertEqual(d['decision'], 'abstain_tied')

    def test_abstain_unsupported_without_a_discriminating_agreement(self):
        d = self.decide({'C0': 'growth'}, visited=('C0',))
        self.assertEqual(d['decision'], 'abstain_unsupported')
        self.assertFalse(d['discriminating_visited'])
        # Agreement at a condition where every hypothesis predicts the same phenotype is not support.
        d = self.decide({'C1': 'no_growth'}, visited=('C0', 'C1'), run=run_of(('C0', 'C1'),
                        contrasts={('H0', 'C1'): 'knockout_no_growth', ('H2', 'C1'): 'knockout_no_growth'}))
        self.assertEqual(d['decision'], 'abstain_unsupported')

    def test_no_consistent_hypothesis(self):
        d = self.decide({'C1': 'no_growth', 'C2': 'no_growth'})
        self.assertEqual(d['decision'], 'no_consistent_hypothesis')
        self.assertEqual(d['contradicted'], ['H0', 'H1', 'H2'])

    def test_inconclusive_for_missing_failed_or_solver_inconclusive_runs(self):
        measured = {'C1': 'no_growth', 'C2': 'growth'}
        self.assertEqual(E.decide(public_problem(), None, TRUTH_A, measured)['reason'], 'missing_run')
        for status in ('agent_failed', 'cell_mismatch', 'invalid_budget', 'transcript_unaudited'):
            d = self.decide(measured, run=run_of(('C0', 'C1', 'C2'), status=status))
            self.assertEqual((d['decision'], d['reason']), ('inconclusive', status))
        d = self.decide(measured, run=run_of(('C0', 'C1'), contrasts={('H2', 'C0'): 'inconclusive'}))
        self.assertEqual((d['decision'], d['reason']), ('inconclusive', 'solver_or_budget_inconclusive'))
        broken = run_of(('C0', 'C1'))
        broken['cells'].pop()
        d = self.decide(measured, run=broken)
        self.assertEqual(d['decision'], 'inconclusive')
        self.assertIn('rectangle', d['reason'])
        outside = run_of(('C0', 'C1'))
        outside['visited_conditions'].append('C9')
        self.assertEqual(self.decide(measured, run=outside)['decision'], 'inconclusive')

    def test_unconfirmed_labels_are_never_negatives(self):
        for status in ('partial', 'unknown', 'unmeasured', 'wt_unverified'):
            d = self.decide({'C1': status, 'C2': status})
            self.assertEqual(d['decision'], 'abstain_unsupported', status)
            self.assertEqual(d['contradicted'], [], status)
        # Partial growth at C1 does not contradict H1's no-growth prediction; C2 alone decides.
        d = self.decide({'C1': 'partial', 'C2': 'no_growth'}, truth=TRUTH_B)
        self.assertEqual((d['decision'], d['called'], d['supporting']), ('correct_call', 'H2', ['C2']))
        self.assertEqual([p['measured'] for p in d['pairs']], ['unmeasured', 'partial', 'no_growth'])

    def test_no_growing_control_cells_are_not_evaluable(self):
        d = self.decide({'C1': 'no_growth', 'C2': 'growth'},
                        run=run_of(('C0', 'C1', 'C2'), contrasts={('H1', 'C2'): 'no_growing_control'}))
        self.assertEqual(d['decision'], 'correct_call')

    def test_measured_class_mapping(self):
        self.assertEqual(E.measured_class(None), 'unmeasured')
        self.assertEqual(E.measured_class(label('growth')), 'growth')
        self.assertEqual(E.measured_class(label('no_growth', wt='stated_by_source')), 'no_growth')
        self.assertEqual(E.measured_class(label('no_growth', stable=False)), 'unknown')
        self.assertEqual(E.measured_class(label('partial')), 'partial')
        self.assertEqual(E.measured_class(label('weird')), 'unknown')
        for wt in ('no_growth', 'unknown', None):
            self.assertEqual(E.measured_class(label('no_growth', wt=wt)), 'wt_unverified')

    def test_sequence_only_reference_and_policy_gate(self):
        d = E.sequence_only(public_problem(), TRUTH_A)
        self.assertEqual((d['decision'], d['score'], d['visited']), ('correct_call', 0.9, []))
        row = dict(d, eligible=True)
        self.assertEqual(E.gate(row, 0.95)['decision'], 'abstain_policy')
        self.assertEqual(E.gate(row, None)['decision'], 'abstain_policy')   # no threshold: the arm abstains
        self.assertEqual(E.gate(row, 0.9)['decision'], 'correct_call')
        null = {'decision': 'null_call', 'null_outcome': 'false_null'}
        self.assertEqual(E.gate(null, None), null)

    def test_capacity_fva_classes_of_called_and_truth_assignments(self):
        d = self.decide({'C1': 'no_growth', 'C2': 'growth'})
        cap = E.capacity(public_problem(), run_of(('C0', 'C1', 'C2')), TRUTH_A, d)
        self.assertEqual(cap, {'called': {'required': 3}, 'truth': {'required': 3}})
        self.assertEqual(E.fva_class({'minimum': 0.2, 'maximum': 1.0}), 'required')
        self.assertEqual(E.fva_class({'minimum': 0.0, 'maximum': 0.0}), 'no_flux_at_floor')
        self.assertEqual(E.fva_class(None), 'unassessed')


def row(case, decision, score=0.9, truth=True, eligible=True, primary=True, arm='evidence', k=3, seed=0, **extra):
    r = {'case': case, 'problem': case + '-v', 'decision': decision, 'score': score, 'truth_in_candidates': truth,
         'eligible': eligible, 'primary_variant': primary, 'arm': arm, 'k': k, 'seed': seed, 'stratum': 'catalytic',
         'species': 'Ec', 'source': 'keio', 'pairs': [], 'n_hypotheses': 3}
    if decision == 'null_call':
        r.setdefault('null_outcome', extra.pop('null_outcome', 'false_null'))
    r.update(extra)
    return r


class Denominators(unittest.TestCase):
    def test_recovery_harm_abstention_and_false_null_use_all_eligible_cases(self):
        rows = [row('A', 'correct_call'), row('B', 'wrong_call'), row('C', 'null_call', null_outcome='false_null'),
                row('D', 'no_admissible_experiment', truth=False), row('E', 'inconclusive'),
                row('F', 'abstain_unsupported'), row('G', 'correct_call', score=0.5)]
        s = E.summarize(rows, 0.9)
        self.assertEqual(s['denominator'], 7)
        self.assertAlmostEqual(s['recovery'], 1 / 7)
        self.assertAlmostEqual(s['harm'], 1 / 7)
        self.assertAlmostEqual(s['false_null'], 1 / 7)
        self.assertAlmostEqual(s['abstention'], 2 / 7)          # unsupported + the gated call
        self.assertEqual(s['counts']['inconclusive'], 1)
        self.assertEqual(s['counts']['no_admissible_experiment'], 1)
        self.assertAlmostEqual(s['candidate_recall'], 6 / 7)
        self.assertEqual(s['conditional_denominator'], 6)
        self.assertAlmostEqual(s['conditional_recovery'], 1 / 6)
        self.assertEqual(s['counts']['abstain_policy'], 1)
        self.assertEqual(sum(s['counts'].values()), 7)

    def test_random_seeds_are_averaged_per_case(self):
        rows = [row('A', 'correct_call' if s < 3 else 'abstain_unsupported', arm='random', seed=s) for s in range(10)]
        rows += [row('B', 'wrong_call', arm='random', seed=s) for s in range(10)]
        s = E.summarize(rows, 0.0)
        self.assertEqual(s['denominator'], 2)
        self.assertAlmostEqual(s['recovery'], 0.15)
        self.assertAlmostEqual(s['harm'], 0.5)

    def test_transport_control_harm_is_counted(self):
        rows = [row('T', 'wrong_call', stratum='transport_only', transport_catalytic_call=True),
                row('U', 'null_call', null_outcome='correct_negative', stratum='transport_only')]
        s = E.summarize(rows, 0.0)
        self.assertEqual(s['transport_catalytic_calls'], 1)
        self.assertEqual(s['counts']['correct_negative'], 1)
        self.assertEqual(s['harm'], 0.5)

    def test_calibration_rule(self):
        items = [row(f'c{i}', 'correct_call', score=0.9) for i in range(8)] + [row('w', 'wrong_call', score=0.5)]
        self.assertEqual(E.calibrate_rows(items)['evidence:k3']['status'], 'no_supported_threshold_abstain')  # 9 < 10
        items.append(row('c9', 'correct_call', score=0.8))
        policy = E.calibrate_rows(items)['evidence:k3']
        self.assertEqual((policy['status'], policy['threshold']), ('selected_on_calibration', 0.5))   # 9/10 >= 0.9
        items.append(row('w2', 'wrong_call', score=0.85))
        policy = E.calibrate_rows(items)['evidence:k3']
        self.assertEqual(policy['threshold'], 0.8)        # 0.5 gives 9/11 < 0.9; 0.8 gives 9/10 with ten decisions
        self.assertEqual(policy['calibration_counts']['wrong_call'], 1)
        ineligible = [dict(r, eligible=False) for r in items]
        self.assertEqual(E.calibrate_rows(ineligible), {})

    def test_phenotype_agreement_excludes_partial_growth_from_binary_pairs(self):
        pairs = [{'condition_id': 'C1', 'measured': 'no_growth', 'predictions': {'H0': 'growth', 'H1': 'no_growth'}},
                 {'condition_id': 'C2', 'measured': 'partial', 'predictions': {'H0': 'growth', 'H1': 'growth'}},
                 {'condition_id': 'C3', 'measured': 'unknown', 'predictions': {'H0': 'growth', 'H1': 'growth'}}]
        r = row('A', 'correct_call', called='H1', pairs=pairs)
        a = E.phenotype_agreement([r], 0.0, [row('A', 'abstain_unsupported', arm='exhaustive', pairs=pairs)])
        self.assertEqual(a['binary_pairs'], 1)
        self.assertEqual(a['partial_pairs'], 1)
        self.assertEqual(a['not_evaluable_pairs'], 1)
        self.assertEqual((a['called_pairs'], a['called_agree']), (1, 1))
        self.assertEqual(a['constant_no_growth_agree'], 1)
        self.assertEqual(a['constant_growth_agree'], 0)
        self.assertEqual(a['common_exhaustive']['binary_pairs'], 1)
        self.assertEqual(a['common_exhaustive']['partial_pairs'], 1)
        self.assertEqual(a['unmasked_original']['status'], 'not_available')

    def test_statistics(self):
        self.assertEqual(E.mcnemar({'a': True, 'b': True}, {'a': True, 'b': True})['p'], 1.0)
        m = E.mcnemar({str(i): True for i in range(6)}, {str(i): False for i in range(6)})
        self.assertEqual((m['only_a'], m['only_b']), (6, 0))
        self.assertAlmostEqual(m['p'], 2 / 64)
        lo, hi = E.wilson(5, 10)
        self.assertTrue(0.2 < lo < 0.5 < hi < 0.8)
        self.assertEqual(E.wilson(0, 0), [None, None])
        self.assertEqual(E.bootstrap_difference({'a': 1, 'b': 1}, {'a': 0, 'b': 0}, n=200), [1.0, 1.0])


# --------------------------------------------------------------------------- end-to-end synthetic world

SCENARIOS = {
    # name: (truth reactions, stratum, {menu index: mutant status}, wt status)
    'A': (['R_A'], 'catalytic', {1: 'no_growth', 2: 'growth'}, 'stated_by_source'),
    'W': (['R_C'], 'catalytic', {1: 'no_growth', 2: 'growth'}, 'stated_by_source'),
    'B': (['R_B'], 'complex_subunit', {1: 'growth', 2: 'no_growth'}, 'growth'),
    'T': (['R_T'], 'transport_only', {1: 'no_growth', 2: 'growth'}, 'growth'),
    'N': (['R_T'], 'transport_only', {1: 'growth', 2: 'growth'}, 'growth'),
    'P': (['R_B'], 'complex_subunit', {1: 'partial', 2: 'no_growth'}, 'growth'),
    'X': (['R_A', 'R_B', 'R_C'], 'ambiguous_truth', {1: 'no_growth'}, 'growth'),
}


class World:
    """A bundle and run directory written directly as JSON (shape of the real pipeline outputs)."""

    def __init__(self, root, calibration='A' * 11 + 'W', test='ABTNPX', unconstructed_test=1, qwen_audited=(3,)):
        self.root = Path(root)
        self.bundle, self.runs = self.root / 'bundle', self.root / 'runs'
        pub = self.bundle / 'public'
        cases, truth, labels, problems = [], {'calibration': [], 'test': []}, {'calibration': [], 'test': []}, {}
        specs = [('calibration', s, 'keio', ['manuscript_regulatory']) for s in calibration]
        specs += [('test', s, 'nutrient', ['printed', 'additional_inorganic']) for s in test]
        specs += [('test', 'A', 'sulfonate', ['m63_sulfur_ceiling_1'])] * unconstructed_test
        for n, (split, scenario, source, variants) in enumerate(specs, 1):
            case = f'W{n:03d}'
            reactions, stratum, measured, wt = SCENARIOS[scenario]
            truth[split].append({'case': case, 'reactions': reactions, 'stratum': stratum})
            for variant in variants:
                name = f'{case}-{variant}'
                primary = variant == variants[0]
                cases.append({'case': case, 'problem': name, 'split': split, 'source': source, 'species': 'Ec',
                              'variant': variant, 'primary_variant': primary})
                constructed = source != 'sulfonate'
                record = {'case': case, 'problem': name, 'status': 'constructed' if constructed else
                          'no_admissible_experiment', 'masked_model': f'models/{case}.xml', 'masked_model_sha256': 'f' * 64}
                if constructed:
                    problem = public_problem()
                    record.update(public_problem=problem, public_problem_sha256=R.canonical_hash(problem),
                                  menu_map={'0': 'C0', '1': 'C1', '2': 'C2'})
                problems[name] = record
                P5.write_json(pub / 'problems' / f'{name}.json', record)
                for index, status in measured.items():
                    labels[split].append({'case': case, 'variant': variant, 'menu_index': index, 'source': source,
                                          'mutant_status': status, 'stable': True, 'wt_status': wt})
        P5.write_jsonl(pub / 'cases.jsonl', cases)
        for split in ('calibration', 'test'):
            P5.write_jsonl(self.bundle / f'private/truth.{split}.jsonl', truth[split])
            P5.write_jsonl(self.bundle / f'private/labels.{split}.jsonl', labels[split])
        F.pass_audit(self.bundle)
        entries = {}
        for name, record in problems.items():
            entry = {'status': record['status'], 'masked_model_sha256': record['masked_model_sha256']}
            if record['status'] == 'constructed':
                entry.update(public_problem_sha256=record['public_problem_sha256'], cost_per_condition=10,
                             n_conditions=3, panels={'minimal_medium': ['C0'], 'fixed_panel': ['C0', 'C1', 'C2']},
                             caps={str(k): {'max_experiments': k, 'max_solver_calls': 10 * k} for k in R.BUDGETS})
            entries[name] = entry
        self.protocol = {'schema': R.SCHEMA + ':protocol', 'public_manifest_sha256': P5.sha_file(pub / 'manifest.json'),
                         'masking_audit_sha256': P5.sha_file(pub / 'masking_audit.json'), 'arms': list(R.ARMS),
                         'budgets': list(R.BUDGETS), 'random_seeds': list(R.RANDOM_SEEDS), 'problems': entries}
        P5.write_json(self.runs / 'protocol.json', self.protocol)
        self.protocol_sha = P5.sha_file(self.runs / 'protocol.json')
        transcripts = {}
        for stage in ('exhaustive', 'deterministic', 'qwen'):
            for name, arm, k, seed in R.plan_runs(self.protocol, stage):
                visited = (['C0', 'C1', 'C2'] if arm == 'exhaustive' else
                           ['C1', 'C2', 'C0'][:k] if arm == 'qwen' else
                           R.planned_selection(problems[name]['public_problem'], arm, k, seed,
                                               entries[name]['panels'].get(arm)))
                run = self.write_run(name, arm, k, seed, visited)
                if arm == 'qwen':
                    tname = R.transcript_name(name, arm, k, seed)
                    P5.write_jsonl(self.runs / tname, [{'request': [{'role': 'user', 'content': '{}'}]}])
                    run.update(transcript=tname, transcript_sha256=P5.sha_file(self.runs / tname))
                    P5.write_json(self.runs / 'runs' / R.run_name(name, arm, k, seed), run)
                    if k in qwen_audited:
                        transcripts[tname] = run['transcript_sha256']
        P5.write_json(self.runs / 'transcript_audit.json', {'status': 'pass', 'protocol_sha256': self.protocol_sha,
                                                            'transcripts_sha256': transcripts})
        self.out = self.root / 'evaluation'

    def write_run(self, name, arm, k, seed, visited, **extra):
        n, cap = (3, 30) if arm == 'exhaustive' else (k, 10 * k)
        run = dict({'schema': R.SCHEMA + ':run', 'problem': name, 'arm': arm, 'k': k, 'seed': seed,
                    'public_problem_sha256': self.protocol['problems'][name]['public_problem_sha256'],
                    'protocol_sha256': self.protocol_sha, 'engine_entry': R.ENGINE_ENTRY, 'status': 'completed',
                    'max_experiments': n, 'max_solver_calls': cap,
                    'config': {'max_experiments': n, 'max_solver_calls': cap},
                    'visited_conditions': visited, 'cells': run_of(visited)['cells'],
                    'solver_calls': [{}] * (10 * len(visited)), 'solver_calls_used': 10 * len(visited),
                    'solver_calls_unused': cap - 10 * len(visited), 'wall_seconds': 0.1,
                    'agent_usage': {'requests': len(visited), 'failures': 0, 'seconds': 1.0} if arm == 'qwen' else None},
                   **extra)
        P5.write_json(self.runs / 'runs' / R.run_name(name, arm, k, seed), run)
        return run

    def private_test_files(self):
        return {(self.bundle / f'private/{k}.test.jsonl').resolve() for k in ('truth', 'labels')}


class EndToEnd(unittest.TestCase):
    def setUp(self):
        self.root = Path(tempfile.mkdtemp(prefix='wp5-eval-', dir=os.environ.get('TMPDIR')))
        self.addCleanup(shutil.rmtree, self.root, True)
        self.world = World(self.root)
        self.w = self.world

    def freeze(self):
        guard = self.w.private_test_files()
        with F.refuse_paths(lambda p: p in guard) as opened:
            E.join(self.w.bundle, self.w.runs, self.w.out)
            sha = E.calibrate(self.w.bundle, self.w.runs, self.w.out)
        self.assertFalse(guard & set(opened))      # test truth and labels are never opened before the policy
        return sha

    def test_join_and_calibrate_never_open_test_files_and_select_a_threshold(self):
        sha = self.freeze()
        policy = P5.read_json(self.w.out / 'policy.json')
        self.assertEqual(sha, P5.sha_file(self.w.out / 'policy.json'))
        self.assertEqual(policy['arms']['fixed_panel:k3']['status'], 'selected_on_calibration')
        self.assertEqual(policy['arms']['fixed_panel:k3']['threshold'], 0.9)
        self.assertEqual(policy['arms']['minimal_medium:k3']['status'], 'no_supported_threshold_abstain')
        self.assertEqual(policy['arms']['sequence_only']['status'], 'selected_on_calibration')
        self.assertIn('design', policy['note'])
        rows = P5.read_jsonl(self.w.out / 'joined.calibration.jsonl')
        self.assertEqual({r['split'] for r in rows}, {'calibration'})

    def test_evaluate_refuses_a_policy_whose_hash_or_bindings_do_not_match(self):
        sha = self.freeze()
        with self.assertRaisesRegex(P5.WP5Error, 'policy hash mismatch'):
            E.evaluate(self.w.bundle, self.w.runs, self.w.out, '0' * 64)
        for bad in (None, 'ABC', sha.upper()):
            with self.assertRaises(P5.WP5Error):
                E.evaluate(self.w.bundle, self.w.runs, self.w.out, bad)
        policy = P5.read_json(self.w.out / 'policy.json')
        P5.write_json(self.w.out / 'policy.json', dict(policy, arms={}))
        with self.assertRaisesRegex(P5.WP5Error, 'policy hash mismatch'):
            E.evaluate(self.w.bundle, self.w.runs, self.w.out, sha)
        P5.write_json(self.w.out / 'policy.json', policy)
        self.assertEqual(P5.sha_file(self.w.out / 'policy.json'), sha)
        name = 'W001-manuscript_regulatory'
        self.w.write_run(name, 'evidence', 3, 0, ['C0'])                # a run changed after the freeze
        with self.assertRaisesRegex(P5.WP5Error, 'run files changed'):
            E.evaluate(self.w.bundle, self.w.runs, self.w.out, sha)
        self.assertFalse((self.w.out / 'joined.test.jsonl').exists())

    def test_changed_calibration_labels_after_the_freeze_are_refused(self):
        sha = self.freeze()
        path = self.w.bundle / 'private/labels.calibration.jsonl'
        path.write_text(path.read_text() + '\n')
        with self.assertRaisesRegex(P5.WP5Error, 'calibration truth or labels changed'):
            E.evaluate(self.w.bundle, self.w.runs, self.w.out, sha)

    def test_evaluate_tables_denominators_transport_harm_and_report(self):
        sha = self.freeze()
        result = E.evaluate(self.w.bundle, self.w.runs, self.w.out, sha)
        t = result['tables']['fixed_panel:k3|primary']
        # Six constructed test cases (A B T N P X) plus one unconstructed; X is ambiguous (reported only).
        self.assertEqual(t['denominator'], 6)
        self.assertEqual(t['reported_only'], {'ambiguous_truth': 1})
        c = t['counts']
        self.assertEqual(c['correct_call'], 1)                  # A (score 0.9 passes the 0.9 threshold)
        self.assertEqual(c['abstain_policy'], 2)                # B and P are correct at score 0.6 but gated
        self.assertEqual(c['wrong_call'], 1)                    # T: catalytic call on a transport control
        self.assertEqual(c['correct_negative'], 1)              # N: null on a transport control
        self.assertEqual(c['no_admissible_experiment'], 1)
        self.assertEqual(t['transport_catalytic_calls'], 1)
        self.assertAlmostEqual(t['recovery'], 1 / 6)
        self.assertAlmostEqual(t['candidate_recall'], 3 / 6)   # A, B, P have truth among candidates
        self.assertEqual(t['strata']['stratum:transport_only']['denominator'], 2)
        self.assertEqual(t['phenotype_agreement']['partial_pairs'], 1)
        mm = result['tables']['minimal_medium:k3|primary']
        self.assertEqual(mm['policy_status'], 'no_supported_threshold_abstain')
        self.assertEqual(mm['recovery'], 0.0)
        sens = result['tables']['fixed_panel:k3|sensitivity']
        self.assertEqual(sens['denominator'], 5)               # nutrient sensitivity arm only
        # Qwen k=3 transcripts were audited, k=6 were not: unaudited runs are inconclusive, not credited.
        self.assertEqual(result['tables']['qwen:k6|primary']['counts']['inconclusive'], 5)
        self.assertEqual(result['tables']['qwen:k6|primary']['recovery'], 0.0)
        self.assertGreater(result['tables']['qwen:k3|primary']['counts']['correct_call'], 0)
        comparators = {(r['arm'], r['comparator']) for r in result['tests']}
        self.assertTrue({('qwen:k3', 'evidence:k3'), ('qwen:k3', 'minimal_medium:k3'), ('qwen:k3', 'sequence_only'),
                         ('evidence:k3', 'minimal_medium:k3')} <= comparators)
        random_rows = [r for r in result['tests'] if 'random' in r['arm'] + r['comparator']]
        self.assertTrue(random_rows and all(not r.get('mcnemar_exact') for r in random_rows))
        E.report(self.w.out)
        text = (self.w.out / 'report.md').read_text()
        self.assertGreaterEqual(text.count(E.KEIO_NOTE), 2)
        self.assertIn('oxygen uptake 20', text)
        self.assertIn('Calibration policy', text)
        with self.assertRaisesRegex(P5.WP5Error, 'immutable'):
            E.evaluate(self.w.bundle, self.w.runs, self.w.out, sha)

    def test_invalid_runs_are_kept_as_inconclusive_failures(self):
        name = 'W013-printed'      # a test case of scenario A
        self.w.write_run(name, 'evidence', 3, 0, ['C0', 'C1', 'C2'], max_solver_calls=31)
        self.w.write_run(name, 'evidence', 6, 0, ['C0', 'C1', 'C2'], public_problem_sha256='1' * 64)
        self.w.write_run(name, 'fixed_panel', 3, 0, ['C0', 'C1', 'C2'], engine_entry='my_own_loop')
        run = self.w.write_run(name, 'fixed_panel', 6, 0, ['C0', 'C1', 'C2'])
        run['cells'][4]['contrast'] = 'both_grow'                # differs from the exhaustive cell
        P5.write_json(self.w.runs / 'runs' / R.run_name(name, 'fixed_panel', 6, 0), run)
        self.w.write_run(name, 'minimal_medium', 3, 0, ['C0'], status='agent_failed')
        self.w.write_run(name, 'minimal_medium', 6, 0, ['C0'], solver_calls_used=0)
        index, _ = E.run_index(self.w.runs, self.w.protocol)
        self.assertEqual(index[(name, 'evidence', 3, 0)]['status'], 'invalid_budget')
        self.assertEqual(index[(name, 'evidence', 6, 0)]['status'], 'invalid_provenance')
        self.assertEqual(index[(name, 'fixed_panel', 3, 0)]['status'], 'invalid_engine')
        self.assertEqual(index[(name, 'fixed_panel', 6, 0)]['status'], 'cell_mismatch')
        self.assertEqual(index[(name, 'minimal_medium', 6, 0)]['status'], 'invalid_solver_ledger')
        sha = self.freeze()
        rows = {(r['problem'], r['arm'], r['k']): r for r in
                E.evaluate(self.w.bundle, self.w.runs, self.w.out, sha) and P5.read_jsonl(self.w.out / 'joined.test.jsonl')}
        for key in (('evidence', 3), ('evidence', 6), ('fixed_panel', 3), ('fixed_panel', 6), ('minimal_medium', 3),
                    ('minimal_medium', 6)):
            self.assertEqual(rows[(name,) + key]['decision'], 'inconclusive', key)

    def test_missing_truth_row_fails_closed(self):
        path = self.w.bundle / 'private/truth.calibration.jsonl'
        P5.write_jsonl(path, P5.read_jsonl(path)[1:])
        with self.assertRaisesRegex(P5.WP5Error, 'truth file does not cover'):
            E.join(self.w.bundle, self.w.runs, self.w.out)

    def test_failed_masking_audit_invalidates_the_evaluation(self):
        audit = self.w.bundle / 'public/masking_audit.json'
        P5.write_json(audit, dict(P5.read_json(audit), status='fail'))
        with self.assertRaisesRegex(P5.WP5Error, 'masking audit'):
            E.join(self.w.bundle, self.w.runs, self.w.out)

    def test_cli_requires_the_policy_hash_for_evaluate(self):
        with self.assertRaises(SystemExit):
            E.main(['evaluate', '--bundle', str(self.w.bundle), '--runs', str(self.w.runs), '--out', str(self.w.out)])


if __name__ == '__main__':
    unittest.main()
