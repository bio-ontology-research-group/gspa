"""Learned-stage report in the results browser: schema-2 ledger, `integrated.tsv.learned.json`, trust checks.

Fixtures have the artifact shape the JVM writers produce (`IntegrationRunAudit.report`,
`IntegrationArtifacts.writeLearnedReport`); a contract test reads the Groovy sources so the fixture
cannot drift from them. Probabilities here are synthetic: nothing asserts or implies accuracy.
"""
import copy
import json
import re
import unittest
from pathlib import Path

from scripts.tests import test_visualize_ledger as base

ROOT = Path(__file__).resolve().parents[2]
SHA = 'a' * 64
OTHER_SHA = 'b' * 64
RUN_ORDER = ['reference_preparation', 'cohort_context', 'evidence_integration', 'dark_matter', 'reaction_context', 'cross_genome',
             'symbolic_enforcement', 'metabolic_refinement', 'agent_completion', 'metabolic_experiments', 'community_metabolism', 'learned_integration',
             'feature_export']


def v2_ledger(mode='learned', sha=SHA, order=None):
    """Schema-2 ledger in the writer's stage order; the learned stage rescored one enzyme (.85 -> .9)."""
    ledger = base.ledger_fixture()
    by = {s['name']: s for s in ledger['stages']}
    learned_change = [base.change('p_hexokinase', 'EC', 'EC:2.7.1.1', 'increased', .85, .9)]
    stages = []
    for name in order or RUN_ORDER:
        s = by.get(name) or base.stage(name, 'not_run')
        if name in ('reference_preparation', 'cohort_context'): s = base.stage(name, 'executed')
        if name == 'learned_integration' and mode == 'learned':
            s = base.stage(name, 'executed', learned_change, model_id='deepgo-test-1', model_sha256=sha, model_kind='logistic',
                           contract_id='gspa.features/1', candidates_scored=6, note='Scores were replaced by the frozen model.')
        if name == 'feature_export': s = base.stage(name, 'executed', output_dir='features')
        stages.append(s)
    ledger['stages'][:] = stages
    ledger['final_changes'] = ledger['final_changes'] + (learned_change if mode == 'learned' else [])
    ledger.update(schema_version=2, scoring_mode=mode,
                  scoring_mode_note='Annotation scores were replaced by a frozen learned model (see learned_integration).' if mode == 'learned'
                  else 'Legacy heuristic scoring (evidence combiner and hand-set priors). No trained model produced these scores.',
                  baseline_definitions={'stage0_raw_max_score': 'Stage 0: strongest raw predictor score per function.',
                                        'legacy_prior_baseline': 'Posterior after the evidence_integration stage.'},
                  changes_measured_from='legacy_prior_baseline')
    return ledger


def cand(key, role, p, margin=0., projected=True, locus=None, null=None, sources=('diamond',), notes=()):
    return {'candidate_key': f'{key}#{role}', 'function_key': key, 'role': role, 'margin': margin, 'probability': p,
            'above_display_threshold': p >= .5, 'projected': projected, 'feature_digest': 'd' * 64, 'locus_probability': locus,
            'locus_null_probability': null, 'notes': list(notes), 'supporting_sources': list(sources)}


def learned_report(sha=SHA, locus_reason='locus layer disabled for this run'):
    hexo, tran, ghost, hyp = 'p_hexokinase|EC|EC:2.7.1.1', 'p_transporter|GO|GO:0005215', 'p_ghost|EC|EC:9.9.9.9', 'p_hypothetical|EC|EC:3.1.1.1'
    return {'schema': 'gspa.learned-stage-report/1',
            'model': {'id': 'deepgo-test-1', 'sha256': sha, 'kind': 'logistic', 'intended_use': 'unit test', 'training': {}},
            'contract': {'id': 'gspa.features/1', 'fingerprint': 'f' * 16},
            'settings': {'retract_below': .1, 'add_context_candidates': True, 'add_unretained_claims': False,
                         'compute_locus_assignments': False, 'structurally_forbidden_go_terms': 0},
            'locus_layer_unavailable_reason': locus_reason,
            'counts': {'candidates': 5, 'rescored': 1, 'added': 0, 'retracted': 0, 'removed_no_projection': 0, 'scored_only': 2},
            'pinned_floors_cleared': [],
            'outcomes': [
                {'function_key': hexo, 'action': 'rescored', 'reason': None, 'projected_role': 'CATALYTIC', 'prior_probability': .85,
                 'learned_probability': .9, 'role_probabilities': {'CATALYTIC': .9, 'ACCESSORY': .3}},
                {'function_key': tran, 'action': 'scored_only', 'reason': 'role is not an enzyme', 'projected_role': None,
                 'prior_probability': .7, 'learned_probability': None, 'role_probabilities': {'TRANSPORTER': .8}},
                {'function_key': hyp, 'action': 'scored_only', 'reason': 'below display threshold', 'projected_role': None,
                 'prior_probability': None, 'learned_probability': None, 'role_probabilities': {'CATALYTIC': .35}},
                {'function_key': ghost, 'action': 'scored_only', 'reason': None, 'projected_role': None, 'prior_probability': None,
                 'learned_probability': None, 'role_probabilities': {'CATALYTIC': .5}}],
            'candidates': [
                cand(hexo, 'CATALYTIC', .9, 2.2, True, locus=.6, null=.2),
                cand(hexo, 'ACCESSORY', .3, -.8, False, locus=.1, null=.2, sources=('diamond', 'interproscan')),
                cand(tran, 'TRANSPORTER', .8, 1.4, False),
                cand(hyp, 'CATALYTIC', .35, -.6, False, sources=(), notes=('context-only candidate',)),
                cand(ghost, 'CATALYTIC', .5, 0., False)],
            'interpretation': 'Scores replace the heuristic posterior; they are calibrated marginals of a frozen model, not a measurement of accuracy.'}


class LearnedPayloadTest(base.LedgerTest):
    def learned_of(self, **kw):
        data = self.render(ledger=kw.pop('ledger', v2_ledger()), learned=kw.pop('learned', learned_report()), **kw)
        return data['ledger']

    def test_stage_order_is_the_ledger_order_with_readable_labels(self):
        ledger = self.learned_of()
        names = [s['name'] for s in ledger['stages']]
        self.assertEqual(names, RUN_ORDER)
        labels = {s['name']: s['label'] for s in ledger['stages']}
        self.assertEqual(labels['cohort_context'], 'Cohort context (two-pass)')
        self.assertEqual(labels['learned_integration'], 'Learned integration (frozen model)')
        self.assertEqual(labels['feature_export'], 'Feature export')
        self.assertEqual(labels['metabolic_experiments'], 'Metabolic experiments')
        shuffled = ['learned_integration'] + [n for n in RUN_ORDER if n != 'learned_integration']
        self.assertEqual([s['name'] for s in self.learned_of(ledger=v2_ledger(order=shuffled))['stages']], shuffled)

    def test_prediction_free_stages_are_not_reported_as_zero_changes(self):
        stages = {s['name']: s for s in self.learned_of()['stages']}
        for name in ('cohort_context', 'feature_export', 'reference_preparation'):
            self.assertTrue(stages[name]['no_prediction_changes'], name)
            self.assertNotIn('summary', stages[name])
        self.assertEqual(stages['learned_integration']['summary']['functions'], 1)

    def test_scoring_mode_learned_legacy_and_unknown(self):
        ledger = self.learned_of()
        self.assertEqual(ledger['scoring']['mode'], 'learned')
        self.assertTrue(ledger['scoring']['consistent'])
        self.assertIn('not a measurement of accuracy', ledger['scoring']['claim'])
        legacy = self.render(ledger=v2_ledger(mode='legacy_heuristic'))['ledger']
        self.assertEqual(legacy['scoring']['mode'], 'legacy_heuristic')
        self.assertEqual(legacy['learned']['status'], 'not_applicable')
        old = self.render()['ledger']          # schema 1: the ledger does not say
        self.assertEqual(old['scoring']['mode'], 'unknown')
        self.assertIsNone(old['scoring']['consistent'])
        self.assertEqual(old['learned']['status'], 'unavailable')
        self.assertNotIn('candidates', old['learned'])

    def test_mode_claiming_learned_without_a_learned_stage_is_inconsistent(self):
        ledger = v2_ledger()
        next(s for s in ledger['stages'] if s['name'] == 'learned_integration').update(status='failed', reason='boom')
        self.assertFalse(self.learned_of(ledger=ledger)['scoring']['consistent'])

    def test_baseline_uses_changes_measured_from_and_old_schema_falls_back(self):
        ledger = self.learned_of()
        self.assertEqual(ledger['baseline']['kind'], 'legacy_prior_baseline')
        self.assertFalse(ledger['baseline']['sequence_only'])
        self.assertIn('not a sequence-only baseline', ledger['baseline']['label'])
        self.assertIn('evidence_integration', ledger['baseline']['definition'])
        self.assertEqual(self.render()['ledger']['baseline']['kind'], 'after_evidence_integration')

    def test_candidates_join_by_canonical_key_and_marginals_stay_separate(self):
        learned = self.learned_of()['learned']
        self.assertEqual(learned['status'], 'available')
        self.assertTrue(learned['hash_verified'])
        hexo = [c for c in learned['candidates'] if c['k'] == 'p_hexokinase|EC|EC:2.7.1.1']
        self.assertEqual({c['role']: c['pr'] for c in hexo}, {'CATALYTIC': .9, 'ACCESSORY': .3})
        self.assertEqual({c['role']: c['lp'] for c in hexo}, {'CATALYTIC': .6, 'ACCESSORY': .1})
        self.assertTrue(all(c['joined'] for c in hexo))
        # No derived probability exists anywhere in the payload: nothing sums the marginals.
        self.assertNotIn(1.2, [c.get('pr') for c in learned['candidates']])
        self.assertEqual(learned['outcomes']['p_hexokinase|EC|EC:2.7.1.1']['roles'], {'CATALYTIC': .9, 'ACCESSORY': .3})

    def test_null_locus_probability_stays_null_with_the_reason(self):
        learned = self.learned_of()['learned']
        tran = next(c for c in learned['candidates'] if c['role'] == 'TRANSPORTER')
        self.assertIsNone(tran['lp'])
        self.assertIsNone(tran['ln'])
        self.assertEqual(learned['locus_unavailable_reason'], 'locus layer disabled for this run')

    def test_unknown_protein_is_listed_as_not_joined(self):
        learned = self.learned_of()['learned']
        self.assertEqual(learned['unjoined_proteins'], [{'p': 'p_ghost', 'n': 1}])
        self.assertEqual(learned['n_unjoined_candidates'], 1)
        ghost = next(c for c in learned['candidates'] if c['p'] == 'p_ghost')
        self.assertFalse(ghost['joined'])

    def test_protein_ids_containing_pipes_are_split_from_the_right(self):
        report = learned_report()
        report['candidates'].append(cand('sp|P00956|SYI_ECOLI|GO|GO:0005215', 'TRANSPORTER', .4))
        learned = self.learned_of(learned=report)['learned']
        c = next(c for c in learned['candidates'] if c['f'] == 'GO:0005215' and c['p'].startswith('sp|'))
        self.assertEqual((c['p'], c['t']), ('sp|P00956|SYI_ECOLI', 'GO'))

    def test_non_enzyme_roles_and_unprojected_candidates_are_kept(self):
        learned = self.learned_of()['learned']
        roles = {(c['k'], c['role']): c['proj'] for c in learned['candidates']}
        self.assertFalse(roles[('p_transporter|GO|GO:0005215', 'TRANSPORTER')])
        self.assertFalse(roles[('p_hexokinase|EC|EC:2.7.1.1', 'ACCESSORY')])     # rejected role is still scored
        self.assertTrue(roles[('p_hexokinase|EC|EC:2.7.1.1', 'CATALYTIC')])
        hyp = next(c for c in learned['candidates'] if c['f'] == 'EC:3.1.1.1')   # scored, never in the final EC table
        self.assertEqual((hyp['proj'], hyp['src']), (False, []))
        self.assertEqual(learned['outcomes'][hyp['k']]['action'], 'scored_only')

    def test_model_hash_mismatch_is_not_trusted_and_withholds_scores(self):
        learned = self.learned_of(learned=learned_report(sha=OTHER_SHA))['learned']
        self.assertEqual(learned['status'], 'not_trusted')
        self.assertIn('mismatch', learned['reason'])
        self.assertNotIn('candidates', learned)

    def test_missing_stage_hash_is_flagged_unverified_not_trusted_silently(self):
        ledger = v2_ledger()
        next(s for s in ledger['stages'] if s['name'] == 'learned_integration').pop('model_sha256')
        learned = self.learned_of(ledger=ledger)['learned']
        self.assertEqual(learned['status'], 'available')
        self.assertFalse(learned['hash_verified'])
        self.assertIn('could not be checked', learned['hash_note'])

    def test_missing_stage_hash_is_explicitly_unverified(self):
        ledger = v2_ledger()
        next(s for s in ledger['stages'] if s['name'] == 'learned_integration').pop('model_sha256')
        learned = self.learned_of(ledger=ledger)['learned']
        self.assertIs(learned['hash_verified'], False)
        self.assertIn('unverified', learned['hash_note'])
        self.assertNotIn('equals', learned['hash_note'])

    def test_failed_or_missing_learned_stage_withholds_scores_even_in_learned_mode(self):
        failed = v2_ledger()
        next(s for s in failed['stages'] if s['name'] == 'learned_integration').update(status='failed', reason='boom')
        missing = v2_ledger()
        missing['stages'][:] = [s for s in missing['stages'] if s['name'] != 'learned_integration']
        for ledger in (failed, missing):
            out = self.learned_of(ledger=ledger)
            self.assertFalse(out['scoring']['learned_ran'])
            self.assertEqual(out['learned']['status'], 'not_trusted')
            self.assertIn('did not execute', out['learned']['reason'])
            for withheld in ('candidates', 'outcomes', 'hash_verified'):
                self.assertNotIn(withheld, out['learned'])

    def test_report_hash_absent_while_the_ledger_expects_one_is_not_trusted(self):
        report = learned_report()
        del report['model']['sha256']
        learned = self.learned_of(learned=report)['learned']
        self.assertEqual(learned['status'], 'not_trusted')
        self.assertIn('no hash', learned['reason'])
        self.assertNotIn('candidates', learned)

    def test_wrong_shaped_report_is_unavailable_and_bad_scalars_do_not_crash(self):
        good = learned_report()
        wrong = [[], 3, None, {**good, 'model': 'deepgo'}, {**good, 'model': None}, {**good, 'schema': 'other/1'}, {**good, 'schema': None},
                 {**good, 'candidates': {}}, {**good, 'outcomes': 'x'}]
        for bad in wrong:
            learned = self.learned_of(learned=json.dumps(bad))['learned']
            self.assertEqual(learned['status'], 'unavailable', bad)
            self.assertNotIn('candidates', learned)
        odd = learned_report()
        odd['model'].update(id=7, kind=['x'])
        odd.update(counts='many', contract=[], settings=3, interpretation=5)
        odd['candidates'][0].update(probability='high', margin=float('nan'), locus_probability=True, supporting_sources='diamond', notes=3, role=4)
        odd['candidates'].append('not a candidate')
        odd['outcomes'][0].update(role_probabilities=[1], prior_probability='x', action=1)
        learned = self.learned_of(learned=json.dumps(odd))['learned']
        self.assertEqual(learned['status'], 'available')
        first = next(c for c in learned['candidates'] if c['k'] == 'p_hexokinase|EC|EC:2.7.1.1' and c['ck'].endswith('#CATALYTIC'))
        self.assertEqual((first['pr'], first['m'], first['lp'], first['src'], first['notes'], first['role']), (None, None, None, [], [], None))
        self.assertEqual((learned['model']['id'], learned['counts'], learned['contract']), (None, None, None))

    def test_absent_unreadable_or_contradicted_report_is_unavailable_not_zero(self):
        absent = self.render(ledger=v2_ledger())['ledger']['learned']
        self.assertEqual(absent['status'], 'unavailable')
        self.assertIn('not zero', absent['reason'])
        self.assertNotIn('counts', absent)
        for bad in ('{not json', '{"candidates": 3}', '{"schema": "gspa.learned-stage-report/1", "candidates": 3}'):
            self.assertEqual(self.learned_of(learned=bad)['learned']['status'], 'unavailable', bad)
        stray = self.render(ledger=v2_ledger(mode='legacy_heuristic'), learned=learned_report())['ledger']['learned']
        self.assertEqual(stray['status'], 'not_trusted')


class FixtureContractTest(unittest.TestCase):
    """The fixtures use exactly the field names the JVM writers emit."""
    def source(self, rel):
        return (ROOT / 'gspa-core/src/main/groovy/gspa/integration' / rel).read_text()

    def test_learned_report_fields_match_the_groovy_writer(self):
        text = self.source('workflow/IntegrationArtifacts.groovy')
        block = text[text.index('static File writeLearnedReport'):]
        report = learned_report()
        for key in list(report) + list(report['candidates'][0]) + list(report['outcomes'][0]) + list(report['model']) + list(report['counts']):
            self.assertRegex(block, rf'\b{re.escape(key)}\s*:' if key != 'schema' else r'schema\s*:', key)
        self.assertIn("'gspa.learned-stage-report/1'", block)

    def test_ledger_fields_match_the_groovy_audit(self):
        text = self.source('IntegrationRunAudit.groovy')
        ledger = v2_ledger()
        for key in ('schema_version', 'scoring_mode', 'scoring_mode_note', 'baseline_definitions', 'changes_measured_from',
                    'baseline_candidates', 'final_candidates', 'final_changes'):
            self.assertIn(key, ledger)
            self.assertIn(key + ':', text)
        self.assertIn('schema_version: 2', text)
        stages = re.search(r"STAGES = \[(.*?)\]\.asImmutable", text, re.S).group(1)
        self.assertEqual(re.findall(r"'(\w+)'", stages), RUN_ORDER)
        flow = self.source('workflow/IntegrationWorkflow.groovy')
        for key in ('model_id', 'model_sha256', 'model_kind'):
            self.assertIn(key + ':', flow)


if __name__ == '__main__':
    unittest.main()
