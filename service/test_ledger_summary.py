"""Job-level ledger summary: no fastapi needed. Absent ledger is unavailable, never zero changes."""
import json
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import unittest

from service import ledger_summary

WEB = Path(__file__).parent / 'web/index.html'


def change(pid, ftype, fid, kind, before, after):
    return {'function_key': f'{pid}|{ftype}|{fid}', 'protein_id': pid, 'function_type': ftype, 'function_id': fid,
            'change': kind, 'before_probability': before, 'after_probability': after}


LEDGER = {'schema_version': 1, 'baseline_candidates': 3, 'final_candidates': 3,
          'stages': [{'name': 'evidence_integration', 'status': 'executed', 'changes': [change('a', 'GO', 'GO:1', 'added', None, .9)]},
                     {'name': 'dark_matter', 'status': 'executed', 'changes': []},
                     {'name': 'cross_genome', 'status': 'skipped', 'reason': 'disabled', 'changes': []},
                     {'name': 'symbolic_enforcement', 'status': 'failed', 'reason': 'checker unavailable', 'changes': []},
                     {'name': 'agent_completion', 'status': 'not_run', 'changes': []}],
          'final_changes': [change('p1', 'EC', 'EC:1.1.1.1', 'removed', .7, None),
                            change('p1', 'GO', 'GO:2', 'increased', .2, .6),
                            change('p2', 'GO', 'GO:2', 'increased', .2, .6)]}


class LedgerSummaryTest(unittest.TestCase):
    def summarize(self, ledger):
        with tempfile.TemporaryDirectory() as tmp:
            if ledger is not None: (Path(tmp) / ledger_summary.LEDGER_NAME).write_text(ledger if isinstance(ledger, str) else json.dumps(ledger))
            return ledger_summary.summarize(tmp)

    def test_absent_ledger_is_unavailable_not_zero(self):
        result = self.summarize(None)
        self.assertFalse(result['available'])
        self.assertNotIn('net', result)
        self.assertIn('not zero', result['reason'])

    def test_unreadable_ledger_is_unavailable(self):
        for bad in ('not json', '{"stages": 3}', '{}'):
            self.assertFalse(self.summarize(bad)['available'], bad)

    def test_counts_are_unique_functions_and_proteins_from_ledger_actions(self):
        result = self.summarize(LEDGER)
        self.assertEqual(result['net']['functions'], 3)
        self.assertEqual(result['net']['proteins'], 2)
        self.assertEqual(result['net']['by_change'], {'added': 0, 'removed': 1, 'increased': 2, 'decreased': 0})
        self.assertEqual(result['status_counts'], {'executed': 2, 'skipped': 1, 'failed': 1, 'not_run': 1, 'unknown': 0})
        self.assertEqual([s['reason'] for s in result['stages'] if s['status'] == 'failed'], ['checker unavailable'])

    def test_baseline_is_not_called_sequence_only_unless_the_ledger_says_so(self):
        self.assertIn('not a sequence-only baseline', self.summarize(LEDGER)['baseline'])
        explicit = {**LEDGER, 'baseline_kind': 'sequence_only'}
        self.assertEqual(self.summarize(explicit)['baseline'], 'Sequence-only evidence')


SHA = 'a' * 64


def learned_ledger(mode='learned', sha=SHA):
    stages = [{'name': n, 'status': 'executed', 'changes': []} for n in ('reference_preparation', 'cohort_context')]
    stages += [{'name': 'learned_integration', 'status': 'executed' if mode == 'learned' else 'not_run', 'changes': [],
                **({'model_sha256': sha} if sha else {})}, {'name': 'feature_export', 'status': 'skipped', 'changes': []}]
    return {'schema_version': 2, 'scoring_mode': mode, 'changes_measured_from': 'legacy_prior_baseline', 'stages': stages,
            'baseline_candidates': 1, 'final_candidates': 1, 'final_changes': []}


class LearnedSummaryTest(LedgerSummaryTest):
    def summarize_with(self, ledger, report=None):
        with tempfile.TemporaryDirectory() as tmp:
            (Path(tmp) / ledger_summary.LEDGER_NAME).write_text(json.dumps(ledger))
            if report is not None: (Path(tmp) / ledger_summary.LEARNED_NAME).write_text(report if isinstance(report, str) else json.dumps(report))
            return ledger_summary.summarize(tmp)

    def report(self, sha=SHA):
        return {'model': {'id': 'm1', 'sha256': sha}, 'candidates': [], 'outcomes': [], 'counts': {'candidates': 0}}

    def test_stage_order_labels_and_baseline_come_from_the_new_schema(self):
        out = self.summarize_with(learned_ledger(), self.report())
        self.assertEqual([s['name'] for s in out['stages']], ['reference_preparation', 'cohort_context', 'learned_integration', 'feature_export'])
        self.assertEqual(out['stages'][1]['label'], 'Cohort context (two-pass)')
        self.assertEqual(out['baseline_kind'], 'legacy_prior_baseline')
        self.assertIn('not a sequence-only baseline', out['baseline'])

    def test_scoring_modes_are_distinct(self):
        self.assertEqual(self.summarize_with(learned_ledger(), self.report())['scoring']['mode'], 'learned')
        legacy = self.summarize_with(learned_ledger('legacy_heuristic'))
        self.assertEqual((legacy['scoring']['mode'], legacy['learned']['status']), ('legacy_heuristic', 'not_applicable'))
        self.assertEqual(self.summarize(LEDGER)['scoring']['mode'], 'unknown')

    def test_learned_report_availability_identity_and_trust(self):
        ok = self.summarize_with(learned_ledger(), self.report())['learned']
        self.assertEqual((ok['status'], ok['model_id'], ok['hash_verified']), ('available', 'm1', True))
        self.assertEqual(ok['model_sha256'], SHA[:12])
        missing = self.summarize_with(learned_ledger())['learned']
        self.assertEqual(missing['status'], 'unavailable')
        self.assertIn('not zero', missing['reason'])
        self.assertEqual(self.summarize_with(learned_ledger(), self.report(sha='b' * 64))['learned']['status'], 'not_trusted')
        self.assertEqual(self.summarize_with(learned_ledger(), 'not json')['learned']['status'], 'unavailable')
        self.assertEqual(self.summarize_with(learned_ledger('legacy_heuristic'), self.report())['learned']['status'], 'not_trusted')
        unverified = self.summarize_with(learned_ledger(sha=None), self.report())['learned']
        self.assertEqual((unverified['status'], unverified['hash_verified']), ('available', False))


@unittest.skipUnless(shutil.which('node'), 'Node is required')
class PortalTextTest(unittest.TestCase):
    def text(self, summary):
        core = re.search(r'// LEDGER-TEXT-BEGIN\n(.*)// LEDGER-TEXT-END', WEB.read_text(), re.S).group(1)
        with tempfile.TemporaryDirectory() as tmp:
            script = Path(tmp) / 'run.js'
            script.write_text(core + f'\nprocess.stdout.write(ledgerText({json.dumps(summary)}));')
            return subprocess.run(['node', str(script)], check=True, capture_output=True, text=True).stdout

    def test_text_reports_unavailable_failed_and_counts(self):
        unavailable = self.text(ledger_summary.summarize(tempfile.gettempdir() + '/nonexistent-gspa-ledger'))
        self.assertIn('Prediction ledger unavailable', unavailable)
        self.assertNotIn('pairs changed', unavailable)   # unavailable must not read as a count
        with tempfile.TemporaryDirectory() as tmp:
            (Path(tmp) / ledger_summary.LEDGER_NAME).write_text(json.dumps(LEDGER))
            text = self.text(ledger_summary.summarize(tmp))
        self.assertIn('2 executed, 1 skipped, 1 failed, 1 not run', text)
        self.assertIn('3 (protein, function) pairs changed across 2 proteins', text)
        self.assertIn('not measurements of accuracy', text)
        self.assertNotIn('improved', text.lower())

    def test_text_states_scoring_mode_and_model_without_accuracy_claims(self):
        with tempfile.TemporaryDirectory() as tmp:
            (Path(tmp) / ledger_summary.LEDGER_NAME).write_text(json.dumps(learned_ledger()))
            unavailable = self.text(ledger_summary.summarize(tmp))
            self.assertIn('Scoring: learned (frozen model)', unavailable)
            self.assertIn('Learned-stage report unavailable', unavailable)
            (Path(tmp) / ledger_summary.LEARNED_NAME).write_text(json.dumps({'model': {'id': 'm1', 'sha256': SHA}, 'candidates': []}))
            text = self.text(ledger_summary.summarize(tmp))
        self.assertIn('Scoring: learned (frozen model m1, sha256 ' + SHA[:12] + ')', text)
        self.assertNotIn('accura', text.replace('not measurements of accuracy', ''))
        legacy = self.text({'available': True, 'baseline': 'X', 'status_counts': {}, 'net': {'functions': 0, 'proteins': 0, 'by_change': {'added': 0, 'removed': 0, 'increased': 0, 'decreased': 0}},
                            'scoring': {'mode': 'legacy_heuristic'}, 'learned': {'status': 'not_applicable'}})
        self.assertIn('Scoring: legacy heuristic (no trained model)', legacy)
        with tempfile.TemporaryDirectory() as tmp:
            (Path(tmp) / ledger_summary.LEDGER_NAME).write_text(json.dumps(LEDGER))
            self.assertIn('Scoring: unknown', self.text(ledger_summary.summarize(tmp)))

    def test_page_script_parses(self):
        script = re.search(r'<script>(.*)</script>', WEB.read_text(), re.S).group(1)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'page.js'; path.write_text(script)
            subprocess.run(['node', '--check', str(path)], check=True, capture_output=True, text=True)


if __name__ == '__main__': unittest.main()
