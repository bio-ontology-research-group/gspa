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

    def test_page_script_parses(self):
        script = re.search(r'<script>(.*)</script>', WEB.read_text(), re.S).group(1)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'page.js'; path.write_text(script)
            subprocess.run(['node', '--check', str(path)], check=True, capture_output=True, text=True)


if __name__ == '__main__': unittest.main()
