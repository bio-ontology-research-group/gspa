"""Real JVM writer -> results browser and portal summary, for the learned stage.

`LearnedArtifactsWriterSpec` runs the shared integration workflow on the `ContextFixtures` reaction-network
genome with a frozen, hand-set, TEST ONLY learned model (`intended_use = mechanism_test`) and writes the
actual `integrated.tsv`, `.stages.json`, `.learned.json` and provenance. This test reads those files through
`make_viz.py` and `service/ledger_summary.py`. It checks that the writer and the readers agree on the file
contract. The model was never trained: nothing here validates a trained model or says anything about accuracy.
(The production loader refuses a mechanism-test model, so this goes through the workflow's test seam, not
through `gspa integrate`.)

Set GSPA_E2E_SKIP_JVM=1 to skip it where Gradle cannot run.
"""
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import unittest

from scripts.tests import test_visualize_ledger as base
from service import ledger_summary
from service import test_ledger_summary

ROOT = Path(__file__).resolve().parents[2]
KEY = ('p3', 'EC', 'EC:2.2.2.2')


@unittest.skipIf(os.environ.get('GSPA_E2E_SKIP_JVM') == '1' or not shutil.which('java'), 'needs a JVM and Gradle')
class RealWriterToBrowserTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.tmp.cleanup)
        cls.root = Path(cls.tmp.name)
        cls.run_dir = cls.root / 'gspa_out'
        spec = 'gspa.integration.workflow.LearnedArtifactsWriterSpec'
        # cleanTest: an up-to-date test task would not rewrite the artifacts
        subprocess.run([str(ROOT / 'gradlew'), '--offline', '-q', ':gspa-core:cleanTest', ':gspa-core:test', '--tests', spec],
                       cwd=ROOT, check=True, capture_output=True, text=True, timeout=900,
                       env={**os.environ, 'GSPA_LEARNED_E2E_DIR': str(cls.run_dir)})

    def stages_json(self):
        return json.loads((self.run_dir / 'integrated.tsv.stages.json').read_text())

    def render(self, run_dir=None):
        run_dir = Path(run_dir or self.run_dir)
        root = run_dir.parent
        ids = [line.split('\t')[0] for line in (run_dir / 'g1.layout.tsv').read_text().splitlines()[1:]]
        (root / 'prokka_out').mkdir(exist_ok=True)
        (root / 'prokka_out/g1.tsv').write_text('locus_tag\tftype\tlength_bp\tproduct\n' + ''.join(f'{p}\tCDS\t300\thypothetical protein\n' for p in ids))
        (root / 'prokka_out/g1.gff').write_text('##gff-version 3\n' + ''.join(
            f'c1\ttest\tCDS\t{i * 1000 + 1}\t{i * 1000 + 900}\t.\t+\t0\tID=gene{i};locus_tag={p}\n' for i, p in enumerate(ids, 1)))
        (root / 'go.obo').write_text(base.OBO)
        (root / 'ec2go').write_text('! none\n')
        out = root / 'browser.html'
        subprocess.run([sys.executable, str(base.SCRIPT)], check=True, capture_output=True, text=True, env={
            **os.environ, 'GSPA_WORKDIR': str(root), 'GSPA_GENOME_ID': 'g1', 'GSPA_OUT': str(out), 'GSPA_RUN_DIR': run_dir.name,
            'GSPA_FASTA': str(root / 'none.fna'), 'GSPA_GO_OBO': str(root / 'go.obo'), 'GSPA_EC2GO': str(root / 'ec2go')})
        return json.loads(re.search(r'const DATA = (.*);', out.read_text()).group(1))['ledger']

    def test_browser_reads_the_real_ledger_and_learned_report(self):
        ledger = self.render()
        self.assertEqual(ledger['status'], 'available')
        self.assertEqual((ledger['scoring']['mode'], ledger['scoring']['consistent'], ledger['scoring']['learned_ran']), ('learned', True, True))
        learned = ledger['learned']
        self.assertEqual((learned['status'], learned['hash_verified']), ('available', True))
        self.assertEqual(learned['model']['sha256'], self.stages_json_stage('learned_integration')['model_sha256'])
        self.assertEqual(learned['model']['intended_use'], 'mechanism_test')
        cand = next(c for c in learned['candidates'] if (c['p'], c['t'], c['f']) == KEY and c['role'] == 'CATALYTIC')
        self.assertTrue(cand['joined'])
        self.assertIsNotNone(cand['pr'])
        self.assertEqual(learned['n_unjoined_candidates'], 0)
        self.assertEqual([s['name'] for s in ledger['stages']], [s['name'] for s in self.stages_json()['stages']])

    def stages_json_stage(self, name):
        return next(s for s in self.stages_json()['stages'] if s['name'] == name)

    def test_portal_summary_reads_the_same_files(self):
        summary = ledger_summary.summarize(self.run_dir)
        self.assertEqual(summary['scoring']['mode'], 'learned')
        self.assertEqual((summary['learned']['status'], summary['learned']['hash_verified']), ('available', True))
        if shutil.which('node'):
            text = test_ledger_summary.PortalTextTest().text(summary)
            self.assertIn('Scoring: learned (frozen model', text)
            self.assertNotIn('unverified', text)
            self.assertNotIn('accura', text.replace('not measurements of accuracy', ''))

    def test_a_ledger_without_the_model_hash_stays_unverified_in_both_readers(self):
        copy = Path(tempfile.mkdtemp(dir=self.root)) / 'gspa_out'
        shutil.copytree(self.run_dir, copy)
        ledger = json.loads((copy / 'integrated.tsv.stages.json').read_text())
        next(s for s in ledger['stages'] if s['name'] == 'learned_integration').pop('model_sha256')
        (copy / 'integrated.tsv.stages.json').write_text(json.dumps(ledger))
        learned = self.render(copy)['learned']
        self.assertEqual((learned['status'], learned['hash_verified']), ('available', False))
        self.assertIn('unverified', learned['hash_note'])
        summary = ledger_summary.summarize(copy)
        self.assertEqual((summary['learned']['status'], summary['learned']['hash_verified']), ('available', False))
        if shutil.which('node'): self.assertIn('model hash unverified', test_ledger_summary.PortalTextTest().text(summary))

    def test_a_failed_learned_stage_withholds_the_real_report(self):
        copy = Path(tempfile.mkdtemp(dir=self.root)) / 'gspa_out'
        shutil.copytree(self.run_dir, copy)
        ledger = json.loads((copy / 'integrated.tsv.stages.json').read_text())
        next(s for s in ledger['stages'] if s['name'] == 'learned_integration').update(status='failed', reason='injected')
        (copy / 'integrated.tsv.stages.json').write_text(json.dumps(ledger))
        learned = self.render(copy)['learned']
        self.assertEqual(learned['status'], 'not_trusted')
        self.assertNotIn('candidates', learned)
        self.assertEqual(ledger_summary.summarize(copy)['learned']['status'], 'not_trusted')


if __name__ == '__main__':
    unittest.main()
