"""Web and Nextflow routes compose the canonical integration workflow configuration.

No tool, network or Java process runs here: the tests check that each route emits only keys of the
typed ``IntegrationWorkflowConfig`` schema (the JVM rejects anything else) and no longer carries a
private flag list of its own.
"""
import json
import re
import unittest
from pathlib import Path

from service import integration_options

ROOT = Path(__file__).resolve().parent.parent
SCHEMA = set(json.loads((ROOT / 'gspa-core/src/main/resources/workflow/integration-workflow-config.keys.json').read_text()))
RESOURCES = Path('/resources')


class WebWorkflowConfigTest(unittest.TestCase):
    def build(self, preset, **kw):
        defaults = dict(resources=RESOURCES, claims=Path('/job/claims.jsonl'), browser=Path('/job/gspa_out'), pathways=True)
        defaults.update(kw)
        return integration_options.workflow_config(preset=preset, **defaults)

    def test_every_key_is_in_the_typed_schema(self):
        full = self.build('bacteria', operons=Path('/job/annot/operons_for_integrate.tsv'), fba_gram='neg')
        self.assertLessEqual(set(full), SCHEMA)
        self.assertLessEqual(set(self.build('virus', pathways=False)), SCHEMA)

    def test_bacteria_request_carries_references_and_model_assumptions(self):
        cfg = self.build('bacteria', operons=Path('/job/annot/operons_for_integrate.tsv'), fba_gram='pos')
        self.assertEqual(cfg['enableFbaRefiner'], True)
        self.assertEqual((cfg['fbaGram'], cfg['fbaMedium']), ('pos', 'minimal'))
        self.assertEqual(cfg['operonsFile'], '/job/annot/operons_for_integrate.tsv')
        self.assertEqual(cfg['pathwaysFile'], '/resources/kegg_pathways.tsv')
        self.assertEqual(cfg['moduleFiles'], ['/resources/kegg_modules.tsv'])
        self.assertEqual(cfg['goOwlFile'], '/resources/go/go.owl')

    def test_disabled_analyses_are_not_requested(self):
        cfg = self.build('eukaryote', pathways=False)
        for key in ('pathwaysFile', 'moduleFiles', 'operonsFile', 'enableFbaRefiner'):
            self.assertNotIn(key, cfg)

    def test_virus_integrates_without_priors_or_enforcement(self):
        cfg = self.build('virus', pathways=False)
        self.assertEqual((cfg['enabledPriors'], cfg['enforceStage2'], cfg['lite']), ('', False, True))
        command = integration_options.integrate_command('gspa', 'virus', Path('/job/c.json'))
        self.assertNotIn('--preset', command)

    def test_command_carries_only_the_config_file_and_preset(self):
        command = integration_options.integrate_command('gspa', 'bacteria', Path('/job/c.json'))
        self.assertEqual(command, ['gspa', 'integrate', '--workflow-config', '/job/c.json', '--preset', 'bacteria'])

    def test_portal_has_no_private_integration_flag_list(self):
        source = (ROOT / 'service/analysis.py').read_text()
        for flag in ("'--enable-priors'", "'--pathways'", "'--enable-fba-refiner'", "'--fba-medium'", "'--operons'"):
            self.assertNotIn(flag, source)
        self.assertIn('integration_options.workflow_config', source)


class NextflowWorkflowConfigTest(unittest.TestCase):
    def process(self):
        text = (ROOT / 'gspa-nf/modules/integrate.nf').read_text()
        return text[text.index('process INTEGRATE'):text.index('process VISUALIZE')]

    def test_integrate_composes_the_same_schema(self):
        body = self.process()
        block = body[body.index('def settings = ['):body.index(']\n', body.index('def settings = ['))]
        keys = set(re.findall(r'^\s*(\w+)\s*:', block, re.M)) | set(re.findall(r'settings\.(\w+)\s*=', body))
        self.assertGreaterEqual(len(keys), 8)
        self.assertLessEqual(keys, SCHEMA, sorted(keys - SCHEMA))

    def test_integrate_delegates_to_the_workflow_config_not_a_flag_list(self):
        body = self.process()
        self.assertIn('--workflow-config', body)
        for flag in ('--enable-priors', '--pathways', '--ec2go', '--go-owl', '--operons'):
            self.assertNotIn(flag, body)

    def test_extra_settings_use_the_same_schema_file(self):
        self.assertIn('params.integration_config', self.process())
        self.assertIn('integration_config', (ROOT / 'gspa-nf/nextflow.config').read_text())


if __name__ == '__main__':
    unittest.main()
