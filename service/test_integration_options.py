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


class WebClaimsTest(unittest.TestCase):
    RAW = [
        # a weak candidate below the display threshold, with its evidence type and metadata
        {'protein_id': 'p1', 'function_type': 'GO', 'function_id': 'GO:0004340', 'source': 'diamond',
         'evidence_type': 'SEQUENCE_SIMILARITY', 'go_aspect': 'MF', 'raw_score': 0.12, 'calibrated_prob': 0.0731,
         'metadata': {'target': 'P12345', 'evalue': 1e-9}},
        {'protein_id': 'p1', 'function_type': 'EC', 'function_id': '2.7.1.1', 'source': 'deepgo-plusplus-light',
         'evidence_type': 'SEQUENCE_DEEPLEARNING', 'go_aspect': None, 'raw_score': 0.31, 'calibrated_prob': 0.2, 'metadata': {}},
        {'protein_id': 'p2', 'function_type': 'GO', 'function_id': 'GO:0006096', 'source': 'deepgo-plusplus-light',
         'evidence_type': 'SEQUENCE_DEEPLEARNING', 'go_aspect': None, 'raw_score': 0.77, 'calibrated_prob': 0.7, 'metadata': {}},
        {'protein_id': 'p2', 'function_type': 'GO', 'function_id': 'GO:0003674', 'source': 'operon',
         'evidence_type': 'GENOMIC_CONTEXT', 'go_aspect': 'MF', 'raw_score': 0.9, 'calibrated_prob': 0.9, 'metadata': {}},
    ]

    def run_build(self, **kw):
        import tempfile
        tmp = Path(tempfile.mkdtemp())
        raw = tmp / 'raw.jsonl'
        raw.write_text(''.join(json.dumps(r) + '\n' for r in self.RAW))
        seen = integration_options.build_claims(raw, tmp / 'claims.jsonl', **kw)
        lines = (tmp / 'claims.jsonl').read_text().splitlines()
        return seen, lines, raw.read_text().splitlines()

    def test_raw_claims_pass_through_verbatim_including_weak_candidates(self):
        _, lines, raw = self.run_build(aspects={})
        self.assertEqual(lines, [raw[0], raw[1], raw[2]])       # operon context is out of the portal's scope
        first = json.loads(lines[0])
        self.assertEqual((first['raw_score'], first['calibrated_prob'], first['evidence_type']),
                         (0.12, 0.0731, 'SEQUENCE_SIMILARITY'))
        self.assertEqual(first['metadata']['target'], 'P12345')

    def test_only_a_missing_go_aspect_is_filled(self):
        _, lines, raw = self.run_build(aspects={'GO:0006096': 'BP', 'GO:0004340': 'CC'})
        self.assertEqual(lines[0], raw[0])                      # already had an aspect: untouched
        filled = json.loads(lines[2])
        self.assertEqual(filled['go_aspect'], 'BP')
        self.assertEqual({k: v for k, v in filled.items() if k != 'go_aspect'},
                         {k: v for k, v in self.RAW[2].items() if k != 'go_aspect'})

    def test_excluded_source_is_dropped_without_touching_others(self):
        _, lines, _ = self.run_build(exclude_sources=('diamond',))
        self.assertEqual([json.loads(l)['source'] for l in lines], ['deepgo-plusplus-light'] * 2)

    def test_added_psortb_claims_carry_no_invented_probability(self):
        extra = [{'protein_id': 'p1', 'type': 'GO', 'value': 'GO:0005737', 'source': 'psortb', 'score': 0.8, 'aspect': 'CC'},
                 {'protein_id': 'p1', 'type': 'GO', 'value': 'GO:0005737', 'source': 'psortb', 'score': 0.5, 'aspect': 'CC'}]
        seen, lines, _ = self.run_build(extra_claims=extra)
        psort = [json.loads(l) for l in lines if json.loads(l)['source'] == 'psortb']
        self.assertEqual(len(psort), 1)
        self.assertEqual(psort[0]['raw_score'], 0.8)
        self.assertNotIn('calibrated_prob', psort[0])           # the integrator calibrates it like any raw score
        self.assertIn(('p1', 'GO', 'GO:0005737', 'psortb'), seen)

    def test_portal_consumes_the_annotator_export(self):
        source = (ROOT / 'service/analysis.py').read_text()
        self.assertIn("'--claims-out'", source)
        self.assertIn('integration_options.build_claims', source)
        self.assertNotIn('calibrated_prob=min(1-1e-6', source)


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

    def test_reference_files_are_staged_not_host_paths(self):
        body = self.process()
        # theta and the overlay are task inputs (content-hashed), never host paths in the script.
        self.assertIn('path theta_file', body)
        self.assertIn('path overlay_file', body)
        self.assertNotIn('params.theta_file', body)
        self.assertNotIn('params.integration_config', body)
        self.assertNotIn('new File(', body)
        main = (ROOT / 'gspa-nf/main.nf').read_text()
        self.assertIn("file(params.integration_config) : empty_file('NO_INTEGRATION_CONFIG')", main)
        self.assertIn("file(params.theta_file) : empty_file('NO_THETA')", main)

    def test_overlay_cannot_rebind_a_sample(self):
        # The overlay is applied by the JVM, which rejects genomeId/claims/out and file-valued keys;
        # the process must not merge raw JSON over its bindings itself.
        body = self.process()
        self.assertIn('--workflow-overlay', body)
        self.assertNotIn('JsonSlurper', body)
        self.assertNotIn('settings <<', body)
        self.assertIn('integration_config', (ROOT / 'gspa-nf/nextflow.config').read_text())


if __name__ == '__main__':
    unittest.main()
