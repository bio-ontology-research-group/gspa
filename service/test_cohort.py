"""The web service's multi-genome cohort job is the CLI batch route, offered only on request.

No Java process runs here. The tests pin the job to the Nextflow cohort module (same manifest columns,
same resource keys, same workflow settings, same engine command), check the bundle validation, prove that
the kind is off by default, and prove that single-genome jobs are untouched.
"""
import hashlib
import io
import json
import re
import tempfile
import unittest
import zipfile
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from fastapi.testclient import TestClient

from service import cohort, jobs, portal

ROOT = Path(__file__).resolve().parent.parent
SCHEMA = set(json.loads((ROOT / 'gspa-core/src/main/resources/workflow/integration-workflow-config.keys.json').read_text()))
MODULE = (ROOT / 'gspa-nf/modules/cohort.nf').read_text()


def groovy_list(source, function):
    body = source[source.index(f'def {function}('):]
    body = body[:body.index('\n}\n')]
    return re.findall(r"'([A-Za-z0-9_]+)'", body)


def bundle_bytes(rows=None, metadata=None, orthogroups=None, extra=None, manifest_header=None):
    rows = rows if rows is not None else [{'genome_id': 'G1', 'claims_path': 'claims/G1.jsonl', 'genome_layout_path': 'layout/G1.tsv'},
                                           {'genome_id': 'G2', 'claims_path': 'claims/G2.jsonl', 'genome_layout_path': 'layout/G2.tsv'}]
    header = manifest_header or ['genome_id', 'claims_path', 'genome_layout_path']
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, 'w') as z:
        z.writestr(cohort.MANIFEST, '\n'.join(['\t'.join(header)] + ['\t'.join(r.get(h, '-') for h in header) for r in rows]) + '\n')
        ids = [r['genome_id'] for r in rows]
        z.writestr(cohort.METADATA, json.dumps(metadata if metadata is not None else {'genomes': {g: {} for g in ids}}))
        z.writestr(cohort.ORTHOGROUPS, orthogroups if orthogroups is not None else 'genome_id\tprotein_id\torthogroup_id\tsubfamily_id\tweight\nG1\tp1\tfam1\t-\t1\n')
        for r in rows:
            for key in ('claims_path', 'genome_layout_path'):
                if r.get(key, '-') not in ('-', '') and not r.get('skip_' + key):
                    z.writestr(r[key], 'protein_id\tcontig\tstart\tend\tstrand\np1\tc1\t1\t300\t+\n' if key.endswith('layout_path') else '{}\n')
        for name, content in (extra or {}).items():
            z.writestr(name, content)
    return buf.getvalue()


class ModuleParityTest(unittest.TestCase):
    def test_columns_equal_the_nextflow_cohort_module(self):
        self.assertEqual(cohort.COLUMNS, groovy_list(MODULE, 'cohortColumns'))

    def test_resource_keys_equal_the_nextflow_cohort_module(self):
        allowed = groovy_list(MODULE, 'cohortResources')
        allowed = allowed[:allowed.index('reasonerCacheDir') + 1]
        self.assertEqual(cohort.RESOURCE_KEYS, allowed)
        self.assertEqual(set(cohort.REQUIRED_RESOURCES) | set(cohort.BUNDLE_KEYS),
                         {'reactionGraphFile', 'reactionEcAliasesFile', 'cohortMetadataFile', 'cohortOrthogroupsFile'})

    def test_settings_equal_the_nextflow_cohort_module_defaults(self):
        settings = MODULE[MODULE.index('process INTEGRATE_COHORT'):]
        nf = (ROOT / 'gspa-nf/cohort.nf').read_text()
        self.assertIn("lite: true", nf)
        self.assertIn("params.essential_profile ?: 'bacteria'", nf)
        self.assertIn("params.enable_priors ?: 'essentiality,coherence,gap_filling,genomic_context'", nf)
        self.assertIn("holdout: params.cohort_holdout ?: 'CLADE'", nf)
        self.assertIn("contextMode: params.cohort_context_mode ?: 'UNIVERSAL'", nf)
        self.assertEqual(cohort.SETTINGS, {'lite': True, 'essentialProfile': 'bacteria',
                                           'enabledPriors': 'essentiality,coherence,gap_filling,genomic_context',
                                           'cohortConfig': {'holdout': 'CLADE', 'contextMode': 'UNIVERSAL'}})
        for key in ('config.cohortInference = true', "config.provenanceFile = 'provenance.json'", "config.suggestionsOut = 'suggestions.tsv'",
                    "config.featureExportDir = 'results/features'"):
            self.assertIn(key, settings)

    def test_workflow_config_uses_only_typed_keys_and_the_module_values(self):
        res = {'reactionGraphFile': '/ref/seed_reactions.tsv', 'reactionEcAliasesFile': '/ref/aliases.tsv', 'goOwlFile': '/ref/go.owl',
               'ec2goFile': '/ref/ec2go', 'learnedModelFile': '/models/m.json', 'moduleFiles': ['/ref/kegg_modules.tsv'], 'overlayFile': '/ref/overlay.json'}
        cfg = cohort.workflow_config(resources=res, metadata=Path('/b/cohort_metadata.json'), orthogroups=Path('/b/cohort_orthogroups.tsv'),
                                     out_dir=Path('/job/gspa_out/cohort'), feature_export=True)
        self.assertLessEqual(set(cfg), SCHEMA, sorted(set(cfg) - SCHEMA))
        self.assertEqual({k: cfg[k] for k in cohort.SETTINGS}, cohort.SETTINGS)
        self.assertEqual((cfg['cohortInference'], cfg['cohortMetadataFile'], cfg['cohortOrthogroupsFile']),
                         (True, '/b/cohort_metadata.json', '/b/cohort_orthogroups.tsv'))
        self.assertEqual((cfg['provenanceFile'], cfg['suggestionsOut'], cfg['featureExportDir']),
                         ('/job/gspa_out/cohort/provenance.json', '/job/gspa_out/cohort/suggestions.tsv', '/job/gspa_out/cohort/features'))
        self.assertEqual(cfg['moduleFiles'], ['/ref/kegg_modules.tsv'])
        self.assertNotIn('overlayFile', cfg)               # the CLI applies the overlay itself
        self.assertNotIn('featureExportDir', cohort.workflow_config(resources=res, metadata=Path('/m'), orthogroups=Path('/o'), out_dir=Path('/j')))

    def test_command_is_the_batch_route_without_private_flags(self):
        cmd = cohort.integrate_command('gspa', Path('/job/c.json'), Path('/job/batch.tsv'), '8')
        self.assertEqual(cmd, ['gspa', 'integrate', '--workflow-config', '/job/c.json', '--batch-manifest', '/job/batch.tsv', '--batch-threads', '8'])
        self.assertEqual(cohort.integrate_command('gspa', Path('/c'), Path('/b'), 4, Path('/o.json'))[-2:], ['--workflow-overlay', '/o.json'])
        self.assertIn("java -jar engine/gspa.jar integrate --workflow-config workflow.json --batch-manifest batch.tsv --batch-threads", MODULE)
        source = (ROOT / 'service/cohort.py').read_text()
        for flag in ("'--cohort'", "'--enable-priors'", "'--learned-model'", "'--cohort-metadata'", "'--cohort-orthogroups'"):
            self.assertNotIn(flag, source)

    def test_batch_manifest_has_the_module_columns_and_service_owned_outputs(self):
        rows = [{'genome_id': 'G1', 'claims_path': Path('/b/c1'), 'genome_layout_path': Path('/b/l1'), 'taxon': 'NCBITaxon_1'}]
        text = cohort.batch_manifest(rows, Path('/job/gspa_out/cohort'))
        header, row = text.rstrip('\n').split('\n')
        self.assertEqual(header.split('\t'), cohort.COLUMNS)
        values = dict(zip(cohort.COLUMNS, row.split('\t')))
        self.assertEqual((values['genome_id'], values['claims_path'], values['out_path'], values['taxon'], values['gaps_path']),
                         ('G1', '/b/c1', '/job/gspa_out/cohort/G1/integrated.tsv', 'NCBITaxon_1', '-'))


class RenderLayoutTest(unittest.TestCase):
    def test_batch_artifacts_are_linked_under_the_single_job_names(self):
        # Regression: the first web cohort job integrated all genomes identically to the CLI batch and then failed at
        # the browser step because the renderer looks for gspa_out/provenance.json, not integrated.tsv.provenance.json.
        with tempfile.TemporaryDirectory() as tmp:
            run = Path(tmp) / 'gspa_out/cohort/G1'
            run.mkdir(parents=True)
            for name in ('integrated.tsv', 'integrated.tsv.stages.json', 'integrated.tsv.learned.json', 'integrated.tsv.provenance.json',
                         'integrated.tsv.enforcement.tsv', 'integrated.tsv.cohort.json'):
                (run / name).write_text(name)
            target = cohort.render_run_dir(run, Path(tmp) / 'render/G1/gspa_out')
            names = sorted(p.name for p in target.iterdir())
            self.assertEqual(names, ['actions.tsv', 'integrated.tsv', 'integrated.tsv.cohort.json', 'integrated.tsv.learned.json',
                                     'integrated.tsv.stages.json', 'provenance.json'])
            self.assertEqual((target / 'provenance.json').read_text(), 'integrated.tsv.provenance.json')
            self.assertTrue(all(p.is_symlink() for p in target.iterdir()))
            cohort.render_run_dir(run, target)      # idempotent
            self.assertEqual(sorted(p.name for p in target.iterdir()), names)


class BundleValidationTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.folder = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def accept(self, data):
        (self.folder / 'cohort.zip').write_bytes(data)
        return cohort.accept(self.folder)

    def test_valid_bundle_is_accepted_sorted(self):
        out = self.accept(bundle_bytes())
        self.assertEqual(out, {'genomes': ['G1', 'G2'], 'proteins': 2})
        rows = cohort.read_bundle(self.folder / 'bundle')['rows']
        self.assertTrue(all(r['claims_path'].is_file() and r['genome_layout_path'].is_file() for r in rows))

    def test_rejects_unsafe_and_unknown_inputs(self):
        cases = {
            'absolute path': bundle_bytes(rows=[{'genome_id': 'G1', 'claims_path': '/etc/passwd', 'genome_layout_path': 'layout/G1.tsv', 'skip_claims_path': True}]),
            'parent path': bundle_bytes(rows=[{'genome_id': 'G1', 'claims_path': '../x.jsonl', 'genome_layout_path': 'layout/G1.tsv', 'skip_claims_path': True}]),
            'unknown column': bundle_bytes(manifest_header=['genome_id', 'claims_path', 'genome_layout_path', 'taxon_lineage'],
                                           rows=[{'genome_id': 'G1', 'claims_path': 'claims/G1.jsonl', 'genome_layout_path': 'layout/G1.tsv', 'taxon_lineage': 'x'}]),
            'duplicate id': bundle_bytes(rows=[{'genome_id': 'G1', 'claims_path': 'claims/G1.jsonl', 'genome_layout_path': 'layout/G1.tsv'}] * 2),
            'unsafe id': bundle_bytes(rows=[{'genome_id': '../G1', 'claims_path': 'claims/G1.jsonl', 'genome_layout_path': 'layout/G1.tsv'}]),
            'missing layout': bundle_bytes(rows=[{'genome_id': 'G1', 'claims_path': 'claims/G1.jsonl', 'genome_layout_path': '-'}]),
            'missing claims file': bundle_bytes(rows=[{'genome_id': 'G1', 'claims_path': 'claims/none.jsonl', 'genome_layout_path': 'layout/G1.tsv', 'skip_claims_path': True}]),
            'metadata mismatch': bundle_bytes(metadata={'genomes': {'G1': {}, 'G3': {}}}),
            'bad orthogroups header': bundle_bytes(orthogroups='protein\tfamily\n'),
            'unsafe member': bundle_bytes(extra={'../escape.txt': 'x'}),
            'not a zip': b'garbage',
        }
        for name, data in cases.items():
            with self.subTest(name):
                with self.assertRaises(cohort.BundleError):
                    self.accept(data)
                for p in self.folder.iterdir():
                    if p.is_dir():
                        import shutil
                        shutil.rmtree(p)

    def test_resources_are_validated_against_the_module_keys(self):
        ref = self.folder / 'ref'
        ref.mkdir()
        for name in ('seed_reactions.tsv', 'aliases.tsv', 'go.owl'):
            (ref / name).write_text('x')
        spec = self.folder / 'cohort_resources.json'
        spec.write_text(json.dumps({'reactionGraphFile': 'ref/seed_reactions.tsv', 'reactionEcAliasesFile': 'ref/aliases.tsv', 'goOwlFile': str(ref / 'go.owl')}))
        res = cohort.resources({'GSPA_COHORT_RESOURCES': str(spec)})
        self.assertEqual(res['reactionGraphFile'], str(ref / 'seed_reactions.tsv'))
        self.assertEqual(res['goOwlFile'], str(ref / 'go.owl'))
        for bad in ({'reactionGraphFile': 'ref/seed_reactions.tsv'},
                    {'reactionGraphFile': 'ref/seed_reactions.tsv', 'reactionEcAliasesFile': 'ref/aliases.tsv', 'cohortMetadataFile': 'ref/go.owl'},
                    {'reactionGraphFile': 'ref/seed_reactions.tsv', 'reactionEcAliasesFile': 'ref/aliases.tsv', 'privateKey': 'x'},
                    {'reactionGraphFile': 'ref/seed_reactions.tsv', 'reactionEcAliasesFile': 'ref/missing.tsv'}):
            spec.write_text(json.dumps(bad))
            with self.assertRaises(ValueError):
                cohort.resources({'GSPA_COHORT_RESOURCES': str(spec)})


class PortalCohortTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = patch.object(jobs, 'ROOT', Path(self.tmp.name)); self.root.start()
        self.disk = patch.object(jobs.shutil, 'disk_usage', return_value=SimpleNamespace(free=100 * 1024**3)); self.disk.start()
        self.client = TestClient(portal.app)

    def tearDown(self):
        self.disk.stop(); self.root.stop(); self.tmp.cleanup()

    def post(self, data, environ):
        with patch.object(jobs, 'ready', return_value=True), patch.dict(cohort.os.environ, environ, clear=False):
            for key in ('GSPA_COHORT_RESOURCES',):
                if key not in environ:
                    cohort.os.environ.pop(key, None)
            return self.client.post('/api/jobs', data={'kind': 'cohort', 'preset': 'bacteria'}, files={'file': ('cohort.zip', data)})

    def test_cohort_kind_is_off_by_default(self):
        r = self.post(bundle_bytes(), {})
        self.assertEqual(r.status_code, 400, r.text)
        self.assertIn('not offered', r.json()['detail'])
        self.assertEqual(list(jobs.ROOT.iterdir()), [])

    def test_cohort_kind_is_accepted_when_offered_and_invalid_bundles_leave_nothing(self):
        env = {'GSPA_COHORT_RESOURCES': '/deployment/cohort_resources.json'}
        r = self.post(bundle_bytes(), env)
        self.assertEqual(r.status_code, 202, r.text)
        job = r.json()
        self.assertEqual((job['kind'], job['genomes'], job['analyses']), ('cohort', ['G1', 'G2'], {}))
        self.assertTrue((jobs.ROOT / job['id'] / 'bundle' / cohort.MANIFEST).is_file())
        before = {p.name for p in jobs.ROOT.iterdir()}
        r = self.post(bundle_bytes(extra={'../escape.txt': 'x'}), env)
        self.assertEqual(r.status_code, 400, r.text)
        self.assertEqual({p.name for p in jobs.ROOT.iterdir()}, before)
        # per-genome artifacts are refused until the job completes; unknown genomes are never resolved
        self.assertEqual(self.client.get(f'/api/jobs/{job["id"]}/genomes/G1/browser').status_code, 404)
        self.assertEqual(self.client.get(f'/api/jobs/{job["id"]}/genomes/G9/browser').status_code, 404)

    def test_status_of_a_complete_cohort_job_has_one_ledger_per_genome(self):
        env = {'GSPA_COHORT_RESOURCES': '/deployment/cohort_resources.json'}
        job = self.post(bundle_bytes(), env).json()
        folder = jobs.ROOT / job['id']
        run = folder / 'gspa_out' / 'cohort' / 'G1'
        run.mkdir(parents=True)
        (run / 'integrated.tsv').write_text('protein_id\n')
        jobs.state(folder, status='complete')
        data = self.client.get('/api/jobs/' + job['id']).json()
        self.assertNotIn('ledger', data)
        self.assertEqual(set(data['genome_ledgers']), {'G1', 'G2'})
        self.assertFalse(data['genome_ledgers']['G1']['available'])        # no ledger written: unavailable, not zero
        self.assertEqual(self.client.get(f'/api/jobs/{job["id"]}/genomes/G1/table').status_code, 200)
        self.assertEqual(self.client.get(f'/api/jobs/{job["id"]}/genomes/G1/ledger').status_code, 404)


class SingleGenomeUnchangedTest(unittest.TestCase):
    """Single-genome jobs are byte-identical in behaviour: same dispatch, same worker module, no cohort import."""

    def test_single_kinds_dispatch_to_run_full_and_never_import_cohort(self):
        import sys
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            calls = []
            fake = SimpleNamespace(run_full=lambda f: calls.append(f))
            for kind in ('gff3', 'genome', 'proteins'):
                jobs.state(folder, kind=kind, created=1.0)
                with patch.dict(sys.modules, {'service.analysis': fake}):
                    sys.modules.pop('service.cohort', None)
                    jobs.run(folder)
                    self.assertNotIn('service.cohort', sys.modules)
            self.assertEqual(calls, [folder] * 3)
            jobs.state(folder, kind='cohort')
            with patch('service.cohort.run_cohort') as rc:
                jobs.run(folder)
                rc.assert_called_once_with(folder)

    def test_single_genome_worker_does_not_know_the_cohort_module(self):
        source = (ROOT / 'service/analysis.py').read_text()
        self.assertNotIn('cohort', source)

    def test_single_genome_worker_source_hash_is_recorded(self):
        # The hash recorded before this change (benchmark/v3/results/routes-learned-cohort-v1/plan/service-before.sha256);
        # analysis.py is the single-genome worker and must be byte-identical.
        recorded = ROOT / 'benchmark/v3/results/routes-learned-cohort-v1/plan/service-before.sha256'
        if not recorded.exists():
            self.skipTest('evidence folder not present')
        expected = dict(line.split()[::-1] for line in recorded.read_text().splitlines() if line.strip())
        self.assertEqual(hashlib.sha256((ROOT / 'service/analysis.py').read_bytes()).hexdigest(), expected['service/analysis.py'])
        self.assertEqual(hashlib.sha256((ROOT / 'service/integration_options.py').read_bytes()).hexdigest(), expected['service/integration_options.py'])

    def test_portal_rejects_cohort_kind_exactly_like_before_when_not_offered(self):
        # The pre-existing kinds and their validation messages are untouched.
        source = (ROOT / 'service/portal.py').read_text()
        self.assertIn("if kind not in ('gff3', 'genome', 'proteins', 'cohort') or preset not in PRESETS or gram not in ('neg', 'pos'):", source)
        self.assertIn("raise HTTPException(400, 'Select a valid input type and biological context.')", source)


if __name__ == '__main__':
    unittest.main()
