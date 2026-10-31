"""Execute real Nextflow and the production JAR; opt in with GSPA_NEXTFLOW_BIN.

Synthetic genomes test route parity and failure/cache contracts, not accuracy.
"""
import csv
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[3]
NEXTFLOW = os.environ.get('GSPA_NEXTFLOW_BIN')
ENGINE = Path(os.environ.get('GSPA_JAR', ROOT / 'gspa-cli/build/libs/gspa-3.0.0.jar'))


@unittest.skipUnless(NEXTFLOW and Path(NEXTFLOW).is_file() and ENGINE.is_file(),
                     'set GSPA_NEXTFLOW_BIN to the official launcher and build the production engine')
class NextflowCohortTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='gspa-nf-cohort-')
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.rows = []
        genomes = {}
        groups = ['genome_id\tprotein_id\torthogroup_id']
        for g in ('g1', 'g2', 'g3', 'g4'):
            directory = self.root / g
            directory.mkdir()
            claims = []
            for pid, ec, score in [('p2', '1.1.1.1', .9), ('p3', '2.2.2.2', .4), ('p4', '3.3.3.3', .9)]:
                claims.append({'protein_id': pid, 'function_type': 'EC', 'function_id': 'EC:' + ec,
                               'source': 'diamond', 'raw_score': score, 'calibrated_prob': score,
                               'evidence_type': 'SEQUENCE_SIMILARITY'})
                groups.append(f'{g}\t{pid}\tOG_{pid}')
            (directory / 'claims.jsonl').write_text(''.join(json.dumps(c) + '\n' for c in claims))
            (directory / 'layout.tsv').write_text('protein_id\tcontig\tstart\tend\tstrand\n' +
                ''.join(f'p{i}\tchr\t{i*1000}\t{i*1000+900}\t+\n' for i in range(1, 9)))
            (directory / 'operons.tsv').write_text('p2\tp3\tp4\n')
            self.rows.append({'genome_id': g, 'claims_path': str(directory / 'claims.jsonl'),
                              'genome_layout_path': str(directory / 'layout.tsv'),
                              'operons_path': str(directory / 'operons.tsv')})
            genomes[g] = {'clade': 'A' if g in ('g1', 'g2') else g, 'strainCluster': g,
                          'assemblyStatus': 'COMPLETE', 'familySearchComplete': True,
                          'contigs': {'chr': {'length': 9000, 'closed': True}}}
        graph = self.write('graph/data.tsv', 'id\tstoichiometry\tec_numbers\n' +
            'r1\t-1:a:0:0;1:b:0:0\t1.1.1.1\n' +
            'r2\t-1:b:0:0;1:c:0:0\t2.2.2.2\n' +
            'r3\t-1:c:0:0;1:d:0:0\t3.3.3.3\n')
        aliases = self.write('aliases/data.tsv', 'MS ID\tOld MS ID\tExternal ID\tSource\n' +
            'r1\t\t1.1.1.1\tt\nr2\t\t2.2.2.2\tt\nr3\t\t3.3.3.3\tt\n')
        metadata = self.write('metadata.json', json.dumps({'genomes': genomes}))
        og = self.write('orthogroups.tsv', '\n'.join(groups) + '\n')
        overlay = self.write("overlay's config.json", json.dumps({'rlcCurrencyPercentile': 100.0}))
        self.resources = {'reactionGraphFile': str(graph), 'reactionEcAliasesFile': str(aliases),
                          'cohortMetadataFile': str(metadata), 'cohortOrthogroupsFile': str(og),
                          'overlayFile': str(overlay),
                          'learnedModelFile': str(ROOT / 'benchmark/v3/results/cohort-evidence-dedup-v1/models/model_logistic.json')}
        self.write('resources.json', json.dumps(self.resources))
        self.write_manifest()
        self.config = self.write('test.config', '''
process.executor = 'local'
process.cpus = 2
trace.overwrite = true
docker.enabled = false
singularity.enabled = false
params.essential_profile = 'bacteria'
params.enable_priors = 'none'
params.cohort_holdout = 'CLADE'
params.feature_export = true
''')
        self.env = dict(os.environ, NXF_OFFLINE='true', NXF_ANSI_LOG='false')

    def write(self, name, text):
        path = self.root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
        return path

    def write_manifest(self):
        with (self.root / 'manifest.tsv').open('w') as fh:
            writer = csv.DictWriter(fh, fieldnames=list(self.rows[0]), delimiter='\t')
            writer.writeheader()
            writer.writerows(self.rows)

    def nf(self, *, resume=False, ok=True):
        cmd = [NEXTFLOW, '-C', str(self.config), 'run', str(ROOT / 'gspa-nf/cohort.nf'),
               '-work-dir', str(self.root / 'work'), '-with-trace', str(self.root / 'trace.tsv'),
               '--cohort_manifest', str(self.root / 'manifest.tsv'), '--cohort_resources', str(self.root / 'resources.json'),
               '--gspa_jar', str(ENGINE), '--outdir', str(self.root / 'published')]
        if resume:
            cmd.append('-resume')
        proc = subprocess.run(cmd, cwd=self.root, env=self.env, text=True, stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT, timeout=180)
        self.assertEqual(proc.returncode == 0, ok, proc.stdout[-15000:])
        return proc.stdout

    def direct(self):
        rows = [{**r, 'out_path': str(self.root / 'direct' / (r['genome_id'] + '.tsv'))} for r in self.rows]
        (self.root / 'direct').mkdir()
        columns = ['genome_id', 'claims_path', 'out_path', 'taxon', 'essential_profile', 'operons_path',
                   'gaps_path', 'max_iter', 'enforce_stage2', 'enforce_mode', 'enable_priors', 'genome_layout_path']
        manifest = self.root / 'direct.tsv'
        with manifest.open('w') as fh:
            writer = csv.DictWriter(fh, fieldnames=columns, delimiter='\t', restval='-')
            writer.writeheader()
            writer.writerows(rows)
        config = {k: v for k, v in self.resources.items() if k != 'overlayFile'}
        config.update(lite=True, enabledPriors='none', essentialProfile='bacteria', cohortInference=True,
                      cohortConfig={'holdout': 'CLADE'}, provenanceFile='enabled', suggestionsOut='enabled',
                      featureExportDir=str(self.root / 'direct/features'), rlcCurrencyPercentile=100.0)
        conf = self.write('direct.json', json.dumps(config))
        proc = subprocess.run(['java', '-jar', str(ENGINE), 'integrate', '--workflow-config', str(conf),
                               '--batch-manifest', str(manifest), '--batch-threads', '2'], cwd=self.root,
                              env=self.env, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=120)
        self.assertEqual(proc.returncode, 0, proc.stdout[-6000:])

    def test_joint_route_matches_direct_cli_scores_context_and_stage_states(self):
        self.nf()
        self.direct()
        published = self.root / 'published/cohort'
        for g in ('g1', 'g2', 'g3', 'g4'):
            nf = published / g / f'{g}_integrated.tsv'
            direct = self.root / 'direct' / f'{g}.tsv'
            self.assertEqual(nf.read_text(), direct.read_text())
            a = json.loads(Path(str(nf) + '.cohort.json').read_text())
            b = json.loads(Path(str(direct) + '.cohort.json').read_text())
            self.assertEqual(a, b)
            self.assertGreater(a['n_peer_genomes'], 0)
            if g == 'g1':
                self.assertEqual(set(a['held_out_genomes']), {'g1', 'g2'})
            for suffix in ('.stages.json', '.learned.json', '.provenance.json'):
                self.assertTrue(Path(str(nf) + suffix).is_file())
            stages = lambda p: [(s['name'], s['status']) for s in json.loads(p.read_text())['stages']]
            self.assertEqual(stages(Path(str(nf)+'.stages.json')), stages(Path(str(direct)+'.stages.json')))
            self.assertIn(('cohort_context', 'executed'), stages(Path(str(nf)+'.stages.json')))
            # The cohort model's observed source is swissprot_diamond_baseline;
            # it cannot treat these synthetic DIAMOND inputs as calibrated coverage.
            self.assertIn(('learned_integration', 'skipped'), stages(Path(str(nf)+'.stages.json')))
            learned = json.loads(Path(str(nf)+'.learned.json').read_text())
            self.assertTrue(learned['candidates'])
            self.assertTrue(all(c['score_status'] == 'unsupported' and c['probability'] is None
                                for c in learned['candidates']))
            self.assertTrue((published / 'features' / g / 'candidate_features.jsonl').is_file())
        with (self.root / 'trace.tsv').open() as fh:
            tasks = list(csv.DictReader(fh, delimiter='\t'))
        self.assertEqual(len(tasks), 1, tasks)

    def test_version_two_model_is_rejected_by_the_actual_route(self):
        self.resources['learnedModelFile'] = str(ROOT / 'benchmark/v3/results/model-input-support-v1/models/model_logistic.json')
        self.write('resources.json', json.dumps(self.resources))
        self.nf(ok=False)
        # Nextflow truncates the displayed stderr tail; inspect the actual task
        # logs so an unrelated failure cannot satisfy this compatibility test.
        output = '\n'.join(p.read_text() for p in (self.root / 'work').rglob('.command.err'))
        self.assertIn('feature contract', output)
        self.assertIn('gspa.context-features/2', output)
        self.assertIn('gspa.context-features/3', output)

    def test_deep_cache_reuses_then_invalidates_same_size_same_mtime_claim_edit(self):
        self.nf()
        self.nf(resume=True)
        with (self.root / 'trace.tsv').open() as fh:
            self.assertEqual(next(csv.DictReader(fh, delimiter='\t'))['status'], 'CACHED')
        path = Path(self.rows[0]['claims_path'])
        old = path.stat()
        path.write_text(path.read_text().replace('0.4', '0.8'))
        os.utime(path, ns=(old.st_atime_ns, old.st_mtime_ns))
        self.nf(resume=True)
        with (self.root / 'trace.tsv').open() as fh:
            self.assertEqual(next(csv.DictReader(fh, delimiter='\t'))['status'], 'COMPLETED')

    def test_duplicate_genomes_fail_instead_of_becoming_separate_singletons(self):
        self.rows.append(dict(self.rows[0]))
        self.write_manifest()
        self.assertIn('distinct, safe genome IDs', self.nf(ok=False))

    def test_swapping_claim_files_between_genomes_invalidates_cache(self):
        path = Path(self.rows[0]['claims_path'])
        path.write_text(path.read_text().replace('0.4', '0.8'))
        self.nf()
        output = self.root / 'published/cohort/g1/g1_integrated.tsv'
        before = output.read_text()
        self.rows[0]['claims_path'], self.rows[1]['claims_path'] = self.rows[1]['claims_path'], self.rows[0]['claims_path']
        self.write_manifest()
        self.nf(resume=True)
        with (self.root / 'trace.tsv').open() as fh:
            self.assertEqual(next(csv.DictReader(fh, delimiter='\t'))['status'], 'COMPLETED')
        self.assertNotEqual(before, output.read_text())

    def test_swapping_resource_roles_does_not_reuse_successful_cache(self):
        self.nf()
        self.resources['overlayFile'], self.resources['cohortMetadataFile'] = (
            self.resources['cohortMetadataFile'], self.resources['overlayFile'])
        self.write('resources.json', json.dumps(self.resources))
        # Same file collection, invalid new roles: the JVM must run and reject the input.
        self.nf(resume=True, ok=False)
        with (self.root / 'trace.tsv').open() as fh:
            self.assertEqual(next(csv.DictReader(fh, delimiter='\t'))['status'], 'FAILED')

    def test_missing_layout_fails_before_inference(self):
        Path(self.rows[0]['genome_layout_path']).unlink()
        self.assertIn('layout.tsv', self.nf(ok=False))

    def test_empty_manifest_fails_instead_of_successfully_doing_nothing(self):
        self.write('manifest.tsv', 'genome_id\tclaims_path\tgenome_layout_path\n')
        self.assertIn('cohort input is empty', self.nf(ok=False))

    def test_runtime_failure_preserves_failed_stage_ledgers_for_every_genome(self):
        path = Path(self.rows[0]['genome_layout_path'])
        path.write_text('\n'.join(line for line in path.read_text().splitlines() if not line.startswith('p3\t')) + '\n')
        self.nf(ok=False)
        ledgers = list((self.root / 'work').rglob('*_integrated.tsv.stages.json'))
        self.assertEqual(len(ledgers), 4)
        for p in ledgers:
            stages = {s['name']: s['status'] for s in json.loads(p.read_text())['stages']}
            self.assertEqual(stages['cohort_context'], 'failed')
        self.assertFalse((self.root / 'published/cohort').exists())


if __name__ == '__main__':
    unittest.main()
