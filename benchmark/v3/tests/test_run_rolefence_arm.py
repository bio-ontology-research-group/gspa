import copy
import csv
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import run_rolefence_arm as M

F = M.F
G = M.G
R = M.R
S = M.S


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')


def write_jsonl(path, rows):
    Path(path).write_text(''.join(json.dumps(row) + '\n' for row in rows))


class TempCase(unittest.TestCase):
    def setUp(self):
        temporary = Path(os.environ.get('TMPDIR', str(Path.home() / 'tmp')))
        temporary.mkdir(parents=True, exist_ok=True)
        self.temp = tempfile.TemporaryDirectory(dir=temporary)
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)


class PrepareTest(TempCase):
    def setUp(self):
        super().setUp()
        self.genomes = ['GCA_%09d.1' % index for index in range(18)]
        self.contract = dict(
            schema='gspa.predictor-contract/1', tool='diamond', tool_version='fixture',
            tool_sha256='a' * 64, reference_sha256='b' * 64,
            parameters_sha256='c' * 64, adapter_sha256='d' * 64,
            score_model='fixture')
        self.claims_dir = self.root / 'claims'
        self.claims_dir.mkdir()
        genome_records = {}
        for gid in self.genomes:
            path = self.claims_dir / (gid + '.claims.jsonl')
            row = dict(source='swissprot_diamond_baseline', protein_id='p',
                       function_type='EC', function_id='EC:1.1.1.1', raw_score=0.5,
                       evidence_type='SEQUENCE_SIMILARITY',
                       metadata=dict(predictor_contract=self.contract, genome_id=gid))
            write_jsonl(path, [row])
            genome_records[gid] = dict(file=str(path), sha256=S.sha(path), records=1)
        self.report = dict(
            schema='gspa.role-genome-claims/1', mode='homology',
            genomes=genome_records, contract=self.contract,
            contract_differs_from_cohort=['reference_sha256'],
            cohort_contract=dict(reference_sha256='0' * 64),
            reference=dict(sequence_sha256='e' * 64, annotation_sha256='f' * 64))
        self.claims_report = self.claims_dir / 'claims_report.json'
        write_json(self.claims_report, self.report)
        self.training = self.root / 'training'
        self.context = self.root / 'context'
        self.context.mkdir()
        supervision = {}
        for name in ('labels', 'splits', 'label_policy', 'registry', 'go_obo',
                     'example_clusters'):
            path = self.root / (name + '.input')
            path.write_text(name + '\n')
            supervision[name] = str(path)
        recipe = dict(supervision=supervision, disjoint_keys=['genome_id'],
                      settings={kind: {} for kind in F.KINDS})
        self.fitting = dict(recipes={mode: copy.deepcopy(recipe) for mode in G.MODES},
                            comparison='comparison')
        self.training.mkdir()
        write_json(self.training / 'registration.json',
                   dict(inputs=F.hashes(supervision.values())))
        manifest = self.root / 'baseline.tsv'
        with manifest.open('w') as stream:
            writer = csv.DictWriter(stream, fieldnames=R.MANIFEST_COLUMNS,
                                    delimiter='\t', lineterminator='\n')
            writer.writeheader()
            for gid in self.genomes:
                writer.writerow(dict(
                    genome_id=gid, claims_path='old.claims', out_path='old.out',
                    taxon='bacteria', essential_profile='bacteria', operons_path='-',
                    gaps_path=str(self.root / (gid + '.gaps')), max_iter='-',
                    enforce_stage2='false', enforce_mode='-', enable_priors='none',
                    genome_layout_path='-', orthogroups_path='-',
                    cluster_consensus_path='-', fba_medium='-', fba_gram='-'))
        command = ['java', '-jar', 'old.jar', 'integrate', '--batch-manifest', str(manifest),
                   '--feature-export-dir', 'old-features', '--batch-stats-out', 'old-stats',
                   '--features-out', 'old-rlc', '--suggestions-out', 'old-suggestions',
                   '--cohort', '--rlc-suggester', '--cohort-context-mode', 'UNIVERSAL',
                   '--cohort-pathways', 'old-pathways', '--cohort-metadata', 'old-metadata',
                   '--enable-priors', 'all']
        self.base = dict(
            genomes=self.genomes,
            gaps={gid: dict(path=str(self.root / (gid + '.gaps'))) for gid in self.genomes},
            plans={mode: dict(observed_gaps=dict(manifest=str(manifest), command=command))
                   for mode in G.MODES})

    def prepare(self, name='out'):
        with patch.object(M.R, 'base_inputs',
                          return_value=(self.fitting, self.base, F.hashes([self.root / 'baseline.tsv']))), \
             patch.object(M.F.V, 'describe_run', return_value={'inputs': {}}):
            return M.prepare(self.training, self.context, self.root / 'gspa.jar',
                             self.claims_dir, self.claims_report, self.root / name)

    def test_prepare_registers_eighteen_claims_and_rolefence_paths(self):
        protocol = self.prepare()
        self.assertEqual(protocol['schema'], M.SCHEMA)
        self.assertEqual(protocol['arms'], [M.ARM])
        self.assertEqual(len(protocol['claims'][M.ARM]), 18)
        self.assertEqual(protocol['rolefence_reference']['predictor_contract'], self.contract)
        for ctx in M.CONTEXTS:
            plan = protocol['plans'][M.ARM][ctx]
            with Path(plan['manifest']).open() as stream:
                rows = list(csv.DictReader(stream, delimiter='\t'))
            self.assertEqual(len(rows), 18)
            for row in rows:
                expected = self.claims_dir / (row['genome_id'] + '.claims.jsonl')
                self.assertEqual(row['claims_path'], str(expected.resolve()))
        self.assertEqual(M.verify(self.root / 'out'), protocol)
        registration = json.loads((self.root / 'out/registration.json').read_text())
        for gid in self.genomes:
            path = self.claims_dir / (gid + '.claims.jsonl')
            self.assertEqual(registration['inputs'][str(path.resolve())], S.sha(path))

    def test_prepare_rejects_wrong_record_contract(self):
        gid = self.genomes[0]
        path = self.claims_dir / (gid + '.claims.jsonl')
        row = json.loads(path.read_text())
        row['metadata']['predictor_contract']['reference_sha256'] = '9' * 64
        write_jsonl(path, [row])
        self.report['genomes'][gid]['sha256'] = S.sha(path)
        write_json(self.claims_report, self.report)
        with self.assertRaisesRegex(ValueError, 'unregistered predictor identity'):
            self.prepare('wrong-contract')

    def test_prepare_rejects_claim_sha_mismatch(self):
        path = self.claims_dir / (self.genomes[0] + '.claims.jsonl')
        path.write_text(path.read_text() + path.read_text())
        with self.assertRaisesRegex(ValueError, 'claim hash differs'):
            self.prepare('wrong-sha')

    def test_verify_detects_claim_change_after_registration(self):
        self.prepare()
        path = self.claims_dir / (self.genomes[0] + '.claims.jsonl')
        path.write_text(path.read_text() + path.read_text())
        with self.assertRaisesRegex(ValueError, 'inputs or outputs changed'):
            M.verify(self.root / 'out')


class StageTest(TempCase):
    def setUp(self):
        super().setUp()
        self.out = self.root / 'run'
        self.arm = self.out / M.ARM
        (self.arm / 'SEQONLY').mkdir(parents=True)
        recipe = dict(
            supervision={name: 'input' for name in
                         ('labels', 'splits', 'label_policy', 'registry', 'go_obo',
                          'example_clusters')},
            original_table='original', disjoint_keys=['genome_id'],
            settings={kind: {} for kind in F.KINDS})
        self.protocol = dict(
            recipes={mode: copy.deepcopy(recipe) for mode in G.MODES},
            pending_stages=['full_reference_ec_endpoint', 'production_replay'],
            limitations='fixture', genomes=['g'], gaps={},
            claims={M.ARM: {'g': {}}}, contracts={M.ARM: []},
            plans={M.ARM: {'SEQONLY': dict(
                command=['x', '--feature-export-dir', str(self.arm / 'SEQONLY/features')],
                manifest='manifest', inference='inference',
                features=str(self.arm / 'SEQONLY/features'))}},
            zero_audit=dict(prefixes=list(R.ZERO_PREFIXES)))

    comparison = dict(shared_examples=5, original_only=3, new_only=2,
                      representation_changes=1)

    def test_one_context_export_train_evaluate_writes_comparison_shaped_files(self):
        rows = [dict(split=split) for split in F.T.SPLITS]

        def batch(command, manifest, inference, sources):
            features = Path(command[command.index('--feature-export-dir') + 1]) / 'g'
            features.mkdir(parents=True)
            (features / 'feature_contract.json').write_text('{}')

        def build_table(*args, **kwargs):
            destination = args[5]
            destination.mkdir()
            (destination / 'training_table.jsonl').write_text('{}\n')
            (destination / 'training_table_manifest.json').write_text('{}')
            return dict(table=str(destination / 'training_table.jsonl'))

        def fit(train, calibration, contract, kind, **kwargs):
            self.assertEqual(train, [dict(split='train')])
            self.assertEqual(calibration, [dict(split='calibration')])
            self.assertFalse((self.out / 'training.json').exists())
            return dict(kind=kind)

        with patch.object(M, 'CONTEXTS', ('SEQONLY',)), \
             patch.object(M, 'verify', return_value=self.protocol), \
             patch.object(M.F.V, 'run_batch', side_effect=batch) as run_batch, \
             patch.object(M.S, 'audit_export', return_value={'g': 1}) as audit, \
             patch.object(M.R, 'zero_audit', return_value=4), \
             patch.object(M.R, 'sequence_stage_audit', return_value={'g': []}), \
             patch.object(M.F.P, 'build_table', side_effect=build_table), \
             patch.object(M.R, 'compare_supervision', return_value=self.comparison), \
             patch.object(M.F.T, 'load_contract', return_value={}), \
             patch.object(M.F.T, 'read_table', return_value=rows), \
             patch.object(M.F.T, 'check_disjoint'), \
             patch.object(M.F.T, 'split_identity', return_value={}), \
             patch.object(M.F.T, 'train', side_effect=fit) as fit_model, \
             patch.object(M.F.T, 'write_frozen',
                          side_effect=lambda artifact, path: write_json(path, artifact)), \
             patch.object(M.F.T, 'evaluate_frozen', return_value={'rows': 1}) as evaluate:
            M.export(self.out)
            self.assertEqual(run_batch.call_count, 1)
            self.assertEqual(audit.call_count, 1)
            M.train(self.out)
            self.assertEqual(fit_model.call_count, 2)
            evaluate.assert_not_called()
            M.evaluate(self.out)
            self.assertEqual(evaluate.call_count, 2)

        context = self.arm / 'SEQONLY'
        for relative in ('features/g/feature_contract.json', 'table/training_table.jsonl',
                         'table/training_table_manifest.json', 'models/model_logistic.json',
                         'models/model_gbdt.json'):
            self.assertTrue((context / relative).is_file(), relative)
        for name in ('export.json', 'training.json', 'evaluation.json'):
            self.assertTrue((self.out / name).is_file(), name)
        training = json.loads((self.out / 'training.json').read_text())
        self.assertEqual(training['counts'][M.ARM]['SEQONLY']['supervision_comparison'],
                         self.comparison)
        self.assertEqual(training['supervision_comparison'],
                         {M.ARM: {'SEQONLY': self.comparison}})
        exported = json.loads((self.out / 'export.json').read_text())
        self.assertEqual(exported['summaries'][M.ARM]['SEQONLY']['candidates'], 1)
        evaluation = json.loads((self.out / 'evaluation.json').read_text())
        self.assertEqual(evaluation['results'][M.ARM]['SEQONLY'],
                         {kind: {'rows': 1} for kind in F.KINDS})
        with self.assertRaises(FileExistsError), \
             patch.object(M, 'CONTEXTS', ('SEQONLY',)), \
             patch.object(M, 'verify', return_value=self.protocol):
            M.export(self.out)

    def test_train_fails_closed_without_shared_examples(self):
        (self.arm / 'SEQONLY/features/g').mkdir(parents=True)
        (self.arm / 'SEQONLY/features/g/feature_contract.json').write_text('{}')
        write_json(self.out / 'export.json', dict(status='verified', outputs={}))

        def build_table(*args, **kwargs):
            destination = args[5]
            destination.mkdir()
            (destination / 'training_table.jsonl').write_text('{}\n')
            return dict(table=str(destination / 'training_table.jsonl'))

        empty = dict(self.comparison, shared_examples=0)
        with patch.object(M, 'CONTEXTS', ('SEQONLY',)), \
             patch.object(M, 'verify', return_value=self.protocol), \
             patch.object(M.F.P, 'build_table', side_effect=build_table), \
             patch.object(M.R, 'compare_supervision', return_value=empty), \
             patch.object(M.F.T, 'train') as fit_model:
            with self.assertRaisesRegex(ValueError, 'no shared supervised examples'):
                M.train(self.out)
        fit_model.assert_not_called()
        self.assertFalse((self.out / 'training.json').exists())


if __name__ == '__main__':
    unittest.main()
