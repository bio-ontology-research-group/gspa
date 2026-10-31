import csv
import copy
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import audit_clean_cohort as A
import observe_clean_cohort as O


class CollectionAuditTest(unittest.TestCase):
    def setUp(self):
        base = Path(os.environ.get('TMPDIR', str(Path.home() / 'tmp')))
        base.mkdir(exist_ok=True)
        self.temp = tempfile.TemporaryDirectory(dir=base)
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.folder = self.root / 'job-1'
        self.folder.mkdir()
        fasta = self.root / 'g.faa'
        fasta.write_text('>p\nMPEPTIDE\n>q\nMPROTEIN\n')
        self.g = dict(genome_id='g', proteome=str(fasta), sha256=A.sha(fasta), proteins=2)
        self.plan = dict(genomes=[self.g], parameters=dict(runtime_packages={'torch': 'fixture'}),
                         python='/python', sidecar='/sidecar',
                         native_inputs={'g': self.g}, unsupported_inputs={'g': []},
                         staged_app_files={'CLEAN_infer_fasta.py': 'a' * 64},
                         predictor_contract={'schema': 'gspa.predictor-contract/1', 'tool': 'fixture'})
        self.plan_sha = 'b' * 64
        p = self.folder / 'predictions'
        p.mkdir()
        native = p / 'g.clean.native.csv'
        prefix = 'gspa_' + '0' * 32 + '_p'
        native.write_text(prefix + '00000000,EC:1.1.1.1/0.0000,EC:2.2.2.2/0.3456\n' +
                          prefix + '00000001,EC:3.3.3.n1/1.0000\n')
        O.C.rewrite(native, p / 'g.clean.tsv', 0, score_kind='gmm_confidence',
                    protein_ids={prefix + '00000000': 'p', prefix + '00000001': 'q'},
                    provenance=dict(driver_sha256='a' * 64, driver_score_kind='gmm_confidence',
                                    input_fasta_sha256=self.g['sha256']))
        summary = O.audit_native(p, self.g, self.plan['predictor_contract'],
                                 self.folder / 'claims.jsonl', driver_sha256='a' * 64)
        O.save(self.folder / 'started.json', dict(genome=self.g, index=0,
               plan_sha256=self.plan_sha, job_id='1'))
        O.save(self.folder / 'runtime.json', dict(cuda='fixture', device='fixture',
               packages=self.plan['parameters']['runtime_packages']))
        O.save(self.folder / 'command.json', ['/python', '/sidecar', '--predictor', 'clean', '--manifest',
               str(self.folder / 'manifest.tsv'), '--model-dir', '/scratch/1/clean-app',
               '--clean-score-kind', 'gmm_confidence', '--min-score', '0', '--seed', '0'])
        (self.folder / 'manifest.tsv').write_text('tag\tfasta_path\toutput_dir\n' +
               'g\t' + str(fasta) + '\t' + str(p) + '\n')
        self.done = dict(status='verified', genome=self.g, plan_sha256=self.plan_sha, audit=summary,
                         native_genome=self.g, unsupported_inputs=[])
        self.reseal()

    def reseal(self):
        self.done['artifacts'] = {str(p.relative_to(self.folder)): A.sha(p)
                                 for p in self.folder.rglob('*')
                                 if p.is_file() and p.name != 'completed.json'}
        (self.folder / 'completed.json').write_text(json.dumps(self.done))

    def verify(self):
        return A.verify_attempt(self.folder, self.g, self.plan, self.plan_sha)

    def test_full_native_join_including_zero_and_provisional(self):
        result = self.verify()
        self.assertEqual((result['proteins'], result['calls'], result['zero_confidence_calls']), (2, 3, 1))

    def test_claim_mutation_rejected_even_with_resealed_artifacts(self):
        p = self.folder / 'claims.jsonl'
        rows = [json.loads(line) for line in p.read_text().splitlines()]
        rows[0]['raw_score'] = .25
        p.write_text(''.join(json.dumps(row) + '\n' for row in rows))
        self.reseal()
        with self.assertRaisesRegex(ValueError, 'claim/native mismatch'):
            self.verify()

    def test_calibration_invention_rejected(self):
        p = self.folder / 'claims.jsonl'
        rows = [json.loads(line) for line in p.read_text().splitlines()]
        rows[0]['calibrated_prob'] = .25
        p.write_text(''.join(json.dumps(row) + '\n' for row in rows))
        self.reseal()
        with self.assertRaisesRegex(ValueError, 'claim/native mismatch'):
            self.verify()

    def test_typed_mutation_rejected_after_receipt_resealed(self):
        p = self.folder / 'predictions/g.clean.tsv'
        with p.open() as stream:
            reader = csv.DictReader(stream, delimiter='\t')
            header = reader.fieldnames
            rows = list(reader)
        rows[0]['score'] = '0.25'
        with p.open('w') as stream:
            writer = csv.DictWriter(stream, fieldnames=header, delimiter='\t')
            writer.writeheader()
            writer.writerows(rows)
        receipt = p.with_suffix('.receipt.json')
        value = A.read(receipt)
        value['output_sha256'] = A.sha(p)
        receipt.write_text(json.dumps(value))
        self.reseal()
        with self.assertRaisesRegex(ValueError, 'typed score mismatch'):
            self.verify()

    def test_unexpected_metadata_rejected_even_when_claim_matches(self):
        p = self.folder / 'predictions/g.clean.tsv'
        with p.open() as stream:
            reader = csv.DictReader(stream, delimiter='\t')
            header, rows = reader.fieldnames, list(reader)
        metadata = json.loads(rows[0]['metadata_json'])
        metadata['calibrated_probability'] = .97
        rows[0]['metadata_json'] = json.dumps(metadata)
        with p.open('w') as stream:
            writer = csv.DictWriter(stream, fieldnames=header, delimiter='\t')
            writer.writeheader()
            writer.writerows(rows)
        receipt = p.with_suffix('.receipt.json')
        value = A.read(receipt)
        value['output_sha256'] = A.sha(p)
        receipt.write_text(json.dumps(value))
        claims = self.folder / 'claims.jsonl'
        data = [json.loads(line) for line in claims.read_text().splitlines()]
        data[0]['metadata']['calibrated_probability'] = .97
        claims.write_text(''.join(json.dumps(row) + '\n' for row in data))
        self.reseal()
        with self.assertRaisesRegex(ValueError, 'unexpected typed metadata'):
            self.verify()

    def test_changed_command_rejected_even_when_resealed(self):
        p = self.folder / 'command.json'
        command = A.read(p)
        command[-1] = '9'
        p.write_text(json.dumps(command))
        self.reseal()
        with self.assertRaisesRegex(ValueError, 'native command mismatch'):
            self.verify()

    def test_boolean_is_not_native_confidence(self):
        p = self.folder / 'claims.jsonl'
        rows = [json.loads(line) for line in p.read_text().splitlines()]
        rows[0]['raw_score'] = False
        p.write_text(''.join(json.dumps(row) + '\n' for row in rows))
        self.reseal()
        with self.assertRaisesRegex(ValueError, 'claim/native mismatch'):
            self.verify()

    def test_duplicate_and_nonfinite_json_rejected(self):
        for text in ('{"a":1,"a":2}', '{"a":NaN}', '{"a":Infinity}'):
            with self.subTest(text=text), self.assertRaises(ValueError):
                A.decode(text)

    def test_unregistered_file_rejected(self):
        (self.folder / 'extra.json').write_text('{}')
        with self.assertRaisesRegex(ValueError, 'artifact inventory/hash mismatch'):
            self.verify()

    def test_wrong_driver_rejected(self):
        self.plan['staged_app_files']['CLEAN_infer_fasta.py'] = 'c' * 64
        with self.assertRaisesRegex(ValueError, 'native receipt mismatch'):
            self.verify()

    def collection_shape(self):
        plan = dict(schema='gspa.clean-cohort/2', evaluation_labels_read=False,
                    genomes=[dict(genome_id='g' + str(i), proteins=2748 if i < 17 else 2747)
                             for i in range(18)])
        path = self.root / 'plan.json'
        path.write_text(json.dumps(plan))
        attempts = self.root / 'attempts'
        attempts.mkdir()
        for i in range(18):
            (attempts / ('g' + str(i))).mkdir()
        return path, attempts

    def test_missing_whole_genome_rejected_before_scoring(self):
        path, attempts = self.collection_shape()
        (attempts / 'g17').rmdir()
        with self.assertRaisesRegex(ValueError, 'missing/unknown genome attempts'):
            A.audit(self.root, A.sha(path))

    def test_multiple_successes_cannot_be_selected_by_score(self):
        path, attempts = self.collection_shape()
        for name in ('job-1', 'job-2'):
            folder = attempts / 'g0' / name
            folder.mkdir()
            (folder / 'completed.json').write_text('{}')
        with self.assertRaisesRegex(ValueError, 'exactly one completed attempt'):
            A.audit(self.root, A.sha(path))

    def test_unfinished_attempt_is_not_ignored(self):
        path, attempts = self.collection_shape()
        (attempts / 'g0' / 'job-1').mkdir()
        with self.assertRaisesRegex(ValueError, 'unfinished attempt'):
            A.audit(self.root, A.sha(path))

    def test_failed_attempt_cannot_also_be_completed(self):
        (self.folder / 'failed.json').write_text('{}')
        with self.assertRaisesRegex(ValueError, 'also marked failed'):
            self.verify()

    def native_plan(self):
        app = self.root / 'app'
        source = self.root / 'source/benchmark'
        home = self.root / 'torch'
        files = {'CLEAN_infer_fasta.py': "infer_maxsep('split100', args.fasta_data, gmm='data/pretrained/gmm_ensumble.pkl')\n",
                 'esm/scripts/extract.py': '# fixture', 'src/CLEAN/infer.py': '# fixture',
                 'data/split100.csv': 'Entry\tSequence\nx\tMPEPTIDE\n',
                 'data/pretrained/split100.pth': 'fixture model',
                 'data/pretrained/100.pt': 'fixture reference',
                 'data/pretrained/gmm_ensumble.pkl': 'fixture gmm'}
        inputs, staged = {}, {}
        for name, value in files.items():
            p = app / name
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_text(value)
            inputs[str(p)] = staged[name] = A.sha(p)
        for name in ('esm1b_t33_650M_UR50S.pt', 'esm1b_t33_650M_UR50S-contact-regression.pt'):
            p = home / 'hub/checkpoints' / name
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_text('fixture weight')
            inputs[str(p)] = A.sha(p)
        adapters = {}
        for name in ('v3/observe_clean_cohort.py', 'neural/clean_contract.py', 'neural/run_neural_predictors.py'):
            p = source / name
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_text('# fixture adapter')
            inputs[str(p)] = adapters[p.name] = A.sha(p)
        parent = self.root / 'cohort.json'
        parent.write_text(json.dumps(dict(schema='gspa.cohort-pathway-search/1', genomes=[self.g])))
        inputs[str(parent)] = A.sha(parent)
        params = dict(mode='gmm_confidence', min_score=0, seed=0, native_selection='maxsep',
                      esm_maximum_residues=1022, runtime_packages={},
                      unsupported_input_policy='abstain_without_sequence_changes/v1')
        tool = {k: v for k, v in staged.items() if k.endswith('.py')}
        refs = {k: v for k, v in staged.items() if not k.endswith('.py')}
        refs.update({Path(k).name: v for k, v in inputs.items() if '/hub/checkpoints/' in k})
        contract = dict(schema='gspa.predictor-contract/1', tool='CLEAN',
            tool_version='installed-source-' + O.digest(tool)[:16], tool_sha256=O.digest(tool),
            reference_sha256=O.digest(refs), parameters_sha256=O.digest(params),
            adapter_sha256=O.digest(adapters), score_model='native_maxsep_gmm_confidence_unfiltered_round4/v1')
        inputs[self.g['proteome']] = self.g['sha256']
        return dict(inputs=inputs, parent_cohort_sha256=A.sha(parent), genomes=[self.g],
                    native_inputs={'g': self.g}, unsupported_inputs={'g': []},
                    native_alphabet=sorted(set('LAGVSERTIDPKQNFYMHWCXBUZO.-')),
                    app=str(app), staged_app_files=staged, torch_home=str(home),
                    sidecar=str(source / 'neural/run_neural_predictors.py'),
                    predictor_contract=contract, parameters=params)

    def test_independent_contract_reconstruction(self):
        plan = self.native_plan()
        A.verify_plan(plan)
        plan['predictor_contract']['reference_sha256'] = 'f' * 64
        with self.assertRaisesRegex(ValueError, 'contract digest mismatch'):
            A.verify_plan(plan)

    def test_parent_population_must_match(self):
        plan = self.native_plan()
        plan['genomes'] = [dict(self.g, proteins=3)]
        with self.assertRaisesRegex(ValueError, 'parent genome population mismatch'):
            A.verify_plan(plan)

    def test_driver_must_explicitly_request_gmm(self):
        plan = self.native_plan()
        (Path(plan['app']) / 'CLEAN_infer_fasta.py').write_text("infer_maxsep('split100', args.fasta_data)\n")
        with self.assertRaisesRegex(ValueError, 'does not request GMM'):
            A.verify_plan(plan)

    def partition_plan(self):
        plan = self.native_plan()
        f = Path(self.g['proteome'])
        f.write_text(f.read_text() + '>unsupported\nMPJPTIDE\n')
        self.g.update(proteins=3, sha256=A.sha(f))
        plan['inputs'][str(f)] = self.g['sha256']
        parent = self.root / 'cohort.json'
        parent.write_text(json.dumps(dict(schema='gspa.cohort-pathway-search/1', genomes=plan['genomes'])))
        plan['inputs'][str(parent)] = plan['parent_cohort_sha256'] = A.sha(parent)
        native, failed = O.native_inputs(plan['genomes'], set(plan['native_alphabet']), self.root / 'native')
        plan['native_inputs'] = native
        plan['unsupported_inputs'] = failed
        plan['inputs'][native['g']['proteome']] = native['g']['sha256']
        return plan

    def test_partition_keeps_denominator_and_requires_abstention_marker(self):
        plan = self.partition_plan()
        A.verify_plan(plan)
        self.assertEqual(plan['genomes'][0]['proteins'], 3)
        self.assertEqual(plan['native_inputs']['g']['proteins'], 2)
        plan['unsupported_inputs']['g'] = []
        with self.assertRaisesRegex(ValueError, 'unsupported input outcomes differ'):
            A.verify_plan(plan)

    def test_partition_cannot_rewrite_valid_sequences_even_with_updated_hashes(self):
        plan = self.partition_plan()
        native = plan['native_inputs']['g']
        f = Path(native['proteome'])
        f.write_text(f.read_text().replace('MPEPTIDE', 'MQEPTIDE'))
        native['sha256'] = plan['inputs'][str(f)] = A.sha(f)
        with self.assertRaisesRegex(ValueError, 'native input sequences changed or omitted'):
            A.verify_plan(plan)

    def test_partitioned_attempt_audits_native_calls_and_original_denominator(self):
        original = self.root / 'full.faa'
        original.write_text(Path(self.g['proteome']).read_text() + '>unsupported\nMPJPTIDE\n')
        self.g = dict(self.g, proteome=str(original), sha256=A.sha(original), proteins=3)
        self.plan['genomes'] = [self.g]
        _, failures = O.native_inputs([self.g], set('LAGVSERTIDPKQNFYMHWCXBUZO.-'), self.root / 'partition')
        self.plan['unsupported_inputs'] = failures
        self.done.update(genome=self.g, unsupported_inputs=failures['g'])
        start = self.folder / 'started.json'
        data = A.read(start)
        data['genome'] = self.g
        start.write_text(json.dumps(data))
        self.reseal()
        result = self.verify()
        self.assertEqual((result['original_proteins'], result['proteins'], len(result['unsupported_inputs'])), (3, 2, 1))
        self.assertEqual(result['original_fasta_sha256'], A.sha(original))
        good = copy.deepcopy(self.done)
        for mutation in ('abstention', 'native', 'count'):
            self.done = copy.deepcopy(good)
            if mutation == 'abstention':self.done['unsupported_inputs'] = []
            if mutation == 'native':self.done['native_genome'] = self.g
            if mutation == 'count':self.done['audit']['proteins'] = 3
            self.reseal()
            with self.subTest(mutation=mutation), self.assertRaises(ValueError):
                self.verify()


if __name__ == '__main__':
    unittest.main()
