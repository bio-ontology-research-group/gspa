import copy
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import cohort_sequence_inputs as C


class MultisourceEvidenceTest(unittest.TestCase):
    def setUp(self):
        base = Path(os.environ.get('TMPDIR', str(Path.home() / 'tmp')))
        base.mkdir(exist_ok=True)
        self.tmp = tempfile.TemporaryDirectory(dir=base)
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.contracts = {}
        for name, char in (('clean', 'c'), ('swissprot_diamond_baseline', 'd')):
            self.contracts[name] = dict(schema='gspa.predictor-contract/1', tool=name,
                tool_version='fixture', tool_sha256=char * 64, reference_sha256=char * 64,
                parameters_sha256=char * 64, adapter_sha256=char * 64, score_model='fixture')
        self.allowed = [dict(source=k, contract=v) for k, v in self.contracts.items()]
        self.diamond = [self.claim('swissprot_diamond_baseline', 'p', 'EC', '1.1.1.1', .4),
                        self.claim('swissprot_diamond_baseline', 'q', 'GO', 'GO:0000001', .6)]
        self.clean = [self.claim('clean', 'p', 'EC', 'EC:1.1.1.1', .9),
                      self.claim('clean', 'r', 'EC', 'EC:3.3.3.n1', 0.0)]
        self.inputs = []
        for name, rows in (('diamond', self.diamond), ('clean', self.clean)):
            p = self.root / (name + '.jsonl')
            self.jsonl(p, rows)
            self.inputs.append(dict(path=str(p), sha256=C.sha(p)))
        self.merged = self.root / 'merged.jsonl'
        self.record = C.merge(self.inputs, 'g', self.allowed, self.merged)
        self.rows = [self.row('p', 'EC', 'EC:1.1.1.1', 'CATALYTIC', {'clean': .9, 'swissprot_diamond_baseline': .4}),
                     self.row('q', 'GO', 'GO:0000001', 'COMPLEX_SUBUNIT', {'swissprot_diamond_baseline': .6}),
                     self.row('r', 'EC', 'EC:3.3.3.n1', 'CATALYTIC', {'clean': 0.0})]
        self.features = self.root / 'features/g'
        self.features.mkdir(parents=True)
        self.export()

    @staticmethod
    def jsonl(path, rows):
        path.write_text(''.join(json.dumps(r) + '\n' for r in rows))

    def claim(self, source, protein, kind, function, score):
        return dict(source=source, protein_id=protein, function_type=kind, function_id=function,
                    raw_score=score, evidence_type='SEQUENCE_DEEPLEARNING' if source == 'clean' else 'SEQUENCE_SIMILARITY',
                    metadata=dict(predictor_contract=self.contracts[source], genome_id='g'))

    def row(self, protein, kind, function, role, sources):
        fields = {'seq.n_sources_reporting': len(sources), 'seq.n_groups_reporting': len(sources),
                  'seq.max_raw_any': max(sources.values())}
        for source, column, group in (('clean', 'clean', 'ml_protein_seq'),
                                      ('swissprot_diamond_baseline', 'other', 'homology')):
            for prefix, name in (('seq.src.', column), ('seq.group.', group)):
                fields[prefix + name + '.reported'] = float(source in sources)
                fields[prefix + name + '.max_raw'] = sources.get(source, 0.0)
        identities = [dict(source=s, contract=self.contracts[s]) for s in sources]
        identities.sort(key=lambda x: json.dumps(x, sort_keys=True, separators=(',', ':')))
        return dict(genome_id='g', protein_id=protein, function_type=kind, function_id=function,
                    function_key='|'.join((protein, kind, function)), candidate_key='|'.join((protein, kind, function, role)),
                    role=role, sequence_sources=sorted(sources), sequence_contracts=identities, features=fields)

    def export(self):
        p = self.features / 'candidate_features.jsonl'
        self.jsonl(p, self.rows)
        (self.features / 'export_manifest.json').write_text(json.dumps(dict(
            schema='gspa.context-feature-export/1', rows_sha256=C.sha(p), rows=len(self.rows), labels_included=False)))

    def audit(self):
        return C.audit_export(self.root / 'features', {'g': self.record}, self.allowed)

    def test_exact_merge_and_candidate_specific_evidence(self):
        self.assertEqual(self.merged.read_bytes(), b''.join(Path(i['path']).read_bytes() for i in self.inputs))
        self.assertEqual((self.record['records'], self.record['candidates']), (4, 3))
        self.assertEqual(self.audit(), {'g': 3})

    def test_real_prefix_difference_is_reported_without_rewriting_inputs(self):
        self.assertEqual(self.record['identity_normalization'], [dict(protein_id='p',function_type='EC',
            function_id='EC:1.1.1.1',raw_function_ids=['1.1.1.1','EC:1.1.1.1'])])
        expected,_=C.profiles(self.merged,'g',self.allowed)
        self.assertEqual(expected[('p','EC','EC:1.1.1.1')]['sources'], {'clean','swissprot_diamond_baseline'})
        self.assertEqual(self.merged.read_bytes(),b''.join(Path(i['path']).read_bytes() for i in self.inputs))

    def test_old_split_ec_candidate_export_is_rejected(self):
        raw=self.row('p','EC','1.1.1.1','CATALYTIC',{'swissprot_diamond_baseline':.4})
        self.rows[0]=self.row('p','EC','EC:1.1.1.1','CATALYTIC',{'clean':.9})
        self.rows.insert(0,raw);self.export()
        with self.assertRaisesRegex(ValueError,'candidate/claim identity differs'):
            self.audit()

    def test_ec_prefix_only_never_resolves_obsolete_or_partial_functions(self):
        for raw,wanted in [('  ec:1.3.5.4  ','EC:1.3.5.4'),('2.8.1.-','EC:2.8.1.-'),
                           ('3.3.3.n1','EC:3.3.3.n1'),('EC:bad','EC:bad'),(' EC:EC:1.1.1.1 ',' EC:EC:1.1.1.1 '),
                           ('\u00a01.1.1.1\u00a0','\u00a01.1.1.1\u00a0')]:
            self.assertEqual(C.function_id('EC',raw),wanted)
        self.assertEqual(C.function_id('GO','1.1.1.1'),'1.1.1.1')

    def test_conflicting_raw_identity_fails_before_jvm_ingestion(self):
        rows=copy.deepcopy(self.diamond)
        rows[0]['metadata']['raw_function_id']='EC:9.9.9.9'
        path=self.root/'conflict.jsonl';self.jsonl(path,rows)
        with self.assertRaisesRegex(ValueError,'contradictory raw_function_id'):
            C.profiles(path,'g',self.allowed)
        rows[0]['metadata']['raw_function_id']='1.1.1.1';self.jsonl(path,rows)
        self.assertEqual(len(C.profiles(path,'g',self.allowed)[0]),2)

    def test_training_export_audit_dispatches_multisource_inputs(self):
        import run_family_excluded_training as F
        self.assertEqual(F.audit_export(self.root / 'features', {'g': self.record}, self.allowed), {'g': 3})
        self.rows[1]['sequence_sources'] = ['clean', 'swissprot_diamond_baseline']
        self.export()
        with self.assertRaisesRegex(ValueError, 'candidate predictor provenance differs'):
            F.audit_export(self.root / 'features', {'g': self.record}, self.allowed)

    def test_zero_score_still_reports_source(self):
        self.rows[2]['features']['seq.src.clean.reported'] = 0.0
        self.export()
        with self.assertRaisesRegex(ValueError, 'raw sequence feature differs'):
            self.audit()

    def test_combined_probability_cannot_replace_raw_score(self):
        self.rows[0]['features']['seq.max_raw_any'] = .94
        self.export()
        with self.assertRaisesRegex(ValueError, 'raw sequence feature differs'):
            self.audit()

    def test_global_source_list_cannot_be_copied_to_single_source_candidate(self):
        self.rows[1]['sequence_contracts'] = copy.deepcopy(self.rows[0]['sequence_contracts'])
        self.rows[1]['sequence_sources'] = self.rows[0]['sequence_sources'][:]
        self.export()
        with self.assertRaisesRegex(ValueError, 'candidate predictor provenance differs'):
            self.audit()

    def test_missing_candidate_is_not_a_correct_abstention(self):
        self.rows.pop()
        self.export()
        with self.assertRaisesRegex(ValueError, 'candidate/claim identities differ'):
            self.audit()

    def test_duplicate_role_rejected(self):
        self.rows.append(copy.deepcopy(self.rows[0]))
        self.export()
        with self.assertRaisesRegex(ValueError, 'duplicate candidate role'):
            self.audit()

    def test_go_membership_not_forced_to_catalytic(self):
        self.rows[1]['role'] = 'CATALYTIC'
        self.rows[1]['candidate_key'] = self.rows[1]['function_key'] + '|CATALYTIC'
        self.export()
        with self.assertRaisesRegex(ValueError, 'unexpected candidate role'):
            self.audit()

    def test_supplied_distinct_roles_remain_separate(self):
        rows = copy.deepcopy(self.diamond + self.clean)
        rows[0]['metadata']['candidate_role'] = 'CATALYTIC'
        rows[2]['metadata']['candidate_role'] = 'ACCESSORY'
        self.jsonl(self.merged, rows)
        self.record['sha256'] = C.sha(self.merged)
        extra = copy.deepcopy(self.rows[0])
        extra['role'] = 'ACCESSORY'
        extra['candidate_key'] = extra['function_key'] + '|ACCESSORY'
        self.rows.append(extra)
        self.export()
        self.assertEqual(self.audit(), {'g': 4})
        self.rows.pop()
        self.export()
        with self.assertRaisesRegex(ValueError, 'supplied candidate role missing'):
            self.audit()

    def test_changed_reference_contract_rejected(self):
        self.rows[0] = copy.deepcopy(self.rows[0])
        self.rows[0]['sequence_contracts'][0]['contract']['reference_sha256'] = 'e' * 64
        self.export()
        with self.assertRaisesRegex(ValueError, 'candidate predictor provenance differs'):
            self.audit()

    def test_claims_cannot_cross_genome_boundaries(self):
        with self.assertRaisesRegex(ValueError, 'claim genome mismatch'):
            C.merge(self.inputs, 'another-genome', self.allowed, self.root / 'bad.jsonl')
        self.assertFalse((self.root / 'bad.jsonl').exists())

    def test_changed_source_rejected_before_publication(self):
        Path(self.inputs[0]['path']).write_text('changed')
        with self.assertRaisesRegex(ValueError, 'source claim hash changed'):
            C.merge(self.inputs, 'g', self.allowed, self.root / 'bad.jsonl')
        self.assertFalse((self.root / 'bad.jsonl').exists())

    def test_previous_output_not_overwritten(self):
        before = self.merged.read_bytes()
        with self.assertRaises(FileExistsError):
            C.merge(self.inputs, 'g', self.allowed, self.merged)
        self.assertEqual(self.merged.read_bytes(), before)

    def test_empty_clean_not_treated_as_successful_observation(self):
        p = self.root / 'empty.jsonl'
        p.write_text('')
        inputs = [self.inputs[0], dict(path=str(p), sha256=C.sha(p))]
        with self.assertRaisesRegex(ValueError, 'CLEAN observation is empty'):
            C.merge(inputs, 'g', self.allowed, self.root / 'bad.jsonl')
        self.assertFalse((self.root / 'bad.jsonl').exists())

    def test_empty_diamond_is_valid_no_hit_observation(self):
        p = self.root / 'empty.jsonl'
        p.write_text('')
        inputs = [dict(path=str(p), sha256=C.sha(p)), self.inputs[1]]
        result = C.merge(inputs, 'g', self.allowed, self.root / 'valid.jsonl')
        self.assertEqual(result['by_source'], {'clean': 2})

    def test_generic_pipe_ids_preserved_but_loader_rewrites_rejected(self):
        for protein in ('custom|p', 'sp|P123|protein'):
            rows = copy.deepcopy(self.clean)
            rows[0]['protein_id'] = protein
            p = self.root / 'ids.jsonl'
            self.jsonl(p, rows)
            inputs = [self.inputs[0], dict(path=str(p), sha256=C.sha(p))]
            if protein.startswith('sp|'):
                with self.assertRaisesRegex(ValueError, 'would be rewritten'):
                    C.merge(inputs, 'g', self.allowed, self.root / 'bad.jsonl')
            else:
                result = C.merge(inputs, 'g', self.allowed, self.root / 'pipe.jsonl')
                self.assertEqual(result['records'], 4)

    def test_candidate_key_cannot_name_a_different_protein(self):
        self.rows[0]['candidate_key'] = 'other|EC|EC:1.1.1.1|CATALYTIC'
        self.export()
        with self.assertRaisesRegex(ValueError, 'candidate key encoding differs'):
            self.audit()

    def test_same_go_term_has_consistent_unsupplied_role(self):
        rows = copy.deepcopy(self.diamond + self.clean)
        extra = copy.deepcopy(self.diamond[1])
        extra['protein_id'] = 'q2'
        rows.append(extra)
        self.jsonl(self.merged, rows)
        self.record['sha256'] = C.sha(self.merged)
        self.rows.append(self.row('q2', 'GO', 'GO:0000001', 'UNSPECIFIED', {'swissprot_diamond_baseline': .6}))
        self.export()
        with self.assertRaisesRegex(ValueError, 'inconsistent GO default role'):
            self.audit()


if __name__ == '__main__':
    unittest.main()
