import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import run_family_excluded_training as F


class FamilyExcludedTrainingTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(dir=os.environ.get('TMPDIR', str(Path.home()/'tmp')))
        self.root = Path(self.temp.name)
        self.addCleanup(self.temp.cleanup)

    def fixture(self):
        claims = self.root/'claims.jsonl'
        claims.write_text(json.dumps(dict(protein_id='p',function_type='EC',function_id='EC:1.1.1.1'))+'\n')
        feature = self.root/'features/g/candidate_features.jsonl';feature.parent.mkdir(parents=True)
        expected = [dict(source=F.SOURCE,contract=dict(tool='diamond',parameters_sha256='family'))]
        row = dict(genome_id='g',protein_id='p',function_type='EC',function_id='EC:1.1.1.1',sequence_contracts=expected)
        feature.write_text(json.dumps(row)+'\n')
        return claims,feature,row,expected

    def test_complete_candidate_identity_and_correct_mode_are_required(self):
        claims, feature, row, expected = self.fixture()
        files = {'g':{'path':str(claims)}}
        self.assertEqual(F.audit_export(self.root/'features', files, expected), {'g':1})
        row['function_id']='EC:2.2.2.2';feature.write_text(json.dumps(row)+'\n')
        with self.assertRaisesRegex(ValueError,'identities differ'):
            F.audit_export(self.root/'features', files, expected)
        row['function_id']='EC:1.1.1.1';row['sequence_contracts']=[dict(source=F.SOURCE,contract=dict(tool='diamond',parameters_sha256='homology'))]
        feature.write_text(json.dumps(row)+'\n')
        with self.assertRaisesRegex(ValueError,'wrong predictor identity'):
            F.audit_export(self.root/'features', files, expected)

    def test_missing_genome_and_duplicate_candidate_do_not_pass_counts(self):
        claims, feature, row, expected = self.fixture()
        with self.assertRaisesRegex(ValueError,'genome coverage'):
            F.audit_export(self.root/'features', {'g':{'path':str(claims)},'missing':{'path':str(claims)}}, expected)
        feature.write_text((json.dumps(row)+'\n')*2)
        with self.assertRaisesRegex(ValueError,'duplicate candidate'):
            F.audit_export(self.root/'features', {'g':{'path':str(claims)}}, expected)

    def test_two_versions_of_one_source_keep_uniform_identity_check(self):
        claims, feature, row, expected = self.fixture()
        expected.append(dict(source=F.SOURCE, contract=dict(tool='diamond', parameters_sha256='another')))
        row['sequence_contracts'] = expected
        feature.write_text(json.dumps(row)+'\n')
        self.assertEqual(F.audit_export(self.root/'features', {'g':{'path':str(claims)}}, expected), {'g':1})
        row['sequence_contracts'] = expected[:1]
        feature.write_text(json.dumps(row)+'\n')
        with self.assertRaisesRegex(ValueError,'wrong predictor identity'):
            F.audit_export(self.root/'features', {'g':{'path':str(claims)}}, expected)

    def test_changed_protocol_or_registered_input_invalidates_later_stages(self):
        source=self.root/'source';source.write_text('original')
        protocol=dict(mode=F.MODE,model_kinds=list(F.KINDS),numpy_version=F.T.np.__version__,python_version=sys.version)
        F.save(self.root/'protocol.json',protocol)
        registration=dict(inputs=F.hashes([source]),protocol_sha256=F.T.sha256_file(self.root/'protocol.json'))
        F.save(self.root/'registration.json',registration)
        self.assertEqual(F.verify_protocol(self.root),protocol)
        source.write_text('changed')
        with self.assertRaises(Exception):F.verify_protocol(self.root)
        source.write_text('original');protocol['mode']='homology';F.save(self.root/'protocol.json',protocol)
        with self.assertRaisesRegex(ValueError,'protocol changed'):F.verify_protocol(self.root)

    def test_settings_copy_selection_rules_not_fitted_or_test_results(self):
        model=dict(training=dict(weighting='group',cv=dict(folds=4,seed=17)),
                   display_threshold=dict(rule='max_f1',value=.2),test=dict(f1=.99),logistic=dict(weights=[99]))
        self.assertEqual(F.training_settings(model),dict(weighting='group',folds=4,seed=17,threshold_rule='max_f1'))
        model['display_threshold']['value']=.9;model['test']['f1']=0
        self.assertEqual(F.training_settings(model)['threshold_rule'],'max_f1')
        self.assertNotIn('test',F.training_settings(model))

    def test_training_receives_no_test_rows_and_freezes_both_models_before_evaluation(self):
        arm=self.root/'default';(arm/'features/g').mkdir(parents=True)
        (arm/'features/g/feature_contract.json').write_text('{}')
        inputs={k:'unused' for k in ('labels','splits','label_policy','registry','go_obo','example_clusters')}
        protocol=dict(supervision=inputs,disjoint_keys=['group'],settings={k:{} for k in F.KINDS})
        rows=[dict(split=s,example_id=s,group=s,genome_id=s) for s in F.T.SPLITS]
        calls=[]
        def table(*args,**kwargs):
            (arm/'table').mkdir();(arm/'table/training_table.jsonl').write_text('{}\n')
            (arm/'table/training_table_manifest.json').write_text('{}')
            return dict(table=str(arm/'table/training_table.jsonl'))
        def fit(train,cal,contract,kind,**kwargs):
            self.assertEqual([r['example_id'] for r in train],['train'])
            self.assertEqual([r['example_id'] for r in cal],['calibration'])
            calls.append(kind);return dict(kind=kind)
        with patch.object(F,'verified_export',return_value=(protocol,arm)), \
             patch.object(F.P,'build_table',side_effect=table), \
             patch.object(F.T,'load_contract',return_value={}),patch.object(F.T,'read_table',return_value=rows), \
             patch.object(F.T,'split_identity',return_value={}),patch.object(F.T,'train',side_effect=fit), \
             patch.object(F.T,'write_frozen',side_effect=lambda a,p:F.save(p,a)), \
             patch.object(F.T,'evaluate_frozen',side_effect=AssertionError('test evaluation during training')):
            F.train(self.root)
            self.assertEqual(calls,list(F.KINDS))
            self.assertEqual(len(list((arm/'models').glob('*.json'))),2)
            self.assertEqual(json.loads((arm/'training.json').read_text())['rows_by_split'],{s:1 for s in F.T.SPLITS})
            with self.assertRaises(FileExistsError):F.train(self.root)


if __name__=='__main__':unittest.main()
