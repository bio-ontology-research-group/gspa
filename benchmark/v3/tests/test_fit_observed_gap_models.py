import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import fit_observed_gap_models as B


class ObservedGapFittingTest(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(dir=os.environ.get('TMPDIR',str(Path.home()/'tmp')))
        self.addCleanup(self.temp.cleanup);self.root=Path(self.temp.name)

    def fixture(self):
        features=self.root/'features/g';features.mkdir(parents=True)
        (features/'feature_contract.json').write_text('{}')
        inputs={k:'fixture' for k in ('labels','splits','label_policy','registry','go_obo','example_clusters')}
        recipes={m:dict(supervision=inputs,features=str(features.parent),original_table='original',
            disjoint_keys=['group'],settings={k:{} for k in B.F.KINDS}) for m in B.G.MODES}
        return dict(recipes=recipes),[dict(split=s,example_id=s,group=s,genome_id=s) for s in B.F.T.SPLITS]

    def test_all_four_models_freeze_without_test_rows_or_evaluation_during_fit(self):
        protocol,rows=self.fixture();calls=[]
        def table(*args,**kwargs):
            target=args[5];target.mkdir(parents=True)
            (target/'training_table.jsonl').write_text('{}\n')
            (target/'training_table_manifest.json').write_text('{}')
            return dict(table=str(target/'training_table.jsonl'))
        def fit(train,cal,contract,kind,**kwargs):
            self.assertEqual([r['split'] for r in train],['train'])
            self.assertEqual([r['split'] for r in cal],['calibration'])
            calls.append(kind);return {'kind':kind}
        with patch.object(B,'verify',return_value=protocol),patch.object(B.F.P,'build_table',side_effect=table), \
             patch.object(B.M,'compare_supervision',return_value=3),patch.object(B.F.T,'load_contract',return_value={}), \
             patch.object(B.F.T,'read_table',return_value=rows),patch.object(B.F.T,'split_identity',return_value={}), \
             patch.object(B.F.T,'train',side_effect=fit),patch.object(B.F.T,'write_frozen',side_effect=lambda a,p:B.F.save(p,a)), \
             patch.object(B.F.T,'evaluate_frozen',side_effect=AssertionError('test evaluation during fit')):
            counts=B.train(self.root)
            self.assertEqual(calls,['logistic','gbdt']*2)
            self.assertEqual(len(list(self.root.glob('*/default/models/model_*.json'))),4)
            self.assertEqual(counts['homology']['rows_by_split'],{s:1 for s in B.F.T.SPLITS})
            with self.assertRaises(FileExistsError):B.train(self.root)
        evaluated=[]
        def evaluate(model,contract,test,*args):
            self.assertEqual([r['split'] for r in test],['test'])
            self.assertEqual(len(list(self.root.glob('*/default/models/model_*.json'))),4)
            evaluated.append(model);return {'fixture':True}
        with patch.object(B,'verify',return_value=protocol),patch.object(B.F.T,'load_contract',return_value={}), \
             patch.object(B.F.T,'read_table',return_value=rows),patch.object(B.F.T,'evaluate_frozen',side_effect=evaluate):
            B.evaluate(self.root);self.assertEqual(len(evaluated),4)
            with self.assertRaises(FileExistsError):B.evaluate(self.root)

    def test_missing_frozen_model_prevents_any_evaluation(self):
        protocol,_=self.fixture();B.F.save(self.root/'training.json',dict(status='verified',outputs={}))
        with patch.object(B,'verify',return_value=protocol), \
             patch.object(B.F.T,'evaluate_frozen',side_effect=AssertionError('premature evaluation')):
            with self.assertRaisesRegex(ValueError,'missing frozen model'):B.evaluate(self.root)

    def test_changed_protocol_inputs_or_runtime_reject_later_training(self):
        source=self.root/'source';source.write_text('original')
        protocol=dict(python_version=sys.version,numpy_version=B.F.T.np.__version__,
            modes=list(B.G.MODES),model_kinds=list(B.F.KINDS))
        B.F.save(self.root/'protocol.json',protocol)
        reg=dict(inputs=B.F.hashes([source]),protocol_sha256=B.F.T.sha256_file(self.root/'protocol.json'))
        B.F.save(self.root/'registration.json',reg);self.assertEqual(B.verify(self.root),protocol)
        source.write_text('changed')
        with self.assertRaisesRegex(ValueError,'changed'):B.verify(self.root)
        source.write_text('original');protocol['python_version']='different';B.F.save(self.root/'protocol.json',protocol)
        with self.assertRaisesRegex(ValueError,'protocol changed'):B.verify(self.root)
        reg['protocol_sha256']=B.F.T.sha256_file(self.root/'protocol.json');B.F.save(self.root/'registration.json',reg)
        with self.assertRaisesRegex(ValueError,'runtime or scope'):B.verify(self.root)

    def test_unverified_or_modified_template_cannot_register_training(self):
        comparison=self.root/'comparison';comparison.mkdir()
        B.F.save(comparison/'verification.json',dict(status='verified',inputs={},outputs={}))
        B.F.save(comparison/'protocol.json',{})
        template=self.root/'template';template.mkdir()
        source=template/'model.json';source.write_text('{}')
        receipt=dict(status='failed',outputs=B.F.hashes([source]))
        B.F.save(template/'training.json',receipt)
        templates={mode:template for mode in B.G.MODES};out=self.root/'registered'
        with patch.object(B.G,'verify',return_value={}):
            with self.assertRaisesRegex(ValueError,'verified template'):B.register(comparison,templates,out)
            receipt['status']='verified';B.F.save(template/'training.json',receipt);source.write_text('changed')
            with self.assertRaisesRegex(ValueError,'changed'):B.register(comparison,templates,out)
        self.assertFalse(out.exists())


if __name__=='__main__':unittest.main()
