import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import check_family_excluded_models as C


class FamilyExcludedModelAuditTest(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(dir=os.environ.get('TMPDIR',str(Path.home()/'tmp')))
        self.addCleanup(self.temp.cleanup);self.root=Path(self.temp.name)
        self.model=self.root/'model.json';self.model.write_text(json.dumps(dict(contract=dict(features=['x']))))
        self.old=self.root/'before.tsv';self.new=self.root/'after.tsv'
        self.header='protein_id\ttype\tfunction_id\tgo_aspect\tposterior_prob\n'
        self.old.write_text(self.header+'p\tEC\tEC:1.1.1.1\t\t0.5\n')
        self.new.write_text(self.header+'p\tEC\tEC:1.1.1.1\t\t0.7\n')
        self.features=self.root/'features.jsonl';self.old_features=self.root/'old.jsonl'
        self.row=dict(candidate_key='p|EC|EC:1.1.1.1|CATALYTIC',genome_id='g',protein_id='p',
            function_type='EC',function_id='EC:1.1.1.1',features={'x':.4},sequence_sources=[],
            sequence_contracts=[],model={'probability':.7})
        self.old_features.write_text(json.dumps(self.row)+'\n')
        self.record=dict(candidate_key=self.row['candidate_key'],unsupported_reason=None,sequence_contracts=[],
                         probability=.7,score_status='scored')
        self.reason=None
        self.outcome=dict(function_key='p|EC|EC:1.1.1.1',action='rescored',learned_probability=.7,reason=None)

    def audit(self):
        self.features.write_text(json.dumps(self.row)+'\n')
        Path(str(self.new)+'.learned.json').write_text(json.dumps(dict(model={'sha256':C.F.T.sha256_file(self.model)},
                                                                      candidates=[self.record],outcomes=[self.outcome])))
        with patch.object(C.F.T,'unsupported_input_reason',return_value=self.reason), \
             patch.object(C.F.T,'predict_proba',return_value=.7):
            return C.audit_genome(self.old,self.new,self.old_features,self.features,self.model)

    def test_complete_valid_scoring_passes(self):
        self.assertEqual(self.audit()['reasons'],{'scored':1})

    def test_nonfinite_probabilities_cannot_pass_tolerance_comparisons(self):
        for target in (self.record,self.row['model'],self.outcome):
            field='learned_probability' if target is self.outcome else 'probability'
            previous=target[field];target[field]=float('nan')
            with self.subTest(field=field),self.assertRaisesRegex(ValueError,'probability'):
                self.audit()
            target[field]=previous
        self.new.write_text(self.header+'p\tEC\tEC:1.1.1.1\t\tNaN\n')
        with self.assertRaisesRegex(ValueError,'probability'):self.audit()

    def test_supported_record_cannot_claim_unsupported_status(self):
        self.record['score_status']='unsupported'
        with self.assertRaisesRegex(ValueError,'status'):self.audit()

    def test_changed_features_or_report_provenance_fail(self):
        self.row['features']['x']=.6
        with self.assertRaisesRegex(ValueError,'changed input features'):self.audit()
        self.row['features']['x']=.4;self.record['sequence_contracts']=['changed']
        with self.assertRaisesRegex(ValueError,'provenance differs'):self.audit()

    def test_final_score_or_missing_prediction_cannot_hide_behind_correct_report(self):
        self.new.write_text(self.header+'p\tEC\tEC:1.1.1.1\t\t0.8\n')
        with self.assertRaisesRegex(ValueError,'final prediction'):self.audit()
        self.new.write_text(self.header)
        with self.assertRaisesRegex(ValueError,'identity coverage'):self.audit()

    def test_taxon_withheld_score_requires_explicit_disposition_and_unchanged_output(self):
        self.outcome.update(action='scored_only',reason='forbidden by the structural taxon constraints for this organism: left as the enforcement stage decided')
        self.new.write_bytes(self.old.read_bytes())
        self.assertEqual(self.audit()['actions'],{'scored_only':1})
        self.new.write_text(self.header+'p\tEC\tEC:1.1.1.1\t\t0.7\n')
        with self.assertRaisesRegex(ValueError,'taxon-constrained prediction changed'):self.audit()
        self.outcome['reason']='unknown'
        with self.assertRaisesRegex(ValueError,'unexpected unprojected'):self.audit()

    def test_unsupported_requires_null_scores_and_unchanged_fallback(self):
        self.reason='unsupported';self.record.update(unsupported_reason=self.reason,score_status='unsupported',
            **{k:None for k in ('probability','margin','above_display_threshold','locus_probability','locus_null_probability')})
        self.row['model']['probability']=None
        with self.assertRaisesRegex(ValueError,'fallback changed'):self.audit()
        self.new.write_bytes(self.old.read_bytes())
        self.assertEqual(self.audit()['reasons'],{'unsupported':1})
        self.record['margin']=0
        with self.assertRaisesRegex(ValueError,'learned probability'):self.audit()


if __name__=='__main__':unittest.main()
