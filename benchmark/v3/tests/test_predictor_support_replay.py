"""Migration invariants; these fixtures do not represent predictor accuracy."""
import copy
import json
from pathlib import Path
import tempfile
import unittest
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import replay_predictor_support as R


class PredictorSupportReplayTests(unittest.TestCase):
    def setUp(self):
        self.directory=tempfile.TemporaryDirectory();self.addCleanup(self.directory.cleanup)
        self.old=Path(self.directory.name)/'old.jsonl';self.new=Path(self.directory.name)/'new.jsonl'
        c=dict(schema='gspa.predictor-contract/1',tool='fixture',tool_version='1',score_model='fixture/1',
               **{k:'a'*64 for k in ('tool_sha256','reference_sha256','parameters_sha256','adapter_sha256')})
        self.identities=[dict(source='fixture',contract=c)]
        self.rows=[dict(candidate_key=str(i),sequence_sources=['fixture'],features={'seq.max':.5},
                        label=i,split='train') for i in range(2)]
        self.current=[dict(copy.deepcopy(r),sequence_contracts=copy.deepcopy(self.identities)) for r in self.rows]

    def compare(self):
        for p,rows in [(self.old,self.rows),(self.new,self.current)]:
            p.write_text(''.join(json.dumps(r)+'\n' for r in rows))
        return R.compare_rows(self.old,self.new,self.identities)

    def test_only_verified_sideband_changes(self):
        self.assertEqual(self.compare(),2)

    def test_historical_implicit_universal_is_the_only_mode_normalization(self):
        for r in self.current:r['cohort_context_mode']='UNIVERSAL'
        self.assertEqual(self.compare(),2)
        self.current[0]['cohort_context_mode']='PATHWAY'
        with self.assertRaisesRegex(ValueError,'changed'):self.compare()

    def test_numeric_labels_splits_and_candidates_cannot_change(self):
        for key,value in [('features',{'seq.max':.6}),('label',1),('split','test'),('candidate_key','other')]:
            with self.subTest(key=key):
                original=self.current[0][key];self.current[0][key]=value
                with self.assertRaisesRegex(ValueError,'changed'):self.compare()
                self.current[0][key]=original

    def test_missing_extra_or_reordered_rows_fail(self):
        original=self.current
        for rows in [original[:1],original+original[:1],list(reversed(original))]:
            self.current=rows
            with self.assertRaisesRegex(ValueError,'changed'):self.compare()

    def test_missing_identity_fails(self):
        self.current[0].pop('sequence_contracts')
        with self.assertRaisesRegex(ValueError,'missing predictor'):self.compare()

    def test_plausible_but_unverified_identity_fails(self):
        self.current[0]['sequence_contracts'][0]['contract']['tool_version']='2'
        with self.assertRaisesRegex(ValueError,'differs from verified'):self.compare()


if __name__=='__main__':unittest.main()
