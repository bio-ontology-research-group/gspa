import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from run_context_modes import compare_rows


class ContextModeComparisonTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(); self.addCleanup(self.tmp.cleanup)
        self.a = Path(self.tmp.name)/'before.jsonl'; self.b = Path(self.tmp.name)/'after.jsonl'
        self.row = dict(candidate_key='p|EC|1|CATALYTIC', genome_id='g', sequence_sources=['s'],
                        features={'seq.x':.3,'cohort.w.ec_direct':1.,'cohort.w.local_anchor_count':.5})

    def pair(self, changes=None, features=None):
        self.a.write_text(json.dumps(self.row)+'\n')
        row = dict(self.row, cohort_context_mode='PATHWAY')
        if changes: row.update(changes)
        row['features'] = dict(self.row['features'], **(features or {}))
        self.b.write_text(json.dumps(row)+'\n')

    def test_context_changes_are_reported_without_candidate_or_sequence_changes(self):
        self.pair(features={'cohort.w.local_anchor_count':0.})
        self.assertEqual(compare_rows(self.a,self.b,'PATHWAY'),(1,{'cohort.w.local_anchor_count':1}))

    def test_missing_mode_and_candidate_or_sequence_changes_fail(self):
        for changes,features in [({'cohort_context_mode':'UNIVERSAL'},None),
                                 ({'candidate_key':'different'},None),
                                 ({'sequence_sources':['different']},None),
                                 (None,{'seq.x':.4})]:
            with self.subTest(changes=changes,features=features), self.assertRaises(ValueError):
                self.pair(changes,features); compare_rows(self.a,self.b,'PATHWAY')

    def test_default_control_requires_exact_numerical_identity(self):
        self.pair({'cohort_context_mode':'UNIVERSAL'})
        self.assertEqual(compare_rows(self.a,self.b,'UNIVERSAL',exact=True),(1,{}))
        self.pair({'cohort_context_mode':'UNIVERSAL'},{'cohort.w.local_anchor_count':0.})
        with self.assertRaises(ValueError): compare_rows(self.a,self.b,'UNIVERSAL',exact=True)
