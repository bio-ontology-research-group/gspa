import csv
import math
import os
from pathlib import Path
import sys
import tempfile
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import compare_gap_predictions as C


class GapPredictionComparisonTest(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(dir=os.environ.get('TMPDIR',str(Path.home()/'tmp')))
        self.addCleanup(self.temp.cleanup);self.root=Path(self.temp.name)
        self.a=self.root/'a.tsv';self.b=self.root/'b.tsv';self.d=self.root/'changes.tsv'

    def write(self,path,scores,changes=None):
        rows=[]
        for i,p in enumerate(scores):
            row=dict(protein_id='p'+str(i),type='EC',function_id='1.1.1.1',go_aspect='',
                likelihood_logodds='0.1',n_supporting='1',posterior_prob=p,
                final_logodds=math.log(p/(1-p)),priors_fired='gap_filling' if path==self.b else '')
            row.update(changes or {});rows.append(row)
        with path.open('w') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)

    def test_both_directions_and_all_threshold_crossings_are_retained(self):
        self.write(self.a,[.2,.8,.4]);self.write(self.b,[.7,.3,.4])
        result=C.compare(self.a,self.b,self.d)
        self.assertEqual(result['counts'],dict(predictions=3,changed=2,increased=1,decreased=1,
            **{'crossed_above_0.5':1,'crossed_below_0.5':1}))
        self.assertAlmostEqual(result['max_absolute_change'],.5)
        with self.d.open() as f:self.assertEqual(len(list(csv.DictReader(f,delimiter='\t'))),2)

    def test_changed_sequence_evidence_or_identity_is_not_a_gap_effect(self):
        self.write(self.a,[.2])
        for field,value in [('likelihood_logodds','2'),('n_supporting','2'),('protein_id','other')]:
            self.write(self.b,[.7],{field:value})
            with self.assertRaises(ValueError):C.compare(self.a,self.b,self.d)

    def test_nonfinite_or_inconsistent_scores_fail(self):
        self.write(self.a,[.2])
        for change in [{'posterior_prob':'nan'},{'final_logodds':'inf'},
                       {'posterior_prob':2},{'final_logodds':0}]:
            self.write(self.b,[.7],change)
            with self.assertRaises(ValueError):C.compare(self.a,self.b,self.d)


if __name__=='__main__':unittest.main()
