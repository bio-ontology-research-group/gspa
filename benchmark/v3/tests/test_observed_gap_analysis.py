from pathlib import Path
import sys
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import analyze_observed_gap_models as A


def record(i,label,scores,split='test',genome='g1'):
    return dict(example_id=str(i),genome_id=genome,label=label,split=split,
        label_basis='positive_truth' if label else 'unconfirmed',scores=dict(zip(A.ARMS,scores)))


class ObservedGapAnalysisTest(unittest.TestCase):
    def test_calibration_only_selection_and_conservative_ties(self):
        rows=[record(0,1,[.9]*4,'calibration'),record(1,0,[.8]*4,'calibration'),
              record(2,0,[.6]*4,'calibration'),record(3,1,[.5]*4,'calibration')]
        cuts=A.select_thresholds(rows)
        self.assertTrue(all(v['value']==.9 for v in cuts.values()))
        rows[0]['split']='test'
        with self.assertRaisesRegex(ValueError,'calibration rows only'):A.select_thresholds(rows)
        with self.assertRaisesRegex(ValueError,'test rows only'):A.evaluate_rows(rows,cuts,10)

    def test_every_supervised_function_is_required_and_ec_normalized(self):
        row=dict(example_id='x',protein_id='p|pipe',function_type='EC',function_id='EC:1.1.1.1',
            label=1,genome_id='g',split='test',label_basis='positive_truth')
        key=('p|pipe','EC','1.1.1.1');scores={a:{'g':{key:.5}} for a in A.ARMS}
        self.assertEqual(A.attach_scores([row],scores)[0]['scores'],{a:.5 for a in A.ARMS})
        with self.assertRaisesRegex(ValueError,'duplicate supervised'):A.attach_scores([row,row],scores)
        del scores['gbdt']['g'][key]
        with self.assertRaisesRegex(ValueError,'missing supervised'):A.attach_scores([row],scores)
        scores['gbdt']['g'][key]=float('nan')
        with self.assertRaisesRegex(ValueError,'probability'):A.attach_scores([row],scores)

    def test_corrections_harms_and_genome_weighting_are_explicit(self):
        rows=[record(0,1,[.2,.9,.9,.9],genome='small'),
              record(1,0,[.1,.9,.9,.9],genome='large'),
              record(2,0,[.1,.9,.9,.9],genome='large')]
        result=A.paired(rows,'rule_observed_gaps','rule_no_gaps',{a:.5 for a in A.ARMS},100)
        self.assertEqual((result['reference_corrected'],result['reference_harmed'],result['changed_calls']),(1,2,3))
        # Equal genome weights: (-.63 + .80) / 2, not the pooled-row mean.
        self.assertAlmostEqual(result['mean_genome_brier_delta']['estimate'],.085)
        self.assertEqual(result['mean_genome_brier_delta']['n_genomes'],2)
        single=A.paired(rows[:1],'logistic','rule_no_gaps',{a:.5 for a in A.ARMS},100)
        self.assertIsNone(single['mean_genome_brier_delta']['ci_low'])

    def test_all_test_rows_and_unconfirmed_basis_survive_evaluation(self):
        rows=[record(0,1,[.8,.8,.8,.8]),record(1,0,[.2,.2,.2,.2])]
        cuts={a:{'value':.5} for a in A.ARMS}
        result=A.evaluate_rows(rows,cuts,10)
        self.assertEqual(result['label_basis'],{'positive_truth':1,'unconfirmed':1})
        for arm in result['policies']['fixed_0.5']['arms'].values():
            self.assertEqual((arm['n'],arm['tp'],arm['tn']),(2,1,1))
            self.assertAlmostEqual(arm['brier'],.04)
            self.assertAlmostEqual(arm['ece_10bin'],.2)
        self.assertEqual(len(result['policies']['fixed_0.5']['comparisons']),6)


if __name__=='__main__':unittest.main()
