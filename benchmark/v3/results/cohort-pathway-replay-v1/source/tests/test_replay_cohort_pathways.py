from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import replay_cohort_pathways as R


def reaction(**values):
    row=dict(pathway='P',rxn='r',ec='1.1.1.1',status='bad_blast',spont='FALSE',
             is_complex='FALSE',stitle='p1 protein description',pident='42.0')
    row.update(values);return row


class CohortPathwayReplayOracle(unittest.TestCase):
    def test_weak_hits_partial_ecs_and_complex_roles_are_retained(self):
        weak=reaction(ec='EC:1.1.1.1/2.4.1')
        complex_hit=reaction(rxn='r2',is_complex='TRUE')
        local=reaction(rxn='r3',ec='1.5.1.M4')
        expected,roles=R.expected_candidates([weak,complex_hit,local,reaction(status='no_blast')])
        self.assertEqual(set(expected),{('p1','EC','EC:1.1.1.1'),('p1','EC','EC:2.4.1'),
                                       ('p1','CUSTOM','reaction-subunit:P:r2'),('p1','CUSTOM','reaction:P:r3')})
        self.assertEqual(roles[('p1','CUSTOM','reaction-subunit:P:r2')],'COMPLEX_SUBUNIT')
        self.assertEqual(expected[('p1','EC','EC:2.4.1')],[weak])

    def test_all_rows_and_pathway_identity_determine_gap_eligibility(self):
        rows=[reaction(),reaction(status='good_blast'),reaction(pathway='Q'),
              reaction(rxn='complex'),reaction(rxn='complex',is_complex='TRUE'),
              reaction(rxn='unknown',status='no_seq_data'),reaction(rxn='spont',spont='TRUE'),
              reaction(rxn='absent',status='no_blast'),reaction(pathway='not-predicted')]
        pathways=[{'ID':'P','Prediction':'TRUE'},{'ID':'Q','Prediction':'TRUE'}]
        self.assertEqual(R.eligible_groups(rows,pathways),{('Q','r'):True,('P','absent'):False})
        self.assertEqual(R.eligible_groups(list(reversed(rows)),pathways),R.eligible_groups(rows,pathways))

    def test_gap_oracle_rejects_missing_invented_and_wrong_guessed_groups(self):
        expected={('P','r'):True}
        good={'pathway_id':'gapsmith:P','reaction_id':'r','gapseq_guessed':True}
        R.verify_gap_groups([good,dict(good,ec_number='2.4.1')],expected)
        for records in [[],[dict(good,reaction_id='other')],[dict(good,gapseq_guessed=False)],
                        [dict(good,pathway_id='MetaCyc:P')]]:
            with self.subTest(records=records),self.assertRaises(ValueError):
                R.verify_gap_groups(records,expected)

    def test_malformed_tables_are_not_silently_repaired(self):
        with tempfile.TemporaryDirectory(prefix='gspa-pathway-table-') as directory:
            p=Path(directory)/'table.tsv'
            for text in ['a\ta\n1\t2\n','a\tb\n1\n','a\tb\n1\t2\t3\n','']:
                p.write_text(text)
                with self.subTest(text=text),self.assertRaises(ValueError): R.rows(p)

    def test_partial_cohort_is_rejected_before_output_without_explicit_opt_in(self):
        with tempfile.TemporaryDirectory(prefix='gspa-pathway-replay-') as directory:
            p=Path(directory)
            with patch.object(R,'audit',return_value={'status':'incomplete'}):
                with self.assertRaisesRegex(ValueError,'complete cohort acquisition required'):
                    R.run([],p/'engine',p/'training',p/'out')
            self.assertFalse((p/'out').exists())


if __name__=='__main__': unittest.main()
