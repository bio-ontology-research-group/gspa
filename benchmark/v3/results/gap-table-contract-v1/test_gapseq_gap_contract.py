"""Reaction-level gap eligibility through the real converter CLI (synthetic tables)."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

SCRIPT = Path(__file__).resolve().parents[2] / 'parse_gapseq_gaps.py'
HEADER = ['pathway', 'rxn', 'ec', 'status', 'spont', 'is_complex', 'complex.status']


class GapContract(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def convert(self, rows, pathways='ID\tPrediction\nP\tTRUE\n', header=HEADER, valid=True):
        p, r, m, o = [self.root / n for n in ['Pathways.tbl', 'Reactions.tbl', 'ec2go', 'gaps.jsonl']]
        p.write_text(pathways)
        r.write_text('\t'.join(header)+'\n'+''.join('\t'.join(row)+'\n' for row in rows))
        m.write_text('EC:1.1.1.1\tGO:0000001\nEC:1.1.1.1\tGO:0000002\n')
        o.write_text('previous output\n')
        result = subprocess.run([sys.executable,str(SCRIPT),'--pathways-tbl',str(p),
            '--reactions-tbl',str(r),'--ec2go',str(m),'--out',str(o),'--tag','synthetic'],
            text=True,capture_output=True,env=os.environ,timeout=10)
        if not valid:
            self.assertNotEqual(result.returncode,0,result.stdout)
            self.assertEqual(o.read_text(),'previous output\n')
            return
        self.assertEqual(result.returncode,0,result.stderr)
        return [json.loads(line) for line in o.read_text().splitlines()]

    @staticmethod
    def row(status='bad_blast', **kw):
        d=dict(zip(HEADER,['P','R','1.1.1.1',status,'FALSE','FALSE','']))
        d.update(kw)
        return [d[k] for k in HEADER]

    def test_strong_hit_suppresses_weak_hit_gap_in_either_order(self):
        rows=[self.row(),self.row('good_blast')]
        self.assertEqual(self.convert(rows),[])
        self.assertEqual(self.convert(rows[::-1]),[])

    def test_no_reference_is_unknown_not_missing_enzyme(self):
        self.assertEqual(self.convert([self.row('no_seq_data')]),[])

    def test_spontaneous_reactions_are_not_missing_enzymes(self):
        self.assertEqual(self.convert([self.row('no_blast',spont='TRUE')]),[])
        self.assertEqual(self.convert([self.row('spontaneous')]),[])

    def test_complex_missing_subunits_are_not_whole_ec_gaps(self):
        self.assertEqual(self.convert([self.row('no_blast',is_complex='TRUE')]),[])

    def test_supported_negative_search_retains_all_go_alternatives_once(self):
        rows=self.convert([self.row('no_blast'),self.row(),self.row()])
        self.assertEqual({r['go_term'] for r in rows},{'GO:0000001','GO:0000002'})
        self.assertEqual(len(rows),2)
        self.assertTrue(all(r['gapseq_guessed'] for r in rows))

    def test_unknown_evidence_suppresses_gap_even_with_weak_hit(self):
        self.assertEqual(self.convert([self.row(),self.row('no_seq_data')]),[])

    def test_inactive_pathway_does_not_emit_gap(self):
        self.assertEqual(self.convert([self.row()],pathways='ID\tPrediction\nP\tFALSE\n'),[])

    def test_bad_header_rejected_before_overwriting_output(self):
        self.convert([self.row()],header=['wrong',*HEADER[1:]],valid=False)

    def test_truncated_row_rejected_before_overwriting_output(self):
        self.convert([self.row()[:-1]],valid=False)

    def test_unknown_status_rejected_before_overwriting_output(self):
        self.convert([self.row('unexpected')],valid=False)

    def test_invalid_boolean_rejected_before_overwriting_output(self):
        self.convert([self.row(spont='maybe')],valid=False)
        self.convert([self.row()],pathways='ID\tPrediction\nP\tmaybe\n',valid=False)

    def test_ec_prefix_normalized_for_lookup(self):
        rows=self.convert([self.row('no_blast',ec='EC:1.1.1.1')])
        self.assertEqual({r['go_term'] for r in rows},{'GO:0000001','GO:0000002'})
        self.assertTrue(all(r['ec_number']=='1.1.1.1' for r in rows))

    def test_alternative_ecs_and_unmapped_go_remain_explicit(self):
        rows=self.convert([self.row('no_blast',ec='1.1.1.1 / 2.2.2.2')])
        self.assertEqual(len(rows),3)
        self.assertIn(('2.2.2.2',None),{(r['ec_number'],r['go_term']) for r in rows})
        self.assertTrue(all(not r['gapseq_guessed'] for r in rows))

    def test_missing_ec_is_not_invented(self):
        rows=self.convert([self.row('no_blast',ec='NA')])
        self.assertEqual(len(rows),1)
        self.assertIsNone(rows[0]['ec_number'])
        self.assertIsNone(rows[0]['go_term'])

    def test_header_only_is_a_valid_empty_search_table(self):
        self.assertEqual(self.convert([]),[])

    def test_conflicting_pathway_prediction_rejected(self):
        self.convert([self.row()],pathways='ID\tPrediction\nP\tTRUE\nP\tFALSE\n',valid=False)

    def test_invalid_ec_rejected_before_overwriting_output(self):
        self.convert([self.row(ec='nonsense')],valid=False)


if __name__ == '__main__':
    unittest.main()
