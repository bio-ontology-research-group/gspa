import tempfile
import unittest
from pathlib import Path
from benchmark.v3.recover_price_coordinates import exact_hits


class ExactCoordinatesTest(unittest.TestCase):
    def run_rows(self, rows):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / 'matches.tsv'
            p.write_text('\n'.join(rows) + '\n')
            return exact_hits(p, {'a.1': 'MAAA', 'b.1': 'MSSS'})

    def test_reverse_strand_and_duplicate_hsps(self):
        row = 'a.1\tchr\t100\t4\t4\t1\t4\t120\t109\t4\t0\t-1'
        x = self.run_rows([row, row])
        self.assertEqual(x['a.1'], [{'contig': 'chr', 'start': 109, 'end': 120, 'strand': '-'}])
        self.assertEqual(x['b.1'], [])

    def test_partial_and_gapped_hits_do_not_supply_a_coordinate(self):
        rows = ['a.1\tchr\t100\t3\t4\t2\t4\t10\t18\t3\t0\t1',
                'a.1\tchr\t100\t4\t4\t1\t4\t10\t24\t4\t1\t1']
        self.assertEqual(self.run_rows(rows)['a.1'], [])

    def test_multiple_exact_loci_remain_multiple(self):
        rows = [f'a.1\tchr\t100\t4\t4\t1\t4\t{s}\t{s+11}\t4\t0\t1' for s in (10, 100)]
        self.assertEqual(len(self.run_rows(rows)['a.1']), 2)
