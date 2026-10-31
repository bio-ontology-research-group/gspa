"""Private truth must never manufacture negatives or select public cases."""
import copy
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import prepare_agent_truth as P
from test_evaluate_stage3_agent import bundle


class PrivateTruthTests(unittest.TestCase):
    def test_only_experimental_positive_ancestors_are_joined(self):
        with tempfile.TemporaryDirectory() as directory:
            d = Path(directory)
            (d / 'go.obo').write_text('format-version: 1.2\n\n[Term]\nid: GO:0000001\n\n[Term]\nid: GO:0000003\nis_a: GO:0000001\n')
            (d / 'labels.tsv').write_text('genome_id\tprotein_id\tfunction_id\tevidence\ttier\tswissprot_accession\tmatch_basis\n'
                'G1\tp1\tGO:0000003\tIDA\texperimental\tP123\tidentity\n'
                'G2\tp2\tGO:0000001\tIEA\tcomputational\tP456\tidentity\n')
            public = bundle()
            before = copy.deepcopy(public)
            labels, provenance = P.prepare(public, d / 'labels.tsv', d / 'go.obo')
            self.assertEqual(public, before)
            self.assertEqual(len(labels['calibration']), 1)
            self.assertEqual(labels['calibration'][0]['status'], 'positive')
            self.assertEqual(labels['calibration'][0]['function_id'], 'GO:0000001')
            self.assertTrue(labels['calibration'][0]['source']['propagated'])
            self.assertEqual(labels['test'], [])
            self.assertEqual(provenance['counts']['input_rows'], 2)


if __name__ == '__main__':
    unittest.main()
