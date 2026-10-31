from datetime import date
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import prepare_reviewed_thrs_locus as R


class ReviewedSequenceTests(unittest.TestCase):
    def record(self, root, sequence='MKL', taxid='273057', version_date='01-OCT-2001', xref='AAK42622.1'):
        p = root / 'entry.txt'
        p.write_text(f'AC   Q97VW8;\nDT   {version_date}, sequence version 1.\n'
                     f'OX   NCBI_TaxID={taxid};\nDR   EMBL; AE006641; {xref}; -; Genomic_DNA.\n'
                     f'SQ   SEQUENCE\n     {sequence}\n//\n')
        return p

    def test_exact_identity_and_version_are_required(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            p = self.record(d)
            r = R.verify_entry(p, 'Q97VW8', 'AAK42622.1', {'AAK42622.1': 'MKL'}, date(2026, 6, 17))
            self.assertEqual(r['length'], 3)
            for kw in [{'sequence': 'MKI'}, {'taxid': '1'}, {'version_date': '18-JUN-2026'}, {'xref': 'AAK42622.2'}]:
                with self.subTest(kw=kw), self.assertRaises(ValueError):
                    R.verify_entry(self.record(d, **kw), 'Q97VW8', 'AAK42622.1', {'AAK42622.1': 'MKL'}, date(2026, 6, 17))

    def test_multiple_or_truncated_records_are_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p = self.record(Path(d))
            original = p.read_text()
            for text in (original * 2, original.replace('//\n', '')):
                p.write_text(text)
                with self.assertRaisesRegex(ValueError, 'one UniProt record'):
                    R.verify_entry(p, 'Q97VW8', 'AAK42622.1', {'AAK42622.1': 'MKL'}, date(2026, 6, 17))


if __name__ == '__main__':
    unittest.main()
