import gzip
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import audit_gaf_negatives as A


def row(relation='NOT|enables', evidence='IDA', **kw):
    fields = ['UniProtKB', 'P12345', 'gene', relation, 'GO:0000001', 'PMID:12345', evidence,
              '', 'F', 'name', '', 'protein', 'taxon:1', '20260615', 'UniProt', '', '']
    for key, index in {'extension': 15, 'form': 16, 'taxon': 12, 'db': 0}.items():
        if key in kw:
            fields[index] = kw[key]
    return ('\t'.join(fields) + '\n').encode()


class GafNegativeTests(unittest.TestCase):
    def test_polarity_relation_references_and_raw_identity_are_preserved(self):
        r = A.negative_record(row('NOT|contributes_to'), 42)
        self.assertEqual(r['relation'], 'contributes_to')
        self.assertEqual(r['references'], ['PMID:12345'])
        self.assertEqual(r['line'], 42)
        self.assertEqual(r['evidence_category'], 'experimental')
        self.assertIsNone(A.negative_record(row('enables'), 1))
        self.assertEqual(A.negative_record(row(evidence='IKR'), 1)['evidence_category'], 'other')

    def test_missing_required_fields_and_ambiguous_qualifiers_fail(self):
        for raw in [row('NOT'), row('NOT|enables|contributes_to'), row().replace(b'PMID:12345', b''),
                    row().rstrip(b'\n') + b'\textra\n', row(taxon='unknown')]:
            with self.subTest(raw=raw), self.assertRaises(ValueError):
                A.negative_record(raw, 1)

    def test_context_forms_and_taxon_mismatch_are_never_silently_transferred(self):
        matches = {'P12345': [{'genome_id': 'g1', 'protein_id': 'p1', 'status': 'matched', 'taxid': '1', 'basis': 'accession_and_sequence'}]}
        clean = A.candidate_joins(A.negative_record(row(), 1), matches)[0]
        self.assertEqual(clean['status'], 'needs_snapshot_and_reference_review')
        r = A.negative_record(row(extension='occurs_in(CL:0000001)', form='UniProtKB:P12345-2', taxon='taxon:2|taxon:3'), 1)
        joined = A.candidate_joins(r, matches)[0]
        self.assertEqual(set(joined['exclusion_reasons']), {'annotation_extension_requires_context', 'specific_gene_product_form', 'taxon_differs', 'interaction_context'})
        self.assertEqual(A.candidate_joins(A.negative_record(row(db='ComplexPortal'), 1), matches), [])

    def test_complete_scan_validates_crc_counts_hashes_and_refuses_overwrite(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            source = d / 'source.gaf.gz'
            source.write_bytes(gzip.compress(b'!gaf-version: 2.2\n' + row('enables') + row() + row(evidence='IEA')))
            result = A.audit(source, d / 'out')
            self.assertEqual(result['negative_records'], 2)
            self.assertEqual(result['experimental_negative_records'], 1)
            self.assertEqual(result['source']['sha256'], A.sha(source))
            records = [json.loads(s) for s in (d/'out/negatives.jsonl').read_text().splitlines()]
            self.assertEqual([r['line'] for r in records], [3, 4])
            with self.assertRaises(FileExistsError):
                A.audit(source, d / 'out')

    def test_corrupt_compressed_source_cannot_certify_partial_scan(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            source = d / 'bad.gz'
            source.write_bytes(gzip.compress(b'!gaf-version: 2.2\n' + row())[:-5])
            with self.assertRaisesRegex(ValueError, 'CRC'):
                A.audit(source, d / 'out')
            self.assertFalse((d/'out/inventory.json').exists())


if __name__ == '__main__':
    unittest.main()
