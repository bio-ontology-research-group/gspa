import tempfile
import unittest
import json
from pathlib import Path

from benchmark.v3.prepare_price_layout import read_coordinates, fasta_sequences, prepare


class CoordinatesTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.path = Path(self.temp.name) / 'genes.gff'
        self.header = '#!genome-build-accession NCBI_Assembly:GCA_1.2\n##sequence-region chr 1 10000\n'

    def row(self, pid, start, end, strand='+', extras=''):
        return f'chr\tNCBI\tCDS\t{start}\t{end}\t.\t{strand}\t0\tprotein_id={pid};{extras}\n'

    def test_single_cds_and_strand_without_function_leakage(self):
        self.path.write_text(self.header + self.row('enzyme', 100, 399, '-', 'product=hexokinase;EC_number=2.7.1.1'))
        first = read_coordinates(self.path, 'GCA_1.2', {'enzyme'})
        self.path.write_text(self.header + self.row('enzyme', 100, 399, '-', 'product=WRONG;EC_number=9.9.9.9'))
        self.assertEqual(first, read_coordinates(self.path, 'GCA_1.2', {'enzyme'}))
        self.assertEqual(first[0], [('enzyme', 'chr', 100, 399, '-')])

    def test_identical_protein_at_two_loci_is_not_one_invented_locus(self):
        self.path.write_text(self.header + self.row('isozyme', 1, 300) + self.row('isozyme', 1000, 1299))
        rows, unavailable, _ = read_coordinates(self.path, 'GCA_1.2', {'isozyme', 'unmapped'})
        self.assertEqual(rows, [])
        self.assertEqual(unavailable, {'isozyme': 'multiple_or_split_CDS_loci', 'unmapped': 'missing_CDS'})

    def test_split_origin_cds_not_wrapped_to_false_long_gene(self):
        self.path.write_text(self.header + self.row('origin_gene', 9800, 10000) + self.row('origin_gene', 1, 100))
        rows, unavailable, _ = read_coordinates(self.path, 'GCA_1.2', {'origin_gene'})
        self.assertFalse(rows)
        self.assertEqual(unavailable['origin_gene'], 'multiple_or_split_CDS_loci')

    def test_wrong_version_refused(self):
        self.path.write_text(self.header)
        with self.assertRaisesRegex(ValueError, 'assembly'):
            read_coordinates(self.path, 'GCA_1.3', set())

    def test_unknown_strand_and_out_of_bounds_are_unavailable(self):
        self.path.write_text(self.header + self.row('a', 1, 300, '?') + self.row('b', 9999, 10010))
        rows, unavailable, _ = read_coordinates(self.path, 'GCA_1.2', {'a', 'b'})
        self.assertFalse(rows)
        self.assertEqual(unavailable, {'a': 'unknown_strand', 'b': 'coordinates_outside_declared_contig'})

    def test_fasta_preserves_versions_and_normalizes_terminal_stop(self):
        self.path.write_text('>WP_1.2 product name\nMAAA*\n>WP_1.3\nMATA\n')
        self.assertEqual(fasta_sequences(self.path), {'WP_1.2': 'MAAA', 'WP_1.3': 'MATA'})
        self.path.write_text('>a\nMA\n>a\nMV\n')
        with self.assertRaisesRegex(ValueError, 'duplicate'):
            fasta_sequences(self.path)

    def test_preparation_includes_genbank_and_refseq_not_zip_files(self):
        root = Path(self.temp.name) / 'genomes'
        root.mkdir()
        for asm in ['GCA_1.2', 'GCF_2.1']:
            folder = root / asm
            folder.mkdir()
            (folder / 'genomic.gff').write_text(self.header.replace('GCA_1.2', asm) + self.row('p.1', 1, 300))
            (folder / 'protein.faa').write_text('>p.1 secret_function\nMAAAA\n')
            (folder / 'assembly_data_report.json').write_text(json.dumps({'accession': asm, 'assemblyInfo': {'assemblyLevel': 'Complete Genome'}}))
            (root / (asm + '.zip')).write_text('archive ignored')
        out = Path(self.temp.name) / 'out'
        report = prepare(root, out)
        self.assertEqual(set(report['assemblies']), {'GCA_1.2', 'GCF_2.1'})
        self.assertNotIn('secret_function', (out / 'all_proteins.faa').read_text())
        self.assertIn('>GCF_2.1__p.1', (out / 'all_proteins.faa').read_text())
        self.assertNotIn('familySearchComplete', (out / 'cohort_metadata.json').read_text())
