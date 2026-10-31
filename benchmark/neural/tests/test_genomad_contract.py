"""Native-output contract fixtures, not geNomad accuracy or live execution tests."""
import argparse
import csv
import importlib.util
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

SOURCE = Path(__file__).resolve().parents[1] / 'run_genomic_predictors.py'
spec = importlib.util.spec_from_file_location('genomic_sidecar_contract', SOURCE)
S = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = S
spec.loader.exec_module(S)

VIRUS = 'seq_name\tlength\ttopology\tcoordinates\tvirus_score\ttaxonomy\n'
PLASMID = 'seq_name\tlength\ttopology\tplasmid_score\n'


class GenomadContract(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.fasta = self.root / 'genome.fna'
        self.fasta.write_text('>ctg|original\n' + 'A' * 100 + '\n')
        self.row = S.ManifestRow('query', self.fasta, None, self.root / 'out')
        self.output = S.output_path(self.row, 'genomad')
        self.args = argparse.Namespace(db_path=self.root / 'db', genomad_sif=None,
                                       threads=2, min_score=0.0)

    def run_native(self, virus, plasmid, returncode=0):
        def execute(command, **kw):
            self.assertIn('--threads', command)
            summary = Path(command[-2]) / 'genome_summary'
            summary.mkdir()
            for kind, body in [('virus', virus), ('plasmid', plasmid)]:
                if body is not None:
                    (summary / f'genome_{kind}_summary.tsv').write_text(body)
            return subprocess.CompletedProcess(command, returncode, 'fixture stdout', 'fixture stderr')
        with patch.object(S.shutil, 'which', return_value='/controlled/genomad'), \
             patch.object(S.subprocess, 'run', side_effect=execute):
            S.run_genomad([self.row], self.args)

    def assert_rejected(self, virus, plasmid, returncode=0):
        self.output.write_text('stale previous success\n')
        with self.assertRaises((ValueError, RuntimeError, SystemExit)):
            self.run_native(virus, plasmid, returncode)
        self.assertFalse(self.output.exists(), 'failed attempt must not publish stale or partial output')

    def test_missing_native_summaries_are_not_zero_hits(self):
        for virus, plasmid in [(None, None), (VIRUS, None), (None, PLASMID)]:
            with self.subTest(virus=virus, plasmid=plasmid):
                self.assert_rejected(virus, plasmid)

    def test_malformed_native_rows_are_not_filtered_into_success(self):
        cases = {
            'empty native file': '',
            'missing required columns': 'seq_name\tlength\n',
            'truncated row': VIRUS + 'ctg|original\t100\n',
            'bad length': VIRUS + 'ctg|original\toops\tDTR\tNA\t0.9\ttaxon\n',
            'NaN': VIRUS + 'ctg|original\t100\tDTR\tNA\tNaN\ttaxon\n',
            'infinite': VIRUS + 'ctg|original\t100\tDTR\tNA\tinf\ttaxon\n',
            'negative': VIRUS + 'ctg|original\t100\tDTR\tNA\t-0.1\ttaxon\n',
            'above one': VIRUS + 'ctg|original\t100\tDTR\tNA\t1.1\ttaxon\n',
            'bad coordinates': VIRUS + 'ctg|original\t10\tProvirus\tnonsense\t0.9\ttaxon\n',
            'mismatched suffix': VIRUS + 'ctg|original|provirus_2_11\t10\tProvirus\t3-12\t0.9\ttaxon\n',
        }
        for name, virus in cases.items():
            with self.subTest(name=name): self.assert_rejected(virus, PLASMID)

    def test_nonzero_native_exit_is_failure_and_clears_stale_output(self):
        self.assert_rejected(VIRUS, PLASMID, returncode=7)

    def test_late_invalid_row_never_publishes_partial_success(self):
        self.assert_rejected(VIRUS + 'ctg|original\t100\tDTR\tNA\t0.9\ttaxon\n' +
                             'ctg|original\t100\tDTR\tNA\tbroken\ttaxon\n', PLASMID)

    def test_both_native_headers_are_a_valid_zero_result(self):
        self.run_native(VIRUS, PLASMID)
        with self.output.open() as f: self.assertEqual(list(csv.reader(f, delimiter='\t')), [S.OUTPUT_HEADER])

    def test_zero_scores_provirus_coordinates_and_plasmids_are_preserved(self):
        self.run_native(VIRUS + 'ctg|original\t100\tDTR\tNA\t0\ttaxon\n' +
                        'ctg|original|provirus_2_11\t10\tProvirus\t2-11\t0.9\ttaxon\n',
                        PLASMID + 'ctg|original\t100\tDTR\t1\n')
        with self.output.open() as f: rows = list(csv.DictReader(f, delimiter='\t'))
        self.assertEqual([(r['contig_id'], r['region_start'], r['region_end'], r['region_type'], r['score']) for r in rows],
                         [('ctg|original','1','100','viral_contig','0.0000'),
                          ('ctg|original','2','11','prophage','0.9000'),
                          ('ctg|original','1','100','plasmid','1.0000')])


if __name__ == '__main__': unittest.main(verbosity=2)
