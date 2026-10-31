import hashlib
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import verify_reported_reference_overlap as V


class ReportedReferenceOverlap(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory(prefix='gspa-reference-overlap-');self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)
        self.a=self.root/'a.faa';self.a.write_text('>p1 a protein\nMa\nAA\n')
        self.b=self.root/'b.faa';self.b.write_text('>p1 another genome\nMAAT\n')
        self.ref=self.root/'ref.fasta';self.ref.write_text('>ref1 descriptive header\nMAAA\n')
        self.plan={'proteomes':{g:{'path':str(p),'sha256':V.digest(p)} for g,p in [('A',self.a),('B',self.b)]},
                   'reference_hashes':{str(self.ref):V.digest(self.ref)},
                   'pairs':[dict(genome_id=g,protein_id='p1',reference_file=str(self.ref),reference_id='ref1',
                                 query_sequence_sha256=hashlib.sha256(seq.encode()).hexdigest())
                            for g,seq in [('A','MAAA'),('B','MAAT')]]}

    def test_full_sequence_identity_is_genome_scoped(self):
        result=V.verify(self.plan)
        self.assertEqual(result['counts']['A']['exact_proteins'],1)
        self.assertEqual(result['counts']['B']['exact_proteins'],0)
        self.assertEqual([r['exact_sequence'] for r in result['results']],[True,False])

    def test_changed_reference_is_rejected(self):
        self.ref.write_text('>ref1\nMAAT\n')
        with self.assertRaisesRegex(ValueError,'reference file changed'): V.verify(self.plan)

    def test_changed_reported_hit_source_is_rejected(self):
        source=self.root/'reported.tsv';source.write_text('original hit\n')
        self.plan['source_inputs']={str(source):V.digest(source)}
        source.write_text('different hit\n')
        with self.assertRaisesRegex(ValueError,'reported-hit source changed'): V.verify(self.plan)

    def test_conflicting_reference_ids_are_not_silently_overwritten(self):
        self.ref.write_text('>ref1\nMAAA\n>ref1\nMAAT\n')
        self.plan['reference_hashes'][str(self.ref)]=V.digest(self.ref)
        with self.assertRaisesRegex(ValueError,'conflicting sequence ID'): V.verify(self.plan)


if __name__=='__main__': unittest.main()
