import csv
import json
import os
from pathlib import Path
import tempfile
import unittest
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import observe_clean_cohort as O


class CleanObservationTest(unittest.TestCase):
    def setUp(self):
        base=Path(os.environ.get('TMPDIR',str(Path.home()/'tmp')));base.mkdir(exist_ok=True)
        self.tmp=tempfile.TemporaryDirectory(dir=base);self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)

    def fixture(self):
        fasta=self.root/'input.faa';fasta.write_text('>p\nMPEPTIDE\n>q\nMPROTEIN\n')
        g=dict(genome_id='genome',proteome=str(fasta),sha256=O.sha(fasta),proteins=2)
        native=self.root/'genome.clean.native.csv';native.write_text('gspa_00000000000000000000000000000000_p00000000,EC:1.1.1.1/0.0000,EC:2.2.2.2/0.3456\ngspa_00000000000000000000000000000000_p00000001,EC:3.3.3.n1/1.0000\n')
        O.C.rewrite(native,self.root/'genome.clean.tsv',0,score_kind='gmm_confidence',
            protein_ids={'gspa_00000000000000000000000000000000_p00000000':'p','gspa_00000000000000000000000000000000_p00000001':'q'},provenance=dict(driver_sha256='a'*64,
                driver_score_kind='gmm_confidence',input_fasta_sha256=g['sha256']))
        return g

    def test_every_native_call_and_metadata_survives(self):
        g=self.fixture();contract={'schema':'gspa.predictor-contract/1','tool':'fixture'}
        out=self.root/'claims.jsonl';result=O.audit_native(self.root,g,contract,out,driver_sha256="a"*64)
        rows=[json.loads(l) for l in out.read_text().splitlines()]
        self.assertEqual((result['proteins'],result['calls'],result['zero_confidence_calls']),(2,3,1))
        self.assertEqual([r['raw_score'] for r in rows],[0,.3456,1])
        self.assertEqual({r['protein_id'] for r in rows},{'p','q'})
        self.assertTrue(all(r['metadata']['predictor_contract']==contract for r in rows))
        self.assertTrue(all(r['metadata']['native_adapter_contract'].startswith('clean:maxsep/') for r in rows))
        self.assertTrue(all('calibrated_prob' not in r for r in rows))

    def test_native_typed_disagreement_even_with_updated_hash(self):
        g=self.fixture();tsv=self.root/'genome.clean.tsv'
        with tsv.open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
        rows[0]['score']='0.9'
        with tsv.open('w') as f:
            w=csv.DictWriter(f,fieldnames=O.C.HEADER,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
        p=self.root/'genome.clean.receipt.json';r=json.loads(p.read_text());r['output_sha256']=O.sha(tsv);p.write_text(json.dumps(r))
        with self.assertRaisesRegex(ValueError,'score differs'):
            O.audit_native(self.root,g,{},self.root/'claims.jsonl',driver_sha256='a'*64)
        self.assertFalse((self.root/'claims.jsonl').exists())

    def test_missing_input_protein_not_success_zero(self):
        g=self.fixture();g['proteins']=3
        with self.assertRaisesRegex(ValueError,'coverage'):
            O.audit_native(self.root,g,{},self.root/'claims.jsonl',driver_sha256='a'*64)

    def test_driver_is_checked_against_plan(self):
        g=self.fixture()
        with self.assertRaisesRegex(ValueError,'contract mismatch'):
            O.audit_native(self.root,g,{},self.root/'claims.jsonl',driver_sha256='b'*64)

    def test_swapped_native_identity_map_rejected(self):
        g=self.fixture();p=self.root/'genome.clean.receipt.json';r=json.loads(p.read_text())
        keys=sorted(r['native_to_input_ids']);r['native_to_input_ids'][keys[0]]='q';r['native_to_input_ids'][keys[1]]='p'
        p.write_text(json.dumps(r))
        with self.assertRaisesRegex(ValueError,'deterministic input order'):
            O.audit_native(self.root,g,{},self.root/'claims.jsonl',driver_sha256='a'*64)

    def test_previous_claims_cannot_be_overwritten(self):
        g=self.fixture();p=self.root/'claims.jsonl';p.write_text('previous result\n')
        with self.assertRaises(FileExistsError):
            O.audit_native(self.root,g,{},p,driver_sha256='a'*64)
        self.assertEqual(p.read_text(),'previous result\n')

    def test_changed_reference_input_rejected(self):
        g=self.fixture();p=Path(g['proteome']);p.write_text(p.read_text()+'M')
        with self.assertRaisesRegex(ValueError,'frozen input changed'):O.check({str(p):g['sha256']})

    def test_overlap_separates_full_and_truncated_without_ec_labels(self):
        f=self.root/'input.faa';f.write_text('>exact\nMPEPTIDE\n>tail\n'+'A'*1022+'T\n>absent\nQQQQ\n')
        ref=self.root/'reference.tsv';ref.write_text('Entry\tSequence\nreference1\tMPEPTIDE\nreference2\t'+'A'*1022+'G\n')
        result=O.reference_overlap(ref,[dict(genome_id='g',proteome=str(f))])
        self.assertEqual((result['exact_proteins'],result['truncated_input_matches']),(1,2))
        self.assertEqual({r['protein_id'] for r in result['rows']},{'exact','tail'})
        self.assertEqual(result['rows'][1]['exact_reference_entries'],[])


if __name__=='__main__':unittest.main()
