import copy
import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import replay_predictor_contract as R


class PredictorContractReplayTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.baseline, self.prepared = self.root/'baseline', self.root/'prepared'
        self.baseline.mkdir(); self.prepared.mkdir()
        (self.baseline/'reference.faa').write_text('>R1\nMKVL\n')
        (self.baseline/'reference_labels.tsv').write_text('accession\tfunction_id\tevidence\nR1\tGO:1\tEXP\nR1\tEC:1.1.1.1\tEXP\n')
        (self.baseline/'hits.tsv').write_text('q1\tR1\t50\t80\t1e-20\t100\nq2\tR1\t20\t90\t1e-6\t50\n')
        (self.prepared/'unique_proteins.faa').write_text('>q1\nMKVL\n>q2\nMKVA\n')
        (self.prepared/'proteins.tsv').write_text('genome_id\tprotein_id\tgenome_protein_id\ng1\tq1\tp1\ng2\tq2\tp2\n')
        self.command = 'diamond blastp -q query -d db -o hits --threads 2 --evalue 1e-5 --max-target-seqs 10 --outfmt 6 qseqid sseqid pident qcovhsp evalue bitscore --sensitive'
        self.report = {'schema':'gspa.homology-baseline/1','diamond':{'command':self.command,
            'version':'diamond version 2.1.11', 'hits_sha256':R.H.sha256_file(self.baseline/'hits.tsv')},
            'reference':{'reference_faa_sha256':R.H.sha256_file(self.baseline/'reference.faa'),
                         'reference_labels_sha256':R.H.sha256_file(self.baseline/'reference_labels.tsv')}, 'modes':{}}
        for mode in ('homology','family_excluded'):
            claims = R.H.claims_from_hits(self.baseline/'hits.tsv',R.H.read_reference_labels(self.baseline/'reference_labels.tsv'),mode)
            files = R.H.write_claims_jsonl(claims,self.prepared/'proteins.tsv',self.baseline/'claims'/mode)
            (self.baseline/f'claims_{mode}.tsv').write_text('bound historical summary\n')
            self.report['modes'][mode] = {'claim_files':files,'claims_sha256':R.H.sha256_file(self.baseline/f'claims_{mode}.tsv')}
        self.write_report()

    def write_report(self):
        R.save(self.baseline/'baseline_report.json',self.report)

    def search(self, query, reference, out, **kwargs):
        shutil.copyfile(self.baseline/'hits.tsv',out/'hits.tsv')
        return dict(self.report['diamond'], hits_sha256=R.H.sha256_file(out/'hits.tsv'))

    def run_replay(self, search=None):
        with patch.object(R.shutil,'which',return_value='/bin/true'), patch.object(R.H,'run_diamond',side_effect=search or self.search):
            return R.replay(self.baseline,self.prepared,self.root/'out')

    def test_reconstruction_and_enrichment_preserve_every_claim_and_separate_modes(self):
        result = self.run_replay()
        self.assertEqual(result['status'],'verified')
        self.assertEqual({m:sum(r['records'] for r in rows.values()) for m,rows in result['output_claims'].items()},
                         {'homology':4,'family_excluded':2})
        contracts=json.loads((self.root/'out/contracts.json').read_text())
        self.assertNotEqual(contracts['homology']['identity_sha256'],contracts['family_excluded']['identity_sha256'])
        self.assertEqual(contracts['homology']['identity']['tool_version'],'diamond version 2.1.11')
        for mode, genomes in result['output_claims'].items():
            for g, entry in genomes.items():
                source=(self.baseline/'claims'/mode/(g+'.claims.jsonl')).read_text()
                restored=[]
                for line in Path(entry['path']).read_text().splitlines():
                    row=json.loads(line)
                    self.assertEqual(row['metadata'].pop('predictor_contract'),contracts[mode]['identity'])
                    restored.append(json.dumps(row,sort_keys=True)+'\n')
                self.assertEqual(''.join(restored),source)
        with self.assertRaisesRegex(ValueError,'already exists'): self.run_replay()

    def test_changed_historical_input_fails_before_starting(self):
        (self.baseline/'hits.tsv').write_text('changed\n')
        with self.assertRaisesRegex(ValueError,'changed historical'): self.run_replay()
        self.assertFalse((self.root/'out').exists())

    def test_incomplete_original_genome_catalogue_cannot_be_certified(self):
        self.report['modes']['homology']['claim_files'].pop('g2')
        self.write_report()
        with self.assertRaisesRegex(ValueError,'genome coverage'): self.run_replay()
        self.assertEqual(json.loads((self.root/'out/run.json').read_text())['status'],'failed')

    def test_changed_fresh_predictions_remain_failed_and_leave_a_receipt(self):
        def changed(*args,**kwargs):
            got=self.search(*args,**kwargs)
            p=args[2]/'hits.tsv'; p.write_text(p.read_text().replace('\t50\t80\t','\t49\t80\t'))
            return dict(got,hits_sha256=R.H.sha256_file(p))
        with self.assertRaisesRegex(ValueError,'reconstructed claims differ'): self.run_replay(changed)
        self.assertEqual(json.loads((self.root/'out/run.json').read_text())['status'],'failed')
        self.assertTrue((self.root/'out/failure.log').is_file())

    def test_failed_search_does_not_produce_contracts(self):
        def failed(*args,**kwargs): raise RuntimeError('solver/tool unavailable')
        with self.assertRaisesRegex(RuntimeError,'unavailable'): self.run_replay(failed)
        self.assertFalse((self.root/'out/contracts.json').exists())
        self.assertEqual(json.loads((self.root/'out/run.json').read_text())['status'],'failed')

    def test_contract_changes_with_tool_reference_parameters_adapter_and_mode(self):
        params=R.search_parameters(self.command)
        args=[self.report['diamond'],'a'*64,{'sequence_sha256':'b'*64,'annotation_sha256':'c'*64},params,'homology','d'*64]
        original=R.contract(*args)[0]
        for index,value in [(0,dict(args[0],version='diamond version 9')), (1,'e'*64),
                            (2,dict(args[2],annotation_sha256='f'*64)), (3,dict(params,max_hits=20)),
                            (4,'family_excluded'), (5,'e'*64)]:
            changed=copy.deepcopy(args);changed[index]=value
            self.assertNotEqual(R.digest(original),R.digest(R.contract(*changed)[0]))

    def test_unrecognized_commands_are_not_silently_normalized(self):
        for command in [self.command+' --very-sensitive',self.command.replace('1e-5','1e-3'),
                        self.command.replace('--threads 2','--threads 0'),self.command.replace('qcovhsp','qcov'),
                        self.command.replace('--sensitive','')]:
            with self.subTest(command=command), self.assertRaisesRegex(ValueError,'unrecognized'):
                R.search_parameters(command)


if __name__=='__main__': unittest.main()
