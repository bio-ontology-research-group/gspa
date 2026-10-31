import csv
import json
from pathlib import Path
import sys
import subprocess
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import clean_contract as C
import run_neural_predictors as N


class CleanContractTest(unittest.TestCase):
    def setUp(self):
        base = Path.home() / 'tmp';base.mkdir(exist_ok=True)
        self.work = tempfile.TemporaryDirectory(dir=base);self.addCleanup(self.work.cleanup)
        self.root = Path(self.work.name)
        self.src, self.dst = self.root/'native.csv', self.root/'query.clean.tsv'
        self.provenance = {'driver_sha256': 'd'*64, 'driver_score_kind': 'gmm_confidence', 'input_fasta_sha256': 'a'*64}

    def rewrite(self, data, kind='gmm_confidence', cutoff=0, ids=None):
        self.src.write_text(data)
        return C.rewrite(self.src,self.dst,cutoff,score_kind=kind,protein_ids=ids or {'p':'original|with/slash'},
                         provenance={**self.provenance,'driver_score_kind':kind})

    def output(self):
        # metadata is raw JSON in the fifth TSV column, not CSV quote-escaped.
        return [line.split('\t') for line in self.dst.read_text().splitlines()[1:]]

    def test_confidence_precision_and_original_identity_are_preserved(self):
        receipt=self.rewrite('p,EC:1.1.1.1/0.00000001,EC:2.2.2.2/0.00000002\n')
        rows=self.output()
        self.assertEqual([float(r[2]) for r in rows],[1e-8,2e-8])
        self.assertEqual(rows[0][0],'original|with/slash')
        self.assertEqual(receipt['emitted_calls'],2)
        self.assertTrue(json.loads(rows[0][4])['score_available'])

    def test_native_provisional_ec_identifiers_are_preserved(self):
        self.rewrite('p,EC:1.1.1.n11/0.6\n')
        self.assertEqual(self.output()[0][1], 'EC:1.1.1.n11')
        with self.assertRaises(ValueError):self.rewrite('p,EC:1.n1.1.11/0.6\n')

    def test_generic_ensemble_does_not_rank_unscored_distances_as_confidence(self):
        import build_ensemble_preds as ensemble
        self.rewrite('p,EC:1.1.1.1/0.01\n','distance')
        self.assertEqual(ensemble.fuse_rank([ensemble.load_one(self.dst)]), {})
        self.rewrite('p,EC:1.1.1.1/0.8\n')
        self.assertEqual(list(ensemble.load_one(self.dst).values()), [.8])

    def test_generic_claim_converter_rejects_clean_before_writing(self):
        script=Path(__file__).resolve().parents[3]/'gspa-cli/src/main/resources/visualize/sidecar_to_claims.py'
        target=self.root/'claims.jsonl'
        self.rewrite('p,EC:1.1.1.1/0.01\n','distance')
        for source in ('clean','other'):
            result=subprocess.run([sys.executable,str(script),'--input',source,str(self.dst),'--output',str(target)],capture_output=True,text=True)
            self.assertNotEqual(result.returncode,0)
            self.assertIn('clean-claims',result.stderr)
            self.assertFalse(target.exists())
        self.dst.write_text('protein_id\tterm\tscore\tannotation_type\np\tEC:1.1.1.1\t8.0\tEC\n')
        result=subprocess.run([sys.executable,str(script),'--input','clean',str(self.dst),'--output',str(target)],capture_output=True,text=True)
        self.assertNotEqual(result.returncode,0)
        self.assertFalse(target.exists())

    def test_distances_are_retained_without_fabricated_confidence(self):
        self.rewrite('p,EC:1.1.1.1/0.01,EC:2.2.2.2/8.0627\n','distance')
        rows=self.output()
        self.assertEqual([float(r[2]) for r in rows],[0.,0.])
        meta=[json.loads(r[4]) for r in rows]
        self.assertEqual([m['native_score'] for m in meta],[.01,8.0627])
        self.assertTrue(all(m['score_available'] is False and m['native_score_direction']=='lower' for m in meta))
        with self.assertRaises(ValueError):self.rewrite('p,EC:1.1.1.1/0.01\n','distance',.1)

    def test_malformed_empty_nonfinite_duplicate_unknown_and_incomplete_outputs_fail(self):
        fixtures=['','p\n','p,not-an-ec/0.5\n','p,EC:1.1.1.1/nan\n','p,EC:1.1.1.1/inf\n',
                  'p,EC:1.1.1.1/-0.1\n','p,EC:1.1.1.1/1.1\n','p,EC:1.1.1.1/bad\n',
                  'q,EC:1.1.1.1/0.1\n','p,EC:1.1.1.1/0.1,EC:1.1.1.1/0.2\n',
                  'p,EC:1.1.1.1/0.1\np,EC:2.2.2.2/0.1\n']
        for text in fixtures:
            with self.subTest(text=text),self.assertRaises(ValueError):self.rewrite(text)
            self.assertFalse(self.dst.exists())
            self.assertFalse(self.dst.with_suffix('.receipt.json').exists())
        with self.assertRaisesRegex(ValueError,'incomplete'):
            self.rewrite('p,EC:1.1.1.1/0.1\n',ids={'p':'one','q':'two'})

    def test_valid_but_filtered_output_has_explicit_coverage_receipt(self):
        receipt=self.rewrite('p,EC:1.1.1.1/0.01\n',cutoff=.5)
        self.assertEqual(self.output(),[])
        self.assertEqual((receipt['input_proteins'],receipt['native_selected_calls'],receipt['filtered_calls']),(1,1,1))

    def test_driver_mode_is_explicit_and_checked(self):
        driver=self.root/'driver.py'
        driver.write_text("infer_maxsep('split100', data, pretrained=True, gmm='./data/gmm.pkl')\n")
        C.validate_mode(driver,'gmm_confidence')
        for kind in (None,'distance'):
            with self.assertRaises(ValueError):C.validate_mode(driver,kind)
        driver.write_text("infer_maxsep('split100', data, gmm=None)\n")
        C.validate_mode(driver,'distance')
        driver.write_text("infer_maxsep('split100', data, gmm=choose_model())\n")
        with self.assertRaises(ValueError):C.validate_mode(driver,'gmm_confidence')

    def setup_run(self):
        app=self.root/'app';app.mkdir()
        (app/'CLEAN_infer_fasta.py').write_text("infer_maxsep('split100', data, gmm='./data/gmm.pkl')\n")
        fasta=self.root/'input.faa';fasta.write_text('>p|a/unsafe description\nMALW\n>p2\nMEKL\n')
        out=self.root/'out';out.mkdir()
        row=N.ManifestRow('query',fasta,out)
        args=SimpleNamespace(model_dir=str(app),clean_score_kind='gmm_confidence',min_score=0.)
        return app,row,args

    def test_invalid_cutoff_fails_before_native_execution(self):
        app,row,args=self.setup_run();args.min_score=float('nan')
        with patch.object(C.subprocess,'run') as native,self.assertRaises(ValueError):
            C.run([row],args)
        native.assert_not_called()

    def test_shared_query_files_are_never_reused(self):
        app,row,args=self.setup_run()
        stale=app/'results/inputs/query_maxsep.csv';stale.parent.mkdir(parents=True)
        stale.write_text('p|a/unsafe,EC:1.1.1.1/0.9\np2,EC:1.1.1.1/0.9\n')
        with patch.object(C.subprocess,'run'),self.assertRaisesRegex(ValueError,'this run'):
            C.run([row],args)
        self.assertTrue(stale.exists())
        self.assertFalse((row.output_dir/'query.clean.tsv').exists())
        self.assertEqual(list((app/'data/inputs').iterdir()),[])

    def test_runs_use_unique_protein_and_file_names_and_restore_ids(self):
        app,row,args=self.setup_run();names=[]
        unrelated=app/'data/esm_data/existing.pt';unrelated.parent.mkdir(parents=True);unrelated.write_text('leave alone')
        def native(cmd,check,cwd):
            self.assertEqual(cmd[0],sys.executable)
            name=cmd[-1];names.append(name)
            staged=app/'data/inputs'/f'{name}.fasta'
            proteins=C.read_proteins(staged)
            self.assertTrue(all(pid.startswith(name+'_p') for pid,_ in proteins))
            result=app/'results/inputs'/f'{name}_maxsep.csv'
            result.write_text(''.join(pid+',EC:1.1.1.1/0.8\n' for pid,_ in proteins))
            for pid,_ in proteins:(app/'data/esm_data'/f'{pid}.pt').write_text('controlled cache fixture')
        with patch.object(C.subprocess,'run',side_effect=native):
            C.run([row],args);C.run([row],args)
        self.assertNotEqual(names[0],names[1])
        rows=(row.output_dir/'query.clean.tsv').read_text().splitlines()[1:]
        self.assertEqual([r.split('\t')[0] for r in rows],['p|a/unsafe','p2'])
        self.assertEqual(list((app/'data/inputs').iterdir()),[])
        self.assertEqual(list((app/'data/esm_data').iterdir()),[unrelated])
        self.assertEqual(list((app/'results/inputs').iterdir()),[])

    def test_fasta_duplicates_and_missing_sequences_fail(self):
        for text in ('>p\nAA\n>p\nBB\n','>p\n','>\nAA\n','AA\n>p\nBB\n'):
            p=self.root/'bad.faa';p.write_text(text)
            with self.assertRaises(ValueError):C.read_proteins(p)

    def test_empty_input_is_explicit_and_does_not_execute_a_model(self):
        app,row,args=self.setup_run()
        row.fasta_path.write_text('')
        with patch.object(C.subprocess,'run') as native:
            C.run([row],args)
        native.assert_not_called()
        receipt=json.loads((row.output_dir/'query.clean.receipt.json').read_text())
        self.assertEqual(receipt['input_state'],'empty')
        self.assertEqual(receipt['emitted_calls'],0)

    def test_defaults_preserve_other_predictors_and_keep_all_clean_selected_calls(self):
        for name,expected in [('clean',0.),('proteinfer',.1)]:
            with patch.object(sys,'argv',['tool','--predictor',name,'--manifest','manifest.tsv']):
                self.assertEqual(N.parse_args().min_score,expected)


if __name__=='__main__':unittest.main()
