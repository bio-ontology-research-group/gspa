import copy
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import run_price_production_collection as P


def fixture():
    jobs=[dict(predictor=p,assembly=g,matrix_rows=3,expected_claims=n+1)
          for n,(p,g) in enumerate([('a','g1'),('a','g2'),('b','g1'),('b','g2')])]
    export=dict(out='/export',jobs=jobs,expected_claims=10)
    parts=[dict(index=i,predictor=j['predictor'],assembly=j['assembly'],matrix_rows=3,
                claims=j['expected_claims'],receipt_sha256=str(i)) for i,j in enumerate(jobs)]
    audit=dict(status='verified',n_parts=4,claims=10,parts=parts,
               hashes={f"/export/parts/{j['predictor']}/{j['assembly']}/job-receipt.json":str(i) for i,j in enumerate(jobs)})
    return export,audit


class ProductionCollectionTests(unittest.TestCase):
    def completed_fixture(self, root):
        cell=dict(index=0,part='/part',claims=1,assembly='g',predictor='topec',matrix_rows=1)
        P.save(root/'plan.json',dict(cells=[cell],arms=['stage0','stage1'],layout='/layout',
                                   jar='/jar',reference='/ref',expected_claims=1))
        P.save(root/'registration.json',dict(inputs={},plan_sha256=P.C.V.sha256_file(str(root/'plan.json'))))
        def engine(part,layout,jar,reference,out,**kwargs):
            out.mkdir();runs={}
            for arm in kwargs['arms']:
                directory=out/arm;directory.mkdir();tsv=directory/'g.tsv'
                tsv.write_text('protein_id\ttype\tfunction_id\taspect\tscore\np\tEC\t1.1.1.1\t\t.5\n')
                P.save(Path(str(tsv)+'.stages.json'),dict(schema_version=2,genome_id='g',stages=[
                    dict(name='evidence_integration',status='executed')]))
                spec=dict(outputs={'g':str(tsv)},required_stages=['evidence_integration'],inputs={},
                          command=['java']+(['--enable-priors','none'] if arm=='stage0' else []))
                receipt=directory/P.C.V.RECEIPT
                P.save(receipt,dict(spec=spec,returncode=0,artifacts=P.C.V.output_hashes(spec)))
                runs[arm]=dict(prediction_rows=1,input_candidates=1,candidate_universe_preserved=True,
                               receipt_sha256=P.C.V.sha256_file(str(receipt)))
            P.save(out/'verification.json',dict(status='verified',runs=runs,
                  assembly='g',predictor='topec',matrix_rows=1))
        with patch.object(P.C,'run',side_effect=engine):P.run(root,0,'job1','24g',100,root/'scratch')
        return root/'attempts/0/job1'

    def test_success_and_finalization_verify_both_arms_and_refuse_overwrite(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);attempt=self.completed_fixture(root)
            self.assertTrue((attempt/'completed.json').exists())
            P.finalize(root)
            result=json.loads((root/'collection-verification.json').read_text())
            self.assertEqual((result['status'],result['cells'],result['arms']),('verified',1,2))
            with self.assertRaises(FileExistsError):P.finalize(root)

    def test_partial_failed_only_and_duplicate_completed_cells_cannot_finalize(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);attempt=self.completed_fixture(root)
            plan=json.loads((root/'plan.json').read_text())
            self.assertEqual(len(P.completed_attempts(root,plan)),1)
            done=attempt/'completed.json';saved=done.read_text();done.unlink()
            P.save(attempt/'failed.json',dict(status='failed'))
            with self.assertRaisesRegex(ValueError,'found 0'):P.completed_attempts(root,plan)
            done.write_text(saved)
            with self.assertRaisesRegex(ValueError,'contradictory'):P.finalize(root)
            (attempt/'failed.json').unlink()
            other=attempt.parent/'job2';other.mkdir();(other/'completed.json').write_text(saved)
            with self.assertRaisesRegex(ValueError,'found 2'):P.completed_attempts(root,plan)
            (other/'completed.json').unlink()
            plan['cells'].append(dict(index=1))
            with self.assertRaisesRegex(ValueError,'cell 1'):P.completed_attempts(root,plan)

    def test_tampered_prediction_receipt_and_verification_are_rejected(self):
        for suffix in ['stage0/g.tsv','stage0/verified_run.json','verification.json']:
            with self.subTest(suffix=suffix),tempfile.TemporaryDirectory() as folder:
                root=Path(folder);attempt=self.completed_fixture(root)
                path=attempt/'production'/suffix;path.write_text(path.read_text()+'\n')
                with self.assertRaisesRegex(ValueError,'changed'):P.finalize(root)
                self.assertFalse((root/'collection-verification.json').exists())

    def test_complete_population_and_bound_cell_counts(self):
        export,audit=fixture();rows=P.cells(export,audit)
        self.assertEqual(len(rows),4)
        self.assertEqual(sum(r['claims'] for r in rows),10)
        self.assertEqual([r['index'] for r in rows],[0,1,2,3])

    def test_missing_duplicate_swapped_and_changed_cells_fail(self):
        for failure in ['missing','duplicate','swapped','claims','unbound','partial_audit']:
            e,a=fixture()
            if failure=='missing':e['jobs'].pop()
            elif failure=='duplicate':e['jobs'][-1]=copy.deepcopy(e['jobs'][0])
            elif failure=='swapped':a['parts'][0],a['parts'][1]=a['parts'][1],a['parts'][0]
            elif failure=='claims':a['parts'][0]['claims']+=1
            elif failure=='unbound':a['hashes'].clear()
            else:a['status']='running'
            with self.subTest(failure=failure),self.assertRaises(ValueError):P.cells(e,a)

    def test_failed_arm_keeps_attempt_and_cannot_be_certified_or_overwritten(self):
        with tempfile.TemporaryDirectory(dir=os.environ.get('TMPDIR',str(Path.home()/'tmp'))) as folder:
            root=Path(folder);cell=dict(index=0,part='/part',claims=10)
            P.save(root/'plan.json',dict(cells=[cell],arms=['stage0','stage1'],layout='/layout',jar='/jar',reference='/ref'))
            P.save(root/'registration.json',dict(inputs={},plan_sha256=P.C.V.sha256_file(str(root/'plan.json'))))
            def incomplete(part,layout,jar,reference,out,**kwargs):
                out.mkdir();P.save(out/'verification.json',dict(status='verified',runs={'stage0':{}}))
            with patch.object(P.C,'run',side_effect=incomplete),self.assertRaisesRegex(ValueError,'incomplete production arms'):
                P.run(root,0,'job1','24g',100,root/'scratch')
            attempt=root/'attempts/0/job1'
            self.assertTrue((attempt/'failed.json').exists())
            self.assertFalse((attempt/'completed.json').exists())
            with self.assertRaises(FileExistsError):P.run(root,0,'job1','24g',100,root/'scratch')

    def test_invalid_attempt_identifier_cannot_escape_output(self):
        for value in ['..','../x','a/b','']:
            with self.subTest(value=value),self.assertRaisesRegex(ValueError,'invalid attempt'):
                P.run(Path('/does-not-exist'),0,value,'24g',100,Path.home()/'tmp')


if __name__=='__main__':unittest.main()
