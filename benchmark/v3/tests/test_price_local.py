import json
import os
from pathlib import Path
import sys
import subprocess
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import run_price_local as L


class LocalPublicationTests(unittest.TestCase):
    def setUp(self):
        scratch=Path(os.environ.get('TMPDIR') or Path.home()/'tmp')
        scratch.mkdir(parents=True,exist_ok=True)
        temporary=tempfile.TemporaryDirectory(dir=scratch)
        self.addCleanup(temporary.cleanup)
        self.root=Path(temporary.name)
        self.local=self.root/'scratch'/'job';self.local.mkdir(parents=True)
        self.destination=self.root/'collection'/'attempts'/'0'/'job'
        self.staging=self.root/'publication-staging'
        self.original=b'{"status":"verified","unchanged":"receipt bytes"}\n'
        (self.local/'completed.json').write_bytes(self.original)
        (self.local/'prediction.tsv').write_text('complete candidate output\n')

    def test_atomic_publication_preserves_all_original_bytes_and_paths(self):
        before=L.tree_hashes(self.local)
        receipt=L.publish(self.local,self.destination,self.staging,dict(returncode=0))
        report=json.loads(receipt.read_text())
        self.assertEqual(report['status'],'verified')
        self.assertEqual(report['published_path'],str(self.destination))
        self.assertEqual(report['artifacts'],before)
        self.assertEqual((self.destination/'completed.json').read_bytes(),self.original)
        self.assertEqual(L.tree_hashes(self.local),before)
        self.assertEqual(list(self.staging.iterdir()),[])
        self.assertEqual(len(list((self.root/'collection/attempts').glob('*/*/completed.json'))),1)

    def test_existing_attempt_cannot_be_replaced(self):
        self.destination.mkdir(parents=True)
        with self.assertRaises(FileExistsError):
            L.publish(self.local,self.destination,self.staging,dict(returncode=0))
        self.assertEqual(list(self.destination.iterdir()),[])

    def test_corruption_during_copy_never_exposes_completion(self):
        original=L.shutil.copytree
        def corrupt(source,target,**kwargs):
            original(source,target,**kwargs)
            (target/'prediction.tsv').write_text('truncated')
        with patch.object(L.shutil,'copytree',side_effect=corrupt),self.assertRaisesRegex(ValueError,'bytes changed'):
            L.publish(self.local,self.destination,self.staging,dict(returncode=0))
        self.assertFalse(self.destination.exists())
        self.assertEqual(list((self.root/'collection/attempts').glob('*/*/completed.json')),[])
        self.assertTrue((self.local/'completed.json').exists())
        self.assertEqual(len(list(self.staging.iterdir())),1)

    def test_source_mutation_during_copy_is_rejected(self):
        original=L.shutil.copytree
        def mutate(source,target,**kwargs):
            original(source,target,**kwargs)
            (source/'prediction.tsv').write_text('changed source')
        with patch.object(L.shutil,'copytree',side_effect=mutate),self.assertRaisesRegex(ValueError,'bytes changed'):
            L.publish(self.local,self.destination,self.staging,dict(returncode=0))
        self.assertFalse(self.destination.exists())

    def test_failed_attempt_is_retained_without_completion(self):
        (self.local/'completed.json').unlink()
        (self.local/'failed.json').write_text('{"status":"failed"}\n')
        report=json.loads(L.publish(self.local,self.destination,self.staging,dict(returncode=1)).read_text())
        self.assertEqual(report['status'],'failed_execution')
        self.assertFalse((self.destination/'completed.json').exists())
        self.assertTrue((self.destination/'failed.json').exists())

    def test_contradictory_or_missing_completion_is_rejected(self):
        with self.assertRaisesRegex(ValueError,'contradictory'):
            L.publish(self.local,self.destination,self.staging,dict(returncode=1))
        (self.local/'failed.json').write_text('{}')
        with self.assertRaisesRegex(ValueError,'contradictory'):
            L.publish(self.local,self.destination,self.staging,dict(returncode=0))
        (self.local/'completed.json').unlink()
        with self.assertRaisesRegex(ValueError,'without a completed'):
            L.publish(self.local,self.destination,self.staging,dict(returncode=0))

    def test_staging_inside_attempt_glob_and_overlapping_source_are_rejected(self):
        for staging in (self.destination.parent/'staged',self.destination.parent.parent/'staged',self.local/'staged'):
            with self.subTest(staging=staging),self.assertRaisesRegex(ValueError,'staging'):
                L.publish(self.local,self.destination,staging,dict(returncode=0))
        with self.assertRaisesRegex(ValueError,'separate'):
            L.publish(self.local,self.local/'nested',self.staging,dict(returncode=0))

    def test_symlink_artifacts_are_not_followed(self):
        (self.local/'escape').symlink_to(self.root)
        with self.assertRaisesRegex(ValueError,'symlink'):
            L.publish(self.local,self.destination,self.staging,dict(returncode=0))
        self.assertFalse(self.destination.exists())

    def test_binding_preserves_logical_production_paths_and_scratch(self):
        out=Path('/data/project/collection');scratch=Path('/scratch/123')
        cmd=L.command(out,7,'job-123','24g',7200,scratch,scratch/'output',Path('/container'),Path('/usr/bin/singularity'))
        self.assertIn('/scratch/123/output:/data/project/collection/attempts/7',cmd)
        self.assertEqual(cmd[cmd.index('--out')+1],str(out))
        self.assertEqual(cmd[cmd.index('--scratch')+1],str(scratch))
        self.assertEqual(cmd[cmd.index('--env')+1],'TMPDIR='+str(scratch))
        with self.assertRaisesRegex(ValueError,'delimiter'):
            L.command(out,7,'job','24g',7200,scratch,Path('/bad,path'),Path('/container'),Path('/bin/singularity'))

    def test_cross_filesystem_publication_is_rejected_before_execution(self):
        L.preflight(self.local,self.destination,self.staging)
        original=os.stat
        def other_device(path,*args,**kwargs):
            record=original(path,*args,**kwargs)
            if Path(path)==self.staging:
                values=list(record);values[2]+=1
                return os.stat_result(values)
            return record
        with patch.object(L.os,'stat',side_effect=other_device),self.assertRaisesRegex(ValueError,'share a filesystem'):
            L.preflight(self.local,self.destination,self.staging)
        self.assertFalse(self.destination.exists())

    def test_sync_failure_retains_bytes_without_publishing_completion(self):
        with patch.object(L,'sync_tree',side_effect=OSError('writeback failed')),self.assertRaisesRegex(OSError,'writeback'):
            L.publish(self.local,self.destination,self.staging,dict(returncode=0))
        self.assertFalse(self.destination.exists())
        self.assertTrue((self.local/'completed.json').exists())
        self.assertEqual(list((self.root/'collection/attempts').glob('*/*/completed.json')),[])

    def test_inactive_bind_fails_before_engine_can_write_to_shared_paths(self):
        script=self.root/'stub.py';marker=self.root/'ran.txt'
        script.write_text('from pathlib import Path\nimport sys\nPath(sys.argv[1]).write_text("ran")\n')
        logical=self.root/'logical';logical.mkdir()
        args=[sys.executable,'-c',L.MOUNT_GUARD,str(self.local),str(logical),str(script),str(marker)]
        failed=subprocess.run(args,capture_output=True,text=True)
        self.assertNotEqual(failed.returncode,0)
        self.assertIn('private output bind is not active',failed.stderr)
        self.assertFalse(marker.exists())
        args[4]=str(self.local)
        passed=subprocess.run(args,capture_output=True,text=True)
        self.assertEqual(passed.returncode,0,passed.stderr)
        self.assertEqual(marker.read_text(),'ran')

    def test_early_engine_failure_keeps_durable_returncode_without_completion(self):
        out=self.root/'collection';out.mkdir()
        L.P.save(out/'plan.json',dict(cells=[{}]))
        L.P.save(out/'registration.json',dict(inputs={},
            plan_sha256=L.P.C.V.sha256_file(str(out/'plan.json'))))
        container=self.root/'container';container.mkdir()
        with patch.object(L.shutil,'which',return_value=sys.executable), \
             patch.object(L.subprocess,'run',return_value=subprocess.CompletedProcess([],7)), \
             self.assertRaisesRegex(ValueError,'missing local attempt'):
            L.run(out,0,'job-bad','24g',10,self.root/'scratch',container)
        result=json.loads((self.root/'local-executions/0-job-bad.result.json').read_text())
        self.assertEqual(result['returncode'],7)
        self.assertEqual(list((out/'attempts').glob('*/*/completed.json')),[])
        self.assertTrue((self.root/'scratch/price-output-0-job-bad').is_dir())


if __name__=='__main__':unittest.main()
