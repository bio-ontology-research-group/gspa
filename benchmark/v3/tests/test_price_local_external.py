"""External (host-executed) launch of a registered Price cell, e.g. on Shaheen.

The wrapper records the engine command inside a container; the host executes it;
the wrapper then publishes with the host return code. The default launch path
and its command are unchanged.
"""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import run_price_local as L

SCRIPT = Path(L.__file__).resolve()


class ExternalLaunchTests(unittest.TestCase):
    def setUp(self):
        scratch=Path(os.environ.get('TMPDIR') or Path.home()/'tmp')
        scratch.mkdir(parents=True,exist_ok=True)
        temporary=tempfile.TemporaryDirectory(dir=scratch)
        self.addCleanup(temporary.cleanup)
        self.root=Path(temporary.name).resolve()
        self.out=self.root/'collection';self.out.mkdir()
        L.P.save(self.out/'plan.json',dict(cells=[{}]))
        L.P.save(self.out/'registration.json',dict(inputs={},
            plan_sha256=L.P.C.V.sha256_file(str(self.out/'plan.json'))))
        self.container=self.root/'image';(self.container/'usr/bin').mkdir(parents=True)
        (self.container/'usr/bin/env').write_bytes(b'binary')
        (self.container/'bin').symlink_to('usr/bin')
        (self.container/'dangling').symlink_to('/nonexistent/target')
        self.launcher=self.root/'singularity';self.launcher.write_bytes(b'launcher')
        self.prefix=self.root/'prefix.json'
        self.prefix.write_text(json.dumps([str(self.launcher),'exec','--no-home','--bind','/a:/b:ro']))
        self.scratch=self.root/'scratch';self.scratch.mkdir()
        self.environment=patch.dict(os.environ,{k:'' for k in ('JAVA_TOOL_OPTIONS','JDK_JAVA_OPTIONS','_JAVA_OPTIONS')})
        self.environment.start();self.addCleanup(self.environment.stop)

    def prepare(self,attempt='job-1'):
        return L.prepare_external(self.out,0,attempt,'24g',10,self.scratch,self.container,self.prefix)

    def engine_completes(self,attempt='job-1',completed=True):
        """What the recorded command does: write the attempt into the private parent."""
        local=self.scratch/('price-output-0-'+attempt)/attempt;local.mkdir()
        (local/'started.json').write_text('{}\n')
        (local/('completed.json' if completed else 'failed.json')).write_text('{"status":"x"}\n')
        return local

    def test_external_command_matches_default_engine_arguments(self):
        out=Path('/data/project/collection');scratch=Path('/scratch/123')
        default=L.command(out,7,'job-123','24g',7200,scratch,scratch/'output',Path('/image'),Path('/usr/bin/singularity'))
        prefix=['/sw/singularity','exec','--no-home','--bind','/x:/y:ro']
        external=L.external_command(prefix,out,7,'job-123','24g',7200,scratch,scratch/'output',Path('/image'))
        self.assertEqual(external[:len(prefix)],prefix)
        # Everything from --pwd onward, minus the host /usr,... bind of the default skeleton launch, is identical.
        tail=default[default.index('--pwd'):]
        host=tail.index('/usr,/lib,/lib64,/bin,/etc,/data');del tail[host-1:host+1]
        self.assertEqual(external[len(prefix):],tail)
        with self.assertRaisesRegex(ValueError,'delimiter'):
            L.external_command(prefix,out,7,'job','24g',7200,scratch,Path('/bad:path'),Path('/image'))

    def test_prefix_validation(self):
        for value in ([],'singularity',['relative/singularity','exec'],[str(self.launcher),''],[str(self.launcher),'a\nb'],[str(self.launcher),3]):
            with self.subTest(value=value):
                self.prefix.write_text(json.dumps(value))
                with self.assertRaisesRegex(ValueError,'launch prefix'):
                    L.load_prefix(self.prefix)

    def test_prepare_records_without_launching_and_publish_binds_host_returncode(self):
        with patch.object(L.subprocess,'run',side_effect=AssertionError('must not launch')):
            locator=self.prepare()
        self.assertEqual(Path(locator),self.root/'local-executions/0-job-1.json')
        execution=json.loads(Path(locator).read_text())
        self.assertEqual(execution['launch'],'external')
        self.assertEqual(execution['command'][:5],json.loads(self.prefix.read_text()))
        self.assertIn(str(self.container),execution['command'])
        self.assertEqual(execution['container_files']['bin'],'symlink:usr/bin')
        self.assertEqual(execution['container_files']['dangling'],'symlink:/nonexistent/target')
        self.assertEqual(set(execution['runtime_files']),{str(SCRIPT),str(self.launcher)})
        self.assertNotIn('returncode',execution)
        self.assertFalse((self.out/'attempts/0/job-1').exists())
        self.engine_completes()
        L.publish_external(self.out,0,'job-1',self.scratch,self.container,0)
        publication=json.loads((self.out/'attempts/0/job-1/publication.json').read_text())
        self.assertEqual(publication['status'],'verified')
        self.assertEqual(publication['execution']['returncode'],0)
        self.assertEqual(publication['execution']['command'],execution['command'])
        result=json.loads((self.root/'local-executions/0-job-1.result.json').read_text())
        self.assertEqual(result,dict(returncode=0,locator_sha256=L.P.C.V.sha256_file(str(locator))))
        self.assertFalse((self.scratch/'price-output-0-job-1').exists())
        self.assertEqual(len(list((self.out/'attempts').glob('*/*/completed.json'))),1)

    def test_failed_host_execution_is_published_as_failure_and_retained(self):
        self.prepare();self.engine_completes(completed=False)
        with self.assertRaisesRegex(RuntimeError,'production failed \\(3\\)'):
            L.publish_external(self.out,0,'job-1',self.scratch,self.container,3)
        publication=json.loads((self.out/'attempts/0/job-1/publication.json').read_text())
        self.assertEqual(publication['status'],'failed_execution')
        self.assertEqual(list((self.out/'attempts').glob('*/*/completed.json')),[])
        self.assertTrue((self.scratch/'price-output-0-job-1/job-1/failed.json').exists())

    def test_success_without_completion_and_changed_runtime_are_rejected(self):
        self.prepare()
        (self.scratch/'price-output-0-job-1/job-1').mkdir()
        with self.assertRaisesRegex(ValueError,'without a completed'):
            L.publish_external(self.out,0,'job-1',self.scratch,self.container,0)
        self.prepare('job-2');self.engine_completes('job-2')
        (self.container/'usr/bin/env').write_bytes(b'changed')
        with self.assertRaisesRegex(ValueError,'runtime changed'):
            L.publish_external(self.out,0,'job-2',self.scratch,self.container,0)
        self.assertFalse((self.out/'attempts/0/job-2').exists())
        self.assertTrue((self.scratch/'price-output-0-job-2/job-2/completed.json').exists())

    def test_publish_requires_matching_external_locator_and_integer_returncode(self):
        with self.assertRaises(FileNotFoundError):
            L.publish_external(self.out,0,'job-1',self.scratch,self.container,0)
        self.prepare()
        for value in (True,'0',None):
            with self.subTest(value=value),self.assertRaisesRegex(ValueError,'integer return code'):
                L.publish_external(self.out,0,'job-1',self.scratch,self.container,value)
        other=self.root/'other-scratch';other.mkdir()
        with self.assertRaisesRegex(ValueError,'does not describe'):
            L.publish_external(self.out,0,'job-1',other,self.container,0)
        locator=self.root/'local-executions/0-job-1.json'
        execution=json.loads(locator.read_text());del execution['launch']
        locator.write_text(json.dumps(execution))
        with self.assertRaisesRegex(ValueError,'does not describe'):
            L.publish_external(self.out,0,'job-1',self.scratch,self.container,0)
        self.assertFalse((self.root/'local-executions/0-job-1.result.json').exists())

    def test_prepare_refuses_reused_attempt_and_missing_launcher(self):
        self.prepare()
        with self.assertRaises(FileExistsError):
            self.prepare()
        self.prefix.write_text(json.dumps([str(self.root/'missing-launcher'),'exec']))
        with self.assertRaisesRegex(ValueError,'launcher visible'):
            self.prepare('job-3')
        self.assertFalse((self.root/'local-executions/0-job-3.json').exists())

    def test_cli_default_dispatch_is_unchanged_and_mixed_options_are_refused(self):
        base=[sys.executable,str(SCRIPT),'--out',str(self.out),'--scratch',str(self.scratch),
              '--container',str(self.container),'--index','0','--attempt','job-9','--heap','24g','--timeout','10']
        env={**os.environ,'PYTHONPATH':str(SCRIPT.parent),'PATH':'/nonexistent'}
        # The default path still requires a Singularity on PATH and launches nothing else.
        default=subprocess.run(base,capture_output=True,text=True,env=env)
        self.assertNotEqual(default.returncode,0)
        self.assertIn('Singularity and a container root are required',default.stderr)
        for extra in (['--phase','prepare'],['--returncode','0'],['--launch-prefix',str(self.prefix)],
                      ['--phase','publish','--launch-prefix',str(self.prefix),'--returncode','0']):
            with self.subTest(extra=extra):
                mixed=subprocess.run(base+extra,capture_output=True,text=True,env=env)
                self.assertEqual(mixed.returncode,2,mixed.stderr)
                self.assertIn('external launch',mixed.stderr)
        prepared=subprocess.run(base[:11]+['job-8']+base[12:]+['--phase','prepare','--launch-prefix',str(self.prefix)],
                                capture_output=True,text=True,env=env)
        self.assertEqual(prepared.returncode,0,prepared.stderr)
        self.assertEqual(prepared.stdout.strip(),str(self.root/'local-executions/0-job-8.json'))


if __name__=='__main__':unittest.main()
