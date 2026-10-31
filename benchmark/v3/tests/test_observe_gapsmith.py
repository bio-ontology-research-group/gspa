"""Synthetic execution-contract checks; not evidence of real tool accuracy."""
import contextlib
import io
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import observe_gapsmith as O


class ObservationContract(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='gspa-observe-')
        self.addCleanup(self.tmp.cleanup)
        self.base = Path(self.tmp.name)
        self.root = self.base/'out'; self.root.mkdir()
        self.scratch = self.base/'scratch'; self.scratch.mkdir()
        self.ref = self.base/'ref'; self.ref.mkdir()
        for name in O.TABLES: (self.ref/name).write_text('synthetic table\n')
        for tax in ['Bacteria','Archaea']:
            for sub in ['rev','unrev','rxn']:
                d = self.ref/'seq'/tax/sub; d.mkdir(parents=True)
                (d/'r.fasta').write_text('>ref\nMAAA\n')
        self.proteome = self.base/'different_genome.faa'
        self.proteome.write_text('>p1\nMAAA\n')
        self.tool = self.base/'tool'
        self.tool.write_text('''#!/usr/bin/env python3
import os,sys
from pathlib import Path
if sys.argv[1:] == ['--version']:
 print('synthetic gapsmith fixture');sys.exit(0)
args=sys.argv[1:]
tax=args[args.index('-t')+1]
ref=Path(args[args.index('--data-dir')+1])
assert tax in ('Bacteria','Archaea')
assert (ref/'seq'/tax).is_dir()
# The installed 0.1.0 resolver requires this marker even for Archaea.
assert (ref/'seq'/'Bacteria').is_dir()
other=ref/'seq'/('Bacteria' if tax=='Archaea' else 'Archaea')
assert not list(other.rglob('*.fasta'))
assert args[args.index('-K')+1] == '2'
if args[args.index('-p')+1] == '^$': sys.exit(0)
if os.environ.get('FIXTURE_FAIL'): sys.exit(7)
out=Path(args[args.index('-o')+1]);out.mkdir()
stem=Path(args[-1]).stem
(out/(stem+'-all-Reactions.tbl')).write_text('status\\ngood_blast\\n')
(out/(stem+'-all-Pathways.tbl')).write_text('Prediction\\nTRUE\\nFALSE\\n')
''')
        self.tool.chmod(0o755)
        self.diamond = self.base/'diamond'
        self.diamond.write_text('#!/bin/sh\necho synthetic-diamond\n');self.diamond.chmod(0o755)
        self.env = patch.dict(os.environ, SLURM_JOB_ID='synthetic', SLURM_TMPDIR=str(self.scratch),
                              SLURM_CPUS_PER_TASK='2')
        self.env.start();self.addCleanup(self.env.stop)
        # GNU time is present on the Slurm workers but optional on developer
        # machines. Only replace that resource wrapper; execute the fake tool
        # as an actual subprocess so argument, copy and failure checks stay real.
        popen = O.subprocess.Popen
        def without_time(command, *args, **kwargs):
            if command[0] == '/usr/bin/time':
                Path(command[3]).write_text('synthetic resource wrapper\n')
                command = command[4:]
            return popen(command, *args, **kwargs)
        wrapper = patch.object(O.subprocess, 'Popen', side_effect=without_time)
        wrapper.start();self.addCleanup(wrapper.stop)

    def execute(self, taxonomy='Archaea', **kwargs):
        with contextlib.redirect_stdout(io.StringIO()):
            O.run(self.root,self.ref,self.tool,self.diamond,O.digest(self.proteome),
                  proteome=self.proteome,taxonomy=taxonomy,threads=2,**kwargs)
        return json.loads((self.root/'observation.json').read_text())

    def test_arbitrary_genome_uses_declared_taxonomy_and_threads(self):
        r=self.execute()
        self.assertEqual(r['status'],'observed')
        self.assertEqual(r['taxonomy'],'Archaea')
        self.assertEqual(r['threads'],2)
        self.assertEqual(r['reaction_statuses'],{'good_blast':1})
        self.assertEqual(r['predicted_pathways'],1)
        self.assertEqual(r['proteome_sha256'],O.digest(self.proteome))
        self.assertIn('outputs/different_genome-all-Reactions.tbl',r['output_hashes'])
        refs=json.loads((self.root/'reference-hashes.json').read_text())
        self.assertIn('seq/Archaea/rev/r.fasta',refs)
        self.assertFalse(any(p.startswith('seq/Bacteria/') for p in refs))

    def test_bacterial_reference_route(self):
        self.assertEqual(self.execute('Bacteria')['status'],'observed')

    def test_failure_is_recorded_and_existing_observation_cannot_be_overwritten(self):
        with patch.dict(os.environ,FIXTURE_FAIL='1'):
            with self.assertRaises(RuntimeError): self.execute()
        receipt=self.root/'observation.json'; before=receipt.read_bytes()
        self.assertEqual(json.loads(before)['status'],'failed')
        self.assertEqual(json.loads(before)['exit_code'],7)
        with self.assertRaises(FileExistsError): self.execute()
        self.assertEqual(receipt.read_bytes(),before)

    def test_input_hash_mismatch_is_recorded(self):
        with self.assertRaisesRegex(ValueError,'proteome changed'):
            O.run(self.root,self.ref,self.tool,self.diamond,'wrong',self.proteome,threads=2)
        self.assertEqual(json.loads((self.root/'observation.json').read_text())['status'],'failed')

    def test_outside_slurm_and_excess_threads_are_rejected(self):
        with patch.dict(os.environ,SLURM_JOB_ID=''):
            with self.assertRaises(RuntimeError): self.execute()
        with self.assertRaisesRegex(ValueError,'allocation'):
            O.run(self.root,self.ref,self.tool,self.diamond,O.digest(self.proteome),threads=3)
        self.assertFalse((self.root/'observation.json').exists())

    def test_user_only_reference_is_not_a_complete_database(self):
        import shutil
        shutil.rmtree(self.ref/'seq/Archaea')
        user=self.ref/'seq/Archaea/user';user.mkdir(parents=True)
        (user/'small.faa').write_text('>ref\nMAAA\n')
        with self.assertRaisesRegex(ValueError,'missing reference category'):
            self.execute()
        self.assertEqual(json.loads((self.root/'observation.json').read_text())['status'],'failed')


if __name__ == '__main__': unittest.main()
