import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import observe_gapsmith_cohort as C


class CohortSearchContract(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='gspa-cohort-observe-')
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.proteome = self.root/'g.faa';self.proteome.write_text('>p\nMAAA\n')
        self.source = self.root/'frozen.py';self.source.write_text('# frozen\n')
        self.manifest = self.root/'manifest.json'
        self.data = dict(schema='gspa.cohort-pathway-search/1', output_root=str(self.root/'outputs'),
                         reference='reference', binary='binary', diamond='diamond', threads=2,
                         sources={str(self.source):C.observer.digest(self.source)},
                         genomes=[dict(genome_id='g', taxonomy='Archaea', proteome=str(self.proteome),
                                       sha256=C.observer.digest(self.proteome))])
        self.write()

    def write(self):
        self.manifest.write_text(json.dumps(self.data))
        self.sha = C.observer.digest(self.manifest)

    def test_bound_manifest_and_input_selection(self):
        self.assertEqual(C.selected(self.manifest,self.sha,0)[1],self.data['genomes'][0])
        for index in [-1,1]:
            with self.assertRaises(ValueError): C.selected(self.manifest,self.sha,index)
        self.data['genomes'] *= 2;self.write()
        with self.assertRaises(ValueError): C.selected(self.manifest,self.sha,0)

    def test_changed_manifest_source_and_proteome_are_rejected(self):
        for target in [self.manifest,self.source,self.proteome]:
            with self.subTest(target=target.name):
                original=target.read_bytes();target.write_bytes(original+b' ')
                with self.assertRaises(ValueError): C.selected(self.manifest,self.sha,0)
                target.write_bytes(original)

    def test_observed_execution_is_bound_and_not_overwritten(self):
        def fake_observer(root,*args):
            (root/'observation.json').write_text('{"status":"observed"}\n')
        with patch.dict(os.environ,SLURM_JOB_ID='synthetic'), patch.object(C.observer,'run',side_effect=fake_observer):
            C.run(self.manifest,self.sha,0)
            with self.assertRaises(FileExistsError): C.run(self.manifest,self.sha,0)
        result=json.loads((self.root/'outputs/g/execution.json').read_text())
        self.assertEqual(result['status'],'observed')
        self.assertEqual(result['manifest_sha256'],self.sha)
        self.assertEqual(result['observation_sha256'],C.observer.digest(self.root/'outputs/g/observation.json'))

    def test_source_mutation_during_execution_is_failed(self):
        def changing(*args): self.source.write_text('# changed\n')
        with patch.dict(os.environ,SLURM_JOB_ID='synthetic'), patch.object(C.observer,'run',side_effect=changing):
            with self.assertRaisesRegex(ValueError,'source changed'): C.run(self.manifest,self.sha,0)
        result=json.loads((self.root/'outputs/g/execution.json').read_text())
        self.assertEqual(result['status'],'failed')


if __name__ == '__main__': unittest.main()
