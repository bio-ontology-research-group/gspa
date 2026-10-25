"""Validation and access controls independent of prediction assets."""
import io
import json
from pathlib import Path
import tempfile
import unittest
from types import SimpleNamespace
from unittest.mock import patch
from fastapi.testclient import TestClient
from service import jobs, portal


class PortalTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = patch.object(jobs, 'ROOT', Path(self.tmp.name)); self.root.start()
        self.disk = patch.object(jobs.shutil, 'disk_usage', return_value=SimpleNamespace(free=100 * 1024**3)); self.disk.start()
        self.client = TestClient(portal.app)
    def tearDown(self):
        self.disk.stop(); self.root.stop(); self.tmp.cleanup()
    def upload(self, text, kind='proteins', preset='eukaryote'):
        with patch.object(jobs, 'ready', return_value=True):
            return self.client.post('/api/jobs', data={'kind':kind, 'preset':preset}, files={'file':('sample.faa', text)})
    def test_persistent_job_and_private_address(self):
        r=self.upload('>p1\nMEEKL\n'); self.assertEqual(r.status_code,202,r.text)
        job=r.json(); self.assertEqual(len(job['id']),32)
        self.assertEqual(self.client.get('/api/jobs/'+job['id']).json()['status'],'queued')
        self.assertEqual(self.client.get('/api/jobs/'+job['id']+'/status.json').status_code,404)
        self.assertEqual(self.client.get('/api/jobs/not-a-job').status_code,404)
    def test_rejects_bad_fasta_without_orphans(self):
        for data in ('garbage','>\nMKT','>p\nMKT\n>p\nMKT','>bad<script>\nMKT','>p\n'):
            r=self.upload(data); self.assertEqual(r.status_code,400,r.text)
        self.assertEqual(list(jobs.ROOT.iterdir()),[])
    def test_type_and_size_limits(self):
        self.assertEqual(self.upload('>p\nMEEKL',kind='genome',preset='eukaryote').status_code,400)
        self.assertEqual(self.upload('>p\nMEEKL',kind='genome',preset='bacteria').status_code,400)
        with patch.object(jobs,'MAX_BYTES',4):
            self.assertEqual(self.upload('>p\nMEEKL').status_code,413)
    def test_queue_bound(self):
        for i in range(10):
            self.assertEqual(self.upload('>p\nMKT').status_code,202)
        self.assertEqual(self.upload('>p\nMKT').status_code,429)
    def test_failed_predictor_is_not_reported_as_zero_hits(self):
        log = jobs.ROOT / 'run.log'
        log.write_text('deepgo-plusplus-light failed: missing asset\nAnnotation complete: 0/10 proteins annotated')
        with self.assertRaises(RuntimeError):
            jobs.verify_prediction(log)
        log.write_text('deepgo-plusplus-light: 0 annotations assigned')
        jobs.verify_prediction(log)

    def test_full_profile_is_default_and_can_be_disabled(self):
        result=self.upload('>p\nMKT').json()
        self.assertTrue(all(result['analyses'].values()))
        with patch.object(jobs, 'ready', return_value=True):
            r=self.client.post('/api/jobs',data={'kind':'proteins','disabled':'names,amr'},files={'file':('x.faa','>p\nMKT')})
        self.assertEqual(r.status_code,202,r.text)
        self.assertFalse(r.json()['analyses']['names'])
        self.assertFalse(r.json()['analyses']['amr'])
        self.assertTrue(r.json()['analyses']['metabolism'])

    def test_gff3_requires_matching_sequence_and_retains_annotation(self):
        gff='##gff-version 3\nc\tRefSeq\tCDS\t1\t9\t.\t+\t0\tID=cds1;protein_id=P1;gene=abc;product=An enzyme\n'
        with patch.object(jobs, 'ready', return_value=True):
            self.assertEqual(self.client.post('/api/jobs',files={'file':('x.gff3',gff)}).status_code,400)
            r=self.client.post('/api/jobs',files={'file':('x.gff3',gff),'genome':('x.fna','>c\nATGAAATAA\n')})
            self.assertEqual(r.status_code,202,r.text)
            folder=jobs.ROOT/r.json()['id']
            self.assertIn('product=An enzyme',(folder/'sample.gff3').read_text())
            r=self.client.post('/api/jobs',files={'file':('x.gff3',gff+'##FASTA\n>c\nATGAAATAA\n')})
            self.assertEqual(r.status_code,202,r.text)
            r=self.client.post('/api/jobs',files={'file':('x.gff3',gff),'genome':('x.fna','>wrong\nATGAAATAA\n')})
            self.assertEqual(r.status_code,400)

    def test_unavailable_assets(self):
        with patch.object(jobs,'ready',return_value=False):
            self.assertEqual(self.client.get('/health').status_code,503)
            self.assertEqual(self.client.post('/api/jobs',files={'file':('x','>p\nMKT')}).status_code,503)

    def test_named_example_only_exposes_allowlisted_completed_artifacts(self):
        job = self.upload('>p\nMKT').json()
        folder = jobs.ROOT/job['id']
        (jobs.ROOT/'examples.json').write_text(json.dumps({'ecoli': job['id']}))
        (folder/'gspa_out').mkdir()
        (folder/'gspa_out/sample_browser.html').write_text('<html>Example</html>')
        self.assertEqual(self.client.get('/examples/ecoli').status_code,409)
        jobs.state(folder,status='complete',example=True)
        self.assertEqual(self.client.get('/examples/ecoli').status_code,200)
        self.assertEqual(self.client.get('/examples/ecoli/status.json').status_code,404)
        self.assertEqual(self.client.get('/examples/missing').status_code,404)

if __name__=='__main__': unittest.main()
