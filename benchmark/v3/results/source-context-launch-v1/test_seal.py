import importlib.util,json,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
spec=importlib.util.spec_from_file_location('seal',Path(__file__).with_name('seal.py'));S=importlib.util.module_from_spec(spec);spec.loader.exec_module(S)
SHA='37aeee2b565951dd3a4e049a23aacd2e8155f29bd8b6216fed7006e4273965e6'
class SealTests(unittest.TestCase):
 def setUp(self):
  self.tmp=tempfile.TemporaryDirectory(dir=Path(__file__).parent);self.addCleanup(self.tmp.cleanup);self.r=Path(self.tmp.name)
  self.probe=self.r/'probe/results/verification.json';self.clean=self.r/'clean/audit.json'
  for p in [self.probe,self.clean]:p.parent.mkdir(parents=True)
  self.p=dict(schema='gspa.source-context-probe/2',status='verified',ec_identity=dict(shared_source_candidates=511,verified_both_source_candidates=511,combined_candidates=29512,raw_claim_files_unchanged=True),inputs={str(S.JAR):SHA})
  self.probe.write_text(json.dumps(self.p));self.clean.write_text(json.dumps(dict(schema='gspa.clean-cohort-audit/2',genomes=18,proteins=49463,native_proteins=49462,unsupported_proteins=1)))
  self.exit=self.probe.parent.parent/'exit-code.txt';self.exit.write_text('0\n');(self.clean.parent/'exit-code.txt').write_text('0\n')
  self.source=self.r/'source';self.source.write_text('frozen source');self.static=self.r/'static.sha256'
  self.static.write_text(SHA+'  '+str(S.JAR)+'\n'+S.digest(self.source)+'  '+str(self.source)+'\n')
  digest=S.digest;self.patch=patch.object(S,'digest',lambda p:SHA if p==S.JAR else digest(p));self.patch.start();self.addCleanup(self.patch.stop)
 def run_seal(self):return S.seal(self.r,self.probe,self.clean)
 def test_complete_and_exclusive(self):
  records=self.run_seal();self.assertIn(str(self.probe),records);self.assertIn(str(self.clean),records)
  before=(self.r/'packages.sha256').read_bytes()
  with self.assertRaises(FileExistsError):self.run_seal()
  self.assertEqual(before,(self.r/'packages.sha256').read_bytes())
 def test_failed_job(self):
  self.exit.write_text('1\n')
  with self.assertRaisesRegex(ValueError,'did not exit'):self.run_seal()
 def test_missing_job_exit(self):
  self.exit.unlink()
  with self.assertRaises(FileNotFoundError):self.run_seal()
 def test_old_split_identity_rejected(self):
  self.p['ec_identity']['combined_candidates']=30023;self.probe.write_text(json.dumps(self.p))
  with self.assertRaisesRegex(ValueError,'EC identity'):self.run_seal()
 def test_changed_frozen_source(self):
  self.source.write_text('changed')
  with self.assertRaisesRegex(ValueError,'static input'):self.run_seal()
 def test_duplicate_static_path(self):
  self.static.write_text(self.static.read_text()*2)
  with self.assertRaisesRegex(ValueError,'duplicate'):self.run_seal()
 def test_smaller_clean_population(self):
  c=json.loads(self.clean.read_text());c['proteins']-=1;self.clean.write_text(json.dumps(c))
  with self.assertRaisesRegex(ValueError,'population'):self.run_seal()
if __name__=='__main__':unittest.main(verbosity=2)
