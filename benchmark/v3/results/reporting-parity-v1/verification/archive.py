from pathlib import Path
import gzip,hashlib,json,shutil,subprocess
repo=Path('/home/leechuck/Public/software/gspa');scratch=Path(__file__).parent
out=repo/'benchmark/v3/results/reporting-parity-v1'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for block in iter(lambda:f.read(1048576),b''):h.update(block)
 return h.hexdigest()
r=json.loads((scratch/'replay-v4/verification.json').read_text())
assert r['status']=='verified' and len(r['cases'])==7
assert {c['name']:c['input_proteins'] for c in r['cases']}=={'mixed-logistic':2994,'mixed-gbdt':2994,'unsupported-logistic':2661,'unsupported-gbdt':2661,'absence-unknown':1,'absence-budget':1,'absence-baseline':1}
for path,digest in r['sources'].items():assert sha(Path(path))==digest,path
assert 'Ran 46 tests' in (scratch/'service-gate-final2.log').read_text()
assert 'Ran 73 tests' in (scratch/'browser-tests.log').read_text()
for name in ['service-gate-final2.log','browser-tests.log']:assert (scratch/name).read_text().rstrip().endswith('OK')
out.mkdir(exist_ok=False);index={}
def add(p,key):
 assert key not in index
 stored=key+'.gz' if p.stat().st_size>524288 else key
 target=out/stored;target.parent.mkdir(parents=True,exist_ok=True)
 if key==stored:shutil.copyfile(p,target)
 else:
  with p.open('rb') as src,target.open('wb') as dst,gzip.GzipFile(filename='',mode='wb',mtime=0,compresslevel=1,fileobj=dst) as z:shutil.copyfileobj(src,z,1048576)
 index[key]=dict(source=str(p),stored=stored,bytes=p.stat().st_size,sha256=sha(p))
for p in sorted(scratch.rglob('*')):
 if p.is_file() and p.name!='archive.log' and '__pycache__' not in p.parts:add(p,'verification/'+str(p.relative_to(scratch)))
for name in ['service/ledger_summary.py','service/web/index.html','service/portal.py','service/jobs.py','service/test_ledger_summary.py','service/test_portal.py','service/test_integration_options.py','service/requirements.txt','gspa-cli/src/main/resources/visualize/make_viz.py','gspa-cli/src/main/groovy/gspa/cli/VisualizeCommand.groovy','scripts/tests/test_visualize_browser.py','scripts/tests/test_visualize_learned.py','scripts/tests/test_visualize_ledger.py','docs/benchmarks/v3-reporting-parity.md','tasks/todo.md','tasks/v3-acceptance.md']:
 add(repo/name,'source/'+name)
for label,p in [('go-basic.obo',Path.home()/'.cache/gspa-v3-data/training-cohort/go-basic.obo'),('ec2go.txt',repo/'benchmark/whole_genome_cafa/data/reference/ec2go.txt')]:add(p,'references/'+label)
(out/'storage-index.json').write_text(json.dumps(index,indent=2,sort_keys=True)+'\n')
for logical,record in index.items():
 p=out/record['stored'];h=hashlib.sha256();n=0
 opener=gzip.open if logical!=record['stored'] else open
 with opener(p,'rb') as f:
  for block in iter(lambda:f.read(1048576),b''):h.update(block);n+=len(block)
 assert h.hexdigest()==record['sha256'] and n==record['bytes'],logical
summary=dict(status='verified',payloads=len(index),original_bytes=sum(v['bytes'] for v in index.values()),stored_bytes=sum((out/v['stored']).stat().st_size for v in index.values()),final_replay='verification/replay-v4',production_replays=7,service_tests=46,visualization_tests=73,whole_item_accepted=False)
(out/'archive-verification.json').write_text(json.dumps(summary,indent=2)+'\n')
(out/'README.md').write_text('''# Reporting parity: evidence and limitations

Three reporting defects are repaired: mixed scoring appeared unknown in the
portal; complete model abstention reports were marked untrusted; nested stage
assumptions and outcomes rendered as `[object Object]` in the results browser.

**Final evidence is `verification/replay-v4/`.** Seven stored production-output
sets pass through the packaged visualization CLI, real FastAPI endpoints and
Chromium. Sources, scores, candidate identities and all input proteins remain
intact. Two models each retain 26,328 mixed-scoring candidates (7,630 unsupported)
and 3,469 wholly unsupported candidates. Three synthetic absence-mechanism
outputs retain their uncertainty and budget/no-work states.

`verification/replay/` and `replay-final/` are intermediate harness runs: their
protein-only adapter omitted input proteins without annotations. `replay-v3/`
stopped on an incorrect assertion key for the browser's compact protein schema.
Do not use their protein totals as final evidence. The final adapter includes all
2,661 input proteins and reports 802 annotated proteins. Original production
prediction files were unchanged throughout.

Final gates pass 46 service tests and 73 visualization tests, including live
Chromium. Failing-before logs, the missing-FastAPI environment failure, package
versions, source snapshots, exact CLI commands, old/new portal summaries and
screenshots are retained. `verification/jar-verification.json` proves that the
new JAR changes only `visualize/make_viz.py`; all other ZIP entries are identical.
The immutable JAR is cached locally by its full SHA-256
`b410dadec00ead589f1f2fd0280ced8eea460755c202c41784f4de3df8a012d5`.

The source report gives reproduction commands and scope. Replay helpers retain
original workspace/cache paths in their receipts; restore or remap those inputs
when running elsewhere. `references/` holds the rendering ontology and EC map.
Large files use deterministic gzip. The storage index gives original sizes and
uncompressed hashes; run `sha256sum -c SHA256SUMS` to verify stored files.

This is reporting correctness, not biological validation. Batch provenance was
unavailable and remains empty. Protein-only inputs lack genomic coordinates.
No public service was deployed. No whole audit item is signed off.
''')
files=sorted(p for p in out.rglob('*') if p.is_file() and p.name!='SHA256SUMS')
(out/'SHA256SUMS').write_text(''.join(sha(p)+'  '+str(p.relative_to(out))+'\n' for p in files))
summary['checksummed_files']=len(files)
print(json.dumps(summary))
