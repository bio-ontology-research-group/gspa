"""Actual stored production artifacts -> packaged CLI browser -> portal HTTP/DOM.

No new biological inference, no outcome labels, and no manufactured provenance.
"""
from pathlib import Path
import csv, hashlib, importlib.util, json, os, re, shutil, subprocess, sys
from unittest.mock import patch
from fastapi.testclient import TestClient
from playwright.sync_api import sync_playwright

REPO=Path('/home/leechuck/Public/software/gspa');sys.path.insert(0,str(REPO))
from service import jobs,portal
CACHE=Path.home()/'.cache/gspa-v3-data';ROOT=Path(__file__).parent
OUT=ROOT/'replay-v3';OUT.mkdir(exist_ok=False)
JAR=REPO/'gspa-cli/build/libs/gspa-3.0.0.jar'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
sources={str(p):sha(p) for p in [Path(__file__),JAR,REPO/'service/ledger_summary.py',REPO/'service/web/index.html',REPO/'gspa-cli/src/main/resources/visualize/make_viz.py']}
spec=importlib.util.spec_from_file_location('old_summary',ROOT/'ledger_summary-before.py');old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
cases=[]
for model in ['logistic','gbdt']:
 cases.append((f'mixed-{model}',CACHE/f'model-input-support-cli-v2-final/homology/{model}/GCA_000007005.1.integrated.tsv','mixed',26328,7630,'genome'))
 cases.append((f'unsupported-{model}',REPO/f'benchmark/v3/results/gapsmith-observed-v2/final-cli-replay/{model}.tsv','legacy_heuristic',3469,3469,'proteins'))
for name in ['unknown','budget','baseline']:
 cases.append((f'absence-{name}',REPO/f'benchmark/v3/results/absence-comparison-v1/cli/{name}.tsv','legacy_heuristic',1,None,'mechanism'))
reports=[]
client=TestClient(portal.app)
with patch.object(jobs,'ROOT',OUT),sync_playwright() as pw:
 browser=pw.chromium.launch(headless=True)
 for i,(name,source,mode,n,unsupported,metadata) in enumerate(cases,1):
  jid=f'{i:032x}';folder=OUT/jid;run=folder/'gspa_out';run.mkdir(parents=True)
  for suffix in ['','.stages.json','.learned.json']:
   src=Path(str(source)+suffix)
   if src.exists():
    sources[str(src)]=sha(src);shutil.copyfile(src,run/('integrated.tsv'+suffix))
  rows=list(csv.DictReader(source.open(),delimiter='\t'));assert len(rows)==n
  pids={r['protein_id'] for r in rows}
  adapter=folder/'prokka_out';adapter.mkdir()
  gid='report'
  if metadata=='genome':
   for ext in ['tsv','gff']:
    src=CACHE/f'model-input-support-browser-v2-final/prokka_out/GCA_000007005.1.{ext}'
    sources[str(src)]=sha(src);shutil.copyfile(src,adapter/f'{gid}.{ext}')
  else:
   meta={}
   if metadata=='proteins':
    src=REPO/'benchmark/v3/results/gapsmith-observed-v2/final-cli-replay/integrated/protein_proteins.tsv'
    sources[str(src)]=sha(src);meta={r['protein_id']:r for r in csv.DictReader(src.open(),delimiter='\t')}
    assert pids<=set(meta)
   with (adapter/f'{gid}.tsv').open('w') as f:
    w=csv.writer(f,delimiter='\t');w.writerow(['locus_tag','ftype','length_bp','product'])
    for pid in sorted(meta if metadata=='proteins' else pids):
     m=meta.get(pid,{});w.writerow([pid,'CDS',3*int(m.get('length_aa',0)),m.get('product','')])
   # Protein-only and synthetic mechanism inputs have no genomic coordinates.
   (adapter/f'{gid}.gff').write_text('##gff-version 3\n')
  (run/'provenance.json').write_text('[]\n')
  go=CACHE/'training-cohort/go-basic.obo';ec=REPO/'benchmark/whole_genome_cafa/data/reference/ec2go.txt'
  for src in [go,ec]:sources[str(src)]=sha(src)
  cmd=['java','-jar',str(JAR),'visualize','--workdir',str(folder),'--run-dir','gspa_out','--genome-id',gid,'--out',str(run/'sample_browser.html'),'--go-obo',str(go),'--ec2go',str(ec),'--python',sys.executable]
  (folder/'command.json').write_text(json.dumps(cmd,indent=2)+'\n')
  with (folder/'render.log').open('w') as log:subprocess.run(cmd,cwd=REPO,env={**os.environ,'TMPDIR':str(Path.home()/'tmp'),'JAVA_TOOL_OPTIONS':'-Djava.io.tmpdir='+str(Path.home()/'tmp')},stdout=log,stderr=subprocess.STDOUT,check=True,timeout=120)
  data=json.loads(re.search(r'const DATA = (.*);',(run/'sample_browser.html').read_text()).group(1))
  ledger=data['ledger'];raw=json.loads((run/'integrated.tsv.stages.json').read_text())
  assert ledger['scoring']['mode']==mode
  if metadata=='proteins':assert {r['id'] for r in data['proteins']}==set(meta)
  assert ledger['reconciliation']['counts_match'] and ledger['reconciliation']['n_probability_mismatches']==0
  assert not ledger['reconciliation']['net_changes_without_stage_action']
  assert [(s['name'],s['status'],s['reason']) for s in ledger['stages']]==[(s['name'],s['status'],s.get('reason')) for s in raw['stages']]
  for stage,original in zip(ledger['stages'],raw['stages']):
   for key,value in original.items():
    if isinstance(value,(dict,list)) and key!='changes':assert stage['notes'][key]==value,(name,key)
  jobs.state(folder,id=jid,status='complete',filename=name,stage='Stored production replay',proteins=len(data['proteins']),annotated_proteins=len(pids),annotations=n,analysis_results={})
  response=client.get('/api/jobs/'+jid);assert response.status_code==200
  summary=response.json()['ledger'];before=old.summarize(run)
  assert summary['scoring']['mode']==ledger['scoring']['mode']
  assert summary['learned']['status']==ledger['learned']['status']
  if unsupported is not None:
   assert summary['learned']['status']=='available' and summary['learned']['hash_verified']
   assert summary['learned']['unsupported_only']==name.startswith('unsupported-')
   assert summary['learned']['counts']['unsupported']==unsupported
   learned=json.loads((run/'integrated.tsv.learned.json').read_text())
   assert len(ledger['learned']['candidates'])==len(learned['candidates'])==n
   assert ledger['learned']['n_unjoined_candidates']==0
   for shown,original in zip(ledger['learned']['candidates'],learned['candidates']):
    assert (shown['ck'],shown['k'],shown['pr'],shown['ss'],shown['src'])==(original['candidate_key'],original['function_key'],original['probability'],original.get('score_status','scored'),original['supporting_sources'])
   artifact=client.get('/api/jobs/'+jid+'/learned');assert artifact.content==(run/'integrated.tsv.learned.json').read_bytes()
  (folder/'summary-before.json').write_text(json.dumps(before,indent=2)+'\n')
  (folder/'summary-after.json').write_text(json.dumps(summary,indent=2)+'\n')
  page=browser.new_page();errors=[];page.on('pageerror',lambda e:errors.append(str(e)))
  page.goto((run/'sample_browser.html').as_uri());page.locator('[data-tab="ledger"]').click()
  assert '[object Object]' not in page.locator('#ledger-body table').first.inner_text()
  if unsupported is not None:
   assert f'{unsupported} unsupported' in page.locator('#ln-count').inner_text()
   assert len(ledger['learned']['candidates'])==n
  if metadata=='mechanism':
   stage=next(s for s in raw['stages'] if s['name']=='metabolic_experiments')
   if stage['status']=='executed':
    row=page.locator('#ledger-body table').first.locator('tr').filter(has_text='Metabolic experiments')
    row.locator('summary').click()
    expected=[v for k,v in stage.items() if k!='changes' and isinstance(v,(list,dict))]
    assert [json.loads(t) for t in row.locator('pre').all_text_contents()]==expected
    if name=='absence-unknown':assert 'UNKNOWN' in row.inner_text()
  page.screenshot(path=str(folder/'results-browser.png'),full_page=False)
  # Route the actual portal page's fetches through its real FastAPI endpoints.
  # This avoids starting a public server or a prediction worker.
  def route(request):
   path=request.request.url.removeprefix('http://gspa.test')
   if path=='/':request.fulfill(status=200,body=(REPO/'service/web/index.html').read_text(),content_type='text/html')
   else:
    r=client.get(path);request.fulfill(status=r.status_code,body=r.content,content_type=r.headers.get('content-type','application/octet-stream'))
  page.route('http://gspa.test/**',route);page.goto('http://gspa.test/#'+jid)
  page.wait_for_function("document.querySelector('#ledger').textContent.includes('Scoring:')")
  text=page.locator('#ledger').inner_text()
  assert 'not measurements of accuracy' in text
  assert f'{len(pids):,} of {len(data["proteins"]):,} proteins' in page.locator('#summary').inner_text()
  if unsupported is not None:
   assert page.locator('#learned-file').is_visible()
   assert summary['learned']['model_id'] in text
   assert f'{unsupported:,} unsupported' in text
   if mode=='mixed':assert 'Scoring: mixed' in text and 'Scoring: unknown' not in text
   else:assert 'model abstained' in text and 'baseline scores retained' in text
  else:assert not page.locator('#learned-file').is_visible()
  (folder/'portal-text.txt').write_text(text+'\n')
  page.screenshot(path=str(folder/'portal.png'),full_page=False)
  assert not errors,errors;page.close()
  result=dict(name=name,id=jid,scoring_mode=mode,rows=n,input_proteins=len(data['proteins']),annotated_proteins=len(pids),unsupported=unsupported,metadata=metadata,stage_statuses_match=True,export_reconciliation=True,page_errors=errors,before_report_status=before['learned']['status'],after_report_status=summary['learned']['status'])
  reports.append(result);print(json.dumps(result),flush=True)
 browser.close()
assert sources=={p:sha(Path(p)) for p in sources}
(OUT/'verification.json').write_text(json.dumps({'status':'verified','cases':reports,'sources':sources,'limits':['Stored production inference was not rerun.','The absence cases use a synthetic mechanism model.','Batch provenance is absent and remains empty; no provenance completeness claim.','Protein-only and mechanism cases have no genomic coordinates.','This is reporting correctness, not biological accuracy.']},indent=2)+'\n')
