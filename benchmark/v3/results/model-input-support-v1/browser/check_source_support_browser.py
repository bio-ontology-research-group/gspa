import csv, hashlib, json, os, re, subprocess, sys
from pathlib import Path
from urllib.parse import unquote
from playwright.sync_api import sync_playwright
root=Path('/home/leechuck/Public/software/gspa')
cache=Path('/home/leechuck/.cache/gspa-v3-data')
out=cache/'model-input-support-browser-v2-final'
out.mkdir(exist_ok=False)
gid='GCA_000007005.1'
source=cache/'model-input-support-cli-v2-final/homology/logistic'
genome=cache/'training-cohort/genomes'/gid
adapter=out/'prokka_out'; adapter.mkdir()
inputs=[Path(__file__),root/'gspa-cli/src/main/resources/visualize/make_viz.py',genome/'genomic.gff',genome/'genomic.fna',cache/'training-cohort/go-basic.obo',root/'benchmark/whole_genome_cafa/data/reference/ec2go.txt']
rows=[]; gff=['##gff-version 3']; seen=set()
for line in (genome/'genomic.gff').read_text().splitlines():
    if line.startswith('#'): continue
    cols=line.split('\t')
    if len(cols)<9 or cols[2]!='CDS': continue
    attrs=dict(p.split('=',1) for p in cols[8].split(';') if '=' in p)
    pid=attrs.get('protein_id')
    if not pid: continue
    assert pid not in seen, pid
    seen.add(pid)
    rows.append([pid,'CDS',int(cols[4])-int(cols[3])+1,unquote(attrs.get('product','')),attrs.get('locus_tag','')])
    attrs['locus_tag']=pid
    cols[8]=';'.join(k+'='+v for k,v in attrs.items());gff.append('\t'.join(cols))
with (adapter/(gid+'.tsv')).open('w') as fh:
    w=csv.writer(fh,delimiter='\t');w.writerow(['locus_tag','ftype','length_bp','product','original_locus_tag']);w.writerows(rows)
(adapter/(gid+'.gff')).write_text('\n'.join(gff)+'\n')
run=out/'gspa_out';run.mkdir()
for suffix in ('','.stages.json','.learned.json'):
    p=source/(gid+'.integrated.tsv'+suffix); inputs.append(p)
    (run/('integrated.tsv'+suffix)).symlink_to(p)
# Batch output has no provenance sidecar. Keep unavailable provenance empty;
# do not synthesize supporting evidence. Learned report retains actual sources.
(run/'provenance.json').write_text('[]\n')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
identity={str(p):sha(p) for p in inputs}
html=out/'browser.html'
env={**os.environ,'GSPA_WORKDIR':str(out),'GSPA_GENOME_ID':gid,'GSPA_OUT':str(html),'GSPA_RUN_DIR':'gspa_out','GSPA_FASTA':str(genome/'genomic.fna'),'GSPA_GO_OBO':str(cache/'training-cohort/go-basic.obo'),'GSPA_EC2GO':str(root/'benchmark/whole_genome_cafa/data/reference/ec2go.txt')}
with (out/'render.log').open('w') as log:
    subprocess.run([sys.executable,str(inputs[1])],env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
data=json.loads(re.search(r'const DATA = (.*);',html.read_text()).group(1))
learned=data['ledger']['learned']; assert learned['status']=='available'
assert len(learned['candidates'])==26328
assert sum(c['ss']=='unsupported' for c in learned['candidates'])==7630
assert data['ledger']['scoring']['mode']=='mixed'
assert not learned['n_unjoined_candidates']
errors=[]
with sync_playwright() as pw:
    browser=pw.chromium.launch(headless=True); page=browser.new_page()
    page.on('pageerror',lambda e: errors.append(str(e)))
    page.goto(html.as_uri());page.locator('[data-tab="ledger"]').click()
    assert 'Mixed' in page.locator('#lg-scoring').inner_text()
    assert '7630 unsupported' in page.locator('#ln-count').inner_text()
    candidate=next(c for c in learned['candidates'] if c['ss']=='unsupported')
    page.locator('#ln-search').fill(candidate['p']); page.locator('#ln-role').select_option(candidate['role'])
    matching=page.locator('#ln-rows tr').filter(has_text=candidate['f']).first
    text=matching.inner_text(); assert 'unsupported model input' in text and candidate['ur'] in text
    assert 'below' not in text.lower()
    matching.locator('[data-lg-protein]').click()
    assert 'unsupported model input' in page.locator('#d-ledger').inner_text()
    page.screenshot(path=str(out/'unsupported-detail.png'),full_page=False)
    browser.close()
assert not errors,errors
assert identity=={str(p):sha(p) for p in inputs}
report={'status':'verified','candidate_rows':26328,'unsupported_rows':7630,'unjoined_candidates':0,'page_errors':errors,'inputs':identity,'artifacts':{p.name:sha(p) for p in [html,out/'unsupported-detail.png']},'scope':'Real mixed production report. NCBI CDS metadata adapted to the browser input schema. Missing batch provenance remains empty; no provenance completeness claim.'}
(out/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k not in ('inputs','artifacts')}))
