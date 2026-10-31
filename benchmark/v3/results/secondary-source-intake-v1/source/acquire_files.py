from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.request import Request,urlopen
from datetime import datetime,timezone
import json,hashlib
root=Path('/home/leechuck/tmp/gspa-secondary-audit-v1/sources');commit='63003f8b01a692dc2b73f2c9e0bbde871a3fe6f0'
tree=json.loads((root/'tree.json').read_text())['tree']
names=[r['path'] for r in tree if r['type']=='blob' and (r['path'].startswith(('ComplementaryData/growth/','ComplementaryData/media/')) or r['path'] in ['ModelFiles/xml/Sco-GEM.xml','README.md','LICENSE.md','ComplementaryData/experiments.yml','ComplementaryData/data/germicidin.csv','ComplementaryData/data/CDW.csv','ComplementaryScripts/rate_estimation.py','ComplementaryScripts/ecModel/reconstruction/preprocessModel.m','ComplementaryScripts/ecModel/simulation/fixFluxes.m','ComplementaryScripts/tests/growth_test.py'])]
def fetch(name):
 url='https://raw.githubusercontent.com/SysBioChalmers/Sco-GEM/'+commit+'/'+name
 with urlopen(Request(url,headers={'User-Agent':'GSPA-research-audit/1'}),timeout=60) as r:data=r.read();content_type=r.headers.get('Content-Type')
 p=root/'repository'/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(data)
 expected=next(e['sha'] for e in tree if e['path']==name)
 blob=hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest();assert blob==expected,name
 return dict(path=name,url=url,sha256=hashlib.sha256(data).hexdigest(),git_blob=blob,bytes=len(data),content_type=content_type,retrieved_at=datetime.now(timezone.utc).isoformat())
r=list(ThreadPoolExecutor(4).map(fetch,names));(root/'repository-receipt.json').write_text(json.dumps({'commit':commit,'tag':'v1.2.1','files':r},indent=2)+'\n');print('Verified',len(r),'Git blobs;',sum(x['bytes'] for x in r),'bytes')
