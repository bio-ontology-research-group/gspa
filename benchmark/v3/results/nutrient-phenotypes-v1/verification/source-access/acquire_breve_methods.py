from pathlib import Path
import urllib.request,json,hashlib,concurrent.futures
b=Path.home()/'.cache/gspa-v3-data/breve-functional-genomics-v1'
commit='076ab75968b37591f5b484d86f7d3a5806cebef6'
root='https://api.bitbucket.org/2.0/repositories/kchuanglab/biffungenpaper/src/'+commit+'/'
def get(path):
 req=urllib.request.Request(root+path,headers={'User-Agent':'GSPA-reproducibility-audit/1'})
 with urllib.request.urlopen(req,timeout=40) as r:
  data=r.read(8*1024*1024+1)
  if len(data)>8*1024*1024:raise ValueError('size cap')
 return data
receipts=json.loads((b/'method-acquisition.json').read_text()) if (b/'method-acquisition.json').exists() else []
for path in ['cellshape_growthcurves/scripts','aromatic_lactic_acid/230117_Metabolomics/scripts']:
 try:
  data=get(path+'?pagelen=100');d=json.loads(data);(b/('tree-'+path.replace('/','--')+'.json')).write_bytes(data)
  print(path,[(x['path'],x.get('size')) for x in d.get('values',[])],flush=True)
  if d.get('next'):print('pagination incomplete',d['next'],flush=True)
  files=[x for x in d.get('values',[]) if x.get('type')=='commit_file' and x.get('size',0)<500000 and not x['path'].endswith('.DS_Store')]
  for x in files:
   target=b/'author-code'/x['path']
   if target.exists():
    assert target.stat().st_size==x['size']
    continue
   payload=get(x['path']);assert len(payload)==x['size']
   target=b/'author-code'/x['path'];target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(payload)
   receipts.append(dict(path=x['path'],url=root+x['path'],bytes=len(payload),sha256=hashlib.sha256(payload).hexdigest(),executed=False))
 except Exception as e:receipts.append(dict(path=path,error=str(e)))
path='essentialgaps/data/name.conversion.xlsx'
try:
 data=get(path);assert len(data)==59851
 target=b/'raw'/path;target.write_bytes(data)
 receipts.append(dict(path=path,url=root+path,bytes=len(data),sha256=hashlib.sha256(data).hexdigest()))
except Exception as e:receipts.append(dict(path=path,error=str(e)))
(b/'method-acquisition.json').write_text(json.dumps(receipts,indent=2)+'\n')
print('acquired',len(receipts),flush=True)
