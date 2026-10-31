from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.request import Request,urlopen
from urllib.error import HTTPError
from datetime import datetime,timezone
import json,hashlib
root=Path('/home/leechuck/tmp/gspa-secondary-audit-v1/sources');root.mkdir(exist_ok=True)
sources=[('tags.json','https://api.github.com/repos/SysBioChalmers/Sco-GEM/tags?per_page=100'),('article.xml','https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7501462/fullTextXML'),('metadata.json','https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:32942174%20AND%20SRC:MED&format=json')]
def get(item):
 name,url=item
 try:
  with urlopen(Request(url,headers={'User-Agent':'GSPA-research-audit/1'}),timeout=60) as r:data=r.read();status=r.status;ctype=r.headers.get('Content-Type');final=r.url
 except HTTPError as e:data=e.read();status=e.code;ctype=e.headers.get('Content-Type');final=e.url
 path=root/name;path.write_bytes(data)
 return dict(file=name,url=url,status=status,content_type=ctype,final_url=final,bytes=len(data),sha256=hashlib.sha256(data).hexdigest(),retrieved_at=datetime.now(timezone.utc).isoformat())
receipts=list(ThreadPoolExecutor(3).map(get,sources));(root/'initial-receipts.json').write_text(json.dumps(receipts,indent=2)+'\n')
print(json.dumps(receipts,indent=2))
for x in json.loads((root/'tags.json').read_text()):print(x['name'],x['commit']['sha'])
