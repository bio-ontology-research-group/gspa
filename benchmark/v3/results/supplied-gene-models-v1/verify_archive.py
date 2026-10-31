#!/usr/bin/env python3
"""Independent portable verification of retained native joins and production outputs."""
import csv,gzip,hashlib,io,json,re,xml.etree.ElementTree as ET
from pathlib import Path
root=Path(__file__).resolve().parent
manifest=json.loads((root/'retained.json').read_text())
raw={}
for key,entry in manifest.items():
 data=gzip.decompress((root/entry['file']).read_bytes())
 assert len(data)==entry['bytes'] and hashlib.sha256(data).hexdigest()==entry['sha256'],key
 raw[key]=data

def txt(key):return raw[key].decode()
def js(key):return json.loads(txt(key))
def table(key):return list(csv.DictReader(io.StringIO(txt(key)),delimiter='\t'))
def sequences(key):
 result={}
 for block in txt(key).split('>'):
  if not block:continue
  lines=block.splitlines();pid=lines[0].split()[0]
  assert pid not in result
  result[pid]=(lines[0],''.join(lines[1:]))
 return result

def native(caller,mode):
 prefix='numeric' if caller=='numeric' else 'genome'
 base='native/'+caller+'/'+prefix+'-'+mode
 records=sequences(base+'.faa');ref={}
 for pid,(header,seq) in records.items():
  m=re.fullmatch(r'\S+\s+#\s+(\d+)\s+#\s+(\d+)\s+#\s+(-?1)\s+#\s+.*',header);assert m
  ref[pid]=(pid.rsplit('_',1)[0],int(m[1]),int(m[2]),'+' if m[3]=='1' else '-',seq)
 gff=[line.split('\t') for line in txt(base+'.gff').splitlines() if line and not line.startswith('#')]
 coords=[(r[0],int(r[3]),int(r[4]),r[6]) for r in gff if r[2]=='CDS']
 assert len(coords)==len(ref)==len(set(coords)) and set(coords)=={r[:4] for r in ref.values()}
 return ref,prefix+'-'+mode

def groups(ref,excluded):
 by_contig={}
 for pid,row in ref.items():
  if pid not in excluded:by_contig.setdefault(row[0],[]).append(pid)
 answer=set()
 for proteins in by_contig.values():
  ordered=sorted(proteins,key=lambda p:ref[p][1]);breaks=[0]
  for i in range(1,len(ordered)):
   a,b=(ref[p] for p in ordered[i-1:i+1])
   if a[3]!=b[3] or b[1]>a[2]+301:breaks.append(i)
  breaks.append(len(ordered))
  for start,end in zip(breaks,breaks[1:]):
   if end-start>=2:answer.add(tuple(ordered[start:end]))
 return answer

source=js('preparation/source-v3.json')
assert len(source)==845
for path,digest in source.items():assert manifest['source/'+path]['sha256']==digest,path
for job in [10523,10524,10526,10527]:assert txt(f'jobs/{job}/exit-code.txt').strip()=='0'
counts={}
for label,prefix in [('before','jobs/10526/before-tests/'),('focused','jobs/10526/focused-tests/'),('full','jobs/10526/final-tests/')]:
 suites=[ET.fromstring(data) for key,data in raw.items() if key.startswith(prefix) and key.endswith('.xml') and '/TEST-' in key]
 counts[label]={k:sum(int(s.attrib[k]) for s in suites) for k in ['tests','failures','errors','skipped']}
assert counts['before']['failures']>0
assert counts['focused']['failures']==counts['focused']['errors']==counts['full']['failures']==counts['full']['errors']==0
assert counts['full']['tests']>1552 and counts['full']['skipped']==6
for job,total in [(10523,4),(10527,10)]:
 runs=js(f'jobs/{job}/probe/runs.json');assert len(runs)==total
 assert all(r['exit_code']==0 and not r.get('timed_out') for r in runs)

summary={}
for name,record in js('jobs/10527/probe/summary.json').items():
 if name.startswith('numeric-'):caller='numeric';mode=name.split('-')[1]
 elif name.startswith('pyrodigal-'):caller='pyrodigal';mode=name.split('-')[1]
 else:caller='prodigal';mode='meta' if name=='prodigal-meta' else 'single'
 ref,prefix=native(caller,mode)
 excluded={next(iter(ref))} if name=='prodigal-header-conflict' else set()
 base='jobs/10527/probe/'+name+'/'
 rows=table(base+prefix+'_proteins.tsv');emitted=sequences(base+prefix+'.faa')
 assert len(rows)==len(ref) and {r['protein_id'] for r in rows}==set(ref)==set(emitted)
 for row in rows:
  pid=row['protein_id'];expected=ref[pid]
  assert emitted[pid][1]==expected[4]
  observed=(row['contig'],int(row['start']),int(row['end']),row['strand'])
  if pid in excluded:assert observed[1:3]==(0,0)
  else:assert observed==expected[:4],(name,pid)
 wanted=groups(ref,excluded)
 got=[tuple(r.split('\t')) for r in txt(base+'operons_for_integrate.tsv').splitlines()]
 assert len(got)==len(set(got)) and set(got)==wanted,name
 members=table(base+'protein_to_operon.tsv');assert len(members)==sum(map(len,wanted))
 assert {r['protein_id'] for r in members}=={p for group in wanted for p in group}
 actual=dict(proteins=len(ref),placed=len(ref)-len(excluded),unplaced=len(excluded),operons=len(wanted),members=len(members))
 assert all(record[k]==v for k,v in actual.items())
 if name.startswith('pyrodigal-'):
  for file in [prefix+s for s in ['.faa','.gff3','_proteins.tsv','_annotations.tsv']]+['operons.tsv','protein_to_operon.tsv','operons_for_integrate.tsv']:
   assert raw[base+file]==raw['jobs/10523/probe/'+name+'/'+file]
 summary[name]=actual

for caller in ['pyrodigal','prodigal']:
 for mode in ['single','meta']:
  name=caller+'-'+mode;ref,prefix=native(caller,mode)
  rows=table('jobs/10523/probe/'+name+'/'+prefix+'_proteins.tsv')
  assert len(rows)==len(ref) and {r['protein_id'] for r in rows}==set(ref)
  for r in rows:
   if caller=='prodigal':assert r['contig']=='unplaced' and int(r['start'])==int(r['end'])==0
   else:assert (r['contig'],int(r['start']),int(r['end']),r['strand'])==ref[r['protein_id']][:4]

available={entry['sha256'] for entry in manifest.values()}
omitted={};retained=0
for key in ['jobs/10523/probe/inputs.sha256.json','jobs/10527/probe/bindings.json']:
 for path,digest in js(key).items():
  if digest in available:retained+=1
  else:omitted[key+"::"+path]=digest
assert all(p.endswith(('gspa.jar','/java')) for p in omitted),omitted
assert set(omitted.values())=={txt('jobs/10526/jar.sha256').split()[0],
 '1ae842e5905eecfff4ecbcb07bfa21d320ff9de0797273e4aeafd2abaccbdb5d',
 '6698f6f10143ed8463b06062281c727152d2e0ae2a3569b1716fbb2c6cedf0ba'},omitted
report=dict(verified=True,retained_files=len(raw),source_files=len(source),tests=counts,cli_runs=14,
 after=summary,retained_runtime_bindings=retained,omitted_runtime_bindings=omitted,
 limitation='Coordinate and neighborhood preservation, not biological accuracy or whole-item acceptance.')
print(json.dumps(report,indent=2))
