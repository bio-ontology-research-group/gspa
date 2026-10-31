import argparse, hashlib, json, os, platform, shutil, subprocess, sys, time, zipfile
from pathlib import Path
p=argparse.ArgumentParser()
p.add_argument('--out', type=Path, required=True)
p.add_argument('--work', type=Path, required=True)
a=p.parse_args(); out=a.out; work=a.work
out.mkdir(parents=True, exist_ok=False); work.mkdir(parents=True, exist_ok=True)
base=Path('/data/hohndor/gspa-v3/operon-identity-v1')
expected={
'pyrodigal':'b18910b02f79f5dc59568e2873f4ab9f0e528c8fe01f253b0df562db3feee4c4',
'archspec':'22463371fe8e9fa92fa144bf320cb445373a99f2ca6a5914c890e0c0b7816f8f'}
bindings={}
libs=work/'libs';libs.mkdir()
for prefix,digest in expected.items():
 files=list((base/'wheels').glob(prefix+'-*.whl'));assert len(files)==1,files
 f=files[0];actual=hashlib.sha256(f.read_bytes()).hexdigest();assert actual==digest
 bindings[str(f)]=actual
 with zipfile.ZipFile(f) as z:z.extractall(libs)
assembly=base/'observed/10484/genome.fna'
assert hashlib.sha256(assembly.read_bytes()).hexdigest()=='76a0087d8d90585696bd21d4c01b0e1f6699c72cce2a9a4682898c0255154ed5'
bindings[str(assembly)]=hashlib.sha256(assembly.read_bytes()).hexdigest()
inputs=out/'inputs';inputs.mkdir();shutil.copyfile(assembly,inputs/'genome.fna')
(inputs/'all-n.fna').write_text('>no_information\n'+'N'*10000+'\n')
(inputs/'short.fna').write_text('>short\nATGAAATAA\n')
(inputs/'empty.fna').write_text('')
env=dict(os.environ,PYTHONPATH=str(libs))
commands=[]
def run(name,cmd):
 started=time.monotonic()
 with (out/(name+'.stdout')).open('wb') as stdout,(out/(name+'.stderr')).open('wb') as stderr:
  proc=subprocess.run(cmd,stdout=stdout,stderr=stderr,env=env,timeout=180)
 rec=dict(name=name,argv=cmd,exit_code=proc.returncode,elapsed_seconds=time.monotonic()-started)
 commands.append(rec)
 (out/'commands.json').write_text(json.dumps(commands,indent=2)+'\n')
 return rec
binary=Path('/data/hohndor/gspa-v3/gene-caller-status-v1/prodigal-2.6.3')
assert hashlib.sha256(binary.read_bytes()).hexdigest()=='80a3d7acb08ba9f63dde6b6ef4e95ec356128ac355f36960ea722f381150f246'
bindings[str(binary)]='80a3d7acb08ba9f63dde6b6ef4e95ec356128ac355f36960ea722f381150f246'
run('prodigal-version',[str(binary),'-v'])
prodigal=str(binary)
cases=[('genome-single','genome','single'),('genome-meta','genome','meta'),('all-n-meta','all-n','meta'),('all-n-single','all-n','single'),('short-single','short','single'),('short-meta','short','meta'),('empty-meta','empty','meta')]
observations=[]
for name,source,mode in cases:
 gff=out/(name+'.gff');faa=out/(name+'.faa')
 rec=run(name,[str(binary),'-i',str(inputs/(source+'.fna')),'-o',str(gff),'-a',str(faa),'-f','gff','-g','11','-p',mode])
 rec=dict(rec)
 for label,f in [('gff',gff),('faa',faa)]:
  rec[label]=dict(exists=f.exists())
  if f.exists():
   data=f.read_bytes();rec[label].update(bytes=len(data),sha256=hashlib.sha256(data).hexdigest())
   lines=data.decode().splitlines()
   if label=='gff':rec[label].update(cds_rows=sum(1 for l in lines if not l.startswith('#') and len(l.split('\t'))==9 and l.split('\t')[2]=='CDS'),comment_lines=sum(l.startswith('#') for l in lines))
   else:rec[label]['protein_headers']=sum(l.startswith('>') for l in lines)
 observations.append(rec)
(out/'observations.json').write_text(json.dumps(dict(python=sys.version,host=platform.node(),native_prodigal_path=prodigal,bindings=bindings,cases=observations),indent=2)+'\n')
for path,digest in bindings.items():assert hashlib.sha256(Path(path).read_bytes()).hexdigest()==digest
print(json.dumps(observations,indent=2))
