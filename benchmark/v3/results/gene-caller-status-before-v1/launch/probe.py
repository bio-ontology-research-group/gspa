import argparse,csv,hashlib,json,os,shlex,shutil,subprocess,zipfile
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--jar',type=Path,required=True);p.add_argument('--observed',type=Path,required=True);p.add_argument('--wheels',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--groovy',type=Path,required=True);a=p.parse_args()
assert hashlib.sha256(a.jar.read_bytes()).hexdigest()=='069263f7b876222815928a637d02106b9c3b31a756edd1cf692432acc2e52303'
a.out.mkdir(parents=True,exist_ok=False)
libs=a.out/'libraries';libs.mkdir()
for f in a.wheels.glob('*.whl'):
 with zipfile.ZipFile(f) as z:z.extractall(libs)
wrapper=a.out/'pyrodigal'
wrapper.write_text('#!/bin/bash\nexport PYTHONPATH='+shlex.quote(str(libs))+'\nexec /usr/bin/python3 -m pyrodigal "$@"\n');wrapper.chmod(0o755)
commands=[]
def run(name,input,exe,mode,extra=()):
 cmd=['java','-cp',str(a.jar),'groovy.ui.GroovyMain',str(a.groovy),str(input),str(exe),mode,str(a.out/name),*map(str,extra)]
 with (a.out/(name+'.log')).open('wb') as log:proc=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,timeout=180)
 commands.append(dict(name=name,argv=cmd,exit_code=proc.returncode));(a.out/'commands.json').write_text(json.dumps(commands,indent=2)+'\n')
 assert proc.returncode==0,(name,proc.returncode)
 return json.loads((a.out/name/'result.json').read_text())
summary={}
for mode in ('single','meta'):
 summary['native-'+mode]=run('native-'+mode,a.observed/'inputs/genome.fna',wrapper,mode)
for mode in ('single','meta'):
 summary['short-'+mode]=run('short-'+mode,a.observed/'inputs/short.fna',wrapper,mode)
# Controlled corruption of captured native output: not a native caller observation.
records=['>'+r for r in (a.observed/'genome-single.faa').read_text().split('>') if r]
for name,selected in [('reordered',list(reversed(records))),('missing-protein',records[1:])]:
 faa=a.out/(name+'.faa');faa.write_text(''.join(selected))
 replay=a.out/(name+'-replay')
 replay.write_text('#!/usr/bin/python3\nimport sys,shutil\na=sys.argv\nshutil.copyfile('+repr(str(a.observed/'genome-single.gff'))+',a[a.index("-o")+1])\nshutil.copyfile('+repr(str(faa))+',a[a.index("-a")+1])\n');replay.chmod(0o755)
 extra=(a.observed/'genome-single.gff',faa) if name=='reordered' else ()
 summary[name]=run(name,a.observed/'inputs/genome.fna',replay,'single',extra)
def rows(path):
 with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def coordinates(r):return (r['contig_id'],r['start'],r['end'],r['strand'])
# Independent reference uses explicit native FASTA header coordinates, not row order.
native={}
for record in records:
 lines=record.splitlines();parts=lines[0][1:].split(' # ');ident=parts[0]
 key=(ident.rsplit('_',1)[0],parts[1],parts[2],'+' if parts[3]=='1' else '-')
 assert key not in native
 native[key]=''.join(lines[1:])
for name in ('native-single','reordered','missing-protein'):
 got=rows(a.out/name/'genes.tsv');summary[name]['sequence_coordinate_mismatches']=sum(native.get(coordinates(r))!=r['sequence'] for r in got)
 strict=a.out/name/'strict-genes.tsv'
 if strict.exists():summary[name]['strict_reader_mismatches']=sum(native.get(coordinates(r))!=r['sequence'] for r in rows(strict))
assert summary['native-single']['protein_count']==4359
assert summary['native-meta']['protein_count']==4291
assert summary['native-single']['sequence_coordinate_mismatches']==0
assert summary['short-meta']['protein_count']==0
assert summary['short-meta']['statuses'][0]['outcome']=='FAILED'
assert summary['short-single']['statuses'][0]['outcome']=='FAILED'
assert summary['reordered']['sequence_coordinate_mismatches']>0
assert summary['reordered']['strict_reader_mismatches']==0
assert summary['missing-protein']['feature_count']==4359
assert summary['missing-protein']['protein_count']==4358
assert summary['missing-protein']['statuses'][0]['outcome']=='SUCCESS'
(a.out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
