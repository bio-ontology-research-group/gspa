#!/usr/bin/env python3
"""Controlled CLEAN processes through real sidecar, CLI and Nextflow; no model inference."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--jar', type=Path, required=True)
p.add_argument('--nextflow', type=Path, required=True)
p.add_argument('--out', type=Path, required=True)
args = p.parse_args()
args.out.mkdir(parents=True, exist_ok=False)
out = args.out.resolve()
commands = []


def write(name, text):
    path = out / name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


def run(name, command, env, succeeds=True):
    result = subprocess.run([str(x) for x in command], cwd=out, env=env,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=600)
    write(name + '.log', result.stdout.decode(errors='replace'))
    commands.append(dict(name=name, command=[str(x) for x in command], exit_code=result.returncode))
    write('commands.json', json.dumps(commands, indent=2) + '\n')
    assert (result.returncode == 0) == succeeds, (name, result.returncode)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


fasta = write('q.faa', '>p_a\nMKKLL\n>p_b\nMALW\n')
genome = write('q.fna', '>chr\n' + 'ATG' * 500 + '\n')
# The wrapper retains actual sidecar outputs before the JVM removes its work directory.
wrapper = write('capture_sidecar.py', '''import csv,os,sys,subprocess,shutil
from pathlib import Path
a=sys.argv[1:]
r=subprocess.run([sys.executable,'''+repr(str(ROOT/'benchmark/neural/run_neural_predictors.py'))+''']+a)
capture=os.environ.get('CLEAN_CAPTURE')
if capture:
 dest=Path(capture);dest.mkdir(parents=True,exist_ok=True)
 manifest=Path(a[a.index('--manifest')+1])
 for row in csv.DictReader(manifest.open(),delimiter='\\t'):
  shutil.copyfile(row['fasta_path'],dest/'input.faa')
  for path in Path(row['output_dir']).glob('*.clean.*'):shutil.copyfile(path,dest/path.name)
sys.exit(r.returncode)
''')
# Resolve container script paths while keeping the production Nextflow process unchanged.
python = write('bin/python3', '#!'+sys.executable+'\nimport os,sys\na=sys.argv[1:]\n'+
    'if a and a[0]=="/opt/gspa/run_neural_predictors.py":a[0]='+repr(str(wrapper))+'\n'+
    'elif a and a[0].startswith("/opt/gspa/"):a[0]='+repr(str(ROOT/'benchmark/neural'))+'+"/"+a[0].rsplit("/",1)[1]\n'+
    'os.execv('+repr(sys.executable)+',['+repr(sys.executable)+']+a)\n')
python.chmod(0o755)
caller = write('bin/pyrodigal', '#!'+sys.executable+'''\nimport sys
from pathlib import Path
a=sys.argv[1:]
Path(a[a.index('-a')+1]).write_text('''+repr(fasta.read_text())+''')
Path(a[a.index('-o')+1]).write_text('##gff-version 3\\nchr\\tfixture\\tCDS\\t100\\t114\\t.\\t+\\t0\\tID=p_a\\nchr\\tfixture\\tCDS\\t130\\t141\\t.\\t+\\t0\\tID=p_b\\n')
''')
caller.chmod(0o755)
env = dict(os.environ, PATH=str(out/'bin')+':'+os.environ['PATH'],
           NXF_OFFLINE='true', NXF_ANSI_LOG='false', NXF_VER='25.10.4',
           JDK_JAVA_OPTIONS='-Djava.io.tmpdir='+os.environ['TMPDIR'])
java = ['java', '-jar', args.jar.resolve()]
for kind in ('gmm_confidence','distance'):
    model = out / ('model-'+kind)
    # Deliberate native process fixture, with one statically inspectable infer_maxsep call.
    driver = write(str(model.relative_to(out))+'/CLEAN_infer_fasta.py', '''import os,sys
from pathlib import Path

def infer_maxsep(training,name,gmm=None):
 ids=[line[1:].split()[0] for line in Path('data/inputs/'+name+'.fasta').read_text().splitlines() if line.startswith('>')]
 fault=os.environ.get('CLEAN_FAULT')
 if fault=='missing':return
 if fault=='unknown':ids[0]='unknown'
 values=('0.00000001','0.8') if gmm else ('0.01','8.0627')
 text=''.join(pid+',EC:1.1.1.1/'+values[0]+',EC:1.1.1.n11/'+values[1]+'\\n' for pid in ids)
 if fault=='malformed':text='malformed\\n'
 Path('results/inputs/'+name+'_maxsep.csv').write_text(text)

infer_maxsep('fixture',sys.argv[-1],gmm='''+("'fixture.pkl'" if kind=='gmm_confidence' else 'None')+''')
''')
    cfg = write('config-'+kind+'.json',json.dumps(dict(predictors=dict(
        enable=['clean'],similarity=dict(enabled=False),domains=dict(enabled=False),operons=dict(enabled=False),
        pathway=dict(enabled=False),neural=dict(pythonExecutable=sys.executable,sidecarScript=str(wrapper),
        clean=dict(enabled=True,modelDir=str(model),scoreKind=kind,minScore=0))),integration=dict(enabled=False))))
    capture = out / ('capture-'+kind)
    claims = out / ('annotate-'+kind+'.jsonl')
    run('annotate-'+kind,java+['annotate','--proteins',fasta,'--config',cfg,'--output',out/('annotation-'+kind),
        '--claims-out',claims,'--no-metrics','--clean-score-kind',kind],dict(env,CLEAN_CAPTURE=str(capture)))
    imported = out / ('import-'+kind+'.jsonl')
    run('import-'+kind,java+['clean-claims','--proteins',capture/'input.faa','--tsv',capture/'query.clean.tsv',
        '--receipt',capture/'query.clean.receipt.json','--score-kind',kind,'--output',imported],env)
    assert claims.read_bytes()==imported.read_bytes(), 'annotation/import claim parity'
    rows=[json.loads(line) for line in claims.read_text().splitlines()]
    assert len(rows)==4 and all(r['metadata']['score_available']==(kind=='gmm_confidence') for r in rows)
    assert {r['metadata']['native_score'] for r in rows}==({1e-8,.8} if kind=='gmm_confidence' else {.01,8.0627})
    if kind=='distance':
        integrated=out/'distance-integrated.tsv'
        run('integrate-distance',java+['integrate','--claims',claims,'--out',integrated,'--lite','--enable-priors','none'],env)
        assert len(integrated.read_text().splitlines())==1,'unscored distance became functional confidence'
    else:
        for fault in ('missing','unknown','malformed'):
            run('annotate-'+fault,java+['annotate','--proteins',fasta,'--config',cfg,'--output',out/('annotation-'+fault),
                '--claims-out',out/(fault+'.jsonl'),'--no-metrics'],dict(env,CLEAN_FAULT=fault))
            # annotate is currently best effort: its status manifest, not exit code, records tool failure.
            status=json.loads((out/('annotation-'+fault)/'q_predictor_status.json').read_text())
            assert status['summary']['FAILED']==1 and status['summary']['SUCCESS_ZERO']==0
            assert (out/(fault+'.jsonl')).read_bytes()==b''

    goa = write('goa.gaf','!gaf-version: 2.2\n')
    owl = write('go.owl','<?xml version="1.0"?><rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:owl="http://www.w3.org/2002/07/owl#"><owl:Ontology rdf:about="http://example.org/clean-route-fixture"/></rdf:RDF>\n')
    ec = write('ec2go','! controlled empty mapping\n')
    pathways = write('pathways.tsv','')
    config = write('nextflow.config', '''process.executor = 'local'
process.cpus = 2
process.memory = '2 GB'
docker.enabled = false
singularity.enabled = false
params.run_diamond = false
params.run_hmmer = false
params.run_barrnap = false
params.run_minced = false
params.run_visualize = false
params.run_integrate = true
params.run_clean = true
params.run_ensemble = true
params.enable_priors = 'none'
''')
    command=[args.nextflow,'-c',config,'run',ROOT/'gspa-nf/main.nf','-work-dir',out/('nf-work-'+kind),
        '-with-trace',out/('trace-'+kind+'.tsv'),'--input',genome,'--clean_model_dir',model,'--clean_score_kind',kind,
        '--gspa_jar',args.jar.resolve(),'--goa',goa,'--go_owl',owl,'--ec2go',ec,'--pathways',pathways,
        '--outdir',out/('nf-output-'+kind)]
    run('nextflow-'+kind,command,env)
    nfout=out/('nf-output-'+kind)/'q'
    merged=nfout/'q_claims.jsonl'
    # Native UUIDs differ per execution. Compare complete import bytes within each execution.
    assert merged.read_bytes()==(nfout/'clean/q_claims_clean.jsonl').read_bytes()
    nfrows=[json.loads(line) for line in merged.read_text().splitlines()]
    assert len(nfrows)==4 and all(r['metadata']['score_available']==(kind=='gmm_confidence') for r in nfrows)
    if kind=='distance':
        assert len((nfout/'q_integrated.tsv').read_text().splitlines())==1
        assert not (nfout/'ensemble/q.ensemble.tsv').exists(),'unscored distance entered confidence ensemble'
    else:
        fault=list(command)
        fault[fault.index('-work-dir')+1]=out/'nf-missing-work'
        fault[fault.index('-with-trace')+1]=out/'missing-trace.tsv'
        fault[fault.index('--outdir')+1]=out/'nf-missing-output'
        run('nextflow-missing',fault,dict(env,CLEAN_FAULT='missing'),succeeds=False)
write('report.json',json.dumps(dict(status='controlled_route_parity_verified',commands=len(commands),
    scope='Controlled native processes; real Python sidecar, CLI and Nextflow; no CLEAN model or biological accuracy claim.',
    jar_sha256=sha(args.jar),score_kinds=['gmm_confidence','distance']),indent=2)+'\n')
print('CLEAN confidence/distance routes, import parity and failure propagation verified.')
