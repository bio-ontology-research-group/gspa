#!/usr/bin/env python3
"""Real Pyrodigal/Nextflow ID transport with controlled function-context probes.

Gene calling uses observed genome sequences and the real pinned Pyrodigal wheel.
DIAMOND and GO/pathway evidence are synthetic route controls, not function truth.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import zipfile

ROOT = Path(__file__).resolve().parents[2]


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def fasta(path):
    rows = {}; headers = {}; key = None
    for line in path.read_text().splitlines():
        if line.startswith('>'):
            key = line[1:].split()[0]
            assert key not in rows
            rows[key] = ''; headers[key] = line[1:]
        elif line.strip():
            assert key is not None
            rows[key] += line.strip()
    return rows, headers


def gff(path):
    rows = {}
    for line in path.read_text().splitlines():
        if not line or line.startswith('#'):
            continue
        fields = line.split('\t')
        assert len(fields) == 9 and fields[2] == 'CDS'
        attrs = dict(p.split('=',1) for p in fields[8].split(';') if '=' in p)
        key = attrs['ID']; assert key not in rows
        rows[key] = (fields[0], int(fields[3]), int(fields[4]), fields[6], fields[7])
    return rows


def predictions(path):
    with path.open() as stream:
        return {(r['protein_id'],r['function_id']):r for r in csv.DictReader(stream,delimiter='\t')}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--observed', type=Path, required=True)
    p.add_argument('--wheels', type=Path, required=True)
    p.add_argument('--jar', type=Path, required=True)
    p.add_argument('--nextflow', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    out = args.out.resolve(); observed = args.observed.resolve()
    commands = []; summaries = {}

    def write(name, text):
        path = out/name; path.parent.mkdir(parents=True,exist_ok=True);path.write_text(text);return path

    def execute(name, command, env):
        run = subprocess.run([str(x) for x in command], cwd=out, env=env,
                             stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=900)
        write(name+'.log',run.stdout.decode(errors='replace'))
        commands.append(dict(name=name,command=[str(x) for x in command],exit_code=run.returncode))
        write('commands.json',json.dumps(commands,indent=2)+'\n')
        assert run.returncode == 0,(name,run.returncode)

    libraries = out/'libraries';libraries.mkdir()
    for wheel in args.wheels.glob('*.whl'):
        with zipfile.ZipFile(wheel) as z:z.extractall(libraries)
    caller = write('bin/pyrodigal', '#!/usr/bin/python3\nimport sys\nsys.path.insert(0,'+repr(str(libraries))+')\nfrom pyrodigal.cli import main\nraise SystemExit(main())\n')
    caller.chmod(0o755)
    probe_registry = {}
    for mode in ('single','meta'):
        sequence,headers = fasta(observed/(mode+'.faa'))
        coordinates = gff(observed/(mode+'.gff'))
        assert sequence.keys() == coordinates.keys()
        groups = [line.split('\t') for line in (observed/(mode+'-operons-before/operons_for_integrate.tsv')).read_text().splitlines()]
        assert all(set(group) <= sequence.keys() for group in groups)
        selected = next(group for group in groups if len(group)>=3)
        decoy = next(pid for pid in sequence if pid not in selected)
        probe_registry[mode] = dict(anchor=selected[0],candidate=selected[1],partner=selected[2],decoy=decoy)
    registry = write('probes.json',json.dumps(probe_registry,indent=2)+'\n')
    diamond = write('bin/diamond', '''#!/usr/bin/python3
import json,sys
from pathlib import Path
a=sys.argv[1:];query=Path(a[a.index('--query')+1]);mode=query.stem
p=json.loads(Path('''+repr(str(registry))+''').read_text())[mode]
ids={line[1:].split()[0] for line in query.read_text().splitlines() if line.startswith('>')}
assert set(p.values()) <= ids
with Path(a[a.index('--out')+1]).open('w') as f:
 for key,accession,identity in [('anchor','CONTROL_A',95),('partner','CONTROL_C',95),('candidate','CONTROL_B',35),('decoy','CONTROL_B',35)]:
  f.write(f"{p[key]}\\tsp|{accession}|controlled\\t{identity}\\t100\\t100\\t100\\t1e-20\\t100\\tcontrolled route observation\\n")
''');diamond.chmod(0o755)
    env = dict(os.environ,PATH=str(out/'bin')+':'+os.environ['PATH'],NXF_OFFLINE='true',
        NXF_ANSI_LOG='false',NXF_VER='25.10.4',NXF_HOME=str(args.nextflow.resolve().parent.parent/'.nextflow'),
        JDK_JAVA_OPTIONS='-Djava.io.tmpdir='+os.environ['TMPDIR'])
    execute('native-version',[caller,'--version'],env)
    assert 'pyrodigal v3.7.1' in (out/'native-version.log').read_text()
    goa = write('control.gaf','!gaf-version: 2.2\n'+''.join(
        f'UniProtKB\t{acc}\tcontrol\t\t{go}\tTEST:control\tEXP\t\tF\tcontrolled probe\n'
        for acc,go in [('CONTROL_A','GO:0004340'),('CONTROL_B','GO:0004332'),('CONTROL_C','GO:0004618')]))
    ontology = write('control.owl','''<?xml version="1.0"?>
<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:owl="http://www.w3.org/2002/07/owl#">
<owl:Ontology rdf:about="http://example.org/controlled-operon-route"/>
</rdf:RDF>
''')
    ec = write('ec2go','EC:2.7.1.1\tGO:0004340\nEC:4.1.2.13\tGO:0004332\nEC:5.3.1.9\tGO:0004618\n')
    pathways = write('pathways.tsv','pathway_id\tpathway_name\tgo_term\treaction_id\tec_number\tdepends_on_reaction\n'+
        'CONTROL\tSynthetic context control\tGO:0006096\tr1\tEC:2.7.1.1\t-\n'+
        'CONTROL\tSynthetic context control\tGO:0006096\tr2\tEC:4.1.2.13\tr1\n'+
        'CONTROL\tSynthetic context control\tGO:0006096\tr3\tEC:5.3.1.9\tr2\n')
    database = write('controlled.dmnd','not a real database; DIAMOND process is a declared control\n')
    java = ['java','-XX:ActiveProcessorCount=3','-jar',args.jar.resolve()]
    for mode in ('single','meta'):
        genome = write(mode+'.fna',(observed/'genome.fna').read_text())
        config = write(mode+'.config', '''process.executor = 'local'
process.cpus = 1
process.memory = '1 GB'
executor.queueSize = 2
docker.enabled = false
singularity.enabled = false
params.run_hmmer = false
params.run_barrnap = false
params.run_minced = false
params.run_visualize = false
params.run_integrate = true
params.run_diamond = true
params.run_antismash = false
params.enable_priors = 'genomic_context'
params.mag = '''+('true' if mode=='meta' else 'false')+'\n')
        published = out/(mode+'-output')
        execute(mode+'-nextflow',[args.nextflow.resolve(),'-c',config,'run',ROOT/'gspa-nf/main.nf',
            '-work-dir',out/(mode+'-work'),'-with-trace',out/(mode+'-trace.tsv'),
            '--input',genome,'--outdir',published,'--diamond_db',database,'--gspa_jar',args.jar.resolve(),
            '--goa',goa,'--go_owl',ontology,'--ec2go',ec,'--pathways',pathways],env)
        sample=published/mode;cds=sample/'gene_calling'/(mode+'.gff');faa=sample/'gene_calling'/(mode+'.faa')
        actual,headers=fasta(faa);expected,_=fasta(observed/(mode+'.faa'))
        assert actual == expected,mode+' protein sequences changed'
        coordinates=gff(cds)
        assert coordinates == gff(observed/(mode+'.gff'))
        assert actual.keys()==coordinates.keys()
        ops=sample/(mode+'_operons_for_integrate.tsv')
        assert ops.read_bytes()==(observed/(mode+'-operons-before/operons_for_integrate.tsv')).read_bytes()
        groups=[line.split('\t') for line in ops.read_text().splitlines()]
        for group in groups:
            assert len(group)>=2 and set(group)<=actual.keys()
            for a,b in zip(group,group[1:]):
                ca,cb=coordinates[a],coordinates[b]
                assert ca[0]==cb[0] and ca[3]==cb[3] and cb[1]>=ca[1] and cb[1]-ca[2]-1<=300
        raw=sample/(mode+'_claims.jsonl');claims=[json.loads(line) for line in raw.read_text().splitlines()]
        assert len(claims)==4 and {r['protein_id'] for r in claims}==set(probe_registry[mode].values())
        common=java+['integrate','--claims',raw,'--lite','--go-owl',ontology,'--ec2go',ec,'--pathways',pathways]
        direct=out/(mode+'-direct.tsv')
        execute(mode+'-direct',common+['--operons',ops,'--enable-priors','genomic_context','--out',direct],env)
        integrated=sample/(mode+'_integrated.tsv')
        assert direct.read_bytes()==integrated.read_bytes(),'Nextflow/direct probability mismatch'
        baseline=out/(mode+'-baseline.tsv')
        execute(mode+'-baseline',common+['--enable-priors','none','--out',baseline],env)
        # Counterfactual: replace full protein IDs with numeric header metadata IDs.
        aliases={pid:re.search(r'(?:^|;)ID=([^;]+)',header.split(' # ')[-1])[1] for pid,header in headers.items()}
        wrong=write(mode+'-numeric-operons.tsv',''.join('\t'.join(aliases[pid] for pid in group)+'\n' for group in groups))
        wrong_result=out/(mode+'-numeric.tsv')
        execute(mode+'-numeric-control',common+['--operons',wrong,'--enable-priors','genomic_context','--out',wrong_result],env)
        live,base,bad=predictions(integrated),predictions(baseline),predictions(wrong_result)
        probe=probe_registry[mode];candidate=(probe['candidate'],'GO:0004332');decoy=(probe['decoy'],'GO:0004332')
        assert live.keys()==base.keys()==bad.keys()
        assert float(live[candidate]['posterior_prob'])>float(base[candidate]['posterior_prob'])
        assert float(live[decoy]['posterior_prob'])==float(base[decoy]['posterior_prob'])
        assert 'genomic_context' in live[candidate]['priors_fired']
        assert all(float(base[k]['posterior_prob'])==float(bad[k]['posterior_prob']) for k in base)
        summaries[mode]=dict(proteins=len(actual),operons=len(groups),operon_members=sum(map(len,groups)),
            candidate=probe['candidate'],decoy=probe['decoy'],
            candidate_baseline=float(base[candidate]['posterior_prob']),candidate_context=float(live[candidate]['posterior_prob']),
            decoy_baseline=float(base[decoy]['posterior_prob']),decoy_context=float(live[decoy]['posterior_prob']),
            full_ids_match=True,numeric_control_loses_context=True,identical_standalone_and_nextflow=True)
    bound={str(p):sha(p) for p in [args.jar,Path(__file__),observed/'input.json',observed/'genome.fna',*args.wheels.glob('*.whl')]}
    write('report.json',json.dumps(dict(status='real_gene_calling_and_operon_route_verified',summary=summaries,
        limitation='Real Pyrodigal3.7.1 on observed genome; synthetic DIAMOND/GO/pathway controls. No measured operon or function accuracy.',
        commands=commands,inputs_sha256=bound),indent=2)+'\n')
    print(json.dumps(summaries,indent=2))


if __name__=='__main__':
    main()
