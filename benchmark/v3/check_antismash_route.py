"""Replay captured antiSMASH output through actual GSPA CLI routes.

The external process is a declared replay stub, not a new antiSMASH execution.
Checks ingestion, coordinates, raw evidence and route parity, not biological accuracy.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def run(jar, before_jar, native, java, out):
    out.mkdir(parents=True, exist_ok=False)
    data = json.loads(native.read_text())
    assert data['schema'] == 4 and data['version'] == '8.0.4'
    fasta, gff, proteins = (out/name for name in ('input.fna', 'input.gff3', 'input.faa'))
    expected = set()
    region_records = []
    n_proteins = 0
    with fasta.open('w') as nt, gff.open('w') as gf, proteins.open('w') as aa:
        gf.write('##gff-version 3\n')
        for record in data['records']:
            cid = record.get('original_id') or record['id']
            nt.write('>'+cid+'\n'+record['seq']['data']+'\n')
            gf.write(f'##sequence-region {cid} 1 {len(record["seq"]["data"])}\n')
            for area in record['areas']:
                assert area['start'] < area['end'], 'this observed replay uses linear regions'
                region_records.append((cid, area['start']+1, area['end']))
            for feature in record['features']:
                if feature['type'] != 'CDS':
                    continue
                match = re.fullmatch(r'\[(\d+):(\d+)\]\(([+-])\)', feature['location'])
                assert match, feature['location']
                start, end = int(match[1]), int(match[2])
                pid = 'GSPA_'+feature['qualifiers']['locus_tag'][0]
                gf.write(f'{cid}\tobserved\tCDS\t{start+1}\t{end}\t.\t{match[3]}\t0\tID={pid};protein_id={pid}\n')
                aa.write('>'+pid+'\n'+feature['qualifiers']['translation'][0]+'\n')
                n_proteins += 1
                for i, area in enumerate(record['areas']):
                    if start < area['end'] and end > area['start']:
                        expected.add((pid, cid+':region'+str(i+1), ', '.join(area['products'])))
    assert len(region_records) == 3 and len(expected) == 70
    tool = out/'antismash-replay'
    tool.write_text('''#!/usr/bin/env python3
import json,sys,shutil
from pathlib import Path
if '--version' in sys.argv:
 print('antiSMASH 8.0.4 captured-output replay');sys.exit(0)
args=sys.argv[1:]
source=Path('''+repr(str(native))+''')
reference=json.loads(source.read_text())
expected={r.get('original_id') or r['id']:r['seq']['data'].upper() for r in reference['records']}
actual={};key=None
for line in Path(args[0]).read_text().splitlines():
 if line.startswith('>'):key=line[1:].split()[0];actual[key]=''
 else:actual[key]+=line.upper()
assert actual==expected, 'replay input differs from observed record sequences'
assert args[args.index('--cpus')+1]=='1'
out=Path(args[args.index('--output-dir')+1]);out.mkdir(parents=True)
shutil.copyfile(source,out/'genome.json')
''')
    tool.chmod(0o755)
    config = out/'config.json'
    config.write_text(json.dumps(dict(input=dict(kingdom='bacteria', requireExistingGenes=True),
        predictors=dict(similarity=dict(enabled=False),domains=dict(enabled=False),pathway=dict(enabled=False),
            operons=dict(enabled=False),localization=dict(enabled=False),enable=['antismash'],
            antismash=dict(executablePath=str(tool),threads=3)),quality=dict(enabled=False)),indent=2)+'\n')
    workflow = out/'workflow.json'
    workflow.write_text(json.dumps(dict(lite=True, enabledPriors='none', genomeId='input'))+'\n')
    bound = {str(p):sha(p) for p in (jar,before_jar,native,java,Path(__file__),fasta,gff,proteins,tool,config,workflow)}
    env = dict(os.environ, TMPDIR=str(out/'tmp'))
    Path(env['TMPDIR']).mkdir()
    runs, results, statuses = [], {}, {}
    def execute(name, command):
        with (out/(name+'.log')).open('w') as log:
            result = subprocess.run(command,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=600)
        runs.append(dict(name=name,command=command,exit_code=result.returncode))
        (out/'runs.json').write_text(json.dumps(runs,indent=2)+'\n')
        assert result.returncode == 0, name
    for name, engine, integrated in [('before',before_jar,False),('annotate',jar,False),('annotate_integrated',jar,True)]:
        folder=out/name
        command=[str(java),'-XX:ActiveProcessorCount=3','-jar',str(engine),'annotate',
            '--genome',str(fasta),'--gff3',str(gff),'--proteins',str(proteins),
            '--config',str(config),'--output',str(folder),'--claims-out',str(folder/'raw.jsonl'),
            '--threads','1','--no-metrics']
        if integrated:command += ['--integration-config',str(workflow)]
        execute(name,command)
        status=json.loads((folder/'input_predictor_status.json').read_text())['predictors']
        rows=[json.loads(line) for line in (folder/'raw.jsonl').read_text().splitlines()]
        statuses[name]=status
        if name=='before':
            assert not rows and not any(s['predictor']=='antismash' for s in status)
            continue
        assert len(status)==1 and status[0]['predictor']=='antismash' and status[0]['outcome']=='SUCCESS',status
        assert status[0]['annotation_count']==3
        assert len(rows)==70
        assert {(r['protein_id'],r['metadata']['region_id'],r['function_id']) for r in rows}==expected
        assert all(r['function_type']=='SECONDARY_METABOLITE' and r.get('calibrated_prob') is None
            and r['metadata']['score_available'] is False and r['metadata']['source_sha256']==sha(native) for r in rows)
        results[name]=rows
        outputs=list(folder.glob('*.gff3'))
        assert len(outputs)==1,outputs
        annotation_text=outputs[0].read_text().split('##FASTA',1)[0]
        records=[line.split('\t') for line in annotation_text.splitlines() if line and not line.startswith('#')]
        assert all(len(row)==9 for row in records), 'invalid GFF feature row'
        genomic=[r for r in records if r[2]=='region']
        assert sorted((r[0],int(r[3]),int(r[4])) for r in genomic)==sorted(region_records)
        assert all(r[5]=='.' for r in genomic)
        assert sum(r[2]=='CDS' for r in records)==n_proteins
    canonical=lambda rows:sorted(json.dumps(r,sort_keys=True) for r in rows)
    assert canonical(results['annotate'])==canonical(results['annotate_integrated'])
    direct=out/'direct.tsv'
    execute('direct',[str(java),'-jar',str(jar),'integrate','--claims',str(out/'annotate/raw.jsonl'),
        '--workflow-config',str(workflow),'--out',str(direct)])
    assert direct.read_bytes()==(out/'annotate_integrated/input_integrated.tsv').read_bytes()
    assert len(direct.read_text().splitlines())==1, 'unscored regions must not become functional calls'
    for p,digest in bound.items():assert sha(p)==digest,p
    report=dict(status='verified',scope=__doc__,observed_version=data['version'],records=len(data['records']),
        input_proteins=n_proteins,regions=3,raw_memberships=70,cli_runs=len(runs),predictor_statuses=statuses,
        raw_route_parity=True,integrated_route_parity=True,invented_function_calls=0,inputs=bound,
        outputs={str(p.relative_to(out)):sha(p) for p in out.rglob('*') if p.is_file()})
    (out/'verification.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:report[k] for k in ('status','records','input_proteins','regions','raw_memberships','cli_runs')}))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('jar','before-jar','native','java','out'):parser.add_argument('--'+name,type=Path,required=True)
    run(**{k:v.resolve() for k,v in vars(parser.parse_args()).items()})
