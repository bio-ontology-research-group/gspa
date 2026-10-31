#!/usr/bin/env python3
"""Observed antiSMASH JSON through actual annotation/import/Nextflow integration.

The antiSMASH and gene-calling processes replay captured data. This checks workflow
transport and coordinate interpretation, not new tool execution or biological accuracy.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess

ROOT = Path(__file__).resolve().parents[2]


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--jar', type=Path, required=True)
    p.add_argument('--before-jar', type=Path, required=True)
    p.add_argument('--native', type=Path, required=True)
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

    def execute(name, command, env, succeeds=True):
        result = subprocess.run([str(x) for x in command], cwd=out, env=env,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=900)
        write(name + '.log', result.stdout.decode(errors='replace'))
        commands.append(dict(name=name, command=[str(x) for x in command], exit_code=result.returncode))
        write('commands.json', json.dumps(commands, indent=2) + '\n')
        assert (result.returncode == 0) == succeeds, (name, result.returncode)

    data = json.loads(args.native.read_text())
    assert data['schema'] == 4 and data['version'] == '8.0.4'
    zero_record = next(r for r in data['records'] if r['id'] == 'QZKU01000103.1')
    assert zero_record['areas'] == []
    samples = {}
    expected = {}
    counts = {}
    for name, payload in [('q', data), ('zero', dict(data, records=[zero_record]))]:
        # Keep the full observed JSON's original bytes and provenance for q.
        native = write(f'{name}.json', args.native.read_text() if name == 'q' else json.dumps(payload))
        nt, aa, gff, annotation_gff = [], [], ['##gff-version 3'], ['##gff-version 3']
        expectations = set()
        count = 0
        for ri, record in enumerate(payload['records'], 1):
            cid = record.get('original_id') or record['id']
            nt += ['>' + cid, record['seq']['data']]
            for pi, cds in enumerate((f for f in record['features'] if f['type'] == 'CDS'), 1):
                match = re.fullmatch(r'\[(\d+):(\d+)\]\(([+-])\)', cds['location'])
                assert match, cds['location']
                start, end, strand = int(match[1]) + 1, int(match[2]), match[3]
                pid, native_id = f'{cid}_{pi}', f'{ri}_{pi}'
                aa += [f'>{pid} # {start} # {end} # {1 if strand == "+" else -1} # ID={native_id};partial=00',
                       cds['qualifiers']['translation'][0]]
                row = f'{cid}\trecorded-CDS\tCDS\t{start}\t{end}\t.\t{strand}\t0\tID='
                gff.append(row + native_id)
                annotation_gff.append(row + pid)
                count += 1
                # Independent oracle uses native region SeqFeatures, not the parser's areas.
                for region in (f for f in record['features'] if f['type'] == 'region'):
                    rm = re.fullmatch(r'\[(\d+):(\d+)\]\(([+-])\)', region['location'])
                    assert rm, region['location']
                    rs, re_ = int(rm[1]) + 1, int(rm[2])
                    if start <= re_ and end >= rs:
                        rid = cid + ':region' + region['qualifiers']['region_number'][0]
                        overlap = 'contained' if start >= rs and end <= re_ else 'partial'
                        expectations.add((pid, rid, rs, re_, overlap))
        samples[name] = dict(native=str(native),
            genome=str(write(name + '.fna', '\n'.join(nt) + '\n')),
            proteins=str(write(name + '.faa', '\n'.join(aa) + '\n')),
            gff=str(write(name + '.gff', '\n'.join(gff) + '\n')),
            annotation_gff=str(write(name + '-annotation.gff', '\n'.join(annotation_gff) + '\n')))
        expected[name], counts[name] = expectations, count
    assert len(expected['q']) == 70 and counts['q'] == 4915
    assert not expected['zero'] and counts['zero'] == 6
    registry = write('samples.json', json.dumps(samples, indent=2))
    tool = write('bin/antismash', '''#!/usr/bin/env python3
import json,os,shutil,sys
from pathlib import Path
if '--version' in sys.argv:print('antiSMASH 8.0.4 captured-output replay');sys.exit(0)
a=sys.argv[1:];registry=json.loads(Path(''' + repr(str(registry)) + ''').read_text())
def fasta(path):
 out={};key=None
 for line in Path(path).read_text().splitlines():
  if line.startswith('>'):key=line[1:].split()[0];out[key]=''
  elif line:out[key]+=line.upper()
 return out
input=fasta(a[0]);matches=[s for s in registry.values() if fasta(s['genome'])==input]
assert len(matches)==1, 'unexpected replay input'
assert a[a.index('--hmmdetection-strictness')+1]=='relaxed'
assert int(a[a.index('--cpus')+1]) in (1,2)
d=Path(a[a.index('--output-dir')+1]);d.mkdir(parents=True)
if os.environ.get('ANTISMASH_REPLAY_FAULT')=='missing':sys.exit(0)
shutil.copyfile(matches[0]['native'],d/'recorded-result.json')
''')
    caller = write('bin/pyrodigal', '''#!/usr/bin/env python3
import json,shutil,sys
from pathlib import Path
a=sys.argv[1:];registry=json.loads(Path(''' + repr(str(registry)) + ''').read_text())
s=registry[Path(a[a.index('-i')+1]).stem]
assert Path(a[a.index('-i')+1]).read_bytes()==Path(s['genome']).read_bytes()
shutil.copyfile(s['gff'],a[a.index('-o')+1]);shutil.copyfile(s['proteins'],a[a.index('-a')+1])
''')
    tool.chmod(0o755); caller.chmod(0o755)
    env = dict(os.environ, PATH=str(out/'bin') + ':' + os.environ['PATH'],
               NXF_OFFLINE='true', NXF_ANSI_LOG='false', NXF_VER='25.10.4',
               NXF_HOME=str(args.nextflow.resolve().parent.parent/'.nextflow'),
               JDK_JAVA_OPTIONS='-Djava.io.tmpdir=' + os.environ['TMPDIR'])
    java = ['java', '-XX:ActiveProcessorCount=3', '-jar', args.jar.resolve()]
    q = samples['q']
    import_args = ['antismash-claims', '--results', q['native'], '--genome', q['genome'],
                   '--gff3', q['gff'], '--proteins', q['proteins']]
    execute('before-import', ['java', '-jar', args.before_jar.resolve()] + import_args +
            ['--output', out/'before.jsonl'], env, succeeds=False)
    assert not (out/'before.jsonl').exists()
    assert 'antismash-claims' in (out/'before-import.log').read_text()
    imported = out/'import.jsonl'
    execute('import', java + import_args + ['--output', imported], env)
    config = write('annotation.json', json.dumps(dict(input=dict(kingdom='bacteria', requireExistingGenes=True),
        predictors=dict(similarity=dict(enabled=False),domains=dict(enabled=False),pathway=dict(enabled=False),
            operons=dict(enabled=False),localization=dict(enabled=False),enable=['antismash'],
            antismash=dict(executablePath=str(tool),threads=1)),quality=dict(enabled=False))))
    raw = out/'annotation.jsonl'
    execute('annotate', java + ['annotate', '--genome', q['genome'], '--gff3', q['annotation_gff'],
        '--proteins', q['proteins'], '--config', config, '--output', out/'annotate', '--claims-out', raw,
        '--threads', '1', '--no-metrics'], env)
    assert imported.read_bytes() == raw.read_bytes(), 'annotation/import raw claims mismatch'
    goa = write('goa.gaf', '!gaf-version: 2.2\n')
    owl = write('go.owl', '''<?xml version="1.0"?>
<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:owl="http://www.w3.org/2002/07/owl#">
<owl:Ontology rdf:about="http://example.org/antismash-route-fixture"/>
</rdf:RDF>
''')
    ec = write('ec2go', '! controlled empty mapping\n')
    pathways = write('pathways.tsv', '')
    database = out/'database'; database.mkdir()
    sheet = write('samples.csv', 'sample_id,fasta\n' + ''.join(f'{n},{s["genome"]}\n' for n,s in samples.items()))
    nfconfig = write('nextflow.config', '''process.executor = 'local'
process.cpus = 1
process.memory = '1 GB'
executor.queueSize = 2
docker.enabled = false
singularity.enabled = false
params.run_diamond = false
params.run_hmmer = false
params.run_barrnap = false
params.run_minced = false
params.run_visualize = false
params.run_integrate = true
params.run_antismash = true
params.enable_priors = 'none'
''')
    command = [args.nextflow.resolve(), '-c', nfconfig, 'run', ROOT/'gspa-nf/main.nf',
        '-work-dir', out/'nf-work', '-with-trace', out/'trace.tsv', '--input', sheet,
        '--antismash_db', database, '--gspa_jar', args.jar.resolve(), '--goa', goa, '--go_owl', owl,
        '--ec2go', ec, '--pathways', pathways, '--outdir', out/'nextflow-output']
    execute('nextflow', command, env)
    for name in samples:
        folder = out/'nextflow-output'/name
        claims = folder/(name+'_claims.jsonl')
        rows = [json.loads(line) for line in claims.read_text().splitlines()]
        assert {(r['protein_id'], r['metadata']['region_id'], r['metadata']['native_start'] + 1,
                 r['metadata']['native_end'], r['metadata']['cdsOverlap']) for r in rows} == expected[name]
        assert len(rows) == len(expected[name])
        assert all(r['function_type']=='SECONDARY_METABOLITE' and r['metadata']['score_available'] is False
                   and r.get('calibrated_prob') is None for r in rows)
        receipt = json.loads((folder/'antismash'/(name+'_claims_antismash.jsonl.import.json')).read_text())
        assert receipt['coverage'] == 'complete'
        assert receipt['protein_count'] == counts[name]
        assert receipt['region_count'] == (3 if name == 'q' else 0)
        assert len((folder/(name+'_integrated.tsv')).read_text().splitlines()) == 1
        if name == 'q':
            assert claims.read_bytes() == imported.read_bytes(), 'Nextflow changed/dropped raw memberships'
    execute('integrate', java + ['integrate', '--claims', imported, '--out', out/'direct.tsv',
                                '--lite', '--enable-priors', 'none'], env)
    assert (out/'direct.tsv').read_bytes() == (out/'nextflow-output/q/q_integrated.tsv').read_bytes()
    fault = list(command)
    for option,value in [('-work-dir', out/'nf-missing-work'), ('-with-trace', out/'missing-trace.tsv'),
                         ('--outdir', out/'missing-output')]:
        fault[fault.index(option)+1] = value
    execute('nextflow-missing', fault, dict(env, ANTISMASH_REPLAY_FAULT='missing'), succeeds=False)
    assert not list((out/'missing-output').rglob('*_integrated.tsv'))
    write('report.json', json.dumps(dict(status='observed_replay_route_parity_verified',
        limitation='Captured antiSMASH and CDS replay; no fresh native inference or biological accuracy claim.',
        commands=commands, jar_sha256=sha(args.jar), native_sha256=sha(args.native),
        checker_sha256=sha(Path(__file__)), protein_counts=counts, memberships={'q':70,'zero':0},
        independent_oracle='native region SeqFeatures, not areas',
        annotated_imported_nextflow_bytes_equal=True, functional_predictions=0), indent=2)+'\n')
    print('Observed annotation/import/Nextflow parity, zero-region sample and missing-output rejection passed.')


if __name__ == '__main__':
    main()
