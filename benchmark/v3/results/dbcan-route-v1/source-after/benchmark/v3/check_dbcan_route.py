#!/usr/bin/env python3
"""Execute controlled dbCAN annotation/import/Nextflow routes, without model inference."""
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


def run(name, command, env=None, succeeds=True):
    result = subprocess.run([str(x) for x in command], cwd=out, env=env,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=600)
    write(name + '.log', result.stdout.decode(errors='replace'))
    commands.append(dict(name=name, command=[str(x) for x in command], exit_code=result.returncode))
    write('commands.json', json.dumps(commands, indent=2) + '\n')
    assert (result.returncode == 0) == succeeds, (name, result.returncode)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


fasta = write('q.faa', '>p_a\nMKKLL\n>p_b\nMKKLL\n>p_c\nMKKLL\n')
genome = write('q.fna', '>chr\n' + 'ATG' * 500 + '\n')
fixture = (ROOT / 'gspa-predictors/src/test/resources/observed-output/dbcan-formatter-v4/overview-standard.txt').read_text()
overview = write('overview.tsv', fixture.replace('HMMER', 'dbCAN_hmm'))
database = out / 'database'
database.mkdir()
native = write('bin/run_dbcan', '''#!/usr/bin/env python3
import json,os,sys
from pathlib import Path
a=sys.argv[1:]
modern=a and a[0]=='CAZyme_annotation'
if '--help' in a:
 print('--input_raw_data --mode --db_dir --threads' if modern else '--db_dir --hmm_cpu --dia_cpu --dbcan_thread');sys.exit(0)
if a==['version']:print('dbCAN version: 5.2.8');sys.exit(0)
if '--version' in a:sys.exit(2)
flag='--output_dir' if modern else '--out_dir'
d=Path(a[a.index(flag)+1]);d.mkdir(parents=True,exist_ok=True)
if os.environ.get('DBCAN_ROUTE_FAULT')=='missing':sys.exit(0)
text=Path(''' + repr(str(overview)) + ''').read_text()
if not modern:text=text.replace('dbCAN_hmm','HMMER')
(d/('overview.tsv' if modern else 'overview.txt')).write_text(text)
''')
native.chmod(0o755)
caller = write('bin/pyrodigal', '''#!/usr/bin/env python3
import sys
from pathlib import Path
a=sys.argv[1:]
Path(a[a.index('-a')+1]).write_text(''' + repr(fasta.read_text()) + ''')
Path(a[a.index('-o')+1]).write_text('##gff-version 3\\n'+''.join('chr\\tfixture\\tCDS\\t%d\\t%d\\t.\\t+\\t0\\tID=p_%s\\n'%(100+i*30,114+i*30,c) for i,c in enumerate('abc')))
''')
caller.chmod(0o755)
env = dict(os.environ, PATH=str(out / 'bin') + ':' + os.environ['PATH'],
           NXF_OFFLINE='true', NXF_ANSI_LOG='false', NXF_VER='25.10.4',
           JDK_JAVA_OPTIONS='-Djava.io.tmpdir=' + os.environ['TMPDIR'])
java = ['java', '-jar', args.jar.resolve()]
for version in ('4', '5'):
    cfg = write('config-' + version + '.json', json.dumps(dict(
        predictors=dict(enable=['dbcan'], similarity=dict(enabled=False), domains=dict(enabled=False),
                        operons=dict(enabled=False), pathway=dict(enabled=False),
                        dbcan=dict(database=str(database), executablePath=str(native), interfaceVersion=version)),
        integration=dict(enabled=False))))
    claims = out / ('annotate-' + version + '.jsonl')
    run('annotate-' + version, java + ['annotate', '--proteins', fasta, '--config', cfg,
        '--output', out / ('annotation-' + version), '--claims-out', claims, '--threads', '2', '--no-metrics'], env)
    imported = out / ('import-' + version + '.jsonl')
    native_output = write('native-' + version + '.tsv', fixture if version == '4' else overview.read_text())
    run('import-' + version, java + ['dbcan-claims', '--proteins', fasta, '--overview', native_output,
                                   '--output', imported], env)
    assert claims.read_bytes() == imported.read_bytes(), 'annotation/import mismatch for interface ' + version
    rows = [json.loads(line) for line in claims.read_text().splitlines()]
    assert len(rows) == 4 and all(r['metadata']['score_available'] is False for r in rows)
    integrated = out / ('integrated-' + version + '.tsv')
    run('integrate-' + version, java + ['integrate', '--claims', imported, '--out', integrated,
                                      '--lite', '--enable-priors', 'none'], env)
    assert 'Loaded 4 claims' in (out / ('integrate-' + version + '.log')).read_text()
    assert all(isinstance(json.loads(r['metadata']['rawFieldsJson']), dict) for r in rows)
    assert len(integrated.read_text().splitlines()) == 1, 'unscored observations became confidence predictions'

goa = write('goa.gaf', '!gaf-version: 2.2\n')
owl = write('go.owl', '''<?xml version="1.0"?>
<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:owl="http://www.w3.org/2002/07/owl#">
<owl:Ontology rdf:about="http://example.org/dbcan-route-fixture"/>
</rdf:RDF>
''')
ec = write('ec2go', '! controlled empty mapping\n')
pathways = write('pathways.tsv', '')
# A component with no observations must not abort unrelated predictors.
empty_results = out / 'empty-builtin-inputs'
empty_results.mkdir()
empty_command = [sys.executable, ROOT / 'benchmark/02b_parse_predictors_to_claims.py',
                 '--results-dir', empty_results, '--goa', goa,
                 '--output', out / 'empty-builtin.jsonl']
run('empty-builtin-default', empty_command, env, succeeds=False)
assert 'ERROR: no claims produced' in (out / 'empty-builtin-default.log').read_text()
run('empty-builtin-allowed', empty_command + ['--allow-empty'], env)
assert (out / 'empty-builtin.jsonl').read_bytes() == b''
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
params.run_dbcan = true
params.enable_priors = 'none'
''')
command = [args.nextflow, '-c', config, 'run', ROOT / 'gspa-nf/main.nf',
    '-work-dir', out / 'nf-work', '-with-trace', out / 'trace.tsv',
    '--input', genome, '--dbcan_db', database, '--gspa_jar', args.jar.resolve(),
    '--goa', goa, '--go_owl', owl, '--ec2go', ec, '--pathways', pathways,
    '--outdir', out / 'nextflow-output']
run('nextflow', command, env)
merged = out / 'nextflow-output/q/q_claims.jsonl'
assert merged.read_bytes() == (out / 'import-5.jsonl').read_bytes(), 'Nextflow dropped or changed dbCAN observations'
assert len((out / 'nextflow-output/q/q_integrated.tsv').read_text().splitlines()) == 1
# A clean native exit without the mandatory overview must fail, not become zero hits.
fault = list(command)
fault[fault.index('-work-dir') + 1] = out / 'nf-missing-work'
fault[fault.index('-with-trace') + 1] = out / 'missing-trace.tsv'
fault[fault.index('--outdir') + 1] = out / 'missing-output'
run('nextflow-missing-output', fault, dict(env, DBCAN_ROUTE_FAULT='missing'), succeeds=False)
write('report.json', json.dumps(dict(status='controlled_route_parity_verified',
    scope='Controlled native processes and real GSPA/Nextflow execution; no dbCAN model or biological accuracy claim.',
    commands=len(commands), native_interfaces=['4', '5'], raw_observations_per_success=4,
    unscored_observations_promoted=0, jar_sha256=sha(args.jar),
    nextflow_merged_claims_sha256=sha(merged)), indent=2) + '\n')
print('All annotation/import/integration comparisons and full Nextflow success/failure checks passed.')
