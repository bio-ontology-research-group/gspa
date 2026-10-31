"""Capture genuine HMMER output on explicitly synthetic, unlabelled sequences.

The GSPA_TEST accession is a fixture identifier, not a Pfam accession or a
biological function claim. This exercises the external-tool format/adapter.
"""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import random
import shutil
import subprocess


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def generate(out, hmmbuild='hmmbuild', hmmsearch='hmmsearch'):
    out.mkdir(parents=True, exist_ok=False)
    tools = {}
    commands = []
    for name, binary in [('hmmbuild', hmmbuild), ('hmmsearch', hmmsearch)]:
        executable = Path(shutil.which(binary) or binary).resolve()
        proc = subprocess.run([str(executable), '-h'], text=True, capture_output=True, check=True, timeout=10)
        version = next(line[2:] for line in proc.stdout.splitlines() if line.startswith('# HMMER '))
        (out/(name+'-help.txt')).write_text(proc.stdout)
        tools[name] = {'path':str(executable),'sha256':digest(executable),'version':version}
    rng = random.Random(381)
    alphabet = 'ACDEFGHIKLMNPQRSTVWY'
    sequence = ''.join(rng.choice(alphabet) for _ in range(80))
    sequences = [sequence]
    for i in range(1, 4):
        letters = list(sequence)
        for j in range(i, len(letters), 19):
            letters[j] = alphabet[(alphabet.index(letters[j])+i) % len(alphabet)]
        sequences.append(''.join(letters))
    (out/'alignment.sto').write_text('# STOCKHOLM 1.0\n#=GF ID Gspa_test_domain\n'
        '#=GF AC GSPA_TEST_001.1\n#=GF DE Synthetic adapter-contract fixture, not a biological family\n'+
        ''.join(f'synthetic_{i} {seq}\n' for i,seq in enumerate(sequences))+'//\n')
    (out/'input.faa').write_text('>target\n'+sequence+'\n>negative\n'+'A'*80+'\n')
    for name, args in [('hmmbuild',['--amino','profile.hmm','alignment.sto']),
                       ('hmmsearch',['--domtblout','hmmer_domtbl.tsv','--noali','-E','1',
                                     '--domE','1e-5','--cpu','1','profile.hmm','input.faa'])]:
        command = [tools[name]['path'],*args]
        proc = subprocess.run(command,cwd=out,text=True,capture_output=True,check=True,timeout=30)
        (out/(name+'.stdout.txt')).write_text(proc.stdout)
        (out/(name+'.stderr.txt')).write_text(proc.stderr)
        commands.append({'argv':command,'cwd':str(out.resolve()),'returncode':proc.returncode})
    report = {'schema':'gspa.observed-tool-fixture/1','tool':'HMMER','tools':tools,'commands':commands,
              'python':platform.python_version(),'generator_sha256':digest(__file__),
              'scope':'Real HMMER execution on synthetic sequences/profile; no biological accuracy or genuine Pfam model claim.',
              'artifacts':{p.name:digest(p) for p in sorted(out.iterdir()) if p.is_file()}}
    (out/'provenance.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'out':str(out),'tool_versions':{k:v['version'] for k,v in tools.items()}}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--hmmbuild',default='hmmbuild')
    parser.add_argument('--hmmsearch',default='hmmsearch')
    args = parser.parse_args()
    generate(args.out.resolve(),args.hmmbuild,args.hmmsearch)
