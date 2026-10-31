"""Observe a real gapsmith run on an exact training-genome proteome, inside Slurm.

Copies reference inputs to node-local scratch, checks their bytes before/after,
retains every output and reports tool execution independently of biological truth.
"""
import argparse
from collections import Counter
import csv
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import signal
import subprocess
import time

TABLES = ('meta_pwy.tbl', 'custom_pwy.tbl', 'exception.tbl', 'complex_subunit_dict.tsv',
          'seed_Enzyme_Class_Reactions_Aliases_unique_edited.tsv', 'altec.csv',
          'mnxref_seed-other.tsv', 'seed_Enzyme_Name_Reactions_Aliases.tsv', 'taxonomy.tbl')


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as fh:
        for block in iter(lambda: fh.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def run(root, reference, binary, diamond, expected_proteome_hash):
    if not os.environ.get('SLURM_JOB_ID'):
        raise RuntimeError('run this observation through Slurm, not on a login node')
    started = time.time()
    scratch = Path(os.environ['SLURM_TMPDIR']) / 'gapsmith-observation'
    scratch.mkdir(exist_ok=False)
    ref = scratch / 'reference'
    ref.mkdir()
    record = {'schema': 'gspa.gapsmith-observation/1', 'status': 'preparing',
              'job_id': os.environ['SLURM_JOB_ID'], 'node': platform.node(),
              'python': platform.python_version(), 'scope': 'Actual tool output on a training genome, not biological accuracy.'}
    receipt = root / 'observation.json'
    save = lambda: receipt.write_text(json.dumps(record, indent=2, sort_keys=True) + '\n')
    save()
    source_files = [reference / n for n in TABLES]
    source_files += sorted(p for p in (reference / 'seq/Bacteria').rglob('*') if p.is_file())
    if len(source_files) <= len(TABLES):
        raise ValueError('no bacterial reference sequences')
    source_hashes = {}
    for i, source in enumerate(source_files):
        h = digest(source)
        relative = source.relative_to(reference)
        target = ref / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        if digest(target) != h:
            raise ValueError('reference changed while copied: ' + str(source))
        source_hashes[str(relative)] = h
        if i % 1000 == 0:
            print(json.dumps({'copied': i, 'files': len(source_files)}), flush=True)
    (root / 'reference-hashes.json').write_text(json.dumps(source_hashes, indent=2, sort_keys=True) + '\n')
    proteome = root / 'GCA_000012525.1.faa'
    if digest(proteome) != expected_proteome_hash:
        raise ValueError('exact acquired proteome changed')
    tool = root / 'gapsmith'
    binary_hash = digest(binary)
    shutil.copy2(binary, tool)
    if digest(tool) != binary_hash:
        raise ValueError('tool changed during copy')
    tool.chmod(0o755)
    initial = {str(p): digest(p) for p in (binary, diamond, proteome, Path(__file__).resolve())}
    command = [str(tool), '--data-dir', str(ref), 'find', '-p', 'all', '-t', 'Bacteria', '-K', '4',
               '-o', str(scratch / 'output'), str(proteome)]
    env = dict(os.environ, PATH=str(diamond.parent) + os.pathsep + os.environ['PATH'])
    record.update(status='running', command=command, inputs=initial, reference_files=len(source_files),
                  tool_version=subprocess.check_output([str(tool), '--version'], text=True).strip(),
                  diamond_version=subprocess.check_output([str(diamond), 'version'], env=env, text=True).strip())
    save()
    with (root / 'tool.log').open('w') as log:
        process = subprocess.Popen(['/usr/bin/time', '-v', '-o', str(root / 'resources.txt'), *command],
                                   cwd=scratch, env=env, stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
        try:
            code = process.wait(timeout=7200)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGTERM)
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
            record.update(status='failed', reason='tool timeout', elapsed_seconds=time.time()-started)
            save()
            raise
    record['exit_code'] = code
    if (scratch / 'output').exists():
        shutil.copytree(scratch / 'output', root / 'outputs')
    if code:
        record.update(status='failed', reason='tool exited nonzero', elapsed_seconds=time.time()-started)
        save()
        raise RuntimeError('gapsmith failed; see tool.log')
    if initial != {p: digest(Path(p)) for p in initial}:
        raise ValueError('bound input changed during observation')
    for relative, h in source_hashes.items():
        if digest(reference / relative) != h or digest(ref / relative) != h:
            raise ValueError('reference changed during observation: ' + relative)
    output = root / 'outputs'
    rxn = output / 'GCA_000012525.1-all-Reactions.tbl'
    pwy = output / 'GCA_000012525.1-all-Pathways.tbl'
    with rxn.open() as fh:
        reader = csv.DictReader(fh, delimiter='\t'); columns = reader.fieldnames
        statuses = Counter(r['status'] for r in reader)
    with pwy.open() as fh:
        pathways = list(csv.DictReader(fh, delimiter='\t'))
    record.update(status='observed', reaction_columns=columns, reaction_statuses=dict(statuses),
                  pathways=len(pathways), predicted_pathways=sum(r['Prediction'].upper() == 'TRUE' for r in pathways),
                  elapsed_seconds=time.time()-started,
                  output_hashes={str(p.relative_to(root)): digest(p) for p in output.rglob('*') if p.is_file()})
    save()
    print(json.dumps({k: v for k, v in record.items() if k not in ('inputs', 'command', 'output_hashes')}), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('root', 'reference', 'binary', 'diamond'):
        parser.add_argument('--' + key, type=Path, required=True)
    parser.add_argument('--proteome-sha256', required=True)
    a = parser.parse_args()
    run(a.root.resolve(), a.reference.resolve(), a.binary.resolve(), a.diamond.resolve(), a.proteome_sha256)
