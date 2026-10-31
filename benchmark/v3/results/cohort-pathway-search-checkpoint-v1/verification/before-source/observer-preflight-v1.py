"""Observe a real gapsmith run on an exact acquired proteome, inside Slurm.

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

TABLES = ('seed_reactions_corrected.tsv', 'meta_pwy.tbl', 'custom_pwy.tbl', 'exception.tbl', 'complex_subunit_dict.tsv',
          'seed_Enzyme_Class_Reactions_Aliases_unique_edited.tsv', 'altec.csv',
          'mnxref_seed-other.tsv', 'seed_Enzyme_Name_Reactions_Aliases.tsv', 'taxonomy.tbl')


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as fh:
        for block in iter(lambda: fh.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def run(root, reference, binary, diamond, expected_proteome_hash,
        proteome=None, taxonomy='Bacteria', threads=4):
    if not os.environ.get('SLURM_JOB_ID'):
        raise RuntimeError('run this observation through Slurm, not on a login node')
    if taxonomy not in ('Bacteria', 'Archaea'):
        raise ValueError('taxonomy must select a supported reference directory')
    if threads < 1 or threads > int(os.environ.get('SLURM_CPUS_PER_TASK', threads)):
        raise ValueError('threads must fit the Slurm allocation')
    receipt = root / 'observation.json'
    if receipt.exists() or (root / 'outputs').exists():
        raise FileExistsError('refusing to overwrite an existing observation')
    proteome = proteome or root / 'GCA_000012525.1.faa'
    try:
        _run(root, reference, binary, diamond, expected_proteome_hash, proteome, taxonomy, threads)
    except Exception as exc:
        record = json.loads(receipt.read_text()) if receipt.exists() else {}
        record.update(status='failed', failure_type=type(exc).__name__, failure=str(exc))
        receipt.write_text(json.dumps(record, indent=2, sort_keys=True)+'\n')
        raise


def _run(root, reference, binary, diamond, expected_proteome_hash, proteome, taxonomy, threads):
    started = time.time()
    scratch = Path(os.environ['SLURM_TMPDIR']) / 'gapsmith-observation'
    scratch.mkdir(exist_ok=False)
    ref = scratch / 'reference'
    ref.mkdir()
    record = {'schema': 'gspa.gapsmith-observation/1', 'status': 'preparing',
              'job_id': os.environ['SLURM_JOB_ID'], 'node': platform.node(),
              'python': platform.python_version(), 'scope': 'Actual tool output on an acquired proteome, not biological accuracy.',
              'proteome': str(proteome), 'proteome_sha256': expected_proteome_hash,
              'taxonomy': taxonomy, 'threads': threads}
    receipt = root / 'observation.json'
    save = lambda: receipt.write_text(json.dumps(record, indent=2, sort_keys=True) + '\n')
    save()
    if digest(proteome) != expected_proteome_hash:
        raise ValueError('exact acquired proteome changed')
    tool = root / 'gapsmith'
    binary_hash = digest(binary)
    shutil.copy2(binary, tool)
    if digest(tool) != binary_hash:
        raise ValueError('tool changed during copy')
    tool.chmod(0o755)
    initial = {str(p): digest(p) for p in (binary, diamond, proteome, Path(__file__).resolve())}
    env = dict(os.environ, PATH=str(diamond.parent) + os.pathsep + os.environ['PATH'])
    (ref / 'seq' / taxonomy).mkdir(parents=True)
    source_files = [reference / n for n in TABLES]
    source_files += sorted(p for p in (reference / 'seq' / taxonomy).rglob('*') if p.is_file())
    if len(source_files) <= len(TABLES):
        raise ValueError('no reference sequences for ' + taxonomy)
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
        if i == len(TABLES) - 1:
            # Resolve the staged data root and parse all tables before copying
            # tens of thousands of sequence files. No pathway matches this regex.
            preflight = [str(tool), '--data-dir', str(ref), 'find', '-p', '^$',
                         '-t', taxonomy, '-K', str(threads), '-o', str(scratch / 'preflight'), str(proteome)]
            with (root / 'preflight.log').open('w') as log:
                checked = subprocess.run(preflight, cwd=scratch, env=env, stdout=log,
                                         stderr=subprocess.STDOUT, timeout=60)
            record.update(preflight_command=preflight, preflight_exit_code=checked.returncode)
            if checked.returncode:
                record.update(status='failed', reason='reference preflight failed',
                              elapsed_seconds=time.time()-started)
                save()
                raise RuntimeError('reference preflight failed; see preflight.log')
            save()
        if i % 1000 == 0:
            print(json.dumps({'copied': i, 'files': len(source_files)}), flush=True)
    (root / 'reference-hashes.json').write_text(json.dumps(source_hashes, indent=2, sort_keys=True) + '\n')
    command = [str(tool), '-v', '--data-dir', str(ref), 'find', '-p', 'all', '-t', taxonomy, '-K', str(threads),
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
    rxn = output / (proteome.stem + '-all-Reactions.tbl')
    pwy = output / (proteome.stem + '-all-Pathways.tbl')
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
    parser.add_argument('--proteome', type=Path)
    parser.add_argument('--taxonomy', choices=['Bacteria','Archaea'], default='Bacteria')
    parser.add_argument('--threads', type=int, default=4)
    a = parser.parse_args()
    run(a.root.resolve(), a.reference.resolve(), a.binary.resolve(), a.diamond.resolve(), a.proteome_sha256,
        a.proteome.resolve() if a.proteome else None, a.taxonomy, a.threads)
