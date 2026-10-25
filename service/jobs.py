"""Persistent, bounded, single-worker prediction queue for the public portal."""
from __future__ import annotations
import csv
from functools import lru_cache
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import threading
import time
import zipfile

ROOT = Path(os.environ.get('GSPA_JOBS', '/var/lib/gspa/jobs'))
ASSETS = Path(os.environ.get('DGPP_ASSETS', '/opt/dgpp_assets'))
RESOURCES = Path(os.environ.get('GSPA_RESOURCES', '/opt/gspa-resources'))
DEFAULT_ANALYSES = ('homology', 'names', 'operons', 'pathways', 'amr', 'bgc', 'localization', 'quality', 'metabolism')
BIN = os.environ.get('GSPA_BIN', 'gspa-cli')
THREADS = os.environ.get('GSPA_THREADS', '4')
TIMEOUT = int(os.environ.get('GSPA_TIMEOUT', '3600'))
MAX_BYTES = int(os.environ.get('GSPA_MAX_UPLOAD_BYTES', str(20 * 1024 * 1024)))
MAX_SEQUENCES = 20000
RETENTION = 7 * 86400
LOCK = threading.Lock()
STOP = threading.Event()
REQUIRED = ('train_db.dmnd', 'train_net_index.tsv', 'train_terms.tsv', 'go-dag.tsv', 'go.obo')


def ready():
    try:
        return (bool(shutil.which(BIN)) and all((ASSETS / f).is_file() and (ASSETS / f).stat().st_size for f in REQUIRED)
                and all((RESOURCES / p).exists() for p in ('diamond/swissprot.dmnd', 'go/go.owl', 'ec2go/ec2go',
                  'kegg_pathways.tsv', 'kegg_modules.tsv', 'reasoner-cache/has_part_pairs.tsv', 'fba_data/fba_matrices_v6RE.pkl',
                  'tool-data/amrfinder/latest', 'tool-data/antismash')))
    except OSError:
        return False


def state(folder, **changes):
    path = folder / 'status.json'
    data = json.loads(path.read_text()) if path.exists() else {}
    if changes:
        data.update(changes, updated=time.time())
        temporary = folder / 'status.tmp'
        temporary.write_text(json.dumps(data))
        temporary.replace(path)
    return data


def fasta(path):
    pid, sequence = None, []
    with path.open() as handle:
        for line in handle:
            if line.startswith('>'):
                if pid is not None:
                    yield pid, ''.join(sequence)
                pid, sequence = line[1:].split()[0], []
            else:
                sequence.append(line.strip())
        if pid is not None:
            yield pid, ''.join(sequence)


def validate(path, kind):
    allowed = set('ACGTRYSWKMBDHVN' if kind == 'genome' else 'ACDEFGHIKLMNPQRSTVWYBXZJUO*')
    ids, count, length = set(), 0, 0
    with path.open() as f:
        if not f.readline().startswith('>'):
            raise ValueError('Use an uncompressed FASTA file starting with a > header.')
    try:
        for pid, seq in fasta(path):
            if not re.fullmatch(r'[A-Za-z0-9_.:|+-]{1,120}', pid) or pid in ids:
                raise ValueError('FASTA identifiers must be unique, short, and contain only letters, numbers, or _ . : | + - .')
            if not seq or set(seq.upper()) - allowed:
                raise ValueError('A sequence is empty or contains symbols incompatible with the selected input type.')
            if kind == 'proteins' and len(seq) > 20000:
                raise ValueError('Individual proteins must be at most 20,000 residues.')
            ids.add(pid)
            count += 1
            length += len(seq)
            if count > MAX_SEQUENCES:
                raise ValueError('Please submit at most 20,000 sequences per job.')
    except (UnicodeError, IndexError):
        raise ValueError('Use an uncompressed text FASTA with a nonempty identifier on each header.')
    if not count:
        raise ValueError('No sequences found.')
    return count, length


def command(folder, stage, argv):
    state(folder, stage=stage)
    with (folder / 'run.log').open('a') as log:
        log.write('\n' + stage + '\n' + repr(argv) + '\n')
        log.flush()
        process = subprocess.Popen(argv, stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
        deadline = time.monotonic() + TIMEOUT
        while process.poll() is None:
            if STOP.wait(.5) or time.monotonic() > deadline:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
                raise RuntimeError('Prediction interrupted or exceeded the one-hour stage limit.')
        if process.returncode:
            raise RuntimeError(f'{stage} failed. Download the run log for details.')


def verify_prediction(log):
    # The generic CLI allows optional predictors to fail and continues. For this
    # service the selected predictor is mandatory, including a valid zero-hit run.
    text = log.read_text(errors='replace')
    if ('deepgo-plusplus-light failed:' in text or
            not re.search(r'deepgo-plusplus-light: \d+ annotations assigned', text)):
        raise RuntimeError('The function predictor did not complete. Download the run log for details.')


@lru_cache(maxsize=1)
def reference_metadata():
    aspects, current = {}, None
    for line in (ASSETS / 'go.obo').read_text().splitlines():
        if line.startswith('id: GO:'):
            current = line[4:]
        elif current and line.startswith('namespace: '):
            aspects[current] = {'biological_process': 'BP', 'molecular_function': 'MF',
                                'cellular_component': 'CC'}.get(line[11:], '')
    hashes = {}
    for name in REQUIRED:
        with (ASSETS / name).open('rb') as handle:
            hashes[name] = hashlib.file_digest(handle, 'sha256').hexdigest()
    return aspects, hashes


def run(folder):
    from service.analysis import run_full
    run_full(folder)


def worker():
    while not STOP.is_set():
        with LOCK:
            folders = sorted(ROOT.glob('*/status.json'), key=lambda p: state(p.parent)['created'])
            pending = []
            for path in folders:
                info = state(path.parent)
                if info['status'] in ('complete', 'failed') and not info.get('example') and time.time()-info['updated'] > RETENTION:
                    shutil.rmtree(path.parent)
                elif info['status'] == 'queued':
                    pending.append(path.parent)
        if pending:
            try:
                run(pending[0])
            except Exception as exc:
                state(pending[0], status='failed', stage='Job stopped', error=str(exc), finished=time.time())
        else:
            STOP.wait(2)
