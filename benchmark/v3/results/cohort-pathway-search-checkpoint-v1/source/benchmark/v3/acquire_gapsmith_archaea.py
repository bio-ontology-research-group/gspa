"""Acquire the pinned archaeal reference in an isolated Slurm staging directory.

The outer and nested MD5 checksums are published by Zenodo record 16908828.
SHA-256 receipts additionally bind every acquired and extracted file locally.
Existing shared bacterial references are read, never modified.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import tarfile
from urllib.request import urlopen

from observe_gapsmith import TABLES, digest

RECORD = '16908828'
OUTER_MD5 = '56bdc5265a1d3eba801936711625e160'
NESTED_MD5 = {'rev':'c6e0d23c40ae8f569686ddff5310fcca',
              'rxn':'d29bfc01e1a91cb6e8ae056a3e886ad7',
              'unrev':'593efe58728d272da9b211d63e9ca9fa'}


def extract(archive, destination):
    destination.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive, 'r:gz') as tar:
        for member in tar.getmembers():
            target=(destination/member.name).resolve()
            if (not target.is_relative_to(destination.resolve()) or
                    not (member.isdir() or member.isfile())):
                raise ValueError('unsafe archive member: '+member.name)
        tar.extractall(destination, filter='data')


def run(root, tables):
    if not os.environ.get('SLURM_JOB_ID'): raise RuntimeError('acquire through Slurm')
    root.mkdir(parents=True, exist_ok=False)
    scratch=Path(os.environ['SLURM_TMPDIR'])/'archaeal-reference'
    scratch.mkdir()
    record=dict(schema='gspa.archaeal-reference-acquisition/1', status='preparing',
                zenodo_record=RECORD, job_id=os.environ['SLURM_JOB_ID'],
                started_utc=datetime.now(timezone.utc).isoformat())
    receipt=root/'acquisition.json'
    def save(): receipt.write_text(json.dumps(record,indent=2,sort_keys=True)+'\n')
    save()
    try:
        acquired={}
        for name in ['Archaea.tar.gz','md5sums.txt']:
            url=f'https://zenodo.org/records/{RECORD}/files/{name}?download=1'
            path=scratch/name
            with urlopen(url,timeout=180) as response, path.open('wb') as output:
                shutil.copyfileobj(response,output)
                final_url=response.url
            acquired[name]=dict(url=url,final_url=final_url,sha256=digest(path),bytes=path.stat().st_size)
            shutil.copy2(path,root/name)
        outer=scratch/'Archaea.tar.gz'
        if hashlib.md5(outer.read_bytes()).hexdigest()!=OUTER_MD5:
            raise ValueError('outer archive checksum mismatch')
        advertised={line.split()[1].removeprefix('./'):line.split()[0]
                    for line in (scratch/'md5sums.txt').read_text().splitlines() if line.strip()}
        for sub,h in NESTED_MD5.items():
            if advertised.get(f'Archaea/{sub}/sequences.tar.gz')!=h:
                raise ValueError('published nested checksum changed')
        reference=scratch/'reference'
        extract(outer,reference/'seq')
        for sub,h in NESTED_MD5.items():
            archive=reference/'seq/Archaea'/sub/'sequences.tar.gz'
            if hashlib.md5(archive.read_bytes()).hexdigest()!=h:
                raise ValueError('nested archive checksum mismatch: '+sub)
            extract(archive,archive.parent)
            if not list(archive.parent.glob('*.fasta')):
                raise ValueError('no extracted FASTA files: '+sub)
        sources={}
        for name in TABLES:
            source=tables/name;sources[str(source)]=digest(source)
            shutil.copy2(source,reference/name)
        user=tables/'seq/Archaea/user'
        for source in user.glob('*.fasta'):
            sources[str(source)]=digest(source)
            target=reference/'seq/Archaea/user'/source.name
            target.parent.mkdir(exist_ok=True)
            shutil.copy2(source,target)
        # Resolver marker only; no bacterial sequence files in this reference.
        (reference/'seq/Bacteria').mkdir()
        if sources!={p:digest(Path(p)) for p in sources}:
            raise ValueError('shared source changed during acquisition')
        hashes={str(p.relative_to(reference)):digest(p) for p in reference.rglob('*') if p.is_file()}
        shutil.copytree(reference,root/'reference')
        for p,h in hashes.items():
            if digest(root/'reference'/p)!=h: raise ValueError('published reference copy changed')
        (root/'reference-hashes.json').write_text(json.dumps(hashes,indent=2,sort_keys=True)+'\n')
        record.update(status='verified',downloads=acquired,source_inputs=sources,
                      reference_files=len(hashes),reference_manifest_sha256=digest(root/'reference-hashes.json'),
                      category_fastas={sub:len(list((reference/'seq/Archaea'/sub).glob('*.fasta')))
                                       for sub in NESTED_MD5})
    except Exception as exc:
        record.update(status='failed',failure_type=type(exc).__name__,failure=str(exc))
        raise
    finally:
        record['finished_utc']=datetime.now(timezone.utc).isoformat();save()


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True)
    p.add_argument('--tables',type=Path,required=True)
    a=p.parse_args();run(a.root.resolve(),a.tables.resolve())
