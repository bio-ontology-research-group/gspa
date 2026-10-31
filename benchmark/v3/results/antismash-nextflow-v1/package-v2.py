import hashlib,json,tarfile
from pathlib import Path
root=Path('/home/leechuck/Public/software/gspa');out=Path('/home/leechuck/tmp/gspa-antismash-nextflow-v1')
prior=json.loads(Path('/home/leechuck/tmp/gspa-antismash-observed-v1/build-source-files-v2.json').read_text())
paths={root/p for p in prior if (root/p).is_file()}
for module in ('gspa-core','gspa-cli','gspa-predictors','gspa-nf'):
 paths.update(p for p in (root/module).rglob('*') if p.is_file() and ('src' in p.parts or p.suffix in ('.nf','.config','.groovy')) and 'build' not in p.parts)
paths.update((root/'benchmark/v3').glob('*.py'));paths.update((root/'benchmark/v3/tests').glob('*.py'))
identity={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(paths)}
with tarfile.open(out/'source-v2.tar.gz','w:gz') as archive:
 for p in sorted(paths):archive.add(p,arcname=str(p).lstrip('/'),recursive=False)
 # Use the corrected launcher that was exercised in the dbCAN validation, not the old local shell launcher.
 nf=Path('/home/leechuck/Public/software/gspa/benchmark/v3/results/dbcan-route-v1/upstream/nextflow-25.10.4')
 archive.add(nf,arcname='home/leechuck/bin/nextflow',recursive=False)
 jar=Path('/home/leechuck/.nextflow/framework/25.10.4/nextflow-25.10.4-one.jar')
 archive.add(jar,arcname=str(jar).lstrip('/'),recursive=False)
 (out/'nextflow-runtime.json').write_text(json.dumps({str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (nf,jar)},indent=2)+'\n')
(out/'source-files-v2.json').write_text(json.dumps(identity,sort_keys=True,indent=2)+'\n')
print(len(paths),'files', (out/'source-v2.tar.gz').stat().st_size,'bytes')
