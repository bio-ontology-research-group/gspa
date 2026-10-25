#!/usr/bin/env python3
"""Copy only build inputs into an otherwise empty tutorial workspace.

Includes working-tree edits so unpublished features can be demonstrated, but no
host JAR, build cache, reference database or virtualenv can enter the image.
"""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

source = Path(__file__).resolve().parents[2]
root = Path(sys.argv[1]).resolve()
root.mkdir(parents=True, exist_ok=True)
if shutil.disk_usage(root).free < 25 * 1024**3:
    raise SystemExit('Need at least 25 GiB free for CPU + optional CUDA tutorial images.')
context = root / 'source'
if context.exists():
    if not (root / 'source-manifest.json').is_file():
        raise SystemExit('Refusing to replace an unmanaged source directory.')
    shutil.rmtree(context)
context.mkdir()
paths = ['gradlew', 'gradle', 'gradle.properties', 'settings.gradle.kts', 'build.gradle.kts',
         'gspa', 'gspa-core', 'gspa-predictors', 'gspa-cli', 'deepgo-plusplus/models',
         'scripts/gspa_resources.py']
paths += [str(p.relative_to(source)) for p in sorted((source / 'benchmark/neural').glob('run_*predictors.py'))]
paths += ['benchmark/neural/build_esm2_centroids.py', 'tutorials/workflows/Dockerfile']
for name in paths:
    src, dst = source / name, context / name
    dst.parent.mkdir(parents=True, exist_ok=True)
    if src.is_dir():
        shutil.copytree(src, dst, dirs_exist_ok=True,
                        ignore=shutil.ignore_patterns('build', '.gradle', '__pycache__', '*.pyc'))
    else:
        shutil.copy2(src, dst)
files = {str(p.relative_to(context)): hashlib.sha256(p.read_bytes()).hexdigest()
         for p in sorted(context.rglob('*')) if p.is_file()}
(root / 'source-manifest.json').write_text(json.dumps({
    'base_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=source, text=True).strip(),
    'description': 'Working-tree source snapshot, including uncommitted usability work',
    'files': files}, indent=2) + '\n')
print('Build inputs copied; no host JAR, databases or Gradle caches included.')
print('Storage check: at least 25 GiB available.')
