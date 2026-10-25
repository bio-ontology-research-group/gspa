#!/usr/bin/env python3
"""Record the actual local tutorial images, package versions and disk footprint."""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

root = Path('.tutorial-work')
images = {}
for tag in ('cpu', 'neural-cpu', 'neural-cuda'):
    name = 'gspa-tutorial:' + tag
    result = subprocess.run(['docker','image','inspect',name], capture_output=True,text=True)
    if result.returncode:
        continue
    image = json.loads(result.stdout)[0]
    images[name] = {'id':image['Id'], 'size_bytes':image['Size'], 'architecture':image['Architecture']}
    if tag.startswith('neural-'):
        images[name]['packages'] = subprocess.check_output(['docker','run','--rm',name,
            '/opt/neural/bin/python','-c',
            'import importlib.metadata as m,json; print(json.dumps({d.metadata["Name"]:d.version for d in m.distributions()},sort_keys=True))'],text=True).strip()
used = sum(p.stat().st_size for p in root.rglob('*') if p.is_file() and not p.is_symlink())
helpers = {p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(Path('tutorials/workflows').iterdir()) if p.is_file()}
record = {'tutorial_files_sha256':helpers, 'images':images, 'workspace_bytes':used, 'free_bytes':shutil.disk_usage(root).free}
(root/'environment.json').write_text(json.dumps(record,indent=2)+'\n')
print(f'Workspace files: {used / 1024**3:.2f} GiB')
print(f'Free space: {record["free_bytes"] / 1024**3:.1f} GiB')
for name, image in images.items(): print(f'{name}: {image["size_bytes"] / 1024**3:.2f} GiB (layers may be shared)')
