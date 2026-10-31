"""Run a registered Price cell on local disk, preserving its published paths.

Singularity binds a private directory over one logical attempt parent. Production
commands and receipts retain their final /data paths. After execution, publication
copies and hashes every file outside collection/attempts, then renames the checked
tree into place. It never edits a production receipt or accepts partial success.
"""
import argparse
import fcntl
import json
import os
from pathlib import Path
import re
import shutil
import socket
import stat
import subprocess
import sys
import uuid

import run_price_production_collection as P

MOUNT_GUARD = '''import os,sys
if not os.path.samefile(sys.argv[1],sys.argv[2]):
    raise RuntimeError("private output bind is not active; inference was not started")
os.execv(sys.executable,[sys.executable]+sys.argv[3:])
'''


def tree_hashes(root):
    result = {}
    for path in sorted(root.rglob('*')):
        if path.is_symlink():
            raise ValueError(f'symlink in publication: {path}')
        if path.is_dir():
            continue
        if not stat.S_ISREG(path.stat().st_mode):
            raise ValueError(f'non-regular publication artifact: {path}')
        result[str(path.relative_to(root))] = P.C.V.sha256_file(str(path))
    return result


def preflight(local, destination, staging_root):
    attempts = destination.parent.parent
    if staging_root == attempts or attempts in staging_root.parents:
        raise ValueError('publication staging must be outside the attempts directory')
    if local == destination or local in destination.parents or destination in local.parents:
        raise ValueError('local and published trees must be separate')
    if staging_root == local or local in staging_root.parents or staging_root in local.parents:
        raise ValueError('local and staging trees must be separate')
    if destination.exists():
        raise FileExistsError(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    staging_root.mkdir(parents=True, exist_ok=True)
    if destination.parent.stat().st_dev != staging_root.stat().st_dev:
        raise ValueError('publication staging and destination must share a filesystem')


def sync_directory(path):
    descriptor = os.open(path, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def sync_tree(root):
    for path in root.rglob('*'):
        if path.is_file():
            with path.open('rb') as stream:
                os.fsync(stream.fileno())
    for path in sorted((p for p in root.rglob('*') if p.is_dir()), reverse=True):
        sync_directory(path)
    sync_directory(root)


def publish(local, destination, staging_root, execution):
    """Publish immutable bytes; retain local/staged files after any failure."""
    local, destination, staging_root = map(lambda p: Path(p).resolve(),
                                         (local, destination, staging_root))
    preflight(local,destination,staging_root)
    if not local.is_dir() or (local/'publication.json').exists():
        raise ValueError('missing local attempt or reserved publication filename')
    completed = (local/'completed.json').exists()
    if completed and (execution['returncode'] != 0 or (local/'failed.json').exists()):
        raise ValueError('contradictory completed attempt')
    if execution['returncode'] == 0 and not completed:
        raise ValueError('successful process without a completed attempt')
    with (destination.parent/'.publication.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if destination.exists():
            raise FileExistsError(destination)
        hashes = tree_hashes(local)
        staged = staging_root/str(uuid.uuid4())
        shutil.copytree(local, staged, symlinks=True)
        if tree_hashes(staged) != hashes or tree_hashes(local) != hashes:
            raise ValueError('publication bytes changed during copy')
        P.save(staged/'publication.json', dict(schema='gspa.local-publication/1',
            status='verified' if completed else 'failed_execution',
            local_path=str(local), published_path=str(destination),
            artifacts=hashes, execution=execution,
            interpretation='Byte-checked storage publication; does not establish biological accuracy.'))
        sync_tree(staged)
        # Staging is outside the glob used by the collection finalizer. Completion
        # becomes visible only with the complete, checked attempt directory.
        if destination.exists():
            raise FileExistsError(destination)
        staged.rename(destination)
        sync_directory(destination.parent)
        sync_directory(staging_root)
    return destination/'publication.json'


def command(out, index, attempt, heap, timeout, scratch, local_parent, container, launcher):
    logical_parent = out/'attempts'/str(index)
    for path in (scratch, local_parent, logical_parent):
        if any(c in str(path) for c in ':,\n'):
            raise ValueError('unsupported bind-path delimiter')
    return [str(launcher), 'exec', '--no-home', '--no-mount', 'tmp', '--pwd', str(out.parent),
        '--env', 'TMPDIR='+str(scratch),
        '--bind', '/usr,/lib,/lib64,/bin,/etc,/data', '--bind', str(scratch)+':'+str(scratch),
        '--bind', str(local_parent)+':'+str(logical_parent), str(container), sys.executable,
        '-c', MOUNT_GUARD, str(local_parent), str(logical_parent),
        str(Path(P.__file__).resolve()), 'run', '--out', str(out), '--index', str(index),
        '--attempt', attempt, '--heap', heap, '--timeout', str(timeout), '--scratch', str(scratch)]


# External launch (e.g. Shaheen): the wrapper itself runs inside a container,
# where the launcher cannot be nested. The wrapper prepares and records the
# engine command; the host executes that recorded command verbatim; the wrapper
# then publishes with the recorded return code. Publication rules are unchanged.

def image_identity(root):
    """Hash a container image tree; symlinks are recorded by target, never followed."""
    result = {}
    for path in sorted(root.rglob('*')):
        relative = str(path.relative_to(root))
        if path.is_symlink():
            result[relative] = 'symlink:'+os.readlink(path)
        elif path.is_dir():
            continue
        elif stat.S_ISREG(path.stat().st_mode):
            result[relative] = P.C.V.sha256_file(str(path))
        else:
            raise ValueError(f'non-regular container image artifact: {path}')
    return result


def load_prefix(path):
    """Host launcher argv up to, excluding, the container image."""
    prefix = json.loads(Path(path).read_text())
    if (not isinstance(prefix, list) or not prefix or
            not all(isinstance(t, str) and t and '\n' not in t for t in prefix)):
        raise ValueError('launch prefix must be a nonempty JSON list of strings')
    if not Path(prefix[0]).is_absolute():
        raise ValueError('launch prefix must start with an absolute launcher path')
    return prefix


def external_command(prefix, out, index, attempt, heap, timeout, scratch, local_parent, container):
    """Same bind contract and engine arguments as command(), after a host-built prefix."""
    logical_parent = out/'attempts'/str(index)
    for path in (scratch, local_parent, logical_parent):
        if any(c in str(path) for c in ':,\n'):
            raise ValueError('unsupported bind-path delimiter')
    return [*prefix, '--pwd', str(out.parent),
        '--env', 'TMPDIR='+str(scratch), '--bind', str(scratch)+':'+str(scratch),
        '--bind', str(local_parent)+':'+str(logical_parent), str(container), sys.executable,
        '-c', MOUNT_GUARD, str(local_parent), str(logical_parent),
        str(Path(P.__file__).resolve()), 'run', '--out', str(out), '--index', str(index),
        '--attempt', attempt, '--heap', heap, '--timeout', str(timeout), '--scratch', str(scratch)]


def locator_path(out, index, attempt):
    return out.parent/'local-executions'/f'{index}-{attempt}.json'


def prepare_external(out, index, attempt, heap, timeout, scratch, container, prefix_path):
    """Check, reserve and record one cell; return the durable locator. Launches nothing."""
    out, scratch, container = (p.resolve() for p in (out, scratch, container))
    prefix = load_prefix(prefix_path)
    logical_parent, destination, local_parent = checked_attempt(out, index, attempt, scratch)
    launcher = Path(prefix[0])
    if not launcher.is_file() or not container.is_dir():
        raise ValueError('a launcher visible to the wrapper and a container root are required')
    cmd = external_command(prefix,out,index,attempt,heap,timeout,scratch,local_parent,container)
    execution = dict(command=cmd, node=socket.gethostname(),
        registration_sha256=P.C.V.sha256_file(str(out/'registration.json')),
        runtime_files=P.C.V.file_hashes([Path(__file__),launcher]),
        container_files=image_identity(container), local_parent=str(local_parent),
        logical_parent=str(logical_parent), launch='external')
    locator = locator_path(out, index, attempt)
    locator.parent.mkdir(parents=True, exist_ok=True)
    with locator.open('x') as stream:
        json.dump(execution,stream,indent=2);stream.write('\n')
    return locator


def publish_external(out, index, attempt, scratch, container, returncode):
    """Bind the host return code to the recorded command, then publish as run() does."""
    out, scratch, container = (p.resolve() for p in (out, scratch, container))
    if isinstance(returncode, bool) or not isinstance(returncode, int):
        raise ValueError('integer return code required')
    locator = locator_path(out, index, attempt)
    execution = json.loads(locator.read_text())
    local_parent = scratch/('price-output-'+str(index)+'-'+attempt)
    logical_parent = out/'attempts'/str(index)
    if (execution.get('launch') != 'external' or execution.get('local_parent') != str(local_parent) or
            execution.get('logical_parent') != str(logical_parent)):
        raise ValueError('execution locator does not describe this external attempt')
    execution['returncode'] = returncode
    with locator.with_suffix('.result.json').open('x') as stream:
        json.dump(dict(returncode=returncode,
            locator_sha256=P.C.V.sha256_file(str(locator))),stream,indent=2)
        stream.write('\n')
    if (image_identity(container) != execution['container_files'] or
            P.C.V.file_hashes(execution['runtime_files']) != execution['runtime_files']):
        raise ValueError('publication runtime changed')
    publication = publish(local_parent/attempt,logical_parent/attempt,out.parent/'publication-staging',execution)
    print(json.dumps(dict(publication=str(publication),returncode=returncode)),flush=True)
    if returncode:
        raise RuntimeError(f'production failed ({returncode}); failed attempt retained')
    shutil.rmtree(local_parent)


def checked_attempt(out, index, attempt, scratch):
    """Validate a registered cell and reserve its private output parent."""
    if not re.fullmatch(r'[A-Za-z0-9_.-]+', attempt) or attempt in ('.', '..'):
        raise ValueError('invalid attempt identifier')
    plan = json.loads((out/'plan.json').read_text())
    registration = json.loads((out/'registration.json').read_text())
    if type(index) is not int or not 0 <= index < len(plan['cells']):
        raise ValueError('unknown collection index')
    if P.C.V.sha256_file(str(out/'plan.json')) != registration['plan_sha256']:
        raise ValueError('production plan changed')
    P.check(registration['inputs'])
    if any(os.environ.get(k) for k in ('JAVA_TOOL_OPTIONS', 'JDK_JAVA_OPTIONS', '_JAVA_OPTIONS')):
        raise ValueError('explicit Java options required without environment overrides')
    logical_parent = out/'attempts'/str(index)
    destination = logical_parent/attempt
    local_parent = scratch/('price-output-'+str(index)+'-'+attempt)
    preflight(local_parent/attempt,destination,out.parent/'publication-staging')
    local_parent.mkdir(parents=True, exist_ok=False)
    return logical_parent, destination, local_parent


def run(out, index, attempt, heap, timeout, scratch, container):
    out, scratch, container = (p.resolve() for p in (out, scratch, container))
    logical_parent, destination, local_parent = checked_attempt(out, index, attempt, scratch)
    launcher = Path(shutil.which('singularity') or '')
    if not launcher.is_file() or not container.is_dir():
        raise ValueError('Singularity and a container root are required')
    cmd = command(out,index,attempt,heap,timeout,scratch,local_parent,container,launcher)
    execution = dict(command=cmd, node=socket.gethostname(),
        registration_sha256=P.C.V.sha256_file(str(out/'registration.json')),
        runtime_files=P.C.V.file_hashes([Path(__file__),launcher]),
        container_files=tree_hashes(container), local_parent=str(local_parent),
        logical_parent=str(logical_parent))
    # This durable locator also survives a killed batch or failed publication.
    locator = out.parent/'local-executions'/f'{index}-{attempt}.json'
    locator.parent.mkdir(parents=True, exist_ok=True)
    with locator.open('x') as stream:
        json.dump(execution,stream,indent=2);stream.write('\n')
    result = subprocess.run(cmd)
    execution['returncode'] = result.returncode
    with locator.with_suffix('.result.json').open('x') as stream:
        json.dump(dict(returncode=result.returncode,
            locator_sha256=P.C.V.sha256_file(str(locator))),stream,indent=2)
        stream.write('\n')
    if (tree_hashes(container) != execution['container_files'] or
            P.C.V.file_hashes(execution['runtime_files']) != execution['runtime_files']):
        raise ValueError('publication runtime changed')
    local = local_parent/attempt
    publication = publish(local,destination,out.parent/'publication-staging',execution)
    print(json.dumps(dict(publication=str(publication),returncode=result.returncode)),flush=True)
    if result.returncode:
        raise RuntimeError(f'production failed ({result.returncode}); failed attempt retained')
    # Only successfully published private output is removed. Any failure leaves
    # its local bytes and durable locator available for recovery.
    shutil.rmtree(local_parent)


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('out','scratch','container'):
        parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--index',type=int,required=True)
    parser.add_argument('--attempt',required=True)
    parser.add_argument('--heap',required=True)
    parser.add_argument('--timeout',type=int,required=True)
    parser.add_argument('--launch-prefix',type=Path,
        help='External launch only: JSON list with the host launcher argv before the container image.')
    parser.add_argument('--phase',choices=['prepare','publish'],
        help='External launch only: prepare records the engine command; publish binds --returncode.')
    parser.add_argument('--returncode',type=int,help='External publish only: host return code of the recorded command.')
    args=parser.parse_args()
    if args.launch_prefix is None and args.phase is None and args.returncode is None:
        run(args.out,args.index,args.attempt,args.heap,args.timeout,args.scratch,args.container)
    elif args.phase=='prepare' and args.launch_prefix is not None and args.returncode is None:
        print(prepare_external(args.out,args.index,args.attempt,args.heap,args.timeout,args.scratch,
                               args.container,args.launch_prefix),flush=True)
    elif args.phase=='publish' and args.returncode is not None and args.launch_prefix is None:
        publish_external(args.out,args.index,args.attempt,args.scratch,args.container,args.returncode)
    else:
        parser.error('external launch: --phase prepare --launch-prefix FILE, or --phase publish --returncode N')
