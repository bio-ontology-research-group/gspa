#!/usr/bin/env python3
"""GSPA resource manager. Python 3.9+, standard library only; also shipped as a zipapp."""
import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
import time
import urllib.error
import urllib.request
import zipfile

REPO = Path(__file__).resolve().parents[1]


def bundled_bytes(name):
    if zipfile.is_zipfile(sys.argv[0]):
        with zipfile.ZipFile(sys.argv[0]) as archive:
            return archive.read(name)
    if name == 'catalog.json':
        return (REPO / 'gspa-core/src/main/resources/downloads/catalog.json').read_bytes()
    if name.startswith('bundled/models/'):
        return (REPO / 'deepgo-plusplus/models' / name.removeprefix('bundled/models/')).read_bytes()
    if name.startswith('bundled/sidecars/'):
        return (REPO / 'benchmark/neural' / name.removeprefix('bundled/sidecars/')).read_bytes()
    return (REPO / 'gspa-core/src/main/resources' / name.removeprefix('bundled/')).read_bytes()


def catalog():
    data = json.loads(bundled_bytes('catalog.json'))
    entries = {r['id']: r for r in data['resources']}
    if len(entries) != len(data['resources']):
        raise ValueError('Duplicate resource IDs')
    for key, resource in entries.items():
        if not re.fullmatch(r'[a-z0-9][a-z0-9.-]*', key):
            raise ValueError('Invalid resource ID: ' + key)
        for dependency in resource.get('dependencies', []):
            if dependency not in entries:
                raise ValueError(f'{key}: unknown dependency {dependency}')
    return data, entries


def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def safe_path(root, name):
    path = PurePosixPath(name)
    if path.is_absolute() or '..' in path.parts or '\\' in name or ':' in name:
        raise ValueError('Unsafe archive/output path: ' + name)
    dest = root.joinpath(*path.parts)
    if not dest.resolve().is_relative_to(root.resolve()):
        raise ValueError('Path escapes resource directory: ' + name)
    return dest


def unpack(archive, target, mode, strip=0):
    """Write files first; permit only relative archive links confined to the resource."""
    def destination(name):
        safe_path(target, name)  # validate before stripping components
        parts = PurePosixPath(name).parts[strip:]
        return safe_path(target, '/'.join(parts)) if parts else None

    if mode == 'gzip':
        target.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(archive, 'rb') as src, target.open('wb') as dst:
            shutil.copyfileobj(src, dst)
        return
    target.mkdir(parents=True, exist_ok=True)
    if mode == 'zip':
        with zipfile.ZipFile(archive) as src:
            for member in src.infolist():
                if (member.external_attr >> 16) & 0o170000 == 0o120000:
                    raise ValueError('Archive symlinks are not supported: ' + member.filename)
                dest = destination(member.filename)
                if dest is None:
                    continue
                if member.is_dir():
                    dest.mkdir(parents=True, exist_ok=True)
                else:
                    dest.parent.mkdir(parents=True, exist_ok=True)
                    with src.open(member) as inp, dest.open('wb') as out:
                        shutil.copyfileobj(inp, out)
    elif mode == 'tar':
        links = []
        with tarfile.open(archive) as src:
            for member in src:
                dest = destination(member.name)
                if member.issym() or member.islnk():
                    if dest is not None:
                        links.append((dest, member))
                    continue
                if not (member.isfile() or member.isdir()):
                    raise ValueError('Archive links/devices are not supported: ' + member.name)
                if dest is None:
                    continue
                if member.isdir():
                    dest.mkdir(parents=True, exist_ok=True)
                else:
                    dest.parent.mkdir(parents=True, exist_ok=True)
                    with src.extractfile(member) as inp, dest.open('wb') as out:
                        shutil.copyfileobj(inp, out)
                    dest.chmod(member.mode & 0o755)
        # No regular file is ever written through an archive-provided symlink.
        for dest, member in links:
            if PurePosixPath(member.linkname).is_absolute() or '\\' in member.linkname:
                raise ValueError('Archive links/devices escape the resource: ' + member.name)
            link_target = (dest.parent / member.linkname) if member.issym() else destination(member.linkname)
            if link_target is None or not link_target.resolve().is_relative_to(target.resolve()):
                raise ValueError('Archive links/devices escape the resource: ' + member.name)
            dest.parent.mkdir(parents=True, exist_ok=True)
            if member.issym():
                dest.symlink_to(member.linkname)
            else:
                os.link(link_target, dest)
    else:
        raise ValueError('Unknown archive type: ' + mode)


def fetch(url, output, expected=None, attempts=3):
    """Restart interrupted transfers, validate before promotion; never cache an error page."""
    if not url.startswith(('https://', 'http://')):
        raise ValueError('Download URL must use HTTP(S)')
    output.parent.mkdir(parents=True, exist_ok=True)
    cache_digest = output.with_name(output.name + '.sha256')
    known_digest = expected or (cache_digest.read_text().strip() if cache_digest.is_file() else None)
    if output.is_file() and known_digest and sha256(output) == known_digest:
        return
    partial = output.with_name(output.name + '.part')
    for attempt in range(attempts):
        try:
            request = urllib.request.Request(url, headers={'User-Agent': 'GSPA-resource-manager/1'})
            with urllib.request.urlopen(request, timeout=60) as response, partial.open('wb') as dest:
                if response.headers.get_content_type() == 'text/html':
                    raise ValueError('Server returned HTML instead of a resource: ' + url)
                length = response.headers.get('Content-Length')
                copied = 0
                while True:
                    chunk = response.read(1024 * 1024)
                    if not chunk:
                        break
                    dest.write(chunk)
                    copied += len(chunk)
                if copied == 0 or (length is not None and copied != int(length)):
                    raise OSError('Empty or truncated download: ' + url)
            actual_digest = sha256(partial)
            if expected and actual_digest != expected:
                raise ValueError('SHA-256 mismatch: ' + url)
            os.replace(partial, output)
            cache_digest.write_text(actual_digest + '\n')
            return
        except (OSError, urllib.error.URLError) as error:
            if isinstance(error, urllib.error.HTTPError) and 400 <= error.code < 500 and error.code not in (408, 429):
                raise
            if attempt + 1 == attempts:
                raise
            time.sleep(attempt + 1)
        finally:
            partial.unlink(missing_ok=True)


def file_inventory(directory):
    return {str(p.relative_to(directory)): {'bytes': p.stat().st_size, 'sha256': sha256(p)}
            for p in sorted(directory.rglob('*')) if p.is_file() and p.name != '.gspa-resource.json'}


def installed(directory, fingerprint):
    receipt = directory / '.gspa-resource.json'
    if not receipt.is_file():
        return False
    try:
        data = json.loads(receipt.read_text())
        return data['fingerprint'] == fingerprint and bool(data['files']) and all(
            safe_path(directory, name).is_file() and safe_path(directory, name).stat().st_size == record['bytes']
            and sha256(safe_path(directory, name)) == record['sha256']
            for name, record in data['files'].items())
    except (ValueError, OSError, KeyError):
        return False


def expand(value, root, target):
    return value.replace('{root}', str(root)).replace('{target}', str(target)).replace('{python}', sys.executable)


def prepare_swissprot(directory):
    """Preserve the GO/EC evidence that GSPA reads from reference hit titles."""
    accession, terms, sequence, in_sequence = None, set(), [], False
    count = 0
    with (directory / 'uniprot_sprot.dat').open() as inp, (directory / 'uniprot_sprot.fasta').open('w') as out:
        for line in inp:
            if line.startswith('AC   ') and accession is None:
                accession = line[5:].split(';')[0].strip()
            elif line.startswith('DR   GO;'):
                terms.update(re.findall(r'GO:\d{7}', line))
            elif line.startswith('DE   '):
                terms.update('EC:' + ec for ec in re.findall(r'EC=([\d-]+\.[\d-]+\.[\d-]+\.[\d-]+)', line))
            elif line.startswith('SQ   '):
                in_sequence = True
            elif line.startswith('//'):
                if not accession or not sequence:
                    raise ValueError('Malformed Swiss-Prot record: missing accession or sequence')
                out.write('>' + accession + ' ' + ' '.join(sorted(terms)) + '\n' + ''.join(sequence) + '\n')
                count += 1
                accession, terms, sequence, in_sequence = None, set(), [], False
            elif in_sequence:
                sequence.append(''.join(line.split()))
    if not count or accession is not None:
        raise ValueError('Empty or truncated Swiss-Prot DAT file')


def install(resource, root, force=False):
    name = resource['id']
    target = root / name
    fingerprint = hashlib.sha256(json.dumps(resource, sort_keys=True).encode()).hexdigest()
    if installed(target, fingerprint) and not force:
        print(f'{name}: verified existing download', flush=True)
        return target
    if resource['method'] == 'manual':
        raise ValueError(resource['notes'] + ' Source: ' + resource['source'])
    for command in resource.get('commands', []):
        executable = expand(command[0], root, target)
        if not shutil.which(executable):
            raise ValueError(f'Requires installed command: {executable}. See {resource["source"]}')
    if target.exists() and not (target / '.gspa-resource.json').is_file():
        raise ValueError(f'{target} exists without a GSPA receipt; choose another --dir (existing data is not overwritten).')
    work = Path(tempfile.mkdtemp(prefix='.' + name + '-', dir=root))
    try:
        for item in resource.get('bundled', []):
            dest = safe_path(work, item['path'])
            dest.parent.mkdir(parents=True, exist_ok=True)
            content = bundled_bytes(item['source'])
            if item.get('sha256') and hashlib.sha256(content).hexdigest() != item['sha256']:
                raise ValueError('Bundled resource checksum mismatch: ' + item['source'])
            dest.write_bytes(content)
        for item in resource.get('files', []):
            # Completed transfers survive a failed extraction/command and are reusable.
            key = hashlib.sha256((item['url'] + (item.get('sha256') or '')).encode()).hexdigest()
            cached = root / '.downloads' / key
            print(f'{name}: downloading {item["url"]}', flush=True)
            if force:
                cached.unlink(missing_ok=True)
            fetch(item['url'], cached, item.get('sha256'))
            dest = safe_path(work, item['path'])
            if item.get('unpack'):
                unpack(cached, dest, item['unpack'], item.get('strip', 0))
            else:
                dest.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(cached, dest)
        if resource.get('prepare') == 'swissprot':
            prepare_swissprot(work)
        for command in resource.get('commands', []):
            argv = [expand(arg, root, work) for arg in command]
            print(f'{name}: running {argv}', flush=True)
            subprocess.run(argv, cwd=work, check=True)
        for item in resource.get('archives', []):
            archive_target = safe_path(work, item['target'])
            unpack(safe_path(work, item['path']), archive_target, item['unpack'], item.get('strip', 0))
            # Some upstream bundles wrap the payload in an extra directory.
            if item.get('root_file') and not (archive_target / item['root_file']).is_file():
                matches = list(archive_target.rglob(item['root_file']))
                if len(matches) != 1:
                    raise ValueError('Cannot locate unique model bundle root: ' + item['root_file'])
                shutil.copytree(matches[0].parent, archive_target, dirs_exist_ok=True)
        for pattern in resource.get('expected', []):
            if not any(p.is_file() and p.stat().st_size > 0 for p in work.glob(pattern)):
                raise ValueError(f'{name}: expected output is missing/empty: {pattern}')
        inventory = file_inventory(work)
        if not inventory:
            raise ValueError(f'{name}: downloader produced no files')
        (work / '.gspa-resource.json').write_text(json.dumps({
            'id': name, 'version': resource.get('version'), 'source': resource['source'],
            'fingerprint': fingerprint, 'files': inventory,
            'note': 'Locally recorded SHA-256 values are integrity receipts, not upstream signatures.'
        }, indent=2) + '\n')
        # Do not remove a working copy until its replacement is complete.
        backup = root / ('.' + name + '-previous')
        if backup.exists():
            raise ValueError(f'Recovery directory exists: {backup}; inspect it before retrying.')
        if target.exists():
            target.rename(backup)
        try:
            work.rename(target)
        except BaseException:
            if backup.exists():
                backup.rename(target)
            raise
        if backup.exists():
            shutil.rmtree(backup)
        print(f'{name}: installed in {target}', flush=True)
        return target
    finally:
        if work.exists():
            shutil.rmtree(work)


def resolve(ids, entries):
    selected, visiting = [], set()
    def visit(name):
        if name not in entries:
            raise ValueError(f'Unknown resource {name!r}; run gspa resources.')
        if name in visiting:
            raise ValueError('Dependency cycle: ' + name)
        if name in selected:
            return
        visiting.add(name)
        for dependency in entries[name].get('dependencies', []):
            visit(dependency)
        visiting.remove(name)
        selected.append(name)
    for name in ids:
        visit(name)
    return selected


def write_paths(root, entries):
    paths, config = {}, {}
    for name, resource in entries.items():
        target = root / name
        receipt = target / '.gspa-resource.json'
        if not receipt.is_file():
            continue
        paths[name] = {'directory': str(target), 'usage': expand(resource.get('usage', ''), root, target)}
        for key, value in resource.get('config', {}).items():
            branch = config
            parts = key.split('.')
            for part in parts[:-1]:
                branch = branch.setdefault(part, {})
            branch[parts[-1]] = expand(value, root, target) if isinstance(value, str) else value
    # JSON is accepted by the YAML config loader; only paths are set, never enable flags.
    for name, value in [('paths.json', paths), ('gspa-paths.yaml', config)]:
        temp = root / ('.' + name + '.tmp')
        temp.write_text(json.dumps(value, indent=2) + '\n')
        temp.replace(root / name)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action', required=True)
    listing = sub.add_parser('list')
    listing.add_argument('id', nargs='?')
    listing.add_argument('--json', action='store_true')
    download = sub.add_parser('download')
    download.add_argument('ids', nargs='*')
    download.add_argument('--all', action='store_true')
    download.add_argument('--group', action='append', default=[])
    download.add_argument('--dir', type=Path, default=Path(os.getenv('GSPA_DATA_DIR', Path.home() / '.local/share/gspa')))
    download.add_argument('--dry-run', action='store_true')
    download.add_argument('--force', action='store_true')
    download.add_argument('--source-url', help='Override one-file resource URL, or supply a custom resource bundle')
    download.add_argument('--sha256', help='Expected SHA-256 for --source-url')
    args = parser.parse_args(argv)
    data, entries = catalog()
    if args.action == 'list':
        if args.id:
            if args.id not in entries:
                parser.error('Unknown resource: ' + args.id)
            print(json.dumps(entries[args.id], indent=2))
        elif args.json:
            print(json.dumps(data, indent=2))
        else:
            print('ID                                    KIND        METHOD    DESCRIPTION')
            for entry in entries.values():
                print(f'{entry["id"]:37} {entry["kind"]:11} {entry["method"]:9} {entry["description"]}')
            print('\nGroups: ' + ', '.join(data['groups']))
            print('Details: gspa resources ID; preview: gspa download --all --dry-run')
            print('manual = no verified public bundle or licensed/custom data; see resource details.')
        return 0
    ids = list(entries) if args.all else list(args.ids)
    for group in args.group:
        if group not in data['groups']:
            parser.error('Unknown group: ' + group)
        ids.extend(data['groups'][group])
    if not ids:
        parser.error('Choose resource IDs, --group NAME, or --all. Run gspa resources.')
    if args.sha256 and (not args.source_url or not re.fullmatch(r'[a-fA-F0-9]{64}', args.sha256)):
        parser.error('--sha256 requires --source-url and a 64-digit hex digest')
    if args.source_url:
        if len(ids) != 1 or args.all or args.group:
            parser.error('--source-url requires exactly one resource ID')
        entry = dict(entries[ids[0]]) if ids[0] in entries else None
        template = entry.get('custom_file') if entry else None
        if template is None and entry and len(entry.get('files', [])) == 1:
            template = entry['files'][0]
        if template is None:
            parser.error('--source-url needs a one-file resource or custom_file metadata; run gspa resources ID')
        entry['method'] = 'files'
        entry['files'] = [dict(template, url=args.source_url, sha256=args.sha256.lower() if args.sha256 else None)]
        entries[ids[0]] = entry
    selected = resolve(ids, entries)
    root = args.dir.expanduser().resolve()
    if args.dry_run:
        for name in selected:
            resource = entries[name]
            print(f'{name} [{resource["method"]}] -> {root / name}')
            for item in resource.get('files', []):
                print('  GET ' + item['url'])
            for command in resource.get('commands', []):
                print('  RUN ' + repr([expand(v, root, root / name) for v in command]))
            if resource.get('notes'):
                print('  ' + resource['notes'])
        return 0
    root.mkdir(parents=True, exist_ok=True)
    lock = root / '.download.lock'
    try:
        lock.mkdir()
    except FileExistsError:
        raise ValueError(f'Another download may be running: {lock}. Remove this directory only after confirming it is stale.')
    failures, successful = {}, set()
    try:
        for name in selected:
            resource = entries[name]
            try:
                failed_deps = [d for d in resource.get('dependencies', []) if d not in successful]
                if failed_deps:
                    raise ValueError('Unavailable dependencies: ' + ', '.join(failed_deps))
                install(resource, root, args.force)
                successful.add(name)
            except (ValueError, OSError, subprocess.CalledProcessError, tarfile.TarError, zipfile.BadZipFile) as error:
                failures[name] = str(error)
                print(f'{name}: FAILED: {error}', file=sys.stderr, flush=True)
        write_paths(root, entries)
        print(f'\nDownloaded/verified: {len(successful)}; unavailable/failed: {len(failures)}. Paths: {root / "paths.json"}')
        print(f'Annotation path configuration: {root / "gspa-paths.yaml"} (does not enable predictors).')
        for name, error in failures.items():
            print(f'  {name}: {error}', file=sys.stderr)
        return 1 if failures else 0
    finally:
        lock.rmdir()


if __name__ == '__main__':
    try:
        sys.exit(main())
    except (ValueError, OSError) as error:
        print('gspa resources: ' + str(error), file=sys.stderr)
        sys.exit(1)
