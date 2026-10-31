"""Check retained evidence and independently recompute native memberships."""
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import sys
import tarfile
import tempfile
import xml.etree.ElementTree as E
sys.dont_write_bytecode = True
from check_native_memberships import check

root = Path(__file__).resolve().parent

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

inventory = {}
for line in (root/'archive.sha256').read_text().splitlines():
    digest, name = line.split('  ', 1)
    assert name not in inventory
    inventory[name] = digest
actual = {str(p.relative_to(root)) for p in root.rglob('*') if p.is_file() and p.name != 'archive.sha256'}
assert actual == set(inventory)
for name, digest in inventory.items():
    assert sha(root/name) == digest, name
for version, count in (('v1', 801), ('v2', 802)):
    manifest = json.loads((root/f'build-source-files-{version}.json').read_text())
    assert len(manifest) == count
    prefix = 'home/leechuck/Public/software/gspa/'
    with tarfile.open(root/f'build-source-{version}.tar.gz') as tar:
        assert {m.name for m in tar.getmembers()} == {prefix+p for p in manifest}
        for name, digest in manifest.items():
            assert hashlib.sha256(tar.extractfile(prefix+name).read()).hexdigest() == digest, name
for label, expected in [('before-tests', (27,27,0,0)), ('focused-tests', (124,0,0,0)),
                        ('final-tests', (1502,0,0,6))]:
    values = dict(tests=0, failures=0, errors=0, skipped=0)
    for path in (root/'job-10439'/label).rglob('TEST-*.xml'):
        doc = E.parse(path).getroot()
        for key in values:
            values[key] += int(doc.get(key, '0'))
    assert tuple(values.values()) == expected, (label, values)
for job, code in ((10437, '1'), (10439, '1'), (10459, '0')):
    assert (root/f'job-{job}/exit-code.txt').read_text().strip() == code
cli = root/'job-10459/cli'
report = json.loads((cli/'verification.json').read_text())
for path in cli.rglob('*'):
    if path.is_file() and path.name != 'verification.json':
        assert sha(path) == report['outputs'][str(path.relative_to(cli))], path
assert report['cli_runs'] == 4 and report['invented_function_calls'] == 0
assert 'b8b28f522732f8873ff2f899d2f6f4119c4ced2ca187f6c79f4d1a33b7fd0d7d' in report['inputs'].values()
assert (cli/'direct.tsv').read_bytes() == (cli/'annotate_integrated/input_integrated.tsv').read_bytes()
scratch = Path.home()/'tmp'
scratch.mkdir(exist_ok=True)
with tempfile.TemporaryDirectory(prefix='gspa-antismash-verify-', dir=scratch) as directory:
    native = Path(directory)/'native.json'
    with gzip.open(root/'GCA_003598085.1.json.gz', 'rb') as source, native.open('wb') as target:
        shutil.copyfileobj(source, target)
    provenance = json.loads((root/'observed-fixture/provenance.json').read_text())
    assert sha(native) == provenance['source_sha256']
    subset = json.loads((root/'observed-fixture/records.json').read_text())
    full = json.loads(native.read_text())
    assert subset['records'] == [full['records'][i] for i in provenance['record_indices']]
    assert sha(root/'observed-fixture/records.json') == provenance['subset_sha256']
    result = check(native, cli)
    assert result == json.loads((root/'job-10459/independent-memberships.json').read_text())
result.update(archive_files=len(inventory), full_jvm_cases=1502, optional_skips=6)
print(json.dumps(result, indent=2))
