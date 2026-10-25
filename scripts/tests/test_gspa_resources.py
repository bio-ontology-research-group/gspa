"""Offline integration tests for downloads, preparation and registry coverage."""
import contextlib
import hashlib
import http.server
import importlib.util
import io
import json
from pathlib import Path
import tarfile
import tempfile
import threading
import unittest
from unittest.mock import patch
import zipfile

spec = importlib.util.spec_from_file_location('resources', Path(__file__).parents[1] / 'gspa_resources.py')
r = importlib.util.module_from_spec(spec)
spec.loader.exec_module(r)


class Handler(http.server.BaseHTTPRequestHandler):
    payload = b'fixture data\n'
    calls = 0

    def do_GET(self):
        type(self).calls += 1
        if self.path == '/missing':
            self.send_error(404)
            return
        if self.path == '/html':
            self.send_response(200)
            self.send_header('Content-Type', 'text/html')
            self.end_headers()
            self.wfile.write(b'<html>login</html>')
            return
        self.send_response(200)
        self.send_header('Content-Length', str(len(self.payload) + (10 if self.path == '/truncated' else 0)))
        self.end_headers()
        self.wfile.write(self.payload)

    def log_message(self, *args):
        pass


class ResourceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.server = http.server.ThreadingHTTPServer(('127.0.0.1', 0), Handler)
        cls.thread = threading.Thread(target=cls.server.serve_forever, daemon=True)
        cls.thread.start()
        cls.url = f'http://127.0.0.1:{cls.server.server_port}'

    @classmethod
    def tearDownClass(cls):
        cls.server.shutdown()
        cls.server.server_close()
        cls.thread.join()

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.resource = dict(id='fixture', source=self.url, method='files',
            files=[dict(url=self.url+'/file', path='data.txt', sha256=hashlib.sha256(Handler.payload).hexdigest())],
            expected=['data.txt'], config={'quality.goOwlFile': '{target}/data.txt'})

    def test_download_reuse_and_corruption_repair(self):
        target = r.install(self.resource, self.root)
        calls = Handler.calls
        r.install(self.resource, self.root)
        self.assertEqual(calls, Handler.calls)
        (target/'data.txt').write_text('bad')
        r.install(self.resource, self.root)
        self.assertEqual((target/'data.txt').read_bytes(), Handler.payload)
        receipt = json.loads((target/'.gspa-resource.json').read_text())
        self.assertEqual(receipt['files']['data.txt']['sha256'], hashlib.sha256(Handler.payload).hexdigest())

    def test_bad_checksum_preserves_previous_installation(self):
        target = r.install(self.resource, self.root)
        self.resource['files'][0]['sha256'] = '0'*64
        with self.assertRaisesRegex(ValueError, 'SHA-256'):
            r.install(self.resource, self.root, force=True)
        self.assertEqual((target/'data.txt').read_bytes(), Handler.payload)
        self.assertTrue((target/'.gspa-resource.json').is_file())
        self.assertFalse(list(self.root.glob('.fixture-*')))

    def test_cache_without_upstream_digest_is_checked_before_reuse(self):
        target = self.root/'cache'
        r.fetch(self.url+'/file', target)
        calls = Handler.calls
        r.fetch(self.url+'/file', target)
        self.assertEqual(Handler.calls, calls)
        target.write_bytes(b'corrupted cache')
        r.fetch(self.url+'/file', target)
        self.assertEqual(target.read_bytes(), Handler.payload)
        self.assertEqual(Handler.calls, calls+1)

    def test_truncated_html_and_missing_files_are_not_promoted(self):
        for route in ['/truncated', '/html', '/missing']:
            with self.subTest(route=route), self.assertRaises((OSError, ValueError)):
                r.fetch(self.url+route, self.root/'file', attempts=1)
            self.assertFalse((self.root/'file').exists())
            self.assertFalse((self.root/'file.part').exists())

    def test_existing_unmanaged_directory_is_never_overwritten(self):
        (self.root/'fixture').mkdir()
        (self.root/'fixture/user-data').write_text('keep')
        with self.assertRaisesRegex(ValueError, 'without a GSPA receipt'):
            r.install(self.resource, self.root, force=True)
        self.assertEqual((self.root/'fixture/user-data').read_text(), 'keep')

    def test_safe_tar_and_zip_extraction(self):
        archive = self.root/'ok.tar.gz'
        with tarfile.open(archive, 'w:gz') as out:
            item = tarfile.TarInfo('release/data.txt')
            item.size = 4
            out.addfile(item, io.BytesIO(b'data'))
        r.unpack(archive, self.root/'extracted', 'tar', strip=1)
        self.assertEqual((self.root/'extracted/data.txt').read_text(), 'data')
        for name in ['../outside', '/absolute', 'x/../../outside', 'C:\\outside']:
            archive = self.root/'bad.zip'
            with zipfile.ZipFile(archive, 'w') as out:
                out.writestr(name, 'bad')
            with self.assertRaises(ValueError):
                r.unpack(archive, self.root/'bad', 'zip')

    def test_archive_links_are_rejected(self):
        archive = self.root/'link.tar'
        with tarfile.open(archive, 'w') as out:
            item = tarfile.TarInfo('link')
            item.type = tarfile.SYMTYPE
            item.linkname = '/tmp'
            out.addfile(item)
        with self.assertRaisesRegex(ValueError, 'links/devices'):
            r.unpack(archive, self.root/'bad', 'tar')

    def test_relative_links_inside_archive_are_supported(self):
        archive = self.root/'link.tar'
        with tarfile.open(archive, 'w') as out:
            link = tarfile.TarInfo('release/alias')
            link.type = tarfile.SYMTYPE
            link.linkname = 'data.txt'
            out.addfile(link)
            item = tarfile.TarInfo('release/data.txt')
            item.size = 4
            out.addfile(item, io.BytesIO(b'data'))
        r.unpack(archive, self.root/'safe', 'tar', strip=1)
        self.assertEqual((self.root/'safe/alias').read_text(), 'data')

    def test_wrong_archive_layout_is_not_marked_installed(self):
        self.resource['expected'] = ['weights.pt']
        with self.assertRaisesRegex(ValueError, 'expected output'):
            r.install(self.resource, self.root)
        self.assertFalse((self.root/'fixture').exists())

    def test_dependency_order_cycles_and_unknown_ids(self):
        entries = {'a': {'dependencies': ['b']}, 'b': {}}
        self.assertEqual(r.resolve(['a','b'], entries), ['b','a'])
        entries['b']['dependencies'] = ['a']
        with self.assertRaisesRegex(ValueError, 'cycle'):
            r.resolve(['a'], entries)
        with self.assertRaisesRegex(ValueError, 'Unknown resource'):
            r.resolve(['typo'], entries)

    def test_bulk_failures_do_not_hide_successful_downloads(self):
        entries = {'fixture': self.resource, 'manual': dict(id='manual', method='manual',
            source=self.url, notes='Unreleased model')}
        with patch.object(r, 'catalog', return_value=({'groups': {}}, entries)):
            self.assertEqual(r.main(['download','--all','--dir',str(self.root)]), 1)
        self.assertTrue((self.root/'fixture/data.txt').exists())
        self.assertFalse((self.root/'manual').exists())
        paths = json.loads((self.root/'gspa-paths.yaml').read_text())
        self.assertEqual(paths['quality']['goOwlFile'], str(self.root/'fixture/data.txt'))
        self.assertFalse((self.root/'.download.lock').exists())

    def test_dry_run_is_offline_and_does_not_create_destination(self):
        dest = self.root/'new'
        with patch.object(r.urllib.request, 'urlopen', side_effect=AssertionError('network')):
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(r.main(['download','--all','--dry-run','--dir',str(dest)]), 0)
        self.assertFalse(dest.exists())

    def test_custom_bundle_download_is_validated(self):
        self.resource.update(method='manual', custom_file={'path':'data.txt'}, notes='No release')
        self.resource.pop('files')
        with patch.object(r, 'catalog', return_value=({'groups':{}}, {'fixture':self.resource})):
            self.assertEqual(r.main(['download','fixture','--source-url',self.url+'/file',
                '--sha256',hashlib.sha256(Handler.payload).hexdigest(),'--dir',str(self.root)]), 0)
        self.assertEqual((self.root/'fixture/data.txt').read_bytes(), Handler.payload)

    def test_swissprot_headers_retain_go_ec_and_sequence(self):
        (self.root/'uniprot_sprot.dat').write_text('''ID   EXAMPLE
AC   P12345; Q00000;
DE   RecName: Full=Enzyme; EC=1.2.3.4;
DR   GO; GO:0003824; F:catalytic activity;
SQ   SEQUENCE   6 AA;
     MA AA CC
//
''')
        r.prepare_swissprot(self.root)
        self.assertEqual((self.root/'uniprot_sprot.fasta').read_text(), '>P12345 EC:1.2.3.4 GO:0003824\nMAAACC\n')

    def test_catalog_groups_models_and_bundled_bytes_are_complete(self):
        data, entries = r.catalog()
        for members in data['groups'].values():
            r.resolve(members, entries)
        shipped = {p.name for p in (r.REPO/'deepgo-plusplus/models').glob('*.json')}
        covered = set()
        for resource in entries.values():
            self.assertIn(resource['method'], ['files','bundled','command','manual'])
            self.assertTrue(resource.get('expected'))
            for item in resource.get('bundled', []):
                content = r.bundled_bytes(item['source'])
                self.assertEqual(hashlib.sha256(content).hexdigest(), item['sha256'])
                if item['source'].startswith('bundled/models/'):
                    covered.add(item['path'])
            for item in resource.get('files', []):
                self.assertTrue(item['url'].startswith('https://'))
                self.assertNotIn('REPLACE', json.dumps(item))
        self.assertEqual(shipped, covered)

    def test_command_failures_do_not_write_success_receipts(self):
        self.resource.update(method='command', commands=[['gspa-nonexistent-test-tool','download']])
        with self.assertRaisesRegex(ValueError, 'Requires installed command'):
            r.install(self.resource, self.root)
        self.assertFalse((self.root/'fixture').exists())


if __name__ == '__main__':
    unittest.main()
