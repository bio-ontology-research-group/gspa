"""Receipt validity and bounded retention, without inference or biological claims."""
import importlib.util
import io
import json
import os
from pathlib import Path
import sys
import tempfile
import tracemalloc
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import verified_batch as V
import v3_runner as R


@unittest.skipUnless(importlib.util.find_spec('ijson'), 'optional streaming parser not installed')
class StreamingLedgerTests(unittest.TestCase):
    def test_parser_binding_includes_bundled_native_library(self):
        import ijson
        bundled=list((Path(ijson.__file__).parent.parent/'ijson.libs').glob('*'))
        self.assertTrue(set(p for p in bundled if p.is_file()) <= set(V.streaming_parser_sources()))

    def parse(self, value):
        return V.stream_ledger_metadata(io.BytesIO(value.encode()))

    def test_fields_after_changes_and_nested_lookalikes(self):
        data = dict(schema_version=2, genome_id='μ', stages=[
            dict(changes=[{'name': 'fake', 'status': 'failed', 'nested': [None, True, 1.2]}],
                 name='evidence_integration', status='executed', model_sha256='abc')],
            final_changes=[{'genome_id': 'spoof'}], **{'stages.item.status': 'failed'})
        expected = dict(schema_version=2, genome_id='μ', stages=[
            dict(name='evidence_integration', status='executed', model_sha256='abc')])
        self.assertEqual(self.parse(json.dumps(data)), expected)

    def test_malformed_skipped_payload_truncation_and_trailing_data_fail(self):
        good = '{"schema_version":2,"genome_id":"g","stages":[],"final_changes":[{"x":1}]}'
        for text in [good[:-1], good+'{}', good+' garbage',
                     good.replace('"x":1', '"x":NaN'), good.replace('"x":1', '"x":'),
                     good.replace('"x":1', '"x":"bad\\x"'),
                     good.replace('"x":1', '"x":01')]:
            with self.subTest(text=text), self.assertRaises(ValueError): self.parse(text)

    def test_duplicate_metadata_and_invalid_record_shapes_fail(self):
        for text in ['[]', '{"stages":{}}', '{"stages":[null]}',
                     '{"genome_id":"g","genome_id":"h"}',
                     '{"stages":[],"stages":[]}',
                     '{"stages":[{"name":"a","name":"b"}]}',
                     '{"stages":[{"status":{}}]}']:
            with self.subTest(text=text), self.assertRaises(ValueError): self.parse(text)

    def test_small_fallback_and_stream_accept_the_same_corpus(self):
        corpus = ['{}', '{"stages":[]}', '{"stages":[{"status":"executed","changes":[{"x":1,"x":2}]}]}',
                  '{"stages":[{"status":"failed","status":"executed"}]}',
                  '{"stages":{}}', '{"genome_id":[]}', '{"stages":[[]]}',
                  '{"stages":[],"ignored":NaN}', '{"stages":[],"ignored":Infinity}',
                  '{"stages":[],"ignored":{}} trailing', '{"stages":[', '[]']
        with tempfile.TemporaryDirectory(dir=os.environ.get('TMPDIR', str(Path.home()/'tmp'))) as folder:
            path=Path(folder)/'ledger.json'
            for text in corpus:
                path.write_text(text)
                def outcome():
                    try: return V.ledger_metadata(path)
                    except ValueError: return 'rejected'
                streamed=outcome()
                with patch.dict(sys.modules, {'ijson': None}): fallback=outcome()
                with self.subTest(text=text): self.assertEqual(streamed,fallback)

    def test_ignored_unicode_and_large_exponents_have_matching_acceptance(self):
        with tempfile.TemporaryDirectory(dir=os.environ.get('TMPDIR', str(Path.home()/'tmp'))) as folder:
            path=Path(folder)/'ledger.json'
            for raw in (b'{"stages":[],"ignored":"\xff"}',
                        b'{"stages":[],"ignored":"\\ud800"}',
                        b'{"stages":[],"ignored":1e999}'):
                path.write_bytes(raw)
                def outcome():
                    try: return V.ledger_metadata(path)
                    except ValueError: return 'rejected'
                streamed=outcome()
                with patch.dict(sys.modules, {'ijson':None}): fallback=outcome()
                with self.subTest(raw=raw): self.assertEqual(streamed,fallback)

    def test_large_changes_do_not_accumulate_in_python_memory(self):
        with tempfile.TemporaryFile(dir=os.environ.get('TMPDIR', str(Path.home()/'tmp'))) as stream:
            stream.write(b'{"stages":[{"changes":[')
            row = b'{"protein_id":"p","after_probability":0.123,"notes":[true,null,"text"]}'
            for i in range(150000): stream.write((b',' if i else b'')+row)
            stream.write(b'],"name":"evidence_integration","status":"executed"}],"schema_version":2,"genome_id":"g"}')
            self.assertGreater(stream.tell(), V.LEDGER_MEMORY_LIMIT)
            stream.seek(0)
            tracemalloc.start()
            try:
                data = V.stream_ledger_metadata(stream)
                _, peak = tracemalloc.get_traced_memory()
            finally: tracemalloc.stop()
            self.assertEqual(data['stages'], [dict(name='evidence_integration',status='executed')])
            self.assertLess(peak, 5 * 1024 * 1024)

    def test_failed_duplicate_and_model_mismatch_cannot_be_certified(self):
        with tempfile.TemporaryDirectory(dir=os.environ.get('TMPDIR', str(Path.home()/'tmp'))) as folder:
            root=Path(folder); out=root/'out.tsv'; model=root/'model.json'
            out.write_text('protein_id\ttype\tfunction_id\taspect\tscore\np\tEC\t1.1.1.1\t\t.5\n')
            model.write_text('{}')
            spec=dict(outputs={'g': str(out)}, required_stages=['learned_integration'],
                      command=['java', '--learned-model', str(model)], inputs={str(model):'expected'})
            Path(str(out)+'.learned.json').write_text('{"model":{"sha256":"expected"}}')
            stage=dict(name='learned_integration',status='executed',model_sha256='expected')
            ledger=Path(str(out)+'.stages.json')
            for records, pattern in [([dict(stage,status='failed')], 'did not execute'),
                                     ([stage,stage], 'duplicate stages'),
                                     ([dict(stage,model_sha256='wrong')], 'model hash mismatch')]:
                ledger.write_text(json.dumps(dict(schema_version=2,genome_id='g',stages=records)))
                with self.subTest(pattern=pattern), self.assertRaisesRegex(ValueError,pattern): V.output_hashes(spec)
            ledger.write_text(json.dumps(dict(schema_version=2,genome_id='g',stages=[stage])))
            self.assertEqual(len(V.output_hashes(spec)),3)


class MissingDependencyTests(unittest.TestCase):
    def test_large_files_refuse_unbounded_fallback(self):
        with tempfile.TemporaryDirectory(dir=os.environ.get('TMPDIR', str(Path.home()/'tmp'))) as folder:
            path=Path(folder)/'ledger.json';path.write_text('{"schema_version":2}')
            with patch.dict(sys.modules, {'ijson':None}):
                self.assertEqual(V.ledger_metadata(path), {'schema_version':2})
                with patch.object(V.L, 'LEDGER_MEMORY_LIMIT', 1), self.assertRaisesRegex(ValueError, 'requires'):
                    V.ledger_metadata(path)


class GenericRunnerLedgerTests(unittest.TestCase):
    def test_generic_runner_uses_bounded_verifier_and_preserves_sidecar_hashes(self):
        with tempfile.TemporaryDirectory() as folder:
            out=Path(folder)/'g.tsv'; ledger=Path(str(out)+'.stages.json')
            ledger.write_text(json.dumps(dict(schema_version=2,genome_id='g',stages=[
                dict(name='reference_preparation',status='executed'),
                dict(name='evidence_integration',status='executed',changes=[dict(protein_id='p')])])) )
            sidecar=Path(str(out)+'.cohort.json');sidecar.write_text('{}')
            expected={str(p):R.hash_record(str(p)) for p in (ledger,sidecar)}
            # Whole-document json.load was the generic runner's memory regression.
            with patch.object(R.json,'load',side_effect=AssertionError('unbounded ledger read')):
                self.assertEqual(R.stage_artifacts('g',str(out),{'gspa':{}},{'args':[]}),expected)

    def test_generic_runner_rejects_corruption_after_stage_metadata(self):
        with tempfile.TemporaryDirectory() as folder:
            out=Path(folder)/'g.tsv';ledger=Path(str(out)+'.stages.json')
            stages=[dict(name=n,status='executed') for n in ('reference_preparation','evidence_integration')]
            valid=json.dumps(dict(schema_version=2,genome_id='g',stages=stages))
            for text in (valid+'{}', valid[:-1]+',"changes":[NaN]}', valid[:-1]):
                ledger.write_text(text)
                with self.subTest(text=text),self.assertRaises(ValueError):
                    R.stage_artifacts('g',str(out),{'gspa':{}},{'args':[]})


if __name__ == '__main__': unittest.main()
