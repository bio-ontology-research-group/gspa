"""Identifier semantics: transfer is not a license to pick a split function."""
import os
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import ec_identity as E


def record(status='active',targets=()):
    return dict(status=status,targets=list(targets),description='fixture')


class EcIdentityTests(unittest.TestCase):
    def test_legacy_duplicates_and_non_string_targets_are_rejected(self):
        with tempfile.TemporaryDirectory(dir=os.environ.get('TMPDIR',str(Path.home()/'tmp'))) as folder:
            p=Path(folder)/'legacy.json'
            for value in ['{"aliases":{"1.1.1.1":"1.1.1.2","1.1.1.1":"1.1.1.3"}}',
                          '{"aliases":{"1.1.1.1":["1.1.1.2"]}}','{"aliases":[]}',
                          '{"aliases":{"unknown":"1.1.1.2"}}']:
                p.write_text(value)
                with self.subTest(value=value),self.assertRaises(ValueError):E.load_legacy(p)

    def test_build_checks_pins_and_refuses_overwriting_an_audit(self):
        with tempfile.TemporaryDirectory(dir=os.environ.get('TMPDIR',str(Path.home()/'tmp'))) as folder:
            root=Path(folder);source=root/'enzyme.dat';legacy=root/'legacy.json'
            source.write_text('ID   1.1.1.1\nDE   Test enzyme.\n//\n')
            legacy.write_text(json.dumps({'aliases':{'1.1.1.1':'1.1.1.1'}}))
            pins=dict(source_sha256=E.sha(source),legacy_sha256=E.sha(legacy))
            with self.assertRaisesRegex(ValueError,'hash mismatch'):
                E.build(source,legacy,root/'wrong',**dict(pins,source_sha256='wrong'))
            self.assertFalse((root/'wrong').exists())
            result=E.build(source,legacy,root/'audit',**pins)
            self.assertEqual(result['legacy_audit']['counts'],{'unique_active_equivalence':1})
            with self.assertRaises(FileExistsError):E.build(source,legacy,root/'audit',**pins)
    def test_single_transfer_chain_and_preliminary_identity(self):
        entries={'1.1.1.n1':record('transferred',['1.1.1.1']),
                 '1.1.1.1':record('transferred',['1.1.1.2']),'1.1.1.2':record(),
                 '2.1.1.n1':record()}
        result=E.resolve(entries,'EC:1.1.1.n1')
        self.assertEqual(result['canonical_id'],'1.1.1.2')
        self.assertEqual(result['path'],['1.1.1.n1','1.1.1.1','1.1.1.2'])
        self.assertFalse(result['preliminary'])
        self.assertTrue(E.resolve(entries,'2.1.1.n1')['preliminary'])

    def test_split_is_ambiguous_even_if_destinations_later_converge(self):
        entries={'1.1.1.5':record('transferred',['1.1.1.303','1.1.1.304']),
                 '1.1.1.303':record('transferred',['1.1.1.304']),'1.1.1.304':record()}
        result=E.resolve(entries,'1.1.1.5')
        self.assertEqual(result['status'],'ambiguous');self.assertIsNone(result['canonical_id'])
        self.assertEqual(E.audit(entries,{'1.1.1.5':'1.1.1.304'})['counts'],{'ambiguous':1})

    def test_transfer_to_deleted_unknown_and_cycle_never_create_equivalence(self):
        entries={'1.6.2.1':record('transferred',['1.6.99.3']), '1.6.99.3':record('deleted'),
                 '1.1.1.1':record('transferred',['1.1.1.2']), '1.1.1.2':record('transferred',['1.1.1.1'])}
        for ec,status in [('1.6.2.1','deleted'),('9.9.9.9','unknown'),('1.1.1.1','cycle')]:
            with self.subTest(ec=ec):
                result=E.resolve(entries,ec);self.assertEqual(result['status'],status)
                self.assertIsNone(result['canonical_id'])

    def test_invalid_branches_are_rejected_even_when_source_is_ambiguous(self):
        entries={'1.1.1.1':record('transferred',['1.1.1.2','1.1.1.3']), '1.1.1.2':record()}
        with self.assertRaisesRegex(ValueError,'missing'):E.validate_transfer_graph(entries)
        entries['1.1.1.3']=record('transferred',['1.1.1.1'])
        with self.assertRaisesRegex(ValueError,'cycle'):E.validate_transfer_graph(entries)
        entries['1.1.1.3']=record();E.validate_transfer_graph(entries)

    def test_distinct_current_reductases_are_not_equated(self):
        entries={'1.5.98.3':record(),'7.1.1.2':record(),'1.6.5.3':record('transferred',['7.1.1.2'])}
        self.assertEqual(E.audit(entries,{'1.6.5.3':'7.1.1.2'})['counts'],{'unique_active_equivalence':1})
        self.assertEqual(E.audit(entries,{'1.5.98.3':'7.1.1.2'})['counts'],{'different_active_function':1})

    def test_parser_continuation_and_malformed_records(self):
        with tempfile.TemporaryDirectory(dir=os.environ.get('TMPDIR',str(Path.home()/'tmp'))) as folder:
            path=Path(folder)/'enzyme.dat'
            valid='CC   header\n//\nID   1.1.1.5\nDE   Transferred entry: 1.1.1.303\nDE   and 1.1.1.304.\n//\n'
            path.write_text(valid);self.assertEqual(E.parse_enzyme(path)['1.1.1.5']['targets'],['1.1.1.303','1.1.1.304'])
            for data in [valid+valid,valid[:-3],valid.replace(' and ',' or '),
                         valid.replace('1.1.1.304','1.1.1.303'), 'DE   Missing ID.\n//\n', '']:
                path.write_text(data)
                with self.subTest(data=data),self.assertRaises(ValueError):E.parse_enzyme(path)


if __name__=='__main__':unittest.main()
