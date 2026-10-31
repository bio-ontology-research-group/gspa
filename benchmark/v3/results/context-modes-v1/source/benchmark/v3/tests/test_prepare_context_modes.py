import copy
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from prepare_context_modes import summarize, load_graph_ids


class ContextModePreparationTests(unittest.TestCase):
    def row(self, **overrides):
        return dict(pathway='P', rxn='G', dbhit='rxn1 rxn2', status='good_blast',
                    spont='FALSE', is_complex='FALSE', **overrides)

    def test_alternatives_and_unknown_mappings_remain_explicit(self):
        r = self.row()
        result = summarize([r, copy.deepcopy(r)], {'rxn1'})
        self.assertEqual(result['groups'], {('P','G'):{'rxn1','rxn2'}})
        self.assertEqual(result['plausible'], {'rxn1'})
        self.assertEqual(result['unknown'], {'rxn2'})

    def test_weak_or_complex_hits_cannot_establish_an_organism_reaction(self):
        for changes in ({'status':'bad_blast'}, {'status':'no_seq_data'}, {'is_complex':'TRUE'}):
            with self.subTest(changes=changes):
                r = self.row(); r.update(changes)
                self.assertEqual(summarize([r], {'rxn1','rxn2'})['plausible'], set())

    def test_spontaneous_reactions_and_empty_observation_are_distinct_from_missing(self):
        r = self.row(); r.update(status='spontaneous', spont='TRUE')
        self.assertEqual(summarize([r], {'rxn1'})['plausible'], {'rxn1'})
        self.assertEqual(summarize([], {'rxn1'})['groups'], {})
        del r['dbhit']
        with self.assertRaises(KeyError): summarize([r], {'rxn1'})

    def test_malformed_mapping_and_status_fail(self):
        for changes in ({'dbhit':'not-a-reaction'}, {'status':'typo'}, {'spont':''}, {'pathway':''}):
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                r = self.row(); r.update(changes)
                summarize([r], {'rxn1'})

    def test_empty_stoichiometry_rows_are_not_production_graph_nodes(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d)/'graph.tsv'
            path.write_text('id\tstoichiometry\tec_numbers\nrxn1\t-1:a;1:b\t1.1.1.1\nrxn2\t\t2.2.2.2\n')
            self.assertEqual(load_graph_ids(path), {'rxn1'})
