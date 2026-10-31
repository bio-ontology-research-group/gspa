import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import propose_negated_functions as P


OBO = '''format-version: 1.2

[Term]
id: GO:0000001
name: fructose-bisphosphate (FBP) aldolase activity
namespace: molecular_function
alt_id: GO:0001001
synonym: "fructose bisphosphate aldolase activity" EXACT []
synonym: "(FBP) aldolase activity" RELATED []
synonym: "FBP" EXACT []

[Term]
id: GO:0000002
name: glutaminase activity
namespace: molecular_function
synonym: "glutaminase" RELATED []

[Term]
id: GO:0000003
name: 2-isopropylmalate synthase activity
namespace: molecular_function
synonym: "2 isopropylmalate synthase" EXACT []

[Term]
id: GO:0000004
name: citramalate synthase activity
namespace: molecular_function
synonym: "citramalate synthase" EXACT []
synonym: "synthase" RELATED []

[Term]
id: GO:0016832
name: aldehyde-lyase activity
namespace: molecular_function
synonym: "aldolase activity" BROAD []

[Term]
id: GO:0097023
name: fructose 6-phosphate aldolase activity
namespace: molecular_function

[Term]
id: GO:0052749
name: F420-dependent glucose-6-phosphate dehydrogenase activity
namespace: molecular_function

[Term]
id: GO:0004345
name: glucose-6-phosphate dehydrogenase activity
namespace: molecular_function

[Term]
id: GO:0000005
name: ATPase activity
namespace: molecular_function
is_obsolete: true

[Term]
id: GO:0000006
name: ATPase activity
namespace: biological_process
'''

ENZYME = '''ID   3.6.1.3
DE   Adenosinetriphosphatase.
AN   ATPase.
//
ID   2.3.3.99
DE   Citramalate synthase.
//
ID   1.2.3.4
DE   Transferred entry: 9.9.9.9.
//
'''


def statement(accession, sentence, protein, *, go_mf=None, entry_ec=None):
    return {
        'kind': 'swissprot_cc', 'taxid': '1', 'accession': accession,
        'sentence': sentence, 'statement_class': 'activity_negative',
        'pubmed_ids': ['123'], 'joins': [{'genome_id': 'g1', 'protein_id': protein}],
        'entry_ec': entry_ec or [], 'entry_go_mf': go_mf or [],
        'review_status': 'pending', 'accepted_training_label': False,
    }


def jsonl(path, rows):
    path.write_text(''.join(json.dumps(row) + '\n' for row in rows))


class ProposeNegatedFunctionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory()
        cls.root = Path(cls.temp.name)
        cls.obo = cls.root / 'go-basic.obo'
        cls.enzyme = cls.root / 'enzyme.dat'
        cls.statements = cls.root / 'statements.jsonl'
        cls.obo.write_text(OBO)
        cls.enzyme.write_text(ENZYME)
        sentences = [
            'Does not possess fructose-bisphosphate (FBP) aldolase activity',
            'Does not exhibit glutaminase activity',
            'Does not display 2-isopropylmalate synthase and citramalate synthase activities since it cannot use 2-oxoisovalerate or pyruvate as substrate',
            'lacks detectable in vitro ATPase activity',
        ]
        rows = [
            statement('P0001', sentences[0], 'p1',
                      go_mf=[{'go_id': 'GO:0000001', 'evidence': 'IDA'}]),
            statement('P0002', sentences[1], 'p2'),
            statement('P0003', sentences[2], 'p3'),
            statement('P0004', sentences[3], 'p4', entry_ec=['EC:3.6.1.3']),
            {
                'kind': 'gaf_not', 'taxid': '1', 'accession': 'P0005',
                'sentence': 'GO:0000002 glutaminase activity',
                'statement_class': 'gaf_not', 'go_id': 'GO:0000002',
                'pubmed_ids': ['456'],
                'joins': [{'genome_id': 'g1', 'protein_id': 'p5'}],
                'entry_ec': [], 'entry_go_mf': [],
            },
        ]
        jsonl(cls.statements, rows)
        cls.proposals = cls.root / 'proposals'
        cls.inventory = P.propose(cls.statements, cls.obo, cls.enzyme, cls.proposals)
        cls.records = [json.loads(line) for line in
                       (cls.proposals / 'proposals.jsonl').read_text().splitlines()]

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def test_matches_spans_normalisation_and_entry_conflicts(self):
        self.assertEqual([p['function_id'] for p in self.records[0]['proposed_functions']],
                         ['GO:0000001'])
        fbp = self.records[0]['proposed_functions'][0]
        self.assertEqual(fbp['source'], 'go_name')
        self.assertEqual(fbp['matched_text'], 'fructose-bisphosphate (FBP) aldolase')
        sentence = self.records[0]['sentence']
        self.assertEqual(sentence[slice(*fbp['span'])], fbp['matched_text'])
        self.assertTrue(fbp['in_negation_scope'])
        self.assertEqual(fbp['also_matched'][0]['source'], 'go_synonym:RELATED')
        self.assertEqual(self.records[0]['proposed_conflicts_with_entry'], ['GO:0000001'])
        self.assertTrue(fbp['also_matched'])

        third = {p['function_id']: p for p in self.records[2]['proposed_functions']}
        self.assertEqual(set(third), {'GO:0000003', 'GO:0000004', 'EC:2.3.3.99'})
        self.assertTrue(third['GO:0000003']['in_negation_scope'])
        self.assertTrue(third['GO:0000004']['in_negation_scope'])
        for proposal in third.values():
            self.assertEqual(self.records[2]['sentence'][slice(*proposal['span'])],
                             proposal['matched_text'])

        atp = self.records[3]['proposed_functions'][0]
        self.assertEqual((atp['function_id'], atp['source']), ('EC:3.6.1.3', 'enzyme_an'))
        self.assertEqual(self.records[3]['proposed_conflicts_with_entry'], ['EC:3.6.1.3'])
        all_ids = {p['function_id'] for r in self.records for p in r['proposed_functions']}
        self.assertNotIn('GO:0000005', all_ids)  # obsolete
        self.assertNotIn('GO:0000006', all_ids)  # biological process
        self.assertNotIn('EC:1.2.3.4', all_ids)  # transferred

        matcher = P.build_matcher(self.obo, self.enzyme, 6)
        spaced = statement('X', 'Does not exhibit fructose bisphosphate aldolase activity', 'x')
        self.assertEqual(P.propose_record(spaced, matcher)['proposed_functions'][0]['function_id'],
                         'GO:0000001')
        self.assertEqual(P.normalize_name('Fructose-bisphosphate, aldolase activity'),
                         'fructose bisphosphate aldolase')

    def test_negation_scope_stop_list_and_gaf_passthrough(self):
        for function_id in ('GO:0000003', 'GO:0000004'):
            proposal = next(p for p in self.records[2]['proposed_functions']
                            if p['function_id'] == function_id)
            self.assertTrue(proposal['in_negation_scope'])
        self.assertIn('Does not display 2-isopropylmalate synthase and citramalate synthase',
                      self.records[2]['negation_scope'][0])

        matcher = P.build_matcher(self.obo, self.enzyme, 6)
        positive_first = statement(
            'X', 'Exhibits glutaminase activity, but does not exhibit ATPase activity.', 'x')
        proposed = P.propose_record(positive_first, matcher)['proposed_functions']
        by_id = {p['function_id']: p for p in proposed}
        self.assertFalse(by_id['GO:0000002']['in_negation_scope'])
        self.assertTrue(by_id['EC:3.6.1.3']['in_negation_scope'])
        stopped = statement('X', 'Does not display synthase activity', 'x')
        self.assertEqual(P.propose_record(stopped, matcher)['proposed_functions'], [])
        punctuated = statement('X', 'Does not exhibit (glutaminase).', 'x')
        self.assertEqual(P.propose_record(punctuated, matcher)['proposed_functions'][0]['span'],
                         [18, 29])

        gaf = self.records[4]
        self.assertEqual(len(gaf['proposed_functions']), 1)
        proposal = gaf['proposed_functions'][0]
        self.assertEqual(proposal['function_id'], 'GO:0000002')
        self.assertEqual(proposal['source'], 'gaf_not')
        self.assertTrue(proposal['in_negation_scope'])
        self.assertEqual(gaf['review_status'], 'pending')
        self.assertFalse(gaf['accepted_training_label'])
        self.assertEqual(gaf['pubmed_ids'], ['456'])
        self.assertEqual(gaf['joins'][0]['protein_id'], 'p5')

    def test_abbreviations_specificity_and_broad_mapping_flag(self):
        matcher = P.build_matcher(self.obo, self.enzyme, 6)
        tal = P.propose_record(
            statement('P19669', 'Does not show fructose-6-P aldolase activity.', 'tal'),
            matcher,
        )['proposed_functions']
        self.assertEqual([item['function_id'] for item in tal], ['GO:0097023'])
        self.assertEqual(tal[0]['matched_text'], 'fructose-6-P aldolase')
        self.assertEqual(tal[0]['span'], [14, 35])
        self.assertIn('fructose 6-phosphate aldolase', tal[0]['normalised_text'])
        self.assertNotIn('mapping_too_broad', tal[0])
        self.assertIn('GO:0016832', {item['function_id'] for item in tal[0]['also_matched']})

        fgd = P.propose_record(statement(
            'P96809',
            'Does not exhibit F420-dependent glucose-6-phosphate dehydrogenase (FGD) activity',
            'fgd2',
        ), matcher)['proposed_functions']
        self.assertEqual([item['function_id'] for item in fgd], ['GO:0052749'])
        self.assertEqual(
            [item['function_id'] for item in fgd[0]['also_matched']], ['GO:0004345'])

        bare = P.propose_record(
            statement('X', 'Does not show aldolase activity.', 'x'), matcher
        )['proposed_functions']
        self.assertEqual([item['function_id'] for item in bare], ['GO:0016832'])
        self.assertTrue(bare[0]['mapping_too_broad'])
        self.assertTrue(bare[0]['in_negation_scope'])

    def test_summary_and_inventory_hashes(self):
        inventory = json.loads((self.proposals / 'inventory.json').read_text())
        self.assertEqual(inventory['schema'], 'gspa.negated-function-proposals/1')
        for key, path in [('statements', self.statements), ('go_obo', self.obo),
                          ('enzyme', self.enzyme)]:
            self.assertEqual(inventory['inputs'][key]['sha256'], P.sha256(path))
        self.assertEqual(inventory['min_name_length'], 6)
        self.assertEqual(inventory['stop_list'], sorted(P.STOP_LIST))
        self.assertEqual(inventory['counts']['statements_selected'], 5)
        self.assertEqual(inventory['counts']['proposed_mapping_too_broad'], 0)
        self.assertEqual(inventory['proposals_sha256'],
                         P.sha256(self.proposals / 'proposals.jsonl'))
        lines = (self.proposals / 'proposals_summary.tsv').read_text().splitlines()
        self.assertEqual(lines[0], 'accession\ttaxid\tclass\tn_proposed\tn_in_scope\tproposed_ids')
        self.assertEqual(len(lines), 6)
        self.assertEqual(lines[3].split('\t')[3:5], ['3', '3'])

    def test_feature_coverage_and_reviewable_loci(self):
        feature_dir = self.root / 'features'
        feature_dir.mkdir()
        rows = [
            {'protein_id': 'p1', 'function_type': 'GO', 'function_id': 'GO:0000001',
             'role': 'x', 'candidate_key': '1'},
            {'protein_id': 'p2', 'function_type': 'GO', 'function_id': 'GO:0000001',
             'role': 'x', 'candidate_key': '2'},
            {'protein_id': 'p2', 'function_type': 'GO', 'function_id': 'GO:0000002',
             'role': 'x', 'candidate_key': '3'},
            {'protein_id': 'p6', 'function_type': 'GO', 'function_id': 'GO:0000002',
             'role': 'x', 'candidate_key': '4'},
            {'protein_id': 'p3', 'function_type': 'GO', 'function_id': 'GO:0000003',
             'role': 'x', 'candidate_key': '5'},
            {'protein_id': 'p3', 'function_type': 'GO', 'function_id': 'GO:0000004',
             'role': 'x', 'candidate_key': '6'},
            {'protein_id': 'p4', 'function_type': 'EC', 'function_id': '3.6.1.3',
             'role': 'x', 'candidate_key': '7'},
            {'protein_id': 'p7', 'function_type': 'EC', 'function_id': 'EC:3.6.1.3',
             'role': 'x', 'candidate_key': '8'},
        ]
        (feature_dir / 'g1').mkdir()
        jsonl(feature_dir / 'g1' / 'candidate_features.jsonl', rows)
        report = P.coverage(self.proposals / 'proposals.jsonl', self.root / 'feature-coverage',
                            features=feature_dir)
        self.assertEqual(report['proposals_total'], 7)
        self.assertEqual(report['proposals_in_scope'], 7)
        coverage_rows = [json.loads(line) for line in
                         (self.root / 'feature-coverage/coverage.jsonl').read_text().splitlines()]
        atp = next(row for row in coverage_rows if row['function_id'] == 'EC:3.6.1.3')
        self.assertTrue(atp['protein_in_universe'])
        self.assertEqual(atp['universe_members'], ['p4', 'p7'])
        self.assertEqual(atp['universe_size'], 2)
        self.assertFalse(atp['truncated'])
        gaf = next(row for row in coverage_rows if row['statement_index'] == 4)
        self.assertFalse(gaf['protein_in_universe'])
        self.assertTrue(any(c['function_id'] == 'EC:3.6.1.3'
                            for c in report['reviewable_locus_candidates']))
        self.assertEqual(report['accepted_training_labels'], 0)

    def test_claims_coverage_uses_accessions_and_max_source_scores(self):
        claims = self.root / 'claims.jsonl'
        rows = [
            {'protein_id': 'P0002', 'function_type': 'GO', 'function_id': 'GO:0000002',
             'source': 'diamond', 'raw_score': 0.2},
            {'protein_id': 'P0002', 'function_type': 'GO', 'function_id': 'GO:0000002',
             'source': 'diamond', 'raw_score': 0.7},
            {'protein_id': 'P0002', 'function_type': 'GO', 'function_id': 'GO:0000002',
             'source': 'deepgo', 'raw_score': 0.4},
            {'protein_id': 'P0006', 'function_type': 'GO', 'function_id': 'GO:0000002',
             'source': 'diamond', 'raw_score': 0.9},
        ]
        jsonl(claims, rows)
        report = P.coverage(self.proposals / 'proposals.jsonl', self.root / 'claims-coverage',
                            claims=claims, genome_id='genome-claims')
        coverage_rows = [json.loads(line) for line in
                         (self.root / 'claims-coverage/coverage.jsonl').read_text().splitlines()]
        row = next(row for row in coverage_rows if row['statement_index'] == 1)
        self.assertEqual(row['genome_id'], 'genome-claims')
        self.assertEqual(row['protein_id'], 'P0002')
        self.assertTrue(row['protein_in_universe'])
        self.assertEqual(row['universe_members'], ['P0002', 'P0006'])
        self.assertEqual(row['member_scores']['P0002'], {'deepgo': 0.4, 'diamond': 0.7})
        self.assertEqual(row['member_scores']['P0006'], {'diamond': 0.9})
        self.assertEqual(report['mode'], 'claims')


if __name__ == '__main__':
    unittest.main()
