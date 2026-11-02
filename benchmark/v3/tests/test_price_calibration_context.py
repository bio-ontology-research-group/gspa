import hashlib
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from price_calibration_context import row1_calibration as R1  # noqa: E402
from price_calibration_context import build_context_inputs as B  # noqa: E402


class SplitRule(unittest.TestCase):
    def test_greedy_hash_order_and_merge(self):
        metadata = {'genomes': {
            'A': {'clade': 'ncbi_genus_286'}, 'B': {'clade': 'ncbi_genus_2901164'},
            'C': {'clade': 'ncbi_genus_22'}, 'D': {'clade': 'ncbi_genus_1'}}}
        reference = {('A', f'a{i}'): frozenset({'1.1.1.1'}) for i in range(50)}
        reference.update({('B', f'b{i}'): frozenset({'1.1.1.1'}) for i in range(10)})
        reference.update({('C', f'c{i}'): frozenset({'1.1.1.1'}) for i in range(45)})
        reference.update({('D', f'd{i}'): frozenset({'1.1.1.1'}) for i in range(5)})
        split = R1.genus_split(metadata, reference)
        # Stutzerimonas merged into the Pseudomonas block
        self.assertEqual(sorted(split['blocks']['ncbi_genus_286']), ['A', 'B'])
        self.assertNotIn('ncbi_genus_2901164', split['blocks'])
        order = sorted(split['blocks'], key=lambda b: hashlib.sha256((R1.SPLIT_SALT + b).encode()).hexdigest())
        self.assertEqual(split['block_order'], order)
        # greedy: the first block goes to calibration, every block lands in the lighter arm
        totals = {'calibration': 0, 'test': 0}
        for b in order:
            arm = 'calibration' if totals['calibration'] <= totals['test'] else 'test'
            self.assertIn(b, split['arms'][arm]['blocks'])
            totals[arm] += sum(1 for g, _ in reference if g in split['blocks'][b])
        self.assertEqual({a: s['reference_proteins'] for a, s in split['arms'].items()}, totals)
        self.assertEqual(set(split['assignment']), {'A', 'B', 'C', 'D'})


class ThresholdRule(unittest.TestCase):
    def rows(self, supported_scores, unsupported_scores):
        rows = [('p', f's{i}', s, True) for i, s in enumerate(supported_scores)]
        rows += [('p', f'u{i}', s, False) for i, s in enumerate(unsupported_scores)]
        return rows

    def test_abstains_below_ten_calls(self):
        pol = R1.choose_threshold(self.rows([0.9] * 9, []))
        self.assertIsNone(pol['threshold'])
        self.assertEqual(pol['status'], 'no_supported_threshold_abstain')

    def test_abstains_below_ninety_percent(self):
        pol = R1.choose_threshold(self.rows([0.9] * 10, [0.95, 0.95]))
        self.assertIsNone(pol['threshold'])

    def test_selects_lowest_admissible_threshold_with_most_supported(self):
        # ten supported at 0.8, one unsupported at 0.7, nine supported at 0.6:
        # t=0.8 -> 10/10; t=0.7 -> 10/11 (0.909); t=0.6 -> 19/20 (0.95) -> most supported wins
        pol = R1.choose_threshold(self.rows([0.8] * 10 + [0.6] * 9, [0.7]))
        self.assertEqual(pol['status'], 'selected_on_calibration')
        self.assertAlmostEqual(pol['threshold'], 0.6)
        self.assertEqual(pol['calibration_calls']['supported'], 19)
        self.assertEqual(pol['calibration_calls']['unsupported'], 1)

    def test_ties_prefer_fewer_unsupported_then_higher_threshold(self):
        # t=0.9 -> 10/10 ; t=0.5 -> 10 supported + 1 unsupported (0.909): same supported, more unsupported -> keep 0.9
        pol = R1.choose_threshold(self.rows([0.9] * 10, [0.5]))
        self.assertAlmostEqual(pol['threshold'], 0.9)

    def test_call_counts(self):
        rows = self.rows([0.9, 0.3], [0.95])
        truth = {('p', 's0'), ('p', 's1'), ('q', 'x')}
        c = R1.call_counts(rows, truth, 0.5)
        self.assertEqual((c['calls'], c['supported'], c['unconfirmed'], c['missed']), (2, 1, 1, 2))


class Bootstrap(unittest.TestCase):
    def test_point_estimate_and_interval_contain_truth(self):
        vectors = {f'g{i}': [i % 2, 1] for i in range(6)}
        b = R1.genome_bootstrap(vectors, R1.ratio(0, 1), n=200, seed=1)
        self.assertAlmostEqual(b['estimate'], 0.5)
        self.assertLessEqual(b['ci_low'], 0.5)
        self.assertGreaterEqual(b['ci_high'], 0.5)


class ObservationConverter(unittest.TestCase):
    REACTION_HEADER = ['pathway', 'rxn', 'ec', 'keyrea', 'spont', 'is_complex', 'qseqid', 'pident', 'evalue',
                       'bitscore', 'qcovs', 'stitle', 'sstart', 'send', 'status']

    def write(self, d, name, header, rows):
        p = Path(d) / name
        p.write_text('\t'.join(header) + '\n' + ''.join('\t'.join(r) + '\n' for r in rows))
        return p

    def test_tables_to_observation(self):
        with tempfile.TemporaryDirectory() as d:
            rx = self.write(d, 'G-all-Reactions.tbl', self.REACTION_HEADER, [
                ['|PWY-1|', 'rxn00001', '1.1.1.1', 'TRUE', 'FALSE', 'FALSE', 'ref1', '80.0', '1e-50', '200', '90', 'prot1 desc', '1', '300', 'good_blast'],
                ['|PWY-1|', 'rxn00002', '2.2.2.2', 'FALSE', 'FALSE', 'FALSE', '', 'NA', 'NA', 'NA', 'NA', '', 'NA', 'NA', 'no_blast']])
            pw = self.write(d, 'G-all-Pathways.tbl', ['ID', 'Name', 'Prediction'], [['|PWY-1|', 'x', 'true']])
            obs = B.observation('G', rx, pw)
            self.assertEqual(obs['schema'], 'gspa.pathway-search-observations/1')
            self.assertEqual((obs['source'], obs['target'], obs['outcome']), ('gapsmith', 'proteome', 'SUCCESS'))
            self.assertEqual(len(obs['reaction_rows']), 2)
            self.assertEqual(obs['pathway_rows'][0]['Prediction'], 'true')
            self.assertEqual(obs['reaction_table_sha256'], hashlib.sha256(rx.read_bytes()).hexdigest())
            self.assertEqual(obs['pathway_table_sha256'], hashlib.sha256(pw.read_bytes()).hexdigest())

    def test_conflicting_pipe_variants_are_dropped_and_recorded(self):
        with tempfile.TemporaryDirectory() as d:
            rx = self.write(d, 'G-all-Reactions.tbl', self.REACTION_HEADER, [])
            pw = self.write(d, 'G-all-Pathways.tbl', ['ID', 'Prediction'], [['P1', 'TRUE'], ['|P1|', 'FALSE'], ['|P2|', 'TRUE']])
            obs = B.observation('G', rx, pw)
            self.assertEqual(obs['conflicting_pathway_ids_dropped'], ['P1'])
            self.assertEqual([r['ID'] for r in obs['pathway_rows']], ['|P2|'])

    def test_rejects_unknown_status_and_missing_column(self):
        with tempfile.TemporaryDirectory() as d:
            rx = self.write(d, 'G-all-Reactions.tbl', self.REACTION_HEADER, [
                ['P', 'r', '1.1.1.1', 'TRUE', 'FALSE', 'FALSE', '', 'NA', 'NA', 'NA', 'NA', '', 'NA', 'NA', 'odd']])
            pw = self.write(d, 'G-all-Pathways.tbl', ['ID', 'Prediction'], [['P', 'true']])
            with self.assertRaises(ValueError):
                B.observation('G', rx, pw)
            rx2 = self.write(d, 'H-all-Reactions.tbl', self.REACTION_HEADER[:-1], [])
            with self.assertRaises(ValueError):
                B.observation('H', rx2, pw)

    def test_plausible_reactions_from_sbml(self):
        sbml = ('<?xml version="1.0"?><sbml xmlns="http://www.sbml.org/sbml/level3/version1/core" level="3" version="1">'
                '<model id="m"><listOfReactions><reaction id="R_rxn00001_c0"/><reaction id="R_rxn00002_c0"/>'
                '<reaction id="R_EX_cpd00001_e0"/><reaction id="R_bio1"/></listOfReactions></model></sbml>')
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / 'G-draft.xml'
            p.write_text(sbml)
            kept, dropped = B.plausible_reactions(p, {'rxn00001', 'rxn00002', 'rxn99999'})
            self.assertEqual(kept, ['rxn00001', 'rxn00002'])
            self.assertEqual(sorted(dropped), ['EX_cpd00001', 'bio1'])

    def test_gap_linked_loci(self):
        gaps = [dict(pathway_id='gapsmith:P', reaction_id='rxn1', ec_number='1.1.1.1', go_term=None, gapseq_guessed=False)]
        linked = B.gap_linked({'p1': ['1.1.1.1'], 'p2': ['2.2.2.2']}, gaps)
        self.assertEqual(linked, {'p1': ['rxn1']})


if __name__ == '__main__':
    unittest.main()


class CellHelpers(unittest.TestCase):
    def test_retain_and_differences(self):
        from price_calibration_context import run_context_cell as C
        header = 'protein_id\ttype\tfunction_id\tposterior_prob\textra\n'
        with tempfile.TemporaryDirectory() as d:
            a = Path(d) / 'a.tsv'
            b = Path(d) / 'b.tsv'
            a.write_text(header + 'p1\tEC\tEC:1.1.1.1\t0.5\tx\np1\tGO\tGO:0001\t0.2\tx\np2\tEC\t2.2.2.2\t0.9\tx\n')
            b.write_text(header + 'p1\tEC\tEC:1.1.1.1\t0.6\tx\np1\tGO\tGO:0001\t0.2\tx\np2\tEC\t2.2.2.2\t0.9\tx\n')
            rows, retained, digest = C.retain_scores(a, {'p1': ['1.1.1.1']})
            self.assertEqual(rows, 3)
            self.assertEqual(retained, {'p1': {'1.1.1.1': 0.5}})
            self.assertEqual(digest, hashlib.sha256(a.read_bytes()).hexdigest())
            self.assertEqual(C.count_differences(a, b), dict(rows=3, differing=1))
            c = Path(d) / 'c.tsv'
            c.write_text(header + 'p9\tEC\tEC:1.1.1.1\t0.6\tx\n')
            with self.assertRaises(ValueError):
                C.count_differences(a, c)
