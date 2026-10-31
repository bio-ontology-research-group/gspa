import copy
import hashlib
import json
import tempfile
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import summarize_source_context_replay as m

GENOMES = ('G1', 'G2')


def ec_entry(i):
    recs = [{'protein_id': f'p{j}'} for j in range(3)]
    return {'emitted_reference_pairs': 3, 'missed_reference': 2, 'policy': {'schema': 'x'}, 'proteins': 3,
            'proteins_with_candidates': 3, 'records': recs, 'reference_matches': 2, 'reference_micro_f1': 0.5 + i / 100,
            'reference_pairs': 4, 'reference_pairs_by_genome': {'G1': 3, 'G2': 1}, 'reference_precision': None if i == 7 else 0.6,
            'reference_recall': 0.5, 'test_genomes': ['G1', 'G2'], 'threshold': 0.25, 'top1_matches': 2, 'top1_tied': 0,
            'top3_matches': 2, 'top5_matches': 3, 'unconfirmed_calls': 1}


def cond_entry(i):
    return {'brier_denominator': 6, 'brier_on_available_scores': 0.2, 'label_basis': {'closed_world_unconfirmed': 2,
            'positive_truth': 4}, 'missed_reference': 1, 'missing_predictions': 0, 'positive': 4, 'reference_f1': 0.7,
            'reference_matched_calls': 3, 'rows': 6, 'threshold': 0.25, 'uncalled_unconfirmed_or_negative': 1,
            'unconfirmed_or_negative_calls': 1}


def cmp_entry(gain=1):
    pg = {'G1': {'proteins': 2, 'reference_calls_gained': gain, 'reference_calls_lost': 0, 'top1_gained': 0,
                 'top1_lost': 0, 'unconfirmed_call_change': -1},
          'G2': {'proteins': 1, 'reference_calls_gained': 0, 'reference_calls_lost': 1, 'top1_gained': 0,
                 'top1_lost': 1, 'unconfirmed_call_change': 0}}
    return {'per_genome': pg, 'proteins': 3, 'reference_calls_gained': gain, 'reference_calls_lost': 1,
            'top1_gained': 0, 'top1_lost': 1, 'unconfirmed_call_change': -1}


def comparison_keys():
    keys = []
    for s in m.SOURCES:
        for c in m.CONTEXTS:
            for k in ('logistic', 'gbdt'):
                keys.append(f'{s}/{c}/{k}_vs_{s}/{c}/rule')
        for c in m.CONTEXTS[1:]:
            for k in m.KINDS:
                keys.append(f'{s}/{c}/{k}_vs_{s}/SEQONLY/{k}')
    keys.append('clean_diamond_homology/UNIVERSAL/rule_vs_clean/UNIVERSAL/rule')
    for a in m.arm_names():
        keys.append(f'{a}_vs_clean_native')
    return keys


def make_results():
    arms = m.arm_names() + [m.NATIVE]
    res = {}
    for op in m.ARM_OPERATING_POINTS:
        res[op] = {
            'full_ec_reference': {a: ec_entry(i) for i, a in enumerate(arms)},
            'conditional_candidate_metrics': {a: cond_entry(i) for i, a in enumerate(arms)},
            'comparisons_to_native_emitted_calls': {a: cmp_entry() for a in arms},
            'comparisons': {k: cmp_entry(2) for k in comparison_keys()}}
    res[m.NATIVE_OPERATING_POINT] = {'full_ec_reference': ec_entry(1), 'conditional_candidate_metrics': cond_entry(1)}
    return {'policy': {'schema': 'gspa.source-context-analysis/1'}, 'results': res}


def replay_cell(ties=None, scored_only=0):
    actions = {'rescored': 8, 'unchanged_unsupported': 2}
    reasons = {'scored': 8 + scored_only, 'profile_lacks_both_classes_in_train_and_calibration': 2 - scored_only}
    if scored_only:
        actions = {'rescored': 8, 'scored_only': scored_only, 'unchanged_unsupported': 2 - scored_only}
    cell = {'actions': actions, 'added': 0, 'removed': 0, 'candidates': 10, 'functions': 10,
            'max_probability_error': 1e-16, 'reasons': reasons}
    if ties is not None:
        cell['display_ties'] = len(ties)
        cell['display_tie_candidates'] = ties
    return cell


def make_verification():
    summ = {s: {c: {k: {g: replay_cell() for g in GENOMES} for k in m.REPLAY_KINDS} for c in m.CONTEXTS}
            for s in m.SOURCES}
    summ['clean']['SEQONLY']['gbdt']['G1'] = replay_cell(ties=['a', 'b'], scored_only=1)
    summ['clean']['SEQONLY']['gbdt']['G2']['max_probability_error'] = 3e-15
    return {'schema': 'gspa.source-context-production/1', 'status': 'verified', 'summaries': summ}


class ReplaySummaryTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(dir=Path.home() / 'tmp' if (Path.home() / 'tmp').is_dir() else None)
        self.root = Path(self.tmp.name)
        self.analysis = make_results()
        self.verification = make_verification()
        self.write()

    def tearDown(self):
        self.tmp.cleanup()

    def write(self):
        (self.root / 'analysis').mkdir(exist_ok=True)
        (self.root / 'replay').mkdir(exist_ok=True)
        (self.root / 'analysis' / 'results.json').write_text(json.dumps(self.analysis))
        (self.root / 'replay' / 'verification.json').write_text(json.dumps(self.verification))

    def run_into(self, name):
        out = self.root / name
        m.run(self.root, out)
        return out

    def test_outputs_are_byte_identical_across_runs(self):
        a, b = self.run_into('o1'), self.run_into('o2')
        for f in (m.RESULTS_TSV, m.INTEGRITY_TSV, m.MD_NAME):
            self.assertEqual(hashlib.sha256((a / f).read_bytes()).hexdigest(),
                             hashlib.sha256((b / f).read_bytes()).hexdigest())

    def test_tsv_rows_and_fields(self):
        out = self.run_into('o')
        lines = (out / m.RESULTS_TSV).read_text().splitlines()
        header = lines[0].split('\t')
        rows = [dict(zip(header, l.split('\t'))) for l in lines[1:]]
        per_op = 61 * 3 + len(comparison_keys())  # ec + conditional + to-native + pairwise comparisons
        self.assertEqual(len(rows), 3 * per_op + 2)
        self.assertEqual(len({(r['operating_point'], r['family'], r['name']) for r in rows}), len(rows))
        ec = next(r for r in rows if r['operating_point'] == 'fixed_0.5' and r['family'] == 'full_ec_reference'
                  and r['name'] == 'clean/ORGANISM/gbdt')
        self.assertEqual(ec['reference_pairs_by_genome.G1'], '3')
        self.assertEqual(ec['proteins'], '3')
        self.assertEqual(ec['test_genomes'], 'G1;G2')
        none_row = [r for r in rows if r['family'] == 'full_ec_reference' and r['reference_precision'] == '']
        self.assertTrue(none_row)  # undefined precision stays an empty cell
        cmp_row = next(r for r in rows if r['name'] == 'clean/ORGANISM/gbdt_vs_clean/ORGANISM/rule')
        self.assertEqual(cmp_row['comparison_type'], 'learned_vs_rule')
        nat_row = next(r for r in rows if r['family'] == 'comparison_to_native_emitted')
        self.assertEqual(nat_row['baseline'], 'native_emitted_calls')
        self.assertEqual(cmp_row['per_genome.G1.reference_calls_gained'], '2')

    def test_comparison_types(self):
        self.assertEqual(m.classify('clean/ORGANISM/rule_vs_clean_native')[2], 'vs_clean_native')
        self.assertEqual(m.classify('clean/PATHWAY/gbdt_vs_clean/SEQONLY/gbdt')[2], 'context_vs_seqonly')
        self.assertEqual(m.classify('clean/PATHWAY/gbdt_vs_clean/PATHWAY/rule')[2], 'learned_vs_rule')
        self.assertEqual(m.classify('clean_diamond_homology/PATHWAY/rule_vs_clean/PATHWAY/rule')[2], 'source_vs_source')
        with self.assertRaises(m.SummaryError):
            m.classify('clean/PATHWAY/rule_vs_clean/UNIVERSAL/gbdt')

    def test_integrity_totals_and_markdown(self):
        out = self.run_into('o')
        rows = m.integrity_rows(self.verification)
        totals = m.source_totals(rows)
        self.assertEqual(totals['clean']['genomes'], 2)
        self.assertEqual(totals['clean']['candidates'], 8 * 2 * 10)
        self.assertEqual(totals['clean']['scored_only'], 1)
        self.assertEqual(totals['clean']['display_ties'], 2)
        self.assertEqual(totals['clean']['unchanged_unsupported'], 8 * 2 * 2 - 1)
        self.assertEqual(totals['clean']['max_probability_error'], 3e-15)
        md = (out / m.MD_NAME).read_text()
        self.assertIn('## Replay integrity per source', md)
        for op in m.OP_ORDER:
            self.assertIn(f'## Operating point `{op}`', md)
        self.assertIn('n/a', md)
        self.assertIn('### Learned versus rule', md)
        self.assertIn('### Arm versus native CLEAN emitted calls', md)

    def test_unverified_replay_fails(self):
        self.verification['status'] = 'failed'
        self.write()
        with self.assertRaises(m.SummaryError):
            self.run_into('o')

    def test_missing_operating_point_fails(self):
        del self.analysis['results']['fixed_0.5']
        self.write()
        with self.assertRaises(m.SummaryError):
            self.run_into('o')

    def test_unexpected_field_fails(self):
        self.analysis['results']['fixed_0.5']['full_ec_reference']['clean/ORGANISM/gbdt']['surprise'] = 1
        self.write()
        with self.assertRaises(m.SummaryError):
            self.run_into('o')

    def test_missing_arm_fails(self):
        del self.analysis['results']['calibration_selected']['conditional_candidate_metrics']['clean/SEQONLY/rule']
        self.write()
        with self.assertRaises(m.SummaryError):
            self.run_into('o')

    def test_inconsistent_comparison_sum_fails(self):
        k = comparison_keys()[0]
        self.analysis['results']['ec_calibration_selected']['comparisons'][k]['reference_calls_gained'] += 1
        self.write()
        with self.assertRaises(m.SummaryError):
            self.run_into('o')

    def test_replay_actions_must_sum(self):
        self.verification['summaries']['clean']['UNIVERSAL']['gbdt']['G1']['candidates'] = 11
        self.write()
        with self.assertRaises(m.SummaryError):
            self.run_into('o')

    def test_tie_count_must_match_list(self):
        c = self.verification['summaries']['clean']['SEQONLY']['gbdt']['G1']
        c['display_ties'] = 5
        self.write()
        with self.assertRaises(m.SummaryError):
            self.run_into('o')

    def test_genome_set_must_match_within_source(self):
        del self.verification['summaries']['clean']['PATHWAY']['logistic']['G2']
        self.write()
        with self.assertRaises(m.SummaryError):
            self.run_into('o')


if __name__ == '__main__':
    unittest.main()
