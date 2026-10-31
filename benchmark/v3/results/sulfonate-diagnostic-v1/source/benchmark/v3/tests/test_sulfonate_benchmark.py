"""Protect real study coverage, ambiguous genotypes, and observation separation."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest

V3 = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(V3))
import prepare_sulfonate_benchmark as S

SOURCES = V3 / 'results/sulfonate-diagnostic-v1/sources'
TABLE = SOURCES / 'publisher-web-extraction.json'
METHODS = SOURCES / 'publisher-methods-web-extraction.json'


class SulfonateExtractionTest(unittest.TestCase):
    def test_complete_population_and_source_coordinates(self):
        inputs, observations, report = S.extract(TABLE, METHODS)
        self.assertEqual(len(inputs), 260)
        self.assertEqual(report['reported_cells'], 247)
        self.assertEqual(report['categories'], {'+++': 68, '−': 141, '+': 28, '++': 10, 'unreported': 13})
        self.assertEqual(len({r['case_id'] for r in inputs}), 260)
        self.assertEqual([r['case_id'] for r in inputs], [r['case_id'] for r in observations])
        cells = {(o['source']['row'], o['source']['column']): o for o in observations}
        self.assertEqual(cells['Taurine', 'tauA']['symbol'], '−')
        self.assertEqual(cells['Butanesulfonate', 'ssuCB tauBC']['symbol'], '+')
        self.assertIsInstance(cells['Taurine', 'tauA']['source']['rendered_line'], int)

    def test_missing_compound_is_not_a_negative(self):
        inputs, observations, _ = S.extract(TABLE, METHODS)
        missing_ids = {r['case_id'] for r in inputs if r['sulfur_source'] == 'Dodecanesulfonate'}
        missing = [o for o in observations if o['case_id'] in missing_ids]
        self.assertEqual(len(missing), 13)
        for o in missing:
            self.assertIsNone(o['symbol'])
            self.assertIsNone(o['ordinal'])
            self.assertIsNone(o['minimum_independent_experiments'])
            self.assertEqual(o['status'], 'not_reported_in_table')

    def test_conflicted_genotype_retains_label_but_no_executable_deletion(self):
        inputs, _, report = S.extract(TABLE, METHODS)
        self.assertEqual(report['uncertain_genotype_cases'], 20)
        for r in inputs:
            if r['genotype_label'] == 'ssuAB tauC':
                self.assertIsNone(r['deleted_genes'])
                self.assertEqual(r['reported_deleted_genes'], ['ssuA', 'ssuB', 'tauC'])
                self.assertEqual(r['genotype_status'], 'table_text_conflict')
        self.assertEqual(S.genes('ssuCB tauBC'), ['ssuB', 'ssuC', 'tauB', 'tauC'])

    def test_normalization_reference_and_measured_control_are_distinct(self):
        _, observations, _ = S.extract(TABLE, METHODS)
        for o in observations:
            if o['source']['column'] == 'None':
                self.assertEqual(o['observation_role'], 'wild_type_normalization_reference')
            elif o['source']['row'] == 'Sulfate':
                self.assertEqual(o['observation_role'], 'sulfate_control')
        self.assertEqual(S.CATEGORIES['+']['relative_yield_interval'], '5–33%')

    def test_outcome_changes_do_not_change_inputs(self):
        inputs, observations, _ = S.extract(TABLE, METHODS)
        text = json.loads(TABLE.read_text())
        line = next(line for line in text.splitlines() if ': Taurine ' in line)
        changed = text.replace(line, line.replace('−', '+++', 1))
        scratch = Path.home() / 'tmp'
        scratch.mkdir(exist_ok=True)
        with tempfile.TemporaryDirectory(dir=scratch) as work:
            path = Path(work) / 'changed.json'
            path.write_text(json.dumps(changed))
            altered_inputs, altered_observations, _ = S.extract(path, METHODS)
        self.assertEqual(inputs, altered_inputs)
        self.assertNotEqual(observations, altered_observations)

    def test_rejects_missing_duplicate_or_corrupt_cells(self):
        lines = S.source_lines(TABLE)
        index = next(i for i, (_, line) in enumerate(lines) if line.split('|')[0].strip() == 'Taurine')
        missing = copy.deepcopy(lines);missing.pop(index)
        duplicate = copy.deepcopy(lines);duplicate.insert(index, lines[index])
        malformed = copy.deepcopy(lines);malformed[index] = (lines[index][0], 'Taurine | +++')
        invalid = copy.deepcopy(lines);invalid[index] = (lines[index][0], lines[index][1].replace('−', '?', 1))
        for fixture in (missing, duplicate, malformed, invalid):
            with self.subTest(fixture=fixture[index]), self.assertRaises(ValueError):
                S.parse_table(fixture)


if __name__ == '__main__':
    unittest.main()
