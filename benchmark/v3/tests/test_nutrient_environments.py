"""Input independence, source integrity and raw-replicate reconstruction."""
import copy
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import prepare_nutrient_environments as P


def workbook(offset=0, strings=False):
    w = P.openpyxl.Workbook()
    s = w.active
    s.title = 'WT'
    s.cell(4, offset + 1, None)
    for col in range(1, 13):
        s.cell(4, offset + col + 1, col)
    for i, letter in enumerate('ABCDEFGH', 5):
        s.cell(i, offset + 1, letter)
        for col in range(1, 13):
            value = .1 + i / 100 + col / 1000
            if strings:
                nitrogen = list(P.NITROGEN)[col - 1] if col <= 10 else 'NoN'
                carbon = list(P.CARBON)[(i - 5) // 2]
                value = nitrogen + ' + ' + carbon
            s.cell(i, offset + col + 1, value)
    return w


def measurement(value, experiment=1, well='A1'):
    return {'species': 'Ec', 'environment': 'Glut + Gluc', 'pair': 'M+M', 'kind': 'monoculture', 'strains': ['M'],
            'nitrogen': 'Glut', 'carbon': 'Gluc', 'experiment': experiment, 'well': well,
            'corrected_od600': value, 'status': 'unavailable' if value is None else 'measured'}


class NutrientDataTests(unittest.TestCase):
    def test_model_preflight_uses_canonical_gene_identity_and_preserves_isozyme_bypass(self):
        sys.path.insert(0, str(Path(__file__).resolve().parents[3] / 'scripts/tests'))
        import register_nutrient_models as R
        from test_metabolic_experiments import write_sbml, libsbml
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            path, inputs = root / 'model.xml', root / 'inputs.jsonl'
            write_sbml(path)
            doc = libsbml.readSBMLFromFile(str(path))
            model = doc.getModel()
            model.getPlugin('fbc').getGeneProduct('g5').setName('metA')
            libsbml.writeSBMLToFile(doc, str(path))
            record = {'species': 'Ec', 'organisms': [{'deleted_genes': ['metA', 'missing_symbol']}]}
            inputs.write_text(json.dumps(record) + '\n')
            first = R.register(inputs, {'Ec': path})['models']['Ec']['genes']
            self.assertEqual(first[0]['blocked_reactions'], ['VITS'])
            self.assertEqual(first[0]['canonical'], 'g5')
            self.assertIsNone(first[1]['blocked_reactions'])
            association = model.getReaction('VITS').getPlugin('fbc').getGeneProductAssociation()
            association.unsetAssociation()
            either = association.createOr()
            either.createGeneProductRef().setGeneProduct('g5')
            either.createGeneProductRef().setGeneProduct('g6')
            libsbml.writeSBMLToFile(doc, str(path))
            second = R.register(inputs, {'Ec': path})['models']['Ec']['genes']
            self.assertEqual(second[0]['blocked_reactions'], [])
            self.assertEqual(second[0]['reactions_with_locus'], ['VITS'])

    def test_plate_locator_accepts_both_real_export_offsets_and_preserves_cells(self):
        for offset in (0, 1):
            grid = P.plate_grid(workbook(offset).active)
            self.assertEqual(len(grid), 96)
            self.assertEqual(grid['A1'][1:], (5, offset + 2))
            self.assertAlmostEqual(grid['H12'][0], .232)

    def test_plate_locator_rejects_duplicate_or_incomplete_layout(self):
        w = workbook()
        w.active['A12'] = 'Z'
        with self.assertRaises(P.DataError):
            P.plate_grid(w.active)
        w = workbook()
        for row in list(w.active.iter_rows(min_row=4, max_row=12)):
            for c in row:
                w.active.cell(c.row + 15, c.column, c.value)
        with self.assertRaises(P.DataError):
            P.plate_grid(w.active)

    def test_sample_map_validates_full_environment_replication_and_blank_locations(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / 'samples.xlsx'
            w = workbook(strings=True)
            w.save(p)
            mapping = P.sample_map(p)
            self.assertEqual(len(mapping), 96)
            self.assertEqual(mapping['G11']['environment'], 'NoN + NoC')
            # Preserve counts but move one blank into a carbon-supplied condition.
            w.active['L5'], w.active['L11'] = w.active['L11'].value, w.active['L5'].value
            w.save(p)
            with self.assertRaisesRegex(P.DataError, 'blank wells'):
                P.sample_map(p)

    def test_pair_identity_separates_supplementation_from_coculture(self):
        self.assertEqual(P.pair_kind('M+Met'), ('supplemented_monoculture', ['M']))
        self.assertEqual(P.pair_kind('M+F'), ('coculture', ['M', 'F']))
        self.assertEqual(P.pair_kind('M+M'), ('monoculture', ['M']))
        self.assertEqual(P.pair_kind('WT'), ('wild_type', []))
        with self.assertRaises(P.DataError):
            P.pair_kind('M+Ser')

    def test_invalid_measurements_are_not_silently_zero(self):
        self.assertIsNone(P.number(None))
        for value in ('0.1', True, float('nan'), float('inf')):
            with self.assertRaises(P.DataError):
                P.number(value)

    def test_pooled_wells_and_experiment_uncertainty_are_separate(self):
        rows = [measurement(.01, 1, 'A1'), measurement(.01, 1, 'B1'),
                measurement(.09, 2, 'A1'), measurement(.09, 2, 'B1')]
        r = next(r for r in P.summarize(rows) if r['case_id'] == 'NUTRIENT2026:Ec:Glut + Gluc:M+M')
        self.assertAlmostEqual(r['mean_od600'], .05)
        self.assertAlmostEqual(r['experiment_mean_sem'], .04)
        self.assertLess(r['pooled_well_sem'], r['experiment_mean_sem'])
        self.assertEqual(r['experiments_missing'], [3])
        self.assertEqual(r['n_wells_missing'], 2)

    def test_entire_missing_assay_stays_unknown_and_catalogue_complete(self):
        rows = P.summarize([])
        self.assertEqual(len(rows), 2464)
        self.assertEqual(sum(r['primary_environment'] for r in rows), 2240)
        self.assertTrue(all(r['mean_od600'] is None for r in rows))
        self.assertTrue(all(r['above_author_monoculture_threshold'] is None for r in rows))

    def test_outcome_and_coverage_changes_cannot_change_model_inputs(self):
        a = P.summarize([measurement(.005)])
        b = P.summarize([measurement(.9, 3)])
        c = P.summarize([])
        self.assertEqual(P.case_inputs(a), P.case_inputs(b))
        self.assertEqual(P.case_inputs(a), P.case_inputs(c))
        self.assertNotEqual(a, b)

    def test_species_specific_genotypes_and_real_supplement_concentrations_are_retained(self):
        rows = P.case_inputs(P.summarize([]))
        by_key = {(r['species'], r['environment'], r['pair']): r for r in rows}
        self.assertEqual(by_key['Ec', 'Glut + Gluc', 'R+R']['organisms'][0]['deleted_genes'], ['argA'])
        self.assertEqual(by_key['Bs', 'Glut + Gluc', 'R+R']['organisms'][0]['deleted_genes'], ['argH'])
        self.assertEqual([o['deleted_genes'] for o in by_key['Ec', 'Glut + Gluc', 'M+F']['organisms']], [['metA'], ['pheA']])
        self.assertEqual(by_key['Ec', 'Glut + Gluc', 'S+Ser']['medium']['supplement_concentration_uM'], 1000)
        self.assertEqual(by_key['Ec', 'Glut + Gluc', 'M+Met']['medium']['supplement_concentration_uM'], 200)
        self.assertIsNone(by_key['Ec', 'Glut + Gluc', 'M+Met']['medium']['uptake_rates'])
        self.assertIn('L-tryptophan', by_key['Ec', 'Glut + NoC', 'M+M']['medium']['base'])
        self.assertIsNone(by_key['Ec', 'Glut + NoC', 'M+M']['medium']['carbon_source'])

    def test_source_mean_in_sem_column_is_detected_without_replacing_source(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            folder = root / 'raw/Generated Tables'
            folder.mkdir(parents=True)
            text = 'Species,Environment,Pair,Mean,SEM\nEc,Glut + Gluc,M+M,0.05,0.05\n'
            for name in ('big_norm_df.csv', 'Supplemented_AA_data.csv'):
                (folder / name).write_text(text)
            observed = P.summarize([measurement(.03, 1), measurement(.07, 2)])
            report = P.audit_author_exports(root, observed)
            for result in report.values():
                self.assertEqual(result['mismatches_by_column'], {'SEM': 1})
            self.assertEqual((folder / 'Supplemented_AA_data.csv').read_text(), text)

    def test_source_content_and_pinned_blob_must_both_validate(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            paper = '<article><article-id>' + P.DOI + '</article-id><article-id>biorxiv;2026.05.20.726445v3</article-id></article>'
            files = {name: name.encode() for name in P.FILES}
            tree_sha = 'a' * 40
            tree = {'sha': tree_sha, 'truncated': False, 'tree': [
                {'type': 'blob', 'path': name, 'size': len(data),
                 'sha': hashlib.sha1(b'blob ' + str(len(data)).encode() + b'\0' + data).hexdigest()}
                for name, data in files.items()]}
            def fetch(url):
                if '/git/commits/' in url:
                    return json.dumps({'sha': P.COMMIT, 'tree': {'sha': tree_sha}}).encode()
                if '/git/trees/' in url:
                    self.assertIn(tree_sha, url)
                    return json.dumps(tree).encode()
                if url == P.PAPER_XML:
                    return paper.encode()
                return files[P.urllib.parse.unquote(url.split(P.COMMIT + '/')[1])]
            with patch.object(P, 'fetch', side_effect=fetch):
                P.acquire(root)
            P.verify_inputs(root)
            p = root / 'raw/README.md'
            p.write_text('changed')
            with self.assertRaisesRegex(P.DataError, 'source content changed'):
                P.verify_inputs(root)
            receipt_path = root / 'acquisition.json'
            receipt = json.loads(receipt_path.read_text())
            record = next(r for r in receipt['files'] if r['path'] == 'README.md')
            record.update(sha256=P.sha(p), bytes=p.stat().st_size)
            receipt_path.write_text(json.dumps(receipt))
            with self.assertRaisesRegex(P.DataError, 'pinned Git tree'):
                P.verify_inputs(root)

    def test_blank_or_measurement_missing_does_not_manufacture_growth(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            map_file = root / 'samples.xlsx'
            workbook(strings=True).save(map_file)
            mapping = P.sample_map(map_file)
            for species in P.GENES:
                for experiment in (1, 2, 3):
                    path = root / f'raw/Data/Crossfeeding Experiments/{species}_data_Ex{experiment}.xlsx'
                    path.parent.mkdir(parents=True, exist_ok=True)
                    w = workbook()
                    # Missing one of the four registered blank wells invalidates correction.
                    w.active['L11'] = None
                    w.save(path)
            rows, coverage = P.read_measurements(root, mapping)
            self.assertTrue(all(r['status'] == 'unavailable' for r in rows))
            self.assertTrue(all(r['corrected_od600'] is None for r in rows))
            self.assertTrue(all('M+M' in c['sheets_missing'] for c in coverage))


if __name__ == '__main__':
    unittest.main()
