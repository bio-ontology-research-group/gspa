"""WP5 step 1 (design §3, §8.1): checklist exclusions, label-flip invariance, gene-independent menus, Tier B gate."""
import json
import os
import shutil
import sys
import tempfile
import unittest
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import wp5_fixtures as F  # noqa: E402
import prepare_wp5_cases as P5  # noqa: E402


class Base(unittest.TestCase):
    def setUp(self):
        self.root = Path(tempfile.mkdtemp(prefix='wp5-cases-', dir=os.environ.get('TMPDIR')))
        self.addCleanup(shutil.rmtree, self.root, True)

    def prepare(self, name='bundle', flip=False, **extra):
        inputs = F.all_inputs(self.root / ('in-flip' if flip else 'in'), flip=flip)
        inputs.update(extra)
        receipt = P5.prepare(self.root / name, **inputs)
        return self.root / name, receipt


class ExclusionReasons(Base):
    def test_every_registered_reason_fires_and_nothing_is_dropped_silently(self):
        bundle, receipt = self.prepare()
        rows = P5.read_jsonl(bundle / 'private/ineligible.jsonl')
        fired = Counter(r['reason'] for r in rows)
        # tier_b_gated is a refusal of the whole Price input (tested below), never a row.
        self.assertEqual(set(fired), set(P5.REASONS) - {'tier_b_gated'})
        self.assertEqual(receipt['counts']['ineligible'], len(rows))
        for r in rows:
            self.assertIn(r['level'], ('source', 'case', 'pair', 'condition'))
            self.assertTrue(r['source_ids'])

    def test_specific_exclusions_land_on_the_expected_inputs(self):
        bundle, _ = self.prepare()
        by_reason = {}
        for r in P5.read_jsonl(bundle / 'private/ineligible.jsonl'):
            by_reason.setdefault(r['reason'], []).append(r)
        self.assertIn('R_EX_arab__L_e', by_reason['strain_background_affected_condition'][0]['detail'])
        self.assertIn('R_EX_xyl__D_e', by_reason['condition_not_expressible'][0]['detail'])
        self.assertEqual(by_reason['gene_ambiguous'][0]['detail'], 'dup')
        self.assertEqual(by_reason['split_overlap'][0]['source'], 'nutrient')
        self.assertTrue(any('Ec:MP' in r['detail'] for r in by_reason['multi_gene_genotype']))
        self.assertTrue(any(r['detail'] == 'tauB?' for r in by_reason['genotype_uncertain']))
        self.assertTrue(any(r['detail'] == 'LB' for r in by_reason['medium_protocol_unregistered']))
        self.assertTrue(any('fum' in r['detail'] for r in by_reason['medium_unrepresented_component']))

    def test_unregistered_reason_is_refused(self):
        builder = P5.Builder(P5.STRAINS)
        with self.assertRaises(P5.WP5Error):
            builder.reject('case', 'keio', ['x'], 'looked_at_the_outcome')

    def test_strain_without_cited_derivation_excludes_the_source(self):
        registry = self.root / 'registry.json'
        registry.write_text(json.dumps({'nutrient:Bs': {'handling': 'exclude_affected_conditions'}}))
        bundle, _ = self.prepare(strain_registry=registry)
        rows = P5.read_jsonl(bundle / 'private/ineligible.jsonl')
        self.assertTrue(any(r['reason'] == 'strain_unregistered' and r['source_ids'] == ['nutrient:Bs'] for r in rows))
        self.assertNotIn('Bs', {r['species'] for r in P5.read_jsonl(bundle / 'public/cases.jsonl')})


class LabelFlipInvariance(Base):
    def test_flipping_every_measured_label_changes_no_public_file_and_no_exclusion(self):
        a, _ = self.prepare('a')
        b, _ = self.prepare('b', flip=True)
        for name in ('public/cases.jsonl', 'public/manifest.json', 'private/ineligible.jsonl', 'private/targets.jsonl',
                     'private/truth.calibration.jsonl', 'private/truth.test.jsonl'):
            self.assertEqual((a / name).read_bytes(), (b / name).read_bytes(), name)
        # The labels themselves do differ, so the flip reached the transcription.
        self.assertNotEqual((a / 'private/labels.calibration.jsonl').read_bytes(),
                            (b / 'private/labels.calibration.jsonl').read_bytes())
        self.assertNotEqual((a / 'private/labels.test.jsonl').read_bytes(), (b / 'private/labels.test.jsonl').read_bytes())

    def test_public_rows_carry_no_label_truth_or_target_identity(self):
        bundle, _ = self.prepare()
        text = (bundle / 'public/cases.jsonl').read_text()
        for forbidden in ('growth', 'no_growth', 'label', 'mask_gene', 'G_b0001', 'pgiX', 'b0001', 'R_PGI', 'stratum'):
            self.assertNotIn(forbidden, text)

    def test_truth_and_labels_are_split_into_separate_private_files(self):
        bundle, _ = self.prepare()
        cal = {r['case'] for r in P5.read_jsonl(bundle / 'private/truth.calibration.jsonl')}
        test = {r['case'] for r in P5.read_jsonl(bundle / 'private/truth.test.jsonl')}
        self.assertTrue(cal and test and not cal & test)
        self.assertTrue({r['case'] for r in P5.read_jsonl(bundle / 'private/labels.test.jsonl')} <= test)
        splits = {r['case']: r['split'] for r in P5.read_jsonl(bundle / 'public/cases.jsonl')}
        self.assertEqual({splits[c] for c in cal}, {'calibration'})
        self.assertEqual({splits[c] for c in test}, {'test'})


class KeioMenu(Base):
    def keio_rows(self, records):
        path = F.write_jsonl(self.root / f'keio-{len(list(self.root.iterdir()))}.jsonl', records)
        model = F.write_sbml(self.root / 'ec.xml')
        out = self.root / f'b{len(list(self.root.iterdir()))}'
        P5.prepare(out, keio_cases=path, keio_model=model)
        return out, [r for r in P5.read_jsonl(out / 'public/cases.jsonl') if r['source'] == 'keio']

    def test_menu_is_the_full_grid_for_every_locus_and_baseline_first(self):
        _, rows = self.keio_rows(F.keio_records())
        menus = {json.dumps(r['menu'], sort_keys=True) for r in rows}
        self.assertEqual(len(menus), 1)
        menu = rows[0]['menu']
        self.assertEqual(menu[0], {'index': 0, 'mods': []})
        opened = sorted(m['mods'][1]['exchange'] for m in menu[1:])
        self.assertEqual(opened, ['R_EX_ac_e', 'R_EX_succ_e'])   # arab (background) and xyl (undeclared) excluded

    def test_menu_does_not_depend_on_the_measured_pair(self):
        records = F.keio_records()
        b1, rows1 = self.keio_rows(records)
        swapped = [dict(r, case=dict(r['case'])) for r in records]
        for r in swapped:     # b0004 measured on acetate instead of succinate; the grid is unchanged
            if r['case']['locus_tag'] == 'b0004' and r['protocol'] == 'manuscript_regulatory':
                r['case']['substrate_exchange'] = 'R_EX_ac_e'
                r['constraints'] = [F.closure('R_EX_ac_e')]
        b2, rows2 = self.keio_rows(swapped)
        self.assertEqual([r['menu'] for r in rows1], [r['menu'] for r in rows2])
        self.assertEqual((b1 / 'public/cases.jsonl').read_bytes(), (b2 / 'public/cases.jsonl').read_bytes())
        labels1 = P5.read_jsonl(b1 / 'private/labels.calibration.jsonl')
        labels2 = P5.read_jsonl(b2 / 'private/labels.calibration.jsonl')
        self.assertNotEqual(labels1, labels2)

    def test_glucose_baseline_is_unmeasured_and_labels_index_the_menu(self):
        bundle, rows = self.keio_rows(F.keio_records())
        labels = P5.read_jsonl(bundle / 'private/labels.calibration.jsonl')
        self.assertTrue(labels)
        self.assertNotIn(0, {r['menu_index'] for r in labels})
        self.assertEqual({r['wt_status'] for r in labels}, {'stated_by_source'})

    def test_keio_note_and_protocol_constants_are_recorded(self):
        _, receipt = self.prepare()
        self.assertIn('design', receipt['keio_note'])
        self.assertEqual(receipt['nutrient_oxygen'], 20.0)
        self.assertEqual(receipt['primary_variant']['nutrient'], 'printed')
        self.assertEqual(receipt['sulfonate_sulfur_ceiling'], 1.0)


class NutrientAndSulfonate(Base):
    def test_nutrient_variants_share_one_gene_independent_menu_per_species(self):
        bundle, _ = self.prepare()
        rows = [r for r in P5.read_jsonl(bundle / 'public/cases.jsonl') if r['source'] == 'nutrient']
        self.assertEqual({r['variant'] for r in rows}, {'printed', 'additional_inorganic'})
        self.assertEqual({r['primary_variant'] for r in rows if r['variant'] == 'printed'}, {True})
        for r in rows:
            self.assertEqual(r['menu'][0]['mods'], [])
            self.assertEqual(r['base']['baseline']['medium']['R_EX_o2_e'], 20.0)
            self.assertEqual(len(r['menu']), 2 * 7)      # two represented environments x (none + six supplements)

    def test_supplemented_wild_type_inherits_the_unsupplemented_control(self):
        bundle, _ = self.prepare()
        labels = [r for r in P5.read_jsonl(bundle / 'private/labels.test.jsonl') if r['source'] == 'nutrient']
        self.assertIn(P5.NUTRIENT_SUPPLEMENTED_WT, {r.get('wt_basis') for r in labels})

    def test_sulfonate_targets_share_a_menu_with_background_knockouts(self):
        bundle, _ = self.prepare()
        rows = [r for r in P5.read_jsonl(bundle / 'public/cases.jsonl') if r['source'] == 'sulfonate']
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0]['menu'], rows[1]['menu'])
        knockouts = [m for e in rows[0]['menu'] for m in e['mods'] if m['op'] == 'knockout']
        self.assertTrue(knockouts)
        self.assertTrue(all(m['genes'] == ['G_b0010'] and 'reactions' not in m for m in knockouts))
        truth = {t['case']: t for t in P5.read_jsonl(bundle / 'private/truth.test.jsonl')}
        self.assertEqual({truth[r['case']]['stratum'] for r in rows}, {'transport_only'})
        labels = [r for r in P5.read_jsonl(bundle / 'private/labels.test.jsonl') if r['source'] == 'sulfonate']
        self.assertIn('partial', {r['mutant_status'] for r in labels})


class TierBGate(Base):
    def manifest(self, columns):
        path = self.root / 'price.tsv'
        path.write_text('\t'.join(columns) + '\nx' + '\tx' * (len(columns) - 1) + '\n')
        return path

    def test_price_placeholder_and_truth_derived_columns_are_rejected(self):
        for column in sorted(P5.PRICE_FORBIDDEN_COLUMNS):
            with self.assertRaisesRegex(P5.WP5Error, column):
                P5.check_tier_b(self.manifest(['locus', column]))

    def test_price_stays_gated_even_with_clean_columns_and_a_complete_gate_file(self):
        gate = self.root / 'gate.json'
        with self.assertRaisesRegex(P5.WP5Error, 'tier_b_gated'):
            P5.check_tier_b(self.manifest(['locus', 'protein']))
        gate.write_text(json.dumps({k: True for k in P5.TIER_B_GATE}))
        with self.assertRaisesRegex(P5.WP5Error, 'tier_b_gated'):
            P5.check_tier_b(self.manifest(['locus', 'protein']), gate)
        with self.assertRaises(P5.WP5Error):
            self.prepare(price_manifest=self.manifest(['locus']))
        self.assertFalse((self.root / 'bundle').exists())


class FailClosed(Base):
    def test_existing_bundle_directory_is_refused(self):
        self.prepare()
        with self.assertRaisesRegex(P5.WP5Error, 'refusing existing'):
            self.prepare()

    def test_reaction_level_condition_mods_and_products_are_refused(self):
        with self.assertRaises(P5.WP5Error):
            P5.check_mods([{'op': 'knockout', 'reactions': ['R_PGI'], 'genes': []}])
        with self.assertRaises(P5.WP5Error):
            P5.check_mods([{'op': 'open', 'exchange': 'R_EX_ac_e', 'product_exchange': 'R_EX_ac_e'}])
        with self.assertRaises(P5.WP5Error):
            P5.check_mods([{'op': 'set_flux', 'reaction': 'R_PGI'}])


if __name__ == '__main__':
    unittest.main()
