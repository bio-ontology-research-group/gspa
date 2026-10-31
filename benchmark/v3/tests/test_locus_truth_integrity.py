"""Incomplete candidate/label sets must never manufacture null-outcome truth."""
import json
from pathlib import Path
import tempfile
import unittest
import importlib.util

from test_prepare_training_data import P, standard, rows_of, tsv


class LocusTruthIntegrity(unittest.TestCase):
    def setup_case(self, path, members=('p1', 'p2'), labels=None):
        policy = {'locus_truth_complete': 'synthetic exhaustive assay panel, fixture only',
                  'datasets': {'ref': {'EC': {'allowed_evidence': 'experimental',
                     'explicit_negative_truth': 'synthetic direct assay, fixture only'}}}}
        fx, kw = standard(path, policy=policy, labels=labels,
                          ec_rows=[('p1', 'EC', 'EC:1.1.1.1', 'CATALYTIC'),
                                   ('p2', 'EC', 'EC:1.1.1.1', 'CATALYTIC')])
        kw['locus_groups'] = fx.write('locus.tsv', tsv(
            ('genome_id', 'locus_group_id', 'candidate_key', 'truth_complete'),
            [(g, 'same_name', f'{p}|EC|EC:1.1.1.1|CATALYTIC', 'true')
             for g in ('g1', 'g2', 'g3') for p in members]))
        return fx, kw

    def test_unknown_member_disqualifies_the_whole_group_not_the_marginal_rows(self):
        with tempfile.TemporaryDirectory() as d:
            _, kw = self.setup_case(Path(d))
            result = P.build_table(**kw)
            self.assertEqual(len(rows_of(result)), 3)
            self.assertTrue(all(r['locus_group_id'] is None for r in rows_of(result)))
            self.assertEqual(result['locus_audit']['eligible_groups'], 0)
            self.assertTrue(all('members_not_labelled_or_masked' in g['reasons'] for g in result['locus_audit']['groups']))

    def test_omitting_an_unlabelled_candidate_from_the_declaration_does_not_make_truth_complete(self):
        with tempfile.TemporaryDirectory() as d:
            _, kw = self.setup_case(Path(d), members=('p1',))
            result = P.build_table(**kw)
            self.assertTrue(all(r['locus_group_id'] is None for r in rows_of(result)))
            self.assertTrue(all('candidate_universe_incomplete' in g['reasons'] for g in result['locus_audit']['groups']))

    def explicit_labels(self, fx, kw, second_label='0'):
        kw['labels_path'] = fx.write('explicit-labels.tsv', tsv(
            ('genome_id', 'protein_id', 'function_id', 'evidence', 'label', 'reference'),
            [(g, p, '1.1.1.1', 'IDA', ('1' if p == 'p1' else second_label), 'fixture:direct-assay')
             for g in ('g1', 'g2', 'g3') for p in ('p1', 'p2')]))

    def test_explicit_negative_assays_are_not_inverted_to_positives_and_complete_groups_are_stamped(self):
        with tempfile.TemporaryDirectory() as d:
            fx, kw = self.setup_case(Path(d))
            self.explicit_labels(fx, kw)
            result = P.build_table(**kw)
            rows = rows_of(result)
            self.assertEqual([r['label'] for r in rows], [1, 0] * 3)
            self.assertEqual(result['locus_audit']['eligible_groups'], 3)
            for r in rows:
                self.assertEqual(len(r['locus_members']), 2)
                self.assertIn(r['example_id'], r['locus_members'])
                self.assertTrue(r['locus_truth_source'])
            self.assertTrue(all(r['label_basis'] == 'explicit_negative_truth' for r in rows if r['label'] == 0))
            spec = importlib.util.spec_from_file_location('locus_trainer', Path(__file__).resolve().parents[3] / 'scripts/train_context_model.py')
            trainer = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(trainer)
            certified = trainer.read_table(Path(result['table']), trainer.load_contract(kw['contract_path']))
            trainer.validate_locus_truth(certified)

    def test_conflicting_positive_and_negative_assay_records_fail(self):
        with tempfile.TemporaryDirectory() as d:
            fx, kw = self.setup_case(Path(d))
            self.explicit_labels(fx, kw)
            with kw['labels_path'].open('a') as fh:
                fh.write('g1\tp1\t1.1.1.1\tIDA\t0\tfixture:conflict\n')
            with self.assertRaisesRegex(P.PreparationError, 'conflicting'):
                P.build_table(**kw)

    def test_duplicate_or_conflicting_group_members_do_not_silently_overwrite(self):
        with tempfile.TemporaryDirectory() as d:
            fx, kw = self.setup_case(Path(d))
            with kw['locus_groups'].open('a') as fh:
                fh.write('g1\tdifferent_group\tp1|EC|EC:1.1.1.1|CATALYTIC\ttrue\n')
            with self.assertRaisesRegex(P.PreparationError, 'duplicate|multiple'):
                P.build_table(**kw)

    def test_closed_world_unconfirmed_negatives_never_supervise_null_even_when_declared_complete(self):
        with tempfile.TemporaryDirectory() as d:
            labels = [(g, p, '2.2.2.2', 'IDA') for g in ('g1', 'g2', 'g3') for p in ('p1', 'p2')]
            _, kw = self.setup_case(Path(d), labels=labels)
            policy = json.loads(kw['policy_path'].read_text())
            policy['datasets']['ref']['EC'].update(negative_policy='closed_world_declared', justification='fixture only')
            kw['policy_path'].write_text(json.dumps(policy))
            result = P.build_table(**kw)
            self.assertTrue(all(r['label'] == 0 and r['locus_group_id'] is None for r in rows_of(result)))
            self.assertTrue(all('unconfirmed_negative_labels' in g['reasons'] for g in result['locus_audit']['groups']))

    def test_negative_truth_requires_a_source_and_does_not_label_other_functions(self):
        for missing in ('policy', 'reference'):
            with self.subTest(missing=missing), tempfile.TemporaryDirectory() as d:
                fx, kw = self.setup_case(Path(d))
                self.explicit_labels(fx, kw)
                if missing == 'policy':
                    policy = json.loads(kw['policy_path'].read_text())
                    del policy['datasets']['ref']['EC']['explicit_negative_truth']
                    kw['policy_path'].write_text(json.dumps(policy))
                else:
                    kw['labels_path'].write_text(kw['labels_path'].read_text().replace('fixture:direct-assay', ''))
                with self.assertRaisesRegex(P.PreparationError, 'explicit negative truth requires'):
                    P.build_table(**kw)

    def test_family_masking_a_positive_does_not_turn_surviving_negative_into_null(self):
        with tempfile.TemporaryDirectory() as d:
            fx, kw = self.setup_case(Path(d))
            self.explicit_labels(fx, kw)
            # Add unrelated positives so ordinary family masking still has usable
            # marginal examples after removing the target positive from train.
            for g in ('g1', 'g2', 'g3'):
                folder = fx.feature_dir / g
                for p in folder.iterdir():
                    p.unlink()
                folder.rmdir()
                fx.export(g, [('p1', 'EC', 'EC:1.1.1.1', 'CATALYTIC'),
                              ('p2', 'EC', 'EC:1.1.1.1', 'CATALYTIC'),
                              ('p3', 'EC', 'EC:9.9.9.9', 'CATALYTIC')])
            with kw['labels_path'].open('a') as fh:
                for g in ('g1', 'g2', 'g3'):
                    fh.write(f'{g}\tp3\t9.9.9.9\tIDA\t1\tfixture:other-assay\n')
            kw['example_clusters'] = fx.write('clusters.tsv', tsv(('genome_id', 'protein_id', 'cluster'),
                [(g, p, 'SHARED' if p == 'p1' and g in ('g1', 'g3') else g + p)
                 for g in ('g1', 'g2', 'g3') for p in ('p1', 'p2', 'p3')]))
            result = P.build_table(**kw)
            negative = next(r for r in rows_of(result) if r['genome_id'] == 'g1' and r['protein_id'] == 'p2')
            self.assertEqual(negative['label'], 0)
            self.assertIsNone(negative['locus_group_id'])
            audit = next(g for g in result['locus_audit']['groups'] if g['genome_id'] == 'g1')
            self.assertIn('members_not_labelled_or_masked', audit['reasons'])

    def test_negative_go_annotation_is_exact_and_never_propagated_to_ancestors(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d)
            policy = {'datasets': {'ref': {'GO': {'allowed_evidence': 'experimental', 'explicit_negative_truth': 'fixture assay'}}}}
            fx, kw = standard(path, policy=policy, ec_rows=[('p1', 'GO', 'GO:0000001', 'UNSPECIFIED'),
                                                            ('p1', 'GO', 'GO:0000002', 'UNSPECIFIED')])
            kw['go_obo'] = fx.write('go.obo', '[Term]\nid: GO:0000002\nis_a: GO:0000001\n\n[Term]\nid: GO:0000001\n')
            kw['labels_path'] = fx.write('go-labels.tsv', tsv(('genome_id', 'protein_id', 'function_id', 'evidence', 'qualifier', 'reference'),
                 [(g, 'p1', 'GO:0000002', 'IDA', 'NOT', 'fixture:negative-assay') for g in ('g1', 'g2', 'g3')]))
            rows = rows_of(P.build_table(**kw))
            self.assertEqual(len(rows), 3)
            self.assertTrue(all(r['function_id'] == 'GO:0000002' and r['label'] == 0 for r in rows))

    def test_mixed_roles_or_polarity_and_boolean_completeness_are_rejected(self):
        for kind in ('role', 'complete', 'source'):
            with self.subTest(kind=kind), tempfile.TemporaryDirectory() as d:
                fx, kw = self.setup_case(Path(d))
                if kind == 'role':
                    kw['locus_groups'].write_text(kw['locus_groups'].read_text().replace('p2|EC|EC:1.1.1.1|CATALYTIC', 'p2|EC|EC:1.1.1.1|ACCESSORY'))
                elif kind == 'complete':
                    kw['locus_groups'].write_text(kw['locus_groups'].read_text().replace('p2|EC|EC:1.1.1.1|CATALYTIC\ttrue', 'p2|EC|EC:1.1.1.1|CATALYTIC\tfalse'))
                else:
                    p = json.loads(kw['policy_path'].read_text())
                    p['locus_truth_complete'] = True
                    kw['policy_path'].write_text(json.dumps(p))
                with self.assertRaises(P.PreparationError):
                    P.build_table(**kw)


if __name__ == '__main__':
    unittest.main()
