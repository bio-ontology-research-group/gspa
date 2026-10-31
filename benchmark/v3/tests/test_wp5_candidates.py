"""WP5 step 2 (design §4, §8.2; ruling §13.4): header stripping, DB exclusion, no GPR reads, construction budget."""
import json
import os
import shutil
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import wp5_fixtures as F  # noqa: E402
import prepare_wp5_cases as P5  # noqa: E402
import build_wp5_candidates as B  # noqa: E402


class Base(unittest.TestCase):
    def setUp(self):
        self.root = Path(tempfile.mkdtemp(prefix='wp5-cand-', dir=os.environ.get('TMPDIR')))
        self.addCleanup(shutil.rmtree, self.root, True)
        self.bundle, self.model = F.keio_bundle(self.root)
        self.paths = B.bundle_paths(self.bundle)
        self.targets = P5.read_jsonl(self.paths['targets'])


class HeaderStripping(Base):
    def test_query_fasta_has_only_opaque_ids_and_header_text_goes_to_private_aliases(self):
        B.build_queries({'Ec': F.write_proteome(self.root / 'ec.faa')}, self.targets, self.paths['work'])
        text = (self.paths['work'] / 'query.Ec.faa').read_text()
        headers = [line for line in text.splitlines() if line.startswith('>')]
        self.assertEqual(headers, ['>Q001', '>Q002', '>Q003'])
        for leak in ('b0001', 'pgiX', 'NP_000001', 'P0A001', 'isomerase', 'locus_tag', '*'):
            self.assertNotIn(leak, text)
        extra = {r['case']: r for r in P5.read_jsonl(self.paths['work'] / 'aliases_extra.jsonl')}
        self.assertIn('NP_000001.1', extra['W001']['aliases'])
        self.assertIn('NP_000001', extra['W001']['aliases'])
        self.assertIn('P0A001', extra['W001']['aliases'])
        self.assertEqual(extra['W001']['names'], ['glucose-6-phosphate isomerase X'])

    def test_target_matching_no_or_several_records_fails_closed(self):
        path = self.root / 'bad.faa'
        path.write_text('>x [locus_tag=b0001]\nMKK\n>y [locus_tag=b0001]\nMKA\n')
        with self.assertRaisesRegex(P5.WP5Error, 'matched 2'):
            B.build_queries({'Ec': path}, self.targets, self.paths['work'])
        with self.assertRaisesRegex(P5.WP5Error, 'no proteome'):
            B.build_queries({}, self.targets, self.paths['work'])

    def test_locus_tags_match_after_underscore_and_case_normalisation(self):
        target = dict(self.targets[0], aliases=['BSU00930'])
        path = self.root / 'bs.faa'
        path.write_text('>lcl|x [locus_tag=BSU_00930] [protein=foo]\nMKK*\n>lcl|y [locus_tag=BSU_00931]\nMKA\n')
        extra = B.build_queries({target['model_key']: path}, [target], self.paths['work'])
        self.assertEqual(len(extra), 1)
        self.assertIn('BSU_00930', extra[0]['aliases'])
        path.write_text('>x [locus_tag=BSU_00930]\nMKK\n>y [locus_tag=bsu00930]\nMKA\n')
        with self.assertRaisesRegex(P5.WP5Error, 'matched 2'):
            B.build_queries({target['model_key']: path}, [target], self.paths['work'])


class DatabaseExclusion(Base):
    def test_reference_drops_target_accessions_and_excluded_genera(self):
        counts = B.build_reference(None, self.targets, self.paths['work'], entries=F.SWISSPROT)
        self.assertEqual(counts['Ec'], {'read': 6, 'dropped_accession': 1, 'dropped_genus': 1, 'no_ec': 1, 'kept': 3})
        faa = (self.paths['work'] / 'reference.Ec.faa').read_text()
        self.assertNotIn('P0A001', faa)
        self.assertNotIn('Q9SHG1', faa)
        self.assertEqual([l for l in faa.splitlines() if l.startswith('>')], ['>A00001', '>A00002', '>A00004'])
        receipt = P5.read_json(self.paths['work'] / 'reference.receipt.json')
        self.assertEqual(receipt['excluded']['Ec']['genera'], ['escherichia', 'shigella'])

    def test_hit_to_an_excluded_accession_or_genus_or_unknown_subject_fails_closed(self):
        B.build_reference(None, self.targets, self.paths['work'], entries=F.SWISSPROT)
        tsv = self.paths['work'] / 'reference.Ec.tsv'
        hits = self.root / 'hits.tsv'
        for subject, message in (('sp|P0A001|PGI_ECOLI', 'DB exclusion'), ('Z99999', 'outside the frozen reference')):
            hits.write_text(f'Q001\t{subject}\t90\t100\t1e-30\t200\n')
            with self.assertRaisesRegex(P5.WP5Error, message):
                B.claims_from_hits(hits, tsv, {'P0A001'}, {'escherichia', 'shigella'})
        with open(tsv, 'a') as stream:   # a planted same-genus reference row is caught at claim time
            stream.write('S00001\tEC:5.3.1.9\tshigella\n')
        hits.write_text('Q001\tS00001\t90\t100\t1e-30\t200\n')
        with self.assertRaisesRegex(P5.WP5Error, 'DB exclusion'):
            B.claims_from_hits(hits, tsv, {'P0A001'}, {'escherichia', 'shigella'})

    def test_claims_are_sequence_only_with_identity_times_coverage(self):
        F.candidate_stages(self.bundle, self.root)
        claims = P5.read_jsonl(self.paths['work'] / 'claims.W001.jsonl')
        self.assertEqual(claims, [{'evidence_type': 'SEQUENCE_SIMILARITY', 'function_id': '5.3.1.9', 'function_type': 'EC',
                                   'metadata': {}, 'protein_id': 'Q001', 'raw_score': 0.8, 'source': B.CLAIM_SOURCE}])
        self.assertEqual(P5.read_jsonl(self.paths['work'] / 'claims.W002.jsonl')[0]['raw_score'], 0.54)

    def test_diamond_receipt_records_commands_parameters_and_db_hash(self):
        calls = []
        B.build_reference(None, self.targets, self.paths['work'], entries=F.SWISSPROT)
        B.build_queries({'Ec': F.write_proteome(self.root / 'ec.faa')}, self.targets, self.paths['work'])
        B.run_diamond(self.paths['work'], ['Ec'], runner=F.fake_diamond(calls=calls), diamond='/opt/diamond')
        self.assertEqual([c[:2] for c in calls], [['/opt/diamond', 'makedb'], ['/opt/diamond', 'blastp'],
                                                  ['/opt/diamond', 'version']])
        receipt = P5.read_json(self.paths['work'] / 'diamond.receipt.json')
        self.assertIn('reference.Ec.dmnd', receipt['outputs'])
        self.assertEqual(receipt['parameters'], B.DIAMOND)

    def test_changed_stage_output_is_refused_downstream(self):
        F.candidate_stages(self.bundle, self.root)
        (self.paths['work'] / 'hits.Ec.tsv').write_text('tampered\n')
        with self.assertRaisesRegex(P5.WP5Error, 'changed after its receipt'):
            B.write_claims(self.bundle)


class NoGprReads(Base):
    def test_reference_query_diamond_and_claims_never_open_the_model_truth_or_labels(self):
        model = Path(P5.read_json(self.paths['models'])['Ec']['path']).resolve()
        forbidden = {model} | {(self.bundle / f'private/{k}.{s}.jsonl').resolve()
                               for k in ('truth', 'labels') for s in ('calibration', 'test')}
        with F.refuse_paths(lambda p: p in forbidden) as opened:
            F.candidate_stages(self.bundle, self.root)
        self.assertFalse(forbidden & set(opened))

    def test_template_carries_no_gene_rule_and_no_gene_annotation(self):
        row = P5.read_jsonl(self.paths['cases'])[0]
        target = {t['case']: t for t in self.targets}[row['case']]
        template = B.template_for(row, target, self.model)
        self.assertEqual(set(template), {'model', 'baseline', 'objectives', 'condition_constraints',
                                         'assignment_recovery', 'condition_menu'})
        self.assertEqual(template['condition_menu'], {})
        self.assertEqual(set(template['assignment_recovery']), {'mask_gene', 'candidate_protein', 'conditions'})
        # The mask gene is the only target identity, and only the masking boundary consumes it.
        rest = json.loads(json.dumps(template))
        self.assertEqual(rest['assignment_recovery'].pop('mask_gene'), 'G_b0001')
        text = json.dumps(rest)
        for leak in ('R_PGI', 'pgiX', 'b0001', 'geneProduct', 'stratum', 'truth'):
            self.assertNotIn(leak, text)
        self.assertEqual(len(template['assignment_recovery']['conditions']), len(row['menu']) - 1)

    def test_construct_one_never_reads_truth_or_labels(self):
        F.candidate_stages(self.bundle, self.root)
        row = P5.read_jsonl(self.paths['cases'])[0]
        target = {t['case']: t for t in self.targets}[row['case']]
        guard = {(self.bundle / f'private/{k}.{s}.jsonl').resolve() for k in ('truth', 'labels')
                 for s in ('calibration', 'test')}
        with F.refuse_paths(lambda p: p in guard):
            status, problem, report = B.construct_one(row, target, self.model, self.paths['work'] / 'claims.W001.jsonl',
                                                      self.root / 'w', F.stub_java())
        self.assertEqual(status, 'constructed')
        self.assertEqual(problem['hypotheses'][1]['assignments'], {'R_PGI': [['MASKED_TARGET']]})

    def test_mask_sbml_removes_target_and_and_clauses_keeps_orphans_and_notes_go(self):
        out = self.root / 'masked.xml'
        B.mask_sbml(self.model, out, 'G_b0002')
        index = P5.index_model(out)
        self.assertNotIn('G_b0002', index['genes'])
        self.assertIn('R_ACK', index['reactions'])                 # orphaned reaction stays in the network
        self.assertEqual(index['reactions']['R_ACK']['dnf'], [])   # whole AND clause (with partner) removed
        self.assertEqual(index['reactions']['R_SUCD']['dnf'], [['G_b0004'], ['G_b0005']])
        text = out.read_text()
        self.assertNotIn('ackX', text)
        self.assertNotIn('free text', text)
        with self.assertRaisesRegex(P5.WP5Error, 'not a declared geneProduct'):
            B.mask_sbml(self.model, out, 'G_nope')

    def test_mask_sbml_keeps_repeated_or_clauses_so_the_production_fingerprint_is_reproduced(self):
        # iML1515 R_3NTD4pp is (G_b4213 or G_b0383 or G_b0383). Production mask_model keeps the repeated
        # clause and fingerprint() hashes the clause multiset; de-duplicating it failed WP5 W002 on the cluster.
        sys.path.insert(0, str(B.ROOT / 'scripts'))
        import metabolic_experiments as ME
        import masked_assignment_experiments as MAE
        reactions = [dict(r, gpr=[['G_b0005'], ['G_b0004'], ['G_b0004']]) if r['id'] == 'R_SUCD' else r
                     for r in F.EC_REACTIONS]
        original = self.root / 'dup.xml'
        F.write_sbml(original, reactions=reactions)
        for gene in ('G_b0001', 'G_b0005'):
            production = MAE.mask_model(ME.load_sbml_model(original), gene)
            self.assertEqual(production.gpr['R_SUCD'].count(['G_b0004']), 2)
            out = self.root / f'dup-masked-{gene}.xml'
            B.mask_sbml(original, out, gene)
            self.assertEqual(P5.index_model(out)['reactions']['R_SUCD']['dnf'].count(['G_b0004']), 2)
            B.verify_fingerprint(out, production.path)


class SbmlNamespaces(Base):
    """The masked SBML keeps the core namespace as default (the JVM reader is namespace-unaware)."""

    def assert_jvm_readable(self, original, out):
        from xml.dom import minidom
        text = out.read_text()
        self.assertTrue(text.split('?>', 1)[1].lstrip().startswith(f'<sbml xmlns="{F.SBML_CORE}"'))
        self.assertNotIn('<sbml:', text)
        self.assertNotIn('<ns0:', text)
        for tag in ('<model ', '<listOfReactions>', '<reaction ', '<species ', '<fbc:geneProduct ',
                    '<fbc:geneProductRef ', '<fbc:geneProductAssociation>', 'fbc:id="G_b0002"'):
            self.assertIn(tag, text)
        # SbmlModelInfo.groovy: namespace-unaware DOM, getElementsByTagName('reaction') etc.
        dom, ref = minidom.parse(str(out)), minidom.parse(str(original))
        for tag in ('sbml', 'model', 'reaction', 'species', 'fbc:geneProduct', 'fbc:listOfObjectives'):
            self.assertEqual(len(dom.getElementsByTagName(tag)), len(ref.getElementsByTagName(tag)) - (
                1 if tag == 'fbc:geneProduct' else 0), tag)
        return text

    def test_core_namespace_declared_default_and_prefixed_in_annotations_stays_default(self):
        original = F.write_cobra_style_sbml(self.root / 'cobra.xml')
        source = original.read_text()
        self.assertIn('<sbml:annotation xmlns:sbml=', source)
        out = self.root / 'cobra-masked.xml'
        B.mask_sbml(original, out, 'G_b0001')
        text = self.assert_jvm_readable(original, out)
        self.assertIn('<rdf:RDF', text)
        self.assertIn('<bqbiol:is', text)
        self.assertIn('<annotation>', text)                 # same core namespace, now via the default
        plain = F.write_sbml(self.root / 'plain.xml')
        plain_out = self.root / 'plain-masked.xml'
        B.mask_sbml(plain, plain_out, 'G_b0001')
        self.assertEqual(P5.index_model(out), P5.index_model(plain_out))
        # Namespace registration is process-global: masking again (and in the other order) is stable.
        again = self.root / 'cobra-masked-again.xml'
        B.mask_sbml(original, again, 'G_b0001')
        self.assertEqual(again.read_bytes(), out.read_bytes())
        self.assert_jvm_readable(plain, plain_out)

    def test_variants_written_for_the_invariance_audit_are_jvm_readable(self):
        original = F.write_cobra_style_sbml(self.root / 'cobra.xml')
        for name, kwargs in (('decoy', {'decoy': {'R_DECOY1'}, 'keep_target_product': True}),
                             ('orphan', {'keep_target_product': True})):
            out = self.root / f'{name}.xml'
            B.mask_sbml(original, out, 'G_b0001', **kwargs)
            text = out.read_text()
            self.assertNotIn('<sbml:', text)
            from xml.dom import minidom
            self.assertEqual(len(minidom.parse(str(out)).getElementsByTagName('reaction')),
                             len(minidom.parse(str(original)).getElementsByTagName('reaction')))

    def test_production_fingerprint_is_reproduced_for_the_cobra_style_fixture(self):
        sys.path.insert(0, str(B.ROOT / 'scripts'))
        import metabolic_experiments as ME
        import masked_assignment_experiments as MAE
        original = F.write_cobra_style_sbml(self.root / 'cobra.xml')
        out = self.root / 'cobra-masked.xml'
        B.mask_sbml(original, out, 'G_b0004')
        B.verify_fingerprint(out, MAE.mask_model(ME.load_sbml_model(original), 'G_b0004').path)


class ConstructionBudget(Base):
    """Ruling §13.4: one experiment and one solver call, recorded; at most one cell accepted."""

    def result(self, cells, calls=None, budget=1, experiments=1):
        return {'max_solver_calls': budget, 'config': {'max_solver_calls': budget, 'max_experiments': experiments},
                'cells': [{}] * cells, 'solver_calls': [{}] * (cells if calls is None else calls),
                'visited_conditions': ['C0'] * min(cells, 1)}

    def test_java_command_requests_the_one_experiment_one_call_budget(self):
        cmd = B.java_command('gspa.jar', 'py', 'side.py', 'c.jsonl', 't.json', 'o.tsv', java='/jdk/bin/java')
        self.assertEqual(cmd[0], '/jdk/bin/java')
        self.assertEqual(cmd[cmd.index('--experiment-max-experiments') + 1], '1')
        self.assertEqual(cmd[cmd.index('--experiment-max-solver-calls') + 1], '1')
        self.assertEqual(cmd[cmd.index('--experiment-max-hypotheses') + 1], '8')
        self.assertEqual(cmd[cmd.index('--experiment-mode') + 1], 'evidence_guided')

    def test_validation_accepts_zero_or_one_cell_and_rejects_more_or_a_different_budget(self):
        self.assertEqual(B.validate_construction_result(self.result(0))['cells'], 0)
        self.assertEqual(B.validate_construction_result(self.result(1))['cells'], 1)
        for bad in (self.result(2), self.result(0, calls=2), self.result(0, budget=2), self.result(0, experiments=3)):
            with self.assertRaises(P5.WP5Error):
                B.validate_construction_result(bad)
        with self.assertRaises(P5.WP5Error):
            B.validate_construction_result({'max_solver_calls': 1})

    def test_construct_records_the_budget_in_the_receipt_and_every_problem(self):
        bundle, _ = F.constructed_bundle(self.root / 'c')
        paths = B.bundle_paths(bundle)
        receipt = P5.read_json(paths['work'] / 'construct.receipt.json')
        self.assertEqual(receipt['jvm_budget'], {'max_experiments': 1, 'max_solver_calls': 1})
        self.assertEqual(set(receipt['statuses'].values()), {'constructed'})
        self.assertEqual({v['cells'] for v in receipt['construction_usage'].values()}, {0})
        for path in paths['problems'].glob('*.json'):
            self.assertEqual(P5.read_json(path)['construction_budget'], B.JVM_BUDGET)

    def test_construction_with_one_cell_is_accepted_and_two_cells_stop_the_stage(self):
        F.candidate_stages(self.bundle, self.root)
        B.construct(self.bundle, F.stub_java(cells=1), verify=None)
        with self.assertRaisesRegex(P5.WP5Error, 'budget'):
            B.construct(self.bundle, F.stub_java(cells=2, budget=2), verify=None)


class PublicProblem(Base):
    def test_menu_map_covers_the_frozen_menu_and_products_bgc_and_extra_fields_are_refused(self):
        bundle, _ = F.constructed_bundle(self.root / 'c')
        record = P5.read_json(B.bundle_paths(bundle)['problems'] / 'W001-manuscript_regulatory.json')
        self.assertEqual(record['menu_map']['0'], 'C0')
        self.assertEqual(len(record['menu_map']), len(P5.read_jsonl(B.bundle_paths(bundle)['cases'])[0]['menu']))
        problem = record['public_problem']
        for mutate in (lambda p: p['conditions'].append({'id': 'C9', 'mods': [], 'product_exchange': 'R_EX_ac_e'}),
                       lambda p: p['hypotheses'][1]['evidence'].update(bgc=0.5),
                       lambda p: p.update(truth=['R_PGI']),
                       lambda p: p['baseline'].update(product_demands=['DM_x'])):
            bad = json.loads(json.dumps(problem))
            mutate(bad)
            with self.assertRaises(P5.WP5Error):
                B.check_public(bad)

    def test_mapping_sources_other_than_sbml_annotations_are_refused(self):
        B.check_mapping_sources({'hypotheses': [{'candidates': [{'mapping_sources': ['sbml_annotation']}]}]})
        with self.assertRaises(P5.WP5Error):
            B.check_mapping_sources({'hypotheses': [{'candidates': [{'mapping_sources': ['reaction_graph_id']}]}]})

    def test_no_admissible_candidate_is_recorded_not_dropped(self):
        F.candidate_stages(self.bundle, self.root, hits=[])
        summary = B.construct(self.bundle, F.stub_java(), verify=None)
        self.assertEqual(set(summary.values()), {'no_admissible_experiment'})
        record = P5.read_json(self.paths['problems'] / 'W001-manuscript_regulatory.json')
        self.assertNotIn('public_problem', record)


if __name__ == '__main__':
    unittest.main()
