"""WP5 step 3 (design §4, §8.3): every leak channel fails the audit; GPR- or truth-reading builders fail invariance."""
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
import audit_wp5_masking as A  # noqa: E402


def truth_reading_java():
    """A builder that reads the truth file below the Python ``open`` guard (as a JVM would) and uses it."""
    honest = F.stub_java()

    def run(out, template_path, claims_path, log):
        honest(out, template_path, claims_path, log)
        truth_path = Path(claims_path).parents[1] / 'truth.calibration.jsonl'
        fd = os.open(truth_path, os.O_RDONLY)
        try:
            truth = [json.loads(line) for line in os.read(fd, 1 << 20).decode().splitlines() if line.strip()]
        finally:
            os.close(fd)
        result_path = Path(str(out) + '.experiments') / 'result.json'
        result = P5.read_json(result_path)
        problem = result['public_problem']
        n = len(problem['hypotheses'])
        problem['hypotheses'].append({'id': f'H{n}', 'assignments': {truth[0]['reactions'][0]: [['MASKED_TARGET']]},
                                      'evidence': {'sequence': 0.1}})
        result['public_problem_sha256'] = B.canonical_hash(problem)
        P5.write_json(result_path, result)
    return run


class Base(unittest.TestCase):
    java = None

    def setUp(self):
        self.root = Path(tempfile.mkdtemp(prefix='wp5-audit-', dir=os.environ.get('TMPDIR')))
        self.addCleanup(shutil.rmtree, self.root, True)
        self.bundle, self.model = F.constructed_bundle(self.root, java=self.java)
        self.work = self.bundle / 'private/candidates'

    def audit(self, java=None, invariance=True):
        return A.audit(self.bundle, java or F.stub_java(), self.root / 'scratch', skip_invariance=not invariance)

    def assert_channel_fails(self, channel):
        report = self.audit(invariance=False)
        self.assertEqual(report['status'], 'fail')
        self.assertEqual(report['checks']['alias_scan']['status'], 'fail')
        self.assertIn(channel, report['checks']['alias_scan']['hits_by_channel'])
        public = (self.bundle / 'public/masking_audit.json').read_text()
        for alias in ('pgiX', 'b0001', 'P0A001'):
            self.assertNotIn(alias, public)          # the public report never names an alias
        detail = P5.read_json(self.bundle / 'private/audit/detail.json')
        self.assertTrue(any(h['channel'] == channel for h in detail['hits']))


class CleanAudit(Base):
    def test_target_blind_builder_passes_every_check_and_binds_the_public_manifest(self):
        report = self.audit()
        self.assertEqual({k: v['status'] for k, v in report['checks'].items()},
                         {'alias_scan': 'pass', 'objectives': 'pass', 'db': 'pass', 'invariance': 'pass'})
        self.assertEqual(report['status'], 'pass')
        self.assertEqual(report['checks']['invariance']['problems_checked'], 3)
        self.assertEqual(report['public_manifest_sha256'], P5.sha_file(self.bundle / 'public/manifest.json'))
        self.assertEqual(set(report['checks']['invariance']['variants']),
                         {'decoy_gpr', 'truth_permuted', 'labels_flipped', 'orphan'})

    def test_skipped_invariance_or_missing_builder_is_never_a_pass(self):
        self.assertEqual(self.audit(invariance=False)['status'], 'fail')
        report = A.audit(self.bundle, None, None)
        self.assertEqual(report['checks']['invariance']['status'], 'fail')

    def test_private_variants_really_permute_truth_and_flip_labels(self):
        truth = A.private_variant(self.bundle, self.root / 't', 'truth')
        labels = A.private_variant(self.bundle, self.root / 'l', 'labels')
        orig = [r['reactions'] for r in P5.read_jsonl(self.bundle / 'private/truth.calibration.jsonl')]
        perm = [r['reactions'] for r in P5.read_jsonl(truth / 'private/truth.calibration.jsonl')]
        self.assertEqual(sorted(map(tuple, orig)), sorted(map(tuple, perm)))
        self.assertTrue(all(a != b for a, b in zip(orig, perm)))
        before = [r['mutant_status'] for r in P5.read_jsonl(self.bundle / 'private/labels.calibration.jsonl')]
        after = [r['mutant_status'] for r in P5.read_jsonl(labels / 'private/labels.calibration.jsonl')]
        self.assertEqual(after, [{'growth': 'no_growth', 'no_growth': 'growth'}[x] for x in before])


class LeakChannels(Base):
    def test_alias_in_a_fasta_header_fails(self):
        path = self.work / 'query.Ec.faa'
        path.write_text(path.read_text().replace('>Q001', '>Q001 b0001', 1))
        self.assert_channel_fails('fasta_header')

    def test_alias_in_claim_metadata_fails(self):
        path = self.work / 'claims.W001.jsonl'
        rows = P5.read_jsonl(path)
        rows[0]['metadata'] = {'note': 'homolog of pgiX'}
        P5.write_jsonl(path, rows)
        self.assert_channel_fails('claims_evidence')

    def test_alias_in_a_condition_mod_fails(self):
        path = self.bundle / 'public/problems/W001-manuscript_regulatory.json'
        record = P5.read_json(path)
        record['public_problem']['conditions'].append({'id': 'C9', 'mods': [{'op': 'knockout', 'genes': ['b0001']}],
                                                       'evidence_channels': []})
        P5.write_json(path, record)
        self.assert_channel_fails('public_problem')

    def test_alias_in_the_masked_sbml_fbc_name_fails(self):
        path = self.bundle / 'public/models/W001.xml'
        text = path.read_text()
        self.assertIn('fbc:name="partX"', text)
        path.write_text(text.replace('fbc:name="partX"', 'fbc:name="pgiX"'))
        self.assert_channel_fails('masked_model')

    def masked_model_hits(self, report):
        detail = P5.read_json(self.bundle / 'private/audit/detail.json')
        return report['checks']['alias_scan'], [h for h in detail['hits'] if h['channel'] == 'masked_model'], detail

    def test_alias_in_another_gene_products_annotation_fails_with_its_gene_product_and_context(self):
        path = self.bundle / 'public/models/W001.xml'
        text = path.read_text()
        anchor = 'http://identifiers.org/uniprot/P0A004"'
        self.assertIn(anchor, text)
        path.write_text(text.replace(anchor, 'http://identifiers.org/refseq_synonym/pgiX"', 1))
        report = self.audit(invariance=False)
        scan, hits, _ = self.masked_model_hits(report)
        self.assertEqual(scan['status'], 'fail')
        self.assertEqual([(h['gene_product'], h['alias']) for h in hits], [('G_b0004', 'pgiX')])
        self.assertIn('refseq_synonym/pgiX', hits[0]['context'])
        self.assertLessEqual(len(hits[0]['context']), A.CONTEXT)
        self.assertNotIn('pgiX', (self.bundle / 'public/masking_audit.json').read_text())

    def test_alias_shared_with_another_gene_in_the_original_model_is_adjudicated_not_gating(self):
        # Real-bundle pattern: an operon/complex name (ptsX on manY/manZ) or a protein-family id carried by
        # other genes in the ORIGINAL model. It is present whether or not the target is masked.
        original = Path(P5.read_json(self.bundle / 'private/models.json')['Ec']['path'])
        anchor = 'http://identifiers.org/uniprot/P0A004"'
        for path in (original, self.bundle / 'public/models/W001.xml'):
            text = path.read_text()
            self.assertIn(anchor, text)
            path.write_text(text.replace(anchor, 'http://identifiers.org/refseq_synonym/pgiX"', 1))
        report = self.audit(invariance=False)
        scan, hits, detail = self.masked_model_hits(report)
        self.assertEqual(hits, [])
        self.assertNotIn('masked_model', scan['hits_by_channel'])
        shared = scan['informational']['masked_model_shared']
        self.assertEqual((shared['gating'], shared['hits'], shared['hits_by_case']), (False, 1, {'W001': 1}))
        [entry] = [i for i in detail['informational'] if i['channel'] == 'masked_model_shared']
        self.assertEqual([(h['gene_product'], h['alias'], h['adjudication']) for h in entry['adjudicated']],
                         [('G_b0004', 'pgiX', 'shared_identifier_in_original_model')])
        self.assertNotIn('pgiX', (self.bundle / 'public/masking_audit.json').read_text())

    def test_shared_alias_in_a_gene_product_ref_still_fails(self):
        self.assertEqual(A.split_model_hits(
            [{'element': 'geneProductRef', 'gene_product': 'G_b0004', 'alias': 'pgiX'}], {'pgiX'}, 'G_b0001'),
            ([{'element': 'geneProductRef', 'gene_product': 'G_b0004', 'alias': 'pgiX'}], []))

    def test_target_id_in_a_gene_product_ref_fails(self):
        path = self.bundle / 'public/models/W001.xml'
        text = path.read_text()
        anchor = 'fbc:geneProduct="G_b0004"'
        self.assertIn(anchor, text)
        path.write_text(text.replace(anchor, 'fbc:geneProduct="G_b0001"', 1))
        scan, hits, _ = self.masked_model_hits(self.audit(invariance=False))
        self.assertEqual(scan['status'], 'fail')
        self.assertEqual({(h['element'], h['alias'], h['reaction']) for h in hits},
                         {('geneProductRef', 'G_b0001', 'R_SUCD')})

    def test_aliases_outside_gene_records_or_in_another_case_are_informational_only(self):
        # Real-bundle collisions: BioCyc ``META:`` vs symbol metA, the word "out" in a reaction name,
        # a BiGG metabolite id equal to a symbol; and a different-case symbol in a gene name.
        path = self.bundle / 'public/models/W001.xml'
        text = path.read_text()
        text = text.replace('<reaction id="R_ACt"', '<reaction name="pgiX b0001 symport" id="R_ACt"', 1)
        text = text.replace('<species id="M_glc__D_e"', '<species name="P0A001" id="M_glc__D_e"', 1)
        self.assertIn('fbc:name="partX"', text)
        text = text.replace('fbc:name="partX"', 'fbc:name="PGIX"', 1)
        path.write_text(text)
        report = self.audit(invariance=False)
        scan, hits, detail = self.masked_model_hits(report)
        self.assertEqual(hits, [])
        self.assertNotIn('masked_model', scan['hits_by_channel'])
        self.assertEqual(scan['status'], 'pass')
        info = scan['informational']['masked_model_text']
        self.assertFalse(info['gating'])
        self.assertEqual(info['hits'], 4)
        self.assertEqual(info['hits_by_case'], {'W001': 4})
        public = (self.bundle / 'public/masking_audit.json').read_text()
        for alias in ('pgiX', 'PGIX', 'b0001', 'P0A001'):
            self.assertNotIn(alias, public)
        (row,) = detail['informational']
        self.assertEqual((row['channel'], row['case'], row['hits']), ('masked_model_text', 'W001', 4))
        self.assertEqual(sorted(row['by_alias'].items()), [('P0A001', 1), ('PGIX', 1), ('b0001', 1), ('pgiX', 1)])
        self.assertEqual(len(row['snippets']), 4)
        self.assertTrue(all(len(s['context']) <= A.CONTEXT for s in row['snippets']))

    def test_informational_snippets_are_capped(self):
        path = self.bundle / 'public/models/W001.xml'
        text = path.read_text()
        path.write_text(text.replace('<listOfSpecies>', '<listOfSpecies>' + ''.join(
            f'<species name="b0001 {i}" id="M_x{i}_c"/>' for i in range(30)), 1))
        report = self.audit(invariance=False)
        self.assertEqual(report['checks']['alias_scan']['status'], 'pass')
        (row,) = P5.read_json(self.bundle / 'private/audit/detail.json')['informational']
        self.assertEqual(row['hits'], 30)
        self.assertEqual(len(row['snippets']), A.INFORMATIONAL_SNIPPETS)

    def test_structural_scan_is_case_sensitive_except_for_numeric_ids(self):
        self.assertEqual(A.alias_pattern({'metA'}, case_sensitive=True).findall('META:CPD metA METa'), ['metA'])
        numeric = A.alias_pattern({'945289'}, case_sensitive=True)
        self.assertEqual(numeric.findall('ncbigene/945289 9452890 945289.5'), ['945289'])
        self.assertEqual(A.alias_pattern({'metA'}).findall('META:CPD'), ['META'])   # informational stays insensitive

    def test_alias_in_the_menu_file_fails(self):
        path = self.bundle / 'public/cases.jsonl'
        rows = P5.read_jsonl(path)
        rows[0]['menu'][1]['comment'] = 'designed for ackX'
        P5.write_jsonl(path, rows)
        self.assert_channel_fails('menus_constraints')

    def test_product_condition_or_reaction_knockout_fails_the_objective_check(self):
        path = self.bundle / 'public/problems/W002-manuscript_regulatory.json'
        for extra in ({'id': 'C9', 'mods': [], 'product_exchange': 'R_EX_ac_e'},
                      {'id': 'C9', 'mods': [{'op': 'knockout', 'reactions': ['R_ACK'], 'genes': ['G_b0005']}]}):
            record = P5.read_json(path)
            record['public_problem']['conditions'].append(extra)
            P5.write_json(path, record)
            report = self.audit(invariance=False)
            self.assertEqual(report['checks']['objectives']['status'], 'fail')
            record['public_problem']['conditions'].pop()
            P5.write_json(path, record)

    def test_reference_entry_of_an_excluded_genus_fails_the_db_check(self):
        with open(self.work / 'reference.Ec.tsv', 'a') as stream:
            stream.write('S00001\tEC:5.3.1.9\tescherichia\n')
        self.assertEqual(self.audit(invariance=False)['checks']['db']['status'], 'fail')

    def test_alias_in_a_qwen_request_fails_the_transcript_audit(self):
        runs = self.root / 'runs'
        P5.write_json(runs / 'protocol.json', {'schema': 'gspa.wp5.recovery/1:protocol', 'problems': {}})
        tdir = runs / 'transcripts/W001-manuscript_regulatory'
        P5.write_jsonl(tdir / 'qwen-k3-s0.jsonl', [{'request': [{'role': 'user', 'content': '{"genes": ["P1"]}'}]}])
        report = A.transcript_audit(self.bundle, runs)
        self.assertEqual(report['status'], 'pass')
        self.assertEqual(set(report['transcripts_sha256']), {'transcripts/W001-manuscript_regulatory/qwen-k3-s0.jsonl'})
        P5.write_jsonl(tdir / 'qwen-k6-s0.jsonl', [{'request': [{'role': 'user', 'content': 'is this PGIX?'}]}])
        self.assertEqual(A.transcript_audit(self.bundle, runs)['status'], 'fail')

    def test_transcript_audit_without_transcripts_is_not_a_pass(self):
        runs = self.root / 'runs'
        P5.write_json(runs / 'protocol.json', {'problems': {}})
        (runs / 'transcripts').mkdir(parents=True)
        self.assertEqual(A.transcript_audit(self.bundle, runs)['status'], 'fail')


class GprReadingBuilder(Base):
    java = staticmethod(F.stub_java(read_gpr=True))

    def test_builder_that_reads_the_target_gpr_fails_invariance(self):
        report = self.audit(java=F.stub_java(read_gpr=True))
        self.assertEqual(report['checks']['alias_scan']['status'], 'pass')
        self.assertEqual(report['checks']['invariance']['status'], 'fail')
        failures = P5.read_json(self.bundle / 'private/audit/detail.json')['invariance']
        self.assertTrue({'decoy_gpr', 'orphan'} <= {f['variant'] for f in failures})
        self.assertFalse({'truth_permuted', 'labels_flipped'} & {f['variant'] for f in failures})
        self.assertEqual(report['status'], 'fail')


class TruthReadingBuilder(Base):
    java = staticmethod(truth_reading_java())

    def test_builder_that_reads_truth_outside_the_python_guard_fails_on_the_permuted_variant(self):
        report = self.audit(java=truth_reading_java())
        self.assertEqual(report['checks']['invariance']['status'], 'fail')
        failures = P5.read_json(self.bundle / 'private/audit/detail.json')['invariance']
        self.assertIn('truth_permuted', {f['variant'] for f in failures})

    def test_python_reads_of_truth_are_refused_by_the_guard(self):
        guard = A.private_files(self.bundle)
        with A.forbid_open(guard):
            for path in guard:
                with self.assertRaises(PermissionError):
                    path.read_text()
                with self.assertRaises(PermissionError):
                    open(path).close()
            (self.bundle / 'public/cases.jsonl').read_text()


if __name__ == '__main__':
    unittest.main()
