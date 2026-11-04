#!/usr/bin/env python3
"""Unit tests for benchmark/v3/row5_masked_context (plan tasks/v3-row5-context-improvement-2026-11-02.md)."""
from __future__ import annotations

import csv
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
V3 = HERE.parent
if str(V3) not in sys.path:
    sys.path.insert(0, str(V3))

from price_calibration_context.build_context_inputs import observation  # noqa: E402
from row5_masked_context import evaluate, registry, run_cell  # noqa: E402

TMP_BASE = Path(os.environ.get('TMPDIR', '/home/leechuck/tmp'))
TMP_BASE.mkdir(parents=True, exist_ok=True)


def _make_temp_dir() -> tempfile.TemporaryDirectory:
    return tempfile.TemporaryDirectory(dir=str(TMP_BASE))


class TestRow5MaskedContextRegistry(unittest.TestCase):
    def _write_synthetic_tables(self, g_dir: Path, genome: str) -> tuple[Path, Path]:
        g_dir.mkdir(parents=True, exist_ok=True)
        pw_tbl = g_dir / f'{genome}-all-Pathways.tbl'
        rx_tbl = g_dir / f'{genome}-all-Reactions.tbl'

        pw_rows = [
            # ID, Name, Prediction, Completeness, Vague, KeyReactions, NrReaction, NrReactionFound
            ['|PWY-MASK|', 'maskable pathway', 'TRUE', '100', '0', 'RXN-OTHER1', '4', '4'],
            ['|PWY-ALRGAP|', 'already gap pathway', 'TRUE', '80', '0', 'RXN-OTHER2', '5', '4'],
            ['|PWY-REDUND|', 'redundant pathway', 'TRUE', '100', '0', '', '3', '3'],
            ['|PWY-NOTPRED|', 'not predicted pathway', 'FALSE', '25', '0', '', '4', '1'],
            ['|PWY-EXCL|', 'excluded rule pathway', 'TRUE', '100', '0', '', '3', '3'],
            ['|PWY-NOCTX|', 'maskable no context pathway', 'TRUE', '75', '0', '', '4', '3'],
        ]
        with pw_tbl.open('w', encoding='utf-8', newline='') as fh:
            w = csv.writer(fh, delimiter='\t')
            w.writerow(['ID', 'Name', 'Prediction', 'Completeness', 'Vague', 'KeyReactions', 'NrReaction', 'NrReactionFound'])
            w.writerows(pw_rows)

        rx_rows = [
            # 1) maskable_with_context (maskable_gap): P_MASK is the sole good_blast on RXN-MASK in PWY-MASK,
            #    while PWY-MASK also has RXN-OTHER1 supported by unlabelled protein P_UNLAB1
            ['RXN-MASK', 'maskable rxn', '1.1.1.1', 'FALSE', 'FALSE', 'P_MASK desc', 'good_blast', '|PWY-MASK|', 'FALSE', '250', '95', '80', '1e-60', 'U1', 'Q1', '100', '500'],
            ['RXN-MASK', 'maskable rxn', '1.1.1.1', 'FALSE', 'FALSE', 'P_WEAK desc', 'bad_blast', '|PWY-MASK|', 'FALSE', '80', '70', '35', '1e-10', 'U2', 'Q2', '600', '900'],
            ['RXN-OTHER1', 'context step 1', '1.1.1.99', 'TRUE', 'FALSE', 'P_UNLAB1 desc', 'good_blast', '|PWY-MASK|', 'FALSE', '300', '98', '85', '1e-70', 'U3', 'Q3', '1000', '1500'],
            ['RXN-OTHER1B', 'context step 2', '1.1.1.98', 'FALSE', 'FALSE', 'P_UNLAB2 desc', 'good_blast', '|PWY-MASK|', 'FALSE', '290', '97', '82', '1e-68', 'U4', 'Q4', '1600', '2100'],

            # 2) already_gap (also maskable_with_context, gap_kind='already_gap'):
            #    no good_blast row on RXN-ALRGAP, while PWY-ALRGAP has other good_blast steps
            ['RXN-ALRGAP', 'already gap rxn', '2.2.2.2', 'FALSE', 'FALSE', 'P_ALRGAP desc', 'bad_blast', '|PWY-ALRGAP|', 'FALSE', '90', '75', '38', '1e-12', 'U5', 'Q5', '2200', '2600'],
            ['RXN-OTHER2', 'context step 2', '2.2.2.99', 'TRUE', 'FALSE', 'P_UNLAB3 desc', 'good_blast', '|PWY-ALRGAP|', 'FALSE', '310', '99', '88', '1e-75', 'U6', 'Q6', '2700', '3200'],

            # 3) redundant: both P_REDUND and P_PARALOG have good_blast on RXN-REDUND
            ['RXN-REDUND', 'redundant rxn', '3.3.3.3', 'FALSE', 'FALSE', 'P_REDUND desc', 'good_blast', '|PWY-REDUND|', 'FALSE', '260', '96', '81', '1e-62', 'U7', 'Q7', '3300', '3800'],
            ['RXN-REDUND', 'redundant rxn', '3.3.3.3', 'FALSE', 'FALSE', 'P_PARALOG desc', 'good_blast', '|PWY-REDUND|', 'FALSE', '255', '95', '80', '1e-61', 'U8', 'Q8', '3900', '4400'],
            ['RXN-REDUND-CTX', 'ctx', '3.3.3.99', 'FALSE', 'FALSE', 'P_UNLAB4 desc', 'good_blast', '|PWY-REDUND|', 'FALSE', '270', '96', '83', '1e-64', 'U9', 'Q9', '4500', '5000'],

            # 4) pathway_not_predicted: PWY-NOTPRED has Prediction=FALSE
            ['RXN-NOTPRED', 'not pred rxn', '4.4.4.4', 'FALSE', 'FALSE', 'P_NOTPRED desc', 'good_blast', '|PWY-NOTPRED|', 'FALSE', '240', '94', '79', '1e-58', 'U10', 'Q10', '5100', '5600'],

            # 5) excluded_rule: is_complex=TRUE causes complex_requires_subunit_inference when P_EXCL is masked
            ['RXN-EXCL', 'complex rxn', '5.5.5.5', 'FALSE', 'FALSE', 'P_EXCL desc', 'good_blast', '|PWY-EXCL|', 'TRUE', '280', '97', '84', '1e-66', 'U11', 'Q11', '5700', '6200'],
            ['RXN-EXCL-CTX', 'ctx', '5.5.5.99', 'FALSE', 'FALSE', 'P_UNLAB5 desc', 'good_blast', '|PWY-EXCL|', 'FALSE', '275', '96', '82', '1e-65', 'U12', 'Q12', '6300', '6800'],

            # 6) maskable_no_context: P_NOCTX is the sole good_blast in PWY-NOCTX (so no other reaction has another protein)
            ['RXN-NOCTX', 'no ctx rxn', '6.6.6.6', 'FALSE', 'FALSE', 'P_NOCTX desc', 'good_blast', '|PWY-NOCTX|', 'FALSE', '245', '95', '80', '1e-59', 'U13', 'Q13', '6900', '7400'],
        ]
        with rx_tbl.open('w', encoding='utf-8', newline='') as fh:
            w = csv.writer(fh, delimiter='\t')
            w.writerow(['rxn', 'name', 'ec', 'keyrea', 'spont', 'stitle', 'status', 'pathway', 'is_complex', 'bitscore', 'qcovs', 'pident', 'evalue', 'dbhit', 'qseqid', 'sstart', 'send'])
            w.writerows(rx_rows)

        return rx_tbl, pw_tbl

    def test_synthetic_classification_and_sham_partner(self):
        with _make_temp_dir() as tmp:
            g_dir = Path(tmp) / 'GCF_SYNTH.1'
            rx_tbl, pw_tbl = self._write_synthetic_tables(g_dir, 'GCF_SYNTH.1')
            gt = registry.load_genome_tables(rx_tbl, pw_tbl)

            labelled_proteins = {
                'P_MASK', 'P_ALRGAP', 'P_REDUND', 'P_NOTPRED', 'P_EXCL', 'P_NOCTX', 'P_ABSENT',
            }

            # 1) maskable_with_context (maskable_gap)
            c_mask = registry.classify_locus_in_genome(gt, 'P_MASK', '1.1.1.1', labelled_proteins, 'lid_mask')
            self.assertEqual(c_mask['cls'], 'maskable_with_context')
            self.assertTrue(c_mask['maskable'])
            self.assertEqual(c_mask['gap_kind'], 'maskable_gap')
            self.assertTrue(c_mask['pathway_robust'])
            self.assertIn(c_mask['sham_protein'], {'P_UNLAB1', 'P_UNLAB2'})
            self.assertNotIn(c_mask['sham_protein'], labelled_proteins)
            self.assertNotEqual(c_mask['sham_reaction'], 'RXN-MASK')

            # 2) already_gap (maskable_with_context, gap_kind == 'already_gap')
            c_alr = registry.classify_locus_in_genome(gt, 'P_ALRGAP', '2.2.2.2', labelled_proteins, 'lid_alr')
            self.assertEqual(c_alr['cls'], 'maskable_with_context')
            self.assertTrue(c_alr['maskable'])
            self.assertEqual(c_alr['gap_kind'], 'already_gap')
            self.assertEqual(c_alr['sham_protein'], 'P_UNLAB3')
            self.assertNotIn(c_alr['sham_protein'], labelled_proteins)
            self.assertEqual(c_alr['sham_reaction'], 'RXN-OTHER2')

            # 3) redundant
            c_red = registry.classify_locus_in_genome(gt, 'P_REDUND', '3.3.3.3', labelled_proteins, 'lid_red')
            self.assertEqual(c_red['cls'], 'redundant')
            self.assertFalse(c_red['maskable'])

            # 4) pathway_not_predicted
            c_np = registry.classify_locus_in_genome(gt, 'P_NOTPRED', '4.4.4.4', labelled_proteins, 'lid_np')
            self.assertEqual(c_np['cls'], 'pathway_not_predicted')
            self.assertFalse(c_np['maskable'])

            # 5) excluded_rule
            c_ex = registry.classify_locus_in_genome(gt, 'P_EXCL', '5.5.5.5', labelled_proteins, 'lid_ex')
            self.assertEqual(c_ex['cls'], 'excluded_rule')
            self.assertFalse(c_ex['maskable'])

            # 6) maskable_no_context
            c_nc = registry.classify_locus_in_genome(gt, 'P_NOCTX', '6.6.6.6', labelled_proteins, 'lid_nc')
            self.assertEqual(c_nc['cls'], 'maskable_no_context')
            self.assertFalse(c_nc['maskable'])

            # 7) ec_not_in_tables
            c_ab = registry.classify_locus_in_genome(gt, 'P_ABSENT', '7.7.7.7', labelled_proteins, 'lid_ab')
            self.assertEqual(c_ab['cls'], 'ec_not_in_tables')
            self.assertFalse(c_ab['maskable'])

    def test_masking_touches_only_target_good_blast_and_validates(self):
        with _make_temp_dir() as tmp:
            g_dir = Path(tmp) / 'GCF_SYNTH.1'
            rx_tbl, pw_tbl = self._write_synthetic_tables(g_dir, 'GCF_SYNTH.1')
            unmasked_obs = observation('GCF_SYNTH.1', rx_tbl, pw_tbl)
            masked_obs = registry.build_masked_observation(
                'GCF_SYNTH.1',
                rx_tbl,
                pw_tbl,
                'P_MASK',
                'lid_mask',
                '1.1.1.1',
                scratch_dir=Path(tmp),
            )
            self.assertEqual(masked_obs['schema'], 'gspa.pathway-search-observations/1')
            self.assertEqual(masked_obs['masking']['rows_masked'], 1)
            self.assertEqual(masked_obs['masking']['good_blast_rows_naming_protein_after'], 0)
            self.assertEqual(len(masked_obs['reaction_rows']), len(unmasked_obs['reaction_rows']))

            for r_before, r_after in zip(unmasked_obs['reaction_rows'], masked_obs['reaction_rows']):
                tok = (r_before.get('stitle') or '').split(' ')[0].strip()
                if tok == 'P_MASK' and r_before['status'] == 'good_blast':
                    self.assertEqual(r_after['status'], 'no_blast')
                else:
                    self.assertEqual(r_after, r_before)

    def test_organism_set_removal_from_sbml(self):
        sbml_text = """<?xml version="1.0" encoding="UTF-8"?>
<sbml xmlns="http://www.sbml.org/sbml/level3/version1/core"
      xmlns:fbc="http://www.sbml.org/sbml/level3/version1/fbc/version2"
      level="3" version="1" fbc:required="false">
  <model id="test_model">
    <fbc:listOfGeneProducts>
      <fbc:geneProduct fbc:id="G_WP_0001_1" fbc:label="WP_0001.1"/>
      <fbc:geneProduct fbc:id="G_WP_0002_1" fbc:label="WP_0002.1"/>
    </fbc:listOfGeneProducts>
    <listOfReactions>
      <reaction id="rxn00001_c0" reversible="false">
        <fbc:geneProductAssociation>
          <fbc:geneProductRef fbc:geneProduct="G_WP_0001_1"/>
        </fbc:geneProductAssociation>
      </reaction>
      <reaction id="rxn00002_c0" reversible="false">
        <fbc:geneProductAssociation>
          <fbc:or>
            <fbc:geneProductRef fbc:geneProduct="G_WP_0001_1"/>
            <fbc:geneProductRef fbc:geneProduct="G_WP_0002_1"/>
          </fbc:or>
        </fbc:geneProductAssociation>
      </reaction>
      <reaction id="rxn00003_c0" reversible="false">
        <fbc:geneProductAssociation>
          <fbc:unknownTag fbc:geneProduct="G_WP_0001_1"/>
        </fbc:geneProductAssociation>
      </reaction>
    </listOfReactions>
  </model>
</sbml>
"""
        with _make_temp_dir() as tmp:
            sbml_path = Path(tmp) / 'draft.xml'
            sbml_path.write_text(sbml_text, encoding='utf-8')
            graph_rxns = {'rxn00001', 'rxn00002', 'rxn00003'}
            res = registry.masked_organism_reactions(sbml_path, graph_rxns, 'WP_0001.1')
            self.assertEqual(res['unmasked'], ['rxn00001', 'rxn00002', 'rxn00003'])
            self.assertEqual(res['removed'], ['rxn00001'])
            self.assertEqual(res['masked'], ['rxn00002', 'rxn00003'])
            self.assertEqual(res['unparsable_associations'], 1)
            self.assertEqual(res['unparsable_reactions'], ['rxn00003'])

    def test_claim_filter_top10_top64_ties_and_go_passthrough(self):
        with _make_temp_dir() as tmp:
            in_claims = Path(tmp) / 'in_claims.jsonl'
            out_claims = Path(tmp) / 'out_claims.jsonl'
            rows = []
            # Protein P1 has 12 ECs: top 9 have scores 0.99..0.91, 10th and 11th tie at 0.50, 12th is 0.10
            for i in range(1, 10):
                rows.append({'protein_id': 'P1', 'function_type': 'EC', 'function_id': f'1.1.1.{i}', 'raw_score': round(1.0 - 0.01 * i, 4)})
            rows.append({'protein_id': 'P1', 'function_type': 'EC', 'function_id': '1.1.1.10', 'raw_score': 0.50})
            rows.append({'protein_id': 'P1', 'function_type': 'EC', 'function_id': '1.1.1.11', 'raw_score': 0.50})
            rows.append({'protein_id': 'P1', 'function_type': 'EC', 'function_id': '1.1.1.12', 'raw_score': 0.10})
            # Also add a GO claim on P1 with low score 0.01 -> must pass through
            rows.append({'protein_id': 'P1', 'function_type': 'GO', 'function_id': 'GO:0003824', 'raw_score': 0.01})
            # And a 13th EC on P1 ('9.9.9.9') with low score 0.05 on P1, but P1 is the ONLY protein for 9.9.9.9 -> kept via top-64 per EC!
            rows.append({'protein_id': 'P1', 'function_type': 'EC', 'function_id': '9.9.9.9', 'raw_score': 0.05})
            # For 1.1.1.12 (where P1 had 0.10), add 64 other proteins with score 0.80 so P1 is NOT in top-64 for 1.1.1.12
            for j in range(1, 65):
                rows.append({'protein_id': f'Q{j}', 'function_type': 'EC', 'function_id': '1.1.1.12', 'raw_score': 0.80})

            with in_claims.open('w', encoding='utf-8') as fh:
                for r in rows:
                    fh.write(json.dumps(r) + '\n')

            rec = registry.filter_claims_file(
                in_claims,
                out_claims,
                top_ecs_per_protein=10,
                top_proteins_per_ec=64,
            )
            kept = [json.loads(line) for line in out_claims.read_text(encoding='utf-8').splitlines() if line.strip()]
            p1_kept = {r['function_id'] for r in kept if r['protein_id'] == 'P1'}

            # Both tied 10th/11th ECs (1.1.1.10 and 1.1.1.11) must be kept
            self.assertIn('1.1.1.10', p1_kept)
            self.assertIn('1.1.1.11', p1_kept)
            # 1.1.1.12 is 13th on P1 and 65th on 1.1.1.12 -> dropped!
            self.assertNotIn('1.1.1.12', p1_kept)
            # 9.9.9.9 is kept because P1 is in top-64 for 9.9.9.9
            self.assertIn('9.9.9.9', p1_kept)
            # GO claim passes through unchanged
            self.assertIn('GO:0003824', p1_kept)
            self.assertEqual(rec['non_ec_rows_in'], 1)
            self.assertEqual(rec['non_ec_rows_out'], 1)
            self.assertEqual(rec['rows_out'], len(rows) - 1)

    def test_claim_sets_assembly_from_fake_plan(self):
        with _make_temp_dir() as tmp:
            base = Path(tmp)
            part_clean = base / 'parts/0'
            part_mapred = base / 'parts/1'
            part_clean.mkdir(parents=True)
            part_mapred.mkdir(parents=True)
            (part_clean / 'claims.jsonl').write_text('{"protein_id":"P1","function_id":"1.1.1.1","raw_score":0.9}\n' * 5, encoding='utf-8')
            (part_mapred / 'claims.jsonl').write_text('{"protein_id":"P1","function_id":"1.1.1.1","raw_score":0.9}\n' * 2, encoding='utf-8')

            filtered_dir = base / 'filtered'
            (filtered_dir / 'clean' / 'GCF_P1.1').mkdir(parents=True)
            filt_clean = filtered_dir / 'clean' / 'GCF_P1.1' / 'claims.jsonl'
            filt_clean.write_text('{"protein_id":"P1","function_id":"1.1.1.1","raw_score":0.9}\n' * 3, encoding='utf-8')

            price_layout = base / 'price_layout'
            price_layout.mkdir()
            (price_layout / 'report.json').write_text(
                json.dumps({'assemblies': {'GCF_P1.1': {'taxid': 562}}}),
                encoding='utf-8',
            )
            (price_layout / 'GCF_P1.1.layout.tsv').write_text('protein_id\tcontig\nP1\tc1\n', encoding='utf-8')

            plan_path = base / 'plan.json'
            plan_path.write_text(
                json.dumps({
                    'layout': str(price_layout),
                    'cells': [
                        {'index': 0, 'predictor': 'CLEAN', 'assembly': 'GCF_P1.1', 'part': str(part_clean), 'claims': 5},
                        {'index': 1, 'predictor': 'MAPred', 'assembly': 'GCF_P1.1', 'part': str(part_mapred), 'claims': 2},
                    ],
                }),
                encoding='utf-8',
            )

            diamond_price = base / 'diamond_price/claims'
            diamond_price.mkdir(parents=True)
            (diamond_price / 'GCF_P1.1.claims.jsonl').write_text('{"protein_id":"P1","function_id":"1.1.1.1","raw_score":0.8}\n', encoding='utf-8')

            diamond_cohort = base / 'diamond_cohort/claims'
            diamond_cohort.mkdir(parents=True)
            (diamond_cohort / 'GCA_C1.1.claims.jsonl').write_text('{"protein_id":"C1","function_id":"2.2.2.2","raw_score":0.7}\n' * 4, encoding='utf-8')

            cohort_manifest = base / 'cohort_manifest.tsv'
            cohort_manifest.write_text(
                'genome_id\ttaxon\tessential_profile\tgenome_layout_path\n'
                f'GCA_C1.1\tNCBITaxon_2287\tarchaea\t{base}/cohort_layout/GCA_C1.1.layout.tsv\n',
                encoding='utf-8',
            )

            out_cs = base / 'claim_sets.json'
            res = registry.build_claim_sets(
                out_path=out_cs,
                rc1_plan=plan_path,
                diamond_price_dir=base / 'diamond_price',
                diamond_cohort_dir=base / 'diamond_cohort',
                price_layout=price_layout,
                cohort_batch_manifest=cohort_manifest,
                filtered_dir=filtered_dir,
            )
            by_key = res['by_key']
            self.assertTrue(by_key['clean:GCF_P1.1']['filtered'])
            self.assertEqual(by_key['clean:GCF_P1.1']['published_claims'], 3)
            self.assertEqual(by_key['clean:GCF_P1.1']['taxon'], 'NCBITaxon_562')
            self.assertFalse(by_key['mapred:GCF_P1.1']['filtered'])
            self.assertEqual(by_key['mapred:GCF_P1.1']['published_claims'], 2)
            self.assertEqual(by_key['diamond:GCF_P1.1']['published_claims'], 1)
            self.assertEqual(by_key['diamond:GCA_C1.1']['published_claims'], 4)
            self.assertEqual(by_key['diamond:GCA_C1.1']['essential_profile'], 'archaea')

    def test_real_tables_reproduce_119_maskable_loci(self):
        tables_dir = Path('/home/leechuck/tmp/row5-design/tables')
        per_locus_path = Path('/home/leechuck/tmp/row5-design/per_locus.json')
        if not tables_dir.is_dir() or not per_locus_path.is_file():
            self.skipTest('Design tables directory /home/leechuck/tmp/row5-design/tables is absent')
        with _make_temp_dir() as tmp:
            res = registry.build_loci(
                models_dir=tables_dir,
                out_dir=Path(tmp),
                check_against=per_locus_path,
            )
            self.assertTrue(res['check_against']['matches_exactly'])
            self.assertEqual(res['summary']['maskable_loci'], 119)
            self.assertEqual(res['summary']['maskable_by_label_set'], {'cohort_swissprot': 87, 'price148': 25, 'price_swissprot': 7})
            self.assertEqual(res['summary']['maskable_by_scope'], {'development': 58, 'heldout': 61})
            self.assertEqual(res['summary']['maskable_by_gap_kind'], {'maskable_gap': 107, 'already_gap': 12})
            self.assertEqual(res['summary']['sham_assigned'], 119)
            self.assertEqual(res['summary']['sham_missing'], 0)


class TestRow5MaskedContextRunCell(unittest.TestCase):
    def test_command_construction_all_arms(self):
        with _make_temp_dir() as tmp:
            base = Path(tmp)
            inputs = base / 'inputs'
            (inputs / 'observations').mkdir(parents=True)
            (inputs / 'organism').mkdir(parents=True)
            (inputs / 'cohort_metadata').mkdir(parents=True)

            lid = 'abc123def456'
            asm = 'GCF_000005845.2'
            for name in (f'{lid}.json', f'unmasked-{asm}.json', f'sham-{lid}.json'):
                (inputs / 'observations' / name).write_text('{}\n', encoding='utf-8')
            (inputs / 'organism' / f'{lid}.txt').write_text('rxn00001\n', encoding='utf-8')
            (inputs / 'cohort_metadata' / f'{lid}.json').write_text(
                json.dumps({'genomes': {asm: {'assemblyStatus': 'COMPLETE', 'contigs': {}, 'clade': 'c1', 'unplacedProteins': []}}}),
                encoding='utf-8',
            )

            locus = {
                'locus_id': lid,
                'genome': asm,
                'group': 'price',
                'protein_id': 'NP_414542.1',
                'ec': '2.7.2.4',
                'domain': 'bacteria',
            }
            claim_entry = {
                'claim_set': 'diamond',
                'genome': asm,
                'claims_path': str(base / 'claims.jsonl'),
                'published_claims': 100,
                'layout_tsv': str(base / f'{asm}.layout.tsv'),
                'taxon': 'NCBITaxon_511145',
                'essential_profile': 'bacteria',
            }
            parser = run_cell.build_arg_parser()
            args = parser.parse_args([
                '--loci', str(base / 'loci.json'),
                '--locus-id', lid,
                '--claim-set', 'diamond',
                '--claim-sets', str(base / 'claim_sets.json'),
                '--inputs', str(inputs),
                '--jar', str(base / 'gspa-3.0.1.jar'),
                '--jar-ctx', str(base / 'gspa-3.1.0-dev.jar'),
                '--reference', str(base / 'ref'),
                '--reaction-graph', str(base / 'seed_reactions.tsv'),
                '--reaction-aliases', str(base / 'aliases.tsv'),
                '--pathway-groups', str(base / 'pathway-groups.tsv'),
                '--orthogroups', str(base / 'orthogroups.tsv'),
                '--model', str(base / 'model.json'),
                '--window-bp', '5000', '10000', '20000',
                '--scratch', str(base / 'scratch'),
                '--out', str(base / 'out'),
            ])

            all_arms = run_cell.default_arms(args.window_bp)
            self.assertEqual(len(all_arms), 23)

            for arm in all_arms:
                arm_dir = base / 'out' / lid / 'diamond' / arm
                cmd, output, _bound = run_cell.command(args, locus, claim_entry, arm, arm_dir, inputs, base / 'scratch')
                self.assertEqual(output, arm_dir / f'{asm}.integrated.tsv')

                # Check JAR selection
                if arm.startswith('ctx_') or arm.startswith('sham_ctx_') or arm in ('gaps_ctxjar', 'gapsfix', 'sham_gapsfix'):
                    self.assertIn(str(base / 'gspa-3.1.0-dev.jar'), cmd)
                    self.assertNotIn(str(base / 'gspa-3.0.1.jar'), cmd)
                else:
                    self.assertIn(str(base / 'gspa-3.0.1.jar'), cmd)
                    self.assertNotIn(str(base / 'gspa-3.1.0-dev.jar'), cmd)

                if arm == 'stage0':
                    self.assertIn('--enable-priors', cmd)
                    self.assertEqual(cmd[cmd.index('--enable-priors') + 1], 'none')
                    self.assertNotIn('--pathway-observations', cmd)
                else:
                    self.assertIn('--pathway-observations', cmd)
                    obs_arg = cmd[cmd.index('--pathway-observations') + 1]
                    if arm == 'unmasked_gaps':
                        self.assertTrue(obs_arg.endswith(f'observations/unmasked-{asm}.json'))
                    elif arm in ('sham_gaps', 'sham_gapsfix') or arm.startswith('sham_ctx_universal_w'):
                        self.assertTrue(obs_arg.endswith(f'observations/sham-{lid}.json'))
                    else:
                        self.assertTrue(obs_arg.endswith(f'observations/{lid}.json'))

                # EC-namespace repair: on for gapsfix, sham_gapsfix and every context arm; off elsewhere
                if arm in ('gapsfix', 'sham_gapsfix') or arm.startswith('ctx_') or arm.startswith('sham_ctx_'):
                    self.assertIn('--gap-prior-ec-namespace-fix', cmd)
                else:
                    self.assertNotIn('--gap-prior-ec-namespace-fix', cmd)

                if arm in ('universal', 'pathway', 'organism'):
                    self.assertIn('--learned-model', cmd)
                    self.assertIn('--cohort', cmd)
                    self.assertEqual(cmd[cmd.index('--cohort-holdout') + 1], 'CLADE')
                    self.assertEqual(cmd[cmd.index('--cohort-context-mode') + 1], arm.upper())
                    self.assertTrue(cmd[cmd.index('--cohort-metadata') + 1].endswith(f'cohort_metadata/{lid}.json'))

                if arm.startswith('ctx_pathway_w'):
                    self.assertEqual(cmd[cmd.index('--gap-prior-context-mode') + 1], 'pathway')
                    self.assertEqual(cmd[cmd.index('--gap-prior-pathway-groups') + 1], str(base / 'pathway-groups.tsv'))

                if arm.startswith('ctx_organism_w'):
                    self.assertEqual(cmd[cmd.index('--gap-prior-context-mode') + 1], 'organism')
                    self.assertTrue(cmd[cmd.index('--gap-prior-organism-reactions') + 1].endswith(f'organism/{lid}.txt'))


class TestRow5MaskedContextEvaluate(unittest.TestCase):
    def test_sign_test_critical_values(self):
        self.assertIsNone(evaluate.sign_test_critical_k(4))
        self.assertEqual(evaluate.sign_test_critical_k(5), 5)
        self.assertEqual(evaluate.sign_test_critical_k(8), 7)
        self.assertEqual(evaluate.sign_test_critical_k(12), 10)
        self.assertEqual(evaluate.sign_test_critical_k(20), 15)
        self.assertEqual(evaluate.sign_test_critical_k(30), 20)
        self.assertEqual(evaluate.sign_test_critical_k(40), 26)

    def test_outcome_logic_including_mixed_and_ties(self):
        # 1) Unique top-1 vs tied top-1
        self.assertTrue(evaluate.compute_protein_top1({'1.1.1.1': 0.8, '1.1.1.2': 0.7}, '1.1.1.1'))
        self.assertFalse(evaluate.compute_protein_top1({'1.1.1.1': 0.8, '1.1.1.2': 0.8}, '1.1.1.1'))

        # 2) Gene for the gap: unique top protein genome-wide vs tie
        arm_unique = {'ec_top_proteins': {'1.1.1.1': [['P1', 0.6], ['P2', 0.4]]}}
        arm_tied = {'ec_top_proteins': {'1.1.1.1': [['P1', 0.6], ['P2', 0.6]]}}
        self.assertTrue(evaluate.compute_gene_for_gap(arm_unique, 'P1', '1.1.1.1'))
        self.assertFalse(evaluate.compute_gene_for_gap(arm_tied, 'P1', '1.1.1.1'))

        # 3) Mixed outcome (call gained, gene_for_gap lost) -> raw_outcome='mixed', outcome='harmed'
        ind_a = {'call': False, 'protein_top1': False, 'gene_for_gap': True}
        ind_b = {'call': True, 'protein_top1': False, 'gene_for_gap': False}
        outcome, raw_outcome, deltas = evaluate.classify_locus_outcome(ind_a, ind_b)
        self.assertEqual(raw_outcome, 'mixed')
        self.assertEqual(outcome, 'harmed')
        self.assertEqual(deltas['call'], 'gained')
        self.assertEqual(deltas['gene_for_gap'], 'lost')

    def test_harm_elsewhere_and_unconfirmed_accounting(self):
        label_universe = {
            'P_TARGET': ['1.1.1.1'],
            'P_OTHER1': ['2.2.2.2'],
            'P_OTHER2': ['3.3.3.3'],
        }
        arm_a = {
            'protein_scores': {
                'P_TARGET': {'1.1.1.1': 0.30},
                'P_OTHER1': {'2.2.2.2': 0.60, '9.9.9.9': 0.20},
                'P_OTHER2': {'3.3.3.3': 0.40, '3.3.3.4': 0.20},
            }
        }
        arm_b = {
            'protein_scores': {
                'P_TARGET': {'1.1.1.1': 0.70},
                # P_OTHER1 loses call on 2.2.2.2 (0.60 -> 0.45) and gains unconfirmed call on 9.9.9.9 (0.20 -> 0.55)
                'P_OTHER1': {'2.2.2.2': 0.45, '9.9.9.9': 0.55},
                # P_OTHER2 keeps call=False on 3.3.3.3, but loses protein_top1 because 3.3.3.4 rises to 0.45
                'P_OTHER2': {'3.3.3.3': 0.40, '3.3.3.4': 0.45},
            }
        }
        he = evaluate.compute_harm_elsewhere_and_unconfirmed(
            arm_a, arm_b, label_universe, 'P_TARGET', '1.1.1.1'
        )
        # Both P_OTHER1 (call lost) and P_OTHER2 (top1 lost) count as harm elsewhere; unconfirmed 9.9.9.9 is reported separately
        self.assertEqual(he['harm_elsewhere_count'], 2)
        self.assertEqual(he['unconfirmed_calls_created_count'], 1)
        self.assertEqual(he['unconfirmed_calls_created'][0]['ec'], '9.9.9.9')

    def test_four_verdicts(self):
        def make_row(i: int, outcome: str, label_set: str = 'price148', claim_set: str = 'diamond',
                     pathway_robust: str = 'robust', gap_kind: str = 'maskable_gap',
                     domain: str = 'bacteria', harm_elsewhere: int = 0) -> dict:
            return {
                'locus_id': f'loc_{i}',
                'genome': f'G_{i}',
                'outcome': outcome,
                'raw_outcome': outcome,
                'label_set': label_set,
                'claim_set': claim_set,
                'pathway_robust': pathway_robust,
                'gap_kind': gap_kind,
                'domain': domain,
                'harm_elsewhere': harm_elsewhere,
                'unconfirmed_calls_created': 0,
                'delta_calls_ge_05': 0,
                'changed_scores': 1,
                'rank_change': 1.0 if outcome == 'improved' else (-1.0 if outcome == 'harmed' else 0.0),
            }

        # 1) Systematic improvement: 5 improved, 0 harmed, consistent strata, harm_elsewhere=0
        rows_sys = [make_row(i, 'improved') for i in range(5)]
        res_sys = evaluate.apply_decision_rule(rows_sys)
        self.assertEqual(res_sys['verdict'], evaluate.VERDICT_SYSTEMATIC)

        # 2) No detectable effect: 4 improved, 0 harmed (< 5 discordant)
        rows_none = [make_row(i, 'improved') for i in range(4)]
        res_none = evaluate.apply_decision_rule(rows_none)
        self.assertEqual(res_none['verdict'], evaluate.VERDICT_NONE)

        # 3) Inconsistent: 7 improved, 1 harmed (p = 9/256 <= 0.05), but harm_elsewhere = 8 > improved (7)
        rows_inc = [make_row(i, 'improved', harm_elsewhere=(8 if i == 0 else 0)) for i in range(7)] + [make_row(7, 'harmed')]
        res_inc = evaluate.apply_decision_rule(rows_inc)
        self.assertEqual(res_inc['verdict'], evaluate.VERDICT_INCONSISTENT)

        # Also test stratum inconsistency: 10 improved on robust, 2 harmed on marginal (harmed > improved + 1 in marginal)
        rows_inc_stratum = (
            [make_row(i, 'improved', pathway_robust='robust') for i in range(10)]
            + [make_row(10 + j, 'harmed', pathway_robust='marginal') for j in range(2)]
        )
        res_inc_stratum = evaluate.apply_decision_rule(rows_inc_stratum)
        self.assertEqual(res_inc_stratum['verdict'], evaluate.VERDICT_INCONSISTENT)

        # 4) Harm: 0 improved, 5 harmed
        rows_harm = [make_row(i, 'harmed') for i in range(5)]
        res_harm = evaluate.apply_decision_rule(rows_harm)
        self.assertEqual(res_harm['verdict'], evaluate.VERDICT_HARM)

    def test_window_choice_with_ties_and_heldout_refusal(self):
        # Tie between 5000 and 10000 -> must pick 10000
        self.assertEqual(evaluate.choose_best_window({5000: 4, 10000: 4, 20000: 2}), 10000)
        # Strict winner 20000 -> picks 20000
        self.assertEqual(evaluate.choose_best_window({5000: 2, 10000: 3, 20000: 5}), 20000)
        # Tie between 5000 and 20000 (both 5, while 10000 has 3) -> equidistant to 10000, smallest W (5000)
        self.assertEqual(evaluate.choose_best_window({5000: 5, 10000: 3, 20000: 5}), 5000)

        with _make_temp_dir() as tmp:
            base = Path(tmp)
            loci_path = base / 'loci.json'
            cells_dir = base / 'cells'
            loci_data = {
                'loci': {
                    'dev_loc1': {
                        'locus_id': 'dev_loc1',
                        'genome': 'G_DEV',
                        'group': 'cohort',
                        'label_set': 'cohort_swissprot',
                        'split': 'train',
                        'scope': 'development',
                        'domain': 'archaea',
                        'protein_id': 'P_DEV',
                        'ec': '1.1.1.1',
                        'maskable': True,
                        'gap_kind': 'maskable_gap',
                        'pathway_robust': True,
                    },
                    'held_loc1': {
                        'locus_id': 'held_loc1',
                        'genome': 'G_HELD',
                        'group': 'price',
                        'label_set': 'price148',
                        'split': 'test',
                        'scope': 'heldout',
                        'domain': 'bacteria',
                        'protein_id': 'P_HELD',
                        'ec': '2.2.2.2',
                        'maskable': True,
                        'gap_kind': 'maskable_gap',
                        'pathway_robust': True,
                    },
                },
                'all_loci': [],
                'label_universe': {
                    'G_DEV': {'P_DEV': ['1.1.1.1']},
                    'G_HELD': {'P_HELD': ['2.2.2.2']},
                },
            }
            loci_data['all_loci'] = list(loci_data['loci'].values())
            loci_path.write_text(json.dumps(loci_data), encoding='utf-8')

            # Write development cell where w5000 and w10000 both improve (tie -> 10000)
            dev_cell_dir = cells_dir / 'dev_loc1' / 'diamond'
            dev_cell_dir.mkdir(parents=True)
            (dev_cell_dir / 'retained.json').write_text(
                json.dumps({
                    'locus_id': 'dev_loc1',
                    'claim_set': 'diamond',
                    'genome': 'G_DEV',
                    'scope': 'development',
                    'scores': {
                        'stage0': {'P_DEV': {'1.1.1.1': 0.30}},
                        'gaps': {'P_DEV': {'1.1.1.1': 0.40}},
                        'ctx_universal_w5000': {'P_DEV': {'1.1.1.1': 0.60}},
                        'ctx_universal_w10000': {'P_DEV': {'1.1.1.1': 0.65}},
                        'ctx_universal_w20000': {'P_DEV': {'1.1.1.1': 0.42}},
                    },
                }),
                encoding='utf-8',
            )

            # Heldout evaluation MUST refuse when parameter_choice.json is missing
            with self.assertRaises(ValueError):
                evaluate.evaluate(
                    cells_dir=cells_dir,
                    loci_path=loci_path,
                    scope='heldout',
                    out_dir=base / 'eval_held_fail',
                    parameter_choice_path=None,
                )

            # Choose window on development -> writes parameter_choice.json with chosen_window_bp == 10000
            pc_path = base / 'eval_dev' / 'parameter_choice.json'
            dev_res = evaluate.evaluate(
                cells_dir=cells_dir,
                loci_path=loci_path,
                scope='development',
                out_dir=base / 'eval_dev',
                parameter_choice_path=pc_path,
                choose_window=True,
            )
            self.assertEqual(dev_res['chosen_window_bp'], 10000)
            self.assertTrue(pc_path.is_file())

            # Also verify inconsistent scope label in loci.json is refused
            bad_loci_path = base / 'bad_loci.json'
            bad_data = json.loads(loci_path.read_text(encoding='utf-8'))
            bad_data['loci']['held_loc1']['scope'] = 'development'  # price/test cannot be development!
            bad_data['all_loci'] = list(bad_data['loci'].values())
            bad_loci_path.write_text(json.dumps(bad_data), encoding='utf-8')
            with self.assertRaises(ValueError):
                evaluate.evaluate(
                    cells_dir=cells_dir,
                    loci_path=bad_loci_path,
                    scope='heldout',
                    out_dir=base / 'eval_bad_scope',
                    parameter_choice_path=pc_path,
                )

    def test_submit_script_pilot_and_packing(self):
        submit_sh = V3 / 'shaheen' / 'row5-masked-context-submit.sh'
        with _make_temp_dir() as tmp:
            base = Path(tmp)
            loci_tsv = base / 'loci.tsv'
            with loci_tsv.open('w', encoding='utf-8', newline='') as fh:
                w = csv.writer(fh, delimiter='\t')
                w.writerow(['locus_id', 'genome', 'group', 'maskable'])
                for i in range(4):
                    w.writerow([f'0{i}0000000000', f'GCA_C{i}', 'cohort', 'true'])
                for i in range(5):
                    w.writerow([f'1{i}0000000000', f'GCF_P{i}', 'price', 'true'])

            tokens_file = base / 'pilot_tokens.txt'
            proc = subprocess.run(
                [
                    'bash', str(submit_sh),
                    '--loci', str(loci_tsv),
                    '--pilot',
                    '--tokens-out', str(tokens_file),
                ],
                capture_output=True,
                text=True,
                check=True,
            )
            tokens = [l.strip() for l in tokens_file.read_text(encoding='utf-8').splitlines() if l.strip()]
            # 2 cohort loci x 1 claim set + 3 Price loci x 7 claim sets = 23 tokens
            self.assertEqual(len(tokens), 23)
            self.assertIn('--exclude=nid00022,nid00497,nid01033', proc.stdout)
            self.assertIn('REPO_ROOT=', proc.stdout)


if __name__ == '__main__':
    unittest.main()
