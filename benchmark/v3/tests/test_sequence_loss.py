import csv
import json
import sys
import tempfile
import unittest
from unittest import mock
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import prepare_sequence_loss as P
import run_sequence_loss_benchmark as E


def fixture(root):
    for d in ('cohort_inputs', 'layout', 'baseline/claims/homology'):
        (root / d).mkdir(parents=True)
    genomes = {}
    og = ['genome_id\tprotein_id\torthogroup_id\tsubfamily_id\tweight']
    for gid in ('g1', 'g2'):
        genome = root / 'genomes' / gid
        genome.mkdir(parents=True)
        P.write_fasta(genome / 'genomic.fna', {'chr': 'ACGT' * 25, 'plasmid': 'TGCA' * 10})
        # Ten complete CDSs on a chromosome, one boundary CDS, a plasmid CDS and an unplaced protein.
        proteins = {f'p{i}': 'MAAA' for i in range(10)} | {'boundary': 'MKK', 'plasmid_p': 'MWW', 'unplaced': 'MNN'}
        P.write_fasta(genome / 'protein.faa', proteins)
        rows = ['protein_id\tcontig\tstart\tend\tstrand']
        rows += [f'p{i}\tchr\t{i*10+1}\t{i*10+9}\t-' for i in range(10)]
        rows += ['boundary\tchr\t9\t12\t+', 'plasmid_p\tplasmid\t1\t9\t+']
        (root / 'layout' / f'{gid}.layout.tsv').write_text('\n'.join(rows) + '\n')
        claims = [{'protein_id': p, 'function_type': 'EC', 'function_id': '1.1.1.1', 'raw_score': 1e-8}
                  for p in sorted(proteins)]
        (root / 'baseline/claims/homology' / f'{gid}.claims.jsonl').write_text(
            ''.join(json.dumps(c) + '\n' for c in claims))
        genomes[gid] = {'assemblyStatus': 'COMPLETE', 'contigs': {'chr': {'length': 100, 'closed': True},
                       'plasmid': {'length': 40, 'closed': True}}, 'unplacedProteins': ['unplaced'],
                       'clade': gid, 'familySearchComplete': True}
        og += [f'{gid}\t{p}\tfamily_A\t-\t1' for p in proteins]
    P.write_json(root / 'cohort_inputs/cohort_metadata.json', {'genomes': genomes})
    (root / 'cohort_inputs/cohort_orthogroups.tsv').write_text('\n'.join(og) + '\n')


class SequenceLossTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.cache = self.root / 'cache'
        fixture(self.cache)

    def test_actual_sequence_claims_and_memberships_removed_together(self):
        out = self.root / 'out'
        r = P.prepare(self.cache, out, block_bp=10)
        with (out / 'orthogroups.tsv').open() as fh:
            og = list(csv.DictReader(fh, delimiter='\t'))
        meta = json.loads((out / 'metadata_sequence_loss.json').read_text())
        coords = json.loads((out / 'metadata_coordinate_loss.json').read_text())
        for gid, g in r['genomes'].items():
            kept, lost = set(g['retained_proteins']), set(g['removed_proteins'])
            self.assertGreater(len(lost), 1)
            self.assertFalse(kept & lost)
            self.assertEqual(len(kept | lost), 13)
            self.assertIn('unplaced', kept)
            self.assertEqual(g['removed_proteins']['boundary'], 'boundary_crossing_CDS')
            self.assertEqual(set(P.fasta_sequences(out / 'proteins' / f'{gid}.faa')), kept)
            self.assertEqual({x['protein_id'] for x in og if x['genome_id'] == gid}, kept)
            claims = [json.loads(l) for l in (out / 'claims' / f'{gid}.claims.jsonl').read_text().splitlines()]
            self.assertEqual({c['protein_id'] for c in claims}, kept)
            self.assertTrue(all(c['raw_score'] == 1e-8 for c in claims))
            dna = P.fasta_sequences(out / 'assemblies' / f'{gid}.fna')
            self.assertEqual(sum(map(len, dna.values())), 70)
            self.assertEqual(meta['genomes'][gid]['completenessFraction'], 0.5)
            self.assertEqual(meta['genomes'][gid]['unplacedProteins'], ['unplaced'])
            self.assertEqual(set(coords['genomes'][gid]['unplacedProteins']), lost | {'unplaced'})
            self.assertEqual(coords['genomes'][gid]['completenessFraction'], 1.0)
            with (out / 'layout' / f'{gid}.layout.tsv').open() as fh:
                for row in csv.DictReader(fh, delimiter='\t'):
                    self.assertIn(row['protein_id'], kept)
                    self.assertLessEqual(int(row['end']), len(dna[row['contig']]))
                    self.assertGreater(int(row['start']), 0)
        P.verify(out)

    def test_selection_is_reproducible_and_does_not_read_labels(self):
        a = P.prepare(self.cache, self.root / 'a', block_bp=10)
        (self.cache / 'labels').mkdir()
        (self.cache / 'labels/labels.tsv').write_text('arbitrary truth that must not be read')
        b = P.prepare(self.cache, self.root / 'b', block_bp=10)
        self.assertEqual(a, b)
        c = P.prepare(self.cache, self.root / 'c', block_bp=10, seed=1)
        self.assertNotEqual(a['genomes'], c['genomes'])
        self.assertFalse(any('/labels/' in p for p in a['inputs']))

    def test_immutable_outputs_and_changed_inputs_are_rejected(self):
        out = self.root / 'out'
        P.prepare(self.cache, out, block_bp=10)
        with self.assertRaises(FileExistsError):
            P.prepare(self.cache, out, block_bp=10)
        path = out / 'claims/g1.claims.jsonl'
        path.write_text(path.read_text() + '{}\n')
        with self.assertRaisesRegex(ValueError, 'changed sequence-loss output'):
            P.verify(out)
        src = self.cache / 'genomes/g1/protein.faa'
        src.write_text(src.read_text() + '>extra\nMMM\n')
        with self.assertRaisesRegex(ValueError, 'changed sequence-loss source'):
            P.verify(out)

    def test_undeclared_claim_cannot_leak_into_control(self):
        path = self.cache / 'baseline/claims/homology/g1.claims.jsonl'
        path.write_text(path.read_text() + json.dumps({'protein_id': 'not_in_fasta'}) + '\n')
        with self.assertRaisesRegex(ValueError, 'claim for undeclared protein'):
            P.prepare(self.cache, self.root / 'out', block_bp=10)

    def test_layout_fasta_and_metadata_must_cover_same_proteome(self):
        src = self.cache / 'genomes/g1/protein.faa'
        src.write_text(src.read_text() + '>extra\nMMM\n')
        with self.assertRaisesRegex(ValueError, 'FASTA must match'):
            P.prepare(self.cache, self.root / 'out', block_bp=10)

    def test_invalid_configuration_rejected(self):
        for kw in ({'block_bp': 0}, {'block_bp': True}, {'drop_fraction': 1.}, {'drop_fraction': float('nan')}):
            with self.assertRaises(ValueError):
                P.prepare(self.cache, self.root / 'out', **kw)


class EvaluationTests(unittest.TestCase):
    def test_actual_route_uses_filtered_claims_and_groups_only_for_sequence_loss(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            cache, inputs, res = root / 'cache', root / 'inputs', root / 'res'
            fixture(cache)
            P.prepare(cache, inputs, block_bp=10)
            res.mkdir()
            (res / 'splits.tsv').write_text('genome_id\tsplit\ng1\tcalibration\ng2\ttest\n')
            with mock.patch.object(E.R, 'CACHE', cache), mock.patch.object(E.R, 'RES', res), \
                 mock.patch.object(E.R, 'PERT_DIR', root / 'work'), \
                 mock.patch.object(E.R, 'run_perturbation_batch') as batch:
                E.run(inputs, root / 'work', root / 'model.json')
            self.assertEqual([c.args[0] for c in batch.call_args_list], list(E.ARMS))
            base, coord, loss = [c.kwargs for c in batch.call_args_list]
            self.assertTrue(all(c['rlc_suggester'] is False for c in (base, coord, loss)))
            self.assertIsNone(base['claims_override'])
            self.assertIsNone(coord['claims_override'])
            self.assertEqual(base['cohort_og_path'], coord['cohort_og_path'])
            self.assertEqual(loss['cohort_og_path'], inputs / 'orthogroups.tsv')
            self.assertEqual(set(loss['claims_override']), {'g1', 'g2'})
            self.assertTrue(all(p.parent == inputs / 'claims' for p in loss['claims_override'].values()))
            self.assertEqual(coord['layout_dir'], loss['layout_dir'])

    def test_identical_sequences_and_protein_ids_do_not_cross_genomes_or_loci(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            for g, p, s in [('g1', 'copy1', '.9'), ('g2', 'copy1', '.2')]:
                (root / f'{g}.integrated.tsv').write_text(
                    'protein_id\ttype\tfunction_id\tposterior_prob\n' + f'{p}\tEC\t1.1.1.1\t{s}\n')
            pred = E.read_predictions(root, ['g1', 'g2'])
            self.assertEqual(pred['g1']['EC'][('copy1', '1.1.1.1')], .9)
            self.assertEqual(pred['g2']['EC'][('copy1', '1.1.1.1')], .2)
            self.assertNotIn(('copy2', '1.1.1.1'), pred['g1']['EC'])

    def test_calibration_never_uses_test_values_and_unknowns_count_against_support(self):
        pred = {'cal': {'GO': {}, 'EC': {(str(i), 'EC:1.1.1.1'): .8 for i in range(10)}},
                'test': {'GO': {}, 'EC': {('test_pid', 'EC:2.2.2.2'): 1.0}}}
        truth = {'cal': {'GO': set(), 'EC': {(str(i), 'EC:1.1.1.1') for i in range(10)}},
                 'test': {'GO': set(), 'EC': set()}}
        policy = E.select_policy(pred, truth, ['cal'])
        self.assertEqual(policy['kinds']['EC']['threshold'], .8)
        pred['test']['EC'] = {}
        truth['test']['EC'] = {('new', 'EC:3.3.3.3')}
        self.assertEqual(E.select_policy(pred, truth, ['cal']), policy)
        for i in range(10):
            pred['cal']['EC'][(str(i), 'EC:9.9.9.9')] = .9
        self.assertIsNone(E.select_policy(pred, truth, ['cal'])['kinds']['EC']['threshold'])

    def test_removed_genes_stay_in_original_denominator_and_retained_scope_is_separate(self):
        truth = {'g1': {'EC': {('copy1', 'EC:1.1.1.1'), ('copy2', 'EC:1.1.1.1')}, 'GO': set()}}
        intact = {'g1': {'EC': {k: .9 for k in truth['g1']['EC']}, 'GO': {}}}
        loss = {'g1': {'EC': {('copy1', 'EC:1.1.1.1'): .9}, 'GO': {}}}
        pred = {'intact': intact, 'coordinate_loss': intact, 'sequence_loss': loss}
        prep = {'genomes': {'g1': {'retained_proteins': ['copy1']}}}
        pol = {'kinds': {'EC': {'threshold': .5}, 'GO': {'threshold': None}}}
        report = E.summarize(pred, truth, prep, pol, ['g1'], n_boot=10)
        all_ = report['scopes']['all_original']['EC']['sequence_loss']
        kept = report['scopes']['retained']['EC']['sequence_loss']
        lost = report['scopes']['removed']['EC']['sequence_loss']
        self.assertEqual(all_['ec_top1_reference_agreement']['estimate'], .5)
        self.assertEqual(all_['lost_supported_calls'], 1)
        self.assertEqual(all_['ec_top1_delta_vs_intact']['estimate'], -.5)
        self.assertEqual(kept['ec_top1_reference_agreement']['estimate'], 1)
        self.assertEqual(kept['ec_top1_delta_vs_intact']['estimate'], 0)
        self.assertEqual(lost['ec_top1_counts']['called'], 0)
        self.assertEqual(lost['reference_pairs'], 1)
        self.assertEqual(lost['supported_calls'], 0)

    def test_unsupported_predictions_are_not_reported_as_experimental_harms(self):
        truth = {'g1': {'EC': {('p', 'EC:1.1.1.1')}, 'GO': set()}}
        pred = {a: {'g1': {'EC': {('p', 'EC:9.9.9.9'): .9}, 'GO': {}}} for a in E.ARMS}
        prep = {'genomes': {'g1': {'retained_proteins': ['p']}}}
        pol = {'kinds': {'EC': {'threshold': .5}, 'GO': {'threshold': None}}}
        r = E.summarize(pred, truth, prep, pol, ['g1'], n_boot=0)
        m = r['scopes']['all_original']['EC']['sequence_loss']
        self.assertEqual(m['unverified_calls'], 1)
        self.assertEqual(m['supported_calls'], 0)
        self.assertNotIn('harmful', m)


if __name__ == '__main__':
    unittest.main()
