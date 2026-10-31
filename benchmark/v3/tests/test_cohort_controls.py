import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import run_cohort_perturbations as R
import prepare_cohort_controls as P
import run_cohort_controls as C
import run_sequence_loss_benchmark as E
from test_sequence_loss import fixture


class CandidateComparison(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.a, self.b = self.root / 'a', self.root / 'b'
        self.a.mkdir(); self.b.mkdir()

    def write(self, folder, rows):
        (folder / 'g.integrated.tsv.learned.json').write_text(json.dumps({'candidates': rows}))

    def test_null_scores_are_unassessed_and_support_transitions_are_visible(self):
        self.write(self.a, [{'candidate_key': 'both', 'probability': None},
                            {'candidate_key': 'changed_support', 'probability': None}])
        self.write(self.b, [{'candidate_key': 'both', 'probability': None},
                            {'candidate_key': 'changed_support', 'probability': .8}])
        result = R.compare_test_cohort_features(self.a, self.b, ['g'])
        self.assertEqual(result['unsupported_both'], 1)
        self.assertEqual(result['support_transitions'], 1)
        self.assertEqual(result['n_test_candidates_compared'], 0)
        self.assertIsNone(result['max_abs_probability_diff'])

    def test_missing_candidates_do_not_disappear_from_the_comparison(self):
        self.write(self.a, [{'candidate_key': 'shared', 'probability': .4},
                            {'candidate_key': 'removed', 'probability': .8}])
        self.write(self.b, [{'candidate_key': 'shared', 'probability': .5},
                            {'candidate_key': 'added', 'probability': .2}])
        result = R.compare_test_cohort_features(self.a, self.b, ['g'])
        self.assertEqual(result['missing_after'], 1)
        self.assertEqual(result['missing_before'], 1)
        self.assertEqual(result['n_test_candidates_compared'], 1)
        self.assertAlmostEqual(result['max_abs_probability_diff'], .1)

    def test_duplicate_candidate_ids_and_nonfinite_scores_fail(self):
        self.write(self.b, [])
        for rows in ([{'candidate_key': 'p', 'probability': .4}] * 2,
                     [{'candidate_key': 'p', 'probability': float('nan')}],
                     [{'candidate_key': 'p', 'probability': True}]):
            self.write(self.a, rows)
            with self.assertRaises(ValueError):
                R.compare_test_cohort_features(self.a, self.b, ['g'])


class Preparation(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.cache = self.root / 'cache'
        fixture(self.cache)
        meta = self.cache / 'cohort_inputs/cohort_metadata.json'
        d = json.loads(meta.read_text())
        d['genomes']['g1']['strainCluster'] = 'source_strain'
        meta.write_text(json.dumps(d))
        self.splits = self.root / 'splits.tsv'
        self.splits.write_text('genome_id\tsplit\ng1\ttrain\ng2\tcalibration\n')

    def prepare(self, out='prepared', **kwargs):
        return P.prepare(self.cache, self.root / out, self.splits, 'g1', **kwargs)

    def test_shuffle_preserves_contig_slots_and_identity_universe_and_duplicates_stay_in_cluster(self):
        import csv
        report = self.prepare()
        out = self.root / 'prepared'
        for gid in ('g1', 'g2'):
            def rows(p):
                with p.open() as fh:
                    return list(csv.DictReader(fh, delimiter='\t'))
            original = rows(self.cache / 'layout' / f'{gid}.layout.tsv')
            shuffled = rows(out / 'layout_shuffled' / f'{gid}.layout.tsv')
            for contig in ('chr', 'plasmid'):
                a, b = [r for r in original if r['contig'] == contig], [r for r in shuffled if r['contig'] == contig]
                self.assertEqual({r['protein_id'] for r in a}, {r['protein_id'] for r in b})
                self.assertEqual(sorted((r['start'], r['end'], r['strand']) for r in a),
                                 sorted((r['start'], r['end'], r['strand']) for r in b))
            self.assertGreater(report['genomes'][gid]['moved_by_shuffle'], 0)
            self.assertEqual(report['genomes'][gid]['unplaced'], 1)
        meta = json.loads((out / 'metadata_duplicated.json').read_text())['genomes']
        og = rows(out / 'orthogroups_duplicated.tsv')
        source = [{k:v for k,v in row.items() if k!='genome_id'} for row in og if row['genome_id']=='g1']
        for gid in report['duplicate_genomes']:
            self.assertEqual(meta[gid], meta['g1'])
            self.assertEqual([{k:v for k,v in row.items() if k!='genome_id'} for row in og if row['genome_id']==gid], source)
        P.verify(out)

    def test_preparation_is_immutable_deterministic_and_does_not_read_function_labels(self):
        a = self.prepare('a')
        (self.cache / 'labels.tsv').write_text('made-up labels that must not be read')
        b = self.prepare('b')
        self.assertEqual(a, b)
        self.assertFalse(any('labels.tsv' in p for p in a['inputs']))
        with self.assertRaises(FileExistsError):
            self.prepare('a')
        p = self.root / 'a/layout_shuffled/g1.layout.tsv'
        p.write_text(p.read_text() + 'corrupt\n')
        with self.assertRaisesRegex(ValueError, 'changed cohort control output'):
            P.verify(self.root / 'a')

    def test_unknown_genes_and_nontraining_duplicate_sources_fail(self):
        with self.assertRaisesRegex(ValueError, 'training genome'):
            P.prepare(self.cache, self.root/'bad', self.splits, 'g2')
        p = self.cache / 'baseline/claims/homology/g1.claims.jsonl'
        p.write_text(p.read_text()+json.dumps({'protein_id':'undeclared'})+'\n')
        with self.assertRaisesRegex(ValueError, 'undeclared protein'):
            self.prepare()

    def test_runner_has_matched_intact_arm_and_only_duplicate_arm_adds_peers(self):
        report = self.prepare()
        with mock.patch.object(R, 'CACHE', self.cache), mock.patch.object(R, 'PERT_DIR', self.root/'runs'), \
             mock.patch.object(R, 'run_perturbation_batch') as batch:
            C.run(self.root/'prepared', self.root/'runs', self.root/'model.json')
        self.assertEqual([c.args[0] for c in batch.call_args_list], list(C.ARMS))
        a, b, c = [v.kwargs for v in batch.call_args_list]
        self.assertEqual(a['genomes'], b['genomes'])
        self.assertEqual(c['genomes'], a['genomes']+report['duplicate_genomes'])
        self.assertEqual(a['layout_dir'], c['layout_dir'])
        self.assertEqual(a['cohort_og_path'], b['cohort_og_path'])
        self.assertTrue(all(v['rlc_suggester'] is False for v in (a,b,c)))

    def test_shared_evaluator_accepts_control_arms_without_sequence_loss_or_duplicate_truth(self):
        truth = {'g': {'GO': set(), 'EC': {('p','1.1.1.1')}}}
        predictions = {a:{'g': {'GO':{},'EC':{('p','1.1.1.1'):.8}}} for a in C.ARMS}
        result = E.summarize(predictions,truth,{'genomes':{'g':{'retained_proteins':['p']}}},
            {'kinds':{'GO':{'threshold':None},'EC':{'threshold':.8}}},['g'],n_boot=10,
            arms=C.ARMS,scopes=('all_original',))
        self.assertEqual(set(result['scopes']),{'all_original'})
        for arm in C.ARMS:
            x=result['scopes']['all_original']['EC'][arm]
            self.assertEqual(x['ec_top1_counts']['labelled'],1)
            self.assertEqual(x['ec_top1_delta_vs_intact']['estimate'],0)


class EvaluationIntegrity(unittest.TestCase):
    def test_extra_duplicate_outputs_are_verified_but_never_added_to_truth_scope(self):
        with tempfile.TemporaryDirectory() as d:
            work = Path(d)
            out = work / 'duplicates/out'
            out.mkdir(parents=True)
            spec = {'inputs': {}, 'outputs': {'g': 'original', 'g_dup1': 'duplicate'}}
            receipt = out / E.V.RECEIPT
            receipt.write_text(json.dumps({'spec': spec}))
            with mock.patch.object(E.V, 'reusable', return_value=True), \
                 mock.patch.object(E, 'read_predictions', return_value={}) as reader:
                E.verified_predictions(work, ['g'], arms=['duplicates'], extra_genomes={'duplicates':['g_dup1']})
                self.assertEqual(reader.call_args.args[1], ['g'])
                with self.assertRaisesRegex(ValueError, 'coverage differs'):
                    E.verified_predictions(work, ['g'], arms=['duplicates'])
                spec['outputs']['undeclared'] = 'unexpected'
                receipt.write_text(json.dumps({'spec':spec}))
                with self.assertRaisesRegex(ValueError, 'coverage differs'):
                    E.verified_predictions(work, ['g'], arms=['duplicates'], extra_genomes={'duplicates':['g_dup1']})

    def test_policy_freezes_calibration_only_and_rejects_changed_sources(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            inputs = root / 'inputs'; inputs.mkdir()
            (inputs/'preparation.json').write_text('{}')
            labels, splits_file, go_file = (root/name for name in ('labels.tsv','splits.tsv','go.obo'))
            labels.write_text('genome_id\tprotein_id\tfunction_id\ttier\n'+
                ''.join(f'cal\tp{i}\tEC:1.1.1.1\texperimental\n' for i in range(10))+
                'test\tp\tEC:2.2.2.2\texperimental\n')
            splits_file.write_text('bound split source')
            go_file.write_text('bound ontology source')
            splits = {'cal':'calibration','test':'test'}
            predictions = {arm:{'cal':{'GO':{},'EC':{(f'p{i}','1.1.1.1'):.8 for i in range(10)}},
                'test':{'GO':{},'EC':{('p','2.2.2.2'):.95}}} for arm in C.ARMS}
            preparation = {'genomes':{g:{'retained_proteins':[f'p{i}' for i in range(10)]+['p']} for g in splits}}
            go = mock.Mock()
            go.propagate_truth.side_effect=lambda v:v
            go.propagate_scores.side_effect=lambda v:v
            policy = root/'policy.json'
            with mock.patch.object(C,'verified',return_value=(preparation,splits,predictions,{})), \
                 mock.patch.object(E.M,'load_obo',return_value=go), \
                 mock.patch.object(R,'compare_test_cohort_features',return_value={}):
                frozen=C.calibrate(inputs,root,policy,labels,splits_file,go_file)
                self.assertEqual(frozen['kinds']['EC']['threshold'],.8)
                self.assertIsNone(frozen['kinds']['GO']['threshold'])
                self.assertEqual(frozen['calibration_genomes'],['cal'])
                result=C.evaluate(inputs,root,policy,labels,splits_file,go_file,n_boot=10)
                self.assertEqual(result['scopes']['all_original']['EC']['intact']['ec_top1_counts']['labelled'],1)
                with self.assertRaises(FileExistsError):
                    C.calibrate(inputs,root,policy,labels,splits_file,go_file)
                labels.write_text(labels.read_text()+'test\tnew\tEC:3.3.3.3\texperimental\n')
                with self.assertRaisesRegex(ValueError,'changed after calibration'):
                    C.evaluate(inputs,root,policy,labels,splits_file,go_file,n_boot=10)
