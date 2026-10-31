import csv
import contextlib
import io
import json
import os
from pathlib import Path
import tempfile
import unittest
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import ec_identity as E
import price_reference_metrics as M
import evaluate_price_reference as R


def entries():
    return {
        '1.1.1.1': dict(status='transferred', targets=['1.1.1.2']),
        '1.1.1.2': dict(status='active', targets=[]),
        '1.1.1.3': dict(status='transferred', targets=['1.1.1.2','1.1.1.4']),
        '1.1.1.4': dict(status='active', targets=[]),
        '1.1.1.5': dict(status='deleted', targets=[]),
        '1.1.1.n1': dict(status='active', targets=[]),
    }


class MetricsTest(unittest.TestCase):
    def record(self, truth, values, mode='exact', integrated=True):
        return M.protein(set(truth), values, entries(), mode, integrated=integrated)

    def test_missing_candidates_keep_every_reference(self):
        r = self.record(['1.1.1.2','1.1.1.4'], {})
        s = M.summary([r])
        self.assertEqual((s['reference_pairs'],s['proteins_without_candidates']), (2,1))
        self.assertEqual(s['fixed_0_5'], dict(matched=0, unconfirmed=0, missed=2,
            reference_precision=None, reference_recall=0, reference_f1=0))

    def test_unique_transfers_merge_max_not_sum(self):
        args = (['1.1.1.1','1.1.1.2'], {'1.1.1.1':.3, '1.1.1.2':.4, '1.1.1.4':.6})
        exact = self.record(*args); current = self.record(*args, mode='current_unambiguous')
        self.assertEqual(len(exact['reference_ecs']), 2)
        self.assertEqual(current['reference_ecs'], ['1.1.1.2'])
        self.assertEqual(current['ranked_top5'][1], dict(ec='1.1.1.2', score=.4))
        self.assertEqual(current['fixed_0_5']['missed'], ['1.1.1.2'])

    def test_split_deleted_unknown_not_dropped_or_guessed(self):
        r = self.record(['1.1.1.3','1.1.1.5','9.9.9.9'],
                        {'1.1.1.2':.8,'1.1.1.5':.7,'9.9.9.9':.6}, 'current_unambiguous')
        self.assertEqual(r['reference_resolution'], dict(ambiguous=1,deleted=1,unknown=1))
        self.assertEqual(r['fixed_0_5']['missed'], ['1.1.1.3'])
        self.assertEqual(r['fixed_0_5']['unconfirmed'], ['1.1.1.2'])
        self.assertEqual(len(r['reference_ecs']), 3)

    def test_ties_and_preliminary_candidates_remain(self):
        r = self.record(['1.1.1.4'], {'1.1.1.2':.8,'1.1.1.4':.8,'1.1.1.n1':.9})
        self.assertFalse(r['top1_match']); self.assertTrue(r['top3_match'])
        self.assertEqual(r['reciprocal_rank'], 1/3)
        r = self.record(['1.1.1.4'], {'1.1.1.2':.8,'1.1.1.4':.8})
        self.assertFalse(r['top1_match']); self.assertTrue(r['top_tie_any_match'])
        self.assertFalse(r['top_tie_all_match'])
        self.assertEqual(r['top1_credit'],.5)
        self.assertEqual((r['reciprocal_rank_low'],r['reciprocal_rank_high']),(.5,1))
        self.assertEqual(M.summary([r])['top1_reference_accuracy'],.5)

    def test_pair_deltas_are_identity_not_score_changes(self):
        a = self.record(['1.1.1.2'], {'1.1.1.2':.3})
        b = self.record(['1.1.1.2'], {'1.1.1.2':.9})
        c = self.record(['1.1.1.2'], {'1.1.1.4':.9})
        self.assertEqual(M.paired({'p':a},{'p':b})['changed_top1'], 0)
        self.assertEqual(M.paired({'p':a},{'p':c})['top1_lost_reference_match'], 1)
        with self.assertRaises(ValueError): M.paired({'p':a},{})

    def test_raw_has_no_fixed_probability_threshold(self):
        r = self.record(['1.1.1.2'], {'1.1.1.2':.000001}, integrated=False)
        self.assertTrue(r['top1_match']); self.assertNotIn('fixed_0_5',M.summary([r]))

    def test_reject_invalid_scores_and_hierarchy(self):
        for value in (float('nan'),float('inf'),-1,1.01,True,'0.5'):
            with self.subTest(value=value), self.assertRaises(ValueError):
                self.record(['1.1.1.2'], {'1.1.1.2':value})
        with self.assertRaises(ValueError): self.record(['1.1.1.2'], {'1.1.-.-':.9})

    def test_paired_bootstrap_zero_and_exact_delta(self):
        before = {'g1':[True,False], 'g2':[False], 'empty':[]}
        zero = M.bootstrap_top1(before,baseline=before,repetitions=100)
        self.assertEqual((zero['estimate'],zero['ci_low'],zero['ci_high']), (0,0,0))
        after = {'g1':[True,True], 'g2':[True], 'empty':[]}
        self.assertEqual(M.bootstrap_top1(after,baseline=before)['estimate'], 2/3)


class ReaderTest(unittest.TestCase):
    def setUp(self):
        root = Path(os.environ.get('TMPDIR', str(Path.home()/'tmp'))); root.mkdir(exist_ok=True)
        self.tmp = tempfile.TemporaryDirectory(dir=root)
        self.root = Path(self.tmp.name)
        self.addCleanup(self.tmp.cleanup)

    def test_hash_validates_unlabelled_rows_and_preserves_absent_targets(self):
        p = self.root/'scores.tsv'
        p.write_text('protein_id\ttype\tfunction_id\tposterior_prob\nq\tEC\t1.1.1.2\t0.8\n')
        r = R.read_scores(p, {'p'},raw=False,predictor='x',expected_rows=1,expected_sha=R.sha(p))
        self.assertEqual(r, {'p':{}})
        old = R.sha(p); p.write_text(p.read_text().replace('0.8','0.7'))
        with self.assertRaises(ValueError): R.read_scores(p,{'p'},raw=False,predictor='x',expected_rows=1,expected_sha=old)
        p.write_text(p.read_text().replace('0.7','nan'))
        with self.assertRaises(ValueError): R.read_scores(p,{'p'},raw=False,predictor='x',expected_rows=1,expected_sha=R.sha(p))

    def test_actual_price_population(self):
        p = Path(__file__).resolve().parents[2]/'manifests/suite6_price149_manifest.tsv'
        with p.open() as f: genomes = {r['host_assembly'] for r in csv.DictReader(f,delimiter='\t')}
        refs, excluded = R.references(p,genomes)
        self.assertEqual((len(refs),sum(map(len,refs.values())),len(excluded)),(148,151,1))

    def test_duplicates_and_mixed_sources_rejected(self):
        p = self.root/'claims.jsonl'
        row = dict(protein_id='p',function_type='EC',function_id='1.1.1.2',raw_score=.2,source='topec')
        p.write_text(json.dumps(row)+'\n'+json.dumps(row)+'\n')
        with self.assertRaises(ValueError): R.read_scores(p,{'p'},raw=True,predictor='topec',expected_rows=2,expected_sha=R.sha(p))
        with self.assertRaises(ValueError): R.read_scores(p,{'p'},raw=True,predictor='mapred',expected_rows=2,expected_sha=R.sha(p))


class CollectionTest(unittest.TestCase):
    """Whole-panel fixture uses the real 149-row reference population."""
    setUp = ReaderTest.setUp
    def fixture(self):
        self.manifest = Path(__file__).resolve().parents[2]/'manifests/suite6_price149_manifest.tsv'
        with self.manifest.open() as f:
            genomes = sorted({r['host_assembly'] for r in csv.DictReader(f,delimiter='\t')})
        refs, _ = R.references(self.manifest, genomes)
        self.collection = self.root/'collection'; self.collection.mkdir()
        self.audit = self.root/'audit.json'
        R.save(self.audit, dict(schema='gspa.ec-identity-audit/1',inputs={'fixture':'fixture'},
            entries=entries(), resolutions={x:E.resolve(entries(),x) for x in entries()}))
        layout = self.root/'layout'; layout.mkdir(); inventories = {}
        for g in genomes:
            p = layout/(g+'.sequence_hashes.tsv')
            ids = [pid for gid,pid in refs if gid == g]+['unlabelled']
            p.write_text('protein_id\tlength\tsequence_sha256\n'+''.join(f'{pid}\t1\tfixture\n' for pid in ids))
            inventories[g] = dict(proteins=len(ids),sequence_hashes_sha256=R.sha(p))
        R.save(layout/'report.json',dict(assemblies=inventories))
        cells = []
        for pred in ('topec','mapred','enzbert','dpz','graphec','clean'):
            for g in genomes:
                count = sum(len(v)+1 for (gid,p),v in refs.items() if gid == g)+1
                cells.append(dict(index=len(cells),predictor=pred,assembly=g,matrix_rows=100,
                    claims=count,part=str(self.root/'parts'/pred/g)))
        R.save(self.collection/'plan.json',dict(arms=['stage0','stage1'],cells=cells,layout=str(layout)))
        R.save(self.collection/'registration.json',dict(inputs={str(layout/'report.json'):R.sha(layout/'report.json')},plan_sha256=R.sha(self.collection/'plan.json')))
        for cell in cells[:22]:
            g = cell['assembly']; part = Path(cell['part']); part.mkdir(parents=True)
            claims = [dict(protein_id=p, function_type='EC', function_id=term, raw_score=.2, source='topec')
                      for (gid,p),terms in refs.items() if gid == g for term in sorted(terms)]
            claims += [dict(protein_id=p,function_type='EC',function_id='9.9.9.9',raw_score=.3,source='topec')
                       for gid,p in refs if gid == g]
            claims += [dict(protein_id='unlabelled',function_type='EC',function_id='9.9.9.9',raw_score=.3,source='topec')]
            raw = part/'claims.jsonl'; raw.write_text(''.join(json.dumps(x)+'\n' for x in claims))
            root = self.collection/'attempts'/str(cell['index'])/'fixture'; root.mkdir(parents=True)
            R.save(root/'started.json',dict(cell=cell,registration_sha256=R.sha(self.collection/'registration.json')))
            runs = {}
            for arm in ('stage0','stage1'):
                d = root/'production'/arm; d.mkdir(parents=True)
                output = d/'output.tsv'
                output.write_text('protein_id\ttype\tfunction_id\tposterior_prob\n'+''.join(
                    f"{x['protein_id']}\tEC\t{x['function_id']}\t"+
                    ('0.7' if x['function_id']=='9.9.9.9' else ('0.6' if arm=='stage0' else '0.8'))+'\n' for x in claims))
                ledger = d/'output.tsv.stages.json'; ledger.write_text('{}\n')
                command = ['integrate']+(['--enable-priors','none'] if arm=='stage0' else [])
                R.save(d/'verified_run.json',dict(returncode=0,spec=dict(command=command,outputs={g:str(output)},inputs={str(raw):R.sha(raw)}),
                    artifacts={str(output):R.sha(output),str(ledger):R.sha(ledger)}))
                runs[arm] = dict(candidate_universe_preserved=True,input_candidates=len(claims),prediction_rows=len(claims),receipt_sha256=R.sha(d/'verified_run.json'))
            verification = root/'production/verification.json'
            R.save(verification,dict(status='verified',assembly=g,predictor='topec',matrix_rows=100,runs=runs))
            R.save(root/'completed.json',dict(status='verified',cell=cell,verification_sha256=R.sha(verification)))

    def test_complete_panel_and_digest_failure(self):
        self.fixture()
        with contextlib.redirect_stdout(io.StringIO()):
            report = R.evaluate(self.collection,self.manifest,self.audit,['topec'],self.root/'success')
        self.assertEqual(report['cells'],22)
        for mode in ('exact','current_unambiguous'):
            panel = report['summaries'][mode]['topec']
            self.assertEqual(panel['arms']['raw_export']['top1_match'],0)
            self.assertEqual(panel['arms']['stage1']['top1_match'],148)
            self.assertEqual(panel['arms']['stage1']['reference_pairs'],151)
            self.assertEqual(panel['paired']['stage0_to_stage1']['top1_gained_reference_match'],148)
        ledger = next(self.collection.glob('attempts/*/*/production/stage0/*.stages.json'))
        ledger.write_text('{"changed":true}\n')
        with contextlib.redirect_stdout(io.StringIO()), self.assertRaisesRegex(ValueError,'bound artifact changed'):
            R.evaluate(self.collection,self.manifest,self.audit,['topec'],self.root/'tampered')
        self.assertFalse((self.root/'tampered/report.json').exists())
        self.assertTrue((self.root/'tampered/failure.log').is_file())

    def test_incomplete_panel_is_not_successful_subset(self):
        self.fixture()
        next(self.collection.glob('attempts/*/*/completed.json')).unlink()
        with self.assertRaisesRegex(ValueError,'not uniquely complete'):
            R.evaluate(self.collection,self.manifest,self.audit,['topec'],self.root/'partial')
        self.assertFalse((self.root/'partial').exists())

    def test_mismatched_canonical_ids_fail_even_with_no_candidates(self):
        self.fixture()
        p = next((self.root/'layout').glob('*.sequence_hashes.tsv'))
        text = p.read_text(); identity = text.splitlines()[1].split('\t')[0]
        p.write_text(text.replace(identity,'wrong_protein'))
        with self.assertRaisesRegex(ValueError,'reference protein missing'):
            R.evaluate(self.collection,self.manifest,self.audit,['topec'],self.root/'identity-failure')

    def test_mismatched_arm_flags_rejected(self):
        left = dict(command=['integrate','--enable-priors','none'],inputs={'a':'digest'})
        right = dict(command=['integrate','--unexpected','yes'],inputs={'a':'digest'})
        with self.assertRaisesRegex(ValueError,'arms differ'):
            R.comparable_specs(left,right,self.root)
        right['command']=['integrate']; right['inputs']={'a':'different'}
        with self.assertRaisesRegex(ValueError,'arms differ'):
            R.comparable_specs(left,right,self.root)


if __name__ == '__main__':
    unittest.main()
