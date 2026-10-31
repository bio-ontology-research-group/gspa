"""Synthetic receipt integrity tests; no LP execution or biological validation."""
import copy
import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import audit_nutrient_community as A

ORDER=['monoculture_0','monoculture_1','permissive','native_unsupplied','no_unsupplied_uptake']


def save(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,sort_keys=True)+'\n')


def rows(path,values):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(''.join(json.dumps(v,sort_keys=True)+'\n' for v in values))


def fixture(root):
    package=root/'package';package.mkdir()
    sources=['benchmark/v3/'+n for n in ['run_nutrient_community.py','run_nutrient_phenotypes.py',
              'check_nutrient_community_cobra.py','check_community_cobra.py']]+['scripts/community_metabolism.py','scripts/metabolic_experiments.py']
    for n in sources+['models/toy.xml']:
        p=package/n;p.parent.mkdir(parents=True,exist_ok=True);p.write_text('synthetic receipt fixture '+n)
    plans=[dict(case_id='case'+str(i),arm='arm',members=[dict(id='a',blocked_reactions=['R']),dict(id='b',blocked_reactions=[])]) for i in range(4)]
    rows(package/'plan/plans.jsonl',plans);save(package/'plan/contexts.json',{})
    implementation={n:A.sha(package/n) for n in sources if 'check_' not in n}
    protocol=dict(solve_order=ORDER,max_solver_calls=5,growth_threshold=1e-6)
    frozen=dict(cases=4,plans=4,plans_sha256=A.sha(package/'plan/plans.jsonl'),contexts_sha256=A.sha(package/'plan/contexts.json'),
                implementation=implementation,models={'toy':{'sha256':A.sha(package/'models/toy.xml')}},protocol=protocol)
    save(package/'plan/protocol.json',frozen)
    files={str(p.relative_to(package)):A.sha(p) for p in package.rglob('*') if p.is_file()}
    (package/'SHA256SUMS').write_text(''.join(value+'  '+name+'\n' for name,value in sorted(files.items())))
    jobs=[]
    for index in range(2):
        job=root/('job'+str(index));shutil.copytree(package,job);jobs.append(job)
        (job/'exit-code.txt').write_text('0\n');(job/'packages.txt').write_text('cobra==fixture\n')
        (job/'input-verification.log').write_text(''.join(n+': OK\n' for n in sorted(files)))
        production=[];comparisons=[]
        for plan in plans[index::2]:
            recs=[];calls=[]
            for i,kind in enumerate(ORDER):
                members=[plan['members'][int(kind[-1])]] if kind.startswith('monoculture_') else plan['members']
                recs.append(dict(id=kind,status='completed',outcome='growth',mu=1.,organisms=[m['id'] for m in members],solver_calls=[i]))
                calls.append(dict(call=i,status='optimal'))
                comparisons.append(dict(case_id=plan['case_id'],arm=plan['arm'],kind=kind,agrees=True,
                    production_status='completed',production_mu=1.,oracle=dict(status='optimal',mu=1.,max_balance_residual=0.,
                    max_bound_violation=0.,blocked_reactions={m['id']:m['blocked_reactions'] for m in members})))
            production.append(dict(case_id=plan['case_id'],arm=plan['arm'],plan_sha256=A.digest(plan),production=recs,
                predictions={k:True for k in ORDER},solver_calls=calls,promotion='none',gene_attribution='none'))
        rows(job/'predictions/predictions.jsonl',production)
        receipt=dict(shard_index=index,shard_count=2,protocol_sha256=files['plan/protocol.json'],plans_sha256=files['plan/plans.jsonl'],
            implementation=implementation,predictions_sha256=A.sha(job/'predictions/predictions.jsonl'),counts={'True':2})
        save(job/'predictions/run.json',receipt);rows(job/'oracle/comparisons.jsonl',comparisons)
        report=dict(status='verified',comparisons=10,agree=10,gpr_comparisons=4,gpr_agree=4,protocol_sha256=files['plan/protocol.json'],
            run_sha256=A.sha(job/'predictions/run.json'),comparisons_sha256=A.sha(job/'oracle/comparisons.jsonl'),
            cobra_version='fixture',tolerances=A.TOLERANCES,sources={'/worker/'+n:files['benchmark/v3/'+n]
                for n in ['check_nutrient_community_cobra.py','check_community_cobra.py']})
        save(job/'oracle/verification.json',report)
    return package,jobs


class CommunityAuditTests(unittest.TestCase):
    def setUp(self):
        self.directory=tempfile.TemporaryDirectory();self.addCleanup(self.directory.cleanup)
        self.package,self.jobs=fixture(Path(self.directory.name))

    def audit(self):return A.audit(self.package,self.jobs,2)

    def change_comparisons(self,change):
        p=self.jobs[0]/'oracle/comparisons.jsonl';values=list(A.read_rows(p));change(values);rows(p,values)
        q=self.jobs[0]/'oracle/verification.json';report=json.loads(q.read_text());report['comparisons_sha256']=A.sha(p);save(q,report)

    def change_production(self,change):
        job=self.jobs[0];p=job/'predictions/predictions.jsonl';values=list(A.read_rows(p));change(values);rows(p,values)
        q=job/'predictions/run.json';receipt=json.loads(q.read_text());receipt['predictions_sha256']=A.sha(p)
        receipt['counts']=dict(A.Counter(str(r['predictions']['permissive']) for r in values));save(q,receipt)
        r=job/'oracle/verification.json';report=json.loads(r.read_text());report['run_sha256']=A.sha(q);save(r,report)

    def test_complete_universe_and_exact_receipts_pass(self):
        result=self.audit()
        self.assertEqual((result['status'],result['audited_comparisons'],result['agree']),('verified',20,20))
        self.assertFalse(result['observation_files_read'])

    def test_missing_and_partial_shards_are_never_complete(self):
        self.assertEqual(A.audit(self.package,self.jobs[:1],2)['missing_shards'],[1])
        (self.jobs[1]/'oracle/verification.json').unlink();(self.jobs[1]/'exit-code.txt').write_text('143\n')
        result=self.audit();self.assertEqual(result['status'],'incomplete')
        self.assertEqual(result['audited_comparisons'],10)
        self.assertEqual(result['shards'][1]['unverified_prefix_lines'],10)

    def test_bound_file_mutation_is_rejected(self):
        (self.jobs[0]/'models/toy.xml').write_text('changed')
        with self.assertRaisesRegex(ValueError,'changed bound file'):self.audit()

    def test_missing_duplicate_and_extra_comparisons_are_rejected_even_if_rehashed(self):
        original=(self.jobs[0]/'oracle/comparisons.jsonl').read_text()
        for change in (lambda x:x.pop(),lambda x:x.insert(0,copy.deepcopy(x[0])),lambda x:x.append(copy.deepcopy(x[-1]))):
            with self.subTest(change=change):
                (self.jobs[0]/'oracle/comparisons.jsonl').write_text(original);self.change_comparisons(change)
                with self.assertRaisesRegex(ValueError,'comparison'):self.audit()

    def test_endpoint_and_solver_budget_are_recomputed(self):
        self.change_production(lambda x:x[0]['predictions'].update(permissive=False))
        with self.assertRaisesRegex(ValueError,'endpoint'):self.audit()
        self.change_production(lambda x:(x[0]['predictions'].update(permissive=True),x[0]['solver_calls'].append(dict(call=5))))
        with self.assertRaisesRegex(ValueError,'budget'):self.audit()

    def test_false_agreement_flag_cannot_hide_numerical_disagreement(self):
        self.change_comparisons(lambda x:x[0]['oracle'].update(mu=2.))
        with self.assertRaisesRegex(ValueError,'reported numerical agreement'):self.audit()

    def test_honest_disagreement_is_retained_and_prevents_numerical_acceptance(self):
        self.change_comparisons(lambda x:(x[0]['oracle'].update(mu=2.),x[0].update(agrees=False)))
        p=self.jobs[0]/'oracle/verification.json';report=json.loads(p.read_text());report.update(status='disagreements',agree=9);save(p,report)
        (self.jobs[0]/'exit-code.txt').write_text('1\n');result=self.audit()
        self.assertEqual((result['status'],result['audited_comparisons'],result['agree']),('disagreements',20,19))
        self.assertEqual(len(result['disagreements']),1)

    def test_invalid_primal_and_gpr_identity_cannot_be_certified(self):
        self.change_comparisons(lambda x:x[0]['oracle'].update(max_balance_residual=1.))
        with self.assertRaisesRegex(ValueError,'primal'):self.audit()
        self.change_comparisons(lambda x:x[0]['oracle'].update(max_balance_residual=0.,blocked_reactions={'a':[]}))
        with self.assertRaisesRegex(ValueError,'GPR'):self.audit()

    def test_oracle_must_bind_the_actual_production_value(self):
        self.change_comparisons(lambda x:x[0].update(production_mu=.5))
        with self.assertRaisesRegex(ValueError,'different production'):self.audit()

    def test_scalar_tolerance_cannot_hide_a_growth_threshold_disagreement(self):
        self.change_production(lambda x:(x[0]['production'][0].update(mu=.8e-6),x[0]['predictions'].update(monoculture_0=False)))
        self.change_comparisons(lambda x:(x[0].update(production_mu=.8e-6),x[0]['oracle'].update(mu=1.2e-6)))
        result=self.audit()
        self.assertEqual((result['status'],result['agree']),('disagreements',20))
        self.assertEqual(result['endpoints'],dict(compared=20,agree=19))
        self.assertEqual(len(result['endpoint_disagreements']),1)

    def test_exit_status_and_independent_source_identity_are_checked(self):
        (self.jobs[0]/'exit-code.txt').write_text('143\n')
        with self.assertRaisesRegex(ValueError,'exit differ'):self.audit()
        (self.jobs[0]/'exit-code.txt').write_text('0\n')
        p=self.jobs[0]/'oracle/verification.json';r=json.loads(p.read_text());r['sources']['/worker/check_community_cobra.py']='0'*64;save(p,r)
        with self.assertRaisesRegex(ValueError,'source bindings'):self.audit()

    def test_duplicate_shards_are_rejected(self):
        with self.assertRaisesRegex(ValueError,'duplicate/invalid shard'):A.audit(self.package,[self.jobs[0]]*2,2)

    def test_input_mutation_during_audit_is_rejected(self):
        original=A.compare_row
        def mutate(*args):
            result=original(*args);(self.jobs[0]/'packages.txt').write_text('changed');return result
        with patch.object(A,'compare_row',side_effect=mutate),self.assertRaisesRegex(ValueError,'input changed during audit'):
            self.audit()

    def add_replay(self):
        job=self.jobs[0];old=list(A.read_rows(job/'oracle/comparisons.jsonl'))[:7]
        rows(job/'previous-oracle/comparisons.jsonl',old);(job/'previous-exit-code.txt').write_text('143\n')
        (job/'oracle-replay.sbatch').write_text('synthetic replay submission')
        record=dict(schema='gspa.nutrient-community-oracle-replay/1',previous_exit_code=143,production_rerun=False,
            oracle_prefix_reused=False,shard_index=0,shard_count=2,submission_sha256=A.sha(job/'oracle-replay.sbatch'),
            copied_inputs={n:A.sha(job/n) for n in ['predictions/run.json','predictions/predictions.jsonl',
                'previous-oracle/comparisons.jsonl','previous-exit-code.txt']})
        save(job/'oracle-replay.json',record)

    def test_interrupted_prefix_is_preserved_and_compared_not_reused_as_new_evidence(self):
        self.add_replay();result=self.audit();replay=result['shards'][0]['replay']
        self.assertEqual((result['status'],replay['prefix_records'],replay['exactly_equal']),('verified',7,7))
        self.assertEqual(result['audited_comparisons'],20)

    def test_fresh_success_cannot_hide_an_earlier_prefix_disagreement(self):
        self.add_replay();job=self.jobs[0];p=job/'previous-oracle/comparisons.jsonl'
        old=list(A.read_rows(p));old[0]['oracle']['mu']=2.;old[0]['agrees']=False;rows(p,old)
        record=json.loads((job/'oracle-replay.json').read_text());record['copied_inputs']['previous-oracle/comparisons.jsonl']=A.sha(p)
        save(job/'oracle-replay.json',record);result=self.audit()
        self.assertEqual((result['status'],result['agree']),('disagreements',20))
        self.assertEqual(len(result['replay_disagreements']),1)

    def test_replay_cannot_claim_production_reruns_or_hide_changed_inputs(self):
        self.add_replay();p=self.jobs[0]/'oracle-replay.json';record=json.loads(p.read_text());record['production_rerun']=True;save(p,record)
        with self.assertRaisesRegex(ValueError,'replay receipt'):self.audit()
        record['production_rerun']=False;save(p,record);(self.jobs[0]/'previous-oracle/comparisons.jsonl').write_text('changed')
        with self.assertRaisesRegex(ValueError,'changed bound file'):self.audit()


if __name__=='__main__':unittest.main()
