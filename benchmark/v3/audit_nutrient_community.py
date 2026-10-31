"""Audit complete community shards and their independent solver receipts.

This checks provenance, coverage and reported numerical comparisons; it does not
solve LPs, read outcomes, certify transport biology or hide disagreements.
"""
import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path

TOLERANCES = dict(objective_absolute=1e-6, objective_relative=1e-6,
                  primal_residual_and_bounds=1e-6, GLPK_primal=1e-9)


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1048576),b''):h.update(block)
    return h.hexdigest()


def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,allow_nan=False).encode()).hexdigest()


def finite(value):
    return type(value) in (int,float) and math.isfinite(value)


def read_rows(path):
    with path.open() as stream:
        for line in stream:
            yield json.loads(line)


def checked(path, inputs, expected=None):
    path=path.resolve();value=sha(path)
    if expected is not None and value!=expected:raise ValueError('changed bound file: '+str(path))
    previous=inputs.setdefault(str(path),value)
    if previous!=value:raise ValueError('input changed while auditing: '+str(path))
    return value


def package_files(root,inputs):
    checked(root/'SHA256SUMS',inputs)
    files={}
    for line in (root/'SHA256SUMS').read_text().splitlines():
        value,name=line.split('  ',1)
        p=(root/name).resolve()
        if not p.is_relative_to(root.resolve()) or name in files:
            raise ValueError('invalid/duplicate package path')
        checked(p,inputs,value);files[name]=value
    return files


def endpoint(record,threshold):
    if record.get('status')!='completed':return None
    if record.get('outcome')=='infeasible_with_maintenance':return False
    mu=record.get('mu')
    return mu>=threshold if finite(mu) else None


def agreement(actual,oracle):
    if actual.get('status')!='completed':return False
    if oracle['status']=='infeasible':return actual.get('outcome')=='infeasible_with_maintenance'
    return (oracle['status']=='optimal' and finite(actual.get('mu')) and finite(oracle.get('mu'))
            and math.isclose(actual['mu'],oracle['mu'],rel_tol=1e-6,abs_tol=1e-6))


def replay_history(job,index,shard_count,threshold,inputs):
    path=job/'oracle-replay.json'
    if not path.exists():return None
    checked(path,inputs);record=json.loads(path.read_text())
    expected={'predictions/run.json','predictions/predictions.jsonl','previous-oracle/comparisons.jsonl','previous-exit-code.txt'}
    if (record.get('schema')!='gspa.nutrient-community-oracle-replay/1' or record.get('previous_exit_code')!=143
            or record.get('production_rerun') is not False or record.get('oracle_prefix_reused') is not False
            or (record.get('shard_index'),record.get('shard_count'))!=(index,shard_count)
            or set(record.get('copied_inputs',{}))!=expected):
        raise ValueError('invalid interrupted-verification replay receipt')
    for name,value in record['copied_inputs'].items():checked(job/name,inputs,value)
    checked(job/'oracle-replay.sbatch',inputs,record['submission_sha256'])
    if (job/'previous-exit-code.txt').read_text().strip()!='143':raise ValueError('previous job was not interrupted')
    old_lines=(job/'previous-oracle/comparisons.jsonl').read_bytes().splitlines(keepends=True)
    count=exact=equivalent=trailing=0;differences=[]
    current=iter(read_rows(job/'oracle/comparisons.jsonl'))
    try:
        for i,line in enumerate(old_lines):
            try:old=json.loads(line)
            except (ValueError,UnicodeDecodeError):
                if i!=len(old_lines)-1 or line.endswith(b'\n'):raise ValueError('corrupt earlier oracle prefix')
                trailing=len(line);break
            new=next(current,None)
            if new is None or tuple(old[k] for k in ('case_id','arm','kind'))!=tuple(new[k] for k in ('case_id','arm','kind')):
                raise ValueError('replayed oracle prefix identity/order differs')
            count+=1;exact+=old==new
            a,b=old['oracle'],new['oracle']
            mu_equal=(a.get('mu') is None and b.get('mu') is None) or (finite(a.get('mu')) and finite(b.get('mu'))
                and math.isclose(a['mu'],b['mu'],rel_tol=1e-6,abs_tol=1e-6)
                and (a['mu']>=threshold)==(b['mu']>=threshold))
            same=(all(old.get(k)==new.get(k) for k in ('agrees','production_status','production_mu'))
                and a['status']==b['status'] and a['blocked_reactions']==b['blocked_reactions'] and mu_equal)
            if a['status']=='optimal':same= same and all(finite(a.get(k)) and 0<=a[k]<=1e-6 for k in ('max_balance_residual','max_bound_violation'))
            equivalent+=same
            if not same:differences.append(dict(case_id=old['case_id'],arm=old['arm'],kind=old['kind'],before=old,after=new))
    finally:current.close()
    return dict(status='verified' if equivalent==count else 'disagreements',record=record,prefix_records=count,
                exactly_equal=exact,numerically_equivalent=equivalent,trailing_fragment_bytes=trailing,differences=differences)


def compare_row(row,plan,comparisons,protocol):
    order=protocol['solve_order']
    if row.get('failure'):
        if row.get('production') or row.get('solver_calls') or row.get('predictions')!={k:None for k in order}:
            raise ValueError('failed production row invented results')
        by_kind={}
    else:
        records=row.get('production',[]);calls=row.get('solver_calls',[])
        if [r['id'] for r in records]!=order or len(calls)>protocol['max_solver_calls']:
            raise ValueError('production solve coverage/order/budget mismatch')
        if [r['call'] for r in calls]!=list(range(len(calls))):raise ValueError('solver call identities differ')
        if sorted(i for r in records for i in r['solver_calls'])!=list(range(len(calls))):
            raise ValueError('solver calls are missing or multiply attributed')
        for r in records:
            expected=[plan['members'][int(r['id'][-1])]['id']] if r['id'].startswith('monoculture_') else [m['id'] for m in plan['members']]
            if r['organisms']!=expected:raise ValueError('production organism identity differs')
        by_kind={r['id']:r for r in records}
        expected={k:endpoint(v,protocol['growth_threshold']) for k,v in by_kind.items()}
        if set(row['predictions'])!=set(expected) or any(row['predictions'][k] is not v for k,v in expected.items()):
            raise ValueError('endpoint differs from production solver result')
        if row.get('promotion')!='none' or row.get('gene_attribution')!='none':
            raise ValueError('unregistered promotion or gene attribution')
    matches=0;disagreements=[];endpoints=Counter();endpoint_disagreements=[]
    for kind in order:
        rec=next(comparisons,None)
        if rec is None or (rec['case_id'],rec['arm'],rec['kind'])!=(plan['case_id'],plan['arm'],kind):
            raise ValueError('independent comparison coverage/order differs')
        actual=by_kind.get(kind,{'status':'missing','mu':None})
        if rec['production_status']!=actual['status'] or rec['production_mu']!=actual.get('mu'):
            raise ValueError('independent comparison is bound to different production result')
        oracle=rec['oracle']
        selected=[plan['members'][int(kind[-1])]] if kind.startswith('monoculture_') else plan['members']
        if oracle['blocked_reactions']!={m['id']:m['blocked_reactions'] for m in selected}:
            raise ValueError('independent GPR identity differs')
        if oracle['status']=='optimal':
            if not finite(oracle.get('mu')) or oracle['mu']<0 or any(not finite(oracle.get(k)) or not 0<=oracle[k]<=1e-6
                    for k in ('max_balance_residual','max_bound_violation')):
                raise ValueError('optimal oracle report lacks a valid primal check')
        if oracle['status']=='infeasible' and oracle.get('mu') is not None:
            raise ValueError('infeasible oracle has numeric growth')
        agrees=agreement(actual,oracle)
        if type(rec.get('agrees')) is not bool or rec['agrees']!=agrees:
            raise ValueError('reported numerical agreement differs')
        matches+=agrees
        if not agrees:disagreements.append(rec)
        actual_endpoint=endpoint(actual,protocol['growth_threshold'])
        oracle_endpoint=(False if oracle['status']=='infeasible' else
            oracle['mu']>=protocol['growth_threshold'] if oracle['status']=='optimal' else None)
        if actual_endpoint is None or oracle_endpoint is None:endpoints['unknown']+=1
        else:
            endpoints['compared']+=1;endpoints['agree']+=actual_endpoint==oracle_endpoint
            if actual_endpoint!=oracle_endpoint:
                endpoint_disagreements.append(dict(case_id=rec['case_id'],arm=rec['arm'],kind=kind,
                    production=actual_endpoint,oracle=oracle_endpoint,production_mu=actual.get('mu'),oracle_mu=oracle.get('mu')))
    return matches,disagreements,endpoints,endpoint_disagreements


def audit(package,jobs,shard_count=6):
    if type(shard_count) is not int or shard_count<1:raise ValueError('invalid shard count')
    inputs={};checked(Path(__file__),inputs)
    files=package_files(package,inputs)
    frozen=json.loads((package/'plan/protocol.json').read_text())
    if frozen['plans_sha256']!=files['plan/plans.jsonl'] or frozen['contexts_sha256']!=files['plan/contexts.json']:
        raise ValueError('frozen plan/contexts binding differs')
    for name,value in frozen['implementation'].items():
        if files.get(name)!=value:raise ValueError('frozen implementation differs')
    if any(v['sha256'] not in files.values() for v in frozen['models'].values()):
        raise ValueError('frozen model missing from package')
    plans=list(read_rows(package/'plan/plans.jsonl'));keys=[(p['case_id'],p['arm']) for p in plans]
    if len(keys)!=len(set(keys)) or len(plans)!=frozen['plans'] or len({p['case_id'] for p in plans})!=frozen['cases']:
        raise ValueError('plan identities/count differ')
    results={};unknown=[];all_disagreements=[];all_endpoint_disagreements=[];all_endpoints=Counter();total=agree=0
    replay_disagreements=[]
    for job in jobs:
        if package_files(job,inputs)!=files:raise ValueError('job package differs')
        checked(job/'exit-code.txt',inputs);exit_code=int((job/'exit-code.txt').read_text())
        checked(job/'input-verification.log',inputs)
        if Counter((job/'input-verification.log').read_text().splitlines())!=Counter(n+': OK' for n in files):
            raise ValueError('job input verification is incomplete')
        p=job/'predictions/run.json'
        if not p.exists():unknown.append(dict(path=str(job),exit_code=exit_code,reason='production receipt missing'));continue
        checked(p,inputs);receipt=json.loads(p.read_text());index=receipt['shard_index']
        if type(index) is not int or not 0<=index<shard_count or receipt['shard_count']!=shard_count or index in results:
            raise ValueError('duplicate/invalid shard identity')
        if receipt['protocol_sha256']!=files['plan/protocol.json'] or receipt['plans_sha256']!=files['plan/plans.jsonl'] or receipt['implementation']!=frozen['implementation']:
            raise ValueError('production receipt provenance differs')
        checked(job/'predictions/predictions.jsonl',inputs,receipt['predictions_sha256'])
        rows=list(read_rows(job/'predictions/predictions.jsonl'));expected=plans[index::shard_count]
        if len(rows)!=len(expected):raise ValueError('production shard coverage differs')
        for row,plan in zip(rows,expected):
            if (row['case_id'],row['arm'],row['plan_sha256'])!=(plan['case_id'],plan['arm'],digest(plan)):
                raise ValueError('production shard plan identity/order differs')
        if dict(Counter(str(r['predictions']['permissive']) for r in rows))!=receipt['counts']:
            raise ValueError('production counts differ')
        result=dict(path=str(job),exit_code=exit_code,plans=len(rows));results[index]=result
        report_path=job/'oracle/verification.json'
        if not report_path.exists():
            result.update(status='incomplete',reason='independent verification receipt missing')
            prefix=job/'oracle/comparisons.jsonl'
            if prefix.exists():
                checked(prefix,inputs)
                with prefix.open() as stream:result['unverified_prefix_lines']=sum(1 for _ in stream)
            continue
        checked(report_path,inputs);report=json.loads(report_path.read_text())
        checked(job/'oracle/comparisons.jsonl',inputs,report['comparisons_sha256'])
        checked(job/'packages.txt',inputs)
        if 'cobra=='+report['cobra_version'] not in (job/'packages.txt').read_text().splitlines():
            raise ValueError('oracle package version differs')
        if report['protocol_sha256']!=files['plan/protocol.json'] or report['run_sha256']!=sha(p) or report['tolerances']!=TOLERANCES:
            raise ValueError('oracle receipt provenance/tolerance differs')
        sources={Path(k).name:v for k,v in report['sources'].items()}
        if len(sources)!=len(report['sources']) or sources!={n:files['benchmark/v3/'+n] for n in ('check_nutrient_community_cobra.py','check_community_cobra.py')}:
            raise ValueError('oracle source bindings differ')
        comparisons=iter(read_rows(job/'oracle/comparisons.jsonl'));matches=0;endpoints=Counter();endpoint_disagreements=[]
        for row,plan in zip(rows,expected):
            n,disagreements,ep,ep_disagree=compare_row(row,plan,comparisons,frozen['protocol'])
            matches+=n;all_disagreements+=disagreements;endpoints.update(ep);endpoint_disagreements+=ep_disagree
        if next(comparisons,None) is not None:raise ValueError('extra independent comparison')
        count=len(rows)*len(frozen['protocol']['solve_order']);gprs=sum(len(p['members']) for p in expected)
        status='verified' if matches==count else 'disagreements'
        if (report['comparisons'],report['agree'],report['gpr_comparisons'],report['gpr_agree'],report['status'],exit_code)!=(count,matches,gprs,gprs,status,0 if status=='verified' else 1):
            raise ValueError('oracle counts/status/exit differ')
        result.update(status='disagreements' if endpoint_disagreements else status,receipt_status=status,
            comparisons=count,agree=matches,endpoints=dict(endpoints),endpoint_disagreements=endpoint_disagreements)
        replay=replay_history(job,index,shard_count,frozen['protocol']['growth_threshold'],inputs)
        if replay:
            result['replay']=replay
            if replay['status']!='verified':result['status']='disagreements';replay_disagreements+=replay['differences']
        total+=count;agree+=matches;all_endpoints.update(endpoints);all_endpoint_disagreements+=endpoint_disagreements
    missing=sorted(set(range(shard_count))-set(results))
    complete=not missing and not unknown and all(r['status'] in ('verified','disagreements') for r in results.values())
    status=('verified' if total==agree and not all_endpoint_disagreements and not replay_disagreements else 'disagreements') if complete else 'incomplete'
    for path,value in inputs.items():
        if sha(path)!=value:raise ValueError('input changed during audit')
    return dict(schema='gspa.nutrient-community-audit/1',status=status,expected_shards=shard_count,
        expected_plans=len(plans),expected_comparisons=len(plans)*len(frozen['protocol']['solve_order']),
        audited_comparisons=total,agree=agree,disagreements=all_disagreements,
        endpoints=dict(all_endpoints),endpoint_disagreements=all_endpoint_disagreements,
        replay_disagreements=replay_disagreements,
        shards=results,missing_shards=missing,unidentified_shards=unknown,inputs=inputs,
        observation_files_read=False,interpretation=__doc__)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--package',type=Path,required=True);p.add_argument('--jobs',nargs='+',type=Path,required=True)
    p.add_argument('--shard-count',type=int,default=6);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    if a.out.exists():raise FileExistsError('audit refuses overwrite')
    result=audit(a.package.resolve(),[j.resolve() for j in a.jobs],a.shard_count)
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:result[k] for k in ('status','expected_comparisons','audited_comparisons','agree')}))
    raise SystemExit(0 if result['status']=='verified' else 1 if result['status']=='disagreements' else 2)
