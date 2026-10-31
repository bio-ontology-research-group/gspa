"""Independent arithmetic check of retained Price reference score tables.

Uses no evaluator implementation. Verifies reduced-score arithmetic, not the
remote extraction from full files; that extraction is bound and rehashed by the
production evaluation. Exact and unique-transfer policies are checked separately.
"""
import argparse
from collections import defaultdict
import csv
import hashlib
import json
from pathlib import Path


def sha(p):
    with Path(p).open('rb') as f:
        return hashlib.file_digest(f,'sha256').hexdigest()


def check(root, audit, manifest):
    report=json.loads((root/'report.json').read_text())
    registration=json.loads((root/'registration.json').read_text())
    assert report['status']=='verified'
    assert sha(root/'registration.json')==report['registration_sha256']
    assert sha(audit) in registration['inputs'].values()
    assert sha(manifest) in registration['inputs'].values()
    source=json.loads(audit.read_text())['resolutions']
    refs={}
    with manifest.open() as f:
        for r in csv.DictReader(f,delimiter='\t'):
            if (r['price_entry'],r['host_assembly'],r['proteome_protein_id'])==('NP_384885','GCA_000006965.1','NC_003047.1'):
                continue
            refs[r['host_assembly'],r['proteome_protein_id']]=set(r['ground_truth_ecs'].split(';'))
    assert len(refs)==148 and sum(map(len,refs.values()))==151
    projected=lambda x,mode: source[x]['canonical_id'] if mode=='current_unambiguous' and source.get(x,{}).get('status')=='resolved' else x
    counts=defaultdict(lambda:defaultdict(float)); proteins=0; panels=set()
    for relative,digest in report['retained_files'].items():
        p=root/relative; assert sha(p)==digest
        x=json.loads(p.read_text()); cell=x['cell']; pred=cell['predictor']; g=cell['assembly']
        assert cell in registration['cells']
        assert (pred,g) not in panels; panels.add((pred,g))
        local={pid:terms for (gid,pid),terms in refs.items() if gid==g}
        assert {p:set(t) for p,t in x['reference'].items()}==local
        for mode in ('exact','current_unambiguous'):
            for arm,values in x['scores'].items():
                assert set(values)==set(local)
                for pid,raw in values.items():
                    truth={projected(t,mode) for t in local[pid]}
                    scores={}
                    for old,v in raw.items():
                        key=projected(old,mode); scores[key]=max(scores.get(key,-1),v)
                    order=sorted(scores,key=lambda t:(-scores[t],t))
                    ties={t for t,v in scores.items() if v==max(scores.values())} if scores else set()
                    credit=len(ties & truth)/len(ties) if ties else 0
                    r=report['records'][mode][pred][arm][g][pid]
                    assert set(r['reference_ecs'])==truth
                    assert set(r['candidate_reference_matches'])==truth & scores.keys()
                    assert set(r['absent_reference_candidates'])==truth-scores.keys()
                    assert r['top1_credit']==credit
                    assert [p['ec'] for p in r['ranked_top5']]==order[:5]
                    for k in (1,3,5): assert r[f'top{k}_match']==bool(truth & set(order[:k]))
                    assert r['top_tie_all_match']==(bool(ties) and ties<=truth)
                    assert r['top_tie_any_match']==bool(ties & truth)
                    c=counts[mode,pred,arm]
                    c['proteins']+=1; c['reference_pairs']+=len(truth)
                    c['candidate_reference_matches']+=len(truth & scores.keys())
                    c['top1_credit']+=credit
                    for k in (1,3,5): c[f'top{k}_match']+=bool(truth & set(order[:k]))
                    c['proteins_without_candidates']+=not scores
                    if arm!='raw_export':
                        called={t for t,v in scores.items() if v>=.5}
                        parts=dict(matched=called & truth,unconfirmed=called-truth,missed=truth-called)
                        assert {k:set(v) for k,v in r['fixed_0_5'].items()}==parts
                        for key,v in parts.items(): c[key]+=len(v)
                    else: assert 'fixed_0_5' not in r
                    proteins+=1
    assert panels=={(p,g) for p in report['predictors'] for g in report['genomes']}
    for (mode,pred,arm),c in counts.items():
        s=report['summaries'][mode][pred]['arms'][arm]
        assert c['proteins']==148
        for k in ('proteins','reference_pairs','candidate_reference_matches','top1_match','top3_match','top5_match','proteins_without_candidates'):
            assert c[k]==s[k], (mode,pred,arm,k)
        assert abs(s['top1_reference_accuracy']-c['top1_credit']/148)<1e-12
        assert abs(s['top1_genome_bootstrap']['estimate']-c['top1_credit']/148)<1e-12
        if arm!='raw_export':
            for k in ('matched','unconfirmed','missed'): assert c[k]==s['fixed_0_5'][k]
            assert abs(s['fixed_0_5']['reference_recall']-c['matched']/c['reference_pairs'])<1e-12
    return dict(status='verified',panels=len(panels),protein_arm_policy_records=proteins,
        original_reference_proteins=148,original_reference_pairs=151,
        scope='Independent reduced-table identity, ranking, tie, coverage and fixed-threshold arithmetic; remote full-file extraction is separately bound by the evaluation registration.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('root','audit','manifest'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();print(json.dumps(check(a.root,a.audit,a.manifest),indent=2))
