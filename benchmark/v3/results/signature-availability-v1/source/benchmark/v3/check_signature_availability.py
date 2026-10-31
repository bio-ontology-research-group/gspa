"""Production CLI regression for untestable reaction context (synthetic data).

Compares frozen engines, keeps the full candidate universe, and checks connected
controls. This validates mechanisms, not biological prediction accuracy.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(before, after, out):
    out.mkdir(parents=True, exist_ok=False)
    java=Path(shutil.which('java')).resolve()
    inputs={str(p):digest(p) for p in [before,after,java,Path(__file__).resolve()]}
    def write(name,value):
        p=out/name;p.parent.mkdir(parents=True,exist_ok=True)
        p.write_text(value if isinstance(value,str) else json.dumps(value,indent=2,sort_keys=True)+'\n')
        inputs[str(p)]=digest(p);return p
    graph=write('graph.tsv','id\tstoichiometry\tec_numbers\n'+
        'r1\t-1:a:0:0;1:b:0:0\t1.1.1.1\n'+
        'r2\t-1:b:0:0;1:c:0:0\t2.2.2.2\n'+
        'r3\t-1:c:0:0;1:d:0:0\t3.3.3.3\n'+
        'r4\t-1:x:0:0;1:y:0:0\t4.4.4.4\n'+
        'r5\t-1:z:0:0;1:w:0:0\t5.5.5.5\n')
    aliases=write('aliases.tsv','MS ID\tOld MS ID\tExternal ID\tSource\n'+''.join(
        f'r{i}\t\t{i}.{i}.{i}.{i}\tindependent_fixture\n' for i in range(1,6)))
    groups=write('groups.tsv','genome_id\tprotein_id\torthogroup_id\n'+''.join(
        f'{g}\tp{i}\tOG{i}\n' for g in ['target','peer'] for i in [2,3,4,6]))
    claims=[]
    for pid,ec,score in [('p2','1.1.1.1',.9),('p3','2.2.2.2',.4),('p4','3.3.3.3',.9),('p6','4.4.4.4',.3)]:
        claims.append(dict(protein_id=pid,function_type='EC',function_id='EC:'+ec,source='diamond',
                           raw_score=score,calibrated_prob=score,evidence_type='SEQUENCE_SIMILARITY'))
    claimfiles={g:write(f'{g}/claims.jsonl',''.join(json.dumps(c)+'\n' for c in claims)) for g in ['target','peer']}
    layouts={g:write(f'{g}/layout.tsv','protein_id\tcontig\tstart\tend\tstrand\n'+''.join(
        f'p{i}\tchr\t{i*1000}\t{i*1000+900}\t+\n' for i in range(1,9))) for g in ['target','peer']}
    manifest=write('peers.tsv',f'genome_id\tclaims_path\tlayout_path\npeer\t{claimfiles["peer"]}\t{layouts["peer"]}\n')
    runs=[];summaries={}
    def features(directory):
        return {r['candidate_key']:r for r in map(json.loads,(directory/'features/candidate_features.jsonl').read_text().splitlines())}
    for arm,edge,minimum,with_peer in [('empty',None,1,True),('supported','r1',1,True),
                ('negative','r5',1,True),('insufficient','r1',2,True),('single',None,1,False)]:
        active=['target','peer'] if with_peer else ['target']
        metadata=write(f'{arm}/metadata.json',{'genomes':{g:dict(clade=g,strainCluster=g,
            assemblyStatus='COMPLETE',familySearchComplete=True,contigs={'chr':{'length':10000,'closed':True}}) for g in active}})
        edges=write(f'{arm}/edges.tsv',f'reaction_a\treaction_b\tsource\tdirected\nr4\t{edge}\tindependent_fixture\tfalse\n') if edge else None
        outputs={}
        for label,jar in [('before',before),('after',after)]:
            directory=out/arm/label;directory.mkdir(parents=True)
            config=dict(lite=True,enabledPriors='none',genomeId='target',cohortInference=True,
                reactionGraphFile=str(graph),reactionEcAliasesFile=str(aliases),rlcCurrencyPercentile=100.,
                cohortMetadataFile=str(metadata),cohortOrthogroupsFile=str(groups),genomeLayoutFile=str(layouts['target']),
                cohortConfig={'minAnchors':minimum,'permutations':20},featureExportDir=str(directory/'features'))
            if with_peer:config['cohortManifestFile']=str(manifest)
            if edges:config['cohortInteractionEdgesFile']=str(edges)
            cfg=write(f'{arm}/{label}/workflow.json',config)
            command=[str(java),'-Xmx2g','-jar',str(jar),'integrate','--claims',str(claimfiles['target']),
                     '--workflow-config',str(cfg),'--out',str(directory/'predictions.tsv')]
            with (directory/'run.log').open('w') as log:
                process=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,timeout=120)
            runs.append(dict(arm=arm,engine=label,command=command,exit_code=process.returncode))
            (out/'runs.json').write_text(json.dumps(runs,indent=2)+'\n')
            if process.returncode:raise ValueError(f'{arm}/{label} CLI failed')
            outputs[label]=directory
        a,b=features(outputs['before']),features(outputs['after'])
        assert a.keys()==b.keys()
        assert (outputs['before']/'predictions.tsv').read_bytes()==(outputs['after']/'predictions.tsv').read_bytes()
        focal=next(k for k,r in b.items() if r['protein_id']=='p6')
        old,new=a[focal]['features'],b[focal]['features']
        unknown=arm in ['empty','insufficient','single']
        assert old['cohort.w.local_available']==1
        assert new['cohort.w.layout_available']==1
        assert new['cohort.w.local_available']==(0 if unknown else 1)
        if with_peer:
            assert new['cohort.w.effective_n']==(0 if unknown else old['cohort.w.effective_n'])
            assert new['cohort.w.frac_unknown_context']==(1 if unknown else 0)
            assert new['cohort.w.frac_negative']==(1 if arm=='negative' else 0)
            assert new['cohort.w.frac_supported']==(1 if arm=='supported' else 0)
        else:assert new['cohort.w.available']==0
        if unknown:assert new['cohort.w.local_enrich_available']==0
        if arm in ['supported','negative']:assert a==b
        else:
            # r3 has only one neighbor too: a minimum of two makes it untestable.
            affected={'p6','p4'} if arm=='insufficient' else {'p6'}
            for k in a:
                if a[k]['protein_id'] not in affected:assert a[k]==b[k],(arm,k)
                else:assert b[k]['features']['cohort.w.local_available']==0,(arm,k)
        summaries[arm]=dict(candidates=len(a),predictions_byte_identical=True,
            changed_candidates=sum(a[k]!=b[k] for k in a),before=old,after=new)
    assert inputs=={p:digest(Path(p)) for p in inputs}
    result=dict(status='verified',scope=__doc__,inputs=inputs,runs=runs,summaries=summaries,
                outputs={str(p.relative_to(out)):digest(p) for p in out.rglob('*') if p.is_file()})
    (out/'verification.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':'verified','cli_runs':len(runs),'arms':list(summaries)}))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['before','after','out']:p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();run(a.before.resolve(),a.after.resolve(),a.out.resolve())
