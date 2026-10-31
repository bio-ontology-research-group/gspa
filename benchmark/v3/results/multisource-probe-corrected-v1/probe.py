"""Remote single-genome production probe; no fitted models or accuracy claims."""
import argparse, copy, csv, json, shutil, subprocess, sys
from pathlib import Path
p=argparse.ArgumentParser()
for name in ('source','fs','out','jar'):p.add_argument('--'+name,type=Path,required=True)
p.add_argument('--clean',action='store_true');args=p.parse_args()
sys.path.insert(0,str(args.source/'benchmark/v3'))
import run_source_arm_context_comparison as R
import check_source_arm_context_models as M
F=R.F;S=R.S
gid='GCA_000007005.1'
def relocate(value):
    if isinstance(value,str) and value.startswith('/home/leechuck/'):
        return str(args.fs/value.lstrip('/'))
    if isinstance(value,list):return [relocate(v) for v in value]
    if isinstance(value,dict):return {k:relocate(v) for k,v in value.items()}
    return value
args.out.mkdir(parents=True,exist_ok=False)
original=args.fs/'home/leechuck/.cache/gspa-v3-data/observed-gap-cohort-v2/protocol.json'
base=relocate(json.loads(original.read_text()))
base['genomes']=[gid]
for transfer in R.G.MODES:
    baseline=base['plans'][transfer]['observed_gaps']
    with Path(baseline['manifest']).open() as stream:
        reader=csv.DictReader(stream,delimiter='\t');fields=reader.fieldnames
        rows=[relocate(r) for r in reader if r['genome_id']==gid]
    assert len(rows)==1
    manifest=args.out/(transfer+'.input.tsv')
    with manifest.open('x') as stream:
        writer=csv.DictWriter(stream,fieldnames=fields,delimiter='\t',lineterminator='\n')
        writer.writeheader();writer.writerows(rows)
    baseline['manifest']=str(manifest)
    command=baseline['command'];command[command.index('-Xmx12g')]='-Xmx6g'
    F.replace(command,'--batch-threads','1')
context_original=args.fs/'home/leechuck/.cache/gspa-v3-data/cohort-context-modes-v2-inputs'
# This is explicitly a one-genome local-context probe. The production cohort
# metadata reader correctly rejects metadata for genomes outside that cohort.
# Keep the selected genome's facts unchanged and bind the full original source.
context=args.out/'one-genome-context';context.mkdir()
metadata=json.loads((context_original/'organism-metadata.json').read_text())
assert set(metadata)=={'genomes'}, 'probe needs an explicit policy for any phylogenetic tree'
F.save(context/'organism-metadata.json',dict(genomes={gid:metadata['genomes'][gid]}))
shutil.copyfile(context_original/'pathway-groups.tsv',context/'pathway-groups.tsv')
arms=['diamond_homology'];clean={}
if args.clean:
    root=Path('/data/hohndor/gspa-v3/clean-cohort-v3/collection')
    expected='2ffa72dca65200eb8d21baeadfb8804da846630385ab87d61b5df37272c309d5'
    assert S.sha(root/'plan.json')==expected
    plan=R.A.read(root/'plan.json');genome=next(g for g in plan['genomes'] if g['genome_id']==gid)
    folder=root/'attempts'/gid/'job-10366'
    # The immutable independent auditor uses Python 3.11's file_digest; the
    # native numerical runtime is Python 3.10. Run that auditor in system Python.
    checker="""import json,sys
from pathlib import Path
sys.path.insert(0,sys.argv[1])
import audit_clean_cohort as A
root=Path(sys.argv[2]);expected=sys.argv[3];gid=sys.argv[4]
assert A.sha(root/'plan.json')==expected
plan=A.read(root/'plan.json');genome=next(g for g in plan['genomes'] if g['genome_id']==gid)
print(json.dumps(A.verify_attempt(root/'attempts'/gid/'job-10366',genome,plan,expected)))
"""
    audit=json.loads(subprocess.check_output(['/usr/bin/python3','-c',checker,str(args.source/'benchmark/v3'),
                                              str(root),expected,gid],text=True))
    F.save(args.out/'native-pilot-audit.json',audit)
    clean=dict(allowed=dict(source='clean',contract=plan['predictor_contract']),
               claims={gid:dict(path=str(folder/'claims.jsonl'),sha256=audit['claims_sha256'])})
    arms+=['clean','clean_diamond_homology']
sources=R.C.sources();plans={};inputs=F.hashes([Path(__file__),original,args.jar,
    context_original/'organism-metadata.json',context_original/'pathway-groups.tsv'])
for arm in arms:
    claims,contracts=R.arm_claims(arm,base,clean,args.out,[gid])
    plans[arm]=dict(claims=claims,contracts=contracts,modes={})
    for mode in ('SEQONLY','UNIVERSAL','WITHOUT_RLC'):
        realmode='UNIVERSAL' if mode=='WITHOUT_RLC' else mode
        destination=args.out/'without-rlc' if mode=='WITHOUT_RLC' else args.out/'runs'
        plan,bound=R.make_plan(arm,realmode,base,claims,args.jar,context,destination,sources)
        if mode=='WITHOUT_RLC':
            plan['command'].remove('--rlc-suggester')
            F.save(Path(plan['arm'])/'command.json',plan['command'])
            bound=F.V.describe_run(plan['command'],plan['manifest'],sources=sources)['inputs']
            bound.update(F.hashes([Path(plan['arm'])/'command.json']))
        plans[arm]['modes'][mode]=plan;inputs.update(bound)
F.save(args.out/'registration.json',dict(scope=__doc__,genomes=[gid],plans=plans,inputs=inputs,
    limitations='One genome only: verifies local context/sequence controls, not cross-genome recurrence or accuracy.'))
summaries={}
for arm,entry in plans.items():
    summaries[arm]={}
    for mode,plan in entry['modes'].items():
        F.V.run_batch(plan['command'],Path(plan['manifest']),Path(plan['inference']),sources=sources)
        counts=S.audit_export(plan['features'],entry['claims'],entry['contracts'])
        summary=dict(candidates=counts)
        if mode=='SEQONLY':
            summary['zero_context_values']=R.zero_audit(plan['features'])
            summary['executed_stages']=R.sequence_stage_audit(plan['inference'],[gid])
        summaries[arm][mode]=summary
        print(json.dumps(dict(arm=arm,mode=mode,summary=summary)),flush=True)
    files={m:Path(plan['features'])/gid/'candidate_features.jsonl' for m,plan in entry['modes'].items()}
    summaries[arm]['sequence_context_comparison']=R.compare_candidates(files['SEQONLY'],files['UNIVERSAL'],
        'UNIVERSAL',lambda n:n.startswith(R.ZERO_PREFIXES))
    assert S.sha(files['UNIVERSAL'])==S.sha(files['WITHOUT_RLC']), 'RLC unexpectedly changed learned feature candidates'
    summaries[arm]['rlc_identical_feature_sha256']=S.sha(files['UNIVERSAL'])
# Namespace spelling must not split equivalent protein/EC assignments.
identity_check={}
if args.clean:
    profiles={arm:S.profiles(entry['claims'][gid]['path'],gid,entry['contracts'])[0]
              for arm,entry in plans.items()}
    shared=set(profiles['diamond_homology']) & set(profiles['clean'])
    assert len(shared)==511, 'real source overlap changed'
    assert len(profiles['clean_diamond_homology'])==29512, 'combined canonical universe changed'
    feature=Path(plans['clean_diamond_homology']['modes']['SEQONLY']['features'])/gid/'candidate_features.jsonl'
    checked=set()
    for line in feature.open():
        row=S.decode(line);key=row['protein_id'],row['function_type'],row['function_id']
        if key in shared:
            assert key not in checked, 'equivalent EC split across candidates'
            assert row['role']=='CATALYTIC'
            assert row['sequence_sources']==['clean','swissprot_diamond_baseline']
            assert row['features']['seq.n_sources_reporting']==2
            assert row['features']['seq.n_groups_reporting']==2
            checked.add(key)
    assert checked==shared
    identity_check=dict(shared_source_candidates=len(shared),verified_both_source_candidates=len(checked),
                        combined_candidates=29512,raw_claim_files_unchanged=True)

# Check the new numerical/role auditor against a real existing complete learned
# genome run. This does not refit or rerun the old model.
prior=Path('/data/hohndor/gspa-v3/context-resume-v1/jobs/10176/results')
oldbase=Path('/data/hohndor/gspa-v3/observed-gap-context-modes-v1-relocation/jobs/10163/results/observed-gap-context-modes-v1')
historical={}
for kind in F.KINDS:
    original_arm=oldbase/'homology/UNIVERSAL'
    learned=prior/'observed-gap-context-model-replay-v1/homology/UNIVERSAL'/kind
    model=prior/'observed-gap-context-modes-v1/homology/UNIVERSAL/models'/('model_'+kind+'.json')
    paths=[original_arm/'inference'/(gid+'.integrated.tsv'),learned/'inference'/(gid+'.integrated.tsv'),
           original_arm/'features'/gid/'candidate_features.jsonl',learned/'features'/gid/'candidate_features.jsonl',model]
    bound=F.hashes(paths+[Path(str(paths[1])+'.learned.json')]);inputs.update(bound)
    historical[kind]=M.audit_genome(*paths)
    F.check_hashes(bound)
F.check_hashes(inputs)
F.save(args.out/'verification.json',dict(schema='gspa.source-context-probe/2',status='verified',
    source_root=str(args.source),results_root=str(args.out),ec_identity=identity_check,
    summaries=summaries,historical_model_audit=historical,
    inputs=inputs,outputs=F.hashes([p for p in args.out.rglob('*') if p.is_file()]),scope=__doc__))
