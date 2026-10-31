"""Source-arm x context-mode comparison on the complete observed-gap cohort.

Five sequence-evidence arms (DIAMOND homology, DIAMOND family-excluded, native
CLEAN, and CLEAN plus each DIAMOND recipe) are exported under four contexts
(SEQONLY, UNIVERSAL, PATHWAY, ORGANISM) and fitted with two model kinds: 40 fits
frozen on train/calibration before test evaluation. All 18 genomes are kept;
the single CLEAN unsupported protein stays an explicit abstention (no invented
negative claim). Supervision files, splits and settings are inherited from the
original recipes; labels/splits are compared on the examples both tables share.

This reuses a previously examined cohort that overlaps the CLEAN reference
(3,730 exact known sequences); it is not unseen-family evidence. The evaluation
stage here reports conditional candidate-table metrics only. The full-reference
primary EC endpoint and production replay are separate later stages, not run here.
"""
import argparse
from collections import Counter
import copy,csv,json,math
from numbers import Real
from pathlib import Path
import sys
import audit_clean_cohort as A
import cohort_sequence_inputs as S
import run_observed_gap_context_modes as C
F=C.F;G=C.G
ARMS=('diamond_homology','diamond_family_excluded','clean','clean_diamond_homology','clean_diamond_family_excluded')
CONTEXTS=('SEQONLY','UNIVERSAL','PATHWAY','ORGANISM')
# diamond: transfer whose DIAMOND claims are used; recipe: transfer whose supervision/settings are inherited.
ARM_SOURCES={'diamond_homology':dict(diamond='homology',clean=False,recipe='homology'),
    'diamond_family_excluded':dict(diamond='family_excluded',clean=False,recipe='family_excluded'),
    'clean':dict(diamond=None,clean=True,recipe='homology'),
    'clean_diamond_homology':dict(diamond='homology',clean=True,recipe='homology'),
    'clean_diamond_family_excluded':dict(diamond='family_excluded',clean=True,recipe='family_excluded')}
EXPECTED=dict(genomes=18,proteins=49463,native=49462,unsupported=1,exact_known=3730)
SHARED_SUPERVISION=('labels','splits','label_policy','go_obo','registry','example_clusters','locus_groups')
ZERO_PREFIXES=('ctx.','cohort.','asm.','fba.')
EXAMPLE_FIELDS=('genome_id','group','split','label','label_basis','dataset','label_evidence',
    'label_evidence_categories','label_role_basis','label_references','locus_group_id')
IDENTITY=('candidate_key','genome_id','protein_id','function_type','function_id','role','sequence_sources',
    'sequence_contracts','quality_evidence')
MANIFEST_COLUMNS=['genome_id','claims_path','out_path','taxon','essential_profile','operons_path','gaps_path','max_iter',
    'enforce_stage2','enforce_mode','enable_priors','genome_layout_path','orthogroups_path','cluster_consensus_path',
    'fba_medium','fba_gram']
SEQONLY_BLANK=('operons_path','gaps_path','genome_layout_path','orthogroups_path','cluster_consensus_path','fba_medium','fba_gram')
# Context routes removed from the SEQONLY command. GO ontology, EC->GO, KEGG tables and the taxon filter stay:
# they define function identity, not organism or cohort context.
CONTEXT_PREFIXES=('--cohort','--rlc-','--experiment-','--community-','--fba-','--agent-','--cg-','--learned-',
    '--gapseq-','--workflow-','--beam-')
CONTEXT_FLAGS={'--rlc-suggester','--features-out','--suggestions-out','--operons','--gaps','--pathway-observations',
    '--orthogroups','--cluster-consensus','--genome-layout','--rxn-locus-catalog','--reaction-graph','--diffusion-mets',
    '--reaction-ec-aliases','--external-context','--enable-fba-refiner','--enable-agent','--iterate-gapseq',
    '--max-gapseq-iter','--intragenome-cluster','--dark-matter','--promotion-strategy','--maxsat-coherence-bonus',
    '--refined-bf','--anomalies-out','--infer-taxon','--enforce-stage2','--enable-priors'}
LIMITS=('Previously examined cohort that overlaps the CLEAN reference (exact known sequences, not unseen families); '
    'development reuse limits apply. Candidate-table metrics are conditional on the exported candidates; no biological '
    'accuracy, calibration or promotion claim. SEQONLY retains the declared organism taxon and essential-profile '
    'metadata and function definitions, but disables priors and context stages.')


def require(condition,message):
    if not condition:raise ValueError(message)


def is_context_flag(token):
    flag=token.split('=',1)[0]
    if flag.startswith('--no-'):flag='--'+flag[5:]
    return flag in CONTEXT_FLAGS or flag.startswith(CONTEXT_PREFIXES)


def strip_context(command):
    """Drop every context route (flag plus value) and disable priors; keep everything else in order."""
    result=[];skip=False
    for i,token in enumerate(command):
        if skip:skip=False;continue
        if token.startswith('--') and is_context_flag(token):
            skip='=' not in token and i+1<len(command) and not command[i+1].startswith('--')
            continue
        result.append(token)
    F.replace(result,'--enable-priors','none')
    return result


def sequence_only_manifest(fields,rows):
    """Keep the batch manifest's positional layout; blank context columns, disable Stage 2."""
    require(fields==MANIFEST_COLUMNS[:len(fields)],'batch manifest positional columns changed')
    out=[]
    for row in rows:
        r={c:row.get(c) or '-' for c in MANIFEST_COLUMNS}
        r.update({c:'-' for c in SEQONLY_BLANK},enforce_stage2='false',enable_priors='none')
        out.append(r)
    return list(MANIFEST_COLUMNS),out


def verify_clean(root,audit_path,genomes,expected=EXPECTED):
    """Bind the independent CLEAN audit to its plan and claim files; the audit itself is not rerun."""
    plan_path=root/'plan.json';report=A.read(audit_path);plan=A.read(plan_path)
    require(report['schema']=='gspa.clean-cohort-audit/2' and plan['schema']=='gspa.clean-cohort/2','CLEAN schema 2 required')
    require(report['plan_sha256']==S.sha(plan_path) and report['auditor_sha256']==S.sha(A.__file__),'audit does not bind this plan/auditor')
    ids=[g['genome_id'] for g in plan['genomes']]
    require(len(ids)==len(set(ids))==expected['genomes'] and set(ids)==set(genomes) and report['original_genomes']==plan['genomes'],
        'CLEAN genome population differs from the original cohort')
    require(report['genomes']==len(report['panels'])==expected['genomes'] and report['proteins']==expected['proteins']
        ==sum(g['proteins'] for g in plan['genomes']) and report['native_proteins']==expected['native']
        and report['unsupported_proteins']==expected['unsupported']==len(report['unsupported_outcomes'])
        and expected['native']+expected['unsupported']==expected['proteins'],'CLEAN protein accounting differs')
    require({p['genome_id'] for p in report['panels']}==set(ids), 'CLEAN panel genome coverage differs')
    unsupported={g:{r['protein_id'] for r in rows} for g,rows in plan['unsupported_inputs'].items()}
    outcomes=report['unsupported_outcomes']
    require({(o['genome_id'],o['protein_id']) for o in outcomes}=={(g,p) for g,ps in unsupported.items() for p in ps}
        and all(o['status']=='unsupported_input' and o['prediction']=='abstain' for o in outcomes),'unsupported outcomes differ from plan')
    overlap=A.read(root/'reference-overlap.json')
    require(S.sha(root/'reference-overlap.json')==plan['reference_overlap_sha256'] and overlap==report['reference_overlap']
        and overlap['exact_proteins']==expected['exact_known']==sum(v['exact_proteins'] for v in overlap['by_genome'].values()),
        'reference overlap differs from the audited acquisition report')
    allowed=dict(source='clean',contract=plan['predictor_contract']);S.identities([allowed])
    claims={};paths=[plan_path,audit_path,root/'reference-overlap.json']
    for panel in report['panels']:
        gid=panel['genome_id'];original=next(g for g in plan['genomes'] if g['genome_id']==gid)
        folder=root/'attempts'/gid/Path(panel['attempt']).name
        claim=folder/'claims.jsonl'
        require(S.sha(claim)==panel['claims_sha256'] and S.sha(folder/'completed.json')==panel['completed_sha256'],'CLEAN panel hash differs')
        require(panel['proteins']+len(unsupported.get(gid,()))==original['proteins']==panel['original_proteins'],'CLEAN genome denominator differs')
        with claim.open() as f:
            for line in f:
                require(A.decode(line)['protein_id'] not in unsupported.get(gid,()),'claim invented for an unsupported protein')
        claims[gid]=dict(path=str(claim.resolve()),sha256=panel['claims_sha256']);paths+=[claim,folder/'completed.json']
    summary={k:v for k,v in overlap.items() if k!='rows'}
    return dict(allowed=allowed,claims=claims,inputs=F.hashes(paths),summary=dict(genomes=ids,proteins=report['proteins'],
        native_proteins=report['native_proteins'],unsupported_proteins=report['unsupported_proteins'],unsupported_outcomes=outcomes,
        original_genomes=plan['genomes'],failed_attempts=len(report['failed_attempts']),reference_overlap=summary,
        plan_sha256=report['plan_sha256'],audit=str(audit_path),root=str(root)))


def supervision_agreement(recipes):
    """Both transfer recipes claim the same truth sources; confirm those bytes really agree."""
    a,b=(recipes[t]['supervision'] for t in G.MODES);agree={}
    for key in sorted(set(a)|set(b)):
        agree[key]=bool(a.get(key) and b.get(key) and F.T.sha256_file(a[key])==F.T.sha256_file(b[key]))
    for key in SHARED_SUPERVISION:
        if a.get(key) or b.get(key):
            require(agree.get(key),'supervision dependency differs between recipes: '+key)
    return agree


def arm_claims(arm,base,clean,out,genomes):
    """Whole-file concatenation of the arm's registered claim files per genome; contracts per source."""
    spec=ARM_SOURCES[arm];allowed=[];files={}
    if spec['diamond']:allowed+=base['contracts'][spec['diamond']]
    if spec['clean']:allowed.append(clean['allowed'])
    directory=out/'claims'/arm;directory.mkdir(parents=True,exist_ok=False)
    for gid in genomes:
        inputs=[]
        if spec['diamond']:inputs.append({k:base['claims'][spec['diamond']][gid][k] for k in ('path','sha256')})
        if spec['clean']:inputs.append(clean['claims'][gid])
        files[gid]=S.merge(inputs,gid,allowed,directory/(gid+'.claims.jsonl'))
    return files,allowed


def make_plan(arm,context,base,claims,jar,context_dir,out,sources):
    baseline=base['plans'][ARM_SOURCES[arm]['recipe']]['observed_gaps']
    command=list(baseline['command'])
    require('--cohort' in command and '--rlc-suggester' in command and '--learned-model' not in command,
        'observed-gap unscored joint-cohort baseline required')
    with Path(baseline['manifest']).open() as f:
        r=csv.DictReader(f,delimiter='\t');fields=r.fieldnames;original=list(r)
    require(len(original)==len(base['genomes']) and {x['genome_id'] for x in original}==set(base['genomes']),'baseline genome coverage differs')
    require(all(x['gaps_path']==base['gaps'][x['genome_id']]['path'] for x in original),'baseline lacks the registered observed gaps')
    arm_dir=out/arm/context;inference=arm_dir/'inference';inference.mkdir(parents=True,exist_ok=False)
    rows=[dict(x,claims_path=claims[x['genome_id']]['path'],out_path=str(inference/(x['genome_id']+'.integrated.tsv'))) for x in original]
    manifest=arm_dir/'manifest.tsv'
    if context=='SEQONLY':
        fields,rows=sequence_only_manifest(fields,rows);command=strip_context(command)
    with manifest.open('w') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    flags=[('-jar',jar),('--batch-manifest',manifest),('--feature-export-dir',arm_dir/'features'),('--batch-stats-out',inference/'batch_stats.tsv')]
    if context!='SEQONLY':
        flags+=[('--features-out',inference/'rlc-features.tsv'),('--suggestions-out',inference/'rlc-suggestions.tsv'),
            ('--cohort-context-mode',context),('--cohort-pathways',context_dir/'pathway-groups.tsv'),
            ('--cohort-metadata',context_dir/'organism-metadata.json')]
    for flag,value in flags:F.replace(command,flag,value)
    F.save(arm_dir/'command.json',command)
    inputs=F.V.describe_run(command,manifest,sources=sources)['inputs'];inputs.update(F.hashes([arm_dir/'command.json']))
    return dict(arm=str(arm_dir),features=str(arm_dir/'features'),inference=str(inference),manifest=str(manifest),command=command),inputs


def base_inputs(training,context):
    """The old runner's verified frozen-input chain, reused unchanged."""
    fitting=C.B.verify(training);frozen=json.loads((training/'training.json').read_text())
    require(frozen.get('status')=='verified','verified reference training required')
    F.check_hashes(frozen['outputs'])
    comparison=Path(fitting['comparison']);base=G.verify(comparison)
    exported=json.loads((comparison/'verification.json').read_text())
    require(exported.get('status')=='verified','complete observed-gap export required')
    F.check_hashes(exported['inputs']);F.check_hashes(exported['outputs'])
    receipt=json.loads((context/'receipt.json').read_text())
    require(receipt.get('status')=='verified' and receipt['genomes']==len(base['genomes']),'complete verified context acquisition required')
    F.check_hashes(receipt['inputs']);F.check_hashes({str(context/k):v for k,v in receipt['outputs'].items()})
    meta=json.loads((context/'organism-metadata.json').read_text())
    require(set(meta['genomes'])==set(base['genomes']),'context metadata genome coverage differs')
    inputs=F.hashes([training/'training.json',training/'protocol.json',training/'registration.json',comparison/'verification.json',context/'receipt.json'])
    for part in (frozen['outputs'],exported['inputs'],exported['outputs'],receipt['inputs'],{str(context/k):v for k,v in receipt['outputs'].items()}):
        inputs.update(part)
    return fitting,base,inputs


def prepare(training,context,jar,clean_root,clean_audit,out,expected=EXPECTED):
    fitting,base,inputs=base_inputs(training,context)
    clean=verify_clean(clean_root,clean_audit,base['genomes'],expected);inputs.update(clean['inputs'])
    recipes=copy.deepcopy(fitting['recipes'])
    for t in G.MODES:recipes[t]['original_table']=str(training/t/'default/table/training_table.jsonl')
    inputs.update(C.supervision_inputs(recipes,json.loads((training/'registration.json').read_text())['inputs']))
    agreement=supervision_agreement(recipes)
    out.mkdir(parents=True,exist_ok=False)
    sources=C.sources();plans={};claims={};contracts={}
    for arm in ARMS:
        claims[arm],contracts[arm]=arm_claims(arm,base,clean,out,base['genomes']);plans[arm]={}
        for ctx in CONTEXTS:
            plans[arm][ctx],found=make_plan(arm,ctx,base,claims[arm],jar,context,out,sources);inputs.update(found)
    protocol=dict(schema='gspa.source-arm-context-comparison/1',scope=__doc__,arms=list(ARMS),contexts=list(CONTEXTS),
        kinds=list(F.KINDS),arm_sources=ARM_SOURCES,plans=plans,claims=claims,contracts=contracts,recipes=recipes,
        supervision_agreement=agreement,clean=clean['summary'],genomes=base['genomes'],gaps=base['gaps'],
        zero_audit=dict(prefixes=list(ZERO_PREFIXES),contexts=['SEQONLY'],required=True),limitations=LIMITS,
        python_version=sys.version,numpy_version=F.T.np.__version__,comparison=str(fitting['comparison']),
        context_inputs=str(context),labels_read_for_export=False,model_promotion='none',
        pending_stages=['full_reference_ec_endpoint','production_replay'])
    F.save(out/'protocol.json',protocol)
    F.save(out/'registration.json',dict(status='registered',registered_at=F.datetime.now(F.timezone.utc).isoformat(),
        protocol_sha256=F.T.sha256_file(out/'protocol.json'),inputs=inputs))
    return verify(out)


def verify(out):
    reg=json.loads((out/'registration.json').read_text());F.check_hashes(reg['inputs'])
    require(F.T.sha256_file(out/'protocol.json')==reg['protocol_sha256'],'registered comparison protocol changed')
    p=json.loads((out/'protocol.json').read_text());C.supervision_inputs(p['recipes'],reg['inputs'])
    require(p['arms']==list(ARMS) and p['contexts']==list(CONTEXTS) and p['kinds']==list(F.KINDS)
        and p['python_version']==sys.version and p['numpy_version']==F.T.np.__version__,'registered runtime or comparison scope changed')
    return p


def zero_audit(features,prefixes=ZERO_PREFIXES):
    """Every context column of every SEQONLY candidate is exactly zero, and its mode is UNIVERSAL."""
    checked=0
    for file in sorted(Path(features).glob('*/candidate_features.jsonl')):
        with file.open() as f:
            for line in f:
                row=S.decode(line)
                require(row['cohort_context_mode']=='UNIVERSAL','SEQONLY candidate not UNIVERSAL')
                for name,value in row['features'].items():
                    if name.startswith(tuple(prefixes)):
                        require(not isinstance(value,bool) and isinstance(value,Real) and math.isfinite(value) and value==0,
                            'SEQONLY context column is nonzero: '+name)
                        checked+=1
    require(checked,'no context columns were audited')
    return checked


def sequence_stage_audit(inference,genomes,learned=False):
    allowed={'reference_preparation','evidence_integration','feature_export'}
    if learned:allowed.add('learned_integration')
    observed={}
    for genome in genomes:
        data=F.V.ledger_metadata(Path(inference)/(genome+'.integrated.tsv.stages.json'))
        ran={s['name'] for s in data['stages'] if s.get('status')=='executed'}
        require(ran<=allowed,'SEQONLY executed a context stage: '+str(sorted(ran-allowed)))
        require({'reference_preparation','evidence_integration','feature_export'}<=ran,'SEQONLY required stage missing')
        observed[genome]=sorted(ran)
    return observed


def read_candidates(path):
    rows={}
    with Path(path).open() as f:
        for line in f:
            row=S.decode(line);require(row['candidate_key'] not in rows,'duplicate candidate')
            rows[row['candidate_key']]=row
    return rows


def compare_candidates(before,after,mode,may_change):
    """Compare by candidate key, never line order; identity/provenance must match, `may_change` features may differ."""
    a,b=read_candidates(before),read_candidates(after)
    require(a.keys()==b.keys(),'candidate universe differs')
    changes=Counter()
    for key,y in b.items():
        x=a[key]
        require(all(x.get(f)==y.get(f) for f in IDENTITY),'candidate identity or provenance changed')
        require(y['cohort_context_mode']==(mode if mode!='UNIVERSAL' and y['features']['cohort.w.ec_direct'] else 'UNIVERSAL'),
            'wrong candidate context mode')
        require(x['features'].keys()==y['features'].keys(),'feature contract changed')
        for name,value in x['features'].items():
            other=y['features'][name]
            require(all(not isinstance(v,bool) and isinstance(v,Real) and math.isfinite(v) for v in (value,other)),'nonfinite or invalid feature')
            if value!=other:
                require(may_change(name),'protected feature changed: '+name)
                changes[name]+=1
    return dict(candidates=len(b),changed_columns=dict(changes))


def supervision_identity(row):
    """Match semantic candidates across the EC-prefix correction, never across roles."""
    fields=('genome_id','protein_id','function_type','function_id','role','candidate_key','example_id')
    require(all(isinstance(row.get(k),str) and row[k] for k in fields),'missing supervised identity')
    parts=row['candidate_key'].rsplit('|',3)
    require(len(parts)==4 and row['example_id']==row['genome_id']+'|'+row['candidate_key'],
            'invalid supervised key encoding')
    fid=S.function_id(row['function_type'],row['function_id'])
    require((parts[0],parts[1],S.function_id(parts[1],parts[2]),parts[3])==
            (row['protein_id'],row['function_type'],fid,row['role']),'supervised identity fields disagree')
    return row['genome_id'],row['protein_id'],row['function_type'],fid,row['role']


def compare_supervision(original,current,context,splits=None):
    """Labels/splits must agree on shared examples; candidate rows may differ because sources differ."""
    old={}
    with Path(original).open() as f:
        for line in f:
            row=json.loads(line);key=supervision_identity(row)
            require(key not in old,'duplicate canonical supervised example')
            old[key]=dict(fields={k:row.get(k) for k in EXAMPLE_FIELDS},example_id=row['example_id'])
    genome_split={v['fields']['genome_id']:v['fields']['split'] for v in old.values()}
    if splits is not None:
        with Path(splits).open() as stream:
            assignments=list(csv.DictReader(stream,delimiter='\t'))
        genome_split={r['genome_id']:r['split'] for r in assignments}
        require(len(genome_split)==len(assignments) and all(s in F.T.SPLITS for s in genome_split.values()),'invalid split registry')
    shared=0;new_only=0;representation_changes=0;seen=set()
    with Path(current).open() as f:
        for line in f:
            row=json.loads(line)
            require(row['cohort_context_mode']==(context if context not in ('SEQONLY','UNIVERSAL') and row['features']['cohort.w.ec_direct'] else 'UNIVERSAL'),
                'wrong supervised context mode')
            require(genome_split.get(row['genome_id'])==row['split'],'genome split changed')
            key=supervision_identity(row)
            require(key not in seen, 'duplicate canonical supervised example')
            seen.add(key)
            before=old.get(key)
            if before is None:new_only+=1;continue
            require(before['fields']=={k:row.get(k) for k in EXAMPLE_FIELDS},'supervision identity, split or label changed')
            representation_changes+=before['example_id']!=row['example_id']
            shared+=1
    # A source may supply no label-matched candidates in a genome. Its proteins
    # still belong to the manifest/full-reference denominator, not this table.
    return dict(shared_examples=shared,original_only=len(old)-shared,new_only=new_only,
                representation_changes=representation_changes)


def export(out):
    p=verify(out)
    if (out/'export.json').exists():raise FileExistsError('export refuses overwrite')
    sources=C.sources();summaries={}
    for arm in ARMS:
        summaries[arm]={};plans=p['plans'][arm]
        for ctx in CONTEXTS:
            plan=plans[ctx]
            F.V.run_batch(plan['command'],Path(plan['manifest']),Path(plan['inference']),sources=sources)
            counts=S.audit_export(Path(plan['features']),p['claims'][arm],p['contracts'][arm])
            require(set(counts)==set(p['genomes']),'export genome coverage differs')
            summary=dict(candidates=sum(counts.values()))
            if ctx=='SEQONLY':
                summary['zero_context_values']=zero_audit(plan['features'],p['zero_audit']['prefixes'])
                summary['executed_stages']=sequence_stage_audit(plan['inference'],p['genomes'])
            else:
                for gid in p['genomes']:
                    expected=G.R.canonical(G.R.jsonl(Path(p['gaps'][gid]['path'])))
                    require(G.R.canonical(G.R.jsonl(Path(plan['inference'])/(gid+'.integrated.tsv.gaps.jsonl')))==expected,'observed gap records changed')
            summaries[arm][ctx]=summary
            print(json.dumps(dict(arm=arm,context=ctx,genomes=len(counts),**summary)),flush=True)
        comparison={}
        for gid in p['genomes']:
            files={c:Path(plans[c]['features'])/gid/'candidate_features.jsonl' for c in CONTEXTS}
            comparison[gid]={m:compare_candidates(files['UNIVERSAL'],files[m],m,lambda n:n.startswith('cohort.')) for m in ('PATHWAY','ORGANISM')}
            comparison[gid]['SEQONLY_vs_UNIVERSAL']=compare_candidates(files['SEQONLY'],files['UNIVERSAL'],'UNIVERSAL',lambda n:n.startswith(ZERO_PREFIXES))
        summaries[arm]['comparison']=comparison;F.save(out/arm/'comparison.json',comparison)
    verify(out)
    F.save(out/'export.json',dict(status='verified',summaries=summaries,outputs=F.hashes([q for q in out.rglob('*') if q.is_file()])))


def train(out):
    p=verify(out);e=json.loads((out/'export.json').read_text())
    require(e.get('status')=='verified','complete verified exports required')
    F.check_hashes(e['outputs'])
    if any((out/a/c/'table').exists() or (out/a/c/'models').exists() for a in ARMS for c in CONTEXTS):raise FileExistsError('training refuses overwrite')
    counts={};outputs=[]
    for arm in ARMS:
        counts[arm]={};recipe=p['recipes'][ARM_SOURCES[arm]['recipe']];inputs=recipe['supervision']
        for ctx in CONTEXTS:
            d=out/arm/ctx;contract_path=next((d/'features').glob('*/feature_contract.json'))
            table=F.P.build_table(d/'features',contract_path,Path(inputs['labels']),Path(inputs['splits']),Path(inputs['label_policy']),
                d/'table',registry_path=Path(inputs['registry']),go_obo=inputs['go_obo'],example_clusters=inputs['example_clusters'],
                locus_groups=inputs.get('locus_groups'),disjoint=tuple(recipe['disjoint_keys']))
            preserved=compare_supervision(Path(recipe['original_table']),Path(table['table']),ctx,inputs['splits'])
            if not ARM_SOURCES[arm]['clean']:
                require(preserved['original_only']==preserved['new_only']==0,'DIAMOND-only supervision universe changed')
            contract=F.T.load_contract(contract_path);rows=F.T.read_table(Path(table['table']),contract);F.T.check_disjoint(rows,['group','genome_id'])
            splits={s:[r for r in rows if r['split']==s] for s in F.T.SPLITS};(d/'models').mkdir()
            for kind in F.KINDS:
                model=d/'models'/f'model_{kind}.json'
                artifact=F.T.train(splits['train'],splits['calibration'],contract,kind,**recipe['settings'][kind],identity=F.T.split_identity(rows,['group','genome_id']))
                F.T.write_frozen(artifact,model);outputs.append(model)
            outputs.extend([Path(table['table']),d/'table/training_table_manifest.json'])
            counts[arm][ctx]=dict(supervision=preserved,rows_by_split={s:len(v) for s,v in splits.items()})
    verify(out);F.check_hashes(e['outputs'])
    F.save(out/'training.json',dict(status='verified',frozen_at=F.datetime.now(F.timezone.utc).isoformat(),counts=counts,outputs=F.hashes(outputs)))


def evaluate(out):
    p=verify(out)
    if (out/'evaluation.json').exists():raise FileExistsError('evaluation refuses overwrite')
    frozen=json.loads((out/'training.json').read_text())
    require(frozen.get('status')=='verified','all forty models must be frozen')
    F.check_hashes(frozen['outputs'])
    for arm in ARMS:
        for ctx in CONTEXTS:
            for kind in F.KINDS:
                require(str((out/arm/ctx/'models'/f'model_{kind}.json').resolve()) in frozen['outputs'],'missing frozen model')
    results={}
    for arm in ARMS:
        results[arm]={}
        for ctx in CONTEXTS:
            d=out/arm/ctx;table=d/'table/training_table.jsonl'
            contract=F.T.load_contract(next((d/'features').glob('*/feature_contract.json')))
            rows=[r for r in F.T.read_table(table,contract) if r['split']=='test'];results[arm][ctx]={}
            for kind in F.KINDS:
                model=d/'models'/f'model_{kind}.json'
                results[arm][ctx][kind]=F.T.evaluate_frozen(model,contract,rows,F.T.sha256_file(model),F.T.sha256_file(table))
    verify(out);F.check_hashes(frozen['outputs'])
    F.save(out/'evaluation.json',dict(scope='Conditional candidate-table metrics on a reused reference-overlapping cohort; not the primary endpoint.',
        limitations=p['limitations'],pending_stages=dict.fromkeys(p['pending_stages'],'not_run'),results=results))


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('action',choices=['prepare','export','train','evaluate']);ap.add_argument('--out',type=Path,required=True)
    for k in ('training','context','jar','clean-root','clean-audit'):ap.add_argument('--'+k,type=Path)
    a=ap.parse_args()
    if a.action=='prepare':
        if not all((a.training,a.context,a.jar,a.clean_root,a.clean_audit)):ap.error('training, context, jar, clean-root and clean-audit required')
        prepare(a.training.resolve(),a.context.resolve(),a.jar.resolve(),a.clean_root.resolve(),a.clean_audit.resolve(),a.out.resolve())
    else:globals()[a.action](a.out.resolve())
