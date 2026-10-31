"""Matched context-mode exports and refits with complete observed gap inputs.

All 18 genomes, both transfer regimes and all 521 observed gaps remain included.
This reused, reference-overlapping cohort has sparse gap-linked supervision; it
cannot establish biological accuracy or unseen-family missing-enzyme recovery.
All twelve models are frozen before any test evaluation. No model is promoted.
"""
import argparse
from collections import Counter
import copy,csv,itertools,json,math
from numbers import Real
from pathlib import Path
import sys
import fit_observed_gap_models as B
F=B.F
G=B.G
MODES=('UNIVERSAL','PATHWAY','ORGANISM')
CONTROLS=('default',)+MODES


def sources():
    return [Path(m.__file__).resolve() for m in tuple(sys.modules.values())
        if getattr(m,'__file__',None) and Path(m.__file__).suffix=='.py'
        and Path(m.__file__).resolve().is_relative_to(F.ROOT)]


def make_plans(protocol,context,jar,out,source_files):
    plans={};inputs={}
    for transfer in G.MODES:
        baseline=protocol['plans'][transfer]['observed_gaps']
        with Path(baseline['manifest']).open() as f:
            r=csv.DictReader(f,delimiter='\t');fields=r.fieldnames;original=list(r)
        if len(original)!=len(protocol['genomes']) or {r['genome_id'] for r in original}!=set(protocol['genomes']):
            raise ValueError('baseline genome coverage differs')
        base=baseline['command']
        if '--cohort' not in base or '--rlc-suggester' not in base or '--learned-model' in base:
            raise ValueError('observed-gap unscored joint-cohort baseline required')
        for row in original:
            if row['gaps_path']!=protocol['gaps'][row['genome_id']]['path']:
                raise ValueError('baseline lacks the registered observed gaps')
        plans[transfer]={}
        for mode in CONTROLS:
            arm=out/transfer/mode;arm.mkdir(parents=True,exist_ok=False)
            inference=arm/'inference';inference.mkdir();manifest=arm/'manifest.tsv'
            rows=[dict(r,out_path=str(inference/(r['genome_id']+'.integrated.tsv'))) for r in original]
            with manifest.open('w') as f:
                w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
            command=list(base)
            for flag,value in [('-jar',jar),('--batch-manifest',manifest),('--feature-export-dir',arm/'features'),
                ('--batch-stats-out',inference/'batch_stats.tsv'),('--features-out',inference/'rlc-features.tsv'),
                ('--suggestions-out',inference/'rlc-suggestions.tsv')]:F.replace(command,flag,value)
            if mode!='default':
                for flag,value in [('--cohort-context-mode',mode),('--cohort-pathways',context/'pathway-groups.tsv'),
                    ('--cohort-metadata',context/'organism-metadata.json')]:F.replace(command,flag,value)
            F.save(arm/'command.json',command)
            inputs.update(F.V.describe_run(command,manifest,sources=source_files)['inputs'])
            inputs.update(F.hashes([arm/'command.json']))
            plans[transfer][mode]=dict(arm=str(arm),features=str(arm/'features'),inference=str(inference),
                manifest=str(manifest),command=command,baseline=baseline)
    return plans,inputs


def compare_features(before,after,mode,exact=False):
    changes=Counter();seen=set()
    with before.open() as a,after.open() as b:
        for left,right in itertools.zip_longest(a,b):
            if left is None or right is None:raise ValueError('candidate coverage changed')
            x,y=json.loads(left),json.loads(right)
            for key in ('candidate_key','genome_id','protein_id','function_type','function_id','role',
                'sequence_sources','sequence_contracts','quality_evidence'):
                if x.get(key)!=y.get(key):raise ValueError('candidate identity or provenance changed: '+key)
            if y['candidate_key'] in seen:raise ValueError('duplicate candidate')
            seen.add(y['candidate_key'])
            expected=mode if y['features']['cohort.w.ec_direct'] else 'UNIVERSAL'
            if y['cohort_context_mode']!=expected:raise ValueError('wrong candidate context mode')
            if x['features'].keys()!=y['features'].keys():raise ValueError('feature contract changed')
            for name,value in x['features'].items():
                if any(isinstance(v,bool) or not isinstance(v,Real) or not math.isfinite(v) for v in (value,y['features'][name])):
                    raise ValueError('nonfinite or invalid context feature')
                if value!=y['features'][name]:
                    if exact or not name.startswith('cohort.'):raise ValueError('non-context feature changed: '+name)
                    changes[name]+=1
            if x['features']['ctx.gap.available']!=1 or y['features']['ctx.gap.available']!=1:
                raise ValueError('observed gaps unavailable')
    return dict(candidates=len(seen),changed_columns=dict(changes))


def compare_supervision(before,after,mode):
    n=0
    with before.open() as a,after.open() as b:
        for left,right in itertools.zip_longest(a,b):
            if left is None or right is None:raise ValueError('supervision coverage changed')
            old,new=json.loads(left),json.loads(right)
            expected=mode if new['features']['cohort.w.ec_direct'] else 'UNIVERSAL'
            if new['cohort_context_mode']!=expected:raise ValueError('wrong supervised context mode')
            for row in (old,new):row.pop('features');row.pop('cohort_context_mode',None)
            if old!=new:raise ValueError('supervision identity, split, label, role or source changed')
            n+=1
    return n


def prepare(training,context,jar,out):
    fitting=B.verify(training);frozen=json.loads((training/'training.json').read_text())
    if frozen.get('status')!='verified':raise ValueError('verified reference training required')
    F.check_hashes(frozen['outputs'])
    comparison=Path(fitting['comparison']);base=G.verify(comparison)
    exported=json.loads((comparison/'verification.json').read_text())
    if exported.get('status')!='verified':raise ValueError('complete observed-gap export required')
    F.check_hashes(exported['inputs']);F.check_hashes(exported['outputs'])
    receipt=json.loads((context/'receipt.json').read_text())
    if receipt.get('status')!='verified' or receipt['genomes']!=len(base['genomes']):raise ValueError('complete verified context acquisition required')
    F.check_hashes(receipt['inputs']);F.check_hashes({str(context/k):v for k,v in receipt['outputs'].items()})
    meta=json.loads((context/'organism-metadata.json').read_text())
    if set(meta['genomes'])!=set(base['genomes']):raise ValueError('context metadata genome coverage differs')
    out.mkdir(parents=True,exist_ok=False)
    plans,inputs=make_plans(base,context,jar,out,sources())
    inputs.update(F.hashes([training/'training.json',training/'protocol.json',training/'registration.json',
        comparison/'verification.json',context/'receipt.json']))
    inputs.update(frozen['outputs']);inputs.update(exported['inputs']);inputs.update(exported['outputs'])
    inputs.update(receipt['inputs']);inputs.update({str(context/k):v for k,v in receipt['outputs'].items()})
    recipes=copy.deepcopy(fitting['recipes'])
    for transfer in G.MODES:recipes[transfer]['original_table']=str(training/transfer/'default/table/training_table.jsonl')
    protocol=dict(schema='gspa.observed-gap-context-modes/1',scope=__doc__,plans=plans,recipes=recipes,
        modes=list(MODES),controls=list(CONTROLS),transfers=list(G.MODES),genomes=base['genomes'],gaps=base['gaps'],
        claims=base['claims'],contracts=base['contracts'],python_version=sys.version,numpy_version=F.T.np.__version__,
        comparison=str(comparison),context_inputs=str(context),labels_read_for_export=False,model_promotion='none')
    F.save(out/'protocol.json',protocol)
    F.save(out/'registration.json',dict(status='registered',registered_at=F.datetime.now(F.timezone.utc).isoformat(),
        protocol_sha256=F.T.sha256_file(out/'protocol.json'),inputs=inputs))
    return verify(out)


def verify(out):
    reg=json.loads((out/'registration.json').read_text());F.check_hashes(reg['inputs'])
    if F.T.sha256_file(out/'protocol.json')!=reg['protocol_sha256']:raise ValueError('registered context protocol changed')
    p=json.loads((out/'protocol.json').read_text())
    if p['python_version']!=sys.version or p['numpy_version']!=F.T.np.__version__ or p['modes']!=list(MODES) or p['controls']!=list(CONTROLS) or p['transfers']!=list(G.MODES):
        raise ValueError('registered runtime or comparison scope changed')
    return p


def export(out):
    p=verify(out)
    if (out/'export.json').exists():raise FileExistsError('export refuses overwrite')
    summaries={}
    for transfer in G.MODES:
        summaries[transfer]={}
        for mode in CONTROLS:
            plan=p['plans'][transfer][mode];arm=Path(plan['arm'])
            F.V.run_batch(plan['command'],Path(plan['manifest']),Path(plan['inference']),sources=sources())
            counts=F.audit_export(Path(plan['features']),p['claims'][transfer],p['contracts'][transfer])
            if set(counts)!=set(p['genomes']):raise ValueError('export genome coverage differs')
            reference=plan['baseline'] if mode=='default' else p['plans'][transfer]['default' if mode=='UNIVERSAL' else 'UNIVERSAL']
            by_genome={}
            for gid in p['genomes']:
                before=Path(reference['inference'])/(gid+'.integrated.tsv');after=Path(plan['inference'])/before.name
                if mode=='default' and F.T.sha256_file(before)!=F.T.sha256_file(after):raise ValueError('new-engine control changed functional predictions')
                r=compare_features(Path(reference['features'])/gid/'candidate_features.jsonl',
                    Path(plan['features'])/gid/'candidate_features.jsonl','UNIVERSAL' if mode=='default' else mode,exact=mode=='default')
                if G.R.canonical(G.R.jsonl(Path(str(after)+'.gaps.jsonl')))!=G.R.canonical(G.R.jsonl(Path(p['gaps'][gid]['path']))):raise ValueError('observed gap records changed')
                r.update(functional_changes=G.compare_predictions(before,after,arm/(gid+'.prediction-changes.tsv')),
                    gaps=p['gaps'][gid]['records'],rlc_rows=len(G.R.rows(Path(str(after)+'.features.tsv'))),
                    suggestions=len(G.R.rows(Path(str(after)+'.suggestions.tsv'))))
                by_genome[gid]=r
            summaries[transfer][mode]=by_genome;F.save(arm/'comparison.json',by_genome)
            print(json.dumps(dict(transfer=transfer,mode=mode,genomes=len(by_genome),candidates=sum(r['candidates'] for r in by_genome.values()))),flush=True)
    verify(out);F.save(out/'export.json',dict(status='verified',summaries=summaries,
        outputs=F.hashes([q for q in out.rglob('*') if q.is_file()])))


def train(out):
    p=verify(out);e=json.loads((out/'export.json').read_text())
    if e.get('status')!='verified':raise ValueError('complete verified context exports required')
    F.check_hashes(e['outputs'])
    if any((out/t/m/'table').exists() or (out/t/m/'models').exists() for t in G.MODES for m in MODES):raise FileExistsError('training refuses overwrite')
    counts={};outputs=[]
    for transfer in G.MODES:
        recipe=p['recipes'][transfer];inputs=recipe['supervision'];counts[transfer]={}
        for mode in MODES:
            arm=out/transfer/mode;features=arm/'features';contract_path=next(features.glob('*/feature_contract.json'))
            table=F.P.build_table(features,contract_path,Path(inputs['labels']),Path(inputs['splits']),
                Path(inputs['label_policy']),arm/'table',registry_path=Path(inputs['registry']),go_obo=inputs['go_obo'],
                example_clusters=inputs['example_clusters'],locus_groups=inputs.get('locus_groups'),disjoint=tuple(recipe['disjoint_keys']))
            preserved=compare_supervision(Path(recipe['original_table']),Path(table['table']),mode)
            contract=F.T.load_contract(contract_path);rows=F.T.read_table(Path(table['table']),contract);F.T.check_disjoint(rows,['group','genome_id'])
            splits={s:[r for r in rows if r['split']==s] for s in F.T.SPLITS};(arm/'models').mkdir()
            for kind in F.KINDS:
                model=arm/'models'/f'model_{kind}.json'
                artifact=F.T.train(splits['train'],splits['calibration'],contract,kind,**recipe['settings'][kind],identity=F.T.split_identity(rows,['group','genome_id']))
                F.T.write_frozen(artifact,model);outputs.append(model)
            outputs.extend([Path(table['table']),arm/'table/training_table_manifest.json'])
            counts[transfer][mode]=dict(supervision_rows_preserved=preserved,rows_by_split={s:len(v) for s,v in splits.items()})
    verify(out);F.check_hashes(e['outputs'])
    F.save(out/'training.json',dict(status='verified',frozen_at=F.datetime.now(F.timezone.utc).isoformat(),counts=counts,outputs=F.hashes(outputs)))


def evaluate(out):
    p=verify(out)
    if (out/'evaluation.json').exists():raise FileExistsError('evaluation refuses overwrite')
    frozen=json.loads((out/'training.json').read_text())
    if frozen.get('status')!='verified':raise ValueError('all twelve models must be frozen')
    F.check_hashes(frozen['outputs'])
    for transfer in G.MODES:
        for mode in MODES:
            for kind in F.KINDS:
                if str((out/transfer/mode/'models'/f'model_{kind}.json').resolve()) not in frozen['outputs']:raise ValueError('missing frozen model')
    results={}
    for transfer in G.MODES:
        results[transfer]={}
        for mode in MODES:
            arm=out/transfer/mode;table=arm/'table/training_table.jsonl'
            contract=F.T.load_contract(next((arm/'features').glob('*/feature_contract.json')))
            rows=[r for r in F.T.read_table(table,contract) if r['split']=='test'];results[transfer][mode]={}
            for kind in F.KINDS:
                model=arm/'models'/f'model_{kind}.json'
                results[transfer][mode][kind]=F.T.evaluate_frozen(model,contract,rows,F.T.sha256_file(model),F.T.sha256_file(table))
    verify(out);F.check_hashes(frozen['outputs']);F.save(out/'evaluation.json',results)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=['prepare','export','train','evaluate']);p.add_argument('--out',type=Path,required=True)
    for k in ('training','context','jar'):p.add_argument('--'+k,type=Path)
    a=p.parse_args()
    if a.action=='prepare':
        if not all((a.training,a.context,a.jar)):p.error('training, context and jar required')
        prepare(a.training.resolve(),a.context.resolve(),a.jar.resolve(),a.out.resolve())
    else:globals()[a.action](a.out.resolve())
