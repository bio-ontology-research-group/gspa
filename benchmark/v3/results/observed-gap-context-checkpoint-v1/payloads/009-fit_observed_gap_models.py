"""Fit both model families on observed-gap exports using frozen supervision.

This is reference-support development on a previously examined cohort. No
independent biological calibration, novel-family or trained-null claim is made.
Both transfer regimes and both model kinds are frozen before test evaluation.
"""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import sys

import run_observed_gap_cohort as G
import migrate_cohort_feature_contract as M
F=G.F


def register(comparison, templates, out):
    protocol=G.verify(comparison)
    verified=json.loads((comparison/'verification.json').read_text())
    if verified.get('status')!='verified':raise ValueError('complete verified gap comparison required')
    F.check_hashes(verified['inputs']);F.check_hashes(verified['outputs'])
    sources=[Path(m.__file__).resolve() for m in tuple(sys.modules.values())
             if getattr(m,'__file__',None) and Path(m.__file__).suffix=='.py'
             and Path(m.__file__).resolve().is_relative_to(F.ROOT)]
    inputs=F.hashes(sources+[comparison/'verification.json',comparison/'protocol.json'])
    inputs.update(verified['inputs']);inputs.update(verified['outputs'])
    recipes={}
    for mode in G.MODES:
        template=templates[mode];manifest=template/'table/training_table_manifest.json'
        prior=json.loads((template/'training.json').read_text())
        if prior.get('status')!='verified':raise ValueError('verified template training required')
        F.check_hashes(prior['outputs'])
        inputs.update(F.hashes([template/'training.json']))
        original=json.loads(manifest.read_text());old=original['inputs'];supervision={}
        if F.T.sha256_file(original['table'])!=original['table_sha256']:
            raise ValueError('original supervision table changed')
        inputs.update(F.hashes([manifest,Path(original['table'])]))
        for key in ('labels','splits','label_policy','registry','go_obo','example_clusters','locus_groups'):
            if old.get(key):
                p=Path(old[key]).resolve()
                if old.get(key+'_sha256') and F.T.sha256_file(p)!=old[key+'_sha256']:
                    raise ValueError('original supervision input changed: '+key)
                supervision[key]=str(p);inputs.update(F.hashes([p]))
        settings={}
        for kind in F.KINDS:
            model=template/'models'/f'model_{kind}.json'
            settings[kind]=F.training_settings(json.loads(model.read_text()))
            inputs.update(F.hashes([model]))
        recipes[mode]=dict(supervision=supervision,settings=settings,original_table=original['table'],
            disjoint_keys=original['split_audit']['disjoint_keys'],
            features=protocol['plans'][mode]['observed_gaps']['features'])
    out.mkdir(parents=True,exist_ok=False)
    recipe=dict(schema='gspa.observed-gap-training/1',scope=__doc__,recipes=recipes,
        model_kinds=list(F.KINDS),modes=list(G.MODES),comparison=str(comparison),
        python_version=sys.version,numpy_version=F.T.np.__version__,
        selection='Preserve original splits and train/calibration selection rules; freeze all four models before evaluation; no test tuning or promotion.')
    F.save(out/'protocol.json',recipe)
    F.save(out/'registration.json',dict(status='registered',registered_at=datetime.now(timezone.utc).isoformat(),
        protocol_sha256=F.T.sha256_file(out/'protocol.json'),inputs=inputs))
    verify(out)
    return recipe


def verify(out):
    reg=json.loads((out/'registration.json').read_text());F.check_hashes(reg['inputs'])
    if F.T.sha256_file(out/'protocol.json')!=reg['protocol_sha256']:
        raise ValueError('registered training protocol changed')
    protocol=json.loads((out/'protocol.json').read_text())
    if (protocol['python_version']!=sys.version or protocol['numpy_version']!=F.T.np.__version__
            or protocol['modes']!=list(G.MODES) or protocol['model_kinds']!=list(F.KINDS)):
        raise ValueError('training runtime or scope changed')
    return protocol


def train(out):
    protocol=verify(out)
    if any((out/mode).exists() for mode in G.MODES):raise FileExistsError('training refuses overwrite')
    counts={};outputs=[]
    for mode in G.MODES:
        recipe=protocol['recipes'][mode];inputs=recipe['supervision'];arm=out/mode/'default'
        features=Path(recipe['features']);contract_path=next(features.glob('*/feature_contract.json'))
        table=F.P.build_table(features,contract_path,Path(inputs['labels']),Path(inputs['splits']),
            Path(inputs['label_policy']),arm/'table',registry_path=Path(inputs['registry']),
            go_obo=inputs['go_obo'],example_clusters=inputs['example_clusters'],locus_groups=inputs.get('locus_groups'),
            disjoint=tuple(recipe['disjoint_keys']))
        preserved=M.compare_supervision(Path(recipe['original_table']),Path(table['table']))
        contract=F.T.load_contract(contract_path);rows=F.T.read_table(Path(table['table']),contract)
        F.T.check_disjoint(rows,['group','genome_id'])
        splits={s:[r for r in rows if r['split']==s] for s in F.T.SPLITS}
        (arm/'models').mkdir()
        for kind in F.KINDS:
            artifact=F.T.train(splits['train'],splits['calibration'],contract,kind,
                **recipe['settings'][kind],identity=F.T.split_identity(rows,['group','genome_id']))
            model=arm/'models'/f'model_{kind}.json';F.T.write_frozen(artifact,model);outputs.append(model)
        outputs.extend([Path(table['table']),arm/'table/training_table_manifest.json'])
        counts[mode]=dict(supervision_rows_preserved=preserved,rows_by_split={s:len(r) for s,r in splits.items()})
    verify(out)
    F.save(out/'training.json',dict(status='verified',frozen_at=datetime.now(timezone.utc).isoformat(),
        counts=counts,outputs=F.hashes(outputs)))
    return counts


def evaluate(out):
    protocol=verify(out)
    if (out/'evaluation.json').exists():raise FileExistsError('evaluation refuses overwrite')
    frozen=json.loads((out/'training.json').read_text());F.check_hashes(frozen['outputs'])
    if frozen['status']!='verified':raise ValueError('all four models must be frozen before evaluation')
    for mode in G.MODES:
        for kind in F.KINDS:
            model=out/mode/'default/models'/f'model_{kind}.json'
            if str(model.resolve()) not in frozen['outputs']:raise ValueError('missing frozen model')
    results={}
    for mode in G.MODES:
        arm=out/mode/'default';table=arm/'table/training_table.jsonl'
        contract=F.T.load_contract(next(Path(protocol['recipes'][mode]['features']).glob('*/feature_contract.json')))
        rows=[r for r in F.T.read_table(table,contract) if r['split']=='test']
        results[mode]={}
        for kind in F.KINDS:
            model=arm/'models'/f'model_{kind}.json'
            results[mode][kind]=F.T.evaluate_frozen(model,contract,rows,F.T.sha256_file(model),F.T.sha256_file(table))
    F.check_hashes(frozen['outputs']);verify(out)
    F.save(out/'evaluation.json',results)
    return results


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=['register','train','evaluate'])
    p.add_argument('--out',type=Path,required=True)
    for name in ('comparison','homology-template','family-template'):p.add_argument('--'+name,type=Path)
    a=p.parse_args()
    if a.action=='register':
        if not all((a.comparison,a.homology_template,a.family_template)):p.error('all input paths required')
        result=register(a.comparison.resolve(),dict(homology=a.homology_template.resolve(),
            family_excluded=a.family_template.resolve()),a.out.resolve())
    else:result=globals()[a.action](a.out.resolve())
    print(json.dumps(result,indent=2))
