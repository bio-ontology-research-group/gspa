"""Migrate verified predictor identities without changing numerical supervision.

Export, refit and evaluation are separate stages. This reuses an examined cohort;
observed score ranges and reference labels are not biological calibration proof.
"""
import argparse
from collections import Counter
import csv
import itertools
import json
from pathlib import Path
import sys

import prepare_training_data as P
import verified_batch as V
from migrate_cohort_feature_contract import hashes, check_hashes

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts'))
import train_context_model as T


def save(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')


def compare_rows(previous, current, expected_contract=None):
    """Add identity and make the historical UNIVERSAL default explicit; preserve all else."""
    count = 0
    with previous.open() as a, current.open() as b:
        for left, right in itertools.zip_longest(a, b):
            if left is None or right is None:
                raise ValueError('candidate/supervision row coverage changed')
            old, new = json.loads(left), json.loads(right)
            old.setdefault('cohort_context_mode', 'UNIVERSAL')
            new.setdefault('cohort_context_mode', 'UNIVERSAL')
            identities = new.pop('sequence_contracts', None)
            if not T.validate_sequence_contracts(identities, new.get('sequence_sources')):
                raise ValueError('missing predictor contract in migrated rows')
            if expected_contract is not None and identities != expected_contract:
                raise ValueError('exported predictor identity differs from verified claims')
            if old != new:
                raise ValueError('existing candidate features, metadata or supervision changed')
            count += 1
    return count


def export(original, claims_replay, out, jar):
    out.mkdir(parents=True, exist_ok=False)
    receipt = json.loads((claims_replay / 'run.json').read_text())
    if receipt.get('status') != 'verified' or not receipt.get('enriched_claims_reversible'):
        raise ValueError('verified reversible predictor reconstruction required')
    for filename, key in [('inputs.json','inputs_sha256'), ('contracts.json','contract_manifest_sha256')]:
        if T.sha256_file(claims_replay / filename) != receipt[key]:
            raise ValueError('predictor reconstruction receipt changed')
    claim_inputs = json.loads((claims_replay / 'inputs.json').read_text())
    check_hashes(claim_inputs)
    enriched = receipt['output_claims']['homology']
    expected = [dict(source='swissprot_diamond_baseline',
        contract=json.loads((claims_replay / 'contracts.json').read_text())['homology']['identity'])]
    old_manifest = original / 'table/training_table_manifest.json'
    old = json.loads(old_manifest.read_text()); inputs = old['inputs']
    if T.sha256_file(old['table']) != old['table_sha256']:
        raise ValueError('original training table changed')
    paths = [old_manifest, Path(old['table']), jar, original/'command.json', original/'manifest.tsv',
             claims_replay/'run.json', claims_replay/'contracts.json', claims_replay/'inputs.json']
    for key in ('contract','labels','splits','label_policy','registry','go_obo','example_clusters','locus_groups'):
        if inputs.get(key):
            p = Path(inputs[key]); paths.append(p)
            if inputs.get(key+'_sha256') and T.sha256_file(p) != inputs[key+'_sha256']:
                raise ValueError('supervision input changed: '+key)
    paths += [p for p in Path(inputs['feature_dir']).rglob('*') if p.is_file()]
    paths += list((original/'models').glob('model_*.json'))
    paths += [Path(m.__file__) for m in tuple(sys.modules.values()) if getattr(m,'__file__',None)
              and Path(m.__file__).suffix == '.py' and Path(m.__file__).resolve().is_relative_to(ROOT)]
    inference, features = out/'inference', out/'features'
    inference.mkdir()
    with (original/'manifest.tsv').open() as stream:
        reader = csv.DictReader(stream, delimiter='\t'); fields = reader.fieldnames; rows = list(reader)
    if {r['genome_id'] for r in rows} != set(enriched):
        raise ValueError('claim and cohort genome coverage differ')
    for row in rows:
        record = enriched[row['genome_id']]; claim = Path(record['path'])
        if T.sha256_file(claim) != record['sha256']:
            raise ValueError('enriched claims changed')
        if T.sha256_file(row['claims_path']) != record['original_sha256']:
            raise ValueError('original claims differ from reconstructed baseline')
        paths.extend([claim, Path(row['claims_path']), Path(row['out_path'])])
        row['claims_path'] = str(claim)
        row['out_path'] = str(inference/(row['genome_id']+'.integrated.tsv'))
    manifest = out/'manifest.tsv'
    with manifest.open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=fields,delimiter='\t',lineterminator='\n')
        writer.writeheader();writer.writerows(rows)
    command = json.loads((original/'command.json').read_text())
    replacements={'-jar':str(jar),'--batch-manifest':str(manifest),'--feature-export-dir':str(features),
                  '--batch-stats-out':str(inference/'batch_stats.tsv')}
    for flag,value in replacements.items():
        if flag in command:command[command.index(flag)+1]=value
    if '--learned-model' in command:raise ValueError('migration export must not score')
    save(out/'command.json',command)
    save(out/'settings.json',dict(original=str(original),claims_replay=str(claims_replay),expected_contract=expected))
    identity=hashes(paths);identity.update(claim_inputs);save(out/'inputs.json',identity)
    V.run_batch(command,manifest,inference,sources=[Path(__file__)])
    previous=sorted(Path(inputs['feature_dir']).glob('*/candidate_features.jsonl'))
    current=sorted(features.glob('*/candidate_features.jsonl'))
    if {p.parent.name for p in previous}!={p.parent.name for p in current}:
        raise ValueError('export genome coverage changed')
    n=sum(compare_rows(p,features/p.parent.name/p.name,expected) for p in previous)
    with (original/'manifest.tsv').open() as stream:
        for row in csv.DictReader(stream,delimiter='\t'):
            if Path(row['out_path']).read_bytes()!=(inference/(row['genome_id']+'.integrated.tsv')).read_bytes():
                raise ValueError('legacy annotation bytes changed')
    check_hashes(identity)
    save(out/'export.json',dict(status='verified',rows=n,numeric_features_unchanged=True,
        legacy_predictions_byte_identical=True,contract=str(current[0].parent/'feature_contract.json'),
        outputs=hashes([p for p in features.rglob('*') if p.is_file()]+[inference/V.RECEIPT,out/'settings.json',out/'command.json'])))


def verified(out):
    check_hashes(json.loads((out/'inputs.json').read_text()))
    report=json.loads((out/'export.json').read_text());check_hashes(report['outputs'])
    receipt=json.loads((out/'inference'/V.RECEIPT).read_text())
    check_hashes(receipt['spec']['inputs'])
    if not V.reusable(out/'inference',receipt['spec']):raise ValueError('stale export receipt')
    return report,json.loads((out/'settings.json').read_text())


def train(out):
    report,settings=verified(out)
    if (out/'table').exists() or (out/'models').exists():raise FileExistsError('refit refuses overwrite')
    original=Path(settings['original'])
    old=json.loads((original/'table/training_table_manifest.json').read_text());inputs=old['inputs']
    table=P.build_table(out/'features',Path(report['contract']),Path(inputs['labels']),Path(inputs['splits']),
        Path(inputs['label_policy']),out/'table',registry_path=Path(inputs['registry']),go_obo=inputs['go_obo'],
        example_clusters=inputs['example_clusters'],locus_groups=inputs.get('locus_groups'),
        disjoint=tuple(old['split_audit']['disjoint_keys']))
    n=compare_rows(Path(old['table']),Path(table['table']),settings['expected_contract'])
    contract=T.load_contract(Path(report['contract']));rows=T.read_table(Path(table['table']),contract)
    T.check_disjoint(rows,['group','genome_id'])
    splits={s:[r for r in rows if r['split']==s] for s in T.SPLITS}
    (out/'models').mkdir();summaries={}
    for kind in ('logistic','gbdt'):
        prior=json.loads((original/'models'/f'model_{kind}.json').read_text())
        art=T.train(splits['train'],splits['calibration'],contract,kind,weighting=prior['training']['weighting'],
            folds=prior['training']['cv']['folds'],seed=prior['training']['cv']['seed'],
            threshold_rule=prior['display_threshold']['rule'],identity=T.split_identity(rows,['group','genome_id']))
        for key in ('kind','normalization','logistic','gbdt','calibration','display_threshold'):
            if prior.get(key)!=art.get(key):raise ValueError('fitted scoring parameters changed: '+key)
        if not art['input_support']['available']:raise ValueError('migrated support unavailable')
        T.write_frozen(art,out/'models'/f'model_{kind}.json')
        summaries[kind]=dict(scoring_parameters_unchanged=True,input_support=art['input_support'],locus_model=art['locus_model'])
    verified(out)
    save(out/'training.json',dict(status='verified',supervision_rows_preserved=n,models=summaries,
        outputs=hashes(list((out/'models').glob('model_*.json'))+[Path(table['table']),out/'table/training_table_manifest.json']),
        limitation='Previously examined cohort. No test-based model, range or policy selection.'))


def evaluate(out):
    report,_=verified(out)
    trained=json.loads((out/'training.json').read_text());check_hashes(trained['outputs'])
    if (out/'evaluation.json').exists():raise FileExistsError('evaluation refuses overwrite')
    contract=T.load_contract(Path(report['contract']));table=out/'table/training_table.jsonl'
    test=[r for r in T.read_table(table,contract) if r['split']=='test'];results={}
    for kind in ('logistic','gbdt'):
        model=out/'models'/f'model_{kind}.json'
        results[kind]=T.evaluate_frozen(model,contract,test,T.sha256_file(model),T.sha256_file(table))
    check_hashes(trained['outputs']);save(out/'evaluation.json',results)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('action',choices=('export','train','evaluate'));p.add_argument('--out',type=Path,required=True)
    for k in ('original','claims-replay','jar'):p.add_argument('--'+k,type=Path)
    a=p.parse_args()
    if a.action=='export':
        if not all((a.original,a.claims_replay,a.jar)):p.error('export requires original, claims replay and JAR')
        export(a.original.resolve(),a.claims_replay.resolve(),a.out.resolve(),a.jar.resolve())
    else:globals()[a.action](a.out.resolve())
