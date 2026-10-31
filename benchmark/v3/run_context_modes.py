"""Matched production context-mode export and frozen-supervision model comparison.

This previously examined, reference-overlapping cohort measures implementation
effects and reference support. It cannot certify independent biological accuracy.
"""
import argparse
from collections import Counter
import itertools
import json
from pathlib import Path
import sys

import migrate_cohort_feature_contract as M
import run_training_cohort_ablations as C
import prepare_training_data as P
import verified_batch as V
from observe_gapsmith import digest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
import train_context_model as T

MODES = ('UNIVERSAL', 'PATHWAY', 'ORGANISM')


def write(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True)+'\n')


def replace(command, flag, value):
    if flag in command:
        command[command.index(flag)+1] = str(value)
    else:
        command.extend([flag, str(value)])


def compare_rows(previous, current, mode, *, exact=False):
    changes = Counter(); count = 0
    with previous.open() as left, current.open() as right:
        for a,b in itertools.zip_longest(left, right):
            if a is None or b is None: raise ValueError('candidate count changed')
            a,b = json.loads(a),json.loads(b)
            for key in ('candidate_key','genome_id','sequence_sources'):
                if a[key] != b[key]: raise ValueError('candidate/source identity changed')
            expected = mode if b['features']['cohort.w.ec_direct'] else 'UNIVERSAL'
            if b['cohort_context_mode'] != expected: raise ValueError('context mode missing from export')
            if a['features'].keys() != b['features'].keys(): raise ValueError('numerical contract changed')
            for key,value in a['features'].items():
                if value != b['features'][key]:
                    if exact or not key.startswith('cohort.'):
                        raise ValueError('unexpected feature change: '+key)
                    changes[key] += 1
            count += 1
    return count, dict(changes)


def export(inputs, baseline, jar, out):
    receipt = json.loads((inputs/'receipt.json').read_text())
    if receipt['status'] != 'verified': raise ValueError('unverified context inputs')
    M.check_hashes(receipt['inputs'])
    M.check_hashes({str(inputs/k):v for k,v in receipt['outputs'].items()})
    out.mkdir(parents=True, exist_ok=False)
    settings = dict(inputs=str(inputs), baseline=str(baseline), jar=str(jar))
    write(out/'settings.json', settings)
    source = [Path(__file__), Path(M.__file__), Path(C.__file__), Path(P.__file__), Path(T.__file__), inputs/'receipt.json']
    for name in ('default',) + MODES:
        arm = out/name; arm.mkdir()
        manifest = arm/'manifest.tsv'; inference = arm/'inference'; features = arm/'features'
        C.make_manifest(C.CACHE/'baseline/claims/homology', inference, manifest)
        command = json.loads((baseline/'command.json').read_text())
        for flag,value in [('-jar',jar),('--batch-manifest',manifest),('--feature-export-dir',features)]:
            replace(command,flag,value)
        if name != 'default':
            replace(command,'--cohort-context-mode',name)
            replace(command,'--cohort-pathways',inputs/'pathway-groups.tsv')
            replace(command,'--cohort-metadata',inputs/'organism-metadata.json')
        write(arm/'command.json',command)
        V.run_batch(command,manifest,inference,sources=source)
        reference = baseline if name == 'default' else out/'UNIVERSAL' if name != 'UNIVERSAL' else out/'default'
        count = 0; changes = Counter()
        old = sorted((reference/'features').glob('*/candidate_features.jsonl'))
        new = sorted(features.glob('*/candidate_features.jsonl'))
        if {p.parent.name for p in old} != {p.parent.name for p in new}: raise ValueError('genome coverage changed')
        for p in old:
            n,changed = compare_rows(p,features/p.parent.name/p.name,'UNIVERSAL' if name == 'default' else name,exact=name=='default')
            count += n; changes.update(changed)
            before = reference/'inference'/(p.parent.name+'.integrated.tsv')
            after = inference/before.name
            if digest(before) != digest(after): raise ValueError('functional predictions changed during feature export')
        write(arm/'export.json',dict(status='verified', candidates=count, changed_columns=dict(changes),
            outputs=M.hashes([p for p in features.rglob('*') if p.is_file()] + [inference/V.RECEIPT,arm/'command.json'])))
    write(out/'export.json',dict(status='verified', modes=list(MODES), default_preserves_baseline=True))


def train(out):
    if json.loads((out/'export.json').read_text())['status'] != 'verified': raise ValueError('unverified exports')
    baseline = Path(json.loads((out/'settings.json').read_text())['baseline'])
    original = json.loads((baseline/'table/training_table_manifest.json').read_text())
    inputs = original['inputs']
    for mode in MODES:
        arm = out/mode
        report = json.loads((arm/'export.json').read_text()); M.check_hashes(report['outputs'])
        record = json.loads((arm/'inference'/V.RECEIPT).read_text())
        if not V.reusable(arm/'inference',record['spec']): raise ValueError('stale export receipt')
        if (arm/'table').exists() or (arm/'models').exists(): raise FileExistsError('refusing training overwrite')
        contract_path = next((arm/'features').glob('*/feature_contract.json'))
        table = P.build_table(arm/'features',contract_path,Path(inputs['labels']),Path(inputs['splits']),
            Path(inputs['label_policy']),arm/'table',registry_path=Path(inputs['registry']),
            go_obo=inputs['go_obo'],example_clusters=inputs['example_clusters'],locus_groups=inputs.get('locus_groups'),
            disjoint=tuple(original['split_audit']['disjoint_keys']))
        n = 0
        with Path(original['table']).open() as a, Path(table['table']).open() as b:
            for old,new in itertools.zip_longest(a,b):
                if old is None or new is None: raise ValueError('supervision row loss')
                old,new=json.loads(old),json.loads(new)
                for r in (old,new): r.pop('features'); r.pop('cohort_context_mode',None)
                if old != new: raise ValueError('supervision/splits/source/locus changed')
                n += 1
        contract = T.load_contract(contract_path); rows = T.read_table(Path(table['table']),contract)
        T.check_disjoint(rows,['group','genome_id'])
        splits = {s:[r for r in rows if r['split']==s] for s in T.SPLITS}
        (arm/'models').mkdir(); scores = {}
        for kind in ('logistic','gbdt'):
            previous=json.loads((baseline/'models'/f'model_{kind}.json').read_text())
            artifact=T.train(splits['train'],splits['calibration'],contract,kind,
                weighting=previous['training']['weighting'],folds=previous['training']['cv']['folds'],
                seed=previous['training']['cv']['seed'],threshold_rule=previous['display_threshold']['rule'],
                identity=T.split_identity(rows,['group','genome_id']))
            model=arm/'models'/f'model_{kind}.json'; T.write_frozen(artifact,model)
            scores[kind]=T.evaluate_frozen(model,contract,splits['test'],digest(model),digest(Path(table['table'])))
        write(arm/'evaluation.json',scores)
        write(arm/'training.json',dict(status='verified',supervision_rows_preserved=n,
            outputs=M.hashes(list((arm/'models').glob('*.json')) + [Path(table['table']),arm/'evaluation.json'])))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('stage',choices=['export','train'])
    for name in ('inputs','baseline','jar'): p.add_argument('--'+name,type=Path)
    p.add_argument('--out',type=Path,required=True); a=p.parse_args()
    if a.stage=='export': export(a.inputs.resolve(),a.baseline.resolve(),a.jar.resolve(),a.out.resolve())
    else: train(a.out.resolve())
