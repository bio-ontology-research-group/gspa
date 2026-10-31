"""Re-export and refit contract v3, preserving supervision and original splits.

Export, train and evaluation are separate immutable stages. This previously
examined cohort is not an independent biological validation set.
"""
import argparse
from collections import Counter
import itertools
import json
from pathlib import Path
import sys

import prepare_training_data as P
import run_training_cohort_ablations as C
import verified_batch as V
from prepare_sequence_loss import digest, write_json

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts'))
import train_context_model as T

REMOVED = {'cohort.available', 'cohort.log_n_genomes', 'cohort.frac_with_orthogroup',
           'cohort.frac_with_function_seq', 'cohort.frac_with_context', 'cohort.n_genomes_incomplete',
           'cohort.w.n_peer_genomes', 'cohort.w.bg_windows'}


def hashes(paths):
    return {str(Path(p).resolve()): digest(p) for p in paths}


def check_hashes(identity):
    if hashes(identity) != identity:
        raise ValueError('migration inputs or outputs changed')


def compare_vectors(previous, current):
    changes, rows = Counter(), 0
    with previous.open() as a, current.open() as b:
        for left, right in itertools.zip_longest(a, b):
            if left is None or right is None:
                raise ValueError('candidate universe changed')
            old, new = json.loads(left), json.loads(right)
            for key in ('candidate_key', 'genome_id', 'sequence_sources'):
                if old[key] != new[key]:
                    raise ValueError('candidate identity or source changed: ' + key)
            x, y = old['features'], new['features']
            if set(x) - set(y) != REMOVED or set(y) - set(x):
                raise ValueError('unexpected contract column changes')
            for key in y:
                if x[key] != y[key]:
                    if not key.startswith('cohort.w.'):
                        raise ValueError('non-cohort feature changed: ' + key)
                    changes[key] += 1
            rows += 1
    return rows, changes


def compare_supervision(previous, current):
    count = 0
    with previous.open() as a, current.open() as b:
        for left, right in itertools.zip_longest(a, b):
            if left is None or right is None:
                raise ValueError('supervised row coverage changed')
            old, new = json.loads(left), json.loads(right)
            old.pop('features'); new.pop('features')
            if old != new:
                raise ValueError('supervision, locus grouping, source identity or split changed')
            count += 1
    return count


def export(original_manifest, original_models, out, jar):
    out.mkdir(parents=True, exist_ok=False)
    old = json.loads(original_manifest.read_text()); inputs = old['inputs']
    paths = [original_manifest, Path(old['table']), jar, Path(__file__), Path(P.__file__),
             Path(T.__file__), Path(C.__file__), Path(V.__file__)]
    paths += [Path(m.__file__) for m in tuple(sys.modules.values())
              if getattr(m, '__file__', None) and Path(m.__file__).suffix == '.py'
              and Path(m.__file__).resolve().is_relative_to(ROOT)]
    for key in ('contract', 'labels', 'splits', 'label_policy', 'registry', 'go_obo', 'example_clusters', 'locus_groups'):
        if inputs.get(key):
            p = Path(inputs[key]); paths.append(p)
            if inputs.get(key + '_sha256') and digest(p) != inputs[key + '_sha256']:
                raise ValueError('original supervision input changed: ' + key)
    if digest(old['table']) != old['table_sha256']:
        raise ValueError('original training table changed')
    paths += [p for p in Path(inputs['feature_dir']).rglob('*') if p.is_file()]
    paths += list(original_models.glob('model_*.json'))
    identity = hashes(paths)
    write_json(out / 'inputs.json', identity)
    write_json(out / 'settings.json', dict(original_manifest=str(original_manifest), original_models=str(original_models)))
    manifest, inference, features = out / 'manifest.tsv', out / 'inference', out / 'features'
    C.make_manifest(C.CACHE / 'baseline/claims/homology', inference, manifest)
    ref = ROOT / 'benchmark/whole_genome_cafa/data/reference'
    command = ['java', '-Xmx12g', '-jar', str(jar), 'integrate', '--batch-manifest', str(manifest),
        '--batch-threads', '4', '--go-owl', str(ref / 'go.owl'), '--ec2go', str(ref / 'ec2go.txt'),
        '--pathways', str(ref / 'kegg_pathways.tsv'), '--modules', str(ref / 'kegg_modules.tsv'),
        '--reasoner-cache', str(C.CACHE / 'reasoner-cache'),
        '--reaction-graph', str(ROOT.parent / 'gapsmith/data/seed_reactions.tsv'),
        '--diffusion-mets', str(ROOT.parent / 'gapsmith/data/diffusion_mets.tsv'),
        '--reaction-ec-aliases', str(ROOT.parent / 'gapsmith/data/seed_Enzyme_Class_Reactions_Aliases_unique.tsv'),
        '--cohort', '--cohort-metadata', str(C.CACHE / 'cohort_inputs/cohort_metadata.json'),
        '--cohort-orthogroups', str(C.CACHE / 'cohort_inputs/cohort_orthogroups.tsv'), '--cohort-holdout', 'CLADE',
        '--feature-export-dir', str(features)]
    write_json(out / 'command.json', command)
    V.run_batch(command, manifest, inference, sources=[Path(__file__), Path(C.__file__)])
    rows, changes = 0, Counter()
    previous = sorted(Path(inputs['feature_dir']).glob('*/candidate_features.jsonl'))
    current = sorted(features.glob('*/candidate_features.jsonl'))
    if {p.parent.name for p in previous} != {p.parent.name for p in current}:
        raise ValueError('export genome coverage changed')
    for a in previous:
        n, changed = compare_vectors(a, features / a.parent.name / a.name)
        rows += n; changes.update(changed)
    contracts = [p.parent / 'feature_contract.json' for p in current]
    if len({digest(p) for p in contracts}) != 1:
        raise ValueError('mixed feature contracts')
    contract = T.load_contract(contracts[0])
    if contract['id'] != 'gspa.context-features/3':
        raise ValueError('expected version 3 export')
    check_hashes(identity)
    write_json(out / 'export.json', dict(status='verified', rows=rows, changed_columns=dict(changes),
        removed_columns=sorted(REMOVED), contract=str(contracts[0]),
        outputs=hashes([p for p in features.rglob('*') if p.is_file()] +
                      [inference / V.RECEIPT, out / 'settings.json', out / 'command.json'])))


def verified(out):
    check_hashes(json.loads((out / 'inputs.json').read_text()))
    report = json.loads((out / 'export.json').read_text())
    check_hashes(report['outputs'])
    receipt = json.loads((out / 'inference' / V.RECEIPT).read_text())
    if V.file_hashes(receipt['spec']['inputs']) != receipt['spec']['inputs'] or not V.reusable(out / 'inference', receipt['spec']):
        raise ValueError('stale production export receipt')
    return report, json.loads((out / 'settings.json').read_text())


def train(out):
    report, settings = verified(out)
    if (out / 'table').exists() or (out / 'models').exists():
        raise FileExistsError('training stage refuses overwrite')
    old = json.loads(Path(settings['original_manifest']).read_text()); inputs = old['inputs']
    table = P.build_table(out / 'features', Path(report['contract']), Path(inputs['labels']), Path(inputs['splits']),
        Path(inputs['label_policy']), out / 'table', registry_path=Path(inputs['registry']),
        go_obo=inputs['go_obo'], example_clusters=inputs['example_clusters'], locus_groups=inputs.get('locus_groups'),
        disjoint=tuple(old['split_audit']['disjoint_keys']))
    n = compare_supervision(Path(old['table']), Path(table['table']))
    contract = T.load_contract(Path(report['contract'])); rows = T.read_table(Path(table['table']), contract)
    T.check_disjoint(rows, ['group', 'genome_id'])
    splits = {s: [r for r in rows if r['split'] == s] for s in T.SPLITS}
    models = out / 'models'; models.mkdir()
    for kind in ('logistic', 'gbdt'):
        old_model = json.loads((Path(settings['original_models']) / f'model_{kind}.json').read_text())
        artifact = T.train(splits['train'], splits['calibration'], contract, kind,
            weighting=old_model['training']['weighting'], folds=old_model['training']['cv']['folds'],
            seed=old_model['training']['cv']['seed'], threshold_rule=old_model['display_threshold']['rule'],
            identity=T.split_identity(rows, ['group', 'genome_id']))
        T.write_frozen(artifact, models / f'model_{kind}.json')
    verified(out)
    write_json(out / 'training.json', dict(status='verified', supervision_rows_preserved=n,
        outputs=hashes(list(models.glob('model_*.json')) + [Path(table['table']), out / 'table/training_table_manifest.json']),
        limitation='same examined development cohort; no test-based model or policy tuning'))


def evaluate(out):
    report, _ = verified(out)
    training = json.loads((out / 'training.json').read_text()); check_hashes(training['outputs'])
    if (out / 'evaluation.json').exists():
        raise FileExistsError('evaluation refuses overwrite')
    contract = T.load_contract(Path(report['contract'])); table = out / 'table/training_table.jsonl'
    test = [r for r in T.read_table(table, contract) if r['split'] == 'test']
    results = {}
    for kind in ('logistic', 'gbdt'):
        model = out / 'models' / f'model_{kind}.json'
        results[kind] = T.evaluate_frozen(model, contract, test, digest(model), digest(table))
    check_hashes(training['outputs'])
    write_json(out / 'evaluation.json', results)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('export', 'train', 'evaluate'))
    parser.add_argument('--out', type=Path, required=True)
    for name in ('original-manifest', 'original-models', 'jar'):
        parser.add_argument('--' + name, type=Path)
    args = parser.parse_args()
    if args.action == 'export':
        if not all((args.original_manifest, args.original_models, args.jar)):
            parser.error('export requires original manifest/models and engine JAR')
        export(args.original_manifest.resolve(), args.original_models.resolve(), args.out.resolve(), args.jar.resolve())
    else:
        globals()[args.action](args.out.resolve())
