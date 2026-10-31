#!/usr/bin/env python3
"""Re-export real features, preserve supervision, and refit input-coverage manifests.

Uses previously examined cohort data. This is a provenance/compatibility replay,
not independent biological validation. Unsupported rows remain explicit abstentions.
"""
import argparse
from collections import Counter
import json
from pathlib import Path
import sys

import prepare_training_data as P
import run_training_cohort_ablations as C
from verified_batch import run_batch
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
import train_context_model as T

ROOT = Path(__file__).resolve().parents[2]


def save(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')


def replay(original_manifest, original_models, out, jar):
    out.mkdir(parents=True, exist_ok=False)
    old = json.loads(original_manifest.read_text()); inputs = old['inputs']
    paths = [original_manifest, Path(old['table']), Path(__file__), Path(P.__file__), Path(T.__file__), jar]
    for key in ('contract', 'labels', 'splits', 'label_policy', 'registry', 'go_obo', 'example_clusters', 'locus_groups'):
        if inputs.get(key):
            p = Path(inputs[key]); paths.append(p)
            if inputs.get(key + '_sha256') and T.sha256_file(p) != inputs[key + '_sha256']:
                raise ValueError('original input changed: ' + key)
    if T.sha256_file(old['table']) != old['table_sha256']:
        raise ValueError('original table changed')
    paths.extend(p for p in Path(inputs['feature_dir']).rglob('*') if p.is_file())
    paths.extend(original_models.glob('model_*.json'))
    identity = {str(p.resolve()): T.sha256_file(p) for p in paths}
    save(out / 'inputs.json', identity)
    manifest = out / 'manifest.tsv'
    inference = out / 'inference'; features = out / 'features'
    C.make_manifest(C.CACHE / 'baseline/claims/homology', inference, manifest)
    reference = ROOT / 'benchmark/whole_genome_cafa/data/reference'
    command = ['java', '-Xmx12g', '-jar', str(jar), 'integrate', '--batch-manifest', str(manifest),
        '--batch-threads', '4', '--go-owl', str(reference / 'go.owl'), '--ec2go', str(reference / 'ec2go.txt'),
        '--pathways', str(reference / 'kegg_pathways.tsv'), '--modules', str(reference / 'kegg_modules.tsv'),
        '--reasoner-cache', str(C.CACHE / 'reasoner-cache'),
        '--reaction-graph', str(ROOT.parent / 'gapsmith/data/seed_reactions.tsv'),
        '--diffusion-mets', str(ROOT.parent / 'gapsmith/data/diffusion_mets.tsv'),
        '--reaction-ec-aliases', str(ROOT.parent / 'gapsmith/data/seed_Enzyme_Class_Reactions_Aliases_unique.tsv'),
        '--cohort', '--cohort-metadata', str(C.CACHE / 'cohort_inputs/cohort_metadata.json'),
        '--cohort-orthogroups', str(C.CACHE / 'cohort_inputs/cohort_orthogroups.tsv'), '--cohort-holdout', 'CLADE',
        '--feature-export-dir', str(features), '--batch-stats-out', str(inference / 'batch_stats.tsv')]
    save(out / 'command.json', command)
    run_batch(command, manifest, inference, sources=[Path(__file__), Path(C.__file__)])
    export_rows = 0
    for previous in sorted(Path(inputs['feature_dir']).glob('*/candidate_features.jsonl')):
        current = features / previous.parent.name / previous.name
        if not current.exists(): raise ValueError('missing re-export: ' + str(current))
        with previous.open() as before, current.open() as after:
            import itertools
            for a, b in itertools.zip_longest(before, after):
                if a is None or b is None: raise ValueError('candidate universe changed')
                old_row, new_row = json.loads(a), json.loads(b)
                if any(old_row[k] != new_row[k] for k in ('candidate_key', 'features', 'genome_id')):
                    raise ValueError('re-export changed numerical features: ' + old_row['candidate_key'])
                if new_row['sequence_sources'] != ['swissprot_diamond_baseline']:
                    raise ValueError('unexpected production source profile')
                export_rows += 1
    table = P.build_table(features, Path(inputs['contract']), Path(inputs['labels']), Path(inputs['splits']),
        Path(inputs['label_policy']), out / 'table', registry_path=Path(inputs['registry']),
        go_obo=inputs['go_obo'], example_clusters=inputs['example_clusters'], locus_groups=inputs.get('locus_groups'),
        disjoint=tuple(old['split_audit']['disjoint_keys']))
    original_rows = [json.loads(s) for s in Path(old['table']).read_text().splitlines()]
    rows_json = [json.loads(s) for s in Path(table['table']).read_text().splitlines()]
    if original_rows != [{k:r[k] for k in old_r} for old_r,r in zip(original_rows, rows_json)] or len(original_rows) != len(rows_json):
        raise ValueError('supervision or features changed during source metadata replay')
    contract = T.load_contract(Path(inputs['contract'])); rows = T.read_table(Path(table['table']), contract)
    T.check_disjoint(rows, ['group', 'genome_id'])
    splits = {s:[r for r in rows if r['split'] == s] for s in T.SPLITS}
    models = out / 'models'; models.mkdir()
    reports = {}
    for kind in ('logistic', 'gbdt'):
        prior = json.loads((original_models / f'model_{kind}.json').read_text())
        art = T.train(splits['train'], splits['calibration'], contract, kind,
            weighting=prior['training']['weighting'], folds=prior['training']['cv']['folds'], seed=prior['training']['cv']['seed'],
            threshold_rule=prior['display_threshold']['rule'], identity=T.split_identity(rows, ['group', 'genome_id']))
        keys = ('kind', 'normalization', 'logistic', 'gbdt', 'calibration', 'display_threshold')
        if any(prior.get(k) != art.get(k) for k in keys): raise ValueError('fitted scoring parameters changed')
        path = models / f'model_{kind}.json'; digest = T.write_frozen(art, path)
        evaluation = T.evaluate_frozen(path, contract, splits['test'], digest, T.sha256_file(table['table']))
        save(models / f'test_report_{kind}.json', evaluation)
        reports[kind] = {'sha256': digest, 'input_support': art['input_support'], 'test_coverage': evaluation['input_coverage'],
                         'scoring_parameters_unchanged': True, 'locus_model': art['locus_model']}
    for p,digest in identity.items():
        if T.sha256_file(p) != digest: raise ValueError('input changed during replay: ' + p)
    result = {'schema':'gspa.model-input-support-replay/1', 'status':'verified', 'export_rows':export_rows,
              'training_rows':len(rows), 'models':reports, 'scope':'Previously examined data; coverage is not a generalization guarantee.'}
    save(out / 'verification.json', result)
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('original-manifest', 'original-models', 'out', 'jar'):
        p.add_argument('--' + name, type=Path, required=True)
    a = p.parse_args()
    r = replay(a.original_manifest.resolve(), a.original_models.resolve(), a.out.resolve(), a.jar.resolve())
    print(json.dumps({k:r[k] for k in ('status','export_rows','training_rows')}))
