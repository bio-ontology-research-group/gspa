#!/usr/bin/env python3
"""Rebuild a prior training table under the repaired truth protocol and refit its models.

This is a regression replay on previously used data, not new biological validation.
No group declarations, negatives, labels or sources are invented for the replay.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import sys

import prepare_training_data as P

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts'))
import train_context_model as T


def sha(p):
    return T.sha256_file(p)


def save(p, value):
    p.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')


def replay(original_manifest, original_models, out):
    out.mkdir(parents=True, exist_ok=False)
    old = json.loads(original_manifest.read_text())
    inputs = old['inputs']
    original_table = Path(old['table'])
    if sha(original_table) != old['table_sha256']:
        raise ValueError('original training table differs from its manifest')
    for key in ('labels', 'splits', 'label_policy', 'go_obo', 'example_clusters', 'locus_groups'):
        if inputs.get(key):
            if inputs.get(key + '_sha256') != sha(inputs[key]):
                raise ValueError('original input hash mismatch: ' + key)
    files = [original_manifest, original_table, Path(__file__).resolve(), Path(P.__file__), Path(T.__file__)]
    files.extend(Path(inputs[k]) for k in ('contract', 'labels', 'splits', 'label_policy', 'registry', 'go_obo', 'example_clusters', 'locus_groups') if inputs.get(k))
    files.extend(p for p in Path(inputs['feature_dir']).rglob('*') if p.is_file())
    files.extend(original_models.glob('model_*.json'))
    identity = {str(p.resolve()): sha(p) for p in files}
    save(out / 'input-sha256.json', identity)
    rebuilt = P.build_table(Path(inputs['feature_dir']), Path(inputs['contract']), Path(inputs['labels']),
        Path(inputs['splits']), Path(inputs['label_policy']), out / 'table', registry_path=Path(inputs['registry']),
        go_obo=inputs.get('go_obo'), example_clusters=inputs.get('example_clusters'), locus_groups=inputs.get('locus_groups'),
        disjoint=tuple(old['split_audit']['disjoint_keys']))
    before = [json.loads(s) for s in original_table.read_text().splitlines() if s]
    after = [json.loads(s) for s in Path(rebuilt['table']).read_text().splitlines() if s]
    # New provenance fields may differ; every original supervision field must survive.
    if len(before) != len(after) or before != [{k: r[k] for k in old_row} for old_row, r in zip(before, after)]:
        raise ValueError('original candidate labels/features or split membership changed')
    contract = T.load_contract(Path(inputs['contract']))
    rows = T.read_table(Path(rebuilt['table']), contract)
    T.check_disjoint(rows, ['group', 'genome_id'])
    splits = {s: [r for r in rows if r['split'] == s] for s in T.SPLITS}
    models = out / 'models'
    models.mkdir()
    reports = {}
    for kind in ('logistic', 'gbdt'):
        previous = json.loads((original_models / f'model_{kind}.json').read_text())
        artifact = T.train(splits['train'], splits['calibration'], contract, kind,
            weighting=previous['training']['weighting'], folds=previous['training']['cv']['folds'],
            seed=previous['training']['cv']['seed'], threshold_rule=previous['display_threshold']['rule'],
            identity=T.split_identity(rows, ['group', 'genome_id']), intended_use='production', synthetic=False)
        path = models / f'model_{kind}.json'
        digest = T.write_frozen(artifact, path)
        # Read test outcomes only after the fitted artifact is frozen.
        evaluation = T.evaluate_frozen(path, contract, splits['test'], digest, sha(rebuilt['table']))
        save(models / f'test_report_{kind}.json', evaluation)
        score_keys = ['kind', 'normalization', 'logistic', 'gbdt', 'calibration', 'display_threshold', 'locus_model']
        unchanged = all(previous.get(k) == artifact.get(k) for k in score_keys)
        if not unchanged:
            raise ValueError(f'{kind}: fitted scoring parameters changed in a truth-only protocol replay')
        reports[kind] = {'model_sha256': digest, 'scoring_parameters_unchanged': unchanged,
                         'locus_model': artifact['locus_model'], 'truth': artifact['training']['truth']}
    for path, digest in identity.items():
        if sha(path) != digest:
            raise ValueError('input changed during replay: ' + path)
    summary = {'schema': 'gspa.locus-truth-replay/1', 'status': 'verified',
               'scope': 'previously used data; engineering regression, not new biological validation',
               'rows': len(rows), 'original_supervision_unchanged': True,
               'label_basis_counts': dict(Counter(r.get('label_basis', 'unspecified') for r in rows)),
               'locus_audit': rebuilt['locus_audit'], 'models': reports}
    save(out / 'summary.json', summary)
    return summary


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--original-manifest', type=Path, required=True)
    ap.add_argument('--original-models', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    a = ap.parse_args()
    print(json.dumps(replay(a.original_manifest.resolve(), a.original_models.resolve(), a.out.resolve()), indent=2))
