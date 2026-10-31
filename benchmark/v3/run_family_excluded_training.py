"""Freeze, export, fit and evaluate models for the verified family-excluded transfer.

This previously examined cohort measures reference support, not untouched or
experimentally calibrated biological accuracy. Training never receives test rows.
"""
import argparse
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
import sys

import prepare_training_data as P
import verified_batch as V
from migrate_cohort_feature_contract import hashes, check_hashes
from run_context_modes import replace

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts'))
import train_context_model as T

MODE = 'family_excluded'
SOURCE = 'swissprot_diamond_baseline'
KINDS = ('logistic', 'gbdt')


def save(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')


def training_settings(model):
    """Copy selection rules, never fitted coefficients or test-derived cutoffs."""
    return dict(weighting=model['training']['weighting'], folds=model['training']['cv']['folds'],
                seed=model['training']['cv']['seed'], threshold_rule=model['display_threshold']['rule'])


def verify_protocol(out):
    registration = json.loads((out / 'registration.json').read_text())
    check_hashes(registration['inputs'])
    if T.sha256_file(out / 'protocol.json') != registration['protocol_sha256']:
        raise ValueError('registered protocol changed')
    protocol = json.loads((out / 'protocol.json').read_text())
    if protocol['numpy_version'] != T.np.__version__ or protocol['python_version'] != sys.version:
        raise ValueError('registered training runtime changed')
    if protocol['mode'] != MODE or protocol['model_kinds'] != list(KINDS):
        raise ValueError('unexpected registered training scope')
    return protocol


def prepare(template, claims_replay, jar, out):
    receipt = json.loads((claims_replay / 'run.json').read_text())
    if receipt.get('status') != 'verified' or receipt.get('enriched_claims_reversible') is not True:
        raise ValueError('verified reversible predictor reconstruction required')
    for filename, key in [('inputs.json', 'inputs_sha256'), ('contracts.json', 'contract_manifest_sha256')]:
        if T.sha256_file(claims_replay / filename) != receipt[key]:
            raise ValueError('predictor receipt changed')
    check_hashes(json.loads((claims_replay / 'inputs.json').read_text()))
    contracts = json.loads((claims_replay / 'contracts.json').read_text())
    contract = contracts[MODE]['identity']
    if contract == contracts['homology']['identity']:
        raise ValueError('family-excluded and homology identities must differ')
    training = json.loads((template / 'training.json').read_text())
    if training['status'] != 'verified':
        raise ValueError('unverified template training')
    check_hashes(training['outputs'])
    original = json.loads((template / 'table/training_table_manifest.json').read_text())
    inputs = original['inputs']
    paths = [template/'training.json', template/'table/training_table_manifest.json',
             template/'manifest.tsv', template/'command.json', claims_replay/'run.json',
             claims_replay/'contracts.json', claims_replay/'inputs.json', Path(sys.executable)]
    supervision = {}
    for key in ('labels', 'splits', 'label_policy', 'registry', 'go_obo', 'example_clusters', 'locus_groups'):
        if inputs.get(key):
            p = Path(inputs[key]).resolve()
            if inputs.get(key + '_sha256') and T.sha256_file(p) != inputs[key + '_sha256']:
                raise ValueError('template supervision changed: ' + key)
            supervision[key] = str(p); paths.append(p)
    settings = {}
    for kind in KINDS:
        p = template/'models'/f'model_{kind}.json'; paths.append(p)
        settings[kind] = training_settings(json.loads(p.read_text()))
    with (template / 'manifest.tsv').open() as f:
        reader = csv.DictReader(f, delimiter='\t'); fields = reader.fieldnames; rows = list(reader)
    enriched = receipt['output_claims'][MODE]
    if {r['genome_id'] for r in rows} != set(enriched):
        raise ValueError('complete cohort genome coverage required')
    for row in rows:
        claim = enriched[row['genome_id']]
        if T.sha256_file(claim['path']) != claim['sha256']:
            raise ValueError('verified claims changed')
    # No source or output from the previous run is overwritten.
    out.mkdir(parents=True, exist_ok=False)
    arm = out/'default'; arm.mkdir(); inference = arm/'inference'; inference.mkdir()
    for row in rows:
        row['claims_path'] = enriched[row['genome_id']]['path']
        row['out_path'] = str(inference/(row['genome_id']+'.integrated.tsv'))
    manifest = arm/'manifest.tsv'
    with manifest.open('w') as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)
    command = json.loads((template/'command.json').read_text())
    if '--learned-model' in command:
        raise ValueError('feature export must not use a learned model')
    for flag, value in [('-jar', jar), ('--batch-manifest', manifest),
                        ('--feature-export-dir', arm/'features'), ('--batch-stats-out', inference/'batch_stats.tsv')]:
        replace(command, flag, value)
    save(arm/'command.json', command)
    modules = [Path(m.__file__).resolve() for m in tuple(sys.modules.values())
               if getattr(m, '__file__', None) and Path(m.__file__).suffix == '.py'
               and Path(m.__file__).resolve().is_relative_to(ROOT)]
    spec = V.describe_run(command, manifest, sources=modules)
    bound = hashes(paths + modules + [arm/'command.json']); bound.update(spec['inputs'])
    protocol = dict(schema='gspa.family-excluded-training/1', mode=MODE, model_kinds=list(KINDS),
        expected_contract=[dict(source=SOURCE, contract=contract)], claim_files=enriched,
        supervision=supervision, disjoint_keys=original['split_audit']['disjoint_keys'], settings=settings,
        evaluation='Frozen marginal test evaluation and full-table JVM/Python parity; no model promotion or test-based selection.',
        limitations='Previously examined reference-overlapping cohort; closed-world negatives are unconfirmed. No biological calibration, independent test, or trained null claim.',
        numpy_version=T.np.__version__, python_version=sys.version)
    save(out/'protocol.json', protocol)
    save(out/'registration.json', dict(status='registered', registered_at=datetime.now(timezone.utc).isoformat(),
        protocol_sha256=T.sha256_file(out/'protocol.json'), inputs=bound))
    verify_protocol(out)
    return protocol


def audit_export(features, claim_files, expected_contract):
    """Require complete claim identities and the correct predictor contract, not just row counts."""
    if {entry['source'] for entry in expected_contract} == {'clean', SOURCE}:
        from cohort_sequence_inputs import audit_export as audit_multisource_export
        return audit_multisource_export(features, claim_files, expected_contract)
    files = sorted(features.glob('*/candidate_features.jsonl'))
    if {p.parent.name for p in files} != set(claim_files):
        raise ValueError('feature genome coverage differs')
    counts = {}
    for file in files:
        expected = set()
        with Path(claim_files[file.parent.name]['path']).open() as f:
            for line in f:
                r = json.loads(line)
                expected.add((r['protein_id'], r['function_type'], r['function_id']))
        seen = set()
        with file.open() as f:
            for line in f:
                row = json.loads(line)
                key = (row['protein_id'], row['function_type'], row['function_id'])
                if row['genome_id'] != file.parent.name or key in seen:
                    raise ValueError('duplicate candidate or wrong genome')
                if row['sequence_contracts'] != expected_contract:
                    raise ValueError('wrong predictor identity in export')
                seen.add(key)
        if seen != expected:
            raise ValueError('candidate/claim identities differ')
        counts[file.parent.name] = len(seen)
    return counts


def export(out):
    protocol = verify_protocol(out); arm = out/'default'
    if (arm/'export.json').exists():
        raise FileExistsError('export already completed')
    command = json.loads((arm/'command.json').read_text())
    V.run_batch(command, arm/'manifest.tsv', arm/'inference', sources=[Path(__file__)])
    counts = audit_export(arm/'features', protocol['claim_files'], protocol['expected_contract'])
    verify_protocol(out)
    save(arm/'export.json', dict(status='verified', counts=counts,
        outputs=hashes([p for p in (arm/'features').rglob('*') if p.is_file()] + [arm/'inference'/V.RECEIPT])))
    return counts


def verified_export(out):
    protocol = verify_protocol(out); arm = out/'default'
    report = json.loads((arm/'export.json').read_text()); check_hashes(report['outputs'])
    record = json.loads((arm/'inference'/V.RECEIPT).read_text())
    check_hashes(record['spec']['inputs'])
    if not V.reusable(arm/'inference', record['spec']):
        raise ValueError('export no longer verifies')
    return protocol, arm


def train(out):
    protocol, arm = verified_export(out)
    if (arm/'table').exists() or (arm/'models').exists():
        raise FileExistsError('training refuses overwrite')
    inputs = protocol['supervision']
    contract_path = next((arm/'features').glob('*/feature_contract.json'))
    table = P.build_table(arm/'features', contract_path, Path(inputs['labels']), Path(inputs['splits']),
        Path(inputs['label_policy']), arm/'table', registry_path=Path(inputs['registry']),
        go_obo=inputs['go_obo'], example_clusters=inputs['example_clusters'], locus_groups=inputs.get('locus_groups'),
        disjoint=tuple(protocol['disjoint_keys']))
    contract = T.load_contract(contract_path); rows = T.read_table(Path(table['table']), contract)
    T.check_disjoint(rows, ['group', 'genome_id'])
    splits = {s: [r for r in rows if r['split'] == s] for s in T.SPLITS}
    (arm/'models').mkdir()
    for kind in KINDS:
        artifact = T.train(splits['train'], splits['calibration'], contract, kind,
            **protocol['settings'][kind], identity=T.split_identity(rows, ['group', 'genome_id']))
        T.write_frozen(artifact, arm/'models'/f'model_{kind}.json')
    verified_export(out)
    save(arm/'training.json', dict(status='verified', frozen_at=datetime.now(timezone.utc).isoformat(),
        rows_by_split={s:len(r) for s,r in splits.items()},
        outputs=hashes(list((arm/'models').glob('*.json')) + [Path(table['table']),arm/'table/training_table_manifest.json'])))


def evaluate(out):
    _, arm = verified_export(out)
    if (arm/'evaluation.json').exists():
        raise FileExistsError('evaluation refuses overwrite')
    training = json.loads((arm/'training.json').read_text()); check_hashes(training['outputs'])
    contract_path = next((arm/'features').glob('*/feature_contract.json')); contract = T.load_contract(contract_path)
    table = arm/'table/training_table.jsonl'
    test = [r for r in T.read_table(table, contract) if r['split'] == 'test']
    results = {}
    for kind in KINDS:
        model = arm/'models'/f'model_{kind}.json'
        results[kind] = T.evaluate_frozen(model, contract, test, T.sha256_file(model), T.sha256_file(table))
    check_hashes(training['outputs']); verify_protocol(out)
    save(arm/'evaluation.json', results)
    return results


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('stage', choices=('prepare','export','train','evaluate'))
    p.add_argument('--out', type=Path, required=True)
    for key in ('template','claims-replay','jar'):
        p.add_argument('--'+key, type=Path)
    a = p.parse_args()
    if a.stage == 'prepare':
        if not all((a.template,a.claims_replay,a.jar)):p.error('prepare requires template, claims-replay and jar')
        r = prepare(a.template.resolve(), a.claims_replay.resolve(), a.jar.resolve(), a.out.resolve())
    else:r = globals()[a.stage](a.out.resolve())
    print(json.dumps(r, indent=2))
