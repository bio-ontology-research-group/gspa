"""Verified production cohort controls with intact-calibration-only thresholds.

Reuse the locus-preserving sequence-loss evaluator. This previously examined
cohort measures reference support and intervention sensitivity, not a blind test
or experimental false-assignment rate.
"""
import argparse
import json
from pathlib import Path

import prepare_cohort_controls as P
import run_cohort_perturbations as R
import run_sequence_loss_benchmark as E
from prepare_sequence_loss import digest, write_json

ARMS = ('intact', 'shuffled_neighborhoods', 'duplicated_strain_5x')


def run(inputs, work, model):
    preparation = P.verify(inputs)
    R.CACHE = Path(preparation['cache'])
    R.PERT_DIR = work
    genomes = sorted(preparation['genomes'])
    common = dict(learned_model_path=model, rlc_suggester=False,
                  sources=(Path(__file__), Path(P.__file__), Path(E.__file__), inputs / 'preparation.json'))
    for arm in ARMS:
        duplicate = arm == 'duplicated_strain_5x'
        R.run_perturbation_batch(arm, genomes=genomes + (preparation['duplicate_genomes'] if duplicate else []),
            layout_dir=inputs / 'layout_shuffled' if arm == 'shuffled_neighborhoods' else R.CACHE / 'layout',
            cohort_meta_path=inputs / 'metadata_duplicated.json' if duplicate else R.CACHE / 'cohort_inputs/cohort_metadata.json',
            cohort_og_path=inputs / 'orthogroups_duplicated.tsv' if duplicate else R.CACHE / 'cohort_inputs/cohort_orthogroups.tsv',
            **common)
    P.verify(inputs)


def verified(inputs, work, splits_path):
    preparation = P.verify(inputs)
    splits = E.read_splits(splits_path)
    if set(splits) != set(preparation['genomes']):
        raise ValueError('control/split genome coverage differs')
    predictions, receipts = E.verified_predictions(work, splits, arms=ARMS,
        extra_genomes={'duplicated_strain_5x': preparation['duplicate_genomes']})
    return preparation, splits, predictions, receipts


def bindings(inputs, labels, splits_path, go_path):
    return {str(p.resolve()): digest(p) for p in (inputs / 'preparation.json', labels, splits_path, go_path,
                                                Path(__file__), Path(E.__file__), Path(P.__file__))}


def calibrate(inputs, work, policy_path, labels, splits_path, go_path):
    if policy_path.exists():
        raise FileExistsError('policy is frozen; refusing overwrite')
    preparation, splits, predictions, receipts = verified(inputs, work, splits_path)
    calibration = sorted(g for g, split in splits.items() if split == 'calibration')
    if not calibration:
        raise ValueError('no calibration genomes')
    go = E.M.load_obo(str(go_path))
    E.propagated(predictions, go)
    truth = E.read_truth(labels, calibration, go)
    policy = E.select_policy(predictions['intact'], truth, calibration)
    policy.update(schema='gspa.cohort-controls-policy/1', calibration_genomes=calibration,
                  receipts=receipts, inputs=bindings(inputs, labels, splits_path, go_path))
    write_json(policy_path, policy)
    return policy


def evaluate(inputs, work, policy_path, labels, splits_path, go_path, n_boot=2000):
    preparation, splits, predictions, receipts = verified(inputs, work, splits_path)
    policy = json.loads(policy_path.read_text())
    calibration = sorted(g for g, split in splits.items() if split == 'calibration')
    if policy.get('schema') != 'gspa.cohort-controls-policy/1' or policy['calibration_genomes'] != calibration:
        raise ValueError('invalid frozen calibration policy')
    if receipts != policy['receipts'] or policy['inputs'] != bindings(inputs, labels, splits_path, go_path):
        raise ValueError('inputs or receipts changed after calibration')
    test = sorted(g for g, split in splits.items() if split == 'test')
    if not test or set(test) & set(calibration):
        raise ValueError('missing or overlapping test genomes')
    go = E.M.load_obo(str(go_path))
    E.propagated(predictions, go)
    truth = E.read_truth(labels, test, go)
    for arm in ARMS:
        for g in splits:
            allowed = set(preparation['genomes'][g]['retained_proteins'])
            for kind in ('GO', 'EC'):
                if any(p not in allowed for p, _ in predictions[arm][g][kind]):
                    raise ValueError('control output contains undeclared protein')
    report = E.summarize(predictions, truth, preparation, policy, test,
                         n_boot=n_boot, arms=ARMS, scopes=('all_original',))
    report.update(schema='gspa.cohort-controls-evaluation/1', receipts=receipts,
        policy_sha256=digest(policy_path), preparation_sha256=digest(inputs / 'preparation.json'),
        duplicated_strain_comparison=R.compare_test_cohort_features(
            work / 'intact/out', work / 'duplicated_strain_5x/out', test),
        limits=['previously examined development cohort, not a new blind test',
                'unmatched calls are unverified, not experimentally disproved functions',
                'layout shuffle preserves sequence/claims/contig membership; it is not sequence loss',
                'exact copies within a known strain cluster; not near-identical reassembled strains',
                'no metabolic gaps supplied; no reaction-context effect claimed',
                'GO ancestors are correlated annotation pairs, not independent observations'])
    return E.clean_json(report)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare', 'run', 'calibrate', 'evaluate'))
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--work', type=Path)
    parser.add_argument('--cache', type=Path, default=R.CACHE)
    parser.add_argument('--splits', type=Path, default=R.RES / 'splits.tsv')
    parser.add_argument('--labels', type=Path, default=R.RES / 'labels.tsv')
    parser.add_argument('--go', type=Path, default=R.CACHE / 'go-basic.obo')
    parser.add_argument('--source-genome')
    parser.add_argument('--model', type=Path, default=R.ROOT / 'benchmark/v3/results/model-input-support-v1/models/model_logistic.json')
    parser.add_argument('--policy', type=Path)
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    if args.action == 'prepare':
        if not args.source_genome:
            parser.error('--source-genome must name a training genome')
        P.prepare(args.cache, args.inputs, args.splits, args.source_genome)
    elif args.work is None:
        parser.error('--work is required')
    elif args.action == 'run':
        run(args.inputs, args.work, args.model)
    elif args.policy is None:
        parser.error('--policy is required')
    elif args.action == 'calibrate':
        print(json.dumps(calibrate(args.inputs, args.work, args.policy, args.labels, args.splits, args.go), indent=2))
    else:
        if args.out is None or args.out.exists():
            parser.error('--out must name a new report')
        write_json(args.out, evaluate(args.inputs, args.work, args.policy, args.labels, args.splits, args.go))


if __name__ == '__main__':
    main()
