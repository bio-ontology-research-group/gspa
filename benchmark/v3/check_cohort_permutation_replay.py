"""Compare frozen engine replays without reading function labels or refitting policy."""
import argparse
import csv
import json
import math
from pathlib import Path

import prepare_cohort_controls as P
import run_cohort_controls as C
import run_cohort_perturbations as R
import verified_batch as V
from prepare_sequence_loss import digest, write_json


def outputs(work, preparation):
    result = {}
    for arm in C.ARMS:
        folder = work / arm / 'out'
        path = folder / V.RECEIPT
        receipt = json.loads(path.read_text())
        spec = receipt['spec']
        expected = set(preparation['genomes'])
        if arm == 'duplicated_strain_5x':
            expected |= set(preparation['duplicate_genomes'])
        if set(spec['outputs']) != expected or V.file_hashes(spec['inputs']) != spec['inputs'] or not V.reusable(folder, spec):
            raise ValueError('invalid replay receipt: ' + arm)
        result[arm] = digest(path)
    return result


def rows(path):
    result = {}
    with path.open() as stream:
        for r in csv.DictReader(stream, delimiter='\t'):
            key = (r['protein_id'], r['type'], r['function_id'])
            p = float(r['posterior_prob'])
            if key in result or not math.isfinite(p) or not 0 <= p <= 1:
                raise ValueError('invalid prediction identity or score')
            result[key] = p
    return result


def compare(inputs, before, after):
    preparation = P.verify(inputs)
    original = sorted(preparation['genomes'])
    receipts = {name: outputs(folder, preparation) for name, folder in [('before', before), ('after', after)]}
    arms = {}
    for arm in C.ARMS:
        report = dict(compared=0, changed=0, increased=0, decreased=0, missing_before=0,
                      missing_after=0, max_absolute_change=0.)
        for gid in original:
            a = rows(before / arm / 'out' / f'{gid}.integrated.tsv')
            b = rows(after / arm / 'out' / f'{gid}.integrated.tsv')
            report['missing_before'] += len(b.keys() - a.keys())
            report['missing_after'] += len(a.keys() - b.keys())
            for key in a.keys() & b.keys():
                delta = b[key] - a[key]
                report['compared'] += 1
                report['max_absolute_change'] = max(report['max_absolute_change'], abs(delta))
                if abs(delta) > 1e-12:
                    report['changed'] += 1
                    report['increased' if delta > 0 else 'decreased'] += 1
        arms[arm] = report
    return dict(schema='gspa.cohort-permutation-replay/1', receipts=receipts,
                preparation_sha256=digest(inputs / 'preparation.json'), source_sha256=digest(__file__),
                original_genomes=original, prediction_changes=arms,
                duplicate_comparison={name: R.compare_test_cohort_features(folder / 'intact/out',
                    folder / 'duplicated_strain_5x/out', original) for name, folder in [('before', before), ('after', after)]},
                limits=['no function labels read; no threshold refitting or biological accuracy claim',
                        'unchanged fitted model; its calibration is not certified for changed cohort features',
                        'raw peer/background counts still reach the model; full duplication invariance remains open'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('inputs', 'before', 'after', 'out'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        parser.error('output must be a new file')
    write_json(args.out, compare(args.inputs, args.before, args.after))
