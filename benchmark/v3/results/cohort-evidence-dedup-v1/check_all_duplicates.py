"""Verify duplicate comparisons over every original genome, without function labels."""
import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import prepare_cohort_controls as P
import run_cohort_perturbations as R
import verified_batch as V
from prepare_sequence_loss import digest


def check(inputs, root):
    preparation = P.verify(inputs)
    originals = sorted(preparation['genomes'])
    result = {}
    for kind in ('logistic', 'gbdt'):
        receipts, models = {}, set()
        for arm in ('intact', 'shuffled_neighborhoods', 'duplicated_strain_5x'):
            out = root / kind / arm / 'out'
            path = out / V.RECEIPT
            spec = json.loads(path.read_text())['spec']
            expected = set(originals)
            if arm == 'duplicated_strain_5x':
                expected |= set(preparation['duplicate_genomes'])
            if set(spec['outputs']) != expected or V.file_hashes(spec['inputs']) != spec['inputs'] or not V.reusable(out, spec):
                raise ValueError('invalid control execution: ' + kind + '/' + arm)
            model = Path(V.option_values(spec['command'])['--learned-model'][-1]).resolve()
            models.add(spec['inputs'][str(model)])
            receipts[arm] = digest(path)
        if len(models) != 1:
            raise ValueError('model differs between arms')
        comparison = R.compare_test_cohort_features(root / kind / 'intact/out',
                                                   root / kind / 'duplicated_strain_5x/out', originals)
        difference = comparison['max_abs_probability_diff']
        invariant = (difference is not None and difference <= 1e-12 and
                     all(comparison[k] == 0 for k in ('support_transitions', 'missing_before', 'missing_after')))
        result[kind] = dict(comparison=comparison, invariant_at_absolute_tolerance_1e_12=invariant,
                            receipts=receipts, model_sha256=next(iter(models)))
    return dict(original_genomes=originals, models=result,
                sources={str(Path(p).resolve()): digest(p) for p in (__file__, P.__file__, R.__file__, V.__file__)},
                scope='all original genomes in this exact-copy control; not heterogeneous-cluster or biological validation')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('inputs', 'root', 'out'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        parser.error('refusing output overwrite')
    args.out.write_text(json.dumps(check(args.inputs, args.root), indent=2, sort_keys=True) + '\n')
