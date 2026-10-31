"""Audit historical Price calibration summaries against an unfitted zero baseline.

This rechecks report arithmetic, not raw inference or experimental negative truth.
No parameters or thresholds are fitted to these test reports.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path

VARIANTS = ('clean', 'dpz', 'enzbert', 'graphec', 'mapred', 'topec', 'ensemble_pool4')
ARMS = ('stage0', 'stage1', 'stage1_rlc', 'stage1_rlc_cohort',
        'learned_logistic', 'learned_gbdt')


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def count(value):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or \
            not math.isfinite(value) or value < 0 or int(value) != value:
        raise ValueError('invalid nonnegative count')
    return int(value)


def audit(root):
    results = []
    for variant in VARIANTS:
        path = root / f'price-{variant}-report.json'
        raw = json.loads(path.read_text())
        ec = raw['evaluation']['kinds']['ec']
        if raw['status'] != 'evaluated' or set(ec['arms']) != set(ARMS):
            raise ValueError('incomplete historical report')
        pairs = sum(count(g['pairs']) for g in ec['universe'].values())
        positive = sum(count(g['true_pairs']) for g in ec['universe'].values())
        proteins = sum(count(g['labelled_proteins']) for g in ec['universe'].values())
        if not 0 < positive <= pairs or proteins != 148 or positive != 151:
            raise ValueError('unexpected Price evaluation population')
        zero = positive / pairs
        arms = {}
        for arm in ARMS:
            result = ec['arms'][arm]
            calibration = result['calibration']
            bins = calibration['reliability']
            if count(calibration['n_pairs']) != pairs or \
                    sum(count(b['n']) for b in bins) != pairs or \
                    not math.isclose(sum(b['n'] * b['frac_true'] for b in bins),
                                     positive, rel_tol=0, abs_tol=1e-8):
                raise ValueError('universe and reliability-bin counts disagree')
            for metric in ('brier', 'ece'):
                value = calibration[metric]
                if isinstance(value, bool) or not math.isfinite(value) or not 0 <= value <= 1:
                    raise ValueError('invalid reported calibration metric')
            top1 = result['ec_top1']
            if count(top1['n_labelled']) != proteins:
                raise ValueError('top-1 population mismatch')
            arms[arm] = {
                'brier': calibration['brier'], 'ece': calibration['ece'],
                'brier_minus_zero': calibration['brier'] - zero,
                'ece_minus_zero': calibration['ece'] - zero,
                'brier_ratio_to_zero': calibration['brier'] / zero,
                'top1_accuracy': top1['top1_accuracy'],
            }
        results.append({'variant': variant, 'report': str(path), 'report_sha256': digest(path),
                        'genomes': len(ec['universe']), 'labelled_proteins': proteins,
                        'reference_positive_pairs': positive, 'evaluated_pairs': pairs,
                        'zero_baseline': {'brier': zero, 'ece': zero,
                                          'top1_accuracy': None,
                                          'reason': 'All scores are tied; no protein/function ranking is provided.'},
                        'arms': arms})
    return {
        'schema': 'gspa.price-report-calibration-baseline-audit/1',
        'status': 'report_arithmetic_verified', 'source_sha256': digest(Path(__file__)),
        'scope': 'All seven stored historical Price reports and all six arms; no inference rerun.',
        'baseline': 'Predict probability zero for every pair, with no fitting or test-selected threshold.',
        'derivation': 'For binary reference labels and score zero, both Brier and equal-width-bin ECE equal positive_pairs / evaluated_pairs.',
        'limitations': [
            'Unlabeled pairs use the reports\' closed-world reference convention; they are not proven experimental negatives.',
            'Report hashes and internal counts are verified; raw predictions, input provenance and biological accuracy are not revalidated here.',
            'These are historical models before current predictor-support gates; this does not measure current-model accuracy.',
            'Zero is a diagnostic comparator, not a recommended annotator or evidence of an absent function.',
            'Different predictors have different pair universes; their absolute Brier/ECE values are not a matched cross-predictor comparison.'
        ],
        'results': results,
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reports', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    report = audit(args.reports)
    with args.out.open('x') as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write('\n')
