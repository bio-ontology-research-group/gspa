#!/usr/bin/env python3
"""Plan step G1: the production anomaly extractor on the contract-bearing DIAMOND claims of the 22 Price genomes.

Calls ``roles_closure/run_price_anomalies.run`` unchanged (arms stage1, learned_logistic, learned_gbdt; its own
genus-blocked split, seed 20261101) with the predictor name ``swissprot_diamond_baseline``; the only reason for this
wrapper is that the original command line restricts ``--predictor`` to the three EC predictors it was written for.
Parts: ``<parts>/swissprot_diamond_baseline/<genome>/part.json`` with the claims path and sha256.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'roles_closure'))
sys.path.insert(0, str(HERE.parent))
import run_price_anomalies as A  # noqa: E402

PREDICTOR = 'swissprot_diamond_baseline'


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('out', 'parts', 'layouts', 'cohort-metadata', 'jar', 'reference', 'models'):
        p.add_argument('--' + name, required=True, type=Path)
    p.add_argument('--batch-threads', type=int, default=4)
    p.add_argument('--heap', default='40g')
    a = p.parse_args(argv)
    source = A.run(PREDICTOR, a.out, a.parts, a.layouts, a.cohort_metadata, a.jar, a.reference, a.models, a.batch_threads, a.heap)
    print(json.dumps({'predictor': PREDICTOR, 'arms': sorted(source['runs']), 'abstained': sorted(source['abstained_arms']),
                      'genomes': len(source['splits'])}))


if __name__ == '__main__':
    os.environ.setdefault('TMPDIR', str(Path.home() / 'tmp'))
    main()
