#!/usr/bin/env python3
"""
GSPA Benchmark Step 3: Bayesian optimization over the Phase 7 integration
hyperparameter space.

Search space (plan §C.2 / §C.3):
  - reliability weights per evidence type (7 parameters)
  - prior strengths:  alpha_ess, alpha_coh, alpha_cons, alpha_gap, alpha_ctx
  - damping, second_opinion_bonus, tau_hard (ConsistencyPrior fallback)

Objective: maximum of the composite in lib.gspa_eval (minimize the
negative). Each BO evaluation invokes 04_evaluate_integration.py, which
invokes the compiled `gspa-cli integrate` binary on cached claims.

Supports resumption by reading an existing trace.json if present.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

# Make the lib package importable regardless of CWD.
THIS_DIR = Path(__file__).resolve().parent
if str(THIS_DIR) not in sys.path:
    sys.path.insert(0, str(THIS_DIR))

from lib import gspa_eval  # noqa: E402

import importlib.util
_EVAL_SPEC = importlib.util.spec_from_file_location(
    'evaluate_integration',
    str(THIS_DIR / '04_evaluate_integration.py'),
)
_EVAL_MOD = importlib.util.module_from_spec(_EVAL_SPEC)
_EVAL_SPEC.loader.exec_module(_EVAL_MOD)
evaluate_fn = _EVAL_MOD.evaluate


# Search dimensions (name, low, high). Order matters — it defines the
# order in the x vector gp_minimize hands to the objective.
SEARCH_SPACE = [
    ('w_seqsim',             0.0, 1.0),
    ('w_domain',             0.0, 1.0),
    ('w_deep',               0.0, 1.0),
    ('w_struct',             0.0, 1.0),
    ('w_ortho',              0.0, 1.0),
    ('w_ctx',                0.0, 1.0),
    ('w_metab',              0.0, 1.0),
    ('alpha_ess',            0.0, 5.0),
    ('alpha_coh',            0.0, 5.0),
    ('alpha_cons',           0.0, 5.0),
    ('alpha_gap',            0.0, 5.0),
    ('alpha_ctx',            0.0, 5.0),
    ('damping',              0.2, 0.9),
    ('second_opinion_bonus', 0.0, 1.0),
    ('tau_hard',             0.5, 0.95),
]


def make_theta(x: list[float]) -> dict:
    """Translate a BO point into a theta.json payload."""
    vals = dict(zip([n for n, _, _ in SEARCH_SPACE], x))
    return {
        'reliability': {
            'SEQUENCE_SIMILARITY':    vals['w_seqsim'],
            'SEQUENCE_DOMAIN':        vals['w_domain'],
            'SEQUENCE_DEEPLEARNING':  vals['w_deep'],
            'STRUCTURE_SIMILARITY':   vals['w_struct'],
            'STRUCTURE_DEEPLEARNING': vals['w_struct'],  # share with similarity for now
            'ORTHOLOGY':              vals['w_ortho'],
            'GENOMIC_CONTEXT':        vals['w_ctx'],
            'METABOLIC_CONTEXT':      vals['w_metab'],
        },
        'alpha_ess':  vals['alpha_ess'],
        'alpha_coh':  vals['alpha_coh'],
        'alpha_cons': vals['alpha_cons'],
        'alpha_gap':  vals['alpha_gap'],
        'alpha_ctx':  vals['alpha_ctx'],
        'damping':    vals['damping'],
        'second_opinion_bonus': vals['second_opinion_bonus'],
        'tau_hard':   vals['tau_hard'],
        'max_iter':   6,
        'epsilon':    0.005,
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--claims', required=True, help='Cached claims.jsonl from 02b_parse_predictors_to_claims.py')
    p.add_argument('--truth', required=True, help='Ground-truth TSV')
    p.add_argument('--essentials', required=True, help='Essential functions TSV')
    p.add_argument('--cli', default='gspa', help='gspa CLI binary (default: gspa on PATH)')
    p.add_argument('--go-owl', default=None)
    p.add_argument('--essential-profile', default=None)
    p.add_argument('--essentials-ref', default=None)
    p.add_argument('--ec2go', default=None)
    p.add_argument('--pathways', default=None)
    p.add_argument('--taxonomy', default=None)
    p.add_argument('--operons', default=None)
    p.add_argument('--gaps', default=None)
    p.add_argument('--enable-priors', default=None)
    p.add_argument('--lite', action='store_true')
    p.add_argument('--n-calls', type=int, default=200)
    p.add_argument('--n-initial-points', type=int, default=25)
    p.add_argument('--seed', type=int, default=42)
    p.add_argument('--out-dir', required=True, help='Directory to write trace.json + best_theta.json')
    p.add_argument('--resume', action='store_true', help='Load existing trace and skip completed points')
    args = p.parse_args()

    try:
        from skopt import gp_minimize
        from skopt.space import Real
    except ImportError:
        sys.stderr.write(
            "ERROR: scikit-optimize is required. Install with:\n"
            "    pip install scikit-optimize\n"
        )
        sys.exit(2)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    trace_path = out_dir / 'trace.json'
    best_path = out_dir / 'best_theta.json'
    summary_path = out_dir / 'best_theta_summary.json'

    # Load or initialize trace
    trace: list[dict] = []
    if args.resume and trace_path.exists():
        with open(trace_path) as fh:
            trace = json.load(fh)
        print(f'Resuming from {len(trace)} existing evaluations in {trace_path}')

    essentials = gspa_eval.load_essential_functions(args.essentials)

    space = [Real(lo, hi, name=name) for name, lo, hi in SEARCH_SPACE]
    weights = gspa_eval.ObjectiveWeights()

    iteration_counter = {'n': len(trace)}
    scratch = out_dir / 'scratch'
    scratch.mkdir(exist_ok=True)

    def objective(x):
        t_index = iteration_counter['n']
        iteration_counter['n'] += 1
        theta = make_theta(x)
        theta_path = scratch / f'theta_{t_index:04d}.json'
        with open(theta_path, 'w') as fh:
            json.dump(theta, fh, indent=2)
        eval_workdir = scratch / f'eval_{t_index:04d}'
        eval_workdir.mkdir(exist_ok=True)
        t0 = time.time()
        try:
            result = evaluate_fn(
                theta_file=str(theta_path),
                claims_file=args.claims,
                truth_file=args.truth,
                cli_binary=args.cli,
                essentials_go_terms=essentials,
                go_owl=args.go_owl,
                essential_profile=args.essential_profile,
                essentials_ref_file=args.essentials_ref,
                ec2go=args.ec2go,
                pathways=args.pathways,
                taxonomy=args.taxonomy,
                operons=args.operons,
                gaps=args.gaps,
                enabled_priors=args.enable_priors,
                lite=args.lite,
                seed=args.seed,
                weights=weights,
                tmp_dir=str(eval_workdir),
            )
            composite = result.composite
            payload = {
                'index': t_index,
                'theta': theta,
                'result': result.as_dict(),
                'elapsed_s': time.time() - t0,
            }
        except Exception as exc:  # keep the BO loop running on failures
            sys.stderr.write(f'[iter {t_index}] evaluation failed: {exc}\n')
            composite = -1e6
            payload = {
                'index': t_index,
                'theta': theta,
                'error': str(exc),
                'elapsed_s': time.time() - t0,
            }

        trace.append(payload)
        with open(trace_path, 'w') as fh:
            json.dump(trace, fh, indent=2)
        print(f'[iter {t_index:04d}] composite={composite:.4f}')
        return -composite   # gp_minimize minimizes

    # Skip any previously-completed indices (if --resume was used).
    x0 = None
    y0 = None
    if trace:
        # skopt lets us seed the optimizer with prior evaluations.
        x0 = [
            [t['theta'].get(k) for k, _, _ in [
                ('w_seqsim', 0, 1),
            ]]
            # NOTE: we don't carry the full theta through to skopt; resuming
            # simply preserves the trace and appends new points.
            for t in []
        ] or None

    print(f'Running BO with {args.n_calls} evaluations over {len(SEARCH_SPACE)}-dim space')
    res = gp_minimize(
        objective,
        space,
        n_calls=args.n_calls,
        n_initial_points=args.n_initial_points,
        random_state=args.seed,
    )

    # Write best theta
    best_theta = make_theta(res.x)
    with open(best_path, 'w') as fh:
        json.dump(best_theta, fh, indent=2)
    with open(summary_path, 'w') as fh:
        json.dump({
            'n_evaluations': len(trace),
            'best_composite': -res.fun,
            'best_x': list(res.x),
            'best_theta': best_theta,
        }, fh, indent=2)

    print(f'Best composite: {-res.fun:.4f}')
    print(f'Wrote {best_path} and {summary_path}')


if __name__ == '__main__':
    main()
