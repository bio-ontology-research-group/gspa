#!/usr/bin/env python3
"""
GSPA Benchmark Step 4: Evaluate a theta.json against the multi-metric
objective described in plan §C.4.

For a given theta, this script:
  1. Ablates 10% of true GO annotations from the cached claims.jsonl.
  2. Injects 50 taxon-violating annotations into the (ablated) claims.
  3. Invokes `gspa-cli integrate --claims <file> --theta <file> --out <tsv>`
     with the reference ontology/metric handles.
  4. Computes F-max (overall + per aspect), essential recovery, gap-fill
     recovery, taxon-violation suppression, and the composite score.
  5. Writes a JSON summary to stdout (or --result-json file).

Called per BO iteration by 03_learn_integration.py. Can also be run
standalone to evaluate a single theta.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path

# Make the lib package importable regardless of CWD.
THIS_DIR = Path(__file__).resolve().parent
if str(THIS_DIR) not in sys.path:
    sys.path.insert(0, str(THIS_DIR))

from lib import gspa_eval  # noqa: E402
from lib import perturb  # noqa: E402


def run_cli_integrate(
    cli_binary: str,
    claims_file: str,
    theta_file: str,
    out_tsv: str,
    go_owl: str | None,
    essential_profile: str | None,
    essentials_file: str | None,
    ec2go: str | None,
    pathways: str | None,
    taxonomy: str | None,
    operons: str | None,
    gaps: str | None,
    enabled_priors: str | None,
    lite: bool,
) -> float:
    """Invoke gspa-cli integrate and return wall-clock seconds."""
    args = [
        cli_binary, 'integrate',
        '--claims', claims_file,
        '--theta', theta_file,
        '--out', out_tsv,
    ]
    if go_owl:            args += ['--go-owl', go_owl]
    if essential_profile: args += ['--essential-profile', essential_profile]
    if essentials_file:   args += ['--essential-functions', essentials_file]
    if ec2go:             args += ['--ec2go', ec2go]
    if pathways:          args += ['--pathways', pathways]
    if taxonomy:          args += ['--taxonomy', taxonomy]
    if operons:           args += ['--operons', operons]
    if gaps:              args += ['--gaps', gaps]
    if enabled_priors:    args += ['--enable-priors', enabled_priors]
    if lite:              args.append('--lite')

    start = time.time()
    result = subprocess.run(args, capture_output=True, text=True)
    elapsed = time.time() - start
    if result.returncode != 0:
        sys.stderr.write(result.stdout)
        sys.stderr.write(result.stderr)
        raise RuntimeError(f"gspa-cli integrate failed (rc={result.returncode})")
    return elapsed


def evaluate(
    theta_file: str,
    claims_file: str,
    truth_file: str,
    *,
    cli_binary: str,
    essentials_go_terms: set[str],
    go_owl: str | None = None,
    essential_profile: str | None = None,
    essentials_ref_file: str | None = None,
    ec2go: str | None = None,
    pathways: str | None = None,
    taxonomy: str | None = None,
    operons: str | None = None,
    gaps: str | None = None,
    enabled_priors: str | None = None,
    ablation_fraction: float = 0.10,
    taxon_n_injections: int = 50,
    lite: bool = False,
    seed: int = 42,
    weights: gspa_eval.ObjectiveWeights | None = None,
    tmp_dir: str | None = None,
) -> gspa_eval.EvaluationResult:
    weights = weights or gspa_eval.ObjectiveWeights()
    truth = gspa_eval.load_ground_truth(truth_file)

    workdir = Path(tmp_dir) if tmp_dir else Path(tempfile.mkdtemp(prefix='gspa-eval-'))
    workdir.mkdir(parents=True, exist_ok=True)

    # --- 1. Ablate claims for gap-fill test ---
    ablated_path = workdir / 'claims.ablated.jsonl'
    ablation = perturb.ablate_claims_for_gapfill(
        claims_path=claims_file,
        truth=truth,
        output_path=str(ablated_path),
        fraction=ablation_fraction,
        seed=seed,
    )

    # --- 2. Inject taxon violations into the ablated claims ---
    perturbed_path = workdir / 'claims.perturbed.jsonl'
    injection = perturb.inject_taxon_violations(
        claims_path=str(ablated_path),
        truth=truth,
        output_path=str(perturbed_path),
        n_injections=taxon_n_injections,
        seed=seed + 1,
    )

    # --- 3. Run gspa-cli integrate ---
    integrated_tsv = workdir / 'integrated.tsv'
    elapsed = run_cli_integrate(
        cli_binary=cli_binary,
        claims_file=str(perturbed_path),
        theta_file=theta_file,
        out_tsv=str(integrated_tsv),
        go_owl=go_owl,
        essential_profile=essential_profile,
        essentials_file=essentials_ref_file,
        ec2go=ec2go,
        pathways=pathways,
        taxonomy=taxonomy,
        operons=operons,
        gaps=gaps,
        enabled_priors=enabled_priors,
        lite=lite,
    )

    # --- 4. Load predictions and compute metrics ---
    predictions = gspa_eval.load_integrated_predictions(str(integrated_tsv))

    fmax_go, t_go = gspa_eval.compute_fmax(predictions, truth)
    fmax_mf, _   = gspa_eval.compute_fmax(predictions, truth, aspect_filter='MF')
    fmax_bp, _   = gspa_eval.compute_fmax(predictions, truth, aspect_filter='BP')
    fmax_cc, _   = gspa_eval.compute_fmax(predictions, truth, aspect_filter='CC')

    essential_recovery = gspa_eval.compute_essential_recovery(
        predictions, essentials_go_terms, threshold=0.5,
    )
    gapfill_recovery = gspa_eval.compute_gapfill_recovery(
        predictions, ablation.removed_pairs, threshold=0.5,
    )
    taxon_suppression = gspa_eval.compute_taxon_suppression(
        predictions, injection.injected_pairs, injection.input_scores,
    )

    result = gspa_eval.EvaluationResult(
        fmax_go=fmax_go,
        fmax_mf=fmax_mf,
        fmax_bp=fmax_bp,
        fmax_cc=fmax_cc,
        essential_recovery=essential_recovery,
        gapfill_recovery=gapfill_recovery,
        taxon_violation_suppression=taxon_suppression,
        runtime_seconds=elapsed,
        details={
            'ablation_n_removed': ablation.n_removed,
            'injection_n_injected': injection.n_injected,
            'fmax_threshold': t_go,
            'integrated_tsv': str(integrated_tsv),
        },
    )
    result.composite = gspa_eval.compute_composite(result, weights)
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--theta', required=True, help='theta.json with hyperparameters')
    p.add_argument('--claims', required=True, help='Cached claims.jsonl')
    p.add_argument('--truth', required=True, help='Ground-truth TSV (test_go_annotations.tsv)')
    p.add_argument('--essentials', required=True,
                   help='Essential functions TSV (for essential_recovery metric)')
    p.add_argument('--cli', default='gspa',
                   help='Path to the gspa CLI binary (default: gspa on PATH)')
    p.add_argument('--go-owl', default=None)
    p.add_argument('--essential-profile', default=None)
    p.add_argument('--essentials-ref', default=None,
                   help='Essentials file passed to gspa-cli integrate --essential-functions')
    p.add_argument('--ec2go', default=None)
    p.add_argument('--pathways', default=None)
    p.add_argument('--taxonomy', default=None)
    p.add_argument('--operons', default=None)
    p.add_argument('--gaps', default=None)
    p.add_argument('--enable-priors', default=None,
                   help='Comma-separated list of priors to enable in gspa-cli integrate')
    p.add_argument('--lite', action='store_true', help='Skip ELK (no process coherence)')
    p.add_argument('--ablation-fraction', type=float, default=0.10)
    p.add_argument('--taxon-n-injections', type=int, default=50)
    p.add_argument('--seed', type=int, default=42)
    p.add_argument('--tmp-dir', default=None)
    p.add_argument('--result-json', default=None,
                   help='Write EvaluationResult as JSON to this file (default: stdout)')
    args = p.parse_args()

    essentials_go_terms = gspa_eval.load_essential_functions(args.essentials)

    result = evaluate(
        theta_file=args.theta,
        claims_file=args.claims,
        truth_file=args.truth,
        cli_binary=args.cli,
        essentials_go_terms=essentials_go_terms,
        go_owl=args.go_owl,
        essential_profile=args.essential_profile,
        essentials_ref_file=args.essentials_ref,
        ec2go=args.ec2go,
        pathways=args.pathways,
        taxonomy=args.taxonomy,
        operons=args.operons,
        gaps=args.gaps,
        enabled_priors=args.enable_priors,
        ablation_fraction=args.ablation_fraction,
        taxon_n_injections=args.taxon_n_injections,
        lite=args.lite,
        seed=args.seed,
        tmp_dir=args.tmp_dir,
    )

    payload = result.as_dict()
    if args.result_json:
        os.makedirs(os.path.dirname(os.path.abspath(args.result_json)) or '.', exist_ok=True)
        with open(args.result_json, 'w') as fh:
            json.dump(payload, fh, indent=2)
    else:
        json.dump(payload, sys.stdout, indent=2)
        sys.stdout.write('\n')


if __name__ == '__main__':
    main()
