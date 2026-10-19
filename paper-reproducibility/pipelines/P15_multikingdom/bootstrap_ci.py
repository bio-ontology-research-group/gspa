#!/usr/bin/env python3
"""Paired bootstrap confidence intervals for the 88-proteome multi-kingdom evaluation.

Reads   results/multikingdom_per_genome_eval.tsv
Writes  results/multikingdom_bootstrap_ci.tsv

For every domain (and all 88 proteomes pooled), every GSPA configuration
(GSPA_Iter1, GSPA_Iter3_Relax) and every F-max column, the per-genome paired
difference d_g = F(GSPA, g) - F(Noisy-OR, g) is formed. The table reports the
mean difference over genomes and a 95% percentile bootstrap interval obtained by
resampling genomes with replacement (10,000 resamples, numpy default_rng with a
fixed seed). The number of genomes with d_g > 0, = 0 and < 0 is also reported.

Usage: python3 bootstrap_ci.py --repo <repo root> [--n-boot 10000] [--seed 20261019]
"""
from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import numpy as np

BASE_CFG = "Noisy-OR"
GSPA_CFGS = ["GSPA_Iter1", "GSPA_Iter3_Relax"]
METRICS = ["fmax_all", "fmax_lt30", "fmax_mf", "fmax_bp", "fmax_cc"]
DOMAINS = ["archaea", "bacteria", "eukaryote", "virus", "ALL"]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--repo", required=True)
    ap.add_argument("--n-boot", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=20261019)
    a = ap.parse_args()
    repo = Path(a.repo).resolve()
    src = repo / "results" / "multikingdom_per_genome_eval.tsv"
    out = repo / "results" / "multikingdom_bootstrap_ci.tsv"

    with open(src, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    by = defaultdict(dict)  # (genome) -> config -> row
    dom = {}
    for r in rows:
        by[r["genome_id"]][r["config"]] = r
        dom[r["genome_id"]] = r["domain"]

    res = []
    for d in DOMAINS:
        genomes = sorted(g for g in by if d == "ALL" or dom[g] == d)
        for cfg in GSPA_CFGS:
            for met in METRICS:
                diffs = np.array([float(by[g][cfg][met]) - float(by[g][BASE_CFG][met])
                                  for g in genomes if cfg in by[g] and BASE_CFG in by[g]])
                n = len(diffs)
                rng = np.random.default_rng(a.seed)
                idx = rng.integers(0, n, size=(a.n_boot, n))
                boots = diffs[idx].mean(axis=1)
                lo, hi = np.percentile(boots, [2.5, 97.5])
                res.append(dict(domain=d, config=cfg, baseline=BASE_CFG, metric=met, n_genomes=n,
                                mean_base=round(float(np.mean([float(by[g][BASE_CFG][met]) for g in genomes])), 4),
                                mean_gspa=round(float(np.mean([float(by[g][cfg][met]) for g in genomes])), 4),
                                mean_delta=round(float(diffs.mean()), 4),
                                ci95_lo=round(float(lo), 4), ci95_hi=round(float(hi), 4),
                                n_pos=int((diffs > 0).sum()), n_zero=int((diffs == 0).sum()),
                                n_neg=int((diffs < 0).sum()), n_boot=a.n_boot, seed=a.seed,
                                method="paired percentile bootstrap over genomes"))
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(res[0].keys()), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(res)
    print(f"wrote {out} ({len(res)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
