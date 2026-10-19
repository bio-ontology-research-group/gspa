#!/usr/bin/env python3
"""Exact McNemar tests on Hit@1 for the E. coli MG1655 strict 256-case LRO benchmark.

Reads the per-case ranking tables in benchmark/manifests/reference_data/lro_data/results_*.tsv
(same files and same 1-indexed `rank_of_p` convention as run_pillar2_all.py
`run_15_genome_stratified_lro`; rank 0 / -1 / blank = not in candidate pool = miss) and writes
results/pillar2_lro_mcnemar.tsv.

For a pair (A, B): b = cases where A has Hit@1 and B does not, c = the reverse;
p = two-sided exact binomial test b ~ Bin(b + c, 0.5) (scipy.stats.binomtest).

Usage:  python3 benchmark/unified/pillar2_fba_lro_price149/lro_stats.py
"""
from __future__ import annotations

import csv
from pathlib import Path

from scipy.stats import binomtest

ROOT = Path(__file__).resolve().parent.parent.parent.parent
LRO_DIR = ROOT / "benchmark" / "manifests" / "reference_data" / "lro_data"
RESULTS_DIR = ROOT / "results"

# Table 5 label -> ranking file (mapping as in run_pillar2_all.py, methods_256)
METHODS_256 = {
    "No target evidence (results_protein.tsv)": "results_protein.tsv",
    "Function prior only (results_function.tsv)": "results_function.tsv",
    "RLGC kernel density alone (results_rlgc_mg1655.tsv)": "results_rlgc_mg1655.tsv",
    "RLGC + coherence, Stage 1 M2 (results_m2_integrated_mg1655.tsv)": "results_m2_integrated_mg1655.tsv",
    "GSPA Track A strict masked (results_m4_trackA_strict_mg1655.tsv)": "results_m4_trackA_strict_mg1655.tsv",
    "ESM-2 centroids only (results_plmonly_mg1655.tsv)": "results_plmonly_mg1655.tsv",
    "GSPA full Track B (results_m4_trackB_mg1655.tsv)": "results_m4_trackB_mg1655.tsv",
}

PAIRS = [
    # the three comparisons requested for Table 5
    ("RLGC + coherence, Stage 1 M2 (results_m2_integrated_mg1655.tsv)", "No target evidence (results_protein.tsv)"),
    ("GSPA full Track B (results_m4_trackB_mg1655.tsv)", "ESM-2 centroids only (results_plmonly_mg1655.tsv)"),
    ("GSPA full Track B (results_m4_trackB_mg1655.tsv)", "RLGC + coherence, Stage 1 M2 (results_m2_integrated_mg1655.tsv)"),
    # supplementary
    ("RLGC kernel density alone (results_rlgc_mg1655.tsv)", "No target evidence (results_protein.tsv)"),
    ("RLGC + coherence, Stage 1 M2 (results_m2_integrated_mg1655.tsv)", "RLGC kernel density alone (results_rlgc_mg1655.tsv)"),
    ("GSPA Track A strict masked (results_m4_trackA_strict_mg1655.tsv)", "RLGC + coherence, Stage 1 M2 (results_m2_integrated_mg1655.tsv)"),
    ("GSPA full Track B (results_m4_trackB_mg1655.tsv)", "GSPA Track A strict masked (results_m4_trackA_strict_mg1655.tsv)"),
]


def load_1idx_ranks(path: Path) -> dict[int, int]:
    """Identical to run_pillar2_all.load_1idx_ranks."""
    out: dict[int, int] = {}
    if not path.exists():
        return out
    with open(path, "r", encoding="utf-8") as f:
        for idx, r in enumerate(csv.DictReader(f, delimiter="\t")):
            cidx = int(r["case_idx"]) if "case_idx" in r and r["case_idx"] != "" else idx
            rk_s = r.get("rank_of_p", "0")
            out[cidx] = int(rk_s) if rk_s not in ("", "-1", "None") else 0
    return out


def exact_mcnemar(hits_a: list[bool], hits_b: list[bool]) -> tuple[int, int, float]:
    b = sum(1 for x, y in zip(hits_a, hits_b) if x and not y)
    c = sum(1 for x, y in zip(hits_a, hits_b) if y and not x)
    p = binomtest(b, b + c, 0.5, alternative="two-sided").pvalue if b + c else 1.0
    return b, c, float(p)


def main() -> int:
    ranks = {m: load_1idx_ranks(LRO_DIR / fn) for m, fn in METHODS_256.items()}
    missing = [METHODS_256[m] for m, r in ranks.items() if not r]
    if missing:
        raise SystemExit(f"missing ranking tables: {missing}")
    keys = sorted(ranks["GSPA full Track B (results_m4_trackB_mg1655.tsv)"].keys())  # same key set as run_pillar2_all
    assert len(keys) == 256, len(keys)
    hits = {m: [r.get(k, 0) == 1 for k in keys] for m, r in ranks.items()}

    out = RESULTS_DIR / "pillar2_lro_mcnemar.tsv"
    with open(out, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["benchmark", "method_a", "method_b", "n", "hit1_a", "hit1_b", "b_a_only", "c_b_only",
                    "p_exact_two_sided", "test"])
        for a, bm in PAIRS:
            b, c, p = exact_mcnemar(hits[a], hits[bm])
            w.writerow(["Ecoli_MG1655_Strict_256", a, bm, len(keys), sum(hits[a]), sum(hits[bm]), b, c,
                        f"{p:.6g}", "exact McNemar (scipy.stats.binomtest, p=0.5); uncovered = miss"])
            print(f"{a}\n   vs {bm}: hit1 {sum(hits[a])} vs {sum(hits[bm])}, b={b} c={c} p={p:.4g}")
    print(f"Wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
