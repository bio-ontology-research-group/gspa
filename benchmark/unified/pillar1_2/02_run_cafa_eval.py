#!/usr/bin/env python3
"""
Pillar 1 & 2 CAFA Evaluator.
Runs benchmark_pgap_v2.py for all 9 genomes across all 9 ablation configurations against
both exp and all truth sets.
"""

import argparse
import concurrent.futures
import json
import os
import subprocess
import sys
import time
from pathlib import Path

PYTHON = "/storage/miniforge3/envs/metagenomics/bin/python3"
BENCH_SCRIPT = "/data/hohndor/gspa/benchmark-py/benchmark_pgap_v2.py"
ROOT = "/data/hohndor/gspa"
TRUTH_DIR = f"{ROOT}/proteomes/truth_dual"
WORK_DIR = f"{ROOT}/proteomes/bench9_unified"
INTEGRATED_DIR = f"{WORK_DIR}/integrated"
CAFA_OUT_DIR = f"{WORK_DIR}/cafa_eval"

GENOMES = [
    "ecoli",
    "hpylori",
    "mgenitalium",
    "mjannaschii",
    "ecolo157",
    "bsubtilis",
    "mtb",
    "synechocystis",
    "paeruginosa",
]

CONFIGS = [
    "A1_noisy_or",
    "A2_essentiality",
    "A3_coherence",
    "A4_gap_filling",
    "A5_genomic_context",
    "A6_consistency_anatomical",
    "A7_stage1_all",
    "A8_stage2_remove",
    "A9_stage2_relax",
]


def run_eval(tag, cfg, n_bootstrap=50, force=False):
    int_file = f"{INTEGRATED_DIR}/{tag}_{cfg}_integrated.tsv"
    out_json = f"{CAFA_OUT_DIR}/{tag}_{cfg}_cafa.json"

    if not os.path.exists(int_file) or os.path.getsize(int_file) < 1000:
        return f"[MISSING] {tag} {cfg} (no integrated file: {int_file})"

    if not force and os.path.exists(out_json) and os.path.getsize(out_json) > 100:
        return f"[SKIP] {tag} {cfg} (already evaluated)"

    truth_exp = f"{TRUTH_DIR}/{tag}_truth_exp.tsv"
    truth_all = f"{TRUTH_DIR}/{tag}_truth_all.tsv"

    cmd = [
        PYTHON, BENCH_SCRIPT,
        "--truth", f"exp:{truth_exp}",
        "--truth", f"all:{truth_all}",
        "--gspa", int_file,
        "--tag", tag,
        "--n-bootstrap", str(n_bootstrap),
    ]

    t0 = time.time()
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    dt = time.time() - t0

    if proc.returncode != 0:
        return f"[FAIL] {tag} {cfg} (rc={proc.returncode}, {dt:.1f}s): {proc.stderr[:200]}"

    try:
        data = json.loads(proc.stdout)
        with open(out_json, "w") as f:
            json.dump(data, f, indent=2)
        fmax_exp = data["by_truth"]["exp"]["results"][0]["fmax_overall"]
        cafa_exp = data["by_truth"]["exp"]["results"][0]["fmax_cafa_overall"]
        return f"[DONE] {tag} {cfg} ({dt:.1f}s) -> fmax_micro(exp)={fmax_exp:.4f}, fmax_cafa(exp)={cafa_exp:.4f}"
    except Exception as e:
        return f"[ERR] {tag} {cfg} parse error: {e}"


def main():
    parser = argparse.ArgumentParser(description="Run CAFA evaluation on ablation grid")
    parser.add_argument("--jobs", "-j", type=int, default=8, help="Parallel worker threads")
    parser.add_argument("--genome", "-g", default=None, help="Filter to specific genome tag")
    parser.add_argument("--config", "-c", default=None, help="Filter to specific config")
    parser.add_argument("--n-bootstrap", type=int, default=50, help="Number of bootstrap iterations")
    parser.add_argument("--force", "-f", action="store_true", help="Force re-evaluation")
    args = parser.parse_args()

    os.makedirs(CAFA_OUT_DIR, exist_ok=True)

    genomes = [g for g in GENOMES if args.genome is None or g == args.genome]
    configs = [c for c in CONFIGS if args.config is None or c == args.config]

    tasks = [(g, c) for g in genomes for c in configs]
    print(f"Total CAFA evaluation jobs: {len(tasks)} (concurrency: {args.jobs})")

    t_start = time.time()
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = {executor.submit(run_eval, g, c, args.n_bootstrap, args.force): (g, c) for g, c in tasks}
        for future in concurrent.futures.as_completed(futures):
            tag, cfg = futures[future]
            try:
                res = future.result()
                print(res, flush=True)
            except Exception as e:
                print(f"[ERROR] {tag} {cfg}: {e}", flush=True)

    total_time = time.time() - t_start
    print(f"\nAll CAFA evaluation tasks completed in {total_time:.1f}s ({total_time/60:.2f}m)")


if __name__ == "__main__":
    main()
