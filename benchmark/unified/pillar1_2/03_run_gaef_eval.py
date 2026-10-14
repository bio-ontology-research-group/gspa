#!/usr/bin/env python3
"""
Pillar 1 & 2 GAEF Quality Evaluator.
For each genome and ablation configuration:
1. Converts integrated annotations + FASTA into synthetic GFF3 + GAF using make_synth_gff_gaf.py
2. Runs gspa-unified.jar evaluate --lite --taxon <taxon> -k <k> ...
3. Extracts completeness, pathway_coherence, complex_coherence, consistent, violation_count, mean_ic, composite_score.
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
SYNTH_SCRIPT = "/data/hohndor/gspa/benchmark-py/make_synth_gff_gaf.py"
JAVA = "/data/hohndor/envs/java21/bin/java"
JAR = "/data/hohndor/gspa/bin/gspa-unified.jar"
ROOT = "/data/hohndor/gspa"
REF = f"{ROOT}/reference"
WORK_DIR = f"{ROOT}/proteomes/bench9_unified"
INTEGRATED_DIR = f"{WORK_DIR}/integrated"
GAEF_DIR = f"{WORK_DIR}/gaef_eval"

GENOMES = [
    {"tag": "ecoli", "taxon": "NCBITaxon_562", "k": "bacteria"},
    {"tag": "hpylori", "taxon": "NCBITaxon_210", "k": "bacteria"},
    {"tag": "mgenitalium", "taxon": "NCBITaxon_2097", "k": "bacteria"},
    {"tag": "mjannaschii", "taxon": "NCBITaxon_2190", "k": "archaea"},
    {"tag": "ecolo157", "taxon": "NCBITaxon_562", "k": "bacteria"},
    {"tag": "bsubtilis", "taxon": "NCBITaxon_1423", "k": "bacteria"},
    {"tag": "mtb", "taxon": "NCBITaxon_1773", "k": "bacteria"},
    {"tag": "synechocystis", "taxon": "NCBITaxon_1148", "k": "bacteria"},
    {"tag": "paeruginosa", "taxon": "NCBITaxon_287", "k": "bacteria"},
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


def run_single_gaef(genome, cfg, force=False):
    tag = genome["tag"]
    taxon = genome["taxon"]
    k = genome["k"]

    int_file = f"{INTEGRATED_DIR}/{tag}_{cfg}_integrated.tsv"
    if not os.path.exists(int_file) or os.path.getsize(int_file) < 1000:
        return f"[MISSING] {tag} {cfg} (no integrated file)"

    fasta_file = f"{ROOT}/proteomes/{tag}.faa"
    gff_file = f"{GAEF_DIR}/{tag}_{cfg}.gff"
    gaf_file = f"{GAEF_DIR}/{tag}_{cfg}.gaf"
    quality_json = f"{GAEF_DIR}/{tag}_{cfg}_quality.json"
    eval_log = f"{GAEF_DIR}/{tag}_{cfg}_eval.log"

    if not force and os.path.exists(quality_json) and os.path.getsize(quality_json) > 100:
        return f"[SKIP] {tag} {cfg} (already evaluated)"

    t0 = time.time()

    # Step 1: make synthetic GFF3 + GAF
    if force or not os.path.exists(gff_file) or not os.path.exists(gaf_file) or os.path.getsize(gaf_file) < 100:
        synth_cmd = [
            PYTHON, SYNTH_SCRIPT,
            "--fasta", fasta_file,
            "--tsv", int_file,
            "--threshold", "0.5",
            "--gff-out", gff_file,
            "--gaf-out", gaf_file,
        ]
        synth_res = subprocess.run(synth_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if synth_res.returncode != 0:
            return f"[FAIL-SYNTH] {tag} {cfg}: {synth_res.stderr[:200]}"

    # Step 2: run gspa evaluate
    eval_cmd = [
        JAVA, "-Xmx8g", "-jar", JAR, "evaluate",
        "--lite",
        "--go-owl", f"{REF}/go.owl",
        "-i", gff_file,
        "-a", gaf_file,
        "-k", k,
        "--taxon", taxon,
        "--pathways", f"{REF}/kegg_pathways.tsv",
        "--modules", f"{REF}/kegg_modules.tsv",
        "--ec2go", f"{REF}/ec2go.txt",
        "-o", quality_json,
    ]

    with open(eval_log, "w") as elf:
        eval_res = subprocess.run(eval_cmd, stdout=elf, stderr=subprocess.STDOUT)

    dt = time.time() - t0
    if eval_res.returncode != 0:
        return f"[FAIL-EVAL] {tag} {cfg} (rc={eval_res.returncode}, {dt:.1f}s) - see {eval_log}"

    try:
        with open(quality_json) as qf:
            d = json.load(qf)
        comp = d.get("completeness", {}).get("score", 0.0) if isinstance(d.get("completeness"), dict) else d.get("completeness", 0.0)
        coh = d.get("coherence", {}) if isinstance(d.get("coherence"), dict) else {}
        p_coh = coh.get("pathway_coherence", 0.0)
        c_coh = coh.get("complex_coherence", 0.0)
        cons = d.get("consistency", {}).get("consistent", True) if isinstance(d.get("consistency"), dict) else d.get("consistent", True)
        viols = d.get("consistency", {}).get("violation_count", 0) if isinstance(d.get("consistency"), dict) else 0
        composite = d.get("summary", {}).get("composite_score") or d.get("composite_score", 0.0)
        return f"[DONE] {tag} {cfg} ({dt:.1f}s) -> comp={comp:.3f}, pw_coh={p_coh:.3f}, cx_coh={c_coh:.3f}, cons={cons}(viols={viols}), comp_score={composite:.3f}"
    except Exception as e:
        return f"[ERR] {tag} {cfg} parse error: {e}"


def main():
    parser = argparse.ArgumentParser(description="Run GAEF quality evaluation on ablation grid")
    parser.add_argument("--jobs", "-j", type=int, default=8, help="Parallel worker threads")
    parser.add_argument("--genome", "-g", default=None, help="Filter to specific genome tag")
    parser.add_argument("--config", "-c", default=None, help="Filter to specific config")
    parser.add_argument("--force", "-f", action="store_true", help="Force re-evaluation")
    args = parser.parse_args()

    os.makedirs(GAEF_DIR, exist_ok=True)

    genomes = [g for g in GENOMES if args.genome is None or g["tag"] == args.genome]
    configs = [c for c in CONFIGS if args.config is None or c == args.config]

    tasks = [(g, c) for g in genomes for c in configs]
    print(f"Total GAEF evaluation jobs: {len(tasks)} (concurrency: {args.jobs})")

    t_start = time.time()
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = {executor.submit(run_single_gaef, g, c, args.force): (g["tag"], c) for g, c in tasks}
        for future in concurrent.futures.as_completed(futures):
            tag, cfg = futures[future]
            try:
                res = future.result()
                print(res, flush=True)
            except Exception as e:
                print(f"[ERROR] {tag} {cfg}: {e}", flush=True)

    total_time = time.time() - t_start
    print(f"\nAll GAEF evaluation tasks completed in {total_time:.1f}s ({total_time/60:.2f}m)")


if __name__ == "__main__":
    main()
