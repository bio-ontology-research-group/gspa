#!/usr/bin/env python3
"""
Pillar 1 & 2 Benchmark Grid Orchestrator.
Runs gspa-unified.jar integrate for 9 reference genomes across 9 ablation configurations:
  A1_noisy_or
  A2_essentiality
  A3_coherence
  A4_gap_filling
  A5_genomic_context
  A6_consistency_anatomical
  A7_stage1_all
  A8_stage2_remove
  A9_stage2_relax
"""

import argparse
import concurrent.futures
import os
import subprocess
import sys
import time
from pathlib import Path

JAVA = "/data/hohndor/envs/java21/bin/java"
JAR = "/data/hohndor/gspa/bin/gspa-unified.jar"
ROOT = "/data/hohndor/gspa"
REF = f"{ROOT}/reference"
BENCH9 = f"{ROOT}/proteomes/bench9"
OPERONS = f"{ROOT}/proteomes/operons"
WORK_DIR = f"{ROOT}/proteomes/bench9_unified"
OUT_DIR = f"{WORK_DIR}/integrated"
LOG_DIR = f"{WORK_DIR}/logs"

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


def build_command(genome, cfg):
    tag = genome["tag"]
    taxon = genome["taxon"]
    k = genome["k"]

    claims_file = f"{BENCH9}/with_iprscan/{tag}_claims_merged.jsonl"
    out_file = f"{OUT_DIR}/{tag}_{cfg}_integrated.tsv"
    log_file = f"{LOG_DIR}/{tag}_{cfg}.log"

    # Gaps arg
    gaps_file = f"{BENCH9}/gapseq/{tag}_gaps.jsonl"
    gaps_args = ["--gaps", gaps_file] if os.path.exists(gaps_file) and os.path.getsize(gaps_file) > 0 else []

    # Operons arg
    operons_file = f"{OPERONS}/{tag}_operons.tsv"
    operons_args = ["--operons", operons_file] if os.path.exists(operons_file) and os.path.getsize(operons_file) > 0 else []

    cmd = [
        JAVA, "-Xmx8g", "-jar", JAR, "integrate",
        "--claims", claims_file,
        "--go-owl", f"{REF}/go.owl",
        "--out", out_file,
        "--lite",
    ]

    if cfg == "A1_noisy_or":
        cmd += ["--enable-priors", "none"]
    elif cfg == "A2_essentiality":
        cmd += ["--enable-priors", "essentiality", "--essential-profile", k]
    elif cfg == "A3_coherence":
        cmd += [
            "--enable-priors", "coherence",
            "--pathways", f"{REF}/kegg_pathways.tsv",
            "--modules", f"{REF}/kegg_modules.tsv",
            "--ec2go", f"{REF}/ec2go.txt",
        ]
    elif cfg == "A4_gap_filling":
        cmd += [
            "--enable-priors", "gap_filling",
            "--pathways", f"{REF}/kegg_pathways.tsv",
            "--modules", f"{REF}/kegg_modules.tsv",
            "--ec2go", f"{REF}/ec2go.txt",
        ] + gaps_args
    elif cfg == "A5_genomic_context":
        cmd += [
            "--enable-priors", "genomic_context",
            "--pathways", f"{REF}/kegg_pathways.tsv",
            "--ec2go", f"{REF}/ec2go.txt",
        ] + operons_args
    elif cfg == "A6_consistency_anatomical":
        cmd += ["--enable-priors", "consistency,anatomical", "--taxon", taxon]
    elif cfg == "A7_stage1_all":
        cmd += [
            "--enable-priors", "essentiality,coherence,consistency,anatomical,gap_filling,genomic_context",
            "--essential-profile", k,
            "--taxon", taxon,
            "--pathways", f"{REF}/kegg_pathways.tsv",
            "--modules", f"{REF}/kegg_modules.tsv",
            "--ec2go", f"{REF}/ec2go.txt",
        ] + operons_args + gaps_args
    elif cfg == "A8_stage2_remove":
        actions_file = f"{OUT_DIR}/{tag}_A8_actions.tsv"
        cmd += [
            "--enable-priors", "essentiality,coherence,consistency,anatomical,gap_filling,genomic_context",
            "--essential-profile", k,
            "--taxon", taxon,
            "--pathways", f"{REF}/kegg_pathways.tsv",
            "--modules", f"{REF}/kegg_modules.tsv",
            "--ec2go", f"{REF}/ec2go.txt",
            "--enforce-stage2",
            "--enforce-consistency-mode", "remove",
            "--enforcement-actions-out", actions_file,
        ] + operons_args + gaps_args
    elif cfg == "A9_stage2_relax":
        actions_file = f"{OUT_DIR}/{tag}_A9_actions.tsv"
        cmd += [
            "--enable-priors", "essentiality,coherence,consistency,anatomical,gap_filling,genomic_context",
            "--essential-profile", k,
            "--taxon", taxon,
            "--pathways", f"{REF}/kegg_pathways.tsv",
            "--modules", f"{REF}/kegg_modules.tsv",
            "--ec2go", f"{REF}/ec2go.txt",
            "--enforce-stage2",
            "--enforce-consistency-mode", "relax",
            "--enforcement-actions-out", actions_file,
        ] + operons_args + gaps_args
    else:
        raise ValueError(f"Unknown config {cfg}")

    return cmd, out_file, log_file


def run_job(genome, cfg, force=False):
    tag = genome["tag"]
    cmd, out_file, log_file = build_command(genome, cfg)

    if not force and os.path.exists(out_file) and os.path.getsize(out_file) > 1000:
        return f"[SKIP] {tag} {cfg} (already exists: {os.path.getsize(out_file)} bytes)"

    t0 = time.time()
    with open(log_file, "w") as lf:
        proc = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT)
    dt = time.time() - t0

    if proc.returncode != 0:
        return f"[FAIL] {tag} {cfg} (rc={proc.returncode}, {dt:.1f}s) - see {log_file}"

    lines = 0
    if os.path.exists(out_file):
        with open(out_file) as f:
            lines = sum(1 for _ in f)

    return f"[DONE] {tag} {cfg} ({lines:,} annotations, {dt:.1f}s)"


def run_batch(tasks, jobs=16):
    manifest_path = f"{WORK_DIR}/ablation_batch_manifest.tsv"
    with open(manifest_path, "w") as mf:
        mf.write("genome_id\tclaims_path\tout_path\ttaxon\tessential_profile\toperons_path\tgaps_path\tmax_iter\tenforce_stage2\tenforce_mode\tenable_priors\n")
        for g, cfg in tasks:
            tag = g["tag"]
            taxon = g["taxon"]
            k = g["k"]
            claims_file = f"{BENCH9}/with_iprscan/{tag}_claims_merged.jsonl"
            out_file = f"{OUT_DIR}/{tag}_{cfg}_integrated.tsv"
            gaps_file = f"{BENCH9}/gapseq/{tag}_gaps.jsonl"
            has_gaps = os.path.exists(gaps_file) and os.path.getsize(gaps_file) > 0
            operons_file = f"{OPERONS}/{tag}_operons.tsv"
            has_ops = os.path.exists(operons_file) and os.path.getsize(operons_file) > 0

            row_tax = "-"
            row_ess = "-"
            row_ops = "-"
            row_gaps = "-"
            row_s2 = "false"
            row_mode = "-"
            row_priors = "none"

            if cfg == "A1_noisy_or":
                row_priors = "none"
            elif cfg == "A2_essentiality":
                row_priors = "essentiality"
                row_ess = k
            elif cfg == "A3_coherence":
                row_priors = "coherence"
            elif cfg == "A4_gap_filling":
                row_priors = "gap_filling"
                if has_gaps:
                    row_gaps = gaps_file
            elif cfg == "A5_genomic_context":
                row_priors = "genomic_context"
                if has_ops:
                    row_ops = operons_file
            elif cfg == "A6_consistency_anatomical":
                row_priors = "consistency,anatomical"
                row_tax = taxon
            elif cfg in ("A7_stage1_all", "A8_stage2_remove", "A9_stage2_relax"):
                row_priors = "essentiality,coherence,consistency,anatomical,gap_filling,genomic_context"
                row_tax = taxon
                row_ess = k
                if has_ops:
                    row_ops = operons_file
                if has_gaps:
                    row_gaps = gaps_file
                if cfg == "A8_stage2_remove":
                    row_s2 = "true"
                    row_mode = "remove"
                elif cfg == "A9_stage2_relax":
                    row_s2 = "true"
                    row_mode = "relax"

            mf.write("\t".join([
                f"{tag}_{cfg}", claims_file, out_file, row_tax, row_ess,
                row_ops, row_gaps, "-", row_s2, row_mode, row_priors
            ]) + "\n")

    cmd = [
        JAVA, "-Xmx48g", "-jar", JAR, "integrate",
        "--batch-manifest", manifest_path,
        "--batch-threads", str(jobs),
        "--enable-priors", "essentiality,coherence,consistency,anatomical,gap_filling,genomic_context",
        "--go-owl", f"{REF}/go.owl",
        "--pathways", f"{REF}/kegg_pathways.tsv",
        "--modules", f"{REF}/kegg_modules.tsv",
        "--ec2go", f"{REF}/ec2go.txt",
        "--lite",
    ]
    print("Running shared-memory batch integration for all tasks...", flush=True)
    proc = subprocess.run(cmd)
    if proc.returncode != 0:
        raise RuntimeError(f"Batch integration failed with exit code {proc.returncode}")


def main():
    parser = argparse.ArgumentParser(description="Run Pillar 1 & 2 ablation grid")
    parser.add_argument("--jobs", "-j", type=int, default=8, help="Parallel worker threads")
    parser.add_argument("--genome", "-g", default=None, help="Filter to specific genome tag")
    parser.add_argument("--config", "-c", default=None, help="Filter to specific config")
    parser.add_argument("--force", "-f", action="store_true", help="Force re-run existing files")
    args = parser.parse_args()

    os.makedirs(OUT_DIR, exist_ok=True)
    os.makedirs(LOG_DIR, exist_ok=True)

    genomes = [g for g in GENOMES if args.genome is None or g["tag"] == args.genome]
    configs = [c for c in CONFIGS if args.config is None or c == args.config]

    tasks = []
    for g in genomes:
        for c in configs:
            out_file = f"{OUT_DIR}/{g['tag']}_{c}_integrated.tsv"
            if not args.force and os.path.exists(out_file) and os.path.getsize(out_file) > 1000:
                continue
            tasks.append((g, c))

    print(f"Total jobs to run: {len(tasks)} (concurrency: {args.jobs})", flush=True)
    if not tasks:
        return

    t_start = time.time()
    run_batch(tasks, jobs=args.jobs)
    total_time = time.time() - t_start
    print(f"\nAll tasks completed in {total_time:.1f}s ({total_time/60:.2f}m)", flush=True)


if __name__ == "__main__":
    main()
