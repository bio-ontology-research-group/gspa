#!/usr/bin/env python3
"""
Pillar 1 & 2 Cross-Clade Taxon & Anatomical Perturbation Experiment.
For all 9 reference genomes:
1. Injects 100 cross-clade false-positive GO claims (50 anatomical, 50 clade-specific taxon-constrained).
2. Runs A1_noisy_or, A6_consistency_anatomical, A8_stage2_remove, A9_stage2_relax.
3. Evaluates:
   - Surviving violations at posterior >= 0.5
   - Suppression rate (%)
   - Impact on fmax_micro and fmax_cafa (via benchmark_pgap_v2.py against exp & all truth).
"""

import argparse
import concurrent.futures
import csv
import json
import os
import random
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path

PYTHON = "/storage/miniforge3/envs/metagenomics/bin/python3"
JAVA = "/data/hohndor/envs/java21/bin/java"
JAR = "/data/hohndor/gspa/bin/gspa-unified.jar"
BENCH_SCRIPT = "/data/hohndor/gspa/benchmark-py/benchmark_pgap_v2.py"
ROOT = "/data/hohndor/gspa"
REF = f"{ROOT}/reference"
BENCH9 = f"{ROOT}/proteomes/bench9"
OPERONS = f"{ROOT}/proteomes/operons"
TRUTH_DIR = f"{ROOT}/proteomes/truth_dual"
WORK_DIR = f"{ROOT}/proteomes/bench9_unified"
PERTURB_DIR = f"{WORK_DIR}/perturbation"
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

PERTURB_CONFIGS = [
    "A1_noisy_or",
    "A6_consistency_anatomical",
    "A7_stage1_all",
    "A8_stage2_remove",
    "A9_stage2_relax",
]


def load_aspect_map():
    aspect_map = {}
    path = f"{REF}/go_aspect_map.tsv"
    if os.path.exists(path):
        with open(path) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    aspect_map[parts[0]] = parts[1]
    return aspect_map


def pick_perturbation_terms():
    """Pick 50 anatomical and 50 taxon-constrained terms."""
    aspect_map = load_aspect_map()

    # 1. Anatomical terms
    # Mandatory terms from prompt:
    mandatory_anat = [
        "GO:0007507",  # heart development
        "GO:0009908",  # flower development
        "GO:0001755",  # neural crest cell migration
        "GO:0030282",  # bone mineralization
        "GO:0048437",  # floral organ development
    ]

    anat_candidates = []
    anat_tsv = f"{WORK_DIR}/go-anatomical-constraints.tsv"
    if not os.path.exists(anat_tsv):
        anat_tsv = "/home/leechuck/Public/software/gspa/gspa-core/src/main/resources/anatomical-constraints/go-anatomical-constraints.tsv"

    with open(anat_tsv) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            go_id = r["GO_ID"]
            anat_id = r["Anatomy_ID"]
            # Metazoan, plant, or specialized cellular structures
            if anat_id.startswith("UBERON:") or anat_id.startswith("PO:") or anat_id.startswith("CL:"):
                if go_id not in mandatory_anat:
                    anat_candidates.append(go_id)

    # Deterministic sample
    rng = random.Random(42)
    unique_anat = sorted(list(set(anat_candidates)))
    rng.shuffle(unique_anat)
    selected_anat = mandatory_anat + unique_anat[:45]

    # 2. Clade-specific taxon-constrained terms
    taxon_candidates = []
    tax_tsv = f"{WORK_DIR}/go-taxon-constraints.tsv"
    if not os.path.exists(tax_tsv):
        tax_tsv = "/home/leechuck/Public/software/gspa/gspa-core/src/main/resources/taxon-constraints/go-taxon-constraints.tsv"

    # Terms only in Eukaryota, Metazoa, Viridiplantae, Fungi, Chordata, Mammalia
    euk_taxa = {
        "NCBITaxon_2759", "NCBITaxon_33208", "NCBITaxon_33090",
        "NCBITaxon_4751", "NCBITaxon_7711", "NCBITaxon_40674",
        "NCBITaxon_89593", "NCBITaxon_9606"
    }
    with open(tax_tsv) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            go_id = r["GO_ID"]
            rel = r["Constraint_Type"]
            tax = r["Taxon_ID"]
            if rel == "only_in_taxon" and tax in euk_taxa:
                if go_id not in selected_anat:
                    taxon_candidates.append(go_id)

    unique_tax = sorted(list(set(taxon_candidates)))
    rng.shuffle(unique_tax)
    selected_tax = unique_tax[:50]

    all_terms = []
    for t in selected_anat:
        all_terms.append({"go_id": t, "category": "anatomical", "aspect": aspect_map.get(t, "BP")})
    for t in selected_tax:
        all_terms.append({"go_id": t, "category": "taxon_constraint", "aspect": aspect_map.get(t, "BP")})

    return all_terms


def inject_claims_for_genome(genome, terms, force=False):
    tag = genome["tag"]
    claims_file = f"{BENCH9}/with_iprscan/{tag}_claims_merged.jsonl"
    out_claims_file = f"{PERTURB_DIR}/{tag}_claims_perturbed.jsonl"
    injected_info_file = f"{PERTURB_DIR}/{tag}_injected_pairs.tsv"

    if not force and os.path.exists(injected_info_file) and os.path.exists(out_claims_file) and os.path.getsize(injected_info_file) > 100:
        records = []
        with open(injected_info_file) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for r in reader:
                records.append({
                    "protein_id": r["protein_id"],
                    "go_id": r["go_id"],
                    "aspect": r["aspect"],
                    "category": r["category"],
                    "raw_score": float(r["raw_score"]),
                })
        if len(records) == len(terms):
            return out_claims_file, records

    # Collect proteins
    proteins = []
    with open(claims_file) as f:
        for line in f:
            try:
                c = json.loads(line)
                proteins.append(c["protein_id"])
            except Exception:
                pass
    unique_proteins = sorted(list(set(proteins)))

    import hashlib
    seed = int(hashlib.md5(tag.encode("utf-8")).hexdigest(), 16) % (2**32)
    rng = random.Random(seed)
    # Pick 100 proteins (one per term)
    chosen_proteins = rng.sample(unique_proteins, len(terms)) if len(unique_proteins) >= len(terms) else [rng.choice(unique_proteins) for _ in range(len(terms))]

    injected_claims = []
    injected_records = []

    for idx_t, (term_info, pid) in enumerate(zip(terms, chosen_proteins)):
        raw_score = round(rng.uniform(0.62, 0.90), 4)
        src_primary = "diamond" if idx_t % 2 == 0 else "deepgo"
        c = {
            "protein_id": pid,
            "function_type": "GO",
            "function_id": term_info["go_id"],
            "go_aspect": term_info["aspect"],
            "source": src_primary,
            "raw_score": raw_score,
            "metadata": {"injected_perturbation": True, "category": term_info["category"]},
        }
        injected_claims.append(c)
        # For 60% of cross-clade false positives, simulate a corroborating domain/PLM hit
        if idx_t % 5 in (0, 1, 2):
            sec_score = round(max(0.52, raw_score - 0.12), 4)
            injected_claims.append({
                "protein_id": pid,
                "function_type": "GO",
                "function_id": term_info["go_id"],
                "go_aspect": term_info["aspect"],
                "source": "pfam" if src_primary == "diamond" else "diamond",
                "raw_score": sec_score,
                "metadata": {"injected_perturbation": True, "category": term_info["category"]},
            })
        injected_records.append({
            "protein_id": pid,
            "go_id": term_info["go_id"],
            "aspect": term_info["aspect"],
            "category": term_info["category"],
            "raw_score": raw_score,
        })

    # Write injected info
    with open(injected_info_file, "w") as f:
        writer = csv.DictWriter(f, delimiter="\t", fieldnames=["protein_id", "go_id", "aspect", "category", "raw_score"])
        writer.writeheader()
        writer.writerows(injected_records)

    # Write combined claims
    with open(out_claims_file, "w") as out_f:
        with open(claims_file) as in_f:
            for line in in_f:
                out_f.write(line)
        for c in injected_claims:
            out_f.write(json.dumps(c) + "\n")

    return out_claims_file, injected_records


def build_perturb_command(genome, cfg, claims_file):
    tag = genome["tag"]
    taxon = genome["taxon"]
    k = genome["k"]

    out_file = f"{PERTURB_DIR}/{tag}_{cfg}_perturbed_integrated.tsv"
    log_file = f"{LOG_DIR}/{tag}_{cfg}_perturbed.log"

    gaps_file = f"{BENCH9}/gapseq/{tag}_gaps.jsonl"
    gaps_args = ["--gaps", gaps_file] if os.path.exists(gaps_file) and os.path.getsize(gaps_file) > 0 else []

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
        actions_file = f"{PERTURB_DIR}/{tag}_A8_perturbed_actions.tsv"
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
        actions_file = f"{PERTURB_DIR}/{tag}_A9_perturbed_actions.tsv"
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

    return cmd, out_file, log_file


def run_perturb_job(genome, cfg, claims_file, force=False):
    tag = genome["tag"]
    cmd, out_file, log_file = build_perturb_command(genome, cfg, claims_file)

    if not force and os.path.exists(out_file) and os.path.getsize(out_file) > 1000:
        return f"[SKIP-INT] {tag} {cfg} perturbed (already exists)"

    t0 = time.time()
    with open(log_file, "w") as lf:
        proc = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT)
    dt = time.time() - t0

    if proc.returncode != 0:
        return f"[FAIL-INT] {tag} {cfg} perturbed (rc={proc.returncode}, {dt:.1f}s) - see {log_file}"

    lines = 0
    if os.path.exists(out_file):
        with open(out_file) as f:
            lines = sum(1 for _ in f)

    return f"[DONE-INT] {tag} {cfg} perturbed ({lines:,} annotations, {dt:.1f}s)"


def evaluate_perturbation(genome, cfg, injected_records):
    tag = genome["tag"]
    int_file = f"{PERTURB_DIR}/{tag}_{cfg}_perturbed_integrated.tsv"
    eval_json = f"{PERTURB_DIR}/{tag}_{cfg}_cafa_eval.json"

    # 1. Surviving injected claims
    injected_set = {(r["protein_id"], r["go_id"]) for r in injected_records}
    surviving = set()
    surviving_anat = 0
    surviving_tax = 0

    anat_set = {r["go_id"] for r in injected_records if r["category"] == "anatomical"}
    tax_set = {r["go_id"] for r in injected_records if r["category"] == "taxon_constraint"}

    if os.path.exists(int_file):
        with open(int_file) as f:
            headers = next(f, "").rstrip("\n").split("\t")
            if "protein_id" in headers and "function_id" in headers and "posterior_prob" in headers:
                idx_pid = headers.index("protein_id")
                idx_func = headers.index("function_id")
                idx_prob = headers.index("posterior_prob")
                for line in f:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) > idx_prob:
                        try:
                            prob = float(parts[idx_prob])
                            if prob >= 0.5:
                                key = (parts[idx_pid], parts[idx_func])
                                if key in injected_set:
                                    surviving.add(key)
                                    if parts[idx_func] in anat_set:
                                        surviving_anat += 1
                                    if parts[idx_func] in tax_set:
                                        surviving_tax += 1
                        except ValueError:
                            pass

    n_injected = len(injected_records)
    n_surviving = len(surviving)
    suppression_rate = (1.0 - n_surviving / n_injected) * 100.0

    # 2. CAFA evaluation using benchmark_pgap_v2.py
    truth_exp = f"{TRUTH_DIR}/{tag}_truth_exp.tsv"
    truth_all = f"{TRUTH_DIR}/{tag}_truth_all.tsv"

    fmax_micro_exp = 0.0
    fmax_cafa_exp = 0.0
    fmax_micro_all = 0.0
    fmax_cafa_all = 0.0

    cmd = [
        PYTHON, BENCH_SCRIPT,
        "--truth", f"exp:{truth_exp}",
        "--truth", f"all:{truth_all}",
        "--gspa", int_file,
        "--tag", tag,
        "--n-bootstrap", "50",
    ]
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if proc.returncode == 0:
        try:
            cafa_data = json.loads(proc.stdout)
            with open(eval_json, "w") as f:
                json.dump(cafa_data, f, indent=2)
            fmax_micro_exp = cafa_data["by_truth"]["exp"]["results"][0]["fmax_overall"]
            fmax_cafa_exp = cafa_data["by_truth"]["exp"]["results"][0]["fmax_cafa_overall"]
            fmax_micro_all = cafa_data["by_truth"]["all"]["results"][0]["fmax_overall"]
            fmax_cafa_all = cafa_data["by_truth"]["all"]["results"][0]["fmax_cafa_overall"]
        except Exception:
            pass

    return {
        "tag": tag,
        "config": cfg,
        "n_injected": n_injected,
        "n_surviving": n_surviving,
        "surviving_anat": surviving_anat,
        "surviving_tax": surviving_tax,
        "suppression_rate_pct": suppression_rate,
        "fmax_micro_exp": fmax_micro_exp,
        "fmax_cafa_exp": fmax_cafa_exp,
        "fmax_micro_all": fmax_micro_all,
        "fmax_cafa_all": fmax_cafa_all,
    }


def main():
    parser = argparse.ArgumentParser(description="Run cross-clade perturbation experiment")
    parser.add_argument("--jobs", "-j", type=int, default=8, help="Parallel worker threads")
    parser.add_argument("--force", "-f", action="store_true", help="Force re-run")
    args = parser.parse_args()

    os.makedirs(PERTURB_DIR, exist_ok=True)
    os.makedirs(LOG_DIR, exist_ok=True)

    terms = pick_perturbation_terms()
    print(f"Selected {len(terms)} perturbation terms (50 anatomical, 50 taxon-constrained)")

    injected_by_genome = {}
    claims_by_genome = {}
    for g in GENOMES:
        claims_file, records = inject_claims_for_genome(g, terms, force=args.force)
        injected_by_genome[g["tag"]] = records
        claims_by_genome[g["tag"]] = claims_file
        print(f"Injected 100 claims into {g['tag']} -> {claims_file}")

    # Run integration via shared-memory batch manifest
    tasks = [(g, c) for g in GENOMES for c in PERTURB_CONFIGS]
    to_run = []
    for g, c in tasks:
        out_file = f"{PERTURB_DIR}/{g['tag']}_{c}_perturbed_integrated.tsv"
        if not args.force and os.path.exists(out_file) and os.path.getsize(out_file) > 1000:
            continue
        to_run.append((g, c))

    if to_run:
        manifest_path = f"{PERTURB_DIR}/perturb_batch_manifest.tsv"
        with open(manifest_path, "w") as mf:
            mf.write("genome_id\tclaims_path\tout_path\ttaxon\tessential_profile\toperons_path\tgaps_path\tmax_iter\tenforce_stage2\tenforce_mode\tenable_priors\n")
            for g, cfg in to_run:
                tag = g["tag"]
                taxon = g["taxon"]
                k = g["k"]
                claims_file = claims_by_genome[tag]
                out_file = f"{PERTURB_DIR}/{tag}_{cfg}_perturbed_integrated.tsv"
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
                elif cfg == "A6_consistency_anatomical":
                    row_priors = "consistency,anatomical"
                    row_tax = taxon
                elif cfg == "A7_stage1_all":
                    row_priors = "essentiality,coherence,consistency,anatomical,gap_filling,genomic_context"
                    row_tax = taxon
                    row_ess = k
                    if has_ops:
                        row_ops = operons_file
                    if has_gaps:
                        row_gaps = gaps_file
                elif cfg in ("A8_stage2_remove", "A9_stage2_relax"):
                    row_priors = "essentiality,coherence,consistency,anatomical,gap_filling,genomic_context"
                    row_tax = taxon
                    row_ess = k
                    if has_ops:
                        row_ops = operons_file
                    if has_gaps:
                        row_gaps = gaps_file
                    row_s2 = "true"
                    row_mode = "remove" if cfg == "A8_stage2_remove" else "relax"

                mf.write("\t".join([
                    f"{tag}_{cfg}_perturb", claims_file, out_file, row_tax, row_ess,
                    row_ops, row_gaps, "-", row_s2, row_mode, row_priors
                ]) + "\n")

        cmd = [
            JAVA, "-Xmx48g", "-jar", JAR, "integrate",
            "--batch-manifest", manifest_path,
            "--batch-threads", str(args.jobs),
            "--enable-priors", "essentiality,coherence,consistency,anatomical,gap_filling,genomic_context",
            "--go-owl", f"{REF}/go.owl",
            "--pathways", f"{REF}/kegg_pathways.tsv",
            "--modules", f"{REF}/kegg_modules.tsv",
            "--ec2go", f"{REF}/ec2go.txt",
            "--lite",
        ]
        print(f"\nRunning {len(to_run)} perturbation integration jobs via shared-memory batch mode...", flush=True)
        proc = subprocess.run(cmd)
        if proc.returncode != 0:
            raise RuntimeError(f"Perturbation batch integration failed with exit code {proc.returncode}")

    # Evaluate all
    print("\nEvaluating surviving violations and CAFA metrics...", flush=True)
    summary_results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as executor:
        eval_futures = {
            executor.submit(evaluate_perturbation, g, c, injected_by_genome[g["tag"]]): (g["tag"], c)
            for g in GENOMES for c in PERTURB_CONFIGS
        }
        for future in concurrent.futures.as_completed(eval_futures):
            res = future.result()
            summary_results.append(res)
            print(f"[EVAL] {res['tag']} {res['config']}: surviving={res['n_surviving']}/100 "
                  f"(suppr={res['suppression_rate_pct']:.1f}%), fmax_micro(exp)={res['fmax_micro_exp']:.4f}, "
                  f"fmax_cafa(exp)={res['fmax_cafa_exp']:.4f}", flush=True)

    # Save summary json
    summary_file = f"{PERTURB_DIR}/perturbation_summary.json"
    with open(summary_file, "w") as f:
        json.dump(summary_results, f, indent=2)
    print(f"\nPerturbation summary saved to {summary_file}", flush=True)


if __name__ == "__main__":
    main()
