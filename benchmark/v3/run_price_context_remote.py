#!/usr/bin/env python3
"""Remote runner for Price-148 context, cohort, and learned-model ablations + locus decoy evaluation."""
from __future__ import annotations

import json
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

BASE = Path("/data/hohndor/gspa-v3/price-context-ablation")
PRICE_IN = Path("/data/hohndor/gspa-v3/price-development-e2303a9/inputs")
LAYOUT_DIR = Path("/data/hohndor/gspa-v3/price-context-preparation/layout")
COHORT_DIR = BASE / "cohort"
REF_DIR = Path("/data/hohndor/gspa/reference")
RXN_DIR = BASE / "rxn_ref"
MODELS_DIR = BASE / "learned_models"
JAR = BASE / "gspa-2.2.2.jar"
OUT_RES = BASE / "results"

VARIANTS = ("clean", "dpz", "enzbert", "graphec", "mapred", "topec", "ensemble_pool4")


def build_arms() -> list[dict]:
    rlc_args = [
        "--reaction-graph", str(RXN_DIR / "seed_reactions.tsv"),
        "--diffusion-mets", str(RXN_DIR / "diffusion_mets.tsv"),
        "--reaction-ec-aliases", str(RXN_DIR / "seed_Enzyme_Class_Reactions_Aliases_unique.tsv"),
        "--rlc-suggester",
    ]
    cohort_args = rlc_args + [
        "--cohort",
        "--cohort-metadata", str(COHORT_DIR / "cohort_metadata.json"),
        "--cohort-orthogroups", str(COHORT_DIR / "cohort_orthogroups.tsv"),
        "--cohort-holdout", "CLADE",
    ]
    return [
        {"name": "stage0", "args": ["--enable-priors", "none"]},
        {"name": "stage1", "args": []},
        {"name": "stage1_rlc", "args": rlc_args},
        {"name": "stage1_rlc_cohort", "args": cohort_args},
        {"name": "learned_logistic", "args": cohort_args + ["--learned-model", str(MODELS_DIR / "model_logistic.json")]},
        {"name": "learned_gbdt", "args": cohort_args + ["--learned-model", str(MODELS_DIR / "model_gbdt.json")]},
    ]


def build_manifest(variant: str) -> Path:
    man_dir = BASE / "manifests"
    man_dir.mkdir(parents=True, exist_ok=True)
    cohort_meta = json.loads((COHORT_DIR / "cohort_metadata.json").read_text())
    genomes = []
    for gid in sorted(cohort_meta["genomes"].keys()):
        genomes.append({
            "genome_id": gid,
            "split": "test",
            "essential_profile": "bacteria",
            "claims": str(PRICE_IN / "claims" / variant / f"{gid}.claims.jsonl"),
            "genome_layout": str(LAYOUT_DIR / f"{gid}.layout.tsv"),
            "labels": {
                "ec": str(PRICE_IN / "labels" / f"{gid}.ec.tsv"),
            },
        })
    manifest = {
        "schema_version": 1,
        "experiment_id": f"v3-price148-{variant}-context-ablation",
        "dataset_status": "real",
        "description": (
            f"Price-148 enzyme attribution and locus decoy benchmark, variant {variant}. "
            "Compares stage0 (sequence evidence only), stage1 (pathway priors), stage1_rlc "
            "(reaction local context suggester + layout), stage1_rlc_cohort (two-pass clade-held-out "
            "cohort recurrence over 35,298 MMseqs2 families), and learned_logistic / learned_gbdt "
            "(trained exclusively on the genus/family-disjoint swissprot_label_rich_v1 cohort)."
        ),
        "declared_controls": {
            "context": "target genomic neighborhoods (+-5 genes) with full coordinates (price-layout-v4) and 22-genome MMseqs2 family cohort (price-cohort-v2, clade holdout)",
            "learned_model_provenance": "trained on swissprot_label_rich_v1 (18 genomes, 0 Price genera); Price labels never used in training or calibration",
            "predictor_training_overlap": "not assessed: the six predictors' training sets may contain Price enzymes",
            "test_only": "all 22 Price genomes evaluated as test; no parameter or threshold fitted on Price",
        },
        "baseline_arm": "stage0",
        "arms": build_arms(),
        "genomes": genomes,
        "gspa": {
            "binary": str(JAR),
            "command": ["java", "-Xmx24g", "-jar", str(JAR)],
            "shared_args": [
                "--ec2go", str(REF_DIR / "ec2go.txt"),
                "--pathways", str(REF_DIR / "kegg_pathways.tsv"),
                "--modules", str(REF_DIR / "kegg_modules.tsv"),
                "--batch-threads", "8",
            ],
            "timeout_seconds": 1800,
        },
        "settings": {
            "bootstrap": 1000,
            "seed": 13,
            "threshold": 0.5,
        },
    }
    out_path = man_dir / f"v3_price_{variant}.json"
    out_path.write_text(json.dumps(manifest, indent=2) + "\n")
    return out_path


def evaluate_locus_decoys(variant: str, run_dir: Path, arms: list[str]) -> dict:
    """Evaluate target vs same-locus (+-5 genes on same contig) decoy ranking for each arm."""
    labels_by_genome = {}
    for lf in sorted((PRICE_IN / "labels").glob("*.ec.tsv")):
        gid = lf.name.replace(".ec.tsv", "")
        t_ec = {}
        for line in open(lf):
            if line.strip() and not line.startswith("#") and not line.startswith("protein_id"):
                p, ec = line.rstrip("\n").split("\t")[:2]
                t_ec[p] = ec.replace("EC:", "")
        labels_by_genome[gid] = t_ec

    layouts = {}
    for gid in labels_by_genome:
        l_idx, l_contig = {}, {}
        with open(LAYOUT_DIR / f"{gid}.layout.tsv") as fh:
            for idx, line in enumerate(fh):
                if not line.startswith("#") and not line.startswith("protein_id"):
                    parts = line.rstrip("\n").split("\t")
                    l_idx[parts[0]] = idx
                    l_contig[parts[0]] = parts[1]
        layouts[gid] = (l_idx, l_contig)

    raw_claims_by_genome = {}
    for gid in labels_by_genome:
        by_prot = defaultdict(dict)
        with open(PRICE_IN / "claims" / variant / f"{gid}.claims.jsonl") as fh:
            for line in fh:
                if line.strip():
                    r = json.loads(line)
                    pid, fid, sc = r["protein_id"], r["function_id"].replace("EC:", ""), float(r["raw_score"])
                    if sc > by_prot[pid].get(fid, 0.0):
                        by_prot[pid][fid] = sc
        raw_claims_by_genome[gid] = by_prot

    arm_summaries = {}
    for arm in arms:
        arm_dir = run_dir / "arms" / arm
        n_loci = 0
        n_retained = 0
        n_protein_top1 = 0
        n_protein_top3 = 0
        n_locus_top1_all = 0
        n_decoy_loci = 0
        n_locus_top1_decoy = 0
        mrr_all_sum = 0.0
        mrr_decoy_sum = 0.0
        margins_decoy = []
        role_counts = Counter()

        for gid, t_ec in sorted(labels_by_genome.items()):
            l_idx, l_contig = layouts[gid]
            raw_prot = raw_claims_by_genome[gid]
            pred_prot = defaultdict(dict)
            tsv = arm_dir / f"{gid}.integrated.tsv"
            if tsv.exists():
                with open(tsv) as fh:
                    for line in fh:
                        if not line.strip() or line.startswith("#") or line.startswith("protein_id\t"):
                            continue
                        parts = line.rstrip("\n").split("\t")
                        if len(parts) < 5 or parts[1] != "EC":
                            continue
                        pid, fid, prob = parts[0], parts[2].replace("EC:", ""), float(parts[4])
                        if prob > pred_prot[pid].get(fid, 0.0):
                            pred_prot[pid][fid] = prob

            learned_json = arm_dir / f"{gid}.integrated.tsv.learned.json"
            if learned_json.exists():
                ld = json.loads(learned_json.read_text())
                for c in ld.get("candidates", []):
                    role_counts[c.get("role", "UNKNOWN")] += 1

            for tp, true_ec in sorted(t_ec.items()):
                n_loci += 1
                p_scores = pred_prot.get(tp, {})
                t_score = p_scores.get(true_ec, 0.0)
                if true_ec in p_scores:
                    n_retained += 1
                if p_scores:
                    higher = sum(1 for ec, sc in p_scores.items() if ec != true_ec and sc >= t_score)
                    if t_score > 0 and higher == 0:
                        n_protein_top1 += 1
                    if t_score > 0 and higher < 3:
                        n_protein_top3 += 1

                t_i = l_idx.get(tp)
                t_c = l_contig.get(tp)
                neigh = [
                    p for p in raw_prot
                    if p != tp and l_contig.get(p) == t_c and t_i is not None and abs(l_idx[p] - t_i) <= 5
                ]
                has_decoy = any(true_ec in raw_prot[np] for np in neigh)
                max_decoy_score = max((pred_prot.get(np, {}).get(true_ec, 0.0) for np in neigh), default=0.0)
                higher_in_locus = sum(
                    1 for np in neigh if pred_prot.get(np, {}).get(true_ec, 0.0) >= t_score
                )
                rank_in_locus = (higher_in_locus + 1) if t_score > 0 else (len(neigh) + 1)
                mrr_all_sum += (1.0 / rank_in_locus) if t_score > 0 else 0.0
                if t_score > 0 and higher_in_locus == 0:
                    n_locus_top1_all += 1
                if has_decoy:
                    n_decoy_loci += 1
                    if t_score > 0 and higher_in_locus == 0:
                        n_locus_top1_decoy += 1
                    mrr_decoy_sum += (1.0 / rank_in_locus) if t_score > 0 else 0.0
                    margins_decoy.append(t_score - max_decoy_score)

        arm_summaries[arm] = {
            "n_loci": n_loci,
            "candidate_recall": n_retained / n_loci if n_loci else 0.0,
            "n_retained": n_retained,
            "protein_ec_top1": n_protein_top1 / n_loci if n_loci else 0.0,
            "n_protein_ec_top1": n_protein_top1,
            "protein_ec_top3": n_protein_top3 / n_loci if n_loci else 0.0,
            "n_protein_ec_top3": n_protein_top3,
            "locus_top1_all_148": n_locus_top1_all / n_loci if n_loci else 0.0,
            "n_locus_top1_all_148": n_locus_top1_all,
            "locus_mrr_all_148": mrr_all_sum / n_loci if n_loci else 0.0,
            "n_same_ec_decoy_loci": n_decoy_loci,
            "locus_top1_decoy_loci": n_locus_top1_decoy / n_decoy_loci if n_decoy_loci else 0.0,
            "n_locus_top1_decoy_loci": n_locus_top1_decoy,
            "locus_mrr_decoy_loci": mrr_decoy_sum / n_decoy_loci if n_decoy_loci else 0.0,
            "mean_margin_over_decoy": (sum(margins_decoy) / len(margins_decoy)) if margins_decoy else 0.0,
            "learned_role_counts": dict(role_counts),
        }
    return arm_summaries


def main():
    OUT_RES.mkdir(parents=True, exist_ok=True)
    arms = [a["name"] for a in build_arms()]
    only_variant = None
    if len(sys.argv) >= 3 and sys.argv[1] == "--variant":
        only_variant = sys.argv[2]

    variants_to_run = (only_variant,) if only_variant else VARIANTS
    all_locus = {"schema": "gspa.price_locus_decoy_evaluation/1", "variants": {}}
    summary_table = {}

    for variant in variants_to_run:
        man_path = build_manifest(variant)
        run_dir = BASE / "runs" / variant
        if not (run_dir / "report.json").exists():
            if run_dir.exists():
                shutil.rmtree(run_dir)
            cmd = [
                sys.executable, str(BASE / "v3_runner.py"),
                "run", str(man_path), "--out", str(run_dir),
            ]
            print(f"[running v3_runner] {variant} ...", flush=True)
            proc = subprocess.run(cmd, capture_output=True, text=True)
            if proc.returncode != 0:
                raise RuntimeError(f"{variant} failed:\n{proc.stdout}\n{proc.stderr}")
        if only_variant:
            print(f"Completed single variant {variant}", flush=True)
            return
        shutil.copy2(run_dir / "report.json", OUT_RES / f"price-{variant}-report.json")
        rep = json.loads((run_dir / "report.json").read_text())
        locus_eval = evaluate_locus_decoys(variant, run_dir, arms)
        all_locus["variants"][variant] = locus_eval
        summary_table[variant] = {}
        for arm in arms:
            ec_m = rep["evaluation"]["kinds"]["ec"]["arms"][arm]
            top1 = ec_m["ec_top1"]["top1_accuracy"]
            fmax = ec_m["global_fmax"]["fmax"]
            ece = ec_m["calibration"]["ece"]
            brier = ec_m["calibration"]["brier"]
            loc_d = locus_eval[arm]["locus_top1_decoy_loci"]
            summary_table[variant][arm] = {
                "ec_top1": top1,
                "ec_fmax": fmax,
                "ece": ece,
                "brier": brier,
                "candidate_recall": locus_eval[arm]["candidate_recall"],
                "locus_top1_all_148": locus_eval[arm]["locus_top1_all_148"],
                "locus_mrr_all_148": locus_eval[arm]["locus_mrr_all_148"],
                "n_same_ec_decoy_loci": locus_eval[arm]["n_same_ec_decoy_loci"],
                "locus_top1_decoy": loc_d,
                "locus_mrr_decoy": locus_eval[arm]["locus_mrr_decoy_loci"],
            }
            print(f"  {variant}/{arm}: EC top1={top1:.4f}, Fmax={fmax:.4f}, ECE={ece:.4f}, decoy_locus_top1={loc_d:.4f}", flush=True)

    (OUT_RES / "locus_decoy_evaluation.json").write_text(json.dumps(all_locus, indent=2, sort_keys=True) + "\n")
    (OUT_RES / "summary.json").write_text(json.dumps(summary_table, indent=2, sort_keys=True) + "\n")
    print("Completed all 7 Price variants ->", OUT_RES, flush=True)


if __name__ == "__main__":
    main()
