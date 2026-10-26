#!/usr/bin/env python3
"""Run production gspa integrate ablation arms on the 18-genome cohort and evaluate on held-out test split."""
from __future__ import annotations

import csv
import hashlib
import json
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path("/home/leechuck/Public/software/gspa")
sys.path.insert(0, str(ROOT / "benchmark" / "v3"))
import homology_baseline as H  # noqa: E402
from v3_metrics import load_obo  # noqa: E402

CACHE = Path("/home/leechuck/.cache/gspa-v3-data/training-cohort")
RES = ROOT / "benchmark" / "v3" / "results" / "training-cohort"
JAR = ROOT / "gspa-cli" / "build" / "libs" / "gspa-2.2.2.jar"


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def make_manifest(claims_dir: Path, out_dir: Path, manifest_path: Path) -> list[dict]:
    out_dir.mkdir(parents=True, exist_ok=True)
    splits = list(csv.DictReader(open(RES / "splits.tsv"), delimiter="\t"))
    prov = json.loads((CACHE / "provenance.json").read_text())
    acc_to_tax = {
        r["assembly"]["used_accession"]: t
        for t, r in prov["organisms"].items()
        if r.get("assembly", {}).get("status") == "resolved"
    }
    sel = json.loads((ROOT / "benchmark/v3/results/label-rich-selection/selection_manifest.json").read_text())
    tax_to_sk = {o["taxid"]: o["superkingdom"].lower() for o in sel["organisms"]}

    with open(manifest_path, "w") as fh:
        fh.write("genome_id\tclaims_path\tout_path\ttaxon\tessential_profile\toperons_path\tgaps_path\tmax_iter\tenforce_stage2\tenforce_mode\tenable_priors\tgenome_layout_path\torthogroups_path\tcluster_consensus_path\tfba_medium\tfba_gram\n")
        for s in splits:
            gid = s["genome_id"]
            taxid = acc_to_tax[gid]
            sk = tax_to_sk[taxid]
            preset = "archaea" if "archaea" in sk else "bacteria"
            claims = claims_dir / f"{gid}.claims.jsonl"
            out = out_dir / f"{gid}.integrated.tsv"
            layout = CACHE / "layout" / f"{gid}.layout.tsv"
            fh.write(f"{gid}\t{claims}\t{out}\tNCBITaxon_{taxid}\t{preset}\t-\t-\t-\t-\t-\t-\t{layout}\t-\t-\t-\t-\n")
    return splits


def run_arm(name: str, claims_mode: str, extra_args: list[str]) -> Path:
    claims_dir = CACHE / "baseline" / "claims" / claims_mode
    out_dir = CACHE / "ablation_runs" / f"{claims_mode}_{name}"
    if out_dir.exists() and len(list(out_dir.glob("*.integrated.tsv"))) == 18:
        print(f"[skip existing] {claims_mode}_{name}")
        return out_dir
    manifest = CACHE / "ablation_runs" / f"{claims_mode}_{name}.manifest.tsv"
    manifest.parent.mkdir(parents=True, exist_ok=True)
    make_manifest(claims_dir, out_dir, manifest)
    cmd = [
        "java", "-Xmx12g", "-jar", str(JAR), "integrate",
        "--batch-manifest", str(manifest),
        "--batch-threads", "4",
        "--go-owl", str(ROOT / "benchmark/whole_genome_cafa/data/reference/go.owl"),
        "--ec2go", str(ROOT / "benchmark/whole_genome_cafa/data/reference/ec2go.txt"),
        "--pathways", str(ROOT / "benchmark/whole_genome_cafa/data/reference/kegg_pathways.tsv"),
        "--modules", str(ROOT / "benchmark/whole_genome_cafa/data/reference/kegg_modules.tsv"),
        "--reasoner-cache", str(CACHE / "reasoner-cache"),
        "--batch-stats-out", str(out_dir / "batch_stats.tsv"),
    ] + extra_args
    print(f"[running] {claims_mode}_{name} ...")
    proc = subprocess.run(cmd, capture_output=True, text=True)
    (out_dir / "run.log").write_text(proc.stdout + proc.stderr)
    if proc.returncode != 0:
        raise RuntimeError(f"{claims_mode}_{name} failed ({proc.returncode}):\n{proc.stderr[-1000:]}")
    return out_dir


def load_integrated_predictions(out_dir: Path, pid_to_seq: dict[tuple[str, str], str]) -> dict[str, dict[str, float]]:
    """Map seq_id -> {function_id -> score} across all 18 genomes."""
    preds: dict[str, dict[str, float]] = defaultdict(dict)
    for tsv in sorted(out_dir.glob("*.integrated.tsv")):
        gid = tsv.name.replace(".integrated.tsv", "")
        with open(tsv) as fh:
            for line in fh:
                if not line.strip() or line.startswith("#") or line.startswith("protein_id\t"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 5:
                    continue
                pid, ftype, fid, prob_s = parts[0], parts[1], parts[2], parts[4]
                try:
                    prob = float(prob_s)
                except ValueError:
                    continue
                sid = pid_to_seq.get((gid, pid))
                if not sid:
                    continue
                if ftype == "EC" and not fid.startswith("EC:"):
                    fid = "EC:" + fid
                if prob > preds[sid].get(fid, 0.0):
                    preds[sid][fid] = prob
    return preds


def main():
    # Copy existing cohort_rlc run if not already copied
    cohort_args = [
        "--reaction-graph", "/home/leechuck/Public/software/gapsmith/data/seed_reactions.tsv",
        "--diffusion-mets", "/home/leechuck/Public/software/gapsmith/data/diffusion_mets.tsv",
        "--reaction-ec-aliases", "/home/leechuck/Public/software/gapsmith/data/seed_Enzyme_Class_Reactions_Aliases_unique.tsv",
        "--cohort",
        "--cohort-metadata", str(CACHE / "cohort_inputs/cohort_metadata.json"),
        "--cohort-orthogroups", str(CACHE / "cohort_inputs/cohort_orthogroups.tsv"),
        "--cohort-holdout", "CLADE",
        "--rlc-suggester",
    ]
    arms = {
        "stage0": ["--enable-priors", "none"],
        "stage1": [],
        "cohort_rlc": cohort_args,
        "learned_logistic": cohort_args + [
            "--learned-model", str(RES / "learned_models/model_logistic.json"),
        ],
        "learned_gbdt": cohort_args + [
            "--learned-model", str(RES / "learned_models/model_gbdt.json"),
        ],
    }

    pid_to_seq = {}
    with open(CACHE / "prepared/proteins.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            pid_to_seq[(r["genome_id"], r["genome_protein_id"])] = r["protein_id"]

    split_of = {}
    with open(RES / "splits.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            split_of[r["genome_id"]] = r["split"]

    go = load_obo(str(CACHE / "go-basic.obo"))
    truth = H.read_truth(RES / "labels.tsv", CACHE / "prepared/proteins.tsv")

    results = {"schema": "gspa.training_cohort_ablation/1", "modes": {}}
    for mode in ("homology", "family_excluded"):
        results["modes"][mode] = {}
        for arm_name, extra in arms.items():
            out_dir = run_arm(arm_name, mode, extra)
            preds = load_integrated_predictions(out_dir, pid_to_seq)
            eval_res = H.evaluate(preds, truth, split_of, go=go, n_boot=200)
            results["modes"][mode][arm_name] = eval_res
            t_go = eval_res["test"]["go_fmax_point"]["fmax"]
            t_ec = eval_res["test"]["ec_top1_accuracy"]["estimate"]
            c_go = eval_res["calibration"]["go_fmax_point"]["fmax"]
            print(f"  => {mode}/{arm_name}: test GO Fmax={t_go:.4f}, test EC top1={t_ec:.4f}, cal GO Fmax={c_go:.4f}")

    out_file = RES / "ablation_evaluation.json"
    out_file.write_text(json.dumps(results, indent=1, sort_keys=True) + "\n")
    print("Wrote", out_file)


if __name__ == "__main__":
    main()
