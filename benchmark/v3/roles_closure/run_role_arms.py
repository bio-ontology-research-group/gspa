#!/usr/bin/env python3
"""Run the role benchmark arms with real predictor scores and summarise with intervals.

Arms (frozen plan A3): L0/L1 use the legacy bench9 claims, F0/F1 the role-fenced
Swiss-Prot DIAMOND claims; ``0`` runs ``integrate`` with priors off, ``1`` with
default priors; both receive the genome layout and the operons derived from it
under the production 300 bp rule. S-L and S-F score the same claim files with
the best single claim per candidate and no integration. Every arm is evaluated
with ``evaluate_role_truth.py`` at threshold 0.5 on accepted and pending_review
cases, and the summary adds Wilson intervals for top-1 rates, a bootstrap
interval for the mean margin, and strata by operon relation and case kind.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import random
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "benchmark" / "v3"))
import evaluate_role_truth as E  # noqa: E402

GENOMES = {  # registry genome id -> (layout key, taxon, fenced assembly)
    "bench9-ecoli-K12-MG1655": ("ecoli", "NCBITaxon_511145", "GCF_000005845.2"),
    "bench9-bsubtilis-168": ("bsubtilis", "NCBITaxon_224308", "GCF_000009045.1"),
    "bench9-mtb-H37Rv": ("mtb", "NCBITaxon_83332", "GCF_000195955.2"),
}
HEADER = ("genome_id\tclaims_path\tout_path\ttaxon\tessential_profile\toperons_path\tgaps_path\tmax_iter"
          "\tenforce_stage2\tenforce_mode\tenable_priors\tgenome_layout_path\torthogroups_path"
          "\tcluster_consensus_path\tfba_medium\tfba_gram")
BOOTSTRAP = 2000
SEED = 20261101


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def wilson(k, n, z=1.959964):
    if n == 0:
        return None
    p = k / n
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return [round(centre - half, 6), round(centre + half, 6)]


def bootstrap_mean(values, seed=SEED, rounds=BOOTSTRAP):
    if not values:
        return None
    rng = random.Random(seed)
    means = sorted(sum(rng.choice(values) for _ in values) / len(values) for _ in range(rounds))
    return [round(means[int(0.025 * rounds)], 6), round(means[int(0.975 * rounds) - 1], 6)]


def write_manifest(path, claims, outdir, layouts, operons, priors):
    lines = [HEADER]
    for genome_id, (key, taxon, _) in GENOMES.items():
        lines.append("\t".join([genome_id, str(claims[genome_id]), str(Path(outdir) / f"{genome_id}.integrated.tsv"),
                                taxon, "bacteria", str(operons[key]), "-", "-", "-", "-", priors,
                                str(layouts[key]), "-", "-", "-", "-"]))
    Path(path).write_text("\n".join(lines) + "\n")


def integrate(jar, reference, manifest, outdir, heap, log):
    cmd = ["java", f"-Xmx{heap}", "-jar", str(jar), "integrate", "--batch-manifest", str(manifest),
           "--batch-threads", "3", "--go-owl", str(reference / "go.owl"), "--ec2go", str(reference / "ec2go.txt"),
           "--pathways", str(reference / "kegg_pathways.tsv"), "--modules", str(reference / "kegg_modules.tsv"),
           "--reasoner-cache", str(reference / "reasoner-cache"),
           "--batch-stats-out", str(Path(outdir) / "batch_stats.tsv")]
    with open(log, "w") as handle:
        proc = subprocess.run(cmd, stdout=handle, stderr=subprocess.STDOUT, cwd=str(ROOT))
    return {"command": cmd, "returncode": proc.returncode}


def concatenate(paths, out):
    with Path(out).open("w") as handle:
        for path in paths:
            with Path(path).open() as source:
                for line in source:
                    handle.write(line)
    return sha256(out)


def merge_integrated(outdir, out):
    """One integrated TSV per arm: protein ids are unique across the three genomes."""
    with Path(out).open("w", newline="") as handle:
        first = True
        for genome_id in GENOMES:
            with (Path(outdir) / f"{genome_id}.integrated.tsv").open() as source:
                header = next(source)
                if first:
                    handle.write(header)
                    first = False
                for line in source:
                    handle.write(line)
    return sha256(out)


def strata(rows, tags):
    kinds = {}
    for row in rows:
        tag = tags.get(row["case_id"], {})
        kinds[row["case_id"]] = {"relation": tag.get("operon_context", {}).get("relation", "unplaced"),
                                 "kind": tag.get("case_kind", "unspecified")}
    groups = {}
    for name, key in (("relation", "relation"), ("kind", "kind")):
        for value in sorted({k[key] for k in kinds.values()}):
            subset = [r for r in rows if kinds[r["case_id"]][key] == value]
            groups[f"{name}:{value}"] = summarise(subset)
    return groups


def summarise(rows):
    scored = [r for r in rows if r["scored"]]
    n = len(scored)
    strict = sum(bool(r["top1_correct_strict"]) for r in scored)
    lenient = sum(bool(r["top1_correct_lenient"]) for r in scored)
    margins = [r["margin"] for r in scored]
    unscored = [r["case_id"] for r in rows if not r["scored"]]
    return {"n_cases": len(rows), "n_scored": n, "unscored_cases": unscored,
            "top1_strict": strict / n if n else None, "top1_strict_wilson95": wilson(strict, n),
            "top1_lenient": lenient / n if n else None, "top1_lenient_wilson95": wilson(lenient, n),
            "mrr": sum(1 / r["active_rank"] for r in scored) / n if n else None,
            "mean_margin": sum(margins) / n if n else None, "mean_margin_bootstrap95": bootstrap_mean(margins),
            "negatives_called": sum(r["n_false_positive_negatives_at_threshold"] for r in rows),
            "negatives_total": sum(len([p for p in r["scores"] if p not in r["active_protein_ids"]]) for r in rows),
            "same_operon_strict_top1": None}


def run(args):
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    layouts = {k: args.layouts / f"{k}.layout.tsv" for k in ("ecoli", "bsubtilis", "mtb")}
    operons = {k: args.operons / f"{k}.operons.tsv" for k in layouts}
    legacy = {g: args.legacy / f"{key}_claims.jsonl" for g, (key, _, _) in GENOMES.items()}
    fenced = {g: args.fenced / f"{asm}.claims.jsonl" for g, (_, _, asm) in GENOMES.items()}
    tags = {}
    for line in Path(args.registry).read_text().splitlines():
        if line.strip():
            case = json.loads(line)
            tags[case["case_id"]] = case
    record = {"schema": "gspa.role-arms/1", "jar_sha256": sha256(args.jar),
              "registry_sha256": sha256(args.registry), "negative_only_sha256": sha256(args.negative_only),
              "layouts": {k: sha256(v) for k, v in layouts.items()}, "operons": {k: sha256(v) for k, v in operons.items()},
              "claims": {"legacy": {g: sha256(p) for g, p in legacy.items()}, "fenced": {g: sha256(p) for g, p in fenced.items()}},
              "threshold": args.threshold, "arms": {}}
    arms = {"L0": (legacy, "none"), "L1": (legacy, "-"), "F0": (fenced, "none"), "F1": (fenced, "-")}
    for name, (claims, priors) in arms.items():
        arm_dir = out / name
        arm_dir.mkdir(exist_ok=True)
        manifest = arm_dir / "batch_manifest.tsv"
        write_manifest(manifest, claims, arm_dir, layouts, operons, priors)
        info = integrate(args.jar, args.reference, manifest, arm_dir, args.heap, arm_dir / "integrate.log")
        if info["returncode"] != 0:
            record["arms"][name] = {"status": "failed", **info}
            continue
        merged = arm_dir / "all.integrated.tsv"
        info["scores_sha256"] = merge_integrated(arm_dir, merged)
        info["stages"] = {g: json.loads((arm_dir / f"{g}.integrated.tsv.stages.json").read_text()).get("stages")
                          for g in GENOMES}
        summary = E.evaluate(args.registry, merged, "integrated_tsv", arm_dir / "evaluation",
                             threshold=args.threshold, negative_only=args.negative_only)
        rows = [json.loads(l) for l in (arm_dir / "evaluation" / "per_case.jsonl").read_text().splitlines()]
        record["arms"][name] = {"status": "executed", **info, "claims": "legacy" if claims is legacy else "fenced",
                                "priors": priors, "overall": summarise(rows), "strata": strata(rows, tags),
                                "negative_only_called": summary["negative_only_called_at_threshold"],
                                "negative_only_total": summary["negative_only_total"]}
    for name, claims in (("S-L", legacy), ("S-F", fenced)):
        arm_dir = out / name
        arm_dir.mkdir(exist_ok=True)
        merged = arm_dir / "all.claims.jsonl"
        digest = concatenate([claims[g] for g in GENOMES], merged)
        summary = E.evaluate(args.registry, merged, "claims_jsonl", arm_dir / "evaluation",
                             threshold=args.threshold, negative_only=args.negative_only)
        rows = [json.loads(l) for l in (arm_dir / "evaluation" / "per_case.jsonl").read_text().splitlines()]
        record["arms"][name] = {"status": "executed", "scores_sha256": digest, "claims": "legacy" if claims is legacy else "fenced",
                                "priors": "not_applicable", "overall": summarise(rows), "strata": strata(rows, tags),
                                "negative_only_called": summary["negative_only_called_at_threshold"],
                                "negative_only_total": summary["negative_only_total"]}
    per_case = {}
    for name in record["arms"]:
        path = out / name / "evaluation" / "per_case.jsonl"
        if path.exists():
            for line in path.read_text().splitlines():
                row = json.loads(line)
                per_case.setdefault(row["case_id"], {})[name] = {
                    "scored": row["scored"], "top1_strict": row["top1_correct_strict"], "rank": row["active_rank"],
                    "margin": row["margin"], "scores": row["scores"]}
    record["per_case"] = {cid: {"relation": tags[cid].get("operon_context", {}).get("relation"),
                                "kind": tags[cid].get("case_kind", "unspecified"), "status": tags[cid]["status"],
                                "arms": arms_} for cid, arms_ in per_case.items()}
    (out / "role_arms.json").write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    print(json.dumps({n: (a["status"], a.get("overall", {}).get("n_scored")) for n, a in record["arms"].items()}))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("registry", "negative-only", "layouts", "operons", "legacy", "fenced", "jar", "reference", "out"):
        parser.add_argument("--" + name, required=True, type=Path)
    parser.add_argument("--heap", default="16g")
    parser.add_argument("--threshold", type=float, default=0.5)
    run(parser.parse_args(argv))


if __name__ == "__main__":
    os.environ.setdefault("TMPDIR", str(Path.home() / "tmp"))
    main()
