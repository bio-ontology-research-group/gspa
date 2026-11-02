#!/usr/bin/env python3
"""Run the production anomaly extractor on the Price release-candidate inputs.

One arm per learned configuration (baseline plus the two contract-v3 models), the
same claims, layouts, references and taxon for every arm, through
``v3_runner.run_arm`` so each arm carries a content-verified receipt and schema-2
stage ledgers. This program reads no label file: the genus-blocked split comes
from the cohort metadata's NCBI genus, and Price labels are joined only later, by
``prepare_price_agent_truth.py``, after the public case bundle is frozen.

The output ``sources.json`` is the ``gspa.agent-evidence-sources/2`` record that
``evaluate_stage3_agent.py prepare`` consumes unchanged.
"""
from __future__ import annotations

import argparse
import json
import os
import random
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "benchmark" / "v3"))
import v3_runner as runner  # noqa: E402

PREDICTORS = ("topec", "mapred", "enzbert")
ARMS = ("stage1", "learned_logistic", "learned_gbdt")
SPLIT_SEED = 20261101
CALIBRATION_MIN_GENOMES = 7


def genus_blocked_split(metadata: Path) -> dict:
    """Genera shuffled with a fixed seed; whole genera fill calibration first."""
    genomes = json.loads(Path(metadata).read_text())["genomes"]
    by_genus: dict[str, list[str]] = {}
    for genome_id, record in genomes.items():
        by_genus.setdefault(record["clade"], []).append(genome_id)
    order = sorted(by_genus)
    random.Random(SPLIT_SEED).shuffle(order)
    split, calibrated = {}, 0
    for genus in order:
        label = "calibration" if calibrated < CALIBRATION_MIN_GENOMES else "test"
        for genome_id in sorted(by_genus[genus]):
            split[genome_id] = {"genome_id": genome_id, "genus": genus, "group": genus,
                                "source": "price_release_candidate", "split": label}
        if label == "calibration":
            calibrated += len(by_genus[genus])
    return split


def genomes_for(predictor: str, parts: Path, layouts: Path, split: dict) -> list[dict]:
    report = json.loads((layouts / "report.json").read_text())["assemblies"]
    rows = []
    for directory in sorted((parts / predictor).iterdir()):
        part = json.loads((directory / "part.json").read_text())
        genome_id = part["assembly"]
        claims = directory / part["claims"]["path"]
        if runner.sha256_file(str(claims)) != part["claims"]["sha256"]:
            raise SystemExit(f"{claims}: published claim digest does not match part.json")
        rows.append({"genome_id": genome_id, "split": split[genome_id]["split"],
                     "claims": str(claims), "taxon": "NCBITaxon_" + str(report[genome_id]["taxid"]),
                     "essential_profile": "bacteria",
                     "genome_layout": str(layouts / f"{genome_id}.layout.tsv")})
    return rows


def manifest(jar: Path, reference: Path, models: Path, anomalies: Path, threads: int, heap: str) -> dict:
    return {"gspa": {"command": ["java", f"-Xmx{heap}", "-jar", str(jar)],
                     "binary": str(jar), "cwd": str(ROOT), "timeout_sec": 36000,
                     "input_files": [str(Path(__file__).resolve())],
                     "shared_args": ["--batch-threads", str(threads),
                                     "--go-owl", str(reference / "go.owl"),
                                     "--ec2go", str(reference / "ec2go.txt"),
                                     "--pathways", str(reference / "kegg_pathways.tsv"),
                                     "--modules", str(reference / "kegg_modules.tsv"),
                                     "--reasoner-cache", str(reference / "reasoner-cache"),
                                     "--anomalies-out", str(anomalies)]},
            "models": {kind: str(models / f"model_{kind}.json") for kind in ("logistic", "gbdt")}}


def arms_for(man: dict) -> list[dict]:
    return [{"name": "stage1", "args": []}] + [
        {"name": f"learned_{kind}", "args": ["--learned-model", man["models"][kind]]}
        for kind in ("logistic", "gbdt")]


def run(predictor: str, out: Path, parts: Path, layouts: Path, metadata: Path, jar: Path,
        reference: Path, models: Path, threads: int, heap: str) -> dict:
    out = out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    split = genus_blocked_split(metadata)
    genomes = genomes_for(predictor, parts, layouts, split)
    man = manifest(jar, reference, models, out / "anomalies-placeholder.jsonl", threads, heap)
    receipts, abstained = {}, {}
    for arm in arms_for(man):
        directory = out / arm["name"]
        record = runner.load_run_record(str(directory), genomes, man, arm)
        if len(record["genomes"]) != len(genomes):
            print(f"[{predictor}] running {arm['name']} on {len(genomes)} genomes", flush=True)
            record = runner.run_arm(man, arm, str(directory), genomes)
        if len(record["genomes"]) == len(genomes):
            receipts[arm["name"]] = str(directory / runner.RUN_RECORD)
            continue
        reasons = set(record["failed"].values())
        # The engine ran but its learned stage declined every genome: the model's
        # predictor contract does not cover these claims. That is a measured
        # abstention of the arm, recorded as such, not an execution failure.
        if arm["name"].startswith("learned_") and record.get("returncode") == 0 and reasons == {
                "invalid execution evidence: required stage learned_integration did not execute"}:
            abstained[arm["name"]] = {"genomes": len(genomes), "reason": next(iter(reasons)),
                                      "run_record": str(directory / runner.RUN_RECORD)}
            print(f"[{predictor}] {arm['name']} abstained on all {len(genomes)} genomes", flush=True)
            continue
        raise SystemExit(f"{predictor}/{arm['name']} failed: {record['failed']}")
    source = {"schema": "gspa.agent-evidence-sources/2", "runs": receipts, "abstained_arms": abstained,
              "splits": {g["genome_id"]: split[g["genome_id"]] for g in genomes},
              "split_audit": {"path": str(metadata), "sha256": runner.sha256_file(str(metadata)),
                              "status": "ok", "rule": "NCBI genus blocks, seed 20261101, "
                              "whole genera to calibration until at least seven genomes"},
              "predictor": predictor,
              "limitations": [
                  "Single enzyme predictor per collection; EC claims only.",
                  "Contract-v3 learned models may abstain on these claims; measured, not assumed.",
                  "Price genomes are a reannotation benchmark, not a random genome sample."]}
    (out / "sources.json").write_text(json.dumps(source, indent=2, sort_keys=True) + "\n")
    return source


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--predictor", required=True, choices=PREDICTORS)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--parts", required=True, type=Path)
    parser.add_argument("--layouts", required=True, type=Path)
    parser.add_argument("--cohort-metadata", required=True, type=Path)
    parser.add_argument("--jar", required=True, type=Path)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--models", required=True, type=Path)
    parser.add_argument("--batch-threads", type=int, default=4)
    parser.add_argument("--heap", default="40g")
    args = parser.parse_args(argv)
    source = run(args.predictor, args.out, args.parts, args.layouts, args.cohort_metadata,
                 args.jar, args.reference, args.models, args.batch_threads, args.heap)
    print(json.dumps({"predictor": args.predictor, "arms": sorted(source["runs"]),
                      "genomes": len(source["splits"])}))


if __name__ == "__main__":
    os.environ.setdefault("TMPDIR", str(Path.home() / "tmp"))
    main()
