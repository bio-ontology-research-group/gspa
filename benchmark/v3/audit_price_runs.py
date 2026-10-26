#!/usr/bin/env python3
"""Compare single-predictor raw rankings with recorded production-engine outputs."""
import argparse
import json
from pathlib import Path

from prepare_claims import file_hash
from prepare_price import PREDICTORS, rank_roundtrip
import v3_metrics as M
import v3_runner as R


def raw_predictions(path):
    predictions = {}
    sources = set()
    with Path(path).open() as stream:
        for line in stream:
            row = json.loads(line)
            sources.add(row["source"])
            if row["function_type"] != "EC" or not M.is_complete_ec(row["function_id"]):
                continue
            key = (row["protein_id"], M.normalize_ec(row["function_id"]))
            predictions[key] = max(predictions.get(key, -1), row["raw_score"])
    if len(sources) != 1:
        raise ValueError("raw ranking comparison requires exactly one predictor")
    return predictions


def audit(manifest, run_report):
    genomes = []
    for genome in manifest["genomes"]:
        if genome["split"] != "test":
            continue
        gid = genome["genome_id"]
        registered = run_report["provenance"]["inputs"][gid]
        claims_hash = file_hash(genome["claims"])
        labels_hash = file_hash(genome["labels"]["ec"])
        if claims_hash != registered["claims"]["sha256"] or labels_hash != registered["labels"]["ec"]["sha256"]:
            raise ValueError(f"claims or labels changed since engine run: {gid}")
        artifact = run_report["runs"]["stage0"]["genomes"][gid]
        if artifact["status"] != "ok" or file_hash(artifact["path"]) != artifact["sha256"]:
            raise ValueError(f"prediction changed or not executed: {gid}")
        truth = R.parse_labels(genome["labels"]["ec"], "ec")
        raw = raw_predictions(genome["claims"])
        transformed = R.parse_integrated(artifact["path"], "ec")
        ro, po = M.ec_top1_outcomes(raw, truth), M.ec_top1_outcomes(transformed, truth)
        genomes.append({
            "genome_id": gid, "claims_sha256": claims_hash,
            "labels_sha256": labels_hash, "prediction_sha256": artifact["sha256"],
            "raw_top1_vector": M.top1_vector(ro), "engine_top1_vector": M.top1_vector(po),
            "labelled_outcomes_changed": sum(ro[p] != po[p] for p in ro),
            "roundtrip": rank_roundtrip(genome["claims"], artifact["path"]),
        })
    raw_total = [sum(g["raw_top1_vector"][i] for g in genomes) for i in range(3)]
    engine_total = [sum(g["engine_top1_vector"][i] for g in genomes) for i in range(3)]
    return {"genomes": genomes, "raw_top1": M.top1_metrics(raw_total),
            "engine_top1": M.top1_metrics(engine_total),
            "rank_preserved_all_genomes": all(g["roundtrip"]["raw_rank_preserved"] for g in genomes),
            "labelled_outcomes_changed": sum(g["labelled_outcomes_changed"] for g in genomes)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--run-root", type=Path, action="append", required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise ValueError("refusing to overwrite an existing audit")
    report = {"schema": "gspa.price_raw_ranking_audit/1", "script_sha256": file_hash(Path(__file__)), "predictors": {}}
    for predictor in PREDICTORS:
        paths = [root / predictor / "report.json" for root in args.run_root if (root / predictor / "report.json").exists()]
        if len(paths) != 1:
            raise ValueError(f"expected exactly one run for {predictor}, found {len(paths)}")
        manifest = R.load_manifest(str(args.inputs / "manifests" / f"v3_price_{predictor}.json"))
        result = audit(manifest, json.loads(paths[0].read_text()))
        result["run_report_sha256"] = file_hash(paths[0])
        report["predictors"][predictor] = result
        print(predictor, result["raw_top1"], "fully_preserved", result["rank_preserved_all_genomes"], flush=True)
    args.out.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
