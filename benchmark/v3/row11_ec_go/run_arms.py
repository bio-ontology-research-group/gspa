#!/usr/bin/env python3
"""Run EC-aware anomaly extraction and integration arms on merged Price claim sets.

Models ``roles_closure/run_price_anomalies.py`` and reuses ``genus_blocked_split``,
``manifest`` and ``arms_for`` by import. Genomes are read from
``<parts>/<set>/<assembly>/part.json`` with SHA-256 verification, each genome row
carries ``operons: <operons>/<assembly>.operons.tsv``, and ``shared_args`` adds
``--anomaly-ec-aware --anomaly-floor <F> --anomaly-max-candidates <K> --anomaly-max-locus <N>``.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "benchmark" / "v3" / "roles_closure"))
sys.path.insert(0, str(ROOT / "benchmark" / "v3"))
import run_price_anomalies as RPA  # noqa: E402
import v3_runner as runner  # noqa: E402
from run_price_anomalies import arms_for, genus_blocked_split  # noqa: E402

SETS = ("ec6", "ec6c")


def _format_floor(floor: float) -> str:
    if round(floor, 2) == floor:
        return f"{floor:.2f}"
    return str(floor)


def genomes_for(
    claim_set: str,
    parts: Path,
    layouts: Path,
    operons: Path,
    split: dict,
) -> list[dict]:
    """Build genome rows for ``claim_set`` with verified claims and derived operons."""
    rows = RPA.genomes_for(claim_set, Path(parts), Path(layouts), split)
    operons_dir = Path(operons)
    for row in rows:
        gid = row["genome_id"]
        row["operons"] = str(operons_dir / f"{gid}.operons.tsv")
    return rows


def manifest(
    jar: Path,
    reference: Path,
    models: Path,
    anomalies: Path,
    threads: int,
    heap: str,
    anomaly_floor: float = 0.20,
    anomaly_max_candidates: int = 10,
    anomaly_max_locus: int = 150,
) -> dict:
    """Build runner manifest with EC-aware anomaly extractor flags and this script in input_files."""
    man = RPA.manifest(Path(jar), Path(reference), Path(models), Path(anomalies), threads, heap)
    man["gspa"]["input_files"] = [str(Path(__file__).resolve())]
    man["gspa"]["shared_args"] = list(man["gspa"]["shared_args"]) + [
        "--anomaly-ec-aware",
        "--anomaly-floor",
        _format_floor(anomaly_floor),
        "--anomaly-max-candidates",
        str(anomaly_max_candidates),
        "--anomaly-max-locus",
        str(anomaly_max_locus),
    ]
    return man


def limitations_for(claim_set: str) -> list[str]:
    if claim_set == "ec6c":
        return [
            "Merged six-predictor enzyme claims (clean, dpz, enzbert, graphec, mapred, topec) after a uniform raw-score floor (raw_score >= 0.05); ec6c adds Swiss-Prot DIAMOND contract claims.",
            "Contract-v3 learned models score only candidates carrying a predictor contract and abstain elsewhere; measured, not assumed.",
            "Price genomes are a reannotation benchmark, not a random genome sample.",
        ]
    return [
        "Merged six-predictor enzyme claims (clean, dpz, enzbert, graphec, mapred, topec) after a uniform raw-score floor (raw_score >= 0.05); EC claims only.",
        "Contract-v3 learned models may abstain on these claims; measured, not assumed.",
        "Price genomes are a reannotation benchmark, not a random genome sample.",
    ]


def run(
    claim_set: str,
    out: Path,
    parts: Path,
    layouts: Path,
    operons: Path,
    metadata: Path,
    jar: Path,
    reference: Path,
    models: Path,
    threads: int = 8,
    heap: str = "48g",
    anomaly_floor: float = 0.20,
    anomaly_max_candidates: int = 10,
    anomaly_max_locus: int = 150,
) -> dict:
    out = Path(out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    split = genus_blocked_split(Path(metadata))
    genomes = genomes_for(claim_set, Path(parts), Path(layouts), Path(operons), split)
    man = manifest(
        Path(jar),
        Path(reference),
        Path(models),
        out / "anomalies-placeholder.jsonl",
        threads,
        heap,
        anomaly_floor=anomaly_floor,
        anomaly_max_candidates=anomaly_max_candidates,
        anomaly_max_locus=anomaly_max_locus,
    )
    receipts, abstained = {}, {}
    for arm in arms_for(man):
        directory = out / arm["name"]
        record = runner.load_run_record(str(directory), genomes, man, arm)
        if len(record["genomes"]) != len(genomes):
            print(f"[{claim_set}] running {arm['name']} on {len(genomes)} genomes", flush=True)
            record = runner.run_arm(man, arm, str(directory), genomes)
        if len(record["genomes"]) == len(genomes):
            receipts[arm["name"]] = str(directory / runner.RUN_RECORD)
            continue
        reasons = set(record["failed"].values())
        if (
            arm["name"].startswith("learned_")
            and record.get("returncode") == 0
            and reasons == {"invalid execution evidence: required stage learned_integration did not execute"}
        ):
            abstained[arm["name"]] = {
                "genomes": len(genomes),
                "reason": next(iter(reasons)),
                "run_record": str(directory / runner.RUN_RECORD),
            }
            print(f"[{claim_set}] {arm['name']} abstained on all {len(genomes)} genomes", flush=True)
            continue
        raise SystemExit(f"{claim_set}/{arm['name']} failed: {record['failed']}")

    source = {
        "schema": "gspa.agent-evidence-sources/2",
        "runs": receipts,
        "abstained_arms": abstained,
        "splits": {g["genome_id"]: split[g["genome_id"]] for g in genomes},
        "split_audit": {
            "path": str(metadata),
            "sha256": runner.sha256_file(str(metadata)),
            "status": "ok",
            "rule": "NCBI genus blocks, seed 20261101, whole genera to calibration until at least seven genomes",
        },
        "predictor": claim_set,
        "limitations": limitations_for(claim_set),
    }
    (out / "sources.json").write_text(json.dumps(source, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return source


def run_single_stage1(
    part_dir: Path,
    out_dir: Path,
    layouts: Path,
    operons: Path,
    metadata: Path,
    jar: Path,
    reference: Path,
    models: Path,
    threads: int = 8,
    heap: str = "168g",
    anomaly_floor: float = 0.20,
    anomaly_max_candidates: int = 10,
    anomaly_max_locus: int = 150,
) -> dict:
    """Run ``stage1`` on a single assembly part folder (used by ``MODE=fullcheck``)."""
    part_dir = Path(part_dir).resolve()
    out_dir = Path(out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    part = json.loads((part_dir / "part.json").read_text(encoding="utf-8"))
    genome_id = part["assembly"]
    claims = part_dir / part["claims"].get("path", "claims.jsonl")
    if runner.sha256_file(str(claims)) != part["claims"]["sha256"]:
        raise SystemExit(f"{claims}: published claim digest does not match part.json")
    report = json.loads((Path(layouts) / "report.json").read_text(encoding="utf-8"))["assemblies"]
    split = genus_blocked_split(Path(metadata))
    genomes = [
        {
            "genome_id": genome_id,
            "split": split[genome_id]["split"],
            "claims": str(claims),
            "taxon": "NCBITaxon_" + str(report[genome_id]["taxid"]),
            "essential_profile": "bacteria",
            "genome_layout": str(Path(layouts) / f"{genome_id}.layout.tsv"),
            "operons": str(Path(operons) / f"{genome_id}.operons.tsv"),
        }
    ]
    man = manifest(
        Path(jar),
        Path(reference),
        Path(models),
        out_dir / "anomalies-placeholder.jsonl",
        threads,
        heap,
        anomaly_floor=anomaly_floor,
        anomaly_max_candidates=anomaly_max_candidates,
        anomaly_max_locus=anomaly_max_locus,
    )
    arm = {"name": "stage1", "args": []}
    stage1_dir = out_dir / "stage1"
    record = runner.load_run_record(str(stage1_dir), genomes, man, arm)
    if len(record["genomes"]) != 1:
        record = runner.run_arm(man, arm, str(stage1_dir), genomes)
    if len(record["genomes"]) != 1:
        raise SystemExit(f"fullcheck stage1 failed for {genome_id}: {record['failed']}")
    return record


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--set", dest="claim_set", choices=SETS, help="Claim set name (ec6 or ec6c)")
    parser.add_argument("--parts", type=Path, help="Merged parts root (<parts>/<set>/<assembly>/part.json)")
    parser.add_argument("--fullcheck-part", type=Path, default=None, help="Single assembly part folder for stage1 fullcheck")
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--layouts", required=True, type=Path)
    parser.add_argument("--operons", required=True, type=Path)
    parser.add_argument("--cohort-metadata", required=True, type=Path)
    parser.add_argument("--jar", required=True, type=Path)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--models", required=True, type=Path)
    parser.add_argument("--batch-threads", type=int, default=8)
    parser.add_argument("--heap", default="48g")
    parser.add_argument("--anomaly-floor", type=float, default=0.20)
    parser.add_argument("--anomaly-max-candidates", type=int, default=10)
    parser.add_argument("--anomaly-max-locus", type=int, default=150)
    args = parser.parse_args(argv)
    if args.fullcheck_part is not None:
        rec = run_single_stage1(
            part_dir=args.fullcheck_part,
            out_dir=args.out,
            layouts=args.layouts,
            operons=args.operons,
            metadata=args.cohort_metadata,
            jar=args.jar,
            reference=args.reference,
            models=args.models,
            threads=args.batch_threads,
            heap=args.heap,
            anomaly_floor=args.anomaly_floor,
            anomaly_max_candidates=args.anomaly_max_candidates,
            anomaly_max_locus=args.anomaly_max_locus,
        )
        print(json.dumps({"fullcheck": str(args.fullcheck_part), "genomes": sorted(rec["genomes"])}, sort_keys=True))
        return 0
    if not args.claim_set or args.parts is None:
        parser.error("--set and --parts are required unless --fullcheck-part is specified")
    source = run(
        claim_set=args.claim_set,
        out=args.out,
        parts=args.parts,
        layouts=args.layouts,
        operons=args.operons,
        metadata=args.cohort_metadata,
        jar=args.jar,
        reference=args.reference,
        models=args.models,
        threads=args.batch_threads,
        heap=args.heap,
        anomaly_floor=args.anomaly_floor,
        anomaly_max_candidates=args.anomaly_max_candidates,
        anomaly_max_locus=args.anomaly_max_locus,
    )
    print(
        json.dumps(
            {
                "set": args.claim_set,
                "arms": sorted(source["runs"]),
                "abstained": sorted(source["abstained_arms"]),
                "genomes": len(source["splits"]),
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    os.environ.setdefault("TMPDIR", str(Path.home() / "tmp"))
    raise SystemExit(main())
