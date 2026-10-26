#!/usr/bin/env python3
"""Run measured knockout cases through the production metabolic experiment engine.

This is phenotype validation, not gene-function recovery. Cases selected by a
model's auxotrophy prediction do not establish population recall or specificity.
Labels are passed only to the evaluator, after each simulation has completed.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import sys
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import metabolic_experiments as ME


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def phenotype(rec: dict | None, threshold: float) -> tuple[str, str]:
    """Absolute biomass flux; never turn a failed or missing solve into no growth."""
    if rec is None:
        return "unknown", "experiment_missing"
    if rec.get("status") != "completed":
        return "unknown", rec.get("status", "status_missing")
    if rec.get("outcome") == "infeasible_with_maintenance":
        # Phenotype interpretation only. This is not evidence of protein function.
        return "no_growth", "infeasible_under_declared_maintenance"
    value = rec.get("measurements", {}).get("max_growth")
    if not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
        return "unknown", "growth_unavailable_or_invalid"
    return ("growth" if value >= threshold else "no_growth"), "optimal"


def score_case(case: dict, result: dict, control_id: str, knockout_id: str) -> dict:
    by_id = {r["experiment_id"]: r for r in result["experiments"]}
    threshold = result["config"]["growth_threshold"]
    control, control_reason = phenotype(by_id.get(control_id), threshold)
    prediction, reason = phenotype(by_id.get(knockout_id), threshold)
    ko = by_id.get(knockout_id, {})
    label = case.get("label", "unknown")
    return {
        "case_id": case["case_id"], "locus_tag": case.get("locus_tag"),
        "substrate": case.get("substrate"), "observed": label,
        "prediction": prediction, "prediction_basis": reason,
        "control_prediction": control, "control_basis": control_reason,
        "control_supports_deletion_comparison": control == "growth",
        "binary_agreement": prediction == label if label in ("growth", "no_growth") and prediction != "unknown" else None,
        "knockout_growth": ko.get("measurements", {}).get("max_growth"),
        "knocked_out_reactions": ko.get("measurements", {}).get("knocked_out_reactions", []),
        "raw_od_audited": case.get("raw_od_reproduces_published_metrics") is True,
        "borderline": case.get("borderline_published"),
        "solver_calls": result["solver"]["used"],
        "model_sha256": result["hashes"]["model_sha256"],
        "manifest_sha256": result["hashes"]["manifest_sha256"],
    }


def summarize(rows: list[dict]) -> dict:
    cells = Counter((r["observed"], r["prediction"]) for r in rows)
    agreed = [r for r in rows if r["binary_agreement"] is not None]
    predicted_aux = [r for r in rows if r["prediction"] == "no_growth"]
    return {
        "n_cases": len(rows),
        "n_executed": sum(r.get("execution_status") == "completed" for r in rows),
        "n_unknown_predictions": sum(r["prediction"] == "unknown" for r in rows),
        "n_wild_type_growing": sum(r["control_prediction"] == "growth" for r in rows),
        "n_binary_evaluable": len(agreed),
        "n_binary_agree": sum(r["binary_agreement"] for r in agreed),
        "binary_agreement_among_evaluable": sum(r["binary_agreement"] for r in agreed) / len(agreed) if agreed else None,
        "n_predicted_auxotrophs": len(predicted_aux),
        "n_confirmed_auxotrophs": sum(r["observed"] == "no_growth" for r in predicted_aux),
        # Partial growth remains in this denominator; it is not confirmed no growth.
        "confirmed_fraction_of_predicted_auxotrophs": sum(r["observed"] == "no_growth" for r in predicted_aux) / len(predicted_aux) if predicted_aux else None,
        "observed_by_predicted_counts": {obs: {pred: cells[obs, pred] for pred in ("growth", "no_growth", "unknown")}
                                          for obs in ("growth", "partial", "no_growth", "unknown")},
    }


def run_benchmark(cases_path: Path, out: Path, model_override: Path | None = None) -> dict:
    records = [json.loads(line) for line in cases_path.read_text().splitlines() if line.strip()]
    ids = [r["case"]["case_id"] for r in records]
    if not records or len(ids) != len(set(ids)):
        raise ValueError("case identifiers must be nonempty and unique")
    if out.exists():
        raise ValueError(f"refusing existing run directory: {out}")
    out.mkdir(parents=True)
    model_cache: dict[tuple, Any] = {}
    rows = []
    for i, record in enumerate(records):
        case, manifest = record["case"], record["manifest"]
        result_path = out / f"case-{i:04d}.json"
        try:
            model_path = model_override or Path(manifest["model"])
            if not model_path.is_absolute():
                model_path = cases_path.resolve().parent / model_path
            cfg = ME.resolve_config(manifest)
            key = (str(model_path), cfg["default_bound"])
            if key not in model_cache:
                model_cache[key] = ME.load_sbml_model(model_path, cfg["default_bound"])
            model = model_cache[key]
            expected = manifest.get("model_sha256_expected")
            if not expected or model.sha256 != expected:
                raise ValueError("model hash missing or different from registered case")
            menu, _ = ME.build_menu(model, manifest)
            selected = [m for m in menu if m["id"] in cfg.get("panel", [])]
            controls = [m for m in selected if m["kind"] == "carbon_source"]
            knockouts = [m for m in selected if m["kind"] == "knockout"]
            if cfg["strategy"] != "fixed" or len(controls) != 1 or len(knockouts) != 1:
                raise ValueError("phenotype cases require one explicit matched carbon control and one knockout")
            if knockouts[0].get("unmapped_candidate_genes"):
                raise ValueError("knockout gene does not map to model")
            target = model.gene_alias.get(case.get("locus_tag"))
            deleted = {g for mod in knockouts[0]["mods"] if mod["op"] == "knockout" for g in mod["genes"]}
            if target is None or deleted != {target}:
                raise ValueError("registered locus and simulated gene deletion differ")
            ko_conditions = [m for m in knockouts[0]["mods"] if m["op"] != "knockout"]
            if controls[0]["mods"] != ko_conditions:
                raise ValueError("wild-type control and knockout conditions differ")
            # No measured label, OD value or case metadata is passed to inference.
            result = ME.run_experiments(manifest, model)
            result_path.write_text(json.dumps(result, indent=2) + "\n")
            row = score_case(case, result, controls[0]["id"], knockouts[0]["id"])
            row.update(execution_status="completed", result_file=result_path.name,
                       result_sha256=sha256(result_path))
        except Exception as exc:
            # Keep every registered case in the denominator, including errors.
            row = {"case_id": case["case_id"], "locus_tag": case.get("locus_tag"),
                   "observed": case.get("label", "unknown"), "prediction": "unknown",
                   "control_prediction": "unknown", "binary_agreement": None,
                   "raw_od_audited": case.get("raw_od_reproduces_published_metrics") is True,
                   "execution_status": "failed", "error": f"{type(exc).__name__}: {exc}"}
        rows.append(row)
        print(f"{i + 1}/{len(records)} {row['case_id']}: {row['execution_status']}", flush=True)
    report = {
        "schema": "gspa.metabolic_phenotype_benchmark/1",
        "role": "measured knockout-growth agreement; not independent function recovery",
        "selection": "cases selected by published model-based auxotrophy design; no general recall or specificity claim",
        "inputs": {"cases_sha256": sha256(cases_path), "runner_sha256": sha256(Path(__file__)),
                   "engine_sha256": sha256(Path(ME.__file__))},
        "all_cases": summarize(rows),
        "sensitivity_raw_od_reproduced": summarize([r for r in rows if r["raw_od_audited"]]),
        "sensitivity_growing_wild_type": summarize([r for r in rows if r["control_prediction"] == "growth"]),
        "sensitivity_optimal_knockout_only": summarize([r for r in rows if r.get("prediction_basis") == "optimal"]),
        "cases": rows,
        "status": ("completed_with_failures" if any(r["execution_status"] != "completed" for r in rows)
                   else "completed_with_unassessed_cases" if any(r["prediction"] == "unknown" for r in rows)
                   else "completed"),
    }
    (out / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cases", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--model", type=Path, help="relocated model; must match registered SHA-256")
    args = parser.parse_args()
    report = run_benchmark(args.cases, args.out, args.model)
    print(json.dumps(report["all_cases"], indent=2))
    return 0 if report["status"] == "completed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
