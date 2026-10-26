#!/usr/bin/env python3
"""Independent COBRApy/GLPK replay of registered matched knockout experiments.

This checks SBML/GPR interpretation and flux results, not biological accuracy.
Measured phenotypes are never passed to COBRApy. All protocols/cases are retained.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path

import cobra
from run_metabolic_benchmark import execution_manifest, sha256


def simulate(model, manifest, knockout):
    hypotheses = manifest["context"]["hypotheses"]
    if len(hypotheses) != 1:
        raise ValueError("oracle requires one registered deletion hypothesis")
    hypothesis = hypotheses[0]
    substrate = hypothesis["condition"]["carbon_exchange"]
    medium = dict(manifest["baseline"]["medium"])
    for exchange in manifest["baseline"].get("carbon_exchanges", []):
        medium.pop(exchange, None)
    medium[substrate] = manifest["config"].get("uptake_rate", 10.0)
    with model:
        # Preserve the SBML reaction identifiers and flux directions. No use of
        # the production parser, stoichiometric matrix, GPR evaluator or LP solver.
        for reaction in model.reactions:
            if len(reaction.metabolites) == 1:
                coefficient = next(iter(reaction.metabolites.values()))
                rate = medium.get(reaction.id, 0.0)
                lower, upper = reaction.bounds
                if coefficient < 0:
                    lower = -rate if reaction.id in medium else max(lower, 0.0)
                else:
                    upper = rate if reaction.id in medium else min(upper, 0.0)
                if lower > upper:
                    return {"status": "infeasible", "growth": None, "deleted_reactions": []}
                reaction.bounds = lower, upper
        maintenance = manifest["objectives"]["maintenance"]
        if maintenance != "none":
            reaction = model.reactions.get_by_id(maintenance["reaction"])
            reaction.lower_bound = max(reaction.lower_bound, maintenance["lower_bound"])
        closed = set()
        for constraint in manifest.get("condition_constraints", []):
            if constraint["substrate_exchange"] in medium:
                closed.update(constraint["closed_reactions"])
        for rid in sorted(closed):
            reaction = model.reactions.get_by_id(rid)
            if reaction.lower_bound > 0 or reaction.upper_bound < 0:
                return {"status": "infeasible", "growth": None, "deleted_reactions": []}
            reaction.bounds = (0.0, 0.0)
        before = {r.id for r in model.reactions if r.bounds == (0.0, 0.0)}
        if knockout:
            gene_id = hypothesis["gene"]
            matches = [g for g in model.genes if g.id in (gene_id, "G_" + gene_id) or g.name == gene_id]
            if len(matches) != 1:
                raise ValueError(f"ambiguous or missing gene: {gene_id}")
            matches[0].knock_out()
        deleted = sorted(r.id for r in model.reactions if r.bounds == (0.0, 0.0) and r.id not in before)
        model.objective = model.reactions.get_by_id(manifest["objectives"]["biomass"])
        model.objective_direction = "max"
        growth = float(model.slim_optimize())
        return {"status": model.solver.status, "growth": growth if math.isfinite(growth) else None,
                "deleted_reactions": deleted}


def compare(cases: Path, results: Path, model_path: Path) -> dict:
    records = [json.loads(line) for line in cases.read_text().splitlines() if line.strip()]
    report = json.loads((results / "report.json").read_text())
    if sha256(cases) != report["inputs"]["cases_sha256"]:
        raise ValueError("cases differ from production run")
    expected = {r["case_id"]: r for r in report["cases"]}
    model = cobra.io.read_sbml_model(model_path, f_replace={})
    model.solver = "glpk"
    model_hash = sha256(model_path)
    rows = []
    for record in records:
        cid = record["case"]["case_id"]
        row = {"case_id": cid}
        try:
            manifest = execution_manifest(record)
            if model_hash != manifest["model_sha256_expected"]:
                raise ValueError("model hash differs from registered model")
            source = results / expected[cid]["result_file"]
            if sha256(source) != expected[cid]["result_sha256"]:
                raise ValueError("production result hash mismatch")
            production = json.loads(source.read_text())
            comparisons = []
            for knockout in (False, True):
                kind = "knockout" if knockout else "carbon_source"
                matches = [r for r in production["experiments"] if r["kind"] == kind]
                if len(matches) != 1:
                    raise ValueError("expected one production result per matched condition")
                prod = matches[0]
                ref = simulate(model, manifest, knockout)
                if prod["outcome"] == "infeasible_with_maintenance":
                    agrees = ref["status"] == "infeasible"
                else:
                    flux = prod.get("measurements", {}).get("max_growth")
                    agrees = (prod["status"] == "completed" and ref["status"] == "optimal" and
                              flux is not None and ref["growth"] is not None and
                              math.isclose(flux, ref["growth"], rel_tol=1e-6, abs_tol=1e-7))
                comparisons.append({"kind": kind, "oracle": ref, "production_outcome": prod["outcome"],
                                    "production_growth": prod.get("measurements", {}).get("max_growth"),
                                    "agrees": agrees})
            row.update(status="completed", comparisons=comparisons, agrees=all(c["agrees"] for c in comparisons))
        except Exception as exc:
            row.update(status="failed", agrees=False, error=f"{type(exc).__name__}: {exc}")
        rows.append(row)
    return {"schema": "gspa.cobra_replay/1", "cobra_version": cobra.__version__, "solver": "glpk",
            "interpretation": "Independent implementation agreement, not phenotype or function accuracy",
            "hashes": {"cases": sha256(cases), "production_report": sha256(results / "report.json"),
                       "model": model_hash, "oracle_script": sha256(Path(__file__))},
            "n_cases": len(rows), "n_agree": sum(r["agrees"] for r in rows), "cases": rows}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cases", required=True, type=Path)
    parser.add_argument("--results", required=True, type=Path)
    parser.add_argument("--model", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    if args.out.exists():
        raise ValueError("refusing to overwrite oracle report")
    result = compare(args.cases, args.results, args.model)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k != "cases"}, indent=2))
    return 0 if result["n_agree"] == result["n_cases"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
