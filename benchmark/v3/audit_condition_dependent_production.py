#!/usr/bin/env python3
"""Audit a fetched primary condition-dependent bioactivity dataset (Westhoff et al., Dryad 10.5061/dryad.wpzgmsbkw).

Purpose: decide, from data and not from assumption, what a condition-dependent claim about secondary
metabolism could rest on. MIBiG (positive characterised-production records, no medium) cannot support an empirical
medium-effect claim. This audit records what the only fetched condition-dependent source does and does not contain.

Nothing here is a model prediction. No negative is invented: 'no inhibition' is reported as a bioactivity observation,
never as 'does not produce'.
"""
from __future__ import annotations

import argparse
import json
import statistics
from pathlib import Path

import openpyxl

SCHEMA = "gspa.v3.condition_dependent_production_audit/1"


def sheet_rows(path: Path) -> dict[str, list[list]]:
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    return {ws.title: [list(r) for r in ws.iter_rows(values_only=True)] for ws in wb.worksheets}


def numeric_cells(rows: list[list]) -> list[float]:
    return [float(c) for r in rows for c in r if isinstance(c, (int, float)) and not isinstance(c, bool)]


def summarize_matrix(path: Path) -> dict:
    out = {}
    for name, rows in sheet_rows(path).items():
        vals = numeric_cells(rows)
        out[name] = {"rows": len(rows), "max_cols": max((len(r) for r in rows), default=0),
                     "header_preview": [c for c in (rows[0] if rows else [])][:6],
                     "n_numeric": len(vals), "n_zero": sum(v == 0 for v in vals), "n_positive": sum(v > 0 for v in vals),
                     "n_negative": sum(v < 0 for v in vals),
                     "median_positive": statistics.median([v for v in vals if v > 0]) if any(v > 0 for v in vals) else None}
    return out


def audit(data_dir: Path) -> dict:
    manifest = json.loads((data_dir / "MANIFEST.json").read_text())
    files = {f["key"]: f for f in manifest["files"]}
    report = {"schema": SCHEMA, "source": {k: manifest[k] for k in ("record", "doi", "title", "license", "publication_date", "creators", "retrieved_utc")},
              "files": {k: {"size": v["size"], "sha256": v.get("sha256")} for k, v in files.items()},
              "matrices": {}}
    for key in ("Inhibition_Data_-_Low_Resource_Conditions.xlsx", "Inhibition_Data_-_High_Resource_Conditions.xlsx",
                "Induction_Suppression_-_Low_Resource_Conditions.xlsx", "Induction_Suppression_-_High_Resource_Conditions.xlsx",
                "Heatmap_biolog_results.xlsx", "BGC_Families.xlsx"):
        report["matrices"][key] = summarize_matrix(data_dir / key)
    report["assessment"] = {
        "condition_dependent": "yes: two resource regimes (low / high) measured for the same Streptomyces strain panel",
        "measured_quantity": "pairwise growth-inhibition (bioactivity) between strains, not a titre of a named compound",
        "compound_identity": "no: inhibition cannot be assigned to a metabolite or a BGC product without further chemistry",
        "negatives": "no-inhibition entries exist but are not evidence of non-production (resistance, diffusion, detection limit)",
        "metabolic_models": "none supplied; the strains are soil isolates without curated genome-scale models",
        "usable_for": ["resource-dependence of antibiotic bioactivity as a descriptive measurement",
                       "a held-out descriptive test that BGC content (BGC_Families) alone does not predict resource dependence"],
        "not_usable_for": ["empirical validation of a flux-based medium effect on a named secondary metabolite",
                           "precision or specificity of a production predictor",
                           "any claim that GSPA or a planner recovers production conditions"],
    }
    return report


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data", type=Path, required=True, help="directory holding the fetched files and MANIFEST.json")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    report = audit(args.data)
    args.out.write_text(json.dumps(report, indent=1, default=str) + "\n")
    print(json.dumps(report["assessment"], indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
