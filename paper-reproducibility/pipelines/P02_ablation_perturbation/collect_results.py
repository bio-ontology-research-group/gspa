#!/usr/bin/env python3
"""Build the two P02 perturbation result tables from the per-run summaries.

  results/pillar1_2_perturbation_designs.tsv  all runs (design, jar, claims_ids columns)
  results/pillar1_2_perturbation_summary.tsv  design v2 + current jar + normalised claims IDs,
      historical column layout plus a leading `design` column

Run summaries written before the claims-ID fix have no claims_ids column; they are
labelled claims_ids=raw.

F-max values are rounded to 4 decimals (as in the historical summary).

Usage:
  collect_results.py --runs DIR/v1_old_jar DIR/v1_current_jar DIR/v2_current_jar \
      DIR/v1_current_jar_idfix DIR/v2_current_jar_idfix --results-dir results/
"""
import argparse
import csv
import os

FMAX = ["fmax_micro_all", "fmax_micro_exp", "fmax_cafa_all", "fmax_cafa_exp"]
DESIGNS_COLS = ["design", "jar", "claims_ids", "genome", "config", "injected_pairs", "injected_claims",
                "surviving_violations", "surviving_anatomical", "surviving_taxon",
                "suppression_rate_pct"] + FMAX
SUMMARY_COLS = ["design", "genome", "config", "injected_false_positives", "surviving_violations",
                "surviving_anatomical", "surviving_taxon", "suppression_rate_pct",
                "fmax_micro_exp", "fmax_cafa_exp", "fmax_micro_all", "fmax_cafa_all"]
CURRENT_JAR = "gspa-unified.jar"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", nargs="+", required=True, help="run directories containing perturbation_summary.tsv")
    ap.add_argument("--results-dir", required=True)
    a = ap.parse_args()

    rows = []
    for d in a.runs:
        with open(os.path.join(d, "perturbation_summary.tsv")) as f:
            for r in csv.DictReader(f, delimiter="\t"):
                for k in FMAX:
                    r[k] = "" if r[k] == "" else f"{round(float(r[k]), 4):g}"
                r["claims_ids"] = r.get("claims_ids") or "raw"
                rows.append(r)

    with open(os.path.join(a.results_dir, "pillar1_2_perturbation_designs.tsv"), "w") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=DESIGNS_COLS, extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    with open(os.path.join(a.results_dir, "pillar1_2_perturbation_summary.tsv"), "w") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=SUMMARY_COLS, extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        for r in rows:
            if r["design"] == "v2" and r["jar"] == CURRENT_JAR and r["claims_ids"] == "normalised":
                w.writerow(dict(r, injected_false_positives=r["injected_pairs"]))


if __name__ == "__main__":
    main()
