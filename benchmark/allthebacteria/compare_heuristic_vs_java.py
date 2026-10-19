#!/usr/bin/env python3
"""Compare the Stage C heuristic (delta_events -> insight3_pathway_and_essential_completeness.tsv)
with GSPA JVM output (insight3_per_genome_java*.tsv) on exactly the same genomes.

Per-genome totals are summed over the genomes present in both tables (inner join on genome_id):
  heuristic rescue events  = n_B6_operon + n_B5_gapfill + n_B2_coherence + n_B1_essential
  heuristic suppressions   = n_B2_singleton + n_B3_taxon + n_B4_anatomical
  Java rescued pairs       = n_rescued_pairs (p0 < 0.5 <= p1)
  Java suppressed pairs    = n_suppressed_pairs (p0 >= 0.5 > p1 or removed)
Note: heuristic events include sub-threshold shifts (delta >= 0.08), so this is an upper bound
on what the heuristic would call a threshold-crossing rescue.
"""
import argparse
import csv
from pathlib import Path


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--heuristic", type=Path, required=True)
    ap.add_argument("--java", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    H = {r["genome_id"]: r for r in csv.DictReader(open(a.heuristic), delimiter="\t")}
    J = {r["genome_id"]: r for r in csv.DictReader(open(a.java), delimiter="\t")}
    common = sorted(set(H) & set(J))
    i = lambda r, k: int(float(r[k] or 0))
    rows = [
        ("n_genomes_in_both", len(common), len(common)),
        ("n_orfs", sum(i(H[g], "n_orfs") for g in common), sum(i(J[g], "n_orfs") for g in common)),
        ("n_dark_lt30", sum(i(H[g], "n_dark_orfs_lt30") for g in common), sum(i(J[g], "n_dark_lt30") for g in common)),
        ("n_twilight_30_50", sum(i(H[g], "n_twilight_orfs_30_50") for g in common), sum(i(J[g], "n_twilight_30_50") for g in common)),
        ("rescue_events_or_pairs", sum(i(H[g], k) for g in common for k in (
            "n_B6_operon_rescues", "n_B5_gapfill_rescues", "n_B2_coherence_rescues", "n_B1_essential_rescues")),
         sum(i(J[g], "n_rescued_pairs") for g in common)),
        ("operon_B6_vs_genomic_context", sum(i(H[g], "n_B6_operon_rescues") for g in common),
         sum(i(J[g], "rescued_pairs__genomic_context") for g in common)),
        ("gapfill_B5_vs_gap_filling", sum(i(H[g], "n_B5_gapfill_rescues") for g in common),
         sum(i(J[g], "rescued_pairs__gap_filling") for g in common)),
        ("coherence_B2_vs_coherence", sum(i(H[g], "n_B2_coherence_rescues") for g in common),
         sum(i(J[g], "rescued_pairs__coherence") for g in common)),
        ("essential_B1_vs_essentiality", sum(i(H[g], "n_B1_essential_rescues") for g in common),
         sum(i(J[g], "rescued_pairs__essentiality") for g in common)),
        ("suppressions", sum(i(H[g], k) for g in common for k in (
            "n_B2_singleton_suppressed", "n_B3_taxon_suppressed", "n_B4_anatomical_suppressed")),
         sum(i(J[g], "n_suppressed_pairs") for g in common)),
        ("singleton_B2_vs_singleton_prune", sum(i(H[g], "n_B2_singleton_suppressed") for g in common),
         sum(i(J[g], "suppressed__complex_singleton_prune") for g in common)),
        ("taxon_B3_vs_taxon_forbidden+sat", sum(i(H[g], "n_B3_taxon_suppressed") for g in common),
         sum(i(J[g], k) for g in common for k in ("suppressed__taxon_forbidden", "suppressed__sat_coannotation_conflict",
                                                 "suppressed__stage2_removed_coannotation_conflict"))),
        ("anatomical_B4_vs_anatomy", sum(i(H[g], "n_B4_anatomical_suppressed") for g in common),
         sum(i(J[g], k) for g in common for k in ("suppressed__anatomy_forbidden", "suppressed__anatomical_prior_other"))),
    ]
    with open(a.out, "w") as fo:
        fo.write("quantity\theuristic_delta_events\tgspa_jvm\n")
        for q, h, j in rows:
            fo.write(f"{q}\t{h}\t{j}\n")
    print(open(a.out).read())


if __name__ == "__main__":
    main()
