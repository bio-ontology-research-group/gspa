#!/usr/bin/env python3
"""Result tables from the fetched per-genome table (results/eq_mags_per_bin.tsv).

Writes results/eq_mags_by_dataset.tsv (totals and medians per group: soil bins by CheckM2 class and
together, the three culture datasets and together, everything) and results/eq_mags_phyla.tsv (GTDB
phylum of each genome's taxon per dataset). Group membership: a soil bin belongs to soil_bin and to its
class (HQ, MQ); a culture genome to its dataset and to cultures; every genome to all. Medians skip
genomes whose value is NA (pathway completeness is NA when no predicted pathway could be scored).
"""
import argparse, csv, statistics
from collections import Counter, defaultdict

GROUPS = ["soil_HQ", "soil_MQ", "soil_bin", "culture_enrichment", "culture_isolate", "culture_thermophilic", "cultures", "all"]
SUM_COLS = ["n_proteins", "n_claims", "stage0_proteins", "stage12_proteins", "score_n_annot_stage0_q050",
            "score_n_annot_stage12_q050", "score_n_gap_rxns", "score_n_gap_fillable", "score_new_fills_q050",
            "score_promoted_fills_q050", "score_gaps_resolved_q050", "rlgc_suggestions", "operon_chains"]
MED_COLS = ["n_proteins", "n_claims", "stage0_proteins", "stage12_proteins", "stage2_actions", "rlgc_suggestions",
            "operon_chains", "score_n_gap_rxns", "score_n_gap_fillable", "score_new_fills_q050",
            "score_promoted_fills_q050", "score_gaps_resolved_q050", "score_pathway_completeness_before",
            "score_pathway_completeness_after", "score_gapseq_pathway_completeness",
            "score_essential_completeness_before", "score_essential_completeness_after"]


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def groups_of(r):
    if r["dataset"] == "soil_bin":
        return ["soil_bin", "soil_" + r["class"], "all"]
    return [r["dataset"], "cultures", "all"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--per-bin", default="results/eq_mags_per_bin.tsv")
    ap.add_argument("--out-dir", default="results")
    a = ap.parse_args()
    rows = list(csv.DictReader(open(a.per_bin), delimiter="\t"))
    by = defaultdict(list)
    for r in rows:
        for g in groups_of(r):
            by[g].append(r)
    out = []
    for g in GROUPS:
        rs = by[g]
        rec = {"group": g, "n_genomes": len(rs),
               "bakta": sum(r["protein_source"] == "bakta" for r in rs),
               "pyrodigal": sum(r["protein_source"] == "pyrodigal" for r in rs),
               "gapsmith_refill": sum(r["gapseq_source"] == "computed_gapsmith_refill" for r in rs),
               "gapsmith_fill_failed": sum(r["gapseq_source"] == "gapsmith_fill_failed" for r in rs),
               "archaea": sum(r["domain"] == "Archaea" for r in rs),
               "taxon_rank_genus_or_species": sum(r["taxon_rank"] in ("g", "s") for r in rs)}
        for c in SUM_COLS:
            rec["total_" + c.replace("score_", "")] = int(sum(num(r[c]) or 0 for r in rs))
        for c in MED_COLS:
            v = [num(r[c]) for r in rs if num(r[c]) is not None]
            rec["median_" + c.replace("score_", "")] = round(statistics.median(v), 4) if v else "NA"
            if c.startswith("score_pathway") or c.startswith("score_essential"):
                rec["n_" + c.replace("score_", "")] = len(v)
        for pre, post in (("score_pathway_completeness_before", "score_pathway_completeness_after"),
                          ("score_essential_completeness_before", "score_essential_completeness_after")):
            k = pre.replace("score_", "").replace("_before", "")
            rec["n_" + k + "_improved"] = sum(1 for r in rs if num(r[post]) is not None and num(r[post]) > num(r[pre]))
        out.append(rec)
    cols = list(out[0].keys())
    with open(f"{a.out_dir}/eq_mags_by_dataset.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(out)
    phyla = defaultdict(Counter)
    for r in rows:
        phyla[r["dataset"]][r["phylum"] or "NA"] += 1
        phyla["all"][r["phylum"] or "NA"] += 1
    names = sorted(phyla["all"], key=lambda p: (-phyla["all"][p], p))
    ds = ["soil_bin", "culture_enrichment", "culture_isolate", "culture_thermophilic", "all"]
    with open(f"{a.out_dir}/eq_mags_phyla.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["phylum"] + ds)
        for p in names:
            w.writerow([p] + [phyla[d][p] for d in ds])
    print(f"{len(rows)} genomes; {len(names)} phyla; wrote {a.out_dir}/eq_mags_by_dataset.tsv and eq_mags_phyla.tsv")


if __name__ == "__main__":
    main()
