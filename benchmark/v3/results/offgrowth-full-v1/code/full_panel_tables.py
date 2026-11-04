#!/usr/bin/env python3
"""Secondary tables for the full-panel off-growth run (tasks/v3-offgrowth-full-panel-2026-11-02.md, §2.7).

Reads the aggregate ``summary.json`` (the runner's own output; evaluations already re-labelled with every
experimental label) and writes descriptive tables. Recovery intervals use the runner's own genome-block bootstrap
(``score.summarize_arms_with_bootstrap``, 2,000 resamples, seed 20261030) on each subset.

* ``subsets.tsv``          recovery per arm for all cases and per dataset / domain / case source (off-growth stratum)
* ``per_genome.tsv``       cases, recoveries per arm, union calls / harm / cost per genome
* ``distinct_decisions.tsv`` labelled gene calls counted once per (genome, protein, reaction) instead of once per case
* ``candidate_rank.tsv``   for each hidden case: is the hidden protein the top candidate for one of its ECs (the
                           ceiling of the top-candidate attribution rule), from the unmasked candidate tables
* ``pilot_comparison.tsv`` the 22 Price genomes: the pilot (six-predictor union) vs this run (DIAMOND baseline), paired
                           by case, exact two-sided McNemar (descriptive)
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

ARMS = ["A0", "A1", "A2", "A3", "A1_A2_A3", "A4@10%", "A4@25%"]


def load_json(path):
    p = Path(path)
    with (gzip.open(p, "rt") if p.suffix == ".gz" else open(p)) as fh:
        return json.load(fh)


def mcnemar_two_sided(b, c):
    n = b + c
    if n == 0:
        return 1.0
    k = min(b, c)
    p = sum(math.comb(n, i) for i in range(k + 1)) / 2.0 ** n
    return min(1.0, 2 * p)


def main(argv=None) -> int:
    from offgrowth.runner import _reconstruct_evaluation
    from offgrowth.score import summarize_arms_with_bootstrap, BOOTSTRAP_SEED
    from offgrowth.panel import build_ec_to_reactions_map, load_seed_reactions

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--summary", required=True)
    ap.add_argument("--genomes", required=True)
    ap.add_argument("--candidates-dir", required=True)
    ap.add_argument("--seed-reactions", required=True)
    ap.add_argument("--ec-aliases", required=True)
    ap.add_argument("--pilot-summary", required=True)
    ap.add_argument("--n-resamples", type=int, default=2000)
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args(argv)
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    genomes = {r["genome_id"]: r for r in csv.DictReader(open(a.genomes), delimiter="\t")}
    summ = load_json(a.summary)
    evals = {arm: [_reconstruct_evaluation(e) for e in summ["evaluations_by_arm"][arm]] for arm in ARMS}

    # ---- subsets
    def subset(pred):
        return {arm: [e for e in evals[arm] if pred(e)] for arm in ARMS}
    groups = [("all", lambda e: True)]
    for ds in ("price_matrix", "recurrence_cohort", "recurrence_peer"):
        groups.append((f"dataset={ds}", lambda e, ds=ds: genomes[e.genome_id]["dataset"] == ds))
    for dom in ("Bacteria", "Archaea"):
        groups.append((f"domain={dom}", lambda e, dom=dom: genomes[e.genome_id]["domain"] == dom))
    for src in ("price148", "swissprot"):
        groups.append((f"source={src}", lambda e, src=src: e.case_source == src))
    with open(out / "subsets.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["subset", "arm", "genomes", "cases", "off_growth_cases", "reaction_readded_off_growth",
                    "recovered_off_growth", "recovery_off_growth", "ci_low", "ci_high", "gene_calls", "supported",
                    "harmful", "unconfirmed", "harm_rate", "harm_ci_low", "harm_ci_high"])
        for name, pred in groups:
            sub = subset(pred)
            if not sub["A0"]:
                continue
            sm = summarize_arms_with_bootstrap(sub, n_resamples=a.n_resamples, seed=BOOTSTRAP_SEED)
            for arm in ARMS:
                s = sm[arm]
                off = [e for e in sub[arm] if e.stratum == "off_growth"]
                w.writerow([name, arm, len({e.genome_id for e in sub[arm]}), s.n_cases_total, s.n_cases_off_growth,
                            sum(1 for e in off if e.reaction_recovered), s.n_recovered_off_growth,
                            f"{s.recovery_rate_off_growth:.3f}", f"{s.recovery_ci_off_growth[0]:.3f}",
                            f"{s.recovery_ci_off_growth[1]:.3f}", s.n_gene_calls, s.n_supported_calls, s.n_harmful_calls,
                            s.n_unconfirmed_calls, f"{s.harm_rate:.3f}", f"{s.harm_ci[0]:.3f}", f"{s.harm_ci[1]:.3f}"])

    # ---- per genome
    with open(out / "per_genome.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["genome_id", "dataset", "domain", "organism", "cases", "off_growth_cases"] +
                   [f"recovered_{a_}" for a_ in ARMS] +
                   ["union_reaction_readded", "union_gene_calls", "union_supported", "union_harmful", "union_unconfirmed",
                    "union_lps", "union_cpu_hours", "a4_25_lps", "a4_25_cpu_hours"])
        by_g = defaultdict(lambda: defaultdict(list))
        for arm in ARMS:
            for e in evals[arm]:
                by_g[e.genome_id][arm].append(e)
        for gid in sorted(by_g, key=lambda g: (genomes[g]["dataset"], g)):
            d = by_g[gid]
            u, b = d["A1_A2_A3"], d["A4@25%"]
            w.writerow([gid, genomes[gid]["dataset"], genomes[gid]["domain"], genomes[gid]["organism"], len(d["A0"]),
                        sum(1 for e in d["A0"] if e.stratum == "off_growth")] +
                       [sum(1 for e in d[a_] if e.recovered) for a_ in ARMS] +
                       [sum(1 for e in u if e.reaction_recovered), sum(e.n_gene_calls for e in u),
                        sum(e.n_supported_calls for e in u), sum(e.n_harmful_calls for e in u),
                        sum(e.n_unconfirmed_calls for e in u), sum(e.n_lps for e in u),
                        f"{sum(e.cpu_seconds for e in u) / 3600:.2f}", sum(e.n_lps for e in b),
                        f"{sum(e.cpu_seconds for e in b) / 3600:.2f}"])

    # ---- distinct labelled decisions
    with open(out / "distinct_decisions.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["arm", "distinct_labelled_decisions", "supported", "harmful", "harm_rate", "distinct_unconfirmed_decisions"])
        for arm in ARMS:
            kinds = {}
            unconf = set()
            for e in evals[arm]:
                for c in e.calls:
                    if not c.protein_id:
                        continue
                    key = (e.genome_id, c.protein_id, c.reaction_id)
                    if c.call_kind in ("supported", "harmful"):
                        kinds.setdefault(key, set()).add(c.call_kind)
                    elif c.call_kind == "unconfirmed":
                        unconf.add(key)
            sup = sum(1 for v in kinds.values() if "supported" in v)
            harm = sum(1 for v in kinds.values() if v == {"harmful"})
            n = sup + harm
            w.writerow([arm, n, sup, harm, f"{harm / n:.3f}" if n else "NA", len(unconf)])

    # ---- candidate rank ceiling
    reactions = load_seed_reactions(Path(a.seed_reactions), Path(a.ec_aliases))
    ec_to_rxn = build_ec_to_reactions_map(reactions, usable_only=False)
    cand_cache = {}

    def cands(gid):
        if gid not in cand_cache:
            per_ec = defaultdict(dict)
            p = Path(a.candidates_dir) / f"{gid}.candidates.tsv"
            for r in csv.DictReader(open(p), delimiter="\t"):
                per_ec[r["ec"]][r["protein_id"]] = float(r["score"])
            cand_cache[gid] = per_ec
        return cand_cache[gid]
    rank_rows = []
    rec_union = {e.case_id: e.recovered for e in evals["A1_A2_A3"]}
    for e in evals["A0"]:
        per_ec = cands(e.genome_id)
        has, top = False, False
        for ec in e.truth_ecs:
            m = per_ec.get(ec, {})
            if e.protein_id in m:
                has = True
                if m[e.protein_id] >= max(m.values()):
                    top = True
        rank_rows.append([e.case_id, e.genome_id, e.protein_id, ";".join(e.truth_ecs), e.case_source, e.stratum,
                          int(has), int(top), int(bool(any(ec_to_rxn.get(ec) for ec in e.truth_ecs))),
                          int(rec_union.get(e.case_id, False))])
    with open(out / "candidate_rank.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["case_id", "genome_id", "protein_id", "truth_ecs", "case_source", "stratum", "has_candidate_for_truth_ec",
                    "top_candidate_for_truth_ec", "truth_ec_maps_to_seed_reaction", "recovered_union"])
        w.writerows(rank_rows)

    # ---- pilot comparison
    pilot = load_json(a.pilot_summary)
    with open(out / "pilot_comparison.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["arm", "paired_cases", "pilot_recovered", "full_recovered", "both", "pilot_only", "full_only",
                    "mcnemar_two_sided_p", "pilot_reaction_readded", "full_reaction_readded"])
        for arm in ARMS:
            pv = {e["case_id"]: e for e in pilot["evaluations_by_arm"][arm]}
            fv = {e.case_id: e for e in evals[arm] if e.case_id in pv}
            ids = sorted(fv)
            pr = {i for i in ids if pv[i]["recovered"]}
            fr = {i for i in ids if fv[i].recovered}
            b, c = len(pr - fr), len(fr - pr)
            w.writerow([arm, len(ids), len(pr), len(fr), len(pr & fr), b, c, f"{mcnemar_two_sided(b, c):.3g}",
                        sum(1 for i in ids if pv[i]["reaction_recovered"]), sum(1 for i in ids if fv[i].reaction_recovered)])
    print("tables written to", out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
