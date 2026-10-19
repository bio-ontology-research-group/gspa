#!/usr/bin/env python3
"""Stage 3 trigger rates over the 88 multikingdom proteomes (no LLM).

The '1.84% trigger rate' quoted for Stage 3 had no producing code. This script
regenerates it from the real anomaly extractor:

  manifest   Build one integrate run per proteome from the two batch manifests that produced
             the 88-proteome results (batch_test/manifest_45.tsv rows for the nine bench9
             genomes, multikingdom_batch/manifest_multikingdom.tsv rows for the other 79), using
             the `iter3_relax` configuration of each row (taxon, essential profile, operons,
             gaps, max_iter=3, Stage 2 consistency enforcement in relax mode, all priors).
             The runs go through integrate_one.sh with P10_DM=0 (the batch configuration has
             no dark-matter suggester) and the cap-lifting driver; the production output
             (AnomalyExtractor cap = 150) is the first 150 records of the uncapped output.
  aggregate  Per proteome: n_proteins (distinct protein ids with at least one claim),
             anomalies (capped and uncapped) per category, n_triggered = distinct proteins
             referenced by the capped anomalies (candidates + complex anchors), pct =
             100 * n_triggered / n_proteins. Writes results/pillar4_v2_trigger_rates.tsv with a
             per-domain and an overall summary row (pooled and median).
"""

from __future__ import annotations

import argparse
import collections
import csv
import json
import statistics
import sys
from pathlib import Path

CAP = 150
CATS = ["PATHWAY_GAP", "INCOMPLETE_COMPLEX", "BORDERLINE_ESSENTIAL", "BORDERLINE_CONFLICT"]
BENCH9 = ["ecoli", "hpylori", "mgenitalium", "mjannaschii", "ecolo157", "bsubtilis", "mtb", "synechocystis", "paeruginosa"]


def cmd_manifest(a: argparse.Namespace) -> None:
    domains = {}
    with open(a.per_genome_eval, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            domains[r["genome_id"]] = r["domain"]
    rows = []
    with open(a.manifest45, encoding="utf-8") as fh:
        for r in csv.reader(fh, delimiter="\t"):
            if r and r[0].endswith("_iter3_relax") and r[0].rsplit("_iter3_relax", 1)[0] in BENCH9:
                rows.append(("b9_" + r[0].rsplit("_iter3_relax", 1)[0], r))
    with open(a.manifest_mk, encoding="utf-8") as fh:
        for r in csv.reader(fh, delimiter="\t"):
            if r and r[0].endswith("__iter3_relax"):
                rows.append((r[0].rsplit("__iter3_relax", 1)[0], r))
    pre = a.path_prefix or ""
    fix = lambda p: p if p in ("-", "") else pre + p  # noqa: E731
    with open(a.runlist, "w", encoding="utf-8") as out, open(a.genome_table, "w", encoding="utf-8") as gt:
        out.write("run_id\ttaxon\tprofile\tclaims\toperons\tgaps\toutdir\n")
        gt.write("genome\tdomain\ttaxon\tessential_profile\tclaims\toperons\tgaps\tmax_iter\tenforce_stage2\tenforce_mode\tenable_priors\n")
        for gid, r in rows:
            assert r[7] == "3" and r[8] == "true" and r[9] == "relax" and r[10] in ("-", ""), r
            out.write("\t".join([f"trig_{gid}", r[3], r[4], fix(r[1]), fix(r[5]), fix(r[6]), f"{a.out_root}/{gid}"]) + "\n")
            gt.write("\t".join([gid, domains.get(gid, "?"), r[3], r[4], r[1], r[5], r[6], r[7], r[8], r[9], r[10]]) + "\n")
    print(f"[manifest] {len(rows)} proteomes", file=sys.stderr)


def cmd_aggregate(a: argparse.Namespace) -> None:
    rows = []
    with open(a.genome_table, encoding="utf-8") as fh:
        genomes = list(csv.DictReader(fh, delimiter="\t"))
    pre = a.path_prefix or ""
    for g in genomes:
        d = a.out_root / g["genome"]
        st = (d / "run_status.tsv").read_text().split("\t") if (d / "run_status.tsv").exists() else ["", "NA"]
        anoms = [json.loads(l) for l in open(d / "anomalies.jsonl", encoding="utf-8")] if (d / "anomalies.jsonl").exists() else None
        prots = set()
        with open(pre + g["claims"], encoding="utf-8") as fh:
            for line in fh:
                if line.strip():
                    prots.add(json.loads(line)["protein_id"])
        row = {"genome": g["genome"], "domain": g["domain"], "taxon": g["taxon"], "n_proteins": len(prots),
               "integrate_rc": st[1].strip()}
        if anoms is None:
            row.update({"n_anomalies": "NA", "n_triggered": "NA", "pct": "NA"})
            rows.append(row)
            continue
        capped = anoms[:CAP]
        trig = set()
        for r in capped:
            trig.update(c["protein_id"] for c in r["candidates"])
            if r.get("context", {}).get("anchor_protein_id"):
                trig.add(r["context"]["anchor_protein_id"])
        cc = collections.Counter(r["anomaly_type"] for r in capped)
        cu = collections.Counter(r["anomaly_type"] for r in anoms)
        trig_u = set()
        for r in anoms:
            trig_u.update(c["protein_id"] for c in r["candidates"])
        row.update({"n_anomalies": len(capped), "n_triggered": len(trig),
                    "pct": round(100.0 * len(trig) / len(prots), 3) if prots else 0.0,
                    **{f"n_{c.lower()}": cc.get(c, 0) for c in CATS},
                    "cap_reached": len(anoms) >= CAP,
                    "n_anomalies_uncapped": len(anoms), "n_triggered_uncapped": len(trig_u),
                    "pct_uncapped": round(100.0 * len(trig_u) / len(prots), 3) if prots else 0.0,
                    **{f"n_{c.lower()}_uncapped": cu.get(c, 0) for c in CATS}})
        rows.append(row)
    cols = ["genome", "domain", "taxon", "n_proteins", "n_triggered", "pct", "n_anomalies"] + \
           [f"n_{c.lower()}" for c in CATS] + ["cap_reached", "n_anomalies_uncapped", "n_triggered_uncapped", "pct_uncapped"] + \
           [f"n_{c.lower()}_uncapped" for c in CATS] + ["integrate_rc"]
    ok = [r for r in rows if r.get("n_triggered") not in ("NA", None)]
    summ = []
    for dom in sorted({r["domain"] for r in ok}) + ["ALL"]:
        rs = [r for r in ok if dom == "ALL" or r["domain"] == dom]
        if not rs:
            continue
        s = {"genome": f"SUMMARY_{dom}", "domain": dom, "taxon": f"n_genomes={len(rs)}",
             "n_proteins": sum(r["n_proteins"] for r in rs), "n_triggered": sum(r["n_triggered"] for r in rs),
             "n_anomalies": sum(r["n_anomalies"] for r in rs)}
        s["pct"] = f"pooled={100.0 * s['n_triggered'] / s['n_proteins']:.3f};median={statistics.median(r['pct'] for r in rs):.3f};" \
                   f"mean={statistics.mean(r['pct'] for r in rs):.3f}"
        for c in CATS:
            s[f"n_{c.lower()}"] = sum(r[f"n_{c.lower()}"] for r in rs)
            s[f"n_{c.lower()}_uncapped"] = sum(r[f"n_{c.lower()}_uncapped"] for r in rs)
        s["cap_reached"] = sum(bool(r["cap_reached"]) for r in rs)
        s["n_anomalies_uncapped"] = sum(r["n_anomalies_uncapped"] for r in rs)
        s["n_triggered_uncapped"] = sum(r["n_triggered_uncapped"] for r in rs)
        s["pct_uncapped"] = f"pooled={100.0 * s['n_triggered_uncapped'] / s['n_proteins']:.3f}"
        summ.append(s)
    with open(a.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows + summ:
            w.writerow({k: r.get(k, "") for k in cols})
    for s in summ:
        print(f"{s['genome']}: {s['taxon']} proteins={s['n_proteins']} triggered={s['n_triggered']} {s['pct']}", file=sys.stderr)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    m = sub.add_parser("manifest")
    m.add_argument("--manifest45", type=Path, required=True)
    m.add_argument("--manifest-mk", type=Path, required=True)
    m.add_argument("--per-genome-eval", type=Path, required=True)
    m.add_argument("--out-root", required=True)
    m.add_argument("--path-prefix", default="")
    m.add_argument("--runlist", type=Path, required=True)
    m.add_argument("--genome-table", type=Path, required=True)
    g = sub.add_parser("aggregate")
    g.add_argument("--genome-table", type=Path, required=True)
    g.add_argument("--out-root", type=Path, required=True)
    g.add_argument("--path-prefix", default="")
    g.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    {"manifest": cmd_manifest, "aggregate": cmd_aggregate}[a.cmd](a)
    return 0


if __name__ == "__main__":
    sys.exit(main())
