#!/usr/bin/env python3
"""Aggregate the Empty Quarter bin annotation: one row per bin and a run summary.

  04_aggregate.py --bins lists/bins.tsv --work work --out-dir results

Per bin: dataset (soil bins, culture genomes), class, taxonomy, protein source and count, claims, annotated proteins and terms after
Stage 0 and after Stage 1+2A+2B, Stage 2 actions, RLGC suggestions, gapseq availability, the P14
score columns when present, and wall times of every step. Summary: counts and medians by class, by dataset and overall.
"""
import argparse
import csv
import glob
import gzip
import os
import statistics


def count_integrated(path):
    """(proteins with >= 1 term, rows) of a GSPA integrate TSV (gzipped or not)."""
    if not os.path.exists(path):
        return "", ""
    prots, rows = set(), 0
    with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        try:
            pi = next(i for i, c in enumerate(hdr) if c in ("protein_id", "protein", "query", "gene"))
        except StopIteration:
            pi = 0
        for line in f:
            rows += 1
            prots.add(line.split("\t")[pi])
    return len(prots), rows


def read_kv_times(path):
    out = {}
    if os.path.exists(path):
        for line in open(path):
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2:
                out[f[0]] = f[1]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bins", required=True)
    ap.add_argument("--work", required=True)
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    bins = list(csv.DictReader(open(a.bins), delimiter="\t"))
    rows = []
    score_cols = []
    for b in bins:
        g = b["mag_id"]
        d = os.path.join(a.work, g)
        r = {k: b.get(k, "") for k in ("mag_id", "dataset", "sample", "bin", "class", "completeness", "contamination",
                                      "domain", "phylum", "genus", "taxon", "taxon_label", "taxon_rank", "gram",
                                      "essential_profile")}
        r["proteins_done"] = int(os.path.exists(os.path.join(d, "DONE_proteins")))
        r["gapseq_done"] = int(os.path.exists(os.path.join(d, "DONE_gapseq")))
        r["integrated"] = int(os.path.exists(os.path.join(d, "DONE")))
        src = os.path.join(d, "annot", "SOURCE")
        r["protein_source"] = open(src).read().strip() if os.path.exists(src) else ""
        dp = os.path.join(d, "DONE_proteins")
        if os.path.exists(dp):
            r["n_proteins"] = open(dp).read().rstrip("\n").split("\t")[2]
        dn = os.path.join(d, "DONE")
        if os.path.exists(dn):
            r["n_claims"] = open(dn).read().rstrip("\n").split("\t")[2]
        dg = os.path.join(d, "DONE_gapseq")
        r["gapseq_source"] = open(dg).read().split("\t")[1] if os.path.exists(dg) else ""
        r["stage0_proteins"], r["stage0_rows"] = count_integrated(os.path.join(d, f"{g}_stage0.tsv.gz"))
        r["stage12_proteins"], r["stage12_rows"] = count_integrated(os.path.join(d, f"{g}_stage12.tsv.gz"))
        act = os.path.join(d, f"{g}_stage2_actions.tsv")
        r["stage2_actions"] = max(0, sum(1 for _ in open(act)) - 1) if os.path.exists(act) else ""
        sug = os.path.join(d, f"{g}_rlgc_suggestions.tsv")
        r["rlgc_suggestions"] = max(0, sum(1 for _ in open(sug)) - 1) if os.path.exists(sug) else ""
        ops = os.path.join(d, f"{g}_operons.tsv")
        r["operon_chains"] = sum(1 for _ in open(ops)) if os.path.exists(ops) else ""
        gaps = os.path.join(d, f"{g}_gaps.jsonl")
        r["gap_records"] = sum(1 for _ in open(gaps)) if os.path.exists(gaps) else ""
        sc = os.path.join(d, f"{g}_score.tsv")
        if os.path.exists(sc):
            srows = list(csv.DictReader(open(sc), delimiter="\t"))
            if srows:
                for k, v in srows[0].items():
                    if k not in r:
                        r[f"score_{k}"] = v
                        if f"score_{k}" not in score_cols:
                            score_cols.append(f"score_{k}")
        for name, fn in (("gapseq", "timing_gapseq.tsv"), ("integrate", "timing.tsv")):
            for k, v in read_kv_times(os.path.join(d, fn)).items():
                r[f"wall_{k}"] = v
        rows.append(r)
    cols = []
    for r in rows:
        for k in r:
            if k not in cols:
                cols.append(k)
    with open(os.path.join(a.out_dir, "eq_mags_per_bin.tsv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    def hms(s):
        p = s.split(":")
        try:
            p = [float(x) for x in p]
        except ValueError:
            return None
        return p[0] * 3600 + p[1] * 60 + p[2] if len(p) == 3 else p[0] * 60 + p[1]

    summ = []
    datasets = sorted({r["dataset"] for r in rows if r["dataset"]})
    groups = [("class", c) for c in ("HQ", "MQ")] + [("dataset", d) for d in datasets] + [("all", "all")]
    for key, val in groups:
        sub = [r for r in rows if key == "all" or r[key] == val]
        done = [r for r in sub if r["integrated"]]
        med = lambda key: statistics.median([float(r[key]) for r in done if r.get(key) not in ("", None)]) if any(r.get(key) not in ("", None) for r in done) else ""
        medt = lambda key: statistics.median([hms(r[key]) for r in done if r.get(key) and hms(r[key]) is not None]) if any(r.get(key) for r in done) else ""
        summ.append({
            "group": val, "n_bins": len(sub), "proteins_done": sum(r["proteins_done"] for r in sub),
            "gapseq_done": sum(r["gapseq_done"] for r in sub), "integrated": len(done),
            "pyrodigal_bins": sum(1 for r in sub if r["protein_source"] == "pyrodigal"),
            "gapsmith_computed": sum(1 for r in sub if r["gapseq_source"] == "computed_gapsmith"),
            "median_proteins": med("n_proteins"), "median_claims": med("n_claims"),
            "median_stage0_proteins": med("stage0_proteins"), "median_stage12_proteins": med("stage12_proteins"),
            "median_stage12_rows": med("stage12_rows"), "median_stage2_actions": med("stage2_actions"),
            "median_rlgc_suggestions": med("rlgc_suggestions"),
            "median_wall_gapsmith_s": medt("wall_gapsmith"), "median_wall_stage12_s": medt("wall_stage12"),
        })
    with open(os.path.join(a.out_dir, "eq_mags_summary.tsv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(summ[0].keys()), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(summ)
    for s in summ:
        print(s)


if __name__ == "__main__":
    main()
