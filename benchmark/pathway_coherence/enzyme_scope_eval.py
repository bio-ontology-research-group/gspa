#!/usr/bin/env python3
"""Evaluate refinement on the functions it addresses: the GO terms of KEGG module steps.

Scope = every GO term a KEGG module step maps to through ec2go, with its descendants. Predictions and the curated
truth (truth_all, propagated, then intersected with the scope) are restricted to that scope; per genome and variant
we report micro precision, recall and F1 at 0.50 and micro F-max over thresholds 0.01..0.99, plus pooled totals.
"""
import argparse, csv, json, os, re, sys
from collections import defaultdict
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from run_bench9 import GENOMES, load_obo, load_truth, REF, TRUTH
from held_out_steps import load_desc

ap = argparse.ArgumentParser(); ap.add_argument("--run-dir", required=True); ap.add_argument("--defs", required=True)
ap.add_argument("--variants", default="nocoh,legacy,module"); ap.add_argument("--label", default="")
a = ap.parse_args()
anc = load_obo(f"{REF}/go.obo"); desc = load_desc(f"{REF}/go.obo")
ec2go = {}
for line in open(f"{REF}/ec2go.txt"):
    m = re.match(r"^EC:(\S+)\s*>\s*GO:[^;]+;\s*(GO:\d{7})", line)
    if m and "-" not in m.group(1): ec2go[m.group(1)] = m.group(2)
scope = set()
for md in json.load(open(a.defs))["modules"].values():
    for o in md["orthology"].values():
        for e in o["ec"]:
            if e in ec2go: scope |= desc(ec2go[e])
print(f"scope: {len(scope)} GO terms", flush=True)
TH = [i / 100 for i in range(1, 100)]
rows = []; pooled = defaultdict(lambda: defaultdict(lambda: [0, 0, 0]))
for g in GENOMES:
    tag = g["tag"]; truth = load_truth(f"{TRUTH}/{tag}_truth_all.tsv")
    T = {p: (set().union(*(anc(t) for t in ts)) & scope) for p, ts in truth.items()}
    T = {p: s for p, s in T.items() if s}
    n_truth = sum(len(s) for s in T.values())
    for v in a.variants.split(","):
        f = os.path.join(a.run_dir, f"{tag}_{v}_integrated.tsv")
        if not os.path.exists(f): continue
        pred = {}
        for x in csv.DictReader(open(f), delimiter="\t"):
            if x["type"] == "GO" and x["function_id"] in scope and x["protein_id"] in T:
                k = (x["protein_id"], x["function_id"]); pred[k] = max(pred.get(k, 0), float(x["posterior_prob"]))
        best = (0, 0)
        for t in TH:
            tp = sum(1 for (p, go), s in pred.items() if s >= t and go in T[p]); n = sum(1 for s in pred.values() if s >= t)
            pooled[v][t][0] += tp; pooled[v][t][1] += n; pooled[v][t][2] += n_truth
            if n and tp:
                pr, rc = tp / n, tp / n_truth; f1 = 2 * pr * rc / (pr + rc)
                if f1 > best[0]: best = (f1, t)
        tp5 = sum(1 for (p, go), s in pred.items() if s >= 0.5 and go in T[p]); n5 = sum(1 for s in pred.values() if s >= 0.5)
        rows.append(dict(genome=tag, variant=v, truth_pairs=n_truth, pred_at_050=n5, tp_at_050=tp5, precision_050=tp5 / max(1, n5),
                         recall_050=tp5 / max(1, n_truth), fmax=best[0], fmax_threshold=best[1]))
out = os.path.join(a.run_dir, f"enzyme_scope{a.label}.tsv")
with open(out, "w") as fh:
    cols = list(rows[0].keys()); fh.write("\t".join(cols) + "\n")
    for r in rows: fh.write("\t".join(f"{r[c]:.4f}" if isinstance(r[c], float) else str(r[c]) for c in cols) + "\n")
    for v in a.variants.split(","):
        if v not in pooled: continue
        tp, n, nt = pooled[v][0.5]; best = max(((2 * (tp_ / n_) * (tp_ / nt_) / ((tp_ / n_) + (tp_ / nt_)), t) for t, (tp_, n_, nt_) in pooled[v].items() if tp_ and n_), default=(0, 0))
        better = sum(1 for r in rows if r["variant"] == v and r["recall_050"] > next(x for x in rows if x["genome"] == r["genome"] and x["variant"] == "nocoh")["recall_050"])
        fh.write(f"# pooled {v}: precision@0.5 {tp / max(1, n):.4f} recall@0.5 {tp / max(1, nt):.4f} TP {tp} predictions {n} truth {nt}; pooled Fmax {best[0]:.4f} at {best[1]}; genomes with higher recall@0.5 than nocoh: {better}\n")
print(open(out).read())
