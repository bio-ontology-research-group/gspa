#!/usr/bin/env python3
"""Split the annotations a coherence variant adds relative to the no-coherence run (posterior >= 0.50, GO) into
module rescues and other coherence effects (process has_part and complex boosts, indirect effects), each judged
against truth_all. Writes added_breakdown.tsv in --out-dir."""
import argparse, csv, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from run_bench9 import GENOMES, load_obo, load_truth, REF, TRUTH

def above(path):
    s = set()
    with open(path) as fh:
        for x in csv.DictReader(fh, delimiter="\t"):
            if x["type"] == "GO" and float(x["posterior_prob"]) >= 0.5: s.add((x["protein_id"], x["function_id"]))
    return s

ap = argparse.ArgumentParser(); ap.add_argument("--out-dir", required=True); a = ap.parse_args()
anc = load_obo(f"{REF}/go.obo"); rows = []
for v in ("legacy", "module"):
    tot = {"rescue": [0, 0, 0], "other": [0, 0, 0]}
    for g in GENOMES:
        base = os.path.join(a.out_dir, f"{g['tag']}_nocoh_integrated.tsv"); f = os.path.join(a.out_dir, f"{g['tag']}_{v}_integrated.tsv")
        if not (os.path.exists(base) and os.path.exists(f)): continue
        truth = load_truth(f"{TRUTH}/{g['tag']}_truth_all.tsv")
        resc = set()
        rp = f + ".rescues.tsv"
        if os.path.exists(rp):
            resc = {(r["protein_id"], r["function_id"]) for r in csv.DictReader(open(rp), delimiter="\t") if r["prior"] == "coherence"}
        for pair in above(f) - above(base):
            k = "rescue" if pair in resc else "other"
            tt = truth.get(pair[0])
            if not tt: tot[k][2] += 1; continue
            ok = pair[1] in set().union(*(anc(t) for t in tt))
            tot[k][0 if ok else 1] += 1
    for k, (c, w, u) in tot.items():
        rows.append((v, k, c, w, u, c / (c + w) if c + w else float("nan")))
with open(os.path.join(a.out_dir, "added_breakdown.tsv"), "w") as fh:
    fh.write("variant\tsource\tcorrect\twrong\tunknown_truth\tprecision\n")
    for r in rows: fh.write("\t".join(map(str, r)) + "\n")
for r in rows: print(r)
