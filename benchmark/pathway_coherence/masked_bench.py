#!/usr/bin/env python3
"""Pathway refinement on the reference genomes without close homologs.

The reference genomes' DIAMOND claims include hits to their own Swiss-Prot entries (identity near 100 %), which
leaves pathway refinement little to do. Here every DIAMOND claim with metadata.pident above --max-identity
(default 50) is removed; InterPro and Pfam claims stay. Nothing is attenuated: the remaining claims are the
evidence an organism without close characterized relatives would have. Each genome is integrated with refinement
off (`nocoh`), the v2.1.0 coherence prior (`legacy`) and module mode (`module`); F-max against the curated
annotation and the annotations each refinement adds relative to `nocoh` are reported (masked_summary.tsv).
"""
import argparse, csv, json, os, subprocess, time, sys
from concurrent.futures import ThreadPoolExecutor
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from run_bench9 import GENOMES, BENCH9, TRUTH, command, load_obo, load_truth, fmax, VARIANTS

ap = argparse.ArgumentParser()
ap.add_argument("--jar", required=True); ap.add_argument("--defs", required=True); ap.add_argument("--out-dir", required=True)
ap.add_argument("--max-identity", type=float, default=50.0); ap.add_argument("--variants", default="nocoh,legacy,module")
ap.add_argument("--claims-from", default="masked", help="masked | original"); ap.add_argument("--jobs", type=int, default=4)
a = ap.parse_args(); os.makedirs(a.out_dir, exist_ok=True); V = a.variants.split(",")
for g in GENOMES:
    out = os.path.join(a.out_dir, f"{g['tag']}_claims_masked.jsonl")
    if os.path.exists(out): continue
    kept = dropped = 0
    with open(f"{BENCH9}/with_iprscan/{g['tag']}_claims_merged.jsonl") as fi, open(out, "w") as fo:
        for line in fi:
            o = json.loads(line)
            pid = (o.get("metadata") or {}).get("pident")
            if o.get("source") == "diamond" and pid is not None and float(pid) > a.max_identity: dropped += 1; continue
            kept += 1; fo.write(line)
    print(f"{g['tag']}: kept {kept}, dropped {dropped} DIAMOND claims above {a.max_identity} % identity", flush=True)
thetas = {}
for v in V:
    if VARIANTS[v] is not None:
        tp = os.path.join(a.out_dir, f"theta_{v}.json"); json.dump(VARIANTS[v], open(tp, "w")); thetas[v] = tp
    else: thetas[v] = None
def run(job):
    g, v = job; out = os.path.join(a.out_dir, f"{g['tag']}_{v}_integrated.tsv")
    if os.path.exists(out) and os.path.getsize(out) > 1000: return g["tag"], v, "skipped"
    cmd = command(a.jar, g, v, thetas[v], a.defs, out); cmd[cmd.index("--claims") + 1] = os.path.join(a.out_dir, f"{g['tag']}_claims_masked.jsonl") if a.claims_from == "masked" else f"{BENCH9}/with_iprscan/{g['tag']}_claims_merged.jsonl"
    with open(out + ".log", "w") as log: rc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT).returncode
    return g["tag"], v, rc
with ThreadPoolExecutor(a.jobs) as ex:
    for r in ex.map(run, [(g, v) for g in GENOMES for v in V]): print("  ", r, flush=True)
anc = load_obo(f"/data/hohndor/gspa/reference/go.obo")
def above(path):
    s = set()
    for x in csv.DictReader(open(path), delimiter="\t"):
        if x["type"] == "GO" and float(x["posterior_prob"]) >= 0.5: s.add((x["protein_id"], x["function_id"]))
    return s
rows = []
for g in GENOMES:
    tag = g["tag"]; truth = load_truth(f"{TRUTH}/{tag}_truth_all.tsv"); T = {p: set().union(*(anc(t) for t in ts)) for p, ts in truth.items()}
    base = above(os.path.join(a.out_dir, f"{tag}_nocoh_integrated.tsv"))
    for v in V:
        f = os.path.join(a.out_dir, f"{tag}_{v}_integrated.tsv"); r = dict(genome=tag, variant=v); r.update(fmax(tag, f, 20))
        add = above(f) - base
        r["added_correct"] = sum(1 for p, t in add if p in T and t in T[p]); r["added_wrong"] = sum(1 for p, t in add if p in T and t not in T[p])
        rows.append(r); print(r, flush=True)
keys = ["fmax_micro_all", "fmax_cafa_all", "fmax_micro_exp", "fmax_cafa_exp"]
with open(os.path.join(a.out_dir, "masked_runs.tsv"), "w") as fh:
    cols = ["genome", "variant"] + keys + ["added_correct", "added_wrong"]; fh.write("\t".join(cols) + "\n")
    for r in rows: fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
with open(os.path.join(a.out_dir, "masked_summary.tsv"), "w") as fh:
    fh.write("variant\t" + "\t".join(keys) + "\tadded_correct\tadded_wrong\tgenomes_micro_all_above_nocoh\n")
    for v in V:
        sel = [r for r in rows if r["variant"] == v]
        better = sum(1 for r in sel if r["fmax_micro_all"] > next(x for x in rows if x["genome"] == r["genome"] and x["variant"] == "nocoh")["fmax_micro_all"])
        fh.write(v + "\t" + "\t".join(f"{sum(r[k] for r in sel) / len(sel):.4f}" for k in keys) +
                 f"\t{sum(r['added_correct'] for r in sel)}\t{sum(r['added_wrong'] for r in sel)}\t{better}\n")
print(open(os.path.join(a.out_dir, "masked_summary.tsv")).read())
