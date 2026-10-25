#!/usr/bin/env python3
"""Held-out module steps: does pathway refinement recover correct functions, and at what cost?

For each reference genome, from the round-4 module-mode run, pick one module step per activated, expected,
assessable module with completeness >= 2/3 whose annotated protein carries the step's function in the curated
annotation (truth_all, propagated) at Stage 0 >= 0.50. In a copy of the genome's claims, every claim of that
protein on the step's GO terms (and their descendants, and EC numbers mapping to them) gets calibrated_prob
drawn from U[0.25, 0.45] (seeded): the correct function is pushed into the borderline zone. The perturbed
claims are integrated without the coherence prior (`nocoh`), with the v2.1.0 coherence prior (`legacy`) and
in module mode (`module`), with the same flags as run_bench9.py, and we count
  * recovered targets: the protein again carries a step term at posterior >= 0.50;
  * other annotations added relative to `nocoh` (non-target pairs at >= 0.50), judged against truth_all.
Writes held_out_targets.tsv, held_out_runs.tsv and held_out_summary.tsv into --out-dir.
"""
import argparse, csv, json, os, random, re, subprocess, time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from run_bench9 import GENOMES, BENCH9, REF, TRUTH, command, load_obo, load_truth, VARIANTS

def load_desc(path):
    children, alt, cur = defaultdict(set), {}, None
    for line in open(path):
        line = line.rstrip("\n")
        if line == "[Term]": cur = None
        elif line.startswith("id: GO:"): cur = line[4:]
        elif cur and line.startswith("is_a: "): children[line[6:].split(" ")[0]].add(cur)
        elif cur and line.startswith("relationship: part_of "): children[line.split(" ")[2]].add(cur)
    cache = {}
    def desc(t):
        if t in cache: return cache[t]
        out, st = {t}, [t]
        while st:
            u = st.pop()
            for c in children.get(u, ()):
                if c not in out: out.add(c); st.append(c)
        cache[t] = out; return out
    return desc

def ec2go_map():
    m = defaultdict(set)
    for line in open(f"{REF}/ec2go.txt"):
        r = re.match(r"^(EC:\S+)\s*>\s*GO:[^;]+;\s*(GO:\d{7})", line)
        if r: m[r.group(1)].add(r.group(2))
    return m

def above(path):
    d = {}
    with open(path) as fh:
        for x in csv.DictReader(fh, delimiter="\t"):
            if x["type"] == "GO": d[(x["protein_id"], x["function_id"])] = float(x["posterior_prob"])
    return d

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jar", required=True); ap.add_argument("--defs", required=True)
    ap.add_argument("--round4", required=True, help="bench9_v4 dir (module-mode runs with modules.json)")
    ap.add_argument("--out-dir", required=True); ap.add_argument("--jobs", type=int, default=4); ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    anc = load_obo(f"{REF}/go.obo"); desc = load_desc(f"{REF}/go.obo"); e2g = ec2go_map()
    targets = {}
    for g in GENOMES:
        tag = g["tag"]; rng = random.Random(f"{a.seed}-{tag}")
        truth = load_truth(f"{TRUTH}/{tag}_truth_all.tsv")
        base = above(os.path.join(a.round4, f"{tag}_nocoh_integrated.tsv"))
        mj = json.load(open(os.path.join(a.round4, f"{tag}_module_integrated.tsv.modules.json")))
        tl, used_prot = [], set()
        for m in mj["modules"]:
            if not m["activated"] or m["signature"] or not m.get("assessable", True) or m.get("not_expected_reason"): continue
            if m["completeness"] is None or m["completeness"] < 2 / 3: continue
            cands = []
            for st in m["steps"]:
                if st["status"] != "present" or not st.get("required", True): continue
                terms = {t["id"] if isinstance(t, dict) else t for t in st["go_terms"]}
                for pr in st.get("proteins", []):
                    p, t = pr["protein"], pr["term"]
                    tt = truth.get(p)
                    if not tt or p in used_prot or base.get((p, t), 0) < 0.5: continue
                    if t in set().union(*(anc(x) for x in tt)):
                        cands.append((st["index"], p, sorted(terms)))
            if cands:
                si, p, terms = rng.choice(sorted(cands))
                used_prot.add(p)
                tl.append(dict(genome=tag, module=m["id"], step=si, protein=p, terms=terms, u=round(rng.uniform(0.25, 0.45), 3)))
        targets[tag] = tl
        print(f"{tag}: {len(tl)} held-out steps", flush=True)
    with open(os.path.join(a.out_dir, "held_out_targets.tsv"), "w") as fh:
        fh.write("genome\tmodule\tstep\tprotein\tterms\tu\n")
        for tl in targets.values():
            for t in tl: fh.write(f"{t['genome']}\t{t['module']}\t{t['step']}\t{t['protein']}\t{','.join(t['terms'])}\t{t['u']}\n")
    # perturbed claims
    for g in GENOMES:
        tag = g["tag"]; tl = targets[tag]
        closure = {}
        for t in tl:
            s = set()
            for x in t["terms"]: s |= desc(x)
            closure[t["protein"]] = (s, t["u"])
        n = 0
        with open(f"{BENCH9}/with_iprscan/{tag}_claims_merged.jsonl") as fi, open(os.path.join(a.out_dir, f"{tag}_claims_heldout.jsonl"), "w") as fo:
            for line in fi:
                o = json.loads(line)
                c = closure.get(o["protein_id"])
                if c:
                    fid = o["function_id"]; hit = fid in c[0] if o["function_type"] == "GO" else bool(e2g.get(fid if fid.startswith("EC:") else "EC:" + fid, set()) & c[0])
                    if hit:
                        o["calibrated_prob"] = c[1]; o["raw_score"] = c[1]; n += 1
                fo.write(json.dumps(o) + "\n")
        print(f"{tag}: {n} claims attenuated", flush=True)
    # integrate
    thetas = {}
    for v in ("nocoh", "legacy", "module"):
        if VARIANTS[v] is not None:
            tp = os.path.join(a.out_dir, f"theta_{v}.json"); json.dump(VARIANTS[v], open(tp, "w")); thetas[v] = tp
        else: thetas[v] = None
    jobs = [(g, v) for g in GENOMES for v in ("nocoh", "legacy", "module")]
    def run(job):
        g, v = job; out = os.path.join(a.out_dir, f"{g['tag']}_{v}_integrated.tsv")
        cmd = command(a.jar, g, v, thetas[v], a.defs, out)
        cmd[cmd.index("--claims") + 1] = os.path.join(a.out_dir, f"{g['tag']}_claims_heldout.jsonl")
        t0 = time.time()
        with open(out + ".log", "w") as log: pr = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
        return g["tag"], v, pr.returncode, time.time() - t0
    with ThreadPoolExecutor(a.jobs) as ex:
        for tag, v, rc, dt in ex.map(run, jobs): print(f"  {tag} {v}: rc={rc} {dt:.0f}s", flush=True)
    # evaluate
    rows = []
    for g in GENOMES:
        tag = g["tag"]; truth = load_truth(f"{TRUTH}/{tag}_truth_all.tsv"); tl = targets[tag]
        tset = {}
        for t in tl:
            s = set()
            for x in t["terms"]: s |= desc(x)
            tset[t["protein"]] = s
        post = {v: above(os.path.join(a.out_dir, f"{tag}_{v}_integrated.tsv")) for v in ("nocoh", "legacy", "module")}
        for v in ("nocoh", "legacy", "module"):
            rec = sum(1 for p, s in tset.items() if any(post[v].get((p, x), 0) >= 0.5 for x in s))
            ok = bad = unk = 0
            if v != "nocoh":
                base = {k for k, x in post["nocoh"].items() if x >= 0.5}
                for (p, gterm) in {k for k, x in post[v].items() if x >= 0.5} - base:
                    if p in tset and gterm in tset[p]: continue
                    tt = truth.get(p)
                    if not tt: unk += 1; continue
                    if gterm in set().union(*(anc(x) for x in tt)): ok += 1
                    else: bad += 1
            rows.append(dict(genome=tag, variant=v, held_out=len(tset), recovered=rec, other_added_correct=ok, other_added_wrong=bad, other_added_unknown=unk))
    with open(os.path.join(a.out_dir, "held_out_runs.tsv"), "w") as fh:
        cols = list(rows[0].keys()); fh.write("\t".join(cols) + "\n")
        for r in rows: fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    with open(os.path.join(a.out_dir, "held_out_summary.tsv"), "w") as fh:
        fh.write("variant\theld_out\trecovered\trecovery\tother_added_correct\tother_added_wrong\n")
        for v in ("nocoh", "legacy", "module"):
            sel = [r for r in rows if r["variant"] == v]; h = sum(r["held_out"] for r in sel); rc = sum(r["recovered"] for r in sel)
            fh.write(f"{v}\t{h}\t{rc}\t{rc / h:.3f}\t{sum(r['other_added_correct'] for r in sel)}\t{sum(r['other_added_wrong'] for r in sel)}\n")
    print(open(os.path.join(a.out_dir, "held_out_summary.tsv")).read())
    for r in rows: print(r)

if __name__ == "__main__":
    main()
