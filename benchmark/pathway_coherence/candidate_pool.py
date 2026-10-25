#!/usr/bin/env python3
"""Are the module rescues better than the candidate pool they were chosen from?

For each genome's module-mode run (--run-dir), the pool is every (protein, GO term) with a Stage 0 likelihood in
[0.16, 0.50) whose term is (or descends from) a GO term of an open required step of an activated module (from the
run's modules.json, before-refinement step status). Each pool pair and each module rescue is judged against
truth_all (propagated). A one-sided binomial test asks whether rescue precision exceeds the pool's.
"""
import argparse, csv, json, math, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from run_bench9 import GENOMES, load_obo, load_truth, REF, TRUTH
from held_out_steps import load_desc

def binom_sf(k, n, p):
    return sum(math.comb(n, i) * p ** i * (1 - p) ** (n - i) for i in range(k, n + 1))

ap = argparse.ArgumentParser(); ap.add_argument("--run-dir", required=True); a = ap.parse_args()
anc = load_obo(f"{REF}/go.obo"); desc = load_desc(f"{REF}/go.obo")
pool_ok = pool_n = res_ok = res_n = 0; rows = []
for g in GENOMES:
    tag = g["tag"]; f = os.path.join(a.run_dir, f"{tag}_module_integrated.tsv")
    if not os.path.exists(f): continue
    truth = load_truth(f"{TRUTH}/{tag}_truth_all.tsv"); T = {p: set().union(*(anc(t) for t in ts)) for p, ts in truth.items()}
    mj = json.load(open(f + ".modules.json"))
    terms = set()
    for m in mj["modules"]:
        if not m.get("activated_before") or m["signature"]: continue
        for st in m["steps"]:
            if st.get("required", True) and st.get("status_before", st["status"]) in ("evidence", "missing", "partial"):
                for t in st["go_terms"]: terms |= desc(t if isinstance(t, str) else t["id"])
    pk = pn = 0
    for x in csv.DictReader(open(f), delimiter="\t"):
        if x["type"] != "GO" or x["function_id"] not in terms or x["protein_id"] not in T: continue
        p0 = 1 / (1 + math.exp(-float(x["likelihood_logodds"])))
        if 0.16 <= p0 < 0.5:
            pn += 1; pk += x["function_id"] in T[x["protein_id"]]
    rk = rn = 0
    if os.path.exists(f + ".rescues.tsv"):
        for r in csv.DictReader(open(f + ".rescues.tsv"), delimiter="\t"):
            if r["prior"] != "coherence" or r["protein_id"] not in T: continue
            rn += 1; rk += r["function_id"] in T[r["protein_id"]]
    rows.append((tag, pn, pk, rn, rk)); pool_n += pn; pool_ok += pk; res_n += rn; res_ok += rk
pp = pool_ok / pool_n; rp = res_ok / res_n
with open(os.path.join(a.run_dir, "candidate_pool.tsv"), "w") as fh:
    fh.write("genome\tpool\tpool_correct\trescues\trescues_correct\n")
    for r in rows: fh.write("\t".join(map(str, r)) + "\n")
    fh.write(f"all\t{pool_n}\t{pool_ok}\t{res_n}\t{res_ok}\n")
    fh.write(f"# pool precision {pp:.3f}; rescue precision {rp:.3f}; one-sided binomial P(rescues >= {res_ok} | n={res_n}, p=pool) = {binom_sf(res_ok, res_n, pp):.2e}\n")
print(open(os.path.join(a.run_dir, "candidate_pool.tsv")).read())
