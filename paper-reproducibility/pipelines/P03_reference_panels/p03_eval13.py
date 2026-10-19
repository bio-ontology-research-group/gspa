#!/usr/bin/env python3
"""P03 step 4: evaluate all methods on the 13-genome panel with p03_eval.

Writes  {REPRO}/eval/panel13_per_genome.tsv   (method x genome x truth)
        {REPRO}/eval/panel13_summary.json     (mean + pooled per method x truth)
"""
import json
import os
import sys
from multiprocessing import Pool

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p03_eval as E
from p03_common import GO_OBO, REPRO, TAGS, METHODS, paths, pred_path

KEYS = ["n_truth_proteins", "n_truth_pairs", "fmax_micro", "fmax_micro_t", "fmax_cafa", "fmax_cafa_t",
        "coverage_05", "ic_recall_05", "smin", "smin_t"]


def run_method(method):
    go = E.GO(GO_OBO)
    rows, arrs = [], {"all": [], "exp": []}
    for tag in TAGS:
        P = paths(tag)
        ta = E.load_truth(P["truth_all"], go)
        te = E.load_truth(P["truth_exp"], go)
        ic = E.compute_ic(ta)
        pp = pred_path(tag, method)
        if not os.path.exists(pp):
            rows.append({"method": method, "tag": tag, "truth": "all", "status": "missing"})
            continue
        preds = E.load_preds(pp, go, set(ta) | set(te))
        npairs = sum(len(v) for v in preds.values())
        for name, tr in (("all", ta), ("exp", te)):
            if not tr:
                rows.append({"method": method, "tag": tag, "truth": name, "status": "empty_truth"})
                continue
            a = E.arrays(preds, tr, ic)
            arrs[name].append(a)
            m = E.metrics(a)
            m.update({"method": method, "tag": tag, "truth": name, "status": "ok",
                      "n_pred_pairs_prop": npairs})
            rows.append(m)
    summ = {"method": method}
    for name in ("all", "exp"):
        ok = [r for r in rows if r.get("truth") == name and r.get("status") == "ok"]
        summ[f"{name}_n_genomes"] = len(ok)
        for k in ("fmax_micro", "fmax_cafa", "coverage_05", "ic_recall_05", "smin"):
            vals = [r[k] for r in ok if not (isinstance(r[k], float) and np.isnan(r[k]))]
            summ[f"{name}_mean_{k}"] = float(np.mean(vals)) if vals else None
        if arrs[name]:
            pm = E.metrics(E.concat(arrs[name]))
            for k, v in pm.items():
                summ[f"{name}_pooled_{k}"] = v
    return rows, summ


def main():
    os.makedirs(f"{REPRO}/eval", exist_ok=True)
    with Pool(len(METHODS)) as pool:
        res = pool.map(run_method, METHODS)
    cols = ["method", "tag", "truth", "status", "n_pred_pairs_prop"] + KEYS
    with open(f"{REPRO}/eval/panel13_per_genome.tsv", "w") as f:
        f.write("\t".join(cols) + "\n")
        for rows, _ in res:
            for r in rows:
                f.write("\t".join("" if r.get(c) is None else (f"{r[c]:.6f}" if isinstance(r.get(c), float) else str(r.get(c, ""))) for c in cols) + "\n")
    with open(f"{REPRO}/eval/panel13_summary.json", "w") as f:
        json.dump([s for _, s in res], f, indent=1)
    for _, s in res:
        print(json.dumps({k: (round(v, 4) if isinstance(v, float) else v) for k, v in s.items()
                          if "mean_fmax" in k or "pooled_fmax_micro" == k[-17:] or k == "method"}), flush=True)


if __name__ == "__main__":
    main()
