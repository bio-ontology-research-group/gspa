#!/usr/bin/env python3
"""P07 per-genome extraction of the 108-genome invariant-audit inputs (DeepProZyme-vanilla), with
GSPA's own implementation of the comparator (scripts/meteor_threshold_repair.py).

Read-only inputs (archived third-party results, never written):
  <meteor-root>/meteor_v8_evw_p2mu3_run/meteor_out/dpz_vanilla/meteor_sol_<GCA>.pkl   METEOR y_vals / v_vals
  DeepProZyme-vanilla prediction pkl: <meteor-root>/meteor_v7_run/reinfer/dpz_vanilla/<GCA>_DPZ.pkl when it
      exists (re-inferred genomes), else <meteor-root>/paperA_2026/baseline_preds/dpz_vanilla_genome_collection/
  <meteor-root>/psb_revision/results/thresh_gapfill_deadends/tgd_<GCA>.json              archived counts (comparison only)
  <data-dir>: fba_matrices_v6RE.pkl, tight_bounds_v6_{neg,pos}.pkl, medium.pkl, seedr2ec.pkl, all_ancestors.txt,
      enzymeobsolete.pkl, universal_metadata.json (the universal model as plain data)

Per genome:
  draft    = reactions with any mapped EC scoring >= tau in the prediction matrix
  core     = draft minus structurally excluded reactions
  repaired = core | Grow-Support(avail = not excluded, gamma_min = 0.1, core = core)   (HiGHS LP)
  profiles (dead ends, mass-unbalanced reactions, FBA growth) of the repaired, draft and METEOR networks
and METEOR's active set (y > 0.5) with its MILP flux vector.  Writes <outdir>/<GCA>.pkl and
<outdir>/<GCA>.sha256.tsv.  Every archived count is stored next to the recomputed one.
"""
import argparse
import glob
import hashlib
import json
import os
import pickle
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
for cand in (os.path.join(HERE, "..", "..", "..", "..", "scripts"), HERE):
    if os.path.exists(os.path.join(cand, "meteor_threshold_repair.py")):
        sys.path.insert(0, cand)
        break
from meteor_threshold_repair import (UniversalData, NetworkProfiler, load_prediction_matrix,  # noqa: E402
                                     threshold_draft, maxbio_active, grow_support, GAMMA_MIN)

ap = argparse.ArgumentParser()
ap.add_argument("--gca", required=True)
ap.add_argument("--gram", required=True, choices=["negative", "positive"])
ap.add_argument("--variant", default="vanilla")
ap.add_argument("--tau", type=float, default=0.5)
ap.add_argument("--outdir", required=True)
ap.add_argument("--data-dir", required=True)
ap.add_argument("--meteor-root", required=True, help="root holding meteor_v8_evw_p2mu3_run/, meteor_v7_run/, paperA_2026/, psb_revision/")
ap.add_argument("--no-profiles", action="store_true")
a = ap.parse_args()
os.makedirs(a.outdir, exist_ok=True)
OUT = os.path.join(a.outdir, f"{a.gca}.pkl")
if os.path.exists(OUT):
    print("done", flush=True)
    sys.exit(0)
t_start = time.time()


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def resolve_prediction(root, variant, gca):
    """Re-inferred prediction first, then the original genome collection, then the BacDive-scale run."""
    for p in (f"{root}/meteor_v7_run/reinfer/dpz_{variant}/{gca}_DPZ.pkl",
              f"{root}/paperA_2026/baseline_preds/dpz_{variant}_genome_collection/{gca}_DPZ.pkl"):
        if os.path.exists(p):
            return p
    if variant == "vanilla":
        hits = glob.glob(f"{root}/dpec2_result/result_bacdive/{gca.split('.')[0]}_*_DeepECv2_t5.pkl")
        if hits:
            return hits[0]
    return ""


U = UniversalData(a.data_dir)
NR = len(U.allrxns)
lb, ub, excluded, open_ex = U.milp_bounds(a.gram)
oi = U.rxn_index[U.biomass_id(a.gram)]
rxn_ecs = U.rxn_ec_mask()

pred_path = resolve_prediction(a.meteor_root, a.variant, a.gca)
if not pred_path:
    sys.exit(f"no DeepProZyme prediction for {a.gca}")
pred = load_prediction_matrix(pred_path, U.ancestors_ec, U.ec_obsolete)
draft, rxn_maxscore = threshold_draft(pred.values, rxn_ecs, a.tau)
avail = np.ones(NR, bool)
avail[list(excluded)] = False
core = draft & avail
draft_growth = maxbio_active(U.S, lb, ub, oi, core)
t0 = time.time()
sup, gs_v, gs_obj = grow_support(U.S, lb, ub, oi, avail, GAMMA_MIN, core=core)
t_repair = time.time() - t0
repaired = core | (sup if sup is not None else np.zeros(NR, bool))
added = repaired & ~core
rep_growth = maxbio_active(U.S, lb, ub, oi, repaired)

sol_path = f"{a.meteor_root}/meteor_v8_evw_p2mu3_run/meteor_out/dpz_{a.variant}/meteor_sol_{a.gca}.pkl"
sol = pickle.load(open(sol_path, "rb"))
yv = np.asarray(sol["y_vals"], dtype=float)
vv = np.asarray(sol["v_vals"], dtype=float)
met_active = yv > 0.5
met_vnz = np.where(np.abs(vv) > 1e-12)[0]

tgd_path = f"{a.meteor_root}/psb_revision/results/thresh_gapfill_deadends/tgd_{a.gca}.json"
tgd = json.load(open(tgd_path)) if os.path.exists(tgd_path) else {}
t1_path = f"{a.meteor_root}/psb_revision/results/table1_v2/table1_{a.gca}.json"

profiles = {}
if not a.no_profiles:
    prof = NetworkProfiler(U)
    for label, sel in (("repaired", repaired), ("raw_threshold", draft), ("meteor", met_active)):
        p = prof.profile(sel, a.gram, lb, ub)
        p["deadends_shared_objects"] = prof.deadends_shared_objects(sel, a.gram, lb, ub)
        profiles[label] = p

idx = lambda m: np.where(m)[0].astype(np.int32)  # noqa: E731
rec = {
    "gca": a.gca, "gram": a.gram, "variant": f"dpz_{a.variant}", "tau": a.tau,
    "pred_path": pred_path, "n_proteins": int(pred.shape[0]),
    "draft_idx": idx(draft), "core_idx": idx(core),
    "draft_maxscore": rxn_maxscore[draft],
    "meteor_active_idx": idx(met_active),
    "meteor_v_idx": met_vnz.astype(np.int32), "meteor_v_val": vv[met_vnz].astype(np.float64),
    "meteor_biomass_flux": float(sol.get("biomass_flux", vv[oi])), "meteor_status": sol.get("status"),
    "meteor_n_repaired": int(sol.get("n_repaired", -1)), "meteor_maxbio_core": float(sol.get("maxbio_core", float("nan"))),
    "meteor_n_active": int(sol.get("n_active", met_active.sum())),
    "gs_repair_found": sup is not None,
    "gs_repaired_idx": idx(repaired), "gs_added_idx": idx(added),
    "gs_support_idx": idx(sup) if sup is not None else np.zeros(0, np.int32),
    "gs_v_val": gs_v[sup].astype(np.float64) if sup is not None else np.zeros(0),
    "gs_lp_objective": float(gs_obj) if gs_obj is not None else float("nan"),
    "gs_draft_growth_milp_medium": float(draft_growth), "gs_repaired_growth_milp_medium": float(rep_growth),
    "gs_repair_sec": round(t_repair, 2),
    "profiles": profiles,
    "tgd": tgd,
    "tgd_check": {k: tgd.get(k) for k in ("n_draft", "n_core", "n_added_by_repair", "repaired_growth_milp_medium", "draft_growth_milp_medium")},
    "implementation": "gspa scripts/meteor_threshold_repair.py (scipy HiGHS)",
    "runtime_sec": round(time.time() - t_start, 1),
}
rec["tgd_match"] = bool(tgd) and tgd.get("n_draft") == int(draft.sum()) and tgd.get("n_core") == int(core.sum()) \
    and tgd.get("n_added_by_repair") == int(added.sum())
rec["tgd_match_draft_core"] = bool(tgd) and tgd.get("n_draft") == int(draft.sum()) and tgd.get("n_core") == int(core.sum())
pickle.dump(rec, open(OUT + ".tmp", "wb"))
os.replace(OUT + ".tmp", OUT)
with open(os.path.join(a.outdir, f"{a.gca}.sha256.tsv"), "w") as f:
    for role, pth in (("meteor_sol", sol_path), ("dpz_pred", pred_path), ("tgd_json", tgd_path), ("table1_v2_json", t1_path)):
        if os.path.exists(pth):
            f.write(f"{a.gca}\t{role}\t{pth}\t{os.path.getsize(pth)}\t{sha256(pth)}\n")
print(f"[p07] {a.gca} draft={int(draft.sum())} core={int(core.sum())} met_active={int(met_active.sum())} "
      f"gs_added={int(added.sum())} tgd_added={tgd.get('n_added_by_repair')} draft_core_match={rec['tgd_match_draft_core']} "
      f"added_match={rec['tgd_match']} rep_growth={rep_growth:.4f} tgd_rep_growth={tgd.get('repaired_growth_milp_medium')} "
      f"t={rec['runtime_sec']}s", flush=True)
