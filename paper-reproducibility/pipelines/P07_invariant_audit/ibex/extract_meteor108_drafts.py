#!/usr/bin/env python3
"""P07 per-genome extraction of the METEOR v8 108-genome invariant-audit inputs (DeepProZyme-vanilla).

Read-only inputs (colleague's tree, never written):
  meteor_v8_evw_p2mu3_run/meteor_out/dpz_vanilla/meteor_sol_<GCA>.pkl  (y_vals, v_vals of the deposited METEOR run)
  DeepProZyme-vanilla prediction pkl resolved by psb_revision/code_snapshot/eval/baseline_io.resolve_baseline_pkl
  psb_revision/results/thresh_gapfill_deadends/tgd_<GCA>.json         (Grow-Support counts, for verification)
  psb_revision/code_snapshot/{src,data}                                 (METEOR code + universal model data)

Per genome it reconstructs, with the archived code paths of psb_revision/eval/thresh_gapfill_deadends.py:
  draft  = reactions with any mapped EC scoring >= tau (0.5) in the DPZ-vanilla prediction matrix
  core   = draft minus structurally excluded reactions
  Grow-Support repair = core | grow_support(avail=~excluded, gamma_min=0.1, core=core)   [same LP, PuLP/CBC]
and records METEOR's active set (y > 0.5) and its MILP flux vector.
Writes <outdir>/<GCA>.pkl (small) and <outdir>/<GCA>.sha256.tsv.
"""
import sys, os, json, pickle, time, hashlib, argparse
R = "/ibex/scratch/projects/c2014/kexin/funcarve/psb_revision"
os.environ.setdefault("METEOR_DATA", f"{R}/code_snapshot/data")
for p in (f"{R}/code_snapshot/src", f"{R}/code_snapshot/eval"):
    if p not in sys.path:
        sys.path.insert(0, p)
import numpy as np
from scipy.sparse import csr_matrix as _csr
import pulp as _pulp

RUNS = "/ibex/scratch/projects/c2014/kexin/funcarve/meteor_v8_evw_p2mu3_run/meteor_out"
TGD = f"{R}/results/thresh_gapfill_deadends"
T1V2 = f"{R}/results/table1_v2"

ap = argparse.ArgumentParser()
ap.add_argument("--gca", required=True)
ap.add_argument("--gram", required=True, choices=["negative", "positive"])
ap.add_argument("--variant", default="vanilla")
ap.add_argument("--tau", type=float, default=0.5)
ap.add_argument("--outdir", required=True)
a = ap.parse_args()
os.makedirs(a.outdir, exist_ok=True)
OUT = os.path.join(a.outdir, f"{a.gca}.pkl")
if os.path.exists(OUT):
    print("done", flush=True); sys.exit(0)
t_start = time.time()


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


from meteor_v8.utils import (load_universal, extract_fba_matrices, load_tight_bounds, apply_media,
                             load_refmapping, load_ec, extract_pred, build_rxn_ec_mask, data_path,
                             data_dir, find_excluded_reactions)
from meteor_v8.repair import maxbio_active
from baseline_io import resolve_baseline_pkl, BASELINE_SUFFIX


def grow_support_with_flux(S, lb, ub, obj_idx, avail, gamma_min, core=None):
    """Verbatim meteor_v8.repair.grow_support (same LP, same CBC call) that additionally returns
    the LP flux vector, so the support can be audited. Support = |v| > 1e-9 as in the original."""
    Sc = _csr(S); Mn, Rn = Sc.shape
    mm = _pulp.LpProblem("grow", _pulp.LpMinimize)
    vpp = {j: _pulp.LpVariable("p%d" % j, 0, float(max(0.0, ub[j])) if avail[j] else 0.0) for j in range(Rn)}
    vnn = {j: _pulp.LpVariable("q%d" % j, 0, float(max(0.0, -lb[j])) if avail[j] else 0.0) for j in range(Rn)}
    for i in range(Mn):
        row = Sc.getrow(i)
        if len(row.indices) == 0:
            continue
        mm += _pulp.lpSum(float(row.data[k]) * (vpp[int(row.indices[k])] - vnn[int(row.indices[k])])
                          for k in range(len(row.indices))) == 0
    mm += (vpp[obj_idx] - vnn[obj_idx]) >= gamma_min
    if core is not None:
        mm += _pulp.lpSum((vpp[j] + vnn[j]) for j in range(Rn) if not core[j])
    else:
        mm += _pulp.lpSum(vpp[j] + vnn[j] for j in range(Rn))
    mm.solve(_pulp.PULP_CBC_CMD(msg=False))
    if _pulp.LpStatus[mm.status] != "Optimal":
        return None, None
    v = np.array([(vpp[j].value() or 0) - (vnn[j].value() or 0) for j in range(Rn)])
    return np.abs(v) > 1e-9, v


universal, allrxns, allmet = load_universal()
biomass_id = "biomass_GmPos" if a.gram == "positive" else "biomass_GmNeg"
S, lb, ub = extract_fba_matrices(universal, allrxns, reversed_trans=True)
lt, ut = load_tight_bounds(data_path("tight_bounds_v6_%s.pkl" % a.gram[:3]))
if lt is not None:
    lb = np.maximum(lb, lt); ub = np.minimum(ub, ut)
oi = allrxns.index(biomass_id)
exc = find_excluded_reactions(S, lb, ub, allrxns, biomass_id)
lb, ub, media_mask, media_rxns = apply_media(["default"], allrxns, lb, ub)
NR = len(allrxns)
seedr2ec, _ = load_refmapping(data_dir()); seedr2ec = {k: v for k, v in seedr2ec.items() if v}
anc = load_ec(data_path("all_ancestors.txt"))
mask = build_rxn_ec_mask(allrxns, seedr2ec, anc)

pred_path = resolve_baseline_pkl("dpz", a.variant, a.gca, BASELINE_SUFFIX["dpz"])
pred = extract_pred(pred_path, anc)
hit = (pred.values >= a.tau).any(axis=0)
draft = np.array([bool((mask[j] == 1).any() and hit[mask[j] == 1].any()) for j in range(NR)])
# per-reaction max DPZ score (noisy-OR-free max over mapped ECs), for audit only
colmax = pred.values.max(axis=0)
rxn_maxscore = np.array([float(colmax[mask[j] == 1].max()) if (mask[j] == 1).any() else 0.0 for j in range(NR)], dtype=np.float32)
avail = np.ones(NR, bool); avail[list(exc)] = False
core = draft & avail
draft_growth = maxbio_active(S, lb, ub, oi, core)
t0 = time.time()
sup, gs_v = grow_support_with_flux(S, lb, ub, oi, avail, 0.1, core=core)
t_repair = time.time() - t0
repaired = core | (sup if sup is not None else np.zeros(NR, bool))
added = repaired & ~core
rep_growth = maxbio_active(S, lb, ub, oi, repaired)

sol_path = f"{RUNS}/dpz_{a.variant}/meteor_sol_{a.gca}.pkl"
sol = pickle.load(open(sol_path, "rb"))
yv = np.asarray(sol["y_vals"], dtype=float); vv = np.asarray(sol["v_vals"], dtype=float)
met_active = yv > 0.5
met_vnz = np.where(np.abs(vv) > 1e-12)[0]

tgd_path = f"{TGD}/tgd_{a.gca}.json"
tgd = json.load(open(tgd_path)) if os.path.exists(tgd_path) else {}
t1_path = f"{T1V2}/table1_{a.gca}.json"

idx = lambda m: np.where(m)[0].astype(np.int32)
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
    "gs_draft_growth_milp_medium": float(draft_growth), "gs_repaired_growth_milp_medium": float(rep_growth),
    "gs_repair_sec": round(t_repair, 2),
    "tgd": tgd,
    "tgd_check": {k: tgd.get(k) for k in ("n_draft", "n_core", "n_added_by_repair", "repaired_growth_milp_medium", "draft_growth_milp_medium")},
    "runtime_sec": round(time.time() - t_start, 1),
}
rec["tgd_match"] = bool(tgd) and tgd.get("n_draft") == int(draft.sum()) and tgd.get("n_core") == int(core.sum()) \
    and tgd.get("n_added_by_repair") == int(added.sum())
pickle.dump(rec, open(OUT + ".tmp", "wb")); os.replace(OUT + ".tmp", OUT)
with open(os.path.join(a.outdir, f"{a.gca}.sha256.tsv"), "w") as f:
    for role, pth in (("meteor_sol", sol_path), ("dpz_pred", pred_path), ("tgd_json", tgd_path), ("table1_v2_json", t1_path)):
        if os.path.exists(pth):
            f.write(f"{a.gca}\t{role}\t{pth}\t{os.path.getsize(pth)}\t{sha256(pth)}\n")
print(f"[p07] {a.gca} draft={int(draft.sum())} core={int(core.sum())} met_active={int(met_active.sum())} "
      f"gs_added={int(added.sum())} tgd_added={tgd.get('n_added_by_repair')} match={rec['tgd_match']} "
      f"rep_growth={rep_growth:.4f} t={rec['runtime_sec']}s", flush=True)
