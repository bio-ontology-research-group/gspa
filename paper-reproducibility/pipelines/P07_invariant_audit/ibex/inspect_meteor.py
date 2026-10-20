"""Inspect one METEOR dpz_vanilla output set + repair_rerun + tgd/table1 JSON (read-only)."""
import os, sys, json, pickle, glob
import numpy as np, pandas as pd
F = "/ibex/scratch/projects/c2014/kexin/funcarve"
RUNS = f"{F}/meteor_v8_evw_p2mu3_run/meteor_out"
G = "GCA_003598085.1"
sol = pickle.load(open(f"{RUNS}/dpz_vanilla/meteor_sol_{G}.pkl", "rb"))
print("SOL keys", {k: (type(v).__name__, getattr(v, "shape", v if np.isscalar(v) or isinstance(v, str) else len(v))) for k, v in sol.items()})
y = np.asarray(sol["y_vals"]); v = np.asarray(sol["v_vals"])
print("n y>0.5", int((y > 0.5).sum()), "n |v|>1e-6", int((np.abs(v) > 1e-6).sum()), "n active&|v|>1e-6", int(((y > 0.5) & (np.abs(v) > 1e-6)).sum()),
      "n active&|v|>1e-9", int(((y > 0.5) & (np.abs(v) > 1e-9)).sum()), "frac y", int(((y > 1e-6) & (y < 1 - 1e-6)).sum()))
pr = pickle.load(open(f"{RUNS}/dpz_vanilla/meteor_preds_{G}.pkl", "rb"))
print("PREDS keys", {k: len(v) for k, v in pr.items()})
ra = f"{RUNS}/dpz_vanilla/recon_active_{G}.pkl"
try:
    r = pickle.load(open(ra, "rb")); print("RECON_ACTIVE", type(r).__name__, (list(r.keys())[:10] if isinstance(r, dict) else getattr(r, "shape", len(r))))
except Exception as e:
    print("RECON_ACTIVE err", e)
df = pd.read_pickle(f"{RUNS}/dpz_vanilla/meteor_df_{G}.pkl")
print("DF", type(df).__name__, df.shape, "index[:3]", list(df.index[:3]), "cols[:5]", list(df.columns[:5]), df.dtypes.value_counts().to_dict())
for sub in sorted(glob.glob(f"{F}/psb_revision/meteor_out/repair_rerun/*")):
    fs = glob.glob(f"{sub}/meteor_sol_*.pkl")
    print("REPAIR_RERUN", os.path.basename(sub), len(fs))
s2 = pickle.load(open(f"{F}/psb_revision/meteor_out/repair_rerun/dpz_filt30/meteor_sol_{G}.pkl", "rb"))
print("RR SOL keys", list(s2.keys()), "n_repaired", s2.get("n_repaired"))
print("TGD", json.load(open(f"{F}/psb_revision/results/thresh_gapfill_deadends/tgd_{G}.json")))
t = json.load(open(f"{F}/psb_revision/results/table1_v2/table1_{G}.json"))
print("T1V2 keys", list(t.keys())); print({k: t[k] for k in t if k.startswith(("baseline_dpz", "meteor_dpz"))})
meta = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "universal_metadata.json")))
print("universal", len(meta["allrxns"]), len(meta["allmet"]), "y len", len(y))
