#!/usr/bin/env python3
"""Smoke test: GSPA's own comparator implementation against the archived P07 values.

For every genome of benchmark/manifests/reference_data/meteor108_drafts_compact.pkl.gz this
recomputes, with scripts/meteor_threshold_repair.py only,
  * the MILP medium bounds and the excluded reaction set (archive media[gram]),
  * the Grow-Support repair from the archived core (gs_added_idx, LP objective, growth),
  * the maximal biomass of core and repaired sets (gs_*_growth_milp_medium),
  * the structural profiles of the repaired, raw-threshold and METEOR networks
    (tgd["repaired" | "raw_threshold" | "meteor"]: n_rxn, deadends, mass_imbal, fba_growth),
and writes one row per genome plus a summary.  The threshold draft itself needs the per-protein
prediction matrices, which live on ibex; ibex/extract_meteor108_drafts.py compares that part.

Usage: python3 smoke_test_reimpl.py [--genomes GCA,GCA,...] [--out results/smoke_test_reimpl.tsv]
"""
from __future__ import annotations

import argparse
import gzip
import json
import pickle
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "scripts"))
from meteor_threshold_repair import UniversalData, NetworkProfiler, grow_support, maxbio_active  # noqa: E402

ARCHIVE = ROOT / "benchmark" / "manifests" / "reference_data" / "meteor108_drafts_compact.pkl.gz"
DATA = ROOT / "benchmark" / "manifests" / "reference_data" / "psb_revision" / "code_snapshot" / "data"


def mask(idx, n):
    m = np.zeros(n, dtype=bool)
    m[np.asarray(idx, dtype=np.int64)] = True
    return m


def jaccard(a: set, b: set) -> float:
    return 1.0 if not a and not b else len(a & b) / len(a | b)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genomes", default="")
    ap.add_argument("--out", default=str(Path(__file__).resolve().parent / "results" / "smoke_test_reimpl.tsv"))
    ap.add_argument("--skip-profiles", action="store_true")
    a = ap.parse_args()

    arc = pickle.load(gzip.open(ARCHIVE))
    U = UniversalData(DATA)
    if list(arc["allrxns"]) != U.allrxns:
        sys.exit("reaction order of the archive differs from universal_metadata.json")
    n = len(U.allrxns)
    prof = None if a.skip_profiles else NetworkProfiler(U)

    media_rows = []
    bounds = {}
    for gram in ("negative", "positive"):
        lb, ub, excl, open_ex = U.milp_bounds(gram)
        am = arc["media"][gram]
        bounds[gram] = (lb, ub, excl)
        media_rows.append({
            "gram": gram,
            "lb_equal": bool(np.allclose(lb, am["lb"])), "ub_equal": bool(np.allclose(ub, am["ub"])),
            "n_excluded": len(excl), "n_excluded_archive": len(am["excluded_idx"]),
            "excluded_equal": set(excl) == set(int(x) for x in am["excluded_idx"]),
            "open_exchanges_equal": open_ex == set(int(x) for x in am["media_exchange_idx"]),
            "biomass_idx_equal": U.rxn_index[U.biomass_id(gram)] == int(am["biomass_idx"]),
        })
    for r in media_rows:
        print("[media]", json.dumps(r), flush=True)

    wanted = [g for g in a.genomes.split(",") if g] or list(arc["genomes"])
    rows = []
    for k, gca in enumerate(wanted, 1):
        g = arc["genomes"][gca]
        gram = g["gram"]
        lb, ub, excl = bounds[gram]
        oi = U.rxn_index[U.biomass_id(gram)]
        avail = np.ones(n, dtype=bool)
        avail[list(excl)] = False
        core = mask(g["core_idx"], n)
        draft = mask(g["draft_idx"], n)
        t0 = time.time()
        core_growth = maxbio_active(U.S, lb, ub, oi, core)
        sup, v, obj = grow_support(U.S, lb, ub, oi, avail, core=core)
        t_lp = time.time() - t0
        repaired = core | (sup if sup is not None else np.zeros(n, bool))
        added = set(np.where(repaired & ~core)[0].tolist())
        rep_growth = maxbio_active(U.S, lb, ub, oi, repaired)
        arc_added = set(int(x) for x in g["gs_added_idx"])
        arc_sup = np.asarray(g["gs_support_idx"], dtype=np.int64)
        arc_obj = float(np.abs(np.asarray(g["gs_v_val"]))[~core[arc_sup]].sum()) if len(arc_sup) else 0.0
        row = {
            "gca": gca, "gram": gram,
            "n_core": int(core.sum()), "repair_found": sup is not None, "repair_found_archive": bool(g["gs_repair_found"]),
            "n_added": len(added), "n_added_archive": len(arc_added), "added_jaccard": round(jaccard(added, arc_added), 4),
            "lp_objective": round(obj, 6) if obj is not None else "", "lp_objective_archive": round(arc_obj, 6),
            "core_growth": round(core_growth, 4), "core_growth_archive": round(float(g["gs_draft_growth_milp_medium"]), 4),
            "repaired_growth": round(rep_growth, 4), "repaired_growth_archive": round(float(g["gs_repaired_growth_milp_medium"]), 4),
            "lp_sec": round(t_lp, 1),
        }
        if prof is not None:
            t0 = time.time()
            # The archived MEMOTE dead-end counts were computed on a cobra sub-model built from the
            # universal model's own reaction and metabolite objects, so each metabolite's reaction list
            # still spanned the whole universal model.  ``deadends_shared_objects`` replicates that
            # (every universal reaction of a network metabolite is consulted) to document the
            # discrepancy; GSPA reports the dead ends of the network itself.
            for label, sel, key in (("repaired", repaired, "repaired"), ("archive_repaired", mask(g["gs_repaired_idx"], n), "repaired"),
                                    ("draft", draft, "raw_threshold"), ("meteor", mask(g["meteor_active_idx"], n), "meteor")):
                p = prof.profile(sel, gram, lb, ub)
                t = g["tgd"].get(key, {})
                row[f"{label}_n_rxn"] = p["n_rxn"]
                row[f"{label}_deadends"] = p["deadends"]
                row[f"{label}_mass_imbal"] = p["mass_imbal"]
                row[f"{label}_fba_growth"] = p["fba_growth"]
                row[f"{label}_deadends_shared_objects"] = prof.deadends_shared_objects(sel, gram, lb, ub)
                if label != "repaired":
                    row[f"{label}_n_rxn_archive"] = t.get("n_rxn", "")
                    row[f"{label}_deadends_archive"] = t.get("deadends", "")
                    row[f"{label}_mass_imbal_archive"] = t.get("mass_imbal", "")
                    row[f"{label}_fba_growth_archive"] = t.get("fba_growth", "")
            row["profile_sec"] = round(time.time() - t0, 1)
        rows.append(row)
        print(f"[{k}/{len(wanted)}] {json.dumps(row)}", flush=True)

    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    cols = sorted({c for r in rows for c in r}, key=lambda c: (c not in ("gca", "gram"), c))
    with open(out, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    with open(out.with_suffix(".media.json"), "w") as f:
        json.dump(media_rows, f, indent=1)

    summ = {
        "n_genomes": len(rows),
        "media_identical": all(r["lb_equal"] and r["ub_equal"] and r["excluded_equal"] and r["open_exchanges_equal"] for r in media_rows),
        "repair_found_agree": rate_key(rows, lambda r: r["repair_found"] == r["repair_found_archive"]),
        "n_added_identical": rate_key(rows, lambda r: r["n_added"] == r["n_added_archive"]),
        "added_set_identical": rate_key(rows, lambda r: r["added_jaccard"] == 1.0),
        "added_jaccard_mean": round(float(np.mean([r["added_jaccard"] for r in rows])), 4),
        "lp_objective_within_1e-4": rate_key(rows, lambda r: r["lp_objective"] != "" and abs(r["lp_objective"] - r["lp_objective_archive"]) <= 1e-4 * max(1.0, abs(r["lp_objective_archive"]))),
        "core_growth_within_1e-3": rate_key(rows, lambda r: abs(r["core_growth"] - r["core_growth_archive"]) <= 1e-3),
        "repaired_growth_within_1e-3": rate_key(rows, lambda r: abs(r["repaired_growth"] - r["repaired_growth_archive"]) <= 1e-3),
    }
    if prof is not None:
        for label in ("archive_repaired", "draft", "meteor"):
            for q in ("n_rxn", "deadends", "mass_imbal"):
                summ[f"{label}_{q}_identical"] = rate_key(rows, lambda r, l=label, q=q: r[f"{l}_{q}"] == r[f"{l}_{q}_archive"])
            summ[f"{label}_deadends_shared_objects_within_1"] = rate_key(rows, lambda r, l=label: r[f"{l}_deadends_archive"] != "" and abs(r[f"{l}_deadends_shared_objects"] - int(r[f"{l}_deadends_archive"])) <= 1)
            summ[f"{label}_fba_growth_within_1e-2"] = rate_key(rows, lambda r, l=label: r[f"{l}_fba_growth_archive"] != "" and abs(r[f"{l}_fba_growth"] - float(r[f"{l}_fba_growth_archive"])) <= 1e-2)
    with open(out.with_suffix(".summary.json"), "w") as f:
        json.dump(summ, f, indent=1)
    print("[summary]", json.dumps(summ, indent=1))


def rate_key(rows, pred):
    return round(sum(1 for r in rows if pred(r)) / max(1, len(rows)), 4)


if __name__ == "__main__":
    main()
