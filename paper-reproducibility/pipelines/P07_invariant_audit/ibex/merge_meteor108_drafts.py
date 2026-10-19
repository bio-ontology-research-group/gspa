#!/usr/bin/env python3
"""P07 merge: per-genome pickles -> meteor108_drafts_compact.pkl.gz + manifest_sha256.tsv."""
import os, sys, glob, gzip, pickle, hashlib, json
R = "/ibex/scratch/projects/c2014/kexin/funcarve/psb_revision"
os.environ.setdefault("METEOR_DATA", f"{R}/code_snapshot/data")
sys.path.insert(0, f"{R}/code_snapshot/src")
import numpy as np
from meteor_v8.utils import (load_universal, extract_fba_matrices, load_tight_bounds, apply_media,
                             data_path, find_excluded_reactions)
D = "/ibex/scratch/projects/c2014/hohndor/meteor108_repro_2026-10-19"


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 22), b""):
            h.update(c)
    return h.hexdigest()


panel = [l.rstrip("\n").split("\t")[:2] for l in open(f"{R}/code_snapshot/panel108_gram.tsv") if l.strip()]
universal, allrxns, _ = load_universal()
S0, lb0, ub0 = extract_fba_matrices(universal, allrxns, reversed_trans=True)
media = {}
for gram in ("negative", "positive"):
    bio = "biomass_GmPos" if gram == "positive" else "biomass_GmNeg"
    lt, ut = load_tight_bounds(data_path("tight_bounds_v6_%s.pkl" % gram[:3]))
    lb = np.maximum(lb0, lt); ub = np.minimum(ub0, ut)
    exc = find_excluded_reactions(S0, lb, ub, allrxns, bio)
    lb, ub, mm, mr = apply_media(["default"], allrxns, lb.copy(), ub.copy())
    media[gram] = {"lb": lb.astype(np.float64), "ub": ub.astype(np.float64), "biomass_id": bio,
                   "biomass_idx": allrxns.index(bio), "excluded_idx": np.array(sorted(exc), np.int32),
                   "media_exchange_idx": np.array(sorted(mr), np.int32)}
genomes, missing = {}, []
for gca, gram in panel:
    p = f"{D}/per_genome/{gca}.pkl"
    if not os.path.exists(p):
        missing.append(gca); continue
    genomes[gca] = pickle.load(open(p, "rb"))
out = {
    "schema_version": 1,
    "description": ("METEOR v8 108-genome invariant-audit inputs, DeepProZyme-vanilla (dpz_vanilla) baseline, tau=0.5. "
                    "All *_idx arrays index `allrxns` (universal.pickle reaction order == fba_matrices_v6RE.pkl columns). "
                    "media[gram].lb/ub are METEOR's own MILP bounds (tight FVA bounds + apply_media(['default'])); "
                    "with reversed_trans S, v>0 on EX_ = uptake, so blocked exchanges (lb=0, ub=1000) are OPEN UPTAKE / no secretion."),
    "allrxns": list(allrxns),
    "S_sha256": hashlib.sha256(S0.tocsr().data.tobytes() + S0.tocsr().indices.tobytes() + S0.tocsr().indptr.tobytes()).hexdigest(),
    "media": media, "genomes": genomes, "missing": missing,
    "panel": panel,
}
dst = f"{D}/meteor108_drafts_compact.pkl.gz"
with gzip.open(dst, "wb", compresslevel=9) as f:
    pickle.dump(out, f, protocol=4)
rows = []
for f in sorted(glob.glob(f"{D}/per_genome/*.sha256.tsv")):
    rows += open(f).read().splitlines()
for nm in ("universal.pickle", "fba_matrices_v6RE.pkl", "tight_bounds_v6_neg.pkl", "tight_bounds_v6_pos.pkl",
           "seedr2ec.pkl", "all_ancestors.txt", "medium.pkl", "enzymeobsolete.pkl"):
    p = data_path(nm)
    rows.append(f"*\tcode_snapshot_data\t{p}\t{os.path.getsize(p)}\t{sha256(p)}")
for nm in ("src/meteor_v8/repair.py", "src/meteor_v8/utils.py", "src/meteor_v8/milp_v8.py", "eval/baseline_io.py", "eval/emit_v8.py", "panel108_gram.tsv"):
    p = f"{R}/code_snapshot/{nm}"
    rows.append(f"*\tcode_snapshot\t{p}\t{os.path.getsize(p)}\t{sha256(p)}")
for nm in ("eval/thresh_gapfill_deadends.py", "results/thresh_gapfill_deadends_summary.json"):
    p = f"{R}/{nm}"
    rows.append(f"*\tpsb_revision\t{p}\t{os.path.getsize(p)}\t{sha256(p)}")
rows.append(f"*\toutput\t{dst}\t{os.path.getsize(dst)}\t{sha256(dst)}")
with open(f"{D}/manifest_sha256.tsv", "w") as f:
    f.write("gca\trole\tpath\tbytes\tsha256\n" + "\n".join(rows) + "\n")
nm = sum(1 for g in genomes.values() if g.get("tgd_match"))
print(f"merged {len(genomes)} genomes, missing={missing}, tgd_match={nm}/{len(genomes)}, size={os.path.getsize(dst)/1e6:.1f} MB")
