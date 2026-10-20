#!/usr/bin/env python3
"""P07 merge: per-genome pickles -> meteor108_drafts_compact.pkl.gz + manifest_sha256.tsv.

Uses GSPA's own universal-model loader (scripts/meteor_threshold_repair.py); no third-party pipeline
code.  Arguments: --work-dir (per_genome/ inside), --data-dir, --panel (wave<TAB>... no: GCA<TAB>gram),
--meteor-root (only to checksum the archived inputs listed in the per-genome sha256 files).
"""
import argparse
import glob
import gzip
import hashlib
import os
import pickle
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
for cand in (os.path.join(HERE, "..", "..", "..", "..", "scripts"), HERE):
    if os.path.exists(os.path.join(cand, "meteor_threshold_repair.py")):
        sys.path.insert(0, cand)
        break
from meteor_threshold_repair import UniversalData  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--work-dir", required=True)
ap.add_argument("--data-dir", required=True)
ap.add_argument("--panel", required=True, help="TSV: GCA<TAB>gram (108 rows)")
ap.add_argument("--out", default=None)
a = ap.parse_args()
D = a.work_dir
dst = a.out or f"{D}/meteor108_drafts_compact.pkl.gz"


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 22), b""):
            h.update(c)
    return h.hexdigest()


panel = [l.rstrip("\n").split("\t")[:2] for l in open(a.panel) if l.strip() and not l.startswith("#")]
U = UniversalData(a.data_dir)
media = {}
for gram in ("negative", "positive"):
    lb, ub, excl, open_ex = U.milp_bounds(gram)
    bio = U.biomass_id(gram)
    media[gram] = {"lb": lb.astype(np.float64), "ub": ub.astype(np.float64), "biomass_id": bio,
                   "biomass_idx": U.rxn_index[bio], "excluded_idx": np.array(sorted(excl), np.int32),
                   "media_exchange_idx": np.array(sorted(open_ex), np.int32)}
genomes, missing = {}, []
for gca, gram in panel:
    p = f"{D}/per_genome/{gca}.pkl"
    if not os.path.exists(p):
        missing.append(gca)
        continue
    genomes[gca] = pickle.load(open(p, "rb"))
Sc = U.S.tocsr()
out = {
    "schema_version": 2,
    "description": ("108-genome invariant-audit inputs, DeepProZyme-vanilla (dpz_vanilla) baseline, tau=0.5, computed with "
                    "GSPA's scripts/meteor_threshold_repair.py (threshold draft, excluded reactions, default medium, "
                    "Grow-Support LP with HiGHS, network profiles). All *_idx arrays index `allrxns` "
                    "(universal_metadata.json order == fba_matrices_v6RE.pkl columns). media[gram].lb/ub are METEOR's "
                    "MILP bounds (tight FVA bounds + default medium); with reversed-sign exchange columns, v>0 on EX_ is "
                    "uptake, so blocked exchanges (lb=0, ub=1000) are open for secretion only. Archived third-party "
                    "values (METEOR solutions, tgd_<GCA>.json counts) are stored next to the recomputed ones."),
    "allrxns": list(U.allrxns),
    "S_sha256": hashlib.sha256(Sc.data.tobytes() + Sc.indices.tobytes() + Sc.indptr.tobytes()).hexdigest(),
    "media": media, "genomes": genomes, "missing": missing, "panel": panel,
}
with gzip.open(dst, "wb", compresslevel=9) as f:
    pickle.dump(out, f, protocol=4)
rows = []
for f in sorted(glob.glob(f"{D}/per_genome/*.sha256.tsv")):
    rows += open(f).read().splitlines()
for nm in ("fba_matrices_v6RE.pkl", "tight_bounds_v6_neg.pkl", "tight_bounds_v6_pos.pkl", "seedr2ec.pkl",
           "all_ancestors.txt", "medium.pkl", "enzymeobsolete.pkl", "universal_metadata.json"):
    p = os.path.join(a.data_dir, nm)
    if not os.path.exists(p) and nm == "universal_metadata.json":
        p = str(U.data_dir / nm)
    if os.path.exists(p):
        rows.append(f"*\tdata\t{p}\t{os.path.getsize(p)}\t{sha256(p)}")
rows.append(f"*\tpanel\t{a.panel}\t{os.path.getsize(a.panel)}\t{sha256(a.panel)}")
rows.append(f"*\toutput\t{dst}\t{os.path.getsize(dst)}\t{sha256(dst)}")
with open(f"{D}/manifest_sha256.tsv", "w") as f:
    f.write("gca\trole\tpath\tbytes\tsha256\n" + "\n".join(rows) + "\n")
nm_dc = sum(1 for g in genomes.values() if g.get("tgd_match_draft_core"))
nm_all = sum(1 for g in genomes.values() if g.get("tgd_match"))
print(f"merged {len(genomes)} genomes, missing={missing}, draft/core match={nm_dc}/{len(genomes)}, "
      f"added-count match={nm_all}/{len(genomes)}, size={os.path.getsize(dst)/1e6:.1f} MB")
