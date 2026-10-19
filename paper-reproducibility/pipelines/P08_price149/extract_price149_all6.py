#!/usr/bin/env python3
"""Reproducible extraction of the six enzyme predictors' raw Price-149 outputs (P08).

Runs on ibex inside an sbatch job (never on the login node). Reads, read-only, the
per-host-genome prediction pickles in Kexin Hu's funcarve tree and writes one
compact pickle that scripts/eval_price149_real_all6.py consumes.

Output (gzip pickle) layout::

    {
      "_meta": {extraction_version, created, ec_vocab (list, 7,433 METEOR ancestor ECs),
                predictors: {key: {label, schema_note, files: {asm: {path, sha256, size}},
                                   raw_stats: {...}}},
                manifest: {path, sha256}, target_score_rule},
      "<pred>": {                       # clean, dpz, enzbert, graphec, mapred, topec
        "<host_assembly>": {
          "n_proteins": int,
          "target_ids": [price_entry, ...],            # rows of target_full
          "target_full": np.ndarray float32 (n_targets x len(ec_vocab)),  # LOSSLESS full vectors
          "target_scores": {price_entry: {ec: score}}, # every EC with score >= 0.01, plus the
                                                       # 100 highest-scoring ECs (lossless for all
                                                       # Top-k / MRR(<=100) / F-max(tau>=0.01) metrics)
          "ec_top_scores": {ec: [(protein_id, score), ...top-8, score >= 0.01]},
          "neighbor_ecs": {price_entry: [(offset, protein_id, [(ec, score) top-5, >= 0.15])]},
          "meteor_active_ecs": [...], "meteor_muted_ecs": [...],
          "raw_target_stats": {price_entry: {n_nonzero_raw, n_nonzero_vocab, max}},
        }
      }
    }

Differences from the v1 extractor (scripts/extract_ibex_compact_matrices.py::extract_price149_all6):
  * per-predictor loaders with explicit format checks (see SCHEMA below) instead of a single
    pd.read_pickle + "first index has three dots" orientation heuristic;
  * target vectors are stored in full (target_full) and target_scores is no longer capped at
    100 ECs nor rounded to 5 decimals (v1 rounded CLEAN's small scores to 0.0);
  * sha256 of every input pickle and of the manifest is recorded in _meta.

Usage (inside sbatch):
    python extract_price149_all6.py --manifest MANIFEST --out OUT.pkl.gz \
        [--predictors clean,dpz,enzbert,graphec,mapred,topec] [--sha256-tsv FILE]
"""

import argparse
import csv
import datetime as _dt
import glob
import gzip
import hashlib
import os
import pickle
import re
import sys
import time
from collections import defaultdict

import numpy as np
import pandas as pd

F = "/ibex/scratch/projects/c2014/kexin/funcarve"
D = f"{F}/meteor_v8/data"
MOP = f"{F}/meteor_v8_evw_p2mu3_run/meteor_out_price"
PA = f"{F}/paperA_2026"
EXTRACTION_VERSION = "v2-2026-10-19"
TOPK_TARGET = 100
MIN_TARGET_SCORE = 0.01

EC_RE = re.compile(r"^(EC:)?\d+\.[\dn\-]+\.[\dn\-]+\.[\dn\-]+$")

LABELS = {
    "clean": "CLEAN (2023)",
    "dpz": "DeepProZyme (2023)",
    "enzbert": "EnzBERT (2023)",
    "graphec": "GraphEC (2024)",
    "mapred": "MAPred (2025)",
    "topec": "TopEC (2025)",
}

# Per-predictor raw-format notes (inspection job 53986007, 2026-10-19; GCA_000006965.1 and
# GCF_000146165.2). All six raw files are pandas DataFrames, rows = proteins (proteome
# FASTA order, NCBI protein accessions), columns = 4-level ECs.
SCHEMA = {
    "clean": "DataFrame protein x 5,242 ECs (bare '1.1.1.1' labels), float64, DENSE: every EC has a "
             "score for every protein (min ~1e-58, max 0.9989); median 0 ECs >= 0.01 per protein, max 8. "
             "CLEAN 'confidence' (GMM-calibrated distance) output. v1 rounded these to 5 decimals.",
    "dpz": "DataFrame protein x 2,802 ECs ('EC:' prefix), float32, DENSE (all 2,802 > 0; ~72 >= 0.01 per "
           "protein). DeepProZyme/DeepECv2 '_t5' file (a sibling '_prXec' file with slightly lower scores "
           "exists for only some genomes; v1 and v2 both use '_t5').",
    "enzbert": "DataFrame protein x ~4,740 ECs (bare labels, incl. preliminary 'n' ECs, all present in the "
               "vocabulary), float64, SPARSE top-100: exactly 100 non-zero ECs per protein "
               "(row sums 0.16-1, i.e. a truncated softmax).",
    "graphec": "DataFrame protein x 7,062 ECs ('EC:' prefix, shared label set with MAPred/TopEC), float32, "
               "DENSE multi-label sigmoid: median 5,071 non-zero and ~100-230 ECs >= 0.01 per protein, "
               "max 0.999; hub labels (4.2.99.18, 3.6.1.31) score highly on many proteins. "
               "graphec_price_new/GCA_000006965.1 is byte-identical (md5) to graphec_price/GCA_000006965.",
    "mapred": "DataFrame protein x 7,062 ECs ('EC:' prefix), float32, SPARSE (scores < 0.01 zeroed): "
              "median 11-13, max ~180 non-zero ECs per protein. mapred_price_new is the float32 copy of "
              "mapred_price (same values, float64 there). Same shape and therefore same pickle byte size "
              "as GraphEC, but different content (md5 and DataFrame.equals differ).",
    "topec": "DataFrame protein x 7,068 ECs ('EC:' prefix), float64, ONE-HOT-WITH-SCORE: exactly one non-zero "
             "EC per protein (the argmax class with its softmax probability, 0.06-1.0). Every copy on ibex "
             "(topec_price/, topec_price_new/, paperA_2026/baseline_preds/price22_TopEC_clean/) has this "
             "form; the full softmax was not kept, so TopEC is a single-call predictor here.",
}


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 24), b""):
            h.update(b)
    return h.hexdigest()


def price_pkl_path(b, gca):
    gns = gca.rsplit(".", 1)[0]
    num = gca.split("_")[1].split(".")[0]
    m = {
        "clean": f"{PA}/baseline_preds/price22_CLEAN_clean/{gca}/{gca}_CLEAN_confidence.pkl",
        "enzbert": f"{F}/tfpc/resultprice_newg/{gca}_enzbert_predictions.pkl",
        "graphec": f"{F}/graphec_price_new/{gca}_GraphEC.pkl",
        "mapred": f"{F}/mapred_price_new/{gca}_MAPred.pkl",
        "topec": f"{F}/topec_price_new/{gca}_TopEC.pkl",
    }
    if b == "dpz":
        e = f"{F}/dpec2_result/result_price/{gns}_DeepECv2_t5.pkl"
        if os.path.exists(e):
            return e
        h = sorted(glob.glob(f"{F}/dpec2_result/result_price/GC?_{num}_DeepECv2_t5.pkl"))
        return h[0] if h else ""
    p = m[b]
    for c in (p, p.replace(gca, gns)):
        if os.path.exists(c):
            return c
    return ""


class Vocab:
    def __init__(self):
        self.anc = open(f"{D}/all_ancestors.txt").read().splitlines()
        self.pos = {e: i for i, e in enumerate(self.anc)}
        self.obsolete = pickle.load(open(f"{D}/enzymeobsolete.pkl", "rb"))

    def resolve(self, ec, max_hops=5):
        cur, seen, hops = ec, set(), 0
        while cur in self.obsolete and cur not in seen and hops < max_hops:
            seen.add(cur)
            cur = self.obsolete[cur]
            hops += 1
        return cur


def _is_ec(x):
    return bool(EC_RE.match(str(x).strip()))


def to_protein_by_ec(obj, b):
    """Return (DataFrame protein x EC (str columns, float), orientation_note)."""
    if isinstance(obj, pd.DataFrame):
        df = obj
        idx_ec = np.mean([_is_ec(x) for x in df.index[:50]]) > 0.9
        col_ec = np.mean([_is_ec(x) for x in df.columns[:50]]) > 0.9
        if idx_ec and not col_ec:
            df, note = df.T, "EC x protein (transposed)"
        elif col_ec:
            note = "protein x EC"
        else:
            raise ValueError(f"{b}: cannot identify EC axis; index[:3]={list(df.index[:3])} cols[:3]={list(df.columns[:3])}")
        df = df.copy()
        df.columns = [str(c).strip().split("EC:")[-1] for c in df.columns]
        df.index = [str(x) for x in df.index]
        return df.apply(pd.to_numeric, errors="coerce").fillna(0.0), note
    if isinstance(obj, dict):
        first = next(iter(obj.values()))
        if isinstance(first, dict):
            df = pd.DataFrame.from_dict(obj, orient="index").fillna(0.0)
            return to_protein_by_ec(df, b)[0], "dict protein -> {ec: score}"
        if isinstance(first, (list, tuple)):
            rec = {p: {str(e).split("EC:")[-1]: float(s) for e, s in v} for p, v in obj.items()}
            df = pd.DataFrame.from_dict(rec, orient="index").fillna(0.0)
            return df, "dict protein -> [(ec, score)]"
    raise TypeError(f"{b}: unsupported pickle type {type(obj)}")


def load_matrix(b, path, vocab):
    raw = pd.read_pickle(path)
    df, orient = to_protein_by_ec(raw, b)
    del raw
    cols = list(df.columns)
    vals = df.values.astype(np.float32)
    prot = list(df.index)
    out = np.zeros((len(prot), len(vocab.anc)), dtype=np.float32)
    n_in_vocab = n_renamed = 0
    dropped = []
    for j, c in enumerate(cols):
        k = vocab.pos.get(c)
        if k is None and c in vocab.obsolete:
            cur = vocab.resolve(c)
            k = vocab.pos.get(cur)
            if k is not None:
                n_renamed += 1
        if k is None:
            dropped.append(j)
            continue
        n_in_vocab += 1
        np.maximum(out[:, k], vals[:, j], out=out[:, k])
    info = {
        "orientation": orient,
        "n_proteins": len(prot),
        "n_ec_columns_raw": len(cols),
        "n_ec_columns_in_vocab": n_in_vocab,
        "n_ec_columns_obsolete_renamed": n_renamed,
        "n_ec_columns_dropped": len(dropped),
        "dropped_examples": [cols[j] for j in dropped[:10]],
        "n_full_4level_raw": sum(1 for c in cols if c.count(".") == 3 and "-" not in c),
        "score_min": float(np.nanmin(vals)) if vals.size else 0.0,
        "score_max": float(np.nanmax(vals)) if vals.size else 0.0,
        "nnz_per_protein_raw_median": float(np.median((vals > 0).sum(axis=1))),
        "nnz_per_protein_raw_max": int((vals > 0).sum(axis=1).max()),
        "nnz_ge001_per_protein_raw_median": float(np.median((vals >= 0.01).sum(axis=1))),
        "row_sum_median": float(np.median(vals.sum(axis=1))),
    }
    dropped_vals = vals[:, dropped] if dropped else None
    return prot, out, info, dropped_vals


def find_protein_idx(idx_map, index_list, pid, nid=None):
    for q in [pid] + ([nid] if nid else []):
        if not q:
            continue
        if q in idx_map:
            return idx_map[q]
        base = str(q).rsplit(".", 1)[0]
        if base in idx_map:
            return idx_map[base]
        for i, x in enumerate(index_list):
            if x == base or x.startswith(base + "."):
                return i
    return -1


def target_dict(row, ec_vocab):
    nz = np.where(row > 0.0)[0]
    keep = set(nz[row[nz] >= MIN_TARGET_SCORE].tolist())
    if len(nz):
        keep |= set(nz[np.argsort(-row[nz], kind="stable")[:TOPK_TARGET]].tolist())
    return {ec_vocab[k]: float(row[k]) for k in sorted(keep)}


def extract(manifest_path, out_path, predictors, sha_tsv=None):
    vocab = Vocab()
    rows = list(csv.DictReader(open(manifest_path), delimiter="\t"))
    assemblies = sorted({r["host_assembly"] for r in rows})
    targets_by_asm = defaultdict(list)
    for r in rows:
        targets_by_asm[r["host_assembly"]].append((r["price_entry"], r["proteome_protein_id"]))

    meta = {
        "extraction_version": EXTRACTION_VERSION,
        "created": _dt.datetime.now().isoformat(timespec="seconds"),
        "host": os.uname().nodename,
        "slurm_job_id": os.environ.get("SLURM_JOB_ID", ""),
        "ec_vocab": vocab.anc,
        "manifest": {"path": manifest_path, "sha256": sha256_file(manifest_path), "n_rows": len(rows)},
        "target_score_rule": f"target_scores = ECs with score >= {MIN_TARGET_SCORE} plus top-{TOPK_TARGET}; "
                             "target_full = full float32 vector over ec_vocab",
        "predictors": {},
    }
    data = {"_meta": meta}
    sha_rows = []
    for b in predictors:
        t0 = time.time()
        pm = {"label": LABELS[b], "schema_note": SCHEMA[b], "files": {}, "raw_stats": {}}
        data[b] = {}
        for asm in assemblies:
            bp = price_pkl_path(b, asm)
            if not bp:
                print(f"  WARN missing {b} {asm}", flush=True)
                continue
            sha = sha256_file(bp)
            st = os.stat(bp)
            pm["files"][asm] = {"path": bp, "sha256": sha, "size": st.st_size,
                                "mtime": _dt.datetime.fromtimestamp(st.st_mtime).isoformat(timespec="seconds")}
            sha_rows.append((b, asm, bp, st.st_size, sha, pm["files"][asm]["mtime"]))
            prot, mat, info, dropped_vals = load_matrix(b, bp, vocab)
            pm["raw_stats"][asm] = info
            idx_map = {x: i for i, x in enumerate(prot)}
            ec_vocab = vocab.anc

            ec_top = {}
            colmax = mat.max(axis=0)
            for c in np.where(colmax >= 0.01)[0]:
                col = mat[:, c]
                k = min(8, len(col))
                top = np.argpartition(col, -k)[-k:]
                top = top[np.argsort(-col[top], kind="stable")]
                pairs = [(prot[i], round(float(col[i]), 5)) for i in top if col[i] >= 0.01]
                if pairs:
                    ec_top[ec_vocab[c]] = pairs

            t_ids, t_rows, t_scores, neigh, rstats = [], [], {}, {}, {}
            for p_ent, p_id in targets_by_asm[asm]:
                r = find_protein_idx(idx_map, prot, p_id, p_ent)
                if r < 0:
                    print(f"  WARN {b} {asm}: target {p_ent}/{p_id} not found", flush=True)
                    continue
                row = mat[r]
                t_ids.append(p_ent)
                t_rows.append(row.copy())
                t_scores[p_ent] = target_dict(row, ec_vocab)
                rstats[p_ent] = {
                    "matched_id": prot[r],
                    "n_nonzero_vocab": int((row > 0).sum()),
                    "n_nonzero_dropped_cols": int((dropped_vals[r] > 0).sum()) if dropped_vals is not None else 0,
                    "max": float(row.max()),
                }
                nl = []
                for n in range(max(0, r - 10), min(len(prot), r + 11)):
                    if n == r:
                        continue
                    nv = mat[n]
                    t5 = np.argpartition(nv, -5)[-5:]
                    t5 = t5[np.argsort(-nv[t5], kind="stable")]
                    n_ecs = [(ec_vocab[c], round(float(nv[c]), 4)) for c in t5 if nv[c] >= 0.15]
                    if n_ecs:
                        nl.append((n - r, prot[n], n_ecs))
                neigh[p_ent] = nl

            m_path = f"{MOP}/{b}/meteor_preds_{asm}.pkl"
            m_active, m_muted = [], []
            if os.path.exists(m_path):
                md = pickle.load(open(m_path, "rb"))
                m_active = sorted(md.get("active_ecs", []))  # sorted: source is a set
                m_muted = sorted(md.get("muted_ecs", []))
            else:
                print(f"  WARN {b} {asm}: no METEOR output {m_path}", flush=True)

            data[b][asm] = {
                "n_proteins": len(prot),
                "target_ids": t_ids,
                "target_full": np.vstack(t_rows).astype(np.float32) if t_rows else np.zeros((0, len(ec_vocab)), np.float32),
                "target_scores": t_scores,
                "ec_top_scores": ec_top,
                "neighbor_ecs": neigh,
                "meteor_active_ecs": m_active,
                "meteor_muted_ecs": m_muted,
                "meteor_source": m_path if os.path.exists(m_path) else "",
                "raw_target_stats": rstats,
            }
            del mat, dropped_vals
            nn = [v["n_nonzero_vocab"] for v in rstats.values()]
            print(f"  {b} {asm}: {info['orientation']} {info['n_proteins']}x{info['n_ec_columns_raw']} "
                  f"(vocab {info['n_ec_columns_in_vocab']}, dropped {info['n_ec_columns_dropped']}) "
                  f"range=[{info['score_min']:.3g},{info['score_max']:.3g}] targets={len(t_ids)} "
                  f"nnz/target median={np.median(nn) if nn else 'NA'}", flush=True)
        meta["predictors"][b] = pm
        print(f"Done {b} in {time.time() - t0:.0f}s ({len(data[b])} assemblies)", flush=True)

    tmp = out_path + ".tmp"
    with gzip.open(tmp, "wb") as f:
        pickle.dump(data, f, protocol=4)
    os.replace(tmp, out_path)
    print(f"Saved {out_path} sha256={sha256_file(out_path)}")
    if sha_tsv:
        with open(sha_tsv, "w") as f:
            f.write("predictor\thost_assembly\tpath\tsize_bytes\tsha256\tmtime\n")
            for r in sha_rows:
                f.write("\t".join(map(str, r)) + "\n")
            f.write(f"manifest\t-\t{manifest_path}\t{os.path.getsize(manifest_path)}\t{meta['manifest']['sha256']}\t-\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", default="/ibex/scratch/projects/c2014/hohndor/suite6_price149_manifest.tsv")
    ap.add_argument("--out", required=True)
    ap.add_argument("--predictors", default="clean,dpz,enzbert,graphec,mapred,topec")
    ap.add_argument("--sha256-tsv", default=None)
    ap.add_argument("--merge", nargs="+", default=None,
                    help="merge per-predictor partial outputs (written with --predictors X) into --out")
    a = ap.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        sys.exit("refusing to run outside a Slurm job (ibex login-node policy)")
    if a.merge:
        merge(a.merge, a.out)
    else:
        extract(a.manifest, a.out, a.predictors.split(","), a.sha256_tsv)


def merge(parts, out_path):
    data = None
    for p in parts:
        with gzip.open(p, "rb") as f:
            d = pickle.load(f)
        if data is None:
            data = d
            continue
        for k, v in d.items():
            if k == "_meta":
                data["_meta"]["predictors"].update(v["predictors"])
                data["_meta"].setdefault("part_jobs", []).append(v.get("slurm_job_id", ""))
            else:
                data[k] = v
    data["_meta"]["merged"] = _dt.datetime.now().isoformat(timespec="seconds")
    tmp = out_path + ".tmp"
    with gzip.open(tmp, "wb") as f:
        pickle.dump(data, f, protocol=4)
    os.replace(tmp, out_path)
    print(f"Merged {len(parts)} parts -> {out_path} sha256={sha256_file(out_path)} "
          f"predictors={[k for k in data if k != '_meta']}")


if __name__ == "__main__":
    main()
