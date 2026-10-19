#!/usr/bin/env python3
"""P03 single evaluator (used for every method, both panels).

Metric definitions are those of benchmark/benchmark_pgap_v2.py (reused by
benchmark/unified/pillar1_2/02_run_cafa_eval.py):
  micro F-max   TP/FP/FN summed over all (protein, term) pairs of truth proteins
  CAFA F-max    precision averaged over truth proteins with >=1 prediction at t,
                recall averaged over all truth proteins
  coverage      fraction of truth proteins with >=1 prediction >= 0.5
  IC recall     sum IC of recovered truth terms / sum IC of truth terms at t=0.5
  S-min         min_t sqrt(ru^2 + mi^2), ru/mi = mean IC of missed / extra terms
  IC(g)         -log2(n_proteins annotated with g / n_truth_proteins), from the
                genome's propagated ALL-evidence truth (used for exp truth too,
                because exp truth sets hold 1-277 rows per genome)
with ONE change required by the paper protocol: ancestor propagation.
Truth and predictions are closed under is_a + part_of over GO
(go.obo releases/2026-03-25; alt_ids mapped to primary ids, obsolete terms
dropped); a propagated prediction carries the max score of its descendants;
the three root terms are removed from truth and predictions. Aspect is taken
from the GO namespace. Thresholds 0.01..0.99 in steps of 0.01 (all metrics,
including S-min); prediction rows with score < 0.01 are discarded.

Panel aggregation: mean over genomes of the per-genome metric ("mean"), and a
pooled evaluation in which all genomes' (protein, term) pairs are evaluated as
one set ("pooled"). For experimental-only truth the pooled value is the
meaningful one (per-genome exp sets are tiny).
"""
import argparse
import json
import math
import os
import sys
from collections import defaultdict

import numpy as np

ROOTS = {"GO:0008150", "GO:0003674", "GO:0005575"}
NS = {"biological_process": "BP", "molecular_function": "MF", "cellular_component": "CC"}
THRESH = np.round(np.arange(1, 100) / 100.0, 2)


class GO:
    def __init__(self, obo):
        self.parents = defaultdict(set)
        self.ns = {}
        self.alt = {}
        self.obsolete = set()
        cur = None
        with open(obo) as f:
            for line in f:
                line = line.rstrip("\n")
                if line == "[Term]":
                    cur = {"id": None}
                    continue
                if line.startswith("[") and line.endswith("]"):
                    cur = None
                    continue
                if cur is None or not line:
                    continue
                k, _, v = line.partition(": ")
                if k == "id":
                    cur["id"] = v
                elif k == "namespace":
                    self.ns[cur["id"]] = NS.get(v, "")
                elif k == "alt_id":
                    self.alt[v] = cur["id"]
                elif k == "is_a":
                    self.parents[cur["id"]].add(v.split()[0])
                elif k == "relationship":
                    p = v.split()
                    if p[0] == "part_of":
                        self.parents[cur["id"]].add(p[1])
                elif k == "is_obsolete" and v == "true":
                    self.obsolete.add(cur["id"])
        self._anc = {}

    def norm(self, t):
        t = self.alt.get(t, t)
        if t in self.obsolete or t not in self.ns:
            return None
        return t

    def anc(self, t):
        a = self._anc.get(t)
        if a is not None:
            return a
        out = {t}
        stack = [t]
        while stack:
            x = stack.pop()
            for p in self.parents.get(x, ()):
                if p not in out:
                    out.add(p)
                    stack.append(p)
        a = frozenset(out - ROOTS)
        self._anc[t] = a
        return a


def load_truth(path, go, idmap=None):
    """Return {pid: set(terms)} propagated, roots removed. Columns: id, aspect, term."""
    raw = defaultdict(set)
    with open(path) as f:
        next(f, None)
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 3:
                continue
            t = go.norm(p[2])
            if t is None:
                continue
            raw[p[0]].add(t)
    out = {}
    for pid, ts in raw.items():
        s = set()
        for t in ts:
            s |= go.anc(t)
        if s:
            out[pid] = s
    return out


def load_ec_truth(path):
    out = defaultdict(set)
    with open(path) as f:
        next(f, None)
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 3 and p[2].startswith("EC:") and "-" not in p[2]:
                out[p[0]].add(p[2])
    return dict(out)


def iter_pred_rows(path, ann_type="GO"):
    """Yield (pid, term, score) from either the GSPA integrated schema or the
    neural-sidecar schema (protein_id, term, score, annotation_type)."""
    with open(path) as f:
        hdr = next(f).rstrip("\n").split("\t")
        if "posterior_prob" in hdr:
            ip, it, ifn, isc = hdr.index("protein_id"), hdr.index("type"), hdr.index("function_id"), hdr.index("posterior_prob")
        else:
            ip, it, ifn, isc = hdr.index("protein_id"), hdr.index("annotation_type"), hdr.index("term"), hdr.index("score")
        n = max(ip, it, ifn, isc)
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) <= n or p[it] != ann_type:
                continue
            try:
                s = float(p[isc])
            except ValueError:
                continue
            yield p[ip], p[ifn], s


def load_preds(path, go, keep, ann_type="GO", min_score=0.01, idmap=None):
    """{pid: {term: score}} for pids in `keep`, propagated (max), roots removed.
    idmap: optional {pred_id: [eval_id, ...]}."""
    raw = defaultdict(dict)
    for pid, term, s in iter_pred_rows(path, ann_type):
        if s < min_score:
            continue
        ids = idmap.get(pid, ()) if idmap is not None else (pid,)
        for q in ids:
            if q not in keep:
                continue
            d = raw[q]
            if s > d.get(term, -1.0):
                d[term] = s
    if ann_type != "GO":
        return raw
    out = {}
    for pid, d in raw.items():
        prop = {}
        for t, s in d.items():
            t2 = go.norm(t)
            if t2 is None:
                continue
            for a in go.anc(t2):
                if s > prop.get(a, -1.0):
                    prop[a] = s
        out[pid] = prop
    return out


def compute_ic(truth):
    n = max(1, len(truth))
    cnt = defaultdict(int)
    for ts in truth.values():
        for t in ts:
            cnt[t] += 1
    return {t: -math.log2(c / n) for t, c in cnt.items()}


def arrays(preds, truth, ic):
    """Flatten into numpy arrays over truth proteins."""
    pids = sorted(truth)
    pidx = {p: i for i, p in enumerate(pids)}
    ntrue = np.array([len(truth[p]) for p in pids], dtype=np.float64)
    ictrue = np.array([sum(ic.get(t, 0.0) for t in truth[p]) for p in pids], dtype=np.float64)
    P, S, T, I = [], [], [], []
    for p, d in preds.items():
        i = pidx.get(p)
        if i is None:
            continue
        tr = truth[p]
        for t, s in d.items():
            P.append(i)
            S.append(s)
            T.append(t in tr)
            I.append(ic.get(t, 0.0))
    return {"n": len(pids), "ntrue": ntrue, "ictrue": ictrue,
            "pid": np.array(P, dtype=np.int64), "score": np.array(S, dtype=np.float64),
            "tp": np.array(T, dtype=bool), "ic": np.array(I, dtype=np.float64)}


def concat(arrs):
    off = 0
    out = {"n": 0, "ntrue": [], "ictrue": [], "pid": [], "score": [], "tp": [], "ic": []}
    for a in arrs:
        out["ntrue"].append(a["ntrue"]); out["ictrue"].append(a["ictrue"])
        out["pid"].append(a["pid"] + off); out["score"].append(a["score"])
        out["tp"].append(a["tp"]); out["ic"].append(a["ic"])
        off += a["n"]
    out["n"] = off
    for k in ("ntrue", "ictrue", "score", "ic"):
        out[k] = np.concatenate(out[k]) if out[k] else np.zeros(0)
    out["pid"] = np.concatenate(out["pid"]).astype(np.int64) if out["pid"] else np.zeros(0, dtype=np.int64)
    out["tp"] = np.concatenate(out["tp"]).astype(bool) if out["tp"] else np.zeros(0, dtype=bool)
    return out


def metrics(a, with_smin=True):
    n = a["n"]
    if n == 0:
        return None
    ntrue_tot = a["ntrue"].sum()
    best_micro = (0.0, 0.0, 0.0, 0.0)
    best_cafa = (0.0, 0.0, 0.0, 0.0)
    best_s = (float("inf"), 0.0, 0.0, 0.0)
    pid, score, tp, ic = a["pid"], a["score"], a["tp"], a["ic"]
    for t in THRESH:
        m = score >= t
        npred = np.bincount(pid[m], minlength=n).astype(np.float64)
        ntp = np.bincount(pid[m & tp], minlength=n).astype(np.float64)
        if True:
            TP = ntp.sum(); FP = npred.sum() - TP; FN = ntrue_tot - TP
            pr = TP / (TP + FP) if TP + FP else 0.0
            rc = TP / (TP + FN) if TP + FN else 0.0
            f = 2 * pr * rc / (pr + rc) if pr + rc else 0.0
            if f > best_micro[0]:
                best_micro = (f, t, pr, rc)
            has = npred > 0
            ap = (ntp[has] / npred[has]).mean() if has.any() else 0.0
            ar = (ntp / a["ntrue"]).mean()
            fc = 2 * ap * ar / (ap + ar) if ap + ar else 0.0
            if fc > best_cafa[0]:
                best_cafa = (fc, t, ap, ar)
        if with_smin:
            ic_tp = np.bincount(pid[m & tp], weights=ic[m & tp], minlength=n)
            ic_fp = np.bincount(pid[m & ~tp], weights=ic[m & ~tp], minlength=n)
            ru = (a["ictrue"] - ic_tp).mean()
            mi = ic_fp.mean()
            s = math.sqrt(ru * ru + mi * mi)
            if s < best_s[0]:
                best_s = (s, t, ru, mi)
    m5 = score >= 0.5
    cov = (np.bincount(pid[m5], minlength=n) > 0).mean()
    ictot = a["ictrue"].sum()
    icr = (ic[m5 & tp].sum() / ictot) if ictot > 0 else float("nan")
    return {"n_truth_proteins": int(n), "n_truth_pairs": int(ntrue_tot),
            "fmax_micro": best_micro[0], "fmax_micro_t": best_micro[1],
            "fmax_cafa": best_cafa[0], "fmax_cafa_t": best_cafa[1],
            "coverage_05": float(cov), "ic_recall_05": float(icr),
            "smin": best_s[0] if with_smin else float("nan"), "smin_t": best_s[1]}


def load_idmap_uniprot_to_refseq(path):
    m = defaultdict(list)
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2 and p[0] != "refseq_id":
                m[p[1]].append(p[0])
    return m
