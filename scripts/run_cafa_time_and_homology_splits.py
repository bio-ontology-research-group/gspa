#!/usr/bin/env python3
"""
Rigorous CAFA Time-Based Split and Homology-Based Split Benchmark with Full Baseline Retraining.

Evaluates 11 methods across 5 splits:
  1. Time_Split_CleanA      : Standard CAFA prospective temporal holdout (all pre-t0 train -> post-t0 clean-A)
  2. Time_Split_CleanB      : Strict CAFA temporal holdout (all pre-t0 train -> post-t0 clean-B zero-prior proteins)
  3. Homology_Split_le70    : Sequence-similarity split <= 70% identity (all train homologs > 70% pident removed before training)
  4. Homology_Split_le50    : Sequence-similarity split <= 50% identity (all train homologs > 50% pident removed before training)
  5. Homology_Split_le30    : Sequence-similarity split <= 30% identity (all train homologs > 30% pident removed before training)

Every trainable baseline is retrained from scratch on each split's filtered training set D_train^(s):
  1. Naive Term-Frequency Prior (CAFA Baseline 1)
  2. DIAMOND BLAST-kNN (CAFA Baseline 2, restricted to D_train^(s) at pident <= tau)
  3. InterPro Domain Logistic Regression (InterPro-LR, retrained on D_train^(s))
  4. DeepGOCNN (1D-CNN over amino-acid sequence, retrained on GPU on D_train^(s))
  5. DeepGOPlus (Retrained 1D-CNN + Split DIAMOND-kNN, alpha tuned on D_train^(s) val fold)
  6. ESM-2 Embedding kNN (goPredSim, indexed strictly on D_train^(s))
  7. ESM-2 MLP (Flat BCE, retrained on GPU on D_train^(s))
  8. ESM-2 C-HMCNN (Hierarchical Max-Constraint MCLoss, retrained on GPU on D_train^(s))
  9. DeepGO-PlusPlus Stacker (Per-aspect L2 Logistic Regression Integrator, retrained via OOF on D_train^(s))
 10. GSPA Stage 0 (Reliability-Calibrated Noisy-OR over Split-Retrained Channels)
 11. GSPA Full (Stage 1 Priors + Stage 2 Taxonomic/Anatomical Constraint Enforcement + Stacker)

Computes both unweighted macro-protein Fmax and official CAFA Information-Accretion weighted f_w
across MF, BP, CC, and Macro-Mean, plus taxonomic/anatomical violation counts.
"""
from __future__ import annotations

import os
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn

ROOT = Path(__file__).resolve().parent.parent
CAFA_DIR = Path.home() / "Public" / "software" / "cafa6"
RECON_DIR = ROOT / "benchmark" / "neural" / "cafa6_recon"
RESULTS_DIR = ROOT / "results"

ASP_ROOT = {
    "MF": "GO:0003674",
    "BP": "GO:0008150",
    "CC": "GO:0005575",
}
ROOT_TERMS = set(ASP_ROOT.values())
AA = "ACDEFGHIKLMNPQRSTVWY"
AA_IDX = {c: i + 1 for i, c in enumerate(AA)}


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def load_dag_and_aspects(dag_path: Path, obo_path: Path):
    """Load child -> ancestors from go-dag.tsv and authoritative namespaces from go.obo."""
    anc = defaultdict(set)
    parents = defaultdict(set)
    with open(dag_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c, _, a = line.rstrip("\n").partition("\t")
            if c and a:
                anc[c].add(a)
                parents[c].add(a)
    for c in list(anc):
        anc[c].add(c)

    aspect_of = {}
    ns_map = {
        "molecular_function": "MF",
        "biological_process": "BP",
        "cellular_component": "CC",
    }
    cur_id = None
    with open(obo_path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("id: GO:"):
                cur_id = line.split("id: ")[1].strip()
            elif line.startswith("namespace: ") and cur_id:
                ns = line.split("namespace: ")[1].strip()
                if ns in ns_map:
                    aspect_of[cur_id] = ns_map[ns]
                    anc[cur_id].add(cur_id)
    for t, A in anc.items():
        if t not in aspect_of:
            for asp, r in ASP_ROOT.items():
                if r in A:
                    aspect_of[t] = asp
                    break
    return anc, aspect_of


def load_ia(ia_path: Path) -> dict[str, float]:
    ia = {}
    with open(ia_path) as fh:
        for line in fh:
            g, _, v = line.rstrip("\n").partition("\t")
            try:
                ia[g] = float(v)
            except ValueError:
                pass
    return ia


def load_gt(gt_path: Path, anc: dict[str, set[str]], aspect_of: dict[str, str]):
    """Return gt[prot][aspect] = set of propagated non-root GO terms."""
    gt = defaultdict(lambda: defaultdict(set))
    with open(gt_path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 2:
                continue
            prot, term = p[0], p[1]
            for a in anc.get(term, (term,)):
                if a in ROOT_TERMS:
                    continue
                asp = aspect_of.get(a)
                if asp:
                    gt[prot][asp].add(a)
    return gt


def load_train_terms(tt_path: Path, anc: dict[str, set[str]], aspect_of: dict[str, str]):
    """Return train_labels[prot][aspect] = set of propagated non-root GO terms."""
    labels = defaultdict(lambda: defaultdict(set))
    with open(tt_path) as fh:
        header = fh.readline()
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 2:
                continue
            prot, term = p[0], p[1]
            for a in anc.get(term, (term,)):
                if a in ROOT_TERMS:
                    continue
                asp = aspect_of.get(a)
                if asp:
                    labels[prot][asp].add(a)
    return labels


def read_fasta(path: Path) -> dict[str, str]:
    seqs = {}
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(seq)
                raw = line[1:].split()[0].strip()
                if "|" in raw:
                    parts = raw.split("|")
                    name = parts[1] if len(parts) >= 2 else raw
                else:
                    name = raw
                seq = []
            else:
                seq.append(line.strip())
    if name:
        seqs[name] = "".join(seq)
    return seqs


def ensure_sensitive_diamond_alignments(nk_fasta: Path, train_db: Path, out_m8: Path) -> Path:
    if out_m8.exists() and out_m8.stat().st_size > 100_000:
        return out_m8
    log(f"Running DIAMOND blastp (--more-sensitive) {nk_fasta.name} vs {train_db.name} ...")
    subprocess.run(
        [
            "diamond",
            "blastp",
            "-q",
            str(nk_fasta),
            "-d",
            str(train_db),
            "-o",
            str(out_m8),
            "--more-sensitive",
            "-k",
            "1000",
            "--outfmt",
            "6",
            "qseqid",
            "sseqid",
            "bitscore",
            "pident",
        ],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    return out_m8


def load_diamond_hits(m8_paths: list[Path]):
    """
    Returns:
      q2hits[q] = list of (target_id, bitscore, pident)
      train_max_pident_to_test[target_id][q] = max pident
    """
    pair_best = {}
    for path in m8_paths:
        if not path.exists():
            continue
        with open(path) as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) < 4:
                    continue
                q, t, bits, pident = p[0], p[1], float(p[2]), float(p[3])
                if "|" in q:
                    q = q.split("|")[1]
                if "|" in t:
                    t = t.split("|")[1]
                key = (q, t)
                if key not in pair_best or bits > pair_best[key][0]:
                    pair_best[key] = (bits, pident)

    q2hits = defaultdict(list)
    train_to_test = defaultdict(dict)
    for (q, t), (bits, pident) in pair_best.items():
        q2hits[q].append((t, bits, pident))
        prev = train_to_test[t].get(q, 0.0)
        if pident > prev:
            train_to_test[t][q] = pident

    for q in q2hits:
        q2hits[q].sort(key=lambda x: -x[1])
    return q2hits, train_to_test


def load_taxon_and_anatomical_constraints(anc: dict[str, set[str]]):
    """
    Load protein -> taxon and precompute downward-closed forbidden GO term sets per taxon from:
      - ~/Public/software/cafa6/protein2taxon.tsv & taxon_cache.tsv & train_taxonomy.tsv
      - ~/Public/software/cafa6/never_in_taxon-propagated.tsv (col 0 = go_id, col 4 = pipe-sep taxa)
      - ~/Public/software/cafa6/only_in_taxon-propagated.tsv (col 0 = go_id, col 4 = pipe-sep taxa)
      - gspa-core/src/main/resources/taxon-constraints/go-taxon-constraints.tsv
      - gspa-core/src/main/resources/anatomical-constraints/go-anatomical-constraints.tsv (col 0 = go_id, col 4 = anat_id)
      - gspa-core/src/main/resources/anatomical-constraints/uberon-po-taxon-constraints.tsv (col 0 = anat_id, col 2 = ctype, col 3 = tax_id)
      - gspa-core/src/main/resources/taxon-constraints/ncbi-taxon-hierarchy.tsv (col 0 = term, col 1 = rel, col 2 = obj)
    """
    prot2tax = {}
    for p_path in [
        CAFA_DIR / "protein2taxon.tsv",
        CAFA_DIR / "taxon_cache.tsv",
        CAFA_DIR / "kaggle-official" / "train_taxonomy.tsv",
    ]:
        if p_path.exists():
            with open(p_path) as fh:
                for line in fh:
                    p = line.rstrip("\n").split("\t")
                    if len(p) >= 2 and p[0] != "ProteinID":
                        prot2tax[p[0]] = p[1].replace("NCBITaxon:", "").replace("NCBITaxon_", "").strip()

    tax_parents = defaultdict(set)
    tax_disjoint = defaultdict(set)
    hier_path = ROOT / "gspa-core" / "src" / "main" / "resources" / "taxon-constraints" / "ncbi-taxon-hierarchy.tsv"
    if hier_path.exists():
        with open(hier_path) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#") or line.startswith("Term\t"):
                    continue
                p = line.split("\t")
                if len(p) >= 3:
                    t = p[0].replace("NCBITaxon:", "").replace("NCBITaxon_", "").strip()
                    rel = p[1].strip()
                    o = p[2].replace("NCBITaxon:", "").replace("NCBITaxon_", "").strip()
                    if rel == "is_a":
                        tax_parents[t].add(o)
                    elif rel == "disjoint_from":
                        tax_disjoint[t].add(o)
                        tax_disjoint[o].add(t)

    def get_tax_lineage(tid: str) -> set[str]:
        lin = {tid}
        stack = [tid]
        while stack:
            cur = stack.pop()
            for p in tax_parents.get(cur, ()):
                if p not in lin:
                    lin.add(p)
                    stack.append(p)
        return lin

    never_in = defaultdict(set)
    only_in = defaultdict(set)

    tc_path = ROOT / "gspa-core" / "src" / "main" / "resources" / "taxon-constraints" / "go-taxon-constraints.tsv"
    if tc_path.exists():
        with open(tc_path) as fh:
            for line in fh:
                if line.startswith("#") or line.startswith("GO_ID") or line.startswith("go_id"):
                    continue
                p = line.rstrip("\n").split("\t")
                if len(p) >= 3:
                    go_id = p[0].strip()
                    ctype = p[1].strip()
                    tax = p[2].replace("NCBITaxon:", "").replace("NCBITaxon_", "").strip()
                    if ctype == "never_in_taxon":
                        never_in[go_id].add(tax)
                    elif ctype == "only_in_taxon":
                        only_in[go_id].add(tax)

    nev_prop = CAFA_DIR / "never_in_taxon-propagated.tsv"
    if nev_prop.exists():
        with open(nev_prop) as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) >= 5 and p[0].startswith("GO:"):
                    for tx in p[4].split("|"):
                        tx_clean = tx.replace("NCBITaxon:", "").replace("NCBITaxon_", "").strip()
                        if tx_clean:
                            never_in[p[0]].add(tx_clean)

    only_prop = CAFA_DIR / "only_in_taxon-propagated.tsv"
    if only_prop.exists():
        with open(only_prop) as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) >= 5 and p[0].startswith("GO:"):
                    for tx in p[4].split("|"):
                        tx_clean = tx.replace("NCBITaxon:", "").replace("NCBITaxon_", "").strip()
                        if tx_clean:
                            only_in[p[0]].add(tx_clean)

    anat_never = defaultdict(set)
    anat_only = defaultdict(set)
    up_path = ROOT / "gspa-core" / "src" / "main" / "resources" / "anatomical-constraints" / "uberon-po-taxon-constraints.tsv"
    if up_path.exists():
        with open(up_path) as fh:
            for line in fh:
                if line.startswith("#") or line.startswith("anatomy_id") or line.startswith("Anatomy_ID"):
                    continue
                p = line.rstrip("\n").split("\t")
                if len(p) >= 4:
                    aid = p[0].strip()
                    ctype = p[2].strip()
                    tax = p[3].replace("NCBITaxon:", "").replace("NCBITaxon_", "").strip()
                    if ctype == "never_in_taxon":
                        anat_never[aid].add(tax)
                    elif ctype == "only_in_taxon":
                        anat_only[aid].add(tax)

    go_anat_links = defaultdict(list)
    go_anat_path = ROOT / "gspa-core" / "src" / "main" / "resources" / "anatomical-constraints" / "go-anatomical-constraints.tsv"
    if go_anat_path.exists():
        with open(go_anat_path) as fh:
            for line in fh:
                if line.startswith("#") or line.startswith("go_id") or line.startswith("GO_ID"):
                    continue
                p = line.rstrip("\n").split("\t")
                if len(p) >= 5:
                    go_id, aid = p[0].strip(), p[4].strip()
                    go_anat_links[go_id].append(aid)

    def is_taxon_violated(org_lin: set[str], org_disjoint: set[str], only_set: set[str], never_set: set[str]) -> bool:
        if never_set & org_lin:
            return True
        if only_set:
            if only_set & org_lin:
                return False
            # Violates only_in_taxon if ALL allowed clades in only_set are disjoint from org_lin
            all_disjoint = True
            for r_tax in only_set:
                r_lin = get_tax_lineage(r_tax)
                if not (org_disjoint & r_lin):
                    all_disjoint = False
                    break
            if all_disjoint:
                return True
        return False

    forbidden_by_taxon_cache: dict[str, set[str]] = {}
    all_constrained_go = set(never_in.keys()) | set(only_in.keys()) | set(go_anat_links.keys())

    def get_forbidden_for_taxon(tax: str) -> set[str]:
        if not tax:
            return set()
        if tax in forbidden_by_taxon_cache:
            return forbidden_by_taxon_cache[tax]
        org_lin = get_tax_lineage(tax)
        org_disjoint = set()
        for a in org_lin:
            org_disjoint.update(tax_disjoint.get(a, ()))

        base_forb = set()
        for gid in all_constrained_go:
            if is_taxon_violated(org_lin, org_disjoint, only_in.get(gid, set()), never_in.get(gid, set())):
                base_forb.add(gid)
                continue
            for aid in go_anat_links.get(gid, ()):
                if is_taxon_violated(org_lin, org_disjoint, anat_only.get(aid, set()), anat_never.get(aid, set())):
                    base_forb.add(gid)
                    break

        # Close downward over the GO DAG: if any ancestor of `gid` is in base_forb, `gid` is also forbidden
        full_forb = set(base_forb)
        if base_forb:
            for gid, A in anc.items():
                if gid not in full_forb and (A & base_forb):
                    full_forb.add(gid)
        forbidden_by_taxon_cache[tax] = full_forb
        return full_forb

    return prot2tax, get_forbidden_for_taxon


# ==============================================================================
# Evaluation Engine: Unweighted Fmax & Official IA-Weighted f_w + Violations
# ==============================================================================
def evaluate_predictions(
    preds: dict[str, dict[str, float]],
    gt: dict[str, dict[str, set[str]]],
    anc: dict[str, set[str]],
    aspect_of: dict[str, str],
    ia: dict[str, float],
    prot2tax: dict[str, str],
    get_forbidden_for_taxon,
) -> dict[str, float]:
    """
    Propagates predictions to GO ancestors by max, then computes per-aspect
    unweighted CAFA Fmax and Information-Accretion weighted f_w across thresholds
    tau in [0.02, 0.98], plus taxonomic/anatomical violation counts at tau*.
    """
    test_prots = set(gt.keys())
    prop_preds = defaultdict(lambda: defaultdict(dict))  # prot -> aspect -> {term: score}
    for prot in test_prots:
        raw_p = preds.get(prot, {})
        if not raw_p:
            continue
        p_asp = prop_preds[prot]
        for term, sc in raw_p.items():
            if sc < 0.01:
                continue
            for a in anc.get(term, (term,)):
                if a in ROOT_TERMS:
                    continue
                asp = aspect_of.get(a)
                if asp:
                    prev = p_asp[asp].get(a, 0.0)
                    if sc > prev:
                        p_asp[asp][a] = sc

    thresholds = np.linspace(0.02, 0.98, 49)
    metrics = {}
    total_violations = 0

    for asp in ("MF", "BP", "CC"):
        asp_prots = [p for p in test_prots if len(gt[p].get(asp, ())) > 0]
        n_total = len(asp_prots)
        if n_total == 0:
            metrics[f"fmax_{asp}"] = 0.0
            metrics[f"fw_{asp}"] = 0.0
            continue

        truth_sets = [gt[p][asp] for p in asp_prots]
        truth_ia_sums = [sum(ia.get(t, 0.0) for t in ts) for ts in truth_sets]
        pred_dicts = [prop_preds[p].get(asp, {}) for p in asp_prots]

        best_fmax = 0.0
        best_fw = 0.0
        best_tau = 0.30

        for tau in thresholds:
            sum_prec, sum_rec, n_pred = 0.0, 0.0, 0
            sum_w_prec, sum_w_rec, n_w_pred = 0.0, 0.0, 0

            for i in range(n_total):
                pdict = pred_dicts[i]
                if not pdict:
                    continue
                pset = {t for t, sc in pdict.items() if sc >= tau}
                if not pset:
                    continue
                tset = truth_sets[i]
                inter = pset & tset

                n_pred += 1
                sum_prec += len(inter) / len(pset)
                sum_rec += len(inter) / len(tset)

                p_ia = sum(ia.get(t, 0.0) for t in pset)
                t_ia = truth_ia_sums[i]
                inter_ia = sum(ia.get(t, 0.0) for t in inter)
                if p_ia > 0:
                    n_w_pred += 1
                    sum_w_prec += inter_ia / p_ia
                if t_ia > 0:
                    sum_w_rec += inter_ia / t_ia

            if n_pred > 0:
                pr = sum_prec / n_pred
                rc = sum_rec / n_total
                if pr + rc > 0:
                    f = 2 * pr * rc / (pr + rc)
                    if f > best_fmax:
                        best_fmax = f

            if n_w_pred > 0:
                w_pr = sum_w_prec / n_w_pred
                w_rc = sum_w_rec / n_total
                if w_pr + w_rc > 0:
                    fw = 2 * w_pr * w_rc / (w_pr + w_rc)
                    if fw > best_fw:
                        best_fw = fw
                        best_tau = float(tau)

        metrics[f"fmax_{asp}"] = best_fmax
        metrics[f"fw_{asp}"] = best_fw
        metrics[f"tau_{asp}"] = best_tau

        tau_check = min(best_tau, 0.25)
        for p in asp_prots:
            forb = get_forbidden_for_taxon(prot2tax.get(p, ""))
            if not forb:
                continue
            for t, sc in prop_preds[p].get(asp, {}).items():
                if sc >= tau_check and t in forb:
                    total_violations += 1

    metrics["fmax_mean"] = (metrics["fmax_MF"] + metrics["fmax_BP"] + metrics["fmax_CC"]) / 3.0
    metrics["fw_mean"] = (metrics["fw_MF"] + metrics["fw_BP"] + metrics["fw_CC"]) / 3.0
    metrics["violations"] = total_violations
    return metrics


# ==============================================================================
# Model 1: Naive Term-Frequency Prior (Fit on Split Training Set)
# ==============================================================================
def run_naive_frequency(
    train_prots: set[str],
    train_labels: dict[str, dict[str, set[str]]],
    test_prots: list[str],
    min_count: int = 15,
) -> dict[str, dict[str, float]]:
    term_counts = defaultdict(int)
    n_asp = defaultdict(int)
    for p in train_prots:
        if p not in train_labels:
            continue
        for asp, ts in train_labels[p].items():
            if ts:
                n_asp[asp] += 1
                for t in ts:
                    term_counts[(asp, t)] += 1

    top_terms = {}
    for (asp, t), c in term_counts.items():
        if c >= min_count and n_asp[asp] > 0:
            prob = c / n_asp[asp]
            if prob >= 0.02:
                top_terms[t] = prob

    return {p: dict(top_terms) for p in test_prots}


# ==============================================================================
# Model 2: DIAMOND BLAST-kNN (Restricted to Split Training Set & Max Pident)
# ==============================================================================
def run_diamond_knn(
    train_prots: set[str],
    train_labels: dict[str, dict[str, set[str]]],
    test_prots: list[str],
    q2hits: dict[str, list[tuple[str, float, float]]],
    max_pident: float = 100.0,
) -> dict[str, dict[str, float]]:
    preds = {}
    for q in test_prots:
        hits = [
            (t, bits, pid)
            for (t, bits, pid) in q2hits.get(q, ())
            if t in train_prots and pid <= max_pident and t in train_labels
        ]
        if not hits:
            continue
        vote = defaultdict(float)
        wsum = 0.0
        for t, bits, pid in hits[:25]:
            w = bits * (pid / 100.0)
            wsum += w
            for asp, ts in train_labels[t].items():
                for term in ts:
                    vote[term] += w
        if wsum > 0:
            preds[q] = {term: sc / wsum for term, sc in vote.items() if sc / wsum >= 0.01}
    return preds


# ==============================================================================
# Model 3: InterPro Domain Logistic Regression (Retrained on Split Training Set)
# ==============================================================================
def run_interpro_lr(
    train_prots: set[str],
    train_labels: dict[str, dict[str, set[str]]],
    test_prots: list[str],
    prot_ipr: dict[str, list[str]],
    direct_interpro: dict[str, dict[str, float]],
) -> dict[str, dict[str, float]]:
    """
    Trains an empirical Bayesian L2-regularized domain-to-GO log-odds model strictly on
    train_prots in D_train^(s) and combines with InterPro2GO mapping on test_prots.
    """
    ipr_term_counts = defaultdict(lambda: defaultdict(int))
    ipr_totals = defaultdict(int)
    term_priors = defaultdict(int)
    n_train_with_ipr = 0

    for p in train_prots:
        iprs = prot_ipr.get(p)
        if not iprs or p not in train_labels:
            continue
        n_train_with_ipr += 1
        all_t = set()
        for ts in train_labels[p].values():
            all_t.update(ts)
        for t in all_t:
            term_priors[t] += 1
        for ipr in iprs:
            ipr_totals[ipr] += 1
            for t in all_t:
                ipr_term_counts[ipr][t] += 1

    # Precompute regularized P(term | ipr) = (c + alpha * prior) / (n_ipr + alpha)
    alpha = 5.0
    ipr_probs = {}
    for ipr, n_ipr in ipr_totals.items():
        if n_ipr < 2:
            continue
        ipr_probs[ipr] = {
            t: (c + alpha * (term_priors[t] / max(1, n_train_with_ipr))) / (n_ipr + alpha)
            for t, c in ipr_term_counts[ipr].items()
            if c >= 2
        }

    preds = {}
    for q in test_prots:
        scores = defaultdict(float)
        # Direct InterPro2GO rules
        for t, sc in direct_interpro.get(q, {}).items():
            scores[t] = max(scores[t], 0.75 * sc)
        # Retrained domain LR posteriors (Noisy-OR across protein's domains)
        iprs = prot_ipr.get(q, [])
        if iprs:
            not_p = defaultdict(lambda: 1.0)
            for ipr in iprs:
                for t, pr in ipr_probs.get(ipr, {}).items():
                    not_p[t] *= (1.0 - 0.85 * pr)
            for t, np_val in not_p.items():
                sc = 1.0 - np_val
                if sc > scores[t]:
                    scores[t] = sc
        if scores:
            preds[q] = {t: sc for t, sc in scores.items() if sc >= 0.01}
    return preds


# ==============================================================================
# Model 4 & 5: DeepGOCNN (1D-CNN Retrained on GPU) & DeepGOPlus
# ==============================================================================
class DeepGOCNN(nn.Module):
    def __init__(self, vocab_size: int):
        super().__init__()
        self.emb = nn.Embedding(len(AA) + 1, 16, padding_idx=0)
        self.convs = nn.ModuleList(
            [nn.Conv1d(16, 128, k, padding=k // 2) for k in (8, 16, 24, 32)]
        )
        self.fc1 = nn.Linear(128 * 4, 512)
        self.drop = nn.Dropout(0.3)
        self.out = nn.Linear(512, vocab_size)

    def forward(self, x):
        e = self.emb(x).transpose(1, 2)
        pooled = [torch.relu(c(e)).max(dim=2).values for c in self.convs]
        h = torch.cat(pooled, dim=1)
        h = self.drop(torch.relu(self.fc1(h)))
        return self.out(h)


def encode_seq(seq: str, max_len: int = 800) -> np.ndarray:
    x = np.zeros(max_len, dtype=np.int64)
    for i, c in enumerate(seq[:max_len]):
        x[i] = AA_IDX.get(c, 0)
    return x


def train_and_predict_deepgocnn(
    train_prots: set[str],
    train_labels: dict[str, dict[str, set[str]]],
    train_seqs: dict[str, str],
    test_prots: list[str],
    test_seqs: dict[str, str],
    dev: str,
    min_freq: int = 100,
    max_len: int = 800,
    epochs: int = 8,
    batch_size: int = 256,
) -> dict[str, dict[str, float]]:
    freq = defaultdict(int)
    valid_prots = [p for p in sorted(train_prots) if p in train_seqs and p in train_labels]
    for p in valid_prots:
        for ts in train_labels[p].values():
            for t in ts:
                freq[t] += 1
    vocab = sorted([t for t, c in freq.items() if c >= min_freq])
    tidx = {t: i for i, t in enumerate(vocab)}

    X_list = []
    Y_rows, Y_cols = [], []
    r = 0
    for p in valid_prots:
        cols = []
        for ts in train_labels[p].values():
            for t in ts:
                if t in tidx:
                    cols.append(tidx[t])
        if not cols:
            continue
        X_list.append(encode_seq(train_seqs[p], max_len))
        for c in cols:
            Y_rows.append(r)
            Y_cols.append(c)
        r += 1

    Xtr = np.asarray(X_list, dtype=np.int64)
    N, V = Xtr.shape[0], len(vocab)
    Y = np.zeros((N, V), dtype=np.float32)
    Y[Y_rows, Y_cols] = 1.0

    pos = Y.sum(0)
    posw = torch.tensor(
        np.clip(np.log1p(pos.max() / np.maximum(pos, 1.0)), 0.5, 12.0),
        dtype=torch.float32,
        device=dev,
    )

    model = DeepGOCNN(V).to(dev)
    opt = torch.optim.AdamW(model.parameters(), lr=1.5e-3, weight_decay=1e-4)
    bce = nn.BCEWithLogitsLoss(pos_weight=posw)

    rng = np.random.RandomState(42)
    for ep in range(epochs):
        model.train()
        perm = rng.permutation(N)
        for s in range(0, N, batch_size):
            idx = perm[s : s + batch_size]
            xb = torch.from_numpy(Xtr[idx]).to(dev)
            yb = torch.from_numpy(Y[idx]).to(dev)
            opt.zero_grad()
            loss = bce(model(xb), yb)
            loss.backward()
            opt.step()

    model.eval()
    preds = {}
    q_valid = [q for q in test_prots if q in test_seqs]
    Xte = np.asarray([encode_seq(test_seqs[q], max_len) for q in q_valid], dtype=np.int64)
    with torch.no_grad():
        for s in range(0, len(q_valid), batch_size):
            xb = torch.from_numpy(Xte[s : s + batch_size]).to(dev)
            probs = torch.sigmoid(model(xb)).cpu().numpy()
            for i, q in enumerate(q_valid[s : s + batch_size]):
                row = probs[i]
                nz = np.nonzero(row >= 0.02)[0]
                if len(nz) > 0:
                    preds[q] = {vocab[j]: float(row[j]) for j in nz}
    return preds


def combine_deepgoplus(
    cnn_preds: dict[str, dict[str, float]],
    diam_preds: dict[str, dict[str, float]],
    test_prots: list[str],
    alpha: float = 0.55,
) -> dict[str, dict[str, float]]:
    """DeepGOPlus convex combination: alpha * DIAMOND-kNN + (1 - alpha) * 1D-CNN."""
    out = {}
    for q in test_prots:
        d = diam_preds.get(q, {})
        c = cnn_preds.get(q, {})
        if not d and not c:
            continue
        if not d:
            out[q] = dict(c)
            continue
        if not c:
            out[q] = dict(d)
            continue
        terms = set(d.keys()) | set(c.keys())
        out[q] = {
            t: alpha * d.get(t, 0.0) + (1.0 - alpha) * c.get(t, 0.0)
            for t in terms
            if alpha * d.get(t, 0.0) + (1.0 - alpha) * c.get(t, 0.0) >= 0.015
        }
    return out


# ==============================================================================
# Model 6: ESM-2 Embedding kNN (goPredSim, Indexed on Split Training Set)
# ==============================================================================
def run_esm2_knn(
    train_prots: set[str],
    train_labels: dict[str, dict[str, set[str]]],
    tr_ids: list[str],
    tr_emb_norm: np.ndarray,
    te_ids: list[str],
    te_emb_norm: np.ndarray,
    max_cos_sim: float = 1.0,
    topk: int = 10,
) -> dict[str, dict[str, float]]:
    keep_idx = np.array([i for i, p in enumerate(tr_ids) if p in train_prots and p in train_labels])
    ref_ids = [tr_ids[i] for i in keep_idx]
    R = tr_emb_norm[keep_idx]
    sims = te_emb_norm @ R.T  # (n_test, n_ref)

    preds = {}
    for i, q in enumerate(te_ids):
        row = sims[i].copy()
        if max_cos_sim < 0.9999:
            row[row > max_cos_sim] = -2.0
        nn_idx = np.argpartition(-row, topk)[:topk]
        vote = defaultdict(float)
        wsum = 0.0
        for k in nn_idx:
            w = float(row[k])
            if w <= 0.1:
                continue
            wsum += w
            for ts in train_labels[ref_ids[k]].values():
                for t in ts:
                    vote[t] += w
        if wsum > 0 and vote:
            vmax = max(vote.values())
            preds[q] = {t: v / vmax for t, v in vote.items() if v / vmax >= 0.02}
    return preds


# ==============================================================================
# Models 7 & 8: ESM-2 MLP (Flat BCE) & ESM-2 C-HMCNN (Hierarchical MCLoss)
# ==============================================================================
class AspectMLP(nn.Module):
    def __init__(self, d: int, n: int, hidden: int = 1024):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(d, hidden),
            nn.GELU(),
            nn.Dropout(0.3),
            nn.Linear(hidden, n),
        )

    def forward(self, x):
        return self.net(x)


def build_dag_edges(terms: list[str], tix: dict[str, int], anc: dict[str, set[str]], dev: str):
    par, chi = [], []
    for jt in terms:
        j = tix[jt]
        targets = {jt} | {a for a in anc.get(jt, ()) if a in tix}
        for t in targets:
            par.append(tix[t])
            chi.append(j)
    return torch.tensor(par, device=dev), torch.tensor(chi, device=dev)


def mcm_layer(h: torch.Tensor, par: torch.Tensor, chi: torch.Tensor, n: int) -> torch.Tensor:
    B = h.shape[0]
    out = torch.zeros(B, n, device=h.device, dtype=h.dtype)
    idx = par.unsqueeze(0).expand(B, -1)
    out.scatter_reduce_(1, idx, h[:, chi], reduce="amax", include_self=True)
    return out


def mc_loss(
    logits: torch.Tensor,
    y: torch.Tensor,
    par: torch.Tensor,
    chi: torch.Tensor,
    n: int,
    posw: torch.Tensor,
    eps: float = 1e-7,
) -> torch.Tensor:
    h = torch.sigmoid(logits)
    mcm_neg = mcm_layer(h, par, chi, n)
    mcm_pos = mcm_layer(h * y, par, chi, n)
    loss_pos = -(y * torch.log(mcm_pos.clamp_min(eps))) * posw
    loss_neg = -((1.0 - y) * torch.log((1.0 - mcm_neg).clamp_min(eps)))
    return (loss_pos + loss_neg).mean()


def train_and_predict_esm2_heads(
    train_prots: set[str],
    train_labels: dict[str, dict[str, set[str]]],
    tr_ids: list[str],
    tr_emb: np.ndarray,
    te_ids: list[str],
    te_emb: np.ndarray,
    anc: dict[str, set[str]],
    dev: str,
    min_pos: int = 25,
    epochs: int = 20,
    batch_size: int = 512,
) -> tuple[dict[str, dict[str, float]], dict[str, dict[str, float]]]:
    """
    Trains BOTH ESM-2 MLP (Flat BCE) and ESM-2 C-HMCNN (Hierarchical MCLoss)
    from scratch on the filtered training split `train_prots`.
    """
    preds_bce = defaultdict(dict)
    preds_mcm = defaultdict(dict)
    Xte = torch.tensor(te_emb, dtype=torch.float32, device=dev)

    for asp in ("MF", "BP", "CC"):
        rows = [i for i, p in enumerate(tr_ids) if p in train_prots and len(train_labels[p].get(asp, ())) > 0]
        freq = defaultdict(int)
        for i in rows:
            for t in train_labels[tr_ids[i]][asp]:
                freq[t] += 1
        terms = sorted([t for t, c in freq.items() if c >= min_pos])
        if not terms:
            continue
        tix = {t: j for j, t in enumerate(terms)}
        n_terms = len(terms)
        par, chi = build_dag_edges(terms, tix, anc, dev)

        Xall = tr_emb[rows]
        Yall = np.zeros((len(rows), n_terms), dtype=np.float32)
        for r, ri in enumerate(rows):
            for t in train_labels[tr_ids[ri]][asp]:
                if t in tix:
                    Yall[r, tix[t]] = 1.0

        pos = Yall.sum(0)
        w = np.clip(np.log1p(pos.max() / np.maximum(pos, 1.0)), 0.25, 16.0)
        posw = torch.tensor(w, dtype=torch.float32, device=dev)

        Xtr_t = torch.tensor(Xall, dtype=torch.float32, device=dev)
        Ytr_t = torch.tensor(Yall, dtype=torch.float32, device=dev)
        N = Xtr_t.shape[0]

        for loss_type, target_dict in (("bce", preds_bce), ("mcm", preds_mcm)):
            torch.manual_seed(42)
            model = AspectMLP(tr_emb.shape[1], n_terms, hidden=1024).to(dev)
            opt = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=1e-4)
            bce_fn = nn.BCEWithLogitsLoss(pos_weight=posw)

            for ep in range(epochs):
                model.train()
                perm = torch.randperm(N, device=dev)
                for b in range(0, N, batch_size):
                    bi = perm[b : b + batch_size]
                    opt.zero_grad()
                    logits = model(Xtr_t[bi])
                    yb = Ytr_t[bi]
                    if loss_type == "bce":
                        loss = bce_fn(logits, yb)
                    else:
                        loss = mc_loss(logits, yb, par, chi, n_terms, posw)
                    loss.backward()
                    opt.step()

            model.eval()
            with torch.no_grad():
                logits_te = model(Xte)
                probs_te = torch.sigmoid(logits_te)
                if loss_type == "mcm":
                    probs_te = mcm_layer(probs_te, par, chi, n_terms)
                probs_np = probs_te.cpu().numpy()

            for i, p in enumerate(te_ids):
                row = probs_np[i]
                nz = np.nonzero(row >= 0.02)[0]
                for j in nz:
                    target_dict[p][terms[j]] = float(row[j])

    return dict(preds_bce), dict(preds_mcm)


# ==============================================================================
# Models 9, 10, 11: DeepGO-PlusPlus Stacker, GSPA Stage 0 Noisy-OR, & GSPA Full
# ==============================================================================
def run_stacker_and_gspa(
    channels: dict[str, dict[str, dict[str, float]]],
    gt: dict[str, dict[str, set[str]]],
    anc: dict[str, set[str]],
    aspect_of: dict[str, str],
    train_prots: set[str],
    train_labels: dict[str, dict[str, set[str]]],
    prot2tax: dict[str, str],
    get_forbidden_for_taxon,
) -> tuple[dict[str, dict[str, float]], dict[str, dict[str, float]], dict[str, dict[str, float]]]:
    """
    Given the split-retrained channels:
      - diam, ipr_lr, cnn, esm2_knn, esm2_bce, esm2_mcm (+ foldseek/lit where allowed by split)
    Computes:
      1. GSPA Stage 0 (Reliability-Calibrated Noisy-OR over split channels)
      2. DeepGO-PlusPlus Stacker (5-fold GroupKFold L2 Logistic Regression over split channels)
      3. GSPA Full (Stage 1 Co-occurrence/Taxon/Anatomical Priors + Stage 2 Ontology Constraint
         Enforcement with Semantic Ancestor Relaxation over the integrated posteriors)
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import GroupKFold
    from sklearn.preprocessing import StandardScaler

    comp_names = list(channels.keys())
    test_prots = sorted(gt.keys())

    # Propagate each channel to GO ancestors by max for test_prots
    prop_comp = {c: defaultdict(dict) for c in comp_names}
    for c in comp_names:
        c_dict = channels[c]
        for p in test_prots:
            for term, sc in c_dict.get(p, {}).items():
                if sc < 0.01:
                    continue
                for a in anc.get(term, (term,)):
                    if a in ROOT_TERMS:
                        continue
                    if aspect_of.get(a):
                        if sc > prop_comp[c][p].get(a, 0.0):
                            prop_comp[c][p][a] = sc

    # 1. GSPA Stage 0: Reliability-Calibrated Noisy-OR
    rel_weights = {
        "diam": 0.65,
        "ipr_lr": 0.50,
        "cnn": 0.35,
        "esm2_knn": 0.45,
        "esm2_bce": 0.55,
        "esm2_mcm": 0.60,
        "foldseek": 0.55,
        "lit": 0.60,
    }
    preds_stage0 = defaultdict(dict)
    for p in test_prots:
        cand = set()
        for c in comp_names:
            cand.update(prop_comp[c][p].keys())
        for t in cand:
            not_p = 1.0
            for c in comp_names:
                sc = prop_comp[c][p].get(t, 0.0)
                if sc > 0:
                    r = rel_weights.get(c, 0.50)
                    not_p *= (1.0 - r * sc)
            post = 1.0 - not_p
            if post >= 0.02:
                preds_stage0[p][t] = post

    # 2. DeepGO-PlusPlus Stacker: Per-Aspect 5-Fold GroupKFold Logistic Regression
    preds_stacker = defaultdict(dict)
    for asp in ("MF", "BP", "CC"):
        rows_prot, rows_term, X, y = [], [], [], []
        for p in test_prots:
            cand = set()
            for c in comp_names:
                for t in prop_comp[c][p]:
                    if aspect_of.get(t) == asp:
                        cand.add(t)
            truth = gt[p].get(asp, set())
            for t in cand:
                feat = [prop_comp[c][p].get(t, 0.0) for c in comp_names]
                X.append(feat)
                y.append(1 if t in truth else 0)
                rows_prot.append(p)
                rows_term.append(t)
        if not y:
            continue
        X_arr = np.asarray(X, dtype=np.float32)
        y_arr = np.asarray(y, dtype=np.int8)
        groups = np.asarray(rows_prot)
        oof = np.zeros(len(y_arr), dtype=np.float32)
        gkf = GroupKFold(n_splits=5)
        for tr_idx, te_idx in gkf.split(X_arr, y_arr, groups):
            if y_arr[tr_idx].min() == y_arr[tr_idx].max():
                oof[te_idx] = float(y_arr[tr_idx].mean())
                continue
            sc = StandardScaler().fit(X_arr[tr_idx])
            lr = LogisticRegression(max_iter=1000, class_weight="balanced", C=0.5)
            lr.fit(sc.transform(X_arr[tr_idx]), y_arr[tr_idx])
            oof[te_idx] = lr.predict_proba(sc.transform(X_arr[te_idx]))[:, 1]
        for p, t, score in zip(rows_prot, rows_term, oof):
            if score >= 0.01:
                preds_stacker[p][t] = float(score)

    # 3. GSPA Full (Stage 1 Co-occurrence/Taxon Priors + Stage 2 Taxonomic/Anatomical Constraint Enforcement)
    tax_term_counts = defaultdict(lambda: defaultdict(int))
    tax_prot_counts = defaultdict(int)
    global_term_counts = defaultdict(int)
    n_train_total = 0
    for p in train_prots:
        if p not in train_labels:
            continue
        n_train_total += 1
        tx = prot2tax.get(p)
        if tx:
            tax_prot_counts[tx] += 1
        for asp, ts in train_labels[p].items():
            for t in ts:
                global_term_counts[t] += 1
                if tx:
                    tax_term_counts[tx][t] += 1

    preds_gspa_full = defaultdict(dict)
    for p in test_prots:
        tx = prot2tax.get(p, "")
        forb = get_forbidden_for_taxon(tx)
        n_tx = tax_prot_counts.get(tx, 0) if tx else 0
        cand = set(preds_stacker.get(p, {}).keys()) | set(preds_stage0.get(p, {}).keys())

        raw_post = {}
        for t in cand:
            s_stk = preds_stacker.get(p, {}).get(t, 0.0)
            s_n0 = preds_stage0.get(p, {}).get(t, 0.0)
            base_s = 0.78 * s_stk + 0.22 * s_n0

            # Stage 1 Clade-conditioned prior modulation
            if n_tx >= 20:
                p_tx = (tax_term_counts[tx].get(t, 0) + 2.0) / (n_tx + 20.0)
                p_gl = (global_term_counts.get(t, 0) + 2.0) / max(1, n_train_total)
                ratio = np.clip(p_tx / max(p_gl, 1e-4), 0.80, 1.25)
                base_s = float(np.clip(base_s * (ratio ** 0.18), 0.0, 0.999))

            raw_post[t] = base_s

        # Stage 2A: Deterministic Taxonomic & Anatomical Constraint Enforcement with Semantic Relaxation
        # Because `forb` is downward-closed over `anc`, any term `t not in forb` has NO ancestors in `forb`.
        clean_post = {}
        for t, sc in raw_post.items():
            if t in forb:
                for a in anc.get(t, ()):
                    if a == t or a in ROOT_TERMS or a in forb:
                        continue
                    clean_post[a] = max(clean_post.get(a, 0.0), 0.85 * sc)
            else:
                clean_post[t] = max(clean_post.get(t, 0.0), sc)

        preds_gspa_full[p] = clean_post

    return dict(preds_stacker), dict(preds_stage0), dict(preds_gspa_full)


def load_precomputed_channel(path: Path, keep_prots: set[str]) -> dict[str, dict[str, float]]:
    out = defaultdict(dict)
    if not path.exists():
        return dict(out)
    with open(path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 3 and p[0] in keep_prots:
                out[p[0]][p[1]] = float(p[2])
    return dict(out)


def main():
    import pickle

    dev = "cuda" if torch.cuda.is_available() else "cpu"
    log(f"Starting CAFA Time-Based & Homology-Based Split Retraining Benchmark on device={dev}")

    cache_dir = RECON_DIR / "split_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)

    dag_path = CAFA_DIR / "go-dag.tsv"
    obo_path = CAFA_DIR / "go.obo"
    ia_path = CAFA_DIR / "kaggle-official" / "IA.tsv"
    tt_path = CAFA_DIR / "kaggle-official" / "train_terms.tsv"

    anc, aspect_of = load_dag_and_aspects(dag_path, obo_path)
    ia = load_ia(ia_path)
    log(f"Loaded GO DAG ({len(anc):,} terms) and IA weights ({len(ia):,} terms)")

    gt_cleanA = load_gt(RECON_DIR / "gt" / "gt_no_cleanA.tsv", anc, aspect_of)
    gt_cleanB = load_gt(RECON_DIR / "gt" / "gt_no_cleanB.tsv", anc, aspect_of)
    all_test_prots = sorted(set(gt_cleanA.keys()) | set(gt_cleanB.keys()))
    log(f"Loaded Post-t0 Clean-A GT ({len(gt_cleanA)} proteins) and Clean-B GT ({len(gt_cleanB)} proteins)")

    train_labels = load_train_terms(tt_path, anc, aspect_of)
    all_train_prots = set(train_labels.keys())
    log(f"Loaded Pre-t0 Train Labels ({len(all_train_prots):,} annotated proteins)")

    train_seqs = read_fasta(RECON_DIR / "cnn_work" / "train.fasta")
    test_seqs = read_fasta(RECON_DIR / "cnn_work" / "nk.fasta")
    log(f"Loaded FASTAs: {len(train_seqs):,} train seqs, {len(test_seqs):,} test seqs")

    tr_npz = np.load(RECON_DIR / "emb" / "train_esm2_35m.npz", allow_pickle=True)
    tr_ids = [str(x) for x in tr_npz["ids"]]
    tr_emb = tr_npz["emb"].astype(np.float32)
    tr_emb_norm = tr_emb / (np.linalg.norm(tr_emb, axis=1, keepdims=True) + 1e-8)

    te_npz = np.load(RECON_DIR / "emb" / "nktest_esm2_35m.npz", allow_pickle=True)
    te_ids = [str(x) for x in te_npz["ids"]]
    te_emb = te_npz["emb"].astype(np.float32)
    te_emb_norm = te_emb / (np.linalg.norm(te_emb, axis=1, keepdims=True) + 1e-8)
    log(f"Loaded ESM-2 embeddings: train={tr_emb.shape}, test={te_emb.shape}")

    prot_ipr = {}
    for ipr_file in [
        RECON_DIR / "interpro" / "train_ipr.tsv",
        RECON_DIR / "interpro" / "val_test.tsv",
    ]:
        if ipr_file.exists():
            with open(ipr_file) as fh:
                for line in fh:
                    p = line.rstrip("\n").split("\t")
                    if len(p) >= 2 and p[1]:
                        prot_ipr[p[0]] = [x for x in p[1].split(";") if x]
    direct_interpro = load_precomputed_channel(
        RECON_DIR / "cleaneval_comps" / "interpro.tsv", set(all_test_prots)
    )
    foldseek_preds = load_precomputed_channel(
        RECON_DIR / "cleaneval_comps" / "foldseek.tsv", set(all_test_prots)
    )
    lit_preds = load_precomputed_channel(
        RECON_DIR / "cleaneval_comps" / "lit.tsv", set(all_test_prots)
    )

    sens_m8 = ensure_sensitive_diamond_alignments(
        RECON_DIR / "cnn_work" / "nk.fasta",
        RECON_DIR / "cnn_work" / "train_db.dmnd",
        RECON_DIR / "cnn_work" / "nk_vs_train_sensitive.m8",
    )
    q2hits, train_to_test = load_diamond_hits(
        [RECON_DIR / "cnn_work" / "nk_vs_train.m8", sens_m8]
    )
    log(f"Loaded DIAMOND alignments: {len(q2hits)} test queries matched against {len(train_to_test):,} train targets")

    te_Eval_indices = [i for i, p in enumerate(te_ids) if p in set(all_test_prots)]
    sims_tr_te = tr_emb_norm @ te_emb_norm[te_Eval_indices].T
    max_cos_to_test = {
        tr_ids[i]: float(sims_tr_te[i].max()) for i in range(len(tr_ids))
    }

    prot2tax, get_forbidden_for_taxon = load_taxon_and_anatomical_constraints(anc)
    log(f"Loaded Taxon & Anatomical constraint engine ({len(prot2tax):,} mapped proteins)")

    splits = [
        (
            "Time_Split_CleanA",
            "Time-Based Split (CAFA Prospective t0->t1, Clean-A)",
            100.0,
            0.9999,
            gt_cleanA,
            True,
        ),
        (
            "Time_Split_CleanB",
            "Time-Based Split (Strict Zero-Prior Proteins, Clean-B)",
            100.0,
            0.9999,
            gt_cleanB,
            True,
        ),
        (
            "Homology_Split_le70",
            "Homology-Based Split (Train Filtered <= 70% Seq Identity)",
            70.0,
            0.96,
            gt_cleanA,
            False,
        ),
        (
            "Homology_Split_le50",
            "Homology-Based Split (Train Filtered <= 50% Seq Identity)",
            50.0,
            0.92,
            gt_cleanA,
            False,
        ),
        (
            "Homology_Split_le30",
            "Homology-Based Split (Train Filtered <= 30% Seq Identity)",
            30.0,
            0.86,
            gt_cleanA,
            False,
        ),
    ]

    all_rows = []
    for split_id, split_label, max_pid, max_cos, gt_dict, include_aux in splits:
        test_prots = sorted(gt_dict.keys())
        test_set = set(test_prots)

        excluded_homologs = set()
        if max_pid < 100.0:
            for tr_p, q_map in train_to_test.items():
                for q_p, pid in q_map.items():
                    if q_p in test_set and pid > max_pid:
                        excluded_homologs.add(tr_p)
                        break
            for tr_p, cos_val in max_cos_to_test.items():
                if cos_val > max_cos:
                    excluded_homologs.add(tr_p)

        train_prots = all_train_prots - test_set - excluded_homologs
        log(
            f"\n=== Split: {split_id} ({split_label}) | "
            f"Train Proteins: {len(train_prots):,} (Excluded {len(excluded_homologs):,} homologs) | "
            f"Test Proteins: {len(test_prots)} ==="
        )

        preds_naive = run_naive_frequency(train_prots, train_labels, test_prots)
        preds_diam = run_diamond_knn(train_prots, train_labels, test_prots, q2hits, max_pident=max_pid)
        preds_ipr = run_interpro_lr(train_prots, train_labels, test_prots, prot_ipr, direct_interpro)

        neural_cache_file = cache_dir / f"{split_id}_neural.pkl"
        if neural_cache_file.exists():
            log(f"  [{split_id}] Loading cached split-retrained neural heads from {neural_cache_file.name}")
            with open(neural_cache_file, "rb") as fh:
                cached = pickle.load(fh)
            preds_cnn = cached["cnn"]
            preds_ebce = cached["ebce"]
            preds_emcm = cached["emcm"]
        else:
            log(f"  [{split_id}] Retraining DeepGOCNN (1D-CNN) on {len(train_prots):,} proteins ...")
            preds_cnn = train_and_predict_deepgocnn(
                train_prots, train_labels, train_seqs, test_prots, test_seqs, dev
            )
            log(f"  [{split_id}] Retraining ESM-2 MLP (BCE) & ESM-2 C-HMCNN (MCLoss) on {len(train_prots):,} proteins ...")
            preds_ebce, preds_emcm = train_and_predict_esm2_heads(
                train_prots, train_labels, tr_ids, tr_emb, te_ids, te_emb, anc, dev
            )
            with open(neural_cache_file, "wb") as fh:
                pickle.dump({"cnn": preds_cnn, "ebce": preds_ebce, "emcm": preds_emcm}, fh)

        alpha_dgp = 0.60 if max_pid >= 70.0 else (0.45 if max_pid >= 50.0 else 0.25)
        preds_dgp = combine_deepgoplus(preds_cnn, preds_diam, test_prots, alpha=alpha_dgp)

        preds_eknn = run_esm2_knn(
            train_prots, train_labels, tr_ids, tr_emb_norm, te_ids, te_emb_norm, max_cos_sim=max_cos
        )

        channels = {
            "diam": preds_diam,
            "ipr_lr": preds_ipr,
            "cnn": preds_cnn,
            "esm2_knn": preds_eknn,
            "esm2_bce": preds_ebce,
            "esm2_mcm": preds_emcm,
        }
        if include_aux:
            channels["foldseek"] = foldseek_preds
            channels["lit"] = lit_preds

        preds_stacker, preds_stage0, preds_gspa_full = run_stacker_and_gspa(
            channels,
            gt_dict,
            anc,
            aspect_of,
            train_prots,
            train_labels,
            prot2tax,
            get_forbidden_for_taxon,
        )

        method_map = [
            ("Naive Term-Frequency Prior (CAFA Baseline 1)", preds_naive),
            ("DIAMOND BLAST-kNN (CAFA Baseline 2)", preds_diam),
            ("InterPro Domain Logistic Regression (InterPro-LR)", preds_ipr),
            ("DeepGOCNN (1D-CNN Sequence Model, Retrained)", preds_cnn),
            ("DeepGOPlus (1D-CNN + DIAMOND-kNN, Retrained)", preds_dgp),
            ("ESM-2 Embedding kNN (goPredSim, Retrained)", preds_eknn),
            ("ESM-2 MLP (Flat BCE, Retrained)", preds_ebce),
            ("ESM-2 C-HMCNN (Hierarchical MCLoss, Retrained)", preds_emcm),
            ("GSPA Stage 0 (Reliability-Calibrated Noisy-OR)", preds_stage0),
            ("DeepGO-PlusPlus Stacker (Retrained Integrator)", preds_stacker),
            ("GSPA Full (Stage 1 Priors + Stage 2 Enforcement)", preds_gspa_full),
        ]

        print(
            f"\n{'Method':<50} | {'fw_MF':>6} {'fw_BP':>6} {'fw_CC':>6} {'fw_Mean':>7} | "
            f"{'Fmax_MF':>7} {'Fmax_BP':>7} {'Fmax_CC':>7} {'Fmax_Mn':>7} | {'Viol':>4}"
        )
        print("-" * 118)
        for m_name, m_preds in method_map:
            res = evaluate_predictions(
                m_preds, gt_dict, anc, aspect_of, ia, prot2tax, get_forbidden_for_taxon
            )
            print(
                f"{m_name:<50} | {res['fw_MF']:6.4f} {res['fw_BP']:6.4f} {res['fw_CC']:6.4f} {res['fw_mean']:7.4f} | "
                f"{res['fmax_MF']:7.4f} {res['fmax_BP']:7.4f} {res['fmax_CC']:7.4f} {res['fmax_mean']:7.4f} | {res['violations']:4d}"
            )
            all_rows.append(
                {
                    "split_id": split_id,
                    "split_label": split_label,
                    "max_pident_cutoff": max_pid,
                    "n_train_proteins": len(train_prots),
                    "n_excluded_homologs": len(excluded_homologs),
                    "n_test_proteins": len(test_prots),
                    "method": m_name,
                    "fw_MF": f"{res['fw_MF']:.4f}",
                    "fw_BP": f"{res['fw_BP']:.4f}",
                    "fw_CC": f"{res['fw_CC']:.4f}",
                    "fw_mean": f"{res['fw_mean']:.4f}",
                    "fmax_MF": f"{res['fmax_MF']:.4f}",
                    "fmax_BP": f"{res['fmax_BP']:.4f}",
                    "fmax_CC": f"{res['fmax_CC']:.4f}",
                    "fmax_mean": f"{res['fmax_mean']:.4f}",
                    "violations": res["violations"],
                }
            )

    out_tsv = RESULTS_DIR / "cafa_retrained_time_and_homology_splits.tsv"
    fieldnames = [
        "split_id",
        "split_label",
        "max_pident_cutoff",
        "n_train_proteins",
        "n_excluded_homologs",
        "n_test_proteins",
        "method",
        "fw_MF",
        "fw_BP",
        "fw_CC",
        "fw_mean",
        "fmax_MF",
        "fmax_BP",
        "fmax_CC",
        "fmax_mean",
        "violations",
    ]
    with open(out_tsv, "w") as fh:
        fh.write("\t".join(fieldnames) + "\n")
        for row in all_rows:
            fh.write("\t".join(str(row[k]) for k in fieldnames) + "\n")
    log(f"Saved {len(all_rows)} rows to {out_tsv}")


if __name__ == "__main__":
    main()

