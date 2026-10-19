#!/usr/bin/env python3
"""
Whole-Genome Withheld CAFA Benchmark across Time-Based and Homology-Based Splits.

Withholds 6 complete prokaryotic reference genomes (20,697 total coding sequences,
3,315 experimentally annotated target proteins, 1,060 strictly post-t0 newly annotated
target proteins). Input files hold 4,199 operon records and 1,938 gapseq gap records (per-genome
counts below); the loader keeps the 3,809 operons with >= 2 members in the proteome and the 1,467
gaps that map to at least one GO term, and prints the per-genome retained counts to the log:
  1. Bacillus subtilis 168            (bsubtilis,     NCBITaxon:1423, 4,288 ORFs,   425 exp, 748 operons, 398 gaps)
  2. Pseudomonas aeruginosa PAO1      (paeruginosa,   NCBITaxon:287,  5,563 ORFs,   916 exp, 1,128 operons, 509 gaps)
  3. Mycobacterium tuberculosis H37Rv (mtb,           NCBITaxon:1773, 3,997 ORFs, 1,575 exp, 863 operons, 354 gaps)
  4. Helicobacter pylori 26695        (hpylori,       NCBITaxon:210,  1,554 ORFs,   114 exp, 272 operons, 195 gaps)
  5. Synechocystis sp. PCC 6803       (synechocystis, NCBITaxon:1148, 3,508 ORFs,   208 exp, 819 operons, 351 gaps)
  6. Methanocaldococcus jannaschii    (mjannaschii,   NCBITaxon:2190, 1,787 ORFs,    77 exp, 369 operons, 131 gaps)

Escherichia coli K-12 MG1655 (ecoli, NCBITaxon:562, 4,403 ORFs) and a held-out validation
pool of 600 pre-t0 Swiss-Prot proteins are kept strictly in the training/calibration pool
so zero hyperparameters are ever tuned on the 6 withheld test genomes.

Evaluates 6 Whole-Genome Splits:
  1. WG_Time_Rollover       : Prospective Time-Based Split (t0 -> t1, pre-t0 train & chromosomal rollover -> 1,060 post-t0 targets)
  2. WG_Time_Withheld       : Whole-Genome-Withheld Time Split (t0 train with all 6 genomes purged -> 1,060 post-t0 targets)
  3. WG_Homology_le100      : Leave-Whole-Genome-Out Split (all 6 genomes & species taxids purged -> 3,315 exp targets)
  4. WG_Homology_le70       : Whole-Genome Homology Split <= 70% Sequence Identity (train homologs > 70% purged -> 3,315 exp targets)
  5. WG_Homology_le50       : Whole-Genome Homology Split <= 50% Sequence Identity (train homologs > 50% purged -> 3,315 exp targets)
  6. WG_Homology_le30       : Whole-Genome Homology Split <= 30% Sequence Identity (train homologs > 30% purged -> 3,315 exp targets)

Evaluates 15 Methods / Whole-Genome Ablations per Split (all baselines retrained from scratch on D_train^(s)
and executed across all 20,697 coding sequences of the 6 withheld genomes):
  1. Naive Term-Frequency Prior (CAFA Baseline 1)
  2. DIAMOND BLAST-kNN (CAFA Baseline 2)
  3. InterPro Domain Logistic Regression (InterPro-LR)
  4. DeepGOCNN (1D-CNN Sequence Model, Retrained)
  5. DeepGOPlus (1D-CNN + DIAMOND-kNN, Retrained)
  6. ESM-2 Embedding kNN (goPredSim, Retrained)
  7. ESM-2 MLP (Flat BCE, Retrained)
  8. ESM-2 C-HMCNN (Hierarchical MCLoss, Retrained)
  9. GSPA-Seq (Ours, Per-Protein Integrator)
 10. GSPA Stage 0 (Per-Protein Reliability-Calibrated Noisy-OR)
 11. GSPA + Operon Genomic Context Only (Stage 0/Base + B6 Operons)
 12. GSPA + Pathway & Complex Coherence Only (Stage 0/Base + B1/B2 Essentiality & Coherence)
 13. GSPA + Metabolic Gap-Filling & FBA Only (Stage 0/Base + B5 Gaps & Stage 2B FBA)
 14. GSPA + Taxon & Anatomical Constraints Only (Stage 0/Base + B3/B4 + Stage 2A Relaxation)
 15. GSPA Full Whole-Genome (Stages 0 + 1 All Priors + 2A Constraints + 2B FBA)
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import pickle
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn

ROOT = Path(os.environ.get("GSPA_ROOT", str(Path(__file__).resolve().parent.parent)))
CAFA_DIR = Path(os.environ.get("GSPA_CAFA_DIR", str(Path.home() / "Public" / "software" / "cafa6")))
RECON_DIR = ROOT / "benchmark" / "neural" / "cafa6_recon"
WG_DIR = ROOT / "benchmark" / "whole_genome_cafa"
WG_DATA = WG_DIR / "data"
CACHE_DIR = WG_DIR / "cache"
RESULTS_DIR = ROOT / "results"

ASP_ROOT = {
    "MF": "GO:0003674",
    "BP": "GO:0008150",
    "CC": "GO:0005575",
}
ROOT_TERMS = set(ASP_ROOT.values())
AA = "ACDEFGHIKLMNPQRSTVWY"
AA_IDX = {c: i + 1 for i, c in enumerate(AA)}

TEST_GENOMES = [
    {"tag": "bsubtilis", "taxon": "1423", "strain_tax": "224308", "domain": "bacteria"},
    {"tag": "paeruginosa", "taxon": "287", "strain_tax": "208964", "domain": "bacteria"},
    {"tag": "mtb", "taxon": "1773", "strain_tax": "83332", "domain": "bacteria"},
    {"tag": "hpylori", "taxon": "210", "strain_tax": "85962", "domain": "bacteria"},
    {"tag": "synechocystis", "taxon": "1148", "strain_tax": "1111708", "domain": "bacteria"},
    {"tag": "mjannaschii", "taxon": "2190", "strain_tax": "243232", "domain": "archaea"},
]
TEST_TAXIDS = {
    g["taxon"] for g in TEST_GENOMES
} | {
    g["strain_tax"] for g in TEST_GENOMES
}


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# ==============================================================================
# 1. Ontology, IA Weights, Cross-Aspect Part-Of Links, & Reference Databases
# ==============================================================================
def load_dag_aspects_and_cross_links(dag_path: Path, obo_path: Path):
    anc = defaultdict(set)
    desc = defaultdict(set)
    with open(dag_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c, _, a = line.rstrip("\n").partition("\t")
            if c and a:
                anc[c].add(a)
                desc[a].add(c)
    for c in list(anc):
        anc[c].add(c)
        desc[c].add(c)

    aspect_of = {}
    ns_map = {
        "molecular_function": "MF",
        "biological_process": "BP",
        "cellular_component": "CC",
    }
    mf_to_bp_partof = defaultdict(set)
    cur_id = None
    with open(obo_path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("id: GO:"):
                cur_id = line.split("id: ")[1].strip()
                anc[cur_id].add(cur_id)
                desc[cur_id].add(cur_id)
            elif line.startswith("namespace: ") and cur_id:
                ns = line.split("namespace: ")[1].strip()
                if ns in ns_map:
                    aspect_of[cur_id] = ns_map[ns]
            elif line.startswith("relationship: part_of GO:") and cur_id:
                tgt = "GO:" + line.split("GO:")[1].split()[0].strip()
                mf_to_bp_partof[cur_id].add(tgt)

    for t, A in anc.items():
        if t not in aspect_of:
            for asp, r in ASP_ROOT.items():
                if r in A:
                    aspect_of[t] = asp
                    break

    # Filter cross-aspect links to MF -> BP (and propagate to BP ancestors)
    mf_bp_links = defaultdict(set)
    for src, tgts in mf_to_bp_partof.items():
        if aspect_of.get(src) == "MF":
            for tgt in tgts:
                if aspect_of.get(tgt) == "BP":
                    for a in anc.get(tgt, (tgt,)):
                        if a not in ROOT_TERMS and aspect_of.get(a) == "BP":
                            mf_bp_links[src].add(a)

    return anc, desc, aspect_of, mf_bp_links


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


def load_reference_pathways_complexes_essentials(ref_dir: Path, anc: dict[str, set[str]], desc: dict[str, set[str]]):
    """
    Loads:
      - ec2go: EC -> set of GO terms
      - pathways: pathway_id -> set of required GO terms (from kegg_pathways.tsv & kegg_modules.tsv)
      - heteromeric_complexes: set of GO terms classified as 'n' (obligate heteromeric) in complex_terms.tsv
      - essential_go: domain -> list of (essential_go_id, set of descendant GO terms)
    """
    ec2go = defaultdict(set)
    ec2go_path = ref_dir / "ec2go.txt"
    if ec2go_path.exists():
        with open(ec2go_path) as fh:
            for line in fh:
                if line.startswith("!") or " ; GO:" not in line:
                    continue
                left, _, right = line.strip().partition(">")
                ec = left.strip().replace("EC:", "").strip()
                go_id = "GO:" + right.split("; GO:")[1].split()[0].strip()
                ec2go[ec].add(go_id)
                ec2go[f"EC:{ec}"].add(go_id)

    pathways = defaultdict(set)
    for pw_file in ("kegg_pathways.tsv", "kegg_modules.tsv"):
        p_path = ref_dir / pw_file
        if p_path.exists():
            with open(p_path) as fh:
                fh.readline()
                for line in fh:
                    p = line.rstrip("\n").split("\t")
                    if len(p) >= 3 and p[0] and p[2].startswith("GO:"):
                        pathways[p[0]].add(p[2].strip())

    heteromeric_complexes = set()
    homomeric_complexes = set()
    cplx_path = ROOT / "gspa-core" / "src" / "main" / "resources" / "complex-terms" / "complex_terms.tsv"
    if cplx_path.exists():
        with open(cplx_path) as fh:
            for line in fh:
                if line.startswith("#") or line.startswith("GO_term"):
                    continue
                p = line.rstrip("\n").split("\t")
                if len(p) >= 2:
                    if p[1].strip() == "n":
                        heteromeric_complexes.add(p[0].strip())
                    elif p[1].strip() == "h":
                        homomeric_complexes.add(p[0].strip())
    for c_desc in desc.get("GO:0032991", ()):
        if c_desc != "GO:0032991" and c_desc not in homomeric_complexes:
            heteromeric_complexes.add(c_desc)

    essential_go = {"bacteria": [], "archaea": []}
    ess_dir = ROOT / "gspa-core" / "src" / "main" / "resources" / "essential-functions"
    for dom in ("bacteria", "archaea"):
        ep = ess_dir / f"{dom}.tsv"
        if ep.exists():
            with open(ep) as fh:
                for line in fh:
                    line = line.strip()
                    if not line or line.startswith("#") or line.startswith("go_id") or line.startswith("GO_ID"):
                        continue
                    gid = line.split("\t")[0].strip()
                    if gid.startswith("GO:"):
                        essential_go[dom].append((gid, set(desc.get(gid, {gid}))))

    return ec2go, pathways, heteromeric_complexes, essential_go


def load_whole_genome_structures(
    wg_data: Path,
    anc: dict[str, set[str]],
    desc: dict[str, set[str]],
    aspect_of: dict[str, str],
    ec2go: dict[str, set[str]],
    pre_t0_labels: dict[str, dict[str, set[str]]],
):
    """
    Loads the 6 withheld test genomes' complete proteomes, operons, gapseq metabolic gaps,
    InterProScan domain & GO features, and experimental GO ground truth (both full exp
    and strictly post-t0 newly annotated target proteins).
    """
    genome_prots = {}
    prot_to_genome = {}
    prot_to_taxon = {}
    genome_operons = {}
    genome_gaps = {}
    gt_all_exp = defaultdict(lambda: defaultdict(set))
    gt_post_t0 = defaultdict(lambda: defaultdict(set))
    wg_ipr = defaultdict(list)
    wg_direct_ipr2go = defaultdict(dict)

    generic_roots = {
        "GO:0003824", "GO:0016491", "GO:0016740", "GO:0016787",
        "GO:0016829", "GO:0016853", "GO:0016874",
    }

    for ginfo in TEST_GENOMES:
        tag = ginfo["tag"]
        tax = ginfo["taxon"]
        prots = []
        with open(wg_data / f"{tag}.faa") as fh:
            for line in fh:
                if line.startswith(">"):
                    pid = line[1:].split()[0].strip()
                    if "|" in pid:
                        pid = pid.split("|")[1]
                    prots.append(pid)
                    prot_to_genome[pid] = tag
                    prot_to_taxon[pid] = tax
        genome_prots[tag] = sorted(set(prots))
        pset = set(genome_prots[tag])

        # Operons
        ops = []
        op_path = wg_data / "operons" / f"{tag}_operons.tsv"
        if op_path.exists():
            with open(op_path) as fh:
                for line in fh:
                    members = [x.strip() for x in line.rstrip("\n").split("\t") if x.strip() in pset]
                    if len(members) >= 2:
                        ops.append(members)
        genome_operons[tag] = ops

        # Metabolic gaps (gapseq)
        gaps = []
        gap_path = wg_data / "bench9" / "gapseq" / f"{tag}_gaps.jsonl"
        if gap_path.exists():
            with open(gap_path) as fh:
                for line in fh:
                    if not line.strip():
                        continue
                    rec = json.loads(line)
                    go_terms = set()
                    if rec.get("go_term"):
                        for g_raw in str(rec["go_term"]).replace(";", ",").split(","):
                            g_id = g_raw.strip()
                            if g_id.startswith("GO:"):
                                go_terms.add(g_id)
                    if rec.get("ec_number"):
                        for ec_raw in str(rec["ec_number"]).replace("/", ",").replace(";", ",").split(","):
                            ec_id = ec_raw.strip().replace("EC:", "").strip()
                            if ec_id in ec2go:
                                go_terms.update(ec2go[ec_id])
                            else:
                                ec_parts = ec_id.split(".")
                                if len(ec_parts) >= 3:
                                    ec3 = f"{ec_parts[0]}.{ec_parts[1]}.{ec_parts[2]}.-"
                                    if ec3 in ec2go:
                                        go_terms.update(ec2go[ec3])
                    expanded_go = set(go_terms)
                    for g_id in list(go_terms):
                        if g_id not in generic_roots:
                            expanded_go.update(desc.get(g_id, ()))
                    if expanded_go:
                        gaps.append({
                            "pathway_id": rec.get("pathway_id", ""),
                            "reaction_id": rec.get("reaction_id", ""),
                            "go_terms": expanded_go,
                            "core_go_terms": go_terms,
                            "guessed": bool(rec.get("gapseq_guessed", False)),
                        })
        genome_gaps[tag] = gaps

        # InterProScan TSV
        ipr_path = wg_data / "interproscan" / f"{tag}.tsv"
        if ipr_path.exists():
            ipr_sets = defaultdict(set)
            with open(ipr_path) as fh:
                for line in fh:
                    p = line.rstrip("\n").split("\t")
                    if len(p) < 12:
                        continue
                    pid = p[0].strip()
                    if "|" in pid:
                        pid = pid.split("|")[1]
                    if pid not in pset:
                        continue
                    if p[11].startswith("IPR"):
                        ipr_sets[pid].add(p[11].strip())
                    if len(p) >= 14 and "GO:" in p[13]:
                        for tok in p[13].split("|"):
                            if tok.startswith("GO:"):
                                gid = tok.split("(")[0].strip()
                                wg_direct_ipr2go[pid][gid] = 0.85
            for pid, s in ipr_sets.items():
                wg_ipr[pid] = sorted(s)

        # Experimental GO ground truth
        exp_path = wg_data / "truth_dual" / f"{tag}_truth_exp.tsv"
        raw_exp_by_p = defaultdict(set)
        with open(exp_path) as fh:
            fh.readline()
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) >= 3 and p[0] in pset and p[2].startswith("GO:"):
                    raw_exp_by_p[p[0]].add(p[2].strip())

        for pid, raw_terms in raw_exp_by_p.items():
            for term in raw_terms:
                for a in anc.get(term, (term,)):
                    if a in ROOT_TERMS:
                        continue
                    asp = aspect_of.get(a)
                    if asp:
                        gt_all_exp[pid][asp].add(a)
            # Strictly post-t0 zero-prior newly annotated proteins
            if pid not in pre_t0_labels:
                for asp, ts in gt_all_exp[pid].items():
                    if ts:
                        gt_post_t0[pid][asp] = set(ts)

    return (
        genome_prots,
        prot_to_genome,
        prot_to_taxon,
        genome_operons,
        genome_gaps,
        gt_all_exp,
        gt_post_t0,
        wg_ipr,
        wg_direct_ipr2go,
    )


# ==============================================================================
# 2. Fast Vectorized CAFA Evaluation Engine (Fmax, IA-Weighted f_w, Violations)
# ==============================================================================
def evaluate_predictions(
    preds: dict[str, dict[str, float]],
    gt: dict[str, dict[str, set[str]]],
    anc: dict[str, set[str]],
    aspect_of: dict[str, str],
    ia: dict[str, float],
    prot2tax: dict[str, str],
    get_forbidden_for_taxon,
    already_propagated: bool = False,
) -> dict[str, float]:
    test_prots = sorted(gt.keys())
    prop_preds = defaultdict(lambda: defaultdict(dict))
    for prot in test_prots:
        raw_p = preds.get(prot, {})
        if not raw_p:
            continue
        p_asp = prop_preds[prot]
        if already_propagated:
            for term, sc in raw_p.items():
                if sc >= 0.01 and term not in ROOT_TERMS:
                    asp = aspect_of.get(term)
                    if asp:
                        p_asp[asp][term] = sc
        else:
            for term, sc in raw_p.items():
                if sc < 0.01:
                    continue
                for a in anc.get(term, (term,)):
                    if a in ROOT_TERMS:
                        continue
                    asp = aspect_of.get(a)
                    if asp:
                        if sc > p_asp[asp].get(a, 0.0):
                            p_asp[asp][a] = sc

    thresholds = np.linspace(0.02, 0.98, 49)
    n_tau = len(thresholds)
    metrics = {}
    total_violations = 0
    total_violations_050 = 0

    for asp in ("MF", "BP", "CC"):
        asp_prots = [p for p in test_prots if len(gt[p].get(asp, ())) > 0]
        n_total = len(asp_prots)
        if n_total == 0:
            metrics[f"fmax_{asp}"] = 0.0
            metrics[f"fw_{asp}"] = 0.0
            metrics[f"smin_{asp}"] = 0.0
            continue

        sum_prec = np.zeros(n_tau, dtype=np.float64)
        sum_rec = np.zeros(n_tau, dtype=np.float64)
        n_pred = np.zeros(n_tau, dtype=np.int64)
        sum_w_prec = np.zeros(n_tau, dtype=np.float64)
        sum_w_rec = np.zeros(n_tau, dtype=np.float64)
        n_w_pred = np.zeros(n_tau, dtype=np.int64)
        sum_ru = np.zeros(n_tau, dtype=np.float64)
        sum_mi = np.zeros(n_tau, dtype=np.float64)

        for p in asp_prots:
            tset = gt[p][asp]
            t_len = len(tset)
            t_ia = sum(ia.get(t, 0.0) for t in tset)
            pdict = prop_preds[p].get(asp, {})
            if not pdict:
                sum_ru += t_ia
                continue

            terms = list(pdict.keys())
            scores = np.array([pdict[t] for t in terms], dtype=np.float32)
            in_truth = np.array([1.0 if t in tset else 0.0 for t in terms], dtype=np.float64)
            term_ias = np.array([ia.get(t, 0.0) for t in terms], dtype=np.float64)
            true_ias = in_truth * term_ias

            mask = scores[:, None] >= thresholds[None, :]
            p_cnt = mask.sum(axis=0)
            tp_cnt = (in_truth[:, None] * mask).sum(axis=0)

            has_p = p_cnt > 0
            n_pred[has_p] += 1
            sum_prec[has_p] += tp_cnt[has_p] / p_cnt[has_p]
            sum_rec[has_p] += tp_cnt[has_p] / t_len

            p_ia_sum = (term_ias[:, None] * mask).sum(axis=0)
            tp_ia_sum = (true_ias[:, None] * mask).sum(axis=0)
            sum_ru += np.maximum(0.0, t_ia - tp_ia_sum)
            sum_mi += np.maximum(0.0, p_ia_sum - tp_ia_sum)
            has_wp = p_ia_sum > 0
            n_w_pred[has_wp] += 1
            sum_w_prec[has_wp] += tp_ia_sum[has_wp] / p_ia_sum[has_wp]
            if t_ia > 0:
                sum_w_rec[has_wp] += tp_ia_sum[has_wp] / t_ia

        best_fmax = 0.0
        best_fw = 0.0
        best_tau = 0.30
        ru_avg = sum_ru / n_total
        mi_avg = sum_mi / n_total
        best_smin = float(np.min(np.sqrt(ru_avg * ru_avg + mi_avg * mi_avg)))
        for ti in range(n_tau):
            if n_pred[ti] > 0:
                pr = sum_prec[ti] / n_pred[ti]
                rc = sum_rec[ti] / n_total
                if pr + rc > 0:
                    f = 2.0 * pr * rc / (pr + rc)
                    if f > best_fmax:
                        best_fmax = float(f)
            if n_w_pred[ti] > 0:
                w_pr = sum_w_prec[ti] / n_w_pred[ti]
                w_rc = sum_w_rec[ti] / n_total
                if w_pr + w_rc > 0:
                    fw = 2.0 * w_pr * w_rc / (w_pr + w_rc)
                    if fw > best_fw:
                        best_fw = float(fw)
                        best_tau = float(thresholds[ti])

        metrics[f"fmax_{asp}"] = best_fmax
        metrics[f"fw_{asp}"] = best_fw
        metrics[f"smin_{asp}"] = best_smin
        metrics[f"tau_{asp}"] = best_tau

        # `violations`: forbidden (taxon-constraint) terms scored at or above the method's
        # f_w-optimal threshold for this aspect, capped at 0.25 (legacy definition).
        # `violations_050`: forbidden terms scored at or above a fixed threshold of 0.50.
        tau_check = min(best_tau, 0.25)
        for p in asp_prots:
            forb = get_forbidden_for_taxon(prot2tax.get(p, ""))
            if not forb:
                continue
            for t, sc in prop_preds[p].get(asp, {}).items():
                if t in forb:
                    if sc >= tau_check:
                        total_violations += 1
                    if sc >= 0.50:
                        total_violations_050 += 1

    metrics["fmax_mean"] = (metrics["fmax_MF"] + metrics["fmax_BP"] + metrics["fmax_CC"]) / 3.0
    metrics["fw_mean"] = (metrics["fw_MF"] + metrics["fw_BP"] + metrics["fw_CC"]) / 3.0
    metrics["smin_mean"] = (metrics["smin_MF"] + metrics["smin_BP"] + metrics["smin_CC"]) / 3.0
    metrics["violations"] = total_violations
    metrics["violations_050"] = total_violations_050
    return metrics


# ==============================================================================
# 3. Import Baseline Models from run_cafa_time_and_homology_splits
# ==============================================================================
from run_cafa_time_and_homology_splits import (
    AspectMLP,
    build_dag_edges,
    combine_deepgoplus,
    load_diamond_hits,
    load_taxon_and_anatomical_constraints,
    load_train_terms,
    mc_loss,
    mcm_layer,
    read_fasta,
    run_diamond_knn,
    run_interpro_lr,
    run_naive_frequency,
    train_and_predict_deepgocnn,
)


def train_and_predict_esm2_heads_batched(
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
    preds_bce = defaultdict(dict)
    preds_mcm = defaultdict(dict)

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
                for s in range(0, len(te_ids), batch_size):
                    xb = torch.tensor(te_emb[s : s + batch_size], dtype=torch.float32, device=dev)
                    probs_b = torch.sigmoid(model(xb))
                    if loss_type == "mcm":
                        probs_b = mcm_layer(probs_b, par, chi, n_terms)
                    probs_np = probs_b.cpu().numpy()
                    for i, p in enumerate(te_ids[s : s + batch_size]):
                        row = probs_np[i]
                        nz = np.nonzero(row >= 0.025)[0]
                        for j in nz:
                            target_dict[p][terms[j]] = float(row[j])

    return dict(preds_bce), dict(preds_mcm)


def run_esm2_knn_fast(
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
    ref_id_to_j = {rid: j for j, rid in enumerate(ref_ids)}
    R = tr_emb_norm[keep_idx]
    sims = te_emb_norm @ R.T  # (n_test, n_ref)

    preds = {}
    for i, q in enumerate(te_ids):
        row = sims[i].copy()
        j_self = ref_id_to_j.get(q)
        if j_self is not None:
            row[j_self] = -2.0
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
# 4. Whole-Genome GSPA Inference Engine & Component Ablations
# ==============================================================================
def run_whole_genome_gspa_ablations(
    channels: dict[str, dict[str, dict[str, float]]],
    val_prots: list[str],
    all_wg_prots: list[str],
    eval_target_prots: set[str],
    genome_prots: dict[str, list[str]],
    genome_operons: dict[str, list[list[str]]],
    genome_gaps: dict[str, list[dict]],
    prot_to_taxon: dict[str, str],
    anc: dict[str, set[str]],
    desc: dict[str, set[str]],
    aspect_of: dict[str, str],
    ia: dict[str, float],
    mf_bp_links: dict[str, set[str]],
    pathways: dict[str, set[str]],
    heteromeric_complexes: set[str],
    essential_go: dict[str, list[tuple[str, set[str]]]],
    sub_train_prots: set[str],
    train_labels: dict[str, dict[str, set[str]]],
    train_prot2tax: dict[str, str],
    get_forbidden_for_taxon,
    max_pid: float = 100.0,
    rollover_known_gt: dict[str, dict[str, set[str]]] | None = None,
    hom_scale_mode: str = "val-tuned",
    fam_mode: str = "family",
    log_prefix: str = "",
):
    """
    Executes per-protein integration (GSPA Stage 0 Noisy-OR & GSPA-Seq Stacker)
    across all 20,697 coding sequences of the 6 withheld genomes, followed by each
    whole-genome Stage 1 + Stage 2 ablation:
      - Stage 0 Noisy-OR
      - GSPA-Seq Stacker
      - GSPA + Operon Genomic Context Only (B6)
      - GSPA + Pathway & Complex Coherence Only (B1/B2)
      - GSPA + Metabolic Gap-Filling & FBA Only (B5 + Stage 2B SLR)
      - GSPA + Taxon & Anatomical Constraints Only (B3/B4 + Stage 2A)
      - GSPA Full Whole-Genome (Stages 0 + 1 All Priors + 2A + 2B)
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler

    comp_names = list(channels.keys())
    full_detail_prots = sorted(set(val_prots) | set(eval_target_prots))
    full_detail_set = set(full_detail_prots)

    prop_raw = {c: defaultdict(dict) for c in comp_names}
    prop_cln = {c: defaultdict(dict) for c in comp_names}
    forb_by_prot = {}
    for p in val_prots:
        forb_by_prot[p] = get_forbidden_for_taxon(train_prot2tax.get(p, ""))
    for p in all_wg_prots:
        forb_by_prot[p] = get_forbidden_for_taxon(prot_to_taxon.get(p, ""))

    for c in comp_names:
        c_dict = channels[c]
        for p in full_detail_prots:
            raw_d = c_dict.get(p, {})
            if not raw_d:
                continue
            forb = forb_by_prot[p]
            dest_r = prop_raw[c][p]
            dest_c = prop_cln[c][p]
            for term, sc in raw_d.items():
                if sc < 0.025:
                    continue
                is_forb = term in forb
                for a in anc.get(term, (term,)):
                    if a in ROOT_TERMS or not aspect_of.get(a):
                        continue
                    if sc > dest_r.get(a, 0.0):
                        dest_r[a] = sc
                    if not is_forb and a not in forb:
                        if sc > dest_c.get(a, 0.0):
                            dest_c[a] = sc
                    elif is_forb and a not in forb:
                        rel_sc = 0.25 * sc
                        if rel_sc >= 0.025 and rel_sc > dest_c.get(a, 0.0):
                            dest_c[a] = rel_sc

        for p in all_wg_prots:
            if p in full_detail_set:
                continue
            raw_d = c_dict.get(p, {})
            if not raw_d:
                continue
            forb = forb_by_prot[p]
            dest_c = prop_cln[c][p]
            for term, sc in raw_d.items():
                if sc < 0.12 or term in forb:
                    continue
                for a in anc.get(term, (term,)):
                    if a not in ROOT_TERMS and a not in forb and aspect_of.get(a) and sc > dest_c.get(a, 0.0):
                        dest_c[a] = sc

    # 2. Calibrate Stage 0 channel reliabilities on val_prots with Evidence-Group Correlation Damping
    #    `base_rel[c][asp]` is the reliability before the MF/BP down-weighting `hom_scale` of the
    #    sequence-model channels (cnn, esm2_bce, esm2_mcm) and before clipping.
    #    hom_scale used to be a constant chosen from the split identity (max_pid). It is now
    #    either fixed at 1.0 or selected on the validation proteins of this split (see
    #    --hom-scale); the split identity is never used ("legacy-split-pid" only reproduces the
    #    pre-fix behaviour for auditing).
    base_rel = {c: {} for c in comp_names}
    for asp in ("MF", "BP", "CC"):
        val_asp_prots = [p for p in val_prots if len(train_labels.get(p, {}).get(asp, ())) > 0]
        for c in comp_names:
            best_f = 0.10
            for tau in (0.15, 0.25, 0.35, 0.45):
                sp_sum, sr_sum, np_cnt = 0.0, 0.0, 0
                for p in val_asp_prots:
                    pset = {t for t, sc in prop_raw[c][p].items() if aspect_of.get(t) == asp and sc >= tau}
                    if not pset:
                        continue
                    tset = train_labels[p][asp]
                    inter = len(pset & tset)
                    np_cnt += 1
                    sp_sum += inter / len(pset)
                    sr_sum += inter / len(tset)
                if np_cnt > 0 and len(val_asp_prots) > 0:
                    pr = sp_sum / np_cnt
                    rc = sr_sum / len(val_asp_prots)
                    if pr + rc > 0:
                        f = 2 * pr * rc / (pr + rc)
                        if f > best_f:
                            best_f = f
            w_c = best_f ** 1.35
            if c in ("esm2_knn", "esm2_bce", "esm2_mcm"):
                w_c *= 0.48
            base_rel[c][asp] = w_c

    def make_rel_weights(h_scale):
        rw = {c: {} for c in comp_names}
        for c in comp_names:
            for asp in ("MF", "BP", "CC"):
                w_c = base_rel[c][asp]
                if c in ("cnn", "esm2_bce", "esm2_mcm") and asp in ("MF", "BP"):
                    w_c *= h_scale
                rw[c][asp] = float(np.clip(w_c, 0.05, 0.88))
        return rw

    def compute_stage0(prop_dict, prot_list, min_sc=0.025, rel_weights=None):
        s0 = defaultdict(dict)
        for p in prot_list:
            cand = set()
            for c in comp_names:
                cand.update(prop_dict[c][p].keys())
            dest = s0[p]
            for t in cand:
                asp = aspect_of.get(t, "MF")
                not_p = 1.0
                for c in comp_names:
                    sc = prop_dict[c][p].get(t, 0.0)
                    if sc > 0:
                        not_p *= (1.0 - rel_weights[c].get(asp, 0.4) * sc)
                post = 1.0 - not_p
                if post >= min_sc:
                    dest[t] = post
        return s0

    hom_grid = (1.0, 0.75, 0.55, 0.38)
    hom_val_fw = {}
    if hom_scale_mode == "fixed-1.0":
        hom_scale = 1.0
    elif hom_scale_mode == "legacy-split-pid":
        hom_scale = 1.0 if max_pid >= 99.0 else (0.75 if max_pid >= 65.0 else (0.55 if max_pid >= 45.0 else 0.38))
    elif hom_scale_mode == "val-tuned":
        # Select hom_scale on this split's validation proteins only: Stage 0 noisy-OR on
        # val_prots (prop_raw), scored by mean IA-weighted f_w against their training labels.
        # Ties (to 1e-6) are broken towards the larger value (less down-weighting).
        val_gt_h = {p: train_labels[p] for p in val_prots if p in train_labels}
        for h_cand in hom_grid:
            s0_val = compute_stage0(prop_raw, val_prots, min_sc=0.025, rel_weights=make_rel_weights(h_cand))
            res_h = evaluate_predictions(
                s0_val, val_gt_h, anc, aspect_of, ia, train_prot2tax, get_forbidden_for_taxon, already_propagated=True
            )
            hom_val_fw[h_cand] = res_h["fw_mean"]
        hom_scale = max(hom_grid, key=lambda h: (round(hom_val_fw[h], 6), h))
    else:
        raise ValueError(f"unknown hom_scale_mode {hom_scale_mode!r}")
    rel_weights = make_rel_weights(hom_scale)
    log(
        f"  {log_prefix}hom_scale mode={hom_scale_mode} selected={hom_scale}"
        + (
            " | validation mean f_w by candidate: "
            + ", ".join(f"{h}: {hom_val_fw[h]:.4f}" for h in hom_grid)
            if hom_val_fw
            else ""
        )
    )
    log(
        f"  {log_prefix}Stage 0 channel reliabilities (MF/BP/CC): "
        + "; ".join(
            f"{c}={rel_weights[c]['MF']:.3f}/{rel_weights[c]['BP']:.3f}/{rel_weights[c]['CC']:.3f}" for c in comp_names
        )
    )

    preds_stage0_raw = compute_stage0(prop_raw, full_detail_prots, min_sc=0.025, rel_weights=rel_weights)
    preds_stage0_cln = compute_stage0(prop_cln, all_wg_prots + val_prots, min_sc=0.025, rel_weights=rel_weights)

    # 3. Fit GSPA-Seq Stacker strictly on val_prots
    preds_stacker_raw = defaultdict(dict)
    preds_stacker_cln = defaultdict(dict)
    for asp in ("MF", "BP", "CC"):
        X_val, y_val = [], []
        for p in val_prots:
            truth = train_labels.get(p, {}).get(asp, set())
            if not truth:
                continue
            cand = {t for t in preds_stage0_raw[p] if aspect_of.get(t) == asp}
            for t in cand:
                feat = [
                    prop_raw[c][p].get(t, 0.0) * (hom_scale if (c in ("cnn", "esm2_bce", "esm2_mcm") and asp in ("MF", "BP")) else 1.0)
                    for c in comp_names
                ] + [preds_stage0_raw[p].get(t, 0.0)]
                X_val.append(feat)
                y_val.append(1 if t in truth else 0)

        if not y_val or min(y_val) == max(y_val):
            continue
        X_val_arr = np.asarray(X_val, dtype=np.float32)
        y_val_arr = np.asarray(y_val, dtype=np.int8)
        sc_model = StandardScaler().fit(X_val_arr)
        lr_model = LogisticRegression(max_iter=1000, class_weight="balanced", C=0.5)
        lr_model.fit(sc_model.transform(X_val_arr), y_val_arr)

        for prop_src, s0_src, out_stk in ((prop_raw, preds_stage0_raw, preds_stacker_raw), (prop_cln, preds_stage0_cln, preds_stacker_cln)):
            rows_prot, rows_term, X_te = [], [], []
            for p in full_detail_prots:
                cand = {t for t in s0_src[p] if aspect_of.get(t) == asp}
                for t in cand:
                    feat = [
                        prop_src[c][p].get(t, 0.0) * (hom_scale if (c in ("cnn", "esm2_bce", "esm2_mcm") and asp in ("MF", "BP")) else 1.0)
                        for c in comp_names
                    ] + [s0_src[p].get(t, 0.0)]
                    X_te.append(feat)
                    rows_prot.append(p)
                    rows_term.append(t)
            if X_te:
                X_te_arr = np.asarray(X_te, dtype=np.float32)
                probs_te = lr_model.predict_proba(sc_model.transform(X_te_arr))[:, 1]
                for p, t, score in zip(rows_prot, rows_term, probs_te):
                    if score >= 0.02:
                        out_stk[p][t] = float(score)

    val_gt_dict = {p: train_labels[p] for p in val_prots if p in train_labels}
    beta_by_asp = {}
    for asp in ("MF", "BP", "CC"):
        best_b, best_fw_b = 0.0, -1.0
        b_grid = (0.0, 0.15, 0.25, 0.35, 0.50) if asp == "CC" else ((0.0, 0.08, 0.15) if asp == "MF" else (0.0,))
        for b_cand in b_grid:
            blend_v = {}
            for p in val_prots:
                blend_v[p] = {}
                for t in set(preds_stacker_raw.get(p, ())) | set(preds_stage0_raw.get(p, ())):
                    if aspect_of.get(t) == asp:
                        blend_v[p][t] = b_cand * preds_stacker_raw.get(p, {}).get(t, 0.0) + (1.0 - b_cand) * preds_stage0_raw.get(p, {}).get(t, 0.0)
            res_v = evaluate_predictions(
                blend_v, val_gt_dict, anc, aspect_of, ia, train_prot2tax, get_forbidden_for_taxon, already_propagated=True
            )
            if res_v[f"fw_{asp}"] > best_fw_b:
                best_fw_b = res_v[f"fw_{asp}"]
                best_b = b_cand
        beta_by_asp[asp] = best_b

    base_raw = defaultdict(dict)
    base_cln = defaultdict(dict)
    for p in eval_target_prots:
        for src_stk, src_s0, dst in ((preds_stacker_raw, preds_stage0_raw, base_raw), (preds_stacker_cln, preds_stage0_cln, base_cln)):
            terms = set(src_stk.get(p, {}).keys()) | set(src_s0.get(p, {}).keys())
            for t in terms:
                asp = aspect_of.get(t, "MF")
                b = beta_by_asp.get(asp, 0.0)
                sc = b * src_stk.get(p, {}).get(t, 0.0) + (1.0 - b) * src_s0.get(p, {}).get(t, 0.0)
                if sc >= 0.02:
                    dst[p][t] = sc
            for t, sc in list(dst[p].items()):
                for a in anc.get(t, (t,)):
                    if a not in ROOT_TERMS and aspect_of.get(a) and sc > dst[p].get(a, 0.0):
                        if dst is base_raw or a not in forb_by_prot[p]:
                            dst[p][a] = sc

    # Chromosomal context state across all 20,697 ORFs
    ctx_wg = defaultdict(dict)
    for p in all_wg_prots:
        if p in eval_target_prots:
            ctx_wg[p] = dict(base_cln[p])
        else:
            ctx_wg[p] = dict(preds_stage0_cln[p])
        if rollover_known_gt and p in rollover_known_gt:
            forb = forb_by_prot[p]
            for asp, ts in rollover_known_gt[p].items():
                for t in ts:
                    if t not in forb:
                        ctx_wg[p][t] = max(ctx_wg[p].get(t, 0.0), 0.92)

    multi_subunit_complexes = set(heteromeric_complexes)
    for t in desc.get("GO:0032991", ()):
        if t != "GO:0032991" and ia.get(t, 0.0) >= 0.5 and aspect_of.get(t) == "CC":
            multi_subunit_complexes.add(t)

    # Precompute split-trained empirical MF -> BP associations on sub_train_prots
    all_pw_mf_terms = set()
    term_to_pathways = defaultdict(list)
    for pw_id, req in pathways.items():
        all_pw_mf_terms.update(req)
        for r in req:
            term_to_pathways[r].append(pw_id)
            for a in anc.get(r, ()):
                if ia.get(a, 0.0) >= 3.0 and aspect_of.get(a) == "MF":
                    term_to_pathways[a].append(pw_id)

    mf_counts = defaultdict(int)
    mf_bp_cooc = defaultdict(lambda: defaultdict(int))
    for p in sub_train_prots:
        if p not in train_labels:
            continue
        mf_ts = {t for t in train_labels[p].get("MF", set()) if ia.get(t, 0.0) >= 1.2 or t in all_pw_mf_terms}
        bp_ts = {t for t in train_labels[p].get("BP", set()) if ia.get(t, 0.0) >= 1.2}
        if mf_ts and bp_ts:
            for m in mf_ts:
                mf_counts[m] += 1
                for b in bp_ts:
                    mf_bp_cooc[m][b] += 1

    mf_to_bp_emp = defaultdict(dict)
    for m, n_m in mf_counts.items():
        if n_m < 4:
            continue
        for b, c_mb in mf_bp_cooc[m].items():
            if c_mb >= 4:
                pr = (c_mb + 1.0) / (n_m + 10.0)
                if pr >= 0.42:
                    mf_to_bp_emp[m][b] = pr
    for m, bps in mf_bp_links.items():
        if ia.get(m, 0.0) >= 1.2 or m in all_pw_mf_terms:
            for b in bps:
                if ia.get(b, 0.0) >= 1.2:
                    mf_to_bp_emp[m][b] = max(mf_to_bp_emp[m].get(b, 0.0), 0.80)

    pw_bp_sig = defaultdict(dict)
    for pw_id, req_mfs in pathways.items():
        bp_acc = defaultdict(list)
        for m in req_mfs:
            for b, pr in mf_to_bp_emp.get(m, {}).items():
                bp_acc[b].append(pr)
        for b, vals in bp_acc.items():
            if len(vals) >= 2:
                pw_bp_sig[pw_id][b] = float(np.mean(sorted(vals, reverse=True)[:3]))

    # ==========================================================================
    # Whole-Genome Prior 1: Operon Genomic Context (B6)
    # ==========================================================================
    # Family expansion of KEGG pathway steps in the operon prior (B6). Formerly switched on
    # only for max_pid >= 99 (split identity); now a single split-independent choice
    # (--fam-candidates). "family" matches the B1/B2 coherence prior, which always expands
    # pathway steps to informative (IA >= 3.0) MF ancestors.
    if fam_mode == "legacy-split-pid":
        use_fam_expansion = max_pid >= 99.0
    elif fam_mode == "family":
        use_fam_expansion = True
    elif fam_mode == "exact":
        use_fam_expansion = False
    else:
        raise ValueError(f"unknown fam_mode {fam_mode!r}")
    operon_boosts = defaultdict(lambda: defaultdict(float))
    for ginfo in TEST_GENOMES:
        tag = ginfo["tag"]
        for operon in genome_operons[tag]:
            m_len = len(operon)
            if m_len < 2:
                continue
            member_pws = {}
            for pid in operon:
                pws = set()
                for t, sc in ctx_wg[pid].items():
                    if sc >= 0.42 and t in term_to_pathways:
                        pws.update(term_to_pathways[t])
                member_pws[pid] = pws

            for i, pid in enumerate(operon):
                if pid not in eval_target_prots:
                    continue
                neighbors = [operon[j] for j in range(m_len) if j != i]
                nb_term_scores = defaultdict(list)
                for nb in neighbors:
                    for t, sc in ctx_wg[nb].items():
                        asp = aspect_of.get(t)
                        if asp == "CC" and sc >= 0.30 and ia.get(t, 0.0) >= 0.4:
                            nb_term_scores[t].append(sc)
                        elif asp == "BP" and sc >= 0.48 and ia.get(t, 0.0) >= 1.8:
                            nb_term_scores[t].append(sc)
                    for t, sc in ctx_wg[nb].items():
                        if sc >= 0.52 and aspect_of.get(t) == "MF" and t in mf_to_bp_emp:
                            for bp_t, cond_p in mf_to_bp_emp[t].items():
                                if sc * cond_p >= 0.42 and ia.get(bp_t, 0.0) >= 1.8:
                                    nb_term_scores[bp_t].append(sc * cond_p)

                for t, sc_list in nb_term_scores.items():
                    asp = aspect_of.get(t)
                    max_nb = max(sc_list)
                    mean_nb = sum(sc_list) / len(neighbors)
                    n_sup = len(sc_list)
                    consensus = 0.55 * max_nb + 0.45 * mean_nb
                    cur_p = base_cln[pid].get(t, 0.0)

                    if asp == "CC":
                        if 0.05 <= cur_p <= 0.70 and (t in multi_subunit_complexes or n_sup >= 2 or max_nb >= 0.50):
                            w_op = 0.28 if t in multi_subunit_complexes else (0.20 if n_sup >= 2 else 0.14)
                            operon_boosts[pid][t] = max(operon_boosts[pid][t], w_op * consensus * min(1.0, cur_p / 0.15))
                    elif asp == "BP":
                        if 0.16 <= cur_p <= 0.58 and (n_sup >= 2 or max_nb >= 0.65):
                            w_op = 0.14 if n_sup >= 2 else 0.08
                            operon_boosts[pid][t] = max(operon_boosts[pid][t], w_op * consensus)

                # Complementary KEGG pathway step completion within the operon
                nb_pw_counts = defaultdict(int)
                for nb in neighbors:
                    for pw_id in member_pws[nb]:
                        nb_pw_counts[pw_id] += 1
                for pw_id, nb_pw_hits in nb_pw_counts.items():
                    pw_frac = nb_pw_hits / len(neighbors)
                    for r in pathways[pw_id]:
                        fam_candidates = (r,) if not use_fam_expansion else ((r,) + tuple(a for a in anc.get(r, ()) if ia.get(a, 0.0) >= 3.2 and aspect_of.get(a) == "MF"))
                        for r_or_fam in fam_candidates:
                            cur_r = base_cln[pid].get(r_or_fam, 0.0)
                            if 0.14 <= cur_r <= 0.58:
                                operon_boosts[pid][r_or_fam] = max(operon_boosts[pid][r_or_fam], 0.26 * (0.5 + 0.5 * pw_frac))
                    for bp_t, w_bp in pw_bp_sig.get(pw_id, {}).items():
                        if ia.get(bp_t, 0.0) >= 1.8:
                            cur_b = base_cln[pid].get(bp_t, 0.0)
                            if 0.16 <= cur_b <= 0.58:
                                operon_boosts[pid][bp_t] = max(
                                    operon_boosts[pid][bp_t], 0.14 * w_bp * (0.5 + 0.5 * pw_frac)
                                )

    # ==========================================================================
    # Whole-Genome Prior 2: Pathway, Complex Coherence & Core Essentiality (B1/B2)
    # ==========================================================================
    coherence_boosts = defaultdict(lambda: defaultdict(float))
    coherence_penalties = defaultdict(lambda: defaultdict(float))

    for ginfo in TEST_GENOMES:
        tag = ginfo["tag"]
        dom = ginfo["domain"]
        g_prots = genome_prots[tag]

        max_by_term = defaultdict(float)
        cands_by_term = defaultdict(list)
        for p in g_prots:
            for t, sc in ctx_wg[p].items():
                if sc > max_by_term[t]:
                    max_by_term[t] = sc
                if sc >= 0.08:
                    cands_by_term[t].append((sc, p))
        for t in cands_by_term:
            cands_by_term[t].sort(reverse=True)

        for pw_id, req_mfs in pathways.items():
            if len(req_mfs) < 2:
                continue
            step_scores = []
            for r in req_mfs:
                s_r = max_by_term.get(r, 0.0)
                if s_r < 0.35:
                    for a in anc.get(r, ()):
                        if ia.get(a, 0.0) >= 2.8 and aspect_of.get(a) == "MF":
                            s_r = max(s_r, 0.85 * max_by_term.get(a, 0.0))
                step_scores.append(min(1.0, s_r / 0.42))
            cov = sum(step_scores) / len(req_mfs)
            if cov >= 0.40:
                for r in req_mfs:
                    for r_or_fam in (r,) + tuple(a for a in anc.get(r, ()) if ia.get(a, 0.0) >= 3.0 and aspect_of.get(a) == "MF"):
                        if max_by_term.get(r_or_fam, 0.0) < 0.58 and r_or_fam in cands_by_term:
                            for sc, p in cands_by_term[r_or_fam][:2]:
                                if p in eval_target_prots and 0.14 <= sc <= 0.58:
                                    coherence_boosts[p][r_or_fam] = max(coherence_boosts[p][r_or_fam], 0.28 * (0.4 + 0.6 * cov))
                    for sc, p in cands_by_term.get(r, ())[:2]:
                        if p in eval_target_prots and sc >= 0.30:
                            for bp_t, w_bp in pw_bp_sig.get(pw_id, {}).items():
                                if ia.get(bp_t, 0.0) >= 1.8:
                                    cur_b = base_cln[p].get(bp_t, 0.0)
                                    if 0.16 <= cur_b <= 0.58:
                                        coherence_boosts[p][bp_t] = max(
                                            coherence_boosts[p][bp_t], 0.16 * cov * w_bp
                                        )

        for c_term in multi_subunit_complexes:
            cands = cands_by_term.get(c_term, [])
            if not cands:
                continue
            strong_subunits = [p for sc, p in cands if sc >= 0.38]
            if len(strong_subunits) >= 1 and len(cands) >= 2:
                for sc, p in cands[:5]:
                    if p in eval_target_prots and 0.08 <= sc <= 0.65:
                        coherence_boosts[p][c_term] = max(coherence_boosts[p][c_term], 0.24 * min(1.0, sc / 0.16))
            elif len(cands) == 1 and cands[0][0] < 0.45 and c_term in heteromeric_complexes:
                if cands[0][1] in eval_target_prots:
                    coherence_penalties[cands[0][1]][c_term] = 0.82

        for ess_root, ess_desc_set in essential_go.get(dom, ()):
            if max_by_term.get(ess_root, 0.0) < 0.48:
                best_sc, best_p, best_t = 0.0, None, None
                for dt in ess_desc_set:
                    if dt in cands_by_term and cands_by_term[dt][0][0] > best_sc:
                        best_sc, best_p = cands_by_term[dt][0]
                        best_t = dt
                if best_p in eval_target_prots and best_t and 0.14 <= best_sc <= 0.52:
                    coherence_boosts[best_p][best_t] = max(coherence_boosts[best_p][best_t], 0.25)

    # ==========================================================================
    # Whole-Genome Prior 3: Metabolic Gap-Filling & Stage 2B FBA (B5 + Stage 2B SLR)
    # ==========================================================================
    fba_gap_boosts = defaultdict(lambda: defaultdict(float))

    for ginfo in TEST_GENOMES:
        tag = ginfo["tag"]
        g_prots = genome_prots[tag]
        gaps = list(genome_gaps[tag])

        prot_operon_mates = defaultdict(set)
        for op in genome_operons[tag]:
            for p in op:
                prot_operon_mates[p].update(x for x in op if x != p)

        max_by_term = defaultdict(float)
        cands_by_term = defaultdict(list)
        for p in g_prots:
            for t, sc in ctx_wg[p].items():
                if sc > max_by_term[t]:
                    max_by_term[t] = sc
                if sc >= 0.08:
                    cands_by_term[t].append((sc, p))
        for t in cands_by_term:
            cands_by_term[t].sort(reverse=True)

        for pw_id, req_mfs in pathways.items():
            if len(req_mfs) < 2:
                continue
            active_steps = sum(1 for r in req_mfs if max_by_term.get(r, 0.0) >= 0.40)
            if active_steps >= 1 and active_steps / len(req_mfs) >= 0.30:
                for r in req_mfs:
                    if 0.10 <= max_by_term.get(r, 0.0) < 0.58:
                        fam_r = {r} | {a for a in anc.get(r, ()) if ia.get(a, 0.0) >= 3.0 and aspect_of.get(a) == "MF"}
                        gaps.append({
                            "pathway_id": pw_id,
                            "reaction_id": r,
                            "go_terms": fam_r,
                            "core_go_terms": {r},
                            "pw_partner_terms": set(req_mfs) - {r},
                            "guessed": False,
                        })

        for gap in gaps:
            w_gap = 0.85 if gap["guessed"] else 1.0
            core_terms = gap["core_go_terms"]
            fam_terms = set(core_terms)
            for cg in core_terms:
                for a in anc.get(cg, ()):
                    if ia.get(a, 0.0) >= 3.0 and aspect_of.get(a) == "MF":
                        fam_terms.add(a)
            all_gap_terms = gap["go_terms"] | fam_terms
            partner_terms = gap.get("pw_partner_terms", fam_terms) | fam_terms

            cand_weights = defaultdict(float)
            cand_best_term = {}
            for gt_id in all_gap_terms:
                is_core_or_fam = gt_id in fam_terms
                for sc, p in cands_by_term.get(gt_id, ())[:3]:
                    has_op_mate = any(
                        max((ctx_wg[nb].get(k, 0.0) for k in partner_terms), default=0.0) >= 0.30
                        for nb in prot_operon_mates.get(p, ())
                    )
                    w_p = (sc ** 1.3) * (1.50 if has_op_mate else 1.0) * (1.25 if is_core_or_fam else 0.85)
                    if w_p > cand_weights[p]:
                        cand_weights[p] = w_p
                        cand_best_term[p] = gt_id

            if not cand_weights:
                continue

            ranked_loci = sorted(cand_weights.items(), key=lambda x: -x[1])[:3]
            w0 = 0.03
            denom = w0 + sum(w for _, w in ranked_loci)
            for p, w_p in ranked_loci[:2]:
                if p not in eval_target_prots:
                    continue
                resp = w_p / denom
                if resp >= 0.25:
                    gt_id = cand_best_term[p]
                    cur_g = base_cln[p].get(gt_id, 0.0)
                    if 0.14 <= cur_g <= 0.58:
                        boost_val = 0.28 * w_gap * resp
                        fba_gap_boosts[p][gt_id] = max(fba_gap_boosts[p][gt_id], boost_val)
                        for cg in fam_terms:
                            cur_c = base_cln[p].get(cg, 0.0)
                            if 0.14 <= cur_c <= 0.58:
                                fba_gap_boosts[p][cg] = max(fba_gap_boosts[p][cg], 0.90 * boost_val)
                            for bp_t, cond_p in mf_to_bp_emp.get(cg, {}).items():
                                if ia.get(bp_t, 0.0) >= 1.8:
                                    cur_b = base_cln[p].get(bp_t, 0.0)
                                    if 0.16 <= cur_b <= 0.58:
                                        fba_gap_boosts[p][bp_t] = max(
                                            fba_gap_boosts[p][bp_t], 0.16 * w_gap * resp * cond_p
                                        )

    def build_ablation_predictions(
        use_operons: bool = False,
        use_coherence: bool = False,
        use_fba_gaps: bool = False,
        use_constraints: bool = False,
    ) -> dict[str, dict[str, float]]:
        out_preds = {}
        src_base = base_cln if use_constraints else base_raw
        for p in sorted(eval_target_prots):
            forb = forb_by_prot[p]
            scores = dict(src_base[p])

            if use_operons and p in operon_boosts:
                for t, b in operon_boosts[p].items():
                    if use_constraints and t in forb:
                        continue
                    cur = scores.get(t, 0.0)
                    scores[t] = 1.0 - (1.0 - cur) * (1.0 - min(0.75, b))

            if use_coherence:
                if p in coherence_boosts:
                    for t, b in coherence_boosts[p].items():
                        if use_constraints and t in forb:
                            continue
                        cur = scores.get(t, 0.0)
                        scores[t] = 1.0 - (1.0 - cur) * (1.0 - min(0.75, b))
                if p in coherence_penalties:
                    for t, mult in coherence_penalties[p].items():
                        if t in scores:
                            scores[t] *= mult

            if use_fba_gaps and p in fba_gap_boosts:
                for t, b in fba_gap_boosts[p].items():
                    if use_constraints and t in forb:
                        continue
                    cur = scores.get(t, 0.0)
                    scores[t] = 1.0 - (1.0 - cur) * (1.0 - min(0.75, b))

            modified_terms = set()
            if use_operons and p in operon_boosts:
                modified_terms.update(operon_boosts[p].keys())
            if use_coherence and p in coherence_boosts:
                modified_terms.update(coherence_boosts[p].keys())
            if use_fba_gaps and p in fba_gap_boosts:
                modified_terms.update(fba_gap_boosts[p].keys())

            for t in modified_terms:
                sc = scores.get(t, 0.0)
                if sc < 0.02:
                    continue
                for a in anc.get(t, (t,)):
                    if a not in ROOT_TERMS and aspect_of.get(a) and sc > scores.get(a, 0.0):
                        if not use_constraints or a not in forb:
                            scores[a] = sc

            out_preds[p] = scores
        return out_preds

    preds_operons = build_ablation_predictions(use_operons=True, use_constraints=False)
    preds_coherence = build_ablation_predictions(use_coherence=True, use_constraints=False)
    preds_fba_gaps = build_ablation_predictions(use_fba_gaps=True, use_constraints=False)
    preds_constraints = build_ablation_predictions(use_constraints=True)
    preds_gspa_full = build_ablation_predictions(
        use_operons=True, use_coherence=True, use_fba_gaps=True, use_constraints=True
    )

    return (
        prop_raw,
        dict(preds_stage0_raw),
        dict(base_raw),
        preds_operons,
        preds_coherence,
        preds_fba_gaps,
        preds_constraints,
        preds_gspa_full,
        {"hom_scale_selected": hom_scale, "hom_scale_val_fw": hom_val_fw, "rel_weights": rel_weights,
         "beta_by_asp": beta_by_asp, "fam_expansion": use_fam_expansion},
    )


# ==============================================================================
# 5. Main Execution Pipeline
# ==============================================================================
def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0] if __doc__ else None)
    ap.add_argument("--gpu-only", action="store_true", help="only build missing split-retrained neural caches")
    ap.add_argument(
        "--hom-scale",
        choices=("val-tuned", "fixed-1.0", "legacy-split-pid"),
        default="val-tuned",
        help="MF/BP down-weighting of the cnn/esm2_bce/esm2_mcm Stage 0 reliabilities: selected per split on "
        "the validation proteins from {1.0, 0.75, 0.55, 0.38} (val-tuned, default), disabled (fixed-1.0), or the "
        "pre-fix constant chosen from the split's identity cutoff (legacy-split-pid; audit only)",
    )
    ap.add_argument(
        "--fam-candidates",
        choices=("family", "exact", "legacy-split-pid"),
        default="family",
        help="operon prior (B6) pathway-step expansion to informative MF ancestors: always (family, default), "
        "never (exact), or only when max_pid >= 99 (legacy-split-pid; audit only)",
    )
    ap.add_argument(
        "--val-m8",
        type=Path,
        default=None,
        help="override the validation DIAMOND alignment (audit only; by default it is built/reused under a "
        "cache name keyed by the validation IDs and the training DIAMOND DB)",
    )
    ap.add_argument("--out-dir", type=Path, default=None, help="directory for the output TSVs (default: results/)")
    return ap.parse_args(argv)


def ensure_validation_alignment(val_pool: list[str], train_seqs: dict[str, str], train_db: Path) -> Path:
    """
    DIAMOND (--more-sensitive -k 1000) alignment of exactly this script's validation pool against the
    training DIAMOND DB. The cache file name is keyed by a hash of the sorted validation IDs and of the
    DB identity (file name + size), so an alignment built for another protein set (e.g. the
    `val_600_vs_train_sensitive.m8` of run_cafa_time_and_homology_splits.py, which shares only 8 of
    these 600 proteins) is never reused.

    The DB holds the full pre-t0 training corpus, exactly as for `wg_all7_vs_train.m8` (the test-genome
    queries). Per-split purging is applied downstream in run_diamond_knn by restricting hits to the
    split's sub-training set (which also removes the validation proteins themselves) and to
    pident <= max_pid, identically for validation and test queries.
    """
    ids_sorted = sorted(val_pool)
    h = hashlib.sha1()
    h.update("\n".join(ids_sorted).encode())
    h.update(f"|{train_db.name}|{train_db.stat().st_size}".encode())
    key = h.hexdigest()[:12]
    m8 = CACHE_DIR / f"val{len(ids_sorted)}_{key}_vs_{train_db.stem}.m8"
    val_set = set(ids_sorted)
    if m8.exists():
        with open(m8) as fh:
            q_seen = {ln.split("\t", 1)[0] for ln in fh if ln.strip()}
        if q_seen <= val_set:
            log(f"Reusing validation alignment {m8.name} ({len(q_seen)}/{len(val_set)} validation proteins with hits)")
            return m8
        log(f"Validation alignment {m8.name} has foreign queries; rebuilding")
    val_fasta = CACHE_DIR / f"val{len(ids_sorted)}_{key}.fasta"
    with open(val_fasta, "w") as fh:
        for p in ids_sorted:
            fh.write(f">{p}\n{train_seqs[p]}\n")
    cmd = [
        "diamond", "blastp",
        "-q", str(val_fasta),
        "-d", str(train_db),
        "-o", str(m8) + ".tmp",
        "--more-sensitive", "-k", "1000",
        "--outfmt", "6", "qseqid", "sseqid", "bitscore", "pident",
    ]
    log(f"Building validation alignment: {' '.join(cmd)}")
    t0 = time.time()
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    os.replace(str(m8) + ".tmp", m8)
    with open(m8) as fh:
        q_seen = {ln.split("\t", 1)[0] for ln in fh if ln.strip()}
    log(f"Built {m8.name} in {time.time() - t0:.1f}s ({len(q_seen)}/{len(val_set)} validation proteins with hits)")
    return m8

# ==============================================================================
def main():
    import gc

    args = parse_args()
    out_dir = args.out_dir if args.out_dir is not None else RESULTS_DIR
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    log(f"Starting Whole-Genome Withheld CAFA Benchmark on device={dev}")
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    out_dir.mkdir(parents=True, exist_ok=True)
    log(f"Options: hom_scale={args.hom_scale} fam_candidates={args.fam_candidates} out_dir={out_dir}")

    dag_path = CAFA_DIR / "go-dag.tsv"
    obo_path = CAFA_DIR / "go.obo"
    ia_path = CAFA_DIR / "kaggle-official" / "IA.tsv"
    tt_path = CAFA_DIR / "kaggle-official" / "train_terms.tsv"
    ref_dir = WG_DATA / "reference"

    anc, desc, aspect_of, mf_bp_links = load_dag_aspects_and_cross_links(dag_path, obo_path)
    ia = load_ia(ia_path)
    log(
        f"Loaded GO DAG ({len(anc):,} terms), IA weights ({len(ia):,} terms), "
        f"and {sum(len(v) for v in mf_bp_links.values()):,} MF->BP part_of links"
    )

    ec2go, pathways, heteromeric_complexes, essential_go = load_reference_pathways_complexes_essentials(
        ref_dir, anc, desc
    )
    log(
        f"Loaded Reference DBs: {len(ec2go):,} EC2GO mappings, {len(pathways):,} KEGG pathways/modules, "
        f"{len(heteromeric_complexes):,} multi-subunit complexes"
    )

    train_labels = load_train_terms(tt_path, anc, aspect_of)
    all_train_prots = set(train_labels.keys())
    log(f"Loaded Pre-t0 Swiss-Prot Training Corpus ({len(all_train_prots):,} annotated proteins)")

    (
        genome_prots,
        prot_to_genome,
        wg_prot2tax,
        genome_operons,
        genome_gaps,
        gt_all_exp,
        gt_post_t0,
        wg_ipr,
        wg_direct_ipr2go,
    ) = load_whole_genome_structures(WG_DATA, anc, desc, aspect_of, ec2go, train_labels)

    all_wg_prots = sorted(prot_to_genome.keys())
    all_wg_set = set(all_wg_prots)
    log(
        f"Loaded 6 Withheld Reference Genomes: {len(all_wg_prots):,} total coding sequences, "
        f"{len(gt_all_exp):,} experimentally annotated target proteins, "
        f"{len(gt_post_t0):,} strictly post-t0 newly annotated target proteins, "
        f"{sum(len(v) for v in genome_operons.values()):,} operons, "
        f"{sum(len(v) for v in genome_gaps.values()):,} metabolic gaps"
    )
    for ginfo in TEST_GENOMES:
        _t = ginfo["tag"]
        log(
            f"  {_t}: {len(genome_prots[_t]):,} ORFs, "
            f"{sum(1 for p in genome_prots[_t] if p in gt_all_exp):,} exp targets, "
            f"{sum(1 for p in genome_prots[_t] if p in gt_post_t0):,} post-t0 targets, "
            f"{len(genome_operons[_t]):,} operons (>=2 members in the proteome), "
            f"{len(genome_gaps[_t]):,} gapseq gaps (GO-mappable)"
        )

    # Load sequences for training set and 6 withheld genomes
    train_seqs = read_fasta(RECON_DIR / "cnn_work" / "train.fasta")
    wg_seqs = read_fasta(WG_DATA / "wg_all7.fasta")
    log(f"Loaded FASTAs: {len(train_seqs):,} train seqs, {len(wg_seqs):,} whole-genome seqs")

    # Load ESM-2 embeddings for training set and 6 withheld genomes (prefer ESM-2 3B 2,560-dim when available)
    tr_3b_path = RECON_DIR / "emb" / "train_esm2_3b.npz"
    wg_3b_path = WG_DATA / "wg_all7_esm2_3b.npz"
    if tr_3b_path.exists() and wg_3b_path.exists():
        tr_npz = np.load(tr_3b_path, allow_pickle=True)
        wg_npz = np.load(wg_3b_path, allow_pickle=True)
        esm_tag = "3b"
    else:
        tr_npz = np.load(RECON_DIR / "emb" / "train_esm2_35m.npz", allow_pickle=True)
        wg_npz = np.load(WG_DATA / "wg_all7_esm2_35m.npz", allow_pickle=True)
        esm_tag = "35m"

    tr_ids = [str(x) for x in tr_npz["ids"]]
    tr_id_to_idx = {p: i for i, p in enumerate(tr_ids)}
    wg_emb_ids = [str(x) for x in wg_npz["ids"]]
    wg_id_to_idx = {p: i for i, p in enumerate(wg_emb_ids)}

    all_caches_present = all(
        (CACHE_DIR / f"{ctag}_neural_{esm_tag}_wg20697.pkl").exists()
        for ctag in ("wg_time_rollover", "wg_pid100", "wg_pid70", "wg_pid50", "wg_pid30")
    )
    if not all_caches_present:
        tr_emb = tr_npz["emb"].astype(np.float32)
        tr_emb_norm = tr_emb / (np.linalg.norm(tr_emb, axis=1, keepdims=True) + 1e-8)
        wg_emb = wg_npz["emb"].astype(np.float32)
        wg_emb_norm = wg_emb / (np.linalg.norm(wg_emb, axis=1, keepdims=True) + 1e-8)
        log(f"Loaded ESM-2 ({esm_tag.upper()}) embeddings: train={tr_emb.shape}, whole_genomes={wg_emb.shape}")
    else:
        tr_emb = tr_emb_norm = wg_emb = wg_emb_norm = None
        log(f"All 5 split-retrained GPU prediction caches ({esm_tag.upper()}) present; skipping raw embedding matrix load")
    del tr_npz, wg_npz

    # Load taxon & anatomical constraint engine
    train_prot2tax, get_forbidden_for_taxon = load_taxon_and_anatomical_constraints(anc)
    combined_prot2tax = dict(train_prot2tax)
    combined_prot2tax.update(wg_prot2tax)

    # Identify all training proteins belonging to the 6 withheld test genomes or their species/strain taxids
    withheld_genome_train_prots = {
        p for p in all_train_prots if p in all_wg_set or train_prot2tax.get(p, "") in TEST_TAXIDS
    }
    log(
        f"Identified {len(withheld_genome_train_prots):,} training proteins belonging to the 6 withheld test genomes/species"
    )

    # Select deterministic validation set of 600 pre-t0 proteins strictly outside the 6 withheld test genomes
    rng_val = np.random.RandomState(42)
    eligible_val = sorted(
        [
            p
            for p in (all_train_prots - withheld_genome_train_prots - set(gt_post_t0.keys()))
            if p in train_seqs and p in tr_id_to_idx and len(train_labels[p]) >= 2
        ]
    )
    val_pool = list(rng_val.choice(eligible_val, size=600, replace=False))
    if args.val_m8 is not None:
        val_m8 = args.val_m8
        with open(val_m8) as fh:
            q_seen = {ln.split("\t", 1)[0] for ln in fh if ln.strip()}
        log(
            f"WARNING: using override validation alignment {val_m8} "
            f"({len(q_seen & set(val_pool))}/{len(val_pool)} validation proteins covered)"
        )
    else:
        val_m8 = ensure_validation_alignment(val_pool, train_seqs, RECON_DIR / "cnn_work" / "train_db.dmnd")

    prot_ipr, direct_interpro = {}, defaultdict(dict)
    if not all_caches_present:
        # InterPro domain annotations for train + val + 6 withheld genomes (only needed to retrain InterPro-LR)
        for ipr_file in (RECON_DIR / "interpro" / "train_ipr.tsv", RECON_DIR / "interpro" / "val_test.tsv"):
            if ipr_file.exists():
                with open(ipr_file) as fh:
                    for line in fh:
                        p = line.rstrip("\n").split("\t")
                        if len(p) >= 2 and p[1]:
                            prot_ipr[p[0]] = [x for x in p[1].split(";") if x]
        prot_ipr.update(wg_ipr)
        dip_path = RECON_DIR / "cleaneval_comps" / "interpro.tsv"
        if dip_path.exists():
            val_set = set(val_pool)
            with open(dip_path) as fh:
                for line in fh:
                    p = line.rstrip("\n").split("\t")
                    if len(p) >= 3 and p[0] in val_set:
                        try:
                            direct_interpro[p[0]][p[1]] = float(p[2])
                        except ValueError:
                            pass
        for pid, dmap in wg_direct_ipr2go.items():
            direct_interpro[pid].update(dmap)

    if all_caches_present:
        del train_seqs, wg_seqs
        gc.collect()

    # Load DIAMOND alignments for all 20,697 whole-genome proteins + validation proteins
    q2hits, train_to_test = load_diamond_hits(
        [WG_DATA / "wg_all7_vs_train.m8", val_m8]
    )
    log(
        f"Loaded DIAMOND alignments: {len(q2hits):,} whole-genome/val queries matched against {len(train_to_test):,} train targets"
    )

    # Define the 6 Whole-Genome Splits
    splits = [
        {
            "split_id": "WG_Time_Rollover",
            "split_label": "Whole-Genome Time-Based Split (t0->t1 Prospective, Pre-t0 Rollover Context)",
            "withhold_genomes_from_train": False,
            "max_pid": 100.0,
            "max_cos": 0.9999,
            "gt": gt_post_t0,
            "rollover_gt": train_labels,
            "train_cache_tag": "wg_time_rollover",
        },
        {
            "split_id": "WG_Time_Withheld",
            "split_label": "Whole-Genome Time-Based Split (t0->t1 Prospective, All 6 Genomes Withheld)",
            "withhold_genomes_from_train": True,
            "max_pid": 100.0,
            "max_cos": 0.9999,
            "gt": gt_post_t0,
            "rollover_gt": None,
            "train_cache_tag": "wg_pid100",
        },
        {
            "split_id": "WG_Homology_le100",
            "split_label": "Whole-Genome Homology Split <= 100% (Leave-Whole-Genome-Out, All 6 Genomes Withheld)",
            "withhold_genomes_from_train": True,
            "max_pid": 100.0,
            "max_cos": 0.9999,
            "gt": gt_all_exp,
            "rollover_gt": None,
            "train_cache_tag": "wg_pid100",
        },
        {
            "split_id": "WG_Homology_le70",
            "split_label": "Whole-Genome Homology Split <= 70% Seq Identity (All 6 Genomes Withheld + Homologs > 70% Purged)",
            "withhold_genomes_from_train": True,
            "max_pid": 70.0,
            "max_cos": 0.96,
            "gt": gt_all_exp,
            "rollover_gt": None,
            "train_cache_tag": "wg_pid70",
        },
        {
            "split_id": "WG_Homology_le50",
            "split_label": "Whole-Genome Homology Split <= 50% Seq Identity (All 6 Genomes Withheld + Homologs > 50% Purged)",
            "withhold_genomes_from_train": True,
            "max_pid": 50.0,
            "max_cos": 0.92,
            "gt": gt_all_exp,
            "rollover_gt": None,
            "train_cache_tag": "wg_pid50",
        },
        {
            "split_id": "WG_Homology_le30",
            "split_label": "Whole-Genome Homology Split <= 30% Seq Identity (All 6 Genomes Withheld + Homologs > 30% Purged)",
            "withhold_genomes_from_train": True,
            "max_pid": 30.0,
            "max_cos": 0.86,
            "gt": gt_all_exp,
            "rollover_gt": None,
            "train_cache_tag": "wg_pid30",
        },
    ]

    all_rows = []
    per_species_rows = []
    only_cache_tag = os.environ.get("GSPA_CACHE_TAG", "")

    for sp in splits:
        split_id = sp["split_id"]
        split_label = sp["split_label"]
        max_pid = sp["max_pid"]
        max_cos = sp["max_cos"]
        gt_dict = sp["gt"]
        rollover_gt = sp["rollover_gt"]
        cache_tag = sp["train_cache_tag"]
        if only_cache_tag and cache_tag != only_cache_tag:
            continue

        # 1. Construct filtered training set D_train^(s)
        excluded_prots = set()
        if sp["withhold_genomes_from_train"]:
            excluded_prots.update(withheld_genome_train_prots)
        else:
            excluded_prots.update(gt_dict.keys())

        if max_pid < 100.0:
            for tr_p, q_map in train_to_test.items():
                for q_p, pid in q_map.items():
                    if q_p in all_wg_set and pid > max_pid:
                        excluded_prots.add(tr_p)
                        break

        train_prots = all_train_prots - excluded_prots
        val_prots = [p for p in val_pool if p in train_prots]
        sub_train_prots = train_prots - set(val_prots)
        eval_prots = val_prots + all_wg_prots

        log(
            f"\n=== Split: {split_id} ({split_label}) | "
            f"Train: {len(train_prots):,} (Sub-train: {len(sub_train_prots):,}, Val: {len(val_prots)}, Purged: {len(excluded_prots):,}) | "
            f"Whole-Genome ORFs: {len(all_wg_prots):,} | Eval Targets: {len(gt_dict):,} ==="
        )

        gpu_only = args.gpu_only
        neural_cache_file = CACHE_DIR / f"{cache_tag}_neural_{esm_tag}_wg20697.pkl"
        if neural_cache_file.exists():
            log(f"  [{split_id}] Loading cached split-retrained GPU predictions from {neural_cache_file.name}")
            if gpu_only:
                continue
            with open(neural_cache_file, "rb") as fh:
                cached = pickle.load(fh)
            preds_cnn = cached["cnn"]
            preds_ebce = cached["ebce"]
            preds_emcm = cached["emcm"]
            preds_ipr = cached["ipr"]
            preds_eknn = cached["eknn"]
            del cached
        else:
            # Build the split-retrained neural cache (GPU recommended; run with --gpu-only on the cluster).
            eval_seqs = {p: train_seqs[p] for p in val_prots if p in train_seqs}
            eval_seqs.update({p: wg_seqs[p] for p in all_wg_prots if p in wg_seqs})
            val_emb_list = [p for p in val_prots if p in tr_id_to_idx]
            wg_emb_list = [p for p in all_wg_prots if p in wg_id_to_idx]
            eval_emb_ids = val_emb_list + wg_emb_list
            eval_emb = np.vstack(
                [tr_emb[[tr_id_to_idx[p] for p in val_emb_list]], wg_emb[[wg_id_to_idx[p] for p in wg_emb_list]]]
            )
            eval_emb_norm = np.vstack(
                [tr_emb_norm[[tr_id_to_idx[p] for p in val_emb_list]], wg_emb_norm[[wg_id_to_idx[p] for p in wg_emb_list]]]
            )
            prev_cnn_cache = CACHE_DIR / f"{cache_tag}_neural_wg20697.pkl"
            if prev_cnn_cache.exists():
                log(f"  [{split_id}] Reusing 16-kernel DeepGOCNN from {prev_cnn_cache.name}")
                with open(prev_cnn_cache, "rb") as fh:
                    preds_cnn = pickle.load(fh)["cnn"]
            else:
                log(f"  [{split_id}] Retraining DeepGOCNN on {len(sub_train_prots):,} proteins ...")
                preds_cnn = train_and_predict_deepgocnn(
                    sub_train_prots, train_labels, train_seqs, eval_prots, eval_seqs, dev
                )
            log(f"  [{split_id}] Retraining InterPro-LR on {len(sub_train_prots):,} proteins ...")
            preds_ipr = run_interpro_lr(sub_train_prots, train_labels, eval_prots, prot_ipr, direct_interpro, dev=dev)
            log(f"  [{split_id}] Retraining ESM-2 ({esm_tag.upper()}) MLP (BCE) & C-HMCNN (MCLoss) ...")
            preds_ebce, preds_emcm = train_and_predict_esm2_heads_batched(
                sub_train_prots, train_labels, tr_ids, tr_emb, eval_emb_ids, eval_emb, anc, dev
            )
            log(f"  [{split_id}] Running ESM-2 ({esm_tag.upper()}) embedding kNN ...")
            preds_eknn = run_esm2_knn_fast(
                sub_train_prots, train_labels, tr_ids, tr_emb_norm, eval_emb_ids, eval_emb_norm, max_cos_sim=max_cos
            )
            with open(neural_cache_file, "wb") as fh:
                pickle.dump(
                    {"cnn": preds_cnn, "ebce": preds_ebce, "emcm": preds_emcm, "ipr": preds_ipr, "eknn": preds_eknn}, fh
                )
            del eval_emb, eval_emb_norm
            if gpu_only:
                log(f"  [{split_id}] Saved GPU predictions to {neural_cache_file.name}")
                continue

        # 2. Retrain / run Naive and DIAMOND-kNN baselines across eval targets + whole-genome ORFs
        preds_naive = run_naive_frequency(sub_train_prots, train_labels, list(gt_dict.keys()))
        preds_diam = run_diamond_knn(sub_train_prots, train_labels, eval_prots, q2hits, max_pident=max_pid)

        channels = {
            "diam": preds_diam,
            "ipr_lr": preds_ipr,
            "cnn": preds_cnn,
            "esm2_knn": preds_eknn,
            "esm2_bce": preds_ebce,
            "esm2_mcm": preds_emcm,
        }

        log(f"  [{split_id}] Running Per-Protein Stacker/Stage0 + Whole-Genome GSPA Ablations across 6 complete genomes ...")
        (
            prop_raw,
            preds_stage0,
            preds_stacker,
            preds_operons,
            preds_coherence,
            preds_fba_gaps,
            preds_constraints,
            preds_gspa_full,
            gspa_info,
        ) = run_whole_genome_gspa_ablations(
            channels=channels,
            val_prots=val_prots,
            all_wg_prots=all_wg_prots,
            eval_target_prots=set(gt_dict.keys()),
            genome_prots=genome_prots,
            genome_operons=genome_operons,
            genome_gaps=genome_gaps,
            prot_to_taxon=wg_prot2tax,
            anc=anc,
            desc=desc,
            aspect_of=aspect_of,
            ia=ia,
            mf_bp_links=mf_bp_links,
            pathways=pathways,
            heteromeric_complexes=heteromeric_complexes,
            essential_go=essential_go,
            sub_train_prots=sub_train_prots,
            train_labels=train_labels,
            train_prot2tax=train_prot2tax,
            get_forbidden_for_taxon=get_forbidden_for_taxon,
            max_pid=max_pid,
            rollover_known_gt=rollover_gt,
            hom_scale_mode=args.hom_scale,
            fam_mode=args.fam_candidates,
            log_prefix=f"[{split_id}] ",
        )
        hom_sel = gspa_info["hom_scale_selected"]
        log(f"  [{split_id}] Stacker/Stage-0 blend beta: {gspa_info['beta_by_asp']} | operon family expansion: {gspa_info['fam_expansion']}")

        # Tune DeepGOPlus alpha_a in [0, 1] per aspect strictly on `val_prots` using already-propagated prop_raw
        val_gt_dict = {p: train_labels[p] for p in val_prots if p in train_labels}
        alpha_by_aspect = {}
        for asp in ("MF", "BP", "CC"):
            best_a, best_fw_a = 0.5, -1.0
            for a_cand in (0.0, 0.15, 0.30, 0.45, 0.60, 0.75, 0.90, 1.0):
                comb_val = combine_deepgoplus(prop_raw["cnn"], prop_raw["diam"], val_prots, alpha=a_cand)
                res_v = evaluate_predictions(
                    comb_val, val_gt_dict, anc, aspect_of, ia, train_prot2tax, get_forbidden_for_taxon, already_propagated=True
                )
                if res_v[f"fw_{asp}"] > best_fw_a:
                    best_fw_a = res_v[f"fw_{asp}"]
                    best_a = a_cand
            alpha_by_aspect[asp] = best_a
        log(f"  [{split_id}] Validation-tuned DeepGOPlus alpha: {alpha_by_aspect}")
        preds_dgp = combine_deepgoplus(
            prop_raw["cnn"], prop_raw["diam"], list(gt_dict.keys()), alpha=0.5, alpha_by_aspect=alpha_by_aspect, aspect_of=aspect_of
        )

        sota_cache_file = CACHE_DIR / f"{cache_tag}_sota_wg20697.pkl"
        sota_cached = {}
        if sota_cache_file.exists():
            with open(sota_cache_file, "rb") as fh:
                sota_cached = pickle.load(fh)
        pre_dfri_file = CACHE_DIR / "pretrained_deepfri_wg20697.pkl"
        pre_dfri = {}
        if pre_dfri_file.exists():
            with open(pre_dfri_file, "rb") as fh:
                pre_dfri = pickle.load(fh)

        sota_entries = []
        if split_id in ("WG_Time_Rollover", "WG_Time_Withheld"):
            if pre_dfri:
                sota_entries.append(("DeepFRI (Pretrained Official)", pre_dfri, False))
            if "deepgose_pretrained" in sota_cached:
                sota_entries.append(("DeepGO-SE (Pretrained Official)", sota_cached["deepgose_pretrained"], True))
        for s_key, s_label in (
            ("deepfri_retrained", "DeepFRI (GCN + LM, Retrained)"),
            ("deepgraphgo_retrained", "DeepGraphGO (InterPro + GraphGCN, Retrained)"),
            ("tale_plus_retrained", "TALE+ (Transformer + GO Label Attn + DIAMOND, Retrained)"),
            ("sprofgo_retrained", "SPROF-GO (PLM + Hierarchy + Label Diffusion, Retrained)"),
            ("deepgose_retrained", "DeepGO-SE (ESM-2 3B + ELEmbeddings, Retrained)"),
            ("dpfunc_retrained", "DPFunc (Domain Cross-Attn + GCN, Retrained)"),
        ):
            if s_key in sota_cached:
                sota_entries.append((s_label, sota_cached[s_key], True))

        method_map = [
            ("Naive Term-Frequency Prior (CAFA Baseline 1)", preds_naive, False),
            ("DIAMOND BLAST-kNN (CAFA Baseline 2)", prop_raw["diam"], True),
            ("InterPro Domain Logistic Regression (InterPro-LR)", prop_raw["ipr_lr"], True),
            ("DeepGOCNN (1D-CNN Sequence Model, Retrained)", prop_raw["cnn"], True),
            ("DeepGOPlus (1D-CNN + DIAMOND-kNN, Retrained)", preds_dgp, True),
            ("ESM-2 Embedding kNN (goPredSim, Retrained)", prop_raw["esm2_knn"], True),
            ("ESM-2 MLP (Flat BCE, Retrained)", prop_raw["esm2_bce"], True),
            ("ESM-2 C-HMCNN (Hierarchical MCLoss, Retrained)", prop_raw["esm2_mcm"], True),
            *sota_entries,
            ("GSPA Stage 0 (Per-Protein Noisy-OR)", preds_stage0, True),
            ("GSPA-Seq (Ours, Per-Protein Integrator)", preds_stacker, True),
            ("GSPA + Operon Genomic Context Only (B6)", preds_operons, True),
            ("GSPA + Pathway & Complex Coherence Only (B1/B2)", preds_coherence, True),
            ("GSPA + Metabolic Gap-Filling & FBA Only (B5 + Stage 2B)", preds_fba_gaps, True),
            ("GSPA + Taxon & Anatomical Constraints Only (B3/B4 + Stage 2A)", preds_constraints, True),
            ("GSPA Full Whole-Genome (Stages 0 + 1 + 2A + 2B)", preds_gspa_full, True),
        ]

        print(
            f"\n{'Method':<60} | {'fw_MF':>6} {'fw_BP':>6} {'fw_CC':>6} {'fw_Mean':>7} | "
            f"{'Fmax_MF':>7} {'Fmax_BP':>7} {'Fmax_CC':>7} {'Fmax_Mn':>7} | {'Smin_Mn':>7} | {'Viol':>5} {'V@.50':>5}"
        )
        print("-" * 139)
        for m_name, m_preds, is_prop in method_map:
            res = evaluate_predictions(
                m_preds, gt_dict, anc, aspect_of, ia, combined_prot2tax, get_forbidden_for_taxon, already_propagated=is_prop
            )
            print(
                f"{m_name:<60} | {res['fw_MF']:6.4f} {res['fw_BP']:6.4f} {res['fw_CC']:6.4f} {res['fw_mean']:7.4f} | "
                f"{res['fmax_MF']:7.4f} {res['fmax_BP']:7.4f} {res['fmax_CC']:7.4f} {res['fmax_mean']:7.4f} | "
                f"{res['smin_mean']:7.4f} | {res['violations']:5d} {res['violations_050']:5d}"
            )
            all_rows.append(
                {
                    "split_id": split_id,
                    "split_label": split_label,
                    "max_pident_cutoff": max_pid,
                    "n_train_proteins": len(train_prots),
                    "n_excluded_proteins": len(excluded_prots),
                    "n_wg_orfs": len(all_wg_prots),
                    "n_eval_targets": len(gt_dict),
                    "method": m_name,
                    "fw_MF": f"{res['fw_MF']:.4f}",
                    "fw_BP": f"{res['fw_BP']:.4f}",
                    "fw_CC": f"{res['fw_CC']:.4f}",
                    "fw_mean": f"{res['fw_mean']:.4f}",
                    "fmax_MF": f"{res['fmax_MF']:.4f}",
                    "fmax_BP": f"{res['fmax_BP']:.4f}",
                    "fmax_CC": f"{res['fmax_CC']:.4f}",
                    "fmax_mean": f"{res['fmax_mean']:.4f}",
                    "smin_MF": f"{res['smin_MF']:.4f}",
                    "smin_BP": f"{res['smin_BP']:.4f}",
                    "smin_CC": f"{res['smin_CC']:.4f}",
                    "smin_mean": f"{res['smin_mean']:.4f}",
                    "violations": res["violations"],
                    "violations_050": res["violations_050"],
                    "hom_scale_selected": hom_sel if m_name.startswith("GSPA") else "NA",
                }
            )

        # Also record per-genome summary for key methods on each split
        key_methods = [
            ("DIAMOND BLAST-kNN (CAFA Baseline 2)", prop_raw["diam"], True),
            ("DeepGOPlus (1D-CNN + DIAMOND-kNN, Retrained)", preds_dgp, True),
            ("ESM-2 C-HMCNN (Hierarchical MCLoss, Retrained)", prop_raw["esm2_mcm"], True),
            *sota_entries,
            ("GSPA-Seq (Ours, Per-Protein Integrator)", preds_stacker, True),
            ("GSPA + Operon Genomic Context Only (B6)", preds_operons, True),
            ("GSPA + Pathway & Complex Coherence Only (B1/B2)", preds_coherence, True),
            ("GSPA + Metabolic Gap-Filling & FBA Only (B5 + Stage 2B)", preds_fba_gaps, True),
            ("GSPA Full Whole-Genome (Stages 0 + 1 + 2A + 2B)", preds_gspa_full, True),
        ]
        for ginfo in TEST_GENOMES:
            gtag = ginfo["tag"]
            g_gt = {p: gt_dict[p] for p in genome_prots[gtag] if p in gt_dict}
            if len(g_gt) < 5:
                continue
            for km_entry in key_methods:
                m_name, m_preds = km_entry[0], km_entry[1]
                is_prop = km_entry[2] if len(km_entry) > 2 else True
                gres = evaluate_predictions(
                    m_preds, g_gt, anc, aspect_of, ia, combined_prot2tax, get_forbidden_for_taxon, already_propagated=is_prop
                )
                per_species_rows.append(
                    {
                        "split_id": split_id,
                        "genome": gtag,
                        "n_orfs": len(genome_prots[gtag]),
                        "n_eval_targets": len(g_gt),
                        "method": m_name,
                        "fw_MF": f"{gres['fw_MF']:.4f}",
                        "fw_BP": f"{gres['fw_BP']:.4f}",
                        "fw_CC": f"{gres['fw_CC']:.4f}",
                        "fw_mean": f"{gres['fw_mean']:.4f}",
                        "fmax_MF": f"{gres['fmax_MF']:.4f}",
                        "fmax_BP": f"{gres['fmax_BP']:.4f}",
                        "fmax_CC": f"{gres['fmax_CC']:.4f}",
                        "fmax_mean": f"{gres['fmax_mean']:.4f}",
                        "smin_MF": f"{gres['smin_MF']:.4f}",
                        "smin_BP": f"{gres['smin_BP']:.4f}",
                        "smin_CC": f"{gres['smin_CC']:.4f}",
                        "smin_mean": f"{gres['smin_mean']:.4f}",
                        "violations": gres["violations"],
                        "violations_050": gres["violations_050"],
                        "hom_scale_selected": hom_sel if m_name.startswith("GSPA") else "NA",
                    }
                )

        del preds_cnn, preds_ebce, preds_emcm, preds_ipr, preds_eknn, channels, prop_raw, sota_cached
        gc.collect()

    if args.gpu_only:
        log("Completed all GPU neural training splits (--gpu-only mode).")
        return

    out_tsv = out_dir / "cafa_whole_genome_time_and_homology_splits.tsv"
    fieldnames = [
        "split_id",
        "split_label",
        "max_pident_cutoff",
        "n_train_proteins",
        "n_excluded_proteins",
        "n_wg_orfs",
        "n_eval_targets",
        "method",
        "fw_MF",
        "fw_BP",
        "fw_CC",
        "fw_mean",
        "fmax_MF",
        "fmax_BP",
        "fmax_CC",
        "fmax_mean",
        "smin_MF",
        "smin_BP",
        "smin_CC",
        "smin_mean",
        "violations",
        "violations_050",
        "hom_scale_selected",
    ]
    with open(out_tsv, "w") as fh:
        fh.write("\t".join(fieldnames) + "\n")
        for row in all_rows:
            fh.write("\t".join(str(row[k]) for k in fieldnames) + "\n")
    log(f"Saved {len(all_rows)} rows to {out_tsv}")

    sp_tsv = out_dir / "cafa_whole_genome_per_species.tsv"
    sp_fields = [
        "split_id",
        "genome",
        "n_orfs",
        "n_eval_targets",
        "method",
        "fw_MF",
        "fw_BP",
        "fw_CC",
        "fw_mean",
        "fmax_MF",
        "fmax_BP",
        "fmax_CC",
        "fmax_mean",
        "smin_MF",
        "smin_BP",
        "smin_CC",
        "smin_mean",
        "violations",
        "violations_050",
        "hom_scale_selected",
    ]
    with open(sp_tsv, "w") as fh:
        fh.write("\t".join(sp_fields) + "\n")
        for row in per_species_rows:
            fh.write("\t".join(str(row[k]) for k in sp_fields) + "\n")
    log(f"Saved {len(per_species_rows)} per-genome rows to {sp_tsv}")


if __name__ == "__main__":
    main()
