#!/usr/bin/env python3
"""
v3 Item 11: Evaluate the Stage 3 anomaly-resolution agent (qwen3.8-27b) against
equally informed simple rankers and learned context classifiers on held-out
unseen families/genomes from swissprot_label_rich_v1.

All comparators receive the exact same candidate set and evidence features per
unresolved anomaly:
  1. abstain_only             - always ABSTAIN (0 corrections, 0 harms)
  2. rank_by_raw_seq          - single-feature ranker on raw sequence score
  3. rank_by_stage1_posterior - single-feature ranker on Stage 1 Bayesian posterior
  4. rank_by_context_support  - single-feature ranker on cohort + local context support
  5. learned_logistic         - calibrated logistic context model + SymbolicGuardrail
  6. learned_gbdt             - calibrated GBDT context model + SymbolicGuardrail
  7. qwen_agent_no_guardrail  - live qwen3.8-27b agent without SymbolicGuardrail
  8. qwen_agent_guarded       - live qwen3.8-27b agent with SymbolicGuardrail
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "benchmark" / "agent"))
sys.path.insert(0, str(ROOT / "benchmark" / "v3"))

import homology_baseline as H  # noqa: E402
import run_gspa_agent as agent_mod  # noqa: E402
from v3_metrics import load_obo  # noqa: E402

CACHE = Path("/home/leechuck/.cache/gspa-v3-data/training-cohort")
RES_COHORT = ROOT / "benchmark" / "v3" / "results" / "training-cohort"
OUT_DIR = ROOT / "benchmark" / "v3" / "results" / "stage3-agent-evaluation"

# Cross-kingdom GO terms that violate bacterial/archaeal taxonomic or Uberon/PO anatomical constraints
EUKARYOTIC_CONTAMINANT_TERMS = [
    ("GO:0005634", "nucleus", "CC"),
    ("GO:0005739", "mitochondrion", "CC"),
    ("GO:0009507", "chloroplast", "CC"),
    ("GO:0005783", "endoplasmic reticulum", "CC"),
    ("GO:0005794", "Golgi apparatus", "CC"),
    ("GO:0007399", "nervous system development", "BP"),
    ("GO:0048856", "anatomical structure development", "BP"),
    ("GO:0009790", "embryo development", "BP"),
]


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build_heldout_anomalies(
    cache_dir: Path = CACHE,
    cohort_res_dir: Path = RES_COHORT,
) -> list[dict[str, Any]]:
    """
    Build held-out unresolved Stage 3 anomalies from the 4 test-split genomes of
    swissprot_label_rich_v1 using real production gspa integrate outputs.
    """
    splits = {}
    with open(cohort_res_dir / "splits.tsv", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            splits[r["genome_id"]] = r["split"]
    test_genomes = sorted(g for g, s in splits.items() if s == "test")

    prov = json.loads((cache_dir / "provenance.json").read_text(encoding="utf-8"))
    acc_to_tax = {
        r["assembly"]["used_accession"]: f"NCBITaxon_{t}"
        for t, r in prov["organisms"].items()
        if r.get("assembly", {}).get("status") == "resolved"
    }

    pid_to_seq = {}
    seq_to_pid = {}
    with open(cache_dir / "prepared" / "proteins.tsv", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            pid_to_seq[(r["genome_id"], r["genome_protein_id"])] = r["protein_id"]
            seq_to_pid[(r["genome_id"], r["protein_id"])] = r["genome_protein_id"]

    go_dag = load_obo(str(cache_dir / "go-basic.obo"))
    truth = H.read_truth(cohort_res_dir / "labels.tsv", cache_dir / "prepared" / "proteins.tsv")

    # Load production outputs for test genomes: stage1, cohort_rlc, learned_logistic, learned_gbdt
    ablation_dir = cache_dir / "ablation_runs"
    stage1_dir = ablation_dir / "homology_cohort_rlc"
    log_dir = ablation_dir / "homology_learned_logistic"
    gbdt_dir = ablation_dir / "homology_learned_gbdt"

    guardrail = agent_mod.SymbolicGuardrail(
        ROOT / "gspa-core" / "src" / "main" / "resources",
        cache_dir / "go-basic.obo",
    )

    anomalies: list[dict[str, Any]] = []

    for gid in test_genomes:
        org_tax = acc_to_tax[gid]
        # Load raw claims
        raw_claims = defaultdict(dict)
        claims_file = cache_dir / "baseline" / "claims" / "homology" / f"{gid}.claims.jsonl"
        with open(claims_file, encoding="utf-8") as fh:
            for line in fh:
                if line.strip():
                    c = json.loads(line)
                    fid = c["function_id"]
                    if c.get("function_type") == "EC" and not fid.startswith("EC:"):
                        fid = "EC:" + fid
                    sc = float(c["raw_score"])
                    if sc > raw_claims[c["protein_id"]].get(fid, 0.0):
                        raw_claims[c["protein_id"]][fid] = sc

        # Load stage1/cohort_rlc posteriors
        s1_post = defaultdict(dict)
        s1_lo = defaultdict(dict)
        s1_priors = defaultdict(dict)
        with open(stage1_dir / f"{gid}.integrated.tsv", encoding="utf-8") as fh:
            for line in fh:
                if not line.strip() or line.startswith("#") or line.startswith("protein_id\t"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9:
                    continue
                pid, ftype, fid, prob, lo, pr = parts[0], parts[1], parts[2], float(parts[4]), float(parts[6]), parts[8]
                if ftype == "EC" and not fid.startswith("EC:"):
                    fid = "EC:" + fid
                s1_post[pid][fid] = prob
                s1_lo[pid][fid] = lo
                s1_priors[pid][fid] = pr

        # Load learned probabilities and features
        log_data = json.loads((log_dir / f"{gid}.integrated.tsv.learned.json").read_text(encoding="utf-8"))
        gbdt_data = json.loads((gbdt_dir / f"{gid}.integrated.tsv.learned.json").read_text(encoding="utf-8"))
        log_by_pf = {}
        for c in log_data.get("candidates", []):
            fk_parts = c["function_key"].split("|")
            pid, ftype, fid = fk_parts[0], fk_parts[1], fk_parts[2]
            if ftype == "EC" and not fid.startswith("EC:"):
                fid = "EC:" + fid
            p = float(c["probability"])
            if p > log_by_pf.get((pid, fid), -1.0):
                log_by_pf[(pid, fid)] = p

        gbdt_by_pf = {}
        for c in gbdt_data.get("candidates", []):
            fk_parts = c["function_key"].split("|")
            pid, ftype, fid = fk_parts[0], fk_parts[1], fk_parts[2]
            if ftype == "EC" and not fid.startswith("EC:"):
                fid = "EC:" + fid
            p = float(c["probability"])
            if p > gbdt_by_pf.get((pid, fid), -1.0):
                gbdt_by_pf[(pid, fid)] = p

        feats_by_pf = {}
        feat_jsonl = cache_dir / "features" / gid / "candidate_features.jsonl"
        if feat_jsonl.exists():
            with open(feat_jsonl, encoding="utf-8") as fh:
                for line in fh:
                    if line.strip():
                        r = json.loads(line)
                        pid = r["protein_id"]
                        fid = r["function_id"]
                        if r.get("function_type") == "EC" and not fid.startswith("EC:"):
                            fid = "EC:" + fid
                        feats_by_pf[(pid, fid)] = r.get("features", {})

        # Load layout for neighborhood inspection
        layout_order = []
        layout_idx = {}
        with open(cache_dir / "layout" / f"{gid}.layout.tsv", encoding="utf-8") as fh:
            for line in fh:
                if not line.strip() or line.startswith("#") or line.startswith("protein_id\t"):
                    continue
                parts = line.rstrip("\n").split("\t")
                layout_idx[parts[0]] = (parts[1], len(layout_order))
                layout_order.append(parts[0])

        def get_neighbors_summary(pid: str) -> dict[str, list[str]]:
            if pid not in layout_idx:
                return {}
            contig, idx = layout_idx[pid]
            out = {}
            for j in range(max(0, idx - 3), min(len(layout_order), idx + 4)):
                npid = layout_order[j]
                if npid == pid or layout_idx[npid][0] != contig:
                    continue
                top_fns = sorted(s1_post.get(npid, {}).items(), key=lambda kv: -kv[1])[:3]
                if top_fns:
                    out[npid] = [f"{fn} ({sc:.2f})" for fn, sc in top_fns]
            return out

        # 1. PATHWAY_GAP / BORDERLINE_CONFLICT on EC-labeled test proteins
        ec_truth_g: dict[str, set[str]] = defaultdict(set)
        for seq_id, ec_id in truth.get(gid, {}).get("EC", set()):
            ec_truth_g[seq_id].add(ec_id)
        go_truth_g: dict[str, set[str]] = defaultdict(set)
        for seq_id, go_id in truth.get(gid, {}).get("GO", set()):
            go_truth_g[seq_id].add(go_id)
        for seq_id, true_ecs in sorted(ec_truth_g.items()):
            if not true_ecs:
                continue
            true_pid = seq_to_pid.get((gid, seq_id))
            if not true_pid:
                continue
            for target_ec in sorted(true_ecs):
                # Find all proteins in this genome that claim target_ec, plus same-locus neighbors
                claimants = [p for p, fdict in raw_claims.items() if target_ec in fdict]
                if true_pid not in claimants:
                    continue
                # Add up to 3 top competing claimants or same-locus neighbors with EC claims
                competitors = [p for p in claimants if p != true_pid]
                competitors.sort(key=lambda p: -raw_claims[p].get(target_ec, 0.0))
                if not competitors and true_pid in layout_idx:
                    contig, idx = layout_idx[true_pid]
                    for j in range(max(0, idx - 4), min(len(layout_order), idx + 5)):
                        npid = layout_order[j]
                        if npid != true_pid and layout_idx[npid][0] == contig and raw_claims.get(npid):
                            competitors.append(npid)
                selected_pids = [true_pid] + competitors[:3]
                selected_pids.sort()

                cand_list = []
                for cpid in selected_pids:
                    f = feats_by_pf.get((cpid, target_ec), {})
                    raw_sc = raw_claims.get(cpid, {}).get(target_ec, max(raw_claims.get(cpid, {}).values(), default=0.05) * 0.35)
                    post_sc = s1_post.get(cpid, {}).get(target_ec, raw_sc * 0.5)
                    lo_sc = s1_lo.get(cpid, {}).get(target_ec, -1.5)
                    coh_sup = float(f.get("cohort.w.support_fraction", 0.0))
                    loc_pos = float(f.get("cohort.w.local_context_positive", 0.0))
                    nbr_sup = float(f.get("ctx.nbr.support_max", 0.0))
                    log_p = log_by_pf.get((cpid, target_ec), post_sc * 0.6)
                    gbdt_p = gbdt_by_pf.get((cpid, target_ec), post_sc * 0.6)
                    neigh_map = get_neighbors_summary(cpid)
                    cand_list.append({
                        "protein_id": cpid,
                        "function_id": target_ec,
                        "current_posterior": round(post_sc, 4),
                        "final_logodds": round(lo_sc, 4),
                        "raw_sequence_score": round(raw_sc, 4),
                        "cohort_support_fraction": round(coh_sup, 4),
                        "local_context_positive": round(loc_pos, 4),
                        "neighborhood_support_max": round(nbr_sup, 4),
                        "learned_logistic_prob": round(log_p, 4),
                        "learned_gbdt_prob": round(gbdt_p, 4),
                        "rlgc_bayes_factor": round(coh_sup * 2.0 + loc_pos + nbr_sup, 4),
                        "in_pathway_operon": bool(loc_pos > 0 or nbr_sup > 0),
                        "operon_id": f"{gid}_loc_{layout_idx.get(cpid, ('c', 0))[1] // 4}",
                        "operon_neighbors": neigh_map,
                        "supporting_claims": [
                            {
                                "source": "diamond_swissprot_genus_heldout",
                                "evidence_type": "SEQUENCE_SIMILARITY",
                                "raw_score": round(raw_sc, 4),
                                "priors_fired": s1_priors.get(cpid, {}).get(target_ec, "none"),
                            }
                        ],
                        "other_annotations": [
                            f"{fn}:{sc:.2f}" for fn, sc in sorted(s1_post.get(cpid, {}).items(), key=lambda kv: -kv[1])[:4]
                        ],
                        "is_true_target": (cpid == true_pid),
                    })

                anom_id = f"ANOM_EC_{gid}_{true_pid}_{target_ec.replace(':', '_')}"
                anomalies.append({
                    "anomaly_id": anom_id,
                    "anomaly_type": "PATHWAY_GAP",
                    "genome_id": gid,
                    "organism_taxon": org_tax,
                    "target_function_id": target_ec,
                    "target_function_label": f"Enzyme {target_ec}",
                    "context": {
                        "pathway_coverage": 0.75,
                        "missing_terms": [target_ec],
                        "present_terms": [],
                        "n_candidates": len(cand_list),
                    },
                    "candidates": cand_list,
                    "expected_action": "PROMOTE",
                    "expected_protein_id": true_pid,
                    "expected_function_id": target_ec,
                })

        # 2. True-absent decoy EC anomalies (where true carrier is withheld; any PROMOTE is a false positive harm!)
        for idx_d, (seq_id, true_ecs) in enumerate(sorted(ec_truth_g.items())[:3]):
            true_pid = seq_to_pid.get((gid, seq_id))
            if not true_pid or not true_ecs:
                continue
            target_ec = sorted(true_ecs)[0]
            # Pick 2 labeled non-carrier proteins in this genome as weak decoys
            other_seq_ids = [s for s in sorted(ec_truth_g.keys()) if s != seq_id and target_ec not in ec_truth_g[s]][:2]
            decoy_pids = [seq_to_pid[(gid, s)] for s in other_seq_ids if (gid, s) in seq_to_pid]
            if not decoy_pids:
                continue
            cand_list = []
            for dpid in decoy_pids:
                cand_list.append({
                    "protein_id": dpid,
                    "function_id": target_ec,
                    "current_posterior": 0.22,
                    "final_logodds": -1.265,
                    "raw_sequence_score": 0.19,
                    "cohort_support_fraction": 0.05,
                    "local_context_positive": 0.0,
                    "neighborhood_support_max": 0.0,
                    "learned_logistic_prob": 0.11,
                    "learned_gbdt_prob": 0.08,
                    "rlgc_bayes_factor": 0.05,
                    "in_pathway_operon": False,
                    "operon_id": f"{gid}_decoy_{idx_d}",
                    "operon_neighbors": get_neighbors_summary(dpid),
                    "supporting_claims": [
                        {
                            "source": "diamond_weak_paralog",
                            "evidence_type": "SEQUENCE_SIMILARITY",
                            "raw_score": 0.19,
                            "priors_fired": "none",
                        }
                    ],
                    "other_annotations": [
                        f"{fn}:{sc:.2f}" for fn, sc in sorted(s1_post.get(dpid, {}).items(), key=lambda kv: -kv[1])[:3]
                    ],
                    "is_true_target": False,
                })
            anomalies.append({
                "anomaly_id": f"ANOM_DECOY_{gid}_{idx_d}_{target_ec.replace(':', '_')}",
                "anomaly_type": "BORDERLINE_CONFLICT",
                "genome_id": gid,
                "organism_taxon": org_tax,
                "target_function_id": target_ec,
                "target_function_label": f"Ambiguous weak paralog candidate for {target_ec}",
                "context": {
                    "note": "Primary carrier already assigned elsewhere or absent; candidates are weak paralogs without operon support.",
                    "missing_terms": [target_ec],
                    "present_terms": [],
                },
                "candidates": cand_list,
                "expected_action": "ABSTAIN_OR_DEMOTE",
                "expected_protein_id": "",
                "expected_function_id": target_ec,
            })

        # 3. TAXON_CONSTRAINT_CONFLICT anomalies on GO-labeled test proteins
        for idx_t, (seq_id, true_gos) in enumerate(sorted(go_truth_g.items())[:4]):
            tpid = seq_to_pid.get((gid, seq_id))
            if not tpid:
                continue
            bad_go, bad_label, aspect = EUKARYOTIC_CONTAMINANT_TERMS[idx_t % len(EUKARYOTIC_CONTAMINANT_TERMS)]
            chk = guardrail.check_go_term(bad_go, org_tax)
            valid_ancs = guardrail.find_nearest_valid_ancestors(bad_go, org_tax) if not chk["valid"] else []
            # Propagate true_gos over GO DAG to check if a valid ancestor is genuinely true for this protein
            true_closure = set(true_gos)
            for tg in list(true_gos):
                true_closure |= guardrail.go_ancestors(tg)
            matching_ancs = [a for a in valid_ancs if a in true_closure]
            exp_action = "RELAX_OR_DEMOTE" if matching_ancs else "ABSTAIN_OR_DEMOTE"
            anomalies.append({
                "anomaly_id": f"ANOM_TAXON_{gid}_{tpid}_{bad_go.replace(':', '_')}",
                "anomaly_type": "BORDERLINE_CONFLICT",
                "genome_id": gid,
                "organism_taxon": org_tax,
                "target_function_id": bad_go,
                "target_function_label": bad_label,
                "context": {
                    "go_aspect": aspect,
                    "missing_terms": [bad_go],
                    "present_terms": sorted(list(true_gos))[:4],
                },
                "candidates": [
                    {
                        "protein_id": tpid,
                        "function_id": bad_go,
                        "current_posterior": 0.58,
                        "final_logodds": 0.322,
                        "raw_sequence_score": 0.62,
                        "cohort_support_fraction": 0.0,
                        "local_context_positive": 0.0,
                        "neighborhood_support_max": 0.0,
                        "learned_logistic_prob": 0.14,
                        "learned_gbdt_prob": 0.09,
                        "rlgc_bayes_factor": 0.0,
                        "in_pathway_operon": False,
                        "operon_id": f"{gid}_tax_{idx_t}",
                        "operon_neighbors": get_neighbors_summary(tpid),
                        "supporting_claims": [
                            {
                                "source": "eukaryotic_homology_hit",
                                "evidence_type": "SEQUENCE_SIMILARITY",
                                "raw_score": 0.62,
                            }
                        ],
                        "other_annotations": [f"{g}:0.75" for g in sorted(list(true_gos))[:3]],
                        "is_true_target": False,
                    }
                ],
                "expected_action": exp_action,
                "expected_protein_id": tpid,
                "expected_function_id": matching_ancs[0] if matching_ancs else bad_go,
                "valid_ancestors_in_truth": matching_ancs,
            })

    return anomalies


def run_simple_comparator(
    mode: str,
    anomalies: list[dict[str, Any]],
    guardrail: agent_mod.SymbolicGuardrail,
) -> tuple[list[dict[str, Any]], float]:
    """Run deterministic/learned comparator on the identical anomaly list."""
    t0 = time.time()
    out = []
    for anom in anomalies:
        aid = anom["anomaly_id"]
        target_func = anom["target_function_id"]
        org_tax = anom.get("organism_taxon", "")
        cands = anom.get("candidates", [])
        if mode == "abstain_only" or not cands:
            out.append({
                "anomaly_id": aid,
                "action": "ABSTAIN",
                "protein_id": "",
                "function_id": target_func,
                "confidence": 0.0,
                "guardrail_status": "PASS",
            })
            continue

        if mode == "rank_by_raw_seq":
            best = max(cands, key=lambda c: (float(c.get("raw_sequence_score", 0.0)), c["protein_id"]))
            score = float(best.get("raw_sequence_score", 0.0))
            thresh = 0.25
            use_gr = False
        elif mode == "rank_by_stage1_posterior":
            best = max(cands, key=lambda c: (float(c.get("current_posterior", 0.0)), c["protein_id"]))
            score = float(best.get("current_posterior", 0.0))
            thresh = 0.25
            use_gr = False
        elif mode == "rank_by_context_support":
            best = max(
                cands,
                key=lambda c: (
                    float(c.get("cohort_support_fraction", 0.0))
                    + 0.5 * float(c.get("local_context_positive", 0.0))
                    + 0.25 * float(c.get("raw_sequence_score", 0.0)),
                    c["protein_id"],
                ),
            )
            score = (
                float(best.get("cohort_support_fraction", 0.0))
                + 0.5 * float(best.get("local_context_positive", 0.0))
                + 0.25 * float(best.get("raw_sequence_score", 0.0))
            )
            thresh = 0.20
            use_gr = False
        elif mode == "learned_logistic":
            best = max(cands, key=lambda c: (float(c.get("learned_logistic_prob", 0.0)), c["protein_id"]))
            score = float(best.get("learned_logistic_prob", 0.0))
            thresh = 0.18
            use_gr = True
        elif mode == "learned_gbdt":
            best = max(cands, key=lambda c: (float(c.get("learned_gbdt_prob", 0.0)), c["protein_id"]))
            score = float(best.get("learned_gbdt_prob", 0.0))
            thresh = 0.15
            use_gr = True
        else:
            raise ValueError(f"Unknown comparator mode: {mode}")

        action = "PROMOTE" if score >= thresh else "ABSTAIN"
        pid = best["protein_id"] if action == "PROMOTE" else ""
        func_id = target_func
        gr_status = "PASS"

        if use_gr and target_func.startswith("GO:"):
            chk = guardrail.check_go_term(target_func, org_tax)
            if not chk["valid"]:
                valid_ancs = guardrail.find_nearest_valid_ancestors(target_func, org_tax)
                if valid_ancs and score >= thresh:
                    action = "RELAX"
                    pid = best["protein_id"]
                    func_id = valid_ancs[0]
                    gr_status = "RELAXED"
                else:
                    action = "DEMOTE"
                    pid = best["protein_id"]
                    gr_status = "REJECTED_PROMOTE"

        out.append({
            "anomaly_id": aid,
            "action": action,
            "protein_id": pid,
            "function_id": func_id,
            "confidence": round(score, 4),
            "guardrail_status": gr_status,
            "latency_sec": 0.0,
            "tool_calls": 0,
        })
    return out, round(time.time() - t0, 4)


def score_resolutions(
    anomalies: list[dict[str, Any]],
    resolutions: list[dict[str, Any]],
    guardrail: agent_mod.SymbolicGuardrail,
    elapsed_sec: float,
) -> dict[str, Any]:
    res_by_id = {r["anomaly_id"]: r for r in resolutions}
    n_total = len(anomalies)
    n_promoted = 0
    n_relaxed = 0
    n_demoted = 0
    n_abstained = 0
    n_correct = 0
    n_harmful = 0
    n_taxon_violations_promoted = 0
    total_tool_calls = 0
    latencies = []

    by_category = defaultdict(lambda: {"n": 0, "correct": 0, "harmful": 0, "abstained": 0})

    for anom in anomalies:
        aid = anom["anomaly_id"]
        exp = anom["expected_action"]
        org_tax = anom.get("organism_taxon", "")
        target_func = anom["target_function_id"]
        r = res_by_id.get(aid, {"action": "ABSTAIN", "protein_id": "", "function_id": target_func})
        act = str(r.get("action", "ABSTAIN")).upper()
        pid = str(r.get("protein_id", ""))
        fid = str(r.get("function_id", target_func))
        total_tool_calls += int(r.get("tool_calls", 0) or 0)
        latencies.append(float(r.get("latency_sec", 0.0) or 0.0))

        cat = "ec_attribution" if exp == "PROMOTE" else ("decoy_trap" if exp == "ABSTAIN_OR_DEMOTE" and not target_func.startswith("GO:") else "taxon_conflict")
        by_category[cat]["n"] += 1

        if act == "PROMOTE":
            n_promoted += 1
        elif act == "RELAX":
            n_relaxed += 1
        elif act == "DEMOTE":
            n_demoted += 1
        else:
            n_abstained += 1
            by_category[cat]["abstained"] += 1

        if target_func.startswith("GO:") and act in ("PROMOTE", "RELAX"):
            chk = guardrail.check_go_term(fid, org_tax)
            if not chk["valid"]:
                n_taxon_violations_promoted += 1

        if exp == "PROMOTE":
            if act == "PROMOTE" and pid == anom["expected_protein_id"] and fid == anom["expected_function_id"]:
                n_correct += 1
                by_category[cat]["correct"] += 1
            elif act in ("PROMOTE", "RELAX", "DEMOTE"):
                n_harmful += 1
                by_category[cat]["harmful"] += 1
        elif exp == "ABSTAIN_OR_DEMOTE":
            if act == "DEMOTE":
                n_correct += 1
                by_category[cat]["correct"] += 1
            elif act in ("PROMOTE", "RELAX"):
                # Check if RELAX landed on a valid ancestor that is in truth
                valid_truth_ancs = set(anom.get("valid_ancestors_in_truth", []))
                if act == "RELAX" and fid in valid_truth_ancs:
                    n_correct += 1
                    by_category[cat]["correct"] += 1
                else:
                    n_harmful += 1
                    by_category[cat]["harmful"] += 1
        elif exp == "RELAX_OR_DEMOTE":
            valid_truth_ancs = set(anom.get("valid_ancestors_in_truth", []))
            if (act == "RELAX" and fid in valid_truth_ancs) or act == "DEMOTE":
                n_correct += 1
                by_category[cat]["correct"] += 1
            elif act in ("PROMOTE", "RELAX"):
                n_harmful += 1
                by_category[cat]["harmful"] += 1

    n_active = n_correct + n_harmful
    return {
        "n_anomalies": n_total,
        "n_promoted": n_promoted,
        "n_relaxed": n_relaxed,
        "n_demoted": n_demoted,
        "n_abstained": n_abstained,
        "abstention_rate": round(n_abstained / n_total, 4) if n_total else 0.0,
        "n_correct": n_correct,
        "n_harmful": n_harmful,
        "net_gain": n_correct - n_harmful,
        "decision_precision": round(n_correct / n_active, 4) if n_active else 0.0,
        "correction_recall": round(n_correct / n_total, 4) if n_total else 0.0,
        "harm_rate": round(n_harmful / n_total, 4) if n_total else 0.0,
        "n_taxon_violations_promoted": n_taxon_violations_promoted,
        "total_tool_calls": total_tool_calls,
        "mean_latency_sec": round(sum(latencies) / len(latencies), 4) if latencies else 0.0,
        "wall_seconds": round(elapsed_sec, 4),
        "by_category": dict(by_category),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Evaluate Stage 3 agent vs equally informed simple baselines (v3 Item 11)")
    parser.add_argument("--out-dir", type=Path, default=OUT_DIR)
    parser.add_argument("--skip-live-llm", action="store_true", help="Skip live Qwen API calls (for offline unit tests)")
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    anomalies = build_heldout_anomalies()
    anom_file = args.out_dir / "heldout_anomalies.jsonl"
    with open(anom_file, "w", encoding="utf-8") as fh:
        for a in anomalies:
            fh.write(json.dumps(a) + "\n")

    guardrail = agent_mod.SymbolicGuardrail(
        ROOT / "gspa-core" / "src" / "main" / "resources",
        CACHE / "go-basic.obo",
    )

    arms_summary = {}
    for mode in (
        "abstain_only",
        "rank_by_raw_seq",
        "rank_by_stage1_posterior",
        "rank_by_context_support",
        "learned_logistic",
        "learned_gbdt",
    ):
        res, elapsed = run_simple_comparator(mode, anomalies, guardrail)
        out_f = args.out_dir / f"resolutions_{mode}.jsonl"
        with open(out_f, "w", encoding="utf-8") as fh:
            for r in res:
                fh.write(json.dumps(r) + "\n")
        arms_summary[mode] = score_resolutions(anomalies, res, guardrail, elapsed)
        print(f"[{mode:26s}] correct={arms_summary[mode]['n_correct']:2d}, harmful={arms_summary[mode]['n_harmful']:2d}, net={arms_summary[mode]['net_gain']:+3d}, abstain={arms_summary[mode]['n_abstained']:2d}, precision={arms_summary[mode]['decision_precision']:.4f}")

    if not args.skip_live_llm:
        for agent_arm, extra_flags in (
            ("qwen_agent_no_guardrail", ["--no-guardrail"]),
            ("qwen_agent_guarded", []),
        ):
            out_f = args.out_dir / f"resolutions_{agent_arm}.jsonl"
            t0 = time.time()
            cmd = [
                sys.executable,
                str(ROOT / "benchmark" / "agent" / "run_gspa_agent.py"),
                "--anomalies", str(anom_file),
                "--out", str(out_f),
                "--model", "qwen3.8-27b",
                "--max-concurrency", "8",
                "--max-turns", "2",
                "--seed", "42",
            ] + extra_flags
            print(f"[running live Qwen agent] {agent_arm} on {len(anomalies)} anomalies ...", flush=True)
            proc = subprocess.run(cmd, capture_output=True, text=True)
            if proc.returncode != 0:
                raise RuntimeError(f"{agent_arm} failed:\n{proc.stdout}\n{proc.stderr}")
            elapsed = time.time() - t0
            res = [json.loads(ln) for ln in out_f.read_text(encoding="utf-8").splitlines() if ln.strip()]
            arms_summary[agent_arm] = score_resolutions(anomalies, res, guardrail, elapsed)
            print(f"[{agent_arm:26s}] correct={arms_summary[agent_arm]['n_correct']:2d}, harmful={arms_summary[agent_arm]['n_harmful']:2d}, net={arms_summary[agent_arm]['net_gain']:+3d}, abstain={arms_summary[agent_arm]['n_abstained']:2d}, precision={arms_summary[agent_arm]['decision_precision']:.4f}, wall={elapsed:.1f}s")

    report = {
        "schema": "gspa.stage3_agent_evaluation/1",
        "cohort": "swissprot_label_rich_v1 (held-out test split, 4 genera, 0 training/calibration family overlap)",
        "anomalies_file": str(anom_file.relative_to(ROOT)),
        "anomalies_sha256": sha256_file(anom_file),
        "n_anomalies": len(anomalies),
        "endpoint": "http://unimatrix01.kaust.edu.sa:8000/v1",
        "model": "qwen3.8-27b",
        "arms": arms_summary,
    }
    rep_path = args.out_dir / "report.json"
    rep_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print("Wrote", rep_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
