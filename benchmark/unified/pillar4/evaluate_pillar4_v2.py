#!/usr/bin/env python3
"""Evaluate the Stage 3 agent and its baselines on benchmark v2 (and on LRO-50).

Inputs (all regenerable):
  results/pillar4_v2_benchmark.jsonl                  build_pillar4_dataset_v2.py assemble
  results/pillar4_v2_resolutions_<cfg>_<seed>.jsonl   run_pillar4_v2_agent.py (c0, c2, c3 x3, c4)
  results/pillar4_lro_50_cases.jsonl                  unchanged LRO-50 cases
  results/pillar4_v2_lro_resolutions_<cfg>_<seed>.jsonl

Outputs:
  results/pillar4_v2_evaluation.tsv       same columns as pillar4_agent_evaluation.tsv + benchmark_version
  results/pillar4_v2_per_category.tsv     per-category correct counts for every method
  results/pillar4_v2_ml_baselines.tsv     supervised / threshold / single-feature baselines (+ LRO-50 columns)
  results/pillar4_v2_multiseed_summary.tsv
  results/pillar4_v2_mcnemar.tsv          exact McNemar, agent (c3 seed 42) vs every other method, both benchmarks
  results/pillar4_v2_lro_reranking.tsv    per-case LRO-50 re-ranking for c3 seed 42
  results/pillar4_v2_pg_tp_hit1.tsv       PATHWAY_GAP TPs: correct protein chosen by argmax features / agents

Scoring (v2; written out because it differs from the v1 evaluator):
  TRUE_POSITIVE_PROMOTABLE   correct = PROMOTE of a held-out protein (any held-out protein that is
                             a candidate) with the target term. PROMOTE of another protein = FP.
  TRUE_NEGATIVE_DECOY        correct = ABSTAIN or DEMOTE. PROMOTE = FP.
  BORDERLINE_CONFLICT_REAL   (injected cross-clade term) correct = DEMOTE, or RELAX to a GO
                             ancestor of the target that passes the guardrail. PROMOTE = FP and a
                             violation. ABSTAIN is counted incorrect (column `accuracy`); the
                             lenient variant that accepts ABSTAIN is reported in per_category.
  BORDERLINE_CONFLICT_TRAP   (real, taxon-valid, truth-supported term) correct = ABSTAIN or
                             PROMOTE of the same term on the same protein (keep). DEMOTE/RELAX
                             = over-correction (incorrect).
  precision = TP promotions / all PROMOTE actions on items where a promotion can be judged
  (TP, decoy, genuine conflict); recall = TP promotions / #TP items; accuracy over all 120.
  violations_accepted = final PROMOTE/RELAX whose term fails SymbolicGuardrail.check_go_term.

Baselines use labels only to fit models inside training folds (5-fold StratifiedKFold over
anomalies, stratified by category x label, seed 42); post-processing is a fixed rule:
  no guardrail:  PROMOTE argmax candidate if p >= 0.5 else ABSTAIN
  + guardrail:   if the target term violates the organism's constraints: RELAX to the nearest
                 valid ancestor when the model would promote (DEMOTE if none), DEMOTE otherwise;
                 else as without guardrail.
"""

from __future__ import annotations

import argparse
import csv
import json
import statistics
import sys
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple

REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT))
from benchmark.agent.run_gspa_agent import SymbolicGuardrail, RESOURCES_DIR  # noqa: E402
from benchmark.unified.pillar4.evaluate_pillar4 import (  # noqa: E402
    _cand_feature_vector, evaluate_lro_reranking, evaluate_lro_single_feature_baselines, exact_mcnemar,
)

RESULTS = REPO_ROOT / "results"
SEED_TEMPS = [(42, 0.0), (43, 0.25), (44, 0.25)]
CFG_LABEL = {
    "c0": "Direct Zero-Shot LLM (No Tools, No Guardrail)",
    "c2": "Stage 3 Agent — Fast 1-Step Mode",
    "c3": "Stage 3 Agent — Qwen 27B + Guardrail",
    "c4": "Stage 3 Agent — Qwen 27B Without Guardrail",
}


def load_jsonl(p: Path) -> List[Dict[str, Any]]:
    return [json.loads(l) for l in open(p, encoding="utf-8") if l.strip()]


def acceptable(a: Dict[str, Any]) -> set:
    s = set(a.get("acceptable_protein_ids") or [])
    if a.get("expected_protein_id"):
        s.add(a["expected_protein_id"])
    return s


# --------------------------------------------------------------------------------------
# scoring
# --------------------------------------------------------------------------------------

def judge(a: Dict[str, Any], r: Dict[str, Any], g: SymbolicGuardrail) -> Dict[str, Any]:
    lab = a["ground_truth_label"]
    act = (r.get("action") or "ABSTAIN").upper()
    pid = r.get("protein_id") or ""
    fid = r.get("function_id") or a["target_function_id"]
    tax = a["organism_taxon"]
    viol = act in ("PROMOTE", "RELAX") and not g.check_go_term(fid, tax)["valid"]
    out = {"correct": False, "tp": False, "fp": False, "promote": act == "PROMOTE", "violation": viol,
           "lenient": False, "action": act}
    if lab == "TRUE_POSITIVE_PROMOTABLE":
        if act == "PROMOTE" and pid in acceptable(a) and fid == a["target_function_id"]:
            out.update(correct=True, tp=True)
        elif act == "PROMOTE":
            out["fp"] = True
    elif lab == "TRUE_NEGATIVE_DECOY":
        if act in ("ABSTAIN", "DEMOTE"):
            out["correct"] = True
        elif act == "PROMOTE":
            out["fp"] = True
    elif lab == "BORDERLINE_CONFLICT_REAL":
        t = a["target_function_id"]
        if act == "DEMOTE":
            out["correct"] = True
        elif act == "RELAX" and fid != t and fid in g.go_ancestors(t) and g.check_go_term(fid, tax)["valid"]:
            out["correct"] = True
        elif act == "PROMOTE":
            out["fp"] = True
        out["lenient"] = out["correct"] or act == "ABSTAIN"
    elif lab == "BORDERLINE_CONFLICT_TRAP":
        if act == "ABSTAIN" or (act == "PROMOTE" and fid == a["target_function_id"] and pid in acceptable(a)):
            out["correct"] = True
    if lab != "BORDERLINE_CONFLICT_REAL":
        out["lenient"] = out["correct"]
    return out


def evaluate(anoms: List[Dict[str, Any]], res: List[Dict[str, Any]], g: SymbolicGuardrail, name: str) -> Dict[str, Any]:
    by = {r["anomaly_id"]: r for r in res}
    n = len(anoms)
    js, per = [], {}
    tot_lat = tot_tc = 0.0
    for a in anoms:
        r = by.get(a["anomaly_id"], {"action": "ABSTAIN", "protein_id": "", "function_id": a["target_function_id"]})
        j = judge(a, r, g)
        js.append(j)
        tot_lat += float(r.get("latency_sec", 0.0) or 0.0)
        tot_tc += float(r.get("tool_calls", 0) or 0)
        k = (a["benchmark_category"], a["ground_truth_label"])
        d = per.setdefault(k, {"n": 0, "correct": 0, "lenient": 0, "promote": 0, "fp": 0})
        d["n"] += 1
        d["correct"] += j["correct"]
        d["lenient"] += j["lenient"]
        d["promote"] += j["promote"]
        d["fp"] += j["fp"]
    tp = sum(j["tp"] for j in js)
    fp = sum(j["fp"] for j in js)
    n_tp_items = sum(1 for a in anoms if a["ground_truth_label"] == "TRUE_POSITIVE_PROMOTABLE")
    prec = tp / (tp + fp) if tp + fp else 0.0
    rec = tp / n_tp_items if n_tp_items else 0.0
    f1 = 2 * prec * rec / (prec + rec) if prec + rec else 0.0
    acts = [j["action"] for j in js]
    viol = sum(j["violation"] for j in js)
    return {
        "config": name, "n_anomalies": n, "tp_promoted": tp, "fp_promoted": fp,
        "abstentions_count": acts.count("ABSTAIN"), "relaxations_count": acts.count("RELAX"),
        "demotions_count": acts.count("DEMOTE"), "promotions_count": acts.count("PROMOTE"),
        "precision": round(prec, 4), "recall": round(rec, 4), "f1": round(f1, 4),
        "accuracy": round(sum(j["correct"] for j in js) / n, 4) if n else 0.0,
        "accuracy_lenient": round(sum(j["lenient"] for j in js) / n, 4) if n else 0.0,
        "violations_accepted": viol, "violation_rate_pct": round(100.0 * viol / n, 2) if n else 0.0,
        "mean_latency_sec": round(tot_lat / n, 3) if n else 0.0, "mean_tool_calls": round(tot_tc / n, 2) if n else 0.0,
        "correct_vector": [j["correct"] for j in js], "per_category": per,
    }


# --------------------------------------------------------------------------------------
# baselines
# --------------------------------------------------------------------------------------

def postprocess(a: Dict[str, Any], best_pid: str, best_p: float, use_guard: bool, g: SymbolicGuardrail,
                thr: float = 0.5) -> Dict[str, Any]:
    t, tax = a["target_function_id"], a["organism_taxon"]
    promote = best_p >= thr and bool(best_pid)
    base = {"anomaly_id": a["anomaly_id"], "protein_id": best_pid if promote else "", "function_id": t,
            "confidence": round(float(best_p), 4), "tool_calls": 0, "latency_sec": 0.0}
    if use_guard and not g.check_go_term(t, tax)["valid"]:
        if promote:
            ancs = g.find_nearest_valid_ancestors(t, tax)
            if ancs:
                return {**base, "action": "RELAX", "function_id": ancs[0]}
        return {**base, "action": "DEMOTE", "protein_id": a["candidates"][0]["protein_id"] if a["candidates"] else ""}
    return {**base, "action": "PROMOTE" if promote else "ABSTAIN"}


def cand_label(a: Dict[str, Any], c: Dict[str, Any]) -> int:
    lab = a["ground_truth_label"]
    if lab in ("TRUE_POSITIVE_PROMOTABLE", "BORDERLINE_CONFLICT_TRAP"):
        return int(c["protein_id"] in acceptable(a))
    return 0


def strat_key(a: Dict[str, Any]) -> str:
    return f"{a['benchmark_category']}|{a['ground_truth_label']}"


def supervised_baselines(anoms: List[Dict[str, Any]], g: SymbolicGuardrail) -> List[Tuple[str, List[Dict[str, Any]]]]:
    import numpy as np
    from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import StratifiedKFold
    from sklearn.neural_network import MLPClassifier
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler

    models = {
        "Logistic Regression": lambda: make_pipeline(StandardScaler(), LogisticRegression(max_iter=1000, random_state=42)),
        "Random Forest": lambda: RandomForestClassifier(n_estimators=200, max_depth=5, random_state=42),
        "GBDT": lambda: GradientBoostingClassifier(n_estimators=100, max_depth=3, random_state=42),
        "MLP": lambda: make_pipeline(StandardScaler(), MLPClassifier(hidden_layer_sizes=(32, 16), max_iter=2000, random_state=42)),
    }
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    y_strat = [strat_key(a) for a in anoms]
    out = []
    for feat_name, pre3 in (("3 Pre-Tool Scalars", True), ("8 Post-Tool Features", False)):
        for mname, factory in models.items():
            scores: Dict[str, List[Tuple[str, float]]] = {}
            for tr, te in skf.split(anoms, y_strat):
                X, y = [], []
                for i in tr:
                    for c in anoms[i]["candidates"]:
                        X.append(_cand_feature_vector(c, anoms[i], pre_tool_3=pre3))
                        y.append(cand_label(anoms[i], c))
                clf = factory().fit(np.asarray(X, float), np.asarray(y, int))
                for i in te:
                    a = anoms[i]
                    pr = clf.predict_proba(np.asarray([_cand_feature_vector(c, a, pre_tool_3=pre3) for c in a["candidates"]], float))[:, 1]
                    scores[a["anomaly_id"]] = [(c["protein_id"], float(p)) for c, p in zip(a["candidates"], pr)]
            for use_guard in (False, True):
                res = []
                for a in anoms:
                    pid, p = max(scores[a["anomaly_id"]], key=lambda x: x[1]) if scores[a["anomaly_id"]] else ("", 0.0)
                    res.append(postprocess(a, pid, p, use_guard, g))
                name = f"{mname} (5-Fold Stratified CV, {feat_name})" + (" + Guardrail" if use_guard else "")
                out.append((name, res))
    return out


SINGLE_FEATURES: Dict[str, Callable[[Dict[str, Any]], float]] = {
    "current_posterior": lambda c: float(c.get("current_posterior", 0.0)),
    "final_logodds": lambda c: float(c.get("final_logodds", 0.0)),
    "suggestion_score": lambda c: float(c.get("suggestion_score", 0.0)),
}


def single_feature_baselines(anoms: List[Dict[str, Any]], g: SymbolicGuardrail) -> List[Tuple[str, List[Dict[str, Any]]]]:
    """Argmax-candidate by one feature; PROMOTE threshold chosen on the training folds (max accuracy)."""
    from sklearn.model_selection import StratifiedKFold
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    y_strat = [strat_key(a) for a in anoms]
    out = []
    for fname, f in SINGLE_FEATURES.items():
        best = {}
        for a in anoms:
            c = max(a["candidates"], key=f) if a["candidates"] else None
            best[a["anomaly_id"]] = (c["protein_id"], f(c)) if c else ("", float("-inf"))
        for use_guard in (False, True):
            res_by: Dict[str, Dict[str, Any]] = {}
            for tr, te in skf.split(anoms, y_strat):
                cands = sorted({best[anoms[i]["anomaly_id"]][1] for i in tr})
                grid = [cands[0] - 1e-9] + [(x + y) / 2 for x, y in zip(cands, cands[1:])] + [cands[-1] + 1e-9]
                def acc(thr):
                    return sum(judge(anoms[i], postprocess(anoms[i], *best[anoms[i]["anomaly_id"]], use_guard, g, thr), g)["correct"] for i in tr)
                thr = max(grid, key=lambda t: (acc(t), -abs(t)))
                for i in te:
                    a = anoms[i]
                    res_by[a["anomaly_id"]] = postprocess(a, *best[a["anomaly_id"]], use_guard, g, thr)
            name = f"Single Feature: {fname} (argmax candidate, CV-tuned threshold)" + (" + Guardrail" if use_guard else "")
            out.append((name, [res_by[a["anomaly_id"]] for a in anoms]))
    return out


def threshold_baseline(anoms: List[Dict[str, Any]], g: SymbolicGuardrail) -> List[Dict[str, Any]]:
    res = []
    for a in anoms:
        c = max(a["candidates"], key=lambda c: float(c["current_posterior"])) if a["candidates"] else None
        res.append(postprocess(a, c["protein_id"] if c else "", float(c["current_posterior"]) if c else 0.0, True, g, 0.35))
    return res


def lro_supervised(lro: List[Dict[str, Any]]) -> Dict[str, Dict[str, float]]:
    import numpy as np
    from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import KFold
    from sklearn.neural_network import MLPClassifier
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    models = {
        "Logistic Regression": lambda: make_pipeline(StandardScaler(), LogisticRegression(max_iter=1000, random_state=42)),
        "Random Forest": lambda: RandomForestClassifier(n_estimators=200, max_depth=5, random_state=42),
        "GBDT": lambda: GradientBoostingClassifier(n_estimators=100, max_depth=3, random_state=42),
        "MLP": lambda: make_pipeline(StandardScaler(), MLPClassifier(hidden_layer_sizes=(32, 16), max_iter=2000, random_state=42)),
    }
    kf = KFold(n_splits=5, shuffle=True, random_state=42)
    out = {}
    for feat_name, pre3 in (("3 Pre-Tool Scalars", True), ("8 Post-Tool Features", False)):
        for mname, factory in models.items():
            res = []
            for tr, te in kf.split(lro):
                X = [_cand_feature_vector(c, lro[i], pre_tool_3=pre3) for i in tr for c in lro[i]["candidates"]]
                y = [int(c["protein_id"] == lro[i]["true_protein_id"]) for i in tr for c in lro[i]["candidates"]]
                clf = factory().fit(np.asarray(X, float), np.asarray(y, int))
                for i in te:
                    pr = clf.predict_proba(np.asarray([_cand_feature_vector(c, lro[i], pre_tool_3=pre3) for c in lro[i]["candidates"]], float))[:, 1]
                    k = int(np.argmax(pr))
                    res.append({"anomaly_id": lro[i]["anomaly_id"], "action": "PROMOTE", "protein_id": lro[i]["candidates"][k]["protein_id"]})
            s = evaluate_lro_reranking(lro, res)
            hv = {r["anomaly_id"]: r["agent_position"] == 1 for r in s["case_rows"]}
            out[f"{mname} (5-Fold CV, {feat_name})"] = {"hit1_pct": s["agent_hit1_pct"], "mrr": s["agent_mrr"],
                                                       "hits": [hv[c["anomaly_id"]] for c in lro]}
    return out


# --------------------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--benchmark", type=Path, default=RESULTS / "pillar4_v2_benchmark.jsonl")
    ap.add_argument("--lro", type=Path, default=RESULTS / "pillar4_lro_50_cases.jsonl")
    ap.add_argument("--results-dir", type=Path, default=RESULTS)
    a = ap.parse_args()
    R = a.results_dir
    g = SymbolicGuardrail(RESOURCES_DIR)
    anoms = load_jsonl(a.benchmark)
    lro = load_jsonl(a.lro)

    def res_file(cfg, seed, lro_=False):
        return R / f"pillar4_v2_{'lro_' if lro_ else ''}resolutions_{cfg}_{seed}.jsonl"

    rows: List[Dict[str, Any]] = []
    rows.append(evaluate(anoms, [{"anomaly_id": x["anomaly_id"], "action": "ABSTAIN"} for x in anoms], g,
                         "Stage 1 + Stage 2 Deterministic Only"))
    agent: Dict[Tuple[str, int], Dict[str, Any]] = {}
    for cfg, seeds in (("c0", [42]), ("c2", [42]), ("c3", [42, 43, 44]), ("c4", [42])):
        for s in seeds:
            f = res_file(cfg, s)
            if not f.exists():
                print(f"[warn] missing {f}", file=sys.stderr)
                continue
            m = evaluate(anoms, load_jsonl(f), g, CFG_LABEL[cfg] + (f" (seed={s})" if cfg == "c3" else ""))
            agent[(cfg, s)] = m
    order = [("c2", 42), ("c3", 42), ("c3", 43), ("c3", 44), ("c4", 42), ("c0", 42)]
    rows += [agent[k] for k in order if k in agent]
    base_rows = [evaluate(anoms, threshold_baseline(anoms, g), g, "Posterior Threshold (>=0.35) + Guardrail")]
    for name, res in supervised_baselines(anoms, g) + single_feature_baselines(anoms, g):
        base_rows.append(evaluate(anoms, res, g, name))
    rows += base_rows

    cols = ["configuration", "n_anomalies", "tp_promoted", "fp_promoted", "abstentions", "relaxations", "demotions",
            "precision", "recall", "f1", "accuracy", "violations_accepted", "violation_rate_pct", "mean_latency_sec",
            "mean_tool_calls", "benchmark_version"]
    with open(R / "pillar4_v2_evaluation.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(cols)
        for m in rows:
            w.writerow([m["config"], m["n_anomalies"], m["tp_promoted"], m["fp_promoted"], m["abstentions_count"],
                        m["relaxations_count"], m["demotions_count"], m["precision"], m["recall"], m["f1"], m["accuracy"],
                        m["violations_accepted"], m["violation_rate_pct"], m["mean_latency_sec"], m["mean_tool_calls"], "v2"])
    with open(R / "pillar4_v2_per_category.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["configuration", "category", "label", "n", "correct", "correct_lenient", "promotions", "false_promotions", "benchmark_version"])
        for m in rows:
            for (cat, lab), d in sorted(m["per_category"].items()):
                w.writerow([m["config"], cat, lab, d["n"], d["correct"], d["lenient"], d["promote"], d["fp"], "v2"])

    # ---- LRO-50 ----
    lro_sum: Dict[Tuple[str, int], Dict[str, Any]] = {}
    for cfg, seeds in (("c0", [42]), ("c2", [42]), ("c3", [42, 43, 44]), ("c4", [42])):
        for s in seeds:
            f = res_file(cfg, s, True)
            if f.exists():
                lro_sum[(cfg, s)] = evaluate_lro_reranking(lro, load_jsonl(f))
    lro_hits = {k: [r["agent_position"] == 1 for r in v["case_rows"]] for k, v in lro_sum.items()}
    feats = evaluate_lro_single_feature_baselines(lro)
    lro_sup = lro_supervised(lro)
    rlgc = evaluate_lro_reranking(lro, [])

    with open(R / "pillar4_v2_ml_baselines.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["method", "anom120_precision", "anom120_recall", "anom120_f1", "anom120_accuracy",
                    "anom120_violation_rate_pct", "lro50_hit1_pct", "lro50_mrr", "benchmark_version"])
        for m in base_rows:
            lr = None
            for k, v in lro_sup.items():
                if m["config"].startswith(k.split(" (")[0]) and k.split(", ")[1] in m["config"] and "Guardrail" not in m["config"]:
                    lr = v
            w.writerow([m["config"], m["precision"], m["recall"], m["f1"], m["accuracy"], m["violation_rate_pct"],
                        lr["hit1_pct"] if lr else "n/a", lr["mrr"] if lr else "n/a", "v2"])
        for (cfg, s), m in agent.items():
            ls = lro_sum.get((cfg, s))
            w.writerow([m["config"], m["precision"], m["recall"], m["f1"], m["accuracy"], m["violation_rate_pct"],
                        ls["agent_hit1_pct"] if ls else "n/a", ls["agent_mrr"] if ls else "n/a", "v2"])
        for f in feats:
            w.writerow([f["method"], "n/a", "n/a", "n/a", "n/a", "n/a", f["hit1_pct"], f["mrr"], "v2"])
        w.writerow(["RLGC ordered baseline (rlgc_ordered_pids)", "n/a", "n/a", "n/a", "n/a", "n/a",
                    rlgc["rlgc_hit1_pct"], rlgc["rlgc_mrr"], "v2"])

    # ---- multi-seed ----
    ms = []
    for s, t in SEED_TEMPS:
        m, ls = agent.get(("c3", s)), lro_sum.get(("c3", s))
        if m:
            ms.append({"seed_label": f"seed_{s}", "temperature": t, "precision": m["precision"], "recall": m["recall"],
                       "f1": m["f1"], "accuracy": m["accuracy"], "violation_rate_pct": m["violation_rate_pct"],
                       "lro_hit1_pct": ls["agent_hit1_pct"] if ls else "", "lro_mrr": ls["agent_mrr"] if ls else ""})
    if len(ms) >= 2:
        def msd(k):
            v = [float(r[k]) for r in ms if r[k] != ""]
            return f"{statistics.mean(v):.4f}+/-{statistics.pstdev(v):.4f}" if v else ""
        ms.append({"seed_label": "MEAN_SD", "temperature": "", **{k: msd(k) for k in
                   ("precision", "recall", "f1", "accuracy", "violation_rate_pct", "lro_hit1_pct", "lro_mrr")}})
    if ms:
        with open(R / "pillar4_v2_multiseed_summary.tsv", "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(ms[0].keys()) + ["benchmark_version"], delimiter="\t")
            w.writeheader()
            for r in ms:
                w.writerow({**r, "benchmark_version": "v2"})

    # ---- McNemar ----
    mc = []
    ref_key = ("c3", 42)
    if ref_key in agent:
        ref = agent[ref_key]["correct_vector"]
        for m in rows:
            if m is agent[ref_key]:
                continue
            r = exact_mcnemar(ref, m["correct_vector"])
            mc.append(["ANOM_v2", "Stage 3 Agent c3 (seed 42)", m["config"], r["n"], r["hits_a"], r["hits_b"],
                       r["b"], r["c"], f"{r['p_value']:.6g}"])
    if ref_key in lro_hits:
        ref = lro_hits[ref_key]
        comps = [(f["method"], f["hit_vector"]) for f in feats]
        comps.append(("RLGC ordered baseline (rlgc_ordered_pids)", [r["rlgc_position"] == 1 for r in rlgc["case_rows"]]))
        comps += [(k, v["hits"]) for k, v in lro_sup.items()]
        comps += [(CFG_LABEL[c] + f" (seed={s})", h) for (c, s), h in lro_hits.items() if (c, s) != ref_key]
        for name, hv in comps:
            r = exact_mcnemar(ref, hv)
            mc.append(["LRO_50", "Stage 3 Agent c3 (seed 42)", name, r["n"], r["hits_a"], r["hits_b"], r["b"], r["c"],
                       f"{r['p_value']:.6g}"])
    with open(R / "pillar4_v2_mcnemar.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["benchmark", "method_a", "method_b", "n", "correct_a", "correct_b", "b_a_only", "c_b_only",
                    "p_exact_two_sided", "test", "benchmark_version"])
        for r in mc:
            w.writerow(r + ["exact McNemar (scipy.stats.binomtest, p=0.5)", "v2"])

    if ref_key in lro_sum:
        with open(R / "pillar4_v2_lro_reranking.tsv", "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh, delimiter="\t")
            keys = list(lro_sum[ref_key]["case_rows"][0].keys())
            w.writerow(keys + ["benchmark_version"])
            for r in lro_sum[ref_key]["case_rows"]:
                w.writerow([r[k] for k in keys] + ["v2"])

    # PATHWAY_GAP true positives: is the argmax candidate (by each single feature) or the agent's
    # promoted protein a correct held-out protein? Ties count as fractional hits (1/k).
    tp_items = [a_ for a_ in anoms if a_["ground_truth_label"] == "TRUE_POSITIVE_PROMOTABLE"]
    with open(R / "pillar4_v2_pg_tp_hit1.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["method", "n_tp_items", "hit1_tie_aware", "hit1_pct", "benchmark_version"])
        for fname, f in SINGLE_FEATURES.items():
            hits = 0.0
            for a_ in tp_items:
                vals = [f(c) for c in a_["candidates"]]
                top = max(vals)
                tied = [c for c, v in zip(a_["candidates"], vals) if v == top]
                hits += sum(c["protein_id"] in acceptable(a_) for c in tied) / len(tied)
            w.writerow([f"argmax {fname}", len(tp_items), round(hits, 3), round(100 * hits / len(tp_items), 1), "v2"])
        def ctx_score(a_, c):  # operon neighbours annotated with a present step of the pathway + in_pathway_operon
            pres = set(a_["context"].get("present_terms", []) or [])
            n = sum(1 for v in (c.get("operon_neighbors") or {}).values() for x in v if x[:10] in pres)
            return n + int(bool(c.get("in_pathway_operon")))
        hits = 0.0
        for a_ in tp_items:
            vals = [ctx_score(a_, c) for c in a_["candidates"]]
            tied = [c for c, v in zip(a_["candidates"], vals) if v == max(vals)]
            hits += sum(c["protein_id"] in acceptable(a_) for c in tied) / len(tied)
        w.writerow(["argmax operon context (neighbours with present pathway steps + in_pathway_operon)", len(tp_items),
                    round(hits, 3), round(100 * hits / len(tp_items), 1), "v2"])
        rnd = sum(len(acceptable(a_) & {c["protein_id"] for c in a_["candidates"]}) / len(a_["candidates"]) for a_ in tp_items)
        w.writerow(["uniform random candidate (expected)", len(tp_items), round(rnd, 3), round(100 * rnd / len(tp_items), 1), "v2"])
        for (cfg, sd), m in agent.items():
            res = {r["anomaly_id"]: r for r in load_jsonl(res_file(cfg, sd))}
            h = sum(1 for a_ in tp_items if res.get(a_["anomaly_id"], {}).get("action") == "PROMOTE"
                    and res[a_["anomaly_id"]].get("protein_id") in acceptable(a_))
            w.writerow([m["config"], len(tp_items), h, round(100 * h / len(tp_items), 1), "v2"])

    # console summary
    for m in rows:
        print(f"{m['config'][:72]:<72} P={m['precision']:.3f} R={m['recall']:.3f} F1={m['f1']:.3f} "
              f"Acc={m['accuracy']:.3f} (len {m['accuracy_lenient']:.3f}) viol={m['violations_accepted']} lat={m['mean_latency_sec']:.1f}")
    print(f"LRO-50 RLGC order: Hit@1 {rlgc['rlgc_hit1_pct']} MRR {rlgc['rlgc_mrr']}")
    for k, v in lro_sum.items():
        print(f"LRO-50 {k}: Hit@1 {v['agent_hit1_pct']} MRR {v['agent_mrr']}")
    for f in feats:
        print(f"LRO-50 {f['feature']}: Hit@1 {f['hit1_pct']} MRR {f['mrr']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
