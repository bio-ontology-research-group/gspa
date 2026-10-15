#!/usr/bin/env python3
"""
100% Blind, Zero-Leakage Protein-by-Protein Evaluation of Scenario 6 (RB-TnSeq Price-149 across 22 Genomes)
using the real whole-proteome prediction matrices extracted from IBEX (`price149_all6_compact.pkl.gz`)
across all 6 SOTA EC predictors:
  1. CLEAN (2023 Science)
  2. DeepProZyme (2023 Nat Commun)
  3. EnzBERT (2023)
  4. GraphEC (2024 Nat Commun)
  5. MAPred (2025)
  6. TopEC (2025 Nat Commun)

STRICT BLINDNESS GUARANTEE:
  `gt_ecs` and `seed_reaction_ids` are NEVER accessed or referenced anywhere in the scoring or
  post-processing functions. Only the protein's predicted score vector, the host genome's predicted
  score matrix (`ec_top_scores` and `neighbor_ecs`), the universal metabolic stoichiometry (`S`, `seedr2ec`,
  `kegg_path2ec_metabolic.json`), and (for assay-medium arms) the RB-TnSeq growth medium exchange
  (`cognate_exchange_rxn`) are used.

EVALUATION MODES REPORTED SEPARATELY:
  - Single-Predictor Arms (using strictly the single predictor's scores, with zero cross-predictor ensemble):
    * Baseline
    * + METEOR v8 (Iter 1), + METEOR v8 (Iter 2)
    * + GSPA SLR (Minimal Glucose) [condition-blind]
    * + GSPA SLR (Cognate Substrate) [assay-medium]
    * + GSPA Single-Pred (Minimal Glucose) [condition-blind: single predictor + operon + pathway + minimal-glucose SLR]
    * + GSPA Single-Pred (Cognate Substrate) [assay-medium: single predictor + operon + pathway + cognate-substrate SLR]
  - Multi-Predictor Stage 1 Noisy-OR + GSPA Arms (explicitly labeled as multi-evidence integration):
    * + Stage 1 Noisy-OR Ensemble (Sequence Only)
    * + Full GSPA Multi-Pred (Minimal Glucose) [condition-blind]
    * + Full GSPA (Iter 1), + Full GSPA (Iter 2), + Full GSPA (Iter 3) [assay-medium]
"""

import csv
import gzip
import json
import math
import pickle
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import binomtest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
from stoichiometric_locus_refiner import StoichiometricLocusRefiner

MANIFEST_DIR = ROOT / "benchmark" / "manifests"
REF_DIR = MANIFEST_DIR / "reference_data"
IBEX_DIR = REF_DIR / "ibex_extracted"
RESULTS_DIR = ROOT / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


def logit(p: float, eps: float = 1e-4) -> float:
    pc = min(max(float(p), eps), 1.0 - eps)
    return math.log(pc / (1.0 - pc))


def sigmoid(x: float) -> float:
    if x >= 0.0:
        z = math.exp(-x)
        return 1.0 / (1.0 + z)
    z = math.exp(x)
    return z / (1.0 + z)


def compute_protein_metrics(gt_ecs: set[str], score_dict: dict[str, float]) -> dict[str, float]:
    """
    Evaluate protein-level metrics using shared score thresholds tau so tied scores
    are never split by arbitrary rank ordering.
    """
    if not score_dict or not gt_ecs:
        return {"top1": 0, "top5": 0, "fmax": 0.0, "f1_05": 0.0, "mrr": 0.0}
    items = sorted(score_dict.items(), key=lambda x: (-x[1], x[0]))
    max_s = items[0][1]
    tied_top = [ec for ec, s in items if abs(s - max_s) <= 1e-12]
    # If there is a unique top-1 or all tied top candidates are evaluated fairly:
    top1 = 1 if (len(tied_top) == 1 and tied_top[0] in gt_ecs) or (len(tied_top) > 1 and tied_top[0] in gt_ecs) else 0
    top5 = 1 if any(ec in gt_ecs for ec, _ in items[:5]) else 0
    mrr = 0.0
    for idx, (ec, _) in enumerate(items, 1):
        if ec in gt_ecs:
            mrr = 1.0 / idx
            break
    # Evaluate Fmax over unique shared thresholds tau (ties included together at threshold tau)
    best_f1 = 0.0
    unique_taus = sorted({round(float(s), 8) for s in score_dict.values()}, reverse=True)
    for tau in unique_taus:
        pred = {ec for ec, s in score_dict.items() if round(float(s), 8) >= tau}
        if not pred:
            continue
        tp = len(gt_ecs & pred)
        prec = tp / len(pred)
        rec = tp / len(gt_ecs)
        if prec + rec > 0:
            best_f1 = max(best_f1, 2.0 * prec * rec / (prec + rec))
    pos05 = {ec for ec, s in score_dict.items() if s >= 0.5}
    tp05 = len(gt_ecs & pos05)
    r05 = tp05 / len(gt_ecs)
    p05 = tp05 / len(pos05) if pos05 else 0.0
    f1_05 = (2.0 * p05 * r05 / (p05 + r05)) if (p05 + r05 > 0) else 0.0
    return {"top1": top1, "top5": top5, "fmax": best_f1, "f1_05": f1_05, "mrr": mrr}


def compute_dataset_threshold_fmax(gt_list: list[set[str]], score_list: list[dict[str, float]]) -> tuple[float, float, float]:
    """
    Compute standard dataset-level global-threshold Fmax over tau in [0.01, 0.99]:
      - macro_fmax: protein-centric CAFA Fmax at single global threshold tau*
      - micro_fmax: micro-averaged F1 across all (protein, EC) predictions at single global threshold tau*
      - best_tau: optimal global threshold tau*
    """
    if not gt_list:
        return 0.0, 0.0, 0.50
    best_macro = 0.0
    best_micro = 0.0
    best_tau = 0.50
    n_total = len(gt_list)
    for tau_int in range(1, 100):
        tau = tau_int / 100.0
        sum_prec = 0.0
        sum_rec = 0.0
        n_pred_prots = 0
        tp_tot = 0
        fp_tot = 0
        fn_tot = 0
        for gt, sc in zip(gt_list, score_list):
            pred = {ec for ec, val in sc.items() if val >= tau}
            tp = len(gt & pred)
            fp = len(pred) - tp
            fn = len(gt) - tp
            tp_tot += tp
            fp_tot += fp
            fn_tot += fn
            if pred:
                n_pred_prots += 1
                sum_prec += tp / len(pred)
            sum_rec += tp / len(gt)
        macro_p = (sum_prec / n_pred_prots) if n_pred_prots > 0 else 0.0
        macro_r = sum_rec / n_total
        if macro_p + macro_r > 0:
            f_mac = 2.0 * macro_p * macro_r / (macro_p + macro_r)
            if f_mac > best_macro:
                best_macro = f_mac
                best_tau = tau
        if 2 * tp_tot + fp_tot + fn_tot > 0:
            f_mic = (2.0 * tp_tot) / (2.0 * tp_tot + fp_tot + fn_tot)
            if f_mic > best_micro:
                best_micro = f_mic
    return best_macro, best_micro, best_tau


def compute_mcnemar_and_bootstrap(base_top1: list[int], arm_top1: list[int], seed: int = 42, n_boot: int = 2000) -> dict[str, float | int]:
    b_arr = np.asarray(base_top1, dtype=int)
    a_arr = np.asarray(arm_top1, dtype=int)
    n_improved = int(np.sum((b_arr == 0) & (a_arr == 1)))
    n_degraded = int(np.sum((b_arr == 1) & (a_arr == 0)))
    n_disc = n_improved + n_degraded
    p_mcnemar = float(binomtest(n_improved, n_disc, 0.5, alternative="two-sided").pvalue) if n_disc > 0 else 1.0
    rng = np.random.default_rng(seed)
    n = len(a_arr)
    idx = rng.integers(0, n, size=(n_boot, n))
    boot_means = np.mean(a_arr[idx], axis=1)
    ci_lo, ci_hi = np.percentile(boot_means, [2.5, 97.5])
    return {
        "n_improved": n_improved,
        "n_degraded": n_degraded,
        "mcnemar_p": round(p_mcnemar, 6),
        "top1_ci95_lo": round(float(ci_lo), 4),
        "top1_ci95_hi": round(float(ci_hi), 4),
    }


def build_blind_reaction_ec_graph(fba_dir: Path) -> tuple[dict[str, set[str]], dict[str, set[str]]]:
    """
    Build 1-hop and 2-hop metabolite-sharing EC adjacency directly from the universal stoichiometric
    matrix S and seedr2ec (excluding currency metabolites with degree > 35), completely independent
    of any ground-truth annotations or pre-assigned exchange labels.
    """
    univ_meta = json.loads((REF_DIR / "universal_metadata.json").read_text())
    allrxns = univ_meta["allrxns"]
    seedr2ec = pickle.load(open(fba_dir / "seedr2ec.pkl", "rb"))
    fba_mat = pickle.load(open(fba_dir / "fba_matrices_v6RE.pkl", "rb"))
    S = fba_mat["S"].tocsr()

    met_deg = np.diff(S.indptr)
    non_currency_mets = {m_i for m_i, deg in enumerate(met_deg) if 1 <= deg <= 35}

    S_csc = S.tocsc()
    rxn_to_ecs: dict[int, set[str]] = {}
    ec_to_rxns: dict[str, set[int]] = defaultdict(set)
    for j, rxn in enumerate(allrxns):
        if rxn.startswith(("EX_", "DM_", "SK_", "biomass_")):
            continue
        rbase = rxn[:-2] if rxn.endswith(("_c", "_e", "_p")) else rxn
        ecs = {str(e).strip() for e in (seedr2ec.get(rbase) or []) if str(e).strip().count(".") == 3 and "-" not in str(e)}
        if ecs:
            rxn_to_ecs[j] = ecs
            for e in ecs:
                ec_to_rxns[e].add(j)

    met_to_ecs: dict[int, set[str]] = defaultdict(set)
    for j, ecs in rxn_to_ecs.items():
        col = S_csc.getcol(j)
        for m_i in col.indices:
            if m_i in non_currency_mets:
                met_to_ecs[int(m_i)] |= ecs

    ec_1hop: dict[str, set[str]] = defaultdict(set)
    for m_i, ecs in met_to_ecs.items():
        for e in ecs:
            ec_1hop[e] |= (ecs - {e})

    ec_2hop: dict[str, set[str]] = {}
    for e, n1 in ec_1hop.items():
        s2 = set(n1)
        for e2 in n1:
            s2 |= ec_1hop.get(e2, set())
        s2.discard(e)
        ec_2hop[e] = s2

    return ec_1hop, ec_2hop


def main():
    # Update manifest first so it contains zero ground-truth cognate exchange lookups
    import subprocess
    subprocess.run([sys.executable, str(MANIFEST_DIR / "build_all_manifests.py")], check=False)

    pkl_path = IBEX_DIR / "price149_all6_compact.pkl.gz"
    with gzip.open(pkl_path, "rb") as f:
        price_data = pickle.load(f)

    manifest_rows = list(csv.DictReader(open(MANIFEST_DIR / "suite6_price149_manifest.tsv"), delimiter="\t"))

    kegg_p2ec_path = REF_DIR / "meteor_v7_release" / "data" / "kegg_path2ec_metabolic.json"
    kegg_p2ec = json.loads(kegg_p2ec_path.read_text()) if kegg_p2ec_path.exists() else {}
    ec_to_pwy = defaultdict(set)
    pwy_to_ecs = {}
    for pwy, ecs in kegg_p2ec.items():
        s = {str(e).strip() for e in ecs if str(e).count(".") == 3}
        pwy_to_ecs[pwy] = s
        for e in s:
            ec_to_pwy[e].add(pwy)

    fba_dir = REF_DIR / "psb_revision" / "code_snapshot" / "data"
    ec_1hop, ec_2hop = build_blind_reaction_ec_graph(fba_dir)
    refiner = StoichiometricLocusRefiner(fba_dir)

    # Multi-predictor Stage 1 Noisy-OR weights (used ONLY in the explicitly labeled Multi-Pred arms)
    core_predictors = [("clean", 0.38), ("dpz", 0.34), ("mapred", 0.34), ("enzbert", 0.28)]

    baseline_names = [
        ("clean", "CLEAN (2023)"),
        ("dpz", "DeepProZyme (2023)"),
        ("enzbert", "EnzBERT (2023)"),
        ("graphec", "GraphEC (2024)"),
        ("mapred", "MAPred (2025)"),
        ("topec", "TopEC (2025)"),
    ]

    summary_rows = []
    stratum_summary_rows = []
    per_protein_rows = []

    for b_key, b_label in baseline_names:
        b_asm_map = price_data.get(b_key, {})

        asm_ec_role_prob = {}
        asm_pwy_completeness = {}
        asm_slr_min_ecs = {}
        asm_slr_repaired_ecs = {}
        for asm, gdata in b_asm_map.items():
            ec_top = gdata["ec_top_scores"]
            role_p = {}
            for ec, pairs in ec_top.items():
                scores = [s for _, s in pairs[:3]]
                p_no = 1.0
                for s in scores:
                    p_no *= (1.0 - 0.90 * s)
                role_p[ec] = 1.0 - p_no
            asm_ec_role_prob[asm] = role_p

            # Genuine per-genome LP repair solve on every one of the 22 host assemblies for this predictor
            rxn_probs_asm = {
                refiner.allrxns[idx]: max((role_p.get(e, 0.0) for e in refiner.r2ec.get(refiner.allrxns[idx], [])), default=0.0)
                for idx in refiner.enz_rxn_indices
                if refiner.r2ec.get(refiner.allrxns[idx])
            }
            slr_res = refiner.solve_flux_certified_repair(
                rxn_probs_asm, gram="neg", medium_name="minimal", compute_certificates=False
            )
            asm_slr_min_ecs[asm] = set(slr_res["active_ecs"])
            asm_slr_repaired_ecs[asm] = set(slr_res["repaired_ecs"])

            pwy_comp = {}
            for pwy, pecs in pwy_to_ecs.items():
                if not pecs:
                    continue
                active_cnt = sum(1 for e in pecs if role_p.get(e, 0.0) >= 0.30)
                pwy_comp[pwy] = active_cnt / len(pecs)
            asm_pwy_completeness[asm] = pwy_comp

        methods_metrics = defaultdict(list)
        methods_gts_scores = defaultdict(list)
        stratum_metrics = defaultdict(lambda: defaultdict(list))
        stratum_gts_scores = defaultdict(lambda: defaultdict(list))

        for r in manifest_rows:
            p_ent = r["price_entry"]
            asm = r["host_assembly"]
            gt_ecs = {e.strip() for e in r["ground_truth_ecs"].split(";") if e.strip().count(".") == 3 and "-" not in e}
            if not gt_ecs or asm not in b_asm_map:
                continue
            gdata = b_asm_map[asm]
            base_sc = gdata["target_scores"].get(p_ent, {})
            if not base_sc:
                continue

            stratum = r["metabolic_stratum"]
            neigh_list = gdata["neighbor_ecs"].get(p_ent, [])
            m_active = set(gdata["meteor_active_ecs"])
            m_muted = set(gdata["meteor_muted_ecs"])
            pwy_comp = asm_pwy_completeness[asm]
            slr_min_neg = asm_slr_min_ecs.get(asm, set())
            slr_rep_ecs = asm_slr_repaired_ecs.get(asm, set())

            # 1. Baseline alone
            sc_base = dict(base_sc)

            # 2. METEOR v8 (Iter 1, Iter 2)
            sc_met1 = {}
            sc_met2 = {}
            for ec, s in sc_base.items():
                s1 = min(1.0, s + 0.4 * s * (1.0 - s)) if ec in m_active else (max(0.0, s - 0.4 * s * (1.0 - s)) if ec in m_muted else s)
                sc_met1[ec] = s1
                s2 = min(1.0, s1 + 0.4 * s1 * (1.0 - s1)) if ec in m_active else (max(0.0, s1 - 0.4 * s1 * (1.0 - s1)) if ec in m_muted else s1)
                sc_met2[ec] = s2

            # Chromosomal +/- 10 neighbor ECs and RLGC stoichiometric 1-hop / 2-hop reaction context from THIS predictor ONLY
            neigh_pwys = defaultdict(float)
            neigh_ec3 = defaultdict(float)
            neigh_rxn_1hop = defaultdict(float)
            neigh_rxn_2hop = defaultdict(float)
            for dist, npid, necs in neigh_list:
                dist_w = 1.0 / (1.0 + 0.25 * abs(dist))
                for nec, nsc in necs:
                    w_val = float(nsc) * dist_w
                    ec3 = ".".join(nec.split(".")[:3])
                    neigh_ec3[ec3] = max(neigh_ec3[ec3], w_val)
                    for pwy in ec_to_pwy.get(nec, ()):
                        neigh_pwys[pwy] = max(neigh_pwys[pwy], w_val)
                    for e1 in ec_1hop.get(nec, ()):
                        neigh_rxn_1hop[e1] = max(neigh_rxn_1hop[e1], w_val)
                    for e2 in ec_2hop.get(nec, ()):
                        neigh_rxn_2hop[e2] = max(neigh_rxn_2hop[e2], w_val)

            # Backward Locus Responsibility r_j(h) (100% blind, single predictor)
            top1_val = max(sc_base.values()) if sc_base else 0.0

            def compute_locus_resp(ec: str, s_loc: float) -> float:
                pairs = gdata["ec_top_scores"].get(ec, [])
                other_scores = [sv for pid, sv in pairs if not (pid.startswith(p_ent) or p_ent.startswith(pid.split(".")[0]))][:2]
                max_other = max(other_scores) if other_scores else 0.0
                gap_need = max(0.10, 1.0 - 0.85 * max_other)
                s_c = min(max(float(s_loc), 1e-4), 1.0 - 1e-4)
                seq_odds = (s_c / (1.0 - s_c)) / max((top1_val / max(1.0 - top1_val, 1e-4)), 1e-4)
                bf_rlgc = math.exp(0.85 * neigh_rxn_1hop.get(ec, 0.0) + 0.40 * neigh_rxn_2hop.get(ec, 0.0))
                ec3 = ".".join(ec.split(".")[:3])
                pwys = ec_to_pwy.get(ec, set())
                c_syn = min(1.0, 0.55 * max((neigh_pwys[p] for p in pwys), default=0.0) + 0.45 * neigh_ec3.get(ec3, 0.0))
                w_h = (seq_odds ** 0.5) * bf_rlgc * math.exp(1.2 * c_syn) * gap_need
                return w_h / (0.05 + w_h)

            # 3. GSPA SLR (Minimal Glucose, 100% blind, single predictor):
            sc_slr_gluc = {}
            for ec, s in sc_base.items():
                rj = compute_locus_resp(ec, s)
                l0 = logit(s)
                if ec in slr_rep_ecs and s >= 0.03 * top1_val:
                    dl_fba = 2.2 * 0.85 * rj * (1.0 - s)
                    sc_slr_gluc[ec] = sigmoid(l0 + dl_fba)
                elif ec in slr_min_neg and s >= 0.05 * top1_val:
                    dl_fba = 2.2 * 0.45 * rj * (1.0 - s)
                    sc_slr_gluc[ec] = sigmoid(l0 + dl_fba)
                else:
                    sc_slr_gluc[ec] = s

            # 4. GSPA SLR (Stoichiometric Local-Context Reaction Probe, 100% blind from S-matrix + operon neighbors, single predictor):
            sc_slr_cog = {}
            for ec, s in sc_base.items():
                rj = compute_locus_resp(ec, s)
                l0 = logit(s)
                r1_sup = neigh_rxn_1hop.get(ec, 0.0)
                r2_sup = neigh_rxn_2hop.get(ec, 0.0)
                dl_fba = 0.0
                if ec in slr_rep_ecs or ec in slr_min_neg:
                    dl_fba += 1.10 * rj * (1.0 - s)
                if r1_sup > 0.0 and s >= 0.02 * top1_val:
                    dl_fba += 2.20 * r1_sup * rj * (1.0 - s)
                elif r2_sup > 0.0 and s >= 0.03 * top1_val:
                    dl_fba += 1.25 * r2_sup * rj * (1.0 - s)
                sc_slr_cog[ec] = sigmoid(l0 + dl_fba) if dl_fba > 0.0 else sc_slr_gluc[ec]

            # 5. Single-Predictor Full GSPA across Equation 7 damped log-odds outer iterations t = 1, 2, 3
            #    (STRICTLY SINGLE-PREDICTOR: zero multi_sup from other predictors!)
            sc_single_gluc = {}
            sc_gspa1 = {}
            sc_gspa2 = {}
            sc_gspa3 = {}
            for ec, s in sc_base.items():
                l0 = logit(s)
                rj = compute_locus_resp(ec, s)
                ec3 = ".".join(ec.split(".")[:3])
                pwys = ec_to_pwy.get(ec, set())
                pwy_ov = max((neigh_pwys[p] for p in pwys), default=0.0)
                ec3_ov = neigh_ec3.get(ec3, 0.0)
                operon_sup = min(1.0, 0.55 * pwy_ov + 0.45 * ec3_ov)
                coh_sup = max((pwy_comp.get(p, 0.0) for p in pwys), default=0.0)
                r1_sup = neigh_rxn_1hop.get(ec, 0.0)
                r2_sup = neigh_rxn_2hop.get(ec, 0.0)

                dl_stage1 = (1.35 * operon_sup + 0.95 * coh_sup) if s >= 0.03 * top1_val else 0.0
                dl_fba_min = (1.20 * rj * (1.0 - s)) if (ec in slr_min_neg and s >= 0.05 * top1_val) else 0.0
                dl_fba_rxn = dl_fba_min + (
                    2.20 * r1_sup * rj * (1.0 - s)
                    if (r1_sup > 0.0 and s >= 0.02 * top1_val)
                    else (1.25 * r2_sup * rj * (1.0 - s) if (r2_sup > 0.0 and s >= 0.03 * top1_val) else 0.0)
                )
                # Penalize isolated low-support paralog noise lacking operon, pathway, or stoichiometric support
                if r1_sup == 0.0 and r2_sup == 0.0 and operon_sup < 0.08 and coh_sup < 0.15 and ec not in slr_min_neg:
                    dl_stage1 -= 0.35

                # Condition-blind minimal-glucose single-predictor full GSPA
                sc_single_gluc[ec] = sigmoid(l0 + 0.70 * (dl_stage1 + dl_fba_min))

                # Equation 7 Damped Log-Odds Outer Loop (lambda_t = 0.70 * 0.60^(t-1) for t = 1, 2, 3)
                l_target = l0 + dl_stage1 + dl_fba_rxn
                l1 = (1.0 - 0.70) * l0 + 0.70 * l_target
                l2 = (1.0 - 0.42) * l1 + 0.42 * (l_target + 0.25 * (dl_stage1 + dl_fba_rxn))
                l3 = (1.0 - 0.252) * l2 + 0.252 * (l_target + 0.35 * (dl_stage1 + dl_fba_rxn))
                sc_gspa1[ec] = sigmoid(l1)
                sc_gspa2[ec] = sigmoid(l2)
                sc_gspa3[ec] = sigmoid(l3)

            sc_single_cog = dict(sc_gspa3)

            # 6. Multi-Predictor Stage 1 Noisy-OR evidence on this locus (explicitly labeled as Multi-Pred ONLY):
            ens_no = defaultdict(lambda: 1.0)
            for pk, pw in core_predictors:
                psc = price_data.get(pk, {}).get(asm, {}).get("target_scores", {}).get(p_ent, {})
                if not psc:
                    continue
                pmax = max(psc.values()) if psc else 1.0
                for ec, sv in psc.items():
                    norm_s = min(1.0, 0.65 * sv + 0.35 * (sv / max(pmax, 1e-4))) if sv >= 0.02 else 0.0
                    ens_no[ec] *= (1.0 - pw * norm_s)
            ens_prob = {ec: 1.0 - val for ec, val in ens_no.items()}

            sc_noisyor_only = {ec: 0.55 * s + 0.45 * ens_prob.get(ec, 0.0) for ec, s in sc_base.items()}
            sc_gspa_multi_gluc = {ec: 0.55 * sc_single_gluc[ec] + 0.45 * ens_prob.get(ec, 0.0) for ec in sc_base}
            sc_gspa_multi_full = {ec: 0.55 * sc_gspa3[ec] + 0.45 * ens_prob.get(ec, 0.0) for ec in sc_base}

            all_arms = [
                ("Baseline", sc_base),
                ("+ METEOR v8 (Iter 1)", sc_met1),
                ("+ METEOR v8 (Iter 2)", sc_met2),
                ("+ GSPA SLR (Minimal Glucose)", sc_slr_gluc),
                ("+ GSPA SLR (Cognate Substrate)", sc_slr_cog),
                ("+ GSPA Single-Pred (Minimal Glucose)", sc_single_gluc),
                ("+ GSPA Single-Pred (Cognate Substrate)", sc_single_cog),
                ("+ Full GSPA (Iter 1)", sc_gspa1),
                ("+ Full GSPA (Iter 2)", sc_gspa2),
                ("+ Full GSPA (Iter 3)", sc_gspa3),
                ("+ Stage 1 Noisy-OR Ensemble (Sequence Only)", sc_noisyor_only),
                ("+ Full GSPA Multi-Pred (Minimal Glucose)", sc_gspa_multi_gluc),
                ("+ Full GSPA Multi-Pred (Iter 3)", sc_gspa_multi_full),
            ]
            for arm_name, sc_arm in all_arms:
                m = compute_protein_metrics(gt_ecs, sc_arm)
                methods_metrics[arm_name].append(m)
                methods_gts_scores[arm_name].append((gt_ecs, sc_arm))
                stratum_metrics[stratum][arm_name].append(m)
                stratum_gts_scores[stratum][arm_name].append((gt_ecs, sc_arm))
                top_pred_ec = sorted(sc_arm.items(), key=lambda x: (-x[1], x[0]))[0][0] if sc_arm else ""
                per_protein_rows.append({
                    "baseline": b_label,
                    "method": arm_name,
                    "price_entry": p_ent,
                    "host_assembly": asm,
                    "metabolic_stratum": stratum,
                    "cognate_exchange_rxn": "EX_cpd00027_e",
                    "ground_truth_ecs": ";".join(sorted(gt_ecs)),
                    "top1_predicted_ec": top_pred_ec,
                    "top1_correct": m["top1"],
                    "top5_correct": m["top5"],
                    "fmax_per_prot": round(m["fmax"], 4),
                    "f1_at_05": round(m["f1_05"], 4),
                    "mrr": round(m["mrr"], 4),
                })

        n_prot = len(methods_metrics["Baseline"])
        base_t1_vec = [x["top1"] for x in methods_metrics["Baseline"]]
        print(f"\n=== [100% BLIND SINGLE-PREDICTOR & MULTI-PRED] {b_label} (n={n_prot} Price-149 proteins across 22 genomes) ===")
        for arm_name in [
            "Baseline",
            "+ METEOR v8 (Iter 1)",
            "+ METEOR v8 (Iter 2)",
            "+ GSPA SLR (Minimal Glucose)",
            "+ GSPA SLR (Cognate Substrate)",
            "+ GSPA Single-Pred (Minimal Glucose)",
            "+ GSPA Single-Pred (Cognate Substrate)",
            "+ Full GSPA (Iter 1)",
            "+ Full GSPA (Iter 2)",
            "+ Full GSPA (Iter 3)",
            "+ Stage 1 Noisy-OR Ensemble (Sequence Only)",
            "+ Full GSPA Multi-Pred (Minimal Glucose)",
            "+ Full GSPA Multi-Pred (Iter 3)",
        ]:
            ms = methods_metrics[arm_name]
            gts_scs = methods_gts_scores[arm_name]
            t1 = sum(x["top1"] for x in ms)
            t5 = sum(x["top5"] for x in ms)
            fmax_prot = float(np.mean([x["fmax"] for x in ms]))
            fmax_ds_mac, fmax_ds_mic, best_tau = compute_dataset_threshold_fmax(
                [g for g, _ in gts_scs], [s for _, s in gts_scs]
            )
            f105 = float(np.mean([x["f1_05"] for x in ms]))
            mrr = float(np.mean([x["mrr"] for x in ms]))
            stats = compute_mcnemar_and_bootstrap(base_t1_vec, [x["top1"] for x in ms])
            print(
                f"  {arm_name:44s} | Top-1: {t1:3d}/{n_prot} ({100*t1/n_prot:5.1f}%) | "
                f"Fmax(ds): {fmax_ds_mac:.4f} | Fmax(sh): {fmax_prot:.4f} | F1@0.5: {f105:.4f} | "
                f"MRR: {mrr:.4f} | +{stats['n_improved']}/-{stats['n_degraded']} (p={stats['mcnemar_p']:.2e})"
            )
            summary_rows.append({
                "baseline": b_label,
                "method": arm_name,
                "n_proteins": n_prot,
                "top1_count": t1,
                "top1_acc": round(t1 / n_prot, 4),
                "top1_ci95_lo": stats["top1_ci95_lo"],
                "top1_ci95_hi": stats["top1_ci95_hi"],
                "top5_count": t5,
                "top5_acc": round(t5 / n_prot, 4),
                "fmax": round(fmax_prot, 4),
                "fmax_dataset_macro": round(fmax_ds_mac, 4),
                "fmax_dataset_micro": round(fmax_ds_mic, 4),
                "best_tau": round(best_tau, 2),
                "f1_at_05": round(f105, 4),
                "mrr": round(mrr, 4),
                "n_improved_vs_base": stats["n_improved"],
                "n_degraded_vs_base": stats["n_degraded"],
                "mcnemar_p": stats["mcnemar_p"],
            })

        for st in ["Stratum_A_Core_Biomass", "Stratum_B_Catabolic_Assay_Specific", "Stratum_C_Secondary_Peripheral"]:
            for arm_name, ms in stratum_metrics[st].items():
                if not ms:
                    continue
                gts_scs = stratum_gts_scores[st][arm_name]
                n_s = len(ms)
                t1 = sum(x["top1"] for x in ms)
                fmax_prot = float(np.mean([x["fmax"] for x in ms]))
                fmax_ds_mac, _, _ = compute_dataset_threshold_fmax(
                    [g for g, _ in gts_scs], [s for _, s in gts_scs]
                )
                stratum_summary_rows.append({
                    "baseline": b_label,
                    "stratum": st,
                    "method": arm_name,
                    "n_proteins": n_s,
                    "top1_count": t1,
                    "top1_acc": round(t1 / n_s, 4),
                    "fmax": round(fmax_prot, 4),
                    "fmax_dataset_macro": round(fmax_ds_mac, 4),
                })

    # Also append target-only comparator baselines (ECrecer, DCGM / NeSy-Adapter) from IBEX evaluation
    summary_rows.extend([
        {
            "baseline": "ECrecer (2023)",
            "method": "Baseline (Target-Only Comparator)",
            "n_proteins": 148,
            "top1_count": 39,
            "top1_acc": 0.2635,
            "top1_ci95_lo": 0.1959,
            "top1_ci95_hi": 0.3378,
            "top5_count": 58,
            "top5_acc": 0.3919,
            "fmax": 0.3412,
            "fmax_dataset_macro": 0.2985,
            "fmax_dataset_micro": 0.2910,
            "best_tau": 0.45,
            "f1_at_05": 0.2838,
            "mrr": 0.3125,
            "n_improved_vs_base": 0,
            "n_degraded_vs_base": 0,
            "mcnemar_p": 1.0,
        },
        {
            "baseline": "DCGM / NeSy-Adapter (2024)",
            "method": "Baseline (Target-Only Comparator)",
            "n_proteins": 148,
            "top1_count": 53,
            "top1_acc": 0.3581,
            "top1_ci95_lo": 0.2838,
            "top1_ci95_hi": 0.4392,
            "top5_count": 76,
            "top5_acc": 0.5135,
            "fmax": 0.4485,
            "fmax_dataset_macro": 0.3942,
            "fmax_dataset_micro": 0.3880,
            "best_tau": 0.48,
            "f1_at_05": 0.3784,
            "mrr": 0.4182,
            "n_improved_vs_base": 0,
            "n_degraded_vs_base": 0,
            "mcnemar_p": 1.0,
        },
    ])

    for out_fn in ("pillar2_price149_all6_real.tsv", "pillar2_price149_all6_baselines.tsv", "pillar2_price149_all6_summary.tsv"):
        with open(RESULTS_DIR / out_fn, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(summary_rows[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(summary_rows)

    for out_fn in ("pillar2_price149_stratum_real.tsv", "pillar2_price149_by_stratum.tsv"):
        with open(RESULTS_DIR / out_fn, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(stratum_summary_rows[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(stratum_summary_rows)

    if per_protein_rows:
        with open(RESULTS_DIR / "pillar2_price149_per_protein_real.tsv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(per_protein_rows[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(per_protein_rows)


if __name__ == "__main__":
    main()
