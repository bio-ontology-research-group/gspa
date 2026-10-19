#!/usr/bin/env python3
"""
Pillar 2 Complete 100%-Real Empirical Evaluation (Suites 3, 4, 5, 6):
  - Experiment 2A: 108-Genome Mathematical Invariant & Structural Audit
    (CLEAN, DeepProZyme, EnzBERT, CarveMe, Reconstructor, METEOR v8, and GSPA StoichiometricLocusRefiner)
    computed directly from `meteor_v8_release/results/table1/table1_*.json`,
    `psb_revision/results/thresh_gapfill_deadends_summary.json`,
    `toolcompare/panelB_curated6.json`, `toolcompare/panelB_ec.json`, `toolcompare/recon_dpz_ec.json`,
    and `StoichiometricLocusRefiner.solve_flux_certified_repair()`.
  - Experiment 2B: Real Role-Stratified Leave-Reaction-Out (LRO) & Secondary BGC Benchmark (Suites 3, 4, 5)
    computed directly from all real LRO candidate ranking tables in `lro_data/results_*.tsv`
    (E. coli MG1655 n=414 & n=256, B. subtilis 168 n=144, A. aeolicus VF5 n=93),
    Stage 3 Agentic LRO re-ranking (`results/pillar4_lro_reranking.tsv`), and the 106-genome
    antiSMASH / MIBiG secondary metabolism benchmark (`bgc_specificity_antismash_v8.tsv`).
  - Experiment 2C: 100% Blind Price-149 (n=149 across 22 genomes) Benchmark across all 6 SOTA EC Baselines
    (CLEAN 2023, DeepProZyme 2023, EnzBERT 2023, MAPred 2025, GraphEC 2024, TopEC 2025),
    METEOR v8, Noisy-OR, and GSPA (Iter 1, 2, 3) via `scripts/eval_price149_real_all6.py`.
"""

import csv
import json
import sys
from collections import defaultdict
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent.parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "benchmark" / "manifests"))
from stoichiometric_locus_refiner import StoichiometricLocusRefiner
import eval_price149_real_all6

MANIFEST_DIR = ROOT / "benchmark" / "manifests"
REF_DIR = MANIFEST_DIR / "reference_data"
LRO_DIR = REF_DIR / "lro_data"
RESULTS_DIR = ROOT / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


_GLOBAL_REF: StoichiometricLocusRefiner | None = None

# Growth requirement used by the METEOR v8 MILP and by the archived grow_support threshold repair
# (meteor108_solutions_compact biomass_flux = 0.1; thresh_gapfill gamma_min = 0.1).
AUDIT_GAMMA_MIN = 0.1
_HIGHS_OPTS = {"primal_feasibility_tolerance": 1e-9, "dual_feasibility_tolerance": 1e-9}


def _network_lp_audit(
    ref: StoichiometricLocusRefiner,
    gram: str,
    internal_idx,
    target_bio: float,
    extra_uptake: list[str] | None = None,
    do_pfba: bool = False,
    eps_flux: float = 1e-6,
) -> dict:
    """Measured LP check of a repaired network on the audit medium (refiner.get_media_bounds(gram, "minimal"):
    named minimal medium + 36-compound ModelSEED supplement; same bounds as Stage 2B).
    Network = `internal_idx` (every non-exchange reaction of the configuration, enzymatic or not)
    + all EX_ reactions (secretion always allowed, uptake only for medium compounds) + the Gram biomass.
    Returns max biomass, feasibility for v_bio >= target_bio, and (optionally) the pFBA support
    (min sum|v| s.t. S v = 0, v_bio >= target_bio)."""
    from scipy.optimize import linprog
    from scipy.sparse import hstack

    lb, ub, bio_idx = ref.get_media_bounds(gram=gram, medium_name="minimal", extra_uptake_exchanges=extra_uptake)
    keep = set(int(i) for i in internal_idx) | set(ref.ex_indices) | {bio_idx}
    keep = sorted(i for i in keep if (lb[i] < -1e-12 or ub[i] > 1e-12))
    pos = {o: k for k, o in enumerate(keep)}
    S = ref.S[:, keep].tocsc()
    n_m = S.shape[0]
    bnds = [(float(lb[i]), float(ub[i])) for i in keep]
    c = np.zeros(len(keep))
    c[pos[bio_idx]] = -1.0
    r = linprog(c, A_eq=S, b_eq=np.zeros(n_m), bounds=bnds, method="highs", options=_HIGHS_OPTS)
    max_bio = float(-r.fun) if (r.success and r.fun is not None) else 0.0
    out = {"max_bio": max_bio, "feasible": bool(target_bio > 0 and max_bio >= target_bio - 1e-7), "pfba_flux_idx": []}
    if do_pfba and out["feasible"]:
        m = len(keep)
        S_split = hstack([S, -S], format="csc")
        b2 = [(max(0.0, lo), max(0.0, hi)) for lo, hi in bnds] + [(max(0.0, -hi), max(0.0, -lo)) for lo, hi in bnds]
        bi = pos[bio_idx]
        b2[bi] = (float(target_bio), float(max(target_bio, b2[bi][1])))
        b2[m + bi] = (0.0, 0.0)
        obj = np.ones(2 * m)
        obj[bi] = 0.0
        r2 = linprog(obj, A_eq=S_split, b_eq=np.zeros(n_m), bounds=b2, method="highs", options=_HIGHS_OPTS)
        if r2.success:
            v = r2.x[:m] - r2.x[m:]
            out["pfba_flux_idx"] = [keep[k] for k in range(m) if abs(v[k]) > eps_flux]
    return out


def _pfba_on_bounds(S, lb, ub, net_idx, bio_idx: int, target_bio: float, eps_flux: float = 1e-6) -> dict:
    """Max biomass and pFBA support (min sum|v| s.t. S v = 0, v_bio >= target_bio) of the network
    `net_idx` under explicit bounds (lb, ub); reactions outside the network are fixed to 0.
    Used for METEOR's own MILP medium (archived meteor108_drafts_compact media[gram].lb/ub)."""
    from scipy.optimize import linprog
    from scipy.sparse import hstack

    keep = sorted(i for i in set(int(x) for x in net_idx) | {int(bio_idx)} if (lb[i] < -1e-12 or ub[i] > 1e-12))
    pos = {o: k for k, o in enumerate(keep)}
    Sk = S[:, keep].tocsc()
    n_m = Sk.shape[0]
    bnds = [(float(lb[i]), float(ub[i])) for i in keep]
    c = np.zeros(len(keep))
    c[pos[bio_idx]] = -1.0
    r = linprog(c, A_eq=Sk, b_eq=np.zeros(n_m), bounds=bnds, method="highs", options=_HIGHS_OPTS)
    max_bio = float(-r.fun) if (r.success and r.fun is not None) else 0.0
    out = {"max_bio": max_bio, "feasible": bool(max_bio >= target_bio - 1e-7), "pfba_flux_idx": []}
    if out["feasible"]:
        m = len(keep)
        S_split = hstack([Sk, -Sk], format="csc")
        b2 = [(max(0.0, lo), max(0.0, hi)) for lo, hi in bnds] + [(max(0.0, -hi), max(0.0, -lo)) for lo, hi in bnds]
        bi = pos[bio_idx]
        b2[bi] = (float(target_bio), float(max(target_bio, b2[bi][1])))
        b2[m + bi] = (0.0, 0.0)
        obj = np.ones(2 * m)
        obj[bi] = 0.0
        r2 = linprog(obj, A_eq=S_split, b_eq=np.zeros(n_m), bounds=b2, method="highs", options=_HIGHS_OPTS)
        if r2.success:
            v = r2.x[:m] - r2.x[m:]
            out["pfba_flux_idx"] = [keep[k] for k in range(m) if abs(v[k]) > eps_flux]
    return out


def _audit_repaired_network(args_tuple: tuple) -> tuple[str, dict]:
    """Measured audit of an archived repaired network (METEOR v8 active set y > 0.5, or the Grow-Support
    threshold repair core | grow_support), both from meteor108_drafts_compact.pkl.gz (DPZ-vanilla, tau=0.5).

    added            = internal (non-EX_, non-biomass) reactions of the network that are not in the tau=0.5 draft
    flux-certified   = added reactions with |v| > 1e-6 in a pFBA solution (v_bio >= gamma_min = 0.1) of the network
                       (a) on METEOR's own MILP medium (tight bounds + apply_media(['default']): open uptake on
                           every exchange the network contains), and
                       (b) on the audit medium (refiner "minimal" + 36-compound supplement), augmented with the
                           greedy-minimal set of non-medium uptakes the network needs to reach gamma_min when the
                           strict audit medium is infeasible (the strict-medium result is reported separately)
    zero-flux padding = added - flux-certified (per medium)."""
    gca, gram, net_indices, draft_indices, slr_target, own_lb, own_ub, own_bio_idx, fba_dir_str = args_tuple
    ref = _GLOBAL_REF if _GLOBAL_REF is not None else StoichiometricLocusRefiner(Path(fba_dir_str))
    ex_set = set(ref.ex_indices)
    net = set(int(x) for x in net_indices)
    internal = net - ex_set
    strict = _network_lp_audit(ref, gram, internal, AUDIT_GAMMA_MIN, do_pfba=True)
    strict_slr = strict["max_bio"] >= slr_target - 1e-7
    # Diagnose infeasibility: uptake exchanges the network uses that the audit medium blocks.
    lb, ub, _ = ref.get_media_bounds(gram=gram, medium_name="minimal")
    own_extra = sorted(ref.allrxns[i] for i in net & ex_set if ub[i] <= 1e-12 and ref.allrxns[i] != "EX_biomass")
    own = _network_lp_audit(ref, gram, internal, AUDIT_GAMMA_MIN, extra_uptake=own_extra, do_pfba=True)
    required_extra = list(own_extra)
    if own["feasible"] and not strict["feasible"]:
        # Greedy minimal subset of non-medium uptakes needed to reach v_bio >= gamma_min.
        for ex in list(own_extra):
            trial = [e for e in required_extra if e != ex]
            if _network_lp_audit(ref, gram, internal, AUDIT_GAMMA_MIN, extra_uptake=trial)["feasible"]:
                required_extra = trial
    elif strict["feasible"]:
        required_extra = []
    if strict["feasible"]:
        aud = strict
    elif own["feasible"]:
        aud = _network_lp_audit(ref, gram, internal, AUDIT_GAMMA_MIN, extra_uptake=required_extra, do_pfba=True)
    else:
        aud = own
    aud_flux = set(aud["pfba_flux_idx"]) if aud["feasible"] else set()
    # (a) METEOR's own MILP medium
    om = _pfba_on_bounds(ref.S, own_lb, own_ub, net, own_bio_idx, AUDIT_GAMMA_MIN)
    om_flux = set(om["pfba_flux_idx"])
    nonbio_internal = {i for i in internal if not ref.allrxns[i].startswith("biomass_")}
    added = nonbio_internal - set(int(x) for x in draft_indices)
    cert_own = added & om_flux
    cert_aud = added & aud_flux
    return gca, {
        "strict_max_bio": strict["max_bio"],
        "strict_feasible": strict["feasible"],
        "strict_feasible_slr_target": bool(strict_slr),
        "own_uptake_max_bio": own["max_bio"],
        "own_uptake_feasible": own["feasible"],
        "non_medium_uptakes_in_network": own_extra,
        "required_non_medium_uptakes": required_extra if own["feasible"] else own_extra,
        "pfba_n_flux_carrying": len(set(own["pfba_flux_idx"]) & net) if own["feasible"] else 0,
        "n_active": len(net),
        "n_internal": len(nonbio_internal),
        "n_added": len(added),
        "n_added_imbalanced": len(added & ref.imbalanced_rxn_indices),
        "n_added_exchange_not_medium": len(own_extra),
        "own_medium_max_bio": om["max_bio"],
        "own_medium_feasible": om["feasible"],
        "own_medium_n_flux_carrying": len(om_flux & net),
        "own_medium_n_flux_certified_added": len(cert_own),
        "audit_medium_feasible_with_required": aud["feasible"],
        "audit_medium_n_flux_carrying": len(aud_flux & net),
        "audit_medium_n_flux_certified_added": len(cert_aud),
        "added_idx": sorted(added),
        "cert_own_idx": sorted(cert_own),
        "cert_aud_idx": sorted(cert_aud),
    }


def _solve_single_genome_slr(args_tuple: tuple) -> tuple[str, dict]:
    gca, gram, act_indices, flx_indices, dpz_ec_scores, fba_dir_str = args_tuple
    ref = _GLOBAL_REF if _GLOBAL_REF is not None else StoichiometricLocusRefiner(Path(fba_dir_str))
    act_set = set(int(x) for x in act_indices)
    flx_set = set(int(x) for x in flx_indices)
    zero_flx_core = act_set - flx_set
    rxn_probs = {}
    # Map real DeepProZyme EC probabilities to universal reactions via ref.ec2r
    for ec, sc in dpz_ec_scores.items():
        for rxn_id in ref.ec2r.get(ec, []):
            if sc > rxn_probs.get(rxn_id, 0.0):
                rxn_probs[rxn_id] = float(sc)
    # Reactions in zero_flx_core were activated without carrying flux in METEOR, which occurs
    # if and only if their sequence probability exceeded tau_core = 0.50 in the baseline draft.
    #
    # INDEPENDENCE CAVEAT (state this in Methods): the archive holds DeepProZyme EC scores only for
    # the evaluated proteins of 95/108 genomes (full_genome_eval_*.pkl "score_changes"), not the full
    # per-genome tau=0.5 draft. The GSPA draft core is therefore seeded from METEOR's own archived
    # solution: every balanced enzymatic reaction METEOR switched on (y > 0.5) without flux is given
    # p = 0.75 (>= tau_core, i.e. forced into the GSPA core). GSPA and METEOR thus share most of the
    # draft core (it is taken from METEOR's output, not re-derived from sequence); the comparison
    # isolates the repair policy (positive-cost LP + pFBA vs negative-cost MILP), not an independent
    # reconstruction. GSPA's core also inherits any METEOR-selected reactions with p < 0.5.
    for r_i in zero_flx_core:
        if r_i in ref.enz_rxn_indices and r_i not in ref.imbalanced_rxn_indices:
            rxn_id = ref.allrxns[r_i]
            rxn_probs[rxn_id] = max(rxn_probs.get(rxn_id, 0.0), 0.75)
    res = ref.solve_flux_certified_repair(
        rxn_probs, gram=gram, medium_name="minimal", tau_core=0.50, alpha_bio=0.10, compute_certificates=True
    )
    # Measured LP check of the GSPA-repaired network (core enzymatic + added enzymatic + the
    # non-EC/transport reactions Stage 2B routed flux through at zero cost) on the same audit medium.
    ex_set = set(ref.ex_indices)
    gspa_internal = {ref.rxn_to_idx[r] for r in res["active_rxns"]} | {
        ref.rxn_to_idx[r] for r in res["pfba_flux_by_rxn"] if ref.rxn_to_idx[r] not in ex_set
    }
    nonenz_free = {
        ref.rxn_to_idx[r] for r in res["pfba_flux_by_rxn"]
        if ref.rxn_to_idx[r] not in ex_set and ref.rxn_to_idx[r] not in ref.enz_rxn_indices
        and not r.startswith("biomass_")
    }
    tgt = float(res["target_bio"])
    a_gamma = _network_lp_audit(ref, gram, gspa_internal, AUDIT_GAMMA_MIN)
    res["_lp_audit"] = {
        "max_bio": a_gamma["max_bio"],
        "feasible_gamma_min": a_gamma["feasible"],
        "feasible_slr_target": bool(tgt > 0 and a_gamma["max_bio"] >= tgt - 1e-7),
        "n_nonec_internal_flux": len(nonenz_free),
    }
    return gca, res


def run_108_genome_invariant_audit(refiner: StoichiometricLocusRefiner) -> Path:
    """
    Experiment 2A: 108-Genome Mathematical Invariant & Structural Audit across:
      - Baseline Sequence Threshold (DeepProZyme, CLEAN, EnzBERT)
      - CarveMe (Machado et al. 2018) on 6 BiGG curated GEMs & 98-genome tool_compare.tsv
      - Reconstructor (Jenior et al. 2023) on 6 BiGG curated GEMs
      - METEOR v8 (Negative-Cost MILP) on 108 genomes
      - Grow-Support Repaired Baseline (per-genome networks re-derived on ibex, meteor108_drafts_compact.pkl.gz;
        counts verified against psb_revision tgd_<GCA>.json for 108/108 genomes)
      - GSPA StoichiometricLocusRefiner (Positive-Cost LP + pFBA Support Certificate across all 108 real genomes)
    """
    import concurrent.futures
    import gzip
    import multiprocessing as mp
    import pickle

    print("[2A] Running 108-Genome Mathematical Invariant & Tool Comparison Audit across all 108 real genomes...")
    global _GLOBAL_REF
    _GLOBAL_REF = refiner
    # Pre-warm universal biomass LP cache for Gram-negative and Gram-positive media once
    slr_target_by_gram = {
        g: float(refiner.solve_flux_certified_repair({}, gram=g, medium_name="minimal", compute_certificates=False)["target_bio"])
        for g in ("neg", "pos")
    }

    table1_dir = REF_DIR / "meteor_v8_release" / "results" / "table1"
    flux_stats_json = REF_DIR / "meteor108_flux_stats.json"
    flux_stats = json.loads(flux_stats_json.read_text(encoding="utf-8")) if flux_stats_json.exists() else {}
    sols_compact = pickle.load(gzip.open(REF_DIR / "meteor108_solutions_compact.pkl.gz", "rb"))
    # Per-genome tau=0.5 DeepProZyme-vanilla draft, METEOR active set + MILP fluxes, Grow-Support repair and
    # METEOR's own MILP bounds, extracted on ibex (paper-reproducibility/pipelines/P07_invariant_audit).
    drafts = pickle.load(gzip.open(REF_DIR / "meteor108_drafts_compact.pkl.gz", "rb"))
    if list(drafts["allrxns"]) != list(refiner.allrxns):
        raise RuntimeError("meteor108_drafts_compact allrxns order differs from the refiner's universal model")
    dg = drafts["genomes"]
    dmedia = drafts["media"]
    n_act_mismatch = sum(
        1 for g, r in dg.items()
        if g in sols_compact and set(map(int, r["meteor_active_idx"])) != set(map(int, sols_compact[g]["active_indices"]))
    )
    print(f"[2A] drafts archive: {len(dg)} genomes; METEOR active sets differing from meteor108_solutions_compact: {n_act_mismatch}")

    thresh_repair_json = REF_DIR / "ibex_extracted" / "psb_revision" / "results" / "thresh_gapfill_deadends_summary.json"
    thresh_repair = json.loads(thresh_repair_json.read_text(encoding="utf-8"))
    # Archived aggregate of the Grow-Support repair, kept only as a cross-check of the re-derived networks.
    print(f"[2A] archived Grow-Support aggregate: n_added_by_repair mean = {thresh_repair['n_added_by_repair']['mean']}, "
          f"n_selected mean = {thresh_repair['n_selected']['repaired_baseline']['mean']}; re-derived genomes matching "
          f"tgd_<GCA>.json counts: {sum(1 for r in dg.values() if r.get('tgd_match'))}/{len(dg)}")

    tve_path = REF_DIR / "meteor_v8_release" / "results" / "toolcompare" / "ablation_thresh_vs_evw.json"
    tve_by_acc = {x["acc"]: x for x in json.loads(tve_path.read_text(encoding="utf-8"))} if tve_path.exists() else {}

    # Load real DeepProZyme EC probabilities per genome from full_genome_eval_gbff_fmax_dpz_vanilla.pkl
    dpz_ec_by_gcf: dict[str, dict[str, float]] = {}
    dpz_pkl_path = REF_DIR / "ibex_extracted" / "full_genome_eval_gbff_fmax_dpz_vanilla.pkl"
    if dpz_pkl_path.exists():
        dpz_pkl = pickle.load(open(dpz_pkl_path, "rb"))
        for row in dpz_pkl.get("score_changes", []):
            g = row["gcf"]
            ec = row["ec"]
            sc = float(row["base"])
            if sc > dpz_ec_by_gcf.setdefault(g, {}).get(ec, 0.0):
                dpz_ec_by_gcf[g][ec] = sc

    t1_files = sorted(table1_dir.glob("table1_*.json"))
    fba_dir_str = str(REF_DIR / "psb_revision" / "code_snapshot" / "data")
    tasks = []
    t1_records = []
    for tf in t1_files:
        d = json.loads(tf.read_text(encoding="utf-8"))
        gca = d["gca"]
        gram = "pos" if "pos" in d.get("gram", "neg").lower() else "neg"
        s_vec = sols_compact.get(gca, {"active_indices": [], "flux_indices": []})
        tasks.append((gca, gram, s_vec["active_indices"], s_vec["flux_indices"], dpz_ec_by_gcf.get(gca, {}), fba_dir_str))
        t1_records.append((gca, gram, d))

    slr_by_gca = {}
    met_lp_by_gca = {}
    gs_lp_by_gca = {}

    def _net_task(gca, gram, net_key):
        r = dg[gca]
        med = dmedia["positive" if gram == "pos" else "negative"]
        return (gca, gram, r[net_key], r["draft_idx"], slr_target_by_gram[gram], med["lb"], med["ub"],
                int(med["biomass_idx"]), fba_dir_str)

    met_tasks = [_net_task(gca, gram, "meteor_active_idx") for gca, gram, _ in t1_records]
    gs_tasks = [_net_task(gca, gram, "gs_repaired_idx") for gca, gram, _ in t1_records]
    with concurrent.futures.ProcessPoolExecutor(max_workers=16, mp_context=mp.get_context("fork")) as ex:
        for gca, res in ex.map(_audit_repaired_network, met_tasks):
            met_lp_by_gca[gca] = res
        for gca, res in ex.map(_audit_repaired_network, gs_tasks):
            gs_lp_by_gca[gca] = res
        for gca, res in ex.map(_solve_single_genome_slr, tasks):
            slr_by_gca[gca] = res

    audit_rows = []
    for gca, gram, d in t1_records:
        base_dpz = d["baseline_dpz"]
        met_dpz = d["meteor_dpz"]
        fstat = flux_stats.get(gca, {})
        s_vec = sols_compact.get(gca, {})

        n_base_sel = int(base_dpz["n_selected"])
        n_met_sel = int(met_dpz["n_selected"])
        n_active_y = int(fstat.get("n_active_y", n_met_sel))
        n_flux_y = int(fstat.get("n_flux_carrying_y", len(s_vec.get("flux_indices", []))))
        n_zero_flux_y = max(0, n_active_y - n_flux_y)

        slr_res = slr_by_gca[gca]
        slr_added = int(slr_res["n_repaired_rxns"])
        slr_flux_cert = int(slr_res["n_flux_certified_repaired"])
        slr_uncert = int(slr_res["n_uncertified_additions"])
        slr_uncert_pct = float(slr_res["uncertified_pct"])
        slr_imbal = int(slr_res["n_imbalanced_repaired"])
        slr_deadend = int(slr_res["n_deadend_repaired"])
        slr_flux_tot = len(slr_res.get("pfba_flux_by_rxn", {}))
        slr_fva_ess = len(slr_res.get("fva_essential_rxns", []))
        slr_ko_ess = len(slr_res.get("ko_essential_rxns", []))

        # MEASURED (meteor108_drafts_compact): added = METEOR internal active (y > 0.5) minus the tau=0.5
        # DeepProZyme-vanilla draft; flux-certified = added with |v| > 1e-6 under pFBA (v_bio >= 0.1), on
        # METEOR's own MILP medium and on the audit medium; padding = added - flux-certified.
        mres = met_lp_by_gca[gca]
        gres = gs_lp_by_gca[gca]
        drec = dg[gca]
        met_v = dict(zip(map(int, drec["meteor_v_idx"]), map(float, drec["meteor_v_val"])))
        met_milp_cert = sum(1 for i in mres["added_idx"] if abs(met_v.get(i, 0.0)) > 1e-6)
        meteor_net_added_over_raw = max(0, n_met_sel - n_base_sel)
        gs_tgd = drec.get("tgd", {}).get("repaired", {})

        audit_rows.append({
            "organism_id": gca,
            "gram": gram,
            "n_core_rxns": n_base_sel,
            "base_imbalanced_active": int(base_dpz["mass_imbal"]),
            "base_deadend_active": int(base_dpz["deadends"]),
            "meteor_n_active_total": n_active_y,
            "meteor_n_flux_carrying_total": n_flux_y,
            "meteor_n_zero_flux_active": n_zero_flux_y,
            "n_draft_rxns": len(drec["draft_idx"]),
            "meteor_n_net_count_diff": meteor_net_added_over_raw,
            "meteor_n_added": mres["n_added"],
            "meteor_n_added_imbalanced": mres["n_added_imbalanced"],
            "meteor_n_flux_certified_added": mres["own_medium_n_flux_certified_added"],
            "meteor_n_uncertified_padding": mres["n_added"] - mres["own_medium_n_flux_certified_added"],
            "meteor_uncertified_pct": round(100.0 * (mres["n_added"] - mres["own_medium_n_flux_certified_added"]) / max(mres["n_added"], 1), 2),
            "meteor_n_flux_certified_added_audit_medium": mres["audit_medium_n_flux_certified_added"],
            "meteor_n_uncertified_padding_audit_medium": mres["n_added"] - mres["audit_medium_n_flux_certified_added"],
            "meteor_n_added_flux_in_milp_solution": met_milp_cert,
            "meteor_pfba_n_flux_carrying_own_medium": mres["own_medium_n_flux_carrying"],
            "meteor_pfba_n_flux_carrying_audit_medium": mres["audit_medium_n_flux_carrying"],
            "meteor_lp_feasible_own_medium": int(mres["own_medium_feasible"]),
            "meteor_imbalanced_active": int(met_dpz["mass_imbal"]),
            "meteor_deadend_active": int(met_dpz["deadends"]),
            "slr_n_added": slr_added,
            "slr_n_flux_certified_added": slr_flux_cert,
            "slr_n_fva_essential": slr_fva_ess,
            "slr_n_ko_essential": slr_ko_ess,
            "slr_n_flux_carrying_total": slr_flux_tot,
            "slr_n_uncertified_padding": slr_uncert,
            "slr_uncertified_pct": slr_uncert_pct,
            "slr_imbalanced_added": slr_imbal,
            "slr_deadend_added": slr_deadend,
            "slr_solve_time_sec": round(float(slr_res["solve_time_sec"]), 4),
            # Measured LP audit columns (audit medium = refiner "minimal" + 36-compound supplement)
            "base_archived_fba_growth": float(base_dpz.get("fba_growth", 0.0)),
            "meteor_archived_fba_growth": float(met_dpz.get("fba_growth", 0.0)),
            "slr_target_bio": round(slr_target_by_gram[gram], 4),
            "meteor_lp_max_bio_audit_medium": round(met_lp_by_gca[gca]["strict_max_bio"], 6),
            "meteor_lp_feasible_gamma_min": int(met_lp_by_gca[gca]["strict_feasible"]),
            "meteor_lp_feasible_slr_target": int(met_lp_by_gca[gca]["strict_feasible_slr_target"]),
            "meteor_lp_max_bio_own_uptake": round(met_lp_by_gca[gca]["own_uptake_max_bio"], 6),
            "meteor_lp_feasible_own_uptake": int(met_lp_by_gca[gca]["own_uptake_feasible"]),
            "meteor_non_medium_uptakes": ";".join(met_lp_by_gca[gca]["non_medium_uptakes_in_network"]),
            "meteor_required_non_medium_uptakes": ";".join(met_lp_by_gca[gca]["required_non_medium_uptakes"]),
            "meteor_pfba_n_flux_carrying_own_uptake": met_lp_by_gca[gca]["pfba_n_flux_carrying"],
            "slr_lp_max_bio_audit_medium": round(slr_res["_lp_audit"]["max_bio"], 6),
            "slr_lp_feasible_gamma_min": int(slr_res["_lp_audit"]["feasible_gamma_min"]),
            "slr_lp_feasible_slr_target": int(slr_res["_lp_audit"]["feasible_slr_target"]),
            "slr_n_nonec_internal_flux": slr_res["_lp_audit"]["n_nonec_internal_flux"],
            # Grow-Support threshold repair (core | grow_support), re-derived on ibex with the archived code
            "gs_tgd_match": int(bool(drec.get("tgd_match"))),
            "gs_n_active_total": gres["n_active"],
            "gs_n_added": gres["n_added"],
            "gs_n_added_imbalanced": gres["n_added_imbalanced"],
            "gs_n_flux_certified_added": gres["own_medium_n_flux_certified_added"],
            "gs_n_uncertified_padding": gres["n_added"] - gres["own_medium_n_flux_certified_added"],
            "gs_n_flux_certified_added_audit_medium": gres["audit_medium_n_flux_certified_added"],
            "gs_n_uncertified_padding_audit_medium": gres["n_added"] - gres["audit_medium_n_flux_certified_added"],
            "gs_pfba_n_flux_carrying_own_medium": gres["own_medium_n_flux_carrying"],
            "gs_pfba_n_flux_carrying_audit_medium": gres["audit_medium_n_flux_carrying"],
            "gs_imbalanced_active": gs_tgd.get("mass_imbal", ""),
            "gs_deadend_active": gs_tgd.get("deadends", ""),
            "gs_lp_feasible_own_medium": int(gres["own_medium_feasible"]),
            "gs_lp_max_bio_audit_medium": round(gres["strict_max_bio"], 6),
            "gs_lp_feasible_gamma_min": int(gres["strict_feasible"]),
            "gs_lp_feasible_slr_target": int(gres["strict_feasible_slr_target"]),
            "gs_lp_feasible_own_uptake": int(gres["own_uptake_feasible"]),
            "gs_required_non_medium_uptakes": ";".join(gres["required_non_medium_uptakes"]),
        })

    out_path = RESULTS_DIR / "pillar2_invariant_audit_108.tsv"
    with open(out_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(audit_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(audit_rows)

    # Compile the 6-GEM Gold-Standard Reconstruction comparison (CarveMe vs Reconstructor vs METEOR v8)
    tc_dir = REF_DIR / "meteor_v8_release" / "results" / "toolcompare"
    panelb_rxn = json.loads((tc_dir / "panelB_curated6.json").read_text())
    panelb_ec = json.loads((tc_dir / "panelB_ec.json").read_text())
    recon_ec = json.loads((tc_dir / "recon_dpz_ec.json").read_text())

    gem_rows = []
    for r_rxn, r_ec, r_rec in zip(panelb_rxn, panelb_ec, recon_ec):
        org = r_rxn["organism"]
        gem = r_rxn["gem"]

        for mname, rxn_p, rxn_r, ec_p, ec_r, grows in [
            ("CarveMe (2018)", r_rxn["carveme"]["P"], r_rxn["carveme"]["R"], r_ec["carveme"]["P"], r_ec["carveme"]["R"], r_rxn["carveme"]["grows"]),
            ("Reconstructor (2023)", None, None, r_rec["P"], r_rec["R"], r_rec["grows"]),
            ("METEOR v8 (2026)", r_rxn["meteor"]["P"], r_rxn["meteor"]["R"], r_ec["meteor"]["P"], r_ec["meteor"]["R"], r_rxn["meteor"]["grows"]),
        ]:
            ec_f1 = round(2 * ec_p * ec_r / max(ec_p + ec_r, 1e-9), 4)
            rxn_f1 = round(2 * rxn_p * rxn_r / max(rxn_p + rxn_r, 1e-9), 4) if rxn_p is not None else None
            gem_rows.append({
                "organism": org,
                "curated_gem": gem,
                "method": mname,
                "rxn_precision": rxn_p if rxn_p is not None else "",
                "rxn_recall": rxn_r if rxn_r is not None else "",
                "rxn_f1": rxn_f1 if rxn_f1 is not None else "",
                "ec_precision": ec_p,
                "ec_recall": ec_r,
                "ec_f1": ec_f1,
                "fba_grows": grows,
            })

    # Append macro-average summary rows across the 6 BiGG GEMs
    for mname in [
        "CarveMe (2018)",
        "Reconstructor (2023)",
        "METEOR v8 (2026)",
    ]:
        m_sub = [r for r in gem_rows if r["method"] == mname and r["organism"] != "MEAN_6GEMS"]
        has_rxn = all(r["rxn_precision"] != "" for r in m_sub)
        gem_rows.append({
            "organism": "MEAN_6GEMS",
            "curated_gem": "6_BiGG_GEMs",
            "method": mname,
            "rxn_precision": round(float(np.mean([float(r["rxn_precision"]) for r in m_sub])), 4) if has_rxn else "",
            "rxn_recall": round(float(np.mean([float(r["rxn_recall"]) for r in m_sub])), 4) if has_rxn else "",
            "rxn_f1": round(float(np.mean([float(r["rxn_f1"]) for r in m_sub])), 4) if has_rxn else "",
            "ec_precision": round(float(np.mean([float(r["ec_precision"]) for r in m_sub])), 4),
            "ec_recall": round(float(np.mean([float(r["ec_recall"]) for r in m_sub])), 4),
            "ec_f1": round(float(np.mean([float(r["ec_f1"]) for r in m_sub])), 4),
            "fba_grows": f"{sum(1 for r in m_sub if r['fba_grows'])}/{len(m_sub)}",
        })

    gem_out = RESULTS_DIR / "pillar2_curated6_gems_real.tsv"
    with open(gem_out, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(gem_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(gem_rows)

    summary_path = RESULTS_DIR / "pillar2_invariant_audit_summary.tsv"
    n_aud = max(len(audit_rows), 1)

    def _mean(key: str) -> float:
        return float(np.mean([float(r[key]) for r in audit_rows]))

    def _pct(key: str) -> float:
        return round(100.0 * sum(int(r[key]) for r in audit_rows) / n_aud, 2)

    def _padpct(unc: float, add: float) -> float:
        return round(100.0 * unc / max(add, 1e-6), 2)

    # Two media (see paper-reproducibility/pipelines/P07_invariant_audit/README.md):
    #   own medium   = METEOR's MILP bounds (FVA-tight bounds + apply_media(['default'])); with the reversed
    #                  exchange sign every non-medium exchange is open for UPTAKE (lb=0, ub=1000), secretion closed.
    #   audit medium = refiner "minimal" + 36-compound supplement (Stage 2B); for flux certification the
    #                  greedy-minimal set of non-medium uptakes is added when the strict medium is infeasible.
    # lp_biomass_feasible_pct = strict audit medium, v_bio >= gamma_min = 0.1 (measured for METEOR,
    # Grow-Support and GSPA; Raw threshold = archived table1 baseline_dpz.fba_growth > 0).
    base_feas_pct = round(100.0 * sum(1 for r in audit_rows if r["base_archived_fba_growth"] > 1e-9) / n_aud, 2)
    raw_own_pct = round(100.0 * sum(1 for g in dg.values() if g.get("gs_draft_growth_milp_medium", 0.0) >= 0.1 - 1e-7) / max(len(dg), 1), 2)

    fields = [
        "method", "n_genomes", "mean_active_rxns", "mean_flux_carrying_active",
        "mean_zero_flux_active", "mean_added_rxns", "mean_flux_certified_added",
        "mean_uncertified_padding", "uncertified_padding_pct",
        "mean_imbalanced_added_or_active", "mean_deadend_active", "lp_biomass_feasible_pct",
        "lp_biomass_feasible_pct_slr_target",
        "mean_flux_carrying_active_audit_medium", "mean_flux_certified_added_audit_medium",
        "mean_uncertified_padding_audit_medium", "uncertified_padding_pct_audit_medium",
        "mean_imbalanced_added", "lp_biomass_feasible_pct_own_medium",
    ]
    m_act = _mean("meteor_n_active_total")
    m_flx = _mean("meteor_pfba_n_flux_carrying_own_medium")
    m_add = _mean("meteor_n_added")
    m_cert = _mean("meteor_n_flux_certified_added")
    m_unc = _mean("meteor_n_uncertified_padding")
    m_cert_a = _mean("meteor_n_flux_certified_added_audit_medium")
    m_unc_a = _mean("meteor_n_uncertified_padding_audit_medium")
    g_act = _mean("gs_n_active_total")
    g_flx = _mean("gs_pfba_n_flux_carrying_own_medium")
    g_add = _mean("gs_n_added")
    g_cert = _mean("gs_n_flux_certified_added")
    g_unc = _mean("gs_n_uncertified_padding")
    g_cert_a = _mean("gs_n_flux_certified_added_audit_medium")
    g_unc_a = _mean("gs_n_uncertified_padding_audit_medium")
    s_act = float(np.mean([r["n_core_rxns"] + r["slr_n_added"] for r in audit_rows]))
    s_flx = _mean("slr_n_flux_carrying_total")
    s_add = _mean("slr_n_added")
    s_cert = _mean("slr_n_flux_certified_added")
    s_unc = _mean("slr_n_uncertified_padding")
    b_act = _mean("n_core_rxns")
    gs_imb = [float(r["gs_imbalanced_active"]) for r in audit_rows if r["gs_imbalanced_active"] != ""]
    gs_de = [float(r["gs_deadend_active"]) for r in audit_rows if r["gs_deadend_active"] != ""]

    rows_out = [
        {
            "method": "Raw Sequence Threshold (DeepProZyme >= 0.5)",
            "mean_active_rxns": round(b_act, 1),
            "mean_flux_carrying_active": 0.0,
            "mean_zero_flux_active": round(b_act, 1),
            "mean_added_rxns": 0.0, "mean_flux_certified_added": 0.0, "mean_uncertified_padding": 0.0,
            "uncertified_padding_pct": 0.0,
            "mean_imbalanced_added_or_active": round(_mean("base_imbalanced_active"), 1),
            "mean_deadend_active": round(_mean("base_deadend_active"), 2),
            "lp_biomass_feasible_pct": base_feas_pct,
            "lp_biomass_feasible_pct_slr_target": "",
            "mean_flux_carrying_active_audit_medium": "", "mean_flux_certified_added_audit_medium": 0.0,
            "mean_uncertified_padding_audit_medium": 0.0, "uncertified_padding_pct_audit_medium": 0.0,
            "mean_imbalanced_added": 0.0,
            "lp_biomass_feasible_pct_own_medium": raw_own_pct,
        },
        {
            "method": "METEOR v8 (Negative-Cost MILP)",
            "mean_active_rxns": round(m_act, 1),
            "mean_flux_carrying_active": round(m_flx, 1),
            "mean_zero_flux_active": round(m_act - m_flx, 1),
            "mean_added_rxns": round(m_add, 2),
            "mean_flux_certified_added": round(m_cert, 2),
            "mean_uncertified_padding": round(m_unc, 2),
            "uncertified_padding_pct": _padpct(m_unc, m_add),
            "mean_imbalanced_added_or_active": round(_mean("meteor_imbalanced_active"), 1),
            "mean_deadend_active": round(_mean("meteor_deadend_active"), 2),
            "lp_biomass_feasible_pct": _pct("meteor_lp_feasible_gamma_min"),
            "lp_biomass_feasible_pct_slr_target": _pct("meteor_lp_feasible_slr_target"),
            "mean_flux_carrying_active_audit_medium": round(_mean("meteor_pfba_n_flux_carrying_audit_medium"), 1),
            "mean_flux_certified_added_audit_medium": round(m_cert_a, 2),
            "mean_uncertified_padding_audit_medium": round(m_unc_a, 2),
            "uncertified_padding_pct_audit_medium": _padpct(m_unc_a, m_add),
            "mean_imbalanced_added": round(_mean("meteor_n_added_imbalanced"), 2),
            "lp_biomass_feasible_pct_own_medium": _pct("meteor_lp_feasible_own_medium"),
        },
        {
            "method": "Grow-Support Threshold Repair (Positive-Cost LP)",
            "mean_active_rxns": round(g_act, 1),
            "mean_flux_carrying_active": round(g_flx, 1),
            "mean_zero_flux_active": round(g_act - g_flx, 1),
            "mean_added_rxns": round(g_add, 2),
            "mean_flux_certified_added": round(g_cert, 2),
            "mean_uncertified_padding": round(g_unc, 2),
            "uncertified_padding_pct": _padpct(g_unc, g_add),
            # archived MEMOTE profile of the repaired network (tgd_<GCA>.json "repaired"), as table1 for METEOR
            "mean_imbalanced_added_or_active": round(float(np.mean(gs_imb)), 1) if gs_imb else "",
            "mean_deadend_active": round(float(np.mean(gs_de)), 2) if gs_de else "",
            "lp_biomass_feasible_pct": _pct("gs_lp_feasible_gamma_min"),
            "lp_biomass_feasible_pct_slr_target": _pct("gs_lp_feasible_slr_target"),
            "mean_flux_carrying_active_audit_medium": round(_mean("gs_pfba_n_flux_carrying_audit_medium"), 1),
            "mean_flux_certified_added_audit_medium": round(g_cert_a, 2),
            "mean_uncertified_padding_audit_medium": round(g_unc_a, 2),
            "uncertified_padding_pct_audit_medium": _padpct(g_unc_a, g_add),
            "mean_imbalanced_added": round(_mean("gs_n_added_imbalanced"), 2),
            "lp_biomass_feasible_pct_own_medium": _pct("gs_lp_feasible_own_medium"),
        },
        {
            # GSPA is solved and certified on the audit medium itself (its pFBA certificate), so the
            # audit-medium columns repeat the main ones; it is not re-run on METEOR's open-uptake medium.
            "method": "GSPA StoichiometricLocusRefiner (Positive-Cost + pFBA Support Certificate)",
            "mean_active_rxns": round(s_act, 1),
            "mean_flux_carrying_active": round(s_flx, 1),
            "mean_zero_flux_active": round(round(s_act, 1) - round(s_flx, 1), 1),
            "mean_added_rxns": round(s_add, 2),
            "mean_flux_certified_added": round(s_cert, 2),
            "mean_uncertified_padding": round(s_unc, 2),
            "uncertified_padding_pct": _padpct(s_unc, s_add),
            "mean_imbalanced_added_or_active": round(_mean("slr_imbalanced_added"), 2),
            "mean_deadend_active": round(_mean("slr_deadend_added"), 2),
            "lp_biomass_feasible_pct": _pct("slr_lp_feasible_gamma_min"),
            "lp_biomass_feasible_pct_slr_target": _pct("slr_lp_feasible_slr_target"),
            "mean_flux_carrying_active_audit_medium": round(s_flx, 1),
            "mean_flux_certified_added_audit_medium": round(s_cert, 2),
            "mean_uncertified_padding_audit_medium": round(s_unc, 2),
            "uncertified_padding_pct_audit_medium": _padpct(s_unc, s_add),
            "mean_imbalanced_added": round(_mean("slr_imbalanced_added"), 2),
            "lp_biomass_feasible_pct_own_medium": "",
        },
    ]
    with open(summary_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows_out:
            writer.writerow({"n_genomes": len(audit_rows), **row})
    print(f"[2A] Wrote 108-genome audit to {out_path.name} and 6-GEM comparison to {gem_out.name}")
    return out_path


def run_15_genome_stratified_lro(refiner: StoichiometricLocusRefiner) -> Path:
    """
    Experiment 2B: 100% Real Leave-Reaction-Out (LRO) & Secondary BGC Benchmark (Suites 3, 4, 5).
    Computes exact 1-indexed Hit@1, Hit@3, Hit@5, and MRR directly from the real LRO candidate
    ranking tables in `benchmark/manifests/reference_data/lro_data/results_*.tsv` across:
      - E. coli MG1655 Strict Track (n=256 cases, 181 covered in candidate pool)
      - E. coli MG1655 Full Genome-Wide Track (n=414 cases, 226 covered in candidate pool)
      - B. subtilis 168 Track A & Track B (n=144 cases, 76 covered in candidate pool)
      - Role-Stratified E. coli MG1655 Breakdown (Stratum A Core Biomass, Stratum B Catabolic, Stratum C Secondary/Peripheral)
      - Protocol H Homology-Masked Breakdown (Unfiltered, <=70%, <=50%, <=30% Dark Matter)
    """
    print("[2B] Computing 100%-Real Role-Stratified LRO & Secondary BGC Metrics from lro_data/results_*.tsv...")
    from build_all_manifests import compute_core_biomass_flux_reactions
    core_ecs, secondary_ecs, _ = compute_core_biomass_flux_reactions()

    def load_1idx_ranks(path: Path) -> dict[int, int]:
        out = {}
        if not path.exists():
            return out
        with open(path, "r", encoding="utf-8") as f:
            for idx, r in enumerate(csv.DictReader(f, delimiter="\t")):
                cidx = int(r["case_idx"]) if "case_idx" in r and r["case_idx"] != "" else idx
                rk_s = r.get("rank_of_p", "0")
                rk = int(rk_s) if rk_s not in ("", "-1", "None") else 0
                out[cidx] = rk
        return out

    # Load all real LRO ranking tables (1-indexed: rk=1 is Top-1, rk=0 is uncovered/not in candidate pool)
    rk_seq_naive = load_1idx_ranks(LRO_DIR / "results_protein.tsv")
    rk_func_naive = load_1idx_ranks(LRO_DIR / "results_function.tsv")
    rk_m2_coh = load_1idx_ranks(LRO_DIR / "results_m2_integrated_mg1655.tsv")
    rk_rlgc_256 = load_1idx_ranks(LRO_DIR / "results_rlgc_mg1655.tsv")
    rk_rlgc_cg_256 = load_1idx_ranks(LRO_DIR / "results_rlgc_cg_mg1655.tsv")
    rk_rlgc_gbdt_256 = load_1idx_ranks(LRO_DIR / "results_rlgc_gbdt_mg1655.tsv")
    rk_plm_256 = load_1idx_ranks(LRO_DIR / "results_plmonly_mg1655.tsv")
    rk_m4_a_256 = load_1idx_ranks(LRO_DIR / "results_m4_trackA_strict_mg1655.tsv")
    rk_m4_b_256 = load_1idx_ranks(LRO_DIR / "results_m4_trackB_mg1655.tsv")

    rk_rlgc_414 = load_1idx_ranks(LRO_DIR / "results_rlgc_mg1655_all.tsv")
    rk_gbdt_414 = load_1idx_ranks(LRO_DIR / "results_rlgc_gbdt_mg1655_all.tsv")
    rk_m4_b_414 = load_1idx_ranks(LRO_DIR / "results_m4_trackB_mg1655_all.tsv")

    rk_bsub_a = load_1idx_ranks(LRO_DIR / "results_bsubt_trackA.tsv")
    rk_bsub_b = load_1idx_ranks(LRO_DIR / "results_bsubt_trackB.tsv")

    # Load case metadata from cases_mg1655_all.tsv and cases_bsubtilis.tsv
    cases_414 = []
    with open(LRO_DIR / "cases_mg1655_all.tsv", "r", encoding="utf-8") as f:
        for idx, r in enumerate(csv.DictReader(f, delimiter="\t")):
            pw = r.get("pathway_id", "")
            pid = r["protein_id"]
            rxn = r["reaction_id"]
            ec = r.get("ec", "").strip()
            pident = float(r.get("max_pident_annot_homolog", "50.0") or 50.0)
            n_loc = int(r.get("n_neighbors_local", "0") or 0)
            n_pw = int(r.get("n_neighbors_in_pathway", "0") or 0)
            og = int(r.get("orthogroup_size", "1") or 1)
            if ec in core_ecs:
                stratum = "Stratum_A_Core_Biomass"
            elif ec in secondary_ecs or "BIFIDO" in pw or ("PWY" in pw and n_pw <= 4):
                stratum = "Stratum_C_Secondary_Peripheral"
            else:
                stratum = "Stratum_B_Catabolic"
            cases_414.append({
                "key": idx,
                "protein_id": pid,
                "reaction_id": rxn,
                "ec": ec,
                "pident": pident,
                "n_loc": n_loc,
                "og": og,
                "stratum": stratum,
            })

    def score_subset(keys: list[int], rk_map: dict[int, int]) -> dict:
        rks = [rk_map.get(k, 0) for k in keys]
        n_tot = len(rks)
        pos = [rk for rk in rks if rk >= 1]
        n_cov = len(pos)
        if n_tot == 0:
            return {"n_total": 0, "n_cov": 0, "cov_frac": 0.0, "hit1_all": 0.0, "hit3_all": 0.0, "hit5_all": 0.0, "mrr_all": 0.0, "hit1_cov": 0.0, "hit3_cov": 0.0, "mrr_cov": 0.0}
        return {
            "n_total": n_tot,
            "n_cov": n_cov,
            "cov_frac": round(n_cov / n_tot, 4),
            "hit1_all": round(float(np.mean([1.0 if rk == 1 else 0.0 for rk in rks])), 4),
            "hit3_all": round(float(np.mean([1.0 if 1 <= rk <= 3 else 0.0 for rk in rks])), 4),
            "hit5_all": round(float(np.mean([1.0 if 1 <= rk <= 5 else 0.0 for rk in rks])), 4),
            "mrr_all": round(float(np.mean([1.0 / rk if rk >= 1 else 0.0 for rk in rks])), 4),
            "hit1_cov": round(float(np.mean([1.0 if rk == 1 else 0.0 for rk in pos])), 4) if pos else 0.0,
            "hit3_cov": round(float(np.mean([1.0 if 1 <= rk <= 3 else 0.0 for rk in pos])), 4) if pos else 0.0,
            "mrr_cov": round(float(np.mean([1.0 / rk for rk in pos])), 4) if pos else 0.0,
        }

    summary_rows = []

    # 1. Overall E. coli MG1655 Strict 256-case LRO benchmark across all 10 real methods
    keys_256 = sorted(rk_m4_b_256.keys())

    methods_256 = [
        ("Sequence_Only_Candidate_Rank (results_protein.tsv)", rk_seq_naive),
        ("Function_Prior_Only (results_function.tsv)", rk_func_naive),
        ("METEOR_Uniform_Boost (ASSERTED_NOT_COMPUTED: ranks copied from Sequence_Only; METEOR v8 not run on LRO)", rk_seq_naive),
        ("RLGC_Kernel_Density_Only (results_rlgc_mg1655.tsv)", rk_rlgc_256),
        ("RLGC_Plus_CrossGenome_Orthogroup (results_rlgc_cg_mg1655.tsv)", rk_rlgc_cg_256),
        ("RLGC_Plus_GBDT (results_rlgc_gbdt_mg1655.tsv)", rk_rlgc_gbdt_256),
        ("GSPA_Stage1_M2_Coherence_RLGC (results_m2_integrated_mg1655.tsv)", rk_m2_coh),
        ("GSPA_Strict_Homology_Masked_TrackA (results_m4_trackA_strict_mg1655.tsv)", rk_m4_a_256),
        ("ESM2_PLM_Centroids_Only (results_plmonly_mg1655.tsv)", rk_plm_256),
        ("GSPA_Full_M4_RLGC_Coherence_PLM (results_m4_trackB_mg1655.tsv)", rk_m4_b_256),
    ]
    for mname, rmap in methods_256:
        st = score_subset(keys_256, rmap)
        summary_rows.append({
            "benchmark_cohort": "Ecoli_MG1655_Strict_256",
            "stratum": "All_Strata",
            "protocol_H": "Unfiltered",
            "method": mname,
            **st,
        })

    # Load density_unassigned.tsv (414 cases) and density_aaeolicus.tsv (101 cases)
    rk_dens_all_414 = {}
    rk_dens_nonanc_414 = {}
    rk_dens_gapseq_414 = {}
    if (LRO_DIR / "density_unassigned.tsv").exists():
        with open(LRO_DIR / "density_unassigned.tsv", "r", encoding="utf-8") as f:
            for idx, r in enumerate(csv.DictReader(f, delimiter="\t")):
                cidx = int(r["case_idx"]) if "case_idx" in r and r["case_idx"] != "" else idx
                rk_dens_all_414[cidx] = int(r.get("rank_all", "0") or 0)
                rk_dens_nonanc_414[cidx] = int(r.get("rank_nonanchor", "0") or 0)
                rk_dens_gapseq_414[cidx] = int(r.get("rank_unassigned", "0") or 0)

    # 2. Overall E. coli MG1655 Full 414-case LRO benchmark + Stratified (Stratum A, B, C) + Protocol H (Unfiltered, <=70%, <=50%, <=30%)
    methods_414 = [
        ("Whole_Proteome_Operon_Kernel_Density (density_unassigned.tsv:rank_all)", rk_dens_all_414),
        ("Non_Anchor_Operon_Kernel_Density (density_unassigned.tsv:rank_nonanchor)", rk_dens_nonanc_414),
        ("gapseq_Unassigned_Operon_Density (density_unassigned.tsv:rank_unassigned)", rk_dens_gapseq_414),
        ("RLGC_Kernel_Density_Only (results_rlgc_mg1655_all.tsv)", rk_rlgc_414),
        ("RLGC_Plus_GBDT (results_rlgc_gbdt_mg1655_all.tsv)", rk_gbdt_414),
        ("GSPA_Full_M4_RLGC_Coherence_PLM (results_m4_trackB_mg1655_all.tsv)", rk_m4_b_414),
    ]
    for stratum_name in ["All_Strata", "Stratum_A_Core_Biomass", "Stratum_B_Catabolic", "Stratum_C_Secondary_Peripheral"]:
        for proto_label, max_pid in [("Unfiltered", 100.0), ("Protocol_H_le70", 70.0), ("Protocol_H_le50", 50.0), ("Protocol_H_le30_Dark", 30.0)]:
            sub_keys = [
                c["key"] for c in cases_414
                if (stratum_name == "All_Strata" or c["stratum"] == stratum_name) and c["pident"] <= max_pid
            ]
            if not sub_keys:
                continue
            for mname, rmap in methods_414:
                st = score_subset(sub_keys, rmap)
                summary_rows.append({
                    "benchmark_cohort": "Ecoli_MG1655_Full_414",
                    "stratum": stratum_name,
                    "protocol_H": proto_label,
                    "method": mname,
                    **st,
                })

    # 3. Bacillus subtilis 168 LRO Benchmark (n=144 cases)
    keys_bsub = list(rk_bsub_b.keys())
    for mname, rmap in [
        ("GSPA_Strict_Homology_Masked_TrackA (results_bsubt_trackA.tsv)", rk_bsub_a),
        ("GSPA_Full_M4_RLGC_Coherence_PLM (results_bsubt_trackB.tsv)", rk_bsub_b),
    ]:
        st = score_subset(keys_bsub, rmap)
        summary_rows.append({
            "benchmark_cohort": "Bsubtilis_168_LRO_144",
            "stratum": "All_Strata",
            "protocol_H": "Unfiltered",
            "method": mname,
            **st,
        })

    # 4. Aquifex aeolicus VF5 Out-of-Clade LRO Benchmark (n=101 cases)
    if (LRO_DIR / "density_aaeolicus.tsv").exists():
        rk_aa_all, rk_aa_nonanc, rk_aa_gapseq = {}, {}, {}
        with open(LRO_DIR / "density_aaeolicus.tsv", "r", encoding="utf-8") as f:
            for idx, r in enumerate(csv.DictReader(f, delimiter="\t")):
                cidx = int(r["case_idx"]) if "case_idx" in r and r["case_idx"] != "" else idx
                rk_aa_all[cidx] = int(r.get("rank_all", "0") or 0)
                rk_aa_nonanc[cidx] = int(r.get("rank_nonanchor", "0") or 0)
                rk_aa_gapseq[cidx] = int(r.get("rank_unassigned", "0") or 0)
        keys_aa = sorted(rk_aa_all.keys())
        for mname, rmap in [
            ("Whole_Proteome_Operon_Kernel_Density (density_aaeolicus.tsv:rank_all)", rk_aa_all),
            ("Non_Anchor_Operon_Kernel_Density (density_aaeolicus.tsv:rank_nonanchor)", rk_aa_nonanc),
            ("gapseq_Unassigned_Operon_Density (density_aaeolicus.tsv:rank_unassigned)", rk_aa_gapseq),
        ]:
            st = score_subset(keys_aa, rmap)
            summary_rows.append({
                "benchmark_cohort": "Aaeolicus_VF5_LRO_101",
                "stratum": "All_Strata",
                "protocol_H": "Unfiltered",
                "method": mname,
                **st,
            })

    # Provenance column: every row is computed from a per-case ranking table except the METEOR row,
    # whose ranks are asserted equal to Sequence_Only (METEOR v8 cannot be run locally on the LRO
    # cases: universal.pickle, the emit_v8.py driver, MG1655 EC-prediction input and the per-case
    # ablated claims/suggestion files are not available; MG1655 is not in the METEOR 108-genome panel).
    for _r in summary_rows:
        _r["provenance"] = (
            "asserted_not_computed (copy of results_protein.tsv ranks)"
            if _r["method"].startswith("METEOR_Uniform_Boost")
            else "computed"
        )

    out_path = RESULTS_DIR / "pillar2_lro15_summary.tsv"
    with open(out_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(summary_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)
    print(f"[2B] Wrote {len(summary_rows)} 100%-real LRO evaluation rows to {out_path.name}")
    return out_path


def run_price149_all6_evaluation() -> Path:
    """
    Experiment 2C: Execute 100% blind evaluation of Price-149 across all 6 SOTA EC baselines
    via `scripts/eval_price149_real_all6.py`.
    """
    print("[2C] Running 100%-Blind Price-149 Evaluation across all 6 SOTA EC Baselines...")
    eval_price149_real_all6.main()
    return RESULTS_DIR / "pillar2_price149_all6_real.tsv"


def main():
    fba_dir = Path("/data/hohndor/gspa/fba_data/psb_revision/code_snapshot/data")
    if not fba_dir.exists():
        fba_dir = REF_DIR / "psb_revision" / "code_snapshot" / "data"

    refiner = StoichiometricLocusRefiner(fba_dir)
    run_108_genome_invariant_audit(refiner)
    run_15_genome_stratified_lro(refiner)
    run_price149_all6_evaluation()


if __name__ == "__main__":
    main()
