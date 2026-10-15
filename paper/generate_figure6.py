#!/usr/bin/env python3
"""
Generate Figure 6 (fig6_pillar2_3_fba_bacdive.{pdf,png}) for the GSPA v2.1.0 manuscript
directly from canonical results/ TSV files:
  (a) 108-Genome Stoichiometric Invariant Audit: StoichiometricLocusRefiner vs. Grow-Support vs. METEOR
  (b) Re-Audited Price-149 (n=148 RB-TnSeq Proteins) Across All 6 Core Enzyme Predictors
  (c) Role-Stratified LRO (Suites 3-5, n=414) Under Protocol H Homology Masking
  (d) BacDive Experimental Phenotype & Enzyme Activity Validation (Suites 7-8, 529 Carbon & 3,017 Enzyme Type Strains)
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Dict, List

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


REPO_ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = REPO_ROOT / "results"
FIG_DIR = REPO_ROOT / "paper" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

PALETTE = {
    "navy": "#1B365D",
    "blue": "#0072B2",
    "sky": "#56B4E9",
    "teal": "#009E73",
    "orange": "#E69F00",
    "vermilion": "#D55E00",
    "purple": "#CC79A7",
    "grey": "#7F8C8D",
    "lightgrey": "#E5E8E8",
    "dark": "#2C3E50",
}

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "font.size": 9.5,
    "axes.titlesize": 10.5,
    "axes.titleweight": "bold",
    "axes.labelsize": 9.5,
    "axes.labelweight": "bold",
    "xtick.labelsize": 8.5,
    "ytick.labelsize": 8.5,
    "legend.fontsize": 7.6,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.9,
    "grid.alpha": 0.25,
    "grid.linestyle": "--",
    "figure.dpi": 300,
})


def read_tsv(path: Path) -> List[Dict[str, str]]:
    with open(path, "r", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def add_panel_label(ax, label: str, x: float = -0.14, y: float = 1.09):
    ax.text(
        x, y, label,
        transform=ax.transAxes,
        fontsize=13,
        fontweight="bold",
        va="top",
        ha="left",
        color="#111111",
    )


def generate_figure6():
    fig, axes = plt.subplots(2, 2, figsize=(13.6, 9.6))
    plt.subplots_adjust(wspace=0.28, hspace=0.38)

    # -----------------------------------------------------------------
    # Panel (a): 108-Genome Stoichiometric Invariant Audit (Read from TSV)
    # -----------------------------------------------------------------
    ax_a = axes[0, 0]
    add_panel_label(ax_a, "a")
    inv_rows = {r["method"]: r for r in read_tsv(RESULTS_DIR / "pillar2_invariant_audit_summary.tsv")}
    r_met = inv_rows["METEOR v8 (Negative-Cost MILP)"]
    r_grw = inv_rows["Grow-Support Threshold Repair (Positive-Cost LP)"]
    r_slr = inv_rows["GSPA StoichiometricLocusRefiner (Positive-Cost + pFBA Support Certificate)"]

    categories = [
        "Added Rxns\n(Mean / Genome)",
        "Flux-Certified\nAdded (|v|>1e-7)",
        "Zero-Flux\nPadding Rxns",
        "Imbalanced\nActive Rxns",
    ]
    meteor_vals = [
        float(r_met["mean_added_rxns"]),
        float(r_met["mean_flux_certified_added"]),
        float(r_met["mean_uncertified_padding"]),
        float(r_met["mean_imbalanced_added_or_active"]),
    ]
    grow_vals = [
        float(r_grw["mean_added_rxns"]),
        float(r_grw["mean_flux_certified_added"]),
        float(r_grw["mean_uncertified_padding"]),
        float(r_grw["mean_imbalanced_added_or_active"]),
    ]
    slr_vals = [
        float(r_slr["mean_added_rxns"]),
        float(r_slr["mean_flux_certified_added"]),
        float(r_slr["mean_uncertified_padding"]),
        float(r_slr["mean_imbalanced_added_or_active"]),
    ]

    x = np.arange(len(categories))
    w = 0.26
    met_pad_pct = float(r_met["uncertified_padding_pct"])
    bars_m = ax_a.bar(x - w, meteor_vals, w, label=f"METEOR v8 (Negative-Cost MILP, {met_pad_pct:.1f}% Padding)", color=PALETTE["vermilion"], edgecolor="#111111")
    bars_g = ax_a.bar(x, grow_vals, w, label="Grow-Support Threshold Repair (Positive-Cost LP)", color=PALETTE["orange"], edgecolor="#111111", hatch="//")
    bars_s = ax_a.bar(x + w, slr_vals, w, label="GSPA StoichiometricLocusRefiner (pFBA, 100% Cert.)", color=PALETTE["teal"], edgecolor="#111111")

    ax_a.set_yscale("symlog", linthresh=1.0)
    ax_a.set_ylim(0, 7500)
    ax_a.set_xticks(x)
    ax_a.set_xticklabels(categories, fontsize=8.0)
    ax_a.set_ylabel("Reactions per Genome (symlog scale)")
    ax_a.set_title("108-Genome Stoichiometric Audit: Flux Certification vs. Padding")
    ax_a.grid(axis="y")

    for b, val in zip(bars_m, meteor_vals):
        ax_a.text(b.get_x() + b.get_width() / 2, val * 1.15, f"{val:.1f}", ha="center", va="bottom", fontsize=7.0, fontweight="bold", color=PALETTE["vermilion"])
    for b, val in zip(bars_g, grow_vals):
        ypos = val * 1.15 if val > 0 else 0.12
        ax_a.text(b.get_x() + b.get_width() / 2, ypos, f"{val:.1f}" if val > 0 else "0.0", ha="center", va="bottom", fontsize=7.0, fontweight="bold", color="#8A5A00")
    for b, val in zip(bars_s, slr_vals):
        ypos = val * 1.15 if val > 0 else 0.12
        ax_a.text(b.get_x() + b.get_width() / 2, ypos, f"{val:.2f}" if val > 0 else "0.0", ha="center", va="bottom", fontsize=7.0, fontweight="bold", color=PALETTE["navy"])
    ax_a.legend(loc="upper left", frameon=True, facecolor="white", framealpha=0.92, fontsize=7.2)

    # -----------------------------------------------------------------
    # Panel (b): Price-149 (n=148) Across All 6 Core Baselines (Read from TSV)
    # -----------------------------------------------------------------
    ax_b = axes[0, 1]
    add_panel_label(ax_b, "b")
    price_rows = read_tsv(RESULTS_DIR / "pillar2_price149_all6_baselines.tsv")
    baselines = [
        ("CLEAN (2023)", "CLEAN\n(2023)"),
        ("MAPred (2025)", "MAPred\n(2025)"),
        ("DeepProZyme (2023)", "DeepProZyme\n(2023/24)"),
        ("EnzBERT (2023)", "EnzBERT\n(2023)"),
        ("GraphEC (2024)", "GraphEC\n(2024)"),
        ("TopEC (2025)", "TopEC\n(2025)"),
    ]
    p_lookup = {(r["baseline"], r["method"]): float(r["fmax"]) for r in price_rows}
    base_fmax = [p_lookup[(b, "Baseline")] for b, _ in baselines]
    meteor_fmax = [p_lookup[(b, "+ METEOR v8 (Iter 2)")] for b, _ in baselines]
    slr_cog_fmax = [p_lookup[(b, "+ GSPA SLR (Cognate Substrate)")] for b, _ in baselines]
    gspa_single_fmax = [p_lookup[(b, "+ Full GSPA (Iter 3)")] for b, _ in baselines]
    gspa_multi_fmax = [p_lookup[(b, "+ Full GSPA Multi-Pred (Iter 3)")] for b, _ in baselines]

    x_b = np.arange(len(baselines))
    wb = 0.165
    ax_b.bar(x_b - 2 * wb, base_fmax, wb, label="Baseline Predictor", color=PALETTE["grey"])
    ax_b.bar(x_b - wb, meteor_fmax, wb, label="+ METEOR v8 (Single-Pred)", color=PALETTE["orange"], hatch="//")
    ax_b.bar(x_b, slr_cog_fmax, wb, label="+ GSPA SLR (Cognate)", color=PALETTE["sky"])
    ax_b.bar(x_b + wb, gspa_single_fmax, wb, label="+ GSPA Single-Pred (Iter 3)", color=PALETTE["teal"])
    bars_g = ax_b.bar(x_b + 2 * wb, gspa_multi_fmax, wb, label="+ GSPA Multi-Pred (Iter 3)", color=PALETTE["navy"])

    for b, v_base, v_full in zip(bars_g, base_fmax, gspa_multi_fmax):
        delta = v_full - v_base
        ax_b.text(b.get_x() + b.get_width() / 2, v_full + 0.012, f"+{delta:.2f}", ha="center", va="bottom", fontsize=6.8, fontweight="bold", color=PALETTE["navy"])

    ax_b.set_ylim(0, 0.78)
    ax_b.set_xticks(x_b)
    ax_b.set_xticklabels([lbl for _, lbl in baselines], fontsize=7.8, fontweight="bold")
    ax_b.set_ylabel("Enzyme EC Fmax (Price-149, n=148)")
    ax_b.set_title("Price-149 RB-TnSeq Benchmark Across 6 Enzyme Predictors")
    ax_b.grid(axis="y")
    ax_b.legend(loc="upper right", ncol=2, frameon=True, facecolor="white", framealpha=0.95, fontsize=6.8)

    # -----------------------------------------------------------------
    # Panel (c): Role-Stratified LRO Under Protocol H (Read from pillar2_lro15_summary.tsv)
    # -----------------------------------------------------------------
    ax_c = axes[1, 0]
    add_panel_label(ax_c, "c")
    lro15 = {(r["stratum"], r["protocol_H"], r["method"].split(" ")[0]): r
             for r in read_tsv(RESULTS_DIR / "pillar2_lro15_summary.tsv")
             if r["benchmark_cohort"] == "Ecoli_MG1655_Full_414"}
    prot_keys = ["Unfiltered", "Protocol_H_le70", "Protocol_H_le50", "Protocol_H_le30_Dark"]
    regimes = ["Unfiltered\n(100%)", "Protocol H\n<=70% ID", "Protocol H\n<=50% ID", "Protocol H\n<=30% Dark"]

    rlgc_all_pool = [float(lro15[("All_Strata", pk, "RLGC_Kernel_Density_Only")]["hit1_cov"]) for pk in prot_keys]
    gbdt_all_pool = [float(lro15[("All_Strata", pk, "RLGC_Plus_GBDT")]["hit1_cov"]) for pk in prot_keys]
    gspa_all_pool = [float(lro15[("All_Strata", pk, "GSPA_Full_M4_RLGC_Coherence_PLM")]["hit1_cov"]) for pk in prot_keys]
    gspa_stratumA_pool = [float(lro15[("Stratum_A_Core_Biomass", pk, "GSPA_Full_M4_RLGC_Coherence_PLM")]["hit1_cov"]) for pk in prot_keys]
    gspa_all_total = [float(lro15[("All_Strata", pk, "GSPA_Full_M4_RLGC_Coherence_PLM")]["hit1_all"]) for pk in prot_keys]

    xc = np.arange(len(regimes))
    ax_c.plot(xc, gbdt_all_pool, "s:", color=PALETTE["orange"], linewidth=1.8, markersize=6, label="RLGC + GBDT (All Strata, Given Pool)")
    ax_c.plot(xc, rlgc_all_pool, "^-.", color=PALETTE["sky"], linewidth=2.0, markersize=6.5, label="RLGC Kernel Only (All Strata, Given Pool)")
    ax_c.plot(xc, gspa_all_total, "D--", color=PALETTE["teal"], linewidth=2.0, markersize=6.5, label="GSPA Full M4 (All Strata, Hit@1 All)")
    ax_c.plot(xc[:3], gspa_stratumA_pool[:3], "v-", color=PALETTE["purple"], linewidth=2.2, markersize=7.0, label="GSPA Full M4 (Stratum A Core, Given Pool)")
    ax_c.plot(xc, gspa_all_pool, "o-", color=PALETTE["navy"], linewidth=2.5, markersize=7.5, label="GSPA Full M4 (All Strata, Hit@1 Given Pool)")

    ax_c.set_ylim(-0.04, 1.34)
    ax_c.set_xticks(xc)
    ax_c.set_xticklabels(regimes, fontsize=8.2)
    ax_c.set_ylabel("Enzyme Gene Attribution Hit@1")
    ax_c.set_title("Role-Stratified LRO (n=414) Under Protocol H Homology Masking")
    ax_c.grid(True)
    ax_c.legend(loc="upper right", frameon=True, facecolor="white", framealpha=0.92, fontsize=6.9)

    # -----------------------------------------------------------------
    # Panel (d): BacDive Phenotype & Enzyme Assays (Read from Track 1 & 2 TSVs)
    # -----------------------------------------------------------------
    ax_d = axes[1, 1]
    add_panel_label(ax_d, "d")
    t1_lookup = {(r["dataset"], r["method"]): r for r in read_tsv(RESULTS_DIR / "pillar3_bacdive_track1_summary.tsv")}
    t2_lookup = {(r["dataset"], r["method"]): r for r in read_tsv(RESULTS_DIR / "pillar3_bacdive_track2_summary.tsv")}

    ds_c1 = "ProTraits_CarbonSource_529Strains"
    ds_e1 = "BacDive_Enzyme_3017TypeStrains"
    ds_c2 = "ProTraits_CarbonSource_HeldOut70pct"
    ds_e2 = "BacDive_Enzyme_HeldOut70pct"

    def get_bd_vals(m_t1_c, m_t1_e, m_t2_c, m_t2_e):
        return [
            float(t1_lookup[(ds_c1, m_t1_c)]["f1"]),
            float(t1_lookup[(ds_e1, m_t1_e)]["f1"]),
            float(t2_lookup[(ds_c2, m_t2_c)]["f1"]),
            float(t2_lookup[(ds_e2, m_t2_e)]["f1"]),
            float(t2_lookup[(ds_e2, m_t2_e)]["roc_auc"]),
        ]

    v_carveme = get_bd_vals(
        "CarveMe (Machado et al. 2018)",
        "CarveMe (Machado et al. 2018)",
        "CarveMe (Machado et al. 2018) [Held-Out 70%]",
        "CarveMe (Machado et al. 2018) [Held-Out 70%]",
    )
    v_gapseq = get_bd_vals(
        "gapseq (Zimmermann et al. 2021)",
        "gapseq (Zimmermann et al. 2021)",
        "gapseq (Zimmermann et al. 2021) [Held-Out 70%]",
        "gapseq (Zimmermann et al. 2021) [Held-Out 70%]",
    )
    v_gspa = get_bd_vals(
        "GSPA Track 1 (100% Blind Multi-Model + Catabolic Scope & Taxon Prior)",
        "GSPA Track 1 (100% Blind Multi-Model + Clade Prior)",
        "GSPA Track 2 (30% Revealed Bayesian Anchor -> Evaluated on Held-Out 70%)",
        "GSPA Track 2 (30% Revealed Bayesian Anchor -> Evaluated on Held-Out 70%)",
    )
    v_diaspora = get_bd_vals(
        "DiASPora GBDT (Koblitz et al. 2025) [Genus-Blocked CV]",
        "DiASPora GBDT (Koblitz et al. 2025) [Genus-Blocked CV]",
        "DiASPora GBDT (Koblitz et al. 2025) [Trained on 30% Revealed -> Held-Out 70%]",
        "DiASPora GBDT (Koblitz et al. 2025) [Trained on 30% Revealed -> Held-Out 70%]",
    )
    v_hybrid = get_bd_vals(
        "GSPA + DiASPora Neuro-Symbolic Hybrid (2026) [Genus-Blocked CV]",
        "GSPA + DiASPora Neuro-Symbolic Hybrid (2026) [Genus-Blocked CV]",
        "GSPA Track 2 + DiASPora Hybrid (30% Revealed -> Evaluated on Held-Out 70%)",
        "GSPA Track 2 + DiASPora Hybrid (30% Revealed -> Evaluated on Held-Out 70%)",
    )

    bd_tasks = [
        "Track 1 Carbon\nGrowth F1\n(n=2,332)",
        "Track 1 Enzyme\nActivity F1\n(n=17,496)",
        "Track 2 Carbon\nHeld-Out 70% F1\n(n=1,641)",
        "Track 2 Enzyme\nHeld-Out 70% F1\n(n=12,316)",
        "Track 2 Enzyme\nHeld-Out AUC\n(n=12,316)",
    ]

    xd = np.arange(len(bd_tasks))
    wd = 0.165
    ax_d.bar(xd - 2 * wd, v_carveme, wd, label="CarveMe (Unsupervised)", color=PALETTE["lightgrey"], edgecolor=PALETTE["grey"])
    ax_d.bar(xd - wd, v_gapseq, wd, label="gapseq (Unsupervised)", color=PALETTE["sky"], edgecolor=PALETTE["blue"])
    ax_d.bar(xd, v_gspa, wd, label="GSPA Track 1 / Track 2 (Bayesian Anchor)", color=PALETTE["teal"], edgecolor=PALETTE["navy"])
    ax_d.bar(xd + wd, v_diaspora, wd, label="DiASPora GBDT (Koblitz 2025, Genus-CV)", color=PALETTE["orange"], edgecolor="#111111", hatch="//")
    bars_hyb = ax_d.bar(xd + 2 * wd, v_hybrid, wd, label="GSPA + DiASPora Hybrid (Genus-CV)", color=PALETTE["navy"])

    for b, val in zip(bars_hyb, v_hybrid):
        ax_d.text(b.get_x() + b.get_width() / 2, val + 0.012, f"{val:.3f}", ha="center", va="bottom", fontsize=6.8, fontweight="bold", color=PALETTE["navy"])

    ax_d.set_ylim(0, 1.28)
    ax_d.set_xticks(xd)
    ax_d.set_xticklabels(bd_tasks, fontsize=7.6)
    ax_d.set_ylabel("F1 / ROC-AUC Score")
    ax_d.set_title("BacDive Validation Across 529 Carbon & 3,017 Enzyme Type Strains")
    ax_d.grid(axis="y")
    ax_d.legend(loc="upper left", ncol=2, frameon=True, facecolor="white", framealpha=0.95, fontsize=6.8)

    pdf_path = FIG_DIR / "fig6_pillar2_3_fba_bacdive.pdf"
    png_path = FIG_DIR / "fig6_pillar2_3_fba_bacdive.png"
    fig.savefig(pdf_path, bbox_inches="tight")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[figure6] Saved {pdf_path} and {png_path}")


if __name__ == "__main__":
    generate_figure6()
