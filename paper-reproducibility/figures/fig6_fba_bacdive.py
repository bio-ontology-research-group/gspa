"""Figure 6: stoichiometric flux certification, fitness-validated enzyme discovery, phenotype prediction.

a  results/pillar2_invariant_audit_summary.tsv  method=METEOR v8 / Grow-Support Threshold Repair / GSPA StoichiometricLocusRefiner
b  results/pillar2_price149_all6_real.tsv       extraction_version=<newest>; method=Baseline / + Stage 1 Noisy-OR Ensemble
                                                (Sequence Only) / pool_without_this_baseline / + Full GSPA Multi-Pred (Iter 3)
c  results/pillar2_lro15_summary.tsv            benchmark_cohort=Ecoli_MG1655_Full_414; stratum=All_Strata;
                                                method=GSPA_Full_M4_RLGC_Coherence_PLM*; protocol_H tiers
d  results/pillar3_bacdive_track1_summary.tsv / pillar3_bacdive_track2_summary.tsv (rows of tab:bacdive, no 10-EC panel)
"""

from __future__ import annotations

import re

import numpy as np
from matplotlib.lines import Line2D

import common as C

AUD = "pillar2_invariant_audit_summary.tsv"
P149 = "pillar2_price149_all6_real.tsv"
LRO = "pillar2_lro15_summary.tsv"
T1 = "pillar3_bacdive_track1_summary.tsv"
T2 = "pillar3_bacdive_track2_summary.tsv"

AUDIT_METHODS = [
    ("METEOR v8 (Negative-Cost MILP)", "METEOR v8 (negative-cost MILP)", C.ORANGE),
    ("Grow-Support Threshold Repair (Positive-Cost LP)", "Threshold repair (positive-cost LP)", C.AQUA),
    ("GSPA StoichiometricLocusRefiner (Positive-Cost + pFBA Support Certificate)",
     "GSPA Stage 2B (lexicographic pFBA)", C.BLUE4),
]
AUDIT_METRICS = [
    ("mean_added_rxns", "Added"),
    ("mean_flux_certified_added", "Flux-\ncertified"),
    ("mean_uncertified_padding", "Zero-flux\npadding"),
    ("mean_imbalanced_added", "Imbalanced\nadded"),
]
PREDICTORS = ["CLEAN (2023)", "MAPred (2025)", "DeepProZyme (2023)", "EnzBERT (2023)", "GraphEC (2024)", "TopEC (2025)"]
ARMS = [
    ("Baseline", "Standalone", C.GREY),
    ("pool_without_this_baseline", "Pool without this predictor", C.BLUE1),
    ("+ Stage 1 Noisy-OR Ensemble (Sequence Only)", "+ Stage 0 pooling", C.BLUE2),
    ("+ Full GSPA Multi-Pred (Iter 3)", "+ GSPA multi-predictor", C.BLUE4),
]
TIERS = [("Unfiltered", "none"), ("Protocol_H_le70", "≤70"), ("Protocol_H_le50", "≤50"),
         ("Protocol_H_le30_Dark", "≤30")]
BACDIVE = [  # (track file, method selector, label as in tab:bacdive, is GSPA)
    (T1, C.prefix("CarveMe (Machado et al. 2018)"), "CarveMe calls (binary)", False),
    (T1, C.prefix("ModelSEED (Henry et al. 2010)"), "ModelSEED calls (binary)", False),
    (T1, C.prefix("gapseq (Zimmermann et al. 2021)"), "gapseq calls (binary, reference)", False),
    (T1, C.startswith("GSPA Track 1 (100% Blind"), "GSPA Track 1 (label-free pooling + coupling)", True),
    (T1, C.prefix("BacDive-AI RF (Schober et al. 2025)"), "Random forest on reconstruction calls", False),
    (T1, C.prefix("DiASPora GBDT (Koblitz et al. 2025)"), "GBDT on reconstruction calls", False),
    (T1, C.prefix("Traitar L1-SVM (Weimann et al. 2016)"), "L1-SVM on reconstruction calls", False),
    (T1, C.prefix("GSPA + DiASPora Neuro-Symbolic Hybrid (2026)"), "GSPA posteriors + GBDT", True),
    (T2, C.prefix("gapseq (Zimmermann et al. 2021)"), "gapseq calls (held-out 70%)", False),
    (T2, C.prefix("DiASPora GBDT (Koblitz et al. 2025)"), "GBDT on reconstruction calls (30% revealed)", False),
    (T2, "GSPA Track 2 (30% Revealed Bayesian Anchor -> Evaluated on Held-Out 70%)", "GSPA Track 2", True),
    (T2, "GSPA Track 2 + DiASPora Hybrid (30% Revealed -> Evaluated on Held-Out 70%)", "GSPA Track 2 + GBDT", True),
]
DATASETS = {  # facet -> (track-1 dataset, track-2 dataset)
    "Carbon-source growth": ("ProTraits_CarbonSource_529Strains", "ProTraits_CarbonSource_HeldOut70pct"),
    "Diagnostic enzymes": ("BacDive_Enzyme_3017TypeStrains", "BacDive_Enzyme_HeldOut70pct"),
}


def _vkey(v):
    mo = re.match(r"v(\d+)(.*)", v or "")
    return (int(mo.group(1)), mo.group(2)) if mo else (-1, v or "")


def newest_extraction():
    return max({r["extraction_version"] for r in C.load(P149)} - {""}, key=_vkey)


def panel_a(sf):
    C.panel("fig6", "a")
    ax, axs = sf.subplots(2, 1, gridspec_kw={"height_ratios": [3.6, 1.0], "hspace": 0.08})
    x = np.arange(len(AUDIT_METRICS))
    w = 0.8 / len(AUDIT_METHODS)
    for k, (m, lab, color) in enumerate(AUDIT_METHODS):
        vals = [C.get(AUD, {"method": m}, col) for col, _ in AUDIT_METRICS]
        xs = x + (k - (len(AUDIT_METHODS) - 1) / 2) * w
        ax.bar(xs, vals, width=w * 0.9, color=color, label=lab, zorder=2)
        for xi, v in zip(xs, vals):
            ax.text(xi, v * 1.12 + 0.05, f"{v:.1f}" if v < 100 else f"{v:.0f}", ha="center", va="bottom",
                    fontsize=C.BASE - 3, color=C.INK2)
    ax.set_yscale("symlog", linthresh=1)
    ax.set_ylim(0, 3000)
    ax.set_yticks([0, 1, 10, 100, 1000])
    ax.set_yticklabels(["0", "1", "10", "100", "1000"])
    ax.set_xticks(x)
    ax.set_xticklabels([lab for _, lab in AUDIT_METRICS])
    n = int(C.get(AUD, {"method": AUDIT_METHODS[0][0]}, "n_genomes"))
    ax.set_ylabel("Reactions per genome")
    C.hgrid(ax)
    C.top_legend(ax, ncol=2)
    C.title(ax, f"Stoichiometric audit of network repair (mean of {n} genomes)", pad=50)
    # feasibility strip: growth target reached on own medium vs defined audit medium
    cols = [("lp_biomass_feasible_pct_own_medium", "Own medium"), ("lp_biomass_feasible_pct", "Defined audit medium")]
    mat = np.full((len(AUDIT_METHODS), len(cols)), np.nan)
    for i, (m, _, _) in enumerate(AUDIT_METHODS):
        for j, (col, _) in enumerate(cols):
            v = C.get(AUD, {"method": m}, col)
            mat[i, j] = np.nan if v is None else v
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v = mat[i, j]
            if np.isnan(v):  # GSPA solves on the defined medium: its own medium is the audit medium
                txt, fc, tc = "= audit medium", "white", C.MUTED
            else:
                fc = C.BLUE1 if v >= 99.95 else "#f6d3c4"
                txt, tc = f"{v:.0f}%", C.INK
            axs.add_patch(__import__("matplotlib.patches", fromlist=["Rectangle"]).Rectangle(
                (j - 0.48, i - 0.42), 0.96, 0.84, facecolor=fc, edgecolor=C.AXIS if fc == "white" else fc, lw=0.8))
            axs.text(j, i, txt, ha="center", va="center", fontsize=C.BASE - 2.5, color=tc)
    axs.set_xlim(-0.5, len(cols) - 0.5)
    axs.set_ylim(len(AUDIT_METHODS) - 0.5, -0.5)
    axs.set_xticks(range(len(cols)))
    axs.set_xticklabels([c[1] for c in cols], fontsize=C.BASE - 2)
    axs.xaxis.tick_top()
    axs.set_yticks(range(len(AUDIT_METHODS)))
    axs.set_yticklabels([lab.split(" (")[0] for _, lab, _ in AUDIT_METHODS], fontsize=C.BASE - 2.5)
    axs.tick_params(length=0)
    for s in axs.spines.values():
        s.set_visible(False)
    axs.set_xlabel("Genomes reaching the growth target (feasibility)", fontsize=C.BASE - 2, color=C.INK2)


def panel_b(sf):
    C.panel("fig6", "b")
    ver = newest_extraction()
    C.note(P149, "all rows", "extraction_version", "newest version used", ver)
    axt, axb = sf.subplots(2, 1, sharex=True, gridspec_kw={"hspace": 0.1})
    x = np.arange(len(PREDICTORS))
    w = 0.8 / len(ARMS)
    for k, (meth, lab, color) in enumerate(ARMS):
        xs = x + (k - (len(ARMS) - 1) / 2) * w
        top1 = [C.get(P149, {"extraction_version": ver, "baseline": p, "method": meth}, "top1_count") for p in PREDICTORS]
        fm = [C.get(P149, {"extraction_version": ver, "baseline": p, "method": meth}, "fmax_dataset_macro")
              for p in PREDICTORS]
        axt.bar(xs, top1, width=w * 0.9, color=color, label=lab, zorder=2)
        axb.bar(xs, fm, width=w * 0.9, color=color, zorder=2)
    n = int(C.get(P149, {"extraction_version": ver, "baseline": PREDICTORS[0], "method": "Baseline"}, "n_proteins"))
    axt.set_ylabel(f"Top-1 correct\n(of {n})")
    axb.set_ylabel(r"Macro $F_{\max}$")
    axb.set_ylim(0, 0.6)
    axb.set_xticks(x)
    axb.set_xticklabels([p.replace(" (", "\n(") for p in PREDICTORS], fontsize=C.BASE - 2)
    axt.tick_params(axis="x", length=0)
    C.hgrid(axt)
    C.hgrid(axb)
    C.top_legend(axt, ncol=2)
    C.title(axt, "Price-149 fitness-validated enzymes", pad=46)


def panel_c(sf):
    C.panel("fig6", "c")
    axl, axr = sf.subplots(1, 2, gridspec_kw={"wspace": 0.1})
    meth = "GSPA_Full_M4_RLGC_Coherence_PLM"
    hit, cov, labs = [], [], []
    for tier, lab in TIERS:
        crit = {"benchmark_cohort": "Ecoli_MG1655_Full_414", "stratum": "All_Strata", "protocol_H": tier,
                "method": C.prefix(meth)}
        hit.append(C.get(LRO, crit, "hit1_cov"))
        cov.append(C.get(LRO, crit, "cov_frac"))
        ncov, ntot = int(C.get(LRO, crit, "n_cov")), int(C.get(LRO, crit, "n_total"))
        labs.append((lab, ncov, ntot))
    x = np.arange(len(TIERS))
    ramp = [C.BLUE4, C.BLUE3, C.BLUE2, C.BLUE1]
    axl.bar(x, hit, width=0.66, color=ramp, zorder=2)
    for xi, v in zip(x, hit):
        axl.text(xi, v + 0.02, f"{v:.2f}", ha="center", va="bottom", fontsize=C.BASE - 2.5, color=C.INK2)
    axr.bar(x, cov, width=0.66, color=ramp, zorder=2)
    for xi, v, (_, nc, nt) in zip(x, cov, labs):
        axr.text(xi, v + 0.02, f"{nc}/{nt}", ha="center", va="bottom", fontsize=C.BASE - 2.5, color=C.INK2)
    for ax, yl in [(axl, "Hit@1 within candidate pool"), (axr, "Reactions with true gene in pool")]:
        ax.set_xticks(x)
        ax.set_xticklabels([l for l, _, _ in labs])
        ax.set_ylim(0, 1.08)
        ax.set_ylabel(yl)
        ax.set_xlabel("Donor identity cap (%)", fontsize=C.BASE - 1)
        C.hgrid(ax)
    n = labs[0][2]
    C.title(axl, f"Donor-homology masking (Protocol H, {n} $\\it{{E.\\ coli}}$ reactions)", pad=10)


def panel_d(sf):
    C.panel("fig6", "d")
    axes = sf.subplots(1, len(DATASETS), sharey=True, gridspec_kw={"wspace": 0.08})
    # row layout with two group headers (Track 1, Track 2)
    labels, ypos, rows = [], [], []
    yy = 0
    for hdr, idx in [("Track 2: 30% of assays revealed, held-out 70%", range(11, 7, -1)),
                     ("Track 1: label-free and genus-blocked CV", range(7, -1, -1))]:
        for i in idx:
            rows.append((yy, BACDIVE[i]))
            labels.append(BACDIVE[i][2])
            ypos.append(yy)
            yy += 1
        labels.append(hdr)
        ypos.append(yy)
        yy += 1.3
    for ax, (facet, (d1, d2)) in zip(axes, DATASETS.items()):
        allv = []
        for yi, (tr, sel, lab, gspa) in rows:
            crit = {"dataset": d1 if tr == T1 else d2, "method": sel}
            if not C.rows(tr, crit):
                ax.text(0.5, yi, "not evaluated", va="center", fontsize=C.BASE - 3, color=C.MUTED, style="italic")
                continue
            f1, auc = C.get(tr, crit, "f1"), C.get(tr, crit, "roc_auc")
            allv += [f1, auc]
            col = C.BLUE4 if gspa else C.GREY_D
            ax.plot(f1, yi + 0.15, "o", color=col, ms=7.5, mec="white", mew=1.0, zorder=3)
            ax.plot(auc, yi - 0.15, "D", color=col, ms=6.5, mfc="white", mec=col, mew=1.6, zorder=3)
        ax.set_xlim(0.4, 1.0)
        ax.set_ylim(-0.7, yy - 1.0)
        C.hgrid(ax, "x")
        ax.set_xlabel(facet, fontsize=C.BASE - 1)
        ax.tick_params(axis="y", length=0)
    axes[0].set_yticks(ypos)
    axes[0].set_yticklabels(labels, fontsize=C.BASE - 2.5)
    gspa_labels = {b[2] for b in BACDIVE if b[3]}
    for t in axes[0].get_yticklabels():
        if t.get_text() in gspa_labels:
            t.set_color(C.INK)
            t.set_fontweight("bold")
        elif t.get_text().startswith("Track"):
            t.set_color(C.INK)
            t.set_fontstyle("italic")
    h = [Line2D([], [], ls="none", marker="o", color=C.GREY_D, ms=7, label="$F_1$"),
         Line2D([], [], ls="none", marker="D", mfc="white", mec=C.GREY_D, mew=1.6, ms=6, label="ROC-AUC")]
    C.top_legend(axes[0], ncol=2, handles=h)
    C.title(axes[0], "BacDive phenotype prediction", pad=30)


def generate():
    C.apply_style()
    fig, sf = C.figure()
    panel_a(sf[0])
    panel_b(sf[1])
    panel_c(sf[2])
    panel_d(sf[3])
    C.letters(sf)
    C.save(fig, "fig6_pillar2_3_fba_bacdive")
