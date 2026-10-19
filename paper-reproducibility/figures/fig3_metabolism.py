"""Figure 3: leave-reaction-out metabolic enzyme attribution.

a  results/pillar2_lro15_summary.tsv   benchmark_cohort=Ecoli_MG1655_Strict_256; stratum=All_Strata; protocol_H=Unfiltered
b  results/pillar2_lro15_summary.tsv   benchmark_cohort=Ecoli_MG1655_Full_414; protocol_H=Unfiltered; stratum=A/B/C
c  results/pillar2_lro15_summary.tsv   benchmark_cohort=Bsubtilis_168_LRO_144 (Track A / Track B)
   results/pillar4_v2_ml_baselines.tsv  method=Single Feature: ... (lro50_hit1_pct)
   results/pillar4_v2_multiseed_summary.tsv seed_label=MEAN_SD (lro_hit1_pct)
d  results/pillar3_dark_matter_mags_v2.tsv (phylum_or_group with n_mags >= 4, plus ALL) and
   results/pillar3_dark_matter_mags_v2_per_mag.tsv (new_fills_q050 per MAG)
"""

from __future__ import annotations

import numpy as np
from matplotlib.lines import Line2D

import common as C

LRO = "pillar2_lro15_summary.tsv"
ML = "pillar4_v2_ml_baselines.tsv"
SEED = "pillar4_v2_multiseed_summary.tsv"
MAG = "pillar3_dark_matter_mags_v2.tsv"
MAGP = "pillar3_dark_matter_mags_v2_per_mag.tsv"

STRICT = [  # (method prefix, label as in tab:metabolism)
    ("Sequence_Only_Candidate_Rank", "No target evidence"),
    ("GSPA_Strict_Homology_Masked_TrackA", "Track A (homology- and embedding-masked)"),
    ("RLGC_Kernel_Density_Only", "RLGC operon kernel density only"),
    ("GSPA_Stage1_M2_Coherence_RLGC", "Stage 1 coherence + RLGC (no PLM)"),
    ("ESM2_PLM_Centroids_Only", "ESM-2 reaction centroids only (unmasked)"),
    ("GSPA_Full_M4_RLGC_Coherence_PLM", "Full GSPA (Track B)"),
]
STRATA = [
    ("Stratum_A_Core_Biomass", "A: core biomass"),
    ("Stratum_B_Catabolic", "B: alternative catabolic"),
    ("Stratum_C_Secondary_Peripheral", "C: peripheral / secondary"),
]
STRATA_METHODS = [
    ("gapseq_Unassigned_Operon_Density", "gapseq unassigned operon density", C.ORANGE),
    ("RLGC_Kernel_Density_Only", "RLGC operon kernel density only", C.BLUE1),
    ("GSPA_Full_M4_RLGC_Coherence_PLM", "Full GSPA (Track B)", C.BLUE4),
]
SINGLE = [
    ("Single Feature: RLGC Suggestion Score (Rank Candidates by suggestion_score)", "Stage 1 suggestion score"),
    ("Single Feature: Stage 1 Final Log-Odds (Rank Candidates by final_logodds)", "Stage 1 final log-odds"),
    ("Single Feature: ESM-2 Reaction-Centroid Cosine (Rank Candidates by esm2_centroid calibrated_prob)",
     "ESM-2 reaction-centroid similarity"),
]


def _crit(cohort, stratum, proto, meth):
    return {"benchmark_cohort": cohort, "stratum": stratum, "protocol_H": proto, "method": C.prefix(meth)}


def panel_a(ax):
    C.panel("fig3", "a")
    metrics = [("hit1_all", "Hit@1, all reactions", C.BLUE4), ("hit1_cov", "Hit@1, within candidate pool", C.BLUE2),
               ("hit3_cov", "Hit@3, within candidate pool", C.BLUE1)]
    y = np.arange(len(STRICT))[::-1]
    h = 0.27
    for k, (col, lab, color) in enumerate(metrics):
        vals = [C.get(LRO, _crit("Ecoli_MG1655_Strict_256", "All_Strata", "Unfiltered", m), col) for m, _ in STRICT]
        ax.barh(y + (1 - k) * h, vals, height=h * 0.9, color=color, label=lab, zorder=2)
    labs = []
    for m, lab in STRICT:
        r = C.one(LRO, _crit("Ecoli_MG1655_Strict_256", "All_Strata", "Unfiltered", m))
        labs.append(f"{lab}\npool n = {int(C.num(r['n_cov']))}")
    n = int(C.get(LRO, _crit("Ecoli_MG1655_Strict_256", "All_Strata", "Unfiltered", STRICT[-1][0]), "n_total"))
    ax.set_yticks(y)
    ax.set_yticklabels(labs)
    ax.tick_params(axis="y", labelsize=C.BASE - 2)
    ax.set_xlim(0, 1)
    ax.set_ylim(y[-1] - 0.55, y[0] + 0.55)
    ax.set_xlabel(f"Recovery of {n} held-out $\\it{{E.\\ coli}}$ K-12 reactions")
    C.hgrid(ax, "x")
    C.top_legend(ax, ncol=2)
    C.title(ax, "Strict leave-reaction-out benchmark", pad=48)


def panel_b(ax):
    C.panel("fig3", "b")
    x = np.arange(len(STRATA))
    w = 0.8 / len(STRATA_METHODS)
    for k, (m, lab, color) in enumerate(STRATA_METHODS):
        vals = [C.get(LRO, _crit("Ecoli_MG1655_Full_414", s, "Unfiltered", m), "hit1_all") for s, _ in STRATA]
        ax.bar(x + (k - (len(STRATA_METHODS) - 1) / 2) * w, vals, width=w * 0.9, color=color, label=lab, zorder=2)
    labs = []
    for s, lab in STRATA:
        n = int(C.get(LRO, _crit("Ecoli_MG1655_Full_414", s, "Unfiltered", STRATA_METHODS[-1][0]), "n_total"))
        labs.append(f"{lab}\nn = {n}")
    ax.set_xticks(x)
    ax.set_xticklabels(labs)
    ax.set_ylim(0, 1)
    ax.set_ylabel("Top-1 recovery, all reactions")
    C.hgrid(ax)
    C.top_legend(ax, ncol=2)
    n = int(C.get(LRO, _crit("Ecoli_MG1655_Full_414", "All_Strata", "Unfiltered", STRATA_METHODS[-1][0]), "n_total"))
    C.title(ax, f"Role-stratified recovery ({n} $\\it{{E.\\ coli}}$ reactions)", pad=48)


def panel_c(ax):
    C.panel("fig3", "c")
    ticks, labels = [], []
    yy = 0.0
    # 50 ambiguous multi-candidate LRO gaps (fractions = percentage / 100), drawn bottom-up
    rows = [(C.get(ML, {"method": m}, "lro50_hit1_pct", "/100") / 100, None, lab, C.GREY) for m, lab in SINGLE]
    mu, sd = C.pm(C.get_raw(SEED, {"seed_label": "MEAN_SD"}, "lro_hit1_pct", "mean +/- population s.d. over 3 seeds; /100"))
    rows.append((mu / 100, sd / 100, "Stage 3 guarded re-ranking (3 seeds)", C.BLUE4))
    for v, e, lab, color in rows[::-1]:
        ax.barh(yy, v, height=0.7, color=color, zorder=2, xerr=e if e else None,
                error_kw=dict(ecolor=C.INK, elinewidth=1.4, capsize=4))
        ax.text(v + (e or 0) + 0.015, yy, f"{v:.2f}" if not e else f"{v:.2f} ± {e:.2f}", va="center",
                fontsize=C.BASE - 2, color=C.INK2)
        ticks.append(yy)
        labels.append(lab)
        yy += 1
    ax.text(0.0, yy - 0.35, "50 ambiguous $\\it{E.\\ coli}$ gaps",
            fontsize=C.BASE - 1, color=C.INK, fontweight="bold", va="bottom").set_in_layout(False)
    yy += 1.0
    # B. subtilis transfer with parameters fixed on E. coli
    cohort = "Bsubtilis_168_LRO_144"
    meths = [("GSPA_Strict_Homology_Masked_TrackA", "Track A (masked)"), ("GSPA_Full_M4_RLGC_Coherence_PLM", "Full GSPA (Track B)")]
    metrics = [("hit1_cov", "Hit@1, candidate pool", C.BLUE2), ("hit1_all", "Hit@1, all reactions", C.BLUE4)]
    for m, lab in meths[::-1]:
        for k, (col, mlab, color) in enumerate(metrics):
            v = C.get(LRO, _crit(cohort, "All_Strata", "Unfiltered", m), col)
            ax.barh(yy + (k - 0.5) * 0.4, v, height=0.38, color=color, zorder=2, label=mlab if m == meths[-1][0] else None)
            ax.text(v + 0.015, yy + (k - 0.5) * 0.4, f"{v:.2f}", va="center", fontsize=C.BASE - 2.5, color=C.INK2)
        ticks.append(yy)
        labels.append(lab)
        yy += 1
    n = int(C.get(LRO, _crit(cohort, "All_Strata", "Unfiltered", meths[-1][0]), "n_total"))
    ax.text(0.0, yy - 0.35, f"$\\it{{B.\\ subtilis}}$ 168 transfer ({n} reactions)", fontsize=C.BASE - 1,
            color=C.INK, fontweight="bold", va="bottom").set_in_layout(False)
    ax.set_yticks(ticks)
    ax.set_yticklabels(labels)
    ax.set_ylim(-0.6, yy + 0.4)
    ax.set_xlim(0, 1.1)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_xlabel("Hit@1")
    C.hgrid(ax, "x")
    h, l = ax.get_legend_handles_labels()
    C.top_legend(ax, ncol=2, handles=h[::-1], labels=l[::-1])
    C.title(ax, "Transfer and ambiguous-gap re-ranking", pad=28)


def panel_d(sf):
    C.panel("fig3", "d")
    min_mags = 4  # phyla with at least four MAGs are shown individually
    groups = [r["phylum_or_group"] for r in C.load(MAG)
              if r["phylum_or_group"] != "ALL" and C.num(r["n_mags"]) >= min_mags]
    per_mag = C.load(MAGP)
    axl, axr = sf.subplots(1, 2, sharey=True, gridspec_kw={"width_ratios": [1.25, 1], "wspace": 0.06})
    rows = ["ALL"] + groups
    y = np.arange(len(rows))[::-1].astype(float)
    y[0] += 0.5  # gap between the all-MAG row and the phyla
    labels = []
    rng = np.random.default_rng(20261019)  # fixed jitter for the per-MAG points
    for yi, g in zip(y, rows):
        n = int(C.get(MAG, {"phylum_or_group": g}, "n_mags"))
        mean = C.get(MAG, {"phylum_or_group": g}, "mean_fills_per_mag_q050")
        pts = [C.num(r["new_fills_q050"]) for r in per_mag if g == "ALL" or r["phylum"] == g]
        C.note(MAGP, "all MAGs" if g == "ALL" else f"phylum={g}", "new_fills_q050", "per-MAG points", len(pts))
        col = C.BLUE4 if g == "ALL" else C.BLUE2
        axl.barh(yi, mean, height=0.7, color=col, zorder=2)
        axl.scatter(pts, yi + rng.uniform(-0.22, 0.22, len(pts)), s=9, color=C.INK, alpha=0.55, lw=0, zorder=3)
        before = C.get(MAG, {"phylum_or_group": g}, "pathway_completeness_before")
        after = C.get(MAG, {"phylum_or_group": g}, "pathway_completeness_after")
        axr.plot([before, after], [yi, yi], color=C.AXIS, lw=2.0, zorder=1, solid_capstyle="round")
        axr.plot(before, yi, "o", color=C.GREY, ms=6.5, mec="white", mew=0.8, zorder=3)
        axr.plot(after, yi, "o", color=col, ms=6.5, mec="white", mew=0.8, zorder=3)
        labels.append(f"All {n} MAGs" if g == "ALL" else f"{g} ({n})")
    axl.set_yticks(y)
    axl.set_yticklabels(labels, fontsize=C.BASE - 3)
    axl.get_yticklabels()[0].set_fontweight("bold")
    axl.get_yticklabels()[0].set_color(C.INK)
    axl.set_ylim(y[-1] - 0.6, y[0] + 0.6)
    axl.set_xlabel("Gap fills per MAG\n(posterior ≥ 0.50)", fontsize=C.BASE - 1)
    axr.set_xlabel("MetaCyc pathway\ncompleteness", fontsize=C.BASE - 1)
    axr.set_xlim(0, axr.get_xlim()[1] * 1.06)
    for ax in (axl, axr):
        C.hgrid(ax, "x")
        ax.tick_params(axis="y", length=0)
    h = [plt_bar_handle(C.BLUE2, "Mean fills per MAG"),
         Line2D([], [], ls="none", marker="o", color=C.INK, alpha=0.55, ms=4, label="Individual MAGs"),
         Line2D([], [], ls="none", marker="o", color=C.GREY, ms=6.5, label="Before"),
         Line2D([], [], ls="none", marker="o", color=C.BLUE2, ms=6.5, label="After GSPA")]
    C.top_legend(axl, ncol=4, handles=h, fontsize=C.BASE - 2.5)
    nall = int(C.get(MAG, {"phylum_or_group": "ALL"}, "n_mags"))
    C.title(axl, f"Pathway gap fills in {nall} Empty Quarter MAGs", pad=26)


def plt_bar_handle(color, label):
    from matplotlib.patches import Patch
    return Patch(facecolor=color, label=label)


def generate():
    C.apply_style()
    fig, sf = C.figure()
    panel_a(sf[0].subplots())
    panel_b(sf[1].subplots())
    panel_c(sf[2].subplots())
    panel_d(sf[3])
    C.letters(sf)
    C.save(fig, "fig3_pillar3_metabolism")
