"""Figure 4: heteromeric complex completion and eukaryotic anatomical constraint enforcement.

a  results/pillar3_complex_lso.tsv        genome=ALL (Weighted Mean); decoy_design in {unrelated_loci, within_operon};
                                          the seven method_id rows of tab:complex_lso; column hit_at_1
b  results/pillar3_complex_lso.tsv        decoy_design=unrelated_loci; per genome; f1_at_05
c  results/pillar3_eukaryote_anatomy.tsv  input_source=real; n_violations_e1..e5
d  results/pillar3_relax_vs_remove.tsv    input_source=real; scope=injected_cohort_all_proteins; cafa_fmax for
                                          E1_noisy_or / E4_stage2_remove / E5_stage2_relax
"""

from __future__ import annotations

import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

import common as C

LSO = "pillar3_complex_lso.tsv"
ANAT = "pillar3_eukaryote_anatomy.tsv"
RVR = "pillar3_relax_vs_remove.tsv"
ALL = "ALL (Weighted Mean)"

LSO_METHODS = [  # (method_id, label as in tab:complex_lso, colour)
    ("input_attenuated", "Attenuated input scores (no priors)", C.GREY),
    ("real_diamond_channel", "Real DIAMOND channel", C.ORANGE),
    ("real_interproscan_channel", "Real InterProScan channel", C.YELLOW),
    ("operon_only", "Operon prior only ($B_6$)", C.BLUE1),
    ("coherence", "Complex coherence prior only ($B_2$)", C.BLUE2),
    ("coherence_operon", "Coherence + operon priors ($B_2, B_6$)", C.BLUE3),
    ("full_stage3", "Full GSPA (Stages 1, 2A, and 3)", C.BLUE4),
]
F1_METHODS = ["input_attenuated", "coherence", "coherence_operon", "full_stage3"]
SPECIES = ["scerevisiae", "celegans", "dmelanogaster", "ppatens", "athaliana", "mmusculus"]
VIOL = [
    ("n_violations_e1_noisy_or", "Stage 0 pooling", C.GREY),
    ("n_violations_e2_go_taxon", "Direct GO taxon rules", C.ORANGE),
    ("n_violations_e3_stage1_anat", "Stage 1 soft anatomical prior", C.BLUE1),
    ("n_violations_e4_stage2_remove", "Stage 2A deletion", C.BLUE3),
    ("n_violations_e5_stage2_relax", "Stage 2A relaxation", C.BLUE4),
]
FMAX = [
    ("E1_noisy_or", "Stage 0", C.GREY, "o"),
    ("E4_stage2_remove", "Stage 2A deletion", C.BLUE2, "s"),
    ("E5_stage2_relax", "Stage 2A relaxation", C.BLUE4, "D"),
]


def panel_a(ax):
    C.panel("fig4", "a")
    designs = [("unrelated_loci", "Unrelated-loci decoys", C.BLUE3), ("within_operon", "Within-operon decoys", C.ORANGE)]
    y = np.arange(len(LSO_METHODS))[::-1]
    h = 0.38
    for k, (d, lab, col) in enumerate(designs):
        vals = [C.get(LSO, {"decoy_design": d, "genome": ALL, "method_id": m}, "hit_at_1") for m, _, _ in LSO_METHODS]
        ax.barh(y + (0.5 - k) * h, vals, height=h * 0.92, color=col, label=lab, zorder=2)
    n = C.get(LSO, {"decoy_design": "unrelated_loci", "genome": ALL, "method_id": "full_stage3"}, "n_cases")
    # chance level: one true subunit among the six candidates of each trial (benchmark design, five decoys)
    ax.axvline(1 / 6, color=C.INK2, lw=1.2, ls=(0, (3, 2)), zorder=3)
    ax.text(1 / 6 + 0.01, y[-1] - 0.75, "chance (1 of 6)", fontsize=C.BASE - 2, color=C.INK2, va="center")
    ax.set_yticks(y)
    ax.set_yticklabels([lab for _, lab, _ in LSO_METHODS])
    ax.set_xlim(0, 1.0)
    ax.set_ylim(y[-1] - 1.1, y[0] + 0.6)
    ax.set_xlabel(f"Hit@1 (n = {int(n)} held-out subunits, six genomes)")
    C.hgrid(ax, "x")
    C.top_legend(ax, ncol=2)
    C.title(ax, "Leave-subunit-out recovery", pad=28)


def panel_b(ax):
    C.panel("fig4", "b")
    genomes = [g for g in ["ecoli", "ecolo157", "bsubtilis", "paeruginosa", "mtb", "synechocystis"]
               if C.rows(LSO, {"decoy_design": "unrelated_loci", "genome": g})]
    cols = {m: (lab, col) for m, lab, col in LSO_METHODS}
    marks = dict(zip(F1_METHODS, ["o", "s", "^", "D"]))
    y = np.arange(len(genomes))[::-1]
    for k, m in enumerate(F1_METHODS):
        vals = [C.get(LSO, {"decoy_design": "unrelated_loci", "genome": g, "method_id": m}, "f1_at_05") for g in genomes]
        off = (1.5 - k) * 0.19
        ax.plot(vals, y + off, ls="none", marker=marks[m], color=cols[m][1], ms=7.5, mec="white", mew=1.0,
                label=cols[m][0], zorder=3)
    for yi in y:
        ax.axhline(yi, color=C.GRID, lw=0.8, zorder=0)
    labs = []
    for g in genomes:
        n = int(C.get(LSO, {"decoy_design": "unrelated_loci", "genome": g, "method_id": "full_stage3"}, "n_cases"))
        labs.append(C.genome_label(g) + f" (n = {n})")
    ax.set_yticks(y)
    ax.set_yticklabels(labs)
    ax.set_xlim(0, 1)
    ax.set_ylim(y[-1] - 0.6, y[0] + 0.6)
    ax.set_xlabel("Subunit recovery $F_1$ at threshold 0.50")
    C.hgrid(ax, "x")
    C.top_legend(ax, ncol=2)
    C.title(ax, "Per-genome recovery (unrelated-loci decoys)", pad=48)


def panel_c(ax):
    C.panel("fig4", "c")
    series = list(VIOL)
    v4 = [C.get(ANAT, {"input_source": "real", "species": s}, VIOL[3][0]) for s in SPECIES]
    v5 = [C.get(ANAT, {"input_source": "real", "species": s}, VIOL[4][0]) for s in SPECIES]
    if v4 == v5:  # identical counts: one bar for Stage 2A
        series = VIOL[:3] + [(VIOL[3][0], "Stage 2A deletion or relaxation", C.BLUE4)]
    y = np.arange(len(SPECIES))[::-1]
    h = 0.8 / len(series)
    vmax = 0
    for k, (col, lab, color) in enumerate(series):
        vals = [C.get(ANAT, {"input_source": "real", "species": s}, col) for s in SPECIES]
        vmax = max(vmax, max(vals))
        ys = y + ((len(series) - 1) / 2 - k) * h
        ax.barh(ys, vals, height=h * 0.9, color=color, label=lab, zorder=2)
        for yi, v in zip(ys, vals):
            if v == 0:
                ax.text(4, yi, "0", ha="left", va="center", fontsize=C.BASE - 2.5, color=color)
    labs = []
    for s in SPECIES:
        inj = int(C.get(ANAT, {"input_source": "real", "species": s}, "n_cross_clade_injected"))
        labs.append(C.genome_label(s) + f" ({inj} inj.)")
    ax.set_yticks(y)
    ax.set_yticklabels(labs)
    ax.set_ylim(y[-1] - 0.55, y[0] + 0.55)
    ax.set_xlabel("Violations (posterior ≥ 0.50), injected cohort")
    C.hgrid(ax, "x")
    C.top_legend(ax, ncol=2)
    C.title(ax, "Taxonomic and anatomical violations", pad=48)


def panel_d(ax):
    C.panel("fig4", "d")
    order = SPECIES + ["MEAN"]
    y = np.arange(len(order))[::-1]
    allv = []
    means = {}
    for k, (mode, lab, color, mk) in enumerate(FMAX):
        vals = [C.get(RVR, {"input_source": "real", "species": s, "mode": mode, "scope": "injected_cohort_all_proteins"},
                      "cafa_fmax") for s in SPECIES]
        means[mode] = sum(vals) / len(vals)
        C.note(RVR, f"input_source=real; mode={mode}; scope=injected_cohort_all_proteins", "cafa_fmax",
               "mean over six species", means[mode])
        vals = vals + [means[mode]]
        allv += vals
        off = (1 - k) * 0.22
        ax.plot(vals, y + off, ls="none", marker=mk, color=color, ms=8, mec="white", mew=1.0, label=lab, zorder=3)
    for yi in y:
        ax.axhline(yi, color=C.GRID, lw=0.8, zorder=0)
    ax.set_yticks(y)
    ax.set_yticklabels([C.genome_label(s) for s in SPECIES] + ["Mean of six species"])
    ax.get_yticklabels()[-1].set_fontweight("bold")
    ax.axhline(y[-1] + 0.5, color=C.AXIS, lw=1.0)
    lo, hi = min(allv), max(allv)
    ax.set_xlim(lo - 0.03, hi + 0.03)
    ax.set_ylim(y[-1] - 0.6, y[0] + 0.6)
    ax.set_xlabel(r"Micro $F_{\max}$ against Swiss-Prot, all cohort proteins")
    C.hgrid(ax, "x")
    C.top_legend(ax, ncol=3)
    C.title(ax, "Deletion versus ancestor relaxation", pad=28)


def generate():
    C.apply_style()
    fig, sf = C.figure()
    panel_a(sf[0].subplots())
    panel_b(sf[1].subplots())
    panel_c(sf[2].subplots())
    panel_d(sf[3].subplots())
    C.letters(sf)
    C.save(fig, "fig4_pillar3_complex_anatomy")
