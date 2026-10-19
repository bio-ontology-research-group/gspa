"""Figure 2: Stage 1 prior and Stage 2A ablation across reference proteomes.

a  results/pillar1_2_ablation_summary.tsv          config=A1..A9 (nine-genome means)
b  results/pillar1_2_ablation_per_genome.tsv       config=A1_noisy_or vs A9_stage2_relax, composite_score
c  results/pillar1_2_perturbation_designs.tsv      jar=gspa-unified.jar; claims_ids=normalised; design=v1|v2
d  results/pillar1_suite2_panel13_21_sota.tsv      panel=Panel13_Prok13 (six methods of tab:ablation)
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

import common as C

ABL = "pillar1_2_ablation_summary.tsv"
ABG = "pillar1_2_ablation_per_genome.tsv"
PERT = "pillar1_2_perturbation_designs.tsv"
P13 = "pillar1_suite2_panel13_21_sota.tsv"

CONFIGS = [
    ("A1_noisy_or", "Stage 0 pooling"),
    ("A2_essentiality", "+ essentiality ($B_1$)"),
    ("A3_coherence", "+ coherence ($B_2$)"),
    ("A4_gap_filling", "+ gap filling ($B_5$)"),
    ("A5_genomic_context", "+ operon context ($B_6$)"),
    ("A6_consistency_anatomical", "+ taxon, anatomy ($B_3, B_4$)"),
    ("A7_stage1_all", "All Stage 1 priors"),
    ("A8_stage2_remove", "+ Stage 2A deletion"),
    ("A9_stage2_relax", "+ Stage 2A relaxation"),
]
METRICS = [
    ("composite_score", "Composite audit score", C.BLUE3, "o"),
    ("completeness", "Essential completeness", C.ORANGE, "s"),
    ("complex_coherence", "Complex coherence", C.AQUA, "D"),
    ("pathway_coherence", "Pathway coherence", C.YELLOW, "^"),
]
PERT_CONFIGS = [
    ("A1_noisy_or", "Stage 0"),
    ("A6_consistency_anatomical", "+ $B_3, B_4$"),
    ("A7_stage1_all", "Stage 1"),
    ("A8_stage2_remove", "2A del."),
    ("A9_stage2_relax", "2A rel."),
]
PANEL13 = [  # (method in TSV, label as in tab:ablation)
    ("NCBI PGAP (Li et al. 2021)", "NCBI PGAP"),
    ("HMMER 3.4 / Pfam-A + pfam2go", "HMMER3 Pfam-A + pfam2go"),
    ("InterProScan 5.72-103.0", "InterProScan 5.72"),
    ("DIAMOND-BLASTp", "DIAMOND-BLASTp vs Swiss-Prot"),
    ("metagenomic-deepFRI (Bezshapkin et al. 2026)", "metagenomic-deepFRI"),
    ("GSPA Stage 1+2A (3-Channel D+P+IPR)", "GSPA Stage 1+2A, 3 channels"),
]


def panel_a(sf):
    C.panel("fig2", "a")
    axes = sf.subplots(len(METRICS), 1, sharex=True, gridspec_kw={"hspace": 0.12})
    x = np.arange(len(CONFIGS))
    for ax, (col, lab, color, mk) in zip(axes, METRICS):
        ys = [C.get(ABL, {"config": cfg}, col) for cfg, _ in CONFIGS]
        ax.axhline(ys[0], color=C.MUTED, lw=1.0, ls=(0, (2, 2)), zorder=1)
        ax.plot(x, ys, ls="none", marker="o", ms=8, color=C.BLUE3, mec="white", mew=1.0, zorder=3)
        lo, hi = min(ys), max(ys)
        pad = max((hi - lo) * 0.25, 0.01)
        ax.set_ylim(lo - pad, hi + pad)
        ax.yaxis.set_major_locator(plt.MaxNLocator(2))
        for xv in (0.5, 5.5, 6.5):
            ax.axvline(xv, color=C.GRID, lw=1.0, zorder=0)
        ax.set_ylabel(lab.replace(" ", "\n", 1), rotation=0, ha="right", va="center", fontsize=C.BASE - 1,
                      labelpad=8)
        ax.tick_params(axis="x", length=0)
    axes[-1].set_xticks(x)
    axes[-1].set_xticklabels([lab for _, lab in CONFIGS], rotation=35, ha="right", rotation_mode="anchor")
    axes[-1].set_xlim(-0.5, len(CONFIGS) - 0.5)
    C.title(axes[0], "Structural audit by configuration")
    axes[-1].set_xlabel("Mean over nine genomes; dotted line marks the Stage 0 level", fontsize=C.BASE - 1.5,
                        color=C.INK2)


def panel_b(ax):
    C.panel("fig2", "b")
    genomes = sorted({r["genome"] for r in C.load(ABG)})
    data = []
    for g in genomes:
        a1 = C.get(ABG, {"genome": g, "config": "A1_noisy_or"}, "composite_score")
        a9 = C.get(ABG, {"genome": g, "config": "A9_stage2_relax"}, "composite_score")
        data.append((g, a1, a9))
    data.sort(key=lambda t: t[2] - t[1])
    y = np.arange(len(data))
    for i, (g, a1, a9) in enumerate(data):
        ax.plot([a1, a9], [i, i], color=C.AXIS, lw=2.2, zorder=1, solid_capstyle="round")
        ax.plot(a1, i, "o", color=C.GREY, ms=8.5, mec="white", mew=1.2, zorder=3)
        ax.plot(a9, i, "o", color=C.BLUE4, ms=8.5, mec="white", mew=1.2, zorder=3)
        ax.text(max(a1, a9) + 0.006, i, f"+{a9 - a1:.3f}" if a9 >= a1 else f"{a9 - a1:.3f}",
                va="center", ha="left", fontsize=C.BASE - 1.5, color=C.INK2)
    ax.set_yticks(y)
    ax.set_yticklabels([C.genome_label(g) for g, _, _ in data])
    lo = min(min(a1, a9) for _, a1, a9 in data)
    hi = max(max(a1, a9) for _, a1, a9 in data)
    ax.set_xlim(lo - 0.01, hi + 0.04)
    ax.set_ylim(-0.6, len(data) - 0.4)
    ax.set_xlabel("Composite audit score")
    C.hgrid(ax, "x")
    h = [Line2D([], [], ls="none", marker="o", color=C.GREY, ms=8.5, label="Stage 0 pooling"),
         Line2D([], [], ls="none", marker="o", color=C.BLUE4, ms=8.5, label="Stage 1 + 2A relaxation")]
    C.top_legend(ax, ncol=2, handles=h)
    C.title(ax, "Per-genome audit gain", pad=30)


def panel_c(ax):
    C.panel("fig2", "c")
    genomes = sorted({r["genome"] for r in C.load(PERT)})
    order = ["ecoli", "ecolo157", "bsubtilis", "paeruginosa", "mtb", "hpylori", "mgenitalium", "synechocystis",
             "mjannaschii"]
    genomes = [g for g in order if g in genomes] + [g for g in genomes if g not in order]
    designs = ["v1", "v2"]
    mat = np.zeros((len(genomes) + 1, len(designs) * len(PERT_CONFIGS)))
    for di, d in enumerate(designs):
        for ci, (cfg, _) in enumerate(PERT_CONFIGS):
            col = di * len(PERT_CONFIGS) + ci
            for gi, g in enumerate(genomes):
                mat[gi, col] = C.get(PERT, {"design": d, "jar": "gspa-unified.jar", "claims_ids": "normalised",
                                            "genome": g, "config": cfg}, "surviving_violations")
            mat[-1, col] = mat[:-1, col].sum()
            C.note(PERT, f"design={d}; jar=gspa-unified.jar; claims_ids=normalised; config={cfg}",
                   "surviving_violations", "sum over nine genomes", mat[-1, col])
    inj = [C.get(PERT, {"design": "v2", "jar": "gspa-unified.jar", "claims_ids": "normalised", "genome": g,
                        "config": "A1_noisy_or"}, "injected_pairs") for g in genomes]
    vmax = mat[:-1].max()
    shown = mat.copy()
    shown[-1] = np.nan
    ax.imshow(shown, cmap=C.SEQ_ORANGE, vmin=0, vmax=vmax, aspect="auto")
    for (i, j), v in np.ndenumerate(mat):
        if i == len(genomes):
            ax.text(j, i, f"{int(v)}", ha="center", va="center", fontsize=C.BASE - 1, color=C.INK, fontweight="bold")
        else:
            rgba = C.SEQ_ORANGE(v / vmax)
            ax.text(j, i, f"{int(v)}", ha="center", va="center", fontsize=C.BASE - 1.5, color=C.text_color_for(rgba))
    ax.set_xticks(range(mat.shape[1]))
    ax.set_xticklabels([lab for _ in designs for _, lab in PERT_CONFIGS], rotation=40, ha="right",
                       rotation_mode="anchor")
    ax.set_yticks(range(len(genomes) + 1))
    ax.set_yticklabels([C.genome_label(g) for g in genomes] + [f"Total (of {int(sum(inj))})"])
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.axvline(len(PERT_CONFIGS) - 0.5, color="white", lw=4)
    ax.axhline(len(genomes) - 0.5, color="white", lw=4)
    ax.set_xticks(np.arange(-0.5, mat.shape[1], 1), minor=True)
    ax.set_yticks(np.arange(-0.5, mat.shape[0], 1), minor=True)
    ax.grid(which="minor", color="white", lw=1.5)
    ax.tick_params(which="minor", length=0)
    for di, d in enumerate(designs):
        ax.text(di * len(PERT_CONFIGS) + (len(PERT_CONFIGS) - 1) / 2, -0.75, f"Injection design {d}",
                ha="center", va="bottom", fontsize=C.BASE - 0.5, color=C.INK, fontweight="bold")
    C.title(ax, "Surviving injected cross-clade claims", pad=26)


def panel_d(sf):
    C.panel("fig2", "d")
    gs = sf.add_gridspec(1, 2, wspace=0.04)
    fig = sf
    facets = [("Experimental annotations", "go_fmax_exp_prop", "go_fmax_cafa_exp_prop"),
              ("All GO annotations", "go_fmax_micro_all", "go_fmax_cafa_all")]
    axes = []
    y = np.arange(len(PANEL13))[::-1]
    for k, (title, micro, cafa) in enumerate(facets):
        ax = fig.add_subplot(gs[0, k])
        axes.append(ax)
        vals = []
        for (meth, lab), yi in zip(PANEL13, y):
            crit = {"panel": "Panel13_Prok13", "method": meth}
            vm, vc = C.get(P13, crit, micro), C.get(P13, crit, cafa)
            vals += [vm, vc]
            is_g = meth.startswith("GSPA")
            col = C.BLUE4 if is_g else C.GREY_D
            ax.plot(vm, yi + 0.17, "o", color=col, ms=8.5, mec="white", mew=1.0, zorder=3)
            ax.plot(vc, yi - 0.17, "D", color=col, ms=7.5, mfc="white", mec=col, mew=1.8, zorder=3)
            if is_g:
                ax.axhspan(yi - 0.45, yi + 0.45, color="#eaf2fc", zorder=0)
        lo, hi = min(vals), max(vals)
        pad = (hi - lo) * 0.12
        ax.set_xlim(lo - pad, hi + pad)
        ax.set_ylim(-0.6, len(PANEL13) - 0.4)
        ax.set_yticks(y)
        if k == 0:
            ax.set_yticklabels([lab for _, lab in PANEL13])
            for t in ax.get_yticklabels():
                if t.get_text().startswith("GSPA"):
                    t.set_fontweight("bold")
                    t.set_color(C.INK)
        else:
            ax.set_yticklabels([])
            ax.tick_params(axis="y", length=0)
        C.hgrid(ax, "x")
        ax.set_xlabel(r"$F_{\max}$")
        ax.text(0.0, 1.01, title, transform=ax.transAxes, fontsize=C.BASE - 0.5, color=C.INK2, va="bottom")
    h = [Line2D([], [], ls="none", marker="o", color=C.GREY_D, ms=8, label=r"Micro $F_{\max}$"),
         Line2D([], [], ls="none", marker="D", mfc="white", mec=C.GREY_D, mew=1.8, ms=7, label=r"CAFA $F_{\max}$")]
    C.top_legend(axes[0], ncol=2, handles=h, y=1.09)
    C.title(axes[0], "13-genome prokaryote panel", pad=48)
    return axes[0]


def generate():
    C.apply_style()
    fig, sf = C.figure()
    panel_a(sf[0])
    panel_b(sf[1].subplots())
    panel_c(sf[2].subplots())
    panel_d(sf[3])
    C.letters(sf)
    C.save(fig, "fig2_pillar1_2_ablation")
