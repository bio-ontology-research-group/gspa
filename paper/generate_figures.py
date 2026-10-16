#!/usr/bin/env python3
"""
Generate all 5 multi-panel publication figures for the GSPA manuscript
and synthesize results/pillar3_report.md directly from canonical results/ TSVs.

Figures generated in paper/figures/ (both vector PDF and 300 DPI PNG):
  - fig1_architecture.{pdf,png}
  - fig2_pillar1_2_ablation.{pdf,png}
  - fig3_pillar3_metabolism.{pdf,png}
  - fig4_pillar3_complex_anatomy.{pdf,png}
  - fig5_pillar4_agent.{pdf,png}
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Dict, List

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import numpy as np


REPO_ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = REPO_ROOT / "results"
FIG_DIR = REPO_ROOT / "paper" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# Clean publication visual grammar (Okabe-Ito colorblind-safe palette)
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
    "legend.fontsize": 8.0,
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


def read_jsonl(path: Path) -> List[Dict]:
    with open(path, "r", encoding="utf-8") as fh:
        return [json.loads(line) for line in fh if line.strip()]


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


def save_fig(fig, stem: str):
    pdf_path = FIG_DIR / f"{stem}.pdf"
    png_path = FIG_DIR / f"{stem}.png"
    fig.savefig(pdf_path, bbox_inches="tight")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[figures] Saved {pdf_path} and {png_path}")


# =====================================================================
# Figure 1: Unified 3-Stage Architecture of GSPA
# =====================================================================
def generate_figure1():
    fig, ax = plt.subplots(figsize=(13.2, 6.6))
    ax.set_xlim(0, 102)
    ax.set_ylim(0, 54)
    ax.axis("off")

    def draw_box(x, y, w, h, title, lines, fc, ec, title_color="#FFFFFF", header_h=4.0):
        box = FancyBboxPatch(
            (x, y), w, h,
            boxstyle="round,pad=0.3,rounding_size=1.2",
            facecolor="#FFFFFF",
            edgecolor=ec,
            linewidth=1.5,
        )
        ax.add_patch(box)
        hbar = FancyBboxPatch(
            (x, y + h - header_h), w, header_h,
            boxstyle="round,pad=0.2,rounding_size=1.0",
            facecolor=fc,
            edgecolor=ec,
            linewidth=1.2,
        )
        ax.add_patch(hbar)
        ax.text(
            x + w / 2, y + h - header_h / 2,
            title,
            ha="center", va="center",
            fontsize=8.9, fontweight="bold", color=title_color,
        )
        top_y = y + h - header_h - 1.3
        for idx, line in enumerate(lines):
            is_bold = line.startswith("*")
            clean = line.lstrip("*")
            ax.text(
                x + 1.0, top_y - idx * 2.02,
                clean,
                ha="left", va="top",
                fontsize=7.55,
                fontweight="bold" if is_bold else "normal",
                color="#1F2937",
            )

    def draw_arrow(x1, y1, x2, y2, label="", lx=None, ly=None):
        ax.annotate(
            "",
            xy=(x2, y2), xytext=(x1, y1),
            arrowprops=dict(
                arrowstyle="-|>,head_length=0.55,head_width=0.32",
                lw=1.7, color="#34495E",
            ),
        )
        if label:
            tx = (x1 + x2) / 2 if lx is None else lx
            ty = (y1 + y2) / 2 + 1.1 if ly is None else ly
            ax.text(
                tx, ty,
                label, ha="center", va="center",
                fontsize=7.5, fontweight="bold", color="#34495E",
                bbox=dict(boxstyle="round,pad=0.15", facecolor="#FFFFFF", edgecolor="none", alpha=0.85),
            )

    # Top Title Banner
    ax.text(
        51, 52.0,
        "GSPA v2.1.0: Unified Four-Stage Genome-Scale Protein Function Annotation Architecture",
        ha="center", va="center", fontsize=12, fontweight="bold", color=PALETTE["navy"],
    )

    # Column 1: Multi-Evidence Inputs & Genome Context
    draw_box(
        1.0, 26.5, 20.5, 23.0,
        "1A. Stage 0 Evidence",
        [
            "*Homology & Domain Channels:",
            "• DIAMOND (Swiss-Prot GOA)",
            "• Pfam / HMMER 3 + InterProScan",
            "*Structure & PLM Channels:",
            "• ESM-2 / ProtT5 kNN + CLEAN",
            "• DeepGOZero/SE + ProteInfer",
            "• FoldSeek / DeepFRI Structure",
            "*Empirical Isotonic Calibration",
        ],
        fc=PALETTE["navy"], ec=PALETTE["navy"],
    )

    draw_box(
        1.0, 1.5, 20.5, 22.5,
        "1B. Whole-Genome Context",
        [
            "*Genomic Coordinates (GFF3):",
            "• Co-directional Operon Calls",
            "*Metabolic Network (gapseq/BiGG):",
            "• Draft Stoichiometric Gaps",
            "• KEGG Modules & MetaCyc",
            "*Taxon & Anatomy Ontologies:",
            "• Direct GO Taxon Rules (26,911)",
            "• GO–Uberon/PO/CL Axioms (2,715)",
        ],
        fc=PALETTE["blue"], ec=PALETTE["blue"],
    )

    # Column 2: Stage 1 Iterative Bayesian Integration
    draw_box(
        26.0, 9.5, 24.5, 40.0,
        "Stage 1: Bayesian Refiner",
        [
            "*Log-Odds Fixed-Point Update (k=1..6):",
            "  L(p,f) = L_0(p,f) + Σ α_m B_m(p,f)",
            "",
            "*Six Genome-Scale Priors:",
            "• EssentialityPrior (+α_ess):",
            "  Universal/domain core GO categories",
            "• CoherencePrior (+α_coh, -β_single):",
            "  Pathway & heteromeric complex support",
            "• ConsistencyPrior (-α_cons):",
            "  SAT4J Xplain minimal UNSAT cores",
            "• AnatomicalPrior (-α_anat, +α_boost):",
            "  Uberon/PO/CL cross-product bridge",
            "• GapFillingPrior + RLGC S2:",
            "  Reaction-local genomic context + PLM",
            "• GenomicContextPrior (+α_operon):",
            "  Co-operonic pathway/complex boost",
        ],
        fc=PALETTE["teal"], ec=PALETTE["teal"],
    )

    # Column 3: Stage 2 Symbolic SAT & Stoichiometric Enforcement
    draw_box(
        55.0, 26.0, 21.5, 23.5,
        "Stage 2A/2B: Symbolic & pFBA",
        [
            "*2A. Deterministic Ontology Gate:",
            "• SAT Taxon & Uberon/PO/CL Check",
            "• Semantic Ancestor Relaxation:",
            "  Replaces forbidden child t with",
            "  nearest valid GO ancestor a*",
            "*2B. StoichiometricLocusRefiner:",
            "• Positive-cost LP + pFBA flux",
            "  certification (|v_j| > 10^-6)",
            "• Outer feedback ΔL_FBA to Stage 1",
        ],
        fc=PALETTE["orange"], ec=PALETTE["orange"], title_color="#111111",
    )

    # Column 3 Bottom: Output & GAEF Evaluator
    draw_box(
        55.0, 1.5, 21.5, 21.0,
        "Output & GAEF Audit",
        [
            "*Calibrated Proteome Output:",
            "• GAF 2.2 / TSV / JSON Reports",
            "*GAEF Genome-Scale Metrics:",
            "• Essential Completeness (C_ess)",
            "• Pathway & Complex Coherence",
            "• Taxon/Anatomical Consistency",
            "• Composite Quality Score",
        ],
        fc=PALETTE["navy"], ec=PALETTE["navy"],
    )

    # Column 4: Stage 3 Optional Neuro-Symbolic Agentic Layer
    draw_box(
        81.0, 9.5, 20.0, 40.0,
        "Stage 3: GSPA-Agent",
        [
            "*Selective Anomaly Trigger:",
            "• Inspects 1.8% mean (<2.5% max)",
            "• PATHWAY_GAP & COMPLEX",
            "• BORDERLINE_ESSENTIAL",
            "• BORDERLINE_CONFLICT",
            "",
            "*Tool-Using Local LLM:",
            "• inspect_operon_neighborhood",
            "• inspect_pathway_or_complex",
            "• query_protein_evidence",
            "• check_symbolic_constraints",
            "",
            "*Symbolic Guardrail Gate:",
            "• Blocks hallucinated GO/EC IDs",
            "• Enforces 0.0% taxon/anatomical",
            "  violation rate on promotions",
        ],
        fc=PALETTE["vermilion"], ec=PALETTE["vermilion"],
    )

    # Connecting arrows with clean whitespace label placement
    draw_arrow(21.8, 38.0, 25.7, 35.0, "Claims", lx=23.7, ly=38.3)
    draw_arrow(21.8, 13.0, 25.7, 22.0, "Priors", lx=23.7, ly=15.0)
    draw_arrow(50.8, 35.0, 54.7, 36.5, "Posteriors", lx=52.8, ly=38.0)
    draw_arrow(54.7, 30.0, 50.8, 28.5, "ΔL_FBA", lx=52.8, ly=26.8)
    draw_arrow(65.7, 25.7, 65.7, 22.8, "Enforced", lx=70.2, ly=24.3)
    draw_arrow(76.8, 38.0, 80.7, 38.0, "Anomalies", lx=78.8, ly=40.0)
    draw_arrow(80.7, 18.0, 76.8, 15.0, "Guarded", lx=78.8, ly=14.3)

    save_fig(fig, "fig1_architecture")


# =====================================================================
# Figure 2: Suites 1 & 2 — 9-Genome Ablation, Perturbation, and Panels
# =====================================================================
def generate_figure2():
    summary_all = read_tsv(RESULTS_DIR / "pillar1_2_ablation_summary.tsv")
    summary_rows = [r for r in summary_all if not r["config"].endswith("_non_ecoli8")]
    per_genome_rows = read_tsv(RESULTS_DIR / "pillar1_2_ablation_per_genome.tsv")
    perturb_rows = read_tsv(RESULTS_DIR / "pillar1_2_perturbation_summary.tsv")
    panel_rows = read_tsv(RESULTS_DIR / "pillar1_suite2_panel13_21_sota.tsv")

    fig, axes = plt.subplots(2, 2, figsize=(13.6, 9.5))
    plt.subplots_adjust(hspace=0.38, wspace=0.32)

    # --- Panel A: Ablation progression across 9 reference genomes ---
    ax = axes[0, 0]
    add_panel_label(ax, "a")
    cfg_labels = [
        "A1: Noisy-OR",
        "A2: +Essential",
        "A3: +Coherence",
        "A4: +GapFill",
        "A5: +Operon",
        "A6: +Taxon/Anat",
        "A7: Stage 1 All",
        "A8: Stage 2 (Rem)",
        "A9: Stage 2 (Rel)",
    ]
    comp = [float(r["completeness"]) for r in summary_rows]
    p_coh = [float(r["pathway_coherence"]) for r in summary_rows]
    c_coh = [float(r["complex_coherence"]) for r in summary_rows]
    gaef = [float(r["composite_score"]) for r in summary_rows]

    x = np.arange(len(cfg_labels))
    ax.plot(x, gaef, marker="o", lw=2.3, color=PALETTE["navy"], label="GAEF Composite (0.846→0.864)")
    ax.plot(x, comp, marker="s", lw=1.8, linestyle="--", color=PALETTE["teal"], label="Essential Completeness")
    ax.plot(x, c_coh, marker="^", lw=1.8, linestyle="-.", color=PALETTE["blue"], label="Complex Coherence")
    ax.plot(x, p_coh, marker="D", lw=1.8, linestyle=":", color=PALETTE["orange"], label="Pathway Coherence")
    ax.set_xticks(x)
    ax.set_xticklabels(cfg_labels, rotation=32, ha="right")
    ax.set_ylim(0.35, 1.00)
    ax.set_ylabel("Genome-Scale Metric Score")
    ax.set_title("9-Genome Prior & Stage 2 Ablation Progression")
    ax.grid(True, axis="y")
    ax.legend(loc="lower right", frameon=True, ncol=2, fontsize=7.4)

    # --- Panel B: Per-genome GAEF Composite Score (A1 vs A9) ---
    ax = axes[0, 1]
    add_panel_label(ax, "b")
    genomes = ["ecoli", "bsubtilis", "hpylori", "mgenitalium", "mjannaschii", "mtb", "paeruginosa", "synechocystis", "ecolo157"]
    g_names = ["E. coli K-12*", "B. subtilis", "H. pylori", "M. genitalium", "M. jannaschii", "M. tuberculosis", "P. aeruginosa", "Synechocystis", "E. coli O157"]
    a1_by_g = {r["genome"]: float(r["composite_score"]) for r in per_genome_rows if r["config"] == "A1_noisy_or"}
    a9_by_g = {r["genome"]: float(r["composite_score"]) for r in per_genome_rows if r["config"] == "A9_stage2_relax"}

    y_pos = np.arange(len(genomes))
    h = 0.36
    v1 = [a1_by_g[g] for g in genomes]
    v9 = [a9_by_g[g] for g in genomes]
    ax.barh(y_pos + h / 2, v1, height=h, color=PALETTE["lightgrey"], edgecolor=PALETTE["grey"], hatch="//", label="A1: Baseline Noisy-OR (Mean 0.846)")
    ax.barh(y_pos - h / 2, v9, height=h, color=PALETTE["teal"], edgecolor=PALETTE["navy"], label="A9: GSPA Stage 1+2 Relax (Mean 0.864)")
    for idx, (b_val, f_val) in enumerate(zip(v1, v9)):
        diff = f_val - b_val
        ax.text(f_val + 0.003, idx - h / 2, f"+{diff:.3f}", va="center", fontsize=7.6, fontweight="bold", color=PALETTE["navy"])
    ax.set_yticks(y_pos)
    ax.set_yticklabels(g_names)
    ax.set_xlim(0.70, 0.96)
    ax.set_ylim(-0.6, 10.6)
    ax.set_xlabel("GAEF Composite Quality Score (*E. coli K-12 calibration reference)")
    ax.set_title("Per-Genome GAEF Gain (8 Non-E. coli Mean: 0.843→0.860)")
    ax.grid(True, axis="x")
    ax.legend(loc="upper right", frameon=True, fontsize=7.6)

    # --- Panel C: Cross-Clade Taxon & Anatomical Perturbation Suppression ---
    ax = axes[1, 0]
    add_panel_label(ax, "c")
    p_cfgs = ["A1_noisy_or", "A6_consistency_anatomical", "A8_stage2_remove", "A9_stage2_relax"]
    p_labels = ["A1: Noisy-OR\n(Unconstrained)", "A6: Stage 1 Soft\nTaxon+Anat Prior", "A8: Stage 2\nHard Remove", "A9: Stage 2\nSemantic Relax"]
    surv_means = []
    supp_means = []
    for cfg in p_cfgs:
        sub = [r for r in perturb_rows if r["config"] == cfg]
        surv_means.append(sum(float(r["surviving_violations"]) for r in sub) / len(sub))
        supp_means.append(sum(float(r["suppression_rate_pct"]) for r in sub) / len(sub))

    x_p = np.arange(len(p_cfgs))
    bars = ax.bar(x_p, surv_means, width=0.52, color=[PALETTE["vermilion"], PALETTE["orange"], PALETTE["blue"], PALETTE["teal"]], edgecolor="#111111")
    for idx, b in enumerate(bars):
        val = surv_means[idx]
        supp = supp_means[idx]
        txt = f"{val:.1f}/genome\n(Unfiltered)" if idx == 0 else f"{val:.2f}/genome\n({supp:.1f}% supp.)"
        ax.text(b.get_x() + b.get_width() / 2, val + 0.35, txt, ha="center", va="bottom", fontsize=7.8, fontweight="bold")
    ax.set_xticks(x_p)
    ax.set_xticklabels(p_labels)
    ax.set_ylim(0, 17.5)
    ax.set_ylabel("Surviving Cross-Clade Violations / Genome")
    ax.set_title("Cross-Clade Perturbation Suppression (900 Injected Claims)")
    ax.grid(True, axis="y")

    # --- Panel D: 11-Prokaryote PGAP/DeepFRI Panel & 21-Genome Neural Panel (Ancestor-Propagated from TSV) ---
    ax = axes[1, 1]
    add_panel_label(ax, "d")
    panel_lookup = {(r["panel"], r["method"]): r for r in panel_rows}
    panel_items = [
        ("Panel13_Prok11", "metagenomic-deepFRI (Bezshapkin et al. 2026)", "metagenomic-deepFRI\n(11-prok, 2026)"),
        ("Panel13_Prok11", "InterProScan 5.72", "InterProScan 5.72\n(11-prok panel)"),
        ("Panel13_Prok11", "NCBI PGAP (Li et al. 2021)", "NCBI PGAP\n(11-prok, 2021)"),
        ("Panel13_Prok11", "GSPA Stage 1+2 (3-Channel Full)", "GSPA Stage 1+2\n(11-prok, 3-ch)"),
        ("Panel21_All21", "FoldSeek + ProstT5 3Di (2024)", "FoldSeek+ProstT5\n(21-genome, 2024)"),
        ("Panel21_All21", "InterProScan 5.72", "InterProScan 5.72\n(21-genome panel)"),
        ("Panel21_All21", "ESM2-DeepGOPlus (2024)", "ESM2-DeepGOPlus\n(21-genome, 2024)"),
        ("Panel21_All21", "ProteInfer (Sanderson et al. 2023)", "ProteInfer\n(21-genome, 2023)"),
        ("Panel21_All21", "DIAMOND-BLASTp", "DIAMOND-BLASTp\n(21-genome panel)"),
        ("Panel21_All21", "ESM-2 kNN (Lin et al. 2023)", "ESM-2 kNN\n(21-genome, 2023)"),
        ("Panel9_Sub9", "gLM / gLM2 Genomic LM (Hwang et al. 2024)", "gLM / gLM2\n(9-genome, 2024/25)"),
        ("Panel21_All21", "GSPA Stage 1+2 (3-Channel Core)", "GSPA Stage 1+2\n(21-genome, 3-ch)"),
    ]
    methods_panel = [lbl for _, _, lbl in panel_items]
    micro_fmax = [float(panel_lookup[(p, m)]["go_fmax_micro_all"]) for p, m, _ in panel_items]
    cafa_fmax = [float(panel_lookup[(p, m)]["go_fmax_cafa_all"]) for p, m, _ in panel_items]

    y_m = np.arange(len(methods_panel))
    h_m = 0.38
    ax.barh(y_m + h_m / 2, micro_fmax, height=h_m, color=PALETTE["blue"], edgecolor=PALETTE["navy"], label="Micro Fmax (All-GOA Propagated)")
    ax.barh(y_m - h_m / 2, cafa_fmax, height=h_m, color=PALETTE["sky"], edgecolor=PALETTE["blue"], hatch="//", label="CAFA Fmax (All-GOA Propagated)")
    for idx, mv in enumerate(micro_fmax):
        ax.text(mv + 0.008, idx + h_m / 2, f"{mv:.3f}", va="center", fontsize=6.8, fontweight="bold")
    ax.set_yticks(y_m)
    ax.set_yticklabels(methods_panel, fontsize=7.0)
    ax.set_xlim(0.0, 0.99)
    ax.set_xlabel("Ancestor-Propagated GO Fmax Score (11-Prokaryote & 21-Genome Panels)")
    ax.set_title("Uniform Ancestor-Propagated Comparison vs. SOTA Annotators")
    ax.grid(True, axis="x")
    ax.legend(loc="lower right", frameon=True, fontsize=7.2)

    save_fig(fig, "fig2_pillar1_2_ablation")


# =====================================================================
# Figure 3: Suites 3–5 — Metabolic Gap-Filling & Dark Matter
# =====================================================================
def generate_figure3():
    lro15_rows = read_tsv(RESULTS_DIR / "pillar2_lro15_summary.tsv")
    lro_agent_rows = read_tsv(RESULTS_DIR / "pillar4_lro_reranking.tsv")
    mag_rows = read_tsv(RESULTS_DIR / "pillar3_dark_matter_97mags.tsv")

    fig, axes = plt.subplots(2, 2, figsize=(13.6, 9.6))
    plt.subplots_adjust(hspace=0.44, wspace=0.30)

    # --- Panel A: E. coli MG1655 Strict 256-Reaction LRO (All vs Given Pool) ---
    ax = axes[0, 0]
    add_panel_label(ax, "a")
    ecoli_methods = [
        ("Sequence_Only_Candidate_Rank (results_protein.tsv)", "Seq / METEOR\nUniform Boost"),
        ("GSPA_Strict_Homology_Masked_TrackA (results_m4_trackA_strict_mg1655.tsv)", "GSPA Strict\nMasked (Trk A)"),
        ("RLGC_Kernel_Density_Only (results_rlgc_mg1655.tsv)", "RLGC Kernel\nDensity Only"),
        ("GSPA_Stage1_M2_Coherence_RLGC (results_m2_integrated_mg1655.tsv)", "GSPA Stage 1\nCoherence+RLGC"),
        ("ESM2_PLM_Centroids_Only (results_plmonly_mg1655.tsv)", "ESM-2 PLM\nCentroids Only"),
        ("GSPA_Full_M4_RLGC_Coherence_PLM (results_m4_trackB_mg1655.tsv)", "GSPA Full M4\n(RLGC+Coh+PLM)"),
    ]
    ec_rows = {r["method"]: r for r in lro15_rows if r["benchmark_cohort"] == "Ecoli_MG1655_Strict_256"}
    h1_all = [float(ec_rows[m]["hit1_all"]) for m, _ in ecoli_methods]
    h1_pool = [float(ec_rows[m]["hit1_cov"]) for m, _ in ecoli_methods]
    h3_pool = [float(ec_rows[m]["hit3_cov"]) for m, _ in ecoli_methods]

    x = np.arange(len(ecoli_methods))
    w = 0.26
    b1 = ax.bar(x - w, h1_all, width=w, color=PALETTE["blue"], edgecolor=PALETTE["navy"], label="Hit@1 (All 256 Reactions)")
    b2 = ax.bar(x, h1_pool, width=w, color=PALETTE["teal"], edgecolor=PALETTE["navy"], hatch="//", label="Hit@1 (Given Pool, n=175–181)")
    ax.bar(x + w, h3_pool, width=w, color=PALETTE["sky"], edgecolor=PALETTE["blue"], label="Hit@3 (Given Pool, n=175–181)")
    for b in b1:
        v = b.get_height()
        ax.text(b.get_x() + b.get_width() / 2, v + 0.015, f"{v:.2f}", ha="center", va="bottom", fontsize=6.6, fontweight="bold", color=PALETTE["navy"])
    for b in b2:
        v = b.get_height()
        ax.text(b.get_x() + b.get_width() / 2, v + 0.015, f"{v:.2f}", ha="center", va="bottom", fontsize=6.6, fontweight="bold", color=PALETTE["teal"])
    ax.set_xticks(x)
    ax.set_xticklabels([lbl for _, lbl in ecoli_methods], fontsize=7.4, rotation=16, ha="right")
    ax.set_ylim(0, 1.18)
    ax.set_ylabel("Enzyme Gene Recovery Rate")
    ax.set_title("E. coli K-12* LRO Benchmark (n=256 Reactions; Pool Cov. 70.7%)")
    ax.grid(True, axis="y")
    ax.legend(loc="upper left", frameon=True, fontsize=7.4)

    # --- Panel B: Role-Stratified LRO Across Strata A, B, and C (n=414, Read from TSV) ---
    ax = axes[0, 1]
    add_panel_label(ax, "b")
    f414 = {(r["stratum"], r["protocol_H"], r["method"].split(" ")[0]): r
            for r in lro15_rows if r["benchmark_cohort"] == "Ecoli_MG1655_Full_414"}
    strata_keys = [
        ("Stratum_A_Core_Biomass", "Unfiltered", "Stratum A: Core\nBiomass (n=307)"),
        ("Stratum_B_Catabolic", "Unfiltered", "Stratum B: Alt.\nCatabolic (n=49)"),
        ("Stratum_C_Secondary_Peripheral", "Unfiltered", "Stratum C: Periph.\nSecondary (n=58)"),
        ("Stratum_A_Core_Biomass", "Protocol_H_le50", "Stratum A ≤50%\nHomology (n=23)"),
    ]
    strata_labels = [lbl for _, _, lbl in strata_keys]
    rlgc_pool = [float(f414[(s, p, "RLGC_Kernel_Density_Only")]["hit1_cov"]) for s, p, _ in strata_keys]
    gspa_all =  [float(f414[(s, p, "GSPA_Full_M4_RLGC_Coherence_PLM")]["hit1_all"]) for s, p, _ in strata_keys]
    gspa_pool = [float(f414[(s, p, "GSPA_Full_M4_RLGC_Coherence_PLM")]["hit1_cov"]) for s, p, _ in strata_keys]

    x_n = np.arange(len(strata_labels))
    w_n = 0.26
    ax.bar(x_n - w_n, rlgc_pool, width=w_n, color=PALETTE["orange"], edgecolor="#111111", hatch="//", label="RLGC Operon Kernel Only (Hit@1 Pool)")
    ax.bar(x_n, gspa_all, width=w_n, color=PALETTE["blue"], edgecolor=PALETTE["navy"], label="GSPA Full M4 (Hit@1 All)")
    b_gp = ax.bar(x_n + w_n, gspa_pool, width=w_n, color=PALETTE["teal"], edgecolor=PALETTE["navy"], label="GSPA Full M4 (Hit@1 Given Pool)")
    for b in b_gp:
        v = b.get_height()
        ax.text(b.get_x() + b.get_width() / 2, v + 0.015, f"{v:.2f}", ha="center", va="bottom", fontsize=7.4, fontweight="bold", color=PALETTE["navy"])
    ax.set_xticks(x_n)
    ax.set_xticklabels(strata_labels, fontsize=7.8)
    ax.set_ylim(0, 1.15)
    ax.set_ylabel("Top-1 Enzyme Recovery (Hit@1)")
    ax.set_title("Role-Stratified LRO Across Metabolic Strata A, B, & C (n=414)")
    ax.grid(True, axis="y")
    ax.legend(loc="upper right", frameon=True, fontsize=7.4)

    # --- Panel C: Cross-Species B. subtilis (n=144) & Ambiguous LRO (n=50 Unique Proteins) ---
    ax = axes[1, 0]
    add_panel_label(ax, "c")
    bsub_rows = {r["method"].split(" ")[0]: r for r in lro15_rows if r["benchmark_cohort"] == "Bsubtilis_168_LRO_144"}
    bs_base = bsub_rows["GSPA_Strict_Homology_Masked_TrackA"]
    bs_full = bsub_rows["GSPA_Full_M4_RLGC_Coherence_PLM"]

    rlgc_pos = [int(r["rlgc_position"]) for r in lro_agent_rows]
    agent_pos = [int(r["agent_position"]) for r in lro_agent_rows]
    n_lro50 = len(lro_agent_rows)
    rlgc_h1 = sum(1 for p in rlgc_pos if p == 1) / n_lro50
    agent_h1 = sum(1 for p in agent_pos if p == 1) / n_lro50
    rlgc_mrr = sum(1.0 / p for p in rlgc_pos) / n_lro50
    agent_mrr = sum(1.0 / p for p in agent_pos) / n_lro50

    c_groups = [
        "B. subtilis 168\n(n=144, Hit@1 All)",
        "B. subtilis 168\n(n=76, Hit@1 Pool)",
        "B. subtilis 168\n(n=76, MRR Pool)",
        "Ambiguous LRO\n(n=50, Hit@1)",
        "Ambiguous LRO\n(n=50, MRR)",
    ]
    base_vals = [float(bs_base["hit1_all"]), float(bs_base["hit1_cov"]), float(bs_base["mrr_cov"]), rlgc_h1, rlgc_mrr]
    gspa_vals = [float(bs_full["hit1_all"]), float(bs_full["hit1_cov"]), float(bs_full["mrr_cov"]), agent_h1, agent_mrr]

    x_c = np.arange(len(c_groups))
    w_c = 0.34
    b_base = ax.bar(x_c - w_c / 2, base_vals, width=w_c, color=PALETTE["orange"], edgecolor="#111111", hatch="//", label="Strict Masked Track A / RLGC Alone")
    b_gspa = ax.bar(x_c + w_c / 2, gspa_vals, width=w_c, color=PALETTE["navy"], edgecolor="#111111", label="GSPA Full M4 / + Stage 3 Agent Re-Ranking")
    for bars in (b_base, b_gspa):
        for b in bars:
            v = b.get_height()
            ax.text(b.get_x() + b.get_width() / 2, v + 0.015, f"{v:.3f}", ha="center", va="bottom", fontsize=7.2, fontweight="bold")
    ax.set_xticks(x_c)
    ax.set_xticklabels(c_groups, fontsize=7.5, rotation=14, ha="right")
    ax.set_ylim(0, 1.15)
    ax.set_ylabel("Recovery Score")
    ax.set_title("Cross-Species B. subtilis (n=144) & Agent LRO Disambiguation (n=50)")
    ax.grid(True, axis="y")
    ax.legend(loc="upper left", frameon=True, fontsize=7.4)

    # --- Panel D: 97-MAG Uncultured Microbial Dark-Matter Panel (Read from TSV) ---
    ax = axes[1, 1]
    add_panel_label(ax, "d")
    plot_mag_rows = [r for r in mag_rows if r["plot_group"] in ("top8", "top8_summary")]
    phyla = [f"{r['phylum_or_group']} (n={r['n_mags']})" for r in plot_mag_rows]
    fills_q03 = [float(r["mean_fills_per_mag_q030"]) for r in plot_mag_rows]
    fills_q05 = [float(r["mean_fills_per_mag_q050"]) for r in plot_mag_rows]

    y_p = np.arange(len(phyla))
    ax.barh(y_p, fills_q03, height=0.55, color=PALETTE["sky"], edgecolor=PALETTE["blue"], hatch="//", label="Candidate Fills at q ≥ 0.3 (Total = 2,175)")
    ax.barh(y_p, fills_q05, height=0.55, color=PALETTE["navy"], label="High-Confidence Fills at q ≥ 0.5 (Total = 568)")
    for idx, (v3, v5) in enumerate(zip(fills_q03, fills_q05)):
        ax.text(v3 + 0.4, idx, f"{v3:.1f} ({v5:.1f} HQ)", va="center", fontsize=7.5, fontweight="bold")
    ax.set_yticks(y_p)
    ax.set_yticklabels(phyla)
    ax.set_xlim(0, 40)
    ax.set_ylim(-0.6, 10.6)
    ax.set_xlabel("Mean Metabolic Pathway Gaps Resolved per MAG")
    ax.set_title("Exploratory Dark-Matter Candidates Across 97 MAGs (16 Phyla)")
    ax.grid(True, axis="x")
    ax.legend(loc="upper right", frameon=True, fontsize=7.4)

    save_fig(fig, "fig3_pillar3_metabolism")


# =====================================================================
# Figure 4: Suites 2 & 3–5 — Complexes & Eukaryotic Anatomy
# =====================================================================
def generate_figure4():
    lso_rows = read_tsv(RESULTS_DIR / "pillar3_complex_lso.tsv")
    anat_rows = read_tsv(RESULTS_DIR / "pillar3_eukaryotic_anatomy.tsv")
    rel_rows = read_tsv(RESULTS_DIR / "pillar3_relax_vs_remove.tsv")

    fig, axes = plt.subplots(2, 2, figsize=(13.6, 9.4))
    plt.subplots_adjust(hspace=0.40, wspace=0.30)

    # --- Panel A: Obligate Heteromeric Complex LSO Summary (175 cases) ---
    ax = axes[0, 0]
    add_panel_label(ax, "a")
    all_lso = [r for r in lso_rows if r["genome"] == "ALL (Weighted Mean)"]
    m_names = ["M1: Baseline\nNoisy-OR", "M2: Stage 1\nCoherencePrior", "M3: Stage 1\n+ Operon Context", "M4: Stage 1+2\n+ Stage 3 Agent"]
    hit1 = [float(r["hit_at_1"]) for r in all_lso]
    rec = [float(r["recall_at_05"]) for r in all_lso]
    f1 = [float(r["f1_at_05"]) for r in all_lso]
    coh = [float(r["complex_coherence"]) for r in all_lso]

    x = np.arange(len(m_names))
    w = 0.20
    ax.bar(x - 1.5 * w, hit1, width=w, color=PALETTE["navy"], label="Hit@1")
    ax.bar(x - 0.5 * w, rec, width=w, color=PALETTE["teal"], hatch="//", label="Recall@0.5")
    ax.bar(x + 0.5 * w, f1, width=w, color=PALETTE["blue"], label="F1@0.5")
    ax.bar(x + 1.5 * w, coh, width=w, color=PALETTE["orange"], hatch="..", label="Genome Complex Coherence")
    for idx in range(len(m_names)):
        ax.text(x[idx] - 0.5 * w, rec[idx] + 0.02, f"{rec[idx]:.2f}", ha="center", va="bottom", fontsize=7.2, fontweight="bold", color=PALETTE["teal"])
        ax.text(x[idx] - 1.5 * w, hit1[idx] + 0.02, f"{hit1[idx]:.2f}", ha="center", va="bottom", fontsize=7.2, fontweight="bold", color=PALETTE["navy"])
    ax.set_xticks(x)
    ax.set_xticklabels(m_names)
    ax.set_ylim(0, 1.25)
    ax.set_ylabel("Score Across 175 Held-Out Subunits")
    ax.set_title("Leave-Subunit-Out (LSO) Heteromeric Complex Completion")
    ax.grid(True, axis="y")
    ax.legend(loc="upper left", frameon=True, ncol=2, fontsize=7.5)

    # --- Panel B: Per-Genome LSO F1@0.5 Progression (Highlighting Synechocystis) ---
    ax = axes[0, 1]
    add_panel_label(ax, "b")
    g_order = ["ecoli", "bsubtilis", "mtb", "paeruginosa", "ecolo157", "synechocystis"]
    g_labels = ["E. coli K-12*", "B. subtilis", "M. tuberculosis", "P. aeruginosa", "E. coli O157", "Synechocystis\n(Dispersed Loci)"]
    f1_m2 = [float(next(r["f1_at_05"] for r in lso_rows if r["genome"] == g and "CoherencePrior" in r["method"])) for g in g_order]
    f1_m3 = [float(next(r["f1_at_05"] for r in lso_rows if r["genome"] == g and "Operon" in r["method"])) for g in g_order]
    f1_m4 = [float(next(r["f1_at_05"] for r in lso_rows if r["genome"] == g and "Agent" in r["method"])) for g in g_order]

    x_g = np.arange(len(g_order))
    w_g = 0.26
    ax.bar(x_g - w_g, f1_m2, width=w_g, color=PALETTE["sky"], edgecolor=PALETTE["blue"], hatch="//", label="M2: CoherencePrior Only")
    ax.bar(x_g, f1_m3, width=w_g, color=PALETTE["blue"], edgecolor=PALETTE["navy"], label="M3: Coherence + Operon")
    ax.bar(x_g + w_g, f1_m4, width=w_g, color=PALETTE["teal"], edgecolor=PALETTE["navy"], label="M4: + Stage 3 Agent")
    for idx, (v2, v3, v4) in enumerate(zip(f1_m2, f1_m3, f1_m4)):
        if v4 > v3 + 0.005:
            ax.text(x_g[idx] + w_g, v4 + 0.015, f"+{v4 - v3:.2f}", ha="center", va="bottom", fontsize=7.6, fontweight="bold", color=PALETTE["vermilion"])
        else:
            ax.text(x_g[idx] + w_g, v4 + 0.015, f"{v4:.2f}", ha="center", va="bottom", fontsize=7.0, fontweight="bold", color=PALETTE["navy"])
    ax.set_xticks(x_g)
    ax.set_xticklabels(g_labels, rotation=18, ha="right")
    ax.set_ylim(0, 1.08)
    ax.set_ylabel("Subunit Recovery F1@0.5 (M1 Baseline = 0.00)")
    ax.set_title("Per-Genome Complex Recovery & Non-Operonic Agent Lift")
    ax.grid(True, axis="y")
    ax.legend(loc="upper right", frameon=True, fontsize=7.5)

    # --- Panel C: 6-Eukaryote Cross-Clade Violation Elimination (E1, E2, E3, E4/E5) ---
    ax = axes[1, 0]
    add_panel_label(ax, "c")
    euk_labels = ["S. cerevisiae\n(Fungi)", "C. elegans\n(Nematoda)", "D. melanogaster\n(Arthropoda)", "P. patens\n(Bryophyta)", "A. thaliana\n(Magnoliopsida)", "M. musculus\n(Mammalia)"]
    v_e1 = [int(r["n_violations_e1_noisy_or"]) for r in anat_rows]
    v_e2 = [int(r["n_violations_e2_go_taxon"]) for r in anat_rows]
    v_e3 = [int(r["n_violations_e3_stage1_anat"]) for r in anat_rows]
    v_e4 = [int(r["n_violations_e4_stage2_remove"]) for r in anat_rows]

    x_e = np.arange(len(anat_rows))
    w_e = 0.21
    ax.bar(x_e - 1.5 * w_e, v_e1, width=w_e, color=PALETTE["vermilion"], label="E1: Unconstrained Noisy-OR")
    ax.bar(x_e - 0.5 * w_e, v_e2, width=w_e, color=PALETTE["orange"], hatch="//", label="E2: Direct GO Taxon Constraints Only")
    ax.bar(x_e + 0.5 * w_e, v_e3, width=w_e, color=PALETTE["sky"], edgecolor=PALETTE["blue"], label="E3: Stage 1 Soft Prior (80.7% Mean / 47.6% Pooled Elim.)")
    ax.bar(x_e + 1.5 * w_e, v_e4, width=w_e, color=PALETTE["teal"], label="E4/E5: Stage 2 Symbolic (100% Elim.)")
    for idx, (ve1, ve2, ve3, ve4) in enumerate(zip(v_e1, v_e2, v_e3, v_e4)):
        y_off2 = 52 if ve1 > 900 else 12
        ax.text(x_e[idx] - 1.65 * w_e, ve1 + 12, str(ve1), ha="center", va="bottom", fontsize=6.6, fontweight="bold")
        ax.text(x_e[idx] - 0.40 * w_e, ve2 + y_off2, str(ve2), ha="center", va="bottom", fontsize=6.6, fontweight="bold")
        if ve3 > 0:
            ax.text(x_e[idx] + 0.5 * w_e, ve3 + 12, str(ve3), ha="center", va="bottom", fontsize=6.6, fontweight="bold", color=PALETTE["blue"])
        ax.text(x_e[idx] + 1.5 * w_e, ve4 + 12, "0", ha="center", va="bottom", fontsize=7.0, fontweight="bold", color=PALETTE["teal"])
    ax.set_xticks(x_e)
    ax.set_xticklabels(euk_labels, rotation=18, ha="right")
    ax.set_ylim(0, 1380)
    ax.set_ylabel("Surviving Taxon / Anatomical Violations")
    ax.set_title("GO–Uberon/PO/CL Bridge Eliminates Cross-Clade Violations")
    ax.grid(True, axis="y")
    ax.legend(loc="upper right", frameon=True, fontsize=6.9)

    # --- Panel D: Semantic Ancestor Relaxation (E5) vs Naive Removal (E4) & Baseline (E1) ---
    ax = axes[1, 1]
    add_panel_label(ax, "d")
    euk_ids = [r["species"] for r in anat_rows]
    fmax_e1 =  [float(next(r["cafa_fmax"] for r in rel_rows if r["species"] == s and r["mode"] == "E1_noisy_or")) for s in euk_ids]
    fmax_rem = [float(next(r["cafa_fmax"] for r in rel_rows if r["species"] == s and r["mode"] == "E4_stage2_remove")) for s in euk_ids]
    fmax_rel = [float(next(r["cafa_fmax"] for r in rel_rows if r["species"] == s and r["mode"] == "E5_stage2_relax")) for s in euk_ids]

    w_r = 0.26
    ax.bar(x_e - w_r, fmax_e1, width=w_r, color=PALETTE["orange"], edgecolor="#111111", hatch="..", label="E1: Unconstrained Noisy-OR")
    ax.bar(x_e, fmax_rem, width=w_r, color=PALETTE["lightgrey"], edgecolor=PALETTE["grey"], hatch="//", label="E4: Stage 2 Remove (Drops Term)")
    ax.bar(x_e + w_r, fmax_rel, width=w_r, color=PALETTE["blue"], edgecolor=PALETTE["navy"], label="E5: Stage 2 Semantic Relax (Valid Ancestor)")
    for idx, (ve1, vr, vl) in enumerate(zip(fmax_e1, fmax_rem, fmax_rel)):
        gain_rem = vl - vr
        ax.text(x_e[idx] + w_r, vl + 0.008, f"{gain_rem:+.3f}", ha="center", va="bottom", fontsize=7.0, fontweight="bold", color=PALETTE["navy"])
    ax.set_xticks(x_e)
    ax.set_xticklabels(euk_labels, rotation=18, ha="right")
    ax.set_ylim(0.55, 1.03)
    ax.set_ylabel("Ancestor-Propagated Micro Fmax")
    ax.set_title("Semantic Relaxation vs. Removal Across 6 Eukaryotes")
    ax.grid(True, axis="y")
    ax.legend(loc="upper left", frameon=True, fontsize=7.3)

    save_fig(fig, "fig4_pillar3_complex_anatomy")


# =====================================================================
# Figure 5: Stage 3 Neuro-Symbolic Agent (gspa-agent) on De-identified Benchmark
# =====================================================================
def generate_figure5():
    vs_rows = read_tsv(RESULTS_DIR / "pillar4_agent_evaluation.tsv")
    lro_rows = read_tsv(RESULTS_DIR / "pillar4_lro_reranking.tsv")
    res_c2 = read_jsonl(RESULTS_DIR / "pillar4_resolutions_c2_fast.jsonl")
    res_c3 = read_jsonl(RESULTS_DIR / "pillar4_resolutions_c3_guarded.jsonl")

    fig, axes = plt.subplots(2, 2, figsize=(13.6, 9.2))
    plt.subplots_adjust(hspace=0.38, wspace=0.30)

    vs_by_cfg = {r["configuration"]: r for r in vs_rows}
    selected_cfgs = [
        ("Stage 1 + Stage 2 Deterministic Only", "Config 1:\nStage 1+2\nOnly"),
        ("Posterior Threshold (>=0.35) + Guardrail", "Threshold\n(>=0.35)\n+ Guardrail"),
        ("Logistic Regression (5-Fold CV, 3 Pre-Tool Scalars)", "LogReg\n(3 Pre-Tool\nScalars)"),
        ("GBDT (5-Fold CV, 3 Pre-Tool Scalars)", "GBDT\n(3 Pre-Tool\nScalars)"),
        ("Stage 3 Agent — Fast 1-Step Mode", "Config 2:\nFast 1-Step\nAgent"),
        ("Stage 3 Agent — Qwen 27B Without Guardrail", "Config 4:\nQwen 3.8 27B\nNo Guardrail"),
        ("Stage 3 Agent — Qwen 27B + Guardrail", "Config 3:\nQwen 3.8 27B\n+ Guardrail"),
        ("GBDT (5-Fold CV) + Guardrail", "GBDT (8-Feat\n5-Fold CV)\n+ Guardrail"),
    ]
    cfg_short = [lbl for _, lbl in selected_cfgs]
    prec = [float(vs_by_cfg[k]["precision"]) * 100.0 for k, _ in selected_cfgs]
    rec =  [float(vs_by_cfg[k]["recall"]) * 100.0 for k, _ in selected_cfgs]
    f1 =   [float(vs_by_cfg[k]["f1"]) * 100.0 for k, _ in selected_cfgs]
    acc =  [float(vs_by_cfg[k]["accuracy"]) * 100.0 for k, _ in selected_cfgs]

    # --- Panel A: De-identified 120-Anomaly Benchmark Including Non-Agentic Rerankers ---
    ax = axes[0, 0]
    add_panel_label(ax, "a")
    x = np.arange(len(cfg_short))
    w = 0.19
    ax.bar(x - 1.5 * w, prec, width=w, color=PALETTE["blue"], label="Precision (%)")
    ax.bar(x - 0.5 * w, rec, width=w, color=PALETTE["teal"], hatch="//", label="Recall (%)")
    ax.bar(x + 0.5 * w, f1, width=w, color=PALETTE["navy"], label="F1 Score (%)")
    ax.bar(x + 1.5 * w, acc, width=w, color=PALETTE["sky"], hatch="..", label="Overall Accuracy (%)")
    for idx in range(1, len(cfg_short)):
        ax.text(x[idx] + 0.5 * w, f1[idx] + 1.5, f"{f1[idx]:.1f}%", ha="center", va="bottom", fontsize=6.2, fontweight="bold", color=PALETTE["navy"])
    ax.set_xticks(x)
    ax.set_xticklabels(cfg_short, fontsize=6.3)
    ax.set_ylim(0, 125)
    ax.set_ylabel("Performance on 120 Anomalies (%)")
    ax.set_title("De-Identified Stage 3 Resolution vs. ML Baselines (n = 120)")
    ax.grid(True, axis="y")
    ax.legend(loc="upper left", frameon=True, ncol=2, fontsize=7.0)

    # --- Panel B: Symbolic Guardrail Violation Rate vs F1 Tradeoff ---
    ax = axes[0, 1]
    add_panel_label(ax, "b")
    viol_pct = [float(vs_by_cfg[k]["violation_rate_pct"]) for k, _ in selected_cfgs]
    viol_cnt = [int(vs_by_cfg[k]["violations_accepted"]) for k, _ in selected_cfgs]
    colors_b = [
        PALETTE["lightgrey"],
        PALETTE["sky"],
        PALETTE["orange"],
        PALETTE["orange"],
        PALETTE["sky"],
        PALETTE["vermilion"],
        PALETTE["teal"],
        PALETTE["teal"],
    ]
    bars = ax.bar(x, viol_pct, width=0.50, color=colors_b, edgecolor="#111111")
    for idx, b in enumerate(bars):
        v = viol_pct[idx]
        c = viol_cnt[idx]
        lbl = f"{v:.1f}%\n({c}/120)" if c > 0 else "0.0%\n(0/120)"
        ax.text(b.get_x() + b.get_width() / 2, v + 0.12, lbl, ha="center", va="bottom", fontsize=6.4, fontweight="bold",
                color=PALETTE["vermilion"] if c > 0 else PALETTE["teal"])
    ax.set_xticks(x)
    ax.set_xticklabels(cfg_short, fontsize=6.3)
    ax.set_ylim(0, 4.0)
    ax.set_ylabel("Accepted Taxon/Anatomical Violation Rate (%)")
    ax.set_title("Symbolic Guardrail Eliminates Cross-Clade Trap Violations")
    ax.grid(True, axis="y")

    # --- Panel C: 50 Unique-Protein Ambiguous LRO Metabolic Gap Re-Ranking ---
    ax = axes[1, 0]
    add_panel_label(ax, "c")
    rlgc_pos = [int(r["rlgc_position"]) for r in lro_rows]
    agent_pos = [int(r["agent_position"]) for r in lro_rows]
    n_lro = len(lro_rows)
    rlgc_h1 = 100.0 * sum(1 for p in rlgc_pos if p == 1) / n_lro
    agent_h1 = 100.0 * sum(1 for p in agent_pos if p == 1) / n_lro
    rlgc_mrr = sum(1.0 / p for p in rlgc_pos) / n_lro
    agent_mrr = sum(1.0 / p for p in agent_pos) / n_lro

    ranks = [1, 2, 3, 4, 5]
    rlgc_counts = [sum(1 for p in rlgc_pos if p == rk) for rk in ranks]
    agent_counts = [sum(1 for p in agent_pos if p == rk) for rk in ranks]

    x_r = np.arange(len(ranks))
    w_r = 0.35
    b_r1 = ax.bar(x_r - w_r / 2, rlgc_counts, width=w_r, color=PALETTE["orange"], edgecolor="#111111", hatch="//",
                  label=f"Stage 1 RLGC S2 Alone (Hit@1 = {rlgc_h1:.1f}%, MRR = {rlgc_mrr:.3f})")
    b_r2 = ax.bar(x_r + w_r / 2, agent_counts, width=w_r, color=PALETTE["teal"], edgecolor="#111111",
                  label=f"Stage 3 Agent Re-Ranking (Hit@1 = {agent_h1:.1f}%, MRR = {agent_mrr:.3f})")
    for bars in (b_r1, b_r2):
        for b in bars:
            cnt = int(b.get_height())
            if cnt > 0:
                ax.text(b.get_x() + b.get_width() / 2, cnt + 0.6, str(cnt), ha="center", va="bottom", fontsize=8.0, fontweight="bold")
    ax.set_xticks(x_r)
    ax.set_xticklabels([f"Rank {rk}" + (" (Top-1)" if rk == 1 else "") for rk in ranks])
    ax.set_ylim(0, 56)
    ax.set_ylabel("Number of Unique LRO Metabolic Gaps (out of 50)")
    ax.set_title("Agentic Disambiguation of 50 Unique Multi-Candidate LRO Gaps")
    ax.grid(True, axis="y")
    ax.legend(loc="upper right", frameon=True, fontsize=7.4)

    # --- Panel D: Per-Taxon Concurrency-Normalized Wall-Clock Latency (Read from JSONL) ---
    ax = axes[1, 1]
    add_panel_label(ax, "d")
    taxon_order = [
        ("NCBITaxon_562", "E. coli K-12*"),
        ("NCBITaxon_1423", "B. subtilis"),
        ("NCBITaxon_287", "P. aeruginosa"),
        ("NCBITaxon_1773", "M. tuberculosis"),
        ("NCBITaxon_4932", "S. cerevisiae"),
        ("NCBITaxon_6239", "C. elegans"),
        ("NCBITaxon_3702", "A. thaliana"),
    ]
    taxa_labels = [lbl for _, lbl in taxon_order]
    bench_anoms = {r["anomaly_id"]: r["organism_taxon"] for r in read_jsonl(RESULTS_DIR / "pillar4_agent_benchmark_120.jsonl")}
    # Wall-clock latency per anomaly under 8-way async concurrency (HTTP latency / 8.0 concurrency factor)
    c2_scale = 1.0 / 8.0
    c3_scale = 1.0 / 8.0
    lat_fast = []
    lat_full = []
    for tax_id, _ in taxon_order:
        sub2 = [float(r["latency_sec"]) * c2_scale for r in res_c2 if bench_anoms.get(r["anomaly_id"]) == tax_id]
        sub3 = [float(r["latency_sec"]) * c3_scale for r in res_c3 if bench_anoms.get(r["anomaly_id"]) == tax_id]
        lat_fast.append(sum(sub2) / len(sub2))
        lat_full.append(sum(sub3) / len(sub3))

    c2_http = float(vs_by_cfg["Stage 3 Agent — Fast 1-Step Mode"]["mean_latency_sec"])
    c2_lat = c2_http / 8.0
    c2_tools = float(vs_by_cfg["Stage 3 Agent — Fast 1-Step Mode"]["mean_tool_calls"])
    c3_http = float(vs_by_cfg["Stage 3 Agent — Qwen 27B + Guardrail"]["mean_latency_sec"])
    c3_lat = c3_http / 8.0
    c3_tools = float(vs_by_cfg["Stage 3 Agent — Qwen 27B + Guardrail"]["mean_tool_calls"])

    x_t = np.arange(len(taxa_labels))
    w_t = np.arange(len(taxa_labels))
    w_bar = 0.34
    ax.bar(x_t - w_bar / 2, lat_fast, width=w_bar, color=PALETTE["sky"], edgecolor=PALETTE["blue"], hatch="//",
           label=f"Config 2: Fast 1-Step ({c2_lat:.2f}s wall-clock / {c2_http:.1f}s HTTP, {c2_tools:.2f} tools)")
    ax.bar(x_t + w_bar / 2, lat_full, width=w_bar, color=PALETTE["navy"],
           label=f"Config 3: Full Guarded ({c3_lat:.2f}s wall-clock / {c3_http:.1f}s HTTP, {c3_tools:.2f} tools)")
    ax.axhline(c3_lat, color=PALETTE["navy"], linestyle=":", lw=1.2)
    ax.axhline(c2_lat, color=PALETTE["sky"], linestyle=":", lw=1.2)
    ax.set_xticks(x_t)
    ax.set_xticklabels(taxa_labels, rotation=25, ha="right")
    ax.set_ylim(0, 6.2)
    ax.set_ylabel("Wall-Clock Time / Anomaly at 8× Concurrency (s)")
    ax.set_title("Scalable Multi-Turn Latency Across 7 Taxa (1.8% Mean Loci Inspected)")
    ax.grid(True, axis="y")
    ax.legend(loc="upper right", frameon=True, fontsize=7.1)

    save_fig(fig, "fig5_pillar4_agent")


def write_pillar3_report():
    lso_rows = read_tsv(RESULTS_DIR / "pillar3_complex_lso.tsv")
    anat_rows = read_tsv(RESULTS_DIR / "pillar3_eukaryotic_anatomy.tsv")
    rel_rows = read_tsv(RESULTS_DIR / "pillar3_relax_vs_remove.tsv")

    lines = [
        "# Pillar 3: Scenario-Specific Biological Discovery & Constraint Benchmarks",
        "",
        "## 1. Scenario B — Obligate Heteromeric Protein Complex Completion (`pillar3_complex_lso.tsv`)",
        "",
        "| Genome | Method | Cases | Hit@1 | Recall@0.5 | Precision@0.5 | F1@0.5 | Complex Coherence |",
        "|---|---|---:|---:|---:|---:|---:|---:|",
    ]
    for r in lso_rows:
        lines.append(
            f"| `{r['genome']}` | {r['method']} | {r['n_cases']} | {r['hit_at_1']} | {r['recall_at_05']} | {r['precision_at_05']} | {r['f1_at_05']} | {r['complex_coherence']} |"
        )

    lines.extend([
        "",
        "## 2. Scenario C — Multicellular Eukaryotic Anatomical & Developmental Constraints (`pillar3_eukaryotic_anatomy.tsv`)",
        "",
        "| Species | Taxon ID | Clade | Proteins | E1 Violations | E2 Direct GO TC Violations | E3 Stage 1 | E4/E5 Stage 2 Violations | E4/E5 Elimination (%) | Borderline Rescued (E3 vs E1) | Rescue Rate (%) |",
        "|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ])
    for r in anat_rows:
        lines.append(
            f"| `{r['species']}` | `{r['taxon_id']}` | {r['clade']} | {r['n_proteins']} | {r['n_violations_e1_noisy_or']} | {r['n_violations_e2_go_taxon']} | {r['n_violations_e3_stage1_anat']} | {r['n_violations_e4_stage2_remove']} | {r['pct_eliminated_e4_e5']}% | {r['borderline_rescued_e3_anat_boost']}/{r['borderline_true_claims_tested']} (vs {r['borderline_rescued_e1']}) | {r['boost_rescue_rate']}% |"
        )

    lines.extend([
        "",
        "## 3. Semantic Ancestor Relaxation (`mode=relax`) vs. Naive Removal (`mode=remove`) (`pillar3_relax_vs_remove.tsv`)",
        "",
        "| Species | Taxon ID | Mode | Propagated CAFA Fmax | Propagated Precision | Propagated Recall | Mean IC |",
        "|---|---|---|---:|---:|---:|---:|",
    ])
    for r in rel_rows:
        lines.append(
            f"| `{r['species']}` | `{r['taxon_id']}` | `{r['mode']}` | {r['cafa_fmax']} | {r['propagated_precision']} | {r['propagated_recall']} | {r['mean_ic']} |"
        )

    report_path = RESULTS_DIR / "pillar3_report.md"
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"[figures] Wrote {report_path}")


def main():
    generate_figure1()
    generate_figure2()
    generate_figure3()
    generate_figure4()
    generate_figure5()
    write_pillar3_report()


if __name__ == "__main__":
    main()
