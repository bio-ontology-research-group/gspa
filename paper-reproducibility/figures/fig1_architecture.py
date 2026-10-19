"""Figure 1: four-stage architecture (text-box diagram).

Only the text of the earlier diagram was corrected (model name, Stage 3 trigger rate);
the redesign of this figure is a separate task.  The trigger rate and the anomaly cap
are read from results/pillar4_v2_trigger_rates.tsv (genome=SUMMARY_ALL).
"""

from __future__ import annotations

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

import common as C

PALETTE = {
    "navy": "#1B365D",
    "blue": "#0072B2",
    "teal": "#009E73",
    "orange": "#E69F00",
    "vermilion": "#D55E00",
}

LEGACY_RC = {
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "pdf.fonttype": 42,
    "figure.dpi": 300,
}


def trigger_text():
    C.panel("fig1", "Stage 3 box")
    pooled = C.kv(C.get_raw("pillar4_v2_trigger_rates.tsv", {"genome": "SUMMARY_ALL"}, "pct",
                            "pooled percentage of proteins triggered under the cap"), "pooled")
    capped = [r for r in C.load("pillar4_v2_trigger_rates.tsv")
              if not r["genome"].startswith("SUMMARY") and r.get("cap_reached") == "True"]
    cap = int(max(C.num(r["n_anomalies"]) for r in capped))
    C.note("pillar4_v2_trigger_rates.tsv", "cap_reached=True (per-genome rows)", "n_anomalies",
           "max over capped genomes = production anomaly cap", cap)
    return f"• Inspects about {round(pooled):d}% of loci", f"  under the {cap}-anomaly cap"


def generate_figure1():
    trig1, trig2 = trigger_text()
    with plt.rc_context(LEGACY_RC):
        _draw(trig1, trig2)


def _draw(trig1, trig2):
    fig, ax = plt.subplots(figsize=(13.8, 6.6))
    ax.set_xlim(0, 109)
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
                x + 0.9, top_y - idx * 2.02,
                clean,
                ha="left", va="top",
                fontsize=7.35,
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
                fontsize=7.3, fontweight="bold", color="#34495E",
                bbox=dict(boxstyle="round,pad=0.10", facecolor="#FFFFFF", edgecolor="none", alpha=0.92),
            )

    # Top Title Banner
    ax.text(
        54.5, 52.0,
        "GSPA: Unified Four-Stage Genome-Scale Protein Function Annotation Architecture",
        ha="center", va="center", fontsize=12, fontweight="bold", color=PALETTE["navy"],
    )

    # Column 1: Multi-Evidence Inputs & Genome Context
    draw_box(
        0.8, 26.5, 20.8, 23.0,
        "1A. Stage 0 Evidence",
        [
            "*Homology & Domain Channels:",
            "• DIAMOND (Swiss-Prot GOA)",
            "• Pfam / HMMER 3 + InterProScan",
            "*Structure & PLM Channels:",
            "• ESM-2 / ProtT5 kNN + CLEAN",
            "• DeepGOZero/SE + ProteInfer",
            "• FoldSeek / ProstT5 Structure",
            "*Platt Logistic Calibration",
        ],
        fc=PALETTE["navy"], ec=PALETTE["navy"],
    )

    draw_box(
        0.8, 1.5, 20.8, 22.5,
        "1B. Whole-Genome Context",
        [
            "*Genomic Coordinates (GFF3):",
            "• Co-directional Operon Calls",
            "*Metabolic Network (gapseq/SEED):",
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
        27.4, 9.5, 24.8, 40.0,
        "Stage 1: Bayesian Refiner",
        [
            "*Log-Odds Fixed-Point Update (k=1..6):",
            "  L(p,f) = L₀(p,f) + Σ w_m B_m(p,f)",
            "",
            "*Six Genome-Scale Priors:",
            "• Essentiality Prior (+α_ess = 1.5):",
            "  Universal/domain core GO categories",
            "• Coherence Prior (+α_coh, −w_sing):",
            "  Pathway & heteromeric complex support",
            "• Taxon Consistency Prior (−α_cons = −0.68):",
            "  SAT minimal unsatisfiable cores",
            "• Anatomical Prior (−α_anat, +α_boost):",
            "  Uberon/PO/CL cross-product bridge",
            "• Metabolic Gap Prior + RLGC Kernel:",
            "  Reaction-local genomic context + PLM",
            "• Operon Context Prior (+α_ctx = 1.0):",
            "  Co-operonic pathway/complex boost",
        ],
        fc=PALETTE["teal"], ec=PALETTE["teal"],
    )

    # Column 3: Stage 2 Symbolic SAT & Stoichiometric Enforcement
    draw_box(
        58.0, 26.0, 21.8, 23.5,
        "Stage 2A/2B: Symbolic & pFBA",
        [
            "*2A. Deterministic Ontology Gate:",
            "• SAT Taxon & Uberon/PO/CL Check",
            "• Semantic Ancestor Relaxation:",
            "  Replaces forbidden child t with",
            "  nearest valid GO ancestor a*",
            "*2B. Stoichiometric Locus Refiner:",
            "• Positive-cost LP + pFBA flux",
            "  certification (|v_j| > 10⁻⁶)",
            "• Outer feedback ΔL_FBA to Stage 1",
        ],
        fc=PALETTE["orange"], ec=PALETTE["orange"], title_color="#111111",
    )

    # Column 3 Bottom: Output & GAEF Evaluator
    draw_box(
        58.0, 1.5, 21.8, 21.0,
        "Output & Structural Audit",
        [
            "*Calibrated Proteome Output:",
            "• GAF 2.2 / TSV / JSON Reports",
            "*Genome-Scale Structural Metrics:",
            "• Essential Completeness (C_ess)",
            "• Pathway & Complex Coherence",
            "• Taxon/Anatomical Consistency",
            "• Composite Structural Score",
        ],
        fc=PALETTE["navy"], ec=PALETTE["navy"],
    )

    # Column 4: Stage 3 Optional Neuro-Symbolic Agentic Layer
    draw_box(
        85.8, 9.5, 22.0, 40.0,
        "Stage 3: Guarded Agent",
        [
            "*Selective Anomaly Trigger:",
            trig1,
            trig2,
            "• Pathway Gap & Complex Subunit",
            "• Borderline Essential Function",
            "• Borderline Taxon/Anat Conflict",
            "",
            "*Tool-Using LLM (Qwen3.8-27B):",
            "• Operon neighborhood inspection",
            "• Pathway & complex gap inspection",
            "• Per-protein evidence lookup",
            "• Symbolic constraint checker",
            "",
            "*Symbolic Guardrail Gate:",
            "• Blocks invalid GO/EC identifiers",
            "• Enforces 0.0% taxon/anatomical",
            "  violation rate on promotions",
        ],
        fc=PALETTE["vermilion"], ec=PALETTE["vermilion"],
    )

    # Connecting arrows with generous whitespace label placement
    draw_arrow(21.9, 38.0, 27.1, 35.5, "Claims", lx=24.5, ly=38.6)
    draw_arrow(21.9, 13.5, 27.1, 21.0, "Priors", lx=24.5, ly=19.8)
    draw_arrow(52.5, 36.0, 57.7, 36.8, "Posteriors", lx=55.1, ly=38.5)
    draw_arrow(57.7, 30.0, 52.5, 29.0, "ΔL_FBA", lx=55.1, ly=27.1)
    draw_arrow(68.9, 25.7, 68.9, 22.8, "Enforced", lx=73.4, ly=24.3)
    draw_arrow(80.1, 38.0, 85.5, 38.0, "Anomalies", lx=82.8, ly=40.0)
    draw_arrow(85.5, 17.5, 80.1, 15.0, "Guarded", lx=82.8, ly=14.0)

    fig.savefig(C.FIG_DIR / "fig1_architecture.pdf", bbox_inches="tight")
    fig.savefig(C.FIG_DIR / "fig1_architecture.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("[figures] wrote paper/figures/fig1_architecture.pdf and .png")
