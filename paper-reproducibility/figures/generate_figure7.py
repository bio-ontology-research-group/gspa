#!/usr/bin/env python3
"""
LEGACY - NOT RUN by generate_figures.py.  Figure 7 (AllTheBacteria) awaits the P11 rerun;
its current PDF/PNG in paper/figures/ were produced by this script from the earlier
results/allthebacteria/ files.  Moved here unchanged except for the repository path.

Generate Figure 7 (fig7_allthebacteria_panphylum.{pdf,png}) for the GSPA manuscript
directly from canonical results/allthebacteria/ TSV files:
  (a) Pan-Phylum Dark-Matter Burden (Zero-Alignment Orphans vs. Remote Homologs) & Rescue Rate Across Major Phyla
  (b) Genome-Size Bifurcation Law Across 72,530 Bacterial Genomes: U-Shaped Sequence Drift vs. 136x Metabolic Gap Surge
  (c) Stoichiometric Selection Across 72,533 Genomes: Rescued Multi-Subunit Complexes vs. Suppressed Domain-Reuse Singletons
  (d) Non-Orthologous Gene Displacement (pident = 0%) vs. Remote Sequence Divergence (0 < pident < 30%) Across Complexes & Pathways
"""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path
from typing import Dict, List

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[2]
ATB_DIR = REPO_ROOT / "results" / "allthebacteria"
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
    "font.size": 8.8,
    "axes.titlesize": 9.8,
    "axes.titleweight": "bold",
    "axes.labelsize": 9.0,
    "axes.labelweight": "bold",
    "xtick.labelsize": 8.0,
    "ytick.labelsize": 8.0,
    "legend.fontsize": 7.2,
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


def add_panel_label(ax, label: str, x: float = -0.16, y: float = 1.18):
    ax.text(
        x, y, label,
        transform=ax.transAxes,
        fontsize=13,
        fontweight="bold",
        va="top",
        ha="left",
        color="#111111",
    )


def generate_figure7() -> None:
    fig, axes = plt.subplots(2, 2, figsize=(15.4, 11.2))
    plt.subplots_adjust(wspace=0.46, hspace=0.48)

    # =========================================================================
    # Panel (a): Pan-Phylum Dark-Matter Burden & Whole-Genome Rescue Rate
    # =========================================================================
    ax_a = axes[0, 0]
    add_panel_label(ax_a, "a", x=-0.28, y=1.18)

    phy_summary = {r["phylum"]: r for r in read_tsv(ATB_DIR / "insight3b_phylum_summary.tsv")}
    phy_nogd = {r["phylum"]: r for r in read_tsv(ATB_DIR / "insight5b_phylum_nogd_vs_divergent.tsv")}

    selected_phyla = [
        ("p__Planctomycetota", "Planctomycetota (n=31)"),
        ("p__Bacteroidota", "Bacteroidota (n=3,341)"),
        ("p__Fibrobacterota", "Fibrobacterota (n=28)"),
        ("p__Methanobacteriota", "Methanobacteriota (Archaea)"),
        ("p__Cyanobacteriota", "Cyanobacteriota (n=115)"),
        ("p__Verrucomicrobiota", "Verrucomicrobiota (n=48)"),
        ("p__Spirochaetota", "Spirochaetota (n=216)"),
        ("p__Myxococcota", "Myxococcota (n=18)"),
        ("p__Acidobacteriota", "Acidobacteriota (n=27)"),
        ("p__Chlamydiota", "Chlamydiota (n=39)"),
        ("p__Fusobacteriota", "Fusobacteriota (n=67)"),
        ("p__Deinococcota", "Deinococcota (n=65)"),
        ("p__Campylobacterota", "Campylobacterota (n=327)"),
        ("p__Desulfobacterota_I", "Desulfobacterota_I (n=72)"),
        ("p__Bacillota", "Bacillota (n=12,675)"),
        ("p__Pseudomonadota", "Pseudomonadota (n=20,734)"),
        ("p__Actinomycetota", "Actinomycetota (n=31,304)"),
    ]
    selected_phyla.reverse()

    labels_a = [lbl for _, lbl in selected_phyla]
    orphan_pcts = []
    remote_pcts = []
    twi_resc_pcts = []
    for pkey, _ in selected_phyla:
        ps = phy_summary[pkey]
        pn = phy_nogd[pkey]
        n_orfs = int(ps["n_orfs"])
        o0 = int(pn["n_zero_align_orphan_orfs"])
        o1 = int(pn["n_remote_divergent_0_30_orfs"])
        orphan_pcts.append(100.0 * o0 / n_orfs)
        remote_pcts.append(100.0 * o1 / n_orfs)
        twi_resc_pcts.append(float(ps["pct_twilight_rescued"]))

    y_a = np.arange(len(labels_a))
    ax_a.barh(
        y_a, orphan_pcts, height=0.68,
        color=PALETTE["navy"], edgecolor="#111111", linewidth=0.5,
        label="Zero-Alignment Orphan Dark Matter (pident = 0%)",
    )
    ax_a.barh(
        y_a, remote_pcts, left=orphan_pcts, height=0.68,
        color=PALETTE["sky"], edgecolor="#111111", linewidth=0.5,
        label="Remote Divergent Dark Matter (0 < pident < 30%)",
    )
    ax_a.plot(
        twi_resc_pcts, y_a, "o-",
        color=PALETTE["teal"], markersize=4.6, linewidth=1.5,
        label="Twilight-Zone ORFs Rescued by GSPA (%)",
    )
    for i, (op, rp, tr) in enumerate(zip(orphan_pcts, remote_pcts, twi_resc_pcts)):
        tot_d = op + rp
        ax_a.text(op / 2.0, i, f"{tot_d:.1f}%", va="center", ha="center", fontsize=6.8, fontweight="bold", color="white")
        ax_a.text(max(tot_d, tr) + 1.1, i, f"{tr:.1f}%", va="center", ha="left", fontsize=6.8, color=PALETTE["teal"], fontweight="bold")

    ax_a.set_yticks(y_a)
    ax_a.set_yticklabels(labels_a, fontsize=7.6)
    ax_a.set_xlim(0, 92)
    ax_a.set_xlabel("Percentage of Proteome ORFs / Rescued Twilight Fraction (%)")
    ax_a.set_title("Pan-Phylum Sequence Dark Matter & Rescue (292.1M ORFs)", pad=38)
    ax_a.grid(axis="x", alpha=0.25)
    ax_a.legend(
        loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=2,
        frameon=True, facecolor="white", framealpha=0.95, fontsize=6.8,
    )

    # =========================================================================
    # Panel (b): Genome-Size Bifurcation Curve Across 72,530 Bacterial Genomes
    # =========================================================================
    ax_b = axes[0, 1]
    add_panel_label(ax_b, "b", x=-0.16, y=1.18)

    bif_rows = read_tsv(ATB_DIR / "insight6_genome_size_bifurcation.tsv")
    x_labels_b = [
        "<1,500\n(n=413)",
        "1,500-2,499\n(n=5,114)",
        "2,500-3,499\n(n=15,283)",
        "3,500-4,499\n(n=36,243)",
        "4,500-5,499\n(n=9,323)",
        ">=5,500\n(n=6,154)",
    ]
    x_b = np.arange(len(bif_rows))
    dark_b = [float(r["mean_dark_pct"]) for r in bif_rows]
    lt50_b = [float(r["mean_dark_pct"]) + float(r["mean_twilight_pct"]) for r in bif_rows]
    b5_gaps = [float(r["mean_B5_metabolic_gaps_per_genome"]) for r in bif_rows]
    b6_per_orf = [float(r["mean_B6_operon_rescues_per_orf"]) * 10.0 for r in bif_rows]

    ax_b2 = ax_b.twinx()
    ax_b2.spines["top"].set_visible(False)

    bars_b5 = ax_b2.bar(
        x_b + 0.20, b5_gaps, width=0.28,
        color=PALETTE["vermilion"], alpha=0.82, edgecolor="#111111", linewidth=0.6,
        label="Metabolic Pathway Gaps / Genome (B5, Right Axis)",
    )
    l1 = ax_b.plot(
        x_b - 0.08, lt50_b, "s--",
        color=PALETTE["blue"], linewidth=1.9, markersize=5.6,
        label="Dark + Twilight (<50% Identity, Left Axis)",
    )
    l2 = ax_b.plot(
        x_b - 0.08, dark_b, "o-",
        color=PALETTE["navy"], linewidth=2.2, markersize=6.0,
        label="Dark Matter (<30% Identity, Left Axis)",
    )
    l3 = ax_b.plot(
        x_b - 0.08, b6_per_orf, "^:",
        color=PALETTE["teal"], linewidth=1.8, markersize=5.4,
        label="Operon Synteny Rescues / 10 ORFs (B6, Left Axis)",
    )

    for xi, (dval, gval) in enumerate(zip(dark_b, b5_gaps)):
        y_off = -4.8 if xi == 0 else 1.8
        ax_b.text(xi - 0.10, dval + y_off, f"{dval:.1f}%", ha="center", va="bottom", fontsize=7.1, fontweight="bold", color=PALETTE["navy"])
        ax_b2.text(xi + 0.20, gval + 0.45, f"{gval:.2f}" if gval < 2 else f"{gval:.1f}", ha="center", va="bottom", fontsize=7.1, fontweight="bold", color=PALETTE["vermilion"])

    ax_b.annotate(
        "Streamlined Drift Regime\n(50.6% Dark, 0.14 Gaps/Genome)",
        xy=(0 - 0.08, dark_b[0]), xytext=(0.35, 55.5),
        arrowprops=dict(arrowstyle="->", color=PALETTE["navy"], lw=1.0),
        fontsize=7.0, fontweight="bold", color=PALETTE["navy"],
    )
    ax_b2.annotate(
        "136.4x Metabolic Gap Surge\n(19.1 vs. 0.14 Gaps/Genome)",
        xy=(5.05, b5_gaps[-1]), xytext=(2.55, 19.6),
        arrowprops=dict(arrowstyle="->", color=PALETTE["vermilion"], lw=1.0),
        fontsize=7.0, fontweight="bold", color=PALETTE["vermilion"],
    )

    ax_b.set_xticks(x_b)
    ax_b.set_xticklabels(x_labels_b, fontsize=7.5)
    ax_b.set_ylim(18, 95)
    ax_b2.set_ylim(0, 24.5)
    ax_b.set_xlabel("Bacterial Genome Size Stratum (Coding Sequences / Genome)")
    ax_b.set_ylabel("Proteome Dark / Twilight (%) or B6 / 10 ORFs", color=PALETTE["navy"])
    ax_b2.set_ylabel("Mean Metabolic Pathway Gaps / Genome (B5)", color=PALETTE["vermilion"])
    ax_b.set_title("Genome-Size Bifurcation: Drift vs. Metabolic Innovation (n=72,530)", pad=38)
    ax_b.grid(axis="y", alpha=0.25)

    handles_b = l1 + l2 + l3 + [bars_b5]
    labels_b_leg = [h.get_label() for h in handles_b]
    ax_b.legend(
        handles_b, labels_b_leg,
        loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=2,
        frameon=True, facecolor="white", framealpha=0.95, fontsize=6.6,
    )

    # =========================================================================
    # Panel (c): Stoichiometric Selection vs. Modular Domain Reuse Singletons
    # =========================================================================
    ax_c = axes[1, 0]
    add_panel_label(ax_c, "c", x=-0.28, y=1.18)

    cpx_rows = read_tsv(ATB_DIR / "insight2_complex_stoichiometry_and_operons.tsv")
    cpx_totals: Dict[str, List[int]] = defaultdict(lambda: [0, 0])
    for r in cpx_rows:
        cat = r["complex_or_locus_group"]
        if cat == "Operon_CoLocalized_Loci":
            continue
        cpx_totals[cat][0] += int(r["coherence_rescued_subunits"])
        cpx_totals[cat][1] += int(r["singleton_hallucinations_suppressed"])

    cpx_order = [
        ("GO:0005840_Ribosome", "70S Ribosome\n(GO:0005840)"),
        ("GO:0000428_RNA_polymerase", "RNA Polymerase\n(GO:0000428)"),
        ("GO:0045259_ATP_synthase", "F1F0 ATP Synthase\n(GO:0045259)"),
        ("GO:0009330_DNA_gyrase", "DNA Gyrase\n(GO:0009330)"),
        ("GO:0005747_NADH_dehydrogenase", "Respiratory Complex I\n(GO:0005747)"),
        ("GO:0031522_Sec_translocon", "Sec Translocon\n(GO:0031522)"),
    ]
    cpx_order.reverse()

    y_c = np.arange(len(cpx_order))
    resc_vals = [max(1, cpx_totals[k][0]) for k, _ in cpx_order]
    supp_vals = [max(1, cpx_totals[k][1]) for k, _ in cpx_order]
    raw_resc = [cpx_totals[k][0] for k, _ in cpx_order]
    raw_supp = [cpx_totals[k][1] for k, _ in cpx_order]

    h_c = 0.36
    ax_c.barh(
        y_c + h_c / 2, resc_vals, height=h_c,
        color=PALETTE["teal"], edgecolor="#111111", linewidth=0.6,
        label="Complete-Complex Subunits Rescued (12.51M Total)",
    )
    ax_c.barh(
        y_c - h_c / 2, supp_vals, height=h_c,
        color=PALETTE["vermilion"], edgecolor="#111111", linewidth=0.6,
        label="Isolated Singleton Hallucinations Pruned (1.43M Total)",
    )
    ax_c.set_xscale("log")
    ax_c.set_xlim(0.8, 8.5e7)

    for i, (rv, sv) in enumerate(zip(raw_resc, raw_supp)):
        r_str = f"{rv:,}" if rv < 1_000_000 else f"{rv/1e6:.2f}M"
        s_str = f"{sv:,}" if sv < 1_000_000 else f"{sv/1e6:.2f}M"
        ax_c.text(max(1.2, rv * 1.18), i + h_c / 2, r_str, va="center", ha="left", fontsize=7.1, fontweight="bold", color=PALETTE["teal"])
        ax_c.text(max(1.2, sv * 1.18), i - h_c / 2, s_str, va="center", ha="left", fontsize=7.1, fontweight="bold", color=PALETTE["vermilion"])

    ax_c.set_yticks(y_c)
    ax_c.set_yticklabels([lbl for _, lbl in cpx_order], fontsize=7.8)
    ax_c.set_xlabel("Subunit Annotations Across 72,533 Genomes (Log10 Scale)")
    ax_c.set_title("Stoichiometric Selection Separates Complexes from Domain Reuse", pad=38)
    ax_c.grid(axis="x", alpha=0.25)
    ax_c.legend(
        loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=1,
        frameon=True, facecolor="white", framealpha=0.95, fontsize=6.8,
    )

    # =========================================================================
    # Panel (d): Non-Orthologous Displacement (pident=0%) vs. Remote Divergence
    # =========================================================================
    ax_d = axes[1, 1]
    add_panel_label(ax_d, "d", x=-0.34, y=1.18)

    nogd_rows = read_tsv(ATB_DIR / "insight5_non_orthologous_displacement.tsv")
    nogd_lookup = {(r["mechanism"], r["pathway_or_complex_slot"], r["go_term"]): r for r in nogd_rows}

    selected_slots = [
        ("B5_pathway_gapfill", "pathway=KEGG:00982", "GO:0004364", "Thiol / Mycothiol S-Transferase (B5)"),
        ("B2_complex_module_coherence", "group=GO:0045259_ATP_synthase", "GO:0045259", "ATP Synthase F0 Membrane Subunits (B2)"),
        ("B2_complex_module_coherence", "group=GO:0000428_RNA_polymerase", "GO:0000428", "RNA Polymerase Accessory Subunits (B2)"),
        ("B2_complex_module_coherence", "group=KEGG:M00859", "GO:0004016", "Adenylate Cyclase Signaling (M00859)"),
        ("B2_complex_module_coherence", "group=KEGG:M00888", "GO:0016758", "Cell-Wall LPS Hexosyltransferase (M00888)"),
        ("B2_complex_module_coherence", "group=KEGG:M00836", "GO:0016829", "Cofactor C-O / C-C Lyase Module (M00836)"),
        ("B2_complex_module_coherence", "group=KEGG:M00951", "GO:0016874", "Purine / Peptide C-N Ligase (M00951)"),
        ("B2_complex_module_coherence", "group=GO:0005840_Ribosome", "GO:0006412", "70S Ribosomal Subunits (GO:0005840)"),
        ("B5_pathway_gapfill", "pathway=KEGG:00982", "GO:0015020", "Xenobiotic Glycosyltransferase (B5)"),
        ("B5_pathway_gapfill", "pathway=KEGG:00982", "GO:0016712", "Flavin/Heme Oxidoreductase (B5)"),
        ("B2_complex_module_coherence", "group=GO:0009330_DNA_gyrase", "GO:0006265", "DNA Gyrase Catalytic Core (GO:0009330)"),
        ("B2_complex_module_coherence", "group=GO:0000428_RNA_polymerase", "GO:0006351", "RNA Polymerase Catalytic Core (B2)"),
        ("B2_complex_module_coherence", "group=GO:0045259_ATP_synthase", "GO:0046933", "ATP Synthase Rotary Catalysis (B2)"),
        ("B5_pathway_gapfill", "pathway=KEGG:00980", "GO:0004090", "NADPH Carbonyl Reductase Core (B5)"),
    ]
    selected_slots.reverse()

    labels_d = [lbl for _, _, _, lbl in selected_slots]
    nogd_pct_d = []
    rem_pct_d = []
    tot_dark_d = []
    n_phy_d = []
    for m, s, g, _ in selected_slots:
        r = nogd_lookup[(m, s, g)]
        p_n = float(r["pct_nogd_of_dark"])
        nogd_pct_d.append(p_n)
        rem_pct_d.append(100.0 - p_n)
        tot_dark_d.append(int(r["n_total_dark_rescued_lt30"]))
        n_phy_d.append(int(r["n_phyla_with_dark_rescue"]))

    y_d = np.arange(len(labels_d))
    ax_d.barh(
        y_d, nogd_pct_d, height=0.68,
        color=PALETTE["purple"], edgecolor="#111111", linewidth=0.5,
        label="Non-Orthologous Replacement / Zero-Alignment (pident = 0%)",
    )
    ax_d.barh(
        y_d, rem_pct_d, left=nogd_pct_d, height=0.68,
        color=PALETTE["orange"], edgecolor="#111111", linewidth=0.5,
        label="Remote Divergent Ortholog (0% < pident < 30%)",
    )

    for i, (npct, td, nph) in enumerate(zip(nogd_pct_d, tot_dark_d, n_phy_d)):
        td_str = f"{td:,}" if td < 100_000 else f"{td/1e6:.2f}M"
        ax_d.text(101.5, i, f"{npct:.1f}% NOGD ({td_str}, {nph}p)", va="center", ha="left", fontsize=6.9, color="#222222")

    ax_d.set_yticks(y_d)
    ax_d.set_yticklabels(labels_d, fontsize=7.4)
    ax_d.set_xlim(0, 142)
    ax_d.set_xticks([0, 25, 50, 75, 100])
    ax_d.set_xlabel("Share of Rescued Dark-Matter Annotations (<30% Identity) (%)")
    ax_d.set_title("Non-Orthologous Displacement vs. Remote Sequence Divergence", pad=38)
    ax_d.grid(axis="x", alpha=0.25)
    ax_d.legend(
        loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=1,
        frameon=True, facecolor="white", framealpha=0.95, fontsize=6.8,
    )

    out_pdf = FIG_DIR / "fig7_allthebacteria_panphylum.pdf"
    out_png = FIG_DIR / "fig7_allthebacteria_panphylum.png"
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out_pdf} and {out_png}")


if __name__ == "__main__":
    generate_figure7()
