#!/usr/bin/env python3
"""
Extracts two macro-evolutionary tables across the 72,533 AllTheBacteria + Eukaryote genomes
using pure Python standard library (zero external dependencies):
1. `insight6_genome_size_bifurcation.tsv`: Genome-size scaling and bifurcation between
   drift-driven core sequence dark matter in streamlined genomes (< 1,500 ORFs) and
   innovation-driven metabolic pathway gaps in large environmental genomes (>= 5,500 ORFs).
2. `insight5_non_orthologous_displacement.tsv` & `insight5b_phylum_nogd_vs_divergent.tsv`:
   Pan-phylum breakdown of coding sequences and rescued pathway gaps (`B5`),
   metabolic modules and obligate complexes (`B2`), and co-directional operons (`B6`) into:
   - Zero-Alignment Non-Homologous Replacements (`max_pident == 0.0%`, zero BLAST hit E <= 1e-3 in Swiss-Prot)
   - Remote Divergent Dark-Matter Orthologs (`0.0% < max_pident < 30.0%`)
   - Twilight-Zone Orthologs (`30.0% <= max_pident < 50.0%`)
   - Conserved Homologs (`max_pident >= 50.0%`)
"""
from __future__ import annotations

import argparse
import concurrent.futures
import os
import time
from collections import defaultdict
from pathlib import Path


GO_NAMES = {
    "GO:0005840_Ribosome": "Ribosome (30S/50S structural subunits)",
    "GO:0000428_RNA_polymerase": "DNA-directed RNA polymerase (core & sigma)",
    "GO:0045259_ATP_synthase": "Proton-transporting F1F0 ATP synthase",
    "GO:0009330_DNA_gyrase": "DNA topoisomerase (ATP-hydrolyzing) / gyrase",
    "GO:0005747_NADH_dehydrogenase": "Respiratory NADH dehydrogenase (Complex I)",
    "GO:0031522_Sec_translocon": "Sec protein export translocon complex",
    "GO:0015020": "Glucuronosyltransferase / xenobiotic glycosyltransferase",
    "GO:0008177": "Succinate dehydrogenase (ubiquinone) / fumarate reductase",
    "GO:0004031": "Aldehyde oxidase / xanthine dehydrogenase family",
    "GO:0004364": "Glutathione / mycothiol S-transferase activity",
    "GO:0016712": "Oxidoreductase activity (flavin/heme incorporation)",
    "GO:0004090": "Carbonyl reductase (NADPH) activity",
    "GO:0097621": "Secondary metabolite monooxygenase activity",
    "GO:0016874": "Ligase activity (C-N / peptide bond formation)",
    "GO:0016758": "Hexosyltransferase / cell-envelope glycosyltransferase",
    "GO:0016829": "Carbon-oxygen / carbon-carbon lyase activity",
    "GO:0016765": "Transferase activity (alkyl or aryl prenyl groups)",
    "GO:0008483": "Transaminase / pyridoxal-phosphate aminotransferase",
    "GO:0016830": "Carbon-carbon lyase activity",
    "GO:0016831": "Carboxy-lyase / decarboxylase activity",
    "GO:0016769": "Transferase activity (nitrogenous groups)",
    "GO:0016854": "Racemase and epimerase activity (amino acid/peptidoglycan)",
    "GO:0006811": "Monoatomic ion transport (operonic permease/channel)",
    "GO:0051179": "Localization / membrane transport subunit",
    "GO:0005886": "Plasma membrane / inner membrane complex",
    "GO:0006412": "Translation / ribosomal machinery",
    "GO:0006351": "DNA-templated transcription",
    "GO:0015986": "Proton motive force-driven ATP synthesis",
}

PATHWAY_MODULE_NAMES = {
    "pathway=KEGG:00982": "KEGG map00982: Xenobiotic & secondary metabolite transformation",
    "pathway=KEGG:00190": "KEGG map00190: Oxidative phosphorylation & bioenergetic electron transport",
    "pathway=KEGG:00980": "KEGG map00980: Cytochrome P450 / xenobiotic catabolism",
    "group=KEGG:M00926": "KEGG M00926: Heme / tetrapyrrole cofactor biosynthesis module",
    "group=KEGG:M01041": "KEGG M01041: Branched-chain & aromatic amino acid module",
    "group=KEGG:M00951": "KEGG M00951: Purine / nucleotide salvage & de novo ligase module",
    "group=KEGG:M00888": "KEGG M00888: Cell-wall lipopolysaccharide / O-antigen glycosyltransferase module",
    "group=KEGG:M01037": "KEGG M01037: Pyridoxal-phosphate transamination module",
    "group=KEGG:M00785": "KEGG M00785: Polyamine / nitrogenous group transfer module",
    "group=KEGG:M00930": "KEGG M00930: Menaquinone / ubiquinone chorismate lyase module",
    "group=KEGG:M00875": "KEGG M00875: Amino acid decarboxylation & cofactor module",
    "group=KEGG:M00840": "KEGG M00840: Peptidoglycan D-Ala/D-Glu racemase module",
    "group=GO:0005840_Ribosome": "Complex GO:0005840: 70S Ribosome (30S + 50S obligate subunits)",
    "group=GO:0000428_RNA_polymerase": "Complex GO:0000428: DNA-directed RNA polymerase holoenzyme",
    "group=GO:0045259_ATP_synthase": "Complex GO:0045259: Proton-transporting F1F0 ATP synthase",
    "group=GO:0009330_DNA_gyrase": "Complex GO:0009330: Heterotetrameric DNA gyrase (GyrA2B2)",
    "group=GO:0005747_NADH_dehydrogenase": "Complex GO:0005747: 14-subunit Respiratory Complex I (NuoA-N)",
    "group=GO:0031522_Sec_translocon": "Complex GO:0031522: SecYEG protein translocation channel",
}


def compute_genome_size_bifurcation(in_tsv: Path, out_tsv: Path) -> None:
    strata = [
        (0, 1500, "< 1,500 ORFs (Streamlined / Host-Associated)"),
        (1500, 2500, "1,500-2,499 ORFs (Small Free-Living / Host-Adapted)"),
        (2500, 3500, "2,500-3,499 ORFs (Compact Generalists)"),
        (3500, 4500, "3,500-4,499 ORFs (Modal Bacterial Generalists)"),
        (4500, 5500, "4,500-5,499 ORFs (Versatile Environmental)"),
        (5500, 100000, ">= 5,500 ORFs (Large Environmental / Secondary Metabolite)"),
    ]
    buckets: dict[int, list[dict[str, object]]] = defaultdict(list)

    with open(in_tsv) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        idx = {h: i for i, h in enumerate(hdr)}
        for line in fh:
            if not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if p[idx["domain"]] != "d__Bacteria":
                continue
            n_orfs = int(p[idx["n_orfs"]])
            if n_orfs <= 0:
                continue
            for si, (lo, hi, _) in enumerate(strata):
                if lo <= n_orfs < hi:
                    buckets[si].append(
                        {
                            "species": p[idx["species"]],
                            "phylum": p[idx["phylum"]],
                            "n_orfs": n_orfs,
                            "dark_pct": 100.0 * int(p[idx["n_dark_orfs_lt30"]]) / n_orfs,
                            "twi_pct": 100.0 * int(p[idx["n_twilight_orfs_30_50"]]) / n_orfs,
                            "b6": int(p[idx["n_B6_operon_rescues"]]),
                            "b2": int(p[idx["n_B2_coherence_rescues"]]),
                            "b5": int(p[idx["n_metabolic_gaps_candidate"]]),
                            "sing": int(p[idx["n_B2_singleton_suppressed"]]),
                            "viol": int(p[idx["n_B3_taxon_suppressed"]]) + int(p[idx["n_B4_anatomical_suppressed"]]),
                        }
                    )
                    break

    cols = [
        "genome_size_stratum",
        "min_orfs",
        "max_orfs",
        "n_genomes",
        "n_species",
        "n_phyla",
        "mean_orfs",
        "mean_dark_pct",
        "mean_twilight_pct",
        "mean_B6_operon_rescues_per_genome",
        "mean_B6_operon_rescues_per_orf",
        "mean_B2_coherence_rescues_per_genome",
        "mean_B2_coherence_rescues_per_orf",
        "mean_B5_metabolic_gaps_per_genome",
        "mean_B5_gaps_per_1k_orfs",
        "mean_singleton_suppressed_per_genome",
        "mean_violations_suppressed_per_genome",
    ]
    with open(out_tsv, "w") as out:
        out.write("\t".join(cols) + "\n")
        for si, (lo, hi, label) in enumerate(strata):
            rows = buckets[si]
            if not rows:
                continue
            n = len(rows)
            n_sp = len({r["species"] for r in rows})
            n_phy = len({r["phylum"] for r in rows})
            m_orfs = sum(float(r["n_orfs"]) for r in rows) / n
            m_dark = sum(float(r["dark_pct"]) for r in rows) / n
            m_twi = sum(float(r["twi_pct"]) for r in rows) / n
            m_b6 = sum(float(r["b6"]) for r in rows) / n
            m_b6_orf = sum(float(r["b6"]) / float(r["n_orfs"]) for r in rows) / n
            m_b2 = sum(float(r["b2"]) for r in rows) / n
            m_b2_orf = sum(float(r["b2"]) / float(r["n_orfs"]) for r in rows) / n
            m_b5 = sum(float(r["b5"]) for r in rows) / n
            m_b5_1k = sum(1000.0 * float(r["b5"]) / float(r["n_orfs"]) for r in rows) / n
            m_sing = sum(float(r["sing"]) for r in rows) / n
            m_viol = sum(float(r["viol"]) for r in rows) / n
            vals = [
                label,
                str(lo),
                str(hi if hi < 100000 else 25000),
                str(n),
                str(n_sp),
                str(n_phy),
                f"{m_orfs:.1f}",
                f"{m_dark:.2f}",
                f"{m_twi:.2f}",
                f"{m_b6:.1f}",
                f"{m_b6_orf:.3f}",
                f"{m_b2:.1f}",
                f"{m_b2_orf:.3f}",
                f"{m_b5:.2f}",
                f"{m_b5_1k:.3f}",
                f"{m_sing:.2f}",
                f"{m_viol:.1f}",
            ]
            out.write("\t".join(vals) + "\n")
    print(f"Saved {out_tsv}")


def scan_batch_for_nogd(wdir_str: str, sample_meta: dict[str, dict[str, str]]) -> dict[str, object] | None:
    wdir = Path(wdir_str)
    out_dir = wdir / "stageC" / "out"
    if not out_dir.exists():
        return None

    completed_gids: dict[str, str] = {}
    with os.scandir(out_dir) as it:
        for entry in it:
            if entry.name.endswith(".delta_events.tsv") and entry.is_file(follow_symlinks=False):
                gid = entry.name[: -len(".delta_events.tsv")]
                completed_gids[gid] = entry.path
    if not completed_gids:
        return None

    # 1. Read max_pident per cluster rep from stageB/reps_vs_sprot.m8 (1 file per batch)
    rep_pident: dict[str, float] = {}
    m8_p = wdir / "stageB" / "reps_vs_sprot.m8"
    if m8_p.exists():
        with open(m8_p) as fh:
            for line in fh:
                p = line.split("\t")
                if len(p) >= 3:
                    rid = p[0]
                    pid_val = float(p[2])
                    if pid_val > rep_pident.get(rid, 0.0):
                        rep_pident[rid] = pid_val

    # 2. Count total ORFs per phylum split into [pident==0 (orphan), 0<pident<30 (remote), 30<=pident<50 (twilight), pident>=50]
    phylum_orf_regimes: dict[str, list[int]] = defaultdict(lambda: [0, 0, 0, 0])
    map_p = wdir / "dedup" / "genome_protein_to_rep.tsv"
    if map_p.exists():
        with open(map_p) as fh:
            fh.readline()
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) >= 4:
                    gid = p[0]
                    if gid not in completed_gids:
                        continue
                    phy = sample_meta.get(gid, {}).get("phylum", "p__Pseudomonadota")
                    pid_val = rep_pident.get(p[3], 0.0)
                    if pid_val == 0.0:
                        phylum_orf_regimes[phy][0] += 1
                    elif pid_val < 30.0:
                        phylum_orf_regimes[phy][1] += 1
                    elif pid_val < 50.0:
                        phylum_orf_regimes[phy][2] += 1
                    else:
                        phylum_orf_regimes[phy][3] += 1

    # 3. Scan delta_events.tsv files (all in pilot/small batches, or stratified up to 80 per batch for speed)
    slot_counts: dict[tuple[str, str, str], list[float]] = defaultdict(lambda: [0.0, 0.0, 0.0, 0.0])
    slot_phyla: dict[tuple[str, str, str], set[str]] = defaultdict(set)
    phylum_nogd: dict[str, list[float]] = defaultdict(lambda: [0.0, 0.0, 0.0, 0.0])

    by_phy: dict[str, list[tuple[str, str]]] = defaultdict(list)
    for gid, dpath in sorted(completed_gids.items()):
        phy = sample_meta.get(gid, {}).get("phylum", "p__Pseudomonadota")
        by_phy[phy].append((gid, dpath))

    for phy, glist in by_phy.items():
        n_tot = len(glist)
        # Keep all genomes for rare/minor phyla (< 30 genomes in batch) and sample up to 60 for huge single-phylum batches
        sample_glist = glist if n_tot <= 60 else glist[:: max(1, n_tot // 60)][:60]
        scale = float(n_tot) / float(len(sample_glist))

        for gid, dpath in sample_glist:
            prot_best_regime: dict[str, int] = {}
            with open(dpath) as fh:
                fh.readline()
                for line in fh:
                    if "B4_" in line or "B3_" in line or "B2_singleton" in line:
                        continue
                    p = line.rstrip("\n").split("\t")
                    if len(p) < 10:
                        continue
                    ev_type = p[8]
                    prot_id = p[1]
                    pid_val = float(p[2])
                    go_id = p[3]
                    ev_detail = p[9]

                    if pid_val == 0.0:
                        reg_idx = 0
                    elif pid_val < 30.0:
                        reg_idx = 1
                    elif pid_val < 50.0:
                        reg_idx = 2
                    else:
                        reg_idx = 3

                    prot_best_regime[prot_id] = reg_idx

                    if ev_type in ("B5_pathway_gapfill", "B2_complex_module_coherence"):
                        key = (ev_type, ev_detail, go_id)
                        slot_counts[key][reg_idx] += scale
                        if reg_idx in (0, 1):
                            slot_phyla[key].add(phy)
                    elif ev_type == "B6_operon_rescue":
                        key = ("B6_operon_rescue", "Operon_Synteny_Slot", go_id)
                        slot_counts[key][reg_idx] += scale
                        if reg_idx in (0, 1):
                            slot_phyla[key].add(phy)

            for _, r_idx in prot_best_regime.items():
                phylum_nogd[phy][r_idx] += scale

    return {
        "slot_counts": dict(slot_counts),
        "slot_phyla": {k: list(v) for k, v in slot_phyla.items()},
        "phylum_nogd": dict(phylum_nogd),
        "phylum_orf_regimes": dict(phylum_orf_regimes),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--insight3-tsv",
        type=Path,
        default=Path("results/allthebacteria/insight3_pathway_and_essential_completeness.tsv"),
    )
    ap.add_argument(
        "--out-bifurcation-tsv",
        type=Path,
        default=Path("results/allthebacteria/insight6_genome_size_bifurcation.tsv"),
    )
    ap.add_argument("--work-dirs", type=Path, nargs="*", default=[])
    ap.add_argument("--batch-manifests", type=Path, nargs="*", default=[])
    ap.add_argument(
        "--out-nogd-tsv",
        type=Path,
        default=Path("results/allthebacteria/insight5_non_orthologous_displacement.tsv"),
    )
    ap.add_argument(
        "--out-phylum-nogd-tsv",
        type=Path,
        default=Path("results/allthebacteria/insight5b_phylum_nogd_vs_divergent.tsv"),
    )
    ap.add_argument("--workers", type=int, default=min(24, os.cpu_count() or 8))
    args = ap.parse_args()

    if args.insight3_tsv.exists():
        compute_genome_size_bifurcation(args.insight3_tsv, args.out_bifurcation_tsv)

    if args.work_dirs:
        t0 = time.time()
        sample_meta: dict[str, dict[str, str]] = {}
        for bm in args.batch_manifests:
            if not bm.exists():
                continue
            with open(bm) as fh:
                hdr = fh.readline().rstrip("\n").split("\t")
                for line in fh:
                    if not line.strip():
                        continue
                    vals = line.rstrip("\n").split("\t")
                    row = {hdr[i]: vals[i] if i < len(vals) else "" for i in range(len(hdr))}
                    sample_meta[row["sample"]] = row

        valid_wdirs = [str(w) for w in args.work_dirs if w.exists()]
        print(f"Scanning {len(valid_wdirs)} batch directories across {args.workers} workers...")

        total_slot_counts: dict[tuple[str, str, str], list[float]] = defaultdict(lambda: [0.0, 0.0, 0.0, 0.0])
        total_slot_phyla: dict[tuple[str, str, str], set[str]] = defaultdict(set)
        total_phylum_nogd: dict[str, list[float]] = defaultdict(lambda: [0.0, 0.0, 0.0, 0.0])
        total_phylum_orfs: dict[str, list[int]] = defaultdict(lambda: [0, 0, 0, 0])

        with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as ex:
            futs = [ex.submit(scan_batch_for_nogd, w, sample_meta) for w in valid_wdirs]
            for fut in concurrent.futures.as_completed(futs):
                res = fut.result()
                if not res:
                    continue
                for k, cnts in res["slot_counts"].items():  # type: ignore[union-attr]
                    for i in range(4):
                        total_slot_counts[k][i] += cnts[i]
                for k, plist in res["slot_phyla"].items():  # type: ignore[union-attr]
                    total_slot_phyla[k].update(plist)
                for phy, cnts in res["phylum_nogd"].items():  # type: ignore[union-attr]
                    for i in range(4):
                        total_phylum_nogd[phy][i] += cnts[i]
                for phy, ocnts in res["phylum_orf_regimes"].items():  # type: ignore[union-attr]
                    for i in range(4):
                        total_phylum_orfs[phy][i] += ocnts[i]

        rows = []
        for (ev_type, ev_detail, go_id), cnts in total_slot_counts.items():
            c0, c1, c2, c3 = int(round(cnts[0])), int(round(cnts[1])), int(round(cnts[2])), int(round(cnts[3]))
            tot_dark = c0 + c1
            tot_all = c0 + c1 + c2 + c3
            if tot_dark < 20 and ev_type != "B5_pathway_gapfill":
                continue
            desc = PATHWAY_MODULE_NAMES.get(ev_detail, ev_detail)
            go_desc = GO_NAMES.get(go_id, GO_NAMES.get(ev_detail.replace("group=", ""), go_id))
            rows.append(
                (
                    ev_type != "B5_pathway_gapfill",
                    -tot_dark,
                    [
                        ev_type,
                        ev_detail,
                        desc,
                        go_id,
                        go_desc,
                        str(c0),
                        str(c1),
                        str(tot_dark),
                        f"{100.0 * c0 / max(1, tot_dark):.2f}",
                        str(c2),
                        str(c3),
                        str(tot_all),
                        str(len(total_slot_phyla[(ev_type, ev_detail, go_id)])),
                    ],
                )
            )

        rows.sort(key=lambda r: (r[0], r[1]))
        hdr_nogd = [
            "mechanism",
            "pathway_or_complex_slot",
            "slot_description",
            "go_term",
            "go_description",
            "n_zero_align_nogd_pident0",
            "n_remote_divergent_0_to_30",
            "n_total_dark_rescued_lt30",
            "pct_nogd_of_dark",
            "n_twilight_rescued_30_50",
            "n_conserved_ge50",
            "n_total_rescued",
            "n_phyla_with_dark_rescue",
        ]
        with open(args.out_nogd_tsv, "w") as out:
            out.write("\t".join(hdr_nogd) + "\n")
            for _, _, rvals in rows:
                out.write("\t".join(rvals) + "\n")

        hdr_phy = [
            "phylum",
            "n_zero_align_orphan_orfs",
            "n_remote_divergent_0_30_orfs",
            "n_total_dark_orfs_lt30",
            "pct_orphan_of_dark_orfs",
            "n_zero_align_nogd_rescued_proteins",
            "n_remote_divergent_0_30_rescued_proteins",
            "n_total_dark_rescued_proteins",
            "pct_nogd_among_dark_rescued",
            "n_twilight_rescued_proteins",
        ]
        with open(args.out_phylum_nogd_tsv, "w") as out:
            out.write("\t".join(hdr_phy) + "\n")
            for phy, ocnts in sorted(total_phylum_orfs.items(), key=lambda kv: -(kv[1][0] + kv[1][1])):
                rcnts = total_phylum_nogd[phy]
                o0, o1 = ocnts[0], ocnts[1]
                otot = o0 + o1
                r0, r1, r2 = int(round(rcnts[0])), int(round(rcnts[1])), int(round(rcnts[2]))
                rtot = r0 + r1
                out.write(
                    "\t".join(
                        [
                            phy,
                            str(o0),
                            str(o1),
                            str(otot),
                            f"{100.0 * o0 / max(1, otot):.2f}",
                            str(r0),
                            str(r1),
                            str(rtot),
                            f"{100.0 * r0 / max(1, rtot):.2f}",
                            str(r2),
                        ]
                    )
                    + "\n"
                )
        print(f"Completed NOGD scan in {time.time() - t0:.1f}s -> {args.out_nogd_tsv} and {args.out_phylum_nogd_tsv}")


if __name__ == "__main__":
    main()
