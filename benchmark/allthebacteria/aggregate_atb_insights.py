#!/usr/bin/env python3
"""
Stage D: Cross-Phylum Biological Insights Aggregator (`aggregate_atb_insights.py`).

Uses multiprocessing across batch directories and reads `dedup/genome_protein_to_rep.tsv` +
`stageB/reps_vs_sprot.m8` per batch (instead of opening thousands of per-genome homology files
over Lustre) to aggregate 50,000+ genomes in seconds.

Produces five publication-ready TSV tables plus a unified JSON summary in `--out-dir`:
  1. `insight1_dark_matter_rescue_by_phylum.tsv`
  2. `insight2_complex_stoichiometry_and_operons.tsv`
  3. `insight3_pathway_and_essential_completeness.tsv`
  3b. `insight3b_phylum_summary.tsv`
  4. `insight4_taxon_anatomical_suppression.tsv`
"""
from __future__ import annotations

import argparse
import concurrent.futures
import json
import os
import time
from collections import defaultdict
from pathlib import Path


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] [StageD-Insights] {msg}", flush=True)


def pident_bin(pident: float) -> str:
    if pident < 30.0:
        return "dark_matter_lt30"
    if pident < 50.0:
        return "twilight_30_50"
    if pident < 70.0:
        return "moderate_50_70"
    return "close_homolog_ge70"


def process_single_batch(
    wdir_str: str,
    sample_meta: dict[str, dict[str, str]],
) -> dict[str, object] | None:
    wdir = Path(wdir_str)
    out_dir = wdir / "stageC" / "out"
    if not out_dir.exists():
        return None

    completed_gids: dict[str, Path] = {}
    with os.scandir(out_dir) as it:
        for entry in it:
            if entry.name.endswith(".delta_events.tsv") and entry.is_file(follow_symlinks=False):
                gid = entry.name[: -len(".delta_events.tsv")]
                completed_gids[gid] = Path(entry.path)
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

    # 2. Read dedup/genome_protein_to_rep.tsv (1 file per batch) to count ORFs per homology bin per genome
    genome_orf_counts: dict[str, list[int]] = defaultdict(lambda: [0, 0, 0, 0])  # [lt30, 30_50, 50_70, ge70]
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
                    rid = p[3]
                    pid_val = rep_pident.get(rid, 0.0)
                    cnts = genome_orf_counts[gid]
                    if pid_val < 30.0:
                        cnts[0] += 1
                    elif pid_val < 50.0:
                        cnts[1] += 1
                    elif pid_val < 70.0:
                        cnts[2] += 1
                    else:
                        cnts[3] += 1

    # 3. Load Java batch_stats.tsv if present
    java_stats: dict[str, dict[str, str]] = {}
    bstats_p = wdir / "stageC" / "batch_stats.tsv"
    if bstats_p.exists():
        with open(bstats_p) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                if not line.strip():
                    continue
                vals = line.rstrip("\n").split("\t")
                row = {hdr[i]: vals[i] if i < len(vals) else "" for i in range(len(hdr))}
                java_stats[row["genome_id"]] = row

    phylum_bin_prots: dict[tuple[str, str, str], int] = defaultdict(int)
    phylum_bin_rescued_count: dict[tuple[str, str, str], int] = defaultdict(int)
    phylum_bin_events: dict[tuple[str, str, str], dict[str, int]] = defaultdict(lambda: defaultdict(int))
    complex_stats: dict[tuple[str, str], dict[str, int]] = defaultdict(lambda: defaultdict(int))
    violation_stats: dict[tuple[str, str], dict[str, int]] = defaultdict(lambda: defaultdict(int))
    genome_rows: list[dict[str, object]] = []

    bins_order = ["dark_matter_lt30", "twilight_30_50", "moderate_50_70", "close_homolog_ge70"]

    for gid, delta_p in completed_gids.items():
        meta = sample_meta.get(gid, {})
        dom = meta.get("domain", "d__Bacteria")
        phy = meta.get("phylum", "p__Pseudomonadota")
        sp = meta.get("species", gid)
        taxid = meta.get("ncbi_taxid", "NCBITaxon:2")

        cnts = genome_orf_counts[gid]
        n_orfs = sum(cnts)
        n_dark_orfs = cnts[0]
        n_twilight_orfs = cnts[1]
        for idx, hb in enumerate(bins_order):
            phylum_bin_prots[(dom, phy, hb)] += cnts[idx]

        ev_counts: dict[str, int] = defaultdict(int)
        rescued_prots_by_bin: dict[str, set[str]] = defaultdict(set)
        with open(delta_p) as fh:
            fh.readline()
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) < 10:
                    continue
                prot_id = p[1]
                pid_val = float(p[2])
                ev_type = p[8]
                ev_detail = p[9]
                hb = pident_bin(pid_val)
                ev_counts[ev_type] += 1
                phylum_bin_events[(dom, phy, hb)][ev_type] += 1

                if ev_type in (
                    "B6_operon_rescue",
                    "B5_pathway_gapfill",
                    "B2_complex_module_coherence",
                    "B1_essential_function_rescue",
                ):
                    rescued_prots_by_bin[hb].add(prot_id)

                if ev_type == "B2_singleton_complex_suppressed":
                    complex_stats[(phy, ev_detail)]["singleton_suppressed"] += 1
                elif ev_type == "B2_complex_module_coherence" and ev_detail.startswith("group=GO:"):
                    cname = ev_detail.replace("group=", "")
                    complex_stats[(phy, cname)]["coherence_rescued"] += 1
                elif ev_type == "B6_operon_rescue":
                    complex_stats[(phy, "Operon_CoLocalized_Loci")]["operon_rescued"] += 1

                if ev_type == "B3_taxon_violation_suppressed":
                    violation_stats[(dom, phy)]["B3_taxon_suppressed"] += 1
                elif ev_type == "B4_anatomical_violation_suppressed":
                    violation_stats[(dom, phy)]["B4_anatomical_suppressed"] += 1

        for hb, pset in rescued_prots_by_bin.items():
            phylum_bin_rescued_count[(dom, phy, hb)] += len(pset)

        violation_stats[(dom, phy)]["n_genomes"] += 1
        violation_stats[(dom, phy)]["n_orfs"] += n_orfs

        js_full = java_stats.get(gid, {})
        js_noisy = java_stats.get(f"{gid}__noisy_or", {})

        genome_rows.append(
            {
                "genome_id": gid,
                "species": sp,
                "domain": dom,
                "phylum": phy,
                "ncbi_taxid": taxid,
                "n_orfs": n_orfs,
                "n_dark_orfs_lt30": n_dark_orfs,
                "n_twilight_orfs_30_50": n_twilight_orfs,
                "n_claims_stage0": int(js_full.get("n_claims", 0)),
                "n_integrated_noisy_or": int(js_noisy.get("n_integrated", 0)),
                "n_integrated_gspa_full": int(js_full.get("n_integrated", 0)),
                "n_provenance_gspa_full": int(js_full.get("n_provenance", 0)),
                "n_metabolic_gaps_candidate": int(ev_counts["B5_pathway_gapfill"]),
                "n_B6_operon_rescues": ev_counts["B6_operon_rescue"],
                "n_B5_gapfill_rescues": ev_counts["B5_pathway_gapfill"],
                "n_B2_coherence_rescues": ev_counts["B2_complex_module_coherence"],
                "n_B1_essential_rescues": ev_counts["B1_essential_function_rescue"],
                "n_B2_singleton_suppressed": ev_counts["B2_singleton_complex_suppressed"],
                "n_B3_taxon_suppressed": ev_counts["B3_taxon_violation_suppressed"],
                "n_B4_anatomical_suppressed": ev_counts["B4_anatomical_violation_suppressed"],
                "java_runtime_sec": float(js_full.get("runtime_sec", 0.0)),
            }
        )

    lin_st = None
    lin_json = wdir / "dedup" / "linclust_stats.json"
    if lin_json.exists():
        with open(lin_json) as fh:
            lin_st = json.load(fh)

    return {
        "lin_st": lin_st,
        "phylum_bin_prots": dict(phylum_bin_prots),
        "phylum_bin_rescued_count": dict(phylum_bin_rescued_count),
        "phylum_bin_events": {k: dict(v) for k, v in phylum_bin_events.items()},
        "complex_stats": {k: dict(v) for k, v in complex_stats.items()},
        "violation_stats": {k: dict(v) for k, v in violation_stats.items()},
        "genome_rows": genome_rows,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description="Stage D: Aggregate AllTheBacteria Whole-Genome GSPA Insights")
    ap.add_argument(
        "--work-dirs",
        type=Path,
        nargs="+",
        required=True,
        help="One or more completed batch work directories",
    )
    ap.add_argument(
        "--batch-manifests",
        type=Path,
        nargs="+",
        required=True,
        help="Corresponding batch manifest TSVs mapping sample -> species, domain, phylum, ncbi_taxid",
    )
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--workers", type=int, default=min(16, os.cpu_count() or 8))
    args = ap.parse_args()

    t0 = time.time()
    out_dir: Path = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    # Load sample -> taxonomy metadata
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
    log(f"Aggregating {len(valid_wdirs)} batch work directories across {args.workers} workers...")

    phylum_bin_prots: dict[tuple[str, str, str], int] = defaultdict(int)
    phylum_bin_rescued_count: dict[tuple[str, str, str], int] = defaultdict(int)
    phylum_bin_events: dict[tuple[str, str, str], dict[str, int]] = defaultdict(lambda: defaultdict(int))
    complex_stats: dict[tuple[str, str], dict[str, int]] = defaultdict(lambda: defaultdict(int))
    violation_stats: dict[tuple[str, str], dict[str, int]] = defaultdict(lambda: defaultdict(int))
    genome_summary_rows: list[dict[str, object]] = []
    seen_genomes: set[str] = set()

    linclust_totals = {
        "n_batches": 0,
        "n_genomes": 0,
        "n_total_proteins": 0,
        "n_unique_hashes": 0,
        "n_linclust_reps": 0,
    }

    with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = [ex.submit(process_single_batch, w, sample_meta) for w in valid_wdirs]
        for fut in concurrent.futures.as_completed(futs):
            res = fut.result()
            if not res:
                continue
            lin_st = res["lin_st"]
            if isinstance(lin_st, dict):
                linclust_totals["n_batches"] += 1
                linclust_totals["n_genomes"] += int(lin_st.get("n_genomes_processed", 0))
                linclust_totals["n_total_proteins"] += int(lin_st.get("n_total_proteins", 0))
                linclust_totals["n_unique_hashes"] += int(lin_st.get("n_unique_exact_hashes", 0))
                linclust_totals["n_linclust_reps"] += int(lin_st.get("n_linclust_representatives", 0))

            for k, v in res["phylum_bin_prots"].items():  # type: ignore[union-attr]
                phylum_bin_prots[k] += v
            for k, v in res["phylum_bin_rescued_count"].items():  # type: ignore[union-attr]
                phylum_bin_rescued_count[k] += v
            for k, evs in res["phylum_bin_events"].items():  # type: ignore[union-attr]
                for ek, ev in evs.items():
                    phylum_bin_events[k][ek] += ev
            for k, cst in res["complex_stats"].items():  # type: ignore[union-attr]
                for ck, cv in cst.items():
                    complex_stats[k][ck] += cv
            for k, vst in res["violation_stats"].items():  # type: ignore[union-attr]
                for vk, vv in vst.items():
                    violation_stats[k][vk] += vv
            for gr in res["genome_rows"]:  # type: ignore[union-attr]
                gid = str(gr["genome_id"])
                if gid not in seen_genomes:
                    seen_genomes.add(gid)
                    genome_summary_rows.append(gr)

    phylum_summary: dict[tuple[str, str], dict[str, object]] = defaultdict(
        lambda: {
            "species_set": set(),
            "n_genomes": 0,
            "n_orfs": 0,
            "n_dark_orfs": 0,
            "n_twilight_orfs": 0,
            "B6_operon": 0,
            "B5_gapfill": 0,
            "B2_coherence": 0,
            "B1_essential": 0,
            "B2_singleton_supp": 0,
            "B3_taxon_supp": 0,
            "B4_anat_supp": 0,
        }
    )
    for r in genome_summary_rows:
        dom = str(r["domain"])
        phy = str(r["phylum"])
        ps = phylum_summary[(dom, phy)]
        ps["species_set"].add(str(r["species"]))  # type: ignore[attr-defined]
        ps["n_genomes"] = int(ps["n_genomes"]) + 1  # type: ignore[arg-type]
        ps["n_orfs"] = int(ps["n_orfs"]) + int(r["n_orfs"])  # type: ignore[arg-type]
        ps["n_dark_orfs"] = int(ps["n_dark_orfs"]) + int(r["n_dark_orfs_lt30"])  # type: ignore[arg-type]
        ps["n_twilight_orfs"] = int(ps["n_twilight_orfs"]) + int(r["n_twilight_orfs_30_50"])  # type: ignore[arg-type]
        ps["B6_operon"] = int(ps["B6_operon"]) + int(r["n_B6_operon_rescues"])  # type: ignore[arg-type]
        ps["B5_gapfill"] = int(ps["B5_gapfill"]) + int(r["n_B5_gapfill_rescues"])  # type: ignore[arg-type]
        ps["B2_coherence"] = int(ps["B2_coherence"]) + int(r["n_B2_coherence_rescues"])  # type: ignore[arg-type]
        ps["B1_essential"] = int(ps["B1_essential"]) + int(r["n_B1_essential_rescues"])  # type: ignore[arg-type]
        ps["B2_singleton_supp"] = int(ps["B2_singleton_supp"]) + int(r["n_B2_singleton_suppressed"])  # type: ignore[arg-type]
        ps["B3_taxon_supp"] = int(ps["B3_taxon_supp"]) + int(r["n_B3_taxon_suppressed"])  # type: ignore[arg-type]
        ps["B4_anat_supp"] = int(ps["B4_anat_supp"]) + int(r["n_B4_anatomical_suppressed"])  # type: ignore[arg-type]

    # 1. Write Insight 1: Dark-Matter & Twilight-Zone Rescue by Phylum
    i1_path = out_dir / "insight1_dark_matter_rescue_by_phylum.tsv"
    with open(i1_path, "w") as fo:
        fo.write(
            "domain\tphylum\thomology_regime\tn_proteins\tn_rescued_proteins\t"
            "B6_operon_rescue_terms\tB5_pathway_gapfill_terms\tB2_complex_coherence_terms\tB1_essential_rescue_terms\n"
        )
        for (dom, phy, hb) in sorted(phylum_bin_prots.keys()):
            n_p = phylum_bin_prots[(dom, phy, hb)]
            n_resc = phylum_bin_rescued_count[(dom, phy, hb)]
            evs = phylum_bin_events[(dom, phy, hb)]
            fo.write(
                f"{dom}\t{phy}\t{hb}\t{n_p}\t{n_resc}\t"
                f"{evs.get('B6_operon_rescue', 0)}\t{evs.get('B5_pathway_gapfill', 0)}\t"
                f"{evs.get('B2_complex_module_coherence', 0)}\t{evs.get('B1_essential_function_rescue', 0)}\n"
            )

    # 2. Write Insight 2: Complex Stoichiometry & Operon Co-Localization vs Singleton Suppression
    i2_path = out_dir / "insight2_complex_stoichiometry_and_operons.tsv"
    with open(i2_path, "w") as fo:
        fo.write("phylum\tcomplex_or_locus_group\tcoherence_rescued_subunits\toperon_rescued_loci\tsingleton_hallucinations_suppressed\n")
        for (phy, cname) in sorted(complex_stats.keys()):
            st = complex_stats[(phy, cname)]
            fo.write(
                f"{phy}\t{cname}\t{st.get('coherence_rescued', 0)}\t{st.get('operon_rescued', 0)}\t{st.get('singleton_suppressed', 0)}\n"
            )

    # 3. Write Insight 3: Per-Genome Summary
    i3_path = out_dir / "insight3_pathway_and_essential_completeness.tsv"
    cols_i3 = [
        "genome_id",
        "species",
        "domain",
        "phylum",
        "ncbi_taxid",
        "n_orfs",
        "n_dark_orfs_lt30",
        "n_twilight_orfs_30_50",
        "n_claims_stage0",
        "n_integrated_noisy_or",
        "n_integrated_gspa_full",
        "n_provenance_gspa_full",
        "n_metabolic_gaps_candidate",
        "n_B6_operon_rescues",
        "n_B5_gapfill_rescues",
        "n_B2_coherence_rescues",
        "n_B1_essential_rescues",
        "n_B2_singleton_suppressed",
        "n_B3_taxon_suppressed",
        "n_B4_anatomical_suppressed",
        "java_runtime_sec",
    ]
    with open(i3_path, "w") as fo:
        fo.write("\t".join(cols_i3) + "\n")
        for r in genome_summary_rows:
            fo.write("\t".join(str(r[c]) for c in cols_i3) + "\n")

    # 3b. Write Insight 3b: Per-Phylum Summary Table
    i3b_path = out_dir / "insight3b_phylum_summary.tsv"
    with open(i3b_path, "w") as fo:
        fo.write(
            "domain\tphylum\tn_genomes\tn_species\tn_orfs\tn_dark_orfs_lt30\tpct_dark_orfs\t"
            "n_dark_rescued_prots\tpct_dark_rescued\tn_twilight_orfs_30_50\tn_twilight_rescued_prots\tpct_twilight_rescued\t"
            "B6_operon_terms\tB5_gapfill_terms\tB2_coherence_terms\tB1_essential_terms\t"
            "B2_singleton_suppressed\tB3_taxon_suppressed\tB4_anatomical_suppressed\n"
        )
        for (dom, phy) in sorted(phylum_summary.keys(), key=lambda x: (x[0], -int(phylum_summary[x]["n_genomes"]))):  # type: ignore[arg-type]
            ps = phylum_summary[(dom, phy)]
            n_g = int(ps["n_genomes"])  # type: ignore[arg-type]
            n_sp = len(ps["species_set"])  # type: ignore[arg-type]
            n_o = int(ps["n_orfs"])  # type: ignore[arg-type]
            n_dk = int(ps["n_dark_orfs"])  # type: ignore[arg-type]
            n_tw = int(ps["n_twilight_orfs"])  # type: ignore[arg-type]
            n_dk_r = phylum_bin_rescued_count[(dom, phy, "dark_matter_lt30")]
            n_tw_r = phylum_bin_rescued_count[(dom, phy, "twilight_30_50")]
            pct_dk = 100.0 * n_dk / max(1, n_o)
            pct_dk_r = 100.0 * n_dk_r / max(1, n_dk)
            pct_tw_r = 100.0 * n_tw_r / max(1, n_tw)
            fo.write(
                f"{dom}\t{phy}\t{n_g}\t{n_sp}\t{n_o}\t{n_dk}\t{pct_dk:.2f}\t"
                f"{n_dk_r}\t{pct_dk_r:.2f}\t{n_tw}\t{n_tw_r}\t{pct_tw_r:.2f}\t"
                f"{ps['B6_operon']}\t{ps['B5_gapfill']}\t{ps['B2_coherence']}\t{ps['B1_essential']}\t"
                f"{ps['B2_singleton_supp']}\t{ps['B3_taxon_supp']}\t{ps['B4_anat_supp']}\n"
            )

    # 4. Write Insight 4: Taxonomic & Anatomical Violation Suppression Across Lineages
    i4_path = out_dir / "insight4_taxon_anatomical_suppression.tsv"
    with open(i4_path, "w") as fo:
        fo.write("domain\tphylum\tn_genomes\tn_orfs\tB3_taxon_violations_suppressed\tB4_anatomical_violations_suppressed\ttotal_violations_suppressed\n")
        for (dom, phy) in sorted(violation_stats.keys()):
            vst = violation_stats[(dom, phy)]
            b3 = vst.get("B3_taxon_suppressed", 0)
            b4 = vst.get("B4_anatomical_suppressed", 0)
            fo.write(f"{dom}\t{phy}\t{vst.get('n_genomes', 0)}\t{vst.get('n_orfs', 0)}\t{b3}\t{b4}\t{b3 + b4}\n")

    total_orfs = sum(int(r["n_orfs"]) for r in genome_summary_rows)  # type: ignore[arg-type]
    total_dark_orfs = sum(int(r["n_dark_orfs_lt30"]) for r in genome_summary_rows)  # type: ignore[arg-type]
    total_twilight_orfs = sum(int(r["n_twilight_orfs_30_50"]) for r in genome_summary_rows)  # type: ignore[arg-type]
    total_dark_rescued = sum(v for (d, p, hb), v in phylum_bin_rescued_count.items() if hb == "dark_matter_lt30")
    total_twilight_rescued = sum(v for (d, p, hb), v in phylum_bin_rescued_count.items() if hb == "twilight_30_50")
    summary_obj = {
        "linclust_totals": linclust_totals,
        "n_genomes_aggregated": len(genome_summary_rows),
        "n_distinct_phyla": len(phylum_summary),
        "n_distinct_species": len({str(r["species"]) for r in genome_summary_rows}),
        "total_orfs": total_orfs,
        "total_dark_orfs_lt30": total_dark_orfs,
        "total_dark_rescued_prots": total_dark_rescued,
        "pct_dark_rescued": round(100.0 * total_dark_rescued / max(1, total_dark_orfs), 2),
        "total_twilight_orfs_30_50": total_twilight_orfs,
        "total_twilight_rescued_prots": total_twilight_rescued,
        "pct_twilight_rescued": round(100.0 * total_twilight_rescued / max(1, total_twilight_orfs), 2),
        "total_B6_operon_rescues": sum(int(r["n_B6_operon_rescues"]) for r in genome_summary_rows),  # type: ignore[arg-type]
        "total_B5_gapfill_rescues": sum(int(r["n_B5_gapfill_rescues"]) for r in genome_summary_rows),  # type: ignore[arg-type]
        "total_B2_coherence_rescues": sum(int(r["n_B2_coherence_rescues"]) for r in genome_summary_rows),  # type: ignore[arg-type]
        "total_B1_essential_rescues": sum(int(r["n_B1_essential_rescues"]) for r in genome_summary_rows),  # type: ignore[arg-type]
        "total_B2_singleton_suppressed": sum(int(r["n_B2_singleton_suppressed"]) for r in genome_summary_rows),  # type: ignore[arg-type]
        "total_B3_taxon_suppressed": sum(int(r["n_B3_taxon_suppressed"]) for r in genome_summary_rows),  # type: ignore[arg-type]
        "total_B4_anatomical_suppressed": sum(int(r["n_B4_anatomical_suppressed"]) for r in genome_summary_rows),  # type: ignore[arg-type]
        "wall_time_sec": round(time.time() - t0, 2),
    }
    with open(out_dir / "atb_insights_summary.json", "w") as fo:
        json.dump(summary_obj, fo, indent=2)

    log(f"Wrote 5 biological insight tables and summary JSON ({len(genome_summary_rows)} genomes in {time.time() - t0:.1f}s) -> {out_dir}")


if __name__ == "__main__":
    main()
