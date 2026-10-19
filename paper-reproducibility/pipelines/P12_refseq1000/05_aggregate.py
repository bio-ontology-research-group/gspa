#!/usr/bin/env python3
"""P12 step 05: aggregate per-genome outputs into compact per-genome TSVs and
the 8-row Suppl. Table S8.

Reads under --root (IBEX run directory):
  results/matrix_rows.json            p12_matrix_rows.py (matrix reproduction, 6 tools)
  eval/{tools,gspa}/suite/*.jsonl      p12_structural_suite.py (one evaluator, 8 rows)
  eval/{tools,gspa}/gaef/*_quality.json  gspa evaluate --lite (GSPA evaluator, 8 rows)
  work/<genome>/times/*.time            /usr/bin/time -v per GSPA step
  work/<genome>/DONE
Optional --truth-summary (04_goa_exp_truth.py exp_truth_summary.tsv).

Writes to --out-dir:
  per_genome_suite.tsv, per_genome_gaef.tsv, per_genome_runtime.tsv,
  matrix_reproduction.tsv, table_s8_regenerated.tsv (8 rows, column `source`
  = computed|unavailable per row and `unavailable_columns`), progress.json
"""
import argparse
import glob
import json
import os
import re
import statistics
from collections import defaultdict

LABEL = {
    "deepgometa": "DeepGO-SE / DeepGOMeta (2024)", "sprofgo": "SPROF-GO (2023)",
    "deepgraphgo": "DeepGraphGO (2021)", "deepfri": "DeepFRI (2021)", "tale": "TALE (2021)",
    "IPscan": "InterProScan 5", "gspa_stage0": "GSPA Stage 0 (Noisy-OR)",
    "gspa_stage12": "GSPA Stage 1+2 (Full)"}
ORDER = list(LABEL)
STAGE0_STEPS = ["diamond", "hmmsearch", "parse_claims", "parse_ipr", "integrate_stage0"]
STAGE12_STEPS = ["diamond", "hmmsearch", "parse_claims", "parse_ipr", "integrate_stage12"]


def parse_time(path):
    wall = rss = None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("Elapsed (wall clock) time"):
                v = line.rsplit(": ", 1)[1]
                parts = [float(x) for x in v.split(":")]
                wall = sum(p * 60 ** i for i, p in enumerate(reversed(parts)))
            elif line.startswith("Maximum resident set size"):
                rss = int(line.rsplit(": ", 1)[1]) / 1024.0
            elif line.startswith("Exit status"):
                if int(line.rsplit(": ", 1)[1]) != 0:
                    wall = None
    return wall, rss


def mean(xs):
    xs = [x for x in xs if x is not None]
    return round(statistics.mean(xs), 4) if xs else None


def median(xs):
    xs = [x for x in xs if x is not None]
    return round(statistics.median(xs), 4) if xs else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default="/ibex/scratch/projects/c2014/hohndor/refseq1000_repro_2026-10-19")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--truth-summary", default=None)
    ap.add_argument("--min-genomes", type=int, default=1000,
                    help="a row/column group is reported only when this many genomes are evaluated; "
                         "otherwise it is NA and source=unavailable (partial values stay in the per-genome TSVs)")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    genomes = [x.strip() for x in open(f"{a.root}/genomes.txt") if x.strip()]

    # --- unified structural suite ---------------------------------------------
    suite = {}
    for fn in glob.glob(f"{a.root}/eval/*/suite/*.jsonl"):
        for line in open(fn):
            d = json.loads(line)
            if d.get("status") == "ok":
                suite[(d["genome"], d["method"])] = d
    keys = ["n_proteins_annotated", "n_rows_propagated", "n_rows_leaf", "essential_core22_present",
            "essential_full32_present", "metacyc_complete", "metacyc_annotated", "metacyc_pct",
            "complex_coherent", "complex_incoherent", "complex_coherence_pct", "has_part_missing_pct",
            "taxon_constrained_rows_propagated", "taxon_violating_rows_propagated",
            "taxon_violating_rows_leaf_direct", "taxon_violating_rows_leaf_implied",
            "taxon_violating_proteins"]
    with open(f"{a.out_dir}/per_genome_suite.tsv", "w") as out:
        out.write("genome\tmethod\t" + "\t".join(keys) + "\n")
        for (g, m), d in sorted(suite.items()):
            out.write(f"{g}\t{m}\t" + "\t".join(str(d.get(k, "")) for k in keys) + "\n")

    # --- GSPA evaluator -------------------------------------------------------
    gaef = {}
    for fn in glob.glob(f"{a.root}/eval/*/gaef/*_quality.json"):
        b = os.path.basename(fn)[:-len("_quality.json")]
        g, m = b.rsplit(".", 1)
        try:
            d = json.load(open(fn))
        except Exception:
            continue
        comp = d.get("completeness", {})
        coh = d.get("coherence", {})
        cons = d.get("consistency", {})
        gaef[(g, m)] = {
            "gaef_completeness_pct": 100.0 * comp.get("score", 0.0),
            "gaef_essential_present": comp.get("present_count"),
            "gaef_complex_coherence_pct": 100.0 * coh.get("complex_coherence", 0.0),
            "gaef_pathway_coherence_pct": 100.0 * coh.get("pathway_coherence", 0.0),
            "gaef_consistency_violations": cons.get("violation_count"),
            "gaef_consistent": cons.get("consistent"),
            "gaef_composite": d.get("summary", {}).get("composite_score"),
            "gaef_mean_ic": d.get("information_content", {}).get("mean_ic"),
        }
    gkeys = ["gaef_completeness_pct", "gaef_essential_present", "gaef_complex_coherence_pct",
             "gaef_pathway_coherence_pct", "gaef_consistency_violations", "gaef_consistent",
             "gaef_composite", "gaef_mean_ic"]
    with open(f"{a.out_dir}/per_genome_gaef.tsv", "w") as out:
        out.write("genome\tmethod\t" + "\t".join(gkeys) + "\n")
        for (g, m), d in sorted(gaef.items()):
            out.write(f"{g}\t{m}\t" + "\t".join(str(d.get(k, "")) for k in gkeys) + "\n")

    # --- runtimes -------------------------------------------------------------
    rt = {}
    with open(f"{a.out_dir}/per_genome_runtime.tsv", "w") as out:
        out.write("genome\tn_proteins\tstep\twall_s\tmax_rss_mb\n")
        for g in genomes:
            done = f"{a.root}/work/{g}/DONE"
            if not os.path.exists(done):
                continue
            nprot = open(done).read().split("\t")[2]
            steps = {}
            for tf in glob.glob(f"{a.root}/work/{g}/times/*.time"):
                s = os.path.basename(tf)[:-5]
                steps[s] = parse_time(tf)
                out.write(f"{g}\t{nprot}\t{s}\t{steps[s][0]}\t{steps[s][1]}\n")
            rt[g] = steps
    def total(steps, names):
        if any(n not in steps or steps[n][0] is None for n in names):
            return None, None
        return sum(steps[n][0] for n in names), max(steps[n][1] for n in names)

    # --- matrix reproduction ----------------------------------------------------
    mr = {}
    mpath = f"{a.root}/results/matrix_rows.json"
    if os.path.exists(mpath):
        mr = json.load(open(mpath))
        cols = sorted({k for v in mr.values() for k in v if not isinstance(v[k], list)})
        with open(f"{a.out_dir}/matrix_reproduction.tsv", "w") as out:
            out.write("method\t" + "\t".join(cols) + "\n")
            for lab, v in mr.items():
                out.write(lab + "\t" + "\t".join(str(v.get(c, "")) for c in cols) + "\n")

    # --- truth ------------------------------------------------------------------
    truth_note = "no truth summary supplied"
    if a.truth_summary and os.path.exists(a.truth_summary):
        rows = [l.rstrip("\n").split("\t") for l in open(a.truth_summary)][1:]
        n_same = sum(1 for r in rows if int(r[1]) > 0)
        n_lin = sum(1 for r in rows if int(r[2]) > 0)
        n_any = sum(1 for r in rows if int(r[3]) > 0)
        truth_note = (f"genomes with >=1 experimentally annotated protein: same assembly taxid {n_same}, "
                      f"lineage<=genus {n_lin}, any taxon {n_any} (of {len(rows)})")

    # --- table -------------------------------------------------------------------
    cols = ["method", "n_genomes", "cafa_macro_fmax", "essential_core22_pct", "essential_full29_pct",
            "complex_coherence_pct", "metacyc_completeness_pct", "has_part_coherence_pct",
            "taxon_constrained_rows_per_genome", "runtime_sec_per_genome", "peak_ram_mb",
            "source", "n_genomes_evaluated",
            "essential_full_n_terms", "taxon_violating_rows_propagated_per_genome",
            "taxon_violating_rows_leaf_per_genome", "taxon_violating_leaf_rows_implied_per_genome",
            "taxon_violating_proteins_per_genome", "smin",
            "gaef_completeness_pct", "gaef_complex_coherence_pct", "gaef_pathway_coherence_pct",
            "gaef_consistency_violations_per_genome", "gaef_n_genomes",
            "runtime_sec_per_genome_median", "peak_ram_mb_max", "unavailable_columns", "notes"]
    progress = {}
    with open(f"{a.out_dir}/table_s8_regenerated.tsv", "w") as out:
        out.write("\t".join(cols) + "\n")
        for m in ORDER:
            S = [suite[(g, m)] for g in genomes if (g, m) in suite]
            G = [gaef[(g, m)] for g in genomes if (g, m) in gaef]
            progress[m] = {"suite": len(S), "gaef": len(G)}
            r = {c: None for c in cols}
            r["method"] = LABEL[m]
            r["n_genomes"] = len(genomes)
            r["n_genomes_evaluated"] = len(S)
            if len(S) >= a.min_genomes:
                r["essential_core22_pct"] = mean([100.0 * d["essential_core22_present"] / 22 for d in S])
                r["essential_full29_pct"] = mean([100.0 * d["essential_full32_present"] / 32 for d in S])
                r["essential_full_n_terms"] = 32
                r["complex_coherence_pct"] = mean([d["complex_coherence_pct"] for d in S])
                r["metacyc_completeness_pct"] = mean([d["metacyc_pct"] for d in S])
                r["has_part_coherence_pct"] = mean([100.0 - d["has_part_missing_pct"] for d in S])
                r["taxon_constrained_rows_per_genome"] = mean([d["taxon_violating_rows_propagated"] for d in S])
                r["taxon_violating_rows_propagated_per_genome"] = r["taxon_constrained_rows_per_genome"]
                r["taxon_violating_rows_leaf_per_genome"] = mean([d["taxon_violating_rows_leaf_direct"] for d in S])
                r["taxon_violating_leaf_rows_implied_per_genome"] = mean([d["taxon_violating_rows_leaf_implied"] for d in S])
                r["taxon_violating_proteins_per_genome"] = mean([d["taxon_violating_proteins"] for d in S])
            r["gaef_n_genomes"] = len(G)
            if len(G) >= a.min_genomes:
                r["gaef_completeness_pct"] = mean([d["gaef_completeness_pct"] for d in G])
                r["gaef_complex_coherence_pct"] = mean([d["gaef_complex_coherence_pct"] for d in G])
                r["gaef_pathway_coherence_pct"] = mean([d["gaef_pathway_coherence_pct"] for d in G])
                r["gaef_consistency_violations_per_genome"] = mean([d["gaef_consistency_violations"] for d in G])
            notes = []
            if m.startswith("gspa"):
                names = STAGE0_STEPS if m == "gspa_stage0" else STAGE12_STEPS
                tots = [total(rt[g], names) for g in rt]
                tots = [t for t in tots if t[0] is not None]
                if len(tots) >= a.min_genomes:
                    r["runtime_sec_per_genome"] = mean([t[0] for t in tots])
                    r["runtime_sec_per_genome_median"] = median([t[0] for t in tots])
                    r["peak_ram_mb"] = mean([t[1] for t in tots])
                    r["peak_ram_mb_max"] = round(max(t[1] for t in tots), 1)
                notes.append(f"runtime n={len(tots)}; wall s of DIAMOND+HMMER+claim parsing+integrate on 4 cores; "
                             "InterProScan reused (not timed); peak_ram_mb = mean over genomes of max RSS over steps")
                notes.append("posterior>=0.50")
            else:
                notes.append("tool-specific score cut-offs of the annot_evaluation audit; runtime/RAM not logged by the audit")
            notes.append("cafa/smin: no experimental truth for most genomes; " + truth_note)
            r["notes"] = "; ".join(notes)
            una = [c for c in ("cafa_macro_fmax", "smin", "runtime_sec_per_genome", "peak_ram_mb")
                   if r[c] is None]
            r["unavailable_columns"] = ",".join(una)
            r["source"] = "computed" if len(S) >= a.min_genomes else "unavailable"
            if len(S) < a.min_genomes:
                r["notes"] += f"; only {len(S)} of {len(genomes)} genomes evaluated at aggregation time"
            una = [c for c in cols[2:11] + ["smin"] if r[c] is None]
            r["unavailable_columns"] = ",".join(una)
            out.write("\t".join("NA" if r[c] is None else str(r[c]) for c in cols) + "\n")
    json.dump({"genomes": len(genomes), "gspa_done": len(rt), "per_method": progress},
              open(f"{a.out_dir}/progress.json", "w"), indent=2)
    print(json.dumps({"gspa_done": len(rt), "per_method": progress}))


if __name__ == "__main__":
    main()
