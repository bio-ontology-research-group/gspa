#!/usr/bin/env python3
"""P03 final step (local): turn fetched raw P03 outputs into the two result tables.

Inputs (fetched from unimatrix01 by run.sh into $P03_FETCH, default ./fetched):
  panel13_summary.json, panel13_per_genome.tsv, panel13_audit.tsv,
  panel13_perturbation.tsv, panel21_summary.json, panel21_per_genome.tsv
Outputs:
  results/pillar1_suite2_panel13_21_sota.tsv
  results/pillar1_panel11_structural_audit.tsv
"""
import csv
import json
import os
import statistics as st
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
FETCH = os.environ.get("P03_FETCH", os.path.join(HERE, "fetched"))

LABEL13 = {
    "pgap": "NCBI PGAP (Li et al. 2021)",
    "pfam": "HMMER 3.4 / Pfam-A + pfam2go",
    "interproscan": "InterProScan 5.72-103.0",
    "diamond": "DIAMOND-BLASTp",
    "mdf": "metagenomic-deepFRI (Bezshapkin et al. 2026)",
    "gspa3": "GSPA Stage 1+2A (3-Channel D+P+IPR)",
}
LABEL21 = {
    "esm2-deepgoplus": "ESM2-DeepGOPlus (2024)",
    "proteinfer": "ProteInfer (Sanderson et al. 2023)",
    "clean": "CLEAN Contrastive EC (Yu et al. 2023)",
    "diamond": "DIAMOND-BLASTp",
    "pfam": "HMMER 3.4 / Pfam-A + pfam2go",
    "gspa2": "GSPA Stage 1+2A (DIAMOND-only Stage 0)",
    "esm2t30-knn": "ESM-2 t30 150M kNN (k=10, leave-one-genome-out over 30-genome panel)",
}
UNAVAILABLE = [
    ("Panel21_All21", "HMMER 3.4 / Pfam-A + pfam2go", "bench_gtdb30/preds/{tag}.pfam.tsv hold ~20 lines (search truncated); 15 GO pairs per genome - not a usable channel"),
    ("Panel21_All21", "InterProScan 5.72-103.0", "no InterProScan output exists for the 21-genome panel"),
    ("Panel21_All21", "FoldSeek + ProstT5 3Di (2024)", "never ran: run logs report the afdb_swissprot database missing; no prediction tables"),
    ("Panel21_All21", "ProtT5 kNN (Elnaggar et al. 2021)", "no ProtT5 embeddings or prediction tables found on unimatrix01"),
    ("Panel21_All21", "ESM-2 650M kNN (Lin et al. 2023)", "no ESM-2 650M embeddings or prediction tables found (only ESM-2 t30 150M embeddings exist; see the t30 row)"),
    ("Panel21_All21", "GSPA Stage 1+2 (3-Channel Core)", "requires InterProScan claims, which do not exist for the 21 genomes (2-channel D+P row computed instead)"),
    ("Panel21_All21", "GSPA Stage 1+2 (7-Channel All-Pooled)", "no 7-channel claim sets or GSPA run on the 21-genome panel were found"),
    ("Panel9_Sub9", "gLM / gLM2 Genomic LM (Hwang et al. 2024)", "gLM/gLM2 outputs on disk (gspa-glm/phase1/preds*) are operon co-membership predictions used as GSPA genomic-context priors, not GO predictions; only GSPA integrations with those priors exist (gspa-glm/full_priors_glm{,2}, 8 genomes, not 9) - no standalone gLM GO annotation to evaluate"),
]
COLS = ["panel", "method", "n_genomes", "go_fmax_micro_all", "go_fmax_cafa_all", "go_fmax_strict_holdout",
        "go_fmax_exp_prop", "go_fmax_strict_leaf", "go_smin_all", "go_smin_exp_prop", "ec_fmax_micro",
        "ec_fmax_cafa", "runtime_sec_per_genome", "peak_ram_mb",
        "go_fmax_cafa_exp_prop", "ic_recall_exp", "coverage_exp", "ic_recall_all", "coverage_all",
        "n_exp_truth_proteins", "go_fmax_micro_all_pooled", "go_fmax_cafa_all_pooled",
        "source", "n_genomes_with_output", "notes"]


def f4(x):
    return "" if x is None else f"{x:.4f}"


def row_from_summary(panel, label, s, ntot, ec_key=None, note=""):
    r = {c: "" for c in COLS}
    n_all = s.get("all_n_genomes", 0)
    r.update({"panel": panel, "method": label, "n_genomes": ntot, "source": "computed",
              "n_genomes_with_output": n_all or s.get(f"{ec_key}_n_genomes", 0) if ec_key else n_all})
    if n_all:
        r["go_fmax_micro_all"] = f4(s.get("all_mean_fmax_micro"))
        r["go_fmax_cafa_all"] = f4(s.get("all_mean_fmax_cafa"))
        r["go_smin_all"] = f"{s['all_mean_smin']:.2f}" if s.get("all_mean_smin") is not None else ""
        r["ic_recall_all"] = f4(s.get("all_mean_ic_recall_05"))
        r["coverage_all"] = f4(s.get("all_mean_coverage_05"))
        r["go_fmax_micro_all_pooled"] = f4(s.get("all_pooled_fmax_micro"))
        r["go_fmax_cafa_all_pooled"] = f4(s.get("all_pooled_fmax_cafa"))
    if s.get("exp_n_genomes"):
        r["go_fmax_exp_prop"] = f4(s.get("exp_pooled_fmax_micro"))
        r["go_fmax_cafa_exp_prop"] = f4(s.get("exp_pooled_fmax_cafa"))
        r["go_smin_exp_prop"] = f"{s['exp_pooled_smin']:.2f}" if s.get("exp_pooled_smin") is not None else ""
        r["ic_recall_exp"] = f4(s.get("exp_pooled_ic_recall_05"))
        r["coverage_exp"] = f4(s.get("exp_pooled_coverage_05"))
        r["n_exp_truth_proteins"] = s.get("exp_pooled_n_truth_proteins", "")
    if ec_key and s.get(f"{ec_key}_n_genomes"):
        r["ec_fmax_micro"] = f4(s.get(f"{ec_key}_pooled_fmax_micro"))
        r["ec_fmax_cafa"] = f4(s.get(f"{ec_key}_pooled_fmax_cafa"))
    r["notes"] = note
    return r


def main():
    s13 = json.load(open(os.path.join(FETCH, "panel13_summary.json")))
    rows = []
    base_note = ("GO all-evidence metrics = mean of per-genome values; *_exp_prop and exp IC recall/coverage = "
                 "pooled over genomes (per-genome experimental truth holds 1-277 rows); ancestor propagation "
                 "(is_a+part_of, roots removed); 13 prokaryotes (no eukaryote in the panel)")
    for s in s13:
        rows.append(row_from_summary("Panel13_Prok13", LABEL13[s["method"]], s, 13, note=base_note))
    p21 = os.path.join(FETCH, "panel21_summary.json")
    if os.path.exists(p21):
        for s in json.load(open(p21)):
            if s["method"] == "pfam":
                continue  # bench_gtdb30 hmmsearch output truncated (15 GO pairs/genome); reported as unavailable
            note = "EC = exact 4-digit match, pooled, vs Swiss-Prot EC of panel proteins; " + base_note.replace(
                "13 prokaryotes (no eukaryote in the panel)", "21 bacteria (no eukaryote in the panel)")
            if s["method"] == "esm2t30-knn":
                note = "ESM-2 t30 (150M), NOT 650M; reference = 29 other genomes of bench_gtdb30 with GOA truth; " + note
            if s["method"] == "gspa2":
                note = ("Stage 0 = bench_gtdb30 D+P claims, whose Pfam part is truncated (15 GO pairs/genome), "
                        "so effectively DIAMOND-only; no InterProScan output for these genomes; " + note)
            rows.append(row_from_summary("Panel21_All21", LABEL21[s["method"]], s, 21, ec_key="ec_sprot", note=note))
    for panel, label, why in UNAVAILABLE:
        r = {c: "" for c in COLS}
        r.update({"panel": panel, "method": label, "n_genomes": 21 if panel.startswith("Panel21") else 9,
                  "source": "unavailable", "n_genomes_with_output": 0, "notes": why})
        rows.append(r)
    out = os.path.join(REPO, "results", "pillar1_suite2_panel13_21_sota.tsv")
    with open(out, "w") as f:
        w = csv.DictWriter(f, fieldnames=COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print("wrote", out, len(rows))

    # ---- structural audit + perturbation (13-genome panel)
    aud_by_var = {v: list(csv.DictReader(open(os.path.join(FETCH, f"panel13_audit_{v}.tsv")), delimiter="\t"))
                  for v in ("t05", "fmaxt")}
    per = list(csv.DictReader(open(os.path.join(FETCH, "panel13_perturbation.tsv")), delimiter="\t"))
    acols = ["panel", "method", "audit_threshold", "n_genomes", "n_genomes_with_output", "essential_completeness",
             "pathway_coherence", "complex_coherence", "violations_mean", "violations_total",
             "n_genomes_consistent", "composite_audit_score", "perturb_injected_per_genome",
             "perturb_surviving_threshold_only_mean", "perturb_surviving_gspa_A9_mean",
             "perturb_surviving_gspa_A9_total", "perturb_suppression_gspa_A9_pct", "source", "notes"]
    arows = []
    for var, m in [(v, m) for v in ("t05", "fmaxt") for m in ["pgap", "pfam", "interproscan", "diamond", "mdf", "gspa3"]]:
        aud = aud_by_var[var]
        a = [x for x in aud if x["method"] == m and x["composite_score"] != ""]
        p = [x for x in per if x["method"] == m]
        mean = lambda k: st.mean(float(x[k]) for x in a) if a else None
        pb = [int(x["surv_gspa_A9"]) for x in p if x["surv_gspa_A9"] != ""]
        pa = [int(x["surv_threshold_only"]) for x in p]
        ninj = int(p[0]["n_injected"]) if p else 0
        arows.append({
            "panel": "Panel13_Prok13", "method": LABEL13[m],
            "audit_threshold": "0.5" if var == "t05" else "per-genome micro F-max threshold (all-evidence truth)",
            "n_genomes": 13, "n_genomes_with_output": len(a),
            "essential_completeness": f4(mean("essential_completeness")),
            "pathway_coherence": f4(mean("pathway_coherence")),
            "complex_coherence": f4(mean("complex_coherence")),
            "violations_mean": f"{mean('violation_count'):.2f}" if a else "",
            "violations_total": sum(int(x["violation_count"]) for x in a) if a else "",
            "n_genomes_consistent": sum(1 for x in a if x["consistent"] == "True"),
            "composite_audit_score": f4(mean("composite_score")),
            "perturb_injected_per_genome": ninj,
            "perturb_surviving_threshold_only_mean": f"{st.mean(pa):.2f}" if pa else "",
            "perturb_surviving_gspa_A9_mean": f"{st.mean(pb):.2f}" if pb else "",
            "perturb_surviving_gspa_A9_total": sum(pb) if pb else "",
            "perturb_suppression_gspa_A9_pct": f"{100 * (1 - sum(pb) / (ninj * len(pb))):.1f}" if pb else "",
            "source": "computed",
            "notes": ("audit: gspa evaluate --lite on annotations with score>=0.5 (make_synth_gff_gaf.py); "
                      "perturbation design v2 (100 cross-clade claims/genome, same pairs for all methods); "
                      "threshold_only = injected max raw score>=0.5 without GSPA; gspa_A9 = posterior>=0.5 after "
                      "gspa integrate A9 (Stage 1 priors + Stage 2A relax) on the method's claims"
                      + ("; GSPA row injection adds a corroborating 2nd-source claim for 60% of pairs as in "
                         "04_run_perturbation.py" if m == "gspa3" else "")),
        })
    out2 = os.path.join(REPO, "results", "pillar1_panel11_structural_audit.tsv")
    with open(out2, "w") as f:
        w = csv.DictWriter(f, fieldnames=acols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(arows)
    print("wrote", out2, len(arows))


if __name__ == "__main__":
    main()
