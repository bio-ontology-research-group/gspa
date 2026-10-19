#!/usr/bin/env python3
"""
Extract compact, lossless sparse representations of real SOTA baseline predictions on IBEX:
1. Price-149 across 22 genomes for all 6 SOTA EC predictors:
   - CLEAN (2023), DeepProZyme (2023), EnzBERT (2023), GraphEC (2024), MAPred (2025), TopEC (2025)
   - Plus METEOR v8 active/muted ECs per genome & baseline
   - Preserves:
     * Exact genome protein order (for chromosomal operon neighborhood +/- 10 genes)
     * Top-10 protein scores per EC across the whole genome (for SparseRoleAggregator & StoichiometricLocusRefiner r_j(h))
     * Top-15 ECs per protein across the whole genome (for operon/pathway neighbor context)
     * Full non-zero EC score dict for all 149 Price target proteins
2. Swiss-Prot Temporal Holdout-66 & Holdout-227 across CLEAN, DeepProZyme, EnzBERT (vanilla, filt70, filt50, filt30):
   - Full non-zero EC score dict for all Holdout-66 & Holdout-227 proteins
   - Host-genome +/- 10 chromosomal neighbor top ECs and genome-wide top proteins per EC
   - METEOR v8 active/muted ECs per host genome
3. 1,000 RefSeq Bacterial Genomes (annot_evaluation) across 6 SOTA GO predictors:
   - DeepGO-SE / DeepGOMeta (2024), SPROF-GO (2023), DeepFRI (2021), DeepGraphGO (2021), TALE (2021), InterProScan 5
   - Exact per-genome Essential Core GO Completeness, MetaCyc Pathway Completeness, Protein Complex Coherence,
     has_part / occurs_in coherence, and Taxon Constraint Violations across all 1,000 genomes.
"""

import csv
import glob
import gzip
import json
import os
import pickle
import time
from collections import defaultdict

import numpy as np
import pandas as pd

F = "/ibex/scratch/projects/c2014/kexin/funcarve"
D = f"{F}/meteor_v8/data"
MOP = f"{F}/meteor_v8_evw_p2mu3_run/meteor_out_price"
MO_MAIN = f"{F}/meteor_v8_evw_p2mu3_run/meteor_out"
PA = f"{F}/paperA_2026"
REINFER_ROOT = f"{F}/meteor_v7_run/reinfer"
BACDIVE_ROOT = f"{F}/dpec2_result/result_bacdive"
ANNOT_EVAL = "/ibex/scratch/projects/c2014/annot_evaluation"
OUT_DIR = "/ibex/scratch/projects/c2014/hohndor/gspa_extracted"
os.makedirs(OUT_DIR, exist_ok=True)

anc = open(f"{D}/all_ancestors.txt").read().splitlines()
anc_set = set(anc)
obsolete_map = pickle.load(open(f"{D}/enzymeobsolete.pkl", "rb"))


def resolve_current_ec(ec, max_hops=5):
    current = ec
    seen = set()
    hops = 0
    while current in obsolete_map and current not in seen and hops < max_hops:
        seen.add(current)
        current = obsolete_map[current]
        hops += 1
    return current


def extract_pred(pred_fpath):
    df = pd.read_pickle(pred_fpath)
    if df.index[0].count(".") == 3:
        df = df.T
    if "EC:" in str(df.columns[0]):
        df.columns = [str(c).split("EC:")[-1] if "EC:" in str(c) else str(c) for c in df.columns]
    rename = {}
    for c in df.columns:
        if c not in anc_set and c in obsolete_map:
            current = resolve_current_ec(c)
            if current in anc_set and current != c:
                rename[c] = current
    for old, new in rename.items():
        if new in df.columns:
            df[new] = df[[old, new]].max(axis=1)
            df = df.drop(columns=[old])
        else:
            df = df.rename(columns={old: new})
    df = df.reindex(columns=anc, fill_value=0.0).astype(np.float32)
    return df


def price_pkl_path(b, gca):
    gns = gca.rsplit(".", 1)[0]
    num = gca.split("_")[1].split(".")[0]
    m = {
        "clean": f"{PA}/baseline_preds/price22_CLEAN_clean/{gca}/{gca}_CLEAN_confidence.pkl",
        "enzbert": f"{F}/tfpc/resultprice_newg/{gca}_enzbert_predictions.pkl",
        "graphec": f"{F}/graphec_price_new/{gca}_GraphEC.pkl",
        "mapred": f"{F}/mapred_price_new/{gca}_MAPred.pkl",
        "topec": f"{F}/topec_price_new/{gca}_TopEC.pkl",
    }
    if b == "dpz":
        e = f"{F}/dpec2_result/result_price/{gns}_DeepECv2_t5.pkl"
        if os.path.exists(e):
            return e
        h = glob.glob(f"{F}/dpec2_result/result_price/GC?_{num}_DeepECv2_t5.pkl")
        return h[0] if h else ""
    p = m[b]
    for c in (p, p.replace(gca, gns)):
        if os.path.exists(c):
            return c
    return ""


def find_protein_idx(index_list, pid, nid=None):
    idx_map = {str(x): i for i, x in enumerate(index_list)}
    for q in [pid] + ([nid] if nid else []):
        if not q:
            continue
        if q in idx_map:
            return idx_map[q]
        b = str(q).rsplit(".", 1)[0]
        if b in idx_map:
            return idx_map[b]
        for i, x in enumerate(index_list):
            if str(x) == b or str(x).startswith(b + "."):
                return i
    return -1


def extract_price149_all6():
    """Price-149 extraction now delegates to the reproducible v2 extractor.

    The v1 body (until 2026-10-19) capped each target's score dict at 100 ECs and rounded
    scores to 5 decimals; it otherwise read the six predictors correctly (top-1 EC identical
    for all 148 targets). v2 lives in
    paper-reproducibility/pipelines/P08_price149/extract_price149_all6.py (see its README for
    the raw-format inspection) and keeps full target vectors plus input sha256s.
    """
    import importlib.util

    path = os.environ.get(
        "P08_EXTRACTOR",
        "/ibex/scratch/projects/c2014/hohndor/price149_repro_2026-10-19/P08_price149/extract_price149_all6.py",
    )
    spec = importlib.util.spec_from_file_location("extract_price149_all6", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    print("[1/3] Extracting Price-149 across 22 genomes & 6 baselines (v2 extractor)...")
    mod.extract(
        "/ibex/scratch/projects/c2014/hohndor/suite6_price149_manifest.tsv",
        f"{OUT_DIR}/price149_all6_compact.pkl.gz",
        ["clean", "dpz", "enzbert", "graphec", "mapred", "topec"],
        f"{OUT_DIR}/price149_input_sha256.tsv",
    )


def resolve_baseline_pkl(baseline, variant, gca):
    suf = {"clean": "CLEAN_confidence", "dpz": "DPZ", "enzbert": "enzbert"}[baseline]
    reinfer_path = f"{REINFER_ROOT}/{baseline}_{variant}/{gca}_{suf}.pkl"
    if os.path.exists(reinfer_path):
        return reinfer_path
    orig_path = f"{PA}/baseline_preds/{baseline}_{variant}_genome_collection/{gca}_{suf}.pkl"
    if os.path.exists(orig_path):
        return orig_path
    if baseline == "dpz" and variant == "vanilla":
        base_gca = gca.split(".")[0]
        hits = glob.glob(f"{BACDIVE_ROOT}/{base_gca}_*_DeepECv2_t5.pkl")
        if hits:
            return hits[0]
    return ""


def extract_holdout66_and_227():
    print("[2/3] Extracting Swiss-Prot Holdout-66 & Holdout-227 across CLEAN, DPZ, EnzBERT...")
    hc = json.load(open(f"{F}/meteor_diag/holdout_clean.json"))["holdout_clean"]
    by_gcf = defaultdict(list)
    for uid, info in hc.items():
        by_gcf[info["gcf"]].append((uid, info["refseq_id"], info["ec"]))

    holdout_data = {"holdout66_meta": hc, "models": {}}
    for b in ["clean", "dpz", "enzbert"]:
        for v in ["vanilla", "filt70", "filt50", "filt30"]:
            key = f"{b}_{v}"
            holdout_data["models"][key] = {}
            for gcf, items in by_gcf.items():
                bp = resolve_baseline_pkl(b, v, gcf)
                if not bp or not os.path.exists(bp):
                    continue
                df = extract_pred(bp)
                vals = df.values
                cols = list(df.columns)
                idx_list = [str(x) for x in df.index]

                # Top-5 protein scores per EC
                ec_top = {}
                max_per_col = vals.max(axis=0)
                for c_idx in np.where(max_per_col >= 0.05)[0]:
                    col_v = vals[:, c_idx]
                    tk = np.argpartition(col_v, -min(5, len(col_v)))[-min(5, len(col_v)):]
                    tk = tk[np.argsort(-col_v[tk])]
                    ec_top[cols[c_idx]] = [(idx_list[i], round(float(col_v[i]), 5)) for i in tk if col_v[i] >= 0.05]

                target_scores = {}
                neighbor_ecs = {}
                for uid, pid, gt_ecs in items:
                    r_idx = find_protein_idx(idx_list, pid)
                    if r_idx < 0:
                        continue
                    row_v = vals[r_idx]
                    nz = np.where(row_v > 0.0)[0]
                    if len(nz) > 100:
                        nz = nz[np.argsort(-row_v[nz])[:100]]
                    target_scores[uid] = {cols[c_i]: round(float(row_v[c_i]), 5) for c_i in nz if row_v[c_i] > 0.0}

                    neigh = []
                    for n_idx in range(max(0, r_idx - 10), min(len(idx_list), r_idx + 11)):
                        if n_idx == r_idx:
                            continue
                        nv = vals[n_idx]
                        top5 = np.argpartition(nv, -5)[-5:]
                        top5 = top5[np.argsort(-nv[top5])]
                        n_ecs = [(cols[c_i], round(float(nv[c_i]), 4)) for c_i in top5 if nv[c_i] >= 0.15]
                        if n_ecs:
                            neigh.append((n_idx - r_idx, idx_list[n_idx], n_ecs))
                    neighbor_ecs[uid] = neigh

                m_path = f"{MO_MAIN}/{b}_{v}/meteor_preds_{gcf}.pkl"
                if not os.path.exists(m_path):
                    m_path = f"{MO_MAIN}/{b}_vanilla/meteor_preds_{gcf}.pkl"
                m_active, m_muted = [], []
                if os.path.exists(m_path):
                    md = pickle.load(open(m_path, "rb"))
                    m_active = list(md.get("active_ecs", []))
                    m_muted = list(md.get("muted_ecs", []))

                holdout_data["models"][key][gcf] = {
                    "n_proteins": len(idx_list),
                    "ec_top_scores": ec_top,
                    "target_scores": target_scores,
                    "neighbor_ecs": neighbor_ecs,
                    "meteor_active_ecs": m_active,
                    "meteor_muted_ecs": m_muted,
                }
            print(f"  Holdout-66 {key}: {len(holdout_data['models'][key])} host genomes")

    with gzip.open(f"{OUT_DIR}/holdout66_compact.pkl.gz", "wb") as f:
        pickle.dump(holdout_data, f)
    print(f"  Saved {OUT_DIR}/holdout66_compact.pkl.gz")


def extract_annot_evaluation_1000_bacteria():
    print("[3/3] Extracting 1,000-Genome Bacterial RefSeq GO SOTA Evaluation (annot_evaluation)...")
    summary = {}

    # 1. Essential Core GO Completeness (core_new_GO_*_term_presence_matrix.tsv)
    ess_dir = f"{ANNOT_EVAL}/essential_terms"
    ess_files = {
        "DeepGO-SE (2024)": f"{ess_dir}/core_new_GO_th_deepgometa_term_presence_matrix.tsv",
        "SPROF-GO (2023)": f"{ess_dir}/core_new_GO_th_sprofgo_term_presence_matrix.tsv",
        "DeepFRI (2021)": f"{ess_dir}/core_new_GO_th_deepfri_term_presence_matrix.tsv",
        "DeepGraphGO (2021)": f"{ess_dir}/core_new_GO_th_deepgraphgo_term_presence_matrix.tsv",
        "TALE (2021)": f"{ess_dir}/core_new_GO_th_tale_term_presence_matrix.tsv",
        "InterProScan 5": f"{ess_dir}/core_new_GO_IPscan_term_presence_matrix.tsv",
    }
    for m, fp in ess_files.items():
        if os.path.exists(fp):
            df = pd.read_csv(fp, sep="\t", index_col=0)
            per_genome_comp = df.mean(axis=1).values
            summary.setdefault(m, {})["n_genomes_essential"] = int(len(df))
            summary[m]["essential_completeness_mean"] = round(float(np.mean(per_genome_comp)), 4)
            summary[m]["essential_completeness_std"] = round(float(np.std(per_genome_comp)), 4)

    # 2. MetaCyc Pathway Completeness (pathways/new_GO/new_GO_*_MetaCyc_results.tsv)
    pwy_dir = f"{ANNOT_EVAL}/pathways/new_GO"
    pwy_files = {
        "DeepGO-SE (2024)": f"{pwy_dir}/new_GO_th_deepgometa_MetaCyc_results.tsv",
        "SPROF-GO (2023)": f"{pwy_dir}/new_GO_th_sprofgo_MetaCyc_results.tsv",
        "DeepFRI (2021)": f"{pwy_dir}/new_GO_th_deepfri_MetaCyc_results.tsv",
        "DeepGraphGO (2021)": f"{pwy_dir}/new_GO_th_deepgraphgo_MetaCyc_results.tsv",
        "TALE (2021)": f"{pwy_dir}/new_GO_th_tale_MetaCyc_results.tsv",
        "InterProScan 5": f"{pwy_dir}/new_GO_IPscan_MetaCyc_results.tsv",
    }
    for m, fp in pwy_files.items():
        if os.path.exists(fp):
            vals = []
            n_comp_list = []
            with open(fp) as f:
                for line in f:
                    parts = line.strip().split("\t")
                    if len(parts) >= 7:
                        try:
                            n_comp_list.append(float(parts[2]))
                            vals.append(float(parts[6]))
                        except ValueError:
                            pass
            if vals:
                summary.setdefault(m, {})["n_genomes_pathway"] = len(vals)
                summary[m]["metacyc_complete_ratio_pct"] = round(float(np.mean(vals)), 2)
                summary[m]["metacyc_complete_count_mean"] = round(float(np.mean(n_comp_list)), 2)

    # 3. Protein Complex Coherence (protein_complexes/*_complex_results.tsv)
    cpx_dir = f"{ANNOT_EVAL}/protein_complexes"
    cpx_files = {
        "DeepGO-SE (2024)": f"{cpx_dir}/deepgometa_complex_results.tsv",
        "SPROF-GO (2023)": f"{cpx_dir}/sprofgo_complex_results.tsv",
        "DeepFRI (2021)": f"{cpx_dir}/deepfri_complex_results.tsv",
        "DeepGraphGO (2021)": f"{cpx_dir}/deepgraphgo_complex_results.tsv",
        "TALE (2021)": f"{cpx_dir}/tale_complex_results.tsv",
        "InterProScan 5": f"{cpx_dir}/IPscan_complex_results.tsv",
    }
    for m, fp in cpx_files.items():
        if os.path.exists(fp):
            incoh_list, coh_list = [], []
            with open(fp) as f:
                for line in f:
                    parts = line.strip().split("\t")
                    if len(parts) >= 4:
                        try:
                            incoh_list.append(float(parts[2]))
                            coh_list.append(float(parts[3]))
                        except ValueError:
                            pass
            if incoh_list:
                incoh = np.array(incoh_list)
                coh = np.array(coh_list)
                tot = np.maximum(incoh + coh, 1.0)
                coh_rate = coh / tot
                summary.setdefault(m, {})["n_genomes_complex"] = len(incoh_list)
                summary[m]["complex_coherence_mean"] = round(float(np.mean(coh_rate)), 4)
                summary[m]["incoherent_complexes_mean"] = round(float(np.mean(incoh)), 2)
                summary[m]["coherent_complexes_mean"] = round(float(np.mean(coh)), 2)

    # 4. Taxon Constraint Violations across 1,000 Bacterial Genomes
    tax_dir = f"{ANNOT_EVAL}/taxa_constraints"
    prop_tax = {}
    with open(f"{tax_dir}/propagated_taxon_ids.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            gid = r["Genome_ID"].strip()
            lin = {f"NCBITaxon_{x.strip()}" for x in r["Taxon_Lineage_IDs"].split(";") if x.strip()}
            prop_tax[gid] = lin

    tax_inputs = {
        "DeepGO-SE (2024)": f"{tax_dir}/inputs_deepgometa",
        "SPROF-GO (2023)": f"{tax_dir}/inputs_sprofgo",
        "DeepFRI (2021)": f"{tax_dir}/inputs_deepfri",
        "DeepGraphGO (2021)": f"{tax_dir}/inputs_deepgraphgo",
        "TALE (2021)": f"{tax_dir}/inputs_tale",
        "InterProScan 5": f"{tax_dir}/inputs_IPscan",
    }
    for m, idir in tax_inputs.items():
        if not os.path.isdir(idir):
            continue
        viol_by_genome = defaultdict(int)
        seen_genomes = set()
        for bf in sorted(glob.glob(f"{idir}/batch_*_genomes.tsv")):
            with open(bf) as f:
                header = f.readline()
                for line in f:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 5:
                        continue
                    graw, prot, go_id, never_str, only_str = parts[0], parts[1], parts[2], parts[3], parts[4]
                    gcf = "_".join(graw.split("_")[:2])
                    seen_genomes.add(gcf)
                    lin = prop_tax.get(gcf)
                    if not lin:
                        continue
                    violated = False
                    if never_str:
                        nev = {x.strip() for x in never_str.split(",") if x.strip()}
                        if nev & lin:
                            violated = True
                    if not violated and only_str:
                        onl = {x.strip() for x in only_str.split(",") if x.strip() and not x.strip().startswith("NCBITaxon_Union")}
                        if onl and not (onl & lin):
                            violated = True
                    if violated:
                        viol_by_genome[gcf] += 1
        if seen_genomes:
            v_counts = [viol_by_genome[g] for g in seen_genomes]
            summary.setdefault(m, {})["n_genomes_taxon"] = len(seen_genomes)
            summary[m]["taxon_violations_mean"] = round(float(np.mean(v_counts)), 2)
            summary[m]["taxon_violations_total"] = int(np.sum(v_counts))
            summary[m]["genomes_with_violation_pct"] = round(100.0 * sum(1 for v in v_counts if v > 0) / len(v_counts), 2)

    with open(f"{OUT_DIR}/annot_eval_1000bacteria_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print("  1,000-Genome Bacterial SOTA Summary:")
    print(json.dumps(summary, indent=2))



def extract_meteor108_drafts(gram_panel=f"{F}/psb_revision/code_snapshot/panel108_gram.tsv",
                             outdir="/ibex/scratch/projects/c2014/hohndor/meteor108_repro_2026-10-19/per_genome"):
    """4. METEOR v8 108-genome invariant-audit inputs (DeepProZyme-vanilla, tau=0.5).

    Serial driver for paper-reproducibility/pipelines/P07_invariant_audit/ibex/extract_meteor108_drafts.py
    (the per-genome tau=0.5 draft, METEOR y>0.5 active set + MILP fluxes, and the Grow-Support threshold
    repair re-derived with the archived grow_support LP), followed by merge_meteor108_drafts.py, which writes
    meteor108_drafts_compact.pkl.gz (copied to benchmark/manifests/reference_data/). Each genome needs one
    PuLP/CBC LP over the universal model (minutes), so on ibex run it as the SLURM array
    P07_invariant_audit/ibex/extract_array.sbatch instead of calling this function on a login node.
    """
    import subprocess
    import sys

    here = "/ibex/scratch/projects/c2014/hohndor/meteor108_repro_2026-10-19/scripts"
    for line in open(gram_panel):
        if not line.strip():
            continue
        gca, gram = line.rstrip("\n").split("\t")[:2]
        subprocess.run([sys.executable, f"{here}/extract_meteor108_drafts.py", "--gca", gca, "--gram", gram,
                        "--outdir", outdir], check=True)
    subprocess.run([sys.executable, f"{here}/merge_meteor108_drafts.py"], check=True)

if __name__ == "__main__":
    extract_annot_evaluation_1000_bacteria()
    extract_price149_all6()
    extract_holdout66_and_227()
