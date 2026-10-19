#!/usr/bin/env python3
"""Check whether claims-ID normalisation + exact-duplicate removal changes the ablation grid.

Re-runs the 01_run_ablation_grid.py batch manifest (same configs A1-A9, same jar, same
reference data) on preprocessed claims (run_perturbation.load_base_claims, 'normalised')
into --out-dir, then compares each (genome, config) against the published outputs in
bench9_unified/integrated/ and bench9_unified/cafa_eval/:
  * annotation rows (protein_id, function_id) only-old / only-new,
  * max |delta posterior| over shared rows, rows whose >=0.5 call flips,
  * fmax micro / CAFA (exp, all) via benchmark_pgap_v2.py (--n-bootstrap 50).
Writes ablation_idfix_compare.tsv and claims_preprocessing.tsv to --out-dir.
Does not touch the published ablation outputs.
"""
import argparse
import concurrent.futures
import csv
import json
import os
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import run_perturbation as rp  # noqa: E402

WORK_DIR = f"{rp.ROOT}/proteomes/bench9_unified"
OLD_INT = f"{WORK_DIR}/integrated"
OLD_CAFA = f"{WORK_DIR}/cafa_eval"
CONFIGS = ["A1_noisy_or", "A2_essentiality", "A3_coherence", "A4_gap_filling", "A5_genomic_context",
           "A6_consistency_anatomical", "A7_stage1_all", "A8_stage2_remove", "A9_stage2_relax"]


def manifest_row(g, cfg, claims, out_file):
    """Identical to 01_run_ablation_grid.run_batch."""
    tag, taxon, k = g["tag"], g["taxon"], g["k"]
    gaps = f"{rp.BENCH9}/gapseq/{tag}_gaps.jsonl"
    ops = f"{rp.OPERONS}/{tag}_operons.tsv"
    has_gaps = os.path.exists(gaps) and os.path.getsize(gaps) > 0
    has_ops = os.path.exists(ops) and os.path.getsize(ops) > 0
    tax = ess = o = gp = mode = "-"
    s2, pri = "false", "none"
    if cfg == "A2_essentiality":
        pri, ess = "essentiality", k
    elif cfg == "A3_coherence":
        pri = "coherence"
    elif cfg == "A4_gap_filling":
        pri = "gap_filling"
        gp = gaps if has_gaps else "-"
    elif cfg == "A5_genomic_context":
        pri = "genomic_context"
        o = ops if has_ops else "-"
    elif cfg == "A6_consistency_anatomical":
        pri, tax = "consistency,anatomical", taxon
    elif cfg in ("A7_stage1_all", "A8_stage2_remove", "A9_stage2_relax"):
        pri, tax, ess = rp.FULL_PRIORS, taxon, k
        o = ops if has_ops else "-"
        gp = gaps if has_gaps else "-"
        if cfg != "A7_stage1_all":
            s2, mode = "true", ("remove" if cfg == "A8_stage2_remove" else "relax")
    return "\t".join([f"{tag}_{cfg}", claims, out_file, tax, ess, o, gp, "-", s2, mode, pri]) + "\n"


def read_int(path):
    d = {}
    with open(path) as f:
        h = next(f).rstrip("\n").split("\t")
        ip, ifn, ipp = h.index("protein_id"), h.index("function_id"), h.index("posterior_prob")
        for line in f:
            p = line.rstrip("\n").split("\t")
            d[(p[ip], p[ifn])] = float(p[ipp])
    return d


def fmax(js):
    t = js["by_truth"]
    return (t["exp"]["results"][0]["fmax_overall"], t["exp"]["results"][0]["fmax_cafa_overall"],
            t["all"]["results"][0]["fmax_overall"], t["all"]["results"][0]["fmax_cafa_overall"])


def compare(g, cfg, out_dir):
    tag = g["tag"]
    new_f = f"{out_dir}/{tag}_{cfg}_integrated.tsv"
    old, new = read_int(f"{OLD_INT}/{tag}_{cfg}_integrated.tsv"), read_int(new_f)
    shared = old.keys() & new.keys()
    maxd = max((abs(old[k] - new[k]) for k in shared), default=0.0)
    flips = sum((old[k] >= 0.5) != (new[k] >= 0.5) for k in shared)
    cmd = [rp.PYTHON, rp.BENCH_SCRIPT, "--truth", f"exp:{rp.TRUTH_DIR}/{tag}_truth_exp.tsv",
           "--truth", f"all:{rp.TRUTH_DIR}/{tag}_truth_all.tsv", "--gspa", new_f, "--tag", tag,
           "--n-bootstrap", "50"]
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
    js = json.loads(proc.stdout)
    with open(f"{out_dir}/{tag}_{cfg}_cafa.json", "w") as f:
        json.dump(js, f, indent=2)
    with open(f"{OLD_CAFA}/{tag}_{cfg}_cafa.json") as f:
        fo = fmax(json.load(f))
    fn = fmax(js)
    r = {"genome": tag, "config": cfg, "rows_old": len(old), "rows_new": len(new),
         "only_old": len(old.keys() - new.keys()), "only_new": len(new.keys() - old.keys()),
         "max_abs_dposterior": round(maxd, 6), "call_flips_at_0.5": flips}
    for name, a, b in zip(["fmax_micro_exp", "fmax_cafa_exp", "fmax_micro_all", "fmax_cafa_all"], fo, fn):
        r[f"{name}_old"], r[f"{name}_new"], r[f"d_{name}"] = round(a, 4), round(b, 4), round(b - a, 4)
    return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--jar", default=rp.DEFAULT_JAR)
    ap.add_argument("--jobs", type=int, default=16)
    ap.add_argument("--genomes", default=",".join(g["tag"] for g in rp.GENOMES))
    a = ap.parse_args()
    out = os.path.abspath(a.out_dir)
    os.makedirs(f"{out}/claims", exist_ok=True)
    genomes = [g for g in rp.GENOMES if g["tag"] in a.genomes.split(",")]

    prep = []
    man = f"{out}/ablation_batch_manifest.tsv"
    with open(man, "w") as mf:
        mf.write("genome_id\tclaims_path\tout_path\ttaxon\tessential_profile\toperons_path\tgaps_path"
                 "\tmax_iter\tenforce_stage2\tenforce_mode\tenable_priors\n")
        for g in genomes:
            lines, st = rp.load_base_claims(f"{rp.BENCH9}/with_iprscan/{g['tag']}_claims_merged.jsonl", "normalised")
            cf = f"{out}/claims/{g['tag']}_claims_normalised.jsonl"
            with open(cf, "w") as f:
                f.writelines(lines)
            prep.append(dict(genome=g["tag"], **st))
            for c in CONFIGS:
                mf.write(manifest_row(g, c, cf, f"{out}/{g['tag']}_{c}_integrated.tsv"))
    with open(f"{out}/claims_preprocessing.tsv", "w") as f:
        w = csv.DictWriter(f, delimiter="\t", lineterminator="\n", fieldnames=list(prep[0].keys()))
        w.writeheader()
        w.writerows(prep)

    cmd = [rp.JAVA, "-Xmx48g", "-jar", a.jar, "integrate", "--batch-manifest", man,
           "--batch-threads", str(a.jobs), "--enable-priors", rp.FULL_PRIORS,
           "--go-owl", f"{rp.REF}/go.owl", "--pathways", f"{rp.REF}/kegg_pathways.tsv",
           "--modules", f"{rp.REF}/kegg_modules.tsv", "--ec2go", f"{rp.REF}/ec2go.txt", "--lite"]
    with open(f"{out}/batch_integrate.log", "w") as lf:
        subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT, check=True)

    rows = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=a.jobs) as ex:
        for fut in [ex.submit(compare, g, c, out) for g in genomes for c in CONFIGS]:
            rows.append(fut.result())
    with open(f"{out}/ablation_idfix_compare.tsv", "w") as f:
        w = csv.DictWriter(f, delimiter="\t", lineterminator="\n", fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        print("\t".join(str(v) for v in r.values()), flush=True)


if __name__ == "__main__":
    main()
