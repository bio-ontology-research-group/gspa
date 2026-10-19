#!/usr/bin/env python3
"""
P02 -- nine-genome cross-clade taxon / anatomical perturbation experiment.

Derived from benchmark/unified/pillar1_2/04_run_perturbation.py (commit 3310b97,
design v2) and its predecessor at commit 326c112 (design v1). One script, two
injection designs selected with --design:

  v1 (326c112): every injected (protein, GO term) pair is ONE claim from source
      "deepgo" with raw_score ~ U(0.55, 0.85).
  v2 (3310b97): primary source alternates diamond / deepgo (even / odd term
      index), raw_score ~ U(0.62, 0.90); for 60% of pairs (term index mod 5 in
      {0,1,2}) a second corroborating claim is added (pfam if the primary is
      diamond, else diamond) with score max(0.52, raw_score - 0.12).

Both designs inject the same 100 GO terms (50 anatomical, 50 eukaryote-only
taxon-constrained; term-selection seed 42) into the same 100 proteins per genome
(protein seed = md5(genome tag) mod 2^32). Because random.uniform(a, b) is
a + (b - a) * random(), the v1 and v2 scores come from the same random draws,
only rescaled.

Configs: A1_noisy_or, A6_consistency_anatomical, A7_stage1_all,
A8_stage2_remove, A9_stage2_relax.

Claims preprocessing (--claims-ids normalised, the default since 2026-10-19):
the merged claims of hpylori, mgenitalium and mjannaschii carry raw FASTA IDs
("sp|P21762|AHPC_HELPY") on every pfam claim (11,325 / 5,739 / 13,773 claims).
`gspa integrate` normalises such IDs internally, but the injection step sampled
proteins from the raw ID set, so 47 / 42 / 55 of the 100 injected pairs landed
on phantom IDs and were never matched by the survivor count. Before sampling,
every claims file (all nine genomes) is therefore
  1. ID-normalised: protein_id matching ^(sp|tr)|ACC|NAME -> ACC;
  2. de-duplicated: exact duplicate records (identical JSON after step 1, keys
     sorted) are dropped, first occurrence kept (6-96 per genome).
The preprocessed base claims are what is written (plus injections) to
{tag}_claims_perturbed.jsonl; counts go to claims_preprocessing.tsv.
--claims-ids raw reproduces the earlier runs (no preprocessing).

Integration runner: the shared-memory `integrate --batch-manifest` mode when the
jar supports it (gspa-unified.jar >= 2.1.0), otherwise one JVM per
(genome, config) exactly as the 326c112 script did (needed for the old
gspa-unified-bench.jar, which predates --batch-manifest). --runner forces one.

Outputs (under --out-dir):
  {tag}_injected_pairs.tsv, {tag}_claims_perturbed.jsonl,
  {tag}_{cfg}_perturbed_integrated.tsv, {tag}_{cfg}_cafa_eval.json,
  {tag}_A8/A9_perturbed_actions.tsv (per-job runner only), logs/,
  perturbation_summary.json, perturbation_summary.tsv, run_info.json
"""

import argparse
import concurrent.futures
import csv
import hashlib
import json
import os
import random
import subprocess
import sys
import time
import re
import zipfile

PYTHON = "/storage/miniforge3/envs/metagenomics/bin/python3"
JAVA = "/data/hohndor/envs/java21/bin/java"
DEFAULT_JAR = "/data/hohndor/gspa/bin/gspa-unified.jar"
BENCH_SCRIPT = "/data/hohndor/gspa/benchmark-py/benchmark_pgap_v2.py"
ROOT = "/data/hohndor/gspa"
REF = f"{ROOT}/reference"
BENCH9 = f"{ROOT}/proteomes/bench9"
OPERONS = f"{ROOT}/proteomes/operons"
TRUTH_DIR = f"{ROOT}/proteomes/truth_dual"
# Constraint TSVs used for term selection (read-only; same files the originals read)
CONSTRAINT_DIR = f"{ROOT}/proteomes/bench9_unified"

GENOMES = [
    {"tag": "ecoli", "taxon": "NCBITaxon_562", "k": "bacteria"},
    {"tag": "hpylori", "taxon": "NCBITaxon_210", "k": "bacteria"},
    {"tag": "mgenitalium", "taxon": "NCBITaxon_2097", "k": "bacteria"},
    {"tag": "mjannaschii", "taxon": "NCBITaxon_2190", "k": "archaea"},
    {"tag": "ecolo157", "taxon": "NCBITaxon_562", "k": "bacteria"},
    {"tag": "bsubtilis", "taxon": "NCBITaxon_1423", "k": "bacteria"},
    {"tag": "mtb", "taxon": "NCBITaxon_1773", "k": "bacteria"},
    {"tag": "synechocystis", "taxon": "NCBITaxon_1148", "k": "bacteria"},
    {"tag": "paeruginosa", "taxon": "NCBITaxon_287", "k": "bacteria"},
]

PERTURB_CONFIGS = [
    "A1_noisy_or",
    "A6_consistency_anatomical",
    "A7_stage1_all",
    "A8_stage2_remove",
    "A9_stage2_relax",
]

FULL_PRIORS = "essentiality,coherence,consistency,anatomical,gap_filling,genomic_context"

RAW_FASTA_ID = re.compile(r"^(?:sp|tr)\|([^|]+)\|")

SUMMARY_COLUMNS = [
    "design", "jar", "claims_ids", "genome", "config", "injected_pairs", "injected_claims",
    "surviving_violations", "surviving_anatomical", "surviving_taxon",
    "suppression_rate_pct", "fmax_micro_all", "fmax_micro_exp",
    "fmax_cafa_all", "fmax_cafa_exp",
]


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def jar_supports_batch(jar):
    with zipfile.ZipFile(jar) as z:
        data = z.read("gspa/cli/IntegrateCommand.class")
    return b"--batch-manifest" in data


def load_aspect_map():
    aspect_map = {}
    path = f"{REF}/go_aspect_map.tsv"
    if os.path.exists(path):
        with open(path) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    aspect_map[parts[0]] = parts[1]
    return aspect_map


def pick_perturbation_terms(constraint_dir):
    """Pick 50 anatomical and 50 taxon-constrained terms (identical to originals)."""
    aspect_map = load_aspect_map()
    mandatory_anat = [
        "GO:0007507",  # heart development
        "GO:0009908",  # flower development
        "GO:0001755",  # neural crest cell migration
        "GO:0030282",  # bone mineralization
        "GO:0048437",  # floral organ development
    ]
    anat_candidates = []
    with open(f"{constraint_dir}/go-anatomical-constraints.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            go_id = r["GO_ID"]
            anat_id = r["Anatomy_ID"]
            if anat_id.startswith(("UBERON:", "PO:", "CL:")):
                if go_id not in mandatory_anat:
                    anat_candidates.append(go_id)

    rng = random.Random(42)
    unique_anat = sorted(set(anat_candidates))
    rng.shuffle(unique_anat)
    selected_anat = mandatory_anat + unique_anat[:45]

    euk_taxa = {
        "NCBITaxon_2759", "NCBITaxon_33208", "NCBITaxon_33090",
        "NCBITaxon_4751", "NCBITaxon_7711", "NCBITaxon_40674",
        "NCBITaxon_89593", "NCBITaxon_9606",
    }
    taxon_candidates = []
    with open(f"{constraint_dir}/go-taxon-constraints.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["Constraint_Type"] == "only_in_taxon" and r["Taxon_ID"] in euk_taxa:
                if r["GO_ID"] not in selected_anat:
                    taxon_candidates.append(r["GO_ID"])
    unique_tax = sorted(set(taxon_candidates))
    rng.shuffle(unique_tax)
    selected_tax = unique_tax[:50]

    terms = [{"go_id": t, "category": "anatomical", "aspect": aspect_map.get(t, "BP")} for t in selected_anat]
    terms += [{"go_id": t, "category": "taxon_constraint", "aspect": aspect_map.get(t, "BP")} for t in selected_tax]
    return terms


def make_injected_claims(design, terms, chosen_proteins, rng):
    injected_claims, injected_records = [], []
    for idx_t, (term_info, pid) in enumerate(zip(terms, chosen_proteins)):
        meta = {"injected_perturbation": True, "category": term_info["category"]}
        if design == "v1":
            raw_score = round(rng.uniform(0.55, 0.85), 4)
            src_primary = "deepgo"
        else:
            raw_score = round(rng.uniform(0.62, 0.90), 4)
            src_primary = "diamond" if idx_t % 2 == 0 else "deepgo"
        injected_claims.append({
            "protein_id": pid, "function_type": "GO", "function_id": term_info["go_id"],
            "go_aspect": term_info["aspect"], "source": src_primary,
            "raw_score": raw_score, "metadata": meta,
        })
        if design == "v2" and idx_t % 5 in (0, 1, 2):
            injected_claims.append({
                "protein_id": pid, "function_type": "GO", "function_id": term_info["go_id"],
                "go_aspect": term_info["aspect"],
                "source": "pfam" if src_primary == "diamond" else "diamond",
                "raw_score": round(max(0.52, raw_score - 0.12), 4),
                "metadata": dict(meta),
            })
        injected_records.append({
            "protein_id": pid, "go_id": term_info["go_id"], "aspect": term_info["aspect"],
            "category": term_info["category"], "raw_score": raw_score,
        })
    return injected_claims, injected_records


def load_base_claims(claims_file, claims_ids):
    """Return (lines, stats). claims_ids='raw': file as is; 'normalised': ID
    normalisation + exact-duplicate removal (see module docstring)."""
    lines, seen = [], set()
    st = {"lines_in": 0, "ids_normalised": 0, "duplicates_removed": 0, "lines_out": 0}
    with open(claims_file) as f:
        for line in f:
            st["lines_in"] += 1
            if claims_ids == "raw":
                lines.append(line if line.endswith("\n") else line + "\n")
                continue
            c = json.loads(line)
            m = RAW_FASTA_ID.match(c.get("protein_id", ""))
            if m:
                c["protein_id"] = m.group(1)
                st["ids_normalised"] += 1
                line = json.dumps(c) + "\n"
            key = json.dumps(c, sort_keys=True)
            if key in seen:
                st["duplicates_removed"] += 1
                continue
            seen.add(key)
            lines.append(line if line.endswith("\n") else line + "\n")
    st["lines_out"] = len(lines)
    return lines, st


def inject_claims_for_genome(genome, terms, design, out_dir, claims_ids="normalised"):
    tag = genome["tag"]
    claims_file = f"{BENCH9}/with_iprscan/{tag}_claims_merged.jsonl"
    out_claims_file = f"{out_dir}/{tag}_claims_perturbed.jsonl"
    injected_info_file = f"{out_dir}/{tag}_injected_pairs.tsv"

    base_lines, stats = load_base_claims(claims_file, claims_ids)
    proteins = set()
    for line in base_lines:
        try:
            proteins.add(json.loads(line)["protein_id"])
        except Exception:
            pass
    unique_proteins = sorted(proteins)
    stats["proteins"] = len(unique_proteins)

    seed = int(hashlib.md5(tag.encode("utf-8")).hexdigest(), 16) % (2**32)
    rng = random.Random(seed)
    if len(unique_proteins) >= len(terms):
        chosen = rng.sample(unique_proteins, len(terms))
    else:
        chosen = [rng.choice(unique_proteins) for _ in range(len(terms))]

    injected_claims, injected_records = make_injected_claims(design, terms, chosen, rng)

    with open(injected_info_file, "w") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=["protein_id", "go_id", "aspect", "category", "raw_score"])
        w.writeheader()
        w.writerows(injected_records)
    with open(out_claims_file, "w") as out_f:
        out_f.writelines(base_lines)
        for c in injected_claims:
            out_f.write(json.dumps(c) + "\n")
    return out_claims_file, injected_records, len(injected_claims), stats


def per_job_command(jar, genome, cfg, claims_file, out_dir, xmx):
    """One JVM per (genome, config); identical flags to the 326c112 / 3310b97 scripts."""
    tag, taxon, k = genome["tag"], genome["taxon"], genome["k"]
    out_file = f"{out_dir}/{tag}_{cfg}_perturbed_integrated.tsv"
    gaps_file = f"{BENCH9}/gapseq/{tag}_gaps.jsonl"
    gaps_args = ["--gaps", gaps_file] if os.path.exists(gaps_file) and os.path.getsize(gaps_file) > 0 else []
    operons_file = f"{OPERONS}/{tag}_operons.tsv"
    operons_args = ["--operons", operons_file] if os.path.exists(operons_file) and os.path.getsize(operons_file) > 0 else []
    cmd = [JAVA, f"-Xmx{xmx}", "-jar", jar, "integrate",
           "--claims", claims_file, "--go-owl", f"{REF}/go.owl", "--out", out_file, "--lite"]
    full = ["--enable-priors", FULL_PRIORS, "--essential-profile", k, "--taxon", taxon,
            "--pathways", f"{REF}/kegg_pathways.tsv", "--modules", f"{REF}/kegg_modules.tsv",
            "--ec2go", f"{REF}/ec2go.txt"]
    if cfg == "A1_noisy_or":
        cmd += ["--enable-priors", "none"]
    elif cfg == "A6_consistency_anatomical":
        cmd += ["--enable-priors", "consistency,anatomical", "--taxon", taxon]
    elif cfg == "A7_stage1_all":
        cmd += full + operons_args + gaps_args
    elif cfg in ("A8_stage2_remove", "A9_stage2_relax"):
        short = "A8" if cfg == "A8_stage2_remove" else "A9"
        mode = "remove" if short == "A8" else "relax"
        cmd += full + ["--enforce-stage2", "--enforce-consistency-mode", mode,
                       "--enforcement-actions-out", f"{out_dir}/{tag}_{short}_perturbed_actions.tsv"]
        cmd += operons_args + gaps_args
    return cmd, out_file


def run_per_job(jar, tasks, claims_by_genome, out_dir, jobs, xmx):
    log_dir = f"{out_dir}/logs"

    def one(g, cfg):
        cmd, out_file = per_job_command(jar, g, cfg, claims_by_genome[g["tag"]], out_dir, xmx)
        t0 = time.time()
        with open(f"{log_dir}/{g['tag']}_{cfg}_perturbed.log", "w") as lf:
            proc = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT)
        return g["tag"], cfg, proc.returncode, time.time() - t0

    failed = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as ex:
        for fut in concurrent.futures.as_completed([ex.submit(one, g, c) for g, c in tasks]):
            tag, cfg, rc, dt = fut.result()
            print(f"[{'DONE' if rc == 0 else 'FAIL'}-INT] {tag} {cfg} rc={rc} {dt:.1f}s", flush=True)
            if rc != 0:
                failed.append((tag, cfg))
    if failed:
        raise RuntimeError(f"per-job integration failed for {failed}")


def run_batch(jar, tasks, claims_by_genome, out_dir, jobs, xmx):
    """Shared-memory batch runner (identical manifest to 3310b97 04_run_perturbation.py)."""
    manifest_path = f"{out_dir}/perturb_batch_manifest.tsv"
    with open(manifest_path, "w") as mf:
        mf.write("genome_id\tclaims_path\tout_path\ttaxon\tessential_profile\toperons_path\tgaps_path"
                 "\tmax_iter\tenforce_stage2\tenforce_mode\tenable_priors\n")
        for g, cfg in tasks:
            tag, taxon, k = g["tag"], g["taxon"], g["k"]
            out_file = f"{out_dir}/{tag}_{cfg}_perturbed_integrated.tsv"
            gaps_file = f"{BENCH9}/gapseq/{tag}_gaps.jsonl"
            has_gaps = os.path.exists(gaps_file) and os.path.getsize(gaps_file) > 0
            operons_file = f"{OPERONS}/{tag}_operons.tsv"
            has_ops = os.path.exists(operons_file) and os.path.getsize(operons_file) > 0
            row_tax = row_ess = row_ops = row_gaps = row_mode = "-"
            row_s2, row_priors = "false", "none"
            if cfg == "A6_consistency_anatomical":
                row_priors, row_tax = "consistency,anatomical", taxon
            elif cfg in ("A7_stage1_all", "A8_stage2_remove", "A9_stage2_relax"):
                row_priors, row_tax, row_ess = FULL_PRIORS, taxon, k
                row_ops = operons_file if has_ops else "-"
                row_gaps = gaps_file if has_gaps else "-"
                if cfg != "A7_stage1_all":
                    row_s2 = "true"
                    row_mode = "remove" if cfg == "A8_stage2_remove" else "relax"
            mf.write("\t".join([f"{tag}_{cfg}_perturb", claims_by_genome[tag], out_file, row_tax, row_ess,
                                row_ops, row_gaps, "-", row_s2, row_mode, row_priors]) + "\n")
    cmd = [JAVA, f"-Xmx{xmx}", "-jar", jar, "integrate",
           "--batch-manifest", manifest_path, "--batch-threads", str(jobs),
           "--enable-priors", FULL_PRIORS, "--go-owl", f"{REF}/go.owl",
           "--pathways", f"{REF}/kegg_pathways.tsv", "--modules", f"{REF}/kegg_modules.tsv",
           "--ec2go", f"{REF}/ec2go.txt", "--lite"]
    print(f"Running {len(tasks)} integration jobs via --batch-manifest", flush=True)
    with open(f"{out_dir}/logs/batch_integrate.log", "w") as lf:
        proc = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT)
    if proc.returncode != 0:
        raise RuntimeError(f"batch integration failed rc={proc.returncode}")


def evaluate(genome, cfg, injected_records, out_dir):
    tag = genome["tag"]
    int_file = f"{out_dir}/{tag}_{cfg}_perturbed_integrated.tsv"
    eval_json = f"{out_dir}/{tag}_{cfg}_cafa_eval.json"
    injected_set = {(r["protein_id"], r["go_id"]) for r in injected_records}
    anat_set = {r["go_id"] for r in injected_records if r["category"] == "anatomical"}
    tax_set = {r["go_id"] for r in injected_records if r["category"] == "taxon_constraint"}
    surviving, s_anat, s_tax = set(), 0, 0
    with open(int_file) as f:
        headers = next(f, "").rstrip("\n").split("\t")
        i_pid, i_fn, i_p = headers.index("protein_id"), headers.index("function_id"), headers.index("posterior_prob")
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= i_p:
                continue
            try:
                prob = float(parts[i_p])
            except ValueError:
                continue
            key = (parts[i_pid], parts[i_fn])
            if prob >= 0.5 and key in injected_set:
                surviving.add(key)
                s_anat += parts[i_fn] in anat_set
                s_tax += parts[i_fn] in tax_set

    res = {"genome": tag, "config": cfg, "injected_pairs": len(injected_records),
           "surviving_violations": len(surviving), "surviving_anatomical": s_anat,
           "surviving_taxon": s_tax,
           "suppression_rate_pct": round((1.0 - len(surviving) / len(injected_records)) * 100.0, 1),
           "fmax_micro_exp": None, "fmax_cafa_exp": None, "fmax_micro_all": None, "fmax_cafa_all": None}
    cmd = [PYTHON, BENCH_SCRIPT,
           "--truth", f"exp:{TRUTH_DIR}/{tag}_truth_exp.tsv",
           "--truth", f"all:{TRUTH_DIR}/{tag}_truth_all.tsv",
           "--gspa", int_file, "--tag", tag, "--n-bootstrap", "50"]
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if proc.returncode == 0:
        d = json.loads(proc.stdout)
        with open(eval_json, "w") as f:
            json.dump(d, f, indent=2)
        res["fmax_micro_exp"] = d["by_truth"]["exp"]["results"][0]["fmax_overall"]
        res["fmax_cafa_exp"] = d["by_truth"]["exp"]["results"][0]["fmax_cafa_overall"]
        res["fmax_micro_all"] = d["by_truth"]["all"]["results"][0]["fmax_overall"]
        res["fmax_cafa_all"] = d["by_truth"]["all"]["results"][0]["fmax_cafa_overall"]
    else:
        print(f"[WARN] cafa eval failed {tag} {cfg}: {proc.stderr[-500:]}", flush=True)
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--design", choices=["v1", "v2"], required=True)
    ap.add_argument("--jar", default=DEFAULT_JAR)
    ap.add_argument("--jar-label", default=None, help="label written to the 'jar' column (default: jar basename)")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--claims-ids", choices=["normalised", "raw"], default="normalised",
                    help="normalised (default): FASTA-ID normalisation + exact-duplicate removal; raw: as before 2026-10-19")
    ap.add_argument("--runner", choices=["auto", "batch", "perjob"], default="auto")
    ap.add_argument("--jobs", "-j", type=int, default=16)
    ap.add_argument("--xmx", default=None, help="JVM heap (default 48g batch, 8g per job)")
    ap.add_argument("--constraint-dir", default=CONSTRAINT_DIR)
    ap.add_argument("--configs", default=",".join(PERTURB_CONFIGS))
    args = ap.parse_args()

    out_dir = os.path.abspath(args.out_dir)
    os.makedirs(f"{out_dir}/logs", exist_ok=True)
    configs = args.configs.split(",")
    runner = args.runner
    if runner == "auto":
        runner = "batch" if jar_supports_batch(args.jar) else "perjob"
    xmx = args.xmx or ("48g" if runner == "batch" else "8g")
    jar_label = args.jar_label or os.path.basename(args.jar)

    t_start = time.time()
    run_info = {"design": args.design, "claims_ids": args.claims_ids, "jar": args.jar, "jar_label": jar_label,
                "jar_sha256": sha256_file(args.jar), "runner": runner, "xmx": xmx,
                "jobs": args.jobs, "configs": configs, "host": os.uname().nodename,
                "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
                "go_owl_sha256": sha256_file(f"{REF}/go.owl"),
                "started": time.strftime("%Y-%m-%dT%H:%M:%S")}
    print(json.dumps(run_info, indent=2), flush=True)

    terms = pick_perturbation_terms(args.constraint_dir)
    assert len(terms) == 100, len(terms)
    injected_by_genome, claims_by_genome, n_claims, prep = {}, {}, {}, []
    for g in GENOMES:
        cf, recs, nc, st = inject_claims_for_genome(g, terms, args.design, out_dir, args.claims_ids)
        injected_by_genome[g["tag"]], claims_by_genome[g["tag"]], n_claims[g["tag"]] = recs, cf, nc
        prep.append(dict(genome=g["tag"], **st))
        print(f"Injected {len(recs)} pairs / {nc} claims into {g['tag']} (preprocessing: {st})", flush=True)
    with open(f"{out_dir}/claims_preprocessing.tsv", "w") as f:
        w = csv.DictWriter(f, delimiter="\t", lineterminator="\n", fieldnames=[
            "genome", "lines_in", "ids_normalised", "duplicates_removed", "lines_out", "proteins"])
        w.writeheader()
        w.writerows(prep)

    tasks = [(g, c) for g in GENOMES for c in configs]
    t_int = time.time()
    if runner == "batch":
        run_batch(args.jar, tasks, claims_by_genome, out_dir, args.jobs, xmx)
    else:
        run_per_job(args.jar, tasks, claims_by_genome, out_dir, args.jobs, xmx)
    run_info["integration_seconds"] = round(time.time() - t_int, 1)

    t_ev = time.time()
    results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as ex:
        futs = [ex.submit(evaluate, g, c, injected_by_genome[g["tag"]], out_dir) for g, c in tasks]
        for fut in concurrent.futures.as_completed(futs):
            r = fut.result()
            r["design"], r["jar"], r["claims_ids"] = args.design, jar_label, args.claims_ids
            r["injected_claims"] = n_claims[r["genome"]]
            results.append(r)
            print(f"[EVAL] {r['genome']} {r['config']}: surviving={r['surviving_violations']}/100 "
                  f"fmax_micro_exp={r['fmax_micro_exp']}", flush=True)
    run_info["eval_seconds"] = round(time.time() - t_ev, 1)

    order = {g["tag"]: i for i, g in enumerate(GENOMES)}
    results.sort(key=lambda r: (order[r["genome"]], configs.index(r["config"])))
    with open(f"{out_dir}/perturbation_summary.json", "w") as f:
        json.dump(results, f, indent=2)
    with open(f"{out_dir}/perturbation_summary.tsv", "w") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=SUMMARY_COLUMNS, extrasaction="ignore")
        w.writeheader()
        for r in results:
            w.writerow({k: ("" if r.get(k) is None else r.get(k)) for k in SUMMARY_COLUMNS})

    run_info["total_seconds"] = round(time.time() - t_start, 1)
    run_info["finished"] = time.strftime("%Y-%m-%dT%H:%M:%S")
    with open(f"{out_dir}/run_info.json", "w") as f:
        json.dump(run_info, f, indent=2)

    print("\nPer-config survivors (9-genome total /900):")
    for c in configs:
        rows = [r for r in results if r["config"] == c]
        print(f"  {c}: {sum(r['surviving_violations'] for r in rows)}  "
              f"anat={sum(r['surviving_anatomical'] for r in rows)} tax={sum(r['surviving_taxon'] for r in rows)}")
    print(json.dumps(run_info, indent=2), flush=True)


if __name__ == "__main__":
    sys.exit(main())
