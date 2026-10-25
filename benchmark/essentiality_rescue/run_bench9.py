#!/usr/bin/env python3
"""Essentiality rescue, legacy vs targeted, on the nine reference genomes of the ablation (P02).

Runs ``gspa integrate`` per genome, configuration (A2 essentiality only, A7 all Stage 1 priors, A9 Stage 1 +
Stage 2A relaxation) and rescue mode (``legacy`` = v2.1.0 behaviour, ``targeted`` = scored candidate pinned to
the rescue target), with the same inputs and flags as ``01_run_ablation_grid.py``; then evaluates every run:

* rescues (``<out>.rescues.tsv``) against the genome's curated annotation (``truth_all``, propagated through
  go.obo): a rescue is correct at term level when the protein truly carries the rescued term, and at category
  level when it truly carries the uncovered essential category;
* essential-function completeness at posterior 0.50 (fraction of the domain's essential categories with at
  least one protein at or above 0.50, propagated), and how many essential-domain claims sit in [0.45, 0.55);
* micro and CAFA F-max against the experimental and the full truth (benchmark_pgap_v2.py).

Writes ``<out-dir>/runs.tsv`` (one row per genome x config x mode) and ``<out-dir>/summary.tsv``.
"""
import argparse, csv, json, os, subprocess, sys, time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor

JAVA = "/data/hohndor/envs/java21/bin/java"
ROOT = "/data/hohndor/gspa"
REF = f"{ROOT}/reference"
BENCH9 = f"{ROOT}/proteomes/bench9"
OPERONS = f"{ROOT}/proteomes/operons"
TRUTH = f"{ROOT}/proteomes/truth_dual"
PYTHON = "/storage/miniforge3/envs/metagenomics/bin/python3"
BENCH_SCRIPT = f"{ROOT}/benchmark-py/benchmark_pgap_v2.py"

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
ALL_PRIORS = "essentiality,coherence,consistency,anatomical,gap_filling,genomic_context"


def build_command(jar, g, cfg, theta, out_file):
    tag, taxon, k = g["tag"], g["taxon"], g["k"]
    claims = f"{BENCH9}/with_iprscan/{tag}_claims_merged.jsonl"
    gaps = f"{BENCH9}/{tag}_gaps.jsonl"
    ops = f"{OPERONS}/{tag}_operons.tsv"
    gaps_args = ["--gaps", gaps] if os.path.exists(gaps) and os.path.getsize(gaps) > 0 else []
    ops_args = ["--operons", ops] if os.path.exists(ops) and os.path.getsize(ops) > 0 else []
    cmd = [JAVA, "-Xmx8g", "-XX:ActiveProcessorCount=2", "-jar", jar, "integrate", "--claims", claims,
           "--go-owl", f"{REF}/go.owl", "--out", out_file, "--lite", "--theta", theta]
    common = ["--pathways", f"{REF}/kegg_pathways.tsv", "--modules", f"{REF}/kegg_modules.tsv", "--ec2go", f"{REF}/ec2go.txt"]
    if cfg == "A2":
        cmd += ["--enable-priors", "essentiality", "--essential-profile", k] + ops_args
    elif cfg == "A7":
        cmd += ["--enable-priors", ALL_PRIORS, "--essential-profile", k, "--taxon", taxon] + common + ops_args + gaps_args
    elif cfg == "A9":
        cmd += ["--enable-priors", ALL_PRIORS, "--essential-profile", k, "--taxon", taxon] + common + \
               ["--enforce-stage2", "--enforce-consistency-mode", "relax"] + ops_args + gaps_args
    else:
        raise SystemExit(f"unknown config {cfg}")
    return cmd


def load_obo(path):
    """is_a + part_of parents per term; alt_ids mapped to the primary id."""
    parents, alt, cur = defaultdict(set), {}, None
    for line in open(path):
        line = line.rstrip("\n")
        if line == "[Term]":
            cur = None
        elif line.startswith("id: GO:"):
            cur = line[4:]
        elif cur and line.startswith("is_a: "):
            parents[cur].add(line[6:].split(" ")[0])
        elif cur and line.startswith("relationship: part_of "):
            parents[cur].add(line.split(" ")[2])
        elif cur and line.startswith("alt_id: "):
            alt[line[8:]] = cur
    anc_cache = {}

    def ancestors(t):
        t = alt.get(t, t)
        if t in anc_cache:
            return anc_cache[t]
        out, st = {t}, [t]
        while st:
            u = st.pop()
            for p in parents.get(u, ()):
                if p not in out:
                    out.add(p); st.append(p)
        anc_cache[t] = out
        return out
    return ancestors


def load_truth(path):
    by_prot = defaultdict(set)
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            by_prot[r["accession"]].add(r["go_term"])
    return by_prot


def load_essential(path):
    return [l.split("\t")[0].strip() for l in open(path) if l.startswith("GO:")]


def evaluate(tag, k, out_file, ancestors, truth, essentials, n_bootstrap):
    res = {}
    # rescues
    rp = out_file + ".rescues.tsv"
    n = n_term = n_cat = n_unknown = 0
    if os.path.exists(rp):
        with open(rp) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["prior"] != "essentiality":
                    continue
                n += 1
                tt = truth.get(r["protein_id"])
                if not tt:
                    n_unknown += 1; continue
                prop = set()
                for t in tt:
                    prop |= ancestors(t)
                if r["function_id"] in prop: n_term += 1
                if r["category"] in prop: n_cat += 1
    res.update(n_rescues=n, rescues_unknown_truth=n_unknown, rescues_correct_term=n_term, rescues_correct_category=n_cat)
    # completeness and straddle from the integrated table
    ess_closure = {e: set() for e in essentials}
    present = set(); straddle = 0; straddle_ess = 0
    with open(out_file) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for r in rd:
            if r["type"] != "GO":
                continue
            p = float(r["posterior_prob"])
            if 0.45 <= p < 0.55:
                straddle += 1
            if p >= 0.5:
                present |= ancestors(r["function_id"])
    # essential-domain claims in the straddle band: a claim whose term is under an essential category
    with open(out_file) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["type"] != "GO":
                continue
            p = float(r["posterior_prob"])
            if 0.45 <= p < 0.55 and any(e in ancestors(r["function_id"]) for e in essentials):
                straddle_ess += 1
    res.update(essential_present=sum(1 for e in essentials if e in present), essential_total=len(essentials),
               claims_in_045_055=straddle, essential_domain_claims_in_045_055=straddle_ess)
    # F-max
    cmd = [PYTHON, BENCH_SCRIPT, "--truth", f"exp:{TRUTH}/{tag}_truth_exp.tsv", "--truth", f"all:{TRUTH}/{tag}_truth_all.tsv",
           "--gspa", out_file, "--tag", tag, "--n-bootstrap", str(n_bootstrap)]
    pr = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if pr.returncode == 0:
        d = json.loads(pr.stdout)
        for key in ("exp", "all"):
            r0 = d["by_truth"][key]["results"][0]
            res[f"fmax_micro_{key}"] = r0["fmax_overall"]; res[f"fmax_cafa_{key}"] = r0["fmax_cafa_overall"]
    else:
        res["fmax_error"] = pr.stderr[-200:].replace("\n", " ")
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jar", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--modes", default="legacy,targeted")
    ap.add_argument("--configs", default="A2,A7,A9")
    ap.add_argument("--target", type=float, default=0.60)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--n-bootstrap", type=int, default=20)
    ap.add_argument("--essential-dir", required=True, help="gspa-core essential-functions/ (bacteria.tsv, archaea.tsv)")
    ap.add_argument("--skip-existing", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    modes = a.modes.split(","); cfgs = a.configs.split(",")
    thetas = {}
    for m in modes:
        tp = os.path.join(a.out_dir, f"theta_{m}.json")
        json.dump({"ess_rescue_mode": m, "ess_rescue_target": a.target}, open(tp, "w"))
        thetas[m] = tp
    jobs = [(g, c, m) for g in GENOMES for c in cfgs for m in modes]

    def run(job):
        g, c, m = job
        out_file = os.path.join(a.out_dir, f"{g['tag']}_{c}_{m}_integrated.tsv")
        if a.skip_existing and os.path.exists(out_file) and os.path.getsize(out_file) > 1000:
            return (job, out_file, "skipped", 0.0)
        cmd = build_command(a.jar, g, c, thetas[m], out_file)
        t0 = time.time()
        log = open(out_file + ".log", "w")
        pr = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, text=True)
        log.close()
        return (job, out_file, "ok" if pr.returncode == 0 else f"rc={pr.returncode}", time.time() - t0)
    print(f"{len(jobs)} integrate runs, {a.jobs} in parallel", flush=True)
    with ThreadPoolExecutor(a.jobs) as ex:
        for (g, c, m), out_file, st, dt in ex.map(run, jobs):
            print(f"  {g['tag']} {c} {m}: {st} {dt:.0f}s", flush=True)

    ancestors = load_obo(f"{REF}/go.obo")
    ess = {k: load_essential(os.path.join(a.essential_dir, f"{k}.tsv")) for k in ("bacteria", "archaea")}
    rows = []
    for g in GENOMES:
        truth = load_truth(f"{TRUTH}/{g['tag']}_truth_all.tsv")
        for c in cfgs:
            for m in modes:
                out_file = os.path.join(a.out_dir, f"{g['tag']}_{c}_{m}_integrated.tsv")
                if not os.path.exists(out_file):
                    continue
                r = dict(genome=g["tag"], config=c, mode=m)
                r.update(evaluate(g["tag"], g["k"], out_file, ancestors, truth, ess[g["k"]], a.n_bootstrap))
                rows.append(r); print(f"  eval {g['tag']} {c} {m}: {r}", flush=True)
    cols = sorted({k for r in rows for k in r}, key=lambda k: (k not in ("genome", "config", "mode"), k))
    with open(os.path.join(a.out_dir, "runs.tsv"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(k, "")) for k in cols) + "\n")
    # summary per config x mode
    with open(os.path.join(a.out_dir, "summary.tsv"), "w") as fh:
        keys = ["n_rescues", "rescues_correct_term", "rescues_correct_category", "rescues_unknown_truth", "essential_present",
                "essential_domain_claims_in_045_055", "fmax_micro_exp", "fmax_cafa_exp", "fmax_micro_all", "fmax_cafa_all"]
        fh.write("config\tmode\tn_genomes\t" + "\t".join(f"sum_{k}" if k.startswith(("n_", "rescues", "essential")) else f"mean_{k}" for k in keys) +
                 "\tprecision_term\tprecision_category\n")
        for c in cfgs:
            for m in modes:
                sel = [r for r in rows if r["config"] == c and r["mode"] == m]
                if not sel:
                    continue
                vals = []
                for k in keys:
                    xs = [r[k] for r in sel if k in r]
                    vals.append(f"{sum(xs)}" if k.startswith(("n_", "rescues", "essential")) else (f"{sum(xs) / len(xs):.4f}" if xs else ""))
                judged = sum(r["n_rescues"] - r["rescues_unknown_truth"] for r in sel)
                pt = sum(r["rescues_correct_term"] for r in sel) / judged if judged else float("nan")
                pc = sum(r["rescues_correct_category"] for r in sel) / judged if judged else float("nan")
                fh.write(f"{c}\t{m}\t{len(sel)}\t" + "\t".join(vals) + f"\t{pt:.3f}\t{pc:.3f}\n")
    print("done", flush=True)


if __name__ == "__main__":
    main()
