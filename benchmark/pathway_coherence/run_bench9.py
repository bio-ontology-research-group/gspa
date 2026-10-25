#!/usr/bin/env python3
"""Pathway coherence and refinement on the nine reference genomes: legacy vs module-based.

Variants per genome (all with the six priors and Stage 2A relaxation, as the A9 ablation configuration):
  legacy   : v2.1.0 pathway handling, no module definitions (the current behaviour)
  assess   : module definitions loaded, pathway refinement in legacy mode (same annotations as `legacy`
             plus the module assessment, so the metric change is separated from the refinement change)
  module   : module definitions loaded, module-based refinement on

Evaluation per run:
  * micro and CAFA F-max against experimental and full curated truth (benchmark_pgap_v2.py);
  * coherence rescues (`<out>.rescues.tsv`, prior=coherence) judged against the curated annotation of the
    rescued protein: correct at term level (the rescued GO term or a descendant is curated) and at step level
    (any enzyme of the step is curated for that protein);
  * module-level agreement with KEGG: for the genome's KEGG organism, a module KEGG lists as complete is the
    reference; we count modules predicted complete / activated-partial / not activated, precision and recall
    of "complete", before and after refinement (`<out>.modules.json`);
  * the old triggered-pathway score for reference.
Writes runs.tsv, modules.tsv (one row per genome x module x variant) and summary.tsv.
"""
import argparse, csv, json, os, re, subprocess, time
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
    {"tag": "ecoli", "taxon": "NCBITaxon_562", "k": "bacteria", "org": "eco"},
    {"tag": "hpylori", "taxon": "NCBITaxon_210", "k": "bacteria", "org": "hpy"},
    {"tag": "mgenitalium", "taxon": "NCBITaxon_2097", "k": "bacteria", "org": "mge"},
    {"tag": "mjannaschii", "taxon": "NCBITaxon_2190", "k": "archaea", "org": "mja"},
    {"tag": "ecolo157", "taxon": "NCBITaxon_562", "k": "bacteria", "org": "ece"},
    {"tag": "bsubtilis", "taxon": "NCBITaxon_1423", "k": "bacteria", "org": "bsu"},
    {"tag": "mtb", "taxon": "NCBITaxon_1773", "k": "bacteria", "org": "mtu"},
    {"tag": "synechocystis", "taxon": "NCBITaxon_1148", "k": "bacteria", "org": "syn"},
    {"tag": "paeruginosa", "taxon": "NCBITaxon_287", "k": "bacteria", "org": "pae"},
]
PRIORS = "essentiality,coherence,consistency,anatomical,gap_filling,genomic_context"
VARIANTS = {"legacy": None, "assess": {"coh_pathway_mode": "legacy"}, "module": {"coh_pathway_mode": "module"}}


def command(jar, g, variant, theta, defs, out_file):
    tag = g["tag"]
    gaps = f"{BENCH9}/{tag}_gaps.jsonl"; ops = f"{OPERONS}/{tag}_operons.tsv"
    cmd = [JAVA, "-Xmx8g", "-XX:ActiveProcessorCount=2", "-jar", jar, "integrate", "--claims", f"{BENCH9}/with_iprscan/{tag}_claims_merged.jsonl",
           "--go-owl", f"{REF}/go.owl", "--out", out_file, "--lite", "--enable-priors", PRIORS, "--essential-profile", g["k"], "--taxon", g["taxon"],
           "--pathways", f"{REF}/kegg_pathways.tsv", "--modules", f"{REF}/kegg_modules.tsv", "--ec2go", f"{REF}/ec2go.txt",
           "--enforce-stage2", "--enforce-consistency-mode", "relax"]
    if os.path.exists(ops) and os.path.getsize(ops) > 0: cmd += ["--operons", ops]
    if os.path.exists(gaps) and os.path.getsize(gaps) > 0: cmd += ["--gaps", gaps]
    if theta: cmd += ["--theta", theta]
    if variant != "legacy": cmd += ["--module-defs", defs]
    return cmd


def load_obo(path):
    parents, alt, cur = defaultdict(set), {}, None
    for line in open(path):
        line = line.rstrip("\n")
        if line == "[Term]": cur = None
        elif line.startswith("id: GO:"): cur = line[4:]
        elif cur and line.startswith("is_a: "): parents[cur].add(line[6:].split(" ")[0])
        elif cur and line.startswith("relationship: part_of "): parents[cur].add(line.split(" ")[2])
        elif cur and line.startswith("alt_id: "): alt[line[8:]] = cur
    cache = {}
    def anc(t):
        t = alt.get(t, t)
        if t in cache: return cache[t]
        out, st = {t}, [t]
        while st:
            u = st.pop()
            for p in parents.get(u, ()):
                if p not in out: out.add(p); st.append(p)
        cache[t] = out; return out
    return anc


def load_truth(path):
    d = defaultdict(set)
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"): d[r["accession"]].add(r["go_term"])
    return d


def fmax(tag, out_file, n_bootstrap):
    cmd = [PYTHON, BENCH_SCRIPT, "--truth", f"exp:{TRUTH}/{tag}_truth_exp.tsv", "--truth", f"all:{TRUTH}/{tag}_truth_all.tsv", "--gspa", out_file, "--tag", tag, "--n-bootstrap", str(n_bootstrap)]
    pr = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if pr.returncode != 0: return {"fmax_error": pr.stderr[-200:].replace("\n", " ")}
    d = json.loads(pr.stdout); out = {}
    for key in ("exp", "all"):
        r0 = d["by_truth"][key]["results"][0]; out[f"fmax_micro_{key}"] = r0["fmax_overall"]; out[f"fmax_cafa_{key}"] = r0["fmax_cafa_overall"]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jar", required=True); ap.add_argument("--defs", required=True, help="kegg_module_defs.json")
    ap.add_argument("--out-dir", required=True); ap.add_argument("--jobs", type=int, default=3); ap.add_argument("--n-bootstrap", type=int, default=20)
    ap.add_argument("--variants", default="legacy,assess,module"); ap.add_argument("--skip-existing", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    variants = a.variants.split(",")
    thetas = {}
    for v in variants:
        if VARIANTS[v] is not None:
            tp = os.path.join(a.out_dir, f"theta_{v}.json"); json.dump(VARIANTS[v], open(tp, "w")); thetas[v] = tp
        else: thetas[v] = None
    jobs = [(g, v) for g in GENOMES for v in variants]
    def run(job):
        g, v = job; out_file = os.path.join(a.out_dir, f"{g['tag']}_{v}_integrated.tsv")
        if a.skip_existing and os.path.exists(out_file) and os.path.getsize(out_file) > 1000: return (job, "skipped", 0.0)
        t0 = time.time()
        with open(out_file + ".log", "w") as log:
            pr = subprocess.run(command(a.jar, g, v, thetas[v], a.defs, out_file), stdout=log, stderr=subprocess.STDOUT, text=True)
        return (job, "ok" if pr.returncode == 0 else f"rc={pr.returncode}", time.time() - t0)
    print(f"{len(jobs)} runs", flush=True)
    with ThreadPoolExecutor(a.jobs) as ex:
        for (g, v), st, dt in ex.map(run, jobs): print(f"  {g['tag']} {v}: {st} {dt:.0f}s", flush=True)

    anc = load_obo(f"{REF}/go.obo")
    defs = json.load(open(a.defs))["modules"]
    rows, mrows = [], []
    for g in GENOMES:
        truth = load_truth(f"{TRUTH}/{g['tag']}_truth_all.tsv")
        kegg_complete = {mid for mid, m in defs.items() if g["org"] in m.get("complete_organisms", [])}
        for v in variants:
            out_file = os.path.join(a.out_dir, f"{g['tag']}_{v}_integrated.tsv")
            if not os.path.exists(out_file): continue
            r = dict(genome=g["tag"], org=g["org"], variant=v, kegg_complete_modules=len(kegg_complete))
            r.update(fmax(g["tag"], out_file, a.n_bootstrap))
            # coherence rescues
            rp = out_file + ".rescues.tsv"; n = ct = cs = unk = 0
            if os.path.exists(rp):
                with open(rp) as fh:
                    for x in csv.DictReader(fh, delimiter="\t"):
                        if x["prior"] != "coherence": continue
                        n += 1; tt = truth.get(x["protein_id"])
                        if not tt: unk += 1; continue
                        prop = set().union(*(anc(t) for t in tt))
                        if x["function_id"] in prop: ct += 1
                        step_terms = re.findall(r"GO:\d{7}", x.get("evidence", "")) or [x["function_id"]]
                        if any(s in prop for s in step_terms) or x["function_id"] in prop: cs += 1
            r.update(coh_rescues=n, coh_rescues_correct_term=ct, coh_rescues_correct_step=cs, coh_rescues_unknown_truth=unk)
            # module-level agreement with KEGG
            mp = out_file + ".modules.json"
            if os.path.exists(mp):
                mj = json.load(open(mp)); s = mj["summary"]
                r.update(module_coherence_before=s["module_coherence_before"], module_coherence_after=s["module_coherence_after"],
                         modules_activated_after=s["modules_activated_after"], modules_complete_before=s["modules_complete_before"], modules_complete_after=s["modules_complete_after"])
                for phase in ("before", "after"):
                    tp = fp = fn = act_tp = act_fp = 0
                    for m in mj["modules"]:
                        if m["signature"]: continue
                        if phase == "after": complete = m["complete"]; activated = m["activated"]
                        else:
                            complete = m["completeness_before"] is not None and m["completeness_before"] >= 0.999; activated = m["activated_before"]
                        ref = m["id"] in kegg_complete
                        tp += complete and ref; fp += complete and not ref; fn += (not complete) and ref
                        act_tp += activated and ref; act_fp += activated and not ref
                        if phase == "after":
                            mrows.append(dict(genome=g["tag"], variant=v, module=m["id"], name=m["name"], kegg_complete=int(ref), activated=int(m["activated"]),
                                              completeness_before=m["completeness_before"], completeness_after=m["completeness"], n_required=m["n_required"],
                                              n_unmappable=m["n_unmappable"], n_evidence=m["n_evidence"], n_changes=len(m.get("changes", []))))
                    r[f"complete_tp_{phase}"] = tp; r[f"complete_fp_{phase}"] = fp; r[f"complete_fn_{phase}"] = fn
                    r[f"complete_precision_{phase}"] = tp / (tp + fp) if tp + fp else float("nan"); r[f"complete_recall_{phase}"] = tp / (tp + fn) if tp + fn else float("nan")
                    r[f"activated_tp_{phase}"] = act_tp; r[f"activated_fp_{phase}"] = act_fp
            rows.append(r); print(f"  eval {g['tag']} {v}: fmax_cafa_all={r.get('fmax_cafa_all')} rescues={r.get('coh_rescues')} complete_after={r.get('modules_complete_after')}", flush=True)
    for name, data in (("runs.tsv", rows), ("modules.tsv", mrows)):
        cols = sorted({k for x in data for k in x}, key=lambda k: (k not in ("genome", "org", "variant", "module"), k))
        with open(os.path.join(a.out_dir, name), "w") as fh:
            fh.write("\t".join(cols) + "\n")
            for x in data: fh.write("\t".join(str(x.get(k, "")) for k in cols) + "\n")
    with open(os.path.join(a.out_dir, "summary.tsv"), "w") as fh:
        keys = ["fmax_micro_exp", "fmax_cafa_exp", "fmax_micro_all", "fmax_cafa_all", "coh_rescues", "coh_rescues_correct_term", "coh_rescues_correct_step", "coh_rescues_unknown_truth",
                "modules_complete_before", "modules_complete_after", "complete_tp_after", "complete_fp_after", "complete_fn_after", "complete_tp_before", "complete_fp_before", "complete_fn_before"]
        fh.write("variant\tn_genomes\t" + "\t".join(keys) + "\tcomplete_precision_before\tcomplete_recall_before\tcomplete_precision_after\tcomplete_recall_after\trescue_precision_term\trescue_precision_step\n")
        for v in variants:
            sel = [x for x in rows if x["variant"] == v]
            if not sel: continue
            vals = []
            for k in keys:
                xs = [x[k] for x in sel if k in x and x[k] is not None]
                vals.append(f"{sum(xs) / len(xs):.4f}" if k.startswith("fmax") and xs else (str(sum(xs)) if xs else ""))
            def pr(phase):
                tp = sum(x.get(f"complete_tp_{phase}", 0) for x in sel); fp = sum(x.get(f"complete_fp_{phase}", 0) for x in sel); fn = sum(x.get(f"complete_fn_{phase}", 0) for x in sel)
                return (f"{tp / (tp + fp):.3f}" if tp + fp else "", f"{tp / (tp + fn):.3f}" if tp + fn else "")
            judged = sum(x["coh_rescues"] - x["coh_rescues_unknown_truth"] for x in sel)
            rt = f"{sum(x['coh_rescues_correct_term'] for x in sel) / judged:.3f}" if judged else ""; rs = f"{sum(x['coh_rescues_correct_step'] for x in sel) / judged:.3f}" if judged else ""
            fh.write(f"{v}\t{len(sel)}\t" + "\t".join(vals) + "\t" + "\t".join(pr("before") + pr("after")) + f"\t{rt}\t{rs}\n")
    print(open(os.path.join(a.out_dir, "summary.tsv")).read(), flush=True)


if __name__ == "__main__":
    main()
