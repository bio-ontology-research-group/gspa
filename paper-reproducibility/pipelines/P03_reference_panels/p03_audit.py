#!/usr/bin/env python3
"""P03 step 3: structural audit (`gspa evaluate`, as in
benchmark/unified/pillar1_2/03_run_gaef_eval.py) on every method's output.

  make_synth_gff_gaf.py --threshold 0.5  (annotations with score >= 0.5)
  java -jar gspa-unified.jar evaluate --lite --go-owl go.owl --pathways/--modules/--ec2go
       --batch-manifest (gff, gaf, out_json, kingdom, taxon)
Then parse each *_quality.json into {REPRO}/audit/panel13_audit.tsv
"""
import argparse
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p03_common import (JAVA, JAR, PYTHON, SYNTH_SCRIPT, REF, REPRO, GINFO, TAGS, METHODS,
                        paths, pred_path, out_root)


VARIANT = "t05"
FMAXT = {}


def base(tag, m):
    return f"{out_root(tag)}/audit/{tag}_{m}" + ("" if VARIANT == "t05" else "_fmaxt")


def thr(tag, m):
    return 0.5 if VARIANT == "t05" else FMAXT[(m, tag)]


def synth(tag, m):
    b = base(tag, m)
    if os.path.exists(b + ".gaf") and os.path.getsize(b + ".gaf") > 0:
        return
    subprocess.run([PYTHON, SYNTH_SCRIPT, "--fasta", paths(tag)["faa"], "--tsv", pred_path(tag, m),
                    "--threshold", str(thr(tag, m)), "--gff-out", b + ".gff", "--gaf-out", b + ".gaf"],
                   check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)


def parse(js):
    d = json.load(open(js))
    comp = d.get("completeness", {})
    coh = d.get("coherence", {})
    cons = d.get("consistency", {})
    return {
        "essential_completeness": comp.get("score"),
        "essential_present": comp.get("present_count"),
        "essential_missing": comp.get("missing_count"),
        "pathway_coherence": coh.get("pathway_coherence"),
        "complex_coherence": coh.get("complex_coherence"),
        "process_coherence": coh.get("process_coherence"),
        "consistent": cons.get("consistent"),
        "violation_count": cons.get("violation_count"),
        "composite_score": d.get("summary", {}).get("composite_score"),
        "annotation_coverage": d.get("summary", {}).get("annotation_coverage"),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--xmx", default="48g")
    ap.add_argument("--parse-only", action="store_true")
    ap.add_argument("--variant", choices=["t05", "fmaxt"], default="t05",
                    help="t05: score>=0.5 (03_run_gaef_eval.py); fmaxt: score>= the method's per-genome "
                         "micro F-max threshold (all-evidence truth) from eval/panel13_per_genome.tsv")
    a = ap.parse_args()
    global VARIANT
    VARIANT = a.variant
    if VARIANT == "fmaxt":
        import csv
        for r in csv.DictReader(open(f"{REPRO}/eval/panel13_per_genome.tsv"), delimiter="\t"):
            if r["truth"] == "all" and r["status"] == "ok":
                FMAXT[(r["method"], r["tag"])] = float(r["fmax_micro_t"])
    tasks = [(t, m) for t in TAGS for m in METHODS if os.path.exists(pred_path(t, m))]
    if not a.parse_only:
        with ThreadPoolExecutor(12) as ex:
            list(ex.map(lambda x: synth(*x), tasks))
        todo = [(t, m) for t, m in tasks if not os.path.exists(base(t, m) + "_quality.json")]
        if todo:
            man = f"{REPRO}/audit_batch_manifest_{VARIANT}.tsv"
            with open(man, "w") as f:
                f.write("gff_path\tgaf_path\tout_json_path\tkingdom\ttaxon\n")
                for t, m in todo:
                    b = base(t, m)
                    f.write(f"{b}.gff\t{b}.gaf\t{b}_quality.json\t{GINFO[t]['k']}\t{GINFO[t]['taxon']}\n")
            cmd = [JAVA, f"-Xmx{a.xmx}", "-jar", JAR, "evaluate", "--lite", "--go-owl", f"{REF}/go.owl",
                   "--pathways", f"{REF}/kegg_pathways.tsv", "--modules", f"{REF}/kegg_modules.tsv",
                   "--ec2go", f"{REF}/ec2go.txt", "--batch-manifest", man]
            t0 = time.time()
            rc = subprocess.run(cmd).returncode
            print(f"evaluate rc={rc} wall={time.time() - t0:.1f}s n={len(todo)}", flush=True)
    cols = ["method", "tag", "variant", "threshold", "essential_completeness", "essential_present", "essential_missing",
            "pathway_coherence", "complex_coherence", "process_coherence", "consistent",
            "violation_count", "composite_score", "annotation_coverage"]
    with open(f"{REPRO}/audit/panel13_audit_{VARIANT}.tsv", "w") as f:
        f.write("\t".join(cols) + "\n")
        for t, m in tasks:
            js = base(t, m) + "_quality.json"
            if not os.path.exists(js):
                f.write(f"{m}\t{t}\t{VARIANT}" + "\t" * (len(cols) - 3) + "\n")
                continue
            r = parse(js)
            r.update({"method": m, "tag": t, "variant": VARIANT, "threshold": thr(t, m)})
            f.write("\t".join("" if r.get(c) is None else str(r[c]) for c in cols) + "\n")
    print("audit table written", flush=True)


if __name__ == "__main__":
    main()
