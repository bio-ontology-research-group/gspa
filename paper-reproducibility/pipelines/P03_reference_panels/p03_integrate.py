#!/usr/bin/env python3
"""P03 step 2: `gspa integrate` (gspa-unified.jar 2.1.0), A9_stage2_relax flags
(benchmark/unified/pillar1_2/01_run_ablation_grid.py):
  --enable-priors essentiality,coherence,consistency,anatomical,gap_filling,genomic_context
  --essential-profile <kingdom> --taxon <taxon> --pathways/--modules/--ec2go (KEGG)
  --enforce-stage2 --enforce-consistency-mode relax  [--operons] [--gaps]  --lite
Jobs (one shared-memory batch manifest):
  gspa3 clean       13  (3-channel claims)
  gspa3 perturbed   13  (3-channel claims + design-v2 injection)
  <m> perturbed     5 x 13 (external method claims + injection, m in pgap/pfam/interproscan/diamond/mdf)
"""
import argparse
import os
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p03_common import (JAVA, JAR, REF, REPRO, GINFO, TAGS, CLAIM_SOURCE, paths, pred_path,
                        perturbed_claims_path, perturbed_out_path, ensure_dirs)

PRIORS = "essentiality,coherence,consistency,anatomical,gap_filling,genomic_context"


def row(gid, claims, out, tag):
    P = paths(tag)
    ops = P["operons"] if P["operons"] and os.path.exists(P["operons"]) and os.path.getsize(P["operons"]) > 0 else "-"
    gaps = P["gaps"] if P["gaps"] and os.path.exists(P["gaps"]) and os.path.getsize(P["gaps"]) > 0 else "-"
    return "\t".join([gid, claims, out, GINFO[tag]["taxon"], GINFO[tag]["k"], ops, gaps, "-",
                      "true", "relax", PRIORS])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--threads", type=int, default=12)
    ap.add_argument("--xmx", default="64g")
    ap.add_argument("--only", default=None, help="clean|perturbed")
    a = ap.parse_args()
    ensure_dirs()
    rows = []
    for tag in TAGS:
        P = paths(tag)
        if a.only in (None, "clean"):
            out = pred_path(tag, "gspa3")
            if not os.path.exists(out):
                rows.append(row(f"{tag}_gspa3_A9", P["claims3"], out, tag))
        if a.only in (None, "perturbed"):
            for m in ["gspa3"] + list(CLAIM_SOURCE):
                out = perturbed_out_path(tag, m)
                if not os.path.exists(out):
                    rows.append(row(f"{tag}_{m}_A9_perturb", perturbed_claims_path(tag, m), out, tag))
    if not rows:
        print("nothing to do")
        return
    man = f"{REPRO}/integrate_batch_manifest_{a.only or 'all'}.tsv"
    with open(man, "w") as f:
        f.write("genome_id\tclaims_path\tout_path\ttaxon\tessential_profile\toperons_path\tgaps_path\t"
                "max_iter\tenforce_stage2\tenforce_mode\tenable_priors\n")
        f.write("\n".join(rows) + "\n")
    cmd = [JAVA, f"-Xmx{a.xmx}", "-jar", JAR, "integrate",
           "--batch-manifest", man, "--batch-threads", str(a.threads),
           "--enable-priors", PRIORS, "--go-owl", f"{REF}/go.owl",
           "--pathways", f"{REF}/kegg_pathways.tsv", "--modules", f"{REF}/kegg_modules.tsv",
           "--ec2go", f"{REF}/ec2go.txt", "--lite"]
    print(" ".join(cmd), flush=True)
    t0 = time.time()
    rc = subprocess.run(cmd).returncode
    print(f"integrate rc={rc} wall={time.time() - t0:.1f}s jobs={len(rows)}", flush=True)
    sys.exit(rc)


if __name__ == "__main__":
    main()
