#!/usr/bin/env python3
"""P03 step 5: count surviving injected cross-clade claims per genome and method.

(a) threshold only (no GSPA): an injected (protein, term) survives if the max
    raw score over the method's Stage 0 claims for that pair is >= 0.5
(b) GSPA Stage 1+2A (A9 relax) applied on top: survives if posterior_prob >= 0.5
    in the integrated output (same rule as 04_run_perturbation.evaluate_perturbation)
Writes {REPRO}/perturbation/panel13_perturbation.tsv
"""
import csv
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p03_common import REPRO, TAGS, CLAIM_SOURCE, injected_path, perturbed_claims_path, perturbed_out_path


def main():
    rows = []
    for tag in TAGS:
        recs = list(csv.DictReader(open(injected_path(tag)), delimiter="\t"))
        inj = {(r["protein_id"], r["go_id"]): r["category"] for r in recs}
        for m in ["gspa3"] + list(CLAIM_SOURCE):
            mx = {}
            with open(perturbed_claims_path(tag, m)) as f:
                for line in f:
                    c = json.loads(line)
                    k = (c["protein_id"], c["function_id"])
                    if k in inj:
                        mx[k] = max(mx.get(k, 0.0), float(c["raw_score"]))
            surv_a = {k for k, s in mx.items() if s >= 0.5}
            out = perturbed_out_path(tag, m)
            surv_b, status = set(), "ok"
            if os.path.exists(out):
                with open(out) as f:
                    h = next(f).rstrip("\n").split("\t")
                    ip, ifn, ipr = h.index("protein_id"), h.index("function_id"), h.index("posterior_prob")
                    for line in f:
                        p = line.rstrip("\n").split("\t")
                        k = (p[ip], p[ifn])
                        if k in inj:
                            try:
                                if float(p[ipr]) >= 0.5:
                                    surv_b.add(k)
                            except ValueError:
                                pass
            else:
                status = "missing_integrated_output"
            rows.append({"tag": tag, "method": m, "n_injected": len(inj),
                         "surv_threshold_only": len(surv_a),
                         "surv_gspa_A9": len(surv_b) if status == "ok" else "",
                         "surv_gspa_A9_anatomical": sum(1 for k in surv_b if inj[k] == "anatomical") if status == "ok" else "",
                         "surv_gspa_A9_taxon": sum(1 for k in surv_b if inj[k] == "taxon_constraint") if status == "ok" else "",
                         "status": status})
    cols = list(rows[0])
    with open(f"{REPRO}/perturbation/panel13_perturbation.tsv", "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")
    print("written", len(rows))


if __name__ == "__main__":
    main()
