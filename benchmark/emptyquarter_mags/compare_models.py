#!/usr/bin/env python3
"""Compare gapsmith and gapseq models of the same bins (validation of the gapsmith substitution).

For every bin directory under --work that holds a gapseq model (gapseq/<mag>.xml from the HQ run or
the gapseq fallback) and a gapsmith model (gapseq_gapsmith/run/<mag>-filled.xml, written by
02_gapsmith.sbatch with FORCE_COMPUTE=1), reports reaction-set sizes and Jaccard indices of the draft
and gap-filled networks, of the gap-filled additions (filled minus draft), and of the predicted
pathways (Pathways.tbl, Prediction == true). Reaction ids are normalised to the SEED id
(leading ``R_`` and compartment suffix removed).

  compare_models.py --work DIR [--gapseq-out DIR] [--out compare.tsv]
"""
import argparse
import csv
import glob
import gzip
import os
import re

RX = re.compile(r'<reaction[^>]*\sid="([^"]+)"')


def reactions(path):
    if not path or not os.path.exists(path):
        return None
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt", errors="replace") as f:
        ids = RX.findall(f.read())
    out = set()
    for i in ids:
        i = re.sub(r"^R_", "", i)
        i = re.sub(r"_[a-z]\d?$", "", i)
        out.add(i)
    return out


def pathways(path):
    if not path or not os.path.exists(path):
        return None
    out = set()
    with open(path, errors="replace") as f:
        hdr = None
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if hdr is None:
                hdr = parts
                continue
            r = dict(zip(hdr, parts))
            if str(r.get("Prediction", "")).lower() == "true":
                out.add(r.get("ID", ""))
    return out


def jac(a, b):
    return 1.0 if not a and not b else round(len(a & b) / len(a | b), 4)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", required=True)
    ap.add_argument("--gapseq-out", default="")
    ap.add_argument("--out", default="")
    ap.add_argument("--tag", default="", help="gapsmith output set suffix (gapseq_gapsmith<tag>)")
    a = ap.parse_args()
    rows = []
    for d in sorted(glob.glob(os.path.join(a.work, "*"))):
        mag = os.path.basename(d)
        gq = os.path.join(d, "gapseq", f"{mag}.xml")
        if not os.path.exists(gq) and a.gapseq_out:
            gq = os.path.join(a.gapseq_out, mag, f"{mag}.xml")
        gq_draft = re.sub(r"\.xml$", "-draft.xml", gq)
        gq_pw = re.sub(r"\.xml$", "-all-Pathways.tbl", gq)
        gm_dir = os.path.join(d, "gapseq_gapsmith" + a.tag, "run")
        gm = os.path.join(gm_dir, f"{mag}-filled.xml")
        gm_draft = os.path.join(gm_dir, f"{mag}-draft.xml")
        gm_pw = os.path.join(gm_dir, f"{mag}-all-Pathways.tbl")
        A, Ad, B, Bd = reactions(gq), reactions(gq_draft), reactions(gm), reactions(gm_draft)
        if A is None or B is None:
            continue
        Ad = Ad or set()
        Bd = Bd or set()
        Pa, Pb = pathways(gq_pw) or set(), pathways(gm_pw) or set()
        rows.append({
            "mag_id": mag,
            "gapseq_draft_rxns": len(Ad), "gapsmith_draft_rxns": len(Bd), "draft_jaccard": jac(Ad, Bd),
            "gapseq_filled_rxns": len(A), "gapsmith_filled_rxns": len(B), "filled_jaccard": jac(A, B),
            "gapseq_added": len(A - Ad), "gapsmith_added": len(B - Bd), "added_jaccard": jac(A - Ad, B - Bd),
            "gapseq_pathways": len(Pa), "gapsmith_pathways": len(Pb), "pathway_jaccard": jac(Pa, Pb),
        })
        print(rows[-1])
    if not rows:
        print("no bin with both models")
        return
    out = a.out or os.path.join(a.work, "..", "compare_models.tsv")
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print("wrote", out)


if __name__ == "__main__":
    main()
