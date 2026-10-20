#!/usr/bin/env python3
"""Add EC claims to a claims JSONL so that Stage 2B (StoichiometricLocusRefiner) and the gap prior
have enzyme-level evidence. GSPA's Stage 2B reads claims of function_type "EC" only.

Two sources, both deterministic from evidence the channels already produced:
  * ec2go: a GO claim whose term is the ec2go image of exactly one complete EC number (four fields, no
    "-") yields an EC claim with the same source and raw score (metadata via=ec2go, go=<term>).
  * Swiss-Prot: a DIAMOND hit to a Swiss-Prot entry with EC numbers yields an EC claim per complete EC,
    source "diamond", raw score = percent identity / 100, max over hits (metadata target=<accession>).
Duplicates (protein, EC, source) keep the highest raw score. EC ids are written as "EC:1.1.1.1".

  add_ec_claims.py --claims claims.jsonl --diamond diamond.tsv[.gz] --ec2go ec2go.txt --sprot-ec sprot_ec.tsv
                   [--out claims_with_ec.jsonl]   (default: append to --claims in place)
"""
import argparse
import gzip
import json
import os
import re
from collections import defaultdict

EC_RE = re.compile(r"^EC:(\d+\.\d+\.\d+\.\d+)$")


def load_go2ec(path):
    """GO term -> EC (complete numbers only, and only when a single EC maps to the term)."""
    go2ecs = defaultdict(set)
    with open(path) as f:
        for line in f:
            if line.startswith("!") or not line.strip():
                continue
            m = re.match(r"^(EC:\S+)\s*>\s*GO:[^;]+;\s*(GO:\d{7})", line)
            if not m:
                continue
            ec, go = m.groups()
            if EC_RE.match(ec):
                go2ecs[go].add(ec)
    return {go: next(iter(ecs)) for go, ecs in go2ecs.items() if len(ecs) == 1}


def load_sprot_ec(path):
    acc2ec = defaultdict(set)
    with open(path) as f:
        for line in f:
            acc, ec = line.rstrip("\n").split("\t")[:2]
            if EC_RE.match(ec if ec.startswith("EC:") else "EC:" + ec):
                acc2ec[acc].add(ec if ec.startswith("EC:") else "EC:" + ec)
    return acc2ec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--claims", required=True)
    ap.add_argument("--diamond", default="")
    ap.add_argument("--ec2go", required=True)
    ap.add_argument("--sprot-ec", default="")
    ap.add_argument("--out", default="")
    a = ap.parse_args()
    go2ec = load_go2ec(a.ec2go)
    acc2ec = load_sprot_ec(a.sprot_ec) if a.sprot_ec and os.path.exists(a.sprot_ec) else {}
    best = {}   # (protein, ec, source) -> (raw, meta)
    n_go = 0
    with open(a.claims) as f:
        for line in f:
            rec = json.loads(line)
            if rec.get("function_type") != "GO":
                continue
            ec = go2ec.get(rec["function_id"])
            if not ec:
                continue
            k = (rec["protein_id"], ec, rec["source"])
            raw = float(rec.get("raw_score", 0.0))
            if k not in best or raw > best[k][0]:
                best[k] = (raw, {"via": "ec2go", "go": rec["function_id"]})
            n_go += 1
    n_dia = 0
    if a.diamond and os.path.exists(a.diamond) and acc2ec:
        op = gzip.open if a.diamond.endswith(".gz") else open
        with op(a.diamond, "rt") as f:
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) < 3:
                    continue
                q, t, pid = p[0], p[1], float(p[2])
                for ec in acc2ec.get(t, ()):
                    k = (q, ec, "diamond")
                    raw = min(1.0, pid / 100.0)
                    if k not in best or raw > best[k][0]:
                        best[k] = (raw, {"via": "swissprot_ec", "target": t, "pident": pid})
                    n_dia += 1
    out_path = a.out or a.claims
    mode = "a" if out_path == a.claims else "w"
    with open(out_path, mode) as out:
        if mode == "w":
            with open(a.claims) as f:
                for line in f:
                    out.write(line)
        for (pid, ec, src), (raw, meta) in sorted(best.items()):
            out.write(json.dumps({"protein_id": pid, "function_type": "EC", "function_id": ec, "go_aspect": "",
                                  "source": src, "raw_score": round(raw, 4), "metadata": meta}, separators=(",", ":")) + "\n")
    print(f"[ec] {len(best)} EC claims added ({n_go} GO claims mapped via ec2go, {n_dia} Swiss-Prot EC transfers; "
          f"{len(go2ec)} GO terms with a unique complete EC, {len(acc2ec)} Swiss-Prot entries with EC)")


if __name__ == "__main__":
    main()
