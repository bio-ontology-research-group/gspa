#!/usr/bin/env python3
"""Round include lists used on 2026-10-19 (deterministic: sha256('42|'+case_id) order per genome).
r1: first 45 PATHWAY_GAP cases per genome + every INCOMPLETE_COMPLEX / BORDERLINE_ESSENTIAL case
r2: next 45 PATHWAY_GAP cases per genome (not in r1)
r3: every PATHWAY_GAP case not in r1
Usage: make_includes.py cases_all.jsonl   (writes r1_include.txt, r2_include.txt, r3_include.txt next to it)"""
import collections, hashlib, json, sys
from pathlib import Path
src = Path(sys.argv[1])
cases = [json.loads(l) for l in open(src) if l.strip()]
key = lambda x: hashlib.sha256(("42|" + x["case_id"]).encode()).hexdigest()
by = collections.defaultdict(list)
for x in cases:
    by[(x["genome"], x["category"])].append(x)
r1, r2 = [], []
for k in sorted(by):
    v = sorted(by[k], key=key)
    if k[1] == "PATHWAY_GAP":
        r1 += [x["case_id"] for x in v[:45]]
        r2 += [x["case_id"] for x in v[45:90]]
    else:
        r1 += [x["case_id"] for x in v]
s1 = set(r1)
r3 = [x["case_id"] for x in cases if x["category"] == "PATHWAY_GAP" and x["case_id"] not in s1]
for name, ids in (("r1", r1), ("r2", r2), ("r3", r3)):
    (src.parent / f"{name}_include.txt").write_text("\n".join(ids) + "\n")
    print(name, len(ids))
