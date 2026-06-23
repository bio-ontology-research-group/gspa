#!/usr/bin/env python3
"""Build a leak-free no-knowledge ground truth (fixes the GAF-date contamination).

THE BUG (found 2026-06-23). `build_groundtruth.py` decides whether a (protein,
aspect) is a *novel post-t0* test target from the GOA `date` column. That column
is the annotation's **last-modified** date, not its **creation** date, so an
annotation that existed *before* t0 but was touched afterwards gets re-dated past
t0 and is counted as a brand-new discovery. Result: 95% of the reconstruction's
"no-knowledge MF" targets already had that protein's MF term in the pre-t0
`train_terms.tsv`. Homology/PPI methods (esp. `net`) then "recover" pre-known
labels and look far stronger than they are on genuinely novel proteins
(net MF 0.803 -> 0.430 once the contamination is removed).

THE FIX. Trust the canonical pre-t0 label set (`train_terms.tsv`, the file every
competitor was given) instead of the unreliable GAF date. A (protein, aspect) is
genuinely no-knowledge iff the protein has **no `train_terms` entry in that
aspect**. This sidesteps GAF dates entirely.

Input:
  gt_no.tsv        the (contaminated) no-knowledge GT (protein<TAB>GO)
  train_terms.tsv  CAFA6 official pre-t0 labels (EntryID<TAB>term<TAB>aspect[F/P/C])
  go.obo           GO ontology — authoritative term namespaces (what cafaeval uses)
Output:
  gt_no_clean.tsv  per-aspect-filtered GT (the faithful no-knowledge target set)
Optionally (--novel-proteins) also writes gt_no_novelprot.tsv = only proteins
entirely absent from train_terms (the strictest, smallest novel set).
"""
from __future__ import annotations

import argparse
import sys
from collections import defaultdict

ASPMAP = {"F": "MF", "P": "BP", "C": "CC"}
NS = {"molecular_function": "MF", "biological_process": "BP", "cellular_component": "CC"}


def load_obo_aspect(path):
    """term -> MF/BP/CC from the OBO `namespace:` field. This is the *authoritative*
    aspect (the one cafaeval uses). NB: deriving aspect from a go-dag closure is WRONG
    for terms with cross-aspect `part_of` ancestors (e.g. GO:0005216 ion channel
    activity, an MF term, reaches the BP root through `part_of` transport) — the
    closure reaches >1 root and the answer becomes order-dependent."""
    asp = {}
    cur = None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line == "[Term]":
                cur = None
            elif line.startswith("id: GO:"):
                cur = line[4:].strip()
            elif line.startswith("namespace:") and cur:
                asp[cur] = NS.get(line.split(":", 1)[1].strip())
    return asp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gt-no", required=True, help="contaminated gt_no.tsv (protein, GO)")
    ap.add_argument("--train-terms", required=True,
                    help="train_terms.tsv (EntryID, term, aspect F/P/C)")
    ap.add_argument("--obo", required=True,
                    help="go.obo / go-basic.obo — authoritative term namespaces")
    ap.add_argument("--out", required=True, help="gt_no_clean.tsv")
    ap.add_argument("--novel-proteins-out",
                    help="also write the strict novel-protein GT (absent from train_terms)")
    args = ap.parse_args()

    term_asp = load_obo_aspect(args.obo)

    def aspect_of(term, _ignored=None):
        return term_asp.get(term)

    train_asp = defaultdict(set)   # protein -> {aspects present pre-t0}
    train_any = set()
    with open(args.train_terms) as fh:
        next(fh, None)             # header
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if not c or not c[0]:
                continue
            train_any.add(c[0])
            if len(c) >= 3 and c[2] in ASPMAP:
                train_asp[c[0]].add(ASPMAP[c[2]])

    rows = [tuple(l.rstrip("\n").split("\t")[:2]) for l in open(args.gt_no) if l.strip()]
    clean = [(p, t) for p, t in rows if aspect_of(t) not in train_asp.get(p, set())]

    with open(args.out, "w") as out:
        out.writelines(f"{p}\t{t}\n" for p, t in clean)

    def per_aspect(rws):
        c = defaultdict(lambda: [0, set()])
        for p, t in rws:
            a = aspect_of(t)
            if a:
                c[a][0] += 1
                c[a][1].add(p)
        return {a: (c[a][0], len(c[a][1])) for a in ("MF", "BP", "CC")}

    print(f"dirty  gt_no       : {len(rows)} rows, {len({p for p,_ in rows})} prot  {per_aspect(rows)}",
          file=sys.stderr)
    print(f"clean  (per-aspect): {len(clean)} rows, {len({p for p,_ in clean})} prot  {per_aspect(clean)}",
          file=sys.stderr)
    print(f"wrote {args.out}", file=sys.stderr)

    if args.novel_proteins_out:
        novel = [(p, t) for p, t in rows if p not in train_any]
        with open(args.novel_proteins_out, "w") as out:
            out.writelines(f"{p}\t{t}\n" for p, t in novel)
        print(f"novel-protein      : {len(novel)} rows, {len({p for p,_ in novel})} prot  {per_aspect(novel)}",
              file=sys.stderr)
        print(f"wrote {args.novel_proteins_out}", file=sys.stderr)


if __name__ == "__main__":
    main()
