#!/usr/bin/env python3
"""Ancestor-propagate a GOA truth TSV via is_a + part_of edges.

Input: TSV with header ``protein_id<TAB>aspect<TAB>function_id``.
Output: same format, but each (protein, term) row is expanded to include
every ancestor of ``term`` through is_a + part_of. Aspect column is
preserved (ancestors across aspects would be dropped if we restricted to
the term's aspect — so we do NOT restrict; ancestors from other aspects
are valid reflexes of the annotation and get their own aspect).

Usage::
    propagate_truth.py --go-obo go-basic.obo --in truth.tsv --out truth_prop.tsv
"""
from __future__ import annotations

import argparse
from collections import defaultdict, deque
from pathlib import Path


def parse_obo(path: Path):
    """Return (parents, aspect) maps.

    parents[child] = set of direct parents via is_a or part_of.
    aspect[term]  = MF / BP / CC.
    obsolete terms and alt_ids are resolved.
    """
    parents: dict[str, set[str]] = defaultdict(set)
    aspect: dict[str, str] = {}
    alt_to_primary: dict[str, str] = {}
    ns_map = {"molecular_function": "MF",
              "biological_process": "BP",
              "cellular_component": "CC"}
    in_term = False
    cid = None
    is_obsolete = False
    replaced_by = None
    buf_parents: set[str] = set()
    buf_alts: set[str] = set()
    buf_asp: str | None = None

    def flush():
        nonlocal cid, is_obsolete, replaced_by, buf_parents, buf_alts, buf_asp
        if cid and not is_obsolete:
            parents[cid] = set(buf_parents)
            if buf_asp:
                aspect[cid] = buf_asp
            for alt in buf_alts:
                alt_to_primary[alt] = cid
        elif cid and is_obsolete and replaced_by:
            alt_to_primary[cid] = replaced_by
        cid = None
        is_obsolete = False
        replaced_by = None
        buf_parents = set()
        buf_alts = set()
        buf_asp = None

    with path.open() as fh:
        for raw in fh:
            line = raw.rstrip()
            if line.startswith("[Term]"):
                flush()
                in_term = True
                continue
            if line.startswith("[") and line != "[Term]":
                flush()
                in_term = False
                continue
            if not in_term:
                continue
            if line.startswith("id: "):
                cid = line[4:].strip()
            elif line.startswith("alt_id: "):
                buf_alts.add(line[8:].strip())
            elif line.startswith("is_obsolete: true"):
                is_obsolete = True
            elif line.startswith("replaced_by: "):
                replaced_by = line[13:].strip()
            elif line.startswith("namespace: "):
                ns = line[11:].strip()
                buf_asp = ns_map.get(ns)
            elif line.startswith("is_a: "):
                tgt = line[6:].split(" ", 1)[0].strip()
                buf_parents.add(tgt)
            elif line.startswith("relationship: part_of "):
                tgt = line[22:].split(" ", 1)[0].strip()
                buf_parents.add(tgt)
        flush()
    return parents, aspect, alt_to_primary


def compute_ancestors(parents: dict[str, set[str]]) -> dict[str, set[str]]:
    cache: dict[str, set[str]] = {}

    def anc(t: str) -> set[str]:
        if t in cache:
            return cache[t]
        seen: set[str] = set()
        stack = deque([t])
        while stack:
            cur = stack.popleft()
            for p in parents.get(cur, ()):
                if p not in seen:
                    seen.add(p)
                    stack.append(p)
        cache[t] = seen
        return seen

    for t in list(parents):
        anc(t)
    return cache


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--go-obo", required=True)
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    parents, aspect, alt_to_primary = parse_obo(Path(args.go_obo))
    ancestors = compute_ancestors(parents)

    out_rows: set[tuple[str, str, str]] = set()
    n_in = 0
    n_unknown = 0
    with open(args.inp) as fh:
        header = next(fh).rstrip("\n").split("\t")
        # header: protein_id aspect function_id
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            pid, asp, term = parts[0], parts[1], parts[2]
            term = alt_to_primary.get(term, term)
            if term not in parents:
                n_unknown += 1
            n_in += 1
            # add self (resolved alt)
            out_rows.add((pid, aspect.get(term, asp), term))
            # add ancestors; each ancestor gets its own aspect
            for a in ancestors.get(term, ()):
                out_rows.add((pid, aspect.get(a, asp), a))

    with open(args.out, "w") as fh:
        fh.write("protein_id\taspect\tfunction_id\n")
        for pid, asp, term in sorted(out_rows):
            fh.write(f"{pid}\t{asp}\t{term}\n")
    print(f"in={n_in} unknown={n_unknown} out={len(out_rows)}")


if __name__ == "__main__":
    main()
