#!/usr/bin/env python3
"""Operon membership of role-truth loci, from the genome layouts already built.

The grouping is the production distance rule (``OperonEnsemble`` "distance"
predictor, ``predictors.operons.maxIntergenicDistance: 300`` in the bundled
bacteria/archaea presets): consecutive loci on one contig and strand join when
the intergenic gap is at most 300 bp, transitively along the chain. The strict
variant (50 bp) is reported beside it and is never the primary tag.

Nothing here reads a predictor score, so an operon tag cannot be chosen to suit
a result. ``tag`` adds an ``operon_context`` record to each registry case;
``operons`` writes the operon TSV that ``gspa integrate --operons`` consumes.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

SCHEMA = "gspa.role-truth-operon-context/1"
DISTANCE_MAX_GAP = 300
STRICT_MAX_GAP = 50


def sha256(path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_layout(path) -> list[dict]:
    with Path(path).open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    for row in rows:
        for field in ("protein_id", "contig", "start", "end", "strand"):
            if not row.get(field):
                raise ValueError(f"{path}: layout row missing {field}")
        row["start"], row["end"] = int(row["start"]), int(row["end"])
    rows.sort(key=lambda r: (r["contig"], r["start"], r["end"], r["protein_id"]))
    return rows


def group_operons(rows, max_gap=DISTANCE_MAX_GAP):
    """Return (protein -> group id, group id -> ordered members, protein -> gap)."""
    groups, members, gaps = {}, {}, {}
    index, previous = 0, None
    for row in rows:
        joined = (previous is not None and previous["contig"] == row["contig"]
                  and previous["strand"] == row["strand"]
                  and row["start"] - previous["end"] - 1 <= max_gap)
        if joined:
            gaps[row["protein_id"]] = row["start"] - previous["end"] - 1
        else:
            index += 1
        name = f"op_{index:05d}"
        groups[row["protein_id"]] = name
        members.setdefault(name, []).append(row["protein_id"])
        previous = row
    return groups, members, gaps


def layout_index(layouts: dict) -> dict:
    """{genome_id: {'path', 'sha256', 'loci', 'distance', 'strict'}} for each layout."""
    index = {}
    for genome_id, path in layouts.items():
        rows = read_layout(path)
        distance = group_operons(rows, DISTANCE_MAX_GAP)
        strict = group_operons(rows, STRICT_MAX_GAP)
        index[genome_id] = {"path": str(path), "sha256": sha256(path), "n_loci": len(rows),
                            "loci": {r["protein_id"]: r for r in rows},
                            "distance": distance, "strict": strict}
    return index


def case_context(case: dict, index: dict) -> dict:
    """Operon membership of every member of one case, and the case-level relation."""
    genome = index.get(case["genome_id"])
    if genome is None:
        return {"schema": SCHEMA, "status": "no_layout", "relation": "unplaced",
                "rule": f"same strand, intergenic gap <= {DISTANCE_MAX_GAP} bp, transitive",
                "members": {}, "unplaced_members": [m["protein_id"] for m in case["members"]]}
    distance_groups, distance_members, gaps = genome["distance"]
    strict_groups = genome["strict"][0]
    members, unplaced = {}, []
    for member in case["members"]:
        pid = member["protein_id"]
        locus = genome["loci"].get(pid)
        if locus is None:
            unplaced.append(pid)
            members[pid] = {"label": member["label"], "placed": False}
            continue
        members[pid] = {"label": member["label"], "placed": True, "contig": locus["contig"],
                        "start": locus["start"], "end": locus["end"], "strand": locus["strand"],
                        "operon": distance_groups[pid], "strict_operon": strict_groups[pid],
                        "operon_size": len(distance_members[distance_groups[pid]]),
                        "gap_to_previous_in_operon": gaps.get(pid)}
    actives = [m["protein_id"] for m in case["members"] if m["label"] == 1]
    negatives = [m["protein_id"] for m in case["members"] if m["label"] == 0]
    pairs = []
    # A negative-only locus has no active: its relation is between its negatives.
    pair_sources = [(a, n) for a in actives for n in negatives] if actives else [
        (negatives[i], negatives[j]) for i in range(len(negatives)) for j in range(i + 1, len(negatives))]
    for active, negative in pair_sources:
        if True:
            both = members[active].get("placed") and members[negative].get("placed")
            relation = "unplaced"
            if both:
                same = members[active]["operon"] == members[negative]["operon"]
                relation = "same_operon" if same else "different_operon"
            pairs.append({"active": active, "negative": negative, "relation": relation,
                          "same_strict_operon": bool(both) and members[active]["strict_operon"]
                          == members[negative]["strict_operon"]})
    if any(p["relation"] == "same_operon" for p in pairs):
        relation = "same_operon"
    elif pairs and all(p["relation"] == "different_operon" for p in pairs):
        relation = "different_operon"
    else:
        relation = "unplaced"
    return {"schema": SCHEMA, "status": "ok", "relation": relation,
            "rule": f"same strand, intergenic gap <= {DISTANCE_MAX_GAP} bp, transitive",
            "layout": {"path": genome["path"], "sha256": genome["sha256"], "n_loci": genome["n_loci"]},
            "members": members, "active_negative_pairs": pairs, "unplaced_members": unplaced}


def read_jsonl(path) -> list[dict]:
    return [json.loads(line) for line in Path(path).read_text().splitlines() if line.strip()]


def write_jsonl(path, rows):
    Path(path).write_text("".join(json.dumps(r, sort_keys=True) + "\n" for r in rows))


def tag(registry, layouts: dict, out_tags, in_place=False, negative_only=False):
    cases = read_jsonl(registry)
    index = layout_index(layouts)
    tags = []
    for case in cases:
        context = case_context(case, index)
        case["operon_context"] = context
        tags.append({"case_id": case["case_id"], "genome_id": case["genome_id"],
                     "function_id": case["function_id"], "status": case["status"],
                     "negative_only": bool(negative_only), "operon_context": context})
    if in_place:
        write_jsonl(registry, cases)
    if out_tags:
        write_jsonl(out_tags, tags)
    return tags


def operons_tsv(layout, out, max_gap=DISTANCE_MAX_GAP, min_size=2):
    rows = read_layout(layout)
    _, members, _ = group_operons(rows, max_gap)
    kept = [m for m in members.values() if len(m) >= min_size]
    Path(out).write_text("".join("\t".join(m) + "\n" for m in kept))
    return {"layout": str(layout), "layout_sha256": sha256(layout), "operons": len(kept),
            "proteins_in_operons": sum(len(m) for m in kept), "n_loci": len(rows),
            "max_gap": max_gap, "out": str(out), "out_sha256": sha256(out)}


def parse_layouts(values) -> dict:
    layouts = {}
    for item in values or []:
        if item.count("=") != 1:
            raise SystemExit("--layout entries must be genome_id=path")
        genome_id, path = item.split("=")
        layouts[genome_id] = Path(path)
    return layouts


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    tagger = sub.add_parser("tag")
    tagger.add_argument("--registry", required=True, type=Path)
    tagger.add_argument("--layout", action="append", required=True)
    tagger.add_argument("--out-tags", type=Path)
    tagger.add_argument("--in-place", action="store_true")
    tagger.add_argument("--negative-only", action="store_true")
    ops = sub.add_parser("operons")
    ops.add_argument("--layout", required=True, type=Path)
    ops.add_argument("--out", required=True, type=Path)
    args = parser.parse_args(argv)
    if args.command == "tag":
        tags = tag(args.registry, parse_layouts(args.layout), args.out_tags,
                   in_place=args.in_place, negative_only=args.negative_only)
        counts = {}
        for row in tags:
            counts[row["operon_context"]["relation"]] = counts.get(row["operon_context"]["relation"], 0) + 1
        print(json.dumps({"cases": len(tags), "relations": counts}, sort_keys=True))
    else:
        print(json.dumps(operons_tsv(args.layout, args.out), sort_keys=True))


if __name__ == "__main__":
    main()
