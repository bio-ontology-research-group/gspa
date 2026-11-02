#!/usr/bin/env python3
"""Join Price experimental EC truth to a frozen public case bundle.

Truth rows come from (i) the Price-148 label tables exported with the
release-candidate inputs (complete ECs only) and (ii) the Swiss-Prot
experimental EC table (ECO:0000269 and descendants, exact-sequence matched).
A (genome, candidate protein, target function) is ``positive`` when the function
is one of the protein's labelled ECs, ``negative`` when it is a complete EC
different from every labelled EC of a labelled protein, and ``unknown``
otherwise. GO targets and partial ECs are never matched: the anomaly extractor
keys on GO terms, so the intersection with EC truth is a measurement of what
the production extractor can and cannot be judged on, not an assumption.

The public bundle's sha256 is read and recorded before any label file opens.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path

SCHEMA = "gspa.price-agent-truth/1"
COMPLETE_EC = re.compile(r"^(EC:)?(\d+\.\d+\.\d+\.\d+)$")


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def normalise(function_id):
    match = COMPLETE_EC.match(function_id or "")
    return match.group(2) if match else None


def price_labels(tables: dict) -> list[dict]:
    rows = []
    for genome_id, path in sorted(tables.items()):
        with Path(path).open(newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                ec = normalise(row["ec"])
                if ec:
                    rows.append({"genome_id": genome_id, "protein_id": row["protein_id"], "function_id": ec,
                                 "swissprot_accession": "", "evidence": "Price2018-RBTnSeq",
                                 "match_basis": "price_label_audit", "tier": "experimental"})
    return rows


def swissprot_labels(path) -> list[dict]:
    rows = []
    with Path(path).open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            ec = normalise(row["function_id"])
            if ec:
                rows.append({"genome_id": row["genome_id"], "protein_id": row["protein_id"], "function_id": ec,
                             "swissprot_accession": row.get("swissprot_accession", ""),
                             "evidence": row["evidence"], "match_basis": "exact_sequence", "tier": "experimental"})
    return rows


def label_index(rows) -> dict:
    index: dict = {}
    for row in rows:
        index.setdefault((row["genome_id"], row["protein_id"]), set()).add(row["function_id"])
    return index


def truth_status(index, genome_id, protein_id, target):
    ec = normalise(target)
    labelled = index.get((genome_id, protein_id))
    if ec is None or not labelled:
        return "unknown"
    return "positive" if ec in labelled else "negative"


def join(cases: list[dict], index: dict) -> dict:
    counts = {"positive": 0, "negative": 0, "unknown": 0}
    truth = {}
    for case in cases:
        genome_id = case["genome_id"]
        for candidate in case.get("candidates", []):
            status = truth_status(index, genome_id, candidate["protein_id"], case["target_function_id"])
            counts[status] += 1
            truth.setdefault(case["case_id"], {})[candidate["protein_id"]] = status
    return {"pairs": counts, "cases": len(cases), "truth": truth}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--public", required=True, type=Path, help="frozen public case bundle (JSON list)")
    parser.add_argument("--price-labels", required=True, type=Path, help="directory of <assembly>.ec.tsv")
    parser.add_argument("--swissprot-ecs", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args(argv)
    public_sha = sha256(args.public)
    cases = json.loads(args.public.read_text())
    cases = cases["cases"] if isinstance(cases, dict) else cases
    tables = {p.name.split(".")[0]: p for p in sorted(args.price_labels.glob("*.ec.tsv"))}
    rows = price_labels(tables) + swissprot_labels(args.swissprot_ecs)
    result = {"schema": SCHEMA, "public_bundle": {"path": str(args.public), "sha256": public_sha},
              "label_sources": {"price": {k: sha256(v) for k, v in tables.items()},
                                "swissprot_experimental": sha256(args.swissprot_ecs)},
              "labels": len(rows), "labelled_proteins": len(label_index(rows)),
              **join(cases, label_index(rows))}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"cases": result["cases"], "pairs": result["pairs"], "labels": result["labels"]}))


if __name__ == "__main__":
    main()
