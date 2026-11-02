#!/usr/bin/env python3
"""Where the experimental EC labels sit in a production integration.

For each labelled (genome, protein, EC) of the Price release-candidate genomes,
report the integrated posterior that the production run assigned to that EC on
that protein: ``called`` (at or above 0.5), ``sub_threshold`` (between the
anomaly extractor's candidate floor 0.05 and 0.5), ``below_floor`` or ``absent``.
This measures how many experimentally labelled functions could in principle
have appeared as unresolved candidates had the extractor keyed on EC terms; it
is read only after the integration finished and it never feeds back into it.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

CANDIDATE_FLOOR = 0.05
ANNOTATED = 0.5


def labels(directory: Path, swissprot: Path | None) -> dict:
    out: dict = {}
    for path in sorted(directory.glob("*.ec.tsv")):
        genome_id = path.name.split(".ec.tsv")[0]
        with path.open(newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                out.setdefault(genome_id, {}).setdefault(row["protein_id"], {})[row["ec"]] = "price"
    if swissprot is not None:
        with swissprot.open(newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                ec = row["function_id"].removeprefix("EC:")
                out.setdefault(row["genome_id"], {}).setdefault(row["protein_id"], {}).setdefault(ec, "swissprot")
    return out


def posteriors(path: Path) -> dict:
    out = {}
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        for row in reader:
            if len(row) < 5 or row[1] != "EC":
                continue
            key = (row[0], row[2].removeprefix("EC:"))
            out[key] = max(out.get(key, 0.0), float(row[4]))
    return out


def status(value):
    if value is None:
        return "absent"
    if value >= ANNOTATED:
        return "called"
    if value >= CANDIDATE_FLOOR:
        return "sub_threshold"
    return "below_floor"


def run(arm_dir: Path, label_dir: Path, swissprot: Path | None):
    table = labels(label_dir, swissprot)
    counts, rows = {}, []
    for genome_id, proteins in sorted(table.items()):
        path = arm_dir / f"{genome_id}.integrated.tsv"
        if not path.exists():
            continue
        post = posteriors(path)
        for protein_id, ecs in sorted(proteins.items()):
            for ec, source in sorted(ecs.items()):
                value = post.get((protein_id, ec))
                kind = status(value)
                counts.setdefault(source, {}).setdefault(kind, 0)
                counts[source][kind] += 1
                rows.append({"genome_id": genome_id, "protein_id": protein_id, "ec": ec, "source": source,
                             "posterior": value, "status": kind})
    return {"schema": "gspa.price-label-candidate-status/1", "arm": str(arm_dir),
            "floor": CANDIDATE_FLOOR, "annotated": ANNOTATED, "counts": counts, "rows": rows}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--arm", required=True, type=Path)
    parser.add_argument("--price-labels", required=True, type=Path)
    parser.add_argument("--swissprot-ecs", type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args(argv)
    result = run(args.arm, args.price_labels, args.swissprot_ecs)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result["counts"], sort_keys=True))


if __name__ == "__main__":
    main()
