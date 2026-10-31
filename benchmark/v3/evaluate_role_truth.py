#!/usr/bin/env python3
"""Evaluate integration scores on reviewed role-truth loci using only the stdlib.

The metric is **within-locus ranking only** (plus threshold calls), never a
genome-wide accuracy or a sufficient null-model evaluation. Missing actives have
no rank; strict ranking counts tied negatives as above the active. Unknown-truth
members are reported separately and never enter ranking or threshold-call counts.
No model is fitted and no bootstrap or network access is used.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from pathlib import Path

SCHEMA = "gspa.role-truth-evaluation/1"
FORMATS = ("integrated_tsv", "training_table_jsonl", "claims_jsonl")
STATUSES = ("accepted", "pending_review", "rejected")


class EvaluationError(ValueError):
    """Scores or registry records cannot support a deterministic evaluation."""


def sha256(path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_jsonl(path):
    rows = []
    with Path(path).open() as fh:
        for line_number, line in enumerate(fh, 1):
            if not line.strip():
                continue
            try:
                row = json.loads(line)
            except json.JSONDecodeError as exc:
                raise EvaluationError(f"{path}:{line_number}: invalid JSON: {exc}") from exc
            if not isinstance(row, dict):
                raise EvaluationError(f"{path}:{line_number}: expected an object")
            rows.append(row)
    return rows


def _score(value, where):
    if value is None or value == "" or isinstance(value, bool):
        raise EvaluationError(f"{where}: missing or non-numeric score")
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise EvaluationError(f"{where}: non-numeric score {value!r}") from exc
    if not math.isfinite(result):
        raise EvaluationError(f"{where}: score must be finite")
    return result


def parse_genome_map(text):
    """Parse ``score_genome=registry_genome`` comma-separated mappings."""
    result = {}
    if not text:
        return result
    for item in text.split(","):
        if item.count("=") != 1:
            raise EvaluationError("--genome-map entries must have the form a=b")
        source, target = (part.strip() for part in item.split("="))
        if not source or not target or source in result:
            raise EvaluationError("--genome-map has an empty or duplicate source")
        result[source] = target
    return result


def _put(scores, key, value):
    # All formats use max for repeated claims/predictions of one candidate.
    scores[key] = max(value, scores.get(key, -math.inf))


def load_scores(path, fmt, genome_map=None, score_field="score"):
    """Return a max-aggregated score map keyed to the requested input format."""
    genome_map = genome_map or {}
    scores = {}
    if fmt == "integrated_tsv":
        with Path(path).open(newline="") as fh:
            for line_number, row in enumerate(csv.reader(fh, delimiter="\t"), 1):
                if not row:
                    continue
                if len(row) < 5:
                    raise EvaluationError(f"{path}:{line_number}: integrated_tsv needs at least five columns")
                if line_number == 1 and row[0] == "protein_id":
                    continue  # production `gspa integrate` output carries a header row
                protein, kind, function_id = row[0], row[1], row[2]
                value = _score(row[4], f"{path}:{line_number}")
                _put(scores, (None, protein, kind, function_id, None), value)
        return scores
    for line_number, row in enumerate(read_jsonl(path), 1):
        where = f"{path}:{line_number}"
        try:
            protein, kind, function_id = row["protein_id"], row["function_type"], row["function_id"]
        except KeyError as exc:
            raise EvaluationError(f"{where}: missing score identity field {exc.args[0]}") from exc
        if fmt == "claims_jsonl":
            value = _score(row.get("raw_score"), where)
            _put(scores, (None, protein, kind, function_id, None), value)
        elif fmt == "training_table_jsonl":
            raw = row.get(score_field)
            if raw is None and score_field == "score":
                raw = row.get("probability")
            value = _score(raw, where)
            genome = genome_map.get(row.get("genome_id"), row.get("genome_id"))
            _put(scores, (genome, protein, kind, function_id, row.get("role")), value)
        else:
            raise EvaluationError(f"unknown score format {fmt!r}")
    return scores


def evaluate_case(case, scores, fmt, threshold):
    """Score one case; the best-scoring active represents a multi-active locus."""
    member_scores, unknown_scores, called = {}, {}, {}
    for member in case["members"]:
        genome = case["genome_id"] if fmt == "training_table_jsonl" else None
        role = case["role"] if fmt == "training_table_jsonl" else None
        key = (genome, member["protein_id"], case["function_type"], case["function_id"], role)
        value = scores.get(key)
        if member["label"] is None:
            unknown_scores[member["protein_id"]] = value
            continue
        member_scores[member["protein_id"]] = value
        called[member["protein_id"]] = value is not None and value >= threshold
    actives = [m for m in case["members"] if m["label"] == 1]
    negatives = [m for m in case["members"] if m["label"] == 0]
    present_actives = [member_scores[m["protein_id"]] for m in actives if member_scores[m["protein_id"]] is not None]
    present_negatives = [member_scores[m["protein_id"]] for m in negatives if member_scores[m["protein_id"]] is not None]
    present_unknown = [value for value in unknown_scores.values() if value is not None]
    best_active = max(present_actives) if present_actives else None
    scored = best_active is not None and bool(present_negatives)
    rank = None
    negatives_above = len(present_negatives)
    negatives_tied = 0
    unknown_above = len(present_unknown)
    margin = None
    if best_active is not None:
        negatives_above = sum(value > best_active for value in present_negatives)
        negatives_tied = sum(value == best_active for value in present_negatives)
        unknown_above = sum(value > best_active for value in present_unknown)
        # Strict competition rank: a tied negative is not below the active.
        rank = 1 + negatives_above + negatives_tied
        if present_negatives:
            margin = best_active - max(present_negatives)
    false_positives = [m["protein_id"] for m in negatives if called[m["protein_id"]]]
    top1_lenient = best_active is not None and negatives_above == 0
    top1_strict = top1_lenient and negatives_tied == 0
    return {
        "case_id": case["case_id"], "status": case["status"], "genome_id": case["genome_id"],
        "organism": case["organism"], "function_type": case["function_type"],
        "function_id": case["function_id"], "role": case["role"],
        "active_protein_ids": [m["protein_id"] for m in actives], "scores": member_scores,
        "unknown_scores": unknown_scores, "scored": scored, "active_rank": rank,
        "top1_correct_strict": top1_strict, "top1_correct_lenient": top1_lenient,
        "n_negatives_above_active": negatives_above, "n_negatives_tied_with_active": negatives_tied,
        "n_unknown_above_active": unknown_above, "margin": margin,
        "called_at_threshold": called, "false_positive_negatives_at_threshold": false_positives,
        "n_false_positive_negatives_at_threshold": len(false_positives),
        "active_called_at_threshold": bool(present_actives) and all(called[m["protein_id"]] for m in actives),
    }


def aggregate(rows):
    """Aggregate ranking over fully scored cases and calls over every member."""
    scored = [r for r in rows if r["scored"]]
    n_scored = len(scored)
    negatives_total = actives_total = 0
    negatives_called = actives_called = 0
    for row in rows:
        for protein in row["active_protein_ids"]:
            actives_total += 1
            actives_called += bool(row["called_at_threshold"][protein])
        active_ids = set(row["active_protein_ids"])
        for protein in row["scores"]:
            if protein not in active_ids:
                negatives_total += 1
                negatives_called += bool(row["called_at_threshold"][protein])
    return {
        "n_cases": len(rows), "n_cases_scored": n_scored,
        "top1_strict": (sum(r["top1_correct_strict"] for r in scored) / n_scored) if n_scored else None,
        "top1_lenient": (sum(r["top1_correct_lenient"] for r in scored) / n_scored) if n_scored else None,
        "mrr": (sum(1.0 / r["active_rank"] for r in scored) / n_scored) if n_scored else None,
        "mean_margin": (sum(r["margin"] for r in scored) / n_scored) if n_scored else None,
        "negatives_called_at_threshold": negatives_called, "negatives_total": negatives_total,
        "actives_called_at_threshold": actives_called, "actives_total": actives_total,
        "negative_call_rate": (negatives_called / negatives_total) if negatives_total else None,
        "active_call_rate": (actives_called / actives_total) if actives_total else None,
    }


def _markdown(summary):
    lines = [
        "# Role-truth within-locus evaluation", "",
        "Reviewed within-locus ranking only; this is not genome-wide performance or a null-model fit.", "",
        "| Scope | Cases | Scored | Top-1 strict | Top-1 lenient | MRR | Mean margin | Negatives called | Actives called |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]

    def fmt(value):
        return "NA" if value is None else f"{value:.6g}"

    def row(name, metrics):
        return (f"| {name} | {metrics['n_cases']} | {metrics['n_cases_scored']} | "
                f"{fmt(metrics['top1_strict'])} | {fmt(metrics['top1_lenient'])} | "
                f"{fmt(metrics['mrr'])} | {fmt(metrics['mean_margin'])} | "
                f"{metrics['negatives_called_at_threshold']}/{metrics['negatives_total']} | "
                f"{metrics['actives_called_at_threshold']}/{metrics['actives_total']} |")

    lines.append(row("Overall", summary))
    for kind, metrics in summary["by_function_type"].items():
        lines.append(row("Function: " + kind, metrics))
    for organism, metrics in summary["by_organism"].items():
        lines.append(row("Organism: " + organism, metrics))
    lines.extend(["", "## Negative-only observations", "",
                  "Usable for false-positive counting only; never treated as locus groups. "
                  "Unknown-truth members are excluded.", "",
                  f"Negatives called at threshold: {summary['negative_only_called_at_threshold']}/"
                  f"{summary['negative_only_total']}"])
    return "\n".join(lines) + "\n"


def evaluate(registry, scores_file, fmt, out, statuses=("accepted", "pending_review"), threshold=0.5,
             genome_map=None, score_field="score", negative_only=None, cli_args=None):
    """Evaluate selected cases and write per_case.jsonl, summary.json, summary.md."""
    if fmt not in FORMATS:
        raise EvaluationError(f"format must be one of {FORMATS}")
    if not math.isfinite(threshold):
        raise EvaluationError("threshold must be finite")
    cases = [case for case in read_jsonl(registry) if case.get("status") in set(statuses)]
    score_map = load_scores(scores_file, fmt, genome_map=genome_map, score_field=score_field)
    rows = [evaluate_case(case, score_map, fmt, threshold) for case in cases]
    negative_only_rows = []
    if negative_only is not None:
        negative_only_cases = [case for case in read_jsonl(negative_only) if case.get("status") in set(statuses)]
        negative_only_rows = [evaluate_case(case, score_map, fmt, threshold) for case in negative_only_cases]
    summary = aggregate(rows)
    summary.update({
        "schema": SCHEMA, "registry_sha256": sha256(registry), "scores_sha256": sha256(scores_file),
        "negative_only_sha256": sha256(negative_only) if negative_only is not None else None,
        "negative_only_cases": len(negative_only_rows),
        "negative_only_called_at_threshold": sum(r["n_false_positive_negatives_at_threshold"]
                                                  for r in negative_only_rows),
        "negative_only_total": sum(sum(m["label"] == 0 for m in case["members"])
                                    for case in (negative_only_cases if negative_only is not None else [])),
        "by_function_type": {kind: aggregate([r for r in rows if r["function_type"] == kind])
                             for kind in sorted({r["function_type"] for r in rows})},
        "by_organism": {organism: aggregate([r for r in rows if r["organism"] == organism])
                        for organism in sorted({r["organism"] for r in rows})},
        "args": cli_args if cli_args is not None else {
            "registry": str(registry), "scores": str(scores_file), "format": fmt,
            "statuses": list(statuses), "threshold": threshold, "genome_map": genome_map or {},
            "score_field": score_field, "negative_only": str(negative_only) if negative_only is not None else None,
            "out": str(out),
        },
    })
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    with (out / "per_case.jsonl").open("w") as fh:
        for row in rows:
            fh.write(json.dumps(row, sort_keys=True, allow_nan=False) + "\n")
    if negative_only is not None:
        with (out / "negative_only_per_case.jsonl").open("w") as fh:
            for row in negative_only_rows:
                fh.write(json.dumps(row, sort_keys=True, allow_nan=False) + "\n")
    (out / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n")
    (out / "summary.md").write_text(_markdown(summary))
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--registry", type=Path, required=True)
    parser.add_argument("--scores", type=Path, required=True)
    parser.add_argument("--format", choices=FORMATS, required=True)
    parser.add_argument("--statuses", default="accepted,pending_review")
    parser.add_argument("--threshold", type=float, default=0.5)
    parser.add_argument("--genome-map", default="")
    parser.add_argument("--score-field", default="score")
    parser.add_argument("--negative-only", type=Path,
                        help="negative-only observations for separate false-positive counting")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        statuses = tuple(value.strip() for value in args.statuses.split(",") if value.strip())
        if not statuses or set(statuses) - set(STATUSES):
            raise EvaluationError("--statuses contains an unknown or empty status set")
        genome_map = parse_genome_map(args.genome_map)
        cli_args = {"registry": str(args.registry), "scores": str(args.scores), "format": args.format,
                    "statuses": list(statuses), "threshold": args.threshold, "genome_map_text": args.genome_map,
                    "genome_map": genome_map, "score_field": args.score_field,
                    "negative_only": str(args.negative_only) if args.negative_only is not None else None,
                    "out": str(args.out)}
        summary = evaluate(args.registry, args.scores, args.format, args.out, statuses=statuses,
                           threshold=args.threshold, genome_map=genome_map, score_field=args.score_field,
                           negative_only=args.negative_only, cli_args=cli_args)
    except (OSError, EvaluationError) as exc:
        parser.exit(2, f"error: {exc}\n")
    print(json.dumps(summary, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
