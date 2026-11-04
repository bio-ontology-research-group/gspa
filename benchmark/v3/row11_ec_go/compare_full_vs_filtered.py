#!/usr/bin/env python3
"""Compare full vs. raw-score-filtered stage-1 integrations on kept (protein, EC) pairs.

Streams ``--claims-full`` using the fast line-oriented raw-score check so only
rows with ``raw_score >= min_raw`` (or without a numeric ``raw_score``) are parsed
as JSON, then compares their posteriors between ``--full`` and ``--filtered``
integrated TSVs.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "benchmark" / "v3"))
from row11_ec_go.merge_claims import _should_keep_line  # noqa: E402

Pair = tuple[str, str]
CALLED_THRESHOLD = 0.5


def normalize_ec(value: str) -> str:
    v = value.strip()
    if v.upper().startswith("EC:"):
        v = v[3:]
    return v


def parse_ec_integrated(path: Path) -> dict[Pair, float]:
    """Read an integrated TSV and return ``(protein_id, bare_ec) -> posterior_prob``."""
    out: dict[Pair, float] = {}
    with Path(path).open(encoding="utf-8") as fh:
        header = next(fh, "").rstrip("\r\n").split("\t")
        if "protein_id" not in header or "function_id" not in header or "posterior_prob" not in header:
            raise ValueError(f"{path}: header lacks protein_id/function_id/posterior_prob")
        ip = header.index("protein_id")
        if_ = header.index("function_id")
        ipr = header.index("posterior_prob")
        it = header.index("type") if "type" in header else None
        for line_no, ln in enumerate(fh, 2):
            if not ln.strip():
                continue
            fields = ln.rstrip("\r\n").split("\t")
            need = max(ip, if_, ipr, it if it is not None else 0)
            if len(fields) <= need:
                raise ValueError(f"{path}:{line_no}: truncated row")
            if it is not None and fields[it].strip().upper() != "EC":
                continue
            score = float(fields[ipr])
            if not (0.0 <= score <= 1.0) or not math.isfinite(score):
                raise ValueError(f"{path}:{line_no}: invalid posterior_prob {score!r}")
            key = (fields[ip].strip(), normalize_ec(fields[if_]))
            if score > out.get(key, -1.0):
                out[key] = score
    return out


def read_kept_pairs(claims_full: Path, min_raw: float) -> set[Pair]:
    """Return ``(protein_id, bare_ec)`` pairs with ``raw_score >= min_raw`` in ``claims_full``."""
    kept: set[Pair] = set()
    with Path(claims_full).open("rb") as fh:
        for raw_line in fh:
            if not raw_line.strip():
                continue
            keep, _ = _should_keep_line(raw_line, min_raw)
            if not keep:
                continue
            rec = json.loads(raw_line.decode("utf-8"))
            ftype = str(rec.get("function_type", "EC")).strip().upper()
            if ftype != "EC":
                continue
            pid = str(rec["protein_id"]).strip()
            ec = normalize_ec(str(rec["function_id"]))
            kept.add((pid, ec))
    return kept


def compare_full_vs_filtered(
    full_tsv: Path,
    filtered_tsv: Path,
    min_raw: float,
    claims_full: Path,
    out_json: Path,
) -> dict:
    """Compare posteriors on kept ``(protein, EC)`` pairs and write ``out_json``."""
    kept_pairs = read_kept_pairs(Path(claims_full), min_raw)
    full_scores = parse_ec_integrated(Path(full_tsv))
    filt_scores = parse_ec_integrated(Path(filtered_tsv))

    max_abs_diff = 0.0
    count_diff_gt_1e_6 = 0
    only_in_full = 0
    only_in_filtered = 0
    in_neither = 0
    compared_pairs = 0
    kept_called_full = 0
    kept_called_filtered = 0
    called_status_changes = 0

    for pair in sorted(kept_pairs):
        in_f = pair in full_scores
        in_r = pair in filt_scores
        if in_f and full_scores[pair] >= CALLED_THRESHOLD:
            kept_called_full += 1
        if in_r and filt_scores[pair] >= CALLED_THRESHOLD:
            kept_called_filtered += 1
        if in_f and in_r:
            compared_pairs += 1
            pf = full_scores[pair]
            pr = filt_scores[pair]
            diff = abs(pf - pr)
            if diff > max_abs_diff:
                max_abs_diff = diff
            if diff > 1e-6:
                count_diff_gt_1e_6 += 1
            if (pf >= CALLED_THRESHOLD) != (pr >= CALLED_THRESHOLD):
                called_status_changes += 1
        elif in_f:
            only_in_full += 1
            if full_scores[pair] >= CALLED_THRESHOLD:
                called_status_changes += 1
        elif in_r:
            only_in_filtered += 1
            if filt_scores[pair] >= CALLED_THRESHOLD:
                called_status_changes += 1
        else:
            in_neither += 1

    called_full = sum(1 for v in full_scores.values() if v >= CALLED_THRESHOLD)
    called_filtered = sum(1 for v in filt_scores.values() if v >= CALLED_THRESHOLD)

    report = {
        "min_raw": min_raw,
        "kept_pairs": len(kept_pairs),
        "full_rows": len(full_scores),
        "filtered_rows": len(filt_scores),
        "compared_pairs": compared_pairs,
        "max_abs_diff": max_abs_diff,
        "count_diff_gt_1e_6": count_diff_gt_1e_6,
        "differing_gt_1e_6": count_diff_gt_1e_6,
        "only_in_full": only_in_full,
        "only_in_filtered": only_in_filtered,
        "rows_only_in_one": only_in_full + only_in_filtered,
        "in_neither": in_neither,
        "called_full": called_full,
        "called_filtered": called_filtered,
        "kept_called_full": kept_called_full,
        "kept_called_filtered": kept_called_filtered,
        "called_status_changes": called_status_changes,
        "status_unchanged": called_status_changes == 0 and only_in_full == 0 and only_in_filtered == 0,
    }
    out_path = Path(out_json)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return report


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--full", required=True, type=Path)
    parser.add_argument("--filtered", required=True, type=Path)
    parser.add_argument("--min-raw", type=float, default=0.05)
    parser.add_argument("--claims-full", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args(argv)
    report = compare_full_vs_filtered(
        full_tsv=args.full,
        filtered_tsv=args.filtered,
        min_raw=args.min_raw,
        claims_full=args.claims_full,
        out_json=args.out,
    )
    print(json.dumps(report, sort_keys=True))
    return 0


if __name__ == "__main__":
    os.environ.setdefault("TMPDIR", str(Path.home() / "tmp"))
    raise SystemExit(main())
