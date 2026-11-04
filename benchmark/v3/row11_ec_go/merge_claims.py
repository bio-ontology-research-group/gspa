#!/usr/bin/env python3
"""Stream and merge per-predictor Price claims with a uniform raw-score floor.

Reads each predictor's ``claims.jsonl`` once in binary mode, computing its SHA-256
digest while streaming and verifying it against ``part.json``. Keeps rows with
``raw_score >= min_raw`` (rows without a numeric ``raw_score`` are kept and counted
separately) in predictor order then file order without re-sorting, optionally
appending an un-filtered extra claims file (such as the Swiss-Prot DIAMOND
contract claims) after verifying its SHA-256.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
from typing import Optional

DEFAULT_PREDICTORS = ("clean", "dpz", "enzbert", "graphec", "mapred", "topec")
_RAW_KEY = b'"raw_score"'
_WS = b" \t"
_DELIMS = b",}\r\n \t"
_NON_NUMERIC_PREFIX = (b'"', b"{", b"[", b"t", b"f", b"n", b"N", b"I", b"-I")


def _parse_raw_score_fast(raw_line: bytes) -> tuple[bool, Optional[float]]:
    """Fast substring parser for top-level numeric ``raw_score``.

    Returns ``(True, score)`` when a single unambiguous finite numeric
    ``"raw_score": <num>`` field is found; otherwise ``(False, None)`` so the
    caller falls back to ``json.loads``.
    """
    idx = raw_line.find(_RAW_KEY)
    if idx == -1 or raw_line.find(_RAW_KEY, idx + len(_RAW_KEY)) != -1:
        return False, None
    n = len(raw_line)
    pos = idx + len(_RAW_KEY)
    while pos < n and raw_line[pos : pos + 1] in _WS:
        pos += 1
    if pos >= n or raw_line[pos : pos + 1] != b":":
        return False, None
    pos += 1
    while pos < n and raw_line[pos : pos + 1] in _WS:
        pos += 1
    end = pos
    while end < n and raw_line[end : end + 1] not in _DELIMS:
        end += 1
    token = raw_line[pos:end]
    if not token or token.startswith(_NON_NUMERIC_PREFIX):
        return False, None
    trail = end
    while trail < n and raw_line[trail : trail + 1] in _WS:
        trail += 1
    if trail >= n or raw_line[trail : trail + 1] not in (b",", b"}"):
        return False, None
    try:
        val = float(token)
    except ValueError:
        return False, None
    if not math.isfinite(val):
        return False, None
    return True, val


def _should_keep_line(raw_line: bytes, min_raw: float) -> tuple[bool, bool]:
    """Return ``(keep, is_no_raw_score)`` for one non-empty JSONL line."""
    ok, val = _parse_raw_score_fast(raw_line)
    if ok and val is not None:
        return (val >= min_raw), False
    rec = json.loads(raw_line.decode("utf-8"))
    if not isinstance(rec, dict):
        raise ValueError("claim record is not a JSON object")
    raw_val = rec.get("raw_score")
    if isinstance(raw_val, (int, float)) and not isinstance(raw_val, bool) and math.isfinite(raw_val):
        return (float(raw_val) >= min_raw), False
    return True, True


def _normalize_predictors(values: list[str]) -> list[str]:
    out: list[str] = []
    for item in values:
        for tok in item.replace(",", " ").split():
            if tok:
                out.append(tok)
    if not out:
        raise ValueError("at least one predictor is required")
    return out


def merge_claims(
    parts_root: Path,
    predictors: list[str],
    assembly: str,
    min_raw: float,
    out_dir: Path,
    extra_claims: Optional[Path] = None,
    extra_sha256: Optional[str] = None,
    extra_name: str = "swissprot_diamond_baseline",
) -> dict:
    """Stream predictor claims and optional extra claims into ``out_dir``."""
    parts_root = Path(parts_root)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    if extra_claims is not None and not extra_sha256:
        raise ValueError("--extra-sha256 is required when --extra-claims is set")

    tmp_out = out_dir / "claims.jsonl.tmp"
    final_out = out_dir / "claims.jsonl"
    part_out = out_dir / "part.json"

    out_hasher = hashlib.sha256()
    total_claims = 0
    total_bytes = 0
    pred_stats: dict[str, dict] = {}
    extra_record: Optional[dict] = None

    try:
        with tmp_out.open("wb") as out_fh:
            for pred in predictors:
                pdir = parts_root / pred / assembly
                part_meta_path = pdir / "part.json"
                if not part_meta_path.is_file():
                    raise FileNotFoundError(f"missing part.json: {part_meta_path}")
                part_meta = json.loads(part_meta_path.read_text(encoding="utf-8"))
                if part_meta.get("assembly") and part_meta["assembly"] != assembly:
                    raise ValueError(
                        f"{part_meta_path}: assembly {part_meta.get('assembly')!r} != {assembly!r}"
                    )
                claims_rel = part_meta["claims"].get("path", "claims.jsonl")
                expected_sha = part_meta["claims"]["sha256"].strip().lower()
                claims_path = pdir / claims_rel
                if not claims_path.is_file():
                    raise FileNotFoundError(f"missing claims file: {claims_path}")

                in_hasher = hashlib.sha256()
                claims_in = 0
                claims_kept = 0
                no_raw_score = 0

                with claims_path.open("rb") as in_fh:
                    for raw_line in in_fh:
                        in_hasher.update(raw_line)
                        if not raw_line.strip():
                            continue
                        claims_in += 1
                        keep, is_missing = _should_keep_line(raw_line, min_raw)
                        if is_missing:
                            no_raw_score += 1
                        if keep:
                            claims_kept += 1
                            out_line = raw_line if raw_line.endswith(b"\n") else (raw_line + b"\n")
                            out_fh.write(out_line)
                            out_hasher.update(out_line)
                            total_bytes += len(out_line)
                            total_claims += 1

                actual_sha = in_hasher.hexdigest()
                if actual_sha != expected_sha:
                    raise ValueError(
                        f"{claims_path}: sha256 mismatch (got {actual_sha}, expected {expected_sha})"
                    )
                pred_stats[pred] = {
                    "claims_in": claims_in,
                    "claims_kept": claims_kept,
                    "no_raw_score": no_raw_score,
                    "sha256_in": actual_sha,
                    "verified": True,
                }

            if extra_claims is not None:
                extra_path = Path(extra_claims)
                if not extra_path.is_file():
                    raise FileNotFoundError(f"missing extra claims file: {extra_path}")
                expected_extra_sha = extra_sha256.strip().lower()
                extra_hasher = hashlib.sha256()
                extra_count = 0
                with extra_path.open("rb") as ex_fh:
                    for raw_line in ex_fh:
                        extra_hasher.update(raw_line)
                        if not raw_line.strip():
                            continue
                        extra_count += 1
                        out_line = raw_line if raw_line.endswith(b"\n") else (raw_line + b"\n")
                        out_fh.write(out_line)
                        out_hasher.update(out_line)
                        total_bytes += len(out_line)
                        total_claims += 1
                actual_extra_sha = extra_hasher.hexdigest()
                if actual_extra_sha != expected_extra_sha:
                    raise ValueError(
                        f"{extra_path}: extra sha256 mismatch (got {actual_extra_sha}, expected {expected_extra_sha})"
                    )
                extra_record = {
                    "name": extra_name,
                    "path": str(extra_path),
                    "sha256": actual_extra_sha,
                    "claims": extra_count,
                }

        os.replace(tmp_out, final_out)
    except Exception:
        if tmp_out.exists():
            tmp_out.unlink()
        raise

    record = {
        "assembly": assembly,
        "min_raw": min_raw,
        "predictors": pred_stats,
        "extra": extra_record,
        "claims": {
            "path": "claims.jsonl",
            "sha256": out_hasher.hexdigest(),
            "claims": total_claims,
            "bytes": total_bytes,
        },
    }
    part_out.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return record


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parts-root", required=True, type=Path)
    parser.add_argument("--predictors", nargs="+", default=list(DEFAULT_PREDICTORS))
    parser.add_argument("--assembly", required=True)
    parser.add_argument("--min-raw", type=float, default=0.05)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--extra-claims", type=Path, default=None)
    parser.add_argument("--extra-sha256", default=None)
    parser.add_argument("--extra-name", default="swissprot_diamond_baseline")
    args = parser.parse_args(argv)
    predictors = _normalize_predictors(args.predictors)
    record = merge_claims(
        parts_root=args.parts_root,
        predictors=predictors,
        assembly=args.assembly,
        min_raw=args.min_raw,
        out_dir=args.out,
        extra_claims=args.extra_claims,
        extra_sha256=args.extra_sha256,
        extra_name=args.extra_name,
    )
    print(
        json.dumps(
            {
                "assembly": record["assembly"],
                "claims": record["claims"]["claims"],
                "sha256": record["claims"]["sha256"],
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    os.environ.setdefault("TMPDIR", str(Path.home() / "tmp"))
    raise SystemExit(main())
