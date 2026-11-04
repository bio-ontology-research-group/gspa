"""Candidate table conversion and loaders for the second off-growth evaluation (plan §4 step 2).

Provides:
- ``clean_sidecar_to_candidates(sidecar_dir, out_tsv)``: converts per-tag
  CLEAN sidecar output (``*.clean.tsv``, native GMM confidence in ``(0, 1]``)
  into a sorted 3-column candidate TSV (``protein_id, ec, score``) and writes
  a companion ``.receipt.json`` with counts and SHA-256 digests.
- ``load_candidates(tsv)``: loads either a 3-column candidate TSV or a
  10-column six-predictor ``candidate_claims.tsv`` into ``{(protein_id, ec): score}``.
- ``best_by_ec(candidates)``: returns ``{ec: [(protein_id, score), ...]}`` sorted
  by ``(-score, protein_id)``.
- ``best_by_protein(candidates)``: returns ``{protein_id: [(ec, score), ...]}``
  sorted by ``(-score, ec)``.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Tuple, Union

V3_DIR = Path(__file__).resolve().parents[1]
if str(V3_DIR) not in sys.path:
    sys.path.insert(0, str(V3_DIR))

from offgrowth.panel import is_complete_ec, normalize_ec  # noqa: E402


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def _discover_sidecar_tsvs(sidecar_dir: Path) -> List[Path]:
    if sidecar_dir.is_file():
        return [sidecar_dir]
    if not sidecar_dir.exists():
        raise FileNotFoundError(f"CLEAN sidecar path does not exist: {sidecar_dir}")
    direct = sorted(p for p in sidecar_dir.glob("*.clean.tsv") if p.is_file())
    if direct:
        return direct
    nested = sorted(p for p in sidecar_dir.rglob("*.clean.tsv") if p.is_file())
    if nested:
        return nested
    any_tsv = sorted(
        p for p in sidecar_dir.glob("*.tsv")
        if p.is_file() and not p.name.endswith(".receipt.tsv")
    )
    if any_tsv:
        return any_tsv
    raise FileNotFoundError(f"No CLEAN sidecar TSV files found under {sidecar_dir}")


def _validate_score(raw_val: Any, context: str = "") -> float:
    try:
        score = float(raw_val)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid non-numeric score {raw_val!r} {context}".strip()) from exc
    if math.isnan(score) or math.isinf(score) or score < 0.0 or score > 1.0:
        raise ValueError(f"Candidate score {score} out of bounds [0.0, 1.0] {context}".strip())
    return score


def clean_sidecar_to_candidates(
    sidecar_dir: Union[Path, str],
    out_tsv: Union[Path, str],
    receipt_json: Optional[Union[Path, str]] = None,
) -> Dict[str, Any]:
    """Convert CLEAN sidecar output into a sorted ``protein_id, ec, score`` TSV and receipt JSON.

    Enforces ``0 < score <= 1`` (raises ``ValueError`` if ``score < 0`` or
    ``score > 1`` or ``NaN``; skips rows with ``score == 0``), strips ``EC:``
    prefixes, retains complete 4-level EC numbers, deduplicates ``(protein_id, ec)``
    by taking the maximum score, and writes sorted rows to ``out_tsv``.
    """
    src_path = Path(sidecar_dir)
    dst_path = Path(out_tsv)
    rec_path = Path(receipt_json) if receipt_json is not None else dst_path.with_suffix(".receipt.json")

    tsv_files = _discover_sidecar_tsvs(src_path)

    best_scores: Dict[Tuple[str, str], float] = {}
    input_rows = 0
    skipped_zero = 0
    skipped_incomplete_ec = 0
    skipped_non_ec = 0
    duplicate_pairs = 0

    for tsv_file in tsv_files:
        with tsv_file.open("r", encoding="utf-8", newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row_idx, row in enumerate(reader, start=2):
                input_rows += 1
                protein_id = (row.get("protein_id") or "").strip()
                if not protein_id:
                    raise ValueError(f"Missing protein_id at {tsv_file}:{row_idx}")

                ann_type = (row.get("annotation_type") or "").strip().upper()
                if ann_type and ann_type != "EC":
                    skipped_non_ec += 1
                    continue

                raw_term = (
                    row.get("term")
                    if row.get("term") is not None
                    else row.get("ec", "")
                )
                ec = normalize_ec(str(raw_term or ""))
                if not ec or not is_complete_ec(ec):
                    skipped_incomplete_ec += 1
                    continue

                score = _validate_score(
                    row.get("score"),
                    context=f"at {tsv_file}:{row_idx} ({protein_id}, {ec})",
                )
                if score <= 0.0:
                    skipped_zero += 1
                    continue

                key = (protein_id, ec)
                if key in best_scores:
                    duplicate_pairs += 1
                    if score > best_scores[key]:
                        best_scores[key] = score
                else:
                    best_scores[key] = score

    dst_path.parent.mkdir(parents=True, exist_ok=True)
    sorted_items = sorted(best_scores.items(), key=lambda kv: (kv[0][0], kv[0][1]))
    with dst_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["protein_id", "ec", "score"])
        for (protein_id, ec), score in sorted_items:
            writer.writerow([protein_id, ec, f"{score:.6f}"])

    unique_proteins = len({p for (p, _) in best_scores})
    unique_ecs = len({e for (_, e) in best_scores})
    receipt: Dict[str, Any] = {
        "schema": "gspa.offgrowth2.clean_candidates_receipt/1",
        "sidecar_dir": str(src_path),
        "source_files": [
            {"path": str(p), "sha256": _sha256_file(p)} for p in tsv_files
        ],
        "out_tsv": str(dst_path),
        "out_tsv_sha256": _sha256_file(dst_path),
        "receipt_path": str(rec_path),
        "input_rows": input_rows,
        "written_rows": len(sorted_items),
        "n_rows": len(sorted_items),
        "n_candidates": len(sorted_items),
        "n_proteins": unique_proteins,
        "unique_proteins": unique_proteins,
        "n_ecs": unique_ecs,
        "unique_ecs": unique_ecs,
        "skipped_zero_score_rows": skipped_zero,
        "skipped_incomplete_ec_rows": skipped_incomplete_ec,
        "skipped_non_ec_rows": skipped_non_ec,
        "duplicate_pairs_merged": duplicate_pairs,
    }
    rec_path.parent.mkdir(parents=True, exist_ok=True)
    rec_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return receipt


def load_candidates(
    tsv: Union[Path, str, Mapping[Tuple[str, str], float], Iterable[Any]]
) -> Dict[Tuple[str, str], float]:
    """Load candidate table into ``{(protein_id, ec): score}``.

    Supports:
    - 3-column ``protein_id, ec, score`` TSV
    - 10-column six-predictor ``candidate_claims.tsv`` (filters to ``claim_type == 'EC'``
      when ``claim_type`` is present)
    - Sidecar TSV with ``protein_id, term, score``
    - Pre-loaded ``Mapping[(protein_id, ec), float]`` or iterable of ``CandidateClaim`` objects.
    """
    if isinstance(tsv, Mapping):
        out: Dict[Tuple[str, str], float] = {}
        for k, raw_s in tsv.items():
            if isinstance(k, tuple) and len(k) == 2:
                p, raw_e = str(k[0]).strip(), str(k[1]).strip()
            elif isinstance(k, str) and ("\t" in k or "|" in k):
                sep = "\t" if "\t" in k else "|"
                parts = k.split(sep, 1)
                p, raw_e = parts[0].strip(), parts[1].strip()
            elif isinstance(raw_s, Mapping):
                p = str(k).strip()
                for sub_e, sub_s in raw_s.items():
                    ec = normalize_ec(str(sub_e))
                    if not p or not ec or not is_complete_ec(ec):
                        continue
                    s = _validate_score(sub_s, context=f"({p}, {ec})")
                    if s > 0.0 and ((p, ec) not in out or s > out[(p, ec)]):
                        out[(p, ec)] = s
                continue
            else:
                continue
            ec = normalize_ec(raw_e)
            if not p or not ec or not is_complete_ec(ec):
                continue
            s = _validate_score(raw_s, context=f"({p}, {ec})")
            if s > 0.0:
                if (p, ec) not in out or s > out[(p, ec)]:
                    out[(p, ec)] = s
        return out

    if not isinstance(tsv, (Path, str)):
        out = {}
        for item in tsv:
            if isinstance(item, Mapping):
                p = str(item.get("protein_id") or "").strip()
                raw_e = item.get("ec") if item.get("ec") is not None else item.get("term", "")
                raw_s = item.get("score", 0.0)
            elif isinstance(item, (tuple, list)) and len(item) >= 3:
                p = str(item[0]).strip()
                raw_e = str(item[1]).strip()
                raw_s = item[2]
            else:
                p = str(getattr(item, "protein_id", "")).strip()
                raw_e = str(getattr(item, "ec", ""))
                raw_s = getattr(item, "score", 0.0)
            e = normalize_ec(str(raw_e))
            if not p or not e or not is_complete_ec(e):
                continue
            s = _validate_score(raw_s, context=f"({p}, {e})")
            if s > 0.0:
                if (p, e) not in out or s > out[(p, e)]:
                    out[(p, e)] = s
        return out

    path = Path(tsv)
    out = {}
    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row_idx, row in enumerate(reader, start=2):
            claim_type = (row.get("claim_type") or row.get("annotation_type") or "").strip().upper()
            if claim_type and claim_type != "EC":
                continue
            protein_id = (row.get("protein_id") or "").strip()
            raw_ec = row.get("ec") if row.get("ec") is not None else row.get("term", "")
            ec = normalize_ec(str(raw_ec or ""))
            if not protein_id or not ec or not is_complete_ec(ec):
                continue
            score = _validate_score(
                row.get("score"),
                context=f"at {path}:{row_idx} ({protein_id}, {ec})",
            )
            if score <= 0.0:
                continue
            key = (protein_id, ec)
            if key not in out or score > out[key]:
                out[key] = score
    return out


def best_by_ec(
    candidates: Union[Path, str, Mapping[Tuple[str, str], float], Iterable[Any]]
) -> Dict[str, List[Tuple[str, float]]]:
    """Group candidates by EC, returning ``{ec: [(protein_id, score), ...]}`` sorted by ``(-score, protein_id)``."""
    cand_map = load_candidates(candidates) if not isinstance(candidates, dict) else candidates
    by_ec: Dict[str, List[Tuple[str, float]]] = {}
    for (protein_id, ec), score in cand_map.items():
        by_ec.setdefault(ec, []).append((protein_id, float(score)))
    for ec in by_ec:
        by_ec[ec].sort(key=lambda item: (-item[1], item[0]))
    return by_ec


def best_by_protein(
    candidates: Union[Path, str, Mapping[Tuple[str, str], float], Iterable[Any]]
) -> Dict[str, List[Tuple[str, float]]]:
    """Group candidates by protein, returning ``{protein_id: [(ec, score), ...]}`` sorted by ``(-score, ec)``."""
    cand_map = load_candidates(candidates) if not isinstance(candidates, dict) else candidates
    by_prot: Dict[str, List[Tuple[str, float]]] = {}
    for (protein_id, ec), score in cand_map.items():
        by_prot.setdefault(protein_id, []).append((ec, float(score)))
    for prot in by_prot:
        by_prot[prot].sort(key=lambda item: (-item[1], item[0]))
    return by_prot


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="Convert CLEAN sidecar output to candidate table TSV.")
    ap.add_argument("--sidecar-dir", type=Path, required=True, help="Directory or file with *.clean.tsv output")
    ap.add_argument("--out-tsv", type=Path, required=True, help="Output 3-column candidate TSV path")
    ap.add_argument("--receipt-json", type=Path, default=None, help="Optional receipt JSON path")
    args = ap.parse_args(argv)
    receipt = clean_sidecar_to_candidates(args.sidecar_dir, args.out_tsv, args.receipt_json)
    print(json.dumps(receipt, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
