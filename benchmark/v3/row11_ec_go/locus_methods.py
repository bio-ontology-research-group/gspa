#!/usr/bin/env python3
"""Simple and learned locus-level comparators and shared guardrail for Universe B."""
from __future__ import annotations

import argparse
import copy
import math
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "benchmark" / "v3"))
import evaluate_stage3_agent as E  # noqa: E402
import v3_runner as runner  # noqa: E402
from row11_ec_go import guardrail as G  # noqa: E402
from row11_ec_go.build_bundles import LOCUS_PUBLIC_SCHEMA, reject_public_leakage  # noqa: E402
from row11_ec_go.guardrail import EngineConsistentGuardrail, load as load_guardrail, normalize_go_check  # noqa: E402

MODES = E.MODES
PROPOSALS_SCHEMA = "gspa.locus-proposals/1"
GUARDRAIL_POLICY = "shared-reject-locus-v1"


def abstain_locus(loc: dict, **extra) -> dict:
    return {
        "anomaly_id": loc["anomaly_id"],
        "action": "ABSTAIN",
        "protein_id": loc.get("protein_id", ""),
        "function_id": "",
        "confidence": 0.0,
        **extra,
    }


def rank_locus(mode: str, loc: dict) -> dict:
    if mode not in MODES:
        raise ValueError(f"unknown comparator: {mode!r}")
    if mode == "abstain_only":
        return abstain_locus(loc)

    def score(c: dict):
        if mode == "rank_by_raw_seq":
            return c.get("raw_sequence_score")
        if mode == "rank_by_stage1_posterior":
            return c.get("baseline_probability")
        if mode == "rank_by_context_support":
            return float(bool(c.get("in_called_pathway")))
        return (c.get("model_probabilities") or {}).get(mode)

    eligible = [(score(c), c["function_id"]) for c in loc.get("candidates", []) if score(c) is not None]
    if not eligible:
        return abstain_locus(loc, reason="no measured score")
    if any(isinstance(s, bool) or not isinstance(s, (int, float)) or not math.isfinite(s) for s, _ in eligible):
        raise ValueError("nonfinite candidate score")
    score_value, fid = sorted(eligible, key=lambda p: (-float(p[0]), p[1]))[0]
    return {
        "anomaly_id": loc["anomaly_id"],
        "action": "PROMOTE",
        "protein_id": loc["protein_id"],
        "function_id": fid,
        "confidence": float(score_value),
    }


def shared_guardrail_locus(loc: dict, proposal: dict, guardrail=None) -> dict:
    """Enforce shared locus constraints for every method and the agent; reject to ABSTAIN on violation."""
    if proposal.get("anomaly_id") != loc["anomaly_id"]:
        raise ValueError("resolution has wrong case ID")
    r = copy.deepcopy(proposal)
    carry = {k: r[k] for k in ("latency_sec", "tool_calls") if k in r}
    action = r.get("action")
    if action == "ABSTAIN":
        r.setdefault("protein_id", loc.get("protein_id", ""))
        r.setdefault("function_id", "")
        r.setdefault("confidence", 0.0)
        return r

    cand_by_fid = {c["function_id"]: c for c in loc.get("candidates", [])}
    fid = r.get("function_id")
    pid = r.get("protein_id")
    if not pid:
        r["protein_id"] = loc.get("protein_id", "")
        pid = r["protein_id"]

    confidence = r.get("confidence")
    reason = None
    constraint_state = None
    if action not in ("PROMOTE", "DEMOTE"):
        reason = "invalid action"
    elif pid != loc.get("protein_id"):
        reason = "outside locus protein"
    elif fid not in cand_by_fid:
        reason = "outside candidate set"
    elif isinstance(confidence, bool) or not isinstance(confidence, (int, float)) or not math.isfinite(confidence):
        reason = "nonfinite confidence"
    elif action == "PROMOTE" and guardrail is not None:
        mapped_go = cand_by_fid[fid].get("mapped_go") or ""
        taxon = loc.get("organism_taxon", "")
        if mapped_go.startswith("GO:") and hasattr(guardrail, "check_go_term"):
            chk = normalize_go_check(guardrail.check_go_term(mapped_go, taxon))
            constraint_state = chk.get("constraint_state")
            if not chk.get("valid", False):
                reason = "symbolic constraint"

    if constraint_state is not None:
        r["constraint_state"] = constraint_state
    if reason:
        extra = dict(carry)
        if constraint_state is not None:
            extra["constraint_state"] = constraint_state
        return abstain_locus(loc, guardrail_rejected=reason, proposal=r, **extra)
    r["confidence"] = float(confidence)
    return r


def public_locus_cases(public: dict) -> list[dict]:
    reject_public_leakage(public)
    if public.get("schema") != LOCUS_PUBLIC_SCHEMA:
        raise ValueError(f"expected public bundle {LOCUS_PUBLIC_SCHEMA!r}")
    cases = public["cases"]
    ids = [c["anomaly_id"] for c in cases]
    if len(set(ids)) != len(ids):
        raise ValueError("duplicate case ID")
    for c in cases:
        if c["locus"].get("anomaly_id") != c["anomaly_id"]:
            raise ValueError("mismatched inner/outer anomaly_id")
    return cases


def suggest_locus(public: dict, mode: str, guardrail=None, solve=None) -> dict:
    cases = public_locus_cases(public)
    start = time.monotonic()
    rows = []
    for case in cases:
        loc = copy.deepcopy(case["locus"])
        proposal = solve(loc) if solve else rank_locus(mode, loc)
        rows.append(shared_guardrail_locus(loc, proposal, guardrail))
    return {
        "schema": PROPOSALS_SCHEMA,
        "public_sha256": E.digest(public),
        "mode": mode,
        "guardrail_policy": GUARDRAIL_POLICY,
        "resolutions": rows,
        "decision_wall_seconds": time.monotonic() - start,
        "shared_production_wall_seconds": public["shared_production_wall_seconds"],
    }


def main(argv=None) -> int:
    argv_list = list(sys.argv[1:] if argv is None else argv)
    if argv_list and argv_list[0] == "suggest":
        argv_list = argv_list[1:]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--public", required=True, type=Path)
    parser.add_argument("--mode", required=True, choices=MODES)
    parser.add_argument("--go-obo", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args(argv_list)

    code_paths = (
        Path(__file__).resolve(),
        Path(G.__file__).resolve(),
        Path(E.__file__).resolve(),
        Path(E.agent.__file__).resolve(),
    )
    code_identity = {str(p.relative_to(ROOT)): runner.sha256_file(str(p)) for p in code_paths}
    resources = ROOT / "gspa-core" / "src" / "main" / "resources"
    guard = load_guardrail(resources, args.go_obo)

    result = suggest_locus(E.read(args.public), args.mode, guard)
    result["decision_code_sha256"] = code_identity
    result["constraint_inputs"] = [runner.hash_record(str(args.go_obo))] + [
        runner.hash_record(str(p))
        for directory in ("taxon-constraints", "anatomical-constraints")
        for p in sorted((resources / directory).glob("*.tsv"))
    ]
    if code_identity != {str(p.relative_to(ROOT)): runner.sha256_file(str(p)) for p in code_paths}:
        raise RuntimeError("decision code changed during execution")
    E.write(args.out, result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
