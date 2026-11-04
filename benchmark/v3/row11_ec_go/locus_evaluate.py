#!/usr/bin/env python3
"""Calibration, evaluation, leak check, and paired bootstrap comparison for Universe B locus cases."""
from __future__ import annotations

import argparse
import math
import random
import sys
from collections import Counter, defaultdict
from pathlib import Path
from scipy import stats as scipy_stats

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "benchmark" / "v3"))
import evaluate_stage3_agent as E  # noqa: E402
from run_planner_benchmark import wilson  # noqa: E402
from row11_ec_go.locus_methods import (  # noqa: E402
    GUARDRAIL_POLICY,
    PROPOSALS_SCHEMA,
    abstain_locus,
    public_locus_cases,
)

POLICY_SCHEMA = "gspa.locus-policy/1"
EVALUATION_SCHEMA = "gspa.locus-evaluation/1"
LEAK_SCHEMA = "gspa.locus-leak-check/1"
COMPARE_SCHEMA = "gspa.locus-compare/1"

PRECISION_NOTE = (
    "Known decisions = labelled cases with a non-ABSTAIN gated action; "
    "precision over known decisions = correct / (correct + harmful + unconfirmed). "
    "Unconfirmed counts against precision because the plan's 'unconfirmed never negative' rule "
    "forbids calling it harm, but a policy that mostly names unconfirmed ECs has not demonstrated "
    "precision. Decoy calls at the threshold are reported as unknown_call_rate and excluded from precision."
)


def checked_locus_predictions(public: dict, proposals: dict) -> dict[str, dict]:
    cases = public_locus_cases(public)
    if (
        proposals.get("schema") != PROPOSALS_SCHEMA
        or proposals.get("public_sha256") != E.digest(public)
        or proposals.get("guardrail_policy") != GUARDRAIL_POLICY
    ):
        raise ValueError("proposals do not match Universe B public bundle or constraint policy")
    rows = proposals["resolutions"]
    by_id = {r["anomaly_id"]: r for r in rows}
    if len(by_id) != len(rows) or set(by_id) != {c["anomaly_id"] for c in cases}:
        raise ValueError("missing, duplicate or extra locus resolutions; failures must be explicit")
    return by_id


def extract_private_cases(private_b: dict, cases: list[dict]) -> dict[str, dict]:
    priv = private_b.get("cases") or {k: v for k, v in private_b.items() if k.startswith("L_")}
    expected = {c["anomaly_id"] for c in cases}
    if not expected.issubset(set(priv)):
        raise ValueError("private_B missing case IDs from public_B")
    return {aid: priv[aid] for aid in expected}


def extract_truth_cases(truth_b: dict, cases: list[dict]) -> dict[str, dict]:
    tcases = truth_b.get("cases") or {k: v for k, v in truth_b.items() if k.startswith("L_")}
    expected = {c["anomaly_id"] for c in cases}
    if not expected.issubset(set(tcases)):
        raise ValueError("truth_B missing case IDs from public_B")
    return {aid: tcases[aid] for aid in expected}


def gated_locus(loc: dict, proposal: dict, threshold: float | None) -> dict:
    if proposal["action"] != "ABSTAIN" and (threshold is None or proposal["confidence"] < threshold):
        return abstain_locus(loc, reason="calibration policy abstention")
    return proposal


def outcome_locus(case: dict, resolution: dict, truth_entry: dict, role: str) -> str:
    """Classify a single Universe B locus decision on a labelled or decoy case."""
    if resolution.get("error") or resolution.get("guardrail_status") == "ERROR":
        return "error"
    action = resolution.get("action", "ABSTAIN")
    if role == "decoy":
        if action == "ABSTAIN":
            return "unknown_abstain"
        return "unknown_call"
    if action == "ABSTAIN":
        return "abstained"
    positives = set(truth_entry.get("positives") or [])
    fid = resolution.get("function_id", "")
    if action == "PROMOTE":
        return "correct" if fid in positives else "unconfirmed"
    if action == "DEMOTE":
        return "harmful" if fid in positives else "unconfirmed"
    return "error"


def calibrate_locus(
    public: dict,
    proposals: dict,
    truth_b: dict,
    private_b: dict,
    precision: float = 0.90,
    min_decisions: int = 10,
) -> dict:
    if not 0 < precision <= 1 or min_decisions < 1:
        raise ValueError("invalid calibration requirements")
    cases = public_locus_cases(public)
    rows = checked_locus_predictions(public, proposals)
    priv_map = extract_private_cases(private_b, cases)
    truth_map = extract_truth_cases(truth_b, cases)

    calibration = [c for c in cases if c["split"] == "calibration"]
    n_cal_decoys = sum(1 for c in calibration if priv_map[c["anomaly_id"]]["role"] == "decoy")
    scores = sorted({
        float(rows[c["anomaly_id"]]["confidence"])
        for c in calibration
        if rows[c["anomaly_id"]]["action"] != "ABSTAIN"
    })

    ungated_counts: Counter = Counter()
    for c in calibration:
        aid = c["anomaly_id"]
        raw = rows[aid]
        status = outcome_locus(c, raw, truth_map[aid], priv_map[aid]["role"])
        if raw.get("error") or raw.get("guardrail_status") == "ERROR":
            status = "error"
        ungated_counts[status] += 1
    known_decisions_ungated = (
        ungated_counts["correct"] + ungated_counts["harmful"] + ungated_counts["unconfirmed"]
    )

    choices = []
    known_by_threshold = []
    prec_at_min_decisions = []
    for threshold in scores:
        counts: Counter = Counter()
        for c in calibration:
            aid = c["anomaly_id"]
            raw = rows[aid]
            res = gated_locus(c["locus"], raw, threshold)
            status = outcome_locus(c, res, truth_map[aid], priv_map[aid]["role"])
            if raw.get("error") or raw.get("guardrail_status") == "ERROR":
                status = "error"
            counts[status] += 1
        n_known = counts["correct"] + counts["harmful"] + counts["unconfirmed"]
        known_by_threshold.append(n_known)
        prec = counts["correct"] / n_known if n_known else 0.0
        unk_rate = counts["unknown_call"] / n_cal_decoys if n_cal_decoys else 0.0
        if n_known >= min_decisions:
            prec_at_min_decisions.append(prec)
        if n_known >= min_decisions and prec >= precision:
            choices.append((
                counts["correct"],
                -counts["harmful"],
                threshold,
                dict(counts),
                n_known,
                prec,
                unk_rate,
            ))

    best = max(choices, default=None)
    max_known_decisions_over_thresholds = max(known_by_threshold, default=0)
    best_precision_at_min_decisions = max(prec_at_min_decisions) if prec_at_min_decisions else None
    cal_truth_slice = {c["anomaly_id"]: truth_map[c["anomaly_id"]] for c in calibration}
    return {
        "schema": POLICY_SCHEMA,
        "mode": proposals["mode"],
        "public_sha256": E.digest(public),
        "calibration_truth_sha256": E.digest(cal_truth_slice),
        "calibration_proposals_sha256": E.digest([rows[c["anomaly_id"]] for c in calibration]),
        "target_precision": precision,
        "min_decisions": min_decisions,
        "threshold": best[2] if best else None,
        "calibration_counts": best[3] if best else {},
        "calibration_known_decisions": best[4] if best else 0,
        "calibration_precision": best[5] if best else None,
        "unknown_call_rate": best[6] if best else 0.0,
        "known_decisions_ungated": known_decisions_ungated,
        "max_known_decisions_over_thresholds": max_known_decisions_over_thresholds,
        "best_precision_at_min_decisions": best_precision_at_min_decisions,
        "status": "selected_on_calibration" if best else "no_supported_threshold_abstain",
        "selection_rule": "max_correct_then_min_harm_then_highest_threshold; unconfirmed counts against precision",
        "precision_note": PRECISION_NOTE,
        "limitation": "empirical calibration precision on locus-anchored calibration split, not a population guarantee",
    }


def evaluate_locus(
    public: dict,
    proposals: dict,
    policy: dict,
    truth_b: dict,
    private_b: dict,
) -> dict:
    if (
        policy.get("schema") != POLICY_SCHEMA
        or policy.get("public_sha256") != E.digest(public)
        or policy.get("mode") != proposals.get("mode")
    ):
        raise ValueError("missing or mismatched frozen calibration policy")
    cases = public_locus_cases(public)
    rows = checked_locus_predictions(public, proposals)
    cal_rows = [rows[c["anomaly_id"]] for c in cases if c["split"] == "calibration"]
    if policy.get("calibration_proposals_sha256") != E.digest(cal_rows):
        raise ValueError("calibration predictions changed after policy freeze")

    priv_map = extract_private_cases(private_b, cases)
    truth_map = extract_truth_cases(truth_b, cases)

    details = []
    counts: Counter = Counter()
    labelled_counts: Counter = Counter({k: 0 for k in ("correct", "harmful", "unconfirmed", "abstained", "error")})
    decoy_counts: Counter = Counter({k: 0 for k in ("unknown_call", "unknown_abstain", "error")})
    by_genome: dict[str, Counter] = defaultdict(Counter)

    n_labelled = 0
    n_labelled_with_positives = 0
    n_decoys = 0

    test_cases = [c for c in cases if c["split"] == "test"]
    for case in test_cases:
        aid = case["anomaly_id"]
        gid = case["genome_id"]
        role = priv_map[aid]["role"]
        t_entry = truth_map[aid]
        raw = rows[aid]
        resolution = gated_locus(case["locus"], raw, policy["threshold"])
        status = outcome_locus(case, resolution, t_entry, role)
        if raw.get("error") or raw.get("guardrail_status") == "ERROR":
            status = "error"

        counts[status] += 1
        by_genome[gid][status] += 1
        if role == "labelled":
            n_labelled += 1
            labelled_counts[status] += 1
            if t_entry.get("positives"):
                n_labelled_with_positives += 1
        else:
            n_decoys += 1
            decoy_counts[status] += 1

        details.append({
            "anomaly_id": aid,
            "genome_id": gid,
            "role": role,
            "outcome": status,
            "resolution": resolution,
        })

    n_known_decisions = labelled_counts["correct"] + labelled_counts["harmful"] + labelled_counts["unconfirmed"]
    correct_rate = labelled_counts["correct"] / n_labelled if n_labelled else None
    correct_ci = list(wilson(labelled_counts["correct"], n_labelled))
    harm_rate = labelled_counts["harmful"] / n_labelled if n_labelled else None
    harm_ci = list(wilson(labelled_counts["harmful"], n_labelled))
    abstention_rate = labelled_counts["abstained"] / n_labelled if n_labelled else None
    abstention_ci = list(wilson(labelled_counts["abstained"], n_labelled))
    recall_truth = labelled_counts["correct"] / n_labelled_with_positives if n_labelled_with_positives else None
    recall_ci = list(wilson(labelled_counts["correct"], n_labelled_with_positives))
    unknown_call_rate = decoy_counts["unknown_call"] / n_decoys if n_decoys else None
    unknown_call_ci = list(wilson(decoy_counts["unknown_call"], n_decoys))

    test_truth_slice = {c["anomaly_id"]: truth_map[c["anomaly_id"]] for c in test_cases}
    cost_block = {
        "decision_wall_seconds": proposals["decision_wall_seconds"],
        "api_usage": proposals.get("api_usage"),
        "shared_production_wall_seconds": public["shared_production_wall_seconds"],
    }
    return {
        "schema": EVALUATION_SCHEMA,
        "mode": proposals["mode"],
        "public_sha256": E.digest(public),
        "policy_sha256": E.digest(policy),
        "test_truth_sha256": E.digest(test_truth_slice),
        "counts": dict(counts),
        "labelled_counts": dict(labelled_counts),
        "decoy_counts": dict(decoy_counts),
        "n_cases": len(details),
        "n_labelled_test_cases": n_labelled,
        "n_labelled_with_positives": n_labelled_with_positives,
        "n_decoy_test_cases": n_decoys,
        "known_decision_precision": labelled_counts["correct"] / n_known_decisions if n_known_decisions else None,
        "correct_rate": correct_rate,
        "correct_rate_ci95": correct_ci,
        "harm_rate": harm_rate,
        "harm_rate_ci95": harm_ci,
        "abstention_rate": abstention_rate,
        "abstention_rate_ci95": abstention_ci,
        "recall_of_truth_ec": recall_truth,
        "recall_of_truth_ec_ci95": recall_ci,
        "unknown_call_rate": unknown_call_rate,
        "unknown_call_rate_ci95": unknown_call_ci,
        "by_genome": {g: dict(c) for g, c in sorted(by_genome.items())},
        "details": details,
        "cost": cost_block,
        "decision_wall_seconds": proposals["decision_wall_seconds"],
        "shared_production_wall_seconds": public["shared_production_wall_seconds"],
        "api_usage": proposals.get("api_usage"),
        "limitations": list(public.get("limitations", [])) + [
            PRECISION_NOTE,
            "Repeated locus decisions within a genome are clustered; genome-blocked bootstrap is reported in compare.",
        ],
    }


def leak_check_locus(public: dict, proposals_list: list[dict], private_b: dict) -> dict:
    """Test whether ungated call rate on labelled vs decoy cases across both splits is indistinguishable."""
    if not proposals_list:
        raise ValueError("at least one proposals object is required for leak-check")
    cases = public_locus_cases(public)
    priv_map = extract_private_cases(private_b, cases)

    labelled_ids = [c["anomaly_id"] for c in cases if priv_map[c["anomaly_id"]]["role"] == "labelled"]
    decoy_ids = [c["anomaly_id"] for c in cases if priv_map[c["anomaly_id"]]["role"] == "decoy"]

    methods = []
    by_mode = {}
    for prop in proposals_list:
        rows = checked_locus_predictions(public, prop)
        mode = prop["mode"]
        lab_called = sum(1 for aid in labelled_ids if rows[aid]["action"] != "ABSTAIN")
        lab_total = len(labelled_ids)
        lab_abstained = lab_total - lab_called
        lab_rate = lab_called / lab_total if lab_total else 0.0

        dec_called = sum(1 for aid in decoy_ids if rows[aid]["action"] != "ABSTAIN")
        dec_total = len(decoy_ids)
        dec_abstained = dec_total - dec_called
        dec_rate = dec_called / dec_total if dec_total else 0.0

        table = [[lab_called, lab_abstained], [dec_called, dec_abstained]]
        if lab_total > 0 and dec_total > 0:
            _odds, p_val = scipy_stats.fisher_exact(table, alternative="two-sided")
            p_val = float(p_val)
        else:
            p_val = 1.0
        diff = lab_rate - dec_rate
        abs_diff = abs(diff)
        passed = bool(p_val > 0.05 and abs_diff < 0.10)
        entry = {
            "mode": mode,
            "labelled_called": lab_called,
            "labelled_abstained": lab_abstained,
            "labelled_total": lab_total,
            "labelled_call_rate": lab_rate,
            "decoy_called": dec_called,
            "decoy_abstained": dec_abstained,
            "decoy_total": dec_total,
            "decoy_call_rate": dec_rate,
            "table_2x2": table,
            "diff": diff,
            "abs_diff": abs_diff,
            "fisher_exact_p": p_val,
            "passed": passed,
            "verdict": "PASS" if passed else "FAIL",
        }
        methods.append(entry)
        by_mode[mode] = entry

    return {
        "schema": LEAK_SCHEMA,
        "public_sha256": E.digest(public),
        "rule": "two-sided Fisher exact p > 0.05 and |labelled_call_rate - decoy_call_rate| < 0.10 over both splits (ungated)",
        "all_passed": all(m["passed"] for m in methods),
        "methods": methods,
        "by_mode": by_mode,
    }


def _percentile_ci(draws: list[float], point: float, alpha: float = 0.05) -> dict:
    clean = sorted(x for x in draws if not math.isnan(x))
    if not clean:
        return {"estimate": point, "ci_low": None, "ci_high": None, "n_boot": 0}
    lo = clean[int(math.floor(alpha / 2 * (len(clean) - 1)))]
    hi = clean[int(math.ceil((1 - alpha / 2) * (len(clean) - 1)))]
    return {"estimate": point, "ci_low": lo, "ci_high": hi, "n_boot": len(clean)}


def compare_locus_reports(
    public: dict,
    reports: list[dict],
    private_b: dict,
    n_boot: int = 1000,
    seed: int = 20261102,
) -> dict:
    """Paired case-level and genome-blocked bootstrap comparisons between reports[0] and reports[1:]."""
    if len(reports) < 2:
        raise ValueError("at least two evaluation reports are required for compare")
    cases = public_locus_cases(public)
    priv_map = extract_private_cases(private_b, cases)
    pub_sha = E.digest(public)

    for r in reports:
        if r.get("schema") != EVALUATION_SCHEMA or r.get("public_sha256") != pub_sha:
            raise ValueError("evaluation report does not match Universe B public bundle")

    labelled_test = sorted(
        (c for c in cases if c["split"] == "test" and priv_map[c["anomaly_id"]]["role"] == "labelled"),
        key=lambda c: c["anomaly_id"],
    )
    aid_list = [c["anomaly_id"] for c in labelled_test]
    gid_by_aid = {c["anomaly_id"]: c["genome_id"] for c in labelled_test}
    gids = sorted({c["genome_id"] for c in labelled_test})
    aids_by_gid: dict[str, list[str]] = defaultdict(list)
    for aid in aid_list:
        aids_by_gid[gid_by_aid[aid]].append(aid)

    def outcomes_of(rep: dict) -> dict[str, str]:
        by_aid = {d["anomaly_id"]: d["outcome"] for d in rep.get("details", [])}
        if not set(aid_list).issubset(set(by_aid)):
            raise ValueError(f"report {rep.get('mode')} missing labelled test case details")
        return by_aid

    ref = reports[0]
    ref_mode = ref["mode"]
    ref_out = outcomes_of(ref)

    contrasts = []
    by_comparator = {}
    n_cases = len(aid_list)

    for other in reports[1:]:
        other_mode = other["mode"]
        other_out = outcomes_of(other)

        diff_corr_vec = [
            (1.0 if ref_out[a] == "correct" else 0.0) - (1.0 if other_out[a] == "correct" else 0.0)
            for a in aid_list
        ]
        diff_harm_vec = [
            (1.0 if ref_out[a] == "harmful" else 0.0) - (1.0 if other_out[a] == "harmful" else 0.0)
            for a in aid_list
        ]

        pt_corr = sum(diff_corr_vec) / n_cases if n_cases else 0.0
        pt_harm = sum(diff_harm_vec) / n_cases if n_cases else 0.0

        rng_case = random.Random(seed)
        corr_draws_case = []
        harm_draws_case = []
        if n_cases > 0:
            for _ in range(n_boot):
                idxs = [rng_case.randrange(n_cases) for _ in range(n_cases)]
                corr_draws_case.append(sum(diff_corr_vec[i] for i in idxs) / n_cases)
                harm_draws_case.append(sum(diff_harm_vec[i] for i in idxs) / n_cases)

        g_stats = {}
        for g in gids:
            g_aids = aids_by_gid[g]
            g_stats[g] = (
                sum((1.0 if ref_out[a] == "correct" else 0.0) - (1.0 if other_out[a] == "correct" else 0.0) for a in g_aids),
                sum((1.0 if ref_out[a] == "harmful" else 0.0) - (1.0 if other_out[a] == "harmful" else 0.0) for a in g_aids),
                len(g_aids),
            )

        rng_genome = random.Random(seed)
        corr_draws_genome = []
        harm_draws_genome = []
        if len(gids) >= 2:
            for _ in range(n_boot):
                sampled_g = [gids[rng_genome.randrange(len(gids))] for _ in gids]
                tot_n = sum(g_stats[g][2] for g in sampled_g)
                if tot_n > 0:
                    corr_draws_genome.append(sum(g_stats[g][0] for g in sampled_g) / tot_n)
                    harm_draws_genome.append(sum(g_stats[g][1] for g in sampled_g) / tot_n)

        entry = {
            "reference_mode": ref_mode,
            "comparator_mode": other_mode,
            "n_labelled_test_cases": n_cases,
            "n_test_genomes": len(gids),
            "correct_rate_diff": {
                "paired_bootstrap": _percentile_ci(corr_draws_case, pt_corr),
                "genome_blocked_bootstrap": {
                    **_percentile_ci(corr_draws_genome, pt_corr),
                    "n_genomes": len(gids),
                },
            },
            "harm_rate_diff": {
                "paired_bootstrap": _percentile_ci(harm_draws_case, pt_harm),
                "genome_blocked_bootstrap": {
                    **_percentile_ci(harm_draws_genome, pt_harm),
                    "n_genomes": len(gids),
                },
            },
        }
        contrasts.append(entry)
        by_comparator[other_mode] = entry

    return {
        "schema": COMPARE_SCHEMA,
        "public_sha256": pub_sha,
        "reference_mode": ref_mode,
        "n_boot": n_boot,
        "seed": seed,
        "n_labelled_test_cases": n_cases,
        "n_test_genomes": len(gids),
        "contrasts": contrasts,
        "by_comparator": by_comparator,
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="cmd", required=True)

    cal = sub.add_parser("calibrate", help="Fit calibration threshold policy on labelled calibration cases")
    cal.add_argument("--public", required=True, type=Path)
    cal.add_argument("--proposals", required=True, type=Path)
    cal.add_argument("--truth", required=True, type=Path)
    cal.add_argument("--private", required=True, type=Path)
    cal.add_argument("--precision", type=float, default=0.90)
    cal.add_argument("--min-decisions", type=int, default=10)
    cal.add_argument("--out", required=True, type=Path)

    ev = sub.add_parser("evaluate", help="Evaluate frozen policy on test split")
    ev.add_argument("--public", required=True, type=Path)
    ev.add_argument("--proposals", required=True, type=Path)
    ev.add_argument("--policy", required=True, type=Path)
    ev.add_argument("--truth", required=True, type=Path)
    ev.add_argument("--private", required=True, type=Path)
    ev.add_argument("--out", required=True, type=Path)

    lk = sub.add_parser("leak-check", help="Check ungated call-rate parity between labelled and decoy cases")
    lk.add_argument("--public", required=True, type=Path)
    lk.add_argument("--proposals", required=True, nargs="+", type=Path)
    lk.add_argument("--private", required=True, type=Path)
    lk.add_argument("--out", required=True, type=Path)

    cmp_p = sub.add_parser("compare", help="Paired bootstrap and genome-blocked bootstrap comparisons")
    cmp_p.add_argument("--public", required=True, type=Path)
    cmp_p.add_argument("--reports", required=True, nargs="+", type=Path)
    cmp_p.add_argument("--private", required=True, type=Path)
    cmp_p.add_argument("--n-boot", type=int, default=1000)
    cmp_p.add_argument("--seed", type=int, default=20261102)
    cmp_p.add_argument("--out", required=True, type=Path)

    args = parser.parse_args(argv)
    if args.cmd == "calibrate":
        res = calibrate_locus(
            public=E.read(args.public),
            proposals=E.read(args.proposals),
            truth_b=E.read(args.truth),
            private_b=E.read(args.private),
            precision=args.precision,
            min_decisions=args.min_decisions,
        )
    elif args.cmd == "evaluate":
        res = evaluate_locus(
            public=E.read(args.public),
            proposals=E.read(args.proposals),
            policy=E.read(args.policy),
            truth_b=E.read(args.truth),
            private_b=E.read(args.private),
        )
    elif args.cmd == "leak-check":
        res = leak_check_locus(
            public=E.read(args.public),
            proposals_list=[E.read(p) for p in args.proposals],
            private_b=E.read(args.private),
        )
    elif args.cmd == "compare":
        res = compare_locus_reports(
            public=E.read(args.public),
            reports=[E.read(p) for p in args.reports],
            private_b=E.read(args.private),
            n_boot=args.n_boot,
            seed=args.seed,
        )
    E.write(args.out, res)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
