#!/usr/bin/env python3
"""Pure evaluation functions for PU candidate-layer and null-layer artifacts (row 4).

Unconfirmed rows are never counted as negatives or false positives anywhere in this module.
Explicit negatives (GOA experimental NOT, reviewed and accepted Swiss-Prot activity negations,
and table rows with ``label_basis == 'explicit_negative_truth'``) are evaluated separately for
harms and specificity.
"""
from __future__ import annotations

import csv
import math
import sys
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.stats import rankdata

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "benchmark" / "v3"))
sys.path.insert(0, str(REPO / "scripts"))
import train_context_model as T  # noqa: E402
from pu_learning import pu_trainer as PT  # noqa: E402


# ----------------------------------------------------------------------------- helpers

def norm_function_id(function_type: str | None, value: str | None) -> str:
    if value is None:
        return ""
    v = str(value).strip().replace("_", ":")
    if not v:
        return ""
    if str(function_type or "").upper() == "EC" and not v.upper().startswith("EC:"):
        return "EC:" + v
    if not v.upper().startswith(("EC:", "GO:")) and v[0].isdigit() and "." in v:
        return "EC:" + v
    return v


def extract_protein_id(row: dict) -> str:
    ck = row.get("candidate_key")
    if isinstance(ck, str) and "|" in ck:
        return ck.split("|", 1)[0].strip()
    pid = row.get("protein_id")
    if isinstance(pid, str) and pid.strip():
        return pid.strip()
    eid = row.get("example_id")
    if isinstance(eid, str) and "|" in eid:
        parts = eid.split("|")
        if len(parts) >= 2:
            return parts[1].strip()
    return str(ck or pid or eid or "").strip()


def _auroc_rank(y_true: np.ndarray, scores: np.ndarray) -> float | None:
    y = np.asarray(y_true, dtype=int)
    s = np.asarray(scores, dtype=float)
    npos = int(np.sum(y == 1))
    nneg = int(np.sum(y == 0))
    if npos == 0 or nneg == 0:
        return None
    ranks = rankdata(s)
    return float(( np.sum(ranks[y == 1]) - npos * (npos + 1) / 2.0 ) / (npos * nneg))


def _spearman_corr(x: np.ndarray, y: np.ndarray) -> float | None:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(x) < 2 or len(x) != len(y):
        return None
    rx = rankdata(x)
    ry = rankdata(y)
    sx = float(np.std(rx))
    sy = float(np.std(ry))
    if sx < 1e-12 or sy < 1e-12:
        return None
    return float(np.corrcoef(rx, ry)[0, 1])


def _is_explicit_negative(row: dict) -> bool:
    return row.get("label_basis") == "explicit_negative_truth" or bool(row.get("is_explicit_negative"))


def _is_positive(row: dict) -> bool:
    return not _is_explicit_negative(row) and int(row["label"]) == 1


def _is_unlabelled(row: dict) -> bool:
    return not _is_explicit_negative(row) and int(row["label"]) == 0


# ----------------------------------------------------------------------------- 1. score

def score(art: dict, rows: list[dict]) -> list[dict]:
    """Score each row with ``art``, recording probability, margin, or ``unsupported_reason``.

    Abstentions are retained with ``probability=None`` and ``unsupported_reason=<reason>``;
    no row is ever silently dropped.
    """
    threshold = float(art.get("display_threshold", {}).get("value", 0.5))
    out: list[dict] = []
    for r in rows:
        x = np.asarray(r["x"], dtype=float)
        seq_sources = r.get("sequence_sources")
        qual_ev = r.get("quality_evidence")
        ctx_mode = r.get("cohort_context_mode", "UNIVERSAL")
        seq_contracts = r.get("sequence_contracts")
        reason = T.unsupported_input_reason(art, x, seq_sources, qual_ev, ctx_mode, seq_contracts)
        if reason is None:
            prob = float(T.predict_proba(art, x, seq_sources, qual_ev, ctx_mode, seq_contracts))
            marg = float(T.predict_margin(art, x))
        else:
            prob = None
            marg = None
        scored = dict(r)
        scored.update(
            probability=prob,
            margin=marg,
            unsupported_reason=reason,
            supported=(reason is None),
            threshold=threshold,
        )
        out.append(scored)
    return out


# ----------------------------------------------------------------------------- 2. candidate_metrics

def candidate_metrics(scored_rows: list[dict], c: float, threshold: float) -> dict:
    """Candidate-layer metrics on supported rows, with abstention accounting over all supplied rows."""
    total_rows = len(scored_rows)
    abstained = [r for r in scored_rows if r.get("unsupported_reason") is not None or r.get("probability") is None]
    sup = [r for r in scored_rows if r.get("unsupported_reason") is None and r.get("probability") is not None]
    unsupported_reasons = dict(sorted(Counter(r["unsupported_reason"] for r in abstained if r.get("unsupported_reason")).items()))

    all_split_positives = sum(1 for r in scored_rows if _is_positive(r))
    pos_rows = [r for r in sup if _is_positive(r)]
    unl_rows = [r for r in sup if _is_unlabelled(r)]
    exp_neg_rows = [r for r in sup if _is_explicit_negative(r)]

    n_sup = len(sup)
    n_pos = len(pos_rows)
    n_unl = len(unl_rows)
    n_exp_neg = len(exp_neg_rows)

    thr = float(threshold)
    calls = [r for r in sup if float(r["probability"]) >= thr]
    n_calls = len(calls)
    tp = sum(1 for r in calls if _is_positive(r))
    unl_called = sum(1 for r in calls if _is_unlabelled(r))
    exp_neg_called = sum(1 for r in calls if _is_explicit_negative(r))

    supported_precision = float(tp / n_calls) if n_calls > 0 else None
    pu_expected_precision = float(np.mean([float(r["probability"]) for r in calls])) if n_calls > 0 else None
    labelled_fraction_identity = (
        float((tp / n_calls) / float(c)) if (n_calls > 0 and c is not None and float(c) > 0.0) else None
    )
    recall = float(tp / all_split_positives) if all_split_positives > 0 else None
    recall_supported = float(tp / n_pos) if n_pos > 0 else None

    pu_sup = [r for r in sup if _is_positive(r) or _is_unlabelled(r)]
    if n_pos > 0 and n_unl > 0:
        p_pu = np.array([float(r["probability"]) for r in pu_sup], dtype=float)
        y_pu = np.array([1 if _is_positive(r) else 0 for r in pu_sup], dtype=int)
        auroc_pos_unl = _auroc_rank(y_pu, p_pu)
    else:
        p_pu = np.array([float(r["probability"]) for r in pu_sup], dtype=float)
        y_pu = np.array([1 if _is_positive(r) else 0 for r in pu_sup], dtype=int)
        auroc_pos_unl = None

    # Per-positive rank among the unlabelled rows of its own competition
    by_comp: dict[tuple, list[dict]] = {}
    for r in sup:
        key = PT.competition_key(r)
        if key is not None:
            by_comp.setdefault(key, []).append(r)

    ranks_first: list[float] = []
    recip_ranks: list[float] = []
    comps_no_unl_excluded = 0
    pos_no_unl_excluded = 0
    for _, members in sorted(by_comp.items()):
        c_pos = [m for m in members if _is_positive(m)]
        c_unl = [m for m in members if _is_unlabelled(m)]
        if not c_pos:
            continue
        if not c_unl:
            comps_no_unl_excluded += 1
            pos_no_unl_excluded += len(c_pos)
            continue
        unl_probs = [float(u["probability"]) for u in c_unl]
        for pr in c_pos:
            pp = float(pr["probability"])
            higher = sum(1 for up in unl_probs if up > pp)
            tied = sum(1 for up in unl_probs if up == pp)
            rk = 1.0 + higher + 0.5 * tied
            ranks_first.append(1.0 if (higher == 0 and tied == 0) else 0.0)
            recip_ranks.append(1.0 / rk)

    fraction_ranked_first = float(np.mean(ranks_first)) if ranks_first else None
    mean_reciprocal_rank = float(np.mean(recip_ranks)) if recip_ranks else None

    c_val = float(c) if c is not None else 1.0
    if len(pu_sup) > 0:
        rel_table, pu_ece, pn_ece_ref = PT.pu_reliability_table(p_pu, y_pu, c_val)
        brier_c = float(np.mean((c_val * p_pu - y_pu) ** 2))
    else:
        rel_table, pu_ece, pn_ece_ref = [], None, None
        brier_c = None

    return {
        "total_rows": int(total_rows),
        "n": int(n_sup),
        "abstained": int(len(abstained)),
        "unsupported_reasons": unsupported_reasons,
        "positives": int(n_pos),
        "all_split_positives": int(all_split_positives),
        "unlabelled": int(n_unl),
        "explicit_negatives": int(n_exp_neg),
        "threshold": thr,
        "c": c_val,
        "n_calls": int(n_calls),
        "labelled_positives_called": int(tp),
        "unconfirmed_called": int(unl_called),
        "explicit_negatives_called": int(exp_neg_called),
        "supported_precision": supported_precision,
        "pu_expected_precision": pu_expected_precision,
        "pu_adjusted_expected_precision": pu_expected_precision,
        "labelled_fraction_identity": labelled_fraction_identity,
        "recall": recall,
        "recall_supported": recall_supported,
        "auroc_positive_vs_unlabelled": auroc_pos_unl,
        "fraction_ranked_first": fraction_ranked_first,
        "mean_reciprocal_rank": mean_reciprocal_rank,
        "competitions_without_unlabelled_excluded": int(comps_no_unl_excluded),
        "positives_without_unlabelled_excluded": int(pos_no_unl_excluded),
        "competition_ranking": {
            "n_positives_evaluated": int(len(ranks_first)),
            "fraction_ranked_first": fraction_ranked_first,
            "mean_reciprocal_rank": mean_reciprocal_rank,
            "competitions_without_unlabelled_excluded": int(comps_no_unl_excluded),
            "positives_without_unlabelled_excluded": int(pos_no_unl_excluded),
        },
        "pu_reliability_table": rel_table,
        "pu_ece": pu_ece,
        "pn_ece_reference": pn_ece_ref,
        "brier_labelled_c_scaled": brier_c,
    }


# ----------------------------------------------------------------------------- 3. pairwise_contrast

def pairwise_contrast(
    scored_a: list[dict],
    scored_b: list[dict],
    threshold_a: float,
    threshold_b: float,
) -> dict:
    """Pairwise contrast on rows supported by both arms.

    - ``corrections`` = positives called only by A (and reverse for B).
    - ``harms`` = explicit negatives called only by A (and reverse for B).
    - ``unconfirmed_called_only_by_a`` / ``unconfirmed_called_only_by_b`` are reported separately
      and NEVER added to harms.
    """
    thr_a = float(threshold_a)
    thr_b = float(threshold_b)
    if len(scored_a) == len(scored_b) and all(
        a.get("example_id") == b.get("example_id") for a, b in zip(scored_a, scored_b)
    ):
        pairs = list(zip(scored_a, scored_b))
    else:
        by_id_b = {r["example_id"]: r for r in scored_b}
        pairs = [(a, by_id_b[a["example_id"]]) for a in scored_a if a["example_id"] in by_id_b]

    both_sup = [
        (a, b)
        for a, b in pairs
        if a.get("unsupported_reason") is None
        and a.get("probability") is not None
        and b.get("unsupported_reason") is None
        and b.get("probability") is not None
    ]

    pos_only_a = 0
    pos_only_b = 0
    pos_both = 0
    exp_neg_only_a = 0
    exp_neg_only_b = 0
    exp_neg_both = 0
    unc_only_a = 0
    unc_only_b = 0
    unc_both = 0

    for a, b in both_sup:
        call_a = float(a["probability"]) >= thr_a
        call_b = float(b["probability"]) >= thr_b
        only_a = call_a and not call_b
        only_b = call_b and not call_a
        both = call_a and call_b

        if _is_explicit_negative(a) or _is_explicit_negative(b):
            if only_a:
                exp_neg_only_a += 1
            elif only_b:
                exp_neg_only_b += 1
            elif both:
                exp_neg_both += 1
        elif _is_positive(a):
            if only_a:
                pos_only_a += 1
            elif only_b:
                pos_only_b += 1
            elif both:
                pos_both += 1
        else:
            if only_a:
                unc_only_a += 1
            elif only_b:
                unc_only_b += 1
            elif both:
                unc_both += 1

    return {
        "n_both_supported": int(len(both_sup)),
        "threshold_a": thr_a,
        "threshold_b": thr_b,
        "corrections": int(pos_only_a),
        "corrections_a_over_b": int(pos_only_a),
        "corrections_b_over_a": int(pos_only_b),
        "harms": int(exp_neg_only_a),
        "harms_a_over_b": int(exp_neg_only_a),
        "harms_b_over_a": int(exp_neg_only_b),
        "unconfirmed_called_only_by_a": int(unc_only_a),
        "unconfirmed_called_only_by_b": int(unc_only_b),
        "both_called_positives": int(pos_both),
        "both_called_explicit_negatives": int(exp_neg_both),
        "both_called_unconfirmed": int(unc_both),
    }


# ----------------------------------------------------------------------------- 4. negative_specificity

def negative_specificity(scored_rows: list[dict], threshold: float | None = None) -> dict:
    """Among explicit negatives, compute the fraction below threshold and list their probabilities."""
    exp_neg = [r for r in scored_rows if _is_explicit_negative(r)]
    if threshold is None:
        threshold = float(scored_rows[0]["threshold"]) if scored_rows and "threshold" in scored_rows[0] else 0.5
    thr = float(threshold)

    sup = [r for r in exp_neg if r.get("unsupported_reason") is None and r.get("probability") is not None]
    abstained = [r for r in exp_neg if r.get("unsupported_reason") is not None or r.get("probability") is None]
    probs = [float(r["probability"]) for r in sup]
    n_below = sum(1 for p in probs if p < thr)
    n_called = sum(1 for p in probs if p >= thr)
    spec = float(n_below / len(sup)) if sup else None

    details = []
    for r in exp_neg:
        p = r.get("probability")
        details.append(
            {
                "example_id": r.get("example_id"),
                "genome_id": r.get("genome_id"),
                "protein_id": extract_protein_id(r),
                "function_id": norm_function_id(r.get("function_type"), r.get("function_id")),
                "candidate_key": r.get("candidate_key"),
                "split": r.get("split"),
                "descriptive_train_genome": (r.get("split") == "train"),
                "probability": float(p) if p is not None else None,
                "unsupported_reason": r.get("unsupported_reason"),
                "below_threshold": (float(p) < thr) if p is not None else None,
            }
        )

    return {
        "threshold": thr,
        "n_explicit_negatives": int(len(exp_neg)),
        "n_supported": int(len(sup)),
        "n_abstained": int(len(abstained)),
        "n_below_threshold": int(n_below),
        "n_called": int(n_called),
        "fraction_below_threshold": spec,
        "specificity": spec,
        "probabilities": probs,
        "rows": details,
    }


# ----------------------------------------------------------------------------- 5. deletion_test

def _assign_locus_for_members(art: dict, members: list[dict]) -> tuple[float, list[float]]:
    vecs = [np.asarray(m["x"], dtype=float) for m in members]
    has_sources = any(m.get("sequence_sources") is not None for m in members)
    has_qual = any(m.get("quality_evidence") is not None for m in members)
    has_contracts = any(m.get("sequence_contracts") is not None for m in members)
    return T.assign_locus(
        art,
        vecs,
        sequence_sources=[m.get("sequence_sources") for m in members] if has_sources else None,
        quality_evidence=[m.get("quality_evidence") for m in members] if has_qual else None,
        cohort_context_modes=[m.get("cohort_context_mode", "UNIVERSAL") for m in members],
        sequence_contracts=[m.get("sequence_contracts") for m in members] if has_contracts else None,
    )


def deletion_test(art: dict, rows: list[dict], *, min_members: int = 2) -> dict:
    """Deletion test on single-positive test competitions with >= min_members supported members.

    Computes P(null) intact and after deleting the positive member (via T.assign_locus when the
    null layer is trained), plus the naive null 1 - max p_hat for any artifact.
    """
    test_rows = [r for r in rows if r.get("split") == "test"] if any(r.get("split") == "test" for r in rows) else list(rows)
    by_comp: dict[tuple, list[dict]] = {}
    for r in test_rows:
        if _is_explicit_negative(r):
            continue
        key = PT.competition_key(r)
        if key is not None:
            by_comp.setdefault(key, []).append(r)

    locus_trained = bool(art.get("locus_model", {}).get("trained"))
    records: list[dict] = []
    n_unsupported = 0

    for key, members in sorted(by_comp.items()):
        if len(members) < min_members:
            continue
        pos_indices = [j for j, m in enumerate(members) if _is_positive(m)]
        if len(pos_indices) != 1:
            continue
        scored_members = score(art, members)
        if any(not sm["supported"] for sm in scored_members):
            n_unsupported += 1
            continue

        pos_idx = pos_indices[0]
        rem_members = [m for j, m in enumerate(members) if j != pos_idx]
        p_intact = [float(sm["probability"]) for sm in scored_members]
        p_deleted = [p for j, p in enumerate(p_intact) if j != pos_idx]

        naive_intact = float(1.0 - max(p_intact))
        naive_deleted = float(1.0 - max(p_deleted))

        if locus_trained:
            p0_intact, _ = _assign_locus_for_members(art, members)
            p0_deleted, _ = _assign_locus_for_members(art, rem_members)
        else:
            p0_intact = None
            p0_deleted = None

        records.append(
            {
                "competition": list(key),
                "genome_id": key[0],
                "function_type": key[1],
                "function_id": key[2],
                "role": key[3],
                "n_members_intact": int(len(members)),
                "n_members_deleted": int(len(rem_members)),
                "positive_example_id": members[pos_idx].get("example_id"),
                "p_null_intact": float(p0_intact) if p0_intact is not None else None,
                "p_null_deleted": float(p0_deleted) if p0_deleted is not None else None,
                "naive_null_intact": naive_intact,
                "naive_null_deleted": naive_deleted,
            }
        )

    n_comp = len(records)
    if n_comp > 0:
        y_del = np.array([0] * n_comp + [1] * n_comp, dtype=int)
        naive_scores = np.array(
            [r["naive_null_intact"] for r in records] + [r["naive_null_deleted"] for r in records],
            dtype=float,
        )
        auroc_naive = _auroc_rank(y_del, naive_scores)
        mean_naive_intact = float(np.mean([r["naive_null_intact"] for r in records]))
        mean_naive_deleted = float(np.mean([r["naive_null_deleted"] for r in records]))
        if locus_trained:
            p_null_scores = np.array(
                [r["p_null_intact"] for r in records] + [r["p_null_deleted"] for r in records],
                dtype=float,
            )
            auroc_p_null = _auroc_rank(y_del, p_null_scores)
            mean_p_null_intact = float(np.mean([r["p_null_intact"] for r in records]))
            mean_p_null_deleted = float(np.mean([r["p_null_deleted"] for r in records]))
        else:
            auroc_p_null = None
            mean_p_null_intact = None
            mean_p_null_deleted = None
    else:
        auroc_naive = None
        mean_naive_intact = None
        mean_naive_deleted = None
        auroc_p_null = None
        mean_p_null_intact = None
        mean_p_null_deleted = None

    return {
        "locus_trained": locus_trained,
        "n_competitions": int(n_comp),
        "n_unsupported_competitions_skipped": int(n_unsupported),
        "mean_p_null_intact": mean_p_null_intact,
        "mean_p_null_deleted": mean_p_null_deleted,
        "auroc_p_null": auroc_p_null,
        "mean_naive_null_intact": mean_naive_intact,
        "mean_naive_null_deleted": mean_naive_deleted,
        "auroc_naive_null": auroc_naive,
        "competitions": records,
    }


# ----------------------------------------------------------------------------- 6. null_layer_metrics

def null_layer_metrics(art: dict, rows: list[dict], split: str | None, c: float) -> dict:
    """Evaluate the null layer on groups with >= 2 members, plus singleton groups and candidate-count checks."""
    sub = [r for r in rows if (split is None or r.get("split") == split)]
    scored = score(art, sub)

    # Identify competitions with any unsupported member
    by_comp_all: dict[tuple, list[dict]] = {}
    for sr in scored:
        key = PT.competition_key(sr)
        if key is not None:
            by_comp_all.setdefault(key, []).append(sr)

    unsupported_comp_keys = {k for k, ms in by_comp_all.items() if any(not m["supported"] for m in ms)}
    sup_rows = [
        sr for sr in scored
        if sr["supported"] and not _is_explicit_negative(sr) and PT.competition_key(sr) not in unsupported_comp_keys
    ]
    margins = [float(sr["margin"]) for sr in sup_rows]
    names = art["contract"]["features"]

    groups_ge2, multi_ge2 = PT.build_pu_groups(sup_rows, margins, names, min_members=2)
    groups_all, _ = PT.build_pu_groups(sup_rows, margins, names, min_members=1)
    singletons = [g for g in groups_all if len(g["margins"]) == 1]

    pos_groups = [g for g in groups_ge2 if int(g["target"]) > 0]
    no_label_groups = [g for g in groups_ge2 if int(g["target"]) == 0]
    pos_singletons = [g for g in singletons if int(g["target"]) > 0]
    no_label_singletons = [g for g in singletons if int(g["target"]) == 0]

    c_val = float(c)
    uniform_nll_all = PT.uniform_pu_locus_nll(groups_ge2, c_val) if groups_ge2 else None
    uniform_nll_pos = PT.uniform_pu_locus_nll(pos_groups, c_val) if pos_groups else None
    uniform_nll_no_label = PT.uniform_pu_locus_nll(no_label_groups, c_val) if no_label_groups else None
    mean_log_n1 = PT.plain_uniform_nll(groups_ge2) if groups_ge2 else None

    lm = art.get("locus_model", {})
    trained = bool(lm.get("trained"))
    if not trained:
        return {
            "split": split,
            "trained": False,
            "reason": lm.get("reason"),
            "c": c_val,
            "n_groups_ge2": int(len(groups_ge2)),
            "n_single_positive_groups": int(len(pos_groups)),
            "n_no_label_groups": int(len(no_label_groups)),
            "n_multi_positive_groups_excluded": int(multi_ge2),
            "n_unsupported_groups_skipped": int(len(unsupported_comp_keys)),
            "marginal_nll": None,
            "marginal_nll_single_positive": None,
            "marginal_nll_no_label": None,
            "uniform_marginal_nll": uniform_nll_all,
            "uniform_marginal_nll_single_positive": uniform_nll_pos,
            "uniform_marginal_nll_no_label": uniform_nll_no_label,
            "mean_log_n_plus_1": mean_log_n1,
            "mean_p_null_single_positive": None,
            "mean_p_null_no_label": None,
            "candidate_count_check": {
                "spearman_p_null_vs_log1p_n": None,
                "spearman_p_null_vs_max_margin": None,
            },
            "singleton_groups": {
                "n_groups": int(len(singletons)),
                "n_single_positive": int(len(pos_singletons)),
                "n_no_label": int(len(no_label_singletons)),
                "mean_p_null_single_positive": None,
                "mean_p_null_no_label": None,
            },
        }

    theta = np.concatenate(([math.log(float(lm["scale"])), float(lm["null_bias"])], [float(v) for v in lm["null_weights"]]))
    marginal_nll_all = PT.pu_locus_nll(theta, groups_ge2, c_val, 0.0) if groups_ge2 else None
    marginal_nll_pos = PT.pu_locus_nll(theta, pos_groups, c_val, 0.0) if pos_groups else None
    marginal_nll_no_label = PT.pu_locus_nll(theta, no_label_groups, c_val, 0.0) if no_label_groups else None

    p_null_pos = [T.locus_probabilities(theta, g["margins"], g["g"])[0] for g in pos_groups]
    p_null_no_label = [T.locus_probabilities(theta, g["margins"], g["g"])[0] for g in no_label_groups]
    p_null_all = [T.locus_probabilities(theta, g["margins"], g["g"])[0] for g in groups_ge2]
    log1p_ns = [float(g["g"][3]) for g in groups_ge2]
    max_margins = [float(g["g"][0]) for g in groups_ge2]

    p_null_sing_pos = [T.locus_probabilities(theta, g["margins"], g["g"])[0] for g in pos_singletons]
    p_null_sing_no_label = [T.locus_probabilities(theta, g["margins"], g["g"])[0] for g in no_label_singletons]

    # Production Stage-1 calls removed when null wins (p_null > max candidate locus probability)
    thr = float(art.get("display_threshold", {}).get("value", 0.5))
    removed_pos = 0
    removed_exp_neg = 0
    removed_unc = 0
    for key, comp_members in sorted(by_comp_all.items()):
        if key in unsupported_comp_keys or not comp_members:
            continue
        p0, pc = _assign_locus_for_members(art, comp_members)
        if p0 > (max(pc) if len(pc) > 0 else 0.0):
            for cm in comp_members:
                if float(cm["probability"]) >= thr:
                    if _is_explicit_negative(cm):
                        removed_exp_neg += 1
                    elif _is_positive(cm):
                        removed_pos += 1
                    else:
                        removed_unc += 1

    return {
        "split": split,
        "trained": True,
        "c": c_val,
        "n_groups_ge2": int(len(groups_ge2)),
        "n_single_positive_groups": int(len(pos_groups)),
        "n_no_label_groups": int(len(no_label_groups)),
        "n_multi_positive_groups_excluded": int(multi_ge2),
        "n_unsupported_groups_skipped": int(len(unsupported_comp_keys)),
        "marginal_nll": marginal_nll_all,
        "marginal_nll_single_positive": marginal_nll_pos,
        "marginal_nll_no_label": marginal_nll_no_label,
        "uniform_marginal_nll": uniform_nll_all,
        "uniform_marginal_nll_single_positive": uniform_nll_pos,
        "uniform_marginal_nll_no_label": uniform_nll_no_label,
        "mean_log_n_plus_1": mean_log_n1,
        "mean_p_null_single_positive": float(np.mean(p_null_pos)) if p_null_pos else None,
        "mean_p_null_no_label": float(np.mean(p_null_no_label)) if p_null_no_label else None,
        "candidate_count_check": {
            "spearman_p_null_vs_log1p_n": _spearman_corr(np.array(p_null_all), np.array(log1p_ns)),
            "spearman_p_null_vs_max_margin": _spearman_corr(np.array(p_null_all), np.array(max_margins)),
        },
        "singleton_groups": {
            "n_groups": int(len(singletons)),
            "n_single_positive": int(len(pos_singletons)),
            "n_no_label": int(len(no_label_singletons)),
            "mean_p_null_single_positive": float(np.mean(p_null_sing_pos)) if p_null_sing_pos else None,
            "mean_p_null_no_label": float(np.mean(p_null_sing_no_label)) if p_null_sing_no_label else None,
        },
        "calls_removed_by_null_winning": {
            "positives_removed": int(removed_pos),
            "explicit_negatives_removed": int(removed_exp_neg),
            "unconfirmed_removed": int(removed_unc),
        },
    }


# ----------------------------------------------------------------------------- 7. negatives join

def load_accepted_negatives(
    registry_dir: Path | str,
    table_rows: list[dict] | None = None,
) -> set[tuple[str, str, str]]:
    """Return set of (genome_id, protein_id, function_id) from negative_reviews.tsv (disposition == accepted)
    plus any supplied table rows with label_basis == explicit_negative_truth.
    """
    reg_path = Path(registry_dir)
    if reg_path.is_file():
        tsv_path = reg_path
    elif (reg_path / "negative_reviews.tsv").exists():
        tsv_path = reg_path / "negative_reviews.tsv"
    elif (reg_path / "registry" / "negative_reviews.tsv").exists():
        tsv_path = reg_path / "registry" / "negative_reviews.tsv"
    else:
        raise FileNotFoundError(f"negative_reviews.tsv not found under {registry_dir}")

    accepted: set[tuple[str, str, str]] = set()
    with open(tsv_path, newline="") as fh:
        for rec in csv.DictReader(fh, delimiter="\t"):
            if str(rec.get("disposition", "")).strip() == "accepted":
                gid = str(rec["genome_id"]).strip()
                pid = str(rec["protein_id"]).strip()
                fid = norm_function_id(None, rec["function_id"])
                accepted.add((gid, pid, fid))

    if table_rows:
        for r in table_rows:
            if r.get("label_basis") == "explicit_negative_truth":
                gid = str(r["genome_id"]).strip()
                pid = extract_protein_id(r)
                fid = norm_function_id(r.get("function_type"), r.get("function_id"))
                accepted.add((gid, pid, fid))
    return accepted


def join_negatives_to_rows(
    accepted_negatives: set[tuple[str, str, str]],
    rows: list[dict],
) -> tuple[list[dict], list[dict]]:
    """Join accepted negatives to table/export rows by (genome_id, protein_id from candidate_key, function_id).

    Returns (joined_rows, unjoined_negatives_list). Unjoined negatives are explicitly listed, never dropped silently.
    """
    matched_keys: set[tuple[str, str, str]] = set()
    joined_rows: list[dict] = []
    for r in rows:
        gid = str(r["genome_id"]).strip()
        pid = extract_protein_id(r)
        fid = norm_function_id(r.get("function_type"), r.get("function_id"))
        key = (gid, pid, fid)
        if key in accepted_negatives or r.get("label_basis") == "explicit_negative_truth":
            matched_keys.add(key)
            nr = dict(r)
            nr["label"] = 0
            nr["label_basis"] = "explicit_negative_truth"
            nr["is_explicit_negative"] = True
            joined_rows.append(nr)

    unjoined = [
        {"genome_id": g, "protein_id": p, "function_id": f}
        for (g, p, f) in sorted(accepted_negatives - matched_keys)
    ]
    return joined_rows, unjoined
