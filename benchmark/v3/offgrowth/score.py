#!/usr/bin/env python3
"""Scoring, strata, call attribution, harm/unconfirmed bookkeeping, bootstrap, H1-H3, and controls (§2, §6, §7, §8)."""

from __future__ import annotations

import copy
import csv
from dataclasses import asdict, dataclass
import math
from pathlib import Path
import random
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple

import cobra

from .arms import (
    DEFAULT_FLUX_TOL,
    ArmRunResult,
    SolverReplayRecord,
    _find_biomass_rxn_id,
    apply_medium,
    map_reactions_to_candidate_genes,
)
from .mask import (
    CandidateClaim,
    TruthCase,
    base_seed_rxn_id,
    mask_candidates,
    protein_aliases,
    resolve_equivalent_ecs_and_reactions,
    strip_version,
)
from .panel import is_complete_ec, normalize_ec


BOOTSTRAP_SEED: int = 20261030
BOOTSTRAP_RESAMPLES: int = 2000
SHUFFLE_CONTROL_SEEDS: Tuple[int, ...] = (1, 2, 3, 4, 5)
HARM_EXCESS_LIMIT: float = 0.10
H3_RETENTION_TARGET: float = 0.80


@dataclass(frozen=True)
class StratumClassification:
    case_id: str
    genome_id: str
    protein_id: str
    stratum: str  # "on_growth" or "off_growth"
    wt_biomass: float
    ko_biomass: float
    ko_fraction: float


@dataclass(frozen=True)
class ReactionGeneCall:
    reaction_id: str
    protein_id: Optional[str]
    ec: Optional[str]
    score: float
    call_kind: str  # "supported", "harmful", "unconfirmed", "reaction_only"


@dataclass(frozen=True)
class CaseArmEvaluation:
    case_id: str
    genome_id: str
    protein_id: str
    truth_ecs: Tuple[str, ...]
    target_reactions: Tuple[str, ...]
    stratum: str
    arm: str
    reaction_recovered: bool
    recovered: bool
    calls: Tuple[ReactionGeneCall, ...]
    n_added_reactions: int
    n_gene_calls: int
    n_reaction_only_calls: int
    n_supported_calls: int
    n_harmful_calls: int
    n_unconfirmed_calls: int
    tests_considered: int
    tests_executed: int
    n_lps: int
    n_gapfills: int
    cpu_seconds: float
    wall_seconds: float
    peak_rss_kb: int
    case_source: str = "price148"


@dataclass(frozen=True)
class ArmSummaryMetrics:
    arm: str
    n_cases_total: int
    n_recovered_total: int
    recovery_rate_total: float
    recovery_ci_total: Tuple[float, float]
    n_cases_on_growth: int
    n_recovered_on_growth: int
    recovery_rate_on_growth: float
    recovery_ci_on_growth: Tuple[float, float]
    n_cases_off_growth: int
    n_recovered_off_growth: int
    recovery_rate_off_growth: float
    recovery_ci_off_growth: Tuple[float, float]
    n_gene_calls: int
    n_reaction_only_calls: int
    n_supported_calls: int
    n_harmful_calls: int
    n_unconfirmed_calls: int
    harm_rate: float
    harm_ci: Tuple[float, float]
    exceeds_a0_harm_gate: bool
    total_tests_considered: int
    total_tests_executed: int
    total_lps: int
    total_gapfills: int
    total_cpu_seconds: float
    total_wall_seconds: float
    max_peak_rss_kb: int
    lps_per_recovered_enzyme: Optional[float]
    n_cases_price148: int = 0
    n_recovered_price148: int = 0
    recovery_rate_price148: float = 0.0
    recovery_ci_price148: Tuple[float, float] = (0.0, 0.0)
    n_cases_price148_off_growth: int = 0
    n_recovered_price148_off_growth: int = 0
    recovery_rate_price148_off_growth: float = 0.0
    recovery_ci_price148_off_growth: Tuple[float, float] = (0.0, 0.0)
    n_cases_swissprot: int = 0
    n_recovered_swissprot: int = 0
    recovery_rate_swissprot: float = 0.0
    recovery_ci_swissprot: Tuple[float, float] = (0.0, 0.0)
    n_cases_swissprot_off_growth: int = 0
    n_recovered_swissprot_off_growth: int = 0
    recovery_rate_swissprot_off_growth: float = 0.0
    recovery_ci_swissprot_off_growth: Tuple[float, float] = (0.0, 0.0)


@dataclass(frozen=True)
class McNemarResult:
    hypothesis_id: str
    arm_treatment: str
    arm_baseline: str
    stratum: str
    n_pairs: int
    both_recovered: int
    treatment_only: int  # b
    baseline_only: int   # c
    neither_recovered: int
    p_value_raw: float
    p_value_holm: float
    significant_at_0_05: bool


@dataclass(frozen=True)
class H3EfficiencyResult:
    arm_budget: str
    arm_full: str
    off_growth_recovered_budget: int
    off_growth_recovered_full: int
    retention_ratio: float
    retention_ci: Tuple[float, float]
    gapfills_budget: int
    gapfills_full: int
    tests_executed_budget: int
    tests_considered_full: int
    test_budget_fraction: float
    passes_80pct_retention: bool


@dataclass(frozen=True)
class ShuffledSeedOutcome:
    seed: int
    arm: str
    n_cases: int
    n_recovered: int
    recovery_rate: float
    n_off_growth_cases: int
    n_off_growth_recovered: int
    off_growth_recovery_rate: float


@dataclass(frozen=True)
class ShuffledControlSummary:
    arm: str
    seeds: Tuple[int, ...]
    per_seed: Tuple[ShuffledSeedOutcome, ...]
    mean_recovery_rate: float
    max_recovery_rate: float
    mean_off_growth_recovery_rate: float


@dataclass(frozen=True)
class SolverReplaySummary:
    n_replayed: int
    n_agreed: int
    agreement_rate: float
    max_relative_obj_diff: float


def classify_case_stratum(
    case: TruthCase,
    unmasked_filled_model: cobra.Model,
    ec_to_reactions: Mapping[str, Sequence[str]],
    default_medium: Mapping[str, float],
    flux_tol: float = DEFAULT_FLUX_TOL,
) -> StratumClassification:
    """Classify `case` as `on_growth` or `off_growth` on the unmasked gap-filled model (§2).

    - `on_growth`: deleting the truth reactions from the unmasked, gap-filled
      model drops biomass below 1% of wild-type on the default medium.
    - `off_growth`: deleting the truth reactions leaves default-medium biomass
      at >= 1% of wild-type.
    """
    model = copy.deepcopy(unmasked_filled_model)
    apply_medium(model, default_medium)
    biomass_id = _find_biomass_rxn_id(model)
    model.objective = biomass_id
    model.objective_direction = "max"

    wt_raw = model.slim_optimize(error_value=0.0)
    wt_bio = 0.0 if wt_raw is None or math.isnan(wt_raw) else max(0.0, float(wt_raw))

    _, target_rxns = resolve_equivalent_ecs_and_reactions(case.truth_ecs, ec_to_reactions)
    target_base_set = {base_seed_rxn_id(r) for r in target_rxns}

    for rxn in model.reactions:
        if base_seed_rxn_id(rxn.id) in target_base_set:
            rxn.lower_bound = 0.0
            rxn.upper_bound = 0.0

    ko_raw = model.slim_optimize(error_value=0.0)
    ko_bio = 0.0 if ko_raw is None or math.isnan(ko_raw) else max(0.0, float(ko_raw))

    if wt_bio > flux_tol:
        ko_fraction = ko_bio / wt_bio
        stratum = "on_growth" if ko_fraction < 0.01 else "off_growth"
    else:
        ko_fraction = 1.0
        stratum = "off_growth"

    return StratumClassification(
        case_id=case.case_id,
        genome_id=case.genome_id,
        protein_id=case.protein_id,
        stratum=stratum,
        wt_biomass=wt_bio,
        ko_biomass=ko_bio,
        ko_fraction=ko_fraction,
    )


def assign_calls_for_added_reactions(
    added_reactions: Sequence[str],
    candidates: Sequence[CandidateClaim],
    ec_to_reactions: Mapping[str, Sequence[str]],
    known_experimental_ecs: Mapping[str, Sequence[str]],
) -> Tuple[ReactionGeneCall, ...]:
    """Assign each added reaction `r` to its top-scoring candidate gene and classify harm (§6).

    - Each reaction `r` is assigned to the masked GSPA candidate gene `g` with
      the highest score for an EC mapped to `r` (ties broken deterministically
      by ascending `protein_id`).
    - Reactions with no candidate gene are emitted with `protein_id=None` and
      `call_kind='reaction_only'`.
    - For a called gene `Q`:
      * If `Q` has no known experimental label in `known_experimental_ecs`,
        `call_kind='unconfirmed'` (never counted as a negative or in harm rate).
      * If `Q` has known experimental ECs $E_{\\text{exp}}(Q)$, the call is
        `'supported'` if `r` is mapped from any EC in $E_{\\text{exp}}(Q)$ and
        `'harmful'` otherwise.
    """
    # Map base reaction -> sorted (score desc, protein_id asc, ec asc)
    rxn_best_claims: Dict[str, List[Tuple[float, str, str]]] = {}
    for claim in candidates:
        ec = normalize_ec(claim.ec)
        for rxn_id in ec_to_reactions.get(ec, ()):
            base_r = base_seed_rxn_id(rxn_id)
            rxn_best_claims.setdefault(base_r, []).append(
                (-claim.score, claim.protein_id, ec)
            )
    for base_r in rxn_best_claims:
        rxn_best_claims[base_r].sort()

    # Normalize known experimental ECs per protein alias
    exp_ecs_by_protein: Dict[str, Set[str]] = {}
    for prot_id, ecs in known_experimental_ecs.items():
        norm_set = {normalize_ec(e) for e in ecs if is_complete_ec(e)}
        if not norm_set:
            continue
        for alias in protein_aliases(prot_id):
            exp_ecs_by_protein.setdefault(alias, set()).update(norm_set)

    calls: List[ReactionGeneCall] = []
    for raw_rxn in sorted({base_seed_rxn_id(r) for r in added_reactions}):
        ranked_claims = rxn_best_claims.get(raw_rxn, [])
        if not ranked_claims:
            calls.append(
                ReactionGeneCall(
                    reaction_id=raw_rxn,
                    protein_id=None,
                    ec=None,
                    score=0.0,
                    call_kind="reaction_only",
                )
            )
            continue

        neg_score, top_protein, top_ec = ranked_claims[0]
        top_score = -neg_score

        known_ecs = exp_ecs_by_protein.get(top_protein) or exp_ecs_by_protein.get(
            strip_version(top_protein)
        )
        if not known_ecs:
            call_kind = "unconfirmed"
        else:
            allowed_rxns: Set[str] = set()
            for kec in known_ecs:
                for kr in ec_to_reactions.get(kec, ()):
                    allowed_rxns.add(base_seed_rxn_id(kr))
            call_kind = "supported" if raw_rxn in allowed_rxns else "harmful"

        calls.append(
            ReactionGeneCall(
                reaction_id=raw_rxn,
                protein_id=top_protein,
                ec=top_ec,
                score=top_score,
                call_kind=call_kind,
            )
        )

    return tuple(calls)


def evaluate_case_arm(
    case: TruthCase,
    stratum: str,
    arm_result: ArmRunResult,
    call_candidates: Sequence[CandidateClaim],
    ec_to_reactions: Mapping[str, Sequence[str]],
    known_experimental_ecs: Mapping[str, Sequence[str]],
) -> CaseArmEvaluation:
    """Score one `(case, arm)` run for recovery, harm, unconfirmed calls, and LP cost (§6).

    Note on masking (§3) vs call attribution (§6):
    - During arm execution (§3-§5), `case.protein_id`'s claims for `case.truth_ecs`
      (and equivalent ECs) are removed from the arm inputs so they cannot influence
      draft GPRs, LP weights $w_r$, or A4 rule rankings.
    - During post-arm call attribution (§6), each reaction added by an arm is
      assigned to the highest-scoring GSPA candidate gene in `call_candidates`
      for an EC mapped to that reaction, testing whether the arm's added reactions
      recover `(P, E)`.
    """
    _, target_rxns = resolve_equivalent_ecs_and_reactions(
        case.truth_ecs, ec_to_reactions
    )
    target_base_set = {base_seed_rxn_id(r) for r in target_rxns}
    truth_aliases = protein_aliases(case.protein_id)

    added_base_set = {base_seed_rxn_id(r) for r in arm_result.added_reactions}
    reaction_recovered = bool(added_base_set & target_base_set)

    calls = assign_calls_for_added_reactions(
        added_reactions=arm_result.added_reactions,
        candidates=call_candidates,
        ec_to_reactions=ec_to_reactions,
        known_experimental_ecs=known_experimental_ecs,
    )

    recovered = False
    for call in calls:
        if call.protein_id is None:
            continue
        p_match = (
            call.protein_id in truth_aliases
            or strip_version(call.protein_id) in truth_aliases
        )
        if p_match and base_seed_rxn_id(call.reaction_id) in target_base_set:
            recovered = True
            break

    n_gene = sum(1 for c in calls if c.protein_id is not None)
    n_rxn_only = sum(1 for c in calls if c.call_kind == "reaction_only")
    n_supp = sum(1 for c in calls if c.call_kind == "supported")
    n_harm = sum(1 for c in calls if c.call_kind == "harmful")
    n_unconf = sum(1 for c in calls if c.call_kind == "unconfirmed")

    case_src = getattr(case, "case_source", None) or (
        "swissprot" if case.case_id.startswith("SP:") else "price148"
    )

    return CaseArmEvaluation(
        case_id=case.case_id,
        genome_id=case.genome_id,
        protein_id=case.protein_id,
        truth_ecs=case.truth_ecs,
        target_reactions=tuple(sorted(target_base_set)),
        stratum=stratum,
        arm=arm_result.arm,
        reaction_recovered=reaction_recovered,
        recovered=recovered,
        calls=calls,
        n_added_reactions=len(arm_result.added_reactions),
        n_gene_calls=n_gene,
        n_reaction_only_calls=n_rxn_only,
        n_supported_calls=n_supp,
        n_harmful_calls=n_harm,
        n_unconfirmed_calls=n_unconf,
        tests_considered=arm_result.tests_considered,
        tests_executed=arm_result.tests_executed,
        n_lps=arm_result.n_lps,
        n_gapfills=arm_result.n_gapfills,
        cpu_seconds=arm_result.cpu_seconds,
        wall_seconds=arm_result.wall_seconds,
        peak_rss_kb=arm_result.peak_rss_kb,
        case_source=case_src,
    )


def exact_one_sided_mcnemar_p(b_treatment_only: int, c_baseline_only: int) -> float:
    """One-sided exact McNemar p-value testing `b > c` ($H_0: p \\le 0.5$)."""
    n = b_treatment_only + c_baseline_only
    if n <= 0:
        return 1.0
    # Exact upper tail P(X >= b) for X ~ Binomial(n, 0.5)
    tail = 0.0
    for k in range(b_treatment_only, n + 1):
        tail += math.comb(n, k) * (0.5 ** n)
    return min(1.0, max(0.0, tail))


def holm_bonferroni_adjust(p_values: Sequence[float]) -> Tuple[float, ...]:
    """Holm-Bonferroni step-down multiple-testing correction."""
    m = len(p_values)
    if m == 0:
        return ()
    indexed = sorted(enumerate(p_values), key=lambda kv: kv[1])
    adjusted = [0.0] * m
    running_max = 0.0
    for rank, (orig_idx, pval) in enumerate(indexed):
        multiplier = m - rank
        adj = min(1.0, pval * multiplier)
        running_max = max(running_max, adj)
        adjusted[orig_idx] = running_max
    return tuple(adjusted)


def _percentile_ci(values: Sequence[float], alpha: float = 0.05) -> Tuple[float, float]:
    if not values:
        return (0.0, 0.0)
    ordered = sorted(values)
    n = len(ordered)
    lo_idx = max(0, min(n - 1, int(math.floor((alpha / 2.0) * (n - 1)))))
    hi_idx = max(0, min(n - 1, int(math.ceil((1.0 - alpha / 2.0) * (n - 1)))))
    return (ordered[lo_idx], ordered[hi_idx])


def summarize_arms_with_bootstrap(
    evaluations_by_arm: Mapping[str, Sequence[CaseArmEvaluation]],
    n_resamples: int = BOOTSTRAP_RESAMPLES,
    seed: int = BOOTSTRAP_SEED,
) -> Dict[str, ArmSummaryMetrics]:
    """Compute overall/stratum recovery, harm rates, costs, and genome-block 95% CIs (§6)."""
    raw_metrics: Dict[str, Dict[str, object]] = {}
    for arm, evals in evaluations_by_arm.items():
        n_total = len(evals)
        n_rec_total = sum(1 for e in evals if e.recovered)
        on_evals = [e for e in evals if e.stratum == "on_growth"]
        off_evals = [e for e in evals if e.stratum == "off_growth"]
        n_on = len(on_evals)
        n_rec_on = sum(1 for e in on_evals if e.recovered)
        n_off = len(off_evals)
        n_rec_off = sum(1 for e in off_evals if e.recovered)

        p_evals = [e for e in evals if getattr(e, "case_source", "price148") == "price148"]
        sp_evals = [e for e in evals if getattr(e, "case_source", "price148") == "swissprot"]
        n_p = len(p_evals)
        n_rec_p = sum(1 for e in p_evals if e.recovered)
        p_off = [e for e in p_evals if e.stratum == "off_growth"]
        n_p_off = len(p_off)
        n_rec_p_off = sum(1 for e in p_off if e.recovered)

        n_sp = len(sp_evals)
        n_rec_sp = sum(1 for e in sp_evals if e.recovered)
        sp_off = [e for e in sp_evals if e.stratum == "off_growth"]
        n_sp_off = len(sp_off)
        n_rec_sp_off = sum(1 for e in sp_off if e.recovered)

        n_gene = sum(e.n_gene_calls for e in evals)
        n_rxn_only = sum(e.n_reaction_only_calls for e in evals)
        n_supp = sum(e.n_supported_calls for e in evals)
        n_harm = sum(e.n_harmful_calls for e in evals)
        n_unconf = sum(e.n_unconfirmed_calls for e in evals)
        n_known = n_supp + n_harm
        harm_rate = (n_harm / float(n_known)) if n_known > 0 else 0.0

        # Genome-block bootstrap
        by_genome: Dict[str, List[CaseArmEvaluation]] = {}
        for e in evals:
            by_genome.setdefault(e.genome_id, []).append(e)
        genomes = sorted(by_genome.keys())

        rng = random.Random(seed)
        boot_tot: List[float] = []
        boot_on: List[float] = []
        boot_off: List[float] = []
        boot_p_tot: List[float] = []
        boot_p_off: List[float] = []
        boot_sp_tot: List[float] = []
        boot_sp_off: List[float] = []
        boot_harm: List[float] = []

        if genomes and n_resamples > 0:
            k_g = len(genomes)
            for _ in range(n_resamples):
                sampled_g = [genomes[rng.randrange(k_g)] for _ in range(k_g)]
                sampled_evals = [e for g in sampled_g for e in by_genome[g]]
                if sampled_evals:
                    boot_tot.append(
                        sum(1 for e in sampled_evals if e.recovered) / float(len(sampled_evals))
                    )
                s_on = [e for e in sampled_evals if e.stratum == "on_growth"]
                if s_on:
                    boot_on.append(sum(1 for e in s_on if e.recovered) / float(len(s_on)))
                s_off = [e for e in sampled_evals if e.stratum == "off_growth"]
                if s_off:
                    boot_off.append(sum(1 for e in s_off if e.recovered) / float(len(s_off)))

                s_p = [e for e in sampled_evals if getattr(e, "case_source", "price148") == "price148"]
                if s_p:
                    boot_p_tot.append(sum(1 for e in s_p if e.recovered) / float(len(s_p)))
                s_p_off = [e for e in s_p if e.stratum == "off_growth"]
                if s_p_off:
                    boot_p_off.append(sum(1 for e in s_p_off if e.recovered) / float(len(s_p_off)))

                s_sp = [e for e in sampled_evals if getattr(e, "case_source", "price148") == "swissprot"]
                if s_sp:
                    boot_sp_tot.append(sum(1 for e in s_sp if e.recovered) / float(len(s_sp)))
                s_sp_off = [e for e in s_sp if e.stratum == "off_growth"]
                if s_sp_off:
                    boot_sp_off.append(sum(1 for e in s_sp_off if e.recovered) / float(len(s_sp_off)))

                s_harm = sum(e.n_harmful_calls for e in sampled_evals)
                s_known = s_harm + sum(e.n_supported_calls for e in sampled_evals)
                boot_harm.append((s_harm / float(s_known)) if s_known > 0 else 0.0)

        ci_p_tot = _percentile_ci(boot_p_tot) if boot_p_tot else (0.0, 0.0)
        ci_p_off = _percentile_ci(boot_p_off) if boot_p_off else (0.0, 0.0)
        ci_sp_tot = _percentile_ci(boot_sp_tot) if boot_sp_tot else (0.0, 0.0)
        ci_sp_off = _percentile_ci(boot_sp_off) if boot_sp_off else (0.0, 0.0)

        tot_lps = sum(e.n_lps for e in evals)
        tot_gfs = sum(e.n_gapfills for e in evals)
        tot_cpu = sum(e.cpu_seconds for e in evals)
        tot_wall = sum(e.wall_seconds for e in evals)
        max_rss = max((e.peak_rss_kb for e in evals), default=0)

        raw_metrics[arm] = {
            "n_cases_total": n_total,
            "n_recovered_total": n_rec_total,
            "recovery_rate_total": (n_rec_total / float(n_total)) if n_total > 0 else 0.0,
            "recovery_ci_total": _percentile_ci(boot_tot),
            "n_cases_on_growth": n_on,
            "n_recovered_on_growth": n_rec_on,
            "recovery_rate_on_growth": (n_rec_on / float(n_on)) if n_on > 0 else 0.0,
            "recovery_ci_on_growth": _percentile_ci(boot_on),
            "n_cases_off_growth": n_off,
            "n_recovered_off_growth": n_rec_off,
            "recovery_rate_off_growth": (n_rec_off / float(n_off)) if n_off > 0 else 0.0,
            "recovery_ci_off_growth": _percentile_ci(boot_off),
            "n_cases_price148": n_p,
            "n_recovered_price148": n_rec_p,
            "recovery_rate_price148": (n_rec_p / float(n_p)) if n_p > 0 else 0.0,
            "recovery_ci_price148": ci_p_tot,
            "n_cases_price148_off_growth": n_p_off,
            "n_recovered_price148_off_growth": n_rec_p_off,
            "recovery_rate_price148_off_growth": (n_rec_p_off / float(n_p_off)) if n_p_off > 0 else 0.0,
            "recovery_ci_price148_off_growth": ci_p_off,
            "n_cases_swissprot": n_sp,
            "n_recovered_swissprot": n_rec_sp,
            "recovery_rate_swissprot": (n_rec_sp / float(n_sp)) if n_sp > 0 else 0.0,
            "recovery_ci_swissprot": ci_sp_tot,
            "n_cases_swissprot_off_growth": n_sp_off,
            "n_recovered_swissprot_off_growth": n_rec_sp_off,
            "recovery_rate_swissprot_off_growth": (n_rec_sp_off / float(n_sp_off)) if n_sp_off > 0 else 0.0,
            "recovery_ci_swissprot_off_growth": ci_sp_off,
            "n_gene_calls": n_gene,
            "n_reaction_only_calls": n_rxn_only,
            "n_supported_calls": n_supp,
            "n_harmful_calls": n_harm,
            "n_unconfirmed_calls": n_unconf,
            "harm_rate": harm_rate,
            "harm_ci": _percentile_ci(boot_harm),
            "total_tests_considered": sum(e.tests_considered for e in evals),
            "total_tests_executed": sum(e.tests_executed for e in evals),
            "total_lps": tot_lps,
            "total_gapfills": tot_gfs,
            "total_cpu_seconds": tot_cpu,
            "total_wall_seconds": tot_wall,
            "max_peak_rss_kb": max_rss,
            "lps_per_recovered_enzyme": (tot_lps / float(n_rec_total)) if n_rec_total > 0 else None,
        }

    a0_harm = float(raw_metrics.get("A0", {}).get("harm_rate", 0.0))
    summaries: Dict[str, ArmSummaryMetrics] = {}
    for arm, m in raw_metrics.items():
        harm_val = float(m["harm_rate"])
        summaries[arm] = ArmSummaryMetrics(
            arm=arm,
            exceeds_a0_harm_gate=(harm_val - a0_harm) > HARM_EXCESS_LIMIT,
            **m,  # type: ignore[arg-type]
        )
    return summaries


def run_preregistered_tests(
    evaluations_by_arm: Mapping[str, Sequence[CaseArmEvaluation]],
    n_resamples: int = BOOTSTRAP_RESAMPLES,
    seed: int = BOOTSTRAP_SEED,
) -> Tuple[McNemarResult, Tuple[McNemarResult, ...], H3EfficiencyResult]:
    """Run pre-registered statistical tests H1, H2 (Holm-Bonferroni), and H3 (§7)."""

    def paired_counts(
        treat_evals: Sequence[CaseArmEvaluation],
        base_evals: Sequence[CaseArmEvaluation],
        stratum: str = "off_growth",
    ) -> Tuple[int, int, int, int, int]:
        base_by_case = {e.case_id: e for e in base_evals if e.stratum == stratum}
        treat_by_case = {e.case_id: e for e in treat_evals if e.stratum == stratum}
        shared_ids = sorted(set(base_by_case) & set(treat_by_case))
        both = 0
        b_treat_only = 0
        c_base_only = 0
        neither = 0
        for cid in shared_ids:
            t_rec = treat_by_case[cid].recovered
            b_rec = base_by_case[cid].recovered
            if t_rec and b_rec:
                both += 1
            elif t_rec and not b_rec:
                b_treat_only += 1
            elif not t_rec and b_rec:
                c_base_only += 1
            else:
                neither += 1
        return len(shared_ids), both, b_treat_only, c_base_only, neither

    a0_evals = evaluations_by_arm.get("A0", ())
    a123_evals = evaluations_by_arm.get("A1_A2_A3", ())

    # H1: A1 U A2 U A3 vs A0 on off_growth
    n_h1, both_h1, b_h1, c_h1, neith_h1 = paired_counts(a123_evals, a0_evals, "off_growth")
    p_h1 = exact_one_sided_mcnemar_p(b_h1, c_h1)
    h1_result = McNemarResult(
        hypothesis_id="H1",
        arm_treatment="A1_A2_A3",
        arm_baseline="A0",
        stratum="off_growth",
        n_pairs=n_h1,
        both_recovered=both_h1,
        treatment_only=b_h1,
        baseline_only=c_h1,
        neither_recovered=neith_h1,
        p_value_raw=p_h1,
        p_value_holm=p_h1,
        significant_at_0_05=p_h1 < 0.05,
    )

    # H2: A1, A2, A3 each vs A0 on off_growth (Holm-Bonferroni corrected across 3 arms)
    h2_arms = ("A1", "A2", "A3")
    h2_raw_stats: List[Tuple[str, int, int, int, int, int, float]] = []
    for arm_name in h2_arms:
        arm_evals = evaluations_by_arm.get(arm_name, ())
        npairs, both, b_cnt, c_cnt, neith = paired_counts(arm_evals, a0_evals, "off_growth")
        pval = exact_one_sided_mcnemar_p(b_cnt, c_cnt)
        h2_raw_stats.append((arm_name, npairs, both, b_cnt, c_cnt, neith, pval))

    holm_pvals = holm_bonferroni_adjust([item[6] for item in h2_raw_stats])
    h2_results: List[McNemarResult] = []
    for idx, (arm_name, npairs, both, b_cnt, c_cnt, neith, pval) in enumerate(h2_raw_stats):
        padj = holm_pvals[idx]
        h2_results.append(
            McNemarResult(
                hypothesis_id=f"H2_{arm_name}",
                arm_treatment=arm_name,
                arm_baseline="A0",
                stratum="off_growth",
                n_pairs=npairs,
                both_recovered=both,
                treatment_only=b_cnt,
                baseline_only=c_cnt,
                neither_recovered=neith,
                p_value_raw=pval,
                p_value_holm=padj,
                significant_at_0_05=padj < 0.05,
            )
        )

    # H3: A4@25% off-growth retention relative to A1_A2_A3
    a4_25_evals = evaluations_by_arm.get("A4@25%", ())
    off_rec_budget = sum(1 for e in a4_25_evals if e.stratum == "off_growth" and e.recovered)
    off_rec_full = sum(1 for e in a123_evals if e.stratum == "off_growth" and e.recovered)
    retention = (off_rec_budget / float(off_rec_full)) if off_rec_full > 0 else 0.0

    # Genome-block bootstrap for H3 retention ratio
    full_by_genome: Dict[str, List[CaseArmEvaluation]] = {}
    budg_by_genome: Dict[str, List[CaseArmEvaluation]] = {}
    for e in a123_evals:
        full_by_genome.setdefault(e.genome_id, []).append(e)
    for e in a4_25_evals:
        budg_by_genome.setdefault(e.genome_id, []).append(e)
    genomes = sorted(set(full_by_genome) & set(budg_by_genome))
    rng = random.Random(seed)
    boot_ret: List[float] = []
    if genomes and n_resamples > 0:
        kg = len(genomes)
        for _ in range(n_resamples):
            sampled_g = [genomes[rng.randrange(kg)] for _ in range(kg)]
            f_off = sum(
                1
                for g in sampled_g
                for e in full_by_genome[g]
                if e.stratum == "off_growth" and e.recovered
            )
            b_off = sum(
                1
                for g in sampled_g
                for e in budg_by_genome[g]
                if e.stratum == "off_growth" and e.recovered
            )
            if f_off > 0:
                boot_ret.append(b_off / float(f_off))

    gf_budget = sum(e.n_gapfills for e in a4_25_evals)
    gf_full = sum(e.n_gapfills for e in a123_evals)
    exec_budget = sum(e.tests_executed for e in a4_25_evals)
    cons_full = sum(e.tests_considered for e in a123_evals)
    test_frac = (exec_budget / float(cons_full)) if cons_full > 0 else 0.0

    h3_result = H3EfficiencyResult(
        arm_budget="A4@25%",
        arm_full="A1_A2_A3",
        off_growth_recovered_budget=off_rec_budget,
        off_growth_recovered_full=off_rec_full,
        retention_ratio=retention,
        retention_ci=_percentile_ci(boot_ret),
        gapfills_budget=gf_budget,
        gapfills_full=gf_full,
        tests_executed_budget=exec_budget,
        tests_considered_full=cons_full,
        test_budget_fraction=test_frac,
        passes_80pct_retention=retention >= H3_RETENTION_TARGET,
    )

    return h1_result, tuple(h2_results), h3_result


def permute_candidates_within_genome(
    candidates: Sequence[CandidateClaim],
    seed: int,
    extra_proteins: Sequence[str] = (),
) -> Tuple[CandidateClaim, ...]:
    """Permute protein IDs of a genome's GSPA candidate table for a fixed seed (§8)."""
    distinct_proteins = sorted(
        {c.protein_id for c in candidates} | {p for p in extra_proteins if p}
    )
    if len(distinct_proteins) <= 1:
        return tuple(candidates)
    rng = random.Random(seed)
    shuffled = list(distinct_proteins)
    rng.shuffle(shuffled)
    mapping = dict(zip(distinct_proteins, shuffled))
    return tuple(
        CandidateClaim(
            protein_id=mapping.get(c.protein_id, c.protein_id),
            ec=c.ec,
            score=c.score,
        )
        for c in candidates
    )


def run_shuffled_candidate_control(
    cases: Sequence[TruthCase],
    strata_by_case: Mapping[str, str],
    arm_results_by_case: Mapping[str, Mapping[str, ArmRunResult]],
    unmasked_candidates_by_genome: Mapping[str, Sequence[CandidateClaim]],
    ec_to_reactions: Mapping[str, Sequence[str]],
    known_experimental_ecs_by_genome: Mapping[str, Mapping[str, Sequence[str]]],
    arm_name: str = "A1_A2_A3",
    seeds: Sequence[int] = SHUFFLE_CONTROL_SEEDS,
    extra_proteins_by_genome: Optional[Mapping[str, Sequence[str]]] = None,
) -> ShuffledControlSummary:
    """Rerun call assignment on `arm_name` added reactions across 5 shuffled seeds (§8)."""
    outcomes: List[ShuffledSeedOutcome] = []
    for seed in seeds:
        permuted_by_genome: Dict[str, Tuple[CandidateClaim, ...]] = {}
        for gid, cands in unmasked_candidates_by_genome.items():
            extras = (extra_proteins_by_genome or {}).get(gid, ())
            permuted_by_genome[gid] = permute_candidates_within_genome(
                cands, seed=seed, extra_proteins=extras
            )

        n_cases = 0
        n_rec = 0
        n_off = 0
        n_off_rec = 0

        for case in cases:
            arm_res = arm_results_by_case.get(case.case_id, {}).get(arm_name)
            if arm_res is None:
                continue
            perm_cands = permuted_by_genome.get(case.genome_id, ())
            known_ecs = known_experimental_ecs_by_genome.get(case.genome_id, {})
            stratum = strata_by_case.get(case.case_id, "off_growth")

            ev = evaluate_case_arm(
                case=case,
                stratum=stratum,
                arm_result=arm_res,
                call_candidates=perm_cands,
                ec_to_reactions=ec_to_reactions,
                known_experimental_ecs=known_ecs,
            )
            n_cases += 1
            if ev.recovered:
                n_rec += 1
            if stratum == "off_growth":
                n_off += 1
                if ev.recovered:
                    n_off_rec += 1

        outcomes.append(
            ShuffledSeedOutcome(
                seed=seed,
                arm=arm_name,
                n_cases=n_cases,
                n_recovered=n_rec,
                recovery_rate=(n_rec / float(n_cases)) if n_cases > 0 else 0.0,
                n_off_growth_cases=n_off,
                n_off_growth_recovered=n_off_rec,
                off_growth_recovery_rate=(n_off_rec / float(n_off)) if n_off > 0 else 0.0,
            )
        )

    mean_rec = (
        sum(o.recovery_rate for o in outcomes) / float(len(outcomes))
        if outcomes
        else 0.0
    )
    max_rec = max((o.recovery_rate for o in outcomes), default=0.0)
    mean_off_rec = (
        sum(o.off_growth_recovery_rate for o in outcomes) / float(len(outcomes))
        if outcomes
        else 0.0
    )
    return ShuffledControlSummary(
        arm=arm_name,
        seeds=tuple(seeds),
        per_seed=tuple(outcomes),
        mean_recovery_rate=mean_rec,
        max_recovery_rate=max_rec,
        mean_off_growth_recovery_rate=mean_off_rec,
    )


def summarize_solver_replays(
    replays: Sequence[SolverReplayRecord],
) -> SolverReplaySummary:
    if not replays:
        return SolverReplaySummary(
            n_replayed=0,
            n_agreed=0,
            agreement_rate=1.0,
            max_relative_obj_diff=0.0,
        )
    n_rep = len(replays)
    n_ok = sum(1 for r in replays if r.agreed)
    max_diff = max(r.relative_obj_diff for r in replays)
    return SolverReplaySummary(
        n_replayed=n_rep,
        n_agreed=n_ok,
        agreement_rate=n_ok / float(n_rep),
        max_relative_obj_diff=max_diff,
    )
