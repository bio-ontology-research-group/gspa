#!/usr/bin/env python3
"""Simulation engine for off-growth round 2 (plan §4, steps 0-2)."""

from __future__ import annotations

import copy
import csv
from dataclasses import dataclass
import multiprocessing as mp
from pathlib import Path
import re
import sys
import time
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple, Union

import cobra

V3_DIR = Path(__file__).resolve().parents[1]
if str(V3_DIR) not in sys.path:
    sys.path.insert(0, str(V3_DIR))

from offgrowth.arms import (  # noqa: E402
    DEFAULT_FLUX_TOL,
    DEFAULT_PANEL_UPTAKE,
    GapFillTestResult,
    WeightedGapFiller,
    _ensure_exchange_reactions,
    _find_biomass_rxn_id,
    _init_single_thread_worker,
    apply_medium,
    build_sole_source_medium,
    build_union_medium,
    compute_reaction_weights,
    identify_dead_end_targets,
    map_reactions_to_candidate_genes,
    peak_rss_kb,
)
from offgrowth.mask import (  # noqa: E402
    CandidateClaim,
    MaskingAuditError,
    MaskingAuditReport,
    ReactionEvidenceHit,
    TransporterHit,
    TruthCase,
    assert_masking_clean,
    audit_case_masking,
    base_seed_rxn_id,
    load_extra_cases,
    load_price148_cases,
    mask_candidates,
    mask_cobra_model,
    mask_reaction_evidence,
    mask_transporters,
    resolve_equivalent_ecs_and_reactions,
)
from offgrowth.panel import (  # noqa: E402
    PanelEntry,
    SeedMetabolite,
    SeedReaction,
    is_complete_ec,
    normalize_ec,
)
from offgrowth.score import (  # noqa: E402
    StratumClassification,
    classify_case_stratum,
)

from offgrowth2.medium import (  # noqa: E402
    minimal_medium,
    sole_source_medium_from_minimal,
)


@dataclass(frozen=True)
class CaseRecord:
    """Wrapper around `TruthCase` carrying the `hedged` flag (0 or 1) without altering `TruthCase`."""

    truth_case: TruthCase
    hedged: int = 0

    @property
    def case_id(self) -> str:
        return self.truth_case.case_id

    @property
    def genome_id(self) -> str:
        return self.truth_case.genome_id

    @property
    def protein_id(self) -> str:
        return self.truth_case.protein_id

    @property
    def truth_ecs(self) -> Tuple[str, ...]:
        return self.truth_case.truth_ecs

    @property
    def case_source(self) -> str:
        return self.truth_case.case_source

    def to_schema_dict(self) -> Dict[str, object]:
        return {
            "case_id": self.truth_case.case_id,
            "protein_id": self.truth_case.protein_id,
            "truth_ecs": sorted(self.truth_case.truth_ecs),
            "case_source": self.truth_case.case_source,
            "hedged": int(bool(self.hedged)),
        }


@dataclass(frozen=True)
class BaselineFill:
    """Result of Step 1: one weighted gap-fill of the draft for biomass on the minimal medium."""

    fill: Tuple[str, ...]
    replay_fill: Optional[Tuple[str, ...]]
    stable_fill: Tuple[str, ...]
    n_lps: int
    cpu_seconds: float
    wall_seconds: float = 0.0
    grows_without_fill: bool = False
    feasible_in_universal: bool = True

    @property
    def r_base(self) -> Tuple[str, ...]:
        return self.fill

    def to_schema_dict(self) -> Dict[str, object]:
        return {
            "fill": list(self.fill),
            "replay_fill": list(self.replay_fill) if self.replay_fill is not None else None,
            "stable_fill": list(self.stable_fill),
            "n_lps": int(self.n_lps),
            "cpu_seconds": float(self.cpu_seconds),
            "grows_without_fill": bool(self.grows_without_fill),
        }


@dataclass(frozen=True)
class TestRecord:
    """Single test outcome in the Step 2 superset (`T1`, `T2`, `T3`, `A1old`, `A0`)."""

    test_id: str
    test_type: str
    target_id: str
    panel_id: Optional[str]
    element: Optional[str]
    gapfill_ran: bool
    feasible_in_universal: bool
    already_satisfied: bool
    fill: Tuple[str, ...]
    replay_fill: Optional[Tuple[str, ...]]
    stable_fill: Tuple[str, ...]
    n_lps: int
    cpu_seconds: float

    def to_schema_dict(self) -> Dict[str, object]:
        return {
            "test_id": self.test_id,
            "test_type": self.test_type,
            "target_id": self.target_id,
            "panel_id": self.panel_id,
            "element": self.element,
            "gapfill_ran": bool(self.gapfill_ran),
            "feasible_in_universal": bool(self.feasible_in_universal),
            "already_satisfied": bool(self.already_satisfied),
            "fill": list(self.fill),
            "replay_fill": list(self.replay_fill) if self.replay_fill is not None else None,
            "stable_fill": list(self.stable_fill),
            "n_lps": int(self.n_lps),
            "cpu_seconds": float(self.cpu_seconds),
        }


@dataclass(frozen=True)
class RunSupersetResult:
    """Complete execution record for one masked case or the unmasked deployment run."""

    run_kind: str  # "masked" or "unmasked"
    case: Optional[CaseRecord]
    audit: Optional[MaskingAuditReport]
    stratum: Optional[StratumClassification]
    target_reactions: Tuple[str, ...]
    draft_reactions: Tuple[str, ...]
    baseline_fill: BaselineFill
    tests: Tuple[TestRecord, ...]
    n_lps: int
    cpu_seconds: float
    wall_seconds: float
    peak_rss_kb: int

    def to_schema_dict(self) -> Dict[str, object]:
        audit_dict: Optional[Dict[str, object]] = None
        if self.audit is not None:
            raw_audit = self.audit.to_dict()
            audit_dict = {
                "case_id": str(raw_audit["case_id"]),
                "protein_id": str(raw_audit["protein_id"]),
                "truth_ecs": sorted(str(x) for x in raw_audit.get("truth_ecs", ())),
                "equivalent_ecs": sorted(str(x) for x in raw_audit.get("equivalent_ecs", ())),
                "target_reactions": sorted(str(x) for x in raw_audit.get("target_reactions", ())),
                "passed": bool(raw_audit["passed"]),
                "violations": sorted(str(x) for x in raw_audit.get("violations", ())),
            }

        stratum_dict: Optional[Dict[str, object]] = None
        if self.stratum is not None:
            stratum_dict = {
                "stratum": self.stratum.stratum,
                "wt_biomass": float(self.stratum.wt_biomass),
                "ko_biomass": float(self.stratum.ko_biomass),
                "ko_fraction": float(self.stratum.ko_fraction),
            }

        return {
            "run_kind": self.run_kind,
            "case": self.case.to_schema_dict() if self.case is not None else None,
            "audit": audit_dict,
            "stratum": stratum_dict,
            "target_reactions": list(self.target_reactions),
            "draft_reactions": list(self.draft_reactions),
            "baseline_fill": self.baseline_fill.to_schema_dict(),
            "tests": [t.to_schema_dict() for t in self.tests],
            "n_lps": int(self.n_lps),
            "cpu_seconds": float(self.cpu_seconds),
            "wall_seconds": float(self.wall_seconds),
            "peak_rss_kb": int(self.peak_rss_kb),
        }


def load_label_cases(
    label_cases_tsv: Path,
    genome_id: Optional[str] = None,
    case_id_prefix: str = "FB:",
    case_source: str = "fitness_browser",
) -> Tuple[CaseRecord, ...]:
    """Load Fitness Browser label cases (one case per `(protein_id, function_id)` with `overlap == 'new'`)."""
    if not label_cases_tsv.exists():
        return ()

    by_pair: Dict[Tuple[str, str, str], int] = {}
    with label_cases_tsv.open("r", encoding="utf-8", errors="replace", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            gid = (row.get("genome_id") or row.get("host_assembly") or "").strip()
            if not gid:
                continue
            if genome_id is not None and gid != genome_id:
                continue
            overlap = (row.get("overlap") or "new").strip().lower()
            if overlap != "new":
                continue
            pid = (row.get("protein_id") or "").strip()
            fid = (row.get("function_id") or row.get("ec") or "").strip()
            if not pid or not fid:
                continue
            norm_ec = normalize_ec(fid)
            if not is_complete_ec(norm_ec):
                continue
            try:
                hedged = int(float((row.get("hedged") or "0").strip() or "0"))
            except ValueError:
                hedged = 0
            key = (gid, pid, norm_ec)
            if key not in by_pair or hedged > by_pair[key]:
                by_pair[key] = 1 if hedged else 0

    out: List[CaseRecord] = []
    for (gid, pid, norm_ec), hedged in sorted(by_pair.items()):
        base_cid = pid if pid.startswith(case_id_prefix) else f"{case_id_prefix}{pid}"
        cid = f"{base_cid}:{norm_ec}"
        tc = TruthCase(
            case_id=cid,
            genome_id=gid,
            protein_id=pid,
            truth_ecs=(norm_ec,),
            case_source=case_source,
        )
        out.append(CaseRecord(truth_case=tc, hedged=hedged))
    return tuple(out)


def load_cases_for_genome(
    genome_id: str,
    price_manifest: Optional[Path] = None,
    extra_cases_tsv: Optional[Path] = None,
    label_cases_tsv: Optional[Path] = None,
    label_cases_prefix: str = "FB:",
) -> Tuple[CaseRecord, ...]:
    """Load all requested truth cases for `genome_id` across Price-148, Swiss-Prot, and Fitness Browser."""
    records: List[CaseRecord] = []
    if price_manifest is not None and price_manifest.exists():
        for tc in load_price148_cases(price_manifest):
            if tc.genome_id == genome_id:
                records.append(CaseRecord(truth_case=tc, hedged=0))
    if extra_cases_tsv is not None and extra_cases_tsv.exists():
        for tc in load_extra_cases(
            extra_cases_tsv,
            genome_id=genome_id,
            case_id_prefix="SP:",
            case_source="swissprot",
        ):
            records.append(CaseRecord(truth_case=tc, hedged=0))
    if label_cases_tsv is not None and label_cases_tsv.exists():
        records.extend(
            load_label_cases(
                label_cases_tsv,
                genome_id=genome_id,
                case_id_prefix=label_cases_prefix,
                case_source="fitness_browser",
            )
        )
    return tuple(records)


def collect_known_labels_for_genome(
    genome_id: str,
    price_manifest: Optional[Path] = None,
    extra_cases_tsv: Optional[Path] = None,
    label_cases_tsv: Optional[Path] = None,
) -> Dict[str, List[str]]:
    """Collect `{protein_id: sorted_ecs}` from every label file given (L0 + L1 or L2/L3)."""
    by_prot: Dict[str, Set[str]] = {}
    if price_manifest is not None and price_manifest.exists():
        for tc in load_price148_cases(price_manifest):
            if tc.genome_id == genome_id:
                by_prot.setdefault(tc.protein_id, set()).update(tc.truth_ecs)

    for tsv_path, require_experimental in (
        (extra_cases_tsv, True),
        (label_cases_tsv, False),
    ):
        if tsv_path is None or not tsv_path.exists():
            continue
        with tsv_path.open("r", encoding="utf-8", errors="replace", newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                gid = (row.get("genome_id") or row.get("host_assembly") or "").strip()
                if gid and gid != genome_id:
                    continue
                if require_experimental:
                    tier = (row.get("tier") or "experimental").strip().lower()
                    if tier and tier not in {"experimental", "curated_phenotype"}:
                        continue
                pid = (row.get("protein_id") or "").strip()
                raw_ec = (row.get("function_id") or row.get("ec") or row.get("ground_truth_ecs") or "").strip()
                if not pid or not raw_ec:
                    continue
                for tok in re.split(r"[;,|\s]+", raw_ec):
                    norm_ec = normalize_ec(tok)
                    if is_complete_ec(norm_ec):
                        by_prot.setdefault(pid, set()).add(norm_ec)

    return {p: sorted(ecs) for p, ecs in sorted(by_prot.items())}


def draft_evidenced_base_reactions(model: cobra.Model) -> Tuple[str, ...]:
    """Return sorted base SEED reaction IDs present in `model` (excluding exchanges and biomass)."""
    biomass_ids = {"bio1", "R_bio1", "biomass"}
    try:
        biomass_ids.add(_find_biomass_rxn_id(model))
    except ValueError:
        pass
    rxns = {
        base_seed_rxn_id(r.id)
        for r in model.reactions
        if not r.id.startswith("EX_")
        and r.id not in biomass_ids
        and bool((r.gene_reaction_rule or "").strip())
    }
    return tuple(sorted(rxns))


def baseline_fill(
    filler: WeightedGapFiller,
    biomass_id: Optional[str] = None,
    minimal_med: Optional[Mapping[str, float]] = None,
) -> BaselineFill:
    """Step 1: one weighted gap-fill of the draft for biomass on the minimal medium (`force_replay=True`)."""
    eff_biomass_id = biomass_id or _find_biomass_rxn_id(filler.draft_model)
    eff_medium = minimal_med if minimal_med is not None else {}

    res = filler.solve_weighted_gapfill(
        test_id="baseline_minimal",
        arm="baseline",
        target_rxn_id=eff_biomass_id,
        medium=eff_medium,
        check_draft_first=True,
        force_replay=True,
    )
    primary_fill = tuple(sorted({base_seed_rxn_id(r) for r in res.added_reactions if r}))
    if res.replay is not None:
        rep_fill: Optional[Tuple[str, ...]] = tuple(
            sorted({base_seed_rxn_id(r) for r in res.replay.replay_added_rxns if r})
        )
        stable = tuple(sorted(set(primary_fill) & set(rep_fill)))
    else:
        rep_fill = None
        stable = primary_fill

    return BaselineFill(
        fill=primary_fill,
        replay_fill=rep_fill,
        stable_fill=stable,
        n_lps=res.n_lps,
        cpu_seconds=res.cpu_seconds,
        wall_seconds=res.wall_seconds,
        grows_without_fill=res.already_satisfied_in_draft,
        feasible_in_universal=res.feasible_in_universal,
    )


def build_base_model(
    draft_model: cobra.Model,
    universal_reactions: Mapping[str, SeedReaction],
    R_base: Iterable[str],
) -> cobra.Model:
    """Build $M_{\\text{base}} = D + R_{\\text{base}}$ from `draft_model` and `universal_reactions` with no genes.

    Reuses the metabolite and reaction construction logic of
    `WeightedGapFiller._build_universal_completed_model` (`benchmark/v3/offgrowth/arms.py:723-794`),
    adding each reaction in `R_base` with its ModelSEED stoichiometry and reversibility bounds.
    """
    base_model = copy.deepcopy(draft_model)
    existing_base_ids: Set[str] = {
        base_seed_rxn_id(r.id)
        for r in base_model.reactions
        if not r.id.startswith("EX_")
    }
    r_base_sorted = sorted({base_seed_rxn_id(r) for r in R_base if r})

    # Reused from benchmark/v3/offgrowth/arms.py:723-743
    existing_mets: Dict[str, cobra.Metabolite] = {m.id: m for m in base_model.metabolites}
    new_mets: List[cobra.Metabolite] = []
    for base_id in r_base_sorted:
        if base_id in existing_base_ids:
            continue
        urxn = universal_reactions.get(base_id)
        if urxn is None or not urxn.stoichiometry:
            continue
        for term in urxn.stoichiometry:
            if term.species_id not in existing_mets:
                met = cobra.Metabolite(
                    id=term.species_id,
                    name=term.species_id,
                    compartment=term.compartment,
                )
                existing_mets[term.species_id] = met
                new_mets.append(met)
    if new_mets:
        base_model.add_metabolites(new_mets)

    new_rxns: List[cobra.Reaction] = []
    for base_id in r_base_sorted:
        if base_id in existing_base_ids:
            continue
        urxn = universal_reactions.get(base_id)
        if urxn is None or not urxn.stoichiometry:
            continue
        lb, ub = urxn.bounds
        rxn_id = f"{base_id}_c0" if f"{base_id}_c0" not in base_model.reactions else base_id
        crxn = cobra.Reaction(
            id=rxn_id,
            name=urxn.name or base_id,
            lower_bound=lb,
            upper_bound=ub,
        )
        crxn.add_metabolites(
            {existing_mets[t.species_id]: t.coeff for t in urxn.stoichiometry}
        )
        crxn.gene_reaction_rule = ""
        new_rxns.append(crxn)
        existing_base_ids.add(base_id)

    if new_rxns:
        base_model.add_reactions(new_rxns)

    # Reused from benchmark/v3/offgrowth/arms.py:788-793
    ex_ids = [
        f"EX_{m.id}"
        for m in base_model.metabolites
        if m.id.endswith("_e0") and f"EX_{m.id}" not in base_model.reactions
    ]
    if ex_ids:
        _ensure_exchange_reactions(base_model, ex_ids)

    return base_model


def gapfill_result_to_test_record(
    res: GapFillTestResult,
    test_id: str,
    test_type: str,
    target_id: str,
    panel_id: Optional[str] = None,
    element: Optional[str] = None,
) -> TestRecord:
    """Convert a pilot `GapFillTestResult` into a round-2 `TestRecord`."""
    primary_fill = tuple(sorted({base_seed_rxn_id(r) for r in res.added_reactions if r}))
    if res.replay is not None:
        rep_fill: Optional[Tuple[str, ...]] = tuple(
            sorted({base_seed_rxn_id(r) for r in res.replay.replay_added_rxns if r})
        )
        stable = tuple(sorted(set(primary_fill) & set(rep_fill)))
    else:
        rep_fill = None
        stable = primary_fill

    return TestRecord(
        test_id=test_id,
        test_type=test_type,
        target_id=target_id,
        panel_id=panel_id,
        element=element,
        gapfill_ran=bool(res.gapfill_ran),
        feasible_in_universal=bool(res.feasible_in_universal),
        already_satisfied=bool(res.already_satisfied_in_draft),
        fill=primary_fill,
        replay_fill=rep_fill,
        stable_fill=stable,
        n_lps=int(res.n_lps),
        cpu_seconds=float(res.cpu_seconds),
    )


# ---------------------------------------------------------------------------
# Forked worker pool state and task functions (pattern from arms.py:1578-1680)
# ---------------------------------------------------------------------------

_WORKER_DRAFT_FILLER: Optional[WeightedGapFiller] = None
_WORKER_BASE_FILLER: Optional[WeightedGapFiller] = None
_WORKER_UNION_MED: Dict[str, float] = {}
_WORKER_FLUX_TOL: float = DEFAULT_FLUX_TOL


def _worker_check_blocked_v2(
    task: Tuple[int, str, bool],
) -> Tuple[int, str, bool, bool, int, float, float]:
    """Check whether `rxn_id` is blocked on `_WORKER_DRAFT_FILLER` under `_WORKER_UNION_MED`.

    Adapted from `benchmark/v3/offgrowth/arms.py:1596-1617`.
    """
    idx, rxn_id, is_rev = task
    assert _WORKER_DRAFT_FILLER is not None
    f_max, lps1, cpu1, wall1 = _WORKER_DRAFT_FILLER.check_draft_flux(
        rxn_id, _WORKER_UNION_MED, sense="max"
    )
    total_lps = lps1
    total_cpu = cpu1
    total_wall = wall1
    unblocked = f_max >= _WORKER_FLUX_TOL
    if not unblocked and is_rev:
        f_min, lps2, cpu2, wall2 = _WORKER_DRAFT_FILLER.check_draft_flux(
            rxn_id, _WORKER_UNION_MED, sense="min"
        )
        total_lps += lps2
        total_cpu += cpu2
        total_wall += wall2
        if f_min <= -_WORKER_FLUX_TOL:
            unblocked = True
    return idx, rxn_id, is_rev, unblocked, total_lps, total_cpu, total_wall


def _worker_solve_superset_task(
    task: Tuple[
        int,
        str,  # filler_kind: "base" or "draft"
        str,  # test_id (recorded in output)
        str,  # internal_test_id (for replay hash compatibility)
        str,  # test_type: "T1", "T2", "T3", "A1old", "A0"
        str,  # solve_target_rxn_id
        str,  # record_target_id
        Optional[str],  # panel_id
        Optional[str],  # element
        Dict[str, float],  # medium
        bool,  # check_draft_first
        bool,  # try_reverse_if_infeasible
        Optional[Tuple[str, Dict[str, float], float, float]],  # temp_spec
        bool,  # force_replay
        float,  # replay_fraction
    ],
) -> Tuple[int, TestRecord]:
    """Execute one gap-fill test in the Step 2 superset.

    Adapted from `benchmark/v3/offgrowth/arms.py:1620-1679`.
    """
    (
        task_idx,
        filler_kind,
        test_id,
        internal_test_id,
        test_type,
        solve_target_rxn_id,
        record_target_id,
        panel_id,
        element,
        med,
        check_draft_first,
        try_reverse_if_infeasible,
        temp_spec,
        force_replay,
        replay_fraction,
    ) = task

    filler = _WORKER_BASE_FILLER if filler_kind == "base" else _WORKER_DRAFT_FILLER
    assert filler is not None

    res = filler.solve_weighted_gapfill(
        test_id=internal_test_id,
        arm=test_type,
        target_rxn_id=solve_target_rxn_id,
        medium=med,
        check_draft_first=check_draft_first,
        target_direction="forward",
        temporary_rxn_spec=temp_spec,
        replay_fraction=replay_fraction,
        force_replay=force_replay,
    )
    # Reverse fallback for reversible blocked reactions in T2 (benchmark/v3/offgrowth/arms.py:1653-1678)
    if test_type == "T2" and not res.feasible_in_universal and try_reverse_if_infeasible:
        res_rev = filler.solve_weighted_gapfill(
            test_id=internal_test_id,
            arm="T2",
            target_rxn_id=solve_target_rxn_id,
            medium=med,
            check_draft_first=False,
            target_direction="reverse",
            replay_fraction=replay_fraction,
            force_replay=force_replay,
        )
        res = GapFillTestResult(
            test_id=internal_test_id,
            arm="T2",
            target_id=solve_target_rxn_id,
            feasible_in_universal=res_rev.feasible_in_universal,
            already_satisfied_in_draft=False,
            gapfill_ran=res_rev.gapfill_ran,
            vmax_universal=res_rev.vmax_universal,
            objective_value=res_rev.objective_value,
            added_reactions=res_rev.added_reactions,
            n_lps=res.n_lps + res_rev.n_lps,
            n_gapfills=res.n_gapfills + res_rev.n_gapfills,
            cpu_seconds=res.cpu_seconds + res_rev.cpu_seconds,
            wall_seconds=res.wall_seconds + res_rev.wall_seconds,
            replay=res_rev.replay,
        )

    rec = gapfill_result_to_test_record(
        res=res,
        test_id=test_id,
        test_type=test_type,
        target_id=record_target_id,
        panel_id=panel_id,
        element=element,
    )
    return task_idx, rec


def run_superset_for_model(
    draft_model: cobra.Model,
    universal_reactions: Mapping[str, SeedReaction],
    metabolites: Mapping[str, SeedMetabolite],
    candidates: Sequence[CandidateClaim],
    ec_to_reactions: Mapping[str, Sequence[str]],
    panel_entries: Sequence[PanelEntry],
    default_medium: Mapping[str, float],
    minimal_med: Optional[Mapping[str, float]] = None,
    case: Optional[CaseRecord] = None,
    audit_report: Optional[MaskingAuditReport] = None,
    stratum_info: Optional[StratumClassification] = None,
    target_reactions: Sequence[str] = (),
    allow_unmasked: bool = False,
    replay_fraction: float = 0.05,
    flux_tol: float = DEFAULT_FLUX_TOL,
    uptake_rate: float = DEFAULT_PANEL_UPTAKE,
    workers: int = 1,
) -> RunSupersetResult:
    """Run Step 1 (baseline minimal fill) and Step 2 (test superset `T1`, `T2`, `T3`, `A1old`, `A0`).

    Enforces `assert_masking_clean(audit_report)` before any LP solve unless `allow_unmasked=True`.
    Uses one `WeightedGapFiller` on the draft (`D`) and one on `M_base = D + R_base`.
    """
    t_cpu0 = time.process_time()
    t_wall0 = time.monotonic()

    if not allow_unmasked:
        if audit_report is None:
            raise MaskingAuditError("Masking audit report is required when allow_unmasked is False.")
        assert_masking_clean(audit_report)

    eff_minimal_med = (
        minimal_med
        if minimal_med is not None
        else minimal_medium(default_medium, metabolites)
    )
    biomass_id = _find_biomass_rxn_id(draft_model)
    union_med = build_union_medium(default_medium, panel_entries, uptake_rate=uptake_rate)

    rxn_to_candidates = map_reactions_to_candidate_genes(candidates, ec_to_reactions)
    draft_base_ids: Set[str] = {
        base_seed_rxn_id(r.id)
        for r in draft_model.reactions
        if not r.id.startswith("EX_")
    }
    draft_rxn_weights = compute_reaction_weights(
        universal_reactions, draft_base_ids, rxn_to_candidates
    )

    # Filler 1: masked (or unmasked) draft model D
    draft_filler = WeightedGapFiller(
        draft_model=draft_model,
        universal_reactions=universal_reactions,
        rxn_weights=draft_rxn_weights,
        flux_tol=flux_tol,
    )

    # Step 1: baseline minimal fill on D
    base_fill = baseline_fill(draft_filler, biomass_id=biomass_id, minimal_med=eff_minimal_med)

    # Build M_base = D + R_base and recompute weights for draft_base_ids U R_base
    m_base = build_base_model(draft_model, universal_reactions, base_fill.fill)
    base_ids_union = draft_base_ids | set(base_fill.fill)
    base_rxn_weights = compute_reaction_weights(
        universal_reactions, base_ids_union, rxn_to_candidates
    )

    # Filler 2: M_base for differential sole-source T1 tests
    base_filler = WeightedGapFiller(
        draft_model=m_base,
        universal_reactions=universal_reactions,
        rxn_weights=base_rxn_weights,
        flux_tol=flux_tol,
    )

    # Pre-warm exchange reactions on both fillers before forking workers
    all_panel_exchanges = {e.exchange_rxn for e in panel_entries}
    warm_union = dict(union_med)
    for ex_id in all_panel_exchanges:
        warm_union.setdefault(ex_id, 0.0)
    apply_medium(draft_filler.draft_model, warm_union)
    apply_medium(draft_filler.vmax_model, warm_union)
    apply_medium(draft_filler.universal_model, warm_union)

    warm_min = dict(eff_minimal_med)
    for ex_id in all_panel_exchanges:
        warm_min.setdefault(ex_id, 0.0)
    apply_medium(base_filler.draft_model, warm_min)
    apply_medium(base_filler.vmax_model, warm_min)
    apply_medium(base_filler.universal_model, warm_min)

    # Re-apply union_med on draft_filler for FVA
    apply_medium(draft_filler.draft_model, union_med)
    apply_medium(draft_filler.vmax_model, union_med)
    apply_medium(draft_filler.universal_model, union_med)

    # Identify T2 blocked-reaction FVA tasks (reused from benchmark/v3/offgrowth/arms.py:1777-1826)
    evidenced_rxns = sorted(
        [
            r
            for r in draft_model.reactions
            if not r.id.startswith("EX_")
            and r.id != biomass_id
            and bool((r.gene_reaction_rule or "").strip())
        ],
        key=lambda r: r.id,
    )
    fva_tasks: List[Tuple[int, str, bool]] = [
        (idx, r.id, r.lower_bound < 0.0) for idx, r in enumerate(evidenced_rxns)
    ]

    global _WORKER_DRAFT_FILLER, _WORKER_BASE_FILLER, _WORKER_UNION_MED, _WORKER_FLUX_TOL
    _WORKER_DRAFT_FILLER = draft_filler
    _WORKER_BASE_FILLER = base_filler
    _WORKER_UNION_MED = dict(union_med)
    _WORKER_FLUX_TOL = flux_tol

    use_fork_pool = (
        workers > 1
        and "fork" in mp.get_all_start_methods()
    )

    fva_lps = 0
    fva_cpu = 0.0
    blocked_t2_rxns: List[Tuple[str, bool]] = []

    try:
        if use_fork_pool and len(fva_tasks) > 1:
            ctx = mp.get_context("fork")
            n_pool = min(workers, len(fva_tasks))
            with ctx.Pool(processes=n_pool, initializer=_init_single_thread_worker) as pool:
                fva_out = list(pool.map(_worker_check_blocked_v2, fva_tasks, chunksize=1))
            fva_out.sort(key=lambda item: item[0])
        else:
            fva_out = [_worker_check_blocked_v2(t) for t in fva_tasks]

        for _, rxn_id, is_rev, unblocked, lps, cpu_s, _ in fva_out:
            fva_lps += lps
            fva_cpu += cpu_s
            if not unblocked:
                blocked_t2_rxns.append((rxn_id, is_rev))

        # Identify T3 dead-end targets on draft_model (benchmark/v3/offgrowth/arms.py:1833)
        dead_ends = identify_dead_end_targets(draft_model)

        # Assemble unified task list in deterministic order: T1, T2, T3, A1old, A0
        solve_tasks: List[
            Tuple[
                int,
                str,
                str,
                str,
                str,
                str,
                str,
                Optional[str],
                Optional[str],
                Dict[str, float],
                bool,
                bool,
                Optional[Tuple[str, Dict[str, float], float, float]],
                bool,
                float,
            ]
        ] = []
        task_idx = 0

        # 1. T1 differential sole-source tests on M_base from minimal_medium (force_replay=True)
        for entry in panel_entries:
            test_id = f"T1_{entry.panel_id}"
            sole_min_med = sole_source_medium_from_minimal(
                eff_minimal_med, entry, uptake_rate, metabolites=metabolites
            )
            solve_tasks.append(
                (
                    task_idx,
                    "base",
                    test_id,
                    test_id,
                    "T1",
                    biomass_id,
                    biomass_id,
                    entry.panel_id,
                    entry.element,
                    sole_min_med,
                    True,
                    False,
                    None,
                    True,
                    replay_fraction,
                )
            )
            task_idx += 1

        # 2. T2 unblock tests on draft_model under union_med (force_replay=True)
        for rxn_id, is_rev in blocked_t2_rxns:
            test_id = f"T2_unblock_{rxn_id}"
            solve_tasks.append(
                (
                    task_idx,
                    "draft",
                    test_id,
                    test_id,
                    "T2",
                    rxn_id,
                    rxn_id,
                    None,
                    None,
                    dict(union_med),
                    False,
                    is_rev,
                    None,
                    True,
                    replay_fraction,
                )
            )
            task_idx += 1

        # 3. T3 dead-end tests on draft_model under union_med (force_replay=True)
        # Reused from benchmark/v3/offgrowth/arms.py:1871-1880
        for met_id, kind in dead_ends:
            test_id = f"T3_deadend_{met_id}_{kind}"
            record_target_id = f"{met_id}:{kind}"
            if kind == "unproduced":
                temp_spec = (f"DM_{met_id}", {met_id: -1.0}, 0.0, 1000.0)
            else:
                temp_spec = (f"SRC_{met_id}", {met_id: 1.0}, 0.0, 1000.0)
            solve_tasks.append(
                (
                    task_idx,
                    "draft",
                    test_id,
                    test_id,
                    "T3",
                    temp_spec[0],
                    record_target_id,
                    None,
                    None,
                    dict(union_med),
                    False,
                    False,
                    temp_spec,
                    True,
                    replay_fraction,
                )
            )
            task_idx += 1

        # 4. A1-old sole-source tests on draft_model from default_medium (force_replay=False)
        for entry in panel_entries:
            test_id = f"A1old_{entry.panel_id}"
            pilot_test_id = f"A1_{entry.panel_id}"
            sole_def_med = build_sole_source_medium(
                default_medium, entry, metabolites, uptake_rate=uptake_rate
            )
            solve_tasks.append(
                (
                    task_idx,
                    "draft",
                    test_id,
                    pilot_test_id,
                    "A1old",
                    biomass_id,
                    biomass_id,
                    entry.panel_id,
                    entry.element,
                    sole_def_med,
                    True,
                    False,
                    None,
                    False,
                    replay_fraction,
                )
            )
            task_idx += 1

        # 5. A0 biomass fill on default_medium (force_replay=False)
        solve_tasks.append(
            (
                task_idx,
                "draft",
                "A0_biomass_default",
                "A0_biomass_default",
                "A0",
                biomass_id,
                biomass_id,
                None,
                None,
                dict(default_medium),
                True,
                False,
                None,
                False,
                replay_fraction,
            )
        )
        task_idx += 1

        if use_fork_pool and len(solve_tasks) > 1:
            ctx = mp.get_context("fork")
            n_pool = min(workers, len(solve_tasks))
            with ctx.Pool(processes=n_pool, initializer=_init_single_thread_worker) as pool:
                solved_out = list(pool.map(_worker_solve_superset_task, solve_tasks, chunksize=1))
            solved_out.sort(key=lambda item: item[0])
        else:
            solved_out = [_worker_solve_superset_task(t) for t in solve_tasks]
    finally:
        _WORKER_DRAFT_FILLER = None
        _WORKER_BASE_FILLER = None

    test_records = tuple(rec for _, rec in solved_out)
    total_lps = int(base_fill.n_lps) + int(fva_lps) + sum(t.n_lps for t in test_records)
    total_cpu = (
        float(base_fill.cpu_seconds)
        + float(fva_cpu)
        + sum(t.cpu_seconds for t in test_records)
        + (time.process_time() - t_cpu0)
    )
    total_wall = time.monotonic() - t_wall0

    return RunSupersetResult(
        run_kind="unmasked" if allow_unmasked else "masked",
        case=case,
        audit=audit_report,
        stratum=stratum_info,
        target_reactions=tuple(sorted({base_seed_rxn_id(r) for r in target_reactions if r})),
        draft_reactions=draft_evidenced_base_reactions(draft_model),
        baseline_fill=base_fill,
        tests=test_records,
        n_lps=total_lps,
        cpu_seconds=total_cpu,
        wall_seconds=total_wall,
        peak_rss_kb=peak_rss_kb(),
    )


def execute_masked_case_superset(
    case_record: Union[CaseRecord, TruthCase],
    draft_model: cobra.Model,
    unmasked_filled_model: cobra.Model,
    unmasked_candidates: Sequence[CandidateClaim],
    universal_reactions: Mapping[str, SeedReaction],
    metabolites: Mapping[str, SeedMetabolite],
    ec_to_reactions: Mapping[str, Sequence[str]],
    panel_entries: Sequence[PanelEntry],
    default_medium: Mapping[str, float],
    minimal_med: Optional[Mapping[str, float]] = None,
    transporters: Sequence[TransporterHit] = (),
    reaction_evidence: Sequence[ReactionEvidenceHit] = (),
    replay_fraction: float = 0.05,
    flux_tol: float = DEFAULT_FLUX_TOL,
    uptake_rate: float = DEFAULT_PANEL_UPTAKE,
    workers: int = 1,
) -> Tuple[Optional[RunSupersetResult], MaskingAuditReport]:
    """Mask `case_record`, enforce the mechanical audit before any LP, and run the Step 1-2 superset.

    If the mechanical masking audit fails, returns `(None, failed_audit_report)` without running any LP.
    """
    wrapped = (
        case_record
        if isinstance(case_record, CaseRecord)
        else CaseRecord(truth_case=case_record, hedged=0)
    )
    tc = wrapped.truth_case

    equiv_ecs, target_rxns = resolve_equivalent_ecs_and_reactions(
        tc.truth_ecs, ec_to_reactions
    )
    masked_model, _, _ = mask_cobra_model(draft_model, tc.protein_id, target_rxns)
    masked_cands = mask_candidates(unmasked_candidates, tc.protein_id, equiv_ecs)
    masked_trans = mask_transporters(transporters, tc.protein_id, target_rxns)
    masked_ev = mask_reaction_evidence(
        reaction_evidence, tc.protein_id, equiv_ecs, target_rxns
    )
    audit = audit_case_masking(
        case=tc,
        masked_model=masked_model,
        masked_candidates=masked_cands,
        equivalent_ecs=equiv_ecs,
        target_reactions=target_rxns,
        masked_transporters=masked_trans,
        masked_reaction_evidence=masked_ev,
    )
    try:
        assert_masking_clean(audit)
    except MaskingAuditError:
        return None, audit

    stratum_info = classify_case_stratum(
        case=tc,
        unmasked_filled_model=unmasked_filled_model,
        ec_to_reactions=ec_to_reactions,
        default_medium=default_medium,
        flux_tol=flux_tol,
    )

    run_res = run_superset_for_model(
        draft_model=masked_model,
        universal_reactions=universal_reactions,
        metabolites=metabolites,
        candidates=masked_cands,
        ec_to_reactions=ec_to_reactions,
        panel_entries=panel_entries,
        default_medium=default_medium,
        minimal_med=minimal_med,
        case=wrapped,
        audit_report=audit,
        stratum_info=stratum_info,
        target_reactions=target_rxns,
        allow_unmasked=False,
        replay_fraction=replay_fraction,
        flux_tol=flux_tol,
        uptake_rate=uptake_rate,
        workers=workers,
    )
    return run_res, audit
