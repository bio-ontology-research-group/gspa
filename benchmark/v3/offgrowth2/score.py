"""Scoring, arms, controls, comparator Q, bootstrap, and McNemar tests (plan §6–§8, spec §5).

Implements:
- All evaluation arms on the stored test superset:
  ``S0``, ``S1``, ``S1b``, ``N1``, ``N2``, ``N3_seed1..5``,
  ``N1_nogate``, ``N1_nomargin``, ``N1_nostable``, ``N1_F1``, ``N1_F5``,
  ``R_seed1..5``.
- Sequence-only comparator ``Q`` (volume-matched per genome to ``N1``'s gene calls)
  and ``Q_fixed`` (at a fixed score threshold).
- Masked-run recovery (gene-level including co-named proteins, reaction-level,
  pathway-level) in the off-growth stratum.
- Unmasked-run distinct labelled decisions ``(protein_id, ec)`` per genome:
  ``supported``, ``harmful``, ``related_substrate``, ``unconfirmed``.
- Abstention rates, call volume, LP/CPU cost, T1 fill-size distribution,
  fill stability rate, additional functions beyond ``Q``, and orphan steps.
- Genome-block bootstrap intervals (2,000 resamples, seed 20261102), paired
  exact McNemar test, and the fixed-precision rule (``P = 0.75``, ``n = 20``).
"""
from __future__ import annotations

import csv
from dataclasses import dataclass, field
import gzip
import json
import math
import os
import random
import sys
from pathlib import Path
from typing import (
    Any,
    Callable,
    Dict,
    Iterable,
    List,
    Mapping,
    Optional,
    Sequence,
    Set,
    Tuple,
    Union,
)

V3_DIR = Path(__file__).resolve().parents[1]
if str(V3_DIR) not in sys.path:
    sys.path.insert(0, str(V3_DIR))

from followups.common import GenomeLayout, base_seed_rxn_id, normalize_ec  # noqa: E402
from offgrowth.mask import (  # noqa: E402
    load_gapsmith_reactions_tbl,
    protein_aliases,
    strip_version,
)
from offgrowth.panel import (  # noqa: E402
    build_ec_to_reactions_map,
    is_complete_ec,
    load_seed_reactions,
)
from offgrowth.score import (  # noqa: E402
    _percentile_ci,
    exact_one_sided_mcnemar_p,
    holm_bonferroni_adjust,
)

from .candidates import best_by_protein, load_candidates  # noqa: E402
from .panels import load_org_id_map, measured_panel, random_panels, transporter_panel  # noqa: E402
from .rules import (  # noqa: E402
    FILL_MAX,
    LOCUS_WINDOW,
    MARGIN_MIN,
    MAX_CONAMED,
    SCORE_MIN,
    TIE_WINDOW,
    CandidateRankEntry,
    GeneCall,
    PathwayCall,
    ReactionCall,
    _build_rxn_to_ecs,
    _extract_excluded_reactions,
    _fast_base_rxn_id,
    _fast_normalize_ec,
    assign_gene,
    budget_cut_old_tests,
    candidate_index,
    clear_candidate_index_cache,
    old_rule_calls,
    pathway_calls,
    reaction_calls,
)

BOOTSTRAP_SEED: int = 20261102
BOOTSTRAP_RESAMPLES: int = 2000
FIXED_PRECISION_TARGET: float = 0.75
FIXED_PRECISION_MIN_DECISIONS: int = 20
CONTROL_SEEDS: Tuple[int, ...] = (1, 2, 3, 4, 5)

PRIMARY_ARMS: Tuple[str, ...] = (
    "S0",
    "S1",
    "S1b",
    "N1",
    "N2",
    "N3_seed1",
    "N3_seed2",
    "N3_seed3",
    "N3_seed4",
    "N3_seed5",
    "N1_nogate",
    "N1_nomargin",
    "N1_nostable",
    "N1_F1",
    "N1_F5",
    "N1_blocksub",
    "R_seed1",
    "R_seed2",
    "R_seed3",
    "R_seed4",
    "R_seed5",
)


def ec_first_three_levels(ec: str) -> str:
    """Return the first three levels of a normalized EC number (e.g. ``1.3.8`` for ``1.3.8.1``)."""
    ne = _fast_normalize_ec(ec)
    parts = ne.split(".")
    if len(parts) >= 3:
        return ".".join(parts[:3])
    return ne


_PERMUTED_CANDIDATES_CACHE: Dict[Tuple[str, int], Tuple[Any, Tuple[str, ...], Dict[Tuple[str, str], float]]] = {}


def permute_candidates_map_within_genome(
    candidates: Mapping[Tuple[str, str], float],
    seed: int,
    extra_proteins: Optional[Iterable[str]] = None,
    genome_id: Optional[str] = None,
) -> Dict[Tuple[str, str], float]:
    """Deterministically permute protein IDs of a candidate table within a genome.

    Matches the pilot's ``permute_candidates_within_genome`` semantics: the
    protein universe is the sorted union of candidate proteins and
    ``extra_proteins`` (e.g. ``model_genes``), shuffled with ``random.Random(seed)``.
    When ``genome_id`` is provided, the permuted map is cached once per ``(genome_id, seed)``.
    """
    extra_tuple = tuple(str(p) for p in extra_proteins if p) if extra_proteins is not None else ()
    if genome_id is not None:
        gkey = (str(genome_id), int(seed))
        ghit = _PERMUTED_CANDIDATES_CACHE.get(gkey)
        if ghit is not None and len(ghit[0]) == len(candidates) and ghit[1] == extra_tuple:
            return ghit[2]

    prot_set: Set[str] = {str(p) for (p, _) in candidates}
    if extra_tuple:
        prot_set.update(extra_tuple)
    proteins = sorted(prot_set)
    shuffled = list(proteins)
    rng = random.Random(int(seed))
    rng.shuffle(shuffled)
    remap = dict(zip(proteins, shuffled))
    permuted = {(remap[str(p)], str(e)): float(s) for (p, e), s in candidates.items()}

    if genome_id is not None:
        if len(_PERMUTED_CANDIDATES_CACHE) > 32:
            _PERMUTED_CANDIDATES_CACHE.clear()
        _PERMUTED_CANDIDATES_CACHE[(str(genome_id), int(seed))] = (candidates, extra_tuple, permuted)
    return permuted


def load_label_tsv(
    tsv_paths: Union[Path, str, Sequence[Union[Path, str]], Mapping[str, Sequence[str]], None]
) -> Dict[str, Dict[str, List[str]]]:
    """Load label TSV(s) into ``{genome_id: {protein_id: [ec, ...]}}``.

    Supports:
    - ``fitness_browser_ec_reannotations_*.tsv`` (``genome_id, protein_id, ..., function_id``)
    - ``swissprot_experimental_ecs_*.tsv`` (``genome_id, protein_id, ..., function_id``)
    - ``suite6_price149_manifest.tsv`` (``host_assembly, proteome_protein_id, ground_truth_ecs``)
    - Generic TSV with ``genome_id, protein_id, ec`` (or ``function_id``)
    - Pre-loaded dict ``{genome_id: {protein_id: [ec, ...]}}`` or ``{protein_id: [ec, ...]}``.
    """
    if tsv_paths is None:
        return {}
    if isinstance(tsv_paths, Mapping):
        out_map: Dict[str, Dict[str, List[str]]] = {}
        for k, v in tsv_paths.items():
            if isinstance(v, Mapping):
                out_map[str(k)] = {
                    str(p): sorted({normalize_ec(str(e)) for e in ecs if is_complete_ec(normalize_ec(str(e)))})
                    for p, ecs in v.items()
                }
            else:
                out_map.setdefault("*", {})[str(k)] = sorted(
                    {normalize_ec(str(e)) for e in (v or ()) if is_complete_ec(normalize_ec(str(e)))}
                )
        return out_map

    if isinstance(tsv_paths, (Path, str)):
        raw_list = [p.strip() for p in str(tsv_paths).split(",") if p.strip()]
        paths = [Path(p) for p in raw_list]
    else:
        paths = [Path(p) for p in tsv_paths]

    by_genome: Dict[str, Dict[str, Set[str]]] = {}
    for path in paths:
        if not path.exists():
            continue
        with path.open("r", encoding="utf-8", newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                gid = (
                    row.get("genome_id")
                    or row.get("host_assembly")
                    or row.get("org_id")
                    or "*"
                ).strip()
                pid = (
                    row.get("protein_id")
                    or row.get("proteome_protein_id")
                    or row.get("locus_id")
                    or ""
                ).strip()
                raw_ecs = (
                    row.get("function_id")
                    or row.get("ec")
                    or row.get("ground_truth_ecs")
                    or ""
                ).strip()
                if not pid or not raw_ecs:
                    continue
                for token in raw_ecs.replace(";", ",").split(","):
                    ne = normalize_ec(token.strip())
                    if ne and is_complete_ec(ne):
                        by_genome.setdefault(gid, {}).setdefault(pid, set()).add(ne)

    return {
        gid: {pid: sorted(ecs) for pid, ecs in pmap.items() if ecs}
        for gid, pmap in by_genome.items()
    }


def labels_for_genome(
    genome_id: str,
    label_store: Optional[Mapping[str, Any]],
    fallback_known_labels: Optional[Mapping[str, Sequence[str]]] = None,
) -> Dict[str, List[str]]:
    """Extract ``{protein_id: [ec, ...]}`` for ``genome_id`` from a label store or genome JSON."""
    merged: Dict[str, Set[str]] = {}
    if label_store:
        # Check if label_store is already {protein_id: [ec, ...]}
        sample_val = next(iter(label_store.values()), None)
        if sample_val is not None and not isinstance(sample_val, Mapping):
            for p, ecs in label_store.items():
                for e in ecs or ():
                    ne = normalize_ec(str(e))
                    if ne and is_complete_ec(ne):
                        merged.setdefault(str(p), set()).add(ne)
        else:
            for key in (genome_id, "*"):
                g_map = label_store.get(key)
                if isinstance(g_map, Mapping):
                    for p, ecs in g_map.items():
                        for e in ecs or ():
                            ne = normalize_ec(str(e))
                            if ne and is_complete_ec(ne):
                                merged.setdefault(str(p), set()).add(ne)
    elif fallback_known_labels:
        for p, ecs in fallback_known_labels.items():
            for e in ecs or ():
                ne = normalize_ec(str(e))
                if ne and is_complete_ec(ne):
                    merged.setdefault(str(p), set()).add(ne)

    return {p: sorted(ecs) for p, ecs in merged.items() if ecs}


@dataclass(frozen=True)
class DistinctDecision:
    genome_id: str
    arm: str
    protein_id: str
    ec: str
    score: float
    reactions: Tuple[str, ...]
    classification: str  # "supported" | "harmful" | "unconfirmed"
    related_substrate: bool
    truth_ecs: Tuple[str, ...] = ()


@dataclass(frozen=True)
class OrphanStep:
    genome_id: str
    rxn_id: str
    ecs: Tuple[str, ...]
    max_genome_score: float
    smallest_fill: int
    n_tests: int
    test_types: Tuple[str, ...]
    calling_tests: Tuple[str, ...]


@dataclass
class ArmRunEvaluation:
    genome_id: str
    arm: str
    run_kind: str
    case_id: Optional[str]
    stratum: Optional[str]
    reaction_calls: Dict[str, ReactionCall]
    pathway_calls: List[PathwayCall]
    n_admitted_tests: int
    n_candidate_reactions: int
    n_called_reactions: int
    n_gene_called_reactions: int
    n_reaction_only_calls: int
    n_pathway_calls: int
    gene_abstention_rate: float
    reaction_abstention_rate: float
    n_lps: int
    cpu_seconds: float
    t1_fill_sizes: Tuple[int, ...] = ()
    t1_stable_fill_sizes: Tuple[int, ...] = ()
    stability_rate: float = 0.0
    measured_panel_fallback: bool = False
    # Masked case recovery fields
    reaction_recovered: bool = False
    gene_recovered: bool = False
    pathway_recovered: bool = False
    recovering_fill_size: Optional[int] = None
    hidden_protein_rank: Optional[int] = None
    # Unmasked decision fields (populated against primary/test label set)
    decisions: Tuple[DistinctDecision, ...] = ()
    n_distinct_gene_calls: int = 0
    n_supported: int = 0
    n_harmful: int = 0
    n_related_substrate: int = 0
    n_unconfirmed: int = 0
    q_threshold: Optional[float] = None


@dataclass
class GenomeScoringContext:
    """Per-genome reference data needed to score all runs and arms."""
    genome_id: str
    candidates: Dict[Tuple[str, str], float]
    ec_to_reactions: Dict[str, List[str]]
    rxn_to_ecs: Dict[str, Set[str]]
    draft_reactions: Set[str]
    model_genes: List[str]
    draft_gpr_assignments: Dict[str, Set[str]]
    transporter_panel_ids: List[str]
    measured_panel_ids: Optional[List[str]]
    random_panel_by_seed: Dict[int, List[str]]
    universal_stoichiometry: Dict[str, Any]
    layout: Any
    known_labels: Dict[str, List[str]] = field(default_factory=dict)
    threshold_labels: Dict[str, List[str]] = field(default_factory=dict)
    test_labels: Dict[str, List[str]] = field(default_factory=dict)
    permuted_candidates_by_seed: Dict[int, Dict[Tuple[str, str], float]] = field(default_factory=dict)
    _assign_gene_cache: Dict[Tuple[Any, ...], Tuple[Tuple[CandidateRankEntry, ...], Optional[GeneCall]]] = field(
        default_factory=dict, repr=False
    )


def build_genome_context(
    genome_payload: Mapping[str, Any],
    candidates: Optional[Union[Path, str, Mapping[Tuple[str, str], float]]] = None,
    ec_to_reactions: Optional[Mapping[str, Iterable[str]]] = None,
    rxn_to_ecs: Optional[Mapping[str, Iterable[str]]] = None,
    transporters: Optional[Any] = None,
    matches_tsv: Optional[Any] = None,
    org_id: Optional[str] = None,
    universal_stoichiometry: Optional[Mapping[str, Any]] = None,
    layout: Optional[Any] = None,
    draft_gpr_assignments: Optional[Mapping[str, Iterable[str]]] = None,
    threshold_labels: Optional[Mapping[str, Sequence[str]]] = None,
    test_labels: Optional[Mapping[str, Sequence[str]]] = None,
    seeds: Sequence[int] = CONTROL_SEEDS,
    org_map_tsv: Optional[Any] = None,
) -> GenomeScoringContext:
    """Build a ``GenomeScoringContext`` from a genome JSON payload and optional external tables."""
    gid = str(genome_payload.get("genome_id") or "")
    inputs_map = genome_payload.get("inputs") or {}
    input_paths = [Path(p) for p in inputs_map.keys() if Path(p).is_file()] if isinstance(inputs_map, Mapping) else []

    # Discover reference paths from payload["inputs"] when not explicitly provided
    seed_rxn_path = next((p for p in input_paths if "seed_reactions" in p.name), None)
    ec_alias_path = next((p for p in input_paths if "EC_alias" in p.name or "ec_alias" in p.name.lower()), None)
    trans_tbl_path = next((p for p in input_paths if p.name.endswith("-Transporter.tbl")), None)
    rxns_tbl_path = next((p for p in input_paths if p.name.endswith("-Reactions.tbl")), None)
    layout_tsv_path = next((p for p in input_paths if p.name.endswith(".layout.tsv") or "layout" in p.name), None)
    cand_tsv_path = next(
        (
            p for p in input_paths
            if "candidate" in p.name.lower() or p.name.endswith(".clean.tsv")
        ),
        None,
    )

    # Candidates
    raw_cands = (
        candidates
        if candidates is not None
        else (genome_payload.get("candidates") if "candidates" in genome_payload else cand_tsv_path)
    )
    cand_map = load_candidates(raw_cands) if raw_cands else {}

    # Reaction <-> EC mappings and universal stoichiometry
    raw_e2r = ec_to_reactions if ec_to_reactions is not None else genome_payload.get("ec_to_reactions", {})
    raw_r2e = rxn_to_ecs if rxn_to_ecs is not None else genome_payload.get("rxn_to_ecs", {})
    auto_stoich: Dict[str, List[str]] = {}
    if not raw_e2r and not raw_r2e and seed_rxn_path is not None:
        seed_rxns = load_seed_reactions(seed_rxn_path, ec_alias_path)  # ECs come from the alias table (pilot panel.py)
        raw_e2r = build_ec_to_reactions_map(seed_rxns)
        for r_id, sr in seed_rxns.items():
            if sr.is_usable:
                auto_stoich[base_seed_rxn_id(r_id)] = [t.cpd_id for t in sr.stoichiometry]

    r2e = _build_rxn_to_ecs(rxn_to_ecs=raw_r2e, ec_to_reactions=raw_e2r)
    e2r: Dict[str, List[str]] = {}
    if raw_e2r:
        for e, rxns in raw_e2r.items():
            ne = normalize_ec(str(e))
            if ne:
                e2r[ne] = sorted({base_seed_rxn_id(str(r)) for r in (rxns or ()) if base_seed_rxn_id(str(r))})
    for br, ecs in r2e.items():
        for e in ecs:
            if br not in e2r.setdefault(e, []):
                e2r[e].append(br)
    for e in e2r:
        e2r[e].sort()

    draft_rxns = {base_seed_rxn_id(str(r)) for r in (genome_payload.get("draft_reactions") or ()) if str(r)}
    model_genes = [str(g) for g in (genome_payload.get("model_genes") or ())]

    # Draft GPR assignments: {protein_id: {ec, ...}}
    raw_gpr = (
        draft_gpr_assignments
        if draft_gpr_assignments is not None
        else genome_payload.get("draft_gpr_assignments")
    )
    gpr_map: Dict[str, Set[str]] = {}
    if isinstance(raw_gpr, Mapping):
        for p, ecs_or_rxns in raw_gpr.items():
            for item in ecs_or_rxns or ():
                s_item = str(item).strip()
                ne = normalize_ec(s_item)
                if ne and is_complete_ec(ne):
                    gpr_map.setdefault(str(p), set()).add(ne)
                else:
                    br = base_seed_rxn_id(s_item)
                    for mapped_e in r2e.get(br, set()):
                        gpr_map.setdefault(str(p), set()).add(mapped_e)
    elif isinstance(genome_payload.get("draft_reaction_genes"), Mapping):
        for r, prots in genome_payload["draft_reaction_genes"].items():
            br = base_seed_rxn_id(str(r))
            mapped_ecs = r2e.get(br, set())
            for p in prots or ():
                if mapped_ecs:
                    gpr_map.setdefault(str(p), set()).update(mapped_ecs)
    elif rxns_tbl_path is not None:
        for hit in load_gapsmith_reactions_tbl(rxns_tbl_path):
            # ReactionEvidenceHit carries `ec_numbers` and `reaction_ids` (pilot mask.py); draft assignments are
            # the complete ECs it lists plus the ECs mapped from its reactions.
            if not hit.protein_id:
                continue
            ecs = {e for e in (getattr(hit, "ec_numbers", ()) or ()) if is_complete_ec(e)}
            for r in getattr(hit, "reaction_ids", ()) or ():
                ecs.update(r2e.get(base_seed_rxn_id(str(r)), set()))
            if ecs:
                gpr_map.setdefault(hit.protein_id, set()).update(ecs)

    # Panels (P_T, P_M, N3 random panels)
    panel_spec = genome_payload.get("panel") or []
    eff_transporters = transporters if transporters is not None else (genome_payload.get("transporters") or trans_tbl_path)
    if "transporter_panel" in genome_payload and transporters is None:
        p_t = [str(x) for x in genome_payload["transporter_panel"]]
    elif eff_transporters is not None:
        p_t = transporter_panel(panel_spec, eff_transporters)
    else:
        flagged = [
            str(entry.get("panel_id"))
            for entry in panel_spec
            if isinstance(entry, Mapping) and entry.get("has_transporter")
        ]
        p_t = flagged if flagged else [
            str(entry["panel_id"]) if isinstance(entry, Mapping) else str(entry)
            for entry in panel_spec
        ]

    if "measured_panel" in genome_payload and matches_tsv is None:
        raw_pm = genome_payload["measured_panel"]
        p_m = [str(x) for x in raw_pm] if raw_pm else None
    elif matches_tsv is not None:
        resolved_org = org_id
        if resolved_org is None and org_map_tsv is not None:
            if isinstance(org_map_tsv, Mapping):
                resolved_org = org_map_tsv.get(gid)
            else:
                resolved_org = load_org_id_map(org_map_tsv).get(gid)
        target_org = resolved_org if resolved_org is not None else str(genome_payload.get("org_id") or gid)
        p_m = measured_panel(panel_spec, matches_tsv, target_org)
    else:
        p_m = None

    rand_panels = random_panels(panel_spec, size=len(p_t), seeds=seeds)

    stoich = (
        dict(universal_stoichiometry)
        if universal_stoichiometry is not None
        else dict(genome_payload.get("universal_stoichiometry") or auto_stoich)
    )
    resolved_layout = layout if layout is not None else genome_payload.get("layout")
    if resolved_layout is None and layout_tsv_path is not None:
        resolved_layout = GenomeLayout.from_tsv(layout_tsv_path)

    known_lbls = labels_for_genome(gid, None, genome_payload.get("known_labels") or {})
    thresh_lbls = labels_for_genome(gid, threshold_labels, known_lbls) if threshold_labels is not None else dict(known_lbls)
    tst_lbls = labels_for_genome(gid, test_labels, known_lbls) if test_labels is not None else dict(known_lbls)

    return GenomeScoringContext(
        genome_id=gid,
        candidates=cand_map,
        ec_to_reactions=e2r,
        rxn_to_ecs=r2e,
        draft_reactions=draft_rxns,
        model_genes=model_genes,
        draft_gpr_assignments=gpr_map,
        transporter_panel_ids=p_t,
        measured_panel_ids=p_m,
        random_panel_by_seed=rand_panels,
        universal_stoichiometry=stoich,
        layout=resolved_layout,
        known_labels=known_lbls,
        threshold_labels=thresh_lbls,
        test_labels=tst_lbls,
    )


def classify_decision(
    protein_id: str,
    ec: str,
    labels_by_protein: Mapping[str, Sequence[str]],
    ec_to_reactions: Mapping[str, Sequence[str]],
) -> Tuple[str, bool, Tuple[str, ...]]:
    """Classify a single ``(protein_id, ec)`` decision against ``labels_by_protein``.

    Returns ``(classification, related_substrate, truth_ecs)`` where
    ``classification`` is ``'supported'``, ``'harmful'``, or ``'unconfirmed'``.
    Uses reaction equivalence through ``ec_to_reactions`` (matching the pilot's
    ``classify_call``) as well as direct EC equality.
    """
    truth_ecs: Optional[Sequence[str]] = None
    for alias in protein_aliases(protein_id):
        if alias in labels_by_protein:
            truth_ecs = labels_by_protein[alias]
            break
    if not truth_ecs:
        return "unconfirmed", False, ()

    norm_truth = tuple(sorted({normalize_ec(str(e)) for e in truth_ecs if normalize_ec(str(e))}))
    if not norm_truth:
        return "unconfirmed", False, ()

    norm_call_ec = normalize_ec(ec)
    if norm_call_ec in norm_truth:
        return "supported", False, norm_truth

    call_rxns = {base_seed_rxn_id(str(r)) for r in ec_to_reactions.get(norm_call_ec, ()) if str(r)}
    truth_rxns: Set[str] = set()
    for te in norm_truth:
        for r in ec_to_reactions.get(te, ()):
            if str(r):
                truth_rxns.add(base_seed_rxn_id(str(r)))

    if call_rxns and (call_rxns & truth_rxns):
        return "supported", False, norm_truth

    call_prefix3 = ec_first_three_levels(norm_call_ec)
    related = any(
        call_prefix3 and len(call_prefix3.split(".")) == 3 and ec_first_three_levels(te) == call_prefix3
        for te in norm_truth
    )
    return "harmful", related, norm_truth


def extract_distinct_decisions(
    genome_id: str,
    arm: str,
    rxn_calls: Mapping[str, ReactionCall],
    labels_by_protein: Mapping[str, Sequence[str]],
    ec_to_reactions: Mapping[str, Sequence[str]],
) -> List[DistinctDecision]:
    """Extract distinct ``(protein_id, ec)`` decisions from a run's reaction calls (once per genome)."""
    by_pair: Dict[Tuple[str, str], Tuple[float, Set[str]]] = {}
    for rxn_id, rc in rxn_calls.items():
        if rc.gene_call is None:
            continue
        pairs = (
            rc.gene_call.protein_ec_pairs
            if rc.gene_call.protein_ec_pairs
            else tuple(zip(rc.gene_call.proteins, rc.gene_call.ecs, rc.gene_call.scores))
        )
        for p, e, s in pairs:
            ne = normalize_ec(str(e))
            key = (str(p), ne)
            cur = by_pair.get(key)
            if cur is None:
                by_pair[key] = (float(s), {rxn_id})
            else:
                by_pair[key] = (max(cur[0], float(s)), cur[1] | {rxn_id})

    decisions: List[DistinctDecision] = []
    for (p, e) in sorted(by_pair):
        score, rxns = by_pair[(p, e)]
        cls, related, truth_ecs = classify_decision(p, e, labels_by_protein, ec_to_reactions)
        decisions.append(
            DistinctDecision(
                genome_id=genome_id,
                arm=arm,
                protein_id=p,
                ec=e,
                score=score,
                reactions=tuple(sorted(rxns)),
                classification=cls,
                related_substrate=related,
                truth_ecs=truth_ecs,
            )
        )
    return decisions


def _panel_test_filter(panel_ids: Iterable[str]) -> Callable[[Mapping[str, Any]], bool]:
    panel_set = {str(x) for x in panel_ids}

    def _filter(test: Mapping[str, Any]) -> bool:
        ttype = str(test.get("test_type") or "")
        if ttype == "T1":
            pid = test.get("panel_id") or test.get("target_id") or ""
            return str(pid) in panel_set
        return ttype in ("T2", "T3", "A2", "A3")

    return _filter


_RUN_COST_CACHE: Dict[Tuple[int, Any, bool, int], Tuple[Any, Tuple[int, int, float, Set[str], Tuple[int, ...], Tuple[int, ...], float]]] = {}
_RUN_RXN_CALLS_CACHE: Dict[Tuple[int, Optional[float], bool, Any, int], Tuple[Any, Dict[str, ReactionCall]]] = {}


def _compute_run_cost_and_candidates(
    run: Mapping[str, Any],
    test_filter: Callable[[Mapping[str, Any]], bool],
    include_baseline: bool,
    draft_reactions: Set[str],
    require_stable: bool,
    cache_key: Optional[Any] = None,
) -> Tuple[int, int, float, Set[str], Tuple[int, ...], Tuple[int, ...], float]:
    """Return ``(n_admitted_tests, n_lps, cpu_seconds, candidate_rxns, t1_fills, t1_stable_fills, stability_rate)``."""
    if cache_key is not None:
        ckey = (id(run), cache_key, bool(include_baseline), id(draft_reactions))
        chit = _RUN_COST_CACHE.get(ckey)
        if chit is not None and chit[0] is run:
            return chit[1]

    excluded = _extract_excluded_reactions(run, draft_reactions=draft_reactions, exclude_baseline=include_baseline)
    n_admitted = 0
    n_lps = 0
    cpu_sec = 0.0
    if include_baseline:
        bf = run.get("baseline_fill")
        if isinstance(bf, Mapping):
            n_lps += int(bf.get("n_lps") or 0)
            cpu_sec += float(bf.get("cpu_seconds") or 0.0)

    cand_rxns: Set[str] = set()
    t1_fills: List[int] = []
    t1_stable_fills: List[int] = []
    total_fill_elems = 0
    total_stable_elems = 0

    for test in run.get("tests", ()) or ():
        if not test_filter(test):
            continue
        n_admitted += 1
        n_lps += int(test.get("n_lps") or 0)
        cpu_sec += float(test.get("cpu_seconds") or 0.0)

        raw_fill = [_fast_base_rxn_id(r) for r in (test.get("fill") or ()) if str(r)]
        raw_stable = [_fast_base_rxn_id(r) for r in (test.get("stable_fill") or ()) if str(r)]
        fill_set = {r for r in raw_fill if r and r not in excluded}

        # Candidate reactions added by any admitted test (before the fill-size/stability gate)
        cand_rxns.update(fill_set)

        if raw_fill:
            u_fill = set(raw_fill)
            u_stable = set(raw_stable)
            total_fill_elems += len(u_fill)
            total_stable_elems += len(u_stable & u_fill)
            if str(test.get("test_type") or "") == "T1":
                t1_fills.append(len(u_fill))
                t1_stable_fills.append(len(u_stable))

    stab_rate = (
        float(total_stable_elems) / float(total_fill_elems)
        if total_fill_elems > 0
        else 0.0
    )
    res = (
        n_admitted,
        n_lps,
        cpu_sec,
        cand_rxns,
        tuple(t1_fills),
        tuple(t1_stable_fills),
        stab_rate,
    )
    if cache_key is not None:
        if len(_RUN_COST_CACHE) > 64:
            _RUN_COST_CACHE.clear()
        _RUN_COST_CACHE[(id(run), cache_key, bool(include_baseline), id(draft_reactions))] = (run, res)
    return res


def _evaluate_masked_recovery(
    run: Mapping[str, Any],
    rxn_calls: Mapping[str, ReactionCall],
    pw_calls: Sequence[PathwayCall],
    ctx: GenomeScoringContext,
    candidates_used: Mapping[Tuple[str, str], float],
    seed_num: Optional[int] = None,
) -> Tuple[bool, bool, bool, Optional[int], Optional[int]]:
    """Compute ``(reaction_recovered, gene_recovered, pathway_recovered, recovering_fill_size, hidden_protein_rank)``."""
    case = run.get("case")
    if not isinstance(case, Mapping):
        return False, False, False, None, None

    hidden_prot = str(case.get("protein_id") or "")
    hidden_aliases = protein_aliases(hidden_prot) if hidden_prot else set()
    truth_ecs = [_fast_normalize_ec(e) for e in (case.get("truth_ecs") or ()) if _fast_normalize_ec(e)]

    target_rxns: Set[str] = {
        _fast_base_rxn_id(r) for r in (run.get("target_reactions") or ()) if str(r)
    }
    for te in truth_ecs:
        for r in ctx.ec_to_reactions.get(te, ()):
            target_rxns.add(_fast_base_rxn_id(r))
    # Also include reactions whose rxn_to_ecs intersects truth_ecs
    if truth_ecs:
        truth_ec_set = set(truth_ecs)
        for r, ecs in ctx.rxn_to_ecs.items():
            if ecs & truth_ec_set:
                target_rxns.add(r)

    called_targets = [rc for r_id, rc in rxn_calls.items() if r_id in target_rxns]
    rxn_rec = len(called_targets) > 0
    rec_fill_size = min((rc.smallest_fill for rc in called_targets), default=None)

    gene_rec = False
    for rc in called_targets:
        if rc.gene_call is not None:
            for named_p in rc.gene_call.proteins:
                if named_p in hidden_aliases or strip_version(named_p) in hidden_aliases:
                    gene_rec = True
                    break
        if gene_rec:
            break

    pw_rec = False
    for pc in pw_calls:
        has_target_rxn = any(r in target_rxns for r in pc.reactions)
        has_hidden_gene = any(
            g in hidden_aliases or strip_version(g) in hidden_aliases
            for g in pc.genes
        )
        if has_target_rxn and has_hidden_gene:
            pw_rec = True
            break

    # Compute 1-based rank of hidden_prot among candidates for truth_ecs (or target_rxns ECs)
    eval_ecs: Set[str] = set(truth_ecs)
    for r in target_rxns:
        eval_ecs.update(ctx.rxn_to_ecs.get(r, set()))

    idx = candidate_index(candidates_used, genome_id=ctx.genome_id, seed=seed_num)
    prot_scores: Dict[str, float] = {}
    for ne in eval_ecs:
        for p, s in idx.ec_to_prot_scores.get(ne, ()):
            if s > prot_scores.get(p, 0.0):
                prot_scores[p] = s

    hidden_rank: Optional[int] = None
    if prot_scores and hidden_aliases:
        ranked_prots = sorted(prot_scores.keys(), key=lambda p: (-prot_scores[p], p))
        for idx_rank, p in enumerate(ranked_prots, start=1):
            if p in hidden_aliases or strip_version(p) in hidden_aliases:
                hidden_rank = idx_rank
                break

    return rxn_rec, gene_rec, pw_rec, rec_fill_size, hidden_rank


_NORM_ROWS_CACHE: Dict[int, Tuple[Any, List[Tuple[str, str, float]]]] = {}
_GENOME_NORM_ROWS_CACHE: Dict[Tuple[str, Optional[int]], Tuple[Any, List[Tuple[str, str, float]]]] = {}
_BEST_BY_PROTEIN_CACHE: Dict[int, Tuple[Any, Any]] = {}
_GENOME_BEST_BY_PROTEIN_CACHE: Dict[Tuple[str, Optional[int]], Tuple[Any, Any]] = {}


def _normalized_candidate_rows(
    candidates: Mapping[Tuple[str, str], float],
    genome_id: Optional[str] = None,
    seed: Optional[int] = None,
) -> List[Tuple[str, str, float]]:
    """`(normalized_ec, protein, score)` rows of a candidate map, cached by `(genome_id, seed)` and map identity."""
    gkey: Optional[Tuple[str, Optional[int]]] = None
    if genome_id is not None:
        gkey = (str(genome_id), int(seed) if seed is not None else None)
        ghit = _GENOME_NORM_ROWS_CACHE.get(gkey)
        if ghit is not None and len(ghit[0]) == len(candidates):
            return ghit[1]
    key = id(candidates)
    hit = _NORM_ROWS_CACHE.get(key)
    if hit is not None and hit[0] is candidates:
        rows = hit[1]
        if gkey is not None:
            _GENOME_NORM_ROWS_CACHE[gkey] = (candidates, rows)
        return rows
    rows = [(_fast_normalize_ec(e), p, float(s)) for (p, e), s in candidates.items()]
    if len(_NORM_ROWS_CACHE) > 64:
        _NORM_ROWS_CACHE.clear()
    if len(_GENOME_NORM_ROWS_CACHE) > 64:
        _GENOME_NORM_ROWS_CACHE.clear()
    _NORM_ROWS_CACHE[key] = (candidates, rows)
    if gkey is not None:
        _GENOME_NORM_ROWS_CACHE[gkey] = (candidates, rows)
    return rows


def cached_best_by_protein(
    candidates: Mapping[Tuple[str, str], float],
    genome_id: Optional[str] = None,
    seed: Optional[int] = None,
):
    """`best_by_protein(candidates)` cached by `(genome_id, seed)` and map identity (speed only)."""
    gkey: Optional[Tuple[str, Optional[int]]] = None
    if genome_id is not None:
        gkey = (str(genome_id), int(seed) if seed is not None else None)
        ghit = _GENOME_BEST_BY_PROTEIN_CACHE.get(gkey)
        if ghit is not None and len(ghit[0]) == len(candidates):
            return ghit[1]
    key = id(candidates)
    hit = _BEST_BY_PROTEIN_CACHE.get(key)
    if hit is not None and hit[0] is candidates:
        res = hit[1]
        if gkey is not None:
            _GENOME_BEST_BY_PROTEIN_CACHE[gkey] = (candidates, res)
        return res
    res = best_by_protein(candidates)
    if len(_BEST_BY_PROTEIN_CACHE) > 64:
        _BEST_BY_PROTEIN_CACHE.clear()
    if len(_GENOME_BEST_BY_PROTEIN_CACHE) > 64:
        _GENOME_BEST_BY_PROTEIN_CACHE.clear()
    _BEST_BY_PROTEIN_CACHE[key] = (candidates, res)
    if gkey is not None:
        _GENOME_BEST_BY_PROTEIN_CACHE[gkey] = (candidates, res)
    return res


def clear_genome_caches(genome_id: Optional[str] = None) -> None:
    """Release per-genome candidate and index caches after a genome is scored."""
    clear_candidate_index_cache(genome_id)
    _RUN_COST_CACHE.clear()
    _RUN_RXN_CALLS_CACHE.clear()
    _RUN_MASK_SETUP_CACHE.clear()
    if genome_id is None:
        _PERMUTED_CANDIDATES_CACHE.clear()
        _NORM_ROWS_CACHE.clear()
        _GENOME_NORM_ROWS_CACHE.clear()
        _BEST_BY_PROTEIN_CACHE.clear()
        _GENOME_BEST_BY_PROTEIN_CACHE.clear()
        return
    gid = str(genome_id)
    for k in [k for k in _PERMUTED_CANDIDATES_CACHE if k[0] == gid]:
        _, _, pmap = _PERMUTED_CANDIDATES_CACHE.pop(k)
        _NORM_ROWS_CACHE.pop(id(pmap), None)
        _BEST_BY_PROTEIN_CACHE.pop(id(pmap), None)
    for k in [k for k in _GENOME_NORM_ROWS_CACHE if k[0] == gid]:
        c_ref, _ = _GENOME_NORM_ROWS_CACHE.pop(k)
        _NORM_ROWS_CACHE.pop(id(c_ref), None)
    for k in [k for k in _GENOME_BEST_BY_PROTEIN_CACHE if k[0] == gid]:
        c_ref, _ = _GENOME_BEST_BY_PROTEIN_CACHE.pop(k)
        _BEST_BY_PROTEIN_CACHE.pop(id(c_ref), None)


BLOCK_MIN_FILL: int = 20      # plan §13 deviation 3: T1 fills at least this large define the common block
BLOCK_MIN_SHARE: float = 0.5  # a reaction is in the block when it appears in at least this share of those fills


def common_t1_block(run: Mapping[str, Any]) -> Set[str]:
    """Label-free common block of a run: reactions in >= BLOCK_MIN_SHARE of its T1 fills with >= BLOCK_MIN_FILL reactions."""
    big = [set(t.get("fill") or ()) for t in (run.get("tests") or ())
           if t.get("test_type") == "T1" and t.get("gapfill_ran") and len(t.get("fill") or ()) >= BLOCK_MIN_FILL]
    if not big:
        return set()
    counts: Dict[str, int] = {}
    for fill in big:
        for r in fill:
            counts[r] = counts.get(r, 0) + 1
    return {r for r, c in counts.items() if c >= BLOCK_MIN_SHARE * len(big)}


def _subtract_common_block(run: Mapping[str, Any]) -> Dict[str, Any]:
    """Copy of `run` whose T1 fills and stable fills exclude the run's common block (other tests unchanged)."""
    block = common_t1_block(run)
    out = dict(run)
    tests = []
    for t in run.get("tests") or ():
        if t.get("test_type") == "T1" and block:
            t2 = dict(t)
            t2["fill"] = [r for r in (t.get("fill") or ()) if r not in block]
            t2["stable_fill"] = [r for r in (t.get("stable_fill") or ()) if r not in block]
            tests.append(t2)
        else:
            tests.append(t)
    out["tests"] = tests
    out["common_block"] = sorted(block)
    return out


_RUN_MASK_SETUP_CACHE: Dict[Tuple[int, int], Tuple[Any, Any, Set[str], Dict[str, Set[str]], Set[str]]] = {}


def _resolve_run_draft_and_gpr(
    run: Mapping[str, Any],
    ctx: GenomeScoringContext,
    run_kind: str,
    case: Any,
) -> Tuple[Set[str], Dict[str, Set[str]], Set[str]]:
    ckey = (id(run), id(ctx))
    chit = _RUN_MASK_SETUP_CACHE.get(ckey)
    if chit is not None and chit[0] is run and chit[1] is ctx:
        return chit[2], chit[3], chit[4]

    target_rxn_set = {_fast_base_rxn_id(r) for r in (run.get("target_reactions") or ()) if str(r)}
    if "draft_reactions" in run and run.get("draft_reactions") is not None:
        run_draft_rxns = {_fast_base_rxn_id(r) for r in (run.get("draft_reactions") or ()) if str(r)}
    elif run_kind == "masked":
        run_draft_rxns = ctx.draft_reactions - target_rxn_set
    else:
        run_draft_rxns = ctx.draft_reactions

    run_gpr_assignments: Dict[str, Set[str]] = ctx.draft_gpr_assignments
    hidden_aliases: Set[str] = set()
    if run_kind == "masked" and isinstance(case, Mapping):
        hidden_p = str(case.get("protein_id") or "")
        hidden_aliases = protein_aliases(hidden_p) if hidden_p else set()
        if ctx.draft_gpr_assignments:
            masked_ecs: Set[str] = {_fast_normalize_ec(e) for e in (case.get("truth_ecs") or ()) if _fast_normalize_ec(e)}
            audit_obj = run.get("audit")
            if isinstance(audit_obj, Mapping):
                masked_ecs.update(
                    _fast_normalize_ec(e) for e in (audit_obj.get("equivalent_ecs") or ()) if _fast_normalize_ec(e)
                )
            for tr in target_rxn_set:
                masked_ecs.update(ctx.rxn_to_ecs.get(tr, set()))
            run_gpr_assignments = {}
            for p, ecs in ctx.draft_gpr_assignments.items():
                if p in hidden_aliases or strip_version(p) in hidden_aliases:
                    rem = set(ecs) - masked_ecs
                    if rem:
                        run_gpr_assignments[p] = rem
                else:
                    run_gpr_assignments[p] = set(ecs)

    if len(_RUN_MASK_SETUP_CACHE) > 32:
        _RUN_MASK_SETUP_CACHE.clear()
    _RUN_MASK_SETUP_CACHE[ckey] = (run, ctx, run_draft_rxns, run_gpr_assignments, hidden_aliases)
    return run_draft_rxns, run_gpr_assignments, hidden_aliases


def evaluate_arm_on_run(
    run: Mapping[str, Any],
    arm: str,
    ctx: GenomeScoringContext,
    labels_by_protein: Optional[Mapping[str, Sequence[str]]] = None,
    n1_lp_budget: Optional[int] = None,
    _skip_decisions: bool = False,
) -> ArmRunEvaluation:
    """Evaluate one arm (plan §6) on one run (masked or unmasked)."""
    run_kind = str(run.get("run_kind") or "unmasked")
    case = run.get("case")
    case_id = str(case.get("case_id")) if isinstance(case, Mapping) and case.get("case_id") else None
    stratum_obj = run.get("stratum")
    if isinstance(stratum_obj, Mapping):
        stratum = str(stratum_obj.get("stratum") or "") or None
    elif isinstance(stratum_obj, str):
        stratum = stratum_obj
    else:
        stratum = None

    active_labels = (
        labels_by_protein
        if labels_by_protein is not None
        else (ctx.test_labels if ctx.test_labels else ctx.known_labels)
    )

    measured_fallback = False
    candidates_used: Mapping[Tuple[str, str], float] = ctx.candidates
    seed_for_cands: Optional[int] = None

    # Determine run-specific draft_reactions and draft_gpr_assignments (respecting case masking)
    run_draft_rxns, run_gpr_assignments, hidden_aliases = _resolve_run_draft_and_gpr(run, ctx, run_kind, case)

    cost_cache_key: Optional[Any] = None

    # Determine test filter and rule parameters for the requested arm
    if arm == "S0":
        t_filter: Callable[[Mapping[str, Any]], bool] = lambda t: str(t.get("test_type") or "") == "A0"
        cost_cache_key = ("S0",)
        rxn_calls = old_rule_calls(
            run,
            test_filter=t_filter,
            candidates=candidates_used,
            rxn_to_ecs=ctx.rxn_to_ecs,
            draft_reactions=run_draft_rxns,
            genome_id=ctx.genome_id,
            seed=None,
        )
        pw_calls: List[PathwayCall] = []
        include_baseline = False
        require_stable = False

    elif arm == "S1":
        t_filter = lambda t: str(t.get("test_type") or "") in ("A1old", "A1", "T2", "A2", "T3", "A3")
        cost_cache_key = ("S1",)
        rxn_calls = old_rule_calls(
            run,
            test_filter=t_filter,
            candidates=candidates_used,
            rxn_to_ecs=ctx.rxn_to_ecs,
            draft_reactions=run_draft_rxns,
            genome_id=ctx.genome_id,
            seed=None,
        )
        pw_calls = []
        include_baseline = False
        require_stable = False

    elif arm == "S1b":
        if n1_lp_budget is None:
            n1_eval = evaluate_arm_on_run(run, "N1", ctx, labels_by_protein=active_labels, _skip_decisions=_skip_decisions)
            n1_lp_budget = n1_eval.n_lps
        cost_cache_key = ("S1b", int(n1_lp_budget))
        kept_tests = budget_cut_old_tests(
            run,
            lp_budget=n1_lp_budget,
            transporter_panel_ids=ctx.transporter_panel_ids,
        )
        kept_ids = {str(t.get("test_id") or "") for t in kept_tests}
        kept_obj_ids = {id(t) for t in kept_tests}
        t_filter = lambda t: (str(t.get("test_id") or "") in kept_ids) if t.get("test_id") else (id(t) in kept_obj_ids)
        rxn_calls = old_rule_calls(
            run,
            test_filter=t_filter,
            candidates=candidates_used,
            rxn_to_ecs=ctx.rxn_to_ecs,
            draft_reactions=run_draft_rxns,
            genome_id=ctx.genome_id,
            seed=None,
        )
        pw_calls = []
        include_baseline = False
        require_stable = False

    else:
        # New-rule arms: N1, N2, N3_seed*, N1_nogate, N1_nomargin, N1_nostable, N1_F1, N1_F5, R_seed*
        fill_max: Optional[float] = float(FILL_MAX)
        require_stable = True
        score_min = SCORE_MIN
        margin_min = MARGIN_MIN
        tie_window = TIE_WINDOW
        require_free_or_mutual = True
        panel_ids = ctx.transporter_panel_ids

        if arm == "N1":
            pass
        elif arm == "N2":
            if ctx.measured_panel_ids:
                panel_ids = ctx.measured_panel_ids
            else:
                panel_ids = ctx.transporter_panel_ids
                measured_fallback = True
        elif arm.startswith("N3_seed"):
            seed_num = int(arm.split("N3_seed", 1)[1])
            if seed_num in ctx.random_panel_by_seed:
                panel_ids = ctx.random_panel_by_seed[seed_num]
            else:
                panel_ids = random_panels(ctx.transporter_panel_ids, len(ctx.transporter_panel_ids), seeds=(seed_num,))[seed_num]
        elif arm == "N1_nogate":
            fill_max = float("inf")
            require_stable = False
        elif arm == "N1_nomargin":
            margin_min = 0.0
            require_free_or_mutual = False
        elif arm == "N1_nostable":
            require_stable = False
        elif arm == "N1_blocksub":
            # Plan §13 deviation 3 (reported variant, not primary): fill sizes after removing the run's common block.
            run = _subtract_common_block(run)
        elif arm == "N1_F1":
            fill_max = 1.0
        elif arm == "N1_F5":
            fill_max = 5.0
        elif arm.startswith("R_seed"):
            seed_num = int(arm.split("R_seed", 1)[1])
            seed_for_cands = seed_num
            cached_perm = ctx.permuted_candidates_by_seed.get(seed_num)
            if cached_perm is None:
                cached_perm = permute_candidates_map_within_genome(
                    ctx.candidates,
                    seed=seed_num,
                    extra_proteins=ctx.model_genes,
                    genome_id=ctx.genome_id,
                )
                ctx.permuted_candidates_by_seed[seed_num] = cached_perm
            candidates_used = cached_perm
        else:
            raise ValueError(f"Unknown evaluation arm: {arm!r}")

        panel_tuple = tuple(panel_ids)
        cost_cache_key = ("panel", panel_tuple)
        t_filter = _panel_test_filter(panel_ids)
        include_baseline = True

        rxn_cache_key = (id(run), fill_max, bool(require_stable), cost_cache_key, id(run_draft_rxns))
        rxn_hit = _RUN_RXN_CALLS_CACHE.get(rxn_cache_key)
        if rxn_hit is not None and rxn_hit[0] is run:
            rxn_calls = {
                k: ReactionCall(
                    rxn_id=v.rxn_id,
                    calling_tests=v.calling_tests,
                    smallest_fill=v.smallest_fill,
                    n_tests=v.n_tests,
                    test_types=v.test_types,
                    all_tests=v.all_tests,
                )
                for k, v in rxn_hit[1].items()
            }
        else:
            base_calls = reaction_calls(
                run,
                fill_max=fill_max,
                require_stable=require_stable,
                test_filter=t_filter,
                draft_reactions=run_draft_rxns,
            )
            if len(_RUN_RXN_CALLS_CACHE) > 64:
                _RUN_RXN_CALLS_CACHE.clear()
            _RUN_RXN_CALLS_CACHE[rxn_cache_key] = (run, base_calls)
            rxn_calls = {
                k: ReactionCall(
                    rxn_id=v.rxn_id,
                    calling_tests=v.calling_tests,
                    smallest_fill=v.smallest_fill,
                    n_tests=v.n_tests,
                    test_types=v.test_types,
                    all_tests=v.all_tests,
                )
                for k, v in base_calls.items()
            }

        idx_cands = candidate_index(candidates_used, genome_id=ctx.genome_id, seed=seed_for_cands)
        prot_ec_score = idx_cands.prot_ec_score
        assign_cache = ctx._assign_gene_cache

        for r_id, rc in rxn_calls.items():
            ecs = ctx.rxn_to_ecs.get(r_id, set())
            touches_hidden = bool(hidden_aliases) and any(
                (alias, e) in prot_ec_score for alias in hidden_aliases for e in ecs
            )
            gpr_key = id(run_gpr_assignments) if touches_hidden else id(ctx.draft_gpr_assignments)
            akey = (
                r_id,
                id(candidates_used),
                gpr_key,
                score_min,
                margin_min,
                tie_window,
                MAX_CONAMED,
                require_free_or_mutual,
            )
            ahit = assign_cache.get(akey)
            if ahit is not None:
                rc.top_candidates = ahit[0]
                rc.gene_call = ahit[1]
            else:
                assign_gene(
                    rc,
                    ec_set=ecs,
                    candidates=candidates_used,
                    draft_gpr_assignments=run_gpr_assignments,
                    score_min=score_min,
                    margin_min=margin_min,
                    tie_window=tie_window,
                    max_conamed=MAX_CONAMED,
                    require_free_or_mutual=require_free_or_mutual,
                    genome_id=ctx.genome_id,
                    seed=seed_for_cands,
                )
                assign_cache[akey] = (rc.top_candidates, rc.gene_call)

        if fill_max is not None and fill_max != float("inf"):
            pw_calls = pathway_calls(
                rxn_calls,
                universal_stoichiometry=ctx.universal_stoichiometry,
                layout=ctx.layout,
                window=LOCUS_WINDOW,
            )
        else:
            # Plan §4 step 5 applies to small specific fills; without the fill gate every reaction is "called" by
            # hundreds of tests and the per-test grouping is meaningless (and quadratic). Reported as no pathway calls.
            pw_calls = []

    (
        n_admitted,
        n_lps,
        cpu_sec,
        cand_rxns,
        t1_fills,
        t1_stable_fills,
        stab_rate,
    ) = _compute_run_cost_and_candidates(
        run=run,
        test_filter=t_filter,
        include_baseline=include_baseline,
        draft_reactions=run_draft_rxns,
        require_stable=require_stable,
        cache_key=cost_cache_key,
    )

    n_called = len(rxn_calls)
    n_gene_rxns = sum(1 for rc in rxn_calls.values() if rc.gene_call is not None)
    n_rxn_only = n_called - n_gene_rxns
    n_cand_rxns = max(len(cand_rxns), n_called)

    gene_abst = float(n_rxn_only) / float(n_called) if n_called > 0 else 0.0
    rxn_abst = float(n_cand_rxns - n_called) / float(n_cand_rxns) if n_cand_rxns > 0 else 0.0

    rxn_rec, gene_rec, pw_rec, rec_fill_size, hidden_rank = _evaluate_masked_recovery(
        run=run,
        rxn_calls=rxn_calls,
        pw_calls=pw_calls,
        ctx=ctx,
        candidates_used=candidates_used,
        seed_num=seed_for_cands,
    )

    if _skip_decisions:
        decisions: Tuple[DistinctDecision, ...] = ()
        n_sup = 0
        n_harm = 0
        n_rel = 0
        n_unc = 0
    else:
        decisions = tuple(
            extract_distinct_decisions(
                genome_id=ctx.genome_id,
                arm=arm,
                rxn_calls=rxn_calls,
                labels_by_protein=active_labels,
                ec_to_reactions=ctx.ec_to_reactions,
            )
        )
        n_sup = sum(1 for d in decisions if d.classification == "supported")
        n_harm = sum(1 for d in decisions if d.classification == "harmful")
        n_rel = sum(1 for d in decisions if d.classification == "harmful" and d.related_substrate)
        n_unc = sum(1 for d in decisions if d.classification == "unconfirmed")

    return ArmRunEvaluation(
        genome_id=ctx.genome_id,
        arm=arm,
        run_kind=run_kind,
        case_id=case_id,
        stratum=stratum,
        reaction_calls=rxn_calls,
        pathway_calls=pw_calls,
        n_admitted_tests=n_admitted,
        n_candidate_reactions=n_cand_rxns,
        n_called_reactions=n_called,
        n_gene_called_reactions=n_gene_rxns,
        n_reaction_only_calls=n_rxn_only,
        n_pathway_calls=len(pw_calls),
        gene_abstention_rate=gene_abst,
        reaction_abstention_rate=rxn_abst,
        n_lps=n_lps,
        cpu_seconds=cpu_sec,
        t1_fill_sizes=t1_fills,
        t1_stable_fill_sizes=t1_stable_fills,
        stability_rate=stab_rate,
        measured_panel_fallback=measured_fallback,
        reaction_recovered=rxn_rec,
        gene_recovered=gene_rec,
        pathway_recovered=pw_rec,
        recovering_fill_size=rec_fill_size,
        hidden_protein_rank=hidden_rank,
        decisions=decisions,
        n_distinct_gene_calls=len(decisions),
        n_supported=n_sup,
        n_harmful=n_harm,
        n_related_substrate=n_rel,
        n_unconfirmed=n_unc,
    )


def evaluate_comparator_q(
    ctx: GenomeScoringContext,
    target_volume: int,
    labels_by_protein: Optional[Mapping[str, Sequence[str]]] = None,
    fixed_threshold: Optional[float] = None,
    arm_name: str = "Q",
) -> ArmRunEvaluation:
    """Evaluate sequence-only comparator ``Q`` (volume-matched or fixed-threshold) on a genome.

    Per protein ``p``, considers only ``p``'s best-scoring EC ``e*`` (ties broken
    by ascending EC).
    - Volume-matched (``fixed_threshold is None``): picks the largest threshold ``t``
      giving at least ``target_volume`` calls in the genome (if ``target_volume == 0``,
      makes 0 calls).
    - Fixed-threshold (``fixed_threshold is not None``): calls ``(p, e*)`` whenever
      ``score(p, e*) >= fixed_threshold``.
    """
    active_labels = (
        labels_by_protein
        if labels_by_protein is not None
        else (ctx.test_labels if ctx.test_labels else ctx.known_labels)
    )
    by_prot = cached_best_by_protein(ctx.candidates, genome_id=ctx.genome_id)
    best_per_prot: List[Tuple[str, str, float]] = []
    for p, ec_list in by_prot.items():
        if ec_list:
            top_ec, top_score = ec_list[0]
            best_per_prot.append((p, top_ec, float(top_score)))
    best_per_prot.sort(key=lambda item: (-item[2], item[0], item[1]))

    if fixed_threshold is not None:
        thresh = float(fixed_threshold)
        selected = [item for item in best_per_prot if item[2] >= thresh - 1e-12]
    else:
        vol = max(0, int(target_volume))
        if vol == 0 or not best_per_prot:
            thresh = 1.0
            selected = []
        else:
            idx = min(vol, len(best_per_prot)) - 1
            thresh = best_per_prot[idx][2]
            selected = [item for item in best_per_prot if item[2] >= thresh - 1e-12]

    decisions: List[DistinctDecision] = []
    for p, ec, s in sorted(selected, key=lambda x: (x[0], x[1])):
        cls, related, truth_ecs = classify_decision(p, ec, active_labels, ctx.ec_to_reactions)
        rxns = tuple(sorted(ctx.ec_to_reactions.get(ec, ())))
        decisions.append(
            DistinctDecision(
                genome_id=ctx.genome_id,
                arm=arm_name,
                protein_id=p,
                ec=ec,
                score=s,
                reactions=rxns,
                classification=cls,
                related_substrate=related,
                truth_ecs=truth_ecs,
            )
        )

    n_sup = sum(1 for d in decisions if d.classification == "supported")
    n_harm = sum(1 for d in decisions if d.classification == "harmful")
    n_rel = sum(1 for d in decisions if d.classification == "harmful" and d.related_substrate)
    n_unc = sum(1 for d in decisions if d.classification == "unconfirmed")

    return ArmRunEvaluation(
        genome_id=ctx.genome_id,
        arm=arm_name,
        run_kind="unmasked",
        case_id=None,
        stratum=None,
        reaction_calls={},
        pathway_calls=[],
        n_admitted_tests=0,
        n_candidate_reactions=0,
        n_called_reactions=len(decisions),
        n_gene_called_reactions=len(decisions),
        n_reaction_only_calls=0,
        n_pathway_calls=0,
        gene_abstention_rate=0.0,
        reaction_abstention_rate=0.0,
        n_lps=0,
        cpu_seconds=0.0,
        decisions=tuple(decisions),
        n_distinct_gene_calls=len(decisions),
        n_supported=n_sup,
        n_harmful=n_harm,
        n_related_substrate=n_rel,
        n_unconfirmed=n_unc,
        q_threshold=thresh,
    )


def find_orphan_steps(
    genome_id: str,
    rxn_calls: Mapping[str, ReactionCall],
    rxn_to_ecs: Mapping[str, Set[str]],
    candidates: Mapping[Tuple[str, str], float],
    score_min: float = SCORE_MIN,
) -> List[OrphanStep]:
    """Return called reactions whose mapped ECs have no candidate ``>= score_min`` anywhere in the genome."""
    idx = candidate_index(candidates, genome_id=genome_id)
    ec_to_prot_scores = idx.ec_to_prot_scores

    orphans: List[OrphanStep] = []
    for rxn_id in sorted(rxn_calls):
        rc = rxn_calls[rxn_id]
        ecs = tuple(sorted(rxn_to_ecs.get(rxn_id, set())))
        max_score = max(
            (ec_to_prot_scores[e][0][1] if e in ec_to_prot_scores and ec_to_prot_scores[e] else 0.0 for e in ecs),
            default=0.0,
        )
        if max_score < score_min - 1e-12:
            orphans.append(
                OrphanStep(
                    genome_id=genome_id,
                    rxn_id=rxn_id,
                    ecs=ecs,
                    max_genome_score=float(max_score),
                    smallest_fill=rc.smallest_fill,
                    n_tests=rc.n_tests,
                    test_types=rc.test_types,
                    calling_tests=rc.calling_tests,
                )
            )
    return orphans


def find_additional_functions(
    n1_decisions: Sequence[DistinctDecision],
    q_decisions: Sequence[DistinctDecision],
) -> List[DistinctDecision]:
    """Return ``N1`` gene calls on labelled proteins whose ``(protein_id, ec)`` is not called by ``Q``."""
    q_pairs = {(d.protein_id, d.ec) for d in q_decisions}
    out: List[DistinctDecision] = []
    for d in n1_decisions:
        if d.classification in ("supported", "harmful") and (d.protein_id, d.ec) not in q_pairs:
            out.append(d)
    return out


def fixed_precision_rule(
    decisions: Sequence[DistinctDecision],
    min_precision: float = FIXED_PRECISION_TARGET,
    min_decisions: int = FIXED_PRECISION_MIN_DECISIONS,
) -> Dict[str, Any]:
    """Evaluate the fixed-precision rule (plan §8): precision >= P on >= n labelled decisions."""
    n_sup = sum(1 for d in decisions if d.classification == "supported")
    n_harm = sum(1 for d in decisions if d.classification == "harmful")
    n_rel = sum(1 for d in decisions if d.classification == "harmful" and d.related_substrate)
    n_unc = sum(1 for d in decisions if d.classification == "unconfirmed")
    n_labelled = n_sup + n_harm
    precision = float(n_sup) / float(n_labelled) if n_labelled > 0 else 0.0
    harm = 1.0 - precision if n_labelled > 0 else 0.0
    passes = (n_labelled >= int(min_decisions)) and (precision >= float(min_precision) - 1e-12)
    return {
        "passes": bool(passes),
        "precision": precision,
        "harm": harm,
        "n_labelled": n_labelled,
        "n_supported": n_sup,
        "n_harmful": n_harm,
        "n_related_substrate": n_rel,
        "n_unconfirmed": n_unc,
        "min_precision": float(min_precision),
        "min_decisions": int(min_decisions),
    }


def _extract_q_threshold_items(
    ctx: GenomeScoringContext,
) -> Tuple[List[Tuple[float, str]], List[Dict[str, Any]]]:
    """Build labelled `(top_score, cls)` pairs for Q threshold search and compact per-protein best-EC records."""
    gid = ctx.genome_id
    pooled_labels: Dict[str, List[str]] = {}
    for src in (ctx.known_labels, ctx.threshold_labels, ctx.test_labels):
        for p, ecs in src.items():
            s_set = set(pooled_labels.get(p, ())) | set(ecs)
            pooled_labels[p] = sorted(s_set)
    active_test_labels = ctx.test_labels if ctx.test_labels else pooled_labels
    active_thresh_labels = ctx.threshold_labels if ctx.threshold_labels else pooled_labels
    thresh_search_labels = ctx.threshold_labels if ctx.threshold_labels else ctx.known_labels

    by_prot = cached_best_by_protein(ctx.candidates, genome_id=gid)
    scored_items: List[Tuple[float, str]] = []
    best_records: List[Dict[str, Any]] = []

    for p in sorted(by_prot.keys()):
        ec_list = by_prot[p]
        if not ec_list:
            continue
        top_ec, top_score = ec_list[0]
        f_score = float(top_score)
        cls_search, _, _ = classify_decision(p, top_ec, thresh_search_labels, ctx.ec_to_reactions)
        if cls_search in ("supported", "harmful"):
            scored_items.append((f_score, cls_search))

        rxns = sorted(ctx.ec_to_reactions.get(top_ec, ()))
        rec: Dict[str, Any] = {
            "p": str(p),
            "ec": str(top_ec),
            "s": f_score,
            "rxns": rxns,
        }
        lbl_info: Dict[str, Any] = {}
        for lset_name, lmap in (
            ("test", active_test_labels),
            ("threshold", active_thresh_labels),
            ("pooled", pooled_labels),
        ):
            cls_l, rel_l, truth_l = classify_decision(p, top_ec, lmap, ctx.ec_to_reactions)
            if cls_l != "unconfirmed" or rel_l or truth_l:
                lbl_info[lset_name] = [cls_l, bool(rel_l), list(truth_l)]
        if lbl_info:
            rec["lbl"] = lbl_info
        best_records.append(rec)

    return scored_items, best_records


def find_q_fixed_precision_threshold(
    contexts: Sequence[Any],
    min_precision: float = FIXED_PRECISION_TARGET,
    min_decisions: int = FIXED_PRECISION_MIN_DECISIONS,
) -> Optional[float]:
    """Find the lowest threshold ``t`` on ``threshold_labels`` (L0) where ``Q(t)`` passes the fixed-precision rule.

    Accepts ``GenomeScoringContext`` objects, ``GenomeScoreResult`` objects, or paths/dicts of phase-1 per-genome files.
    """
    scored_items: List[Tuple[float, str]] = []
    for item in contexts:
        if isinstance(item, GenomeScoringContext):
            lbls = item.threshold_labels if item.threshold_labels else item.known_labels
            by_prot = cached_best_by_protein(item.candidates, genome_id=item.genome_id)
            for p, ec_list in by_prot.items():
                if not ec_list:
                    continue
                top_ec, top_score = ec_list[0]
                cls, _, _ = classify_decision(p, top_ec, lbls, item.ec_to_reactions)
                if cls in ("supported", "harmful"):
                    scored_items.append((float(top_score), cls))
        elif isinstance(item, GenomeScoreResult):
            for s, cls in item.q_threshold_items:
                if cls in ("supported", "harmful"):
                    scored_items.append((float(s), str(cls)))
        else:
            data = load_genome_results_json(item) if isinstance(item, (Path, str)) else item
            if isinstance(data, Mapping):
                for entry in data.get("q_threshold_items") or ():
                    if isinstance(entry, (list, tuple)) and len(entry) >= 2:
                        s, cls = float(entry[0]), str(entry[1])
                        if cls in ("supported", "harmful"):
                            scored_items.append((s, cls))

    if len(scored_items) < min_decisions:
        return None

    thresholds = sorted({s for s, _ in scored_items})
    for t in thresholds:
        subset = [cls for s, cls in scored_items if s >= t - 1e-12]
        n_lab = len(subset)
        if n_lab < min_decisions:
            break
        n_sup = sum(1 for cls in subset if cls == "supported")
        if float(n_sup) / float(n_lab) >= min_precision - 1e-12:
            return float(t)
    return None


def genome_block_bootstrap_recovery(
    case_evals: Sequence[ArmRunEvaluation],
    metric: str = "gene_recovered",
    n_resamples: int = BOOTSTRAP_RESAMPLES,
    seed: int = BOOTSTRAP_SEED,
) -> Dict[str, Any]:
    """Genome-block bootstrap (2,000 resamples, seed 20261102) for masked-case recovery."""
    by_genome: Dict[str, List[int]] = {}
    for ev in case_evals:
        val = 1 if getattr(ev, metric) else 0
        by_genome.setdefault(ev.genome_id, []).append(val)

    total_cases = sum(len(v) for v in by_genome.values())
    total_hits = sum(sum(v) for v in by_genome.values())
    if total_cases == 0:
        return {
            "point": 0.0,
            "ci_low": 0.0,
            "ci_high": 0.0,
            "numerator": 0,
            "denominator": 0,
            "n_genomes": 0,
        }

    point = float(total_hits) / float(total_cases)
    genome_ids = sorted(by_genome.keys())
    k = len(genome_ids)
    if k <= 1:
        return {
            "point": point,
            "ci_low": point,
            "ci_high": point,
            "numerator": total_hits,
            "denominator": total_cases,
            "n_genomes": k,
        }

    rng = random.Random(int(seed))
    samples: List[float] = []
    for _ in range(int(n_resamples)):
        drawn = rng.choices(genome_ids, k=k)
        num = sum(sum(by_genome[g]) for g in drawn)
        den = sum(len(by_genome[g]) for g in drawn)
        if den > 0:
            samples.append(float(num) / float(den))
    samples.sort()
    if not samples:
        return {
            "point": point,
            "ci_low": point,
            "ci_high": point,
            "numerator": total_hits,
            "denominator": total_cases,
            "n_genomes": k,
        }
    ci_low, ci_high = _percentile_ci(samples)
    return {
        "point": point,
        "ci_low": ci_low,
        "ci_high": ci_high,
        "numerator": total_hits,
        "denominator": total_cases,
        "n_genomes": k,
    }


def genome_block_bootstrap_precision(
    decisions: Sequence[DistinctDecision],
    all_genome_ids: Optional[Sequence[str]] = None,
    n_resamples: int = BOOTSTRAP_RESAMPLES,
    seed: int = BOOTSTRAP_SEED,
) -> Dict[str, Any]:
    """Genome-block bootstrap (2,000 resamples, seed 20261102) for distinct-decision precision and harm."""
    by_genome: Dict[str, List[int]] = {}
    if all_genome_ids is not None:
        for g in all_genome_ids:
            by_genome.setdefault(str(g), [])

    for d in decisions:
        if d.classification == "supported":
            by_genome.setdefault(d.genome_id, []).append(1)
        elif d.classification == "harmful":
            by_genome.setdefault(d.genome_id, []).append(0)
        else:
            by_genome.setdefault(d.genome_id, [])

    total_sup = sum(sum(v) for v in by_genome.values())
    total_lab = sum(len(v) for v in by_genome.values())
    total_harm = total_lab - total_sup

    if total_lab == 0:
        return {
            "precision": 0.0,
            "precision_ci_low": 0.0,
            "precision_ci_high": 0.0,
            "harm": 0.0,
            "harm_ci_low": 0.0,
            "harm_ci_high": 0.0,
            "n_supported": 0,
            "n_harmful": 0,
            "n_labelled": 0,
            "n_genomes": len(by_genome),
        }

    point_p = float(total_sup) / float(total_lab)
    point_h = 1.0 - point_p
    genome_ids = sorted(by_genome.keys())
    k = len(genome_ids)
    if k <= 1:
        return {
            "precision": point_p,
            "precision_ci_low": point_p,
            "precision_ci_high": point_p,
            "harm": point_h,
            "harm_ci_low": point_h,
            "harm_ci_high": point_h,
            "n_supported": total_sup,
            "n_harmful": total_harm,
            "n_labelled": total_lab,
            "n_genomes": k,
        }

    rng = random.Random(int(seed))
    p_samples: List[float] = []
    for _ in range(int(n_resamples)):
        drawn = rng.choices(genome_ids, k=k)
        sup = sum(sum(by_genome[g]) for g in drawn)
        lab = sum(len(by_genome[g]) for g in drawn)
        if lab > 0:
            p_samples.append(float(sup) / float(lab))
    p_samples.sort()
    if not p_samples:
        p_low, p_high = point_p, point_p
    else:
        p_low, p_high = _percentile_ci(p_samples)

    return {
        "precision": point_p,
        "precision_ci_low": p_low,
        "precision_ci_high": p_high,
        "harm": point_h,
        "harm_ci_low": 1.0 - p_high,
        "harm_ci_high": 1.0 - p_low,
        "n_supported": total_sup,
        "n_harmful": total_harm,
        "n_labelled": total_lab,
        "n_genomes": k,
    }


def paired_mcnemar_comparison(
    evals_better: Sequence[ArmRunEvaluation],
    evals_baseline: Sequence[ArmRunEvaluation],
    metric: str = "gene_recovered",
) -> Dict[str, Any]:
    """Paired exact one-sided McNemar test between two arms across matched masked cases."""
    map_better = {(e.genome_id, e.case_id): bool(getattr(e, metric)) for e in evals_better if e.case_id}
    map_base = {(e.genome_id, e.case_id): bool(getattr(e, metric)) for e in evals_baseline if e.case_id}
    common_keys = sorted(set(map_better.keys()) & set(map_base.keys()))

    b_better_only = 0
    c_base_only = 0
    both_true = 0
    both_false = 0
    for k in common_keys:
        v_bet = map_better[k]
        v_bas = map_base[k]
        if v_bet and not v_bas:
            b_better_only += 1
        elif v_bas and not v_bet:
            c_base_only += 1
        elif v_bet and v_bas:
            both_true += 1
        else:
            both_false += 1

    p_val = exact_one_sided_mcnemar_p(b_better_only, c_base_only)
    return {
        "n_cases": len(common_keys),
        "b_better_only": b_better_only,
        "c_baseline_only": c_base_only,
        "both_recovered": both_true,
        "neither_recovered": both_false,
        "p_value": p_val,
    }


@dataclass
class GenomeScoreResult:
    genome_id: str
    candidate_table_kind: str
    context: Optional[GenomeScoringContext] = None
    masked_evals_by_arm: Dict[str, List[ArmRunEvaluation]] = field(default_factory=dict)
    unmasked_eval_by_arm: Dict[str, ArmRunEvaluation] = field(default_factory=dict)
    threshold_decisions_by_arm: Dict[str, List[DistinctDecision]] = field(default_factory=dict)
    test_decisions_by_arm: Dict[str, List[DistinctDecision]] = field(default_factory=dict)
    pooled_decisions_by_arm: Dict[str, List[DistinctDecision]] = field(default_factory=dict)
    orphan_steps: List[OrphanStep] = field(default_factory=list)
    additional_functions: List[DistinctDecision] = field(default_factory=list)
    q_threshold_items: List[Tuple[float, str]] = field(default_factory=list)
    q_best_per_protein: List[Dict[str, Any]] = field(default_factory=list)


def load_genome_results_json(path: Union[Path, str]) -> Dict[str, Any]:
    """Load ``genome_results2.json.gz`` or ``.json``."""
    p = Path(path)
    if p.name.endswith(".gz"):
        with gzip.open(p, "rt", encoding="utf-8") as fh:
            return json.load(fh)
    with p.open("r", encoding="utf-8") as fh:
        return json.load(fh)


def _decision_to_dict(d: DistinctDecision) -> Dict[str, Any]:
    return {
        "genome_id": d.genome_id,
        "arm": d.arm,
        "protein_id": d.protein_id,
        "ec": d.ec,
        "score": float(d.score),
        "reactions": list(d.reactions),
        "classification": d.classification,
        "related_substrate": bool(d.related_substrate),
        "truth_ecs": list(d.truth_ecs),
    }


def _decision_from_dict(d: Mapping[str, Any]) -> DistinctDecision:
    return DistinctDecision(
        genome_id=str(d["genome_id"]),
        arm=str(d["arm"]),
        protein_id=str(d["protein_id"]),
        ec=str(d["ec"]),
        score=float(d["score"]),
        reactions=tuple(str(r) for r in (d.get("reactions") or ())),
        classification=str(d["classification"]),
        related_substrate=bool(d["related_substrate"]),
        truth_ecs=tuple(str(e) for e in (d.get("truth_ecs") or ())),
    )


def _pathway_call_to_dict(pc: PathwayCall) -> Dict[str, Any]:
    return {
        "test_id": pc.test_id,
        "reactions": list(pc.reactions),
        "genes": list(pc.genes),
        "cluster_span": int(pc.cluster_span),
        "contig": pc.contig,
    }


def _pathway_call_from_dict(d: Mapping[str, Any]) -> PathwayCall:
    return PathwayCall(
        test_id=str(d["test_id"]),
        reactions=tuple(str(r) for r in (d.get("reactions") or ())),
        genes=tuple(str(g) for g in (d.get("genes") or ())),
        cluster_span=int(d["cluster_span"]),
        contig=str(d.get("contig") or ""),
    )


def _orphan_step_to_dict(os_entry: OrphanStep) -> Dict[str, Any]:
    return {
        "genome_id": os_entry.genome_id,
        "rxn_id": os_entry.rxn_id,
        "ecs": list(os_entry.ecs),
        "max_genome_score": float(os_entry.max_genome_score),
        "smallest_fill": int(os_entry.smallest_fill),
        "n_tests": int(os_entry.n_tests),
        "test_types": list(os_entry.test_types),
        "calling_tests": list(os_entry.calling_tests),
    }


def _orphan_step_from_dict(d: Mapping[str, Any]) -> OrphanStep:
    return OrphanStep(
        genome_id=str(d["genome_id"]),
        rxn_id=str(d["rxn_id"]),
        ecs=tuple(str(e) for e in (d.get("ecs") or ())),
        max_genome_score=float(d["max_genome_score"]),
        smallest_fill=int(d["smallest_fill"]),
        n_tests=int(d["n_tests"]),
        test_types=tuple(str(t) for t in (d.get("test_types") or ())),
        calling_tests=tuple(str(t) for t in (d.get("calling_tests") or ())),
    )


def _arm_eval_to_dict(ev: ArmRunEvaluation) -> Dict[str, Any]:
    return {
        "genome_id": ev.genome_id,
        "arm": ev.arm,
        "run_kind": ev.run_kind,
        "case_id": ev.case_id,
        "stratum": ev.stratum,
        "pathway_calls": [_pathway_call_to_dict(pc) for pc in ev.pathway_calls],
        "n_admitted_tests": int(ev.n_admitted_tests),
        "n_candidate_reactions": int(ev.n_candidate_reactions),
        "n_called_reactions": int(ev.n_called_reactions),
        "n_gene_called_reactions": int(ev.n_gene_called_reactions),
        "n_reaction_only_calls": int(ev.n_reaction_only_calls),
        "n_pathway_calls": int(ev.n_pathway_calls),
        "gene_abstention_rate": float(ev.gene_abstention_rate),
        "reaction_abstention_rate": float(ev.reaction_abstention_rate),
        "n_lps": int(ev.n_lps),
        "cpu_seconds": float(ev.cpu_seconds),
        "t1_fill_sizes": [int(x) for x in ev.t1_fill_sizes],
        "t1_stable_fill_sizes": [int(x) for x in ev.t1_stable_fill_sizes],
        "stability_rate": float(ev.stability_rate),
        "measured_panel_fallback": bool(ev.measured_panel_fallback),
        "reaction_recovered": bool(ev.reaction_recovered),
        "gene_recovered": bool(ev.gene_recovered),
        "pathway_recovered": bool(ev.pathway_recovered),
        "recovering_fill_size": None if ev.recovering_fill_size is None else int(ev.recovering_fill_size),
        "hidden_protein_rank": None if ev.hidden_protein_rank is None else int(ev.hidden_protein_rank),
        "n_distinct_gene_calls": int(ev.n_distinct_gene_calls),
        "n_supported": int(ev.n_supported),
        "n_harmful": int(ev.n_harmful),
        "n_related_substrate": int(ev.n_related_substrate),
        "n_unconfirmed": int(ev.n_unconfirmed),
        "q_threshold": None if ev.q_threshold is None else float(ev.q_threshold),
    }


def _arm_eval_from_dict(
    d: Mapping[str, Any],
    decisions: Tuple[DistinctDecision, ...] = (),
) -> ArmRunEvaluation:
    return ArmRunEvaluation(
        genome_id=str(d["genome_id"]),
        arm=str(d["arm"]),
        run_kind=str(d["run_kind"]),
        case_id=d.get("case_id"),
        stratum=d.get("stratum"),
        reaction_calls={},
        pathway_calls=[_pathway_call_from_dict(pc) for pc in (d.get("pathway_calls") or ())],
        n_admitted_tests=int(d["n_admitted_tests"]),
        n_candidate_reactions=int(d["n_candidate_reactions"]),
        n_called_reactions=int(d["n_called_reactions"]),
        n_gene_called_reactions=int(d["n_gene_called_reactions"]),
        n_reaction_only_calls=int(d["n_reaction_only_calls"]),
        n_pathway_calls=int(d["n_pathway_calls"]),
        gene_abstention_rate=float(d["gene_abstention_rate"]),
        reaction_abstention_rate=float(d["reaction_abstention_rate"]),
        n_lps=int(d["n_lps"]),
        cpu_seconds=float(d["cpu_seconds"]),
        t1_fill_sizes=tuple(int(x) for x in (d.get("t1_fill_sizes") or ())),
        t1_stable_fill_sizes=tuple(int(x) for x in (d.get("t1_stable_fill_sizes") or ())),
        stability_rate=float(d.get("stability_rate") or 0.0),
        measured_panel_fallback=bool(d.get("measured_panel_fallback", False)),
        reaction_recovered=bool(d.get("reaction_recovered", False)),
        gene_recovered=bool(d.get("gene_recovered", False)),
        pathway_recovered=bool(d.get("pathway_recovered", False)),
        recovering_fill_size=None if d.get("recovering_fill_size") is None else int(d["recovering_fill_size"]),
        hidden_protein_rank=None if d.get("hidden_protein_rank") is None else int(d["hidden_protein_rank"]),
        decisions=decisions,
        n_distinct_gene_calls=int(d.get("n_distinct_gene_calls") or 0),
        n_supported=int(d.get("n_supported") or 0),
        n_harmful=int(d.get("n_harmful") or 0),
        n_related_substrate=int(d.get("n_related_substrate") or 0),
        n_unconfirmed=int(d.get("n_unconfirmed") or 0),
        q_threshold=None if d.get("q_threshold") is None else float(d["q_threshold"]),
    )


def genome_score_result_to_dict(gr: GenomeScoreResult) -> Dict[str, Any]:
    """Serialize a `GenomeScoreResult` into a compact JSON-compatible dictionary."""
    return {
        "schema": "gspa.offgrowth2.per_genome_score/1",
        "genome_id": gr.genome_id,
        "candidate_table_kind": gr.candidate_table_kind,
        "masked_evals_by_arm": {
            arm: [_arm_eval_to_dict(ev) for ev in evs]
            for arm, evs in gr.masked_evals_by_arm.items()
        },
        "unmasked_eval_by_arm": {
            arm: _arm_eval_to_dict(ev)
            for arm, ev in gr.unmasked_eval_by_arm.items()
        },
        "threshold_decisions_by_arm": {
            arm: [_decision_to_dict(d) for d in decs]
            for arm, decs in gr.threshold_decisions_by_arm.items()
        },
        "test_decisions_by_arm": {
            arm: [_decision_to_dict(d) for d in decs]
            for arm, decs in gr.test_decisions_by_arm.items()
        },
        "pooled_decisions_by_arm": {
            arm: [_decision_to_dict(d) for d in decs]
            for arm, decs in gr.pooled_decisions_by_arm.items()
        },
        "orphan_steps": [_orphan_step_to_dict(o) for o in gr.orphan_steps],
        "additional_functions": [_decision_to_dict(d) for d in gr.additional_functions],
        "q_threshold_items": [[float(s), str(c)] for s, c in gr.q_threshold_items],
        "q_best_per_protein": gr.q_best_per_protein,
    }


def genome_score_result_from_dict(data: Mapping[str, Any]) -> GenomeScoreResult:
    """Deserialize a `GenomeScoreResult` from its compact JSON dictionary."""
    thresh_decs = {
        arm: [_decision_from_dict(d) for d in decs]
        for arm, decs in (data.get("threshold_decisions_by_arm") or {}).items()
    }
    test_decs = {
        arm: [_decision_from_dict(d) for d in decs]
        for arm, decs in (data.get("test_decisions_by_arm") or {}).items()
    }
    pooled_decs = {
        arm: [_decision_from_dict(d) for d in decs]
        for arm, decs in (data.get("pooled_decisions_by_arm") or {}).items()
    }
    masked_evals = {
        arm: [_arm_eval_from_dict(ev) for ev in evs]
        for arm, evs in (data.get("masked_evals_by_arm") or {}).items()
    }
    unmasked_evals = {
        arm: _arm_eval_from_dict(ev, decisions=tuple(test_decs.get(arm, ())))
        for arm, ev in (data.get("unmasked_eval_by_arm") or {}).items()
    }
    return GenomeScoreResult(
        genome_id=str(data["genome_id"]),
        candidate_table_kind=str(data.get("candidate_table_kind") or "six"),
        context=None,
        masked_evals_by_arm=masked_evals,
        unmasked_eval_by_arm=unmasked_evals,
        threshold_decisions_by_arm=thresh_decs,
        test_decisions_by_arm=test_decs,
        pooled_decisions_by_arm=pooled_decs,
        orphan_steps=[_orphan_step_from_dict(o) for o in (data.get("orphan_steps") or ())],
        additional_functions=[_decision_from_dict(d) for d in (data.get("additional_functions") or ())],
        q_threshold_items=[(float(x[0]), str(x[1])) for x in (data.get("q_threshold_items") or ())],
        q_best_per_protein=list(data.get("q_best_per_protein") or ()),
    )


def save_genome_score_result(gr: GenomeScoreResult, path: Union[Path, str]) -> Path:
    """Atomically write a compact `GenomeScoreResult` to `<out-dir>/per-genome/<genome>.json.gz`."""
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = target.parent / f".{target.name}.tmp.{os.getpid()}"
    payload = genome_score_result_to_dict(gr)
    if target.name.endswith(".gz"):
        with gzip.open(tmp_path, "wt", encoding="utf-8") as fh:
            json.dump(payload, fh, sort_keys=True)
    else:
        with tmp_path.open("w", encoding="utf-8") as fh:
            json.dump(payload, fh, sort_keys=True)
    os.replace(tmp_path, target)
    return target


def load_genome_score_result(path: Union[Path, str]) -> GenomeScoreResult:
    """Load a compact `GenomeScoreResult` from `<out-dir>/per-genome/<genome>.json.gz`."""
    return genome_score_result_from_dict(load_genome_results_json(path))


def apply_q_fixed_threshold_to_result(
    gr: GenomeScoreResult,
    q_fixed_threshold: Optional[float],
) -> None:
    """Populate or update `Q_fixed` and `additional_functions` on `gr` from its stored `q_best_per_protein`."""
    if q_fixed_threshold is None:
        gr.unmasked_eval_by_arm.pop("Q_fixed", None)
        gr.test_decisions_by_arm.pop("Q_fixed", None)
        gr.threshold_decisions_by_arm.pop("Q_fixed", None)
        gr.pooled_decisions_by_arm.pop("Q_fixed", None)
        if "N1" in gr.unmasked_eval_by_arm:
            gr.additional_functions = find_additional_functions(
                gr.pooled_decisions_by_arm.get("N1", []),
                gr.pooled_decisions_by_arm.get("Q", []),
            )
        return

    thresh = float(q_fixed_threshold)
    existing = gr.unmasked_eval_by_arm.get("Q_fixed")
    if existing is not None and existing.q_threshold is not None and abs(existing.q_threshold - thresh) <= 1e-12:
        return

    if not gr.q_best_per_protein:
        return

    gid = gr.genome_id
    selected = [rec for rec in gr.q_best_per_protein if float(rec["s"]) >= thresh - 1e-12]
    selected.sort(key=lambda rec: (str(rec["p"]), str(rec["ec"])))

    decs_by_set: Dict[str, List[DistinctDecision]] = {"test": [], "threshold": [], "pooled": []}
    for lset_name in ("test", "threshold", "pooled"):
        for rec in selected:
            lbl = (rec.get("lbl") or {}).get(lset_name)
            if lbl is not None:
                cls_l = str(lbl[0])
                rel_l = bool(lbl[1])
                truth_l = tuple(str(e) for e in (lbl[2] or ()))
            else:
                cls_l = "unconfirmed"
                rel_l = False
                truth_l = ()
            decs_by_set[lset_name].append(
                DistinctDecision(
                    genome_id=gid,
                    arm="Q_fixed",
                    protein_id=str(rec["p"]),
                    ec=str(rec["ec"]),
                    score=float(rec["s"]),
                    reactions=tuple(str(r) for r in (rec.get("rxns") or ())),
                    classification=cls_l,
                    related_substrate=rel_l,
                    truth_ecs=truth_l,
                )
            )

    test_decs = decs_by_set["test"]
    n_sup = sum(1 for d in test_decs if d.classification == "supported")
    n_harm = sum(1 for d in test_decs if d.classification == "harmful")
    n_rel = sum(1 for d in test_decs if d.classification == "harmful" and d.related_substrate)
    n_unc = sum(1 for d in test_decs if d.classification == "unconfirmed")

    gr.unmasked_eval_by_arm["Q_fixed"] = ArmRunEvaluation(
        genome_id=gid,
        arm="Q_fixed",
        run_kind="unmasked",
        case_id=None,
        stratum=None,
        reaction_calls={},
        pathway_calls=[],
        n_admitted_tests=0,
        n_candidate_reactions=0,
        n_called_reactions=len(test_decs),
        n_gene_called_reactions=len(test_decs),
        n_reaction_only_calls=0,
        n_pathway_calls=0,
        gene_abstention_rate=0.0,
        reaction_abstention_rate=0.0,
        n_lps=0,
        cpu_seconds=0.0,
        decisions=tuple(test_decs),
        n_distinct_gene_calls=len(test_decs),
        n_supported=n_sup,
        n_harmful=n_harm,
        n_related_substrate=n_rel,
        n_unconfirmed=n_unc,
        q_threshold=thresh,
    )
    gr.test_decisions_by_arm["Q_fixed"] = test_decs
    gr.threshold_decisions_by_arm["Q_fixed"] = decs_by_set["threshold"]
    gr.pooled_decisions_by_arm["Q_fixed"] = decs_by_set["pooled"]
    if "N1" in gr.unmasked_eval_by_arm:
        gr.additional_functions = find_additional_functions(
            gr.pooled_decisions_by_arm.get("N1", []),
            gr.pooled_decisions_by_arm["Q_fixed"],
        )


def score_genome(
    genome_payload: Mapping[str, Any],
    ctx: Optional[GenomeScoringContext] = None,
    arms: Sequence[str] = PRIMARY_ARMS,
    q_fixed_threshold: Optional[float] = None,
) -> GenomeScoreResult:
    """Score all masked and unmasked runs of a single genome across all requested arms and ``Q``."""
    if ctx is None:
        ctx = build_genome_context(genome_payload)

    gid = ctx.genome_id
    cand_kind = str(genome_payload.get("candidate_table_kind") or "six")

    pooled_labels: Dict[str, List[str]] = {}
    for src in (ctx.known_labels, ctx.threshold_labels, ctx.test_labels):
        for p, ecs in src.items():
            s_set = set(pooled_labels.get(p, ())) | set(ecs)
            pooled_labels[p] = sorted(s_set)

    active_test_labels = ctx.test_labels if ctx.test_labels else pooled_labels
    active_thresh_labels = ctx.threshold_labels if ctx.threshold_labels else pooled_labels

    masked_by_arm: Dict[str, List[ArmRunEvaluation]] = {a: [] for a in arms}
    unmasked_by_arm: Dict[str, ArmRunEvaluation] = {}
    thresh_dec_by_arm: Dict[str, List[DistinctDecision]] = {}
    test_dec_by_arm: Dict[str, List[DistinctDecision]] = {}
    pooled_dec_by_arm: Dict[str, List[DistinctDecision]] = {}
    n1_unmasked_rxn_calls: Dict[str, ReactionCall] = {}

    runs = list(genome_payload.get("runs") or ())
    for run in runs:
        r_kind = str(run.get("run_kind") or "unmasked")
        skip_dec = r_kind == "masked"
        # Evaluate N1 first so S1b has N1's exact LP budget for this run
        n1_eval: Optional[ArmRunEvaluation] = None
        if "N1" in arms or "S1b" in arms:
            n1_eval = evaluate_arm_on_run(
                run, "N1", ctx, labels_by_protein=active_test_labels, _skip_decisions=skip_dec
            )

        for arm in arms:
            if arm == "N1" and n1_eval is not None:
                ev = n1_eval
            else:
                ev = evaluate_arm_on_run(
                    run,
                    arm,
                    ctx,
                    labels_by_protein=active_test_labels,
                    n1_lp_budget=(n1_eval.n_lps if n1_eval is not None else None),
                    _skip_decisions=skip_dec,
                )
            if r_kind == "masked":
                ev.reaction_calls = {}
                masked_by_arm[arm].append(ev)
            else:
                if arm == "N1":
                    n1_unmasked_rxn_calls = ev.reaction_calls
                thresh_dec_by_arm[arm] = extract_distinct_decisions(
                    gid, arm, ev.reaction_calls, active_thresh_labels, ctx.ec_to_reactions
                )
                test_dec_by_arm[arm] = list(ev.decisions)
                pooled_dec_by_arm[arm] = extract_distinct_decisions(
                    gid, arm, ev.reaction_calls, pooled_labels, ctx.ec_to_reactions
                )
                ev.reaction_calls = {}
                unmasked_by_arm[arm] = ev

    # Evaluate comparator Q (volume-matched to N1's distinct gene calls on the unmasked run)
    n1_unmasked = unmasked_by_arm.get("N1")
    target_vol = n1_unmasked.n_distinct_gene_calls if n1_unmasked is not None else 0

    q_eval = evaluate_comparator_q(
        ctx,
        target_volume=target_vol,
        labels_by_protein=active_test_labels,
        fixed_threshold=None,
        arm_name="Q",
    )
    unmasked_by_arm["Q"] = q_eval
    test_dec_by_arm["Q"] = list(q_eval.decisions)
    thresh_dec_by_arm["Q"] = list(
        evaluate_comparator_q(
            ctx, target_volume=target_vol, labels_by_protein=active_thresh_labels, fixed_threshold=None, arm_name="Q"
        ).decisions
    )
    pooled_dec_by_arm["Q"] = list(
        evaluate_comparator_q(
            ctx, target_volume=target_vol, labels_by_protein=pooled_labels, fixed_threshold=None, arm_name="Q"
        ).decisions
    )

    if q_fixed_threshold is not None:
        qf_eval = evaluate_comparator_q(
            ctx,
            target_volume=target_vol,
            labels_by_protein=active_test_labels,
            fixed_threshold=q_fixed_threshold,
            arm_name="Q_fixed",
        )
        unmasked_by_arm["Q_fixed"] = qf_eval
        test_dec_by_arm["Q_fixed"] = list(qf_eval.decisions)
        thresh_dec_by_arm["Q_fixed"] = list(
            evaluate_comparator_q(
                ctx,
                target_volume=target_vol,
                labels_by_protein=active_thresh_labels,
                fixed_threshold=q_fixed_threshold,
                arm_name="Q_fixed",
            ).decisions
        )
        pooled_dec_by_arm["Q_fixed"] = list(
            evaluate_comparator_q(
                ctx,
                target_volume=target_vol,
                labels_by_protein=pooled_labels,
                fixed_threshold=q_fixed_threshold,
                arm_name="Q_fixed",
            ).decisions
        )

    # Orphan steps and additional functions from unmasked N1
    if n1_unmasked is not None:
        orphans = find_orphan_steps(
            gid,
            n1_unmasked_rxn_calls,
            ctx.rxn_to_ecs,
            ctx.candidates,
            score_min=SCORE_MIN,
        )
        q_ref_decisions = (
            pooled_dec_by_arm.get("Q_fixed")
            if "Q_fixed" in pooled_dec_by_arm
            else pooled_dec_by_arm.get("Q", [])
        )
        add_funcs = find_additional_functions(
            pooled_dec_by_arm.get("N1", []),
            q_ref_decisions,
        )
    else:
        orphans = []
        add_funcs = []

    q_thresh_items, q_best_records = _extract_q_threshold_items(ctx)
    ctx._assign_gene_cache.clear()

    return GenomeScoreResult(
        genome_id=gid,
        candidate_table_kind=cand_kind,
        context=None,
        masked_evals_by_arm=masked_by_arm,
        unmasked_eval_by_arm=unmasked_by_arm,
        threshold_decisions_by_arm=thresh_dec_by_arm,
        test_decisions_by_arm=test_dec_by_arm,
        pooled_decisions_by_arm=pooled_dec_by_arm,
        orphan_steps=orphans,
        additional_functions=add_funcs,
        q_threshold_items=q_thresh_items,
        q_best_per_protein=q_best_records,
    )
