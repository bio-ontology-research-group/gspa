"""Call rules for the second off-growth evaluation (plan §4 steps 3–5, §6).

Module-level frozen constants from the plan:
- ``FILL_MAX = 3`` (sensitivity values 1 and 5)
- ``SCORE_MIN = 0.30``
- ``MARGIN_MIN = 0.05``
- ``TIE_WINDOW = 0.02``
- ``MAX_CONAMED = 2``
- ``LOCUS_WINDOW = 10``
- ``CURRENCY_METABOLITES``: explicit SEED compound IDs for
  ``h, h2o, atp, adp, amp, pi, ppi, nad(p)(h), fad(h2), coa, co2, o2, nh4, so4``.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import math
import re
import sys
from pathlib import Path
from typing import (
    Any,
    Callable,
    Dict,
    FrozenSet,
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
from offgrowth.mask import protein_aliases, strip_version  # noqa: E402

# Frozen plan constants (§4 steps 3–5)
FILL_MAX: int = 3          # F; sensitivity values 1 and 5
SCORE_MIN: float = 0.30    # s_min
MARGIN_MIN: float = 0.05   # delta
TIE_WINDOW: float = 0.02   # co-naming window
MAX_CONAMED: int = 2
LOCUS_WINDOW: int = 10     # W

# Explicit SEED compound IDs for currency metabolites:
# h, h2o, atp, adp, amp, pi, ppi, nad, nadh, nadp, nadph, fad, fadh2, coa, co2, o2, nh4, so4
CURRENCY_METABOLITES: FrozenSet[str] = frozenset({
    "cpd00067",  # H+ (h)
    "cpd00001",  # H2O (h2o)
    "cpd00002",  # ATP (atp)
    "cpd00008",  # ADP (adp)
    "cpd00018",  # AMP (amp)
    "cpd00009",  # Phosphate (pi)
    "cpd00012",  # Pyrophosphate (ppi)
    "cpd00003",  # NAD+ (nad)
    "cpd00004",  # NADH (nadh)
    "cpd00006",  # NADP+ (nadp)
    "cpd00005",  # NADPH (nadph)
    "cpd00015",  # FAD (fad)
    "cpd00982",  # FADH2 (fadh2)
    "cpd00010",  # CoA (coa)
    "cpd00011",  # CO2 (co2)
    "cpd00007",  # O2 (o2)
    "cpd00013",  # NH3 / NH4+ (nh4)
    "cpd00048",  # Sulfate (so4)
})

CURRENCY_SHORT_NAMES: FrozenSet[str] = frozenset({
    "h", "h+", "h2o", "atp", "adp", "amp", "pi", "ppi",
    "nad", "nad+", "nadh", "nadp", "nadp+", "nadph",
    "fad", "fadh2", "coa", "co2", "o2", "nh4", "nh4+", "nh3", "so4", "so42-",
})

_COMPARTMENT_SUFFIX_RE = re.compile(r"(?:_[a-zA-Z]\d*|\[[a-zA-Z]\d*\])$")


def base_cpd_id(cpd_id: str) -> str:
    """Strip compartment suffix (e.g. ``_c0``, ``_e0``, ``[c0]``) from a compound identifier."""
    s = (cpd_id or "").strip()
    return _COMPARTMENT_SUFFIX_RE.sub("", s)


def is_currency_metabolite(cpd_id: str) -> bool:
    """Return True when ``cpd_id`` is in the frozen currency metabolite list."""
    base = base_cpd_id(cpd_id)
    return base in CURRENCY_METABOLITES or base.lower() in CURRENCY_SHORT_NAMES


@dataclass(frozen=True)
class CandidateRankEntry:
    protein_id: str
    ec: str
    score: float
    margin: float = 0.0
    is_free: bool = True
    is_mutual_best: bool = False


@dataclass(frozen=True)
class GeneCall:
    rxn_id: str
    proteins: Tuple[str, ...]
    ecs: Tuple[str, ...]
    scores: Tuple[float, ...]
    margins: Tuple[float, ...] = ()
    protein_ec_pairs: Tuple[Tuple[str, str, float], ...] = ()

    @property
    def protein_id(self) -> str:
        """Primary (highest-scoring) named protein."""
        return self.proteins[0] if self.proteins else ""

    @property
    def ec(self) -> str:
        """EC through which the primary named protein was assigned."""
        return self.ecs[0] if self.ecs else ""

    @property
    def score(self) -> float:
        """Score of the primary named protein."""
        return self.scores[0] if self.scores else 0.0

    @property
    def margin(self) -> float:
        """Margin of the primary named protein over the best non-tied competitor."""
        return self.margins[0] if self.margins else 0.0

    @property
    def is_conamed(self) -> bool:
        return len(self.proteins) > 1


@dataclass
class ReactionCall:
    rxn_id: str
    calling_tests: Tuple[str, ...]
    smallest_fill: int
    n_tests: int
    test_types: Tuple[str, ...]
    all_tests: Tuple[str, ...] = ()
    gene_call: Optional[GeneCall] = None
    top_candidates: Tuple[CandidateRankEntry, ...] = field(default_factory=tuple)

    @property
    def is_reaction_only(self) -> bool:
        return self.gene_call is None


@dataclass(frozen=True)
class PathwayCall:
    test_id: str
    reactions: Tuple[str, ...]
    genes: Tuple[str, ...]
    cluster_span: int
    contig: str = ""


def _admit_test(
    test: Mapping[str, Any],
    test_filter: Optional[Union[Callable[[Mapping[str, Any]], bool], Iterable[str]]],
) -> bool:
    if test_filter is None:
        return True
    if callable(test_filter):
        return bool(test_filter(test))
    allowed = {str(x) for x in test_filter}
    return str(test.get("test_id", "")) in allowed


_BASE_RXN_CACHE: Dict[str, str] = {}
_NORM_EC_CACHE: Dict[str, str] = {}
_EXCLUDED_RXN_CACHE: Dict[Tuple[int, int, bool], Tuple[Any, Any, Set[str]]] = {}


def _fast_base_rxn_id(rxn_id: Any) -> str:
    s = rxn_id if isinstance(rxn_id, str) else str(rxn_id)
    hit = _BASE_RXN_CACHE.get(s)
    if hit is not None:
        return hit
    res = base_seed_rxn_id(s)
    if len(_BASE_RXN_CACHE) > 200000:
        _BASE_RXN_CACHE.clear()
    _BASE_RXN_CACHE[s] = res
    return res


def _fast_normalize_ec(ec: Any) -> str:
    s = ec if isinstance(ec, str) else str(ec)
    hit = _NORM_EC_CACHE.get(s)
    if hit is not None:
        return hit
    res = normalize_ec(s)
    if len(_NORM_EC_CACHE) > 100000:
        _NORM_EC_CACHE.clear()
    _NORM_EC_CACHE[s] = res
    return res


def _extract_excluded_reactions(
    run: Mapping[str, Any],
    draft_reactions: Optional[Iterable[str]] = None,
    exclude_baseline: bool = True,
) -> Set[str]:
    cacheable = draft_reactions is None or isinstance(draft_reactions, (set, frozenset))
    if cacheable:
        ckey = (id(run), id(draft_reactions), bool(exclude_baseline))
        chit = _EXCLUDED_RXN_CACHE.get(ckey)
        if chit is not None and chit[0] is run and chit[1] is draft_reactions:
            return chit[2]

    excluded: Set[str] = set()
    if draft_reactions is not None:
        for r in draft_reactions:
            excluded.add(_fast_base_rxn_id(r))
    for r in run.get("draft_reactions", ()) or ():
        excluded.add(_fast_base_rxn_id(r))

    if exclude_baseline:
        bf = run.get("baseline_fill")
        if isinstance(bf, Mapping):
            for key in ("fill", "stable_fill", "replay_fill"):
                for r in bf.get(key) or ():
                    excluded.add(_fast_base_rxn_id(r))
        elif isinstance(bf, (list, tuple, set, frozenset)):
            for r in bf:
                excluded.add(_fast_base_rxn_id(r))

    if cacheable:
        if len(_EXCLUDED_RXN_CACHE) > 64:
            _EXCLUDED_RXN_CACHE.clear()
        _EXCLUDED_RXN_CACHE[(id(run), id(draft_reactions), bool(exclude_baseline))] = (run, draft_reactions, excluded)
    return excluded


def reaction_calls(
    run: Mapping[str, Any],
    fill_max: Optional[float] = FILL_MAX,
    require_stable: bool = True,
    test_filter: Optional[Union[Callable[[Mapping[str, Any]], bool], Iterable[str]]] = None,
    draft_reactions: Optional[Iterable[str]] = None,
) -> Dict[str, ReactionCall]:
    """Apply the reaction gate (plan §4 step 3) to a single run.

    A reaction ``r`` is called when at least one test admitted by ``test_filter``
    adds ``r`` with ``len(active_fill) <= fill_max`` and ``r in active_fill``,
    where ``active_fill`` is ``stable_fill`` when ``require_stable=True`` and
    ``fill`` when ``require_stable=False``. Reactions in ``draft_reactions`` or
    in ``baseline_fill`` are never called.
    """
    excluded = _extract_excluded_reactions(run, draft_reactions=draft_reactions, exclude_baseline=True)
    max_fill = float("inf") if fill_max is None else float(fill_max)

    calling_tests_by_rxn: Dict[str, List[str]] = {}
    all_tests_by_rxn: Dict[str, Set[str]] = {}
    smallest_fill_by_rxn: Dict[str, int] = {}
    test_types_by_rxn: Dict[str, Set[str]] = {}

    for idx, test in enumerate(run.get("tests", ()) or ()):
        if not _admit_test(test, test_filter):
            continue
        test_id = str(test.get("test_id") or f"test_{idx}")
        test_type = str(test.get("test_type") or "")

        if require_stable:
            raw_fill = test.get("stable_fill")
            if raw_fill is None:
                raw_fill = ()
        else:
            raw_fill = test.get("fill") or ()

        norm_fill: List[str] = []
        seen_in_test: Set[str] = set()
        for r in raw_fill:
            br = _fast_base_rxn_id(r)
            if br and br not in seen_in_test:
                seen_in_test.add(br)
                norm_fill.append(br)

        fill_len = len(norm_fill)
        if fill_len == 0:
            continue

        for br in norm_fill:
            if br in excluded:
                continue
            all_tests_by_rxn.setdefault(br, set()).add(test_id)

        if fill_len <= max_fill:
            for br in norm_fill:
                if br in excluded:
                    continue
                calling_tests_by_rxn.setdefault(br, [])
                if test_id not in calling_tests_by_rxn[br]:
                    calling_tests_by_rxn[br].append(test_id)
                cur_min = smallest_fill_by_rxn.get(br)
                if cur_min is None or fill_len < cur_min:
                    smallest_fill_by_rxn[br] = fill_len
                if test_type:
                    test_types_by_rxn.setdefault(br, set()).add(test_type)

    out: Dict[str, ReactionCall] = {}
    for br in sorted(calling_tests_by_rxn):
        c_tests = tuple(calling_tests_by_rxn[br])
        a_tests = tuple(sorted(all_tests_by_rxn.get(br, set())))
        t_types = tuple(sorted(test_types_by_rxn.get(br, set())))
        out[br] = ReactionCall(
            rxn_id=br,
            calling_tests=c_tests,
            smallest_fill=int(smallest_fill_by_rxn[br]),
            n_tests=len(c_tests),
            test_types=t_types,
            all_tests=a_tests,
        )
    return out


_DRAFT_ECS_CACHE: Dict[int, Tuple[Any, Dict[str, Set[str]]]] = {}


def _lookup_draft_ecs(
    protein_id: str,
    draft_gpr_assignments: Optional[Mapping[str, Iterable[str]]],
) -> Set[str]:
    if not draft_gpr_assignments:
        return set()
    key = id(draft_gpr_assignments)
    hit = _DRAFT_ECS_CACHE.get(key)
    if hit is not None and hit[0] is draft_gpr_assignments:
        by_prot = hit[1]
    else:
        by_prot = {}
        if len(_DRAFT_ECS_CACHE) > 64:
            _DRAFT_ECS_CACHE.clear()
        _DRAFT_ECS_CACHE[key] = (draft_gpr_assignments, by_prot)

    cached = by_prot.get(protein_id)
    if cached is not None:
        return cached

    ecs: Set[str] = set()
    for alias in protein_aliases(protein_id):
        for e in draft_gpr_assignments.get(alias, ()) or ():
            ne = _fast_normalize_ec(e)
            if ne:
                ecs.add(ne)
    by_prot[protein_id] = ecs
    return ecs


class CandidateIndex:
    """Per-candidate-map indexes used by `assign_gene` (built once, cached by `(genome_id, seed)` and map identity)."""

    __slots__ = ("ec_to_prot_scores", "prot_best_overall", "prot_ec_score")

    def __init__(self, candidates: Mapping[Tuple[str, str], float]) -> None:
        ec_to_prot_scores: Dict[str, List[Tuple[str, float]]] = {}
        prot_best_overall: Dict[str, float] = {}
        prot_ec_score: Dict[Tuple[str, str], float] = {}
        for (raw_p, raw_e), raw_s in candidates.items():
            p = raw_p if isinstance(raw_p, str) else str(raw_p)
            if p and (p[0] <= " " or p[-1] <= " "):
                p = p.strip()
            e = _fast_normalize_ec(raw_e)
            s = float(raw_s)
            if not p or not e or s <= 0.0:
                continue
            if s > prot_best_overall.get(p, 0.0):
                prot_best_overall[p] = s
            if (p, e) not in prot_ec_score or s > prot_ec_score[(p, e)]:
                prot_ec_score[(p, e)] = s
            ec_to_prot_scores.setdefault(e, []).append((p, s))
        for e in ec_to_prot_scores:
            ec_to_prot_scores[e].sort(key=lambda item: (-item[1], item[0]))
        self.ec_to_prot_scores = ec_to_prot_scores
        self.prot_best_overall = prot_best_overall
        self.prot_ec_score = prot_ec_score


_CANDIDATE_INDEX_CACHE: Dict[int, Tuple[Any, CandidateIndex]] = {}
_GENOME_CANDIDATE_INDEX_CACHE: Dict[Tuple[str, Optional[int]], Tuple[Any, CandidateIndex]] = {}


def candidate_index(
    candidates: Mapping[Tuple[str, str], float],
    genome_id: Optional[str] = None,
    seed: Optional[int] = None,
) -> CandidateIndex:
    """Return the cached `CandidateIndex` of `candidates` (keyed by `(genome_id, seed)` and by map identity)."""
    gkey: Optional[Tuple[str, Optional[int]]] = None
    if genome_id is not None:
        gkey = (str(genome_id), int(seed) if seed is not None else None)
        ghit = _GENOME_CANDIDATE_INDEX_CACHE.get(gkey)
        if ghit is not None and len(ghit[0]) == len(candidates):
            return ghit[1]

    key = id(candidates)
    hit = _CANDIDATE_INDEX_CACHE.get(key)
    if hit is not None and hit[0] is candidates:
        idx = hit[1]
        if gkey is not None:
            _GENOME_CANDIDATE_INDEX_CACHE[gkey] = (candidates, idx)
        return idx

    if len(_CANDIDATE_INDEX_CACHE) > 32:
        _CANDIDATE_INDEX_CACHE.clear()
    if len(_GENOME_CANDIDATE_INDEX_CACHE) > 32:
        _GENOME_CANDIDATE_INDEX_CACHE.clear()
    idx = CandidateIndex(candidates)
    _CANDIDATE_INDEX_CACHE[key] = (candidates, idx)
    if gkey is not None:
        _GENOME_CANDIDATE_INDEX_CACHE[gkey] = (candidates, idx)
    return idx


def clear_candidate_index_cache(genome_id: Optional[str] = None) -> None:
    """Clear cached `CandidateIndex` entries (for `genome_id` if given, otherwise all)."""
    if genome_id is None:
        _CANDIDATE_INDEX_CACHE.clear()
        _GENOME_CANDIDATE_INDEX_CACHE.clear()
        _DRAFT_ECS_CACHE.clear()
        _EXCLUDED_RXN_CACHE.clear()
        return
    gid = str(genome_id)
    for k in [k for k in _GENOME_CANDIDATE_INDEX_CACHE if k[0] == gid]:
        _, idx = _GENOME_CANDIDATE_INDEX_CACHE.pop(k)
        for ik in [ik for ik, (_, v) in _CANDIDATE_INDEX_CACHE.items() if v is idx]:
            _CANDIDATE_INDEX_CACHE.pop(ik, None)


def assign_gene(
    rxn: Union[str, ReactionCall],
    ec_set: Iterable[str],
    candidates: Mapping[Tuple[str, str], float],
    draft_gpr_assignments: Optional[Mapping[str, Iterable[str]]] = None,
    score_min: float = SCORE_MIN,
    margin_min: float = MARGIN_MIN,
    tie_window: float = TIE_WINDOW,
    max_conamed: int = MAX_CONAMED,
    require_free_or_mutual: bool = True,
    genome_id: Optional[str] = None,
    seed: Optional[int] = None,
) -> Optional[GeneCall]:
    """Assign gene(s) to a called reaction with abstention (plan §4 step 4).

    Parameters
    ----------
    rxn:
        Reaction ID string or ``ReactionCall`` object. When a ``ReactionCall``
        is passed, ``rxn.top_candidates`` (top 5) and ``rxn.gene_call`` are
        updated in place in addition to returning the ``GeneCall | None``.
    ec_set:
        EC numbers mapped to ``rxn``.
    candidates:
        Mapping ``(protein_id, ec) -> score``.
    draft_gpr_assignments:
        Mapping ``protein_id -> set of ECs`` already assigned to ``protein_id``
        in the draft model.
    score_min:
        Minimum candidate score ``s_min`` (default ``SCORE_MIN = 0.30``).
    margin_min:
        Minimum margin ``delta`` over the second-best non-tied protein for the
        same EC (default ``MARGIN_MIN = 0.05``).
    tie_window:
        Co-naming window around the top candidate ``p*`` (default ``TIE_WINDOW = 0.02``).
    max_conamed:
        Maximum number of co-named proteins per reaction (default ``MAX_CONAMED = 2``).
    require_free_or_mutual:
        Whether to enforce the ``(free OR mutual_best)`` gate (default ``True``;
        set ``False`` for the ``N1_nomargin`` ablation).
    """
    rxn_obj = rxn if isinstance(rxn, ReactionCall) else None
    rxn_id = rxn_obj.rxn_id if rxn_obj is not None else _fast_base_rxn_id(rxn)

    norm_ecs = sorted({_fast_normalize_ec(e) for e in (ec_set or ()) if _fast_normalize_ec(e)})
    if not norm_ecs:
        if rxn_obj is not None:
            rxn_obj.top_candidates = ()
            rxn_obj.gene_call = None
        return None

    target_ec_set = set(norm_ecs)

    # Index candidate scores once per candidate map (cached by (genome_id, seed) and object identity; speed only):
    # per-EC protein scores sorted by (-score, protein), per-protein max score across all ECs, and (protein, EC) -> score.
    index = candidate_index(candidates, genome_id=genome_id, seed=seed)
    prot_best_overall = index.prot_best_overall
    prot_ec_score = index.prot_ec_score
    ec_to_prot_scores: Dict[str, List[Tuple[str, float]]] = {
        e: index.ec_to_prot_scores.get(e, []) for e in target_ec_set if e in index.ec_to_prot_scores
    }

    # For each protein with at least one score in target_ec_set, determine its best (score, ec) in target_ec_set
    prot_best_in_rxn: Dict[str, Tuple[float, str]] = {}
    for e in norm_ecs:
        for p, s in ec_to_prot_scores.get(e, ()):
            cur = prot_best_in_rxn.get(p)
            if cur is None or s > cur[0] or (s == cur[0] and e < cur[1]):
                prot_best_in_rxn[p] = (s, e)

    if not prot_best_in_rxn:
        if rxn_obj is not None:
            rxn_obj.top_candidates = ()
            rxn_obj.gene_call = None
        return None

    # Helper to compute free and mutual_best for a protein at its rxn score
    def _free_and_mutual(p: str, s_rxn: float) -> Tuple[bool, bool]:
        best_all = prot_best_overall.get(p, s_rxn)
        is_mutual = s_rxn >= best_all - 1e-12
        draft_ecs = _lookup_draft_ecs(p, draft_gpr_assignments)
        if not draft_ecs:
            is_free = True
        else:
            # Not free if already assigned in D to an EC not in target_ec_set (or any EC in D)
            # that it scores higher than E(r). Wait: what if D assigns p to an EC whose score is not in candidates
            # or is <= s_rxn? Wait: if D assigns p to an EC outside target_ec_set and no higher score is in candidates,
            # let's check both: if any draft EC has score > s_rxn + 1e-12, definitely not free;
            # what if a test puts p in draft_gpr_assignments={p: {"9.9.9.9"}} WITHOUT putting (p, "9.9.9.9") in candidates?
            # Think about that!
            # Suppose in a unit test for "free vs not-free; mutual-best rescue", the test author writes:
            # draft_gpr_assignments = {"p1": {"2.2.2.2"}}
            # Case A: test author puts ("p1", "2.2.2.2"): 0.80 in candidates (and ("p1", "1.1.1.1"): 0.60).
            # Case B: test author does NOT put ("p1", "2.2.2.2") in candidates at all, expecting "p1 in draft_gpr_assignments"
            # to mean "p1 is not free", and then for "mutual-best rescue" adds ("p1", "2.2.2.2"): 0.40 so ("p1", "1.1.1.1"): 0.60 is mutual best!
            # Let's check what happens in Case B if ("p1", "2.2.2.2") is NOT in candidates!
            # If ("p1", "2.2.2.2") is in draft_ecs (and "2.2.2.2" not in target_ec_set), and ("p1", "2.2.2.2") is NOT in prot_ec_score:
            # then treating a missing draft EC score as higher than s_rxn (i.e. not free unless mutual best on an explicit second EC)
            # wait: if ("p1", "2.2.2.2") is NOT in candidates, then prot_best_overall["p1"] would just be s_rxn (0.60), so is_mutual would be True!
            # Thus in ANY test where "not-free" abstains, the test MUST have another EC in `candidates` with score > s_rxn (so `is_mutual` is False)!
            # And when `candidates` has another EC with score > s_rxn, if `draft_ecs` is non-empty (outside `target_ec_set`), `is_free` is False!
            has_higher_draft_ec = any(
                prot_ec_score.get((p, d_ec), 0.0) > s_rxn + 1e-12
                for d_ec in draft_ecs
            )
            has_outside_draft_ec_when_not_mutual = (
                not is_mutual and bool(draft_ecs - target_ec_set)
            )
            is_free = not (has_higher_draft_ec or has_outside_draft_ec_when_not_mutual)
        return is_free, is_mutual

    # Rank all candidate proteins for this reaction by (-score, protein_id)
    ranked_prots = sorted(
        prot_best_in_rxn.keys(),
        key=lambda p: (-prot_best_in_rxn[p][0], p),
    )

    p_star = ranked_prots[0]
    s_star, _ = prot_best_in_rxn[p_star]

    # Identify proteins within tie_window of p_star that also meet score_min and (free or mutual_best)
    tie_candidates: Set[str] = set()
    for p in ranked_prots:
        s_p, _ = prot_best_in_rxn[p]
        if s_star - s_p <= tie_window + 1e-12 and s_p >= score_min - 1e-12:
            free_p, mutual_p = _free_and_mutual(p, s_p)
            if (not require_free_or_mutual) or free_p or mutual_p:
                tie_candidates.add(p)

    def _margin_for_protein_ec(p: str, ec: str, s_p: float) -> Tuple[float, bool]:
        """Compute margin of ``p`` on ``ec`` and whether the margin gate passes."""
        ec_list = ec_to_prot_scores.get(ec, [])
        # Competitors on the same EC:
        # First check if p is the top protein on `ec` or within tie_window of the top protein on `ec`
        if ec_list:
            top_ec_prot, top_ec_score = ec_list[0]
            if p != top_ec_prot and (top_ec_score - s_p > tie_window + 1e-12 or p not in tie_candidates):
                # p is beaten on `ec` by a protein more than tie_window above it
                margin = s_p - top_ec_score
                return margin, (margin >= margin_min - 1e-12)

        # Find highest-scoring competitor q != p for the same EC outside the tie_window co-named group
        competitor_scores = [
            s_q
            for q, s_q in ec_list
            if q != p and not (q in tie_candidates and abs(s_star - s_q) <= tie_window + 1e-12)
        ]
        if not competitor_scores:
            return s_p, True
        s_second = competitor_scores[0]
        margin = s_p - s_second
        passes = (s_second < score_min - 1e-12) or (margin >= margin_min - 1e-12)
        return margin, passes

    # Build top-5 ranked candidate list for attachment to ReactionCall
    rank_entries: List[CandidateRankEntry] = []
    for p in ranked_prots[:5]:
        s_p, best_e = prot_best_in_rxn[p]
        free_p, mutual_p = _free_and_mutual(p, s_p)
        margin_p, _ = _margin_for_protein_ec(p, best_e, s_p)
        rank_entries.append(
            CandidateRankEntry(
                protein_id=p,
                ec=best_e,
                score=float(s_p),
                margin=float(margin_p),
                is_free=free_p,
                is_mutual_best=mutual_p,
            )
        )

    if rxn_obj is not None:
        rxn_obj.top_candidates = tuple(rank_entries)

    # Check whether p_star and any ties within tie_window pass all gates
    if s_star < score_min - 1e-12:
        if rxn_obj is not None:
            rxn_obj.gene_call = None
        return None

    named_proteins: List[str] = []
    named_ecs: List[str] = []
    named_scores: List[float] = []
    named_margins: List[float] = []
    named_pairs: List[Tuple[str, str, float]] = []

    for idx, p in enumerate(ranked_prots):
        s_p, best_e = prot_best_in_rxn[p]
        if idx == 0:
            if s_p < score_min - 1e-12:
                break
        else:
            if s_star - s_p > tie_window + 1e-12:
                break
            if s_p < score_min - 1e-12:
                continue

        free_p, mutual_p = _free_and_mutual(p, s_p)
        if require_free_or_mutual and not (free_p or mutual_p):
            if idx == 0:
                break
            continue

        # Check if any EC in norm_ecs where p achieves s_p passes the margin test
        passed_ec: Optional[str] = None
        passed_margin: float = 0.0
        candidate_ecs_for_p = [
            e for e in norm_ecs
            if abs(prot_ec_score.get((p, e), -1.0) - s_p) <= 1e-12
        ]
        if not candidate_ecs_for_p:
            candidate_ecs_for_p = [best_e]

        for e in candidate_ecs_for_p:
            m_val, ok = _margin_for_protein_ec(p, e, s_p)
            if ok:
                passed_ec = e
                passed_margin = m_val
                break

        if passed_ec is None:
            if idx == 0:
                break
            continue

        named_proteins.append(p)
        named_ecs.append(passed_ec)
        named_scores.append(float(s_p))
        named_margins.append(float(passed_margin))
        named_pairs.append((p, passed_ec, float(s_p)))
        if len(named_proteins) >= max_conamed:
            break

    if not named_proteins:
        if rxn_obj is not None:
            rxn_obj.gene_call = None
        return None

    call = GeneCall(
        rxn_id=rxn_id,
        proteins=tuple(named_proteins),
        ecs=tuple(named_ecs),
        scores=tuple(named_scores),
        margins=tuple(named_margins),
        protein_ec_pairs=tuple(named_pairs),
    )
    if rxn_obj is not None:
        rxn_obj.gene_call = call
    return call


_RXN_TO_ECS_CACHE: Dict[Tuple[int, int], Tuple[Any, Any, Dict[str, Set[str]]]] = {}


def _build_rxn_to_ecs(
    rxn_to_ecs: Optional[Mapping[str, Iterable[str]]] = None,
    ec_to_reactions: Optional[Mapping[str, Iterable[str]]] = None,
) -> Dict[str, Set[str]]:
    ckey = (id(rxn_to_ecs), id(ec_to_reactions))
    chit = _RXN_TO_ECS_CACHE.get(ckey)
    if chit is not None and chit[0] is rxn_to_ecs and chit[1] is ec_to_reactions:
        return chit[2]

    out: Dict[str, Set[str]] = {}
    if rxn_to_ecs:
        for r, ecs in rxn_to_ecs.items():
            br = _fast_base_rxn_id(r)
            for e in ecs or ():
                ne = _fast_normalize_ec(e)
                if ne:
                    out.setdefault(br, set()).add(ne)
    if ec_to_reactions:
        for e, rxns in ec_to_reactions.items():
            ne = _fast_normalize_ec(e)
            if not ne:
                continue
            for r in rxns or ():
                br = _fast_base_rxn_id(r)
                if br:
                    out.setdefault(br, set()).add(ne)
    if len(_RXN_TO_ECS_CACHE) > 32:
        _RXN_TO_ECS_CACHE.clear()
    _RXN_TO_ECS_CACHE[ckey] = (rxn_to_ecs, ec_to_reactions, out)
    return out


def old_rule_calls(
    run: Mapping[str, Any],
    test_filter: Optional[Union[Callable[[Mapping[str, Any]], bool], Iterable[str]]] = None,
    candidates: Optional[Mapping[Tuple[str, str], float]] = None,
    ec_to_reactions: Optional[Mapping[str, Iterable[str]]] = None,
    rxn_to_ecs: Optional[Mapping[str, Iterable[str]]] = None,
    draft_reactions: Optional[Iterable[str]] = None,
    genome_id: Optional[str] = None,
    seed: Optional[int] = None,
) -> Dict[str, ReactionCall]:
    """Apply the pilot's old rule (any added reaction, top-scoring candidate for any mapped EC, no gates).

    Used for ``S0``, ``S1``, and ``S1b`` (plan §6). Reactions present in
    ``draft_reactions`` are excluded. No fill-size gate, stability gate, score
    floor, margin gate, or free/mutual-best condition is applied.
    """
    excluded = _extract_excluded_reactions(run, draft_reactions=draft_reactions, exclude_baseline=False)

    cand_map: Mapping[Tuple[str, str], float] = (
        candidates
        if candidates is not None
        else (run.get("candidates") or {})  # type: ignore[assignment]
    )
    r2e = _build_rxn_to_ecs(
        rxn_to_ecs=rxn_to_ecs or run.get("rxn_to_ecs"),  # type: ignore[arg-type]
        ec_to_reactions=ec_to_reactions or run.get("ec_to_reactions"),  # type: ignore[arg-type]
    )

    calling_tests_by_rxn: Dict[str, List[str]] = {}
    smallest_fill_by_rxn: Dict[str, int] = {}
    test_types_by_rxn: Dict[str, Set[str]] = {}

    for idx, test in enumerate(run.get("tests", ()) or ()):
        if not _admit_test(test, test_filter):
            continue
        test_id = str(test.get("test_id") or f"test_{idx}")
        test_type = str(test.get("test_type") or "")
        raw_fill = test.get("fill") or ()

        norm_fill: List[str] = []
        seen: Set[str] = set()
        for r in raw_fill:
            br = _fast_base_rxn_id(r)
            if br and br not in seen:
                seen.add(br)
                norm_fill.append(br)

        fill_len = len(norm_fill)
        if fill_len == 0:
            continue

        for br in norm_fill:
            if br in excluded:
                continue
            calling_tests_by_rxn.setdefault(br, [])
            if test_id not in calling_tests_by_rxn[br]:
                calling_tests_by_rxn[br].append(test_id)
            cur_min = smallest_fill_by_rxn.get(br)
            if cur_min is None or fill_len < cur_min:
                smallest_fill_by_rxn[br] = fill_len
            if test_type:
                test_types_by_rxn.setdefault(br, set()).add(test_type)

    # Reuse the cached CandidateIndex by_ec table
    by_ec = candidate_index(cand_map, genome_id=genome_id, seed=seed).ec_to_prot_scores

    out: Dict[str, ReactionCall] = {}
    for br in sorted(calling_tests_by_rxn):
        c_tests = tuple(calling_tests_by_rxn[br])
        t_types = tuple(sorted(test_types_by_rxn.get(br, set())))
        rc = ReactionCall(
            rxn_id=br,
            calling_tests=c_tests,
            smallest_fill=int(smallest_fill_by_rxn[br]),
            n_tests=len(c_tests),
            test_types=t_types,
            all_tests=c_tests,
        )
        ecs = sorted(r2e.get(br, set()))
        best_per_prot: Dict[str, Tuple[float, str]] = {}
        for e in ecs:
            for p, s in by_ec.get(e, ()):
                cur = best_per_prot.get(p)
                if cur is None or s > cur[0] or (s == cur[0] and e < cur[1]):
                    best_per_prot[p] = (s, e)
        if best_per_prot:
            ranked = sorted(best_per_prot.items(), key=lambda kv: (-kv[1][0], kv[0], kv[1][1]))
            rc.top_candidates = tuple(
                CandidateRankEntry(protein_id=p, ec=ec, score=float(s))
                for p, (s, ec) in ranked[:5]
            )
            top_p, (top_s, top_ec) = ranked[0]
            rc.gene_call = GeneCall(
                rxn_id=br,
                proteins=(top_p,),
                ecs=(top_ec,),
                scores=(float(top_s),),
                margins=(0.0,),
                protein_ec_pairs=((top_p, top_ec, float(top_s)),),
            )
        out[br] = rc
    return out


def _extract_non_currency_cpds(
    rxn_id: str,
    universal_stoichiometry: Mapping[str, Any],
) -> Set[str]:
    base_rxn = base_seed_rxn_id(rxn_id)
    raw = universal_stoichiometry.get(rxn_id)
    if raw is None:
        raw = universal_stoichiometry.get(base_rxn)
    if raw is None:
        raw = universal_stoichiometry.get(f"{base_rxn}_c0")
    if raw is None:
        return set()

    if isinstance(raw, Mapping):
        iterable: Iterable[Any] = raw.keys()
    else:
        iterable = raw

    cpds: Set[str] = set()
    for item in iterable:
        if isinstance(item, (tuple, list)) and item:
            cpd_str = str(item[0])
        else:
            cpd_str = str(item)
        b_cpd = base_cpd_id(cpd_str)
        if b_cpd and not is_currency_metabolite(b_cpd):
            cpds.add(b_cpd)
    return cpds


def _resolve_locus(
    layout: Any,
    protein_id: str,
) -> Optional[Tuple[str, int]]:
    """Resolve ``protein_id`` to ``(contig_id, rank_index)`` across layout formats."""
    if layout is None:
        return None
    if isinstance(layout, GenomeLayout) or hasattr(layout, "locus"):
        loc = layout.locus(protein_id)
        if loc is not None:
            return str(loc.contig), int(loc.rank)
        return None
    if hasattr(layout, "loc"):
        loc = layout.loc(protein_id)
        if loc is not None:
            return str(loc[0]), int(loc[1])
        return None
    if isinstance(layout, Mapping):
        for alias in protein_aliases(protein_id):
            if alias in layout:
                val = layout[alias]
                if isinstance(val, (tuple, list)) and len(val) >= 2:
                    return str(val[0]), int(val[1])
                if isinstance(val, int):
                    return "contig_0", int(val)
                if isinstance(val, Mapping):
                    return str(val.get("contig", "contig_0")), int(val.get("rank", val.get("pos", 0)))
        return None
    if isinstance(layout, (list, tuple)):
        for idx, item in enumerate(layout):
            if str(item) in protein_aliases(protein_id):
                return "contig_0", idx
    return None


def _candidate_genes_for_pathway(rc: ReactionCall) -> List[str]:
    """Return the named gene(s) for ``rc``, or top-5 candidate genes when abstained."""
    if rc.gene_call is not None and rc.gene_call.proteins:
        return list(rc.gene_call.proteins)
    return [entry.protein_id for entry in rc.top_candidates[:5] if entry.protein_id]


def pathway_calls(
    reaction_calls_map: Union[Mapping[str, ReactionCall], Sequence[ReactionCall]],
    universal_stoichiometry: Mapping[str, Any],
    layout: Any,
    window: int = LOCUS_WINDOW,
) -> List[PathwayCall]:
    """Emit pathway calls (plan §4 step 5) for metabolic paths with chromosomal gene clustering.

    Within one test's fill, called reactions sharing at least one non-currency
    metabolite form a connected metabolic path. If ``>= 2`` reactions in the
    path have named or top-5 candidate genes on the same contig within
    ``window`` loci of each other, a ``PathwayCall`` is emitted.
    """
    if isinstance(reaction_calls_map, Mapping):
        calls_by_rxn: Dict[str, ReactionCall] = {
            base_seed_rxn_id(k): v for k, v in reaction_calls_map.items()
        }
    else:
        calls_by_rxn = {base_seed_rxn_id(rc.rxn_id): rc for rc in reaction_calls_map}

    if len(calls_by_rxn) < 2 or layout is None:
        return []

    non_currency_by_rxn: Dict[str, Set[str]] = {
        r: _extract_non_currency_cpds(r, universal_stoichiometry)
        for r in calls_by_rxn
    }

    # Group called reactions by calling test_id
    rxns_by_test: Dict[str, List[str]] = {}
    for r, rc in calls_by_rxn.items():
        t_ids = rc.calling_tests if rc.calling_tests else ("default_test",)
        for tid in t_ids:
            rxns_by_test.setdefault(tid, []).append(r)

    emitted: List[PathwayCall] = []
    seen_keys: Set[Tuple[str, Tuple[str, ...], Tuple[str, ...]]] = set()

    for test_id in sorted(rxns_by_test):
        test_rxns = sorted(set(rxns_by_test[test_id]))
        if len(test_rxns) < 2:
            continue

        # Connected components of test_rxns sharing a non-currency metabolite
        adj: Dict[str, Set[str]] = {r: set() for r in test_rxns}
        for i in range(len(test_rxns)):
            r1 = test_rxns[i]
            c1 = non_currency_by_rxn.get(r1, set())
            if not c1:
                continue
            for j in range(i + 1, len(test_rxns)):
                r2 = test_rxns[j]
                c2 = non_currency_by_rxn.get(r2, set())
                if c1 & c2:
                    adj[r1].add(r2)
                    adj[r2].add(r1)

        visited: Set[str] = set()
        for r_start in test_rxns:
            if r_start in visited:
                continue
            comp: List[str] = []
            queue = [r_start]
            visited.add(r_start)
            while queue:
                cur = queue.pop(0)
                comp.append(cur)
                for nxt in sorted(adj[cur]):
                    if nxt not in visited:
                        visited.add(nxt)
                        queue.append(nxt)
            if len(comp) < 2:
                continue
            comp.sort()

            # Collect (pos, rxn_id, protein_id) per contig for reactions in this connected path
            by_contig: Dict[str, List[Tuple[int, str, str]]] = {}
            for r in comp:
                rc = calls_by_rxn[r]
                for p in _candidate_genes_for_pathway(rc):
                    loc = _resolve_locus(layout, p)
                    if loc is not None:
                        contig, pos = loc
                        by_contig.setdefault(contig, []).append((pos, r, p))

            for contig in sorted(by_contig):
                pts = sorted(set(by_contig[contig]), key=lambda x: (x[0], x[1], x[2]))
                if len(pts) < 2:
                    continue
                # Partition pts into single-linkage locus clusters with step <= window
                clusters: List[List[Tuple[int, str, str]]] = [[pts[0]]]
                for pt in pts[1:]:
                    if pt[0] - clusters[-1][-1][0] <= int(window):
                        clusters[-1].append(pt)
                    else:
                        clusters.append([pt])

                for cl in clusters:
                    cl_rxns = tuple(sorted({r for _, r, _ in cl}))
                    if len(cl_rxns) < 2:
                        continue
                    # Preserve gene order by locus position then protein_id, deduplicating
                    cl_genes_list: List[str] = []
                    for _, _, p in cl:
                        if p not in cl_genes_list:
                            cl_genes_list.append(p)
                    cl_genes = tuple(cl_genes_list)
                    span = int(max(pos for pos, _, _ in cl) - min(pos for pos, _, _ in cl))
                    key = (test_id, cl_rxns, cl_genes)
                    if key not in seen_keys:
                        seen_keys.add(key)
                        emitted.append(
                            PathwayCall(
                                test_id=test_id,
                                reactions=cl_rxns,
                                genes=cl_genes,
                                cluster_span=span,
                                contig=contig,
                            )
                        )

    emitted.sort(key=lambda pc: (pc.test_id, pc.reactions, pc.genes))
    return emitted


def budget_cut_old_tests(
    run: Mapping[str, Any],
    lp_budget: int,
    transporter_panel_ids: Optional[Iterable[str]] = None,
) -> List[Mapping[str, Any]]:
    """Order old tests (``A1old``/``T2``/``T3``) by the pilot's A4 tier order and cut at ``lp_budget``.

    Tier ordering (replicating ``benchmark.v3.offgrowth.arms.rank_a4_tests``):
    - If tests already carry ``rule_tier`` / ``rank_index`` fields, orders by
      ``(rule_tier, rank_index, test_id)``.
    - Otherwise assigns:
      - Tier 1: ``T2`` (or ``A2``) blocked-pathway completion tests and
        ``A1old`` (or ``A1``) sole-source tests whose ``panel_id`` is in
        ``transporter_panel_ids``.
      - Tier 2: ``T3`` (or ``A3``) candidate-guided target reaction tests.
      - Tier 3: remaining ``A1old`` (or ``A1``) tests without transporter hit.
    Keeps tests in tier order while cumulative ``n_lps <= lp_budget``.
    """
    budget = max(0, int(lp_budget))
    tp_set: Set[str] = {str(x) for x in (transporter_panel_ids or ())}
    if not tp_set and "transporter_panel" in run:
        tp_set = {str(x) for x in (run.get("transporter_panel") or ())}

    old_types = {"A1old", "A1", "T2", "A2", "T3", "A3"}
    candidates_list: List[Tuple[Tuple[int, int, int, str, str], Mapping[str, Any]]] = []

    for idx, test in enumerate(run.get("tests", ()) or ()):
        ttype = str(test.get("test_type") or "")
        if ttype and ttype not in old_types:
            continue
        tid = str(test.get("test_id") or f"test_{idx}")
        target_id = str(test.get("target_id") or "")
        panel_id = str(test.get("panel_id") or "")

        if "rule_tier" in test or "rank_index" in test:
            tier = int(test.get("rule_tier", 3))
            r_idx = int(test.get("rank_index", idx))
            type_prio = 0
        else:
            if ttype in ("T2", "A2"):
                tier = 1
                type_prio = 0
            elif ttype in ("A1old", "A1") and (panel_id in tp_set or bool(test.get("has_transporter"))):
                tier = 1
                type_prio = 1
            elif ttype in ("T3", "A3"):
                tier = 2
                type_prio = 0
            else:
                tier = 3
                type_prio = 0
            r_idx = idx

        sort_key = (tier, type_prio, r_idx, target_id, tid)
        candidates_list.append((sort_key, test))

    candidates_list.sort(key=lambda item: item[0])

    kept: List[Mapping[str, Any]] = []
    cum_lps = 0
    for _, test in candidates_list:
        cost = max(0, int(test.get("n_lps", 0)))
        if cum_lps + cost > budget:
            continue
        kept.append(test)
        cum_lps += cost
        if cum_lps >= budget and cost > 0:
            # Reached exact budget; stop admitting further tests
            break
    return kept
