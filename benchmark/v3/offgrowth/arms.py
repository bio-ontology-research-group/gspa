#!/usr/bin/env python3
"""Weighted gap-filling LP (§4) and deterministic arms A0-A4 (§5, §5.2).

No language model or stochastic search is used anywhere in this module.
Every LP solve and every gap-fill solve is counted with CPU time, wall time,
and peak RSS memory. A deterministic 5% subset of gap-fill LPs is re-solved
with a secondary LP solver (GLPK or SciPy/HiGHS) to verify solver agreement.
"""

from __future__ import annotations

import copy
import csv
from dataclasses import asdict, dataclass, field
import hashlib
import heapq
import math
import multiprocessing as mp
import os
from pathlib import Path
import resource
import time
from typing import Dict, Iterable, List, Mapping, NamedTuple, Optional, Sequence, Set, Tuple

import cobra
import numpy as np
import optlang
import scipy
import scipy.optimize
import scipy.sparse as sp

try:
    from swiglpk import (
        doubleArray,
        glp_get_mat_row,
        glp_get_num_cols,
        glp_get_obj_coef,
        intArray,
    )

    _HAS_SWIGLPK = True
except ImportError:  # pragma: no cover
    _HAS_SWIGLPK = False

from .mask import (
    CandidateClaim,
    MaskingAuditReport,
    TransporterHit,
    assert_masking_clean,
    base_seed_rxn_id,
)
from .panel import (
    PanelEntry,
    SeedMetabolite,
    SeedReaction,
    is_complete_ec,
    normalize_ec,
)


WEIGHT_CANDIDATE_SUPPORTED: float = 1.0
WEIGHT_NO_CANDIDATE: float = 100.0
DEFAULT_FLUX_TOL: float = 1e-6
DEFAULT_EPSILON_FRAC: float = 0.01
DEFAULT_PANEL_UPTAKE: float = 10.0


@dataclass(frozen=True)
class SolverInfo:
    primary_solver: str
    replay_solver: str
    highspy_available: bool
    available_solvers: Tuple[str, ...]


@dataclass(frozen=True)
class SolverReplayRecord:
    test_id: str
    primary_solver: str
    replay_solver: str
    primary_objective: float
    replay_objective: float
    relative_obj_diff: float
    primary_added_rxns: Tuple[str, ...]
    replay_added_rxns: Tuple[str, ...]
    agreed: bool
    scipy_version: str = field(
        default_factory=lambda: str(getattr(scipy, "__version__", ""))
    )


@dataclass(frozen=True)
class GapFillTestResult:
    test_id: str
    arm: str  # "A0", "A1", "A2", "A3"
    target_id: str
    feasible_in_universal: bool
    already_satisfied_in_draft: bool
    gapfill_ran: bool
    vmax_universal: float
    objective_value: float
    added_reactions: Tuple[str, ...]
    n_lps: int
    n_gapfills: int
    cpu_seconds: float
    wall_seconds: float
    replay: Optional[SolverReplayRecord] = None


@dataclass(frozen=True)
class RankedTestEntry:
    test_id: str
    arm: str
    target_id: str
    rule_tier: int  # 1, 2, 3, or 4
    supporting_genes: Tuple[str, ...]
    n_supporting_genes: int
    rank_index: int  # 0-based deterministic rank


@dataclass
class ArmRunResult:
    arm: str  # "A0", "A1", "A2", "A3", "A1_A2_A3", "A4@10%", "A4@25%"
    tests_considered: int
    tests_executed: int
    n_lps: int
    n_gapfills: int
    cpu_seconds: float
    wall_seconds: float
    peak_rss_kb: int
    added_reactions: Tuple[str, ...]
    test_results: Tuple[GapFillTestResult, ...]
    ranked_entries: Tuple[RankedTestEntry, ...] = ()
    replays: Tuple[SolverReplayRecord, ...] = ()


@dataclass(frozen=True)
class KeggModuleDef:
    module_id: str
    name: str
    ec_numbers: Tuple[str, ...]
    reaction_ids: Tuple[str, ...]


def detect_solvers() -> SolverInfo:
    """Detect LP solvers that COBRApy itself accepts.

    optlang can report HiGHS as installed while COBRApy has no interface for it (optlang reaches HiGHS only through
    the hybrid interface, which also needs osqp). Only names in cobra's own solver table are usable. Preference:
    HiGHS (``hybrid``) primary with GLPK replay; otherwise GLPK primary with the SciPy interface (scipy.optimize.linprog,
    HiGHS backend) as the independent replay solver; the replay solver always differs from the primary when possible.
    """
    from cobra.util import solver as cobra_solver

    usable = tuple(sorted(cobra_solver.solvers))
    avail = tuple(sorted(k.lower() for k, v in optlang.available_solvers.items() if v))
    highs_ok = "hybrid" in usable
    if highs_ok:
        primary = "hybrid"
        replay = "glpk" if "glpk" in usable else ("scipy" if "scipy" in usable else "hybrid")
    elif "glpk" in usable:
        primary = "glpk"
        replay = "scipy" if "scipy" in usable else ("glpk_exact" if "glpk_exact" in usable else "glpk")
    elif "scipy" in usable:
        primary = "scipy"
        replay = "scipy"
    else:
        raise RuntimeError(f"No usable LP solver in cobra.util.solver.solvers: {usable} (optlang reports {avail})")

    return SolverInfo(
        primary_solver=primary,
        replay_solver=replay,
        highspy_available=highs_ok,
        available_solvers=usable,
    )


def peak_rss_kb() -> int:
    """Return peak resident set size in KiB from `resource.getrusage`."""
    return int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)


def load_gapsmith_medium_csv(path: Path) -> Dict[str, float]:
    """Load `<mag>-medium.csv` into `{ 'EX_cpdXXXXX_e0': max_flux }`."""
    medium: Dict[str, float] = {}
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line_idx, raw_line in enumerate(fh):
            line = raw_line.strip()
            if not line:
                continue
            if line_idx == 0 and line.lower().startswith("compounds,"):
                continue
            parts = line.split(",")
            if len(parts) < 3:
                continue
            cpd_id = parts[0].strip()
            flux_str = parts[-1].strip()
            if not cpd_id.startswith("cpd"):
                continue
            try:
                max_flux = float(flux_str)
            except ValueError:
                continue
            medium[f"EX_{cpd_id}_e0"] = abs(max_flux)
    return medium


def load_kegg_modules(path: Optional[Path]) -> Tuple[KeggModuleDef, ...]:
    """Load KEGG module definitions from `kegg_modules.tsv` or `kegg_module_defs.json`."""
    if path is None or not path.exists():
        return ()
    if path.suffix.lower() == ".json":
        import json

        data = json.loads(path.read_text(encoding="utf-8"))
        defs: List[KeggModuleDef] = []
        if isinstance(data, dict):
            for mod_id, val in sorted(data.items()):
                ecs_raw = val.get("ecs") or val.get("ec_numbers") or []
                rxns_raw = val.get("reactions") or val.get("reaction_ids") or []
                ecs = tuple(sorted({normalize_ec(e) for e in ecs_raw if is_complete_ec(e)}))
                rxns = tuple(sorted({base_seed_rxn_id(r) for r in rxns_raw if r}))
                defs.append(
                    KeggModuleDef(
                        module_id=str(mod_id),
                        name=str(val.get("name") or mod_id),
                        ec_numbers=ecs,
                        reaction_ids=rxns,
                    )
                )
        return tuple(defs)

    by_mod: Dict[str, Tuple[str, Set[str], Set[str]]] = {}
    with path.open("r", encoding="utf-8", errors="replace", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            mod_id = (row.get("pathway_id") or row.get("module_id") or "").strip()
            if not mod_id:
                continue
            name = (row.get("pathway_name") or row.get("name") or mod_id).strip()
            ec = normalize_ec(row.get("ec_number") or row.get("ec") or "")
            rxn = (row.get("reaction_id") or row.get("rxn") or "").strip()
            if mod_id not in by_mod:
                by_mod[mod_id] = (name, set(), set())
            if is_complete_ec(ec):
                by_mod[mod_id][1].add(ec)
            if rxn:
                by_mod[mod_id][2].add(base_seed_rxn_id(rxn))

    return tuple(
        KeggModuleDef(
            module_id=mod_id,
            name=name,
            ec_numbers=tuple(sorted(ecs)),
            reaction_ids=tuple(sorted(rxns)),
        )
        for mod_id, (name, ecs, rxns) in sorted(by_mod.items())
    )


def map_reactions_to_candidate_genes(
    candidates: Sequence[CandidateClaim],
    ec_to_reactions: Mapping[str, Sequence[str]],
) -> Dict[str, Tuple[Tuple[str, float], ...]]:
    """Map each base SEED reaction ID to sorted `(protein_id, max_score)` pairs."""
    rxn_gene_scores: Dict[str, Dict[str, float]] = {}
    for claim in candidates:
        ec = normalize_ec(claim.ec)
        for rxn_id in ec_to_reactions.get(ec, ()):
            base_id = base_seed_rxn_id(rxn_id)
            by_gene = rxn_gene_scores.setdefault(base_id, {})
            if claim.protein_id not in by_gene or claim.score > by_gene[claim.protein_id]:
                by_gene[claim.protein_id] = claim.score

    result: Dict[str, Tuple[Tuple[str, float], ...]] = {}
    for base_id, by_gene in rxn_gene_scores.items():
        # Sort by descending score, then ascending protein_id for deterministic tie-breaking
        ordered = sorted(by_gene.items(), key=lambda kv: (-kv[1], kv[0]))
        result[base_id] = tuple(ordered)
    return result


def compute_reaction_weights(
    universal_reactions: Mapping[str, SeedReaction],
    draft_base_rxn_ids: Set[str],
    rxn_to_candidates: Mapping[str, Sequence[Tuple[str, float]]],
) -> Dict[str, float]:
    """Compute gap-filling penalty weight $w_r \\in \\{1.0, 100.0\\}$ for $r \\in U \\setminus D$ (§4)."""
    weights: Dict[str, float] = {}
    for rxn_id, rxn in universal_reactions.items():
        if not rxn.is_usable:
            continue
        base_id = base_seed_rxn_id(rxn_id)
        if base_id in draft_base_rxn_ids:
            continue
        has_candidate = bool(rxn_to_candidates.get(base_id))
        weights[base_id] = WEIGHT_CANDIDATE_SUPPORTED if has_candidate else WEIGHT_NO_CANDIDATE
    return weights


def _should_replay_test(
    test_id: str,
    replay_fraction: float,
    force_replay: bool = False,
) -> bool:
    if force_replay:
        return True
    if replay_fraction <= 0.0:
        return False
    if replay_fraction >= 1.0:
        return True
    digest = hashlib.sha256(test_id.encode("utf-8")).digest()
    bucket = int.from_bytes(digest[:4], "big") / float(2**32)
    return bucket < replay_fraction


def _find_biomass_rxn_id(model: cobra.Model) -> str:
    for candidate_id in ("bio1", "R_bio1", "biomass"):
        if candidate_id in model.reactions:
            return candidate_id
    for rxn in model.reactions:
        if "bio" in rxn.id.lower() or "biomass" in rxn.name.lower():
            return rxn.id
    raise ValueError("Could not locate biomass reaction ('bio1') in model.")


def _ensure_exchange_reactions(
    model: cobra.Model,
    exchange_ids: Iterable[str],
) -> None:
    """Ensure `EX_<cpd>_e0` exists in `model` for every requested exchange."""
    new_mets: Dict[str, cobra.Metabolite] = {}
    new_rxns: List[cobra.Reaction] = []
    for ex_id in sorted(set(exchange_ids)):
        if ex_id in model.reactions:
            continue
        if not ex_id.startswith("EX_"):
            continue
        species_id = ex_id[len("EX_") :]
        if species_id in model.metabolites:
            met = model.metabolites.get_by_id(species_id)
        elif species_id in new_mets:
            met = new_mets[species_id]
        else:
            comp = species_id.rsplit("_", 1)[-1] if "_" in species_id else "e0"
            met = cobra.Metabolite(id=species_id, name=species_id, compartment=comp)
            new_mets[species_id] = met
        ex_rxn = cobra.Reaction(
            id=ex_id,
            name=f"Exchange {species_id}",
            lower_bound=0.0,
            upper_bound=1000.0,
        )
        ex_rxn.add_metabolites({met: -1.0})
        new_rxns.append(ex_rxn)
    if new_mets:
        model.add_metabolites(list(new_mets.values()))
    if new_rxns:
        model.add_reactions(new_rxns)


def apply_medium(
    model: cobra.Model,
    medium: Mapping[str, float],
) -> None:
    """Apply uptake bounds `-abs(flux)` for exchanges in `medium` and `0.0` for all other `EX_*`."""
    _ensure_exchange_reactions(model, medium.keys())
    for rxn in model.reactions:
        if rxn.id.startswith("EX_"):
            uptake = medium.get(rxn.id, 0.0)
            new_lb = -abs(float(uptake))
            if rxn.lower_bound != new_lb:
                rxn.lower_bound = new_lb
            if rxn.upper_bound < 1000.0:
                rxn.upper_bound = 1000.0


def build_sole_source_medium(
    default_medium: Mapping[str, float],
    panel_entry: PanelEntry,
    metabolites: Mapping[str, SeedMetabolite],
    uptake_rate: float = DEFAULT_PANEL_UPTAKE,
) -> Dict[str, float]:
    """Replace all sources of `panel_entry.element` in `default_medium` with `panel_entry.exchange_rxn`."""
    role = panel_entry.element.upper()
    new_medium: Dict[str, float] = {}
    for ex_id, flux in default_medium.items():
        if flux <= 0.0:
            continue
        cpd_id = ""
        if ex_id.startswith("EX_") and ex_id.endswith("_e0"):
            cpd_id = ex_id[len("EX_") : -len("_e0")]
        met = metabolites.get(cpd_id)
        if met is not None and role in met.elements:
            # Strip existing source of element X
            continue
        new_medium[ex_id] = flux
    new_medium[panel_entry.exchange_rxn] = abs(uptake_rate)
    return new_medium


def build_union_medium(
    default_medium: Mapping[str, float],
    panel_entries: Sequence[PanelEntry],
    uptake_rate: float = DEFAULT_PANEL_UPTAKE,
) -> Dict[str, float]:
    """Union of default medium and all panel exchange reactions (§5 A2)."""
    union_med = dict(default_medium)
    for entry in panel_entries:
        union_med[entry.exchange_rxn] = max(
            union_med.get(entry.exchange_rxn, 0.0), abs(uptake_rate)
        )
    return union_med


class LPConstraintMatrices(NamedTuple):
    """Sparse CSR representation matching `cobra.util.array.constraint_matrices` (`cobra/util/array.py:203-246`)."""

    equalities: sp.csr_matrix
    b: np.ndarray
    inequalities: sp.csr_matrix
    bounds: np.ndarray
    variable_fixed: np.ndarray
    variable_bounds: np.ndarray


def extract_lp_constraint_matrices(
    model: cobra.Model,
    zero_tol: float = DEFAULT_FLUX_TOL,
) -> LPConstraintMatrices:
    """Extract LP constraint matrices in CSR format with the same fields as `cobra.util.array.constraint_matrices`.

    Avoids two defects in COBRApy 0.32's `cobra.util.array.constraint_matrices`:
      1. `array_type="dok"` raises `ValueError: shape must have length in (2,). Got new_shape=(0,)` whenever
         `inequality_rows` is empty (`n_ineq == 0`).
      2. `cobra/util/array.py:226-227` calls `const.get_linear_coefficients(model.variables)` and builds a dense
         `n_constraints x n_variables` Python list-of-lists before converting to sparse (`O(n_constraints * n_variables)`
         time and memory).
    Here non-zero entries are read in `O(nnz)` time directly from the active solver rows (or reaction stoichiometry
    plus non-metabolite constraint expressions on non-GLPK backends) and assembled once into `scipy.sparse.csr_matrix`.
    """
    model.solver.update()
    variables = list(model.variables)
    n_vars = len(variables)

    var_bounds = np.empty((n_vars, 2), dtype=float)
    for j, v in enumerate(variables):
        var_bounds[j, 0] = -np.inf if v.lb is None else float(v.lb)
        var_bounds[j, 1] = np.inf if v.ub is None else float(v.ub)
    variable_fixed = (var_bounds[:, 1] - var_bounds[:, 0]) < zero_tol

    eq_rows: List[int] = []
    eq_cols: List[int] = []
    eq_vals: List[float] = []
    b_list: List[float] = []

    ineq_rows: List[int] = []
    ineq_cols: List[int] = []
    ineq_vals: List[float] = []
    ineq_bounds: List[List[float]] = []

    glp_prob = getattr(model.solver, "problem", None)
    use_glpk_direct = (
        _HAS_SWIGLPK
        and glp_prob is not None
        and (n_vars == 0 or hasattr(variables[0], "_index"))
    )

    if use_glpk_direct:
        num_cols = int(glp_get_num_cols(glp_prob))
        col_to_var_idx = np.full(num_cols + 1, -1, dtype=np.int64)
        for j, v in enumerate(variables):
            c_idx = int(v._index)
            if 0 <= c_idx <= num_cols:
                col_to_var_idx[c_idx] = j
        ia = intArray(num_cols + 1)
        da = doubleArray(num_cols + 1)
        for const in model.constraints:
            lb = -np.inf if const.lb is None else float(const.lb)
            ub = np.inf if const.ub is None else float(const.ub)
            equality = (ub - lb) < zero_tol
            nnz = int(glp_get_mat_row(glp_prob, int(const._index), ia, da))
            if equality:
                r_idx = len(b_list)
                b_list.append(lb if abs(lb) > zero_tol else 0.0)
                for k in range(1, nnz + 1):
                    v_idx = int(col_to_var_idx[ia[k]])
                    if v_idx >= 0:
                        eq_rows.append(r_idx)
                        eq_cols.append(v_idx)
                        eq_vals.append(float(da[k]))
            else:
                r_idx = len(ineq_bounds)
                ineq_bounds.append([lb, ub])
                for k in range(1, nnz + 1):
                    v_idx = int(col_to_var_idx[ia[k]])
                    if v_idx >= 0:
                        ineq_rows.append(r_idx)
                        ineq_cols.append(v_idx)
                        ineq_vals.append(float(da[k]))
    else:
        var_idx_by_name: Dict[str, int] = {v.name: j for j, v in enumerate(variables)}
        met_row_terms: Dict[str, List[Tuple[int, float]]] = {
            m.id: [] for m in model.metabolites
        }
        for rxn in model.reactions:
            f_idx = var_idx_by_name[rxn.forward_variable.name]
            r_idx = var_idx_by_name[rxn.reverse_variable.name]
            for met, coeff in rxn._metabolites.items():
                c_val = float(coeff)
                if c_val != 0.0:
                    terms = met_row_terms.setdefault(met.id, [])
                    terms.append((f_idx, c_val))
                    terms.append((r_idx, -c_val))

        for const in model.constraints:
            lb = -np.inf if const.lb is None else float(const.lb)
            ub = np.inf if const.ub is None else float(const.ub)
            equality = (ub - lb) < zero_tol
            if const.name in met_row_terms:
                row_entries = met_row_terms[const.name]
            else:
                coefs = const.get_linear_coefficients(variables)
                row_entries = [
                    (var_idx_by_name[v.name], float(val))
                    for v, val in coefs.items()
                    if val != 0.0
                ]
            if equality:
                r_idx = len(b_list)
                b_list.append(lb if abs(lb) > zero_tol else 0.0)
                for v_idx, val in row_entries:
                    eq_rows.append(r_idx)
                    eq_cols.append(v_idx)
                    eq_vals.append(val)
            else:
                r_idx = len(ineq_bounds)
                ineq_bounds.append([lb, ub])
                for v_idx, val in row_entries:
                    ineq_rows.append(r_idx)
                    ineq_cols.append(v_idx)
                    ineq_vals.append(val)

    n_eq = len(b_list)
    n_ineq = len(ineq_bounds)
    equalities = sp.csr_matrix((eq_vals, (eq_rows, eq_cols)), shape=(n_eq, n_vars), dtype=float)
    inequalities = sp.csr_matrix(
        (ineq_vals, (ineq_rows, ineq_cols)), shape=(n_ineq, n_vars), dtype=float
    )
    b_arr = np.asarray(b_list, dtype=float)
    bounds_arr = (
        np.asarray(ineq_bounds, dtype=float)
        if n_ineq > 0
        else np.empty((0, 2), dtype=float)
    )
    return LPConstraintMatrices(
        equalities=equalities,
        b=b_arr,
        inequalities=inequalities,
        bounds=bounds_arr,
        variable_fixed=variable_fixed,
        variable_bounds=var_bounds,
    )


def _extract_lp_objective(
    model: cobra.Model,
    variables: Sequence[object],
) -> Tuple[np.ndarray, str]:
    """Extract linear objective coefficient vector `c` and direction (`'max'` or `'min'`) from `model`."""
    model.solver.update()
    direction = str(model.objective_direction)
    n_vars = len(variables)
    glp_prob = getattr(model.solver, "problem", None)
    if (
        _HAS_SWIGLPK
        and glp_prob is not None
        and (n_vars == 0 or hasattr(variables[0], "_index"))
    ):
        c = np.array(
            [float(glp_get_obj_coef(glp_prob, int(v._index))) for v in variables],
            dtype=float,
        )
    else:
        obj_coefs = model.solver.objective.get_linear_coefficients(variables)
        c = np.array([float(obj_coefs.get(v, 0.0)) for v in variables], dtype=float)
    return c, direction


def _solve_model_lp_with_highs(
    model: cobra.Model,
    zero_tol: float = DEFAULT_FLUX_TOL,
) -> Tuple[bool, float, np.ndarray, Dict[str, int]]:
    """Extract the current LP from `model` in CSR form and solve with `scipy.optimize.linprog(method='highs')`."""
    variables = list(model.variables)
    var_idx_by_name: Dict[str, int] = {v.name: j for j, v in enumerate(variables)}
    problem = extract_lp_constraint_matrices(model, zero_tol=zero_tol)
    c, direction = _extract_lp_objective(model, variables)

    c_min = -c if direction == "max" else c

    A_eq = sp.csr_matrix(problem.equalities, dtype=float) if problem.equalities.shape[0] > 0 else None
    b_eq = np.asarray(problem.b, dtype=float) if problem.equalities.shape[0] > 0 else None

    A_ub: Optional[sp.csr_matrix] = None
    b_ub: Optional[np.ndarray] = None
    if problem.inequalities.shape[0] > 0:
        ineq_csr = sp.csr_matrix(problem.inequalities, dtype=float)
        ineq_bnds = (
            problem.bounds.toarray()
            if sp.issparse(problem.bounds)
            else np.asarray(problem.bounds, dtype=float)
        )
        ub_mask = np.isfinite(ineq_bnds[:, 1])
        lb_mask = np.isfinite(ineq_bnds[:, 0])
        ub_mats: List[sp.csr_matrix] = []
        ub_vecs: List[np.ndarray] = []
        if np.any(ub_mask):
            ub_mats.append(ineq_csr[ub_mask])
            ub_vecs.append(ineq_bnds[ub_mask, 1])
        if np.any(lb_mask):
            ub_mats.append(-ineq_csr[lb_mask])
            ub_vecs.append(-ineq_bnds[lb_mask, 0])
        if ub_mats:
            A_ub = sp.vstack(ub_mats, format="csr")
            b_ub = np.concatenate(ub_vecs)

    var_bounds = (
        problem.variable_bounds.toarray().copy()
        if sp.issparse(problem.variable_bounds)
        else np.array(problem.variable_bounds, dtype=float, copy=True)
    )
    if np.any(problem.variable_fixed):
        var_bounds[problem.variable_fixed, 1] = var_bounds[problem.variable_fixed, 0]

    res = scipy.optimize.linprog(
        c=c_min,
        A_ub=A_ub,
        b_ub=b_ub,
        A_eq=A_eq,
        b_eq=b_eq,
        bounds=var_bounds,
        method="highs",
    )
    if not res.success or res.x is None or res.fun is None or math.isnan(float(res.fun)):
        return False, 0.0, np.zeros(len(variables), dtype=float), var_idx_by_name
    opt_obj = -float(res.fun) if direction == "max" else float(res.fun)
    return True, opt_obj, np.asarray(res.x, dtype=float), var_idx_by_name


class WeightedGapFiller:
    """Universal-completed model and LP solver wrapper implementing §4."""

    def __init__(
        self,
        draft_model: cobra.Model,
        universal_reactions: Mapping[str, SeedReaction],
        rxn_weights: Mapping[str, float],
        solver_info: Optional[SolverInfo] = None,
        flux_tol: float = DEFAULT_FLUX_TOL,
        epsilon_frac: float = DEFAULT_EPSILON_FRAC,
        replay_weight_overrides: Optional[Mapping[str, float]] = None,
    ) -> None:
        self.solver_info = solver_info or detect_solvers()
        self.flux_tol = flux_tol
        self.epsilon_frac = epsilon_frac
        self.universal_reactions = universal_reactions
        self.rxn_weights = dict(rxn_weights)
        self.replay_weight_overrides: Dict[str, float] = (
            dict(replay_weight_overrides) if replay_weight_overrides is not None else {}
        )

        self.draft_model = copy.deepcopy(draft_model)
        self.draft_model.solver = self.solver_info.primary_solver
        self.draft_model.objective = self.draft_model.problem.Objective(0.0, direction="max")
        self._draft_active_vars: List[object] = []

        self.draft_base_ids: Set[str] = {
            base_seed_rxn_id(r.id)
            for r in self.draft_model.reactions
            if not r.id.startswith("EX_")
        }

        # Build universal-completed model with split non-negative candidate reactions
        self.universal_model, self.candidate_split_map = self._build_universal_completed_model(
            self.draft_model, self.solver_info.primary_solver
        )
        # Separate copy for Step 1 (v_t^max) so Step 2's 13k-term penalty objective is set once
        self.vmax_model = copy.deepcopy(self.universal_model)
        self.vmax_model.objective = self.vmax_model.problem.Objective(0.0, direction="max")
        self._vmax_active_vars: List[object] = []

        self._cand_var_entries = self._configure_penalty_objective(self.universal_model)

    def _configure_penalty_objective(
        self,
        u_model: cobra.Model,
    ) -> List[Tuple[str, float, Tuple[object, ...]]]:
        """Set the fixed penalty objective -sum_{r in U\\D} w_r v_r once on `u_model`."""
        u_model.objective = u_model.problem.Objective(0.0, direction="max")
        coeffs: Dict[object, float] = {}
        entries: List[Tuple[str, float, Tuple[object, ...]]] = []
        for base_id, split_ids in self.candidate_split_map.items():
            w_r = self.rxn_weights.get(base_id, WEIGHT_NO_CANDIDATE)
            fvars: List[object] = []
            for s_id in split_ids:
                rxn_obj = u_model.reactions.get_by_id(s_id)
                coeffs[rxn_obj.forward_variable] = -w_r
                fvars.append(rxn_obj.forward_variable)
            entries.append((base_id, w_r, tuple(fvars)))
        if coeffs:
            u_model.solver.objective.set_linear_coefficients(coeffs)
        return entries

    def _build_universal_completed_model(
        self,
        base_model: cobra.Model,
        solver_name: str,
    ) -> Tuple[cobra.Model, Dict[str, Tuple[str, ...]]]:
        u_model = copy.deepcopy(base_model)
        u_model.solver = solver_name

        existing_mets: Dict[str, cobra.Metabolite] = {m.id: m for m in u_model.metabolites}
        new_mets: List[cobra.Metabolite] = []
        for rxn_id in sorted(self.universal_reactions):
            urxn = self.universal_reactions[rxn_id]
            if not urxn.is_usable:
                continue
            base_id = base_seed_rxn_id(rxn_id)
            if base_id in self.draft_base_ids:
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
            u_model.add_metabolites(new_mets)

        split_map: Dict[str, Tuple[str, ...]] = {}
        new_rxns: List[cobra.Reaction] = []
        for rxn_id in sorted(self.universal_reactions):
            urxn = self.universal_reactions[rxn_id]
            if not urxn.is_usable:
                continue
            base_id = base_seed_rxn_id(rxn_id)
            if base_id in self.draft_base_ids or not urxn.stoichiometry:
                continue
            lb, ub = urxn.bounds
            split_ids: List[str] = []
            if ub > 0.0:
                fwd_id = f"{base_id}_c0__cand_fwd"
                fwd_rxn = cobra.Reaction(
                    id=fwd_id,
                    name=f"{urxn.name} (candidate fwd)",
                    lower_bound=0.0,
                    upper_bound=ub,
                )
                fwd_rxn.add_metabolites(
                    {existing_mets[t.species_id]: t.coeff for t in urxn.stoichiometry}
                )
                new_rxns.append(fwd_rxn)
                split_ids.append(fwd_id)
            if lb < 0.0:
                rev_id = f"{base_id}_c0__cand_rev"
                rev_rxn = cobra.Reaction(
                    id=rev_id,
                    name=f"{urxn.name} (candidate rev)",
                    lower_bound=0.0,
                    upper_bound=abs(lb),
                )
                rev_rxn.add_metabolites(
                    {existing_mets[t.species_id]: -t.coeff for t in urxn.stoichiometry}
                )
                new_rxns.append(rev_rxn)
                split_ids.append(rev_id)
            if split_ids:
                split_map[base_id] = tuple(split_ids)

        if new_rxns:
            u_model.add_reactions(new_rxns)

        # Ensure every extracellular metabolite has an exchange reaction (closed by default)
        ex_ids = [
            f"EX_{m.id}"
            for m in u_model.metabolites
            if m.id.endswith("_e0") and f"EX_{m.id}" not in u_model.reactions
        ]
        _ensure_exchange_reactions(u_model, ex_ids)

        return u_model, split_map

    def check_draft_flux(
        self,
        target_rxn_id: str,
        medium: Mapping[str, float],
        sense: str = "max",
    ) -> Tuple[float, int, float, float]:
        """Optimize `target_rxn_id` on `self.draft_model` under `medium`.

        Returns `(opt_flux, n_lps, cpu_s, wall_s)`.
        """
        t_cpu0 = time.process_time()
        t_wall0 = time.monotonic()
        apply_medium(self.draft_model, medium)
        if target_rxn_id not in self.draft_model.reactions:
            return 0.0, 0, time.process_time() - t_cpu0, time.monotonic() - t_wall0
        rxn = self.draft_model.reactions.get_by_id(target_rxn_id)
        orig_dir = self.draft_model.objective_direction
        try:
            coeffs: Dict[object, float] = {v: 0.0 for v in self._draft_active_vars}
            coeffs[rxn.forward_variable] = 1.0
            coeffs[rxn.reverse_variable] = -1.0
            self.draft_model.solver.objective.set_linear_coefficients(coeffs)
            self._draft_active_vars = [rxn.forward_variable, rxn.reverse_variable]
            if self.draft_model.objective_direction != sense:
                self.draft_model.objective_direction = sense
            sol = self.draft_model.slim_optimize(error_value=0.0)
            flux = 0.0 if sol is None or math.isnan(sol) else float(sol)
            return flux, 1, time.process_time() - t_cpu0, time.monotonic() - t_wall0
        finally:
            self.draft_model.solver.objective.set_linear_coefficients(
                {rxn.forward_variable: 0.0, rxn.reverse_variable: 0.0}
            )
            self._draft_active_vars.clear()
            if self.draft_model.objective_direction != orig_dir:
                self.draft_model.objective_direction = orig_dir

    def solve_weighted_gapfill(
        self,
        test_id: str,
        arm: str,
        target_rxn_id: str,
        medium: Mapping[str, float],
        check_draft_first: bool = True,
        target_direction: str = "forward",
        temporary_rxn_spec: Optional[Tuple[str, Dict[str, float], float, float]] = None,
        replay_fraction: float = 0.05,
        force_replay: bool = False,
    ) -> GapFillTestResult:
        """Run the two-step weighted gap-filling LP (§4) for a single target."""
        t_cpu0 = time.process_time()
        t_wall0 = time.monotonic()
        n_lps = 0
        n_gapfills = 0

        if check_draft_first and temporary_rxn_spec is None:
            sense = "max" if target_direction == "forward" else "min"
            d_flux, dlps, _, _ = self.check_draft_flux(target_rxn_id, medium, sense=sense)
            n_lps += dlps
            if (target_direction == "forward" and d_flux >= self.flux_tol) or (
                target_direction == "reverse" and d_flux <= -self.flux_tol
            ):
                return GapFillTestResult(
                    test_id=test_id,
                    arm=arm,
                    target_id=target_rxn_id,
                    feasible_in_universal=True,
                    already_satisfied_in_draft=True,
                    gapfill_ran=False,
                    vmax_universal=abs(d_flux),
                    objective_value=0.0,
                    added_reactions=(),
                    n_lps=n_lps,
                    n_gapfills=0,
                    cpu_seconds=time.process_time() - t_cpu0,
                    wall_seconds=time.monotonic() - t_wall0,
                    replay=None,
                )

        should_replay = _should_replay_test(
            test_id, replay_fraction, force_replay=force_replay
        )
        vmax, obj_val, added_rxns, ulps, ugfs, replay_record = self._run_lp_on_models(
            vmax_model=self.vmax_model,
            vmax_active_vars=self._vmax_active_vars,
            u_model=self.universal_model,
            cand_var_entries=self._cand_var_entries,
            target_rxn_id=target_rxn_id,
            medium=medium,
            target_direction=target_direction,
            temporary_rxn_spec=temporary_rxn_spec,
            test_id=test_id,
            should_replay=should_replay,
        )
        n_lps += ulps
        n_gapfills += ugfs

        return GapFillTestResult(
            test_id=test_id,
            arm=arm,
            target_id=target_rxn_id,
            feasible_in_universal=vmax >= self.flux_tol,
            already_satisfied_in_draft=False,
            gapfill_ran=ugfs > 0,
            vmax_universal=vmax,
            objective_value=obj_val,
            added_reactions=added_rxns,
            n_lps=n_lps,
            n_gapfills=n_gapfills,
            cpu_seconds=time.process_time() - t_cpu0,
            wall_seconds=time.monotonic() - t_wall0,
            replay=replay_record,
        )

    def _run_lp_on_models(
        self,
        vmax_model: cobra.Model,
        vmax_active_vars: List[object],
        u_model: cobra.Model,
        cand_var_entries: Sequence[Tuple[str, float, Tuple[object, ...]]],
        target_rxn_id: str,
        medium: Mapping[str, float],
        target_direction: str,
        temporary_rxn_spec: Optional[Tuple[str, Dict[str, float], float, float]],
        test_id: str = "",
        should_replay: bool = False,
    ) -> Tuple[float, float, Tuple[str, ...], int, int, Optional[SolverReplayRecord]]:
        """Execute Step 1 ($v_t^{\\max}$ on `vmax_model`) and Step 2 ($\\min \\sum w_r |v_r|$ on `u_model`)."""
        apply_medium(vmax_model, medium)

        v_temp_rxn: Optional[cobra.Reaction] = None
        if temporary_rxn_spec is not None:
            t_id, stoich_ids, t_lb, t_ub = temporary_rxn_spec
            if not all(m_id in vmax_model.metabolites for m_id in stoich_ids):
                return 0.0, 0.0, (), 0, 0, None
            v_temp_rxn = cobra.Reaction(id=t_id, name=t_id, lower_bound=t_lb, upper_bound=t_ub)
            v_temp_rxn.add_metabolites(
                {vmax_model.metabolites.get_by_id(m_id): coeff for m_id, coeff in stoich_ids.items()}
            )
            vmax_model.add_reactions([v_temp_rxn])
            target_rxn_id = t_id

        if target_rxn_id not in vmax_model.reactions:
            if v_temp_rxn is not None:
                vmax_model.remove_reactions([v_temp_rxn], remove_orphans=False)
            return 0.0, 0.0, (), 0, 0, None

        v_target = vmax_model.reactions.get_by_id(target_rxn_id)
        orig_v_dir = vmax_model.objective_direction
        rep_vmax = 0.0
        try:
            # Step 1: Maximize flux through target in the universal-completed model
            v_coeffs: Dict[object, float] = {v: 0.0 for v in vmax_active_vars}
            v_coeffs[v_target.forward_variable] = 1.0
            v_coeffs[v_target.reverse_variable] = -1.0
            vmax_model.solver.objective.set_linear_coefficients(v_coeffs)
            vmax_active_vars.clear()
            vmax_active_vars.extend([v_target.forward_variable, v_target.reverse_variable])

            want_dir = "max" if target_direction == "forward" else "min"
            if vmax_model.objective_direction != want_dir:
                vmax_model.objective_direction = want_dir

            raw_opt = vmax_model.slim_optimize(error_value=0.0)
            n_lps = 1
            if raw_opt is None or math.isnan(raw_opt):
                vmax = 0.0
            else:
                vmax = float(raw_opt) if target_direction == "forward" else float(-raw_opt)

            if should_replay and vmax >= self.flux_tol:
                rep_ok, rep_raw_opt, _, _ = _solve_model_lp_with_highs(
                    vmax_model, zero_tol=self.flux_tol
                )
                n_lps += 1
                if not rep_ok or math.isnan(rep_raw_opt):
                    rep_vmax = 0.0
                else:
                    rep_vmax = (
                        float(rep_raw_opt)
                        if target_direction == "forward"
                        else float(-rep_raw_opt)
                    )
        finally:
            if v_temp_rxn is not None:
                for v in (v_temp_rxn.forward_variable, v_temp_rxn.reverse_variable):
                    if v in vmax_active_vars:
                        vmax_active_vars.remove(v)
                vmax_model.remove_reactions([v_temp_rxn], remove_orphans=False)
            else:
                vmax_model.solver.objective.set_linear_coefficients(
                    {v_target.forward_variable: 0.0, v_target.reverse_variable: 0.0}
                )
                vmax_active_vars.clear()
            if vmax_model.objective_direction != orig_v_dir:
                vmax_model.objective_direction = orig_v_dir

        if vmax < self.flux_tol:
            return 0.0, 0.0, (), n_lps, 0, None

        # Step 2: Constrain v_t >= epsilon and minimize sum_{r in U\D} w_r |v_r| on u_model
        apply_medium(u_model, medium)
        u_temp_rxn: Optional[cobra.Reaction] = None
        if temporary_rxn_spec is not None:
            t_id, stoich_ids, t_lb, t_ub = temporary_rxn_spec
            u_temp_rxn = cobra.Reaction(id=t_id, name=t_id, lower_bound=t_lb, upper_bound=t_ub)
            u_temp_rxn.add_metabolites(
                {u_model.metabolites.get_by_id(m_id): coeff for m_id, coeff in stoich_ids.items()}
            )
            u_model.add_reactions([u_temp_rxn])
            target_rxn_id = t_id

        u_target = u_model.reactions.get_by_id(target_rxn_id)
        orig_lb, orig_ub = u_target.lower_bound, u_target.upper_bound
        try:
            epsilon = max(self.flux_tol, self.epsilon_frac * vmax)
            if target_direction == "forward":
                u_target.lower_bound = max(orig_lb, epsilon)
            else:
                u_target.upper_bound = min(orig_ub, -epsilon)

            opt_val = u_model.slim_optimize(error_value=float("nan"))
            n_lps += 1
            if u_model.solver.status != "optimal" and (opt_val is None or math.isnan(opt_val)):
                added_tuple: Tuple[str, ...] = ()
                total_penalty = 0.0
            else:
                added: List[str] = []
                total_penalty = 0.0
                for base_id, w_r, fvars in cand_var_entries:
                    flux_sum = 0.0
                    for fvar in fvars:
                        val = float(fvar.primal or 0.0)
                        if val > 0.0:
                            flux_sum += val
                    if flux_sum > self.flux_tol:
                        added.append(base_id)
                        total_penalty += w_r * flux_sum
                added_tuple = tuple(sorted(added))

            replay_record: Optional[SolverReplayRecord] = None
            if should_replay:
                u_target.lower_bound = orig_lb
                u_target.upper_bound = orig_ub
                if rep_vmax >= self.flux_tol:
                    rep_epsilon = max(self.flux_tol, self.epsilon_frac * rep_vmax)
                    if target_direction == "forward":
                        u_target.lower_bound = max(orig_lb, rep_epsilon)
                    else:
                        u_target.upper_bound = min(orig_ub, -rep_epsilon)

                    override_restore: Dict[object, float] = {}
                    if self.replay_weight_overrides:
                        override_apply: Dict[object, float] = {}
                        for base_id, w_r, fvars in cand_var_entries:
                            if base_id in self.replay_weight_overrides:
                                new_w = float(self.replay_weight_overrides[base_id])
                                for fvar in fvars:
                                    override_apply[fvar] = -new_w
                                    override_restore[fvar] = -w_r
                        if override_apply:
                            u_model.solver.objective.set_linear_coefficients(override_apply)
                    try:
                        rep_ok2, _, rep_x, var_idx_by_name = _solve_model_lp_with_highs(
                            u_model, zero_tol=self.flux_tol
                        )
                        n_lps += 1
                    finally:
                        if override_restore:
                            u_model.solver.objective.set_linear_coefficients(override_restore)

                    if not rep_ok2:
                        rep_obj = 0.0
                        rep_added: Tuple[str, ...] = ()
                    else:
                        rep_added_list: List[str] = []
                        rep_obj = 0.0
                        for base_id, w_r, fvars in cand_var_entries:
                            flux_sum = 0.0
                            for fvar in fvars:
                                v_idx = var_idx_by_name.get(fvar.name)
                                if v_idx is not None:
                                    val = float(rep_x[v_idx])
                                    if val > 0.0:
                                        flux_sum += val
                            if flux_sum > self.flux_tol:
                                rep_added_list.append(base_id)
                                eff_w = (
                                    float(self.replay_weight_overrides[base_id])
                                    if base_id in self.replay_weight_overrides
                                    else w_r
                                )
                                rep_obj += eff_w * flux_sum
                        rep_added = tuple(sorted(rep_added_list))
                else:
                    rep_obj = 0.0
                    rep_added = ()

                rel_diff = abs(total_penalty - rep_obj) / max(1.0, abs(total_penalty))
                agreed = (added_tuple == rep_added) and (rel_diff <= 1e-6)
                replay_record = SolverReplayRecord(
                    test_id=test_id,
                    primary_solver=self.solver_info.primary_solver,
                    replay_solver="scipy-linprog-highs",
                    primary_objective=total_penalty,
                    replay_objective=rep_obj,
                    relative_obj_diff=rel_diff,
                    primary_added_rxns=added_tuple,
                    replay_added_rxns=rep_added,
                    agreed=agreed,
                    scipy_version=str(getattr(scipy, "__version__", "")),
                )

            return vmax, total_penalty, added_tuple, n_lps, 1, replay_record
        finally:
            u_target.lower_bound = orig_lb
            u_target.upper_bound = orig_ub
            if u_temp_rxn is not None:
                u_model.remove_reactions([u_temp_rxn], remove_orphans=False)


def identify_dead_end_targets(
    draft_model: cobra.Model,
) -> Tuple[Tuple[str, str], ...]:
    """Identify intracellular dead-end metabolites in gene-evidenced reactions (§5 A3).

    Returns sorted `(metabolite_id, dead_end_kind)` pairs where `dead_end_kind` is:
      - `"unproduced"`: consumed by >=1 gene-evidenced reaction, produced by 0 -> gets sink target
      - `"unconsumed"`: produced by >=1 gene-evidenced reaction, consumed by 0 -> gets source target
    """
    evidenced_rxns = [
        r
        for r in draft_model.reactions
        if not r.id.startswith("EX_")
        and r.id not in {"bio1", "R_bio1"}
        and bool((r.gene_reaction_rule or "").strip())
    ]
    produced_by: Dict[str, int] = {}
    consumed_by: Dict[str, int] = {}
    candidate_mets: Set[str] = set()

    for rxn in evidenced_rxns:
        for met, coeff in rxn.metabolites.items():
            if met.id.endswith("_e0"):
                continue
            candidate_mets.add(met.id)
            can_produce = (coeff > 0 and rxn.upper_bound > 0) or (
                coeff < 0 and rxn.lower_bound < 0
            )
            can_consume = (coeff < 0 and rxn.upper_bound > 0) or (
                coeff > 0 and rxn.lower_bound < 0
            )
            if can_produce:
                produced_by[met.id] = produced_by.get(met.id, 0) + 1
            if can_consume:
                consumed_by[met.id] = consumed_by.get(met.id, 0) + 1

    targets: List[Tuple[str, str]] = []
    for met_id in sorted(candidate_mets):
        n_prod = produced_by.get(met_id, 0)
        n_cons = consumed_by.get(met_id, 0)
        if n_cons > 0 and n_prod == 0:
            targets.append((met_id, "unproduced"))
        elif n_prod > 0 and n_cons == 0:
            targets.append((met_id, "unconsumed"))
    return tuple(targets)


def _build_fill_path_index(
    draft_model: cobra.Model,
    universal_reactions: Mapping[str, SeedReaction],
    rxn_weights: Mapping[str, float],
    union_medium: Mapping[str, float],
) -> Tuple[
    Set[str],
    Set[str],
    Dict[str, List[Tuple[str, float, str]]],
    Dict[str, List[str]],
    Dict[str, float],
    Dict[str, Tuple[str, str]],
]:
    """Pre-build directed metabolite graph and multi-source forward Dijkstra once per case."""
    draft_base_ids = {
        base_seed_rxn_id(r.id)
        for r in draft_model.reactions
        if not r.id.startswith("EX_") and bool((r.gene_reaction_rule or "").strip())
    }
    adj: Dict[str, List[Tuple[str, float, str]]] = {}
    species_to_cand_rxns: Dict[str, List[str]] = {}

    def add_directed_edges(
        substrates: Sequence[str],
        products: Sequence[str],
        cost: float,
        base_id: str,
    ) -> None:
        for s in substrates:
            bucket = adj.setdefault(s, [])
            for p in products:
                if s != p:
                    bucket.append((p, cost, base_id))

    for rxn in draft_model.reactions:
        if rxn.id.startswith("EX_") or not (rxn.gene_reaction_rule or "").strip():
            continue
        subs = [m.id for m, c in rxn.metabolites.items() if c < 0]
        prods = [m.id for m, c in rxn.metabolites.items() if c > 0]
        base_id = base_seed_rxn_id(rxn.id)
        if rxn.upper_bound > 0:
            add_directed_edges(subs, prods, 0.0, base_id)
        if rxn.lower_bound < 0:
            add_directed_edges(prods, subs, 0.0, base_id)

    for rxn_id in sorted(universal_reactions):
        urxn = universal_reactions[rxn_id]
        if not urxn.is_usable:
            continue
        base_id = base_seed_rxn_id(rxn_id)
        if base_id in draft_base_ids:
            continue
        cost = rxn_weights.get(base_id, WEIGHT_NO_CANDIDATE)
        subs = [t.species_id for t in urxn.stoichiometry if t.coeff < 0]
        prods = [t.species_id for t in urxn.stoichiometry if t.coeff > 0]
        for t in urxn.stoichiometry:
            species_to_cand_rxns.setdefault(t.species_id, []).append(base_id)
        lb, ub = urxn.bounds
        if ub > 0:
            add_directed_edges(subs, prods, cost, base_id)
        if lb < 0:
            add_directed_edges(prods, subs, cost, base_id)

    for u_node, edges in adj.items():
        edges.sort()

    evidenced_mets: Set[str] = set()
    for ex_id, flux in union_medium.items():
        if flux > 0 and ex_id.startswith("EX_"):
            evidenced_mets.add(ex_id[len("EX_") :])
    for rxn in draft_model.reactions:
        if not rxn.id.startswith("EX_") and (rxn.gene_reaction_rule or "").strip():
            for m in rxn.metabolites:
                evidenced_mets.add(m.id)

    # Precompute multi-source forward Dijkstra from all evidenced_mets once
    prod_dist: Dict[str, float] = {}
    prod_parent: Dict[str, Tuple[str, str]] = {}
    pq: List[Tuple[float, str]] = []
    for s in sorted(evidenced_mets):
        prod_dist[s] = 0.0
        heapq.heappush(pq, (0.0, s))

    while pq:
        d, u = heapq.heappop(pq)
        if d > prod_dist.get(u, float("inf")):
            continue
        for v, edge_cost, rxn_base in adj.get(u, ()):
            nd = d + edge_cost
            if nd < prod_dist.get(v, float("inf")):
                prod_dist[v] = nd
                prod_parent[v] = (u, rxn_base)
                heapq.heappush(pq, (nd, v))

    return (
        draft_base_ids,
        evidenced_mets,
        adj,
        species_to_cand_rxns,
        prod_dist,
        prod_parent,
    )


def _cheapest_fill_path_genes(
    draft_model: cobra.Model,
    universal_reactions: Mapping[str, SeedReaction],
    rxn_weights: Mapping[str, float],
    rxn_to_candidates: Mapping[str, Sequence[Tuple[str, float]]],
    union_medium: Mapping[str, float],
    target_species: Iterable[str],
    direction: str = "produce",
    prebuilt_index: Optional[
        Tuple[
            Set[str],
            Set[str],
            Dict[str, List[Tuple[str, float, str]]],
            Dict[str, List[str]],
            Dict[str, float],
            Dict[str, Tuple[str, str]],
        ]
    ] = None,
) -> Tuple[str, ...]:
    """Deterministic 0-LP Dijkstra shortest path on the metabolite-reaction graph (§5.2 Rule 3).

    Finds the minimum-weight reaction path in the universal-completed network
    either producing `target_species` from evidenced/medium metabolites
    (`direction='produce'`) or consuming `target_species` into evidenced/sink
    metabolites (`direction='consume'`), and returns sorted distinct candidate
    `protein_id`s supporting candidate reactions ($w_r = 1$) on that path.
    """
    targets = sorted(set(target_species))
    if not targets:
        return ()

    if prebuilt_index is None:
        prebuilt_index = _build_fill_path_index(
            draft_model, universal_reactions, rxn_weights, union_medium
        )
    (
        _draft_base_ids,
        evidenced_mets,
        adj,
        species_to_cand_rxns,
        prod_dist,
        prod_parent,
    ) = prebuilt_index

    tset = set(targets)
    reached_goal: Optional[str] = None
    parent_rxn: Dict[str, Tuple[str, str]] = {}

    if direction == "produce" and not (tset & evidenced_mets):
        # Fast lookup from precomputed multi-source forward Dijkstra
        best_t: Optional[str] = None
        best_d = float("inf")
        for t in targets:
            d = prod_dist.get(t, float("inf"))
            if 0.0 < d < best_d:
                best_d = d
                best_t = t
        reached_goal = best_t
        parent_rxn = prod_parent
    else:
        if direction == "produce":
            start_nodes = sorted(evidenced_mets - tset) or sorted(evidenced_mets)
            goal_set = tset
        else:
            start_nodes = targets
            goal_set = (evidenced_mets - tset) | {
                m.id for m in draft_model.metabolites if m.id.endswith("_e0")
            }

        dist: Dict[str, float] = {}
        pq: List[Tuple[float, str]] = []
        for s in start_nodes:
            dist[s] = 0.0
            heapq.heappush(pq, (0.0, s))

        while pq:
            d, u = heapq.heappop(pq)
            if d > dist.get(u, float("inf")):
                continue
            if u in goal_set and d > 0.0:
                reached_goal = u
                break
            for v, edge_cost, rxn_base in adj.get(u, ()):
                nd = d + edge_cost
                if nd < dist.get(v, float("inf")):
                    dist[v] = nd
                    parent_rxn[v] = (u, rxn_base)
                    heapq.heappush(pq, (nd, v))

    if reached_goal is None:
        direct_genes: Set[str] = set()
        for sp in targets:
            for base_id in species_to_cand_rxns.get(sp, ()):
                for gene_id, _ in rxn_to_candidates.get(base_id, ()):
                    direct_genes.add(gene_id)
        return tuple(sorted(direct_genes))

    genes: Set[str] = set()
    cur = reached_goal
    visited: Set[str] = set()
    while cur in parent_rxn and cur not in visited:
        visited.add(cur)
        prev, rxn_base = parent_rxn[cur]
        for gene_id, _ in rxn_to_candidates.get(rxn_base, ()):
            genes.add(gene_id)
        cur = prev
    return tuple(sorted(genes))


def rank_a4_tests(
    test_descriptors: Sequence[Tuple[str, str, str, Optional[PanelEntry]]],
    draft_model: cobra.Model,
    universal_reactions: Mapping[str, SeedReaction],
    rxn_weights: Mapping[str, float],
    rxn_to_candidates: Mapping[str, Sequence[Tuple[str, float]]],
    transporters: Sequence[TransporterHit],
    kegg_modules: Sequence[KeggModuleDef],
    union_medium: Mapping[str, float],
) -> Tuple[RankedTestEntry, ...]:
    """Deterministically rank all A1-A3 tests of a genome under fixed rules 1-4 (§5.2).

    `test_descriptors` is a sequence of `(test_id, arm, target_id, panel_entry)`.
    Rule order:
      1. A1 source compound named by a genome transporter annotation (GSPA or Transporter.tbl).
      2. Target lies in a KEGG module that is >=50% evidenced in the masked genome.
      3. A2/A3 target with >=1 masked GSPA candidate gene on the cheapest fill path.
      4. All remaining tests.
    Within each rule tier, rank by descending count of distinct supporting
    candidate genes, breaking ties by ascending `test_id`.
    """
    draft_evidenced_base_rxns: Set[str] = {
        base_seed_rxn_id(r.id)
        for r in draft_model.reactions
        if not r.id.startswith("EX_") and bool((r.gene_reaction_rule or "").strip())
    }
    draft_evidenced_ecs: Set[str] = set()
    for base_id in draft_evidenced_base_rxns:
        urxn = universal_reactions.get(base_id)
        if urxn is not None:
            draft_evidenced_ecs.update(urxn.ec_numbers)

    ec_to_univ_rxns: Dict[str, List[str]] = {}
    for r_id, urxn in universal_reactions.items():
        b_id = base_seed_rxn_id(r_id)
        for ec in urxn.ec_numbers:
            ec_to_univ_rxns.setdefault(ec, []).append(b_id)

    # Identify KEGG modules that are >= 50% evidenced in the masked genome
    half_evidenced_modules: List[KeggModuleDef] = []
    module_rxns_map: Dict[str, Set[str]] = {}
    for mod in kegg_modules:
        mod_rxns: Set[str] = set(mod.reaction_ids)
        for ec in mod.ec_numbers:
            mod_rxns.update(ec_to_univ_rxns.get(ec, ()))
        module_rxns_map[mod.module_id] = mod_rxns

        if mod.ec_numbers:
            ev_count = sum(1 for ec in mod.ec_numbers if ec in draft_evidenced_ecs)
            frac = ev_count / float(len(mod.ec_numbers))
        elif mod_rxns:
            ev_count = sum(1 for r in mod_rxns if r in draft_evidenced_base_rxns)
            frac = ev_count / float(len(mod_rxns))
        else:
            frac = 0.0
        if frac >= 0.5:
            half_evidenced_modules.append(mod)

    # Pre-index species participating in >=50% evidenced KEGG modules
    species_to_half_mod_rxns: Dict[str, Set[str]] = {}
    rxn_to_half_mod_rxns: Dict[str, Set[str]] = {}
    for mod in half_evidenced_modules:
        m_rxns = module_rxns_map.get(mod.module_id, set())
        for r_base in m_rxns:
            rxn_to_half_mod_rxns.setdefault(r_base, set()).update(m_rxns)
            urxn = universal_reactions.get(r_base)
            if urxn is not None:
                for term in urxn.stoichiometry:
                    species_to_half_mod_rxns.setdefault(term.species_id, set()).update(m_rxns)

    # Pre-index transporter hits by cpd_id and exchange_rxn
    trans_genes_by_cpd: Dict[str, Set[str]] = {}
    for hit in transporters:
        if hit.seed_id:
            trans_genes_by_cpd.setdefault(hit.seed_id, set()).add(hit.protein_id)
        for r_id in hit.reaction_ids:
            base_r = base_seed_rxn_id(r_id)
            urxn = universal_reactions.get(base_r)
            if urxn is not None:
                for term in urxn.stoichiometry:
                    if term.compartment == "e0":
                        trans_genes_by_cpd.setdefault(term.cpd_id, set()).add(hit.protein_id)

    prebuilt_path_idx = _build_fill_path_index(
        draft_model, universal_reactions, rxn_weights, union_medium
    )

    unranked: List[Tuple[int, int, str, str, str, Tuple[str, ...]]] = []

    for test_id, arm, target_id, panel_entry in test_descriptors:
        rule_tier = 4
        supporting_genes: Set[str] = set()

        if arm == "A1" and panel_entry is not None:
            cpd_id = panel_entry.seed_id
            t_genes = set(trans_genes_by_cpd.get(cpd_id, set()))
            for tr_id in panel_entry.transport_rxns:
                for g_id, _ in rxn_to_candidates.get(base_seed_rxn_id(tr_id), ()):
                    t_genes.add(g_id)
            if t_genes:
                rule_tier = 1
                supporting_genes.update(t_genes)
            else:
                c0_species = f"{cpd_id}_c0"
                mod_rxns = species_to_half_mod_rxns.get(c0_species, set())
                if mod_rxns:
                    rule_tier = 2
                    for mr in mod_rxns:
                        if mr not in draft_evidenced_base_rxns:
                            for g_id, _ in rxn_to_candidates.get(mr, ()):
                                supporting_genes.add(g_id)

        elif arm == "A2":
            base_target = base_seed_rxn_id(target_id)
            target_species: List[str] = []
            if target_id in draft_model.reactions:
                drxn = draft_model.reactions.get_by_id(target_id)
                target_species = [m.id for m in drxn.metabolites if not m.id.endswith("_e0")]

            mod_rxns = set(rxn_to_half_mod_rxns.get(base_target, set()))
            for sp in target_species:
                mod_rxns.update(species_to_half_mod_rxns.get(sp, set()))

            path_genes = set(
                _cheapest_fill_path_genes(
                    draft_model=draft_model,
                    universal_reactions=universal_reactions,
                    rxn_weights=rxn_weights,
                    rxn_to_candidates=rxn_to_candidates,
                    union_medium=union_medium,
                    target_species=target_species,
                    direction="produce",
                    prebuilt_index=prebuilt_path_idx,
                )
            )

            if mod_rxns:
                rule_tier = 2
                for mr in mod_rxns:
                    if mr not in draft_evidenced_base_rxns:
                        for g_id, _ in rxn_to_candidates.get(mr, ()):
                            supporting_genes.add(g_id)
                supporting_genes.update(path_genes)
            elif path_genes:
                rule_tier = 3
                supporting_genes.update(path_genes)

        elif arm == "A3":
            # target_id is formatted as "<species_id>:<unproduced|unconsumed>" or "<species_id>"
            if ":" in target_id:
                sp_id, kind = target_id.split(":", 1)
            else:
                sp_id, kind = target_id, "unproduced"
            direction = "produce" if kind == "unproduced" else "consume"

            mod_rxns = set(species_to_half_mod_rxns.get(sp_id, set()))
            path_genes = set(
                _cheapest_fill_path_genes(
                    draft_model=draft_model,
                    universal_reactions=universal_reactions,
                    rxn_weights=rxn_weights,
                    rxn_to_candidates=rxn_to_candidates,
                    union_medium=union_medium,
                    target_species=[sp_id],
                    direction=direction,
                    prebuilt_index=prebuilt_path_idx,
                )
            )
            if mod_rxns:
                rule_tier = 2
                for mr in mod_rxns:
                    if mr not in draft_evidenced_base_rxns:
                        for g_id, _ in rxn_to_candidates.get(mr, ()):
                            supporting_genes.add(g_id)
                supporting_genes.update(path_genes)
            elif path_genes:
                rule_tier = 3
                supporting_genes.update(path_genes)

        sorted_genes = tuple(sorted(supporting_genes))
        unranked.append(
            (rule_tier, -len(sorted_genes), test_id, arm, target_id, sorted_genes)
        )

    unranked.sort(key=lambda item: (item[0], item[1], item[2]))
    ranked: List[RankedTestEntry] = []
    for idx, (tier, neg_cnt, test_id, arm, target_id, genes) in enumerate(unranked):
        ranked.append(
            RankedTestEntry(
                test_id=test_id,
                arm=arm,
                target_id=target_id,
                rule_tier=tier,
                supporting_genes=genes,
                n_supporting_genes=-neg_cnt,
                rank_index=idx,
            )
        )
    return tuple(ranked)


_WORKER_FILLER: Optional[WeightedGapFiller] = None
_WORKER_UNION_MED: Dict[str, float] = {}
_WORKER_REPLAY_FRAC: float = 0.05
_WORKER_FLUX_TOL: float = DEFAULT_FLUX_TOL


def _init_single_thread_worker() -> None:
    """Pin BLAS/OpenMP thread counts to 1 inside each forked worker process."""
    for var in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
    ):
        os.environ[var] = "1"


def _worker_check_blocked(
    task: Tuple[int, str, bool],
) -> Tuple[int, str, bool, bool, int, float, float]:
    idx, rxn_id, is_rev = task
    assert _WORKER_FILLER is not None
    f_max, lps1, cpu1, wall1 = _WORKER_FILLER.check_draft_flux(
        rxn_id, _WORKER_UNION_MED, sense="max"
    )
    total_lps = lps1
    total_cpu = cpu1
    total_wall = wall1
    unblocked = f_max >= _WORKER_FLUX_TOL
    if not unblocked and is_rev:
        f_min, lps2, cpu2, wall2 = _WORKER_FILLER.check_draft_flux(
            rxn_id, _WORKER_UNION_MED, sense="min"
        )
        total_lps += lps2
        total_cpu += cpu2
        total_wall += wall2
        if f_min <= -_WORKER_FLUX_TOL:
            unblocked = True
    return idx, rxn_id, is_rev, unblocked, total_lps, total_cpu, total_wall


def _worker_solve_test(
    task: Tuple[
        int,
        str,
        str,
        str,
        Dict[str, float],
        bool,
        bool,
        Optional[Tuple[str, Dict[str, float], float, float]],
    ],
) -> Tuple[int, str, GapFillTestResult]:
    (
        task_idx,
        test_id,
        arm,
        target_rxn_id,
        med,
        check_draft_first,
        try_reverse_if_infeasible,
        temp_spec,
    ) = task
    assert _WORKER_FILLER is not None
    res = _WORKER_FILLER.solve_weighted_gapfill(
        test_id=test_id,
        arm=arm,
        target_rxn_id=target_rxn_id,
        medium=med,
        check_draft_first=check_draft_first,
        target_direction="forward",
        temporary_rxn_spec=temp_spec,
        replay_fraction=_WORKER_REPLAY_FRAC,
    )
    if arm == "A2" and not res.feasible_in_universal and try_reverse_if_infeasible:
        res_rev = _WORKER_FILLER.solve_weighted_gapfill(
            test_id=test_id,
            arm="A2",
            target_rxn_id=target_rxn_id,
            medium=med,
            check_draft_first=False,
            target_direction="reverse",
            replay_fraction=_WORKER_REPLAY_FRAC,
        )
        res = GapFillTestResult(
            test_id=test_id,
            arm="A2",
            target_id=target_rxn_id,
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
    return task_idx, arm, res


def _ordered_with_progress(pool, func, tasks, label: str, every: int = 50):
    """pool.map with the same deterministic order, printing progress lines (count, elapsed) to stdout."""
    import sys as _sys
    t0 = time.perf_counter()
    out = []
    n = len(tasks)
    print(f"[progress] {label}: start {n}", flush=True)
    for i, res in enumerate(pool.imap(func, tasks, chunksize=1), 1):
        out.append(res)
        if i % every == 0 or i == n:
            print(f"[progress] {label}: {i}/{n} done, {time.perf_counter() - t0:.1f} s", flush=True)
    _sys.stdout.flush()
    return out


def run_all_arms_for_case(
    masked_model: cobra.Model,
    audit_report: Optional[MaskingAuditReport],
    universal_reactions: Mapping[str, SeedReaction],
    metabolites: Mapping[str, SeedMetabolite],
    candidates: Sequence[CandidateClaim],
    ec_to_reactions: Mapping[str, Sequence[str]],
    panel_entries: Sequence[PanelEntry],
    default_medium: Mapping[str, float],
    transporters: Sequence[TransporterHit] = (),
    kegg_modules: Sequence[KeggModuleDef] = (),
    replay_fraction: float = 0.05,
    flux_tol: float = DEFAULT_FLUX_TOL,
    allow_unmasked: bool = False,
    workers: int = 1,
    max_tests_per_arm: Optional[int] = None,
) -> Dict[str, ArmRunResult]:
    """Run A0, A1, A2, A3, A1_A2_A3, A4@10%, and A4@25% on a case (§4, §5).

    Enforces `assert_masking_clean(audit_report)` before any LP is solved unless
    `allow_unmasked=True` (used only for the §8 unmasked control run).
    When `workers > 1`, dispatches FVA and gap-fill tests across a forked process
    pool with single-threaded LP per worker and deterministic output order.
    """
    if not allow_unmasked:
        if audit_report is None:
            raise ValueError("audit_report is required when allow_unmasked is False.")
        assert_masking_clean(audit_report)

    rxn_to_candidates = map_reactions_to_candidate_genes(candidates, ec_to_reactions)
    draft_base_ids = {
        base_seed_rxn_id(r.id)
        for r in masked_model.reactions
        if not r.id.startswith("EX_")
    }
    rxn_weights = compute_reaction_weights(
        universal_reactions, draft_base_ids, rxn_to_candidates
    )
    filler = WeightedGapFiller(
        draft_model=masked_model,
        universal_reactions=universal_reactions,
        rxn_weights=rxn_weights,
        flux_tol=flux_tol,
    )
    biomass_id = _find_biomass_rxn_id(masked_model)
    union_med = build_union_medium(default_medium, panel_entries)

    # Pre-warm union medium before any worker fork
    apply_medium(filler.draft_model, union_med)
    apply_medium(filler.vmax_model, union_med)
    apply_medium(filler.universal_model, union_med)

    # --- Arm A0: default medium biomass gap-fill ---
    a0_res = filler.solve_weighted_gapfill(
        test_id="A0_biomass_default",
        arm="A0",
        target_rxn_id=biomass_id,
        medium=default_medium,
        check_draft_first=True,
        replay_fraction=replay_fraction,
    )
    arm_a0 = _summarize_arm("A0", [a0_res])

    # Re-apply union_med after A0 default_medium check
    apply_medium(filler.draft_model, union_med)
    apply_medium(filler.vmax_model, union_med)
    apply_medium(filler.universal_model, union_med)

    # Prepare A1 descriptors
    eff_panel = (
        panel_entries[:max_tests_per_arm]
        if (max_tests_per_arm is not None and max_tests_per_arm > 0)
        else panel_entries
    )
    a1_descriptors: List[Tuple[str, str, str, Optional[PanelEntry]]] = []
    for entry in eff_panel:
        test_id = f"A1_{entry.panel_id}"
        a1_descriptors.append((test_id, "A1", biomass_id, entry))

    # Prepare A2 FVA tasks on gene-evidenced reactions
    evidenced_rxns = sorted(
        [
            r
            for r in masked_model.reactions
            if not r.id.startswith("EX_")
            and r.id != biomass_id
            and bool((r.gene_reaction_rule or "").strip())
        ],
        key=lambda r: r.id,
    )
    if max_tests_per_arm is not None and max_tests_per_arm > 0:
        evidenced_rxns = evidenced_rxns[:max_tests_per_arm]

    fva_tasks: List[Tuple[int, str, bool]] = [
        (idx, r.id, r.lower_bound < 0.0) for idx, r in enumerate(evidenced_rxns)
    ]

    global _WORKER_FILLER, _WORKER_UNION_MED, _WORKER_REPLAY_FRAC, _WORKER_FLUX_TOL
    _WORKER_FILLER = filler
    _WORKER_UNION_MED = dict(union_med)
    _WORKER_REPLAY_FRAC = replay_fraction
    _WORKER_FLUX_TOL = flux_tol

    use_fork_pool = (
        workers > 1
        and "fork" in mp.get_all_start_methods()
        and (len(fva_tasks) + len(a1_descriptors)) > 1
    )

    fva_lps = 0
    fva_cpu = 0.0
    fva_wall = 0.0
    blocked_a2_rxns: List[Tuple[str, bool]] = []

    if use_fork_pool and len(fva_tasks) > 1:
        ctx = mp.get_context("fork")
        n_pool = min(workers, len(fva_tasks))
        with ctx.Pool(processes=n_pool, initializer=_init_single_thread_worker) as pool:
            fva_out = _ordered_with_progress(pool, _worker_check_blocked, fva_tasks, "flux-variability checks")
        fva_out.sort(key=lambda item: item[0])
    else:
        fva_out = [_worker_check_blocked(t) for t in fva_tasks]

    for _, rxn_id, is_rev, unblocked, lps, cpu_s, wall_s in fva_out:
        fva_lps += lps
        fva_cpu += cpu_s
        fva_wall += wall_s
        if not unblocked:
            blocked_a2_rxns.append((rxn_id, is_rev))

    a2_descriptors: List[Tuple[str, str, str, Optional[PanelEntry]]] = [
        (f"A2_unblock_{rxn_id}", "A2", rxn_id, None)
        for rxn_id, _ in blocked_a2_rxns
    ]

    # Prepare A3 dead-end descriptors
    dead_ends = identify_dead_end_targets(masked_model)
    if max_tests_per_arm is not None and max_tests_per_arm > 0:
        dead_ends = dead_ends[:max_tests_per_arm]
    a3_descriptors: List[Tuple[str, str, str, Optional[PanelEntry]]] = []
    for met_id, kind in dead_ends:
        target_label = f"{met_id}:{kind}"
        test_id = f"A3_deadend_{met_id}_{kind}"
        a3_descriptors.append((test_id, "A3", target_label, None))

    # Build unified task list across A1, A2, A3 in deterministic order
    solve_tasks: List[
        Tuple[
            int,
            str,
            str,
            str,
            Dict[str, float],
            bool,
            bool,
            Optional[Tuple[str, Dict[str, float], float, float]],
        ]
    ] = []
    task_idx = 0
    for entry in eff_panel:
        test_id = f"A1_{entry.panel_id}"
        sole_med = build_sole_source_medium(default_medium, entry, metabolites)
        solve_tasks.append(
            (task_idx, test_id, "A1", biomass_id, sole_med, True, False, None)
        )
        task_idx += 1

    for rxn_id, is_rev in blocked_a2_rxns:
        test_id = f"A2_unblock_{rxn_id}"
        solve_tasks.append(
            (task_idx, test_id, "A2", rxn_id, union_med, False, is_rev, None)
        )
        task_idx += 1

    for met_id, kind in dead_ends:
        test_id = f"A3_deadend_{met_id}_{kind}"
        if kind == "unproduced":
            temp_spec = (f"DM_{met_id}", {met_id: -1.0}, 0.0, 1000.0)
        else:
            temp_spec = (f"SRC_{met_id}", {met_id: 1.0}, 0.0, 1000.0)
        solve_tasks.append(
            (task_idx, test_id, "A3", temp_spec[0], union_med, False, False, temp_spec)
        )
        task_idx += 1

    if use_fork_pool and len(solve_tasks) > 1:
        ctx = mp.get_context("fork")
        n_pool = min(workers, len(solve_tasks))
        with ctx.Pool(processes=n_pool, initializer=_init_single_thread_worker) as pool:
            solved_out = _ordered_with_progress(pool, _worker_solve_test, solve_tasks, "gap-fill tests")
        solved_out.sort(key=lambda item: item[0])
    else:
        solved_out = [_worker_solve_test(t) for t in solve_tasks]

    _WORKER_FILLER = None

    a1_results: List[GapFillTestResult] = []
    a2_results: List[GapFillTestResult] = []
    a3_results: List[GapFillTestResult] = []
    for _, arm, res in solved_out:
        if arm == "A1":
            a1_results.append(res)
        elif arm == "A2":
            a2_results.append(res)
        elif arm == "A3":
            a3_results.append(res)

    arm_a1 = _summarize_arm("A1", a1_results)
    arm_a2 = _summarize_arm(
        "A2",
        a2_results,
        extra_lps=fva_lps,
        extra_cpu=fva_cpu,
        extra_wall=fva_wall,
    )
    arm_a3 = _summarize_arm("A3", a3_results)

    # --- Combined A1 U A2 U A3 ---
    all_a123_results = a1_results + a2_results + a3_results
    arm_a123 = _summarize_arm(
        "A1_A2_A3",
        all_a123_results,
        extra_lps=fva_lps,
        extra_cpu=fva_cpu,
        extra_wall=fva_wall,
    )

    # --- Arm A4: rule-ranked subsets at B = 10% and B = 25% ---
    all_descriptors = a1_descriptors + a2_descriptors + a3_descriptors
    ranked_entries = rank_a4_tests(
        test_descriptors=all_descriptors,
        draft_model=masked_model,
        universal_reactions=universal_reactions,
        rxn_weights=rxn_weights,
        rxn_to_candidates=rxn_to_candidates,
        transporters=transporters,
        kegg_modules=kegg_modules,
        union_medium=union_med,
    )
    result_by_test_id: Dict[str, GapFillTestResult] = {
        r.test_id: r for r in all_a123_results
    }
    n_total_tests = len(ranked_entries)
    arm_a4_10 = _build_a4_budget_result(
        "A4@10%", 0.10, ranked_entries, result_by_test_id, n_total_tests
    )
    arm_a4_25 = _build_a4_budget_result(
        "A4@25%", 0.25, ranked_entries, result_by_test_id, n_total_tests
    )

    return {
        "A0": arm_a0,
        "A1": arm_a1,
        "A2": arm_a2,
        "A3": arm_a3,
        "A1_A2_A3": arm_a123,
        "A4@10%": arm_a4_10,
        "A4@25%": arm_a4_25,
    }


def run_all_arms_unmasked(
    unmasked_draft_model: cobra.Model,
    universal_reactions: Mapping[str, SeedReaction],
    metabolites: Mapping[str, SeedMetabolite],
    unmasked_candidates: Sequence[CandidateClaim],
    ec_to_reactions: Mapping[str, Sequence[str]],
    panel_entries: Sequence[PanelEntry],
    default_medium: Mapping[str, float],
    transporters: Sequence[TransporterHit] = (),
    kegg_modules: Sequence[KeggModuleDef] = (),
    replay_fraction: float = 0.05,
    flux_tol: float = DEFAULT_FLUX_TOL,
    workers: int = 1,
    max_tests_per_arm: Optional[int] = None,
) -> Dict[str, ArmRunResult]:
    """Run every arm (A0-A4) on the unmasked draft model with unmasked candidates (§8 control)."""
    return run_all_arms_for_case(
        masked_model=unmasked_draft_model,
        audit_report=None,
        universal_reactions=universal_reactions,
        metabolites=metabolites,
        candidates=unmasked_candidates,
        ec_to_reactions=ec_to_reactions,
        panel_entries=panel_entries,
        default_medium=default_medium,
        transporters=transporters,
        kegg_modules=kegg_modules,
        replay_fraction=replay_fraction,
        flux_tol=flux_tol,
        allow_unmasked=True,
        workers=workers,
        max_tests_per_arm=max_tests_per_arm,
    )


def _summarize_arm(
    arm_name: str,
    results: Sequence[GapFillTestResult],
    extra_lps: int = 0,
    extra_cpu: float = 0.0,
    extra_wall: float = 0.0,
    ranked_entries: Tuple[RankedTestEntry, ...] = (),
    tests_considered: Optional[int] = None,
) -> ArmRunResult:
    added: Set[str] = set()
    replays: List[SolverReplayRecord] = []
    n_lps = extra_lps
    n_gfs = 0
    cpu_s = extra_cpu
    wall_s = extra_wall

    for res in results:
        added.update(res.added_reactions)
        n_lps += res.n_lps
        n_gfs += res.n_gapfills
        cpu_s += res.cpu_seconds
        wall_s += res.wall_seconds
        if res.replay is not None:
            replays.append(res.replay)

    return ArmRunResult(
        arm=arm_name,
        tests_considered=tests_considered if tests_considered is not None else len(results),
        tests_executed=len(results),
        n_lps=n_lps,
        n_gapfills=n_gfs,
        cpu_seconds=cpu_s,
        wall_seconds=wall_s,
        peak_rss_kb=peak_rss_kb(),
        added_reactions=tuple(sorted(added)),
        test_results=tuple(results),
        ranked_entries=ranked_entries,
        replays=tuple(replays),
    )


def _build_a4_budget_result(
    arm_name: str,
    budget_fraction: float,
    ranked_entries: Sequence[RankedTestEntry],
    result_by_test_id: Mapping[str, GapFillTestResult],
    n_total_tests: int,
) -> ArmRunResult:
    if n_total_tests == 0:
        return _summarize_arm(arm_name, (), ranked_entries=(), tests_considered=0)
    k_budget = max(1, int(math.ceil(budget_fraction * n_total_tests)))
    selected_entries = tuple(ranked_entries[:k_budget])
    selected_results = [
        result_by_test_id[e.test_id]
        for e in selected_entries
        if e.test_id in result_by_test_id
    ]
    return _summarize_arm(
        arm_name,
        selected_results,
        ranked_entries=selected_entries,
        tests_considered=n_total_tests,
    )
