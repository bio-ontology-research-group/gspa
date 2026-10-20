#!/usr/bin/env python3
"""Threshold draft, Grow-Support repair and network profile on the ModelSEED universal model.

This is GSPA's own implementation of the comparator used in the 108-genome invariant audit
(paper-reproducibility/pipelines/P07_invariant_audit): the tau-threshold draft of a per-protein
EC prediction matrix, the structurally excluded reactions, METEOR's default MILP medium, the
Grow-Support LP repair (minimal flux through non-draft reactions such that biomass >= gamma_min),
the maximal biomass of a reaction set, and the structural profile of a network (dead-end
metabolites, mass-unbalanced reactions, FBA growth) in the sense of the MEMOTE checks.

Inputs are plain data files (numpy/scipy pickles, a JSON description of the universal model,
TSV/text EC tables); no third-party pipeline code is imported.  LPs are solved with scipy's
HiGHS interface.  The smoke test in P07 compares every quantity with the archived values.

Definitions (indices refer to the reaction order ``allrxns`` of ``universal_metadata.json``,
which is the column order of ``fba_matrices_v6RE.pkl``):

* excluded reactions: the other Gram's biomass reaction; reactions of metabolites carried by
  one reaction only, or by two of which one is a sink (``SNK_``); reactions of metabolites whose
  stoichiometric coefficients all have one sign and which touch no reversible reaction.
* default medium: exchanges of the named medium plus a fixed set of inorganic and amino-acid
  uptakes are open in both directions (bounds -100..100); every other exchange is blocked for
  uptake (0..1000).  Because the matrix stores exchanges with reversed sign, v > 0 on an
  exchange is uptake.
* threshold draft: a reaction is in the draft when any mapped EC of any protein scores >= tau.
* core = draft minus excluded; Grow-Support repair = core plus the support of the LP
  ``min sum_{j not in core} |v_j|  s.t.  S v = 0,  v_biomass >= gamma_min,  lb <= v <= ub``
  over the non-excluded reactions, support = |v| > 1e-9.
* network profile: the sub-network of the selected reactions plus the universal model's
  boundary reactions of their metabolites plus the Gram biomass.  Dead ends are metabolites
  whose every reaction is irreversible, not an exchange, and produces the metabolite; a reaction
  is mass-unbalanced when a metabolite has no parsable formula or the element sums differ from
  zero (exchanges included, as in the archived profiles).
"""
from __future__ import annotations

import json
import pickle
import re
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import csr_matrix, csc_matrix, hstack

SUPPORT_EPS = 1e-9
GAMMA_MIN = 0.1
PROFILE_BOUND = 100.0
PROFILE_UPTAKE_LB = -10.0
_HIGHS_OPTS = {"presolve": True, "primal_feasibility_tolerance": 1e-9, "dual_feasibility_tolerance": 1e-9}

# Uptakes METEOR's apply_media always opens in addition to the named medium (amino acids,
# nucleobases, vitamins, inorganic nutrients).  Part of the comparator's definition.
MINIMAL_UPTAKE = {
    "EX_cpd00035_e", "EX_cpd00051_e", "EX_cpd00132_e", "EX_cpd00041_e", "EX_cpd00084_e",
    "EX_cpd00053_e", "EX_cpd00023_e", "EX_cpd00033_e", "EX_cpd00119_e", "EX_cpd00322_e",
    "EX_cpd00107_e", "EX_cpd00039_e", "EX_cpd00060_e", "EX_cpd00066_e", "EX_cpd00129_e",
    "EX_cpd00054_e", "EX_cpd00161_e", "EX_cpd00065_e", "EX_cpd00069_e", "EX_cpd00156_e",
    "EX_cpd00027_e", "EX_cpd00149_e", "EX_cpd00030_e", "EX_cpd00254_e", "EX_cpd00971_e",
    "EX_cpd00063_e", "EX_cpd10515_e", "EX_cpd00205_e", "EX_cpd00099_e",
    "EX_cpd00009_e", "EX_cpd00001_e", "EX_cpd00067_e", "EX_cpd00007_e", "EX_cpd00013_e",
    "EX_cpd00048_e", "EX_cpd00011_e",
}


# ----------------------------------------------------------------------------- data loading
class UniversalData:
    """The universal model as plain arrays plus the EC and medium tables."""

    def __init__(self, data_dir: Path | str, univ_meta_path: Path | str | None = None):
        self.data_dir = Path(data_dir)
        meta_path = Path(univ_meta_path) if univ_meta_path else self.data_dir / "universal_metadata.json"
        if not meta_path.exists():
            repo = Path(__file__).resolve().parents[1] / "benchmark" / "manifests" / "reference_data" / "universal_metadata.json"
            if repo.exists():
                meta_path = repo
        with open(meta_path, encoding="utf-8") as f:
            self.univ: dict[str, Any] = json.load(f)
        self.allrxns: list[str] = list(self.univ["allrxns"])
        self.allmet: list[str] = list(self.univ["allmet"])
        self.rxn_meta: dict[str, dict] = self.univ["rxn_meta"]
        self.met_meta: dict[str, dict] = self.univ["met_meta"]
        self.rxn_index = {r: i for i, r in enumerate(self.allrxns)}
        self.met_index = {m: i for i, m in enumerate(self.allmet)}

        with open(self.data_dir / "fba_matrices_v6RE.pkl", "rb") as f:
            fba = pickle.load(f)
        self.S: csr_matrix = csr_matrix(fba["S"])
        self.lb0 = np.asarray(fba["lb"], dtype=np.float64)
        self.ub0 = np.asarray(fba["ub"], dtype=np.float64)
        if self.S.shape != (len(self.allmet), len(self.allrxns)):
            raise ValueError(f"S is {self.S.shape}, metadata has {len(self.allmet)} x {len(self.allrxns)}")

        self.tight: dict[str, tuple[np.ndarray, np.ndarray]] = {}
        for gram, tag in (("negative", "neg"), ("positive", "pos")):
            with open(self.data_dir / f"tight_bounds_v6_{tag}.pkl", "rb") as f:
                d = pickle.load(f)
            self.tight[gram] = (np.asarray(d["lb_tight"], np.float64), np.asarray(d["ub_tight"], np.float64))

        with open(self.data_dir / "medium.pkl", "rb") as f:
            self.media: dict[str, list[str]] = pickle.load(f)
        with open(self.data_dir / "seedr2ec.pkl", "rb") as f:
            raw = pickle.load(f)
        self.seedr2ec: dict[str, list[str]] = {k: v for k, v in raw.items() if v}
        self.ancestors_ec: list[str] = (self.data_dir / "all_ancestors.txt").read_text().splitlines()
        with open(self.data_dir / "enzymeobsolete.pkl", "rb") as f:
            self.ec_obsolete: dict[str, str] = pickle.load(f)

    # -- per-Gram medium
    def biomass_id(self, gram: str) -> str:
        return "biomass_GmPos" if gram == "positive" else "biomass_GmNeg"

    def milp_bounds(self, gram: str, medium: str = "default") -> tuple[np.ndarray, np.ndarray, set[int], set[int]]:
        """METEOR's MILP bounds for a Gram: universal bounds tightened by FVA, then the medium.
        Returns lb, ub, excluded reaction indices (computed before the medium, as the comparator
        does) and the open exchange indices."""
        lt, ut = self.tight[gram]
        lb = np.maximum(self.lb0, lt)
        ub = np.minimum(self.ub0, ut)
        excluded = find_excluded_reactions(self.S, lb, ub, self.allrxns, self.biomass_id(gram))
        lb, ub, open_ex = apply_default_medium(self.allrxns, lb, ub, self.media, medium)
        return lb, ub, excluded, open_ex

    def rxn_ec_mask(self) -> list[np.ndarray]:
        """For every reaction, the indices (into ancestors_ec) of its mapped ECs."""
        ec_to_idx = {ec: i for i, ec in enumerate(self.ancestors_ec)}
        out = []
        for r in self.allrxns:
            base = r[:-2] if r.endswith(("_c", "_e", "_p")) else r
            out.append(np.array(sorted({ec_to_idx[e] for e in self.seedr2ec.get(base, []) if e in ec_to_idx}), dtype=np.int64))
        return out


def find_excluded_reactions(S, lb, ub, allrxns: list[str], biomass_id: str) -> set[int]:
    Sr = csr_matrix(S)
    excl: set[int] = set()
    if "GmNeg" in biomass_id and "biomass_GmPos" in allrxns:
        excl.add(allrxns.index("biomass_GmPos"))
    elif "GmPos" in biomass_id and "biomass_GmNeg" in allrxns:
        excl.add(allrxns.index("biomass_GmNeg"))
    indptr, indices, data = Sr.indptr, Sr.indices, Sr.data
    for i in range(Sr.shape[0]):
        cols = indices[indptr[i]:indptr[i + 1]]
        vals = data[indptr[i]:indptr[i + 1]]
        n = len(cols)
        if n == 0:
            continue
        if n == 1:
            excl.add(int(cols[0]))
        elif n == 2 and any(allrxns[k].startswith("SNK_") for k in cols):
            excl.add(int(cols[0]))  # the comparator adds the first carrier, whichever it is
        if np.all(vals > 0) or np.all(vals < 0):
            if not any(lb[k] < 0 and ub[k] > 0 for k in cols):
                excl.update(int(k) for k in cols)
    return excl


def apply_default_medium(allrxns: list[str], lb, ub, media: dict[str, list[str]], medium: str = "default"):
    lb = np.array(lb, dtype=np.float64, copy=True)
    ub = np.array(ub, dtype=np.float64, copy=True)
    named = media.get(medium, [medium] if medium not in media else [])
    allowed = {"EX_" + c + "_e" for c in named} | MINIMAL_UPTAKE
    open_ex: set[int] = set()
    for i, r in enumerate(allrxns):
        if not r.startswith("EX_"):
            continue
        if r in allowed:
            lb[i], ub[i] = -100.0, 100.0
            open_ex.add(i)
        else:
            lb[i], ub[i] = 0.0, 1000.0
    return lb, ub, open_ex


# ----------------------------------------------------------------------------- predictions
def current_ec(ec: str, obsolete: dict[str, str], max_hops: int = 5) -> str:
    cur, seen, hops = ec, set(), 0
    while cur in obsolete and cur not in seen and hops < max_hops:
        seen.add(cur)
        cur = obsolete[cur]
        hops += 1
    return cur


def load_prediction_matrix(path: Path | str, ancestors_ec: list[str], obsolete: dict[str, str]):
    """Per-protein EC score matrix (proteins x ancestors_ec, float32) from a pickled DataFrame.
    A matrix stored EC-by-protein is transposed; ``EC:`` prefixes are stripped; a column whose
    EC number is obsolete is renamed to its current number (taking the maximum when both exist);
    ECs absent from the file score 0."""
    import pandas as pd

    df = pd.read_pickle(path)
    if str(df.index[0]).count(".") == 3:
        df = df.T
    if "EC:" in str(df.columns[0]):
        df.columns = [c.split("EC:")[-1] if "EC:" in c else c for c in df.columns]
    anc = set(ancestors_ec)
    rename = {}
    for c in df.columns:
        if c not in anc and c in obsolete:
            cur = current_ec(c, obsolete)
            if cur in anc and cur != c:
                rename[c] = cur
    for old, new in rename.items():
        if new in df.columns:
            df[new] = df[[old, new]].max(axis=1)
            df = df.drop(columns=[old])
        else:
            df = df.rename(columns={old: new})
    return df.reindex(columns=ancestors_ec, fill_value=0.0).astype(np.float32)


def threshold_draft(pred_values: np.ndarray, rxn_ecs: list[np.ndarray], tau: float):
    """Boolean draft mask and per-reaction maximal score over the mapped ECs."""
    hit = (pred_values >= tau).any(axis=0)
    colmax = pred_values.max(axis=0) if pred_values.shape[0] else np.zeros(pred_values.shape[1], np.float32)
    n = len(rxn_ecs)
    draft = np.zeros(n, dtype=bool)
    maxscore = np.zeros(n, dtype=np.float32)
    for j, ecs in enumerate(rxn_ecs):
        if len(ecs):
            draft[j] = bool(hit[ecs].any())
            maxscore[j] = float(colmax[ecs].max())
    return draft, maxscore


# ----------------------------------------------------------------------------- LPs
def _split_bounds(lb, ub, allowed):
    """Bounds of the split variables p (forward) and q (backward); a disallowed reaction is fixed at 0."""
    p_ub = np.where(allowed, np.maximum(0.0, ub), 0.0)
    q_ub = np.where(allowed, np.maximum(0.0, -lb), 0.0)
    return p_ub, q_ub


def _solve(S, lb, ub, allowed, c_p, c_q, obj_idx, gamma_min=None, maximize=False):
    Sc = csc_matrix(S)
    n = Sc.shape[1]
    A_eq = hstack([Sc, -Sc], format="csc")
    b_eq = np.zeros(Sc.shape[0])
    p_ub, q_ub = _split_bounds(lb, ub, allowed)
    bounds = np.column_stack([np.zeros(2 * n), np.concatenate([p_ub, q_ub])])
    c = np.concatenate([c_p, c_q]).astype(np.float64)
    if maximize:
        c = -c
    A_ub = b_ub = None
    if gamma_min is not None:
        row = np.zeros(2 * n)
        row[obj_idx] = -1.0
        row[n + obj_idx] = 1.0
        A_ub = csr_matrix(row.reshape(1, -1))
        b_ub = np.array([-gamma_min])
    res = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, bounds=bounds, method="highs", options=_HIGHS_OPTS)
    if res.status != 0:
        return None, None
    v = res.x[:n] - res.x[n:]
    return v, float(res.fun if not maximize else -res.fun)


def maxbio_active(S, lb, ub, obj_idx: int, active: np.ndarray) -> float:
    """Maximal biomass flux using only the ``active`` reactions; 0 when the set cannot grow."""
    n = S.shape[1]
    c = np.zeros(n)
    c[obj_idx] = 1.0
    v, val = _solve(S, lb, ub, active, c, -c, obj_idx, maximize=True)
    return 0.0 if v is None else float(v[obj_idx])


def grow_support(S, lb, ub, obj_idx: int, avail: np.ndarray, gamma_min: float = GAMMA_MIN, core: np.ndarray | None = None):
    """Support of a minimal-non-core-flux growing solution within ``avail``; (None, None, None)
    when infeasible.  Returns (support mask, flux vector, objective value)."""
    n = S.shape[1]
    cost = np.ones(n) if core is None else (~core).astype(np.float64)
    v, val = _solve(S, lb, ub, avail, cost, cost, obj_idx, gamma_min=gamma_min)
    if v is None:
        return None, None, None
    return np.abs(v) > SUPPORT_EPS, v, val


# ----------------------------------------------------------------------------- network profile
_ELEMENT_RE = re.compile(r"([A-Z][a-z]?)([0-9.]+[0-9.]?|(?=[A-Z])?)")


def formula_elements(formula: str | None) -> dict[str, float] | None:
    """Element counts of a formula in the sense of cobra's ``Metabolite.elements``: ``None``
    for an unparsable formula (parentheses), ``{}`` for a missing one."""
    if formula is None:
        return {}
    f = str(formula)
    if "*" in f:
        f = f.replace("*", "")
    if "(" in f or ")" in f:
        return None
    out: dict[str, float] = {}
    for element, count in _ELEMENT_RE.findall(f):
        if count == "":
            n: float = 1.0
        else:
            try:
                n = float(count)
            except ValueError:
                return None
            if n.is_integer():
                n = int(n)
        out[element] = out.get(element, 0) + n
    return out


class NetworkProfiler:
    """Structural profile of a reaction subset of the universal model."""

    def __init__(self, U: UniversalData):
        self.U = U
        self.elements = {m: formula_elements(U.met_meta.get(m, {}).get("formula")) for m in U.allmet}
        # boundary reactions in cobra's sense: exactly one metabolite (EX_, DM_, SNK_); the
        # metadata's "boundary" flag is also set on the biomass reactions and is not used here
        self.boundary_rxns = [r for r in U.allrxns if len(U.rxn_meta[r].get("stoich", {})) == 1]
        self.balanced = {r: self._is_mass_balanced(r) for r in U.allrxns}

    def _is_mass_balanced(self, rid: str) -> bool:
        balance: dict[str, float] = {}
        for m, coef in self.U.rxn_meta[rid].get("stoich", {}).items():
            el = self.elements.get(m)
            if el is None or len(el) == 0:
                return False
            for e, amount in el.items():
                balance[e] = balance.get(e, 0) + coef * amount
        return all(a == 0 for a in balance.values())

    def submodel(self, keep: np.ndarray, biomass_id: str) -> list[str]:
        """Selected reactions, plus boundary reactions of their metabolites, plus the biomass."""
        rxns = [r for r, k in zip(self.U.allrxns, keep) if k]
        present = set(rxns)
        mets = {m for r in rxns for m in self.U.rxn_meta[r].get("stoich", {})}
        for r in self.boundary_rxns:
            if r not in present:
                st = self.U.rxn_meta[r].get("stoich", {})
                if st and next(iter(st)) in mets:
                    rxns.append(r)
                    present.add(r)
        if biomass_id not in present:
            rxns.append(biomass_id)
        return rxns

    def deadends(self, rxns: list[str], lb: dict[str, float], ub: dict[str, float]) -> int:
        by_met: dict[str, list[tuple[str, float]]] = {}
        for r in rxns:
            for m, coef in self.U.rxn_meta[r].get("stoich", {}).items():
                by_met.setdefault(m, []).append((r, coef))
        n = 0
        for m, lst in by_met.items():
            if all((not (lb[r] < 0 < ub[r])) and (not r.startswith("EX_")) and coef > 0 for r, coef in lst):
                n += 1
        return n

    def _profile_bounds(self, rxns: list[str], lb_milp: np.ndarray, ub_milp: np.ndarray, bio: str):
        """Bounds of the archived profile protocol, after the exchange-uptake discovery step."""
        U = self.U
        lb: dict[str, float] = {}
        ub: dict[str, float] = {}
        for r in rxns:
            meta = U.rxn_meta[r]
            if r.startswith(("EX_", "DM_", "SK_")) or "biomass" in r.lower():
                lb[r], ub[r] = float(meta["lb"]), float(meta["ub"])
            else:
                i = U.rxn_index.get(r)
                lb[r] = max(-PROFILE_BOUND, float(lb_milp[i])) if i is not None else -PROFILE_BOUND
                ub[r] = min(PROFILE_BOUND, float(ub_milp[i])) if i is not None else PROFILE_BOUND
        exs = [r for r in rxns if r.startswith("EX_") and r != "EX_biomass"]
        for r in exs:
            lb[r], ub[r] = -1000.0, 1000.0
        _, flux = self._fba(rxns, lb, ub, bio)
        uptake = {r for r in exs if flux is not None and flux.get(r, 0.0) < -1e-6}
        for r in exs:
            lb[r] = PROFILE_UPTAKE_LB if r in uptake else 0.0
            ub[r] = 1000.0
        return lb, ub

    def deadends_shared_objects(self, keep: np.ndarray, gram: str, lb_milp: np.ndarray, ub_milp: np.ndarray) -> int:
        """Dead-end count as the archived profiles computed it: for every metabolite of the network,
        all reactions of the *universal* model touching it are consulted (network reactions with the
        profile bounds, the others with their universal bounds).  Kept only to document the archive."""
        U = self.U
        bio = U.biomass_id(gram)
        rxns = self.submodel(keep, bio)
        lb, ub = self._profile_bounds(rxns, lb_milp, ub_milp, bio)
        if not hasattr(self, "_umet"):
            umet: dict[str, list[tuple[str, float]]] = {}
            for r in U.allrxns:
                for m, coef in U.rxn_meta[r].get("stoich", {}).items():
                    umet.setdefault(m, []).append((r, coef))
            self._umet = umet
        mets = {m for r in rxns for m in U.rxn_meta[r].get("stoich", {})}
        n = 0
        for m in mets:
            lst = self._umet.get(m, [])
            ok = bool(lst)
            for r, coef in lst:
                l = lb.get(r, float(U.rxn_meta[r]["lb"]))
                u = ub.get(r, float(U.rxn_meta[r]["ub"]))
                if (l < 0 < u) or r.startswith("EX_") or coef <= 0:
                    ok = False
                    break
            n += ok
        return n

    def profile(self, keep: np.ndarray, gram: str, lb_milp: np.ndarray, ub_milp: np.ndarray) -> dict[str, Any]:
        """Dead ends, mass-unbalanced reactions and FBA growth of the sub-network built from ``keep``,
        with the bound protocol of the archived profiles: internal reactions clipped to +-100 of
        their MILP bounds, exchanges first unbounded to find the uptakes an optimum uses, then
        those uptakes limited to 10 and every other exchange closed for uptake."""
        U = self.U
        bio = U.biomass_id(gram)
        rxns = self.submodel(keep, bio)
        lb, ub = self._profile_bounds(rxns, lb_milp, ub_milp, bio)
        growth, _ = self._fba(rxns, lb, ub, bio)
        de = self.deadends(rxns, lb, ub)
        mi = sum(1 for r in rxns if not self.balanced[r])
        internal = [r for r in rxns if len(U.rxn_meta[r].get("stoich", {})) > 1 and not r.startswith("biomass")]
        mi_int = sum(1 for r in internal if not self.balanced[r])
        return {"n_selected": int(keep.sum()), "n_rxn": len(rxns), "fba_growth": round(float(growth), 4),
                "deadends": de, "mass_imbal": mi, "mi_frac": round(mi / max(1, len(rxns)), 4),
                "n_internal": len(internal), "mass_imbal_internal": mi_int}

    def _fba(self, rxns: list[str], lb: dict[str, float], ub: dict[str, float], bio: str):
        """Maximise the biomass reaction over the sub-network (original, non-reversed stoichiometry)."""
        mets: dict[str, int] = {}
        rows, cols, vals = [], [], []
        for j, r in enumerate(rxns):
            for m, coef in self.U.rxn_meta[r].get("stoich", {}).items():
                rows.append(mets.setdefault(m, len(mets)))
                cols.append(j)
                vals.append(float(coef))
        S = csc_matrix((vals, (rows, cols)), shape=(len(mets), len(rxns)))
        c = np.zeros(len(rxns))
        c[rxns.index(bio)] = -1.0
        bounds = [(lb[r], ub[r]) for r in rxns]
        res = linprog(c, A_eq=S, b_eq=np.zeros(len(mets)), bounds=bounds, method="highs", options=_HIGHS_OPTS)
        if res.status != 0:
            return 0.0, None
        return -res.fun, dict(zip(rxns, res.x))
