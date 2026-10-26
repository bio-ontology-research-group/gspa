#!/usr/bin/env python3
"""
Community cross-feeding on separate organism SBML models (LP, budgeted, auditable).

Each organism keeps its own stoichiometric block: intracellular and organism-extracellular
species are namespaced ``<organism>::<species>`` and are never merged.  The only coupling is
through *explicit ports*: an organism exchange reaction is mapped to a shared pool metabolite
(``pool::<P>``), where

    sum_o  alpha_o * scale * secretion_{o,P}  +  supply_P  -  removal_P  =  0

``alpha_o`` is the fixed biomass fraction of organism o (gDW_o per gDW community), fluxes are
mmol/gDW_o/h, ``supply_P`` / ``removal_P`` are the only environment connections (mmol per gDW
community per h, declared and finite; absent means closed).  Organism exchange reactions are
therefore ports to the pool, never free access to the environment.

Growth is balanced: every organism grows at the same rate mu (v_biomass_o = mu) and the
community objective is max mu, so a solution in which one organism does not grow while another
does is not admissible.  Matched monocultures (alpha = 1, same pool, same ports, same supply and
removal) are solved with the same builder.

Everything reported is a model-conditional prediction.  A model-dependent rescue says the
network can grow if the declared transfer exists; it assigns no enzyme or gene, does not show
that the transporter exists, and carries no probability.  Solver failure / infeasibility are
reported as such and are never read as evidence about a transfer.  Mechanistic fixtures in the
tests are not empirical validation; measured co-culture data remain outstanding.
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass, field
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
from scipy.sparse import coo_matrix, csc_matrix

sys.path.insert(0, str(Path(__file__).resolve().parent))
from metabolic_experiments import (  # noqa: E402
    BudgetedSolver, BudgetExhausted, ManifestError, Model, canonical_hash, close_uptake,
    load_sbml_model, resolve_objective,
)

SCHEMA = "gspa.community_metabolism/1"
FLUX_UNITS = "mmol/gDW/h"
BIOMASS_UNITS = "1/h"
SUPPORT_VALUES = ("curated", "predicted", "unknown")

DEFAULT_CONFIG: dict[str, Any] = {
    "max_solver_calls": 40,
    "growth_threshold": 1e-6,
    "flux_tol": 1e-7,
    "residual_tol": 1e-6,
    "relative_tol": 0.05,
    "growth_fraction": 0.95,                         # mu floor for FVA of transfers
    "open_ports_with_support": ["curated", "predicted"],   # ports with other support are closed
    "auto_port_closure": True,
    "max_fva_candidates": 4,
    "default_bound": 1000.0,
}

BUDGET_POLICY = {
    "order": "monocultures (1 LP each) -> community (1 LP) -> user experiments (1 LP each) -> "
             "automatic producer-port closures (1 LP per candidate transfer x producer) -> "
             "FVA of candidate transfers (2 LPs each, at most config.max_fva_candidates)",
    "skip_rule": "an experiment is skipped (skipped_budget) when remaining calls < its cost",
    "cap": "config.max_solver_calls is a hard cap over all LPs",
}

DISCLAIMER = (
    "Model-conditional predictions on supplied models, ports and environment. Not observed "
    "co-culture phenotypes or exchanged metabolites; no enzyme or gene is assigned; transporter "
    "existence is only as supported by the declared port evidence; no calibrated probabilities."
)


@dataclass
class Organism:
    id: str
    model: Model
    sha256: str
    biomass_col: int
    maintenance: dict[str, Any] | str
    ports: dict[str, dict[str, Any]]
    closed: set[str]
    knockouts: set[str]
    ext: set[str]
    ports_closed_by_policy: list[str] = field(default_factory=list)


@dataclass
class CommunityLP:
    S: csc_matrix
    lb: np.ndarray
    ub: np.ndarray
    mu_col: int
    row_names: list[str]
    col_names: list[str]
    bio_cols: dict[str, int]
    port_cols: list[dict[str, Any]]
    env_cols: dict[tuple[str, str], int]
    pool_rows: dict[str, int]
    organisms: list[str]


# --------------------------------------------------------------------------- loading / validation


def _finite_nonneg(d: dict[str, Any], what: str) -> dict[str, float]:
    out = {}
    for k, v in (d or {}).items():
        if not isinstance(v, (int, float)) or not np.isfinite(v) or v < 0:
            raise ManifestError(f"{what}[{k}] must be a finite number >= 0")
        out[k] = float(v)
    return out


def validate_pool(manifest: dict[str, Any]) -> tuple[dict[str, Any], dict[str, float], dict[str, float]]:
    pool = manifest.get("pool", {})
    mets = pool.get("metabolites")
    if not mets:
        raise ManifestError("pool.metabolites is required")
    supply = _finite_nonneg(pool.get("supply", {}), "pool.supply")
    removal = _finite_nonneg(pool.get("removal", {}), "pool.removal")
    for k in list(supply) + list(removal):
        if k not in mets:
            raise ManifestError(f"pool supply/removal names unknown pool metabolite: {k}")
    return mets, supply, removal


def build_organism(oid: str, model: Model, spec: dict[str, Any], pool_mets: dict[str, Any],
                   cfg: dict[str, Any]) -> Organism:
    """Validate one organism spec against an already parsed model (hidden open boundaries rejected).

    Reusable so that many organisms (e.g. gene-deletion variants) can share one parsed model."""
    bio_id, _ = resolve_objective(model, {"objectives": {"biomass": spec.get("biomass")}})
    maint = spec.get("maintenance")
    if maint != "none" and not (isinstance(maint, dict) and {"reaction", "lower_bound"} <= set(maint)):
        raise ManifestError(f"{oid}: maintenance must be {{reaction, lower_bound}} or 'none' (explicit)")
    if maint != "none" and maint["reaction"] not in model.rxn_index:
        raise ManifestError(f"{oid}: maintenance reaction not in model")
    ext_comp = set(spec.get("extracellular_compartments", ["e"]))
    ports = spec.get("ports", {})
    closed = set(spec.get("closed_ports", []))
    ext = {r for r in model.rxn_ids if model.is_single_metabolite(r)
           and model.met_compartment[model.exchange_metabolite(r)] in ext_comp}
    missing = sorted(ext - set(ports) - closed)
    if missing:
        raise ManifestError(f"{oid}: exchanges neither mapped to the pool nor closed (hidden open "
                            f"boundary rejected): {missing}")
    for rid in set(ports) | closed:
        if rid not in ext:
            raise ManifestError(f"{oid}: {rid} is not a single-metabolite exchange on an "
                                f"extracellular species")
    if set(ports) & closed:
        raise ManifestError(f"{oid}: ports both mapped and closed: {sorted(set(ports) & closed)}")
    by_policy = []
    for rid, p in ports.items():
        if p.get("pool") not in pool_mets:
            raise ManifestError(f"{oid}.{rid}: unknown pool metabolite {p.get('pool')}")
        if p.get("support", "unknown") not in SUPPORT_VALUES:
            raise ManifestError(f"{oid}.{rid}: support must be one of {SUPPORT_VALUES}")
        scale = p.get("scale", 1.0)
        if not isinstance(scale, (int, float)) or not np.isfinite(scale) or scale <= 0:
            raise ManifestError(f"{oid}.{rid}: scale must be a finite number > 0")
        met = model.exchange_metabolite(rid)
        if p.get("metabolite") and p["metabolite"] != met:
            raise ManifestError(f"{oid}.{rid}: declared metabolite {p['metabolite']} but the exchange "
                                f"acts on {met}")
        pf, mf = pool_mets[p["pool"]].get("formula"), model.met_formula.get(met)
        if pf and mf and pf != mf:
            raise ManifestError(f"{oid}.{rid}: formula {mf} of {met} differs from pool "
                                f"{p['pool']} formula {pf}")
        if p.get("support", "unknown") not in cfg["open_ports_with_support"]:
            by_policy.append(rid)
    ko = set(spec.get("knockouts", []))
    for rid in ko:
        if rid not in model.rxn_index:
            raise ManifestError(f"{oid}: knockout reaction not in model: {rid}")
    return Organism(oid, model, model.sha256, model.rxn_index[bio_id], maint, ports, closed, ko,
                    ext, sorted(by_policy))


def load_organisms(manifest: dict[str, Any], base_dir: Path, cfg: dict[str, Any]) -> list[Organism]:
    pool_mets, _, _ = validate_pool(manifest)
    specs = manifest.get("organisms", [])
    if not specs:
        raise ManifestError("manifest.organisms is empty")
    ids = [s.get("id") for s in specs]
    if len(set(ids)) != len(ids) or any(not i or "::" in i for i in ids):
        raise ManifestError("organism ids must be unique, non-empty and must not contain '::'")
    orgs = []
    for spec in specs:
        oid = spec["id"]
        path = Path(spec["model"])
        path = path if path.is_absolute() else base_dir / path
        model = load_sbml_model(path, cfg["default_bound"])
        if spec.get("sha256") and spec["sha256"] != model.sha256:
            raise ManifestError(f"{oid}: model sha256 mismatch (manifest {spec['sha256'][:12]}..., "
                                f"file {model.sha256[:12]}...)")
        if spec.get("flux_units", FLUX_UNITS) != FLUX_UNITS:
            raise ManifestError(f"{oid}: flux_units must be {FLUX_UNITS}")
        if spec.get("biomass_units") != BIOMASS_UNITS:
            raise ManifestError(f"{oid}: declare biomass_units '{BIOMASS_UNITS}' "
                                "(biomass reaction normalised to 1 gDW)")
        orgs.append(build_organism(oid, model, spec, pool_mets, cfg))
    return orgs


def community_ratios(manifest: dict[str, Any], orgs: list[Organism]) -> dict[str, float]:
    ratio = manifest.get("community", {}).get("biomass_ratio")
    if not ratio or set(ratio) != {o.id for o in orgs}:
        raise ManifestError("community.biomass_ratio must give a ratio for every organism")
    for k, v in ratio.items():
        if not isinstance(v, (int, float)) or not np.isfinite(v) or v <= 0:
            raise ManifestError(f"community.biomass_ratio[{k}] must be a finite number > 0")
    return {k: float(v) for k, v in ratio.items()}


# --------------------------------------------------------------------------- LP construction


def _direction_close(model: Model, lb: np.ndarray, ub: np.ndarray, rid: str, direction: str) -> None:
    col, coeff = model.exchange(rid)
    secretes_positive = coeff < 0
    if direction == "both":
        lb[col] = ub[col] = 0.0
    elif direction == "uptake":
        close_uptake(model, lb, ub, rid)
    elif direction == "secretion":
        if secretes_positive:
            ub[col] = min(ub[col], 0.0)
        else:
            lb[col] = max(lb[col], 0.0)
    else:
        raise ManifestError(f"unknown port direction: {direction}")


_COO_CACHE: dict[int, tuple[Model, tuple[np.ndarray, np.ndarray, np.ndarray]]] = {}


def _block_coo(m: Model) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Row, column, value triplets of a model's stoichiometry (cached; the matrix is never mutated)."""
    hit = _COO_CACHE.get(id(m))
    if hit is None or hit[0] is not m:
        coo = m.S.tocoo()
        hit = (m, (coo.row.astype(np.int64), coo.col.astype(np.int64), coo.data.astype(float)))
        _COO_CACHE[id(m)] = hit
    return hit[1]


def build_lp(orgs: list[Organism], alphas: dict[str, float], supply: dict[str, float],
             removal: dict[str, float], closures: list[tuple[str, str, str]] | None = None,
             knockouts: dict[str, set[str]] | None = None) -> CommunityLP:
    """Block-diagonal organism models + explicit port couplings + environment + common growth."""
    total = sum(alphas[o.id] for o in orgs)
    frac = {o.id: alphas[o.id] / total for o in orgs}
    rows: dict[str, int] = {}
    entries: list[tuple[int, int, float]] = []
    lbs: list[float] = []
    ubs: list[float] = []
    col_names: list[str] = []
    bio_cols: dict[str, int] = {}
    port_cols: list[dict[str, Any]] = []
    closures = closures or []
    knockouts = knockouts or {}

    def row(name: str) -> int:
        return rows.setdefault(name, len(rows))

    for o in orgs:
        m = o.model
        offset = len(col_names)
        row0 = len(rows)
        for r in m.met_ids:
            row(f"{o.id}::{r}")
        lb, ub = m.lb.copy(), m.ub.copy()
        # organism metabolite rows are contiguous from row0; the COO triplets are cached per model
        coo = _block_coo(m)
        entries.extend(zip((coo[0] + row0).tolist(), (coo[1] + offset).tolist(), coo[2].tolist()))
        # single-metabolite reactions that are not ports (sinks, demands, drains): no uptake
        single = np.flatnonzero(np.diff(m.S.indptr) == 1)       # same test as Model.is_single_metabolite
        for j in single:
            rid = m.rxn_ids[int(j)]
            if rid not in o.ext:
                close_uptake(m, lb, ub, rid)
        for rid in o.closed | set(o.ports_closed_by_policy):
            lb[m.rxn_index[rid]] = ub[m.rxn_index[rid]] = 0.0
        for oid, rid, direction in closures:
            if oid == o.id:
                _direction_close(m, lb, ub, rid, direction)
        if o.maintenance != "none":
            c = m.rxn_index[o.maintenance["reaction"]]
            lb[c] = max(lb[c], float(o.maintenance["lower_bound"]))
        for rid in o.knockouts | knockouts.get(o.id, set()):
            lb[m.rxn_index[rid]] = ub[m.rxn_index[rid]] = 0.0
        for rid, p in o.ports.items():
            if rid in o.closed or rid in o.ports_closed_by_policy:
                continue
            col = offset + m.rxn_index[rid]
            coeff = m.exchange(rid)[1]
            # organism loses |coeff| per unit flux (secretion if coeff < 0): the pool gains the same
            entries.append((row(f"pool::{p['pool']}"), col,
                            -coeff * float(p.get("scale", 1.0)) * frac[o.id]))
            port_cols.append({"organism": o.id, "exchange": rid, "pool": p["pool"], "col": col,
                              "support": p.get("support", "unknown")})
        col_names += [f"{o.id}::{r}" for r in m.rxn_ids]
        lbs += list(lb)
        ubs += list(ub)
        bio_cols[o.id] = offset + o.biomass_col

    pool_rows = {k[6:]: v for k, v in rows.items() if k.startswith("pool::")}
    env_cols: dict[tuple[str, str], int] = {}
    for pname, prow in pool_rows.items():
        for kind, sign, limits in (("supply", 1.0, supply), ("removal", -1.0, removal)):
            env_cols[(kind, pname)] = len(col_names)
            entries.append((prow, len(col_names), sign))
            col_names.append(f"ENV::{kind}::{pname}")
            lbs.append(0.0)
            ubs.append(float(limits.get(pname, 0.0)))
    mu_col = len(col_names)
    col_names.append("mu")
    lbs.append(0.0)
    ubs.append(1000.0)
    for o in orgs:
        r = row(f"growth::{o.id}")
        entries += [(r, bio_cols[o.id], 1.0), (r, mu_col, -1.0)]
    ii, jj, vv = zip(*entries)
    S = coo_matrix((vv, (ii, jj)), shape=(len(rows), len(col_names))).tocsc()
    names = [k for k, _ in sorted(rows.items(), key=lambda kv: kv[1])]
    return CommunityLP(S, np.array(lbs), np.array(ubs), mu_col, names, col_names, bio_cols, port_cols,
                       env_cols, pool_rows, [o.id for o in orgs])


# --------------------------------------------------------------------------- analysis


def transfers(lp: CommunityLP, x: np.ndarray, orgs: dict[str, Organism], tol: float) -> tuple[
        dict[str, dict[str, float]], list[dict[str, Any]], dict[str, dict[str, float]]]:
    """Net pool exchange per metabolite and organism (positive: organism adds to the pool)."""
    net: dict[str, dict[str, float]] = {}
    for pc in lp.port_cols:
        v = float(lp.S[lp.pool_rows[pc["pool"]], pc["col"]] * x[pc["col"]])
        d = net.setdefault(pc["pool"], {})
        d[pc["organism"]] = d.get(pc["organism"], 0.0) + v
    env = {p: {"supply": float(x[lp.env_cols[("supply", p)]]), "removal": float(x[lp.env_cols[("removal", p)]])}
           for p in lp.pool_rows}
    cross = []
    for p, d in sorted(net.items()):
        prod = {o: v for o, v in d.items() if v > tol}
        cons = {o: -v for o, v in d.items() if v < -tol}
        if any(s != c for s in prod for c in cons):
            cross.append({"metabolite": p, "producers": prod, "consumers": cons,
                          "environment_supply": env[p]["supply"],
                          "transport_support": {o: sorted({pc["support"] for pc in lp.port_cols
                                                           if pc["organism"] == o and pc["pool"] == p})
                                                for o in set(prod) | set(cons)},
                          "units": "mmol/gDW_community/h"})
    return net, cross, env


def solve_common_growth(solver: BudgetedSolver, lp: CommunityLP, label: str,
                        growth_threshold: float) -> dict[str, Any]:
    """Maximise the common growth rate mu of an LP built by ``build_lp`` under a shared solver budget.

    Same semantics as ``Analysis.growth``: infeasibility under the fixed maintenance is reported as
    ``infeasible_with_maintenance`` (never silently as zero growth); solver errors and an exhausted
    budget are ``inconclusive``/``skipped_budget`` with no growth value."""
    solver.model = SimpleNamespace(S=lp.S)
    if solver.remaining < 1:
        return {"status": "skipped_budget", "outcome": "not_run", "mu": None, "call": None,
                "note": f"solver budget of {solver.max_calls} calls exhausted"}
    c = np.zeros(len(lp.col_names))
    c[lp.mu_col] = -1.0
    r = solver.optimize(c, lp.lb, lp.ub, label)
    if r.status == "infeasible":
        return {"status": "completed", "outcome": "infeasible_with_maintenance", "mu": None, "call": r.call,
                "note": "no steady state satisfies maintenance and balanced growth constraints"}
    if r.status != "optimal":
        return {"status": "inconclusive", "outcome": f"solver_{r.status}", "mu": None, "call": r.call,
                "note": r.note}
    mu = max(0.0, -r.objective)
    return {"status": "completed", "outcome": "growth" if mu >= growth_threshold else "no_growth",
            "mu": mu, "call": r.call, "note": ""}


class Analysis:
    def __init__(self, manifest: dict[str, Any], orgs: list[Organism], cfg: dict[str, Any]):
        self.manifest, self.cfg = manifest, cfg
        self.orgs = {o.id: o for o in orgs}
        self.order = [o.id for o in orgs]
        _, self.supply, self.removal = validate_pool(manifest)
        self.ratio = community_ratios(manifest, orgs)
        self.solver = BudgetedSolver(SimpleNamespace(S=None), cfg["max_solver_calls"], cfg["residual_tol"])
        self.records: list[dict[str, Any]] = []

    def lp(self, ids: list[str], closures=None, knockouts=None, supply=None) -> CommunityLP:
        sel = [self.orgs[i] for i in ids]
        return build_lp(sel, self.ratio, supply if supply is not None else self.supply, self.removal,
                        closures, knockouts)

    def solve(self, lp: CommunityLP, c: np.ndarray, lb: np.ndarray, label: str):
        self.solver.model = SimpleNamespace(S=lp.S)
        return self.solver.optimize(c, lb, lp.ub, label)

    def growth(self, eid: str, kind: str, ids: list[str], lp: CommunityLP, mods: dict[str, Any] | None = None,
               cost: int = 1) -> dict[str, Any]:
        cfg = self.cfg
        rec: dict[str, Any] = {"id": eid, "kind": kind, "organisms": ids, "modifications": mods or {},
                               "status": "completed", "outcome": None, "mu": None, "solver_calls": []}
        if self.solver.remaining < cost:
            rec.update(status="skipped_budget", outcome="not_run",
                       reason=f"{self.solver.remaining} solver calls left, {cost} needed")
            self.records.append(rec)
            return rec
        first = len(self.solver.log)
        try:
            c = np.zeros(len(lp.col_names))
            c[lp.mu_col] = -1.0
            r = self.solve(lp, c, lp.lb, f"{eid}:max_mu")
        except BudgetExhausted as exc:
            rec.update(status="inconclusive", outcome="budget_exhausted", note=str(exc))
            self.records.append(rec)
            return rec
        rec["solver_calls"] = list(range(first, len(self.solver.log)))
        mono = kind.startswith("monoculture")
        if r.status == "infeasible":
            rec["outcome"] = "infeasible_with_maintenance"
            rec["note"] = "no steady state satisfies maintenance and balanced growth constraints"
        elif r.status != "optimal":
            rec.update(status="inconclusive", outcome=f"solver_{r.status}", note=r.note)
        else:
            mu = max(0.0, -r.objective)
            rec["mu"] = mu
            grows = mu >= cfg["growth_threshold"]
            rec["outcome"] = ("growth" if mono else "common_growth") if grows else \
                             ("no_growth" if mono else "no_common_growth")
            if grows:
                net, cross, env = transfers(lp, r.x, self.orgs, cfg["flux_tol"])
                rec["net_pool_exchange"] = net
                rec["environment"] = env
                rec["cross_fed"] = cross
                rec["per_organism_growth"] = {o: float(r.x[lp.bio_cols[o]]) for o in ids}
            rec["_x"] = r.x
            rec["_lp"] = lp
        self.records.append(rec)
        return rec

    def fva_transfer(self, rec: dict[str, Any], cand: dict[str, Any]) -> dict[str, Any]:
        lp: CommunityLP = rec["_lp"]
        cols = [(pc["col"]) for pc in lp.port_cols
                if pc["pool"] == cand["metabolite"] and pc["organism"] in cand["producers"]]
        c = np.zeros(len(lp.col_names))
        for j in cols:
            c[j] = lp.S[lp.pool_rows[cand["metabolite"]], j]
        lb = lp.lb.copy()
        lb[lp.mu_col] = max(lb[lp.mu_col], self.cfg["growth_fraction"] * rec["mu"])
        out: dict[str, Any] = {"metabolite": cand["metabolite"], "producers": sorted(cand["producers"]),
                               "mu_floor": float(lb[lp.mu_col]), "status": "completed"}
        if self.solver.remaining < 2:
            return {**out, "status": "skipped_budget"}
        # c . v is the producers' net addition to the pool, so min(-c.v) is the maximum transfer
        hi = self.solve(lp, -c, lb, f"fva:{cand['metabolite']}:max_transfer")
        lo = self.solve(lp, c, lb, f"fva:{cand['metabolite']}:min_transfer")
        if hi.status == lo.status == "optimal":
            out.update(min_transfer=float(lo.objective), max_transfer=float(-hi.objective))
            out["class"] = ("transfer_unavoidable_at_mu_floor_in_model" if lo.objective > self.cfg["flux_tol"]
                            else "transfer_optional_in_model")
        else:
            out.update(status="inconclusive", **{"class": "inconclusive"},
                       lp_status={"max": hi.status, "min": lo.status})
        out["calls"] = [hi.call, lo.call]
        return out


def run_community(manifest: dict[str, Any], base_dir: Path, overrides: dict[str, Any] | None = None) -> dict[str, Any]:
    cfg = {**DEFAULT_CONFIG, **manifest.get("config", {}), **{k: v for k, v in (overrides or {}).items() if v is not None}}
    orgs = load_organisms(manifest, base_dir, cfg)
    A = Analysis(manifest, orgs, cfg)
    tol = cfg["growth_threshold"]

    mono: dict[str, dict[str, Any]] = {}
    for o in orgs:
        mono[o.id] = A.growth(f"monoculture:{o.id}", "monoculture", [o.id], A.lp([o.id]))
    comm = A.growth("community", "community", A.order, A.lp(A.order))

    for exp in manifest.get("experiments", []):
        ids = exp.get("organisms") or A.order
        if any(i not in A.orgs for i in ids):
            raise ManifestError(f"experiment {exp.get('id')}: unknown organism in {ids}")
        closures = [(p["organism"], p["exchange"], p.get("direction", "both")) for p in exp.get("closed_ports", [])]
        for oid, rid, _ in closures:
            if oid not in A.orgs or rid not in A.orgs[oid].ports:
                raise ManifestError(f"experiment {exp.get('id')}: {oid}.{rid} is not a mapped port")
        kos = {k["organism"]: set(k.get("reactions", [])) | A.orgs[k["organism"]].model.blocked_by_genes(
            set(A.orgs[k["organism"]].model.gene_alias.get(g, g) for g in k.get("genes", [])))
            for k in exp.get("knockouts", [])}
        supply = {**A.supply, **_finite_nonneg(exp.get("supply", {}), "experiment.supply")}
        A.growth(f"user:{exp['id']}", "community" if len(ids) > 1 else "monoculture", ids,
                 A.lp(ids, closures, kos, supply),
                 {"closed_ports": exp.get("closed_ports", []), "knockouts": exp.get("knockouts", []),
                  "supply": exp.get("supply", {})})

    closure_results, fva_results = [], []
    if cfg["auto_port_closure"] and comm.get("cross_fed"):
        for cand in comm["cross_fed"]:
            for prod in sorted(cand["producers"]):
                ports = [pc["exchange"] for pc in comm["_lp"].port_cols
                         if pc["organism"] == prod and pc["pool"] == cand["metabolite"]]
                cl = [(prod, rid, "secretion") for rid in ports]
                rec = A.growth(f"closure:{prod}:{cand['metabolite']}:secretion", "closure", A.order,
                               A.lp(A.order, cl), {"closed_secretion": {"organism": prod, "ports": ports}})
                if rec["status"] == "completed" and rec["outcome"] in ("common_growth", "no_common_growth"):
                    rec["transfer_assessment"] = (
                        "transfer_required_in_model" if rec["outcome"] == "no_common_growth"
                        else "transfer_dispensable_or_replaced_in_model")
                    rec["mu_ratio_vs_community"] = (rec["mu"] / comm["mu"]) if comm["mu"] else None
                else:
                    rec["transfer_assessment"] = "no_result"
                closure_results.append(rec)
        for cand in comm["cross_fed"][: cfg["max_fva_candidates"]]:
            fva_results.append(A.fva_transfer(comm, cand))

    # matched mono/co-culture comparison per organism
    summary = []
    for o in orgs:
        m, c = mono[o.id], comm
        if m["status"] != "completed" or c["status"] != "completed" or m["outcome"] == "infeasible_with_maintenance" \
                or c["outcome"] == "infeasible_with_maintenance":
            cls = "undetermined"
        elif m["outcome"] == "no_growth" and c["outcome"] == "common_growth":
            cls = "model_dependent_rescue"
        elif m["outcome"] == "growth" and c["outcome"] == "common_growth":
            cls = "grows_alone_and_in_community"
        elif m["outcome"] == "no_growth":
            cls = "not_rescued_in_model"
        else:
            cls = "growth_lost_in_community"
        cands = [x for x in (c.get("cross_fed") or []) if o.id in x["consumers"]] if cls == "model_dependent_rescue" else []
        summary.append({
            "organism": o.id, "monoculture_mu": m["mu"], "community_mu": c["mu"], "classification": cls,
            "candidate_exchanged_metabolites": [
                {k: v for k, v in x.items() if k != "consumers"} | {"consumed_by_this_organism": x["consumers"][o.id]}
                for x in cands if x["environment_supply"] <= cfg["flux_tol"]],
            "ports_closed_by_policy": o.ports_closed_by_policy,
        })
    conditional = []
    for s in summary:
        unresolved = sorted({sup for x in s["candidate_exchanged_metabolites"]
                             for v in x["transport_support"].values() for sup in v if sup != "curated"})
        conditional.append({
            "evidence_type": "community_crossfeeding_conditional", "organism": s["organism"],
            "bearing": s["classification"], "model_conditional": True,
            "candidate_exchanged_metabolites": [x["metabolite"] for x in s["candidate_exchanged_metabolites"]],
            "transport_support_not_curated": unresolved,
            "enzyme_assignment": "none",
            "interpretation": ("growth without an organism-internal route is feasible in this model if the "
                               "listed metabolite is supplied by a partner through a transporter; this explains "
                               "a missing biosynthetic reaction without assigning any enzyme, and does not show "
                               "the transporter exists or the exchange occurs"),
            "calibrated_probability": None, "promotion": "none"})

    for rec in A.records:
        rec.pop("_x", None)
        rec.pop("_lp", None)
    unresolved_ports = [{"organism": o.id, "exchange": rid, "support": p.get("support", "unknown"),
                         "evidence": p.get("evidence"), "closed_by_policy": rid in o.ports_closed_by_policy}
                        for o in orgs for rid, p in o.ports.items() if p.get("support", "unknown") != "curated"]
    return {
        "schema": SCHEMA, "disclaimer": DISCLAIMER, "budget_policy": BUDGET_POLICY,
        "config": cfg,
        "hashes": {"config_sha256": canonical_hash(cfg), "manifest_sha256": canonical_hash(manifest),
                   "models": {o.id: o.sha256 for o in orgs}},
        "organisms": [{"id": o.id, "model_sha256": o.sha256, "biomass_reaction": o.model.rxn_ids[o.biomass_col],
                       "maintenance": o.maintenance, "ports": o.ports, "closed_ports": sorted(o.closed),
                       "knockouts": sorted(o.knockouts), "ports_closed_by_policy": o.ports_closed_by_policy}
                      for o in orgs],
        "pool": {"metabolites": manifest["pool"]["metabolites"], "supply": A.supply, "removal": A.removal,
                 "units": "mmol/gDW_community/h"},
        "community": {"biomass_ratio": A.ratio,
                      "fractions": {k: v / sum(A.ratio.values()) for k, v in A.ratio.items()}},
        "assumptions": assumptions(A, orgs),
        "unresolved_transport": unresolved_ports,
        "experiments": A.records, "closure_tests": [r["id"] for r in closure_results],
        "transfer_fva": fva_results, "organism_summary": summary, "conditional_evidence": conditional,
        "solver": {"max_calls": A.solver.max_calls, "used": A.solver.used, "calls": A.solver.log},
    }


def assumptions(A: Analysis, orgs: list[Organism]) -> list[str]:
    return [
        "organism blocks are separate; intracellular and organism-extracellular species are namespaced and never merged",
        "organisms interact only through declared ports mapped to pool metabolites; organism exchanges are not "
        "access to the environment",
        f"environment supply (closed unless listed): {A.supply}; removal (closed unless listed, so undeclared "
        f"secretion into the void is impossible): {A.removal}; units mmol/gDW_community/h",
        f"balanced growth: every organism grows at the common rate mu; biomass fractions fixed at "
        f"{ {k: v / sum(A.ratio.values()) for k, v in A.ratio.items()} }; objective max mu",
        "biomass reactions are assumed normalised to 1 gDW (flux = growth rate, 1/h), declared per organism",
        f"ports opened only with support in {A.cfg['open_ports_with_support']}; others closed (listed in "
        "unresolved_transport)",
        "monocultures use fraction 1 and the same pool, ports, supply and removal as the community",
        "metabolite identity checked by declared metabolite id and, where both are declared, chemical formula; "
        "otherwise identity is as declared and unverified",
        "model-conditional only: no gene or enzyme assignment, no probabilities",
    ]


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Community cross-feeding on separate SBML organism models")
    p.add_argument("--manifest", type=Path, required=True)
    p.add_argument("--out", type=Path, default=None)
    p.add_argument("--max-solver-calls", type=int, default=None)
    args = p.parse_args(argv)
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    result = run_community(manifest, args.manifest.resolve().parent, {"max_solver_calls": args.max_solver_calls})
    text = json.dumps(result, indent=2, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text, encoding="utf-8")
    else:
        print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
