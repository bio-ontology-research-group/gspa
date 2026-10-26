#!/usr/bin/env python3
"""
Evidence-guided computational metabolic experiments on a real SBML (FBC) model.

Given a frozen manifest (SBML model, baseline medium, explicit maintenance
constraint, context records with competing candidate hypotheses), this module
builds a small *admissible menu* of in-silico experiments (carbon / nitrogen
source swaps, oxygen removal, supplement removal, knockouts, rescue, product
checks), lets a selector pick from that menu (fixed panel, seeded random,
evidence-guided, or Qwen agent-guided), runs them as scipy/HiGHS LPs under a
strict solver-call budget and emits an auditable JSON record.

What the records mean (and do not mean):
  * Every outcome is conditional on the supplied model, medium and bounds.
    A planner-chosen objective is never independent proof of gene function.
  * Product checks report two model properties at a growth floor: the maximum
    secretion flux (conditional feasibility) and the minimum secretion flux
    (secretion the model cannot avoid when growing).  Neither is observed
    production or expression; a nonzero minimum only says the stoichiometry
    couples the exchange to growth.  Optional flux is not required flux.
  * Solver failure, infeasibility and budget exhaustion are reported as such
    and are never read as "not producible" or as evidence about a function.
  * Boundaries are explicit: only exchanges the manifest declares as external
    (extracellular species) can import nutrients; all other single-metabolite
    reactions (sinks, demands, biomass drains) have uptake closed.
  * No probabilities are produced and no gap-fill is promoted; every record
    carries ``promotion: "none"`` and ``calibrated_probability: null``.
  * Qwen only chooses an experiment id from the menu; it never sees or edits
    model matrices and nothing it returns is executed.

Exchange sign convention is read from the stoichiometry of each exchange
reaction: if the single metabolite has coefficient < 0 (``M_e ->``) then
v < 0 is uptake; if > 0 (``-> M_e``, "reversed" as in the refiner's matrices)
then v > 0 is uptake.  Uptake limits and secretion directions follow.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import random
import re
import sys
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix, csc_matrix

SCHEMA = "gspa.metabolic_experiments/1"

DEFAULT_CONFIG: dict[str, Any] = {
    "strategy": "evidence",        # fixed | random | evidence | qwen
    "seed": 0,
    "max_experiments": 8,          # experiments run after the baseline
    "max_solver_calls": 60,        # strict cap on LP solves, baseline included
    "uptake_rate": 10.0,           # uptake limit for opened exchanges
    "default_bound": 1000.0,       # replaces infinite SBML bounds
    "growth_threshold": 1e-6,      # biomass flux counted as growth
    "flux_tol": 1e-7,
    "residual_tol": 1e-6,          # max |S v| accepted for a solution
    "relative_tol": 0.05,          # growth ratio band for "unchanged"
    "growth_fraction": 0.95,       # growth floor for FVA / product checks
    "product_requirement": "growth_fraction",  # growth_fraction | maintenance_only
    "max_focus_reactions": 2,      # reactions given FVA + knockout checks
    "panel": None,                 # explicit experiment ids for strategy=fixed
    "qwen": {"model": "qwen3.8-27b", "api_base": None, "max_invalid": 2, "max_tokens": 300},
}

KIND_ORDER = [
    "carbon_source", "nitrogen_source", "oxygen_removal", "supplement_removal",
    "knockout", "ko_rescue", "product",
]

DISCLAIMER = (
    "Conditional model predictions under the stated medium, bounds and maintenance "
    "constraint. Not observed phenotypes, not independent proof of gene function, "
    "no calibrated probabilities, no automatic gap-fill promotion."
)


class BudgetExhausted(RuntimeError):
    pass


class ManifestError(ValueError):
    pass


# --------------------------------------------------------------------------- model


@dataclass
class Model:
    path: str
    sha256: str
    met_ids: list[str]
    rxn_ids: list[str]
    S: csc_matrix
    lb: np.ndarray
    ub: np.ndarray
    objective_terms: list[tuple[str, float]]  # active SBML objective as (reaction, coefficient)
    objective_sense: str | None               # "maximize" | "minimize" | None
    gpr: dict[str, list[list[str]]]          # reaction -> DNF over gene product ids
    gene_alias: dict[str, str]               # id or label -> gene product id
    met_compartment: dict[str, str] = field(default_factory=dict)
    n_infinite_bounds: int = 0
    met_formula: dict[str, str] = field(default_factory=dict)   # FBC chemical formulas where declared
    rxn_index: dict[str, int] = field(default_factory=dict)

    def __post_init__(self):
        self.rxn_index = {r: i for i, r in enumerate(self.rxn_ids)}

    def exchange(self, rid: str) -> tuple[int, float]:
        """Return (column, coefficient) of a single-metabolite exchange reaction."""
        if rid not in self.rxn_index:
            raise ManifestError(f"unknown reaction/exchange: {rid}")
        col = self.rxn_index[rid]
        coeff = self.S[:, col].data
        if len(coeff) != 1:
            raise ManifestError(f"not a single-metabolite exchange: {rid}")
        return col, float(coeff[0])

    def exchange_metabolite(self, rid: str) -> str:
        col, _ = self.exchange(rid)
        return self.met_ids[int(self.S[:, col].indices[0])]

    def is_single_metabolite(self, rid: str) -> bool:
        return len(self.S[:, self.rxn_index[rid]].data) == 1

    def secretion_sign(self, rid: str) -> float:
        return -float(np.sign(self.exchange(rid)[1]))

    def blocked_by_genes(self, genes: set[str]) -> set[str]:
        out = set()
        for rid, dnf in self.gpr.items():
            if dnf and not any(all(g not in genes for g in clause) for clause in dnf):
                out.add(rid)
        return out


def gpr_clauses(association) -> list[list[str]]:
    """Read FBC gene identities directly, preserving AND subunits and OR isozymes.

    libSBML's toInfix() uses display labels by default. Labels can differ from
    gene-product IDs and contain spaces; parsing that string loses identities.
    """
    kind = association.getElementName()
    if kind == "geneProductRef":
        return [[association.getGeneProduct()]]
    if kind not in ("and", "or"):
        raise ValueError(f"unsupported FBC association: {kind}")
    children = [gpr_clauses(association.getAssociation(i)) for i in range(association.getNumAssociations())]
    if not children:
        raise ValueError(f"empty FBC {kind} association")
    if kind == "or":
        return [clause for child in children for clause in child]
    result = [[]]
    for child in children:
        result = [list(dict.fromkeys(left + right)) for left in result for right in child]
    return result


def load_sbml_model(path: Path, default_bound: float = 1000.0) -> Model:
    import libsbml

    raw = Path(path).read_bytes()
    document = libsbml.readSBMLFromString(raw.decode("utf-8"))
    log = document.getErrorLog()
    if log.getNumFailsWithSeverity(libsbml.LIBSBML_SEV_ERROR) or log.getNumFailsWithSeverity(libsbml.LIBSBML_SEV_FATAL):
        raise ValueError(log.toString())
    model = document.getModel()
    if model is None or model.getPlugin("fbc") is None:
        raise ValueError("SBML model must use the FBC package (flux bounds are required)")
    species = [s for s in model.getListOfSpecies() if not s.getBoundaryCondition()]
    mets = {s.getId(): i for i, s in enumerate(species)}
    reactions = list(model.getListOfReactions())
    rows, cols, vals, lb, ub = [], [], [], [], []
    gpr: dict[str, list[list[str]]] = {}
    for col, r in enumerate(reactions):
        fbc = r.getPlugin("fbc")
        lo, hi = fbc.getLowerFluxBound(), fbc.getUpperFluxBound()
        if not lo or not hi or model.getParameter(lo) is None or model.getParameter(hi) is None:
            raise ValueError(f"reaction {r.getId()} lacks FBC flux bounds")
        lb.append(model.getParameter(lo).getValue())
        ub.append(model.getParameter(hi).getValue())
        for sign, refs in ((-1.0, r.getListOfReactants()), (1.0, r.getListOfProducts())):
            for ref in refs:
                if ref.getSpecies() in mets:  # boundary species are not mass balanced
                    rows.append(mets[ref.getSpecies()])
                    cols.append(col)
                    vals.append(sign * ref.getStoichiometry())
        assoc = fbc.getGeneProductAssociation()
        if assoc is not None and assoc.isSetAssociation():
            gpr[r.getId()] = gpr_clauses(assoc.getAssociation())
    S = coo_matrix((vals, (rows, cols)), shape=(len(mets), len(reactions))).tocsc()
    lb_a, ub_a = np.array(lb, dtype=float), np.array(ub, dtype=float)
    if np.any(np.isnan(lb_a)) or np.any(np.isnan(ub_a)):
        raise ValueError("model has NaN flux bounds")
    # Only infinities are replaced; finite bounds are kept exactly as written.
    n_inf = int(np.isinf(lb_a).sum() + np.isinf(ub_a).sum())
    lb_a[np.isneginf(lb_a)] = -default_bound
    ub_a[np.isposinf(ub_a)] = default_bound
    if np.any(np.isinf(lb_a)) or np.any(np.isinf(ub_a)):
        raise ValueError("lower bound +inf or upper bound -inf is not a valid flux bound")
    bad = np.flatnonzero(lb_a > ub_a)
    if bad.size:
        raise ValueError("lower bound exceeds upper bound for: " + ", ".join(reactions[i].getId() for i in bad[:5]))
    terms: list[tuple[str, float]] = []
    sense = None
    active = model.getPlugin("fbc").getActiveObjective()
    if active is not None:
        sense = active.getType()
        terms = [(t.getReaction(), float(t.getCoefficient())) for t in active.getListOfFluxObjectives()]
    alias: dict[str, str] = {}
    for gp in model.getPlugin("fbc").getListOfGeneProducts():
        alias[gp.getId()] = gp.getId()
        if gp.isSetLabel():
            alias[gp.getLabel()] = gp.getId()
    return Model(str(path), hashlib.sha256(raw).hexdigest(), list(mets), [r.getId() for r in reactions],
                 S, lb_a, ub_a, terms, sense, gpr, alias,
                 {s.getId(): s.getCompartment() for s in species}, n_inf,
                 {s.getId(): s.getPlugin("fbc").getChemicalFormula() for s in species
                  if s.getPlugin("fbc") is not None and s.getPlugin("fbc").isSetChemicalFormula()})


# --------------------------------------------------------------------------- solver


@dataclass
class LPResult:
    status: str                  # optimal | infeasible | unbounded | error
    objective: float | None
    x: np.ndarray | None
    call: int | None             # index into the solver-call log (None: no solver call made)
    note: str = ""


class BudgetedSolver:
    """HiGHS LP wrapper with a strict call cap, finiteness and mass-balance checks."""

    def __init__(self, model: Model, max_calls: int, residual_tol: float = 1e-6):
        self.model, self.max_calls, self.residual_tol = model, int(max_calls), residual_tol
        self.log: list[dict[str, Any]] = []

    @property
    def used(self) -> int:
        return len(self.log)

    @property
    def remaining(self) -> int:
        return self.max_calls - self.used

    def optimize(self, c: np.ndarray, lb: np.ndarray, ub: np.ndarray, label: str) -> LPResult:
        """Minimize c.v subject to S v = 0 and lb <= v <= ub."""
        if np.any(lb > ub + 1e-12):
            return LPResult("infeasible", None, None, None, "inconsistent bounds (no solver call)")
        if self.remaining <= 0:
            raise BudgetExhausted(f"solver budget of {self.max_calls} calls exhausted")
        res = linprog(c, A_eq=self.model.S, b_eq=np.zeros(self.model.S.shape[0]),
                      bounds=np.column_stack([lb, ub]), method="highs")
        idx = len(self.log)
        status, note = {0: ("optimal", ""), 2: ("infeasible", ""), 3: ("unbounded", "")}.get(
            res.status, ("error", f"solver status {res.status}"))
        x, obj = None, None
        if status == "optimal":
            x = np.asarray(res.x, dtype=float)
            resid = float(np.max(np.abs(self.model.S @ x))) if x.size else 0.0
            viol = float(max(np.max(lb - x), np.max(x - ub), 0.0))
            if not (np.all(np.isfinite(x)) and np.isfinite(res.fun)):
                status, note, x = "error", "non-finite solution", None
            elif resid > self.residual_tol or viol > 1e-6:
                status, note, x = "error", f"solution failed checks (residual {resid:.2e}, bound violation {viol:.2e})", None
            else:
                obj = float(res.fun)
        self.log.append({"call": idx, "label": label, "status": status, "objective": obj, "note": note})
        return LPResult(status, obj, x, idx, note)


# --------------------------------------------------------------------------- conditions


def close_uptake(model: Model, lb: np.ndarray, ub: np.ndarray, rid: str) -> None:
    col, coeff = model.exchange(rid)
    if coeff > 0:
        ub[col] = min(ub[col], 0.0)
    else:
        lb[col] = max(lb[col], 0.0)


def open_uptake(model: Model, lb: np.ndarray, ub: np.ndarray, rid: str, rate: float) -> None:
    col, coeff = model.exchange(rid)
    if coeff > 0:
        ub[col] = max(ub[col], rate)
    else:
        lb[col] = min(lb[col], -rate)


def validate_boundary(model: Model, manifest: dict[str, Any]) -> set[str]:
    """Return the trusted external exchanges declared in ``baseline.exchanges``.

    Each must be a single-metabolite reaction on an extracellular species
    (``baseline.extracellular_compartments``, default ["e"]).  Medium entries and the
    carbon/nitrogen/oxygen lists must all be declared exchanges.
    """
    base = manifest.get("baseline", {})
    declared = base.get("exchanges")
    if not declared:
        raise ManifestError("baseline.exchanges (explicit external exchange reactions) is required")
    ext_comp = set(base.get("extracellular_compartments", ["e"]))
    for rid in declared:
        model.exchange(rid)  # raises if unknown or not a single-metabolite reaction
        met = model.exchange_metabolite(rid)
        if model.met_compartment[met] not in ext_comp:
            raise ManifestError(f"declared exchange {rid} acts on {met} in non-extracellular compartment "
                                f"{model.met_compartment[met]}")
    external = set(declared)
    refs = list(base.get("medium", {})) + list(base.get("carbon_exchanges", [])) + \
        list(base.get("nitrogen_exchanges", [])) + _as_list(base.get("oxygen_exchange"))
    for rid in refs:
        if rid not in external:
            raise ManifestError(f"{rid} is used in the baseline but not declared in baseline.exchanges")
    for rid, rate in base.get("medium", {}).items():
        if not isinstance(rate, (int, float)) or not np.isfinite(rate) or rate < 0:
            raise ManifestError(f"medium uptake limit for {rid} must be a finite number >= 0")
    return external


def baseline_bounds(model: Model, manifest: dict[str, Any], cfg: dict[str, Any],
                    external: set[str] | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Closed-by-default medium: every external exchange not in the medium has uptake closed,
    and every undeclared single-metabolite reaction (sink / demand / biomass drain) has its
    uptake direction closed so it cannot supply metabolites from nothing."""
    external = validate_boundary(model, manifest) if external is None else external
    lb, ub = model.lb.copy(), model.ub.copy()
    medium = manifest.get("baseline", {}).get("medium", {})
    for rid in model.rxn_ids:
        if model.is_single_metabolite(rid) and rid not in medium:
            close_uptake(model, lb, ub, rid)
    for rid, rate in medium.items():
        col, coeff = model.exchange(rid)
        if coeff > 0:
            ub[col] = float(rate)
        else:
            lb[col] = -float(rate)
    m = manifest["objectives"]["maintenance"]
    if m != "none":
        if m["reaction"] not in model.rxn_index:
            raise ManifestError(f"maintenance reaction not in model: {m['reaction']}")
        col = model.rxn_index[m["reaction"]]
        lb[col] = max(lb[col], float(m["lower_bound"]))
    return lb, ub


def apply_mods(model: Model, manifest: dict[str, Any], cfg: dict[str, Any],
               lb: np.ndarray, ub: np.ndarray, mods: list[dict[str, Any]],
               external: set[str] | None = None) -> list[str]:
    """Apply admissible modifications in place; return the reactions knocked out.
    Uptake can only be opened or closed on declared external exchanges."""
    base = manifest.get("baseline", {})
    external = validate_boundary(model, manifest) if external is None else external
    killed: set[str] = set()
    for mod in mods:
        op = mod["op"]
        if "exchange" in mod and mod["exchange"] not in external:
            raise ManifestError(f"{mod['exchange']} is not a declared external exchange")
        if op in ("carbon_swap", "nitrogen_swap"):
            key = "carbon_exchanges" if op == "carbon_swap" else "nitrogen_exchanges"
            for rid in base.get(key, []):
                close_uptake(model, lb, ub, rid)
            open_uptake(model, lb, ub, mod["exchange"], cfg["uptake_rate"])
        elif op == "oxygen_removal":
            close_uptake(model, lb, ub, base["oxygen_exchange"])
        elif op == "close":
            close_uptake(model, lb, ub, mod["exchange"])
        elif op == "open":
            open_uptake(model, lb, ub, mod["exchange"], mod.get("rate", cfg["uptake_rate"]))
        elif op == "knockout":
            killed |= set(mod["reactions"])
            killed |= model.blocked_by_genes(set(mod.get("genes", [])))
        else:
            raise ManifestError(f"unknown modification: {op}")
    for rid in killed:
        col = model.rxn_index[rid]
        lb[col] = ub[col] = 0.0
    return sorted(killed)


# --------------------------------------------------------------------------- menu


def _as_list(v: Any) -> list[str]:
    if v is None:
        return []
    return [v] if isinstance(v, str) else list(v)


def build_menu(model: Model, manifest: dict[str, Any], max_focus: int = DEFAULT_CONFIG["max_focus_reactions"],
               external: set[str] | None = None) -> tuple[list[dict[str, Any]], list[dict[str, str]]]:
    """Derive the admissible experiment menu from context records. Returns (menu, rejections).

    A hypothesis may name candidate genes, model reactions, or both.  Reactions given directly
    are used as-is; genes are mapped through the model GPR when present.  A candidate gene that
    is absent from the model GPR (e.g. masked target assignments) is recorded as unmapped and
    does not remove a reaction-level hypothesis from the menu.
    """
    external = validate_boundary(model, manifest) if external is None else external
    ctx = manifest.get("context", {})
    evidence = ctx.get("evidence", [])
    hyps = {h["id"]: h for h in ctx.get("hypotheses", [])}
    base = manifest.get("baseline", {})
    menu: dict[str, dict[str, Any]] = {}
    rejected: list[dict[str, str]] = []

    ev_by_hyp: dict[str, list[dict[str, Any]]] = {h: [] for h in hyps}
    for ev in evidence:
        hid = ev.get("hypothesis_id")
        if hid in ev_by_hyp:
            ev_by_hyp[hid].append(ev)

    def hyp_reactions(h: dict[str, Any]) -> tuple[list[str], list[str]]:
        genes = [model.gene_alias[g] for g in _as_list(h.get("gene")) + _as_list(h.get("genes"))
                 if g in model.gene_alias]
        rxns = [r for r in _as_list(h.get("reactions")) if r in model.rxn_index]
        return sorted(set(rxns) | model.blocked_by_genes(set(genes))), genes

    def unmapped_genes(h: dict[str, Any]) -> list[str]:
        return sorted(g for g in _as_list(h.get("gene")) + _as_list(h.get("genes")) if g not in model.gene_alias)

    def peers(hid: str) -> list[str]:
        grp = hyps[hid].get("group")
        return sorted(o for o, h in hyps.items() if o != hid and grp is not None and h.get("group") == grp)

    def add(eid: str, kind: str, mods: list[dict[str, Any]], desc: str, ev_ids: list[str],
            hyp_ids: list[str], focus: list[str], product: str | None = None) -> None:
        if eid in menu:
            item = menu[eid]
            item["evidence_ids"] = sorted(set(item["evidence_ids"]) | set(ev_ids))
            item["hypothesis_ids"] = sorted(set(item["hypothesis_ids"]) | set(hyp_ids))
            item["focus_reactions"] = sorted(set(item["focus_reactions"]) | set(focus))
            return
        menu[eid] = {"id": eid, "kind": kind, "mods": mods, "description": desc,
                     "evidence_ids": sorted(set(ev_ids)), "hypothesis_ids": sorted(set(hyp_ids)),
                     "focus_reactions": sorted(set(focus)), "product_exchange": product,
                     "competing_hypotheses": sorted({p for h in hyp_ids for p in peers(h)})}

    def exch_ok(rid: str, why: str) -> bool:
        try:
            model.exchange(rid)
        except ManifestError as exc:
            rejected.append({"item": why, "reason": str(exc)})
            return False
        if rid not in external:
            rejected.append({"item": why, "reason": f"{rid} is not a declared external exchange"})
            return False
        return True

    def linked_focus(ev: dict[str, Any]) -> tuple[list[str], list[str]]:
        hid = ev.get("hypothesis_id")
        if hid in hyps:
            return hyp_reactions(hyps[hid])[0], [hid]
        return [], []

    for ev in evidence:
        etype, eid = ev.get("type"), ev.get("id", "?")
        focus, hids = linked_focus(ev)
        if etype == "transport" and ev.get("substrate_exchange"):
            ex, role = ev["substrate_exchange"], ev.get("role", "carbon")
            if role not in ("carbon", "nitrogen"):
                rejected.append({"item": eid, "reason": f"unsupported transport role: {role}"})
                continue
            if exch_ok(ex, eid):
                op = f"{role}_swap"
                add(f"{role}_source:{ex}", f"{role}_source", [{"op": op, "exchange": ex}],
                    f"replace baseline {role} sources by {ex}", [eid], hids, focus)
        elif etype == "bgc" and ev.get("product_exchange"):
            ex = ev["product_exchange"]
            if exch_ok(ex, eid):
                add(f"product:{ex}", "product", [], f"conditional feasibility vs. forced production of {ex}",
                    [eid], hids, focus, product=ex)
        hint = ev.get("condition_hint")
        if hint == "oxygen_removal" and base.get("oxygen_exchange") and exch_ok(base["oxygen_exchange"], eid):
            add("oxygen_removal", "oxygen_removal", [{"op": "oxygen_removal"}], "remove oxygen uptake",
                [eid], hids, focus)
        elif hint == "supplement_removal" and ev.get("exchange") and exch_ok(ev["exchange"], eid):
            add(f"supplement_removal:{ev['exchange']}", "supplement_removal",
                [{"op": "close", "exchange": ev["exchange"]}], f"remove supplement uptake {ev['exchange']}",
                [eid], hids, focus)

    for hid, h in hyps.items():
        rxns, genes = hyp_reactions(h)
        ev_ids = [e.get("id", "?") for e in ev_by_hyp[hid]]
        if not rxns and not genes:
            rejected.append({"item": hid, "reason": "hypothesis resolves to no model reaction"})
            continue
        ko = {"op": "knockout", "reactions": rxns, "genes": genes}
        cond: list[dict[str, Any]] = []
        c = h.get("condition") or {}
        for key, op in (("carbon_exchange", "carbon_swap"), ("nitrogen_exchange", "nitrogen_swap")):
            if c.get(key) and exch_ok(c[key], hid):
                cond.append({"op": op, "exchange": c[key]})
        suffix = "" if not cond else "@" + ",".join(m["exchange"] for m in cond)
        add(f"knockout:{hid}{suffix}", "knockout", cond + [ko], f"knock out {','.join(genes or rxns)}"
            + (f" under {suffix[1:]}" if cond else ""), ev_ids, [hid], rxns)
        resc = h.get("rescue_exchange")
        if resc and exch_ok(resc, hid):
            add(f"ko_rescue:{hid}:unrescued", "ko_rescue", cond + [ko, {"op": "close", "exchange": resc}],
                f"knock out {','.join(rxns)} with {resc} closed", ev_ids, [hid], rxns)
            add(f"ko_rescue:{hid}:rescued", "ko_rescue", cond + [ko, {"op": "open", "exchange": resc}],
                f"knock out {','.join(rxns)} with {resc} supplied", ev_ids, [hid], rxns)

    for item in menu.values():
        hs = [hyps[h] for h in item["hypothesis_ids"]]
        item["candidate_genes"] = sorted({g for h in hs for g in _as_list(h.get("gene")) + _as_list(h.get("genes"))})
        item["unmapped_candidate_genes"] = sorted({g for h in hs for g in unmapped_genes(h)})
        item["candidate_reactions"] = sorted({r for h in hs for r in hyp_reactions(h)[0]})
        item["min_calls"], item["max_calls"] = call_bounds(item, max_focus)
    return sorted(menu.values(), key=lambda m: (KIND_ORDER.index(m["kind"]), m["id"])), rejected


def call_bounds(item: dict[str, Any], max_focus: int) -> tuple[int, int]:
    """Solver-call cost policy, identical for every selection strategy.

    mandatory = 1 growth LP (+2 product LPs when the item has a product exchange);
    optional focus checks = 3 LPs (FVA min, FVA max, knockout growth) per focus reaction not
    knocked out by the item itself, at most ``max_focus`` reactions.
    """
    mandatory = 1 + (2 if item.get("product_exchange") else 0)
    killed = {r for m in item["mods"] if m["op"] == "knockout" for r in m["reactions"]}
    n_focus = min(len([r for r in item["focus_reactions"] if r not in killed]), max_focus)
    return mandatory, mandatory + 3 * n_focus


# --------------------------------------------------------------------------- execution


def _growth_class(ratio: float, cfg: dict[str, Any], growth: float) -> str:
    if growth < cfg["growth_threshold"]:
        return "growth_lost"
    if ratio < 1 - cfg["relative_tol"]:
        return "growth_reduced"
    if ratio > 1 + cfg["relative_tol"]:
        return "growth_increased"
    return "growth_unchanged"


def resolve_objective(model: Model, manifest: dict[str, Any]) -> tuple[str, str]:
    """Biomass reaction id and where it came from.  Never guesses.

    ``objectives.biomass`` is preferred.  Otherwise the SBML active objective is accepted only
    if it is a single maximised term with positive coefficient; anything else is rejected.
    """
    explicit = manifest["objectives"].get("biomass")
    if explicit:
        if explicit not in model.rxn_index:
            raise ManifestError(f"biomass reaction not in model: {explicit}")
        return explicit, "manifest"
    terms = model.objective_terms
    if len(terms) != 1 or terms[0][1] <= 0 or (model.objective_sense or "").lower() != "maximize":
        raise ManifestError("SBML objective is not a single maximised positive-coefficient term "
                            f"({terms}, {model.objective_sense}); set objectives.biomass explicitly")
    if terms[0][0] not in model.rxn_index:
        raise ManifestError(f"SBML objective reaction not in model: {terms[0][0]}")
    return terms[0][0], "sbml_active_objective"


def validate_condition_constraints(model: Model, manifest: dict[str, Any], external: set[str]) -> list[dict]:
    """Validate independently specified regulatory closures, never infer them from outcomes."""
    constraints = manifest.get("condition_constraints", [])
    if not isinstance(constraints, list):
        raise ManifestError("condition_constraints must be a list")
    for con in constraints:
        if not isinstance(con, dict) or con.get("schema") != "gspa.v3.condition_constraint/1" or \
                con.get("type") != "regulatory_reaction_closure":
            raise ManifestError("unsupported condition constraint schema or type")
        if con.get("status") != "fully_specified":
            raise ManifestError("condition constraint is not fully specified")
        if con.get("bounds") != {"lower": 0.0, "upper": 0.0, "units": "mmol/gDW/h"}:
            raise ManifestError("regulatory closure requires zero bounds in mmol/gDW/h")
        if con.get("substrate_exchange") not in external:
            raise ManifestError("regulatory substrate is not a declared external exchange")
        if con.get("applies_to") != ["wild_type_control", "knockout"]:
            raise ManifestError("regulatory closure must apply to both control and knockout")
        reactions = con.get("closed_reactions")
        if not isinstance(reactions, list) or any(not isinstance(r, str) for r in reactions) or \
                len(reactions) != len(set(reactions)) or any(r not in model.rxn_index for r in reactions):
            raise ManifestError("regulatory closures require unique mapped reaction identifiers")
        if not isinstance(con.get("source"), dict) or not con["source"]:
            raise ManifestError("regulatory closure requires source provenance")
    return constraints


def apply_condition_constraints(model: Model, constraints: list[dict], lb, ub) -> list[dict]:
    """Apply closures only when their substrate can be taken up in this experiment.

    Intersect with existing constraints: a closure cannot erase an ATP maintenance
    or biomass lower bound. Activation is evaluated before applying any closure.
    """
    active = []
    for con in constraints:
        col, coeff = model.exchange(con["substrate_exchange"])
        if (lb[col] < 0 if coeff < 0 else ub[col] > 0):
            active.append(con)
    for con in active:
        for rid in con["closed_reactions"]:
            col = model.rxn_index[rid]
            lb[col], ub[col] = max(lb[col], 0.0), min(ub[col], 0.0)
    return active


class Runner:
    def __init__(self, model: Model, manifest: dict[str, Any], cfg: dict[str, Any]):
        self.model, self.manifest, self.cfg = model, manifest, cfg
        if cfg.get("product_requirement", "growth_fraction") not in ("growth_fraction", "maintenance_only"):
            raise ManifestError("product_requirement must be growth_fraction or maintenance_only")
        fraction = cfg["growth_fraction"]
        if isinstance(fraction, bool) or not isinstance(fraction, (int, float)) or not math.isfinite(fraction) or not 0 < fraction <= 1:
            raise ManifestError("growth_fraction must be finite and in (0, 1]")
        self.solver = BudgetedSolver(model, cfg["max_solver_calls"], cfg["residual_tol"])
        obj, self.objective_source = resolve_objective(model, manifest)
        self.bio = model.rxn_index[obj]
        self.external = validate_boundary(model, manifest)
        self.condition_constraints = validate_condition_constraints(model, manifest, self.external)
        self.base_lb, self.base_ub = baseline_bounds(model, manifest, cfg, self.external)
        self.baseline_growth: float | None = None
        self.baseline_outcome: str | None = None

    def _unit(self, col: int, sign: float) -> np.ndarray:
        c = np.zeros(len(self.model.rxn_ids))
        c[col] = sign
        return c

    def _growth(self, lb, ub, label) -> tuple[LPResult, float]:
        r = self.solver.optimize(self._unit(self.bio, -1.0), lb, ub, label)
        return r, (max(0.0, -r.objective) if r.status == "optimal" else 0.0)

    def run(self, item: dict[str, Any]) -> dict[str, Any]:
        cfg, model = self.cfg, self.model
        rec: dict[str, Any] = {
            "experiment_id": item["id"], "kind": item["kind"], "description": item["description"],
            "mods": item["mods"], "evidence_ids": item["evidence_ids"],
            "hypothesis_ids": item["hypothesis_ids"],
            "competing_hypotheses": item.get("competing_hypotheses", []),
            "candidate_genes": item.get("candidate_genes", []),
            "unmapped_candidate_genes": item.get("unmapped_candidate_genes", []),
            "candidate_reactions": item.get("candidate_reactions", []),
            "status": "completed", "outcome": None, "measurements": {}, "solver_calls": [],
        }
        product = item.get("product_exchange")
        needed = item["min_calls"] if "min_calls" in item else call_bounds(item, cfg["max_focus_reactions"])[0]
        reserve = needed - 1  # product LPs are held back so focus checks cannot starve them
        if self.solver.remaining < needed:
            rec.update(status="skipped_budget", outcome="not_run",
                       reason=f"{self.solver.remaining} solver calls left, {needed} needed")
            return self._finish(rec)
        lb, ub = self.base_lb.copy(), self.base_ub.copy()
        first = len(self.solver.log)
        try:
            killed = apply_mods(model, self.manifest, cfg, lb, ub, item["mods"], self.external)
            rec["condition_constraints"] = apply_condition_constraints(model, self.condition_constraints, lb, ub)
            rec["measurements"]["knocked_out_reactions"] = killed
            r, growth = self._growth(lb, ub, f"{item['id']}:max_growth")
            m = rec["measurements"]
            if r.status == "infeasible":
                rec["outcome"] = "infeasible_with_maintenance"
                if item["kind"] == "baseline":
                    self.baseline_outcome = rec["outcome"]
                rec["note"] = "no steady state satisfies the maintenance constraint under this condition"
                m["max_growth"] = None
                return self._finish(rec, first)
            if r.status != "optimal":
                rec.update(status="inconclusive", outcome=f"solver_{r.status}", note=r.note)
                return self._finish(rec, first)
            m["max_growth"] = growth
            if item["kind"] == "baseline":
                self.baseline_growth = growth
                rec["outcome"] = "growth" if growth >= cfg["growth_threshold"] else "no_growth_feasible"
                self.baseline_outcome = rec["outcome"]
            else:
                bg = self.baseline_growth
                m["baseline_growth"] = bg
                m["baseline_outcome"] = self.baseline_outcome
                m["growth_ratio"] = growth / bg if bg is not None and bg >= cfg["growth_threshold"] else None
                if m["growth_ratio"] is not None:
                    rec["outcome"] = _growth_class(m["growth_ratio"], cfg, growth)
                elif self.baseline_outcome in ("no_growth_feasible", "infeasible_with_maintenance"):
                    rec["outcome"] = "growth_recovered" if growth >= cfg["growth_threshold"] else "no_growth_feasible"
                else:
                    # A failed or unrun reference is unknown, not zero growth.
                    rec["outcome"] = "growth_without_baseline" if growth >= cfg["growth_threshold"] else "no_growth_feasible"
            if item["kind"] in ("knockout", "ko_rescue") and r.x is not None:
                alt = [p for p in self._peer_reactions(item) if p not in killed
                       and abs(r.x[model.rxn_index[p]]) > cfg["flux_tol"]]
                if alt and rec["outcome"] != "growth_lost":
                    m["alternative_route_via"] = alt
            if growth >= cfg["growth_threshold"]:
                self._focus_checks(rec, item, lb, ub, growth, reserve)
            if product:
                if growth >= cfg["growth_threshold"] or cfg.get("product_requirement") == "maintenance_only":
                    self._product_checks(rec, product, lb, ub, growth)
                else:
                    rec["measurements"]["product"] = {"status": "not_assessed_no_growth"}
        except BudgetExhausted as exc:
            rec.update(status="inconclusive", outcome="budget_exhausted_mid_experiment", note=str(exc))
        return self._finish(rec, first)

    def _finish(self, rec: dict[str, Any], first: int | None = None) -> dict[str, Any]:
        if first is not None:
            rec["solver_calls"] = list(range(first, len(self.solver.log)))
        rec["conditional_evidence"] = conditional_record(rec)
        return rec

    def _fixed_growth_lb(self, lb: np.ndarray, growth: float) -> np.ndarray:
        out = lb.copy()
        out[self.bio] = max(out[self.bio], self.cfg["growth_fraction"] * growth)
        return out

    def _focus_checks(self, rec, item, lb, ub, growth, reserve) -> None:
        cfg, model = self.cfg, self.model
        killed = set(rec["measurements"].get("knocked_out_reactions", []))
        focus = [r for r in item["focus_reactions"] if r not in killed][: cfg["max_focus_reactions"]]
        if not focus:
            return
        afford = max(0, (self.solver.remaining - reserve) // 3)
        checks, floor_lb = [], self._fixed_growth_lb(lb, growth)
        if afford < len(focus):
            rec["measurements"]["focus_checks_truncated_by_budget"] = len(focus) - afford
        for rid in focus[:afford]:
            col = model.rxn_index[rid]
            lo = self.solver.optimize(self._unit(col, 1.0), floor_lb, ub, f"{item['id']}:fva_min:{rid}")
            hi = self.solver.optimize(self._unit(col, -1.0), floor_lb, ub, f"{item['id']}:fva_max:{rid}")
            ko_lb, ko_ub = lb.copy(), ub.copy()
            ko_lb[col] = ko_ub[col] = 0.0
            kr, kg = self._growth(ko_lb, ko_ub, f"{item['id']}:ko:{rid}")
            entry: dict[str, Any] = {"reaction": rid}
            if lo.status == hi.status == "optimal":
                vmin, vmax = lo.objective, -hi.objective
                entry.update(fva_min=vmin, fva_max=vmax)
                tol = cfg["flux_tol"]
                if abs(vmin) <= tol and abs(vmax) <= tol:
                    entry["flux_class"] = "carries_no_flux_at_growth_floor"
                elif vmin > tol or vmax < -tol:
                    entry["flux_class"] = "required_flux"
                else:
                    entry["flux_class"] = "optional_flux"
            else:
                entry["flux_class"] = "inconclusive"
            if kr.status == "optimal":
                entry["knockout_growth"] = kg
                entry["knockout_growth_ratio"] = kg / growth
                entry["knockout_class"] = ("knockout_dispensable" if kg >= growth * (1 - cfg["relative_tol"])
                                           else "knockout_reduces_growth" if kg >= cfg["growth_threshold"]
                                           else "knockout_abolishes_growth")
                if entry["knockout_class"] == "knockout_dispensable" and kr.x is not None:
                    peers = [r for r in self._peer_reactions(item) if r != rid
                             and abs(kr.x[model.rxn_index[r]]) > cfg["flux_tol"]]
                    if peers:
                        entry["alternative_route_via"] = peers
            elif kr.status == "infeasible":
                entry["knockout_class"] = "knockout_infeasible_with_maintenance"
            else:
                entry["knockout_class"] = "inconclusive"
            checks.append(entry)
        rec["measurements"]["focus_checks"] = checks

    def _peer_reactions(self, item) -> list[str]:
        hyps = {h["id"]: h for h in self.manifest.get("context", {}).get("hypotheses", [])}
        out: set[str] = set()
        for pid in item.get("competing_hypotheses", []):
            h = hyps[pid]
            out |= {r for r in _as_list(h.get("reactions")) if r in self.model.rxn_index}
            out |= self.model.blocked_by_genes({self.model.gene_alias[g] for g in
                                                _as_list(h.get("gene")) + _as_list(h.get("genes"))
                                                if g in self.model.gene_alias})
        return sorted(out)

    def _product_checks(self, rec, ex: str, lb, ub, growth) -> None:
        cfg, model = self.cfg, self.model
        col, sec = model.exchange(ex)[0], model.secretion_sign(ex)
        requirement = cfg.get("product_requirement", "growth_fraction")
        floor_lb = self._fixed_growth_lb(lb, growth) if requirement == "growth_fraction" else lb.copy()
        hi = self.solver.optimize(self._unit(col, -sec), floor_lb, ub, f"{rec['experiment_id']}:product_max")
        lo = self.solver.optimize(self._unit(col, sec), floor_lb, ub, f"{rec['experiment_id']}:product_min")
        info: dict[str, Any] = {"exchange": ex, "growth_floor": float(floor_lb[self.bio]),
                                "requirement": requirement, "maintenance": self.manifest["objectives"]["maintenance"]}
        if hi.status == lo.status == "optimal":
            # the two LPs minimise -secretion and +secretion respectively
            pmax, pmin = max(0.0, -hi.objective), max(0.0, lo.objective)
            info.update(max_secretion=pmax, min_secretion_at_required_state=pmin)
            if requirement == "growth_fraction":
                info["min_secretion_at_growth_floor"] = pmin
            tol = cfg["flux_tol"]
            info["conditional_feasibility"] = pmax > tol
            info["secretion_unavoidable_in_model"] = pmin > tol
            info["class"] = (("no_secretion_possible_in_model_at_growth_floor" if requirement == "growth_fraction"
                              else "no_secretion_possible_in_model_at_maintenance") if pmax <= tol
                             else "secretion_unavoidable_in_model" if pmin > tol
                             else "secretion_feasible_not_unavoidable")
        else:
            # solver failure / infeasibility is not a statement about the product
            info["class"] = "inconclusive"
            info["lp_status"] = {"max": hi.status, "min": lo.status}
        rec["measurements"]["product"] = info


def baseline_item() -> dict[str, Any]:
    return {"id": "baseline", "kind": "baseline", "mods": [], "description": "baseline medium with maintenance",
            "evidence_ids": [], "hypothesis_ids": [], "focus_reactions": [], "product_exchange": None,
            "competing_hypotheses": []}


def conditional_record(rec: dict[str, Any]) -> dict[str, Any]:
    """Auditable, deliberately non-committal statement derived from one experiment."""
    m = rec["measurements"]
    bearing = "uninformative"
    for fc in m.get("focus_checks", []):
        if fc.get("knockout_class") in ("knockout_reduces_growth", "knockout_abolishes_growth") or \
                fc.get("flux_class") == "required_flux":
            bearing = "model_requires_reaction_in_condition"
        elif bearing == "uninformative" and (fc.get("knockout_class") == "knockout_dispensable"
                                            or fc.get("flux_class") == "optional_flux"):
            bearing = "model_does_not_require_reaction_in_condition"
    if rec["kind"] == "knockout" and rec.get("outcome") in ("growth_lost", "growth_reduced"):
        bearing = "model_requires_reaction_in_condition"
    elif rec["kind"] == "knockout" and rec.get("outcome") == "growth_unchanged":
        bearing = "model_does_not_require_reaction_in_condition"
    # infeasible / failed / skipped runs say nothing about the candidate function
    if rec["status"] != "completed" or rec.get("outcome") == "infeasible_with_maintenance":
        bearing = "no_result"
    return {"evidence_type": "metabolic_experiment_conditional",
            "bearing": bearing, "model_conditional": True, "independent_proof_of_function": False,
            "calibrated_probability": None, "promotion": "none",
            "status": rec["status"], "outcome": rec.get("outcome"),
            "hypothesis_ids": rec.get("hypothesis_ids", []), "evidence_ids": rec.get("evidence_ids", []),
            "candidate_genes": rec.get("candidate_genes", []),
            "unmapped_candidate_genes": rec.get("unmapped_candidate_genes", []),
            "candidate_reactions": rec.get("candidate_reactions", []),
            "note": "Model-conditional; not observed phenotype, expression or production. "
                    "A required/dispensable call reflects this model's network, which may omit alternatives."}


# --------------------------------------------------------------------------- selection


def evidence_key(item: dict[str, Any], evidence: list[dict[str, Any]]) -> tuple:
    types = {e.get("type") for e in evidence if e.get("id") in item["evidence_ids"]}
    return (-len(types), -len(item["evidence_ids"]), -len(item["competing_hypotheses"]), item["id"])


def order_menu(strategy: str, menu: list[dict[str, Any]], manifest: dict[str, Any],
               cfg: dict[str, Any]) -> tuple[list[dict[str, Any]], list[dict[str, str]]]:
    rejected: list[dict[str, str]] = []
    if strategy == "fixed":
        panel = cfg.get("panel")
        if not panel:
            return list(menu), rejected
        by_id = {m["id"]: m for m in menu}
        out = []
        for pid in panel:
            if pid in by_id:
                out.append(by_id[pid])
            else:
                rejected.append({"item": pid, "reason": "not in admissible menu"})
        return out, rejected
    if strategy == "random":
        out = list(menu)
        random.Random(cfg["seed"]).shuffle(out)
        return out, rejected
    if strategy == "evidence":
        ev = manifest.get("context", {}).get("evidence", [])
        return sorted(menu, key=lambda m: evidence_key(m, ev)), rejected
    raise ManifestError(f"unknown strategy: {strategy}")


def load_api_key() -> str:
    """Same convention as benchmark/agent/run_gspa_agent.py; the key is never printed."""
    env_key = os.environ.get("BORG_LLM_API_KEY") or os.environ.get("OPENAI_API_KEY")
    if env_key:
        return env_key.strip()
    key_file = Path.home() / ".config" / "borg-llm" / "api-key"
    if key_file.exists():
        return key_file.read_text(encoding="utf-8").strip()
    return "EMPTY"


def make_qwen_chat(cfg: dict[str, Any]) -> Callable[[list[dict[str, str]]], str]:
    q = cfg["qwen"]
    base = (q.get("api_base") or os.environ.get("GSPA_LLM_URL", "http://unimatrix01.kaust.edu.sa:8000/v1")).rstrip("/")
    key = load_api_key()

    def chat(messages: list[dict[str, str]]) -> str:
        payload = {"model": q["model"], "messages": messages, "temperature": 0.0, "seed": cfg["seed"],
                   "max_tokens": q.get("max_tokens", 300), "chat_template_kwargs": {"enable_thinking": False}}
        req = urllib.request.Request(f"{base}/chat/completions", data=json.dumps(payload).encode(),
                                     headers={"Content-Type": "application/json", "Authorization": f"Bearer {key}"},
                                     method="POST")
        with urllib.request.urlopen(req, timeout=180) as resp:
            return json.loads(resp.read().decode())["choices"][0]["message"]["content"]

    return chat


def parse_choice(text: str, admissible: set[str]) -> tuple[str | None, str, str]:
    """Return (experiment_id | None for stop, rationale, error)."""
    m = re.search(r"\{.*\}", text or "", re.S)
    if not m:
        return None, "", "no JSON object in response"
    try:
        obj = json.loads(m.group(0))
    except json.JSONDecodeError:
        return None, "", "malformed JSON"
    if not isinstance(obj, dict):
        return None, "", "response is not an object"
    why = str(obj.get("reason", ""))[:300]
    if obj.get("stop") is True:
        return None, why, ""
    eid = obj.get("experiment_id")
    if not isinstance(eid, str) or eid not in admissible:
        return None, why, f"experiment_id not in admissible menu: {str(eid)[:80]}"
    return eid, why, ""


def qwen_prompt(remaining: list[dict[str, Any]], results: list[dict[str, Any]], budget_left: int) -> list[dict[str, str]]:
    menu = [{"experiment_id": m["id"], "kind": m["kind"], "description": m["description"],
             "evidence_ids": m["evidence_ids"], "competing_hypotheses": m["competing_hypotheses"],
             "min_calls": m["min_calls"], "max_calls": m["max_calls"]}
            for m in remaining]
    prior = [{"experiment_id": r["experiment_id"], "status": r["status"], "outcome": r["outcome"],
              "growth_ratio": r["measurements"].get("growth_ratio")} for r in results]
    system = ("You choose the next in-silico metabolic experiment from a fixed menu to help discriminate "
              "competing gene-function hypotheses. Results are model-conditional, not proof. Reply with one JSON "
              'object: {"experiment_id": "<id from menu>", "reason": "<short>"} or {"stop": true, "reason": "<short>"}.')
    user = json.dumps({"admissible_menu": menu, "completed_results": prior, "solver_calls_left": budget_left})
    return [{"role": "system", "content": system}, {"role": "user", "content": user}]


# --------------------------------------------------------------------------- driver


BUDGET_POLICY = {
    "applies_to": "all strategies (fixed, random, evidence, qwen) identically",
    "baseline": "1 LP (maximise biomass); always run first when the menu is non-empty",
    "per_experiment_mandatory": "1 growth LP, plus 2 LPs (max and min secretion) for product experiments",
    "per_focus_reaction": "3 LPs (FVA min, FVA max, knockout growth) for up to config.max_focus_reactions "
                          "reactions of the linked hypotheses that the experiment does not itself knock out",
    "skip_rule": "an experiment is skipped (skipped_budget) if remaining calls < its mandatory calls",
    "focus_truncation": "focus checks run for floor((remaining - product_LPs) / 3) reactions after the growth "
                        "LP; any shortfall is recorded as focus_checks_truncated_by_budget",
    "no_planner_reserve": "no calls are held back for any planner; qwen requests do not consume solver calls",
    "experiment_limit": "config.max_experiments counts every attempted menu item, skipped_budget included",
}


def assumptions(model: Model, manifest: dict[str, Any], cfg: dict[str, Any], runner: Runner) -> list[str]:
    base = manifest.get("baseline", {})
    internal = sorted(r for r in model.rxn_ids if model.is_single_metabolite(r) and r not in runner.external)
    return [
        f"external exchanges (trusted, declared): {sorted(runner.external)}",
        f"medium uptake limits: {base.get('medium', {})}; all other declared exchanges closed to uptake",
        f"undeclared single-metabolite reactions with uptake closed: {internal}",
        f"opened exchanges use uptake_rate={cfg['uptake_rate']}",
        f"maintenance constraint: {manifest['objectives']['maintenance']}",
        f"biomass objective {model.rxn_ids[runner.bio]} (source: {runner.objective_source}); maximised, "
        "growth counted above growth_threshold",
        f"product checks require {cfg.get('product_requirement', 'growth_fraction')}; "
        "maintenance-only checks retain all original model bounds but impose no additional biomass floor",
        f"{model.n_infinite_bounds} infinite SBML bounds replaced by +-{cfg['default_bound']}; finite bounds unchanged",
        "all outcomes are model-conditional and not observations of growth, expression or production",
    ]


def canonical_hash(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, default=str).encode()).hexdigest()


def resolve_config(manifest: dict[str, Any], overrides: dict[str, Any] | None = None) -> dict[str, Any]:
    cfg = {**DEFAULT_CONFIG, **manifest.get("config", {}), **{k: v for k, v in (overrides or {}).items() if v is not None}}
    cfg["qwen"] = {**DEFAULT_CONFIG["qwen"], **(manifest.get("config", {}).get("qwen") or {})}
    return cfg


def validate_manifest(manifest: dict[str, Any]) -> None:
    for key in ("model", "objectives"):
        if key not in manifest:
            raise ManifestError(f"manifest lacks '{key}'")
    m = manifest["objectives"].get("maintenance")
    if m != "none" and not (isinstance(m, dict) and "reaction" in m and "lower_bound" in m):
        raise ManifestError("objectives.maintenance must be {reaction, lower_bound} or the string 'none' (explicit)")


def run_experiments(manifest: dict[str, Any], model: Model, overrides: dict[str, Any] | None = None,
                    chat: Callable[[list[dict[str, str]]], str] | None = None) -> dict[str, Any]:
    validate_manifest(manifest)
    cfg = resolve_config(manifest, overrides)
    strategy = cfg["strategy"]
    runner = Runner(model, manifest, cfg)
    menu, menu_rejections = build_menu(model, manifest, cfg["max_focus_reactions"], runner.external)
    out: dict[str, Any] = {
        "schema": SCHEMA, "disclaimer": DISCLAIMER, "strategy": strategy, "seed": cfg["seed"],
        "model": {"path": model.path, "sha256": model.sha256, "n_metabolites": len(model.met_ids),
                  "n_reactions": len(model.rxn_ids), "biomass_reaction": model.rxn_ids[runner.bio],
                  "maintenance": manifest["objectives"]["maintenance"],
                  "objective_source": runner.objective_source, "n_infinite_bounds_replaced": model.n_infinite_bounds},
        "budget_policy": BUDGET_POLICY,
        "assumptions": assumptions(model, manifest, cfg, runner),
        "config": {k: v for k, v in cfg.items()},
        "hashes": {"model_sha256": model.sha256, "config_sha256": canonical_hash(cfg),
                   "context_sha256": canonical_hash(manifest.get("context", {})),
                   "manifest_sha256": canonical_hash(manifest)},
        "menu": [{k: v for k, v in m.items() if k != "mods"} | {"mods": m["mods"]} for m in menu],
        "menu_rejections": menu_rejections, "selection_log": [], "experiments": [],
    }
    results: list[dict[str, Any]] = out["experiments"]
    if not menu:
        out.update(outcome="no_candidate_experiments",
                   reason="context records yield no admissible experiment; nothing was solved")
        return _finalise(out, runner, manifest)

    base = runner.run(baseline_item())
    results.append(base)
    # Alternative media are especially useful when the reference cannot grow.
    # Keep absolute condition results even without a usable relative baseline.

    limit = int(cfg["max_experiments"])
    if strategy == "qwen":
        chat = chat or make_qwen_chat(cfg)
        remaining, invalid, fallback = list(menu), 0, False
        fb_order, _ = order_menu("evidence", menu, manifest, cfg)
        while remaining and len(results) - 1 < limit:
            if fallback:
                item, entry = next(m for m in fb_order if m in remaining), {"source": "evidence_fallback"}
            else:
                try:
                    text = chat(qwen_prompt(remaining, results, runner.solver.remaining))
                    eid, why, err = parse_choice(text, {m["id"] for m in remaining})
                except Exception as exc:  # network/endpoint failure; never echo headers or keys
                    eid, why, err = None, "", f"endpoint error: {type(exc).__name__}"
                entry = {"source": "qwen", "rationale": why}
                if err:
                    invalid += 1
                    entry.update(valid=False, error=err)
                    out["selection_log"].append(entry)
                    if invalid >= cfg["qwen"]["max_invalid"]:
                        fallback = True
                        out["selection_log"].append({"source": "evidence_fallback", "note": "qwen output invalid; deterministic evidence ordering used for remaining steps"})
                    continue
                if eid is None:
                    entry.update(valid=True, stop=True)
                    out["selection_log"].append(entry)
                    break
                item = next(m for m in remaining if m["id"] == eid)
                entry["valid"] = True
            entry["experiment_id"] = item["id"]
            out["selection_log"].append(entry)
            remaining.remove(item)
            results.append(runner.run(item))
        for m in remaining:
            results.append(_unrun(m, "not_selected"))
    else:
        ordered, rej = order_menu(strategy, menu, manifest, cfg)
        out["menu_rejections"] += rej
        for item in ordered:
            if len(results) - 1 >= limit:
                break
            out["selection_log"].append({"source": strategy, "experiment_id": item["id"]})
            results.append(runner.run(item))
        ran = {r["experiment_id"] for r in results}
        for m in menu:
            if m["id"] not in ran:
                results.append(_unrun(m, "not_selected"))
    out["outcome"] = "completed" if all(
        r["status"] in ("completed", "skipped_not_selected") for r in results) else "completed_with_skips_or_inconclusive"
    return _finalise(out, runner, manifest)


def _unrun(m: dict[str, Any], why: str) -> dict[str, Any]:
    rec = {"experiment_id": m["id"], "kind": m["kind"], "status": f"skipped_{why}", "outcome": "not_run",
           "evidence_ids": m["evidence_ids"], "hypothesis_ids": m["hypothesis_ids"], "measurements": {},
           "solver_calls": [], "candidate_genes": m.get("candidate_genes", []),
           "unmapped_candidate_genes": m.get("unmapped_candidate_genes", []),
           "candidate_reactions": m.get("candidate_reactions", [])}
    rec["conditional_evidence"] = conditional_record({**rec, "mods": []})
    return rec


def _finalise(out: dict[str, Any], runner: Runner, manifest: dict[str, Any]) -> dict[str, Any]:
    out["solver"] = {"max_calls": runner.solver.max_calls, "used": runner.solver.used, "calls": runner.solver.log}
    all_ev = [e.get("id") for e in manifest.get("context", {}).get("evidence", [])]
    used = sorted({e for r in out["experiments"] if r["status"] == "completed" for e in r["evidence_ids"]})
    out["evidence"] = {"used": used, "not_used": sorted(set(all_ev) - set(used))}
    # integration-ready records: provenance only, no confidence value
    out["conditional_evidence"] = [
        dict(r["conditional_evidence"], experiment_id=r["experiment_id"], kind=r["kind"],
             model_sha256=out["hashes"]["model_sha256"], config_sha256=out["hashes"]["config_sha256"],
             solver_calls=r["solver_calls"])
        for r in out["experiments"] if r["experiment_id"] != "baseline"]
    out["promotion"] = "none"
    return out


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Evidence-guided computational metabolic experiments (SBML + LP)")
    p.add_argument("--manifest", type=Path, required=True, help="frozen manifest JSON (model path relative to it)")
    p.add_argument("--out", type=Path, default=None, help="output JSON (stdout if omitted)")
    p.add_argument("--strategy", choices=["fixed", "random", "evidence", "qwen"], default=None)
    p.add_argument("--seed", type=int, default=None)
    p.add_argument("--max-experiments", type=int, default=None)
    p.add_argument("--max-solver-calls", type=int, default=None)
    p.add_argument("--qwen-model", default=None)
    p.add_argument("--api-base", default=None, help="OpenAI-compatible base URL (default: GSPA_LLM_URL)")
    args = p.parse_args(argv)
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    validate_manifest(manifest)
    cfg_over = {"strategy": args.strategy, "seed": args.seed, "max_experiments": args.max_experiments,
                "max_solver_calls": args.max_solver_calls}
    if args.qwen_model or args.api_base:
        manifest.setdefault("config", {})["qwen"] = {
            **(manifest.get("config", {}).get("qwen") or {}),
            **({"model": args.qwen_model} if args.qwen_model else {}),
            **({"api_base": args.api_base} if args.api_base else {})}
    model_path = Path(manifest["model"])
    if not model_path.is_absolute():
        model_path = args.manifest.resolve().parent / model_path
    cfg = resolve_config(manifest, cfg_over)
    model = load_sbml_model(model_path, cfg["default_bound"])
    result = run_experiments(manifest, model, cfg_over)
    text = json.dumps(result, indent=2, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text, encoding="utf-8")
    else:
        print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
