#!/usr/bin/env python3
"""Independent COBRApy/GLPK check of the v3 community benchmark (shared-pool pairs).

Reproducible oracle, not a claim from scratch: it rebuilds the balanced-growth shared-pool formulation of
``scripts/community_metabolism.py`` with COBRApy objects only (no import of the production parser, LP builder,
GPR evaluator or solver), keeps the SBML reaction/metabolite/gene identifiers, and re-solves a fixed,
declared set of diagnostic pairs and controls. Measured labels are never read: the only case data used are
``inputs`` (members, locus tags, medium).

It checks implementation agreement (same LP semantics, same numbers within declared tolerances) and attributes
a rescue to the transfer and to the model's own transport reactions. It does NOT show that a transporter
exists, that the biology is right, or anything about the 1,035-pair benchmark beyond the selected pairs.

Inputs
  --model        iML1515 SBML (hash must equal the case inputs)
  --cases        wintermute2010_cases.jsonl (full file or a subset made with --write-subset-cases)
  --selection    selection manifest (pairs, controls, arms, protocol values, tolerances; declared before outcomes)
  --out          new oracle report (refuses to overwrite)
Production side (optional; validation of completed revision-2 results)
  --production-report   report.json of a full or bounded production run (its directory holds cases_<arm>.jsonl)
  --arm-cases           one or more per-arm cases_<arm>.jsonl files (no report; provenance is then unverified)
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import sys
from pathlib import Path
from typing import Any

import cobra
from cobra import Metabolite, Model, Reaction

SCHEMA = "gspa.community_cobra_oracle/1"
SELECTION_SCHEMA = "gspa.community_oracle_selection/1"
SIMS = ("mono_a", "mono_b", "community", "control")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def canonical(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


# --------------------------------------------------------------------------- selection manifest


def load_selection(path: Path) -> dict[str, Any]:
    sel = json.loads(path.read_text(encoding="utf-8"))
    if sel.get("schema") != SELECTION_SCHEMA:
        raise ValueError(f"selection schema must be {SELECTION_SCHEMA}")
    decl = sel.get("declaration", {})
    if decl.get("role") != "development" or decl.get("labels_used") is not False:
        raise ValueError("selection must declare role 'development' and labels_used false")
    need = {"glucose_exchange", "glucose_uptake", "ceiling", "water_exchange", "proton_exchange", "co2_exchange",
            "oxygen_exchange", "maintenance_reaction", "biomass_reaction", "growth_threshold", "oxygen_arms",
            "trace_arms"}
    miss = need - set(sel.get("protocol", {}))
    if miss:
        raise ValueError(f"selection.protocol lacks {sorted(miss)}")
    tol = sel.get("tolerances", {})
    if not {"mu_rel", "mu_abs", "cross_flux"} <= set(tol) or not all(tol[k] > 0 for k in ("mu_rel", "mu_abs", "cross_flux")):
        raise ValueError("selection.tolerances needs positive mu_rel, mu_abs and cross_flux")
    pairs = sel.get("pairs", [])
    ids = [(p["a"], p["b"]) for p in pairs]
    if not pairs or len(set(map(frozenset, ids))) != len(ids) or any(a == b for a, b in ids):
        raise ValueError("selection.pairs must be non-empty, distinct unordered pairs of different strains")
    if not any(p["role"] == "viable_monoculture_control" for p in pairs):
        raise ValueError("selection needs a viable-monoculture pair control")
    return sel


def arms_of(sel: dict[str, Any]) -> list[dict[str, Any]]:
    p = sel["protocol"]
    return [{"id": f"o2_{int(o2):02d}_trace_{'present' if tr else 'absent'}", "o2": float(o2), "trace": tr}
            for o2 in p["oxygen_arms"] for tr in ([True, False] if p["trace_arms"] == "both" else [True])]


def arm_conditions(medium: dict[str, Any], arm: dict[str, Any], proto: dict[str, Any]) -> dict[str, Any]:
    """Exchange-id -> uptake rate from the case medium and the declared protocol values only."""
    supply = {c["exchange"]: (proto["glucose_uptake"] if c["exchange"] == proto["glucose_exchange"] else proto["ceiling"])
              for c in medium["printed_components"]}
    if proto["glucose_exchange"] not in supply:
        raise ValueError("medium lacks glucose")
    supply[proto["water_exchange"]] = supply[proto["proton_exchange"]] = proto["ceiling"]
    if arm["trace"]:
        for ex in medium["assumed_not_in_paper"]["trace_exchanges"]:
            supply[ex] = proto["ceiling"]
    supply[proto["oxygen_exchange"]] = arm["o2"]
    supply = {k: float(v) for k, v in supply.items() if v > 0}
    nutrients = set(supply) | {proto["water_exchange"], proto["proton_exchange"], proto["co2_exchange"],
                               proto["oxygen_exchange"]}
    return {"supply": supply, "nutrients": nutrients}


# --------------------------------------------------------------------------- COBRA formulation


class Oracle:
    def __init__(self, model_path: Path, proto: dict[str, Any], tol: dict[str, float]):
        cobra.Configuration().solver = "glpk"
        # f_replace={} keeps the SBML identifiers (R_EX_..., M_..._e, G_b....) untouched
        self.base = cobra.io.read_sbml_model(str(model_path), f_replace={})
        self.proto, self.tol = proto, tol
        b = self.base
        self.exchanges = [r for r in b.reactions if len(r.metabolites) == 1
                          and next(iter(r.metabolites)).compartment == "e"]
        if len({next(iter(r.metabolites)).id for r in self.exchanges}) != len(self.exchanges):
            raise ValueError("an extracellular species has several exchanges")
        self.single = {r.id for r in b.reactions if len(r.metabolites) == 1}
        self.port_ids = {r.id for r in self.exchanges}
        objective = [r.id for r in b.reactions if r.objective_coefficient]
        if objective != [proto["biomass_reaction"]]:
            raise ValueError(f"active objective {objective} is not the declared biomass reaction")
        self.native_maintenance = b.reactions.get_by_id(proto["maintenance_reaction"]).lower_bound
        self._ko: dict[str, list[str]] = {}
        self._models: dict[tuple, Model] = {}

    def knocked_reactions(self, locus: str) -> list[str]:
        """Reactions COBRApy's GPR evaluation blocks when the gene product of ``locus`` is removed."""
        if locus not in self._ko:
            hit = [g for g in self.base.genes if g.id in (locus, "G_" + locus) or g.name == locus]
            if len(hit) != 1:
                raise ValueError(f"locus {locus}: {len(hit)} matching genes (no guessing)")
            native = {r.id: r.bounds for r in self.base.reactions}
            with self.base:
                hit[0].knock_out()
                self._ko[locus] = sorted(r.id for r in self.base.reactions
                                         if r.bounds == (0.0, 0.0) and native[r.id] != (0.0, 0.0))
        return self._ko[locus]

    def native_bounds(self, locus: str | None) -> dict[str, tuple[float, float]]:
        """Organism bounds before any port opening: sinks/demands without uptake, deletion applied last."""
        out: dict[str, tuple[float, float]] = {}
        for r in self.base.reactions:
            lb, ub = r.bounds
            if r.id in self.single and r.id not in self.port_ids:      # sinks / demands: no uptake
                c = next(iter(r.metabolites.values()))
                lb, ub = (max(lb, 0.0), ub) if c < 0 else (lb, min(ub, 0.0))
            if r.id == self.proto["maintenance_reaction"]:
                lb = max(lb, self.native_maintenance)
            out[r.id] = (lb, ub)
        for rid in (self.knocked_reactions(locus) if locus else []):
            out[rid] = (0.0, 0.0)
        return out

    def pool_model(self, members: list[tuple[str, str | None]]) -> Model:
        """Block-diagonal organism copies + one pool metabolite per extracellular exchange + common growth.

        Built once per member list with the model's native exchange bounds; ``configure`` varies bounds only."""
        key = tuple(members)
        if key in self._models:
            return self._models[key]
        frac = 1.0 / len(members)
        comp = Model("pool")
        pool = {next(iter(r.metabolites)).id: Metabolite("pool::" + next(iter(r.metabolites)).id, compartment="pool")
                for r in self.exchanges}
        bio_rxns = []
        for strain, locus in members:
            bounds = self.native_bounds(locus)
            mets = {m.id: Metabolite(f"{strain}::{m.id}", formula=m.formula, compartment=m.compartment)
                    for m in self.base.metabolites}
            new = [Reaction(f"{strain}::{r.id}", lower_bound=bounds[r.id][0], upper_bound=bounds[r.id][1])
                   for r in self.base.reactions]
            comp.add_reactions(new)
            for r, nr in zip(self.base.reactions, new):
                stoich = {mets[m.id]: c for m, c in r.metabolites.items()}
                if r.id in self.port_ids:
                    met, c = next(iter(r.metabolites.items()))
                    stoich[pool[met.id]] = -c * frac               # secretion into the pool is positive flux
                nr.add_metabolites(stoich)
            bio_rxns.append(comp.reactions.get_by_id(f"{strain}::{self.proto['biomass_reaction']}"))
        mu = Reaction("mu", lower_bound=0.0, upper_bound=1000.0)
        comp.add_reactions([mu])
        mu_stoich = {}
        for (strain, _), bio in zip(members, bio_rxns):
            g = Metabolite(f"growth::{strain}")
            bio.add_metabolites({g: 1.0})
            mu_stoich[g] = -1.0
        mu.add_metabolites(mu_stoich)
        env = []
        for r in self.exchanges:
            sid = next(iter(r.metabolites)).id
            env += [Reaction(f"ENV::supply::{sid}", lower_bound=0.0, upper_bound=0.0),
                    Reaction(f"ENV::removal::{sid}", lower_bound=0.0, upper_bound=float(self.proto["ceiling"]))]
        comp.add_reactions(env)
        for r in self.exchanges:
            sid = next(iter(r.metabolites)).id
            comp.reactions.get_by_id(f"ENV::supply::{sid}").add_metabolites({pool[sid]: 1.0})
            comp.reactions.get_by_id(f"ENV::removal::{sid}").add_metabolites({pool[sid]: -1.0})
        comp.objective = mu
        comp.objective_direction = "max"
        comp.solver = "glpk"
        self._models[key] = comp
        return comp

    def configure(self, comp: Model, members: list[tuple[str, str | None]], supply_ex: dict[str, float],
                  open_rate: float | None, forbid: set[str] = frozenset(),
                  blocked: dict[str, set[str]] | None = None) -> None:
        """Set supply, port opening, uptake closures and extra blocked reactions (call inside ``with comp``).

        Order as in the formulation: open uptake (at least ``open_rate``; a more permissive native bound is kept),
        then closures win, then blocked reactions."""
        for r in self.exchanges:
            sid = next(iter(r.metabolites)).id
            comp.reactions.get_by_id(f"ENV::supply::{sid}").upper_bound = float(supply_ex.get(r.id, 0.0))
            c = next(iter(r.metabolites.values()))
            for strain, _ in members:
                rx = comp.reactions.get_by_id(f"{strain}::{r.id}")
                lb, ub = r.bounds
                if open_rate is not None:
                    lb, ub = (min(lb, -open_rate), ub) if c < 0 else (lb, max(ub, open_rate))
                if r.id in forbid:
                    lb, ub = (max(lb, 0.0), ub) if c < 0 else (lb, min(ub, 0.0))
                rx.bounds = (lb, ub)
        for strain, rids in (blocked or {}).items():
            for rid in rids:
                comp.reactions.get_by_id(f"{strain}::{rid}").bounds = (0.0, 0.0)

    def run(self, members, supply_ex, open_rate, forbid=frozenset(), blocked=None, fluxes=False) -> dict[str, Any]:
        comp = self.pool_model(members)
        with comp:
            self.configure(comp, members, supply_ex, open_rate, forbid, blocked)
            return self.solve(comp, with_fluxes=fluxes)

    def solve(self, comp: Model, with_fluxes: bool = False) -> dict[str, Any]:
        sol = comp.optimize(raise_error=False)
        if sol.status == "infeasible":
            res = {"outcome": "infeasible_with_maintenance", "mu": None}
        elif sol.status != "optimal":
            res = {"outcome": f"solver_{sol.status}", "mu": None}
        else:
            mu = max(0.0, float(sol.objective_value))
            res = {"outcome": "growth" if mu >= self.proto["growth_threshold"] else "no_growth", "mu": mu}
        if with_fluxes and sol.status == "optimal":
            res["_fluxes"] = sol.fluxes
        return res

    # ---- simulations of one pair under one arm
    def simulate_pair(self, members: list[tuple[str, str]], cond: dict[str, Any]) -> dict[str, Any]:
        rate = float(self.proto["ceiling"])
        supply = cond["supply"]
        forbid = {r.id for r in self.exchanges if r.id not in cond["nutrients"]}
        out = {"mono_a": self.run(members[:1], supply, rate), "mono_b": self.run(members[1:], supply, rate),
               "community": self.run(members, supply, rate, fluxes=True),
               "control": self.run(members, supply, rate, forbid)}
        flux = out["community"].pop("_fluxes", None)
        # same pair with the model's native exchange bounds (legacy port semantics, no open_uptake)
        out["native_bounds_community"] = self.run(members, supply, None)
        out["transfers"] = self.attribute(members, supply, rate, flux, out["community"]) \
            if flux is not None and out["community"]["outcome"] == "growth" else None
        return out

    def attribute(self, members, supply, rate, flux, community) -> dict[str, Any]:
        """Cross-fed species of one optimal community solution, and which of them growth depends on.

        For each species the consumer's reactions that touch the extracellular species (its porins / transporters,
        i.e. the model's internal transport bounds, kept as in the model) are closed. If growth is then lost the
        rescue rests on that model transport step in addition to the opened exchange bound: opening the exchange does
        not bypass membrane transport. Species whose closure leaves growth are only counted. Alternative optima
        exist; the candidate list comes from one optimal solution."""
        tol = self.tol["cross_flux"]
        frac = 1.0 / len(members)
        net: dict[str, dict[str, float]] = {}
        for strain, _ in members:
            for r in self.exchanges:
                met, c = next(iter(r.metabolites.items()))
                net.setdefault(met.id, {})[strain] = -c * frac * float(flux[f"{strain}::{r.id}"])
        rows, tested = [], 0
        for sid, d in sorted(net.items()):
            ex = next(r for r in self.exchanges if next(iter(r.metabolites)).id == sid)
            if supply.get(ex.id, 0.0) > 0:
                continue                                           # medium nutrient, not a partner product
            prod = {o: v for o, v in d.items() if v > tol}
            cons = {o: -v for o, v in d.items() if v < -tol}
            if not prod or not cons or not any(p != c for p in prod for c in cons):
                continue
            entry, closed_mu = [], {}
            for consumer in cons:
                path = sorted((r for r in self.base.metabolites.get_by_id(sid).reactions if r.id != ex.id),
                              key=lambda r: r.id)
                entry.append({"consumer": consumer, "reactions": [
                    {"id": r.id, "bounds": list(r.bounds), "gpr": r.gene_reaction_rule} for r in path]})
                closed_mu[consumer] = self.run(members, supply, rate, blocked={consumer: {r.id for r in path}})
            tested += 1
            if all(v["outcome"] == "growth" for v in closed_mu.values()):
                continue                                           # dispensable given the other transfers
            rows.append({"species": sid, "producers": prod, "consumers": cons, "consumer_entry": entry,
                         "growth_with_consumer_entry_reactions_closed": closed_mu,
                         "exchange_native_bounds": list(ex.bounds)})
        return {"candidates_tested": tested, "necessary_given_other_transfers": rows,
                "n_dispensable": tested - len(rows)}

    def simulate_single(self, strain: str, locus: str | None, cond: dict[str, Any], supplement: tuple[str, float] | None):
        supply = dict(cond["supply"])
        if supplement:
            supply[supplement[0]] = supplement[1]
        return self.run([(strain, locus)], supply, float(self.proto["ceiling"]))


# --------------------------------------------------------------------------- comparison with production rows


def agree(oracle: dict[str, Any], prod: dict[str, Any], tol: dict[str, float]) -> dict[str, Any]:
    same = oracle["outcome"] == prod.get("outcome")
    mu_o, mu_p = oracle["mu"], prod.get("mu")
    if mu_o is None or mu_p is None:
        close, diff = mu_o is None and mu_p is None, None
    else:
        diff = abs(mu_o - mu_p)
        close = math.isclose(mu_o, mu_p, rel_tol=tol["mu_rel"], abs_tol=tol["mu_abs"])
    return {"outcome_agrees": same, "mu_agrees": close, "abs_diff": diff, "agrees": same and close,
            "production": {"outcome": prod.get("outcome"), "mu": mu_p}}


def load_production(args, cases_sha: str, model_file_sha: str, sel: dict[str, Any]):
    """Return ({arm: {pair frozenset: row}}, provenance). Only rows of selected pairs are read into comparison."""
    files: dict[str, Path] = {}
    prov: dict[str, Any] = {}
    if args.production_report:
        rep_path = args.production_report
        rep = json.loads(rep_path.read_text(encoding="utf-8"))
        inp = rep["inputs"]
        if inp["cases_sha256"] != cases_sha:
            raise ValueError("cases file differs from the production run")
        if inp["model_file_sha256"] != model_file_sha:
            raise ValueError("model file differs from the production run")
        proto, pp = sel["protocol"], rep["protocol"]
        pairs = {"glucose_uptake": pp["glucose_uptake_mmol_gdw_h"], "ceiling": pp["ceiling_mmol_gdw_h"],
                 "growth_threshold": pp["growth_threshold_mu"]}
        bad = {k: (proto[k], v) for k, v in pairs.items() if proto[k] != v}
        if not set(proto["oxygen_arms"]) <= set(pp["oxygen_arms_mmol_gdw_h"]):
            bad["oxygen_arms"] = (proto["oxygen_arms"], pp["oxygen_arms_mmol_gdw_h"])
        if bad:
            raise ValueError(f"declared protocol differs from the production protocol: {bad}")
        prov = {"report_sha256": sha256(rep_path), "report_status": rep["status"], "protocol_revision": pp.get("protocol_revision"),
                "protocol_sha256": inp["protocol_sha256"], "n_registered_cases": rep["n_registered_cases"],
                "n_cases_run": rep["n_cases_run"], "partial_reasons": rep["partial_reasons"],
                "production_hashes": {k: inp[k] for k in ("runner_sha256", "engine_sha256", "community_module_sha256",
                                                          "git_commit")}}
        for arm in rep["arms"]:
            files[arm] = rep_path.parent / f"cases_{arm}.jsonl"
    for f in args.arm_cases or []:
        first = json.loads(f.read_text(encoding="utf-8").splitlines()[0])
        files[first["arm"]] = f
    rows: dict[str, dict[frozenset, dict[str, Any]]] = {}
    for arm, f in files.items():
        rows[arm] = {frozenset((r["strain_a"], r["strain_b"])): r
                     for r in map(json.loads, f.read_text(encoding="utf-8").splitlines()) if r.get("arm") == arm}
    prov["cases_files"] = {arm: {"file": f.name, "sha256": sha256(f)} for arm, f in sorted(files.items())}
    prov["provenance"] = "verified_against_production_report" if args.production_report else "per_arm_files_only"
    return rows, prov


# --------------------------------------------------------------------------- driver


def run(args) -> dict[str, Any]:
    sel = load_selection(args.selection)
    proto, tol = sel["protocol"], sel["tolerances"]
    model_file_sha, cases_sha = sha256(args.model), sha256(args.cases)
    records = [json.loads(l) for l in args.cases.read_text(encoding="utf-8").splitlines() if l.strip()]
    by_pair = {frozenset(m["strain"] for m in r["inputs"]["members"]): r for r in records}
    loci = {m["strain"]: m["locus_tag"] for r in records for m in r["inputs"]["members"]}
    oracle = Oracle(args.model, proto, tol)
    prod_rows, prov = ({}, None)
    if args.production_report or args.arm_cases:
        prod_rows, prov = load_production(args, cases_sha, model_file_sha, sel)
    arms = arms_of(sel)
    report_arms: dict[str, Any] = {}
    n_cmp = n_agree = 0
    max_abs = 0.0
    for arm in arms:
        pairs_out, controls_out = [], []
        medium = None
        for p in sel["pairs"]:
            rec = by_pair.get(frozenset((p["a"], p["b"])))
            if rec is None:
                raise ValueError(f"selected pair {p['a']}+{p['b']} is not in the cases file")
            inputs = rec["inputs"]                      # labels (rec['observed']) are never touched
            if inputs["model_sha256"] != model_file_sha:
                raise ValueError("model hash differs from the registered case")
            medium = medium or inputs["medium"]
            if inputs["medium"] != medium:
                raise ValueError("medium differs between selected cases")
            cond = arm_conditions(medium, arm, proto)
            members = [(m["strain"], m["locus_tag"]) for m in inputs["members"]]
            res = oracle.simulate_pair(members, cond)
            row = {"pair_id": inputs["pair_id"], "role": p["role"], "strains": [m[0] for m in members],
                   "oracle": {k: res[k] for k in SIMS + ("native_bounds_community",)},
                   "transfers": res["transfers"]}
            both_dead = res["mono_a"]["outcome"] != "growth" and res["mono_b"]["outcome"] != "growth"
            row["rescue_in_oracle"] = bool(both_dead and res["community"]["outcome"] == "growth")
            prod = prod_rows.get(arm["id"], {}).get(frozenset((p["a"], p["b"])))
            if prod_rows:
                if prod is None:
                    row["production_comparison"] = {"status": "missing_row", "agrees": False}
                    n_cmp += 1
                else:
                    cmp = {k: agree(res[k], prod[k], tol) for k in SIMS}
                    ok = all(c["agrees"] for c in cmp.values())
                    row["production_comparison"] = {"status": "compared", "sims": cmp, "agrees": ok}
                    max_abs = max([max_abs] + [c["abs_diff"] for c in cmp.values() if c["abs_diff"] is not None])
                    n_cmp += 1
                n_agree += row["production_comparison"]["agrees"]
            pairs_out.append(row)
        for c in sel.get("single_controls", []):
            locus = loci.get(c["strain"])
            if locus is None:
                raise ValueError(f"control strain {c['strain']} has no locus in the cases file")
            cond = arm_conditions(medium, arm, proto)
            controls_out.append({
                "strain": c["strain"], "locus_tag": locus, "supplement": c["supplement_exchange"],
                "wild_type": oracle.simulate_single(c["strain"], None, cond, None),
                "deletion": oracle.simulate_single(c["strain"], locus, cond, None),
                "deletion_plus_supplement": oracle.simulate_single(
                    c["strain"], locus, cond, (c["supplement_exchange"], c["supplement_rate"]))})
        report_arms[arm["id"]] = {"arm": arm, "pairs": pairs_out, "single_controls": controls_out}
    comp = oracle.base
    solver_cfg = comp.solver.configuration
    return {
        "schema": SCHEMA,
        "interpretation": "Independent implementation agreement and mechanism attribution on a declared development "
                          "selection; not phenotype accuracy, not transporter evidence, not the full benchmark",
        "status": "completed",
        "hashes": {"model_file_sha256": model_file_sha, "cases_sha256": cases_sha,
                   "selection_sha256": sha256(args.selection), "oracle_script_sha256": sha256(Path(__file__)),
                   "selection_canonical_sha256": canonical(sel)},
        "production": prov,
        "versions": {"python": platform.python_version(), "cobra": cobra.__version__,
                     "optlang": __import__("optlang").__version__, "solver": "glpk",
                     "swiglpk": getattr(__import__("swiglpk"), "__version__", "unknown"),
                     "glpk": _glpk_version()},
        "solver_parameters": {"interface": comp.solver.interface.__name__, "presolve": str(solver_cfg.presolve),
                              "lp_method": str(getattr(solver_cfg, "lp_method", None)),
                              "timeout": solver_cfg.timeout, "cobra_tolerance": cobra.Configuration().tolerance,
                              "note": "library defaults; nothing tuned"},
        "tolerances": tol,
        "declaration": sel["declaration"],
        "protocol": proto,
        "arms": report_arms,
        "agreement": None if not prod_rows else {
            "n_pair_arm_compared": n_cmp, "n_agree": n_agree, "max_abs_mu_diff": max_abs,
            "all_agree": n_cmp == n_agree and n_cmp > 0},
        "limits": ["one optimal solution per community LP; alternative optima can change which species look cross-fed",
                   "growth classification near the threshold would be solver-sensitive; mu values are reported",
                   "agreement between two implementations is not evidence that a transporter exists"],
    }


def _glpk_version() -> str:
    try:
        import swiglpk
        return str(swiglpk.glp_version())
    except Exception:
        return "unknown"


def write_subset(args) -> int:
    sel = load_selection(args.selection)
    wanted = [frozenset((p["a"], p["b"])) for p in sel["pairs"]]
    lines = {}
    for line in args.cases.read_text(encoding="utf-8").splitlines():
        if line.strip():
            rec = json.loads(line)
            lines[frozenset(m["strain"] for m in rec["inputs"]["members"])] = line     # untouched text
    missing = [sorted(w) for w in wanted if w not in lines]
    if missing:
        raise ValueError(f"selected pairs missing from cases: {missing}")
    if args.write_subset_cases.exists():
        raise ValueError("refusing to overwrite subset cases file")
    args.write_subset_cases.write_text("\n".join(lines[w] for w in wanted) + "\n", encoding="utf-8")
    print(f"wrote {len(wanted)} cases, sha256 {sha256(args.write_subset_cases)}")
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", type=Path)
    ap.add_argument("--cases", type=Path, required=True)
    ap.add_argument("--selection", type=Path, required=True)
    ap.add_argument("--out", type=Path)
    ap.add_argument("--production-report", type=Path)
    ap.add_argument("--arm-cases", type=Path, nargs="+")
    ap.add_argument("--write-subset-cases", type=Path, help="write the selected pairs' case lines and exit")
    args = ap.parse_args(argv)
    if args.write_subset_cases:
        return write_subset(args)
    if not args.model or not args.out:
        ap.error("--model and --out are required")
    if args.out.exists():
        raise ValueError("refusing to overwrite oracle report")
    report = run(args)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    ag = report["agreement"]
    print(json.dumps({"status": report["status"], "agreement": ag, "versions": report["versions"]}, indent=2))
    return 0 if ag is None or ag["all_agree"] else 1


if __name__ == "__main__":
    sys.exit(main())
