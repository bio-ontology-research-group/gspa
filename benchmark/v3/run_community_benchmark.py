#!/usr/bin/env python3
"""Empirical community benchmark: Wintermute & Silver 2010 Keio auxotroph pairs on the production engine.

Question tested: does a *static balanced-growth* model of two gene-deletion strains that share a
medium pool predict which of the 1,035 measured pairs grow together in M9 glucose, where neither
strain grows alone?  The model is iML1515 (development transportability baseline; the authors used
iAF1260).  This is phenotype validation of model-conditional feasibility.  It is NOT a 4-day fold
change, NOT competitive fitness, and NOT gene or function recovery.

Protocol (``PROTOCOL`` below) is fixed in code before any outcome is read and hashed into every
report.  Oxygen supply 0, 2 and 20 mmol/gDW/h are three equal arms; none is selected by outcome.
Labels are used only by ``evaluate``; ``simulate_pair`` receives ``inputs`` alone.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
import resource
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import community_metabolism as CM  # noqa: E402
import metabolic_experiments as ME  # noqa: E402

SCHEMA = "gspa.community_benchmark/1"
PROTOCOL_SCHEMA = "gspa.community_benchmark_protocol/1"

GLUCOSE_EXCHANGE = "R_EX_glc__D_e"
OXYGEN_EXCHANGE = "R_EX_o2_e"
WATER_EXCHANGE = "R_EX_h2o_e"
PROTON_EXCHANGE = "R_EX_h_e"
CO2_EXCHANGE = "R_EX_co2_e"
MAINTENANCE_REACTION = "R_ATPM"

# Fixed before any outcome is read.  Every value is a modelling assumption unless noted.
PROTOCOL: dict[str, Any] = {
    "schema": PROTOCOL_SCHEMA,
    "dataset": "Wintermute & Silver 2010 (10.1038/msb.2010.66), 1,035 unordered pairs of 46 Keio auxotrophs",
    "model": "iML1515 (exact acquired file; hash taken from the case inputs); authors modelled iAF1260",
    "model_role": "development transportability baseline, not the authors' model",
    "glucose_uptake_mmol_gdw_h": 10.0,
    "glucose_note": "normalisation assumption (as in the auxotroph protocols); NOT inferred from 20 g/l",
    "medium_components_from_input": "printed M9 exchanges other than glucose supplied at the ceiling",
    "ceiling_mmol_gdw_h": 1000.0,
    "ceiling_applies_to": ["ammonium", "phosphate", "other printed inorganic ions", "H2O", "H+",
                           "assumed trace metals (arm: trace present)"],
    "trace_metals": "taken from inputs.medium.assumed_not_in_paper.trace_exchanges; assumed required, not measured",
    "oxygen_arms_mmol_gdw_h": [0.0, 2.0, 20.0],
    "oxygen_selection": "all three are reported; none is chosen by outcome",
    "trace_arms": {"present": "primary", "absent": "sensitivity"},
    "organic_supplements": "none beyond glucose",
    "maintenance": f"model-native {MAINTENANCE_REACTION} lower bound, read from the model, not tuned",
    "biomass": "the model's single active-objective biomass reaction, normalised to 1 gDW",
    "community": "balanced growth, common growth rate mu, equal biomass fractions 0.5/0.5, objective max mu",
    "monocultures": "same pool, ports, supply and removal, fraction 1, identical gene deletion",
    "pool": "each model-native extracellular species (unique SBML species id) is one pool metabolite; "
            "every native external exchange of every organism is a port with support 'predicted'",
    "port_support": "predicted (model-derived, not experiment-curated)",
    "removal": "every pool metabolite may be removed up to the ceiling (waste accumulates in the medium); "
               "removal never supplies anything",
    "internal_bounds": "model transport and reaction bounds kept; non-port single-metabolite reactions get no uptake",
    "gene_deletion": "locus tag -> model gene product -> Model.blocked_by_genes (GPR-aware); no guessed genes",
    "growth_threshold_mu": 1e-6,
    "control": "same pair and deletions with uptake from the pool forbidden for every metabolite that is not a "
               "medium nutrient or one of H2O, H+, CO2, O2 (no cross-organism transfer of other metabolites)",
    "prediction": "positive iff community grows AND both matched monocultures do not grow; "
                  "one or two viable monocultures are flagged separately and are not positive",
    "infeasibility": "infeasible under fixed maintenance = predicted no growth (primary); "
                     "strict variant treats it as unknown; solver errors and exhausted budget are unknown",
    "labels": "paper thresholds on mean day-4 fold growth: synergy > 50, no_synergy < 8, intermediate kept separate",
    "scope": "feasibility of mutual growth in static balanced growth only; not fold change or fitness",
}
PROTOCOL_SHA256 = hashlib.sha256(json.dumps(PROTOCOL, sort_keys=True, separators=(",", ":")).encode()).hexdigest()

ARMS = [{"id": f"o2_{int(o2):02d}_trace_{tr}", "o2": o2, "trace": tr == "present"}
        for o2 in PROTOCOL["oxygen_arms_mmol_gdw_h"] for tr in ("present", "absent")]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def canonical_digest(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


# --------------------------------------------------------------------------- context (model + organisms)


class Context:
    """Parsed model, pool, ports and gene-deletion variants.  Built once; label-free."""

    def __init__(self, model_path: Path, loci: dict[str, str]):
        self.cfg = {**CM.DEFAULT_CONFIG, "growth_threshold": PROTOCOL["growth_threshold_mu"]}
        self.model = m = ME.load_sbml_model(model_path, self.cfg["default_bound"])
        self.model_path = model_path
        ext = [r for r in m.rxn_ids if m.is_single_metabolite(r) and m.met_compartment[m.exchange_metabolite(r)] == "e"]
        by_species: dict[str, list[str]] = {}
        for r in ext:
            by_species.setdefault(m.exchange_metabolite(r), []).append(r)
        dup = {s: rs for s, rs in by_species.items() if len(rs) > 1}
        if dup:
            raise ME.ManifestError(f"extracellular species with several exchanges (pool mapping not unique): {dup}")
        self.exchanges = ext
        self.species_of = {r: m.exchange_metabolite(r) for r in ext}
        self.pool_mets = {s: ({"formula": m.met_formula[s]} if s in m.met_formula else {}) for s in by_species}
        self.ports = {r: {"pool": self.species_of[r], "support": "predicted", "metabolite": self.species_of[r]}
                      for r in ext}
        bio_id, _ = ME.resolve_objective(m, {"objectives": {"biomass": None}})
        self.biomass = bio_id
        if MAINTENANCE_REACTION not in m.rxn_index:
            raise ME.ManifestError(f"{MAINTENANCE_REACTION} not in model")
        atpm = float(m.lb[m.rxn_index[MAINTENANCE_REACTION]])
        if not atpm > 0:
            raise ME.ManifestError("native maintenance lower bound is not positive")
        self.maintenance = {"reaction": MAINTENANCE_REACTION, "lower_bound": atpm}
        self.strains: dict[str, dict[str, Any]] = {}
        for strain, locus in sorted(loci.items()):
            gene = m.gene_alias.get(locus)
            blocked = sorted(m.blocked_by_genes({gene})) if gene else []
            self.strains[strain] = {"locus_tag": locus, "gene_product": gene, "blocked_reactions": blocked,
                                    "mapped": gene is not None}
        self._organisms: dict[str, CM.Organism] = {}

    def organism(self, strain: str) -> CM.Organism:
        if strain not in self._organisms:
            info = self.strains[strain]
            if not info["mapped"]:
                raise ME.ManifestError(f"{strain}: locus {info['locus_tag']} does not map to a model gene")
            spec = {"id": strain, "biomass": self.biomass, "maintenance": self.maintenance,
                    "ports": self.ports, "knockouts": info["blocked_reactions"]}
            self._organisms[strain] = CM.build_organism(strain, self.model, spec, self.pool_mets, self.cfg)
        return self._organisms[strain]

    def audit(self) -> dict[str, Any]:
        m = self.model
        return {"model_sha256": m.sha256, "n_reactions": len(m.rxn_ids), "n_metabolites": len(m.met_ids),
                "n_external_exchanges": len(self.exchanges), "biomass_reaction": self.biomass,
                "maintenance": self.maintenance, "n_infinite_bounds_replaced": m.n_infinite_bounds,
                "pool_metabolites": len(self.pool_mets),
                "strains": self.strains,
                "n_strains_mapped": sum(s["mapped"] for s in self.strains.values()),
                "n_strains_blocking_no_reaction": sorted(k for k, s in self.strains.items()
                                                         if s["mapped"] and not s["blocked_reactions"])}


def arm_conditions(ctx: Context, medium: dict[str, Any], arm: dict[str, Any]) -> dict[str, Any]:
    """Supply, removal and control nutrient set for one arm, from the case medium and PROTOCOL only."""
    ceiling = PROTOCOL["ceiling_mmol_gdw_h"]
    supply_ex: dict[str, float] = {}
    for comp in medium["printed_components"]:
        supply_ex[comp["exchange"]] = (PROTOCOL["glucose_uptake_mmol_gdw_h"] if comp["exchange"] == GLUCOSE_EXCHANGE
                                       else ceiling)
    if GLUCOSE_EXCHANGE not in supply_ex:
        raise ME.ManifestError("medium lacks glucose")
    supply_ex[WATER_EXCHANGE] = supply_ex[PROTON_EXCHANGE] = ceiling
    if arm["trace"]:
        for ex in medium["assumed_not_in_paper"]["trace_exchanges"]:
            supply_ex[ex] = ceiling
    supply_ex[OXYGEN_EXCHANGE] = arm["o2"]
    missing = sorted(e for e in supply_ex if e not in ctx.species_of)
    if missing:
        raise ME.ManifestError(f"medium exchanges absent from model: {missing}")
    supply = {ctx.species_of[e]: float(v) for e, v in supply_ex.items() if v > 0}
    removal = {s: ceiling for s in ctx.pool_mets}
    nutrient_ex = {e for e, v in supply_ex.items() if v > 0} | {WATER_EXCHANGE, PROTON_EXCHANGE, CO2_EXCHANGE,
                                                                OXYGEN_EXCHANGE}
    return {"supply": supply, "removal": removal, "nutrient_exchanges": sorted(nutrient_ex),
            "supply_exchanges": dict(sorted(supply_ex.items()))}


# --------------------------------------------------------------------------- simulation (inputs only)


class Simulator:
    """Solves matched monocultures, the community and the transfer-forbidden control with caching.

    Caches key on a hash of the exact condition (model, deletions as blocked reactions, supply, removal,
    closures, fractions, threshold), so a hit returns the result the LP would give.  LP semantics are
    those of ``community_metabolism.build_lp`` unchanged.
    """

    def __init__(self, ctx: Context, conditions: dict[str, Any], max_solver_calls: int):
        self.ctx, self.cond = ctx, conditions
        self.solver = ME.BudgetedSolver(None, max_solver_calls, self.ctx.cfg["residual_tol"])
        self.cache: dict[str, dict[str, Any]] = {}
        self.hits = Counter()
        self.misses = Counter()
        self.lp_seconds = 0.0
        self.build_seconds = 0.0

    def _key(self, strains: list[str], kind: str) -> str:
        kos = sorted(self.ctx.strains[s]["blocked_reactions"] for s in strains)
        return CM_hash({"model": self.ctx.model.sha256, "deletions": kos, "kind": kind,
                        "supply": self.cond["supply"], "removal": self.cond["removal"],
                        "fractions": "equal", "threshold": self.ctx.cfg["growth_threshold"],
                        "maintenance": self.ctx.maintenance})

    def solve(self, strains: list[str], kind: str) -> dict[str, Any]:
        key = self._key(strains, kind)
        if key in self.cache and self.cache[key]["status"] in ("completed",):
            self.hits[kind] += 1
            return {**self.cache[key], "cached": True, "condition_sha256": key}
        self.misses[kind] += 1
        t0 = time.perf_counter()
        orgs = [self.ctx.organism(s) for s in strains]
        closures = []
        if kind == "control":
            for o in orgs:
                closures += [(o.id, rid, "uptake") for rid in o.ports
                             if rid not in self.cond["nutrient_exchanges"]]
        lp = CM.build_lp(orgs, {o.id: 1.0 for o in orgs}, self.cond["supply"], self.cond["removal"], closures)
        t1 = time.perf_counter()
        res = CM.solve_common_growth(self.solver, lp, f"{kind}:{'+'.join(strains)}",
                                     self.ctx.cfg["growth_threshold"])
        self.build_seconds += t1 - t0
        self.lp_seconds += time.perf_counter() - t1
        self.cache[key] = res
        return {**res, "cached": False, "condition_sha256": key}


def CM_hash(obj: Any) -> str:
    return ME.canonical_hash(obj)


def simulate_pair(sim: Simulator, inputs: dict[str, Any]) -> dict[str, Any]:
    """Run the four matched simulations for one pair.  ``inputs`` is the only case data used."""
    a, b = (m["strain"] for m in inputs["members"])
    calls0 = sim.solver.used
    out = {"mono_a": sim.solve([a], "monoculture"), "mono_b": sim.solve([b], "monoculture"),
           "community": sim.solve([a, b], "community"), "control": sim.solve([a, b], "control")}
    out["new_solver_calls"] = sim.solver.used - calls0
    return out


def growth_state(rec: dict[str, Any], strict: bool) -> str:
    if rec["status"] != "completed":
        return "unknown"
    if rec["outcome"] == "infeasible_with_maintenance":
        return "unknown" if strict else "no_growth"
    return rec["outcome"]


def predict(sim_result: dict[str, Any], strict: bool = False) -> dict[str, Any]:
    st = {k: growth_state(sim_result[k], strict) for k in ("mono_a", "mono_b", "community", "control")}
    if "unknown" in (st["community"], st["mono_a"], st["mono_b"]):
        return {"prediction": "unknown", "basis": "a required simulation did not give a growth state", "states": st}
    viable = (st["mono_a"] == "growth") + (st["mono_b"] == "growth")
    if st["community"] == "no_growth":
        return {"prediction": "negative", "basis": "no_community_growth", "states": st, "viable_monocultures": viable}
    if viable == 0:
        return {"prediction": "positive", "basis": "growth_requires_partner_in_model", "states": st,
                "viable_monocultures": 0, "control_grows": st["control"] == "growth"}
    return {"prediction": "negative", "basis": "growth_without_partner" if viable == 2 else "one_viable_monoculture",
            "states": st, "viable_monocultures": viable}


# --------------------------------------------------------------------------- evaluation (labels only here)


def wilson(k: int, n: int, z: float = 1.96) -> list[float] | None:
    if n == 0:
        return None
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return [max(0.0, c - h), min(1.0, c + h)]


def _ratio(k: int, n: int) -> dict[str, Any]:
    return {"k": k, "n": n, "value": k / n if n else None, "wilson95": wilson(k, n)}


def evaluate(rows: list[dict[str, Any]], key: str) -> dict[str, Any]:
    """Confusion of the registered labels against ``row[key]`` (positive/negative/unknown)."""
    cells = Counter((r["observed"], r[key]) for r in rows)
    tp, fn = cells["synergy", "positive"], cells["synergy", "negative"]
    fp, tn = cells["no_synergy", "positive"], cells["no_synergy", "negative"]
    denom = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    n_syn = sum(r["observed"] == "synergy" for r in rows)
    n_no = sum(r["observed"] == "no_synergy" for r in rows)
    return {
        "n_cases": len(rows),
        "observed_counts": dict(Counter(r["observed"] for r in rows)),
        "predicted_counts": dict(Counter(r[key] for r in rows)),
        "confusion_observed_by_predicted": {o: {p: cells[o, p] for p in ("positive", "negative", "unknown")}
                                            for o in ("synergy", "intermediate", "no_synergy")},
        "binary_synergy_vs_no_synergy": {
            "tp": tp, "fn": fn, "fp": fp, "tn": tn,
            "unknown_synergy": cells["synergy", "unknown"], "unknown_no_synergy": cells["no_synergy", "unknown"],
            "precision": _ratio(tp, tp + fp), "recall": _ratio(tp, tp + fn), "specificity": _ratio(tn, tn + fp),
            "recall_counting_unknown_as_miss": _ratio(tp, n_syn),
            "specificity_counting_unknown_as_miss": _ratio(tn, n_no),
            "balanced_accuracy": ((tp / (tp + fn) + tn / (tn + fp)) / 2) if (tp + fn) and (tn + fp) else None,
            "mcc": ((tp * tn - fp * fn) / denom) if denom else None,
            "base_rate_synergy": n_syn / (n_syn + n_no) if n_syn + n_no else None,
        },
        "intermediate_predicted_positive": _ratio(cells["intermediate", "positive"],
                                                  sum(cells["intermediate", p] for p in ("positive", "negative"))),
        "observed_synergy_rate_given_predicted": {
            p: _ratio(sum(r["observed"] == "synergy" for r in rows if r[key] == p),
                      sum(r[key] == p for r in rows)) for p in ("positive", "negative")},
    }


def stratum_rows(rows: list[dict[str, Any]], pred: str = "prediction_primary") -> dict[str, list[dict[str, Any]]]:
    return {
        "all_cases": rows,
        "both_monocultures_not_growing": [r for r in rows if r.get("viable_monocultures") == 0],
        "one_viable_monoculture": [r for r in rows if r.get("viable_monocultures") == 1],
        "both_monocultures_viable": [r for r in rows if r.get("viable_monocultures") == 2],
        "failed_or_unknown_monoculture_state": [r for r in rows if r.get("viable_monocultures") is None],
    }


def arm_summary(arm: dict[str, Any], rows: list[dict[str, Any]], sim: Simulator | None,
                mono_by_strain: dict[str, Any]) -> dict[str, Any]:
    ok = [r for r in rows if r["execution_status"] == "completed"]
    def count(f):
        return sum(1 for r in ok if f(r))
    summary = {
        "arm": arm, "n_cases": len(rows), "n_executed": len(ok), "n_failed": len(rows) - len(ok),
        "coverage": {
            "n_registered_cases": len(rows),
            "n_with_prediction_primary": count(lambda r: r["prediction_primary"] != "unknown"),
            "n_with_prediction_strict": count(lambda r: r["prediction_strict"] != "unknown"),
            "n_unmapped_or_failed": len(rows) - len(ok),
        },
        "monoculture_growth": {
            "strains_growing_alone_in_model": sorted(s for s, v in mono_by_strain.items() if v == "growth"),
            "n_strains": len(mono_by_strain),
            "pairs_both_not_growing": count(lambda r: r["viable_monocultures"] == 0),
            "pairs_one_viable": count(lambda r: r["viable_monocultures"] == 1),
            "pairs_both_viable": count(lambda r: r["viable_monocultures"] == 2),
        },
        "model_outcomes": {
            "community_growth": count(lambda r: r["community"]["outcome"] == "growth"),
            "community_infeasible_with_maintenance": count(
                lambda r: r["community"]["outcome"] == "infeasible_with_maintenance"),
            "community_solver_problem": count(lambda r: r["community"]["status"] != "completed"),
            "control_growth_while_both_monocultures_not_growing": count(
                lambda r: r["viable_monocultures"] == 0 and r["control"]["outcome"] == "growth"),
            "positive_with_control_growth": count(lambda r: r["prediction_primary"] == "positive"
                                                  and r["control"]["outcome"] == "growth"),
        },
        "evaluation_primary": {name: evaluate(rs, "prediction_primary") for name, rs in stratum_rows(rows).items()},
        "evaluation_strict_optimal_only": evaluate(rows, "prediction_strict"),
    }
    if sim is not None:
        summary["solver"] = {"max_calls": sim.solver.max_calls, "used": sim.solver.used,
                             "cache_hits": dict(sim.hits), "cache_misses": dict(sim.misses),
                             "lp_build_seconds": round(sim.build_seconds, 2), "lp_solve_seconds": round(sim.lp_seconds, 2),
                             "n_infeasible": sum(e["status"] == "infeasible" for e in sim.solver.log),
                             "n_error": sum(e["status"] == "error" for e in sim.solver.log)}
    return summary


# --------------------------------------------------------------------------- orchestration


def load_cases(path: Path) -> list[dict[str, Any]]:
    recs = [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    ids = [r["inputs"]["pair_id"] for r in recs]
    if not recs or len(ids) != len(set(ids)):
        raise ValueError("pair identifiers must be nonempty and unique")
    return recs


def strain_loci(recs: list[dict[str, Any]]) -> dict[str, str]:
    loci: dict[str, str] = {}
    for r in recs:
        for m in r["inputs"]["members"]:
            if loci.setdefault(m["strain"], m["locus_tag"]) != m["locus_tag"]:
                raise ValueError(f"strain {m['strain']} maps to several loci")
    return loci


def git_commit() -> str | None:
    try:
        return subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True, text=True,
                              timeout=10, check=True).stdout.strip()
    except Exception:
        return None


def library_versions() -> dict[str, str]:
    out = {"python": platform.python_version()}
    for name in ("numpy", "scipy", "libsbml", "highspy"):
        try:
            mod = __import__(name)
            out[name] = str(getattr(mod, "__version__", None) or
                            (mod.getLibSBMLDottedVersion() if name == "libsbml" else "unknown"))
        except Exception:
            out[name] = "unavailable"
    return out


def run_arm(ctx: Context, recs: list[dict[str, Any]], arm: dict[str, Any], max_solver_calls: int,
            out: Path, log) -> dict[str, Any]:
    first_medium = recs[0]["inputs"]["medium"]
    cond = arm_conditions(ctx, first_medium, arm)
    sim = Simulator(ctx, cond, max_solver_calls)
    rows: list[dict[str, Any]] = []
    t0 = time.perf_counter()
    for i, rec in enumerate(recs):
        inputs, label = rec["inputs"], rec["observed"]
        row: dict[str, Any] = {"case_id": inputs["pair_id"], "arm": arm["id"],
                               "strain_a": inputs["members"][0]["strain"], "strain_b": inputs["members"][1]["strain"],
                               "locus_a": inputs["members"][0]["locus_tag"], "locus_b": inputs["members"][1]["locus_tag"]}
        try:
            if canonical_digest(inputs) != rec.get("inputs_sha256"):
                raise ValueError("inputs digest mismatch")
            if inputs["model_sha256"] != ctx.model.sha256:
                raise ValueError("model hash differs from the registered case")
            if inputs["medium"] != first_medium:
                raise ValueError("medium differs between cases (one medium per run)")
            res = simulate_pair(sim, inputs)
            p, s = predict(res), predict(res, strict=True)
            row.update(execution_status="completed", prediction_primary=p["prediction"], prediction_basis=p["basis"],
                       prediction_strict=s["prediction"], viable_monocultures=p.get("viable_monocultures"),
                       **{k: {f: res[k][f] for f in ("status", "outcome", "mu", "cached", "condition_sha256", "note")}
                          for k in ("mono_a", "mono_b", "community", "control")},
                       new_solver_calls=res["new_solver_calls"])
        except Exception as exc:  # every registered case stays in the denominator
            empty = {"status": "failed", "outcome": "not_run", "mu": None, "cached": False,
                     "condition_sha256": None, "note": ""}
            row.update(execution_status="failed", error=f"{type(exc).__name__}: {exc}",
                       prediction_primary="unknown", prediction_basis="execution_failed", prediction_strict="unknown",
                       viable_monocultures=None, mono_a=empty, mono_b=empty, community=empty, control=empty,
                       new_solver_calls=0)
        rows.append(row)
        if (i + 1) % 100 == 0 or i + 1 == len(recs):
            log(f"[{arm['id']}] {i + 1}/{len(recs)} solver calls {sim.solver.used}/{sim.solver.max_calls} "
                f"({time.perf_counter() - t0:.0f}s)")
    # labels join here, after every simulation of this arm has finished
    for row, rec in zip(rows, recs):
        row["observed"] = rec["observed"]["label"]
        row["observed_day4_fold_growth_mean"] = rec["observed"].get("day4_fold_growth_mean")
    mono_by_strain = {}
    for r in rows:
        if r["execution_status"] == "completed":
            mono_by_strain[r["strain_a"]] = growth_state(r["mono_a"], strict=False)
            mono_by_strain[r["strain_b"]] = growth_state(r["mono_b"], strict=False)
    summary = arm_summary(arm, rows, sim, mono_by_strain)
    summary["conditions"] = {k: cond[k] for k in ("supply_exchanges", "nutrient_exchanges")}
    summary["monoculture_by_strain"] = dict(sorted(mono_by_strain.items()))
    summary["wall_seconds"] = round(time.perf_counter() - t0, 1)
    with open(out / f"cases_{arm['id']}.jsonl", "w", encoding="utf-8") as fh:
        for r in rows:
            fh.write(json.dumps(r) + "\n")
    (out / f"summary_{arm['id']}.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary


def run_benchmark(cases_path: Path, model_path: Path, out: Path, arm_ids: list[str] | None = None,
                  max_solver_calls: int | None = None, max_cases: int | None = None, log=print) -> dict[str, Any]:
    recs = load_cases(cases_path)
    n_registered = len(recs)
    if max_cases is not None:
        recs = recs[:max_cases]
    if out.exists():
        raise ValueError(f"refusing existing run directory: {out}")
    out.mkdir(parents=True)
    arms = [a for a in ARMS if arm_ids is None or a["id"] in arm_ids]
    if arm_ids and len(arms) != len(set(arm_ids)):
        raise ValueError(f"unknown arm in {arm_ids}; known: {[a['id'] for a in ARMS]}")
    t_start = time.perf_counter()
    ctx = Context(model_path, strain_loci(load_cases(cases_path)))
    audit = ctx.audit()
    (out / "model_audit.json").write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    (out / "protocol.json").write_text(json.dumps({"protocol": PROTOCOL, "sha256": PROTOCOL_SHA256}, indent=2) + "\n")
    # one LP per distinct condition: monoculture strains + community and control per pair, plus headroom
    per_arm_budget = max_solver_calls if max_solver_calls is not None else 2 * len(recs) + len(ctx.strains) + 100
    summaries = [run_arm(ctx, recs, a, per_arm_budget, out, log) for a in arms]
    complete = max_cases is None and max_solver_calls is None and arm_ids is None
    report = {
        "schema": SCHEMA,
        "role": "static balanced-growth feasibility of measured Keio auxotroph pairs; not fold change, competitive "
                "fitness, or gene/function recovery; iML1515 transportability baseline for an iAF1260 study",
        "status": "complete" if complete else "partial_run",
        "partial_reasons": [r for r, c in (("max_cases", max_cases is not None), ("max_solver_calls", max_solver_calls is not None),
                                           ("arm_subset", arm_ids is not None)) if c],
        "n_registered_cases": n_registered, "n_cases_run": len(recs),
        "inputs": {"cases_sha256": sha256_file(cases_path), "model_sha256": ctx.model.sha256,
                   "model_file_sha256": sha256_file(model_path), "protocol_sha256": PROTOCOL_SHA256,
                   "runner_sha256": sha256_file(Path(__file__)),
                   "engine_sha256": sha256_file(Path(ME.__file__)),
                   "community_module_sha256": sha256_file(Path(CM.__file__)), "git_commit": git_commit()},
        "protocol": PROTOCOL,
        "model_audit_summary": {k: audit[k] for k in ("n_reactions", "n_metabolites", "n_external_exchanges",
                                                       "biomass_reaction", "maintenance", "n_strains_mapped",
                                                       "n_strains_blocking_no_reaction")},
        "arms": {s["arm"]["id"]: s for s in summaries},
        "resources": {"wall_seconds": round(time.perf_counter() - t_start, 1),
                      "peak_rss_mb": round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 1),
                      "cpu_count": os.cpu_count(), "platform": platform.platform(),
                      "libraries": library_versions(),
                      "slurm": {k: v for k, v in os.environ.items() if k in (
                          "SLURM_JOB_ID", "SLURM_CPUS_PER_TASK", "SLURM_MEM_PER_NODE", "SLURM_JOB_NODELIST")}},
        "limitations": [
            "iML1515 is a transportability baseline; the authors used iAF1260 (not run here)",
            "static balanced growth cannot reproduce 4-day fold change, lag, density or competitive fitness",
            "Trace metals, oxygen, glucose uptake and unlimited waste removal are assumptions, not measurements",
            "labels are measured but the 'intermediate' class (41%) is not binary and is reported separately",
            "port support is model-predicted, not experiment-curated; a predicted rescue assigns no gene or function",
            "development set only: protocol fixed before outcomes, no condition tuned to labels, no held-out split",
        ],
    }
    (out / "report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return report


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cases", type=Path, required=True, help="wintermute2010_cases.jsonl")
    ap.add_argument("--model", type=Path, required=True, help="iML1515 SBML; hash must equal the case inputs")
    ap.add_argument("--out", type=Path, required=True, help="new directory")
    ap.add_argument("--arms", nargs="*", help=f"subset of {[a['id'] for a in ARMS]} (default all; marks run partial)")
    ap.add_argument("--max-solver-calls", type=int, help="per-arm cap (default: planned need); lower marks run partial")
    ap.add_argument("--max-cases", type=int, help="smoke test only; marks run partial")
    args = ap.parse_args()
    rep = run_benchmark(args.cases, args.model, args.out, args.arms, args.max_solver_calls, args.max_cases,
                        log=lambda m: print(m, flush=True))
    for aid, s in rep["arms"].items():
        b = s["evaluation_primary"]["all_cases"]["binary_synergy_vs_no_synergy"]
        print(aid, "executed", s["n_executed"], "/", s["n_cases"], "TP/FP/TN/FN", b["tp"], b["fp"], b["tn"], b["fn"])
    return 0 if all(s["n_failed"] == 0 for s in rep["arms"].values()) else 1


if __name__ == "__main__":
    raise SystemExit(main())
