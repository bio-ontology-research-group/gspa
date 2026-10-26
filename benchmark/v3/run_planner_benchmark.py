#!/usr/bin/env python3
"""Equal-budget condition-planner benchmark on the production metabolic experiment engine.

Question: given one gene and a label-free menu of carbon-source conditions, which planner
finds a condition where the model predicts the deletion abolishes growth while the matched
wild-type control grows, within a fixed experiment / LP-call budget?

Two things are reported and never mixed:

* Experiment A (model-derived condition discovery, DEVELOPMENT). Success is defined by the
  same model's exhaustive sweep. It says whether a planner finds model-predicted conditional
  essentiality cheaply. It is circular with respect to function: it does not measure whether any
  gene function was recovered.
* Experiment B (measured phenotype comparison on visited conditions). For (gene, substrate)
  pairs that a planner happened to visit and that have a measured Keio growth label, the model's
  knockout prediction is compared with the measurement. The 61 measured pairs were selected by a
  model-based auxotrophy design (iML1515), so this is precision of a model-predicted auxotrophy
  on the visited subset, not recall, and not independent function recovery.

Truth independence. Menus, evidence, hypotheses and prompts are built from the gene list, the
substrate panel and the model only (``public_view``). Measured labels are joined after every run.
A regression test flips every label and checks that menus, evidence and Qwen prompts are
byte-identical.

All arms drive the production engine (``scripts/metabolic_experiments.py``): fixed / random /
evidence arms through ``run_experiments`` strategies (paired arms as an explicit ``panel``),
Qwen through the engine's own selection loop (invalid-choice counting, fallback, budget) with
only the prompt text and chat transport replaced.
"""
from __future__ import annotations

import argparse
import collections
import contextlib
import copy
import gzip
import hashlib
import json
import math
import os
import random
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any, Callable
from unittest import mock

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(HERE))
import metabolic_experiments as ME  # noqa: E402
import run_metabolic_benchmark as B  # noqa: E402

SCHEMA = "gspa.v3.planner_benchmark/1"

# Everything below is frozen before any planner is run; protocol.json records it with the input hashes.
FROZEN: dict[str, Any] = {
    "protocol": "manuscript_regulatory",
    "budgets_experiments": [8, 16],             # engine max_experiments; LP cap = budget + 1 (baseline) for every arm
    "lp_calls_per_experiment": 1,               # max_focus_reactions = 0: growth LP only, no FVA, no product LPs
    "random_seeds": list(range(10)),
    "qwen": {"model": "qwen3.8-27b", "seed": 0, "max_invalid": 2, "max_tokens": 300, "timeout_s": 600},
    "evidence": {"currency_degree_gt": 25, "near_hops_le": 4, "reachable_hops_le": 10},
    "bootstrap": {"seed": 20261026, "n": 10000},
    "arms": ["fixed_paired", "random_item", "random_paired", "evidence_item", "evidence_paired", "qwen"],
}
FORBIDDEN_EVIDENCE_TYPES = {"bgc", "product"}


def sha256_file(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


# --------------------------------------------------------------------------- label-free inputs


def public_view(records: list[dict]) -> dict[str, Any]:
    """The only facts the planners may use from the case file: loci, substrate panel, closures, base manifest.

    Pairings of locus and substrate are discarded on purpose: which pairs were measured is a property of the
    model-based design, and the full locus x substrate grid is offered to every arm.
    """
    loci = sorted({r["case"]["locus_tag"] for r in records})
    substrates = sorted({r["case"]["substrate_exchange"] for r in records})
    constraints: dict[str, dict] = {}
    base = None
    for r in records:
        manifest = B.execution_manifest(r)
        for con in manifest.get("condition_constraints", []):
            prior = constraints.setdefault(con["substrate_exchange"], con)
            if prior != con:
                raise ValueError(f"substrate {con['substrate_exchange']} has differing closures across cases")
        trimmed = {k: v for k, v in manifest.items() if k not in ("context", "config", "condition_constraints", "benchmark")}
        if base is None:
            base = trimmed
        elif trimmed != base:
            raise ValueError("cases differ in baseline, objectives or model: one shared base is required")
    with_constraints = sum(1 for r in records if B.execution_manifest(r).get("condition_constraints"))
    if with_constraints not in (0, len(records)) or (with_constraints and not set(substrates) <= set(constraints)):
        raise ValueError("a case or substrate lacks a registered closure set")
    return {"loci": loci, "substrates": substrates, "constraints": [constraints[s] for s in substrates if s in constraints], "base": base}


class Network:
    """Undirected metabolite-reaction graph used only for label-free proximity evidence (no LP)."""

    def __init__(self, model: ME.Model, currency_degree_gt: int):
        S = model.S.tocsc()
        self.model = model
        self.rxn_mets = [S.indices[S.indptr[j]:S.indptr[j + 1]].tolist() for j in range(S.shape[1])]
        mets_rxns: dict[int, list[int]] = collections.defaultdict(list)
        for j, ms in enumerate(self.rxn_mets):
            for m in ms:
                mets_rxns[m].append(j)
        self.met_rxns = mets_rxns
        self.currency = {m for m, rs in mets_rxns.items() if len(rs) > currency_degree_gt}
        self._cache: dict[int, dict[int, int]] = {}

    def hops_from(self, met: int) -> dict[int, int]:
        if met not in self._cache:
            dist, frontier, seen_r = {met: 0}, [met], set()
            while frontier:
                nxt = []
                for m in frontier:
                    for r in self.met_rxns[m]:
                        if r in seen_r:
                            continue
                        seen_r.add(r)
                        for m2 in self.rxn_mets[r]:
                            if m2 not in dist and m2 not in self.currency:
                                dist[m2] = dist[m] + 1
                                nxt.append(m2)
                frontier = nxt
            self._cache[met] = dist
        return self._cache[met]

    def distance(self, substrate_exchange: str, reactions: list[str]) -> int | None:
        model = self.model
        src = model.met_ids.index(model.exchange_metabolite(substrate_exchange))
        dist = self.hops_from(src)
        targets = {m for rid in reactions for m in self.rxn_mets[model.rxn_index[rid]] if m not in self.currency}
        hops = [dist[m] for m in targets if m in dist]
        return min(hops) if hops else None


def hypothesis_id(locus: str, substrate: str) -> str:
    return f"H:{locus}:{substrate}"


def gene_manifest(view: dict, locus: str, network: Network, model: ME.Model, budget: int,
                  strategy: str = "fixed", seed: int = 0, panel: list[str] | None = None) -> dict:
    """Label-free manifest for one gene: 24 candidate conditions, each with a control and a knockout item."""
    ev_cfg = FROZEN["evidence"]
    reactions = sorted(ME.hypothesis_targets(model, {"gene": locus})[0])
    hyps, evidence = [], []
    for s in view["substrates"]:
        hid = hypothesis_id(locus, s)
        hyps.append({"id": hid, "gene": locus, "condition": {"carbon_exchange": s}})
        hops = network.distance(s, reactions) if reactions else None
        common = {"hypothesis_id": hid, "condition_hint": "carbon_swap", "exchange": s}
        evidence.append({"id": f"E:{locus}:{s}:candidate", "type": "transport", **common})
        if hops is not None and hops <= ev_cfg["reachable_hops_le"]:
            evidence.append({"id": f"E:{locus}:{s}:reachable", "type": "network_reachable", "hops": hops, **common})
        if hops is not None and hops <= ev_cfg["near_hops_le"]:
            evidence.append({"id": f"E:{locus}:{s}:near", "type": "network_near", "hops": hops, **common})
    manifest = copy.deepcopy(view["base"])
    manifest["context"] = {"hypotheses": hyps, "evidence": evidence}
    manifest["condition_constraints"] = copy.deepcopy(view["constraints"])
    manifest["config"] = {
        "strategy": strategy, "seed": seed, "max_experiments": budget,
        "max_solver_calls": budget + 1, "uptake_rate": 10.0,
        "max_focus_reactions": 0, "max_baseline_focus_reactions": 0,
        "panel": panel, "qwen": dict(FROZEN["qwen"]),
    }
    assert_no_forced_objective(manifest)
    return manifest


def assert_no_forced_objective(manifest: dict) -> None:
    """Planner menus must not contain product objectives or forced-flux checks (circular as function evidence)."""
    cfg = manifest.get("config", {})
    bad_types = {e.get("type") for e in manifest["context"]["evidence"]} & FORBIDDEN_EVIDENCE_TYPES
    if bad_types or any(e.get("product_exchange") for e in manifest["context"]["evidence"]):
        raise ValueError(f"forced-production evidence is not admissible in planner menus: {sorted(bad_types)}")
    if cfg.get("max_focus_reactions") != 0 or cfg.get("max_baseline_focus_reactions") != 0:
        raise ValueError("focus checks impose a growth floor (a forced flux) and are disabled in this benchmark")


def condition_ids(manifest: dict) -> list[tuple[str, str, str, str]]:
    """(substrate, control item id, knockout item id, hypothesis id), alphabetical by substrate exchange."""
    return [(h["condition"]["carbon_exchange"], f"carbon_source:{h['condition']['carbon_exchange']}",
             f"knockout:{h['id']}@{h['condition']['carbon_exchange']}", h["id"])
            for h in sorted(manifest["context"]["hypotheses"], key=lambda h: h["condition"]["carbon_exchange"])]


def paired_panel(order: list[tuple], budget: int) -> list[str]:
    panel = []
    for _, ctl, ko, _ in order:
        panel += [ctl, ko]
    return panel[:budget]


def arm_panel(arm: str, manifest: dict, budget: int, seed: int) -> list[str] | None:
    conds = condition_ids(manifest)
    evidence = manifest["context"]["evidence"]
    if arm == "fixed_paired":
        return paired_panel(conds, budget)
    if arm == "random_paired":
        shuffled = list(conds)
        random.Random(seed).shuffle(shuffled)
        return paired_panel(shuffled, budget)
    if arm == "evidence_paired":
        # same evidence_key the production evidence strategy uses, applied to the knockout item of each condition
        ev_ids = {}
        for e in evidence:
            ev_ids.setdefault(e["hypothesis_id"], []).append(e["id"])
        def key(c):
            return ME.evidence_key({"evidence_ids": ev_ids[c[3]], "competing_hypotheses": [], "id": c[2]}, evidence)
        return paired_panel(sorted(conds, key=key), budget)
    return None


ENGINE_STRATEGY = {"fixed_paired": "fixed", "random_paired": "fixed", "evidence_paired": "fixed",
                   "random_item": "random", "evidence_item": "evidence", "qwen": "qwen"}


# --------------------------------------------------------------------------- Qwen transport and prompt


class UsageChat:
    """OpenAI-compatible chat call identical in payload to ME.make_qwen_chat, recording latency and token usage.

    The key comes from ME.load_api_key() (environment or the user's key file) and is only placed in the request
    header; it is never stored, logged or echoed. Error records keep the exception class and HTTP status only.
    """

    def __init__(self, cfg: dict, base: str | None = None):
        q = cfg["qwen"]
        self.cfg, self.q = cfg, q
        self.base = (base or q.get("api_base") or os.environ.get("GSPA_LLM_URL", "http://unimatrix01.kaust.edu.sa:8000/v1")).rstrip("/")
        self._key = ME.load_api_key()
        self.log: list[dict] = []
        self.transcript: list[dict] = []

    def __call__(self, messages: list[dict[str, str]]) -> str:
        q = self.q
        payload = {"model": q["model"], "messages": messages, "temperature": 0.0, "seed": self.cfg["seed"],
                   "max_tokens": q.get("max_tokens", 300), "chat_template_kwargs": {"enable_thinking": False}}
        req = urllib.request.Request(f"{self.base}/chat/completions", data=json.dumps(payload).encode(),
                                     headers={"Content-Type": "application/json", "Authorization": f"Bearer {self._key}"},
                                     method="POST")
        t0, entry = time.perf_counter(), {"ok": False}
        try:
            with urllib.request.urlopen(req, timeout=q.get("timeout_s", 180)) as resp:
                body = json.loads(resp.read().decode())
            text = body["choices"][0]["message"]["content"]
            usage = body.get("usage") or {}
            entry.update(ok=True, prompt_tokens=usage.get("prompt_tokens"), completion_tokens=usage.get("completion_tokens"))
            self.transcript.append({"messages": messages, "response": text})
            return text
        except urllib.error.HTTPError as exc:
            entry["error"] = f"HTTPError {exc.code}"
            raise
        except Exception as exc:
            entry["error"] = type(exc).__name__
            raise
        finally:
            entry["seconds"] = round(time.perf_counter() - t0, 3)
            self.log.append(entry)

    def usage(self) -> dict:
        ok = [e for e in self.log if e["ok"]]
        return {"requests": len(self.log), "ok": len(ok), "failed": len(self.log) - len(ok),
                "errors": dict(collections.Counter(e["error"] for e in self.log if not e["ok"])),
                "prompt_tokens": sum(e.get("prompt_tokens") or 0 for e in ok),
                "completion_tokens": sum(e.get("completion_tokens") or 0 for e in ok),
                "tokens_reported": all(e.get("prompt_tokens") is not None for e in ok) if ok else None,
                "seconds": round(sum(e["seconds"] for e in self.log), 3)}


def evidence_view(manifest: dict) -> dict[str, dict]:
    return {e["id"]: {"type": e["type"], **({"hops": e["hops"]} if "hops" in e else {})}
            for e in manifest["context"]["evidence"]}


def make_prompt(manifest: dict, locus: str) -> Callable:
    """Qwen prompt: the engine's fields plus the content of each cited evidence record and the planning goal.

    The deterministic evidence arm orders by the same records, so no arm sees evidence another cannot.
    """
    ev = evidence_view(manifest)
    limit = manifest["config"]["max_experiments"]

    def prompt(remaining: list[dict], results: list[dict], budget_left: int) -> list[dict[str, str]]:
        menu = [{"experiment_id": m["id"], "kind": m["kind"], "description": m["description"],
                 "evidence": [dict(ev[i], id=i) for i in m["evidence_ids"]],
                 "competing_hypotheses": m["competing_hypotheses"],
                 "min_calls": m["min_calls"], "max_calls": m["max_calls"]} for m in remaining]
        prior = [{"experiment_id": r["experiment_id"], "status": r["status"], "outcome": r["outcome"],
                  "max_growth": r["measurements"].get("max_growth"),
                  "growth_ratio": r["measurements"].get("growth_ratio")} for r in results]
        system = (f"You choose the next in-silico metabolic experiment from a fixed menu. Goal: for the deletion of gene "
                  f"{locus}, find carbon-source conditions in which the model predicts the deletion abolishes growth while "
                  "the matched wild-type control (the carbon_source item of the same exchange) grows. A knockout result is "
                  "interpretable only together with its matched carbon_source result. Results are model-conditional, not "
                  "proof. Reply with one JSON object: "
                  '{"experiment_id": "<id from menu>", "reason": "<short>"} or {"stop": true, "reason": "<short>"}.')
        user = json.dumps({"admissible_menu": menu, "completed_results": prior, "solver_calls_left": budget_left,
                           "experiments_left": limit - (len(results) - 1)})
        return [{"role": "system", "content": system}, {"role": "user", "content": user}]

    return prompt


def prompt_parity(manifest: dict, menu: list[dict], first_prompt: list[dict[str, str]]) -> dict:
    """Does the first Qwen prompt expose exactly the admissible menu and the evidence the evidence arm orders by?"""
    shown = json.loads(first_prompt[1]["content"])["admissible_menu"]
    ev = evidence_view(manifest)
    ids_ok = [m["experiment_id"] for m in shown] == [m["id"] for m in menu]
    content_ok = all({e["id"]: {k: v for k, v in e.items() if k != "id"} for e in s["evidence"]} ==
                     {i: ev[i] for i in m["evidence_ids"]} for s, m in zip(shown, menu))
    return {"menu_ids_identical": ids_ok, "evidence_content_identical": content_ok}


# --------------------------------------------------------------------------- execution and scoring


def run_arm(model: ME.Model, view: dict, network: Network, locus: str, arm: str, budget: int, seed: int,
            chat: Callable | None = None, base: str | None = None) -> dict:
    strategy = ENGINE_STRATEGY[arm]
    probe = gene_manifest(view, locus, network, model, budget)
    panel = arm_panel(arm, probe, budget, seed)
    manifest = gene_manifest(view, locus, network, model, budget, strategy=strategy, seed=seed, panel=panel)
    usage_chat = None
    t0 = time.perf_counter()
    if arm == "qwen":
        usage_chat = chat or UsageChat(ME.resolve_config(manifest), base)
        first: list = []
        prompt_fn = make_prompt(manifest, locus)

        def logged(remaining, results, left):
            msgs = prompt_fn(remaining, results, left)
            if not first:
                first.append(msgs)
            return msgs

        with mock.patch.object(ME, "qwen_prompt", logged):
            result = ME.run_experiments(manifest, model, chat=usage_chat)
    else:
        result = ME.run_experiments(manifest, model)
    wall = time.perf_counter() - t0
    rec = compact(result, manifest, locus, arm, budget, seed, wall)
    if arm == "qwen":
        rec["qwen"] = qwen_audit(result, usage_chat)
        if first:
            rec["qwen"]["prompt_parity"] = prompt_parity(manifest, result["menu"], first[0])
        rec["_transcript"] = getattr(usage_chat, "transcript", [])
    return rec


def compact(result: dict, manifest: dict, locus: str, arm: str, budget: int, seed: int, wall: float) -> dict:
    source = {e["experiment_id"]: e["source"] for e in result["selection_log"] if "experiment_id" in e}
    items = {}
    for r in result["experiments"]:
        items[r["experiment_id"]] = {
            "kind": r["kind"], "status": r["status"], "outcome": r["outcome"],
            "max_growth": r["measurements"].get("max_growth"), "source": source.get(r["experiment_id"]),
            "solver_calls": len(r["solver_calls"])}
    return {"gene": locus, "arm": arm, "budget": budget, "seed": seed,
            "lp_calls_used": result["solver"]["used"], "lp_calls_cap": result["solver"]["max_calls"],
            "experiments_attempted": sum(1 for r in result["experiments"] if r["experiment_id"] != "baseline"
                                         and r["status"] != "skipped_not_selected"),
            "wall_seconds": round(wall, 4),
            "menu_sha256": ME.canonical_hash(result["menu"]), "context_sha256": result["hashes"]["context_sha256"],
            "growth_threshold": result["config"]["growth_threshold"],
            "selection_log": result["selection_log"], "items": items, "promotion": result.get("promotion"),
            "conditional_evidence_independent": any(c.get("independent_proof_of_function") for c in result["conditional_evidence"])}


def qwen_audit(result: dict, chat: Any) -> dict:
    log = result["selection_log"]
    invalid = [e for e in log if e.get("valid") is False]
    return {"valid_choices": sum(1 for e in log if e.get("source") == "qwen" and e.get("valid") and "experiment_id" in e),
            "invalid_choices": len(invalid),
            "endpoint_errors": sum(1 for e in invalid if str(e.get("error", "")).startswith("endpoint error")),
            "errors": [e["error"] for e in invalid],
            "fallback_activated": any(e.get("source") == "evidence_fallback" and "experiment_id" in e for e in log),
            "fallback_steps": sum(1 for e in log if e.get("source") == "evidence_fallback" and "experiment_id" in e),
            "stopped_early": any(e.get("stop") for e in log),
            "usage": chat.usage() if hasattr(chat, "usage") else None}


def phenotype_of(item: dict | None, threshold: float) -> str:
    if item is None:
        return "unknown"
    return B.phenotype({"status": item["status"], "outcome": item["outcome"],
                        "measurements": {"max_growth": item["max_growth"]}}, threshold)[0]


def discovered(rec: dict, sources: tuple[str, ...] | None = None) -> list[str]:
    """Substrates whose control grew and whose knockout did not, both run by this arm (optionally by given sources)."""
    thr = rec["growth_threshold"]
    out = []
    for item_id, item in rec["items"].items():
        if not item_id.startswith("carbon_source:"):
            continue
        sub = item_id.split(":", 1)[1]
        ko_id = next((k for k in rec["items"] if k.startswith("knockout:") and k.endswith("@" + sub)), None)
        ko = rec["items"].get(ko_id)
        if sources and (item["source"] not in sources or ko is None or ko["source"] not in sources):
            continue
        if phenotype_of(item, thr) == "growth" and phenotype_of(ko, thr) == "no_growth":
            out.append(sub)
    return sorted(out)


def sweep_gene(model: ME.Model, view: dict, network: Network, locus: str) -> dict:
    """Evaluator-side exhaustive run of all 24 control/knockout pairs (not a planner; LPs are not planner budget)."""
    probe = gene_manifest(view, locus, network, model, 1)
    n = len(probe["context"]["hypotheses"]) * 2
    panel = paired_panel(condition_ids(probe), n)
    manifest = gene_manifest(view, locus, network, model, n, strategy="fixed", panel=panel)
    result = ME.run_experiments(manifest, model)
    rec = compact(result, manifest, locus, "sweep", n, 0, 0.0)
    thr = rec["growth_threshold"]
    cond = {}
    for sub, ctl, ko, _ in condition_ids(probe):
        c, k = rec["items"].get(ctl), rec["items"].get(ko)
        cond[sub] = {"control": phenotype_of(c, thr), "knockout": phenotype_of(k, thr),
                     "control_growth": c and c["max_growth"], "knockout_growth": k and k["max_growth"],
                     "control_outcome": c and c["outcome"], "knockout_outcome": k and k["outcome"]}
    essential = sorted(s for s, v in cond.items() if v["control"] == "growth" and v["knockout"] == "no_growth")
    return {"gene": locus, "knockout_reactions": sorted(ME.hypothesis_targets(model, {"gene": locus})[0]),
            "mapped_to_model": locus in model.gene_alias, "lp_calls_used": rec["lp_calls_used"],
            "model_essential_conditions": essential, "conditions": cond,
            "n_incomplete": sum(1 for v in cond.values() if "unknown" in (v["control"], v["knockout"]))}


# --------------------------------------------------------------------------- statistics


def wilson(k: int, n: int, z: float = 1.959964) -> tuple[float | None, float | None]:
    if n == 0:
        return None, None
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(0.0, c - h), min(1.0, c + h)


def bootstrap_mean(values: list[float], n: int, seed: int) -> tuple[float, float, float]:
    rng = random.Random(seed)
    m = len(values)
    means = sorted(sum(values[rng.randrange(m)] for _ in range(m)) / m for _ in range(n))
    return sum(values) / m, means[int(0.025 * n)], means[int(0.975 * n) - 1]


def exact_mcnemar(b: int, c: int) -> float:
    n = b + c
    if n == 0:
        return 1.0
    k = min(b, c)
    return min(1.0, 2 * sum(math.comb(n, i) for i in range(k + 1)) / 2 ** n)


# --------------------------------------------------------------------------- analysis


def measured_index(records: list[dict]) -> dict[tuple[str, str], dict]:
    out = {}
    for r in records:
        c = r["case"]
        out[(c["locus_tag"], c["substrate_exchange"])] = {"case_id": c["case_id"], "label": c.get("label", "unknown"),
                                                          "raw_od_audited": c.get("raw_od_reproduces_published_metrics") is True}
    return out


def measured_visits(rec: dict, measured: dict) -> list[dict]:
    thr = rec["growth_threshold"]
    rows = []
    for item_id, item in rec["items"].items():
        if not item_id.startswith("knockout:"):
            continue
        sub = item_id.rsplit("@", 1)[1]
        meas = measured.get((rec["gene"], sub))
        if meas is None or item["status"] != "completed":
            continue
        ctl = rec["items"].get(f"carbon_source:{sub}")
        rows.append({"gene": rec["gene"], "substrate": sub, "case_id": meas["case_id"], "observed": meas["label"],
                     "prediction": phenotype_of(item, thr), "control": phenotype_of(ctl, thr)})
    return rows


def summarize_measured(rows: list[dict]) -> dict:
    binary = [r for r in rows if r["observed"] in ("growth", "no_growth") and r["prediction"] != "unknown"]
    agree = sum(r["observed"] == r["prediction"] for r in binary)
    pred_aux = [r for r in rows if r["prediction"] == "no_growth"]
    conf = sum(r["observed"] == "no_growth" for r in pred_aux)
    lo, hi = wilson(agree, len(binary))
    plo, phi = wilson(conf, len(pred_aux))
    return {"n_visited_measured_pairs": len(rows),
            "n_with_growing_control": sum(r["control"] == "growth" for r in rows),
            "observed_counts": dict(collections.Counter(r["observed"] for r in rows)),
            "n_binary_evaluable": len(binary), "n_binary_agree": agree,
            "binary_agreement": agree / len(binary) if binary else None, "binary_agreement_ci95": [lo, hi],
            "n_predicted_auxotroph": len(pred_aux), "n_confirmed_no_growth": conf,
            "confirmed_fraction": conf / len(pred_aux) if pred_aux else None, "confirmed_fraction_ci95": [plo, phi]}


def analyse(sweeps: dict[str, dict], runs: list[dict], measured: dict, budgets: list[int]) -> dict:
    eligible = sorted(g for g, s in sweeps.items() if s["model_essential_conditions"])
    out: dict[str, Any] = {"genes": {"total": len(sweeps), "with_model_essential_condition": len(eligible),
                                      "without": sorted(g for g in sweeps if g not in eligible)},
                           "budgets": {}}
    boot = FROZEN["bootstrap"]
    for b in budgets:
        per_arm: dict[str, dict[str, list[float]]] = {}
        stats: dict[str, Any] = {}
        for arm in FROZEN["arms"]:
            rs = [r for r in runs if r["arm"] == arm and r["budget"] == b]
            if not rs:
                continue
            seeds = sorted({r["seed"] for r in rs})
            gene_vals, yield_vals, lps, wall = {}, {}, [], []
            for g in eligible:
                succ, yld = [], []
                for sd in seeds:
                    r = next((x for x in rs if x["gene"] == g and x["seed"] == sd), None)
                    if r is None:
                        continue
                    d = discovered(r, ("qwen",) if arm == "qwen" else None)
                    succ.append(1.0 if d else 0.0)
                    yld.append(len(d) / min(len(sweeps[g]["model_essential_conditions"]), b // 2))
                if succ:
                    gene_vals[g], yield_vals[g] = sum(succ) / len(succ), sum(yld) / len(yld)
            per_arm[arm] = gene_vals
            vals = [gene_vals[g] for g in eligible if g in gene_vals]
            mean, lo, hi = bootstrap_mean(vals, boot["n"], boot["seed"]) if vals else (None, None, None)
            seed_rates = []
            for sd in seeds:
                sv = [1.0 if discovered(next(x for x in rs if x["gene"] == g and x["seed"] == sd), ("qwen",) if arm == "qwen" else None) else 0.0
                      for g in eligible if any(x["gene"] == g and x["seed"] == sd for x in rs)]
                if sv:
                    seed_rates.append(sum(sv) / len(sv))
            st = {"n_runs": len(rs), "n_seeds": len(seeds), "genes_scored": len(vals),
                  "success_rate": mean, "success_ci95_bootstrap_over_genes": [lo, hi],
                  "per_seed_success_rate": seed_rates if len(seeds) > 1 else None,
                  "mean_yield_fraction": sum(yield_vals.values()) / len(yield_vals) if yield_vals else None,
                  "lp_calls_mean": sum(r["lp_calls_used"] for r in rs) / len(rs),
                  "lp_calls_max": max(r["lp_calls_used"] for r in rs), "lp_cap": rs[0]["lp_calls_cap"],
                  "wall_seconds_total": sum(r["wall_seconds"] for r in rs)}
            if arm == "qwen":
                st["qwen"] = qwen_totals(rs, sweeps, eligible)
            # measured comparison for visited conditions (pooled over seeds; per-seed pooled counts reported as mean)
            rows_by_seed = [[row for r in rs if r["seed"] == sd for row in measured_visits(r, measured)] for sd in seeds]
            st["measured_visited"] = summarize_measured([row for rows in rows_by_seed for row in rows]) if len(seeds) == 1 else \
                {"per_seed": [summarize_measured(rows) for rows in rows_by_seed],
                 "mean_distinct_visited_pairs_per_seed": sum(len({(x["gene"], x["substrate"]) for x in rows}) for rows in rows_by_seed) / len(seeds)}
            stats[arm] = st
        comparisons = {}
        if "qwen" in per_arm:
            for arm, vals in per_arm.items():
                if arm == "qwen":
                    continue
                genes = [g for g in eligible if g in vals and g in per_arm["qwen"]]
                diffs = [per_arm["qwen"][g] - vals[g] for g in genes]
                if not diffs:
                    continue
                mean, lo, hi = bootstrap_mean(diffs, boot["n"], boot["seed"] + 1)
                binary = all(v in (0.0, 1.0) for g in genes for v in (vals[g], per_arm["qwen"][g]))
                comp = {"genes": len(genes), "qwen_minus_arm": mean, "ci95_bootstrap_over_genes": [lo, hi]}
                if binary:
                    bq = sum(1 for g in genes if per_arm["qwen"][g] == 1.0 and vals[g] == 0.0)
                    ba = sum(1 for g in genes if per_arm["qwen"][g] == 0.0 and vals[g] == 1.0)
                    comp.update(qwen_only_success=bq, arm_only_success=ba, exact_mcnemar_p=exact_mcnemar(bq, ba))
                comparisons[arm] = comp
        out["budgets"][str(b)] = {"arms": stats, "qwen_vs": comparisons, "per_gene_success": per_arm}
    return out


def qwen_totals(rs: list[dict], sweeps: dict, eligible: list[str]) -> dict:
    q = [r["qwen"] for r in rs]
    all_found = [bool(discovered(r)) for r in rs if r["gene"] in eligible]
    attributable = [bool(discovered(r, ("qwen",))) for r in rs if r["gene"] in eligible]
    clean = [r for r in rs if r["gene"] in eligible and not r["qwen"]["invalid_choices"] and not r["qwen"]["fallback_activated"]]
    usage = [x["usage"] for x in q if x["usage"]]
    return {"genes_run": len(rs), "genes_with_any_invalid_choice": sum(1 for x in q if x["invalid_choices"]),
            "genes_with_fallback": sum(1 for x in q if x["fallback_activated"]),
            "invalid_choices_total": sum(x["invalid_choices"] for x in q),
            "endpoint_errors_total": sum(x["endpoint_errors"] for x in q),
            "fallback_steps_total": sum(x["fallback_steps"] for x in q),
            "stopped_early_genes": sum(1 for x in q if x["stopped_early"]),
            "success_any_source_hybrid": sum(all_found) / len(all_found) if all_found else None,
            "success_qwen_attributable_only": sum(attributable) / len(attributable) if attributable else None,
            "genes_clean_no_invalid_no_fallback": len(clean),
            "success_among_clean_genes": (sum(1 for r in clean if discovered(r)) / len(clean)) if clean else None,
            "api_requests": sum(u["requests"] for u in usage), "api_failed": sum(u["failed"] for u in usage),
            "prompt_tokens": sum(u["prompt_tokens"] for u in usage), "completion_tokens": sum(u["completion_tokens"] for u in usage),
            "api_seconds": round(sum(u["seconds"] for u in usage), 1),
            "prompt_parity_all_identical": all(x.get("prompt_parity", {}).get("menu_ids_identical") and
                                               x.get("prompt_parity", {}).get("evidence_content_identical") for x in q)}


def parity_audit(runs: list[dict]) -> dict:
    """Every arm must have faced the identical menu and context for a gene (and budget, which only changes the cap)."""
    seen: dict[tuple, set] = collections.defaultdict(set)
    for r in runs:
        seen[(r["gene"], r["budget"])].add((r["menu_sha256"], r["context_sha256"]))
    bad = [list(k) for k, v in seen.items() if len(v) != 1]
    caps = {(r["budget"], r["lp_calls_cap"]) for r in runs}
    return {"genes_budgets_checked": len(seen), "mismatched": bad, "lp_caps_by_budget": sorted(caps),
            "equal_budget_everywhere": all(cap == b + 1 for b, cap in caps),
            "over_cap_runs": [[r["arm"], r["gene"], r["budget"]] for r in runs if r["lp_calls_used"] > r["lp_calls_cap"]]}


def measured_accounting(records: list[dict], sweeps: dict) -> list[dict]:
    """Every registered measured case with its eligibility for the model-vs-measurement comparison."""
    rows = []
    for r in records:
        c = r["case"]
        sw = sweeps.get(c["locus_tag"])
        cond = sw["conditions"].get(c["substrate_exchange"]) if sw else None
        if sw is None:
            status, why = "excluded", "gene_not_in_planner_gene_list"
        elif cond is None or "unknown" in (cond["control"], cond["knockout"]):
            status, why = "excluded", "sweep_solve_incomplete"
        elif cond["control"] != "growth":
            status, why = "excluded", f"wild_type_control_{cond['control']}_on_substrate"
        else:
            status, why = "eligible", ""
        rows.append({"case_id": c["case_id"], "locus_tag": c["locus_tag"], "substrate_exchange": c["substrate_exchange"],
                     "observed": c.get("label", "unknown"), "status": status, "reason": why,
                     "note": ("gene deletion blocks no reaction in the model (redundant isozyme); the engine runs it "
                              "as a no-effect knockout" if sw and not sw["knockout_reactions"] else ""),
                     "model_control": cond and cond["control"], "model_knockout": cond and cond["knockout"],
                     "model_knockout_growth": cond and cond["knockout_growth"]})
    return rows


# --------------------------------------------------------------------------- driver


def load_records(path: Path) -> list[dict]:
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def append_jsonl(path: Path, rec: dict) -> None:
    with path.open("a") as fh:
        fh.write(json.dumps(rec, sort_keys=True) + "\n")


def read_jsonl(path: Path) -> list[dict]:
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()] if path.exists() else []


def prepare(args) -> tuple[ME.Model, dict, Network, list[dict]]:
    records = load_records(args.cases)
    view = public_view(records)
    cfg = ME.resolve_config(view["base"])
    model = ME.load_sbml_model(args.model or Path(view["base"]["model"]), cfg["default_bound"])
    if model.sha256 != view["base"].get("model_sha256_expected"):
        raise ValueError("model hash differs from the registered case file")
    return model, view, Network(model, FROZEN["evidence"]["currency_degree_gt"]), records


def write_protocol(args, out: Path, view: dict) -> dict:
    proto = {"schema": SCHEMA, "frozen": FROZEN,
             "inputs": {"cases_sha256": sha256_file(args.cases), "runner_sha256": sha256_file(Path(__file__)),
                        "engine_sha256": sha256_file(Path(ME.__file__)), "model_sha256": view["base"]["model_sha256_expected"]},
             "n_genes": len(view["loci"]), "n_substrates": len(view["substrates"]), "substrates": view["substrates"],
             "roles": {"experiment_A": "model-derived condition discovery (development); not measured function recovery",
                       "experiment_B": "measured phenotype comparison on visited conditions; selected cases, no recall claim"}}
    path = out / "protocol.json"
    if path.exists():
        prior = json.loads(path.read_text())
        if prior["inputs"]["cases_sha256"] != proto["inputs"]["cases_sha256"] or \
                prior["inputs"]["engine_sha256"] != proto["inputs"]["engine_sha256"]:
            raise ValueError("protocol.json was written for different inputs; refusing to continue this run directory")
        if prior["frozen"] != FROZEN:
            # Only the Qwen transport timeout may be revised, and only before any Qwen request of this directory was
            # evaluated; the revision stays in the protocol file.
            strip = lambda f: {**f, "qwen": {k: v for k, v in f["qwen"].items() if k != "timeout_s"}}
            if strip(prior["frozen"]) != strip(FROZEN) or any(path.parent.glob("qwen*.jsonl")):
                raise ValueError("protocol.json differs from the frozen protocol; refusing to continue this run directory")
            prior.setdefault("revisions", []).append(
                {"changed": "qwen.timeout_s", "from": prior["frozen"]["qwen"]["timeout_s"], "to": FROZEN["qwen"]["timeout_s"],
                 "reason": "endpoint latency of 87-115 s per request observed in a smoke run before any evaluated Qwen request",
                 "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())})
            prior["frozen"] = FROZEN
            path.write_text(json.dumps(prior, indent=2) + "\n")
    else:
        path.write_text(json.dumps(proto, indent=2) + "\n")
    return proto


def stage_sweep(args, out, model, view, network, records) -> dict:
    path = out / "sweep.json"
    if path.exists():
        return json.loads(path.read_text())
    sweeps = {}
    for i, locus in enumerate(view["loci"]):
        sweeps[locus] = sweep_gene(model, view, network, locus)
        print(f"sweep {i + 1}/{len(view['loci'])} {locus}: {len(sweeps[locus]['model_essential_conditions'])} model-essential", flush=True)
    path.write_text(json.dumps(sweeps, indent=1, sort_keys=True) + "\n")
    return sweeps


def stage_runs(args, out, model, view, network, arms: list[str], name: str) -> None:
    loci = view["loci"][: args.limit_genes] if args.limit_genes else view["loci"]
    if args.shard:
        index, count = (int(x) for x in args.shard.split("/"))
        loci, name = loci[index::count], f"{name}.shard{index}of{count}"
    path = out / f"{name}.jsonl"
    done = {(r["gene"], r["arm"], r["budget"], r["seed"]) for r in read_jsonl(path)}
    transcripts = out / f"{name}_transcripts.jsonl.gz"
    for budget in args.budgets:
        for arm in arms:
            seeds = FROZEN["random_seeds"] if arm.startswith("random") else [FROZEN["qwen"]["seed"] if arm == "qwen" else 0]
            for seed in seeds:
                for locus in loci:
                    if (locus, arm, budget, seed) in done:
                        continue
                    rec = run_arm(model, view, network, locus, arm, budget, seed, base=args.api_base)
                    transcript = rec.pop("_transcript", None)
                    if transcript:
                        with gzip.open(transcripts, "at") as fh:
                            fh.write(json.dumps({"gene": locus, "budget": budget, "turns": transcript}) + "\n")
                    append_jsonl(path, rec)
                    tag = f" qwen_invalid={rec['qwen']['invalid_choices']}" if arm == "qwen" else ""
                    print(f"{name} B={budget} {arm} seed={seed} {locus}: lp={rec['lp_calls_used']}{tag}", flush=True)


def stage_report(args, out, view, records, sweeps) -> dict:
    runs = read_jsonl(out / "deterministic.jsonl") + [r for p in sorted(out.glob("qwen*.jsonl")) for r in read_jsonl(p)]
    measured = measured_index(records)
    analysis = analyse(sweeps, runs, measured, args.budgets)
    report = {"schema": SCHEMA + "/report",
              "roles": json.loads((out / "protocol.json").read_text())["roles"],
              "parity_audit": parity_audit(runs),
              "analysis": analysis,
              "sweep_reference": {"genes": len(sweeps),
                                  "model_essential_conditions_per_gene": {g: len(s["model_essential_conditions"]) for g, s in sweeps.items()},
                                  "sweep_lp_calls": sum(s["lp_calls_used"] for s in sweeps.values()),
                                  "genes_whose_deletion_blocks_no_reaction": sorted(g for g, s in sweeps.items() if not s["knockout_reactions"]),
                                  "genes_with_incomplete_sweep_solves": sorted(g for g, s in sweeps.items() if s["n_incomplete"])},
              "measured_cases": measured_accounting(records, sweeps)}
    acc = report["measured_cases"]
    report["measured_reference_all_registered"] = {
        "n_registered": len(acc), "n_eligible": sum(a["status"] == "eligible" for a in acc),
        "excluded": [a for a in acc if a["status"] != "eligible"],
        "sweep_model_vs_measured_eligible": summarize_measured([
            {"observed": a["observed"], "prediction": a["model_knockout"], "control": a["model_control"]}
            for a in acc if a["status"] == "eligible"])}
    if args.phenotype_reference:
        report["engine_consistency_with_phenotype_runner"] = consistency(args.phenotype_reference, sweeps, records)
    (out / "report.json").write_text(json.dumps(report, indent=1, sort_keys=True) + "\n")
    return report


def consistency(reference: Path, sweeps: dict, records: list[dict]) -> dict:
    """Sweep growth vs. the registered phenotype runner (same engine, one matched pair per case)."""
    ref = json.loads(reference.read_text())
    exchange = {r["case"]["case_id"]: (r["case"]["locus_tag"], r["case"]["substrate_exchange"]) for r in records}
    diffs, mismatched, missing = [], [], []
    for row in ref["cases"]:
        gene, sub = exchange[row["case_id"]]
        cond = sweeps.get(gene, {}).get("conditions", {}).get(sub)
        if cond is None or row.get("knockout_growth") is None or cond["knockout_growth"] is None:
            missing.append(row["case_id"])
            continue
        diffs.append(abs(cond["knockout_growth"] - row["knockout_growth"]))
        if cond["knockout"] != row["prediction"] or cond["control"] != row["control_prediction"]:
            mismatched.append(row["case_id"])
    return {"reference_sha256": sha256_file(reference), "n_reference_cases": len(ref["cases"]),
            "n_compared": len(diffs), "max_abs_knockout_growth_difference": max(diffs) if diffs else None,
            "phenotype_calls_differing": mismatched, "not_comparable": missing,
            "reference_all_cases": ref["all_cases"]}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cases", type=Path, required=True, help="cases_manuscript_regulatory.jsonl")
    ap.add_argument("--model", type=Path, help="relocated iML1515; must match the registered SHA-256")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--stages", default="sweep,deterministic,report")
    ap.add_argument("--budgets", default=",".join(map(str, FROZEN["budgets_experiments"])))
    ap.add_argument("--limit-genes", type=int, default=0, help="debug only; a limited run is not the registered evaluation")
    ap.add_argument("--api-base", default=None)
    ap.add_argument("--shard", default=None, help="i/n: run genes i, i+n, ... (parallel Qwen workers; results merge in the report)")
    ap.add_argument("--phenotype-reference", type=Path, help="report.json of run_metabolic_benchmark.py on the same cases")
    ap.add_argument("--cpus", type=int, default=2)
    args = ap.parse_args()
    args.budgets = [int(b) for b in args.budgets.split(",")]
    if hasattr(os, "sched_setaffinity"):
        allowed = sorted(os.sched_getaffinity(0))
        os.sched_setaffinity(0, set(allowed[: args.cpus]))
    args.out.mkdir(parents=True, exist_ok=True)
    model, view, network, records = prepare(args)
    write_protocol(args, args.out, view)
    stages = args.stages.split(",")
    sweeps = stage_sweep(args, args.out, model, view, network, records)
    if "deterministic" in stages:
        stage_runs(args, args.out, model, view, network, [a for a in FROZEN["arms"] if a != "qwen"], "deterministic")
    if "qwen" in stages:
        stage_runs(args, args.out, model, view, network, ["qwen"], "qwen")
    if "report" in stages:
        report = stage_report(args, args.out, view, records, sweeps)
        print(json.dumps(report["parity_audit"], indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
