#!/usr/bin/env python3
"""Competitive gene assignment for reactions added by targeted tests (design §1) and the shared responsibility rule.

Rules (the added reactions are identical in all rules):
  R0  baseline: top-scoring candidate (today's rule, offgrowth.score.assign_calls_for_added_reactions)
  R1  responsibility: w(h) = odds(s_h)^beta1 * exp(beta3 * C_syn(h, r)); argmax w when w_max > tau_null, else no gene
  R2  R1 with tau_null = 1.0
  R3  R1 plus paralog competition: a gene wins at most one reaction per enzyme number
  R4  synteny only: R0 restricted to candidates with C_syn = 1 when any exists
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Dict, List, Mapping, Optional, Sequence, Set, Tuple

from .common import (
    BETA1,
    BETA3,
    NEIGHBOUR_MIN_SCORE,
    SYNTENY_WINDOW,
    TAU_NULL,
    GenomeLayout,
    base_seed_rxn_id,
    normalize_ec,
)
from offgrowth.mask import CandidateClaim, protein_aliases, strip_version
from offgrowth.score import ReactionGeneCall

RULES: Tuple[str, ...] = ("R0", "R1", "R2", "R3", "R4")


@dataclass(frozen=True)
class RuleConfig:
    rule: str
    tau_null: float = TAU_NULL
    beta1: float = BETA1
    beta3: float = BETA3
    window: int = SYNTENY_WINDOW
    neighbour_min_score: float = NEIGHBOUR_MIN_SCORE


def rule_config(rule: str) -> RuleConfig:
    if rule == "R2":
        return RuleConfig(rule="R2", tau_null=1.0)
    return RuleConfig(rule=rule)


def odds(score: float, beta1: float = BETA1) -> float:
    s = min(max(float(score), 1e-4), 1.0 - 1e-4)
    return (s / (1.0 - s)) ** beta1


class SyntenyContext:
    """Label-free genome context for C_syn: layout, candidate claims per protein, KEGG modules, transporter hits."""

    def __init__(
        self,
        layout: Optional[GenomeLayout],
        candidates: Sequence[CandidateClaim],
        ec_modules: Mapping[str, Set[str]],
        transporter_proteins_by_cpd: Mapping[str, Set[str]],
        neighbour_min_score: float = NEIGHBOUR_MIN_SCORE,
        window: int = SYNTENY_WINDOW,
    ):
        self.layout = layout
        self.ec_modules = ec_modules
        self.transporters = {cpd: {a for p in prots for a in protein_aliases(p)} for cpd, prots in transporter_proteins_by_cpd.items()}
        self.window = window
        self._cache: Dict[Tuple[str, Tuple[str, ...], Tuple[str, ...]], int] = {}
        self.claims_by_protein: Dict[str, Dict[str, float]] = {}
        for c in candidates:
            if c.score >= neighbour_min_score:
                ec = normalize_ec(c.ec)
                if ec in ec_modules:
                    d = self.claims_by_protein.setdefault(c.protein_id, {})
                    if c.score > d.get(ec, 0.0):
                        d[ec] = c.score

    def c_syn(self, protein_id: str, rxn_ecs: Set[str], substrate_cpds: Set[str]) -> int:
        if self.layout is None:
            return 0
        key = (protein_id, tuple(sorted(rxn_ecs)), tuple(sorted(substrate_cpds)))
        cached = self._cache.get(key)
        if cached is None:
            cached = self._c_syn(protein_id, rxn_ecs, substrate_cpds)
            self._cache[key] = cached
        return cached

    def _c_syn(self, protein_id: str, rxn_ecs: Set[str], substrate_cpds: Set[str]) -> int:
        rxn_modules: Set[str] = set()
        for ec in rxn_ecs:
            rxn_modules |= self.ec_modules.get(ec, set())
        for g in self.layout.neighbours(protein_id, self.window):
            if substrate_cpds:
                g_aliases = protein_aliases(g)
                if any(g_aliases & self.transporters.get(cpd, set()) for cpd in substrate_cpds):
                    return 1
            if rxn_modules:
                claims = self.claims_by_protein.get(g) or self.claims_by_protein.get(strip_version(g)) or {}
                for ec2 in claims:
                    if ec2 not in rxn_ecs and (self.ec_modules.get(ec2, set()) & rxn_modules):
                        return 1
        return 0


def candidates_by_reaction(
    candidates: Sequence[CandidateClaim],
    ec_to_reactions: Mapping[str, Sequence[str]],
) -> Dict[str, Dict[str, Tuple[float, str]]]:
    """{base reaction: {protein: (best score, EC through which it maps)}}."""
    out: Dict[str, Dict[str, Tuple[float, str]]] = {}
    for c in candidates:
        ec = normalize_ec(c.ec)
        for r in ec_to_reactions.get(ec, ()):
            base = base_seed_rxn_id(r)
            d = out.setdefault(base, {})
            cur = d.get(c.protein_id)
            if cur is None or c.score > cur[0] or (c.score == cur[0] and ec < cur[1]):
                d[c.protein_id] = (c.score, ec)
    return out


def classify_call(
    protein_id: str,
    base_rxn: str,
    allowed_rxns_by_protein: Mapping[str, Set[str]],
) -> str:
    rx = allowed_rxns_by_protein.get(protein_id)
    if rx is None:
        rx = allowed_rxns_by_protein.get(strip_version(protein_id))
    if rx is None:
        return "unconfirmed"
    return "supported" if base_rxn in rx else "harmful"


def allowed_reactions_by_protein(
    known_experimental_ecs: Mapping[str, Sequence[str]],
    ec_to_reactions: Mapping[str, Sequence[str]],
) -> Dict[str, Set[str]]:
    allowed: Dict[str, Set[str]] = {}
    for prot, ecs in known_experimental_ecs.items():
        rx: Set[str] = set()
        for e in ecs:
            for r in ec_to_reactions.get(normalize_ec(e), ()):
                rx.add(base_seed_rxn_id(r))
        for alias in protein_aliases(prot):
            allowed.setdefault(alias, set()).update(rx)
    return allowed


def responsibility_weights(
    cands: Mapping[str, Tuple[float, str]],
    rxn_ecs: Set[str],
    substrate_cpds: Set[str],
    ctx: Optional[SyntenyContext],
    config: RuleConfig,
) -> Dict[str, Tuple[float, int]]:
    """{protein: (w, C_syn)} for the candidates of one reaction."""
    weights: Dict[str, Tuple[float, int]] = {}
    for h, (s, _ec) in cands.items():
        c_syn = ctx.c_syn(h, rxn_ecs, substrate_cpds) if ctx is not None else 0
        weights[h] = (odds(s, config.beta1) * math.exp(config.beta3 * c_syn), c_syn)
    return weights


def assign_reactions(
    added: Mapping[str, Set[str]],
    rxn_candidates: Mapping[str, Mapping[str, Tuple[float, str]]],
    rxn_ecs: Mapping[str, Set[str]],
    substrate_by_test: Mapping[str, str],
    ctx: Optional[SyntenyContext],
    config: RuleConfig,
    allowed_rxns_by_protein: Mapping[str, Set[str]],
) -> Tuple[ReactionGeneCall, ...]:
    """Assign every added reaction under `config.rule`; returns calls in increasing reaction-id order."""
    calls: List[ReactionGeneCall] = []
    won: Set[Tuple[str, str]] = set()  # (protein, ec) pairs already won (R3)
    for base_rxn in sorted(added):
        cands = rxn_candidates.get(base_rxn, {})
        if not cands:
            calls.append(ReactionGeneCall(base_rxn, None, None, 0.0, "reaction_only"))
            continue
        ecs = rxn_ecs.get(base_rxn, set())
        substrates = {substrate_by_test[t] for t in added[base_rxn] if t in substrate_by_test}
        ordered: List[Tuple[float, str, float, str, int]]  # (-w, protein, score, ec, c_syn)
        if config.rule == "R0":
            ordered = [(-s, h, s, ec, 0) for h, (s, ec) in cands.items()]
            ordered.sort()
            chosen = ordered[0]
        elif config.rule == "R4":
            weights = responsibility_weights(cands, ecs, substrates, ctx, config)
            syn = [(-cands[h][0], h, cands[h][0], cands[h][1], c) for h, (w, c) in weights.items() if c == 1]
            if syn:
                syn.sort()
                chosen = syn[0]
            else:
                ordered = [(-s, h, s, ec, 0) for h, (s, ec) in cands.items()]
                ordered.sort()
                chosen = ordered[0]
        else:
            weights = responsibility_weights(cands, ecs, substrates, ctx, config)
            ordered = [(-w, h, cands[h][0], cands[h][1], c) for h, (w, c) in weights.items()]
            ordered.sort()
            chosen = None
            for neg_w, h, s, ec, c in ordered:
                if -neg_w <= config.tau_null:
                    break  # the null gene is at least as plausible as every remaining candidate
                if config.rule == "R3" and (h, ec) in won:
                    continue
                chosen = (neg_w, h, s, ec, c)
                break
            if chosen is None:
                calls.append(ReactionGeneCall(base_rxn, None, None, 0.0, "no_gene"))
                continue
        _, h, s, ec, c = chosen
        won.add((h, ec))
        calls.append(ReactionGeneCall(base_rxn, h, ec, float(s), classify_call(h, base_rxn, allowed_rxns_by_protein)))
    return tuple(calls)


def case_outcome(
    calls: Sequence[ReactionGeneCall],
    truth_protein: str,
    target_reactions: Set[str],
) -> Dict[str, object]:
    aliases = protein_aliases(truth_protein)
    recovered = any(
        c.protein_id is not None
        and (c.protein_id in aliases or strip_version(c.protein_id) in aliases)
        and c.reaction_id in target_reactions
        for c in calls
    )
    reaction_recovered = any(c.reaction_id in target_reactions for c in calls)
    return {
        "recovered": recovered,
        "reaction_recovered": reaction_recovered,
        "n_gene_calls": sum(1 for c in calls if c.protein_id is not None),
        "n_no_gene": sum(1 for c in calls if c.call_kind == "no_gene"),
        "n_reaction_only": sum(1 for c in calls if c.call_kind == "reaction_only"),
        "n_supported": sum(1 for c in calls if c.call_kind == "supported"),
        "n_harmful": sum(1 for c in calls if c.call_kind == "harmful"),
        "n_unconfirmed": sum(1 for c in calls if c.call_kind == "unconfirmed"),
    }
