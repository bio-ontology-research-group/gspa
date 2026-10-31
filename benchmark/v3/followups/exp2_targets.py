#!/usr/bin/env python3
"""Experiment 2: unblocking (A2) and dead-end (A3) tests as a gated evidence source for the Stage 2B update (design §2).

For each genome the stored unmasked A2/A3 test results (label-free) give the added reactions. For each predictor cell
the stage-0 output supplies the candidate scores. For every added reaction r and candidate h (score > 1e-5 for an
enzyme number mapped to r): responsibility r(h) by rule R1 (odds * exp(1.5 C_syn), null tau 0.05); gate G1: r(h)
above the null share; gate G2: r shares a non-currency metabolite with a reaction the draft model assigns to a gene
within +/- 10 loci of h. Update: delta logit = 2.2 * 1 * r(h) * (1 - p_j), p_j = max candidate score for r's enzyme
numbers; written as a sparse update record per cell and arm for `price_score`.

Arms: A2, A3, A2A3, A2A3_nogate (no G2), A2A3_shuffled_seed1..5 (candidate gene identities permuted).
"""

from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
import csv
import json
import math
import os
from pathlib import Path
import random
import tempfile
from typing import Dict, List, Mapping, Optional, Sequence, Set, Tuple

from .assignment import SyntenyContext, odds
from .common import (
    BETA3,
    CANDIDATE_FLOOR,
    GAMMA_FBA,
    SHUFFLE_SEEDS,
    SYNTENY_WINDOW,
    TAU_NULL,
    GenomeLayout,
    added_reactions_with_tests,
    design_receipt,
    environment_receipt,
    read_genome_results,
    stored_tests_by_arm,
    write_json,
)
from .inputs import ReferenceTables, extract_model_member
from .price_score import locate_cell
from offgrowth.mask import CandidateClaim, base_seed_rxn_id, protein_aliases, strip_version
from offgrowth.panel import normalize_ec

N_CURRENCY = 20
EVIDENCE_ARMS = ("A2", "A3", "A2A3", "A2A3_nogate")


def logit(p: float, eps: float = 1e-6) -> float:
    pc = min(max(float(p), eps), 1.0 - eps)
    return math.log(pc / (1.0 - pc))


def sigmoid(x: float) -> float:
    if x >= 0:
        z = math.exp(-x)
        return 1.0 / (1.0 + z)
    z = math.exp(x)
    return z / (1.0 + z)


def currency_metabolites(tables: ReferenceTables, n: int = N_CURRENCY) -> List[str]:
    """The n most connected compounds of the usable universal SEED table (compartment ignored)."""
    counts: Counter = Counter()
    for rxn in tables.reactions.values():
        if not rxn.is_usable:
            continue
        for cpd in {t.cpd_id for t in rxn.stoichiometry}:
            counts[cpd] += 1
    return [c for c, _ in counts.most_common(n)]


def reaction_compounds_index(tables: ReferenceTables, currency: Set[str]) -> Dict[str, Set[str]]:
    """{base reaction id: non-currency compounds} for every reaction of the universal SEED table."""
    out: Dict[str, Set[str]] = {}
    for key, rxn in tables.reactions.items():
        out[base_seed_rxn_id(key)] = {t.cpd_id for t in rxn.stoichiometry if t.cpd_id not in currency}
    return out


def draft_gene_reactions(draft_xml: Path) -> Dict[str, Set[str]]:
    """{protein alias: {base reaction ids the draft model assigns to it}} from the gapsmith draft SBML."""
    import cobra

    model = cobra.io.read_sbml_model(str(draft_xml))
    out: Dict[str, Set[str]] = {}
    for rxn in model.reactions:
        if rxn.id.startswith("EX_") or not (rxn.gene_reaction_rule or "").strip():
            continue
        base = base_seed_rxn_id(rxn.id)
        for g in rxn.genes:
            for alias in protein_aliases(g.id):
                out.setdefault(alias, set()).add(base)
    return out


def stream_stage0_candidates(path: Path, relevant_ecs: Set[str], floor: float = CANDIDATE_FLOOR) -> Tuple[List[CandidateClaim], Dict[str, Dict[str, float]]]:
    """Claims of a stage-0 output whose enzyme number is relevant and whose score exceeds the floor (streamed)."""
    claims: List[CandidateClaim] = []
    by_protein: Dict[str, Dict[str, float]] = {}
    with Path(path).open("r", encoding="utf-8", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        i_pid, i_type, i_fid, i_prob = (header.index(c) for c in ("protein_id", "type", "function_id", "posterior_prob"))
        for row in reader:
            if row[i_type] != "EC":
                continue
            ec = normalize_ec(row[i_fid])
            if ec not in relevant_ecs:
                continue
            s = float(row[i_prob])
            if s > floor:
                claims.append(CandidateClaim(row[i_pid], ec, s))
                by_protein.setdefault(row[i_pid], {})[ec] = s
    return claims, by_protein


def compute_updates(
    added: Mapping[str, Set[str]],
    claims: Sequence[CandidateClaim],
    tables: ReferenceTables,
    layout: GenomeLayout,
    draft_genes: Mapping[str, Set[str]],
    rxn_cpds: Mapping[str, Set[str]],
    gate_g2: bool,
) -> Tuple[Dict[str, Dict[str, float]], Dict[str, int]]:
    """Sparse log-odds updates {protein: {ec: new_prob}} and counters."""
    ctx = SyntenyContext(layout, claims, tables.ec_modules, {})
    by_rxn: Dict[str, Dict[str, Tuple[float, str]]] = {}
    scores_by_protein: Dict[str, Dict[str, float]] = {}
    for c in claims:
        scores_by_protein.setdefault(c.protein_id, {})[c.ec] = c.score
        for r in tables.ec_to_reactions.get(c.ec, ()):
            base = base_seed_rxn_id(r)
            if base in added:
                cur = by_rxn.setdefault(base, {}).get(c.protein_id)
                if cur is None or c.score > cur[0]:
                    by_rxn[base][c.protein_id] = (c.score, c.ec)
    delta: Dict[Tuple[str, str], float] = {}
    counters = Counter()
    for base_rxn in sorted(added):
        cands = by_rxn.get(base_rxn)
        if not cands:
            counters["reactions_without_candidates"] += 1
            continue
        counters["reactions_with_candidates"] += 1
        rxn_ecs = tables.rxn_ecs.get(base_rxn, set())
        p_j = max(s for s, _ in cands.values())
        delta_repair = max(0.0, 1.0 - p_j)
        weights = {h: (odds(s) * math.exp(BETA3 * ctx.c_syn(h, rxn_ecs, set())), ec) for h, (s, ec) in cands.items()}
        denom = TAU_NULL + sum(w for w, _ in weights.values())
        null_share = TAU_NULL / denom
        r_cpds = rxn_cpds.get(base_rxn, set()) if gate_g2 else set()
        for h, (w, ec) in weights.items():
            r_h = w / denom
            if r_h <= null_share:
                counters["blocked_by_null"] += 1
                continue
            if gate_g2:
                ok = False
                for g in layout.neighbours(h, SYNTENY_WINDOW):
                    g_rxns = draft_genes.get(g) or draft_genes.get(strip_version(g)) or set()
                    for gr in g_rxns:
                        if rxn_cpds.get(gr, set()) & r_cpds:
                            ok = True
                            break
                    if ok:
                        break
                if not ok:
                    counters["blocked_by_coherence_gate"] += 1
                    continue
            counters["gated_updates"] += 1
            delta[(h, ec)] = delta.get((h, ec), 0.0) + GAMMA_FBA * 1.0 * r_h * delta_repair
    updates: Dict[str, Dict[str, float]] = {}
    for (h, ec), d in delta.items():
        if d <= 0.0:
            continue
        s_old = scores_by_protein[h][ec]
        s_new = sigmoid(logit(s_old) + d)
        if abs(s_new - s_old) > 1e-5:
            updates.setdefault(h, {})[f"EC:{ec}"] = s_new
    counters["changed_scores"] = sum(len(m) for m in updates.values())
    return updates, dict(counters)


def label_outcomes(updates: Mapping[str, Mapping[str, float]], labels: Mapping[str, Sequence[str]]) -> Dict[str, int]:
    """Among updated (protein, ec) pairs on labelled proteins: supported (ec is a label) vs harmful."""
    lab: Dict[str, Set[str]] = {}
    for prot, ecs in labels.items():
        for alias in protein_aliases(prot):
            lab.setdefault(alias, set()).update(normalize_ec(e) for e in ecs)
    supported = harmful = 0
    for h, m in updates.items():
        known = lab.get(h) or lab.get(strip_version(h))
        if not known:
            continue
        for ec in m:
            if normalize_ec(ec) in known:
                supported += 1
            else:
                harmful += 1
    return {"labelled_supported": supported, "labelled_harmful": harmful}


def run_genome(genome_id: str, a: Mapping[str, str]) -> Dict[str, object]:
    tables = ReferenceTables.load(Path(a["seed_reactions"]), Path(a["ec_aliases"]) if a.get("ec_aliases") else None, None, Path(a["kegg_modules"]), Path(a["panel_tsv"]))
    currency = set(currency_metabolites(tables))
    rxn_cpds = reaction_compounds_index(tables, currency)
    payload = read_genome_results(Path(a["genomes_dir"]) / f"{genome_id}.outputs.tar")
    unmasked = payload["unmasked_arm_results"]
    tests = stored_tests_by_arm(unmasked)
    added_by_arm = {"A2": added_reactions_with_tests(tests.get("A2", ())), "A3": added_reactions_with_tests(tests.get("A3", ()))}
    both: Dict[str, Set[str]] = {}
    for arm in ("A2", "A3"):
        for r, t in added_by_arm[arm].items():
            both.setdefault(r, set()).update(t)
    added_by_arm["A2A3"] = both
    added_by_arm["A2A3_nogate"] = both
    relevant_ecs: Set[str] = set()
    for r in both:
        relevant_ecs |= tables.rxn_ecs.get(r, set())
    layout = GenomeLayout.from_tsv(Path(a["layout_dir"]) / f"{genome_id}.layout.tsv")
    with tempfile.TemporaryDirectory(prefix=f"exp2-{genome_id}-", dir=a.get("tmp_dir") or None) as td:
        draft = extract_model_member(Path(a["models_tar"]), genome_id, "-draft.xml", Path(td))
        draft_genes = draft_gene_reactions(draft)
    labels: Dict[str, List[str]] = {k: list(v) for k, v in (payload.get("known_experimental_ecs") or {}).items()}
    sp = Path(a["swissprot_tsv"]) if a.get("swissprot_tsv") else None
    if sp is not None and sp.exists():
        with sp.open("r", encoding="utf-8", newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                if row.get("genome_id") == genome_id and (row.get("tier") or "experimental") == "experimental":
                    labels.setdefault(row["protein_id"], []).append(row["function_id"])
    plan = json.loads((Path(a["collection"]) / "plan.json").read_text())
    cells = [c for c in plan["cells"] if c["assembly"] == genome_id]
    out_dir = Path(a["out_dir"])
    summary: Dict[str, object] = {"genome_id": genome_id, "n_added_A2": len(added_by_arm["A2"]), "n_added_A3": len(added_by_arm["A3"]), "n_added_union": len(both),
                                  "currency": sorted(currency), "cells": {}}
    for cell in cells:
        located = locate_cell(Path(a["collection"]), cell)
        claims, _ = stream_stage0_candidates(Path(located["outputs"]["stage0"]["path"]), relevant_ecs)
        cell_summary: Dict[str, object] = {"n_relevant_claims": len(claims), "arms": {}}
        for arm in EVIDENCE_ARMS:
            updates, counters = compute_updates(added_by_arm[arm], claims, tables, layout, draft_genes, rxn_cpds, gate_g2=(arm != "A2A3_nogate"))
            lo = label_outcomes(updates, labels)
            write_json(out_dir / f"{cell['index']}.{arm}.json", {"cell": dict(cell), "arm": arm, "updates": updates, "diagnostics": {**counters, **lo}})
            cell_summary["arms"][arm] = {**counters, **lo}
        model_genes = [str(g) for g in payload.get("model_genes", ())]
        for seed in SHUFFLE_SEEDS:
            ids = sorted({c.protein_id for c in claims} | set(model_genes))
            rng = random.Random(seed)
            shuffled = list(ids)
            rng.shuffle(shuffled)
            mapping = dict(zip(ids, shuffled))
            perm = [CandidateClaim(mapping.get(c.protein_id, c.protein_id), c.ec, c.score) for c in claims]
            updates, counters = compute_updates(both, perm, tables, layout, draft_genes, rxn_cpds, gate_g2=True)
            lo = label_outcomes(updates, labels)
            arm = f"A2A3_shuffled_seed{seed}"
            write_json(out_dir / f"{cell['index']}.{arm}.json", {"cell": dict(cell), "arm": arm, "updates": updates, "diagnostics": {**counters, **lo}})
            cell_summary["arms"][arm] = {**counters, **lo}
        summary["cells"][str(cell["index"])] = {"predictor": cell["predictor"], **cell_summary}
        print(json.dumps({"genome": genome_id, "cell": cell["index"], "predictor": cell["predictor"], "A2A3_changed": cell_summary["arms"]["A2A3"]["changed_scores"]}), flush=True)
    write_json(out_dir / f"{genome_id}.summary.json", summary)
    return summary


def _worker(task: Tuple[str, Dict[str, str]]) -> Dict[str, object]:
    genome_id, a = task
    try:
        return run_genome(genome_id, a)
    except Exception as exc:
        err = {"genome_id": genome_id, "status": "failed", "error": repr(exc)}
        write_json(Path(a["out_dir"]) / f"{genome_id}.failed.json", err)
        print(json.dumps(err), flush=True)
        return err


def main(argv: Optional[Sequence[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--collection", type=Path, required=True, help="price-stage2b-rc1/collection")
    p.add_argument("--genomes-dir", type=Path, required=True)
    p.add_argument("--layout-dir", type=Path, required=True)
    p.add_argument("--models-tar", type=Path, required=True)
    p.add_argument("--seed-reactions", type=Path, required=True)
    p.add_argument("--ec-aliases", type=Path, default=None)
    p.add_argument("--kegg-modules", type=Path, required=True)
    p.add_argument("--panel-tsv", type=Path, required=True)
    p.add_argument("--swissprot-tsv", type=Path, default=None)
    p.add_argument("--genome-ids", nargs="*", default=None)
    p.add_argument("--workers", type=int, default=max(1, min(22, int(os.environ.get("SLURM_CPUS_PER_TASK", "4")))))
    p.add_argument("--tmp-dir", type=Path, default=None)
    p.add_argument("--out-dir", type=Path, required=True)
    a = p.parse_args(argv)
    plan = json.loads((a.collection / "plan.json").read_text())
    genome_ids = a.genome_ids or sorted({c["assembly"] for c in plan["cells"]})
    a.out_dir.mkdir(parents=True, exist_ok=True)
    args = {"collection": str(a.collection), "genomes_dir": str(a.genomes_dir), "layout_dir": str(a.layout_dir), "models_tar": str(a.models_tar),
            "seed_reactions": str(a.seed_reactions), "ec_aliases": str(a.ec_aliases) if a.ec_aliases else "", "kegg_modules": str(a.kegg_modules),
            "panel_tsv": str(a.panel_tsv), "swissprot_tsv": str(a.swissprot_tsv) if a.swissprot_tsv else "", "tmp_dir": str(a.tmp_dir) if a.tmp_dir else "", "out_dir": str(a.out_dir)}
    write_json(a.out_dir / "receipt.json", {**environment_receipt({"kegg_modules": a.kegg_modules, "seed_reactions": a.seed_reactions, "panel_tsv": a.panel_tsv}), "genomes": genome_ids,
                                            "constants": {"beta3": BETA3, "tau_null": TAU_NULL, "gamma_fba": GAMMA_FBA, "window": SYNTENY_WINDOW, "candidate_floor": CANDIDATE_FLOOR, "n_currency": N_CURRENCY}})
    tasks = [(g, args) for g in genome_ids]
    if a.workers > 1 and len(tasks) > 1:
        with ProcessPoolExecutor(max_workers=min(a.workers, len(tasks))) as pool:
            results = list(pool.map(_worker, tasks))
    else:
        results = [_worker(t) for t in tasks]
    failed = [r for r in results if r.get("status") == "failed"]
    write_json(a.out_dir / "exp2_genome_summaries.json", {"genomes": results, "failed": len(failed)})
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
