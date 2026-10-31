#!/usr/bin/env python3
"""Experiment 6: knockout growth on the fixed panel against Fitness Browser specific phenotypes (design §6).

Per genome: model variants M0 (gapsmith gap-filled model), M1 (plus the unmasked union A1 ∪ A2 ∪ A3 added reactions
assigned by the baseline rule R0), M2 (same reactions assigned by rule R1), M3_seed1..5 (same reactions assigned to
random model genes). For every panel well matched to a Fitness Browser experiment: sole-source medium, wild-type
growth, single-gene knockouts of every model gene; prediction = 1 - growth(knockout)/growth(wild type). Truth: the
gene has a specific phenotype with negative fitness in that experiment. Metrics: AUC over (gene, experiment) pairs
(model genes; added-reaction genes), precision at the top 50, label-permutation control (seeds 1–5).
"""

from __future__ import annotations

import argparse
import copy
import json
import math
import multiprocessing as mp
import os
from pathlib import Path
import random
import time
from typing import Dict, List, Mapping, Optional, Sequence, Set, Tuple

import cobra

from .assignment import SyntenyContext, assign_reactions, candidates_by_reaction, rule_config
from .common import SHUFFLE_SEEDS, added_reactions_with_tests, environment_receipt, read_genome_results, stored_tests_by_arm, write_json
from .inputs import GenomeInputs, ReferenceTables
from offgrowth.arms import _find_biomass_rxn_id, apply_medium, build_sole_source_medium, load_gapsmith_medium_csv
from offgrowth.mask import base_seed_rxn_id, protein_aliases, strip_version

FLUX_TOL = 1e-6
TOP_K_PRECISION = 50
VARIANTS = ("M0", "M1", "M2") + tuple(f"M3_seed{s}" for s in SHUFFLE_SEEDS)


def auc(scores: Sequence[float], labels: Sequence[int], weights: Optional[Sequence[float]] = None) -> Optional[float]:
    """Rank-based AUC with ties (Mann–Whitney), optionally with per-item integer weights (bootstrap multiplicities).

    Returns None when a class is empty. Implemented with numpy so that bootstrap resampling stays fast.
    """
    import numpy as np

    sc = np.asarray(scores, dtype=float)
    lb = np.asarray(labels, dtype=float)
    w = np.ones(len(sc)) if weights is None else np.asarray(weights, dtype=float)
    if len(sc) == 0:
        return None
    keep = w > 0
    sc, lb, w = sc[keep], lb[keep], w[keep]
    w_pos = float((w * lb).sum())
    w_neg = float((w * (1.0 - lb)).sum())
    if w_pos <= 0 or w_neg <= 0:
        return None
    order = np.argsort(sc, kind="stable")
    s_sorted, w_sorted, l_sorted = sc[order], w[order], lb[order]
    # tie blocks: average rank of a block of total weight W starting after cumulative weight C is C + (W + 1) / 2
    _, starts = np.unique(s_sorted, return_index=True)
    block_id = np.cumsum(np.isin(np.arange(len(s_sorted)), starts)) - 1
    block_w = np.bincount(block_id, weights=w_sorted)
    cum_before = np.concatenate(([0.0], np.cumsum(block_w)[:-1]))
    rank_of_block = cum_before + (block_w + 1.0) / 2.0
    ranks = rank_of_block[block_id]
    r_pos = float((ranks * w_sorted * l_sorted).sum())
    return (r_pos - w_pos * (w_pos + 1.0) / 2.0) / (w_pos * w_neg)


def precision_at_k(scores: Sequence[float], labels: Sequence[int], k: int = TOP_K_PRECISION) -> Optional[float]:
    if not scores:
        return None
    top = sorted(range(len(scores)), key=lambda i: (-scores[i], i))[:k]
    return sum(labels[i] for i in top) / float(len(top))


def model_gene_index(model: cobra.Model) -> Dict[str, str]:
    idx: Dict[str, str] = {}
    for g in model.genes:
        for alias in protein_aliases(g.id):
            idx.setdefault(alias, g.id)
    return idx


def add_assigned_reactions(
    base: cobra.Model,
    tables: ReferenceTables,
    assignments: Mapping[str, Optional[str]],
) -> Tuple[cobra.Model, Dict[str, str]]:
    """Copy of `base` with the added reactions present and their gene–reaction rule set to the assigned gene.

    A reaction already in the model keeps its stoichiometry; its rule is set only when it has none. Returns the model
    and {base reaction: model gene id} of the rules actually set.
    """
    model = copy.deepcopy(base)
    gidx = model_gene_index(model)
    existing = {base_seed_rxn_id(r.id): r for r in model.reactions if not r.id.startswith("EX_")}
    mets = {m.id: m for m in model.metabolites}
    new_mets: List[cobra.Metabolite] = []
    new_rxns: List[cobra.Reaction] = []
    set_rules: Dict[str, str] = {}
    for base_id in sorted(assignments):
        gene = assignments[base_id]
        gene_id = (gidx.get(gene) or gidx.get(strip_version(gene)) or gene) if gene else None
        rxn = existing.get(base_id)
        if rxn is None:
            urxn = tables.reactions.get(base_id)
            if urxn is None:
                for key, val in tables.reactions.items():
                    if base_seed_rxn_id(key) == base_id:
                        urxn = val
                        break
            if urxn is None or not urxn.stoichiometry:
                continue
            for t in urxn.stoichiometry:
                if t.species_id not in mets:
                    m = cobra.Metabolite(id=t.species_id, name=t.species_id, compartment=t.compartment)
                    mets[t.species_id] = m
                    new_mets.append(m)
            lb, ub = urxn.bounds
            rxn = cobra.Reaction(id=f"{base_id}_c0", name=urxn.name, lower_bound=lb, upper_bound=ub)
            rxn.add_metabolites({mets[t.species_id]: t.coeff for t in urxn.stoichiometry})
            new_rxns.append(rxn)
            existing[base_id] = rxn
        if gene_id and not (rxn.gene_reaction_rule or "").strip():
            rxn.gene_reaction_rule = gene_id
            set_rules[base_id] = gene_id
    if new_mets:
        model.add_metabolites(new_mets)
    if new_rxns:
        model.add_reactions(new_rxns)
    for base_id, gene_id in set_rules.items():
        existing[base_id].gene_reaction_rule = gene_id  # rules on new reactions need the reaction in the model
    ex_ids = [f"EX_{m.id}" for m in model.metabolites if m.id.endswith("_e0") and f"EX_{m.id}" not in model.reactions]
    if ex_ids:
        from offgrowth.arms import _ensure_exchange_reactions

        _ensure_exchange_reactions(model, ex_ids)
    return model, set_rules


_WORKER_MODELS: Dict[str, cobra.Model] = {}
_WORKER_MEDIA: Dict[str, Dict[str, float]] = {}
_WORKER_BIOMASS: Dict[str, str] = {}


def _init_worker(models: Dict[str, cobra.Model], media: Dict[str, Dict[str, float]], biomass: Dict[str, str]) -> None:
    global _WORKER_MODELS, _WORKER_MEDIA, _WORKER_BIOMASS
    _WORKER_MODELS, _WORKER_MEDIA, _WORKER_BIOMASS = models, media, biomass
    os.environ["OMP_NUM_THREADS"] = "1"


def _knockout_task(task: Tuple[str, str]) -> Tuple[str, str, float, Dict[str, float], int, float]:
    """(variant, panel id) -> wild-type growth and knockout growth per gene."""
    variant, panel_id = task
    t0 = time.process_time()
    model = _WORKER_MODELS[variant]
    apply_medium(model, _WORKER_MEDIA[panel_id])
    model.objective = _WORKER_BIOMASS[variant]
    model.objective_direction = "max"
    wt = model.slim_optimize(error_value=0.0)
    wt = 0.0 if wt is None or math.isnan(wt) else max(0.0, float(wt))
    ko: Dict[str, float] = {}
    n_lps = 1
    if wt > FLUX_TOL:
        for gene in model.genes:
            if not gene.reactions:
                continue
            with model:
                gene.knock_out()
                v = model.slim_optimize(error_value=0.0)
            n_lps += 1
            v = 0.0 if v is None or math.isnan(v) else max(0.0, float(v))
            ko[gene.id] = v / wt
    return variant, panel_id, wt, ko, n_lps, time.process_time() - t0


def run_genome(a: argparse.Namespace) -> Dict[str, object]:
    t_start = time.time()
    fitness = json.loads(Path(a.fitness_json).read_text())
    tables = ReferenceTables.load(a.seed_reactions, a.ec_aliases, a.seed_metabolites, a.kegg_modules, a.panel_tsv)
    gi = GenomeInputs.load(a.genome_id, a.layout_tsv, a.candidates_tsv, a.genome_model_dir / f"{a.genome_id}-Transporter.tbl")
    filled = cobra.io.read_sbml_model(str(a.genome_model_dir / f"{a.genome_id}-filled.xml"))
    filled.solver = "glpk"
    default_medium = load_gapsmith_medium_csv(a.genome_model_dir / f"{a.genome_id}-medium.csv")

    # Added reactions of the unmasked union run and their assignments
    payload = read_genome_results(a.genomes_dir / f"{a.genome_id}.outputs.tar")
    unmasked = payload["unmasked_arm_results"]
    tests = stored_tests_by_arm(unmasked)
    added = added_reactions_with_tests(t for arm in ("A1", "A2", "A3") for t in tests.get(arm, ()))
    union_added = {base_seed_rxn_id(str(r)) for r in unmasked.get("A1_A2_A3", {}).get("added_reactions", ())}
    added = {r: added.get(r, set()) for r in union_added}
    rxn_cands = candidates_by_reaction(gi.candidates, tables.ec_to_reactions)
    ctx = SyntenyContext(gi.layout, gi.candidates, tables.ec_modules, gi.transporter_proteins_by_cpd)
    assignments: Dict[str, Dict[str, Optional[str]]] = {}
    for rule, label in (("R0", "M1"), ("R1", "M2")):
        calls = assign_reactions(added, rxn_cands, tables.rxn_ecs, tables.substrate_by_test(), ctx, rule_config(rule), {})
        assignments[label] = {c.reaction_id: c.protein_id for c in calls}
    model_genes = sorted(g.id for g in filled.genes if g.reactions)
    for seed in SHUFFLE_SEEDS:
        rng = random.Random(seed)
        assignments[f"M3_seed{seed}"] = {r: rng.choice(model_genes) for r in sorted(added)}
    models: Dict[str, cobra.Model] = {"M0": filled}
    rules_set: Dict[str, Dict[str, str]] = {"M0": {}}
    for label in VARIANTS[1:]:
        models[label], rules_set[label] = add_assigned_reactions(filled, tables, assignments[label])
        models[label].solver = "glpk"
    biomass = {label: _find_biomass_rxn_id(m) for label, m in models.items()}

    # Conditions: panel wells matched to at least one experiment
    wells = sorted({pid for e in fitness["experiments"].values() for pid in e["panel_ids"]})
    media = {pid: build_sole_source_medium(default_medium, tables.panel_by_id[pid], tables.metabolites) for pid in wells if pid in tables.panel_by_id}
    for m in models.values():
        apply_medium(m, dict(default_medium))
    tasks = [(v, pid) for v in VARIANTS for pid in sorted(media)]
    workers = max(1, int(a.workers))
    ctx_mp = mp.get_context("fork")
    results: List[Tuple[str, str, float, Dict[str, float], int, float]] = []
    with ctx_mp.Pool(processes=min(workers, max(1, len(tasks))), initializer=_init_worker, initargs=(models, media, biomass)) as pool:
        for k, res in enumerate(pool.imap_unordered(_knockout_task, tasks), 1):
            results.append(res)
            if k % 25 == 0 or k == len(tasks):
                print(json.dumps({"progress": k, "of": len(tasks)}), flush=True)
    ko_by = {(v, pid): (wt, ko) for v, pid, wt, ko, _, _ in results}
    n_lps = sum(r[4] for r in results)
    cpu = sum(r[5] for r in results)

    # Scoring against fitness
    gene_map: Dict[str, str] = fitness["gene_map"]  # locus -> protein id
    gidx_by_variant = {v: model_gene_index(m) for v, m in models.items()}
    added_genes = {g for lab in ("M1", "M2") for g in rules_set[lab].values()}
    pairs: Dict[str, Dict[str, Dict[str, List]]] = {v: {"model_genes": {"scores": [], "labels": [], "keys": []}, "added_genes": {"scores": [], "labels": [], "keys": []}} for v in VARIANTS}
    unscorable: Dict[str, int] = {v: 0 for v in VARIANTS}
    n_exp_scored: Dict[str, int] = {v: 0 for v in VARIANTS}
    for exp_name, exp in fitness["experiments"].items():
        fit = fitness["fitness"].get(exp_name, {})
        specific = set(fitness["specific_phenotypes"].get(exp_name, ()))
        for v in VARIANTS:
            gidx = gidx_by_variant[v]
            scored_any = False
            for pid in exp["panel_ids"]:
                entry = ko_by.get((v, pid))
                if entry is None or entry[0] <= FLUX_TOL:
                    continue
                wt, ko = entry
                for locus, (f, t) in fit.items():
                    prot = gene_map.get(locus)
                    if prot is None:
                        continue
                    mg = gidx.get(prot) or gidx.get(strip_version(prot))
                    if mg is None or mg not in ko:
                        continue
                    label = 1 if (locus in specific and f < 0) else 0
                    score = 1.0 - ko[mg]
                    key = f"{locus}|{exp_name}|{pid}"
                    pairs[v]["model_genes"]["scores"].append(score)
                    pairs[v]["model_genes"]["labels"].append(label)
                    pairs[v]["model_genes"]["keys"].append(key)
                    if mg in added_genes:
                        pairs[v]["added_genes"]["scores"].append(score)
                        pairs[v]["added_genes"]["labels"].append(label)
                        pairs[v]["added_genes"]["keys"].append(key)
                    scored_any = True
            if scored_any:
                n_exp_scored[v] += 1
            else:
                unscorable[v] += 1
    metrics: Dict[str, object] = {}
    for v in VARIANTS:
        metrics[v] = {}
        for subset, d in pairs[v].items():
            sc, lb = d["scores"], d["labels"]
            perm = []
            for seed in SHUFFLE_SEEDS:
                rng = random.Random(seed)
                pl = list(lb)
                rng.shuffle(pl)
                perm.append(auc(sc, pl))
            metrics[v][subset] = {"n_pairs": len(sc), "n_positive": sum(lb), "auc": auc(sc, lb), "precision_at_50": precision_at_k(sc, lb), "permuted_label_auc": perm}
    # Gene-level coverage
    model_gene_ids = set(model_genes)
    matched_model_genes = sum(1 for prot in gene_map.values() if (gidx_by_variant["M0"].get(prot) or gidx_by_variant["M0"].get(strip_version(prot))) in model_gene_ids)
    out = {
        "genome_id": a.genome_id, "org_id": fitness["org_id"], "variants": list(VARIANTS),
        "n_model_genes": len(model_genes), "n_model_genes_with_fitness_gene": matched_model_genes,
        "model_gene_match_fraction": matched_model_genes / float(len(model_genes)) if model_genes else 0.0,
        "n_added_reactions": len(added), "rules_set": {v: len(r) for v, r in rules_set.items()}, "added_genes": sorted(added_genes),
        "n_wells": len(media), "wells_wt_growing": {v: sum(1 for (vv, pid), (wt, _) in ko_by.items() if vv == v and wt > FLUX_TOL) for v in VARIANTS},
        "experiments_scored": n_exp_scored, "experiments_unscorable": unscorable,
        "n_lps": n_lps, "cpu_seconds": cpu, "wall_seconds": time.time() - t_start,
        "metrics": metrics,
        "pairs": {v: {s: {"scores": [round(x, 6) for x in d["scores"]], "labels": d["labels"]} for s, d in pairs[v].items()} for v in VARIANTS},
        "assignments_M1_M2_differ": sum(1 for r in added if assignments["M1"].get(r) != assignments["M2"].get(r)),
        "receipt": environment_receipt({"fitness_json": a.fitness_json, "panel_tsv": a.panel_tsv, "filled_xml": a.genome_model_dir / f"{a.genome_id}-filled.xml", "candidates_tsv": a.candidates_tsv}),
    }
    a.out_dir.mkdir(parents=True, exist_ok=True)
    write_json(a.out_dir / f"{a.genome_id}.exp6.json", out)
    print(json.dumps({"genome": a.genome_id, "lps": n_lps, "cpu_s": round(cpu), "auc_M0": metrics["M0"]["model_genes"]["auc"], "auc_M1_added": metrics["M1"]["added_genes"]["auc"], "auc_M2_added": metrics["M2"]["added_genes"]["auc"]}), flush=True)
    return out


def main(argv: Optional[Sequence[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--genome-id", required=True)
    p.add_argument("--fitness-json", type=Path, required=True)
    p.add_argument("--genome-model-dir", type=Path, required=True)
    p.add_argument("--genomes-dir", type=Path, required=True, help="off-growth outputs (unmasked union added reactions)")
    p.add_argument("--candidates-tsv", type=Path, required=True)
    p.add_argument("--layout-tsv", type=Path, required=True)
    p.add_argument("--seed-reactions", type=Path, required=True)
    p.add_argument("--ec-aliases", type=Path, default=None)
    p.add_argument("--seed-metabolites", type=Path, required=True)
    p.add_argument("--kegg-modules", type=Path, required=True)
    p.add_argument("--panel-tsv", type=Path, required=True)
    p.add_argument("--workers", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", "4")))
    p.add_argument("--out-dir", type=Path, required=True)
    a = p.parse_args(argv)
    run_genome(a)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
