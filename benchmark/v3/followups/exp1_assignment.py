#!/usr/bin/env python3
"""Experiment 1: competitive gene assignment, re-scored on the stored off-growth outputs (design §1).

No new linear programs. For every case the stored A1 ∪ A2 ∪ A3 test results (added reactions and the test that
added them) are re-assigned under rules R0–R4 with the unmasked candidate table, the genome layout, the KEGG module
table and the gapsmith transporter table of the genome. Controls: shuffled candidate genes (seeds 1–5) for R1 and
R3; shuffled layout (seeds 1–5) for R1.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict
import json
import os
from pathlib import Path
import tempfile
from typing import Dict, List, Mapping, Optional, Sequence, Set, Tuple

from .assignment import (
    RULES,
    SyntenyContext,
    allowed_reactions_by_protein,
    assign_reactions,
    candidates_by_reaction,
    case_outcome,
    rule_config,
)
from .common import (
    SHUFFLE_SEEDS,
    added_reactions_with_tests,
    environment_receipt,
    load_known_experimental_ecs,
    read_genome_results,
    stored_tests_by_arm,
    write_json,
    write_tsv,
)
from .inputs import GenomeInputs, ReferenceTables, extract_model_member
from offgrowth.mask import TruthCase, base_seed_rxn_id, resolve_equivalent_ecs_and_reactions
from offgrowth.score import (
    BOOTSTRAP_RESAMPLES,
    BOOTSTRAP_SEED,
    CaseArmEvaluation,
    exact_one_sided_mcnemar_p,
    permute_candidates_within_genome,
    summarize_arms_with_bootstrap,
)

ARM_OF_RECORD = "A1_A2_A3"
CRITERION_MIN_CORRECT = 35
CRITERION_MAX_HARM = 0.85
CRITERION_MAX_SHUFFLED = 2


def _evaluation(case: TruthCase, stratum: str, arm: str, calls, outcome: Mapping[str, object], target: Set[str]) -> CaseArmEvaluation:
    return CaseArmEvaluation(
        case_id=case.case_id,
        genome_id=case.genome_id,
        protein_id=case.protein_id,
        truth_ecs=tuple(case.truth_ecs),
        target_reactions=tuple(sorted(target)),
        stratum=stratum,
        arm=arm,
        reaction_recovered=bool(outcome["reaction_recovered"]),
        recovered=bool(outcome["recovered"]),
        calls=tuple(calls),
        n_added_reactions=len(calls),
        n_gene_calls=int(outcome["n_gene_calls"]),
        n_reaction_only_calls=int(outcome["n_reaction_only"]),
        n_supported_calls=int(outcome["n_supported"]),
        n_harmful_calls=int(outcome["n_harmful"]),
        n_unconfirmed_calls=int(outcome["n_unconfirmed"]),
        tests_considered=0,
        tests_executed=0,
        n_lps=0,
        n_gapfills=0,
        cpu_seconds=0.0,
        wall_seconds=0.0,
        peak_rss_kb=0,
        case_source=case.case_source,
    )


def rescore_genome(
    genome_id: str,
    payloads: Sequence[Mapping[str, object]],
    inputs: GenomeInputs,
    tables: ReferenceTables,
    swissprot_tsv: Optional[Path],
    arm: str = ARM_OF_RECORD,
    seeds: Sequence[int] = SHUFFLE_SEEDS,
) -> Dict[str, List[CaseArmEvaluation]]:
    """Re-assign every stored case of one genome under all rules and controls; returns evaluations per arm label."""
    substrate_by_test = tables.substrate_by_test()
    rxn_cands = candidates_by_reaction(inputs.candidates, tables.ec_to_reactions)
    ctx = SyntenyContext(inputs.layout, inputs.candidates, tables.ec_modules, inputs.transporter_proteins_by_cpd)
    shuffled_layout_ctx = {
        seed: SyntenyContext(inputs.layout.shuffled(seed) if inputs.layout else None, inputs.candidates, tables.ec_modules, inputs.transporter_proteins_by_cpd)
        for seed in seeds
    }
    model_genes: List[str] = []
    for p in payloads:
        model_genes.extend(str(g) for g in p.get("model_genes", ()))
    shuffled_cands = {seed: permute_candidates_within_genome(inputs.candidates, seed=seed, extra_proteins=model_genes) for seed in seeds}
    shuffled_rxn_cands = {seed: candidates_by_reaction(c, tables.ec_to_reactions) for seed, c in shuffled_cands.items()}
    shuffled_ctx = {
        seed: SyntenyContext(inputs.layout, shuffled_cands[seed], tables.ec_modules, inputs.transporter_proteins_by_cpd) for seed in seeds
    }

    out: Dict[str, List[CaseArmEvaluation]] = {}
    for payload in payloads:
        known = load_known_experimental_ecs(payload, swissprot_tsv, genome_id)
        allowed = allowed_reactions_by_protein(known, tables.ec_to_reactions)
        for cb in payload.get("cases", ()):  # type: ignore[union-attr]
            cd = cb["case"]
            case = TruthCase(
                case_id=str(cd["case_id"]),
                genome_id=str(cd["genome_id"]),
                protein_id=str(cd["protein_id"]),
                truth_ecs=tuple(str(x) for x in cd["truth_ecs"]),
                case_source=str(cd.get("case_source") or ("swissprot" if str(cd["case_id"]).startswith("SP:") else "price148")),
            )
            stratum = str(cb["stratum_info"]["stratum"])
            _, target_rxns = resolve_equivalent_ecs_and_reactions(case.truth_ecs, tables.ec_to_reactions)
            target = {base_seed_rxn_id(r) for r in target_rxns}
            tests = stored_tests_by_arm(cb["arm_results"])
            added = added_reactions_with_tests(t for a in ("A1", "A2", "A3") for t in tests.get(a, ()))
            # The stored union arm lists the same reactions; use it as the authoritative set when present.
            union_added = {base_seed_rxn_id(str(r)) for r in cb["arm_results"].get(arm, {}).get("added_reactions", ())}
            for r in union_added:
                added.setdefault(r, set())
            added = {r: t for r, t in added.items() if r in union_added} if union_added else added

            for rule in RULES:
                cfg = rule_config(rule)
                calls = assign_reactions(added, rxn_cands, tables.rxn_ecs, substrate_by_test, ctx, cfg, allowed)
                oc = case_outcome(calls, case.protein_id, target)
                out.setdefault(rule, []).append(_evaluation(case, stratum, rule, calls, oc, target))
            for seed in seeds:
                for rule in ("R1", "R3"):
                    cfg = rule_config(rule)
                    calls = assign_reactions(added, shuffled_rxn_cands[seed], tables.rxn_ecs, substrate_by_test, shuffled_ctx[seed], cfg, allowed)
                    oc = case_outcome(calls, case.protein_id, target)
                    label = f"{rule}_shuffled_genes_seed{seed}"
                    out.setdefault(label, []).append(_evaluation(case, stratum, label, calls, oc, target))
                cfg = rule_config("R1")
                calls = assign_reactions(added, rxn_cands, tables.rxn_ecs, substrate_by_test, shuffled_layout_ctx[seed], cfg, allowed)
                oc = case_outcome(calls, case.protein_id, target)
                label = f"R1_shuffled_layout_seed{seed}"
                out.setdefault(label, []).append(_evaluation(case, stratum, label, calls, oc, target))
    return out


def _genome_worker(task: Tuple[str, Dict[str, str]]) -> Tuple[str, Dict[str, List[Dict[str, object]]], Dict[str, object]]:
    genome_id, a = task
    import time as _time
    _t0 = _time.time()
    print(json.dumps({"genome_id": genome_id, "status": "start"}), flush=True)
    tables = ReferenceTables.load(
        Path(a["seed_reactions"]), Path(a["ec_aliases"]) if a.get("ec_aliases") else None, None, Path(a["kegg_modules"]), Path(a["panel_tsv"])
    )
    payloads = []
    for tar in (Path(a["genomes_dir"]) / f"{genome_id}.outputs.tar", Path(a["genomes_dir"]) / f"{genome_id}.swissprot.outputs.tar"):
        if tar.exists():
            payloads.append(read_genome_results(tar))
    with tempfile.TemporaryDirectory(prefix=f"exp1-{genome_id}-", dir=a.get("tmp_dir") or None) as td:
        trans = extract_model_member(Path(a["models_tar"]), genome_id, "-Transporter.tbl", Path(td))
        inputs = GenomeInputs.load(
            genome_id,
            Path(a["layout_dir"]) / f"{genome_id}.layout.tsv",
            Path(a["candidates_dir"]) / f"{genome_id}.candidates.tsv",
            trans,
        )
        evals = rescore_genome(genome_id, payloads, inputs, tables, Path(a["swissprot_tsv"]) if a.get("swissprot_tsv") else None)
    print(json.dumps({"genome_id": genome_id, "status": "done", "seconds": round(_time.time() - _t0, 1)}), flush=True)
    info = {
        "genome_id": genome_id,
        "n_payloads": len(payloads),
        "n_cases": len(evals.get("R0", ())),
        "layout_proteins": inputs.layout.n_proteins() if inputs.layout else 0,
        "candidates": len(inputs.candidates),
        "transporter_substrates": len(inputs.transporter_proteins_by_cpd),
    }
    return genome_id, {k: [asdict(e) for e in v] for k, v in evals.items()}, info


def summarize(evaluations_by_arm: Mapping[str, Sequence[CaseArmEvaluation]], n_resamples: int) -> Dict[str, object]:
    summaries = summarize_arms_with_bootstrap(evaluations_by_arm, n_resamples=n_resamples, seed=BOOTSTRAP_SEED)
    r0 = {e.case_id: e for e in evaluations_by_arm["R0"]}
    readded = {cid for cid, e in r0.items() if e.reaction_recovered}
    per_rule: Dict[str, object] = {}
    for label, evals in evaluations_by_arm.items():
        by_case = {e.case_id: e for e in evals}
        correct_among_readded = sum(1 for cid in readded if by_case[cid].recovered)
        b = sum(1 for cid, e in by_case.items() if e.recovered and not r0[cid].recovered)
        c = sum(1 for cid, e in by_case.items() if not e.recovered and r0[cid].recovered)
        s = summaries[label]
        per_rule[label] = {
            "n_cases": s.n_cases_total,
            "n_reaction_readded": len(readded),
            "correct_among_readded": correct_among_readded,
            "n_recovered": s.n_recovered_total,
            "recovery_ci": list(s.recovery_ci_total),
            "n_recovered_off_growth": s.n_recovered_off_growth,
            "n_cases_off_growth": s.n_cases_off_growth,
            "n_recovered_price148": s.n_recovered_price148,
            "n_recovered_swissprot": s.n_recovered_swissprot,
            "n_gene_calls": s.n_gene_calls,
            "n_no_gene_calls": sum(1 for e in evals for cl in e.calls if cl.call_kind == "no_gene"),
            "n_reaction_only_calls": s.n_reaction_only_calls,
            "n_supported": s.n_supported_calls,
            "n_harmful": s.n_harmful_calls,
            "n_unconfirmed": s.n_unconfirmed_calls,
            "harm_rate": s.harm_rate,
            "harm_ci": list(s.harm_ci),
            "vs_R0_recovered_only_here": b,
            "vs_R0_lost": c,
            "vs_R0_mcnemar_one_sided_p": exact_one_sided_mcnemar_p(b, c),
        }
    shuffled_max = {
        rule: max(per_rule[f"{rule}_shuffled_genes_seed{s}"]["n_recovered"] for s in SHUFFLE_SEEDS if f"{rule}_shuffled_genes_seed{s}" in per_rule)
        for rule in ("R1", "R3")
        if any(f"{rule}_shuffled_genes_seed{s}" in per_rule for s in SHUFFLE_SEEDS)
    }
    layout_ctrl = [per_rule[f"R1_shuffled_layout_seed{s}"]["n_recovered"] for s in SHUFFLE_SEEDS if f"R1_shuffled_layout_seed{s}" in per_rule]
    criterion = {}
    for rule in ("R1", "R3"):
        pr = per_rule.get(rule)
        if pr is None:
            continue
        criterion[rule] = {
            "correct_among_readded": pr["correct_among_readded"],
            "required": CRITERION_MIN_CORRECT,
            "harm_rate": pr["harm_rate"],
            "harm_limit": CRITERION_MAX_HARM,
            "shuffled_genes_max_recovered": shuffled_max.get(rule),
            "shuffled_limit": CRITERION_MAX_SHUFFLED,
            "met": bool(
                pr["correct_among_readded"] >= CRITERION_MIN_CORRECT
                and pr["harm_rate"] <= CRITERION_MAX_HARM
                and (shuffled_max.get(rule, 0) or 0) <= CRITERION_MAX_SHUFFLED
            ),
        }
    return {
        "per_rule": per_rule,
        "criterion": criterion,
        "controls": {"shuffled_genes_max_recovered": shuffled_max, "shuffled_layout_R1_recovered": layout_ctrl, "R0_recovered": per_rule["R0"]["n_recovered"]},
        "arm_summaries": {k: asdict(v) for k, v in summaries.items()},
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--genomes-dir", type=Path, required=True, help="offgrowth-v1/genomes with <genome>.outputs.tar")
    p.add_argument("--candidates-dir", type=Path, required=True)
    p.add_argument("--layout-dir", type=Path, required=True)
    p.add_argument("--models-tar", type=Path, required=True)
    p.add_argument("--seed-reactions", type=Path, required=True)
    p.add_argument("--ec-aliases", type=Path, default=None)
    p.add_argument("--kegg-modules", type=Path, required=True)
    p.add_argument("--panel-tsv", type=Path, required=True)
    p.add_argument("--swissprot-tsv", type=Path, default=None)
    p.add_argument("--genome-ids", nargs="*", default=None)
    p.add_argument("--workers", type=int, default=max(1, min(22, int(os.environ.get("SLURM_CPUS_PER_TASK", "4")))))
    p.add_argument("--n-resamples", type=int, default=BOOTSTRAP_RESAMPLES)
    p.add_argument("--tmp-dir", type=Path, default=None)
    p.add_argument("--out-dir", type=Path, required=True)
    a = p.parse_args(argv)

    genome_ids = a.genome_ids or sorted({t.name.split(".outputs.tar")[0].replace(".swissprot", "") for t in a.genomes_dir.glob("*.outputs.tar")})
    args = {
        "genomes_dir": str(a.genomes_dir), "candidates_dir": str(a.candidates_dir), "layout_dir": str(a.layout_dir),
        "models_tar": str(a.models_tar), "seed_reactions": str(a.seed_reactions), "ec_aliases": str(a.ec_aliases) if a.ec_aliases else "",
        "kegg_modules": str(a.kegg_modules), "panel_tsv": str(a.panel_tsv), "swissprot_tsv": str(a.swissprot_tsv) if a.swissprot_tsv else "",
        "tmp_dir": str(a.tmp_dir) if a.tmp_dir else "",
    }
    tasks = [(g, args) for g in genome_ids]
    evaluations_by_arm: Dict[str, List[CaseArmEvaluation]] = {}
    infos = []
    if a.workers > 1 and len(tasks) > 1:
        with ProcessPoolExecutor(max_workers=min(a.workers, len(tasks))) as pool:
            results = list(pool.map(_genome_worker, tasks))
    else:
        results = [_genome_worker(t) for t in tasks]
    from offgrowth.runner import _reconstruct_evaluation

    for genome_id, evals, info in results:
        infos.append(info)
        print(json.dumps(info), flush=True)
        for label, items in evals.items():
            evaluations_by_arm.setdefault(label, []).extend(_reconstruct_evaluation(d) for d in items)
    result = summarize(evaluations_by_arm, a.n_resamples)
    result["genomes"] = infos
    result["receipt"] = environment_receipt({
        "panel_tsv": a.panel_tsv, "kegg_modules": a.kegg_modules, "seed_reactions": a.seed_reactions,
        **({"swissprot_tsv": a.swissprot_tsv} if a.swissprot_tsv else {}),
    })
    a.out_dir.mkdir(parents=True, exist_ok=True)
    write_json(a.out_dir / "exp1_results.json", result)
    rows = []
    for label, evals in evaluations_by_arm.items():
        for e in evals:
            rows.append({
                "arm": label, "case_id": e.case_id, "genome_id": e.genome_id, "protein_id": e.protein_id, "case_source": e.case_source,
                "stratum": e.stratum, "reaction_recovered": int(e.reaction_recovered), "recovered": int(e.recovered),
                "n_gene_calls": e.n_gene_calls, "n_supported": e.n_supported_calls, "n_harmful": e.n_harmful_calls, "n_unconfirmed": e.n_unconfirmed_calls,
                "n_no_gene": sum(1 for c in e.calls if c.call_kind == "no_gene"),
                "truth_calls": ";".join(f"{c.reaction_id}:{c.protein_id}:{c.score:.3f}" for c in e.calls if c.reaction_id in e.target_reactions),
            })
    write_tsv(a.out_dir / "exp1_cases.tsv", rows, list(rows[0].keys()) if rows else ["arm"])
    print(json.dumps({"criterion": result["criterion"], "controls": result["controls"]}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
