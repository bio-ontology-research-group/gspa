#!/usr/bin/env python3
"""Experiment 3: per-genome condition panels from transporter calls, re-scored on the stored sole-source tests (design §3).

Arms: P0 = full Biolog panel (the pilot's A1 tests); P1 = wells whose SEED compound has a gapsmith transporter hit in
the genome (genomes with fewer than 10 such wells keep all wells and are flagged); P2 = size-matched random panels
(seeds 1–5). Halves by well-column parity are reported for P0 and P1 as a stability check. Panels are written and
hashed before scoring. Calls use the baseline rule R0 (today's rule). No new linear programs.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict
import json
import os
from pathlib import Path
import random
import tempfile
from typing import Dict, List, Mapping, Optional, Sequence, Set, Tuple

from .assignment import case_outcome
from .common import (
    SHUFFLE_SEEDS,
    environment_receipt,
    load_known_experimental_ecs,
    read_genome_results,
    sha256_text,
    stored_tests_by_arm,
    write_json,
    write_tsv,
)
from .exp1_assignment import _evaluation
from .inputs import GenomeInputs, ReferenceTables, extract_model_member
from offgrowth.mask import TruthCase, base_seed_rxn_id, resolve_equivalent_ecs_and_reactions
from offgrowth.score import (
    BOOTSTRAP_RESAMPLES,
    BOOTSTRAP_SEED,
    CaseArmEvaluation,
    assign_calls_for_added_reactions,
    exact_one_sided_mcnemar_p,
    permute_candidates_within_genome,
    summarize_arms_with_bootstrap,
)

MIN_PANEL_WELLS = 10
CRITERION_MAX_TEST_FRACTION = 0.25
CRITERION_MIN_PRICE_RECOVERED = 9


def well_column(well: str) -> int:
    digits = "".join(ch for ch in well if ch.isdigit())
    return int(digits) if digits else 0


def build_panels(tables: ReferenceTables, transporter_cpds: Set[str], seeds: Sequence[int] = SHUFFLE_SEEDS) -> Dict[str, object]:
    """Panels as sets of panel ids; built from label-free inputs only."""
    all_ids = [p.panel_id for p in tables.panel]
    p1 = [p.panel_id for p in tables.panel if p.seed_id in transporter_cpds]
    flagged = len(p1) < MIN_PANEL_WELLS
    if flagged:
        p1 = list(all_ids)
    panels: Dict[str, List[str]] = {"P0": all_ids, "P1": p1}
    for seed in seeds:
        rng = random.Random(seed)
        panels[f"P2_seed{seed}"] = sorted(rng.sample(all_ids, len(p1)))
    for label in ("P0", "P1"):
        panels[f"{label}_odd"] = [pid for pid in panels[label] if well_column(tables.panel_by_id[pid].well) % 2 == 1]
        panels[f"{label}_even"] = [pid for pid in panels[label] if well_column(tables.panel_by_id[pid].well) % 2 == 0]
    return {"panels": {k: sorted(v) for k, v in panels.items()}, "p1_flagged_small": flagged, "n_transporter_cpds": len(transporter_cpds)}


def rescore_genome(
    genome_id: str,
    payloads: Sequence[Mapping[str, object]],
    inputs: GenomeInputs,
    tables: ReferenceTables,
    swissprot_tsv: Optional[Path],
    seeds: Sequence[int] = SHUFFLE_SEEDS,
) -> Tuple[Dict[str, List[CaseArmEvaluation]], Dict[str, object], Dict[str, int]]:
    panel_info = build_panels(tables, set(inputs.transporter_proteins_by_cpd), seeds)
    panels: Dict[str, List[str]] = panel_info["panels"]  # type: ignore[assignment]
    model_genes = [str(g) for p in payloads for g in p.get("model_genes", ())]
    shuffled_cands = {seed: permute_candidates_within_genome(inputs.candidates, seed=seed, extra_proteins=model_genes) for seed in seeds}
    out: Dict[str, List[CaseArmEvaluation]] = {}
    tests_per_arm: Dict[str, int] = {}
    for payload in payloads:
        known = load_known_experimental_ecs(payload, swissprot_tsv, genome_id)
        for cb in payload.get("cases", ()):  # type: ignore[union-attr]
            cd = cb["case"]
            case = TruthCase(
                case_id=str(cd["case_id"]), genome_id=str(cd["genome_id"]), protein_id=str(cd["protein_id"]),
                truth_ecs=tuple(str(x) for x in cd["truth_ecs"]),
                case_source=str(cd.get("case_source") or ("swissprot" if str(cd["case_id"]).startswith("SP:") else "price148")),
            )
            stratum = str(cb["stratum_info"]["stratum"])
            _, target_rxns = resolve_equivalent_ecs_and_reactions(case.truth_ecs, tables.ec_to_reactions)
            target = {base_seed_rxn_id(r) for r in target_rxns}
            a1_tests = stored_tests_by_arm(cb["arm_results"]).get("A1", [])
            by_panel_id = {t["test_id"][len("A1_"):]: t for t in a1_tests}
            for label, ids in panels.items():
                selected = [by_panel_id[pid] for pid in ids if pid in by_panel_id]
                tests_per_arm[label] = tests_per_arm.get(label, 0) + len(selected)
                added = sorted({base_seed_rxn_id(str(r)) for t in selected for r in t["added_reactions"]})
                calls = assign_calls_for_added_reactions(added, inputs.candidates, tables.ec_to_reactions, known)
                oc = case_outcome(calls, case.protein_id, target)
                out.setdefault(label, []).append(_evaluation(case, stratum, label, calls, oc, target))
                if label == "P1":
                    for seed in seeds:
                        calls_s = assign_calls_for_added_reactions(added, shuffled_cands[seed], tables.ec_to_reactions, known)
                        oc_s = case_outcome(calls_s, case.protein_id, target)
                        lab = f"P1_shuffled_genes_seed{seed}"
                        out.setdefault(lab, []).append(_evaluation(case, stratum, lab, calls_s, oc_s, target))
    return out, panel_info, tests_per_arm


def _genome_worker(task: Tuple[str, Dict[str, str]]):
    genome_id, a = task
    import time as _time
    _t0 = _time.time()
    print(json.dumps({"genome_id": genome_id, "status": "start"}), flush=True)
    tables = ReferenceTables.load(Path(a["seed_reactions"]), Path(a["ec_aliases"]) if a.get("ec_aliases") else None, None, Path(a["kegg_modules"]), Path(a["panel_tsv"]))
    payloads = []
    for tar in (Path(a["genomes_dir"]) / f"{genome_id}.outputs.tar", Path(a["genomes_dir"]) / f"{genome_id}.swissprot.outputs.tar"):
        if tar.exists():
            payloads.append(read_genome_results(tar))
    with tempfile.TemporaryDirectory(prefix=f"exp3-{genome_id}-", dir=a.get("tmp_dir") or None) as td:
        trans = extract_model_member(Path(a["models_tar"]), genome_id, "-Transporter.tbl", Path(td))
        inputs = GenomeInputs.load(genome_id, None, Path(a["candidates_dir"]) / f"{genome_id}.candidates.tsv", trans)
    evals, panel_info, tests_per_arm = rescore_genome(genome_id, payloads, inputs, tables, Path(a["swissprot_tsv"]) if a.get("swissprot_tsv") else None)
    print(json.dumps({"genome_id": genome_id, "status": "done", "seconds": round(_time.time() - _t0, 1)}), flush=True)
    return genome_id, {k: [asdict(e) for e in v] for k, v in evals.items()}, panel_info, tests_per_arm


def summarize(evaluations_by_arm: Mapping[str, Sequence[CaseArmEvaluation]], tests_per_arm: Mapping[str, int], n_resamples: int) -> Dict[str, object]:
    summaries = summarize_arms_with_bootstrap(evaluations_by_arm, n_resamples=n_resamples, seed=BOOTSTRAP_SEED)
    p0 = {e.case_id: e for e in evaluations_by_arm["P0"]}
    per_arm: Dict[str, object] = {}
    for label, evals in evaluations_by_arm.items():
        s = summaries[label]
        by_case = {e.case_id: e for e in evals}
        b = sum(1 for cid, e in by_case.items() if e.recovered and not p0[cid].recovered)
        c = sum(1 for cid, e in by_case.items() if not e.recovered and p0[cid].recovered)
        per_arm[label] = {
            "n_cases": s.n_cases_total, "n_recovered": s.n_recovered_total, "recovery_ci": list(s.recovery_ci_total),
            "n_reaction_readded": sum(1 for e in evals if e.reaction_recovered),
            "n_recovered_off_growth": s.n_recovered_off_growth, "n_cases_off_growth": s.n_cases_off_growth,
            "n_recovered_price148": s.n_recovered_price148, "n_recovered_swissprot": s.n_recovered_swissprot,
            "n_gene_calls": s.n_gene_calls, "n_supported": s.n_supported_calls, "n_harmful": s.n_harmful_calls,
            "n_unconfirmed": s.n_unconfirmed_calls, "harm_rate": s.harm_rate, "harm_ci": list(s.harm_ci),
            "tests": tests_per_arm.get(label, 0), "test_fraction_of_P0": (tests_per_arm.get(label, 0) / tests_per_arm["P0"]) if tests_per_arm.get("P0") else 0.0,
            "vs_P0_only_here": b, "vs_P0_lost": c,
            "vs_P0_mcnemar_two_sided_p": min(1.0, 2.0 * min(exact_one_sided_mcnemar_p(b, c), exact_one_sided_mcnemar_p(c, b))),
        }
    p2_mean = sum(per_arm[f"P2_seed{s}"]["n_recovered"] for s in SHUFFLE_SEEDS) / float(len(SHUFFLE_SEEDS))
    shuffled_max = max(per_arm[f"P1_shuffled_genes_seed{s}"]["n_recovered"] for s in SHUFFLE_SEEDS)
    p1, p0s = per_arm["P1"], per_arm["P0"]
    criterion = {
        "P1_recovered": p1["n_recovered"], "P0_recovered": p0s["n_recovered"], "P1_test_fraction": p1["test_fraction_of_P0"],
        "P1_price148_recovered": p1["n_recovered_price148"], "required_min_price148": CRITERION_MIN_PRICE_RECOVERED,
        "max_test_fraction": CRITERION_MAX_TEST_FRACTION,
        "P2_random_mean_recovered": p2_mean, "P1_shuffled_genes_max_recovered": shuffled_max,
        "met": bool(p1["n_recovered"] >= p0s["n_recovered"] and p1["test_fraction_of_P0"] <= CRITERION_MAX_TEST_FRACTION and p1["n_recovered_price148"] >= CRITERION_MIN_PRICE_RECOVERED),
        "P1_beats_random_control": bool(p1["n_recovered"] > p2_mean),
    }
    return {"per_arm": per_arm, "criterion": criterion, "arm_summaries": {k: asdict(v) for k, v in summaries.items()}}


def main(argv: Optional[Sequence[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--genomes-dir", type=Path, required=True)
    p.add_argument("--candidates-dir", type=Path, required=True)
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
        "genomes_dir": str(a.genomes_dir), "candidates_dir": str(a.candidates_dir), "models_tar": str(a.models_tar),
        "seed_reactions": str(a.seed_reactions), "ec_aliases": str(a.ec_aliases) if a.ec_aliases else "", "kegg_modules": str(a.kegg_modules),
        "panel_tsv": str(a.panel_tsv), "swissprot_tsv": str(a.swissprot_tsv) if a.swissprot_tsv else "", "tmp_dir": str(a.tmp_dir) if a.tmp_dir else "",
    }
    tasks = [(g, args) for g in genome_ids]
    if a.workers > 1 and len(tasks) > 1:
        with ProcessPoolExecutor(max_workers=min(a.workers, len(tasks))) as pool:
            results = list(pool.map(_genome_worker, tasks))
    else:
        results = [_genome_worker(t) for t in tasks]
    from offgrowth.runner import _reconstruct_evaluation

    a.out_dir.mkdir(parents=True, exist_ok=True)
    panels_by_genome = {g: info for g, _, info, _ in results}
    panels_text = json.dumps(panels_by_genome, indent=2, sort_keys=True) + "\n"
    (a.out_dir / "panels.json").write_text(panels_text, encoding="utf-8")
    panels_sha = sha256_text(panels_text)
    evaluations_by_arm: Dict[str, List[CaseArmEvaluation]] = {}
    tests_per_arm: Dict[str, int] = {}
    for g, evals, info, tpa in results:
        print(json.dumps({"genome_id": g, "n_P1_wells": len(info["panels"]["P1"]), "flagged": info["p1_flagged_small"]}), flush=True)
        for label, items in evals.items():
            evaluations_by_arm.setdefault(label, []).extend(_reconstruct_evaluation(d) for d in items)
        for label, n in tpa.items():
            tests_per_arm[label] = tests_per_arm.get(label, 0) + n
    result = summarize(evaluations_by_arm, tests_per_arm, a.n_resamples)
    result["panels_sha256"] = panels_sha
    result["panel_sizes"] = {g: {k: len(v) for k, v in info["panels"].items()} for g, info in panels_by_genome.items()}
    result["receipt"] = environment_receipt({"panel_tsv": a.panel_tsv, "kegg_modules": a.kegg_modules, "seed_reactions": a.seed_reactions})
    write_json(a.out_dir / "exp3_results.json", result)
    rows = [
        {"arm": label, "case_id": e.case_id, "genome_id": e.genome_id, "protein_id": e.protein_id, "case_source": e.case_source, "stratum": e.stratum,
         "reaction_recovered": int(e.reaction_recovered), "recovered": int(e.recovered), "n_gene_calls": e.n_gene_calls,
         "n_supported": e.n_supported_calls, "n_harmful": e.n_harmful_calls, "n_unconfirmed": e.n_unconfirmed_calls}
        for label, evals in evaluations_by_arm.items() for e in evals
    ]
    write_tsv(a.out_dir / "exp3_cases.tsv", rows, list(rows[0].keys()) if rows else ["arm"])
    print(json.dumps(result["criterion"], indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
