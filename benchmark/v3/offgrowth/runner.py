#!/usr/bin/env python3
"""Per-genome and cohort runner for the pre-registered off-growth enzyme discovery experiment.

Supports:
1. Running all truth cases for a single genome (`run_genome_cases`) and packing
   results into `<genome_id>.outputs.tar` to respect Shaheen scratch inode quotas.
2. Aggregating per-genome results across the 22 Price-148 genomes (`aggregate_cohort`),
   running the genome-block bootstrap (`seed=20261030`, `2,000` resamples),
   pre-registered tests H1-H3, shuffled-candidate control (`seeds=1..5`),
   unmasked reachability control, and solver replay summary, and writing the
   plain-word Markdown report.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import csv
from dataclasses import asdict, dataclass
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import sys
import tarfile
from typing import Dict, List, Mapping, Optional, Sequence, Set, Tuple

import cobra
import scipy

from .arms import (
    DEFAULT_FLUX_TOL,
    ArmRunResult,
    GapFillTestResult,
    KeggModuleDef,
    RankedTestEntry,
    SolverInfo,
    SolverReplayRecord,
    detect_solvers,
    load_gapsmith_medium_csv,
    load_kegg_modules,
    run_all_arms_for_case,
    run_all_arms_unmasked,
)
from .mask import (
    CandidateClaim,
    MaskingAuditReport,
    TransporterHit,
    TruthCase,
    base_seed_rxn_id,
    protein_aliases,
    strip_version,
    load_gapsmith_reactions_tbl,
    load_gapsmith_transporter_tbl,
    load_extra_cases,
    load_price148_cases,
    load_stage1_candidates_from_tsv,
    merge_stage1_candidates_to_tsv,
    prepare_masked_case,
)
from .panel import (
    normalize_ec,
    PanelEntry,
    SeedMetabolite,
    SeedReaction,
    UnmappedPanelEntry,
    build_ec_to_reactions_map,
    load_panel_tsv,
    load_seed_metabolites,
    load_seed_reactions,
    sha256_file,
)
from .report import write_experiment_outputs
from .score import (
    BOOTSTRAP_RESAMPLES,
    BOOTSTRAP_SEED,
    SHUFFLE_CONTROL_SEEDS,
    ArmSummaryMetrics,
    CaseArmEvaluation,
    H3EfficiencyResult,
    McNemarResult,
    ReactionGeneCall,
    ShuffledControlSummary,
    SolverReplaySummary,
    StratumClassification,
    classify_case_stratum,
    evaluate_case_arm,
    run_preregistered_tests,
    run_shuffled_candidate_control,
    summarize_arms_with_bootstrap,
    summarize_solver_replays,
)


@dataclass
class CaseExecutionBundle:
    case: TruthCase
    stratum_info: StratumClassification
    audit: MaskingAuditReport
    arm_results: Dict[str, ArmRunResult]
    evaluations: Dict[str, CaseArmEvaluation]
    masked_candidate_evaluations: Dict[str, CaseArmEvaluation]
    unmasked_arm_results: Dict[str, ArmRunResult]
    unmasked_evaluations: Dict[str, CaseArmEvaluation]


def execute_single_case(
    case: TruthCase,
    draft_model: cobra.Model,
    unmasked_filled_model: cobra.Model,
    unmasked_candidates: Sequence[CandidateClaim],
    universal_reactions: Mapping[str, SeedReaction],
    metabolites: Mapping[str, SeedMetabolite],
    ec_to_reactions: Mapping[str, Sequence[str]],
    panel_entries: Sequence[PanelEntry],
    default_medium: Mapping[str, float],
    known_experimental_ecs: Mapping[str, Sequence[str]],
    transporters: Sequence[TransporterHit] = (),
    kegg_modules: Sequence[KeggModuleDef] = (),
    replay_fraction: float = 0.05,
    flux_tol: float = DEFAULT_FLUX_TOL,
    precomputed_unmasked_arm_results: Optional[Mapping[str, ArmRunResult]] = None,
    workers: int = 1,
    max_tests_per_arm: Optional[int] = None,
) -> CaseExecutionBundle:
    """Mask, audit, classify stratum, run A0-A4 (masked + unmasked control), and score `(P, E)`."""
    stratum_info = classify_case_stratum(
        case=case,
        unmasked_filled_model=unmasked_filled_model,
        ec_to_reactions=ec_to_reactions,
        default_medium=default_medium,
        flux_tol=flux_tol,
    )

    masked_bundle = prepare_masked_case(
        case=case,
        draft_model=draft_model,
        candidates=unmasked_candidates,
        ec_to_reactions=ec_to_reactions,
        transporters=transporters,
    )

    arm_results = run_all_arms_for_case(
        masked_model=masked_bundle.model,
        audit_report=masked_bundle.audit,
        universal_reactions=universal_reactions,
        metabolites=metabolites,
        candidates=masked_bundle.candidates,
        ec_to_reactions=ec_to_reactions,
        panel_entries=panel_entries,
        default_medium=default_medium,
        transporters=masked_bundle.transporters,
        kegg_modules=kegg_modules,
        replay_fraction=replay_fraction,
        flux_tol=flux_tol,
        workers=workers,
        max_tests_per_arm=max_tests_per_arm,
    )

    if precomputed_unmasked_arm_results is not None:
        unmasked_arm_results = dict(precomputed_unmasked_arm_results)
    else:
        unmasked_arm_results = run_all_arms_unmasked(
            unmasked_draft_model=draft_model,
            universal_reactions=universal_reactions,
            metabolites=metabolites,
            unmasked_candidates=unmasked_candidates,
            ec_to_reactions=ec_to_reactions,
            panel_entries=panel_entries,
            default_medium=default_medium,
            transporters=transporters,
            kegg_modules=kegg_modules,
            replay_fraction=replay_fraction,
            flux_tol=flux_tol,
            workers=workers,
            max_tests_per_arm=max_tests_per_arm,
        )

    evaluations: Dict[str, CaseArmEvaluation] = {}
    masked_candidate_evaluations: Dict[str, CaseArmEvaluation] = {}
    unmasked_evaluations: Dict[str, CaseArmEvaluation] = {}
    for arm_name, arm_res in arm_results.items():
        evaluations[arm_name] = evaluate_case_arm(
            case=case,
            stratum=stratum_info.stratum,
            arm_result=arm_res,
            call_candidates=unmasked_candidates,
            ec_to_reactions=ec_to_reactions,
            known_experimental_ecs=known_experimental_ecs,
        )
        masked_candidate_evaluations[arm_name] = evaluate_case_arm(
            case=case,
            stratum=stratum_info.stratum,
            arm_result=arm_res,
            call_candidates=masked_bundle.candidates,
            ec_to_reactions=ec_to_reactions,
            known_experimental_ecs=known_experimental_ecs,
        )
    for arm_name, u_arm_res in unmasked_arm_results.items():
        unmasked_evaluations[arm_name] = evaluate_case_arm(
            case=case,
            stratum=stratum_info.stratum,
            arm_result=u_arm_res,
            call_candidates=unmasked_candidates,
            ec_to_reactions=ec_to_reactions,
            known_experimental_ecs=known_experimental_ecs,
        )

    return CaseExecutionBundle(
        case=case,
        stratum_info=stratum_info,
        audit=masked_bundle.audit,
        arm_results=arm_results,
        evaluations=evaluations,
        masked_candidate_evaluations=masked_candidate_evaluations,
        unmasked_arm_results=unmasked_arm_results,
        unmasked_evaluations=unmasked_evaluations,
    )


def run_cohort_in_memory(
    cases: Sequence[TruthCase],
    draft_models_by_genome: Mapping[str, cobra.Model],
    filled_models_by_genome: Mapping[str, cobra.Model],
    candidates_by_genome: Mapping[str, Sequence[CandidateClaim]],
    universal_reactions: Mapping[str, SeedReaction],
    metabolites: Mapping[str, SeedMetabolite],
    ec_to_reactions: Mapping[str, Sequence[str]],
    panel_entries: Sequence[PanelEntry],
    default_media_by_genome: Mapping[str, Mapping[str, float]],
    known_experimental_ecs_by_genome: Mapping[str, Mapping[str, Sequence[str]]],
    transporters_by_genome: Optional[Mapping[str, Sequence[TransporterHit]]] = None,
    kegg_modules: Sequence[KeggModuleDef] = (),
    unmapped_panel_entries: Sequence[UnmappedPanelEntry] = (),
    replay_fraction: float = 0.05,
    n_resamples: int = BOOTSTRAP_RESAMPLES,
    bootstrap_seed: int = BOOTSTRAP_SEED,
    shuffle_seeds: Sequence[int] = SHUFFLE_CONTROL_SEEDS,
    out_dir: Optional[Path] = None,
    workers: int = 1,
    max_tests_per_arm: Optional[int] = None,
) -> Dict[str, object]:
    """Run all cases across genomes, compute all metrics/tests/controls, and optionally write outputs."""
    solver_info = detect_solvers()
    evaluations_by_arm: Dict[str, List[CaseArmEvaluation]] = {}
    masked_cand_by_arm: Dict[str, List[CaseArmEvaluation]] = {}
    unmasked_by_arm: Dict[str, List[CaseArmEvaluation]] = {}
    arm_results_by_case: Dict[str, Dict[str, ArmRunResult]] = {}
    strata_by_case: Dict[str, str] = {}
    all_replays: List[SolverReplayRecord] = []

    # Precompute the §8 unmasked arm run once per genome
    unmasked_runs_by_genome: Dict[str, Dict[str, ArmRunResult]] = {}
    for gid in sorted({c.genome_id for c in cases}):
        unmasked_runs_by_genome[gid] = run_all_arms_unmasked(
            unmasked_draft_model=draft_models_by_genome[gid],
            universal_reactions=universal_reactions,
            metabolites=metabolites,
            unmasked_candidates=candidates_by_genome.get(gid, ()),
            ec_to_reactions=ec_to_reactions,
            panel_entries=panel_entries,
            default_medium=default_media_by_genome[gid],
            transporters=(transporters_by_genome or {}).get(gid, ()),
            kegg_modules=kegg_modules,
            replay_fraction=replay_fraction,
            workers=workers,
            max_tests_per_arm=max_tests_per_arm,
        )

    for case in cases:
        gid = case.genome_id
        bundle = execute_single_case(
            case=case,
            draft_model=draft_models_by_genome[gid],
            unmasked_filled_model=filled_models_by_genome[gid],
            unmasked_candidates=candidates_by_genome.get(gid, ()),
            universal_reactions=universal_reactions,
            metabolites=metabolites,
            ec_to_reactions=ec_to_reactions,
            panel_entries=panel_entries,
            default_medium=default_media_by_genome[gid],
            known_experimental_ecs=known_experimental_ecs_by_genome.get(gid, {}),
            transporters=(transporters_by_genome or {}).get(gid, ()),
            kegg_modules=kegg_modules,
            replay_fraction=replay_fraction,
            precomputed_unmasked_arm_results=unmasked_runs_by_genome[gid],
            workers=workers,
            max_tests_per_arm=max_tests_per_arm,
        )
        strata_by_case[case.case_id] = bundle.stratum_info.stratum
        arm_results_by_case[case.case_id] = bundle.arm_results

        for arm_name, ev in bundle.evaluations.items():
            evaluations_by_arm.setdefault(arm_name, []).append(ev)
        for arm_name, mev in bundle.masked_candidate_evaluations.items():
            masked_cand_by_arm.setdefault(arm_name, []).append(mev)
        for arm_name, uev in bundle.unmasked_evaluations.items():
            unmasked_by_arm.setdefault(arm_name, []).append(uev)
        for arm_name in ("A0", "A1", "A2", "A3"):
            if arm_name in bundle.arm_results:
                all_replays.extend(bundle.arm_results[arm_name].replays)

    arm_summaries = summarize_arms_with_bootstrap(
        evaluations_by_arm, n_resamples=n_resamples, seed=bootstrap_seed
    )
    masked_candidate_summaries = summarize_arms_with_bootstrap(
        masked_cand_by_arm, n_resamples=0, seed=bootstrap_seed
    )
    unmasked_summaries = summarize_arms_with_bootstrap(
        unmasked_by_arm, n_resamples=0, seed=bootstrap_seed
    )
    h1_res, h2_res, h3_res = run_preregistered_tests(
        evaluations_by_arm, n_resamples=n_resamples, seed=bootstrap_seed
    )

    extra_proteins: Dict[str, List[str]] = {}
    for gid, dmodel in draft_models_by_genome.items():
        extra_proteins[gid] = [g.id for g in dmodel.genes]

    shuffled_ctrl = run_shuffled_candidate_control(
        cases=cases,
        strata_by_case=strata_by_case,
        arm_results_by_case=arm_results_by_case,
        unmasked_candidates_by_genome=candidates_by_genome,
        ec_to_reactions=ec_to_reactions,
        known_experimental_ecs_by_genome=known_experimental_ecs_by_genome,
        arm_name="A1_A2_A3",
        seeds=shuffle_seeds,
        extra_proteins_by_genome=extra_proteins,
    )
    replay_summary = summarize_solver_replays(all_replays)

    report_path: Optional[Path] = None
    if out_dir is not None:
        report_path = write_experiment_outputs(
            out_dir=out_dir,
            arm_summaries=arm_summaries,
            evaluations_by_arm=evaluations_by_arm,
            h1_result=h1_res,
            h2_results=h2_res,
            h3_result=h3_res,
            shuffled_control=shuffled_ctrl,
            replay_summary=replay_summary,
            solver_info=solver_info,
            unmapped_panel_entries=unmapped_panel_entries,
            unmasked_summary=unmasked_summaries,
            masked_candidate_summary=masked_candidate_summaries,
        )

    return {
        "solver_info": solver_info,
        "arm_summaries": arm_summaries,
        "masked_candidate_summaries": masked_candidate_summaries,
        "unmasked_summaries": unmasked_summaries,
        "evaluations_by_arm": evaluations_by_arm,
        "masked_candidate_evaluations_by_arm": masked_cand_by_arm,
        "unmasked_evaluations_by_arm": unmasked_by_arm,
        "h1": h1_res,
        "h2": h2_res,
        "h3": h3_res,
        "shuffled_control": shuffled_ctrl,
        "solver_replay": replay_summary,
        "report_path": report_path,
    }


def build_environment_receipt(input_files: Mapping[str, Path]) -> Dict[str, object]:
    """Build reproducibility receipt with input SHA-256 hashes and software versions."""
    solver_info = detect_solvers()
    hashes: Dict[str, str] = {}
    for label, path in sorted(input_files.items()):
        if path.exists() and path.is_file():
            hashes[label] = sha256_file(path)
    return {
        "python_version": sys.version,
        "platform": platform.platform(),
        "cobra_version": getattr(cobra, "__version__", "unknown"),
        "scipy_version": getattr(scipy, "__version__", "unknown"),
        "solver_info": asdict(solver_info),
        "input_sha256": hashes,
    }


def pack_directory_to_tar(source_dir: Path, tar_path: Path) -> None:
    """Pack `source_dir` into `tar_path` and remove `source_dir` to conserve scratch inodes."""
    tar_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_tar = tar_path.with_suffix(tar_path.suffix + ".partial")
    with tarfile.open(tmp_tar, "w") as tf:
        tf.add(source_dir, arcname=source_dir.name)
    tmp_tar.replace(tar_path)
    shutil.rmtree(source_dir)


def load_swissprot_experimental_ecs(
    labels_tsv: Optional[Path],
    genome_id: str,
) -> Dict[str, List[str]]:
    """Load Swiss-Prot experimental EC annotations for `genome_id` if present."""
    by_prot: Dict[str, List[str]] = {}
    if labels_tsv is None or not labels_tsv.exists():
        return by_prot
    import csv

    with labels_tsv.open("r", encoding="utf-8", errors="replace", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            gid = (row.get("genome_id") or row.get("host_assembly") or "").strip()
            if gid and gid != genome_id:
                continue
            tier = (row.get("tier") or "experimental").strip().lower()
            if tier and tier != "experimental":
                continue
            pid = (row.get("protein_id") or "").strip()
            fid = (row.get("function_id") or row.get("ec") or "").strip()
            if pid and fid:
                by_prot.setdefault(pid, []).append(fid)
    return by_prot


def resolve_collection_stage1_tsvs(
    collection_dir: Path,
    genome_id: str,
) -> Tuple[Tuple[Path, ...], Tuple[Path, ...]]:
    """Resolve the 6 predictors' `stage1` and `stage0` `.integrated.tsv` paths from `plan.json`."""
    plan_path = collection_dir / "plan.json"
    if not plan_path.exists():
        raise FileNotFoundError(f"Collection plan.json not found at '{plan_path}'")
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    cells = [
        c for c in plan.get("cells", []) if str(c.get("assembly", "")).strip() == genome_id
    ]
    cells.sort(key=lambda c: (str(c.get("predictor", "")), int(c.get("index", 0))))

    stage1_paths: List[Path] = []
    stage0_paths: List[Path] = []
    for c in cells:
        idx = c["index"]
        attempt_dir = collection_dir / "attempts" / str(idx)
        s1_matches = sorted(attempt_dir.glob(f"job-*/production/stage1/{genome_id}.integrated.tsv"))
        if s1_matches:
            stage1_paths.append(s1_matches[-1])
        s0_matches = sorted(attempt_dir.glob(f"job-*/production/stage0/{genome_id}.integrated.tsv"))
        if s0_matches:
            stage0_paths.append(s0_matches[-1])

    if not stage1_paths:
        # Fallback for flat stage1 directory (e.g. local extracted samples)
        flat_matches = sorted(collection_dir.glob(f"*{genome_id}.integrated.tsv"))
        stage1_paths.extend(flat_matches)

    if not stage1_paths:
        raise FileNotFoundError(
            f"No stage1 .integrated.tsv files found for genome '{genome_id}' under '{collection_dir}'"
        )
    return tuple(stage1_paths), tuple(stage0_paths)


def _prepass_single_genome_worker(
    task: Tuple[str, Path, Path, Optional[int], Tuple[str, ...]],
) -> Dict[str, object]:
    genome_id, collection_dir, out_dir, max_per_ec, keep_proteins = task
    s1_paths, s0_paths = resolve_collection_stage1_tsvs(collection_dir, genome_id)
    out_tsv = out_dir / f"{genome_id}.candidates.tsv"
    out_receipt = out_dir / f"{genome_id}.candidates.receipt.json"
    receipt = merge_stage1_candidates_to_tsv(
        tsv_paths=s1_paths,
        out_tsv=out_tsv,
        max_per_ec=max_per_ec,
        keep_proteins=keep_proteins,
        stage0_paths=s0_paths,
        out_receipt_json=out_receipt,
    )
    receipt["genome_id"] = genome_id
    return receipt


def run_prepass_cli(args: argparse.Namespace) -> int:
    """Stream predictor stage1 TSVs and write compact candidate TSV(s) per genome."""
    args.out_dir.mkdir(parents=True, exist_ok=True)
    all_cases = load_price148_cases(args.manifest) if args.manifest else ()
    keep_by_genome: Dict[str, List[str]] = {}
    for c in all_cases:
        # Supervisor fix 2026-10-30: truth proteins must not get extra candidates (label-dependent input).
        pass

    if args.stage1_tsvs:
        if not args.genome_id:
            raise ValueError("--genome-id is required when --stage1-tsvs is passed directly.")
        out_tsv = args.out_dir / f"{args.genome_id}.candidates.tsv"
        out_receipt = args.out_dir / f"{args.genome_id}.candidates.receipt.json"
        merge_stage1_candidates_to_tsv(
            tsv_paths=args.stage1_tsvs,
            out_tsv=out_tsv,
            max_per_ec=args.max_per_ec,
            keep_proteins=tuple(keep_by_genome.get(args.genome_id, ())),
            out_receipt_json=out_receipt,
        )
        return 0

    if not args.collection_dir:
        raise ValueError("Provide either --collection-dir or --stage1-tsvs.")

    if args.genome_id:
        genome_ids = [args.genome_id]
    elif all_cases:
        genome_ids = sorted({c.genome_id for c in all_cases})
    else:
        plan = json.loads((args.collection_dir / "plan.json").read_text(encoding="utf-8"))
        genome_ids = sorted({str(c["assembly"]) for c in plan.get("cells", [])})

    tasks = [
        (
            gid,
            args.collection_dir,
            args.out_dir,
            args.max_per_ec,
            tuple(sorted(set(keep_by_genome.get(gid, ())))),
        )
        for gid in genome_ids
    ]

    workers = max(1, int(args.workers))
    if workers > 1 and len(tasks) > 1:
        with ProcessPoolExecutor(max_workers=min(workers, len(tasks))) as pool:
            receipts = list(pool.map(_prepass_single_genome_worker, tasks))
    else:
        receipts = [_prepass_single_genome_worker(t) for t in tasks]

    summary_path = args.out_dir / "prepass_summary.json"
    summary_path.write_text(
        json.dumps(
            {
                "n_genomes": len(receipts),
                "max_per_ec": args.max_per_ec,
                "genomes": receipts,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    return 0


def write_timing_artifacts(
    work_dir: Path,
    genome_id: str,
    bundles: Sequence[CaseExecutionBundle],
    unmasked_arm_results: Mapping[str, ArmRunResult],
) -> Tuple[Path, Path]:
    """Write `timing_per_test.tsv` and `timing_summary.json` for pilot sizing."""
    tsv_path = work_dir / "timing_per_test.tsv"
    json_path = work_dir / "timing_summary.json"

    fieldnames = [
        "genome_id",
        "case_id",
        "run_kind",
        "arm",
        "test_id",
        "target_id",
        "feasible_in_universal",
        "already_satisfied_in_draft",
        "gapfill_ran",
        "n_lps",
        "n_gapfills",
        "cpu_seconds",
        "wall_seconds",
        "n_added_reactions",
    ]
    rows: List[Dict[str, object]] = []
    per_arm_totals: Dict[str, Dict[str, float]] = {}

    def _record_arm_map(case_id: str, run_kind: str, arm_map: Mapping[str, ArmRunResult]) -> None:
        for arm_name, arm_res in sorted(arm_map.items()):
            agg = per_arm_totals.setdefault(
                f"{run_kind}:{arm_name}",
                {
                    "n_runs": 0,
                    "tests_considered": 0,
                    "tests_executed": 0,
                    "n_lps": 0,
                    "n_gapfills": 0,
                    "cpu_seconds": 0.0,
                    "wall_seconds": 0.0,
                },
            )
            agg["n_runs"] += 1
            agg["tests_considered"] += arm_res.tests_considered
            agg["tests_executed"] += arm_res.tests_executed
            agg["n_lps"] += arm_res.n_lps
            agg["n_gapfills"] += arm_res.n_gapfills
            agg["cpu_seconds"] += arm_res.cpu_seconds
            agg["wall_seconds"] += arm_res.wall_seconds

            if arm_name in ("A0", "A1", "A2", "A3"):
                for tr in arm_res.test_results:
                    rows.append(
                        {
                            "genome_id": genome_id,
                            "case_id": case_id,
                            "run_kind": run_kind,
                            "arm": arm_name,
                            "test_id": tr.test_id,
                            "target_id": tr.target_id,
                            "feasible_in_universal": int(tr.feasible_in_universal),
                            "already_satisfied_in_draft": int(tr.already_satisfied_in_draft),
                            "gapfill_ran": int(tr.gapfill_ran),
                            "n_lps": tr.n_lps,
                            "n_gapfills": tr.n_gapfills,
                            "cpu_seconds": f"{tr.cpu_seconds:.6f}",
                            "wall_seconds": f"{tr.wall_seconds:.6f}",
                            "n_added_reactions": len(tr.added_reactions),
                        }
                    )

    _record_arm_map("UNMASKED_CONTROL", "unmasked", unmasked_arm_results)
    for b in bundles:
        _record_arm_map(b.case.case_id, "masked", b.arm_results)

    with tsv_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for r in rows:
            writer.writerow(r)

    total_lps = sum(
        int(v["n_lps"])
        for k, v in per_arm_totals.items()
        if k.endswith(":A0") or k.endswith(":A1_A2_A3")
    )
    total_cpu = sum(
        float(v["cpu_seconds"])
        for k, v in per_arm_totals.items()
        if k.endswith(":A0") or k.endswith(":A1_A2_A3")
    )
    total_wall = sum(
        float(v["wall_seconds"])
        for k, v in per_arm_totals.items()
        if k.endswith(":A0") or k.endswith(":A1_A2_A3")
    )

    json_path.write_text(
        json.dumps(
            {
                "genome_id": genome_id,
                "n_cases_timed": len(bundles),
                "n_test_rows": len(rows),
                "total_lps": total_lps,
                "total_cpu_seconds": total_cpu,
                "total_wall_seconds": total_wall,
                "per_arm_totals": per_arm_totals,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    return tsv_path, json_path


def run_genome_cli(args: argparse.Namespace) -> int:
    """Execute truth cases for a single genome on one Shaheen node."""
    price_cases = [c for c in load_price148_cases(args.manifest) if c.genome_id == args.genome_id]
    extra_cases: List[TruthCase] = []
    if getattr(args, "extra_cases_tsv", None) is not None and args.extra_cases_tsv.exists():
        extra_cases = list(load_extra_cases(args.extra_cases_tsv, genome_id=args.genome_id))

    cases_only = (getattr(args, "cases_only", None) or "all").lower()
    if cases_only in ("extra", "swissprot"):
        genome_cases = extra_cases
        if not genome_cases:
            raise ValueError(f"No extra truth cases found for genome_id='{args.genome_id}'")
    elif cases_only in ("price", "price148"):
        genome_cases = price_cases
        if not genome_cases:
            raise ValueError(f"No Price-148 cases found for genome_id='{args.genome_id}'")
    else:
        genome_cases = price_cases + extra_cases
        if not genome_cases:
            raise ValueError(f"No truth cases found for genome_id='{args.genome_id}'")

    effective_max_cases = args.max_cases
    if effective_max_cases is None and args.run_mode == "timing":
        effective_max_cases = 3
    if effective_max_cases is not None and effective_max_cases > 0:
        selected_cases = genome_cases[:effective_max_cases]
    else:
        selected_cases = genome_cases

    metabolites = load_seed_metabolites(args.seed_metabolites)
    reactions = load_seed_reactions(args.seed_reactions, args.ec_aliases)
    ec_to_rxns = build_ec_to_reactions_map(reactions, usable_only=False)
    panel_entries = load_panel_tsv(args.panel_tsv)
    kegg_modules = load_kegg_modules(args.kegg_modules)

    draft_xml = args.genome_model_dir / f"{args.genome_id}-draft.xml"
    filled_xml = args.genome_model_dir / f"{args.genome_id}-filled.xml"
    medium_csv = args.genome_model_dir / f"{args.genome_id}-medium.csv"
    trans_tbl = args.genome_model_dir / f"{args.genome_id}-Transporter.tbl"

    draft_model = cobra.io.read_sbml_model(str(draft_xml))
    filled_model = cobra.io.read_sbml_model(str(filled_xml))
    default_medium = load_gapsmith_medium_csv(medium_csv)
    transporters = load_gapsmith_transporter_tbl(trans_tbl)

    if args.candidates_tsv is not None:
        stage1_inputs: Sequence[Path] = [args.candidates_tsv]
    elif args.stage1_tsvs:
        stage1_inputs = args.stage1_tsvs
    elif args.collection_dir is not None:
        stage1_inputs, _ = resolve_collection_stage1_tsvs(args.collection_dir, args.genome_id)
    else:
        raise ValueError("Provide --candidates-tsv, --stage1-tsvs, or --collection-dir.")

    # Supervisor fix 2026-10-30: never pin truth proteins' claims; the candidate table must not depend on labels.
    keep_prots: Tuple[str, ...] = ()
    candidates = load_stage1_candidates_from_tsv(
        stage1_inputs,
        keep_proteins=keep_prots,
    )

    known_ecs: Dict[str, List[str]] = {}
    if getattr(args, "swissprot_labels", None) is not None:
        for p, ecs in load_swissprot_experimental_ecs(args.swissprot_labels, args.genome_id).items():
            known_ecs.setdefault(p, []).extend(ecs)
    if getattr(args, "extra_cases_tsv", None) is not None and args.extra_cases_tsv.exists():
        for p, ecs in load_swissprot_experimental_ecs(args.extra_cases_tsv, args.genome_id).items():
            known_ecs.setdefault(p, []).extend(ecs)
    for c in price_cases:
        known_ecs.setdefault(c.protein_id, []).extend(c.truth_ecs)
    for c in extra_cases:
        known_ecs.setdefault(c.protein_id, []).extend(c.truth_ecs)
    known_ecs = {p: sorted(set(ecs)) for p, ecs in known_ecs.items()}

    work_dir = args.out_dir / args.genome_id
    work_dir.mkdir(parents=True, exist_ok=True)

    receipt_inputs = {
        "manifest": args.manifest,
        "seed_metabolites": args.seed_metabolites,
        "seed_reactions": args.seed_reactions,
        "panel_tsv": args.panel_tsv,
        "draft_xml": draft_xml,
        "filled_xml": filled_xml,
        "medium_csv": medium_csv,
    }
    if getattr(args, "extra_cases_tsv", None) is not None and args.extra_cases_tsv.exists():
        receipt_inputs["extra_cases_tsv"] = args.extra_cases_tsv
    for idx_p, sp in enumerate(stage1_inputs):
        receipt_inputs[f"candidates_{idx_p}_{sp.name}"] = sp

    receipt = build_environment_receipt(receipt_inputs)
    receipt["run_mode"] = args.run_mode
    receipt["cases_only"] = cases_only
    receipt["workers"] = args.workers
    receipt["max_cases"] = effective_max_cases
    receipt["max_tests_per_arm"] = args.max_tests_per_arm
    (work_dir / "receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    unmasked_arm_results = run_all_arms_unmasked(
        unmasked_draft_model=draft_model,
        universal_reactions=reactions,
        metabolites=metabolites,
        unmasked_candidates=candidates,
        ec_to_reactions=ec_to_rxns,
        panel_entries=panel_entries,
        default_medium=default_medium,
        transporters=transporters,
        kegg_modules=kegg_modules,
        replay_fraction=args.replay_fraction,
        workers=args.workers,
        max_tests_per_arm=args.max_tests_per_arm,
    )

    raw_bundles: List[CaseExecutionBundle] = []
    bundles: List[Dict[str, object]] = []
    for case in selected_cases:
        b = execute_single_case(
            case=case,
            draft_model=draft_model,
            unmasked_filled_model=filled_model,
            unmasked_candidates=candidates,
            universal_reactions=reactions,
            metabolites=metabolites,
            ec_to_reactions=ec_to_rxns,
            panel_entries=panel_entries,
            default_medium=default_medium,
            known_experimental_ecs=known_ecs,
            transporters=transporters,
            kegg_modules=kegg_modules,
            replay_fraction=args.replay_fraction,
            precomputed_unmasked_arm_results=unmasked_arm_results,
            workers=args.workers,
            max_tests_per_arm=args.max_tests_per_arm,
        )
        raw_bundles.append(b)
        bundles.append(
            {
                "case": asdict(b.case),
                "stratum_info": asdict(b.stratum_info),
                "audit": asdict(b.audit),
                "arm_results": {k: asdict(v) for k, v in b.arm_results.items()},
                "evaluations": {k: asdict(v) for k, v in b.evaluations.items()},
                "masked_candidate_evaluations": {
                    k: asdict(v) for k, v in b.masked_candidate_evaluations.items()
                },
                "unmasked_evaluations": {
                    k: asdict(v) for k, v in b.unmasked_evaluations.items()
                },
            }
        )

    timing_tsv, timing_json = write_timing_artifacts(
        work_dir=work_dir,
        genome_id=args.genome_id,
        bundles=raw_bundles,
        unmasked_arm_results=unmasked_arm_results,
    )
    if args.run_mode == "timing":
        shutil.copy2(timing_tsv, args.out_dir / f"{args.genome_id}.timing_per_test.tsv")
        shutil.copy2(timing_json, args.out_dir / f"{args.genome_id}.timing_summary.json")

    (work_dir / "genome_results.json").write_text(
        json.dumps(
            {
                "genome_id": args.genome_id,
                "run_mode": args.run_mode,
                "n_cases": len(bundles),
                "model_genes": sorted(g.id for g in draft_model.genes),
                "known_experimental_ecs": known_ecs,
                "unmasked_arm_results": {
                    k: asdict(v) for k, v in unmasked_arm_results.items()
                },
                "cases": bundles,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    if args.pack_tar:
        tar_name = getattr(args, "tar_name", None) or (
            f"{args.genome_id}.swissprot.outputs.tar"
            if cases_only in ("extra", "swissprot")
            else f"{args.genome_id}.outputs.tar"
        )
        tar_dest = args.out_dir / tar_name
        pack_directory_to_tar(work_dir, tar_dest)
    return 0


def _reconstruct_evaluation(d: Mapping[str, object]) -> CaseArmEvaluation:
    calls = tuple(
        ReactionGeneCall(
            reaction_id=str(c["reaction_id"]),
            protein_id=str(c["protein_id"]) if c.get("protein_id") is not None else None,
            ec=str(c["ec"]) if c.get("ec") is not None else None,
            score=float(c.get("score", 0.0)),
            call_kind=str(c["call_kind"]),
        )
        for c in d.get("calls", ())  # type: ignore[union-attr]
    )
    case_id = str(d["case_id"])
    case_source = str(
        d.get("case_source")
        or ("swissprot" if case_id.startswith("SP:") else "price148")
    )
    return CaseArmEvaluation(
        case_id=case_id,
        genome_id=str(d["genome_id"]),
        protein_id=str(d["protein_id"]),
        truth_ecs=tuple(str(x) for x in d.get("truth_ecs", ())),  # type: ignore[union-attr]
        target_reactions=tuple(str(x) for x in d.get("target_reactions", ())),  # type: ignore[union-attr]
        stratum=str(d["stratum"]),
        arm=str(d["arm"]),
        recovered=bool(d["recovered"]),
        reaction_recovered=bool(d["reaction_recovered"]),
        n_added_reactions=int(d["n_added_reactions"]),  # type: ignore[arg-type]
        n_gene_calls=int(d["n_gene_calls"]),  # type: ignore[arg-type]
        n_reaction_only_calls=int(d["n_reaction_only_calls"]),  # type: ignore[arg-type]
        n_supported_calls=int(d["n_supported_calls"]),  # type: ignore[arg-type]
        n_harmful_calls=int(d["n_harmful_calls"]),  # type: ignore[arg-type]
        n_unconfirmed_calls=int(d["n_unconfirmed_calls"]),  # type: ignore[arg-type]
        tests_considered=int(d["tests_considered"]),  # type: ignore[arg-type]
        tests_executed=int(d["tests_executed"]),  # type: ignore[arg-type]
        n_lps=int(d["n_lps"]),  # type: ignore[arg-type]
        n_gapfills=int(d["n_gapfills"]),  # type: ignore[arg-type]
        cpu_seconds=float(d["cpu_seconds"]),  # type: ignore[arg-type]
        wall_seconds=float(d["wall_seconds"]),  # type: ignore[arg-type]
        peak_rss_kb=int(d["peak_rss_kb"]),  # type: ignore[arg-type]
        calls=calls,
        case_source=case_source,
    )


def reclassify_calls_with_labels(
    ev: CaseArmEvaluation,
    known_ecs: Mapping[str, Sequence[str]],
    ec_to_reactions: Mapping[str, Sequence[str]],
) -> CaseArmEvaluation:
    """Re-derive supported/harmful/unconfirmed for every gene call from experimental labels (design §6).

    A call on a protein with experimental ECs is supported if its reaction maps from one of them, harmful otherwise;
    a call on a protein without experimental ECs is unconfirmed (never a negative). Reaction-only calls are unchanged.
    """
    import dataclasses

    allowed: Dict[str, Set[str]] = {}
    for prot, ecs in known_ecs.items():
        rx: Set[str] = set()
        for e in ecs:
            for r in ec_to_reactions.get(normalize_ec(e), ()):
                rx.add(base_seed_rxn_id(r))
        for alias in protein_aliases(prot):
            allowed.setdefault(alias, set()).update(rx)
    new_calls = []
    for c in ev.calls:
        if c.protein_id is None:
            new_calls.append(c)
            continue
        rx = allowed.get(c.protein_id) or allowed.get(strip_version(c.protein_id))
        if rx is None:
            kind = "unconfirmed"
        else:
            kind = "supported" if base_seed_rxn_id(c.reaction_id) in rx else "harmful"
        new_calls.append(dataclasses.replace(c, call_kind=kind))
    return dataclasses.replace(
        ev,
        calls=tuple(new_calls),
        n_supported_calls=sum(1 for c in new_calls if c.call_kind == "supported"),
        n_harmful_calls=sum(1 for c in new_calls if c.call_kind == "harmful"),
        n_unconfirmed_calls=sum(1 for c in new_calls if c.call_kind == "unconfirmed"),
    )


def _reconstruct_arm_run_result(d: Mapping[str, object]) -> ArmRunResult:
    replays = tuple(
        SolverReplayRecord(
            test_id=str(r["test_id"]),
            primary_solver=str(r["primary_solver"]),
            replay_solver=str(r["replay_solver"]),
            primary_objective=float(r["primary_objective"]),
            replay_objective=float(r["replay_objective"]),
            relative_obj_diff=float(r["relative_obj_diff"]),
            primary_added_rxns=tuple(str(x) for x in r.get("primary_added_rxns", ())),
            replay_added_rxns=tuple(str(x) for x in r.get("replay_added_rxns", ())),
            agreed=bool(r["agreed"]),
            scipy_version=str(r.get("scipy_version", getattr(scipy, "__version__", ""))),
        )
        for r in d.get("replays", ())  # type: ignore[union-attr]
    )
    return ArmRunResult(
        arm=str(d["arm"]),
        tests_considered=int(d["tests_considered"]),  # type: ignore[arg-type]
        tests_executed=int(d["tests_executed"]),  # type: ignore[arg-type]
        n_lps=int(d["n_lps"]),  # type: ignore[arg-type]
        n_gapfills=int(d["n_gapfills"]),  # type: ignore[arg-type]
        cpu_seconds=float(d["cpu_seconds"]),  # type: ignore[arg-type]
        wall_seconds=float(d["wall_seconds"]),  # type: ignore[arg-type]
        peak_rss_kb=int(d["peak_rss_kb"]),  # type: ignore[arg-type]
        added_reactions=tuple(str(x) for x in d.get("added_reactions", ())),  # type: ignore[union-attr]
        test_results=(),
        ranked_entries=(),
        replays=replays,
    )


def run_aggregate_cli(args: argparse.Namespace) -> int:
    """Merge per-genome `.outputs.tar` archives, run cohort statistics/controls, and write report."""
    solver_info = detect_solvers()
    reactions = load_seed_reactions(args.seed_reactions, args.ec_aliases)
    ec_to_rxns = build_ec_to_reactions_map(reactions, usable_only=False)

    unmapped_entries: List[UnmappedPanelEntry] = []
    if args.unmapped_tsv is not None and args.unmapped_tsv.exists():
        with args.unmapped_tsv.open("r", encoding="utf-8", newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                unmapped_entries.append(
                    UnmappedPanelEntry(
                        plate=row["plate"],
                        well=row["well"],
                        element=row["element"],
                        compound_name=row["compound_name"],
                        seed_id=row.get("seed_id") or "",
                        reason=row["reason"],
                        detail=row.get("detail") or "",
                    )
                )

    cases_only = (getattr(args, "cases_only", None) or "all").lower()
    filter_source: Optional[str] = None
    if cases_only in ("price", "price148"):
        filter_source = "price148"
    elif cases_only in ("extra", "swissprot"):
        filter_source = "swissprot"

    tar_pattern = getattr(args, "tar_pattern", "*.outputs.tar") or "*.outputs.tar"
    genome_payloads: List[Dict[str, object]] = []
    for tar_path in sorted(args.results_dir.glob(tar_pattern)):
        with tarfile.open(tar_path, "r") as tf:
            for member in tf.getmembers():
                if member.name.endswith("genome_results.json"):
                    extracted = tf.extractfile(member)
                    if extracted is not None:
                        genome_payloads.append(json.loads(extracted.read().decode("utf-8")))
    for json_path in sorted(args.results_dir.glob("*/genome_results.json")):
        payload = json.loads(json_path.read_text(encoding="utf-8"))
        if not any(
            p.get("genome_id") == payload.get("genome_id")
            and len(p.get("cases", ())) == len(payload.get("cases", ()))
            for p in genome_payloads
        ):
            genome_payloads.append(payload)

    if not genome_payloads:
        raise FileNotFoundError(
            f"No matching .outputs.tar or */genome_results.json found in '{args.results_dir}'"
        )

    evaluations_by_arm: Dict[str, List[CaseArmEvaluation]] = {}
    masked_cand_by_arm: Dict[str, List[CaseArmEvaluation]] = {}
    unmasked_by_arm: Dict[str, List[CaseArmEvaluation]] = {}
    arm_results_by_case: Dict[str, Dict[str, ArmRunResult]] = {}
    strata_by_case: Dict[str, str] = {}
    cases_list: List[TruthCase] = []
    all_replays: List[SolverReplayRecord] = []
    extra_proteins: Dict[str, List[str]] = {}
    known_ecs_by_genome: Dict[str, Dict[str, Sequence[str]]] = {}
    candidates_by_genome: Dict[str, Tuple[CandidateClaim, ...]] = {}

    import dataclasses
    seen_case_ids: Set[str] = set()

    for gpay in sorted(genome_payloads, key=lambda x: str(x.get("genome_id", ""))):
        gid = str(gpay["genome_id"])
        extra_proteins.setdefault(gid, []).extend([str(g) for g in gpay.get("model_genes", ())])
        gid_known = known_ecs_by_genome.setdefault(gid, {})
        for k, v in (gpay.get("known_experimental_ecs") or {}).items():
            existing = set(gid_known.get(str(k), ()))
            existing.update(str(e) for e in v)
            gid_known[str(k)] = tuple(sorted(existing))
        if getattr(args, "swissprot_labels", None) is not None:
            merged_known: Dict[str, List[str]] = {k: list(v) for k, v in known_ecs_by_genome[gid].items()}
            for prot, ecs in load_swissprot_experimental_ecs(args.swissprot_labels, gid).items():
                merged_known.setdefault(prot, []).extend(ecs)
            known_ecs_by_genome[gid] = {k: tuple(sorted(set(v))) for k, v in merged_known.items()}
        _relabel = (
            (lambda ev, _k=known_ecs_by_genome[gid]: reclassify_calls_with_labels(ev, _k, ec_to_rxns))
            if (getattr(args, "swissprot_labels", None) is not None or filter_source != "price148")
            else (lambda ev: ev)
        )
        if args.candidates_dir is not None:
            cand_tsv = args.candidates_dir / f"{gid}.candidates.tsv"
            if cand_tsv.exists():
                candidates_by_genome[gid] = load_stage1_candidates_from_tsv([cand_tsv])

        for cb in gpay.get("cases", ()):  # type: ignore[union-attr]
            cd = cb["case"]
            cid = str(cd["case_id"])
            c_src = str(cd.get("case_source") or ("swissprot" if cid.startswith("SP:") else "price148"))
            if filter_source == "price148" and c_src != "price148":
                continue
            if filter_source == "swissprot" and c_src != "swissprot":
                continue
            if cid in seen_case_ids:
                continue
            seen_case_ids.add(cid)

            tc = TruthCase(
                case_id=cid,
                genome_id=str(cd["genome_id"]),
                protein_id=str(cd["protein_id"]),
                truth_ecs=tuple(str(x) for x in cd["truth_ecs"]),
                case_source=c_src,
            )
            cases_list.append(tc)
            strata_by_case[tc.case_id] = str(cb["stratum_info"]["stratum"])
            arm_map: Dict[str, ArmRunResult] = {}
            for arm_name, ar_dict in cb.get("arm_results", {}).items():
                ar_obj = _reconstruct_arm_run_result(ar_dict)
                arm_map[arm_name] = ar_obj
                if arm_name in ("A0", "A1", "A2", "A3"):
                    all_replays.extend(ar_obj.replays)
            arm_results_by_case[tc.case_id] = arm_map

            for arm_name, ev_dict in cb.get("evaluations", {}).items():
                ev = _relabel(_reconstruct_evaluation(ev_dict))
                if getattr(ev, "case_source", None) != c_src:
                    ev = dataclasses.replace(ev, case_source=c_src)
                evaluations_by_arm.setdefault(arm_name, []).append(ev)
            for arm_name, mev_dict in cb.get("masked_candidate_evaluations", {}).items():
                mev = _relabel(_reconstruct_evaluation(mev_dict))
                if getattr(mev, "case_source", None) != c_src:
                    mev = dataclasses.replace(mev, case_source=c_src)
                masked_cand_by_arm.setdefault(arm_name, []).append(mev)
            for arm_name, uev_dict in cb.get("unmasked_evaluations", {}).items():
                uev = _relabel(_reconstruct_evaluation(uev_dict))
                if getattr(uev, "case_source", None) != c_src:
                    uev = dataclasses.replace(uev, case_source=c_src)
                unmasked_by_arm.setdefault(arm_name, []).append(uev)

    arm_summaries = summarize_arms_with_bootstrap(
        evaluations_by_arm, n_resamples=args.n_resamples, seed=BOOTSTRAP_SEED
    )
    masked_candidate_summaries = (
        summarize_arms_with_bootstrap(masked_cand_by_arm, n_resamples=0, seed=BOOTSTRAP_SEED)
        if masked_cand_by_arm
        else None
    )
    unmasked_summaries = (
        summarize_arms_with_bootstrap(unmasked_by_arm, n_resamples=0, seed=BOOTSTRAP_SEED)
        if unmasked_by_arm
        else None
    )
    h1_res, h2_res, h3_res = run_preregistered_tests(
        evaluations_by_arm, n_resamples=args.n_resamples, seed=BOOTSTRAP_SEED
    )
    shuffled_ctrl = run_shuffled_candidate_control(
        cases=cases_list,
        strata_by_case=strata_by_case,
        arm_results_by_case=arm_results_by_case,
        unmasked_candidates_by_genome=candidates_by_genome,
        ec_to_reactions=ec_to_rxns,
        known_experimental_ecs_by_genome=known_ecs_by_genome,
        arm_name="A1_A2_A3",
        seeds=SHUFFLE_CONTROL_SEEDS,
        extra_proteins_by_genome=extra_proteins,
    )
    replay_summary = summarize_solver_replays(all_replays)

    write_experiment_outputs(
        out_dir=args.out_dir,
        arm_summaries=arm_summaries,
        evaluations_by_arm=evaluations_by_arm,
        h1_result=h1_res,
        h2_results=h2_res,
        h3_result=h3_res,
        shuffled_control=shuffled_ctrl,
        replay_summary=replay_summary,
        solver_info=solver_info,
        unmapped_panel_entries=unmapped_entries,
        unmasked_summary=unmasked_summaries,
        masked_candidate_summary=masked_candidate_summaries,
    )
    return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    default_workers = max(1, int(os.environ.get("SLURM_CPUS_ON_NODE", "1")))
    parser = argparse.ArgumentParser(
        description="Pre-pass, per-genome runner, and cohort aggregator for off-growth enzyme discovery."
    )
    sub = parser.add_subparsers(dest="mode", required=True)

    p_pre = sub.add_parser(
        "prepass",
        help="Stream 6 predictors' stage1 TSVs into compact per-genome candidate TSVs.",
    )
    p_pre.add_argument("--collection-dir", type=Path, default=None)
    p_pre.add_argument("--stage1-tsvs", type=Path, nargs="*", default=None)
    p_pre.add_argument("--genome-id", default=None)
    p_pre.add_argument("--manifest", type=Path, default=None)
    p_pre.add_argument("--max-per-ec", type=int, default=64)
    p_pre.add_argument("--workers", type=int, default=default_workers)
    p_pre.add_argument("--out-dir", type=Path, required=True)

    p_gen = sub.add_parser("genome", help="Run truth cases for one genome.")
    p_gen.add_argument("--genome-id", required=True)
    p_gen.add_argument("--manifest", type=Path, required=True)
    p_gen.add_argument("--genome-model-dir", type=Path, required=True)
    p_gen.add_argument("--candidates-tsv", type=Path, default=None)
    p_gen.add_argument("--stage1-tsvs", type=Path, nargs="*", default=None)
    p_gen.add_argument("--collection-dir", type=Path, default=None)
    p_gen.add_argument("--seed-metabolites", type=Path, required=True)
    p_gen.add_argument("--seed-reactions", type=Path, required=True)
    p_gen.add_argument("--ec-aliases", type=Path, default=None)
    p_gen.add_argument("--panel-tsv", type=Path, required=True)
    p_gen.add_argument("--kegg-modules", type=Path, default=None)
    p_gen.add_argument("--swissprot-labels", type=Path, default=None)
    p_gen.add_argument("--replay-fraction", type=float, default=0.05)
    p_gen.add_argument(
        "--run-mode",
        choices=("pilot", "timing"),
        default="pilot",
        help="'pilot' runs all cases of the genome; 'timing' runs a fixed subset (default 3 cases) and writes timing tables.",
    )
    p_gen.add_argument("--max-cases", type=int, default=None)
    p_gen.add_argument("--max-tests-per-arm", type=int, default=None)
    p_gen.add_argument("--extra-cases-tsv", type=Path, default=None)
    p_gen.add_argument(
        "--cases-only",
        choices=("all", "price", "price148", "extra", "swissprot"),
        default="all",
        help="Run 'all' truth cases (default), 'price' (Price-148 only), or 'extra'/'swissprot' (extra cases only).",
    )
    p_gen.add_argument("--tar-name", type=str, default=None, help="Custom output tar name when --pack-tar is used.")
    p_gen.add_argument("--workers", type=int, default=default_workers)
    p_gen.add_argument("--out-dir", type=Path, required=True)
    p_gen.add_argument("--pack-tar", action="store_true", default=False)

    p_agg = sub.add_parser(
        "aggregate",
        help="Aggregate per-genome outputs across the 22 Price genomes and write report.md + summary.json.",
    )
    p_agg.add_argument("--results-dir", type=Path, required=True)
    p_agg.add_argument(
        "--cases-only",
        "--case-source",
        dest="cases_only",
        choices=("all", "combined", "price", "price148", "extra", "swissprot"),
        default="all",
        help="Filter aggregated cases by source: 'all'/'combined' (default), 'price148', or 'swissprot'.",
    )
    p_agg.add_argument(
        "--tar-pattern",
        default="*.outputs.tar",
        help="Glob pattern for tar archives in results-dir (default: '*.outputs.tar').",
    )
    p_agg.add_argument("--candidates-dir", type=Path, default=None)
    p_agg.add_argument("--seed-reactions", type=Path, required=True)
    p_agg.add_argument("--ec-aliases", type=Path, default=None)
    p_agg.add_argument("--unmapped-tsv", type=Path, default=None)
    p_agg.add_argument(
        "--swissprot-labels",
        type=Path,
        default=None,
        help="Swiss-Prot experimental EC table; calls are reclassified with Price-148 plus these labels (design §6).",
    )
    p_agg.add_argument("--n-resamples", type=int, default=BOOTSTRAP_RESAMPLES)
    p_agg.add_argument("--out-dir", type=Path, required=True)

    args = parser.parse_args(argv)
    if args.mode == "prepass":
        return run_prepass_cli(args)
    if args.mode == "genome":
        return run_genome_cli(args)
    if args.mode == "aggregate":
        return run_aggregate_cli(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

