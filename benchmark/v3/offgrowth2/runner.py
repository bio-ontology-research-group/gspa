#!/usr/bin/env python3
"""Per-genome CLI and restartable runner for off-growth round 2 (plan §4, §10)."""

from __future__ import annotations

import argparse
import collections
import csv
from dataclasses import asdict
import datetime
import gzip
import json
import os
from pathlib import Path
import platform
import re
import sys
import time
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple

import cobra
import numpy as np
import optlang
import scipy

V3_DIR = Path(__file__).resolve().parents[1]
if str(V3_DIR) not in sys.path:
    sys.path.insert(0, str(V3_DIR))

from offgrowth.arms import (  # noqa: E402
    DEFAULT_FLUX_TOL,
    DEFAULT_PANEL_UPTAKE,
    detect_solvers,
    load_gapsmith_medium_csv,
    load_kegg_modules,
)
from offgrowth.mask import (  # noqa: E402
    load_gapsmith_reactions_tbl,
    load_gapsmith_transporter_tbl,
    load_stage1_candidates_from_tsv,
)
from offgrowth.panel import (  # noqa: E402
    build_ec_to_reactions_map,
    load_panel_tsv,
    load_seed_metabolites,
    load_seed_reactions,
    sha256_file,
)
from offgrowth.runner import pack_directory_to_tar  # noqa: E402

from offgrowth2.engine import (  # noqa: E402
    CaseRecord,
    collect_known_labels_for_genome,
    draft_evidenced_base_reactions,
    execute_masked_case_superset,
    load_cases_for_genome,
    run_superset_for_model,
)
from offgrowth2.medium import minimal_medium  # noqa: E402

SCHEMA_VERSION: str = "gspa.offgrowth2.genome/1"


def _safe_case_filename(case_id: str) -> str:
    """Convert `case_id` into a safe single-component filename."""
    cleaned = re.sub(r"[^A-Za-z0-9_.:-]+", "_", case_id.strip())
    return cleaned or "case"


def write_gzip_json(path: Path, payload: Mapping[str, object]) -> None:
    """Atomically write `payload` to a `.json.gz` file with sorted keys."""
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = path.with_suffix(path.suffix + ".tmp")
    raw_bytes = (json.dumps(payload, indent=2, sort_keys=True) + "\n").encode("utf-8")
    with gzip.open(tmp_path, "wb", compresslevel=6) as fh:
        fh.write(raw_bytes)
    tmp_path.replace(path)


def read_gzip_json(path: Path) -> Dict[str, object]:
    """Read and parse a `.json.gz` file."""
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        return json.load(fh)


def _compute_code_hashes() -> Dict[str, str]:
    """Compute SHA-256 hashes of every file in `offgrowth2/` and imported pilot/followup files."""
    offgrowth2_dir = Path(__file__).resolve().parent
    offgrowth_dir = V3_DIR / "offgrowth"
    followups_dir = V3_DIR / "followups"

    hashes: Dict[str, str] = {}
    for p in sorted(offgrowth2_dir.iterdir()):
        if p.is_file() and p.suffix in {".py", ".md"}:
            hashes[f"offgrowth2/{p.name}"] = sha256_file(p)
    for name in ("arms.py", "mask.py", "panel.py", "score.py", "runner.py"):
        p = offgrowth_dir / name
        if p.is_file():
            hashes[f"offgrowth/{name}"] = sha256_file(p)
    common_p = followups_dir / "common.py"
    if common_p.is_file():
        hashes["followups/common.py"] = sha256_file(common_p)
    return hashes


def _filter_cases_by_ids(
    cases: Sequence[CaseRecord],
    case_ids_csv: Optional[str],
    label_prefix: str = "FB:",
) -> Tuple[CaseRecord, ...]:
    if not case_ids_csv:
        return tuple(cases)
    requested: Set[str] = {tok.strip() for tok in case_ids_csv.split(",") if tok.strip()}
    if not requested:
        return tuple(cases)
    kept: List[CaseRecord] = []
    for c in cases:
        prefixed_pid = (
            c.protein_id
            if c.protein_id.startswith(label_prefix)
            else f"{label_prefix}{c.protein_id}"
        )
        if (
            c.case_id in requested
            or c.protein_id in requested
            or prefixed_pid in requested
            or f"SP:{c.protein_id}" in requested
        ):
            kept.append(c)
    return tuple(kept)


def _write_timing_and_fill_tables(
    work_dir: Path,
    genome_id: str,
    runs: Sequence[Mapping[str, object]],
) -> Tuple[Path, Path]:
    """Write `timing_per_test.tsv` and `fill_sizes.tsv` from the assembled run dicts."""
    timing_path = work_dir / "timing_per_test.tsv"
    fill_sizes_path = work_dir / "fill_sizes.tsv"

    timing_fields = [
        "genome_id",
        "run_kind",
        "case_id",
        "test_id",
        "test_type",
        "target_id",
        "panel_id",
        "element",
        "gapfill_ran",
        "feasible_in_universal",
        "already_satisfied",
        "n_lps",
        "cpu_seconds",
        "fill_size",
        "stable_fill_size",
    ]

    fill_counts: Dict[Tuple[str, int], int] = collections.Counter()

    with timing_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=timing_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for run in runs:
            run_kind = str(run.get("run_kind", ""))
            case_obj = run.get("case")
            case_id = (
                str(case_obj["case_id"])
                if isinstance(case_obj, Mapping) and case_obj.get("case_id") is not None
                else "UNMASKED"
            )
            tests = run.get("tests") or ()
            if isinstance(tests, Sequence):
                for t in tests:
                    if not isinstance(t, Mapping):
                        continue
                    t_type = str(t.get("test_type", ""))
                    fill_list = t.get("fill") or ()
                    stable_list = t.get("stable_fill") or ()
                    f_size = len(fill_list) if isinstance(fill_list, Sequence) else 0
                    s_size = len(stable_list) if isinstance(stable_list, Sequence) else 0
                    fill_counts[(t_type, f_size)] += 1
                    writer.writerow(
                        {
                            "genome_id": genome_id,
                            "run_kind": run_kind,
                            "case_id": case_id,
                            "test_id": str(t.get("test_id", "")),
                            "test_type": t_type,
                            "target_id": str(t.get("target_id", "")),
                            "panel_id": "" if t.get("panel_id") is None else str(t.get("panel_id")),
                            "element": "" if t.get("element") is None else str(t.get("element")),
                            "gapfill_ran": int(bool(t.get("gapfill_ran", False))),
                            "feasible_in_universal": int(bool(t.get("feasible_in_universal", False))),
                            "already_satisfied": int(bool(t.get("already_satisfied", False))),
                            "n_lps": int(t.get("n_lps", 0)),  # type: ignore[arg-type]
                            "cpu_seconds": f"{float(t.get('cpu_seconds', 0.0)):.6f}",  # type: ignore[arg-type]
                            "fill_size": f_size,
                            "stable_fill_size": s_size,
                        }
                    )

    with fill_sizes_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["test_type", "fill_size", "count"])
        for (t_type, f_size), count in sorted(fill_counts.items(), key=lambda kv: (kv[0][0], kv[0][1])):
            writer.writerow([t_type, f_size, count])

    return timing_path, fill_sizes_path


def _load_resumed_state(
    work_dir: Path,
    cases_dir: Path,
) -> Tuple[Optional[Dict[str, object]], Dict[str, Dict[str, object]], Dict[str, Dict[str, object]]]:
    """Load finished unmasked run, masked case runs, and skipped cases for `--resume`."""
    resumed_unmasked: Optional[Dict[str, object]] = None
    resumed_cases: Dict[str, Dict[str, object]] = {}
    resumed_skipped: Dict[str, Dict[str, object]] = {}

    results_gz = work_dir / "genome_results2.json.gz"
    if results_gz.exists():
        try:
            existing = read_gzip_json(results_gz)
            for r in existing.get("runs", ()):  # type: ignore[union-attr]
                if not isinstance(r, Mapping):
                    continue
                if r.get("run_kind") == "unmasked":
                    resumed_unmasked = dict(r)
                elif r.get("run_kind") == "masked" and isinstance(r.get("case"), Mapping):
                    cid = str(r["case"]["case_id"])
                    resumed_cases[cid] = dict(r)
            for sk in existing.get("skipped_cases", ()):  # type: ignore[union-attr]
                if isinstance(sk, Mapping) and isinstance(sk.get("case"), Mapping):
                    cid = str(sk["case"]["case_id"])
                    resumed_skipped[cid] = dict(sk)
        except Exception:
            pass

    unmasked_part = cases_dir / "UNMASKED.json.gz"
    if resumed_unmasked is None and unmasked_part.exists():
        try:
            loaded = read_gzip_json(unmasked_part)
            if loaded.get("run_kind") == "unmasked":
                resumed_unmasked = loaded
        except Exception:
            pass

    if cases_dir.exists():
        for part_file in sorted(cases_dir.glob("*.json.gz")):
            if part_file.name == "UNMASKED.json.gz":
                continue
            try:
                loaded = read_gzip_json(part_file)
                if loaded.get("skipped") and isinstance(loaded.get("case"), Mapping):
                    cid = str(loaded["case"]["case_id"])
                    resumed_skipped.setdefault(cid, loaded)
                elif loaded.get("run_kind") == "masked" and isinstance(loaded.get("case"), Mapping):
                    cid = str(loaded["case"]["case_id"])
                    resumed_cases.setdefault(cid, loaded)
            except Exception:
                continue

    return resumed_unmasked, resumed_cases, resumed_skipped


def run_genome(
    genome_id: str,
    genome_model_dir: Path,
    candidates_tsv: Path,
    candidate_table_kind: str,
    seed_reactions_path: Path,
    seed_metabolites_path: Path,
    panel_tsv_path: Path,
    out_dir: Path,
    ec_aliases_path: Optional[Path] = None,
    kegg_modules_path: Optional[Path] = None,
    price_manifest: Optional[Path] = None,
    extra_cases_tsv: Optional[Path] = None,
    label_cases_tsv: Optional[Path] = None,
    label_cases_prefix: str = "FB:",
    layout_tsv: Optional[Path] = None,
    max_cases: Optional[int] = None,
    case_ids: Optional[str] = None,
    skip_unmasked: bool = False,
    resume: bool = False,
    workers: int = 1,
    pack_tar: bool = False,
    replay_fraction: float = 0.05,
    flux_tol: float = DEFAULT_FLUX_TOL,
    uptake_rate: float = DEFAULT_PANEL_UPTAKE,
) -> Dict[str, object]:
    """Execute the round-2 superset for `genome_id` and write all required artifacts."""
    t_wall0 = time.monotonic()
    work_dir = Path(out_dir) / genome_id
    cases_dir = work_dir / "cases"
    work_dir.mkdir(parents=True, exist_ok=True)
    cases_dir.mkdir(parents=True, exist_ok=True)
    log_path = work_dir / "runner.log"

    def _log(msg: str) -> None:
        stamp = datetime.datetime.now().isoformat(timespec="seconds")
        line = f"[{stamp}] {msg}\n"
        with log_path.open("a", encoding="utf-8") as lfh:
            lfh.write(line)

    _log(
        f"Starting offgrowth2 runner for genome={genome_id} "
        f"table_kind={candidate_table_kind} workers={workers} resume={resume}"
    )

    draft_xml = Path(genome_model_dir) / f"{genome_id}-draft.xml"
    filled_xml = Path(genome_model_dir) / f"{genome_id}-filled.xml"
    medium_csv = Path(genome_model_dir) / f"{genome_id}-medium.csv"
    trans_tbl = Path(genome_model_dir) / f"{genome_id}-Transporter.tbl"
    rxns_tbl = Path(genome_model_dir) / f"{genome_id}-all-Reactions.tbl"

    draft_model = cobra.io.read_sbml_model(str(draft_xml))
    filled_model = (
        cobra.io.read_sbml_model(str(filled_xml))
        if filled_xml.exists()
        else cobra.io.read_sbml_model(str(draft_xml))
    )
    default_med = load_gapsmith_medium_csv(medium_csv)
    transporters = load_gapsmith_transporter_tbl(trans_tbl) if trans_tbl.exists() else ()
    reaction_evidence = load_gapsmith_reactions_tbl(rxns_tbl) if rxns_tbl.exists() else ()

    metabolites = load_seed_metabolites(seed_metabolites_path)
    reactions = load_seed_reactions(seed_reactions_path, ec_aliases_path)
    ec_to_rxns = build_ec_to_reactions_map(reactions, usable_only=False)
    panel_entries = load_panel_tsv(panel_tsv_path)
    _ = load_kegg_modules(kegg_modules_path)
    candidates = load_stage1_candidates_from_tsv([Path(candidates_tsv)])

    min_med = minimal_medium(default_med, metabolites)

    all_genome_cases = load_cases_for_genome(
        genome_id=genome_id,
        price_manifest=price_manifest,
        extra_cases_tsv=extra_cases_tsv,
        label_cases_tsv=label_cases_tsv,
        label_cases_prefix=label_cases_prefix,
    )
    filtered_cases = _filter_cases_by_ids(
        all_genome_cases, case_ids, label_prefix=label_cases_prefix
    )
    if max_cases is not None and max_cases > 0:
        selected_cases = filtered_cases[:max_cases]
    else:
        selected_cases = filtered_cases

    known_labels = collect_known_labels_for_genome(
        genome_id=genome_id,
        price_manifest=price_manifest,
        extra_cases_tsv=extra_cases_tsv,
        label_cases_tsv=label_cases_tsv,
    )

    _log(
        f"Loaded inputs: panel={len(panel_entries)} candidates={len(candidates)} "
        f"cases_selected={len(selected_cases)} known_labelled_proteins={len(known_labels)}"
    )

    resumed_unmasked: Optional[Dict[str, object]] = None
    resumed_cases: Dict[str, Dict[str, object]] = {}
    resumed_skipped: Dict[str, Dict[str, object]] = {}
    if resume:
        resumed_unmasked, resumed_cases, resumed_skipped = _load_resumed_state(work_dir, cases_dir)
        _log(
            f"Resume state: unmasked_done={resumed_unmasked is not None} "
            f"finished_cases={len(resumed_cases)} skipped_cases={len(resumed_skipped)}"
        )

    runs_list: List[Dict[str, object]] = []
    skipped_cases_list: List[Dict[str, object]] = []

    # 1. Unmasked deployment run (unless --skip-unmasked)
    if not skip_unmasked:
        if resume and resumed_unmasked is not None:
            _log("Skipping unmasked run (loaded from resume state)")
            runs_list.append(resumed_unmasked)
        else:
            _log("Running unmasked superset")
            unmasked_res = run_superset_for_model(
                draft_model=draft_model,
                universal_reactions=reactions,
                metabolites=metabolites,
                candidates=candidates,
                ec_to_reactions=ec_to_rxns,
                panel_entries=panel_entries,
                default_medium=default_med,
                minimal_med=min_med,
                case=None,
                audit_report=None,
                stratum_info=None,
                target_reactions=(),
                allow_unmasked=True,
                replay_fraction=replay_fraction,
                flux_tol=flux_tol,
                uptake_rate=uptake_rate,
                workers=workers,
            )
            unmasked_dict = unmasked_res.to_schema_dict()
            write_gzip_json(cases_dir / "UNMASKED.json.gz", unmasked_dict)
            runs_list.append(unmasked_dict)
            _log(
                f"Finished unmasked run: n_lps={unmasked_res.n_lps} "
                f"wall={unmasked_res.wall_seconds:.2f}s"
            )

    # 2. Masked case runs
    for idx, case_rec in enumerate(selected_cases, 1):
        cid = case_rec.case_id
        part_path = cases_dir / f"{_safe_case_filename(cid)}.json.gz"
        if resume and cid in resumed_cases:
            _log(f"[{idx}/{len(selected_cases)}] Skipping finished case {cid} (--resume)")
            runs_list.append(resumed_cases[cid])
            continue
        if resume and cid in resumed_skipped:
            _log(f"[{idx}/{len(selected_cases)}] Skipping previously failed-audit case {cid} (--resume)")
            skipped_cases_list.append(resumed_skipped[cid])
            continue

        _log(f"[{idx}/{len(selected_cases)}] Running masked case {cid} (protein={case_rec.protein_id})")
        case_run, audit = execute_masked_case_superset(
            case_record=case_rec,
            draft_model=draft_model,
            unmasked_filled_model=filled_model,
            unmasked_candidates=candidates,
            universal_reactions=reactions,
            metabolites=metabolites,
            ec_to_reactions=ec_to_rxns,
            panel_entries=panel_entries,
            default_medium=default_med,
            minimal_med=min_med,
            transporters=transporters,
            reaction_evidence=reaction_evidence,
            replay_fraction=replay_fraction,
            flux_tol=flux_tol,
            uptake_rate=uptake_rate,
            workers=workers,
        )
        if case_run is None:
            raw_audit = audit.to_dict()
            skip_entry: Dict[str, object] = {
                "skipped": True,
                "case": case_rec.to_schema_dict(),
                "reason": "masking_audit_failed",
                "audit": {
                    "case_id": str(raw_audit["case_id"]),
                    "protein_id": str(raw_audit["protein_id"]),
                    "truth_ecs": sorted(str(x) for x in raw_audit.get("truth_ecs", ())),
                    "equivalent_ecs": sorted(str(x) for x in raw_audit.get("equivalent_ecs", ())),
                    "target_reactions": sorted(str(x) for x in raw_audit.get("target_reactions", ())),
                    "passed": bool(raw_audit["passed"]),
                    "violations": sorted(str(x) for x in raw_audit.get("violations", ())),
                },
            }
            write_gzip_json(part_path, skip_entry)
            skipped_cases_list.append(skip_entry)
            _log(
                f"Masking audit FAILED for case {cid}; case skipped: "
                f"violations={audit.violations}"
            )
            continue

        run_dict = case_run.to_schema_dict()
        write_gzip_json(part_path, run_dict)
        runs_list.append(run_dict)
        _log(
            f"[{idx}/{len(selected_cases)}] Finished case {cid}: "
            f"n_lps={case_run.n_lps} wall={case_run.wall_seconds:.2f}s"
        )

    # Collect input file hashes
    candidate_inputs: List[Tuple[str, Optional[Path]]] = [
        ("draft_xml", draft_xml),
        ("filled_xml", filled_xml),
        ("medium_csv", medium_csv),
        ("transporter_tbl", trans_tbl),
        ("reactions_tbl", rxns_tbl),
        ("candidates_tsv", Path(candidates_tsv)),
        ("seed_reactions", Path(seed_reactions_path)),
        ("seed_metabolites", Path(seed_metabolites_path)),
        ("ec_aliases", Path(ec_aliases_path) if ec_aliases_path is not None else None),
        ("panel_tsv", Path(panel_tsv_path)),
        ("kegg_modules", Path(kegg_modules_path) if kegg_modules_path is not None else None),
        ("price_manifest", Path(price_manifest) if price_manifest is not None else None),
        ("extra_cases_tsv", Path(extra_cases_tsv) if extra_cases_tsv is not None else None),
        ("label_cases_tsv", Path(label_cases_tsv) if label_cases_tsv is not None else None),
        ("layout_tsv", Path(layout_tsv) if layout_tsv is not None else None),
    ]
    input_hashes_by_path: Dict[str, str] = {}
    input_hashes_by_label: Dict[str, str] = {}
    for label, p in candidate_inputs:
        if p is not None and p.exists() and p.is_file():
            digest = sha256_file(p)
            input_hashes_by_path[str(p)] = digest
            input_hashes_by_label[label] = digest

    panel_schema_list = [
        {
            "panel_id": e.panel_id,
            "element": e.element,
            "seed_id": e.seed_id,
            "exchange_rxn": e.exchange_rxn,
        }
        for e in panel_entries
    ]

    genome_payload: Dict[str, object] = {
        "schema": SCHEMA_VERSION,
        "genome_id": genome_id,
        "candidate_table_kind": candidate_table_kind,
        "inputs": input_hashes_by_path,
        "panel": panel_schema_list,
        "minimal_medium": {k: float(v) for k, v in sorted(min_med.items())},
        "default_medium": {k: float(v) for k, v in sorted(default_med.items())},
        "model_genes": sorted({g.id for g in draft_model.genes}),
        "draft_reactions": list(draft_evidenced_base_reactions(draft_model)),
        "known_labels": known_labels,
        "skipped_cases": skipped_cases_list,
        "runs": runs_list,
    }

    results_gz_path = work_dir / "genome_results2.json.gz"
    write_gzip_json(results_gz_path, genome_payload)
    _write_timing_and_fill_tables(work_dir, genome_id, runs_list)

    wall_total = time.monotonic() - t_wall0
    solver_info = detect_solvers()
    receipt: Dict[str, object] = {
        "schema": SCHEMA_VERSION,
        "genome_id": genome_id,
        "candidate_table_kind": candidate_table_kind,
        "hostname": platform.node(),
        "platform": platform.platform(),
        "python_version": sys.version,
        "package_versions": {
            "cobra": getattr(cobra, "__version__", "unknown"),
            "scipy": getattr(scipy, "__version__", "unknown"),
            "numpy": getattr(np, "__version__", "unknown"),
            "optlang": getattr(optlang, "__version__", "unknown"),
        },
        "solver_info": asdict(solver_info),
        "input_sha256": input_hashes_by_label,
        "input_paths_sha256": input_hashes_by_path,
        "code_sha256": _compute_code_hashes(),
        "n_runs": len(runs_list),
        "n_skipped_cases": len(skipped_cases_list),
        "skipped_cases": skipped_cases_list,
        "wall_seconds": wall_total,
    }
    (work_dir / "receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    _log(f"Completed genome={genome_id} in {wall_total:.2f}s -> {results_gz_path}")

    if pack_tar:
        tar_dest = Path(out_dir) / f"{genome_id}.outputs.tar"
        pack_directory_to_tar(work_dir, tar_dest)

    return genome_payload


def _add_genome_arguments(parser: argparse.ArgumentParser, default_workers: int) -> None:
    parser.add_argument("--genome-id", required=True, help="Genome identifier.")
    parser.add_argument(
        "--genome-model-dir",
        type=Path,
        required=True,
        help="Directory containing <genome>-draft.xml, <genome>-filled.xml, <genome>-medium.csv, <genome>-Transporter.tbl.",
    )
    parser.add_argument(
        "--candidates-tsv",
        type=Path,
        required=True,
        help="Path to candidate table TSV (protein_id, ec, score).",
    )
    parser.add_argument(
        "--candidate-table-kind",
        choices=("six", "clean"),
        default="six",
        help="Candidate table source kind ('six' or 'clean').",
    )
    parser.add_argument("--seed-reactions", type=Path, required=True)
    parser.add_argument("--seed-metabolites", type=Path, required=True)
    parser.add_argument("--ec-aliases", type=Path, default=None)
    parser.add_argument("--panel-tsv", type=Path, required=True)
    parser.add_argument("--kegg-modules", type=Path, default=None)
    parser.add_argument("--price-manifest", "--manifest", dest="price_manifest", type=Path, default=None)
    parser.add_argument("--extra-cases-tsv", type=Path, default=None)
    parser.add_argument("--label-cases-tsv", type=Path, default=None)
    parser.add_argument("--label-cases-prefix", type=str, default="FB:")
    parser.add_argument("--layout-tsv", type=Path, default=None)
    parser.add_argument("--max-cases", type=int, default=None)
    parser.add_argument("--case-ids", type=str, default=None)
    parser.add_argument("--skip-unmasked", action="store_true", default=False)
    parser.add_argument("--resume", action="store_true", default=False)
    parser.add_argument("--workers", type=int, default=default_workers)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--pack-tar", action="store_true", default=False)
    parser.add_argument("--replay-fraction", type=float, default=0.05)


def run_genome_cli(args: argparse.Namespace) -> int:
    run_genome(
        genome_id=args.genome_id,
        genome_model_dir=args.genome_model_dir,
        candidates_tsv=args.candidates_tsv,
        candidate_table_kind=args.candidate_table_kind,
        seed_reactions_path=args.seed_reactions,
        seed_metabolites_path=args.seed_metabolites,
        panel_tsv_path=args.panel_tsv,
        out_dir=args.out_dir,
        ec_aliases_path=args.ec_aliases,
        kegg_modules_path=args.kegg_modules,
        price_manifest=args.price_manifest,
        extra_cases_tsv=args.extra_cases_tsv,
        label_cases_tsv=args.label_cases_tsv,
        label_cases_prefix=args.label_cases_prefix,
        layout_tsv=args.layout_tsv,
        max_cases=args.max_cases,
        case_ids=args.case_ids,
        skip_unmasked=args.skip_unmasked,
        resume=args.resume,
        workers=args.workers,
        pack_tar=args.pack_tar,
        replay_fraction=args.replay_fraction,
    )
    return 0


def build_parser() -> argparse.ArgumentParser:
    default_workers = max(1, int(os.environ.get("SLURM_CPUS_ON_NODE", "1")))
    parser = argparse.ArgumentParser(
        description="Per-genome simulation engine runner for off-growth round 2."
    )
    sub = parser.add_subparsers(dest="mode")
    p_gen = sub.add_parser("genome", help="Run round-2 test superset for one genome.")
    _add_genome_arguments(p_gen, default_workers)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    raw_argv = list(sys.argv[1:] if argv is None else argv)
    if raw_argv and raw_argv[0].startswith("-") and raw_argv[0] not in ("-h", "--help"):
        raw_argv = ["genome"] + raw_argv
    parser = build_parser()
    args = parser.parse_args(raw_argv)
    if args.mode == "genome":
        return run_genome_cli(args)
    parser.print_help()
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
