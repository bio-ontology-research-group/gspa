#!/usr/bin/env python3
"""Stage 2B (flux-certified refiner) Price collection runner (`stage0` vs `stage2b`).

Reuses functions from the frozen Price engine modules (`run_price_production_collection`,
`run_price_matrix_capacity`, `run_price_local`, `verified_batch`, `verify_candidate_universe`,
`stage_ledger`, `v3_runner`) by import without modifying any frozen file.

Arms:
  * `stage0`: unchanged stage0 command (`--enable-priors none`).
  * `stage2b`: stage1 command (default priors) + `--enable-fba-refiner --fba-data-dir <dir>
    --fba-gram neg --fba-medium <medium> --fba-python-bin <py>` (and optional `--fba-sidecar`).
"""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import re
import shutil
import sys
import tempfile
import time
import traceback
from typing import Any

V3_DIR = Path(__file__).resolve().parents[1]
REPO_ROOT = V3_DIR.parents[1]
if str(V3_DIR) not in sys.path:
    sys.path.insert(0, str(V3_DIR))

import run_price_local as L  # noqa: E402
import run_price_matrix_capacity as C  # noqa: E402
import run_price_production_collection as P  # noqa: E402
import stage_ledger as SL  # noqa: E402
from v3_runner import BATCH_HEADER  # noqa: E402
import verified_batch as V  # noqa: E402
import verify_candidate_universe as U  # noqa: E402

save = P.save
check = P.check
cells = P.cells
completed_attempts = P.completed_attempts

ARMS = ("stage0", "stage2b")

FBA_REQUIRED_DATA_FILES = (
    "universal_metadata.json",
    "fba_matrices_v6RE.pkl",
    "tight_bounds_v6_neg.pkl",
    "tight_bounds_v6_pos.pkl",
    "seedr2ec.pkl",
    "seedec2r.pkl",
    "medium.pkl",
)

# All 22 Price genomes with their NCBI taxids (from layout/report.json and
# assembly_data_report.json), NCBI taxonomy lineage (phylum Pseudomonadota,
# taxid 1224; classes Alphaproteobacteria, Betaproteobacteria, Gammaproteobacteria),
# and Gram-negative stain status.
PRICE22_TAXONOMY: dict[str, dict[str, Any]] = {
    "GCA_000006965.1": {
        "taxid": 266834,
        "organism_name": "Sinorhizobium meliloti 1021",
        "phylum": "Pseudomonadota",
        "phylum_taxid": 1224,
        "class_name": "Alphaproteobacteria",
        "order": "Hyphomicrobiales",
        "family": "Rhizobiaceae",
        "gram": "neg",
    },
    "GCA_003443575.1": {
        "taxid": 221822,
        "organism_name": "Phaeobacter inhibens",
        "phylum": "Pseudomonadota",
        "phylum_taxid": 1224,
        "class_name": "Alphaproteobacteria",
        "order": "Rhodobacterales",
        "family": "Roseobacteraceae",
        "gram": "neg",
    },
    "GCA_004402975.1": {
        "taxid": 1211383,
        "organism_name": "Paraburkholderia dipogonis",
        "phylum": "Pseudomonadota",
        "phylum_taxid": 1224,
        "class_name": "Betaproteobacteria",
        "order": "Burkholderiales",
        "family": "Burkholderiaceae",
        "gram": "neg",
    },
    "GCF_000015245.1": {
        "taxid": 326297,
        "organism_name": "Shewanella amazonensis SB2B",
        "phylum": "Pseudomonadota",
        "phylum_taxid": 1224,
        "class_name": "Gammaproteobacteria",
        "order": "Alteromonadales",
        "family": "Shewanellaceae",
        "gram": "neg",
    },
    "GCF_000016065.1": {
        "taxid": 323850,
        "organism_name": "Shewanella loihica PV-4",
        "phylum": "Pseudomonadota",
        "phylum_taxid": 1224,
        "class_name": "Gammaproteobacteria",
        "order": "Alteromonadales",
        "family": "Shewanellaceae",
        "gram": "neg",
    },
    "GCF_000018145.1": {
        "taxid": 398580,
        "organism_name": "Dinoroseobacter shibae DFL 12 = DSM 16493",
        "phylum": "Pseudomonadota",
        "phylum_taxid": 1224,
        "class_name": "Alphaproteobacteria",
        "order": "Rhodobacterales",
        "family": "Roseobacteraceae",
        "gram": "neg",
    },
    "GCF_000020125.1": {
        "taxid": 398527,
        "organism_name": "Paraburkholderia phytofirmans PsJN",
        "phylum": "Pseudomonadota",
        "phylum_taxid": 1224,
        "class_name": "Betaproteobacteria",
        "order": "Burkholderiales",
        "family": "Burkholderiaceae",
        "gram": "neg",
    },
    "GCF_000146165.2": {
        "taxid": 211586,
        "organism_name": "Shewanella oneidensis MR-1",
        "phylum": "Pseudomonadota",
        "phylum_taxid": 1224,
        "class_name": "Gammaproteobacteria",
        "order": "Alteromonadales",
        "family": "Shewanellaceae",
        "gram": "neg",
    },
    "GCF_000166295.1": {
        "taxid": 225937,
        "organism_name": "Marinobacter adhaerens HP15",
        "phylum": "Pseudomonadota",
        "phylum_taxid": 1224,
        "class_name": "Gammaproteobacteria",
        "order": "Alteromonadales",
        "family": "Marinobacteraceae",
        "gram": "neg",
    },
    "GCF_000203935.1": {
        "taxid": 94122,
        "organism_name": "Shewanella sp. ANA-3",
        "phylum": "Pseudomonadota",
        "phylum_taxid": 1224,
        "class_name": "Gammaproteobacteria",
        "order": "Alteromonadales",
        "family": "Shewanellaceae",
        "gram": "neg",
    },
    "GCF_000236665.1": {
        "taxid": 640081,
        "organism_name": "Azospira oryzae PS",
        "phylum": "Pseudomonadota",
        "phylum_taxid": 1224,
        "class_name": "Betaproteobacteria",
        "order": "Rhodocyclales",
        "family": "Rhodocyclaceae",
        "gram": "neg",
    },
    "GCF_000237365.1": {
        "taxid": 1064539,
        "organism_name": "Azospirillum baldaniorum",
        "phylum": "Pseudomonadota",
        "phylum_taxid": 1224,
        "class_name": "Alphaproteobacteria",
        "order": "Rhodospirillales",
        "family": "Azospirillaceae",
        "gram": "neg",
    },
    "GCF_000327065.1": {
        "taxid": 644801,
        "organism_name": "Stutzerimonas stutzeri RCH2",
        "phylum": "Pseudomonadota",
        "phylum_taxid": 1224,
        "class_name": "Gammaproteobacteria",
        "order": "Pseudomonadales",
        "family": "Pseudomonadaceae",
        "gram": "neg",
    },
    "GCF_000698265.1": {
        "taxid": 321846,
        "organism_name": "Pseudomonas simiae",
        "phylum": "Pseudomonadota",
        "phylum_taxid": 1224,
        "class_name": "Gammaproteobacteria",
        "order": "Pseudomonadales",
        "family": "Pseudomonadaceae",
        "gram": "neg",
    },
    "GCF_000832305.1": {
        "taxid": 68895,
        "organism_name": "Cupriavidus basilensis",
        "phylum": "Pseudomonadota",
        "phylum_taxid": 1224,
        "class_name": "Betaproteobacteria",
        "order": "Burkholderiales",
        "family": "Burkholderiaceae",
        "gram": "neg",
    },
    "GCF_001307155.1": {
        "taxid": 294,
        "organism_name": "Pseudomonas fluorescens",
        "phylum": "Pseudomonadota",
        "phylum_taxid": 1224,
        "class_name": "Gammaproteobacteria",
        "order": "Pseudomonadales",
        "family": "Pseudomonadaceae",
        "gram": "neg",
    },
    "GCF_001307275.1": {
        "taxid": 294,
        "organism_name": "Pseudomonas fluorescens",
        "phylum": "Pseudomonadota",
        "phylum_taxid": 1224,
        "class_name": "Gammaproteobacteria",
        "order": "Pseudomonadales",
        "family": "Pseudomonadaceae",
        "gram": "neg",
    },
    "GCF_001440325.1": {
        "taxid": 294,
        "organism_name": "Pseudomonas fluorescens",
        "phylum": "Pseudomonadota",
        "phylum_taxid": 1224,
        "class_name": "Gammaproteobacteria",
        "order": "Pseudomonadales",
        "family": "Pseudomonadaceae",
        "gram": "neg",
    },
    "GCF_001623525.1": {
        "taxid": 294,
        "organism_name": "Pseudomonas fluorescens",
        "phylum": "Pseudomonadota",
        "phylum_taxid": 1224,
        "class_name": "Gammaproteobacteria",
        "order": "Pseudomonadales",
        "family": "Pseudomonadaceae",
        "gram": "neg",
    },
    "GCF_001625455.1": {
        "taxid": 294,
        "organism_name": "Pseudomonas fluorescens",
        "phylum": "Pseudomonadota",
        "phylum_taxid": 1224,
        "class_name": "Gammaproteobacteria",
        "order": "Pseudomonadales",
        "family": "Pseudomonadaceae",
        "gram": "neg",
    },
    "GCF_001633105.1": {
        "taxid": 1813946,
        "organism_name": "Acidovorax sp. GW101-3H11",
        "phylum": "Pseudomonadota",
        "phylum_taxid": 1224,
        "class_name": "Betaproteobacteria",
        "order": "Burkholderiales",
        "family": "Comamonadaceae",
        "gram": "neg",
    },
    "GCF_002797435.1": {
        "taxid": 93064,
        "organism_name": "Sphingomonas koreensis",
        "phylum": "Pseudomonadota",
        "phylum_taxid": 1224,
        "class_name": "Alphaproteobacteria",
        "order": "Sphingomonadales",
        "family": "Sphingomonadaceae",
        "gram": "neg",
    },
}


def verify_price22_gram_negative_pseudomonadota(
    layout_report: dict[str, Any],
    *,
    allow_synthetic_subset: bool = False,
) -> dict[str, Any]:
    """Verify from layout/report.json and taxonomy table that all Price genomes are Gram-negative Pseudomonadota."""
    assemblies = layout_report.get("assemblies", {})
    if not assemblies:
        raise ValueError("layout/report.json contains no assemblies")

    if not allow_synthetic_subset and set(assemblies.keys()) != set(PRICE22_TAXONOMY.keys()):
        raise ValueError(
            f"layout/report.json assemblies do not match the 22 Price genomes: "
            f"got {sorted(assemblies.keys())}"
        )

    verified_genomes: dict[str, dict[str, Any]] = {}
    for asm, rec in sorted(assemblies.items()):
        taxid = int(rec["taxid"])
        if asm in PRICE22_TAXONOMY:
            expected = PRICE22_TAXONOMY[asm]
            if taxid != expected["taxid"]:
                raise ValueError(
                    f"Assembly {asm} taxid mismatch in layout/report.json: got {taxid}, expected {expected['taxid']}"
                )
            if expected["phylum"] != "Pseudomonadota" or expected["gram"] != "neg":
                raise ValueError(f"Assembly {asm} is not Gram-negative Pseudomonadota: {expected}")
            verified_genomes[asm] = dict(expected)
        elif allow_synthetic_subset:
            verified_genomes[asm] = {
                "taxid": taxid,
                "organism_name": f"synthetic_{asm}",
                "phylum": "Pseudomonadota",
                "phylum_taxid": 1224,
                "class_name": "Gammaproteobacteria",
                "gram": "neg",
            }
        else:
            raise ValueError(f"Unknown Price assembly {asm!r}")

    return {
        "all_gram_negative": True,
        "all_pseudomonadota": True,
        "phylum": "Pseudomonadota",
        "phylum_taxid": 1224,
        "gram": "neg",
        "n_genomes": len(verified_genomes),
        "genomes": verified_genomes,
    }


def stage2b_extra_flags(
    fba_data_dir: Path | str,
    fba_gram: str,
    fba_medium: str,
    fba_python_bin: Path | str,
    fba_sidecar: Path | str | None = None,
) -> list[str]:
    flags = [
        "--enable-fba-refiner",
        "--fba-data-dir",
        str(fba_data_dir),
        "--fba-gram",
        str(fba_gram),
        "--fba-medium",
        str(fba_medium),
        "--fba-python-bin",
        str(fba_python_bin),
    ]
    if fba_sidecar is not None:
        flags += ["--fba-sidecar", str(fba_sidecar)]
    return flags


def verify_arm_command_pair(
    stage0_cmd: list[str],
    stage2b_cmd: list[str],
    *,
    fba_data_dir: str,
    fba_gram: str,
    fba_medium: str,
    fba_python_bin: str,
    fba_sidecar: str | None = None,
) -> None:
    """Enforce that stage0 and stage2b differ ONLY in directory paths and the documented Stage 2B flag set."""
    if len(stage0_cmd) < 2 or stage0_cmd[-2:] != ["--enable-priors", "none"]:
        raise ValueError(f"stage0 command must end with ['--enable-priors', 'none']: {stage0_cmd}")

    expected_extra = stage2b_extra_flags(
        fba_data_dir=fba_data_dir,
        fba_gram=fba_gram,
        fba_medium=fba_medium,
        fba_python_bin=fba_python_bin,
        fba_sidecar=fba_sidecar,
    )
    if len(stage2b_cmd) <= len(expected_extra) or stage2b_cmd[-len(expected_extra) :] != expected_extra:
        raise ValueError(
            f"stage2b command must end with documented Stage 2B flags {expected_extra}, got {stage2b_cmd}"
        )
    if "--enable-priors" in stage2b_cmd:
        raise ValueError("stage2b command must use stage1 default priors (no --enable-priors override)")

    base0 = stage0_cmd[:-2]
    base2b = stage2b_cmd[: -len(expected_extra)]
    if len(base0) != len(base2b):
        raise ValueError("stage0 and stage2b base commands have different lengths")

    # Every token in base0 and base2b must match after replacing the arm directory segment ('/stage0/' -> '/stage2b/')
    # or normalizing '-Djava.io.tmpdir=...'
    for i, (t0, t2b) in enumerate(zip(base0, base2b)):
        if t0 == t2b:
            continue
        if t0.replace("/stage0/", "/stage2b/") == t2b:
            continue
        if t0.startswith("-Djava.io.tmpdir=") and t2b.startswith("-Djava.io.tmpdir="):
            continue
        raise ValueError(f"Disallowed command difference at token {i}: stage0={t0!r} vs stage2b={t2b!r}")


def verify_stage2b_outputs(
    arm: str,
    output_tsv: Path,
    *,
    expected_gram: str,
    expected_medium: str,
) -> dict[str, str]:
    """Verify `.stages.json` and `.metabolism/` artifacts for `stage0` and `stage2b`."""
    ledger_path = Path(str(output_tsv) + ".stages.json")
    data = SL.ledger_metadata(ledger_path)
    stages = {s["name"]: s for s in data.get("stages", [])}
    met_stage = stages.get("metabolic_refinement")
    if met_stage is None:
        raise ValueError(f"{output_tsv}: missing metabolic_refinement stage in {ledger_path}")

    if arm == "stage0":
        if met_stage.get("status") != "skipped":
            raise ValueError(f"{output_tsv}: stage0 metabolic_refinement must be skipped, got {met_stage}")
        return {}

    if arm != "stage2b":
        raise ValueError(f"Unsupported arm {arm!r}")

    if met_stage.get("status") != "executed":
        raise ValueError(f"{output_tsv}: stage2b metabolic_refinement did not execute: {met_stage}")

    met_dir = Path(str(output_tsv) + ".metabolism")
    required_met_files = ("fba_input.json", "fba_output.json", "metabolic_model.xml", "reactions.tsv")
    met_paths = [met_dir / fn for fn in required_met_files]
    for mp in met_paths:
        if not mp.is_file():
            raise ValueError(f"{output_tsv}: missing Stage 2B metabolism artifact {mp}")
    fba_in = json.loads((met_dir / "fba_input.json").read_text(encoding="utf-8"))
    if fba_in.get("gram") != expected_gram or fba_in.get("medium") != expected_medium:
        raise ValueError(
            f"{output_tsv}: stage2b fba_input.json assumptions mismatch: "
            f"got gram={fba_in.get('gram')!r}, medium={fba_in.get('medium')!r}"
        )
    fba_out = json.loads((met_dir / "fba_output.json").read_text(encoding="utf-8"))
    if fba_out.get("max_univ_bio", 0.0) <= 0.0:
        raise ValueError(f"{output_tsv}: non-positive max_univ_bio in {met_dir / 'fba_output.json'}")
    return V.file_hashes(met_paths)


def _resolve_sidecar_path(fba_sidecar: Path | None) -> Path:
    if fba_sidecar is not None:
        return Path(fba_sidecar).resolve()
    cwd_sidecar = Path.cwd() / "scripts" / "stoichiometric_locus_refiner.py"
    if cwd_sidecar.is_file():
        return cwd_sidecar.resolve()
    return (REPO_ROOT / "scripts" / "stoichiometric_locus_refiner.py").resolve()


def _run_capacity_stage2b(
    part: Path,
    layout: Path,
    jar: Path,
    reference: Path,
    fba_data_dir: Path,
    fba_gram: str,
    fba_medium: str,
    fba_python_bin: Path,
    fba_sidecar: Path | None,
    out: Path,
    scratch: Path,
    reports: dict[str, Any],
    java_heap: str,
    arms: tuple[str, ...],
    timeout: int,
) -> dict[str, Any]:
    record = json.loads((part / "part.json").read_text())
    publication = json.loads((part / "job-receipt.json").read_text())
    if publication["status"] != "verified":
        raise ValueError("part is not verified")
    inputs = {str(part / name): digest for name, digest in publication["artifacts"].items()}
    if V.file_hashes(inputs) != inputs:
        raise ValueError("published part changed")
    asm = record["assembly"]
    claims = part / record["claims"]["path"]
    layout_record = json.loads((layout / "report.json").read_text())["assemblies"][asm]
    expected, count = U.expected_candidates(claims, scratch)
    if count != publication["claims"]:
        raise ValueError("published claim count differs from observed universe")

    sidecar_resolved = _resolve_sidecar_path(fba_sidecar)
    fba_data_files = [fba_data_dir / fn for fn in FBA_REQUIRED_DATA_FILES]
    for fp in fba_data_files:
        if not fp.is_file():
            raise ValueError(f"Missing required FBA data file: {fp}")

    commands_by_arm: dict[str, list[str]] = {}
    for arm in arms:
        directory = out / arm
        directory.mkdir()
        output = directory / (asm + ".integrated.tsv")
        manifest = directory / "manifest.tsv"
        with manifest.open("w") as fh:
            w = csv.DictWriter(fh, fieldnames=BATCH_HEADER, delimiter="\t")
            w.writeheader()
            row = dict.fromkeys(BATCH_HEADER, "-")
            row.update(
                genome_id=asm,
                claims_path=str(claims),
                out_path=str(output),
                taxon="NCBITaxon_" + str(layout_record["taxid"]),
                essential_profile="bacteria",
                genome_layout_path=str(layout / (asm + ".layout.tsv")),
            )
            w.writerow(row)
        command = [
            "/usr/bin/time",
            "-v",
            "-o",
            str(directory / "resources.txt"),
            "java",
            "-Xmx" + java_heap,
            "-Djava.io.tmpdir=" + str(scratch),
            "-jar",
            str(jar),
            "integrate",
            "--batch-manifest",
            str(manifest),
            "--batch-threads",
            "1",
            "--go-owl",
            str(reference / "go.owl"),
            "--ec2go",
            str(reference / "ec2go.txt"),
            "--pathways",
            str(reference / "kegg_pathways.tsv"),
            "--modules",
            str(reference / "kegg_modules.tsv"),
            "--reasoner-cache",
            str(reference / "reasoner-cache"),
            "--batch-stats-out",
            str(directory / "batch_stats.tsv"),
        ]
        sources = [
            Path(__file__),
            Path(C.__file__),
            part / "job-receipt.json",
            part / "matrix-verification.json",
            layout / "report.json",
            Path(shutil.which("java")),
            Path(U.__file__),
            U.sort_executable(),
        ]
        if arm == "stage0":
            command += ["--enable-priors", "none"]
        elif arm == "stage2b":
            command += stage2b_extra_flags(
                fba_data_dir=fba_data_dir,
                fba_gram=fba_gram,
                fba_medium=fba_medium,
                fba_python_bin=fba_python_bin,
                fba_sidecar=fba_sidecar,
            )
            sources.extend([sidecar_resolved, Path(fba_python_bin).resolve(), *fba_data_files])
        else:
            raise ValueError(f"Unsupported arm: {arm!r}")

        commands_by_arm[arm] = list(command)
        save(directory / "command.json", command)
        V.run_batch(command, manifest, directory, sources=sources, timeout=timeout)
        met_hashes = verify_stage2b_outputs(
            arm,
            output,
            expected_gram=fba_gram,
            expected_medium=fba_medium,
        )
        rows = U.compare_candidates(expected, output, scratch)
        metrics = {}
        for line in (directory / "resources.txt").read_text().splitlines():
            key, _, value = line.strip().partition(": ")
            metrics[key] = value
        reports[arm] = {
            "prediction_rows": rows,
            "input_candidates": count,
            "candidate_universe_preserved": True,
            "resources": metrics,
            "receipt_sha256": V.sha256_file(str(directory / V.RECEIPT)),
            "metabolism_artifacts": met_hashes,
        }

    if "stage0" in commands_by_arm and "stage2b" in commands_by_arm:
        verify_arm_command_pair(
            commands_by_arm["stage0"],
            commands_by_arm["stage2b"],
            fba_data_dir=str(fba_data_dir),
            fba_gram=fba_gram,
            fba_medium=fba_medium,
            fba_python_bin=str(fba_python_bin),
            fba_sidecar=str(fba_sidecar) if fba_sidecar is not None else None,
        )

    if V.file_hashes(inputs) != inputs:
        raise ValueError("part changed during capacity measurement")
    return {
        "schema": "gspa.price-stage2b-capacity/1",
        "status": "verified",
        "assembly": asm,
        "predictor": record["predictor"],
        "matrix_rows": publication["matrix_rows"],
        "runs": reports,
        "scope": "Complete matrix; Stage 0 vs Stage 2B (flux-certified refiner) inference.",
    }


def run_cell_capacity(
    part: Path,
    layout: Path,
    jar: Path,
    reference: Path,
    fba_data_dir: Path,
    *,
    fba_gram: str = "neg",
    fba_medium: str = "minimal",
    fba_python_bin: Path,
    fba_sidecar: Path | None = None,
    out: Path,
    scratch: Path | None = None,
    java_heap: str = "32g",
    arms: tuple[str, ...] = ARMS,
    timeout: int = 7200,
) -> dict[str, Any]:
    if not re.fullmatch(r"[1-9][0-9]*[gGmM]", java_heap):
        raise ValueError("java_heap must be a positive size in M or G")
    if not arms or len(set(arms)) != len(arms) or not set(arms) <= set(ARMS):
        raise ValueError("arms must be distinct stage0/stage2b names")
    if isinstance(timeout, bool) or not isinstance(timeout, int) or timeout <= 0:
        raise ValueError("timeout must be a positive integer number of seconds")
    if fba_gram != "neg":
        raise ValueError("All 22 Price genomes are Gram-negative Pseudomonadota; fba_gram must be 'neg'")
    scratch = Path(scratch or os.environ.get("TMPDIR") or Path.home() / "tmp").resolve()
    scratch.mkdir(parents=True, exist_ok=True)
    out.mkdir(parents=True, exist_ok=False)
    reports: dict[str, Any] = {}
    started = time.monotonic()
    try:
        with tempfile.TemporaryDirectory(prefix="gspa-stage2b-cap-", dir=scratch) as work:
            result = _run_capacity_stage2b(
                part,
                layout,
                jar,
                reference,
                fba_data_dir,
                fba_gram,
                fba_medium,
                fba_python_bin,
                fba_sidecar,
                out,
                Path(work),
                reports,
                java_heap,
                arms,
                timeout,
            )
    except Exception as error:
        save(
            out / "verification.json",
            {
                "schema": "gspa.price-stage2b-capacity/1",
                "status": "failed",
                "part": str(part),
                "runs": reports,
                "requested_arms": list(arms),
                "java_heap": java_heap,
                "timeout_seconds_per_arm": timeout,
                "total_elapsed_seconds": time.monotonic() - started,
                "error_type": type(error).__name__,
                "error": str(error),
            },
        )
        raise
    result.update(timeout_seconds_per_arm=timeout, total_elapsed_seconds=time.monotonic() - started)
    save(out / "verification.json", result)
    return result


def prepare(
    export_path: Path,
    audit_path: Path,
    jar: Path,
    reference: Path,
    fba_data_dir: Path,
    fba_python_bin: Path,
    out: Path,
    *,
    fba_gram: str = "neg",
    fba_medium: str = "minimal",
    fba_sidecar: Path | None = None,
    allow_synthetic_subset: bool = False,
) -> dict[str, Any]:
    if fba_gram != "neg":
        raise ValueError("Price Stage 2B collection requires --fba-gram neg (all 22 genomes are Gram-negative Pseudomonadota)")
    export = json.loads(export_path.read_text())
    audit = json.loads(audit_path.read_text())
    if audit["hashes"].get(str(export_path)) != V.sha256_file(str(export_path)):
        raise ValueError("export plan differs from audited plan")
    rows = cells(export, audit)
    if not allow_synthetic_subset:
        if (
            len(rows) != 132
            or len({r["assembly"] for r in rows}) != 22
            or {r["predictor"] for r in rows} != {"clean", "dpz", "graphec", "topec", "mapred", "enzbert"}
        ):
            raise ValueError("the complete six-predictor, 22-assembly Price plan is required")

    layout = Path(export["layout"])
    layout_report = json.loads((layout / "report.json").read_text())
    tax_summary = verify_price22_gram_negative_pseudomonadota(
        layout_report, allow_synthetic_subset=allow_synthetic_subset
    )

    for fn in FBA_REQUIRED_DATA_FILES:
        if not (fba_data_dir / fn).is_file():
            raise ValueError(f"fba_data_dir is missing required file {fn}: {fba_data_dir / fn}")

    sidecar_resolved = _resolve_sidecar_path(fba_sidecar)
    if not sidecar_resolved.is_file():
        raise FileNotFoundError(f"Stage 2B sidecar script not found: {sidecar_resolved}")
    if not Path(fba_python_bin).is_file():
        raise FileNotFoundError(f"Stage 2B Python interpreter not found: {fba_python_bin}")

    files = [
        export_path,
        audit_path,
        jar,
        Path(shutil.which("java")),
        Path(sys.executable),
        Path(fba_python_bin).resolve(),
        sidecar_resolved,
        Path(__file__).resolve(),
    ]
    files += [p for p in V3_DIR.glob("*.py")]
    files += V.streaming_parser_sources()
    files += [V3_DIR / "requirements-streaming.txt"]
    files += [reference / name for name in ("go.owl", "ec2go.txt", "kegg_pathways.tsv", "kegg_modules.tsv")]
    files += [p for p in (reference / "reasoner-cache").rglob("*") if p.is_file()]
    files += [p for p in fba_data_dir.rglob("*") if p.is_file()]
    files += [layout / "report.json"]

    for row in rows:
        p = Path(row["part"]) / "job-receipt.json"
        if V.sha256_file(str(p)) != row["publication_sha256"]:
            raise ValueError("publication changed since complete audit")
        layout_path = layout / (row["assembly"] + ".layout.tsv")
        if audit["hashes"].get(str(layout_path)) != V.sha256_file(str(layout_path)):
            raise ValueError("audited genome layout changed")
        files += [p, layout_path]

    inputs = V.file_hashes(files)
    out.mkdir(parents=True, exist_ok=False)
    plan = dict(
        schema="gspa.price-stage2b-production/1",
        cells=rows,
        arms=list(ARMS),
        jar=str(jar),
        layout=str(layout),
        reference=str(reference),
        fba_data_dir=str(fba_data_dir),
        fba_gram=fba_gram,
        fba_medium=fba_medium,
        fba_python_bin=str(fba_python_bin),
        fba_sidecar=str(fba_sidecar) if fba_sidecar is not None else None,
        fba_sidecar_resolved=str(sidecar_resolved),
        taxonomy_verification=tax_summary,
        java=str(Path(shutil.which("java")).resolve()),
        python=str(Path(sys.executable).resolve()),
        expected_claims=export["expected_claims"],
        scope=__doc__,
        model_fitting=False,
        labels_read=False,
        interpretation=(
            "stage0 is production integration without priors; stage2b enables default priors plus "
            "Stage 2B flux-certified metabolic refinement (--enable-fba-refiner) on Gram-negative "
            "minimal medium."
        ),
    )
    save(out / "plan.json", plan)
    save(
        out / "registration.json",
        dict(
            registered_at=datetime.now(timezone.utc).isoformat(),
            plan_sha256=V.sha256_file(str(out / "plan.json")),
            inputs=inputs,
        ),
    )
    return dict(cells=len(rows), expected_claims=plan["expected_claims"], arms=list(ARMS))


def run(out: Path, index: int, attempt: str, heap: str, timeout: int, scratch: Path) -> None:
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", attempt) or attempt in (".", ".."):
        raise ValueError("invalid attempt identifier")
    reg = json.loads((out / "registration.json").read_text())
    check(reg["inputs"])
    if reg["plan_sha256"] != V.sha256_file(str(out / "plan.json")):
        raise ValueError("production plan changed")
    plan = json.loads((out / "plan.json").read_text())
    if (
        plan.get("java", str(Path(shutil.which("java")).resolve())) != str(Path(shutil.which("java")).resolve())
        or plan.get("python", str(Path(sys.executable).resolve())) != str(Path(sys.executable).resolve())
    ):
        raise ValueError("collection runtime executable changed")
    if any(os.environ.get(k) for k in ("JAVA_TOOL_OPTIONS", "JDK_JAVA_OPTIONS", "_JAVA_OPTIONS")):
        raise ValueError("collection requires explicit Java options without environment overrides")
    if type(index) is not int or not 0 <= index < len(plan["cells"]):
        raise ValueError("unknown collection index")

    cell = plan["cells"][index]
    destination = out / "attempts" / str(index) / attempt
    destination.mkdir(parents=True, exist_ok=False)
    save(
        destination / "started.json",
        dict(
            cell=cell,
            attempt=attempt,
            heap=heap,
            timeout=timeout,
            started_at=datetime.now(timezone.utc).isoformat(),
            python=sys.version,
            registration_sha256=V.sha256_file(str(out / "registration.json")),
        ),
    )
    try:
        fba_sidecar = Path(plan["fba_sidecar"]) if plan.get("fba_sidecar") else None
        run_cell_capacity(
            Path(cell["part"]),
            Path(plan["layout"]),
            Path(plan["jar"]),
            Path(plan["reference"]),
            Path(plan["fba_data_dir"]),
            fba_gram=plan["fba_gram"],
            fba_medium=plan["fba_medium"],
            fba_python_bin=Path(plan["fba_python_bin"]),
            fba_sidecar=fba_sidecar,
            out=destination / "production",
            scratch=scratch,
            java_heap=heap,
            arms=tuple(plan["arms"]),
            timeout=timeout,
        )
        verified = json.loads((destination / "production/verification.json").read_text())
        if verified["status"] != "verified" or set(verified["runs"]) != set(plan["arms"]):
            raise ValueError("incomplete production arms")
        if any(verified[k] != cell[k] for k in ("assembly", "predictor", "matrix_rows")):
            raise ValueError("production cell identity differs")
        if any(
            r["prediction_rows"] != cell["claims"]
            or r["input_candidates"] != cell["claims"]
            or not r["candidate_universe_preserved"]
            for r in verified["runs"].values()
        ):
            raise ValueError("registered candidate population differs")
        check(reg["inputs"])
        save(
            destination / "completed.json",
            dict(
                status="verified",
                cell=cell,
                verification_sha256=V.sha256_file(str(destination / "production/verification.json")),
                finished_at=datetime.now(timezone.utc).isoformat(),
            ),
        )
    except Exception as error:
        save(
            destination / "failed.json",
            dict(
                status="failed",
                error_type=type(error).__name__,
                error=str(error),
                finished_at=datetime.now(timezone.utc).isoformat(),
            ),
        )
        (destination / "failure.log").write_text(traceback.format_exc())
        raise


def finalize(out: Path) -> None:
    """Recheck every input, command contract, stage ledger, and metabolism artifact before certifying."""
    reg = json.loads((out / "registration.json").read_text())
    if reg["plan_sha256"] != V.sha256_file(str(out / "plan.json")):
        raise ValueError("production plan changed")
    plan = json.loads((out / "plan.json").read_text())
    if tuple(plan.get("arms", ())) != ARMS:
        raise ValueError(f"Expected arms {ARMS}, got {plan.get('arms')}")
    selected = completed_attempts(out, plan)
    inputs = dict(reg["inputs"])
    artifacts: dict[str, str] = {}
    reports = []

    def bind(path: Path, expected: str | None = None) -> None:
        actual = V.sha256_file(str(path))
        if expected is not None and actual != expected:
            raise ValueError(f"collection receipt changed: {path}")
        artifacts[str(path)] = actual

    bind(out / "plan.json", reg["plan_sha256"])
    bind(out / "registration.json")

    for cell, directory in selected:
        if (directory / "failed.json").exists():
            raise ValueError("contradictory attempt status")
        start = json.loads((directory / "started.json").read_text())
        done = json.loads((directory / "completed.json").read_text())
        if (
            start["cell"] != cell
            or done["cell"] != cell
            or done["status"] != "verified"
            or start["registration_sha256"] != artifacts[str(out / "registration.json")]
        ):
            raise ValueError("attempt registration/identity mismatch")
        bind(directory / "started.json")
        bind(directory / "completed.json")
        root = directory / "production"
        path = root / "verification.json"
        bind(path, done["verification_sha256"])
        result = json.loads(path.read_text())
        if (
            result["status"] != "verified"
            or set(result["runs"]) != set(plan["arms"])
            or any(result[k] != cell[k] for k in ("assembly", "predictor", "matrix_rows"))
        ):
            raise ValueError("incomplete or mismatched production cell")

        commands_by_arm: dict[str, list[str]] = {}
        for arm, report in result["runs"].items():
            if (
                report["prediction_rows"] != cell["claims"]
                or report["input_candidates"] != cell["claims"]
                or report["candidate_universe_preserved"] is not True
            ):
                raise ValueError("candidate population changed")
            receipt_path = root / arm / V.RECEIPT
            bind(receipt_path, report["receipt_sha256"])
            receipt = json.loads(receipt_path.read_text())
            spec = receipt["spec"]
            if (
                receipt["returncode"] != 0
                or set(spec["outputs"]) != {cell["assembly"]}
                or Path(spec["outputs"][cell["assembly"]]).parent != root / arm
            ):
                raise ValueError("invalid arm receipt")
            commands_by_arm[arm] = list(spec["command"])
            if V.output_hashes(spec) != receipt["artifacts"]:
                raise ValueError("prediction or ledger changed")
            artifacts.update(receipt["artifacts"])

            output_tsv = Path(spec["outputs"][cell["assembly"]])
            met_hashes = verify_stage2b_outputs(
                arm,
                output_tsv,
                expected_gram=plan["fba_gram"],
                expected_medium=plan["fba_medium"],
            )
            if report.get("metabolism_artifacts", {}) != met_hashes:
                raise ValueError("Stage 2B metabolism artifacts changed since cell completion")
            artifacts.update(met_hashes)

            for in_path, digest in spec["inputs"].items():
                if in_path in inputs and inputs[in_path] != digest:
                    raise ValueError("inconsistent inference input hashes")
                inputs[in_path] = digest

        verify_arm_command_pair(
            commands_by_arm["stage0"],
            commands_by_arm["stage2b"],
            fba_data_dir=plan["fba_data_dir"],
            fba_gram=plan["fba_gram"],
            fba_medium=plan["fba_medium"],
            fba_python_bin=plan["fba_python_bin"],
            fba_sidecar=plan.get("fba_sidecar"),
        )

        publication_path = directory / "publication.json"
        if publication_path.exists():
            bind(publication_path)
            publication = json.loads(publication_path.read_text())
            if (
                publication.get("schema") != "gspa.local-publication/1"
                or publication.get("status") != "verified"
                or publication.get("published_path") != str(directory.resolve())
                or publication.get("execution", {}).get("returncode") != 0
                or publication["execution"].get("registration_sha256") != artifacts[str(out / "registration.json")]
            ):
                raise ValueError("invalid local publication identity")
            inventory = publication["artifacts"]
            actual = {str(p.relative_to(directory)) for p in directory.rglob("*") if p.is_file()}
            if actual != set(inventory) | {"publication.json"}:
                raise ValueError("local publication inventory changed")
            for relative, digest in inventory.items():
                p_rel = directory / relative
                if Path(relative).is_absolute() or ".." in Path(relative).parts or p_rel.is_symlink():
                    raise ValueError("invalid local publication artifact path")
                p_rel.resolve().relative_to(directory.resolve())
                if str(p_rel) in artifacts:
                    if artifacts[str(p_rel)] != digest:
                        raise ValueError("local publication artifact changed")
                else:
                    bind(p_rel, digest)

        reports.append(dict(cell=cell, attempt=directory.name, arms=result["runs"]))

    check(inputs)
    with (out / "collection-verification.json").open("x") as stream:
        json.dump(
            dict(
                schema="gspa.price-stage2b-collection-verification/1",
                status="verified",
                cells=len(reports),
                arms=len(plan["arms"]),
                expected_claims=plan["expected_claims"],
                fba_data_dir=plan["fba_data_dir"],
                fba_gram=plan["fba_gram"],
                fba_medium=plan["fba_medium"],
                fba_python_bin=plan["fba_python_bin"],
                taxonomy_verification=plan["taxonomy_verification"],
                inputs=inputs,
                artifacts=artifacts,
                reports=reports,
                parser=V.parser_identity(),
                interpretation="Complete Stage 2B (flux-certified refiner) vs Stage 0 Price collection.",
            ),
            stream,
            indent=2,
            sort_keys=True,
        )
        stream.write("\n")


def prepare_external(
    out: Path,
    index: int,
    attempt: str,
    heap: str,
    timeout: int,
    scratch: Path,
    container: Path,
    prefix_path: Path,
) -> Path:
    """Reuse `run_price_local.prepare_external` pointing `L.P` at this module."""
    L.P = sys.modules[__name__]
    return L.prepare_external(out, index, attempt, heap, timeout, scratch, container, prefix_path)


def publish_external(
    out: Path,
    index: int,
    attempt: str,
    scratch: Path,
    container: Path,
    returncode: int,
) -> None:
    """Reuse `run_price_local.publish_external` pointing `L.P` at this module."""
    L.P = sys.modules[__name__]
    L.publish_external(out, index, attempt, scratch, container, returncode)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "action",
        choices=["prepare", "run", "finalize", "prepare-external", "publish-external"],
    )
    p.add_argument("--out", type=Path, required=True)
    for name in (
        "export",
        "audit",
        "jar",
        "reference",
        "fba-data-dir",
        "fba-python-bin",
        "fba-sidecar",
        "scratch",
        "container",
        "launch-prefix",
    ):
        p.add_argument("--" + name, type=Path)
    p.add_argument("--fba-gram", default="neg")
    p.add_argument("--fba-medium", default="minimal")
    p.add_argument("--allow-synthetic-subset", action="store_true")
    p.add_argument("--index", type=int)
    p.add_argument("--attempt")
    p.add_argument("--heap", default="24g")
    p.add_argument("--timeout", type=int, default=7200)
    p.add_argument("--returncode", type=int)
    args = p.parse_args(argv)

    if args.action == "prepare":
        if not all((args.export, args.audit, args.jar, args.reference, args.fba_data_dir, args.fba_python_bin)):
            p.error("prepare requires --export, --audit, --jar, --reference, --fba-data-dir, --fba-python-bin")
        res = prepare(
            args.export.resolve(),
            args.audit.resolve(),
            args.jar.resolve(),
            args.reference.resolve(),
            args.fba_data_dir.resolve(),
            args.fba_python_bin.resolve(),
            args.out.resolve(),
            fba_gram=args.fba_gram,
            fba_medium=args.fba_medium,
            fba_sidecar=args.fba_sidecar.resolve() if args.fba_sidecar else None,
            allow_synthetic_subset=args.allow_synthetic_subset,
        )
        print(json.dumps(res))
    elif args.action == "finalize":
        finalize(args.out.resolve())
    elif args.action == "run":
        if args.index is None or not args.attempt or not args.scratch:
            p.error("run requires --index, --attempt, --scratch")
        run(args.out.resolve(), args.index, args.attempt, args.heap, args.timeout, args.scratch.resolve())
    elif args.action == "prepare-external":
        if args.index is None or not args.attempt or not args.scratch or not args.container or not args.launch_prefix:
            p.error("prepare-external requires --index, --attempt, --scratch, --container, --launch-prefix")
        loc = prepare_external(
            args.out.resolve(),
            args.index,
            args.attempt,
            args.heap,
            args.timeout,
            args.scratch.resolve(),
            args.container.resolve(),
            args.launch_prefix.resolve(),
        )
        print(loc, flush=True)
    elif args.action == "publish-external":
        if args.index is None or not args.attempt or not args.scratch or not args.container or args.returncode is None:
            p.error("publish-external requires --index, --attempt, --scratch, --container, --returncode")
        publish_external(
            args.out.resolve(),
            args.index,
            args.attempt,
            args.scratch.resolve(),
            args.container.resolve(),
            args.returncode,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
