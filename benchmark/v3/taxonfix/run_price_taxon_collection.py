#!/usr/bin/env python3
"""Price 132-cell collection runner with the taxon/EC consistency fix (`stage0` vs `stage1_taxon`).

Reuses functions from the frozen Price engine modules (`run_price_production_collection`,
`run_price_matrix_capacity`, `run_price_local`, `verified_batch`, `verify_candidate_universe`,
`stage_ledger`, `v3_runner`) by import without modifying any frozen file.

Arms:
  * `stage0`: unchanged production command (`--enable-priors none`).
  * `stage1_taxon`: default priors + `--enforce-stage2 --enforce-consistency-mode remove --taxon-lineage <lineage_tsv>`.
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

ARMS = ("stage0", "stage1_taxon")
DEFAULT_LINEAGE_TSV = Path(__file__).resolve().parent / "price22_taxon_lineages.tsv"
PRICE22_LINEAGE_SHA256 = "f3eeedb204d5403aeb394df1832f6eb6ac4a2905078ebe70d5390b9fa0415ab2"

# All 22 Price genomes with their NCBI taxids (18 distinct strain/species taxids),
# all belonging to Bacteria (NCBITaxon_2) and Pseudomonadota (NCBITaxon_1224).
PRICE22_TAXONOMY: dict[str, dict[str, Any]] = {
    "GCA_000006965.1": {
        "taxid": 266834,
        "organism_name": "Sinorhizobium meliloti 1021",
        "phylum": "Pseudomonadota",
        "phylum_taxid": 1224,
        "class_name": "Alphaproteobacteria",
        "order": "Hyphomicrobiales",
        "family": "Rhizobiaceae",
    },
    "GCA_003443575.1": {
        "taxid": 221822,
        "organism_name": "Phaeobacter inhibens",
        "phylum": "Pseudomonadota",
        "phylum_taxid": 1224,
        "class_name": "Alphaproteobacteria",
        "order": "Rhodobacterales",
        "family": "Roseobacteraceae",
    },
    "GCA_004402975.1": {
        "taxid": 1211383,
        "organism_name": "Paraburkholderia dipogonis",
        "phylum": "Pseudomonadota",
        "phylum_taxid": 1224,
        "class_name": "Betaproteobacteria",
        "order": "Burkholderiales",
        "family": "Burkholderiaceae",
    },
    "GCF_000015245.1": {
        "taxid": 326297,
        "organism_name": "Shewanella amazonensis SB2B",
        "phylum": "Pseudomonadota",
        "phylum_taxid": 1224,
        "class_name": "Gammaproteobacteria",
        "order": "Alteromonadales",
        "family": "Shewanellaceae",
    },
    "GCF_000016065.1": {
        "taxid": 323850,
        "organism_name": "Shewanella loihica PV-4",
        "phylum": "Pseudomonadota",
        "phylum_taxid": 1224,
        "class_name": "Gammaproteobacteria",
        "order": "Alteromonadales",
        "family": "Shewanellaceae",
    },
    "GCF_000018145.1": {
        "taxid": 398580,
        "organism_name": "Dinoroseobacter shibae DFL 12 = DSM 16493",
        "phylum": "Pseudomonadota",
        "phylum_taxid": 1224,
        "class_name": "Alphaproteobacteria",
        "order": "Rhodobacterales",
        "family": "Roseobacteraceae",
    },
    "GCF_000020125.1": {
        "taxid": 398527,
        "organism_name": "Paraburkholderia phytofirmans PsJN",
        "phylum": "Pseudomonadota",
        "phylum_taxid": 1224,
        "class_name": "Betaproteobacteria",
        "order": "Burkholderiales",
        "family": "Burkholderiaceae",
    },
    "GCF_000146165.2": {
        "taxid": 211586,
        "organism_name": "Shewanella oneidensis MR-1",
        "phylum": "Pseudomonadota",
        "phylum_taxid": 1224,
        "class_name": "Gammaproteobacteria",
        "order": "Alteromonadales",
        "family": "Shewanellaceae",
    },
    "GCF_000166295.1": {
        "taxid": 225937,
        "organism_name": "Marinobacter adhaerens HP15",
        "phylum": "Pseudomonadota",
        "phylum_taxid": 1224,
        "class_name": "Gammaproteobacteria",
        "order": "Alteromonadales",
        "family": "Marinobacteraceae",
    },
    "GCF_000203935.1": {
        "taxid": 94122,
        "organism_name": "Shewanella sp. ANA-3",
        "phylum": "Pseudomonadota",
        "phylum_taxid": 1224,
        "class_name": "Gammaproteobacteria",
        "order": "Alteromonadales",
        "family": "Shewanellaceae",
    },
    "GCF_000236665.1": {
        "taxid": 640081,
        "organism_name": "Azospira oryzae PS",
        "phylum": "Pseudomonadota",
        "phylum_taxid": 1224,
        "class_name": "Betaproteobacteria",
        "order": "Rhodocyclales",
        "family": "Rhodocyclaceae",
    },
    "GCF_000237365.1": {
        "taxid": 1064539,
        "organism_name": "Azospirillum baldaniorum",
        "phylum": "Pseudomonadota",
        "phylum_taxid": 1224,
        "class_name": "Alphaproteobacteria",
        "order": "Rhodospirillales",
        "family": "Azospirillaceae",
    },
    "GCF_000327065.1": {
        "taxid": 644801,
        "organism_name": "Stutzerimonas stutzeri RCH2",
        "phylum": "Pseudomonadota",
        "phylum_taxid": 1224,
        "class_name": "Gammaproteobacteria",
        "order": "Pseudomonadales",
        "family": "Pseudomonadaceae",
    },
    "GCF_000698265.1": {
        "taxid": 321846,
        "organism_name": "Pseudomonas simiae",
        "phylum": "Pseudomonadota",
        "phylum_taxid": 1224,
        "class_name": "Gammaproteobacteria",
        "order": "Pseudomonadales",
        "family": "Pseudomonadaceae",
    },
    "GCF_000832305.1": {
        "taxid": 68895,
        "organism_name": "Cupriavidus basilensis",
        "phylum": "Pseudomonadota",
        "phylum_taxid": 1224,
        "class_name": "Betaproteobacteria",
        "order": "Burkholderiales",
        "family": "Burkholderiaceae",
    },
    "GCF_001307155.1": {
        "taxid": 294,
        "organism_name": "Pseudomonas fluorescens",
        "phylum": "Pseudomonadota",
        "phylum_taxid": 1224,
        "class_name": "Gammaproteobacteria",
        "order": "Pseudomonadales",
        "family": "Pseudomonadaceae",
    },
    "GCF_001307275.1": {
        "taxid": 294,
        "organism_name": "Pseudomonas fluorescens",
        "phylum": "Pseudomonadota",
        "phylum_taxid": 1224,
        "class_name": "Gammaproteobacteria",
        "order": "Pseudomonadales",
        "family": "Pseudomonadaceae",
    },
    "GCF_001440325.1": {
        "taxid": 294,
        "organism_name": "Pseudomonas fluorescens",
        "phylum": "Pseudomonadota",
        "phylum_taxid": 1224,
        "class_name": "Gammaproteobacteria",
        "order": "Pseudomonadales",
        "family": "Pseudomonadaceae",
    },
    "GCF_001623525.1": {
        "taxid": 294,
        "organism_name": "Pseudomonas fluorescens",
        "phylum": "Pseudomonadota",
        "phylum_taxid": 1224,
        "class_name": "Gammaproteobacteria",
        "order": "Pseudomonadales",
        "family": "Pseudomonadaceae",
    },
    "GCF_001625455.1": {
        "taxid": 294,
        "organism_name": "Pseudomonas fluorescens",
        "phylum": "Pseudomonadota",
        "phylum_taxid": 1224,
        "class_name": "Gammaproteobacteria",
        "order": "Pseudomonadales",
        "family": "Pseudomonadaceae",
    },
    "GCF_001633105.1": {
        "taxid": 1813946,
        "organism_name": "Acidovorax sp. GW101-3H11",
        "phylum": "Pseudomonadota",
        "phylum_taxid": 1224,
        "class_name": "Betaproteobacteria",
        "order": "Burkholderiales",
        "family": "Comamonadaceae",
    },
    "GCF_002797435.1": {
        "taxid": 93064,
        "organism_name": "Sphingomonas koreensis",
        "phylum": "Pseudomonadota",
        "phylum_taxid": 1224,
        "class_name": "Alphaproteobacteria",
        "order": "Sphingomonadales",
        "family": "Sphingomonadaceae",
    },
}


def _normalize_taxon_id(raw: str) -> str:
    t = raw.strip()
    if t.startswith("NCBITaxon:"):
        return "NCBITaxon_" + t.split(":", 1)[1]
    if t.isdigit():
        return "NCBITaxon_" + t
    return t


def parse_taxon_lineage_tsv(lineage_path: Path) -> dict[str, dict[str, Any]]:
    """Parse and validate a `--taxon-lineage` TSV file."""
    path = Path(lineage_path)
    if not path.is_file():
        raise FileNotFoundError(f"Taxon lineage TSV not found: {path}")
    records: dict[str, dict[str, Any]] = {}
    with path.open("r", encoding="utf-8") as fh:
        for line_no, raw_line in enumerate(fh, start=1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            cols = [c.strip() for c in line.split("\t")]
            if cols[0].lower() in ("taxon_id", "taxon", "taxid"):
                continue
            if len(cols) < 2 or not cols[0] or not cols[1]:
                raise ValueError(f"{path}:{line_no}: expected taxon_id<TAB>lineage_ids[<TAB>organism_name]")
            taxon_id = _normalize_taxon_id(cols[0])
            lineage_ids = [_normalize_taxon_id(x) for x in cols[1].split(",") if x.strip()]
            if not lineage_ids:
                raise ValueError(f"{path}:{line_no}: empty lineage_ids for {taxon_id}")
            if taxon_id in records:
                raise ValueError(f"{path}:{line_no}: duplicate taxon_id {taxon_id}")
            records[taxon_id] = {
                "taxon_id": taxon_id,
                "lineage_ids": lineage_ids,
                "organism_name": cols[2] if len(cols) >= 3 else "",
            }
    if not records:
        raise ValueError(f"{path}: taxon lineage TSV contains no taxon rows")
    return records


def verify_price22_taxon_lineages(
    layout_report: dict[str, Any],
    lineage_path: Path,
    *,
    allow_synthetic_subset: bool = False,
) -> dict[str, Any]:
    """Verify that every Price assembly's taxid is covered in `lineage_path` and links to Bacteria/Pseudomonadota."""
    assemblies = layout_report.get("assemblies", {})
    if not assemblies:
        raise ValueError("layout/report.json contains no assemblies")

    if not allow_synthetic_subset and set(assemblies.keys()) != set(PRICE22_TAXONOMY.keys()):
        raise ValueError(
            f"layout/report.json assemblies do not match the 22 Price genomes: "
            f"got {sorted(assemblies.keys())}"
        )

    lineage_table = parse_taxon_lineage_tsv(lineage_path)
    lineage_sha256 = V.sha256_file(str(lineage_path))
    if not allow_synthetic_subset and lineage_sha256 != PRICE22_LINEAGE_SHA256:
        raise ValueError(
            f"price22_taxon_lineages.tsv sha256 mismatch: got {lineage_sha256}, "
            f"expected {PRICE22_LINEAGE_SHA256}"
        )

    verified_genomes: dict[str, dict[str, Any]] = {}
    distinct_taxids: set[str] = set()
    for asm, rec in sorted(assemblies.items()):
        taxid = int(rec["taxid"])
        taxon_curie = f"NCBITaxon_{taxid}"
        distinct_taxids.add(taxon_curie)
        if asm in PRICE22_TAXONOMY:
            expected = PRICE22_TAXONOMY[asm]
            if taxid != expected["taxid"]:
                raise ValueError(
                    f"Assembly {asm} taxid mismatch in layout/report.json: "
                    f"got {taxid}, expected {expected['taxid']}"
                )
            meta = dict(expected)
        elif allow_synthetic_subset:
            meta = {
                "taxid": taxid,
                "organism_name": f"synthetic_{asm}",
                "phylum": "Pseudomonadota",
                "phylum_taxid": 1224,
                "class_name": "Gammaproteobacteria",
            }
        else:
            raise ValueError(f"Unknown Price assembly {asm!r}")

        if taxon_curie not in lineage_table:
            raise ValueError(f"Assembly {asm} taxon {taxon_curie} is missing from {lineage_path}")
        ancestors = lineage_table[taxon_curie]["lineage_ids"]
        if "NCBITaxon_2" not in ancestors or "NCBITaxon_1224" not in ancestors:
            raise ValueError(
                f"Assembly {asm} taxon {taxon_curie} lineage does not link to "
                f"NCBITaxon_2 (Bacteria) and NCBITaxon_1224 (Pseudomonadota): {ancestors}"
            )
        meta["taxon_id"] = taxon_curie
        meta["lineage_ids"] = list(ancestors)
        verified_genomes[asm] = meta

    return {
        "all_in_lineage_table": True,
        "all_bacteria": True,
        "all_pseudomonadota": True,
        "phylum": "Pseudomonadota",
        "phylum_taxid": 1224,
        "taxon_lineage_path": str(Path(lineage_path).resolve()),
        "taxon_lineage_sha256": lineage_sha256,
        "n_genomes": len(verified_genomes),
        "n_distinct_taxids": len(distinct_taxids),
        "n_lineage_rows": len(lineage_table),
        "genomes": verified_genomes,
    }


def stage1_taxon_extra_flags(taxon_lineage: Path | str) -> list[str]:
    return [
        "--enforce-stage2",
        "--enforce-consistency-mode",
        "remove",
        "--taxon-lineage",
        str(taxon_lineage),
    ]


def verify_arm_command_pair(
    stage0_cmd: list[str],
    stage1_taxon_cmd: list[str],
    *,
    taxon_lineage: str,
) -> None:
    """Enforce that stage0 and stage1_taxon differ ONLY in arm paths and the documented taxon-fix flags."""
    if len(stage0_cmd) < 2 or stage0_cmd[-2:] != ["--enable-priors", "none"]:
        raise ValueError(f"stage0 command must end with ['--enable-priors', 'none']: {stage0_cmd}")

    expected_extra = stage1_taxon_extra_flags(taxon_lineage)
    if (
        len(stage1_taxon_cmd) <= len(expected_extra)
        or stage1_taxon_cmd[-len(expected_extra) :] != expected_extra
    ):
        raise ValueError(
            f"stage1_taxon command must end with documented taxon-fix flags {expected_extra}, "
            f"got {stage1_taxon_cmd}"
        )
    if "--enable-priors" in stage1_taxon_cmd:
        raise ValueError("stage1_taxon command must use stage1 default priors (no --enable-priors override)")

    base0 = stage0_cmd[:-2]
    base1 = stage1_taxon_cmd[: -len(expected_extra)]
    if len(base0) != len(base1):
        raise ValueError("stage0 and stage1_taxon base commands have different lengths")

    for i, (t0, t1) in enumerate(zip(base0, base1)):
        if t0 == t1:
            continue
        if t0.replace("/stage0/", "/stage1_taxon/") == t1:
            continue
        if t0.startswith("-Djava.io.tmpdir=") and t1.startswith("-Djava.io.tmpdir="):
            continue
        raise ValueError(f"Disallowed command difference at token {i}: stage0={t0!r} vs stage1_taxon={t1!r}")


def _read_stage_ledger_full(ledger_path: Path) -> dict[str, Any]:
    """Read the stage records of `.stages.json` with every SCALAR field, streaming.

    Dense cells write ledgers of hundreds of MB (per-candidate change lists), so the file is never loaded whole:
    `SL.ledger_metadata` validates the complete JSON first, then a second streaming pass keeps each stage's scalar
    fields (status, taxon state, counts) and skips lists and nested objects. Small ledgers without ijson fall back to
    json.loads.
    """
    SL.ledger_metadata(ledger_path)
    try:
        import ijson
    except ImportError:
        if ledger_path.stat().st_size > SL.LEDGER_MEMORY_LIMIT:
            raise ValueError(f"{ledger_path}: large stage ledger requires ijson")
        data = json.loads(ledger_path.read_text(encoding="utf-8"))
        return {"stages": [{k: v for k, v in st.items() if not isinstance(v, (list, dict))}
                           for st in data.get("stages", [])]}
    stages: list[dict[str, Any]] = []
    with ledger_path.open("rb") as stream:
        events = iter(ijson.basic_parse(stream))
        def skip(first):
            depth = int(first[0] in ("start_map", "start_array"))
            while depth:
                ev, _ = next(events)
                if ev in ("start_map", "start_array"):
                    depth += 1
                elif ev in ("end_map", "end_array"):
                    depth -= 1
        def stage_obj(first):
            if first[0] != "start_map":
                raise ValueError(f"{ledger_path}: stage record is not an object")
            rec: dict[str, Any] = {}
            while True:
                ev, key = next(events)
                if ev == "end_map":
                    return rec
                val = next(events)
                if val[0] in ("string", "number", "boolean", "null"):
                    rec[key] = float(val[1]) if hasattr(val[1], "as_tuple") else val[1]
                else:
                    skip(val)
        first = next(events)
        if first[0] != "start_map":
            raise ValueError(f"{ledger_path}: ledger is not an object")
        while True:
            ev, key = next(events)
            if ev == "end_map":
                break
            val = next(events)
            if key == "stages" and val[0] == "start_array":
                nxt = next(events)
                while nxt[0] != "end_array":
                    stages.append(stage_obj(nxt))
                    nxt = next(events)
            else:
                skip(val)
    for st in stages:
        for k, v in list(st.items()):
            if isinstance(v, float) and v.is_integer():
                st[k] = int(v)
    return {"stages": stages}


def verify_taxon_outputs(
    arm: str,
    output_tsv: Path,
    *,
    expected_taxon: str | None = None,
) -> dict[str, Any]:
    """Verify `.stages.json` and `batch_stats.tsv` for `stage0` and `stage1_taxon`."""
    ledger_path = Path(str(output_tsv) + ".stages.json")
    data = _read_stage_ledger_full(ledger_path)
    stages = {s["name"]: s for s in data.get("stages", [])}
    ref_stage = stages.get("reference_preparation")
    ev_stage = stages.get("evidence_integration")
    sym_stage = stages.get("symbolic_enforcement")
    if ref_stage is None or ref_stage.get("status") != "executed":
        raise ValueError(f"{output_tsv}: reference_preparation did not execute in {ledger_path}")
    if ev_stage is None or ev_stage.get("status") != "executed":
        raise ValueError(f"{output_tsv}: evidence_integration did not execute in {ledger_path}")
    if sym_stage is None:
        raise ValueError(f"{output_tsv}: missing symbolic_enforcement stage in {ledger_path}")

    if arm == "stage0":
        if sym_stage.get("status") != "skipped":
            raise ValueError(f"{output_tsv}: stage0 symbolic_enforcement must be skipped, got {sym_stage}")
        return {
            "symbolic_enforcement_status": "skipped",
            "removed_candidates": 0,
            "ec_violations": 0,
            "violations": 0,
        }

    if arm != "stage1_taxon":
        raise ValueError(f"Unsupported arm {arm!r}")

    if sym_stage.get("status") != "executed":
        raise ValueError(f"{output_tsv}: stage1_taxon symbolic_enforcement did not execute: {sym_stage}")

    if expected_taxon is not None and ref_stage.get("organism_taxon") not in (None, expected_taxon):
        raise ValueError(
            f"{output_tsv}: reference_preparation organism_taxon mismatch: "
            f"got {ref_stage.get('organism_taxon')!r}, expected {expected_taxon!r}"
        )
    if ref_stage.get("taxon_in_hierarchy") is not True:
        raise ValueError(
            f"{output_tsv}: reference_preparation taxon_in_hierarchy must be True, got {ref_stage}"
        )
    if ref_stage.get("taxon_consistency_state") != "applied":
        raise ValueError(
            f"{output_tsv}: reference_preparation taxon_consistency_state must be 'applied', got {ref_stage}"
        )
    if ref_stage.get("taxon_not_in_hierarchy_count", 0) != 0:
        raise ValueError(
            f"{output_tsv}: reference_preparation taxon_not_in_hierarchy_count must be 0, got {ref_stage}"
        )

    if sym_stage.get("taxon_in_hierarchy") is not True:
        raise ValueError(f"{output_tsv}: symbolic_enforcement taxon_in_hierarchy must be True: {sym_stage}")
    if sym_stage.get("constraints_applied") is not True:
        raise ValueError(f"{output_tsv}: symbolic_enforcement constraints_applied must be True: {sym_stage}")
    if sym_stage.get("taxon_consistency_state") != "applied":
        raise ValueError(
            f"{output_tsv}: symbolic_enforcement taxon_consistency_state must be 'applied': {sym_stage}"
        )
    if sym_stage.get("taxon_not_in_hierarchy_count", 0) != 0:
        raise ValueError(
            f"{output_tsv}: symbolic_enforcement taxon_not_in_hierarchy_count must be 0: {sym_stage}"
        )
    if sym_stage.get("annotations_relaxed", 0) != 0:
        raise ValueError(
            f"{output_tsv}: remove mode must have annotations_relaxed == 0, got {sym_stage}"
        )

    removed_candidates = sym_stage.get("annotations_affected")
    ec_violations = sym_stage.get("ec_violations", 0)
    violations = sym_stage.get("violations", 0)
    ec_claims_checked = sym_stage.get("ec_claims_checked", 0)
    ec_claims_not_checkable = sym_stage.get("ec_claims_not_checkable", 0)
    if type(removed_candidates) is not int or removed_candidates < 0:
        raise ValueError(f"{output_tsv}: invalid annotations_affected in {sym_stage}")
    if type(ec_violations) is not int or ec_violations < 0 or ec_violations > removed_candidates:
        raise ValueError(f"{output_tsv}: invalid ec_violations in {sym_stage}")

    batch_stats_path = output_tsv.parent / "batch_stats.tsv"
    if batch_stats_path.is_file():
        with batch_stats_path.open("r", encoding="utf-8") as fh:
            stats_rows = list(csv.DictReader(fh, delimiter="\t"))
        if len(stats_rows) != 1:
            raise ValueError(f"{batch_stats_path}: expected 1 row, got {len(stats_rows)}")
        srow = stats_rows[0]
        if srow.get("enforce_stage2") not in (None, "true"):
            raise ValueError(f"{batch_stats_path}: expected enforce_stage2='true', got {srow.get('enforce_stage2')!r}")
        if srow.get("enforce_mode") not in (None, "remove"):
            raise ValueError(f"{batch_stats_path}: expected enforce_mode='remove', got {srow.get('enforce_mode')!r}")
        if srow.get("taxon_in_hierarchy") not in (None, "true"):
            raise ValueError(
                f"{batch_stats_path}: expected taxon_in_hierarchy='true', got {srow.get('taxon_in_hierarchy')!r}"
            )
        if srow.get("taxon_consistency_state") not in (None, "applied"):
            raise ValueError(
                f"{batch_stats_path}: expected taxon_consistency_state='applied', "
                f"got {srow.get('taxon_consistency_state')!r}"
            )
        if srow.get("stage2_annotations_affected") is not None:
            if int(srow["stage2_annotations_affected"]) != removed_candidates:
                raise ValueError(
                    f"{batch_stats_path}: stage2_annotations_affected ({srow['stage2_annotations_affected']}) "
                    f"differs from ledger annotations_affected ({removed_candidates})"
                )

    return {
        "symbolic_enforcement_status": "executed",
        "taxon_consistency_state": "applied",
        "taxon_in_hierarchy": True,
        "constraints_applied": True,
        "removed_candidates": removed_candidates,
        "ec_violations": ec_violations,
        "violations": int(violations),
        "ec_claims_checked": int(ec_claims_checked),
        "ec_claims_not_checkable": int(ec_claims_not_checkable),
    }


def verify_stage1_taxon_candidates(
    expected: Path,
    prediction: Path,
    scratch: Path,
    *,
    expected_removed: int,
) -> dict[str, int]:
    """Verify bounded-memory candidate subset identity under `--enforce-consistency-mode remove`.

    Rejects any duplicate or invented candidate key, verifies that only GO/EC keys are
    removed, and requires `removed_candidates == expected_removed` and
    `prediction_rows + removed_candidates == input_candidates`.
    """
    def keys():
        with Path(prediction).open() as stream:
            reader = csv.DictReader(stream, delimiter="\t")
            required = {"protein_id", "type", "function_id"}
            if not required <= set(reader.fieldnames or []):
                raise ValueError("prediction header lacks candidate identity")
            for row in reader:
                yield U._key(row["protein_id"], row["type"], row["function_id"])

    actual = Path(scratch) / "actual.keys"
    U._sort_keys(keys(), actual)
    retained_rows = 0
    removed_rows = 0
    removed_ec_rows = 0
    removed_go_rows = 0
    try:
        exp_iter = iter(U._unique_lines(expected))
        act_iter = iter(U._unique_lines(actual))
        exp_line = next(exp_iter, None)
        act_line = next(act_iter, None)
        while exp_line is not None or act_line is not None:
            if exp_line is not None and act_line == exp_line:
                retained_rows += 1
                exp_line = next(exp_iter, None)
                act_line = next(act_iter, None)
            elif act_line is None or (exp_line is not None and exp_line < act_line):
                parsed = json.loads(exp_line)
                kind = parsed[1]
                if kind == "EC":
                    removed_ec_rows += 1
                elif kind == "GO":
                    removed_go_rows += 1
                else:
                    raise ValueError(f"unexpected candidate type removed in stage1_taxon: {kind!r}")
                removed_rows += 1
                exp_line = next(exp_iter, None)
            else:
                raise ValueError(
                    f"stage1_taxon output contains invented candidate not in full matrix input: "
                    f"{act_line.strip()}"
                )
        if removed_rows != expected_removed:
            raise ValueError(
                f"stage1_taxon removed candidate count ({removed_rows}) differs from "
                f"symbolic_enforcement annotations_affected ({expected_removed})"
            )
        return {
            "prediction_rows": retained_rows,
            "removed_candidates": removed_rows,
            "removed_ec_candidates": removed_ec_rows,
            "removed_go_candidates": removed_go_rows,
        }
    finally:
        actual.unlink(missing_ok=True)


def _run_capacity_taxon(
    part: Path,
    layout: Path,
    jar: Path,
    reference: Path,
    taxon_lineage: Path,
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
    expected_taxon = "NCBITaxon_" + str(layout_record["taxid"])
    expected, count = U.expected_candidates(claims, scratch)
    if count != publication["claims"]:
        raise ValueError("published claim count differs from observed universe")

    taxon_lineage_resolved = Path(taxon_lineage).resolve()
    if not taxon_lineage_resolved.is_file():
        raise FileNotFoundError(f"Missing taxon lineage TSV: {taxon_lineage_resolved}")

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
                taxon=expected_taxon,
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
        elif arm == "stage1_taxon":
            command += stage1_taxon_extra_flags(taxon_lineage_resolved)
            sources.append(taxon_lineage_resolved)
        else:
            raise ValueError(f"Unsupported arm: {arm!r}")

        commands_by_arm[arm] = list(command)
        save(directory / "command.json", command)
        V.run_batch(command, manifest, directory, sources=sources, timeout=timeout)
        taxon_meta = verify_taxon_outputs(
            arm,
            output,
            expected_taxon=expected_taxon,
        )
        if arm == "stage0":
            rows = U.compare_candidates(expected, output, scratch)
            removed_candidates = 0
            removed_ec_candidates = 0
            removed_go_candidates = 0
            universe_preserved = True
        else:
            cand_check = verify_stage1_taxon_candidates(
                expected,
                output,
                scratch,
                expected_removed=taxon_meta["removed_candidates"],
            )
            rows = cand_check["prediction_rows"]
            removed_candidates = cand_check["removed_candidates"]
            removed_ec_candidates = cand_check["removed_ec_candidates"]
            removed_go_candidates = cand_check["removed_go_candidates"]
            universe_preserved = removed_candidates == 0

        metrics = {}
        for line in (directory / "resources.txt").read_text().splitlines():
            key, _, value = line.strip().partition(": ")
            metrics[key] = value
        reports[arm] = {
            "prediction_rows": rows,
            "input_candidates": count,
            "removed_candidates": removed_candidates,
            "removed_ec_candidates": removed_ec_candidates,
            "removed_go_candidates": removed_go_candidates,
            "candidate_universe_preserved": universe_preserved,
            "candidate_subset_verified": True,
            "taxon_enforcement": taxon_meta,
            "resources": metrics,
            "receipt_sha256": V.sha256_file(str(directory / V.RECEIPT)),
        }

    if "stage0" in commands_by_arm and "stage1_taxon" in commands_by_arm:
        verify_arm_command_pair(
            commands_by_arm["stage0"],
            commands_by_arm["stage1_taxon"],
            taxon_lineage=str(taxon_lineage_resolved),
        )

    if V.file_hashes(inputs) != inputs:
        raise ValueError("part changed during capacity measurement")
    return {
        "schema": "gspa.price-taxon-capacity/1",
        "status": "verified",
        "assembly": asm,
        "predictor": record["predictor"],
        "matrix_rows": publication["matrix_rows"],
        "runs": reports,
        "scope": "Complete matrix; Stage 0 vs Stage 1 + Stage 2 taxon/EC consistency enforcement (remove mode).",
    }


def run_cell_capacity(
    part: Path,
    layout: Path,
    jar: Path,
    reference: Path,
    taxon_lineage: Path,
    *,
    out: Path,
    scratch: Path | None = None,
    java_heap: str = "32g",
    arms: tuple[str, ...] = ARMS,
    timeout: int = 7200,
) -> dict[str, Any]:
    if not re.fullmatch(r"[1-9][0-9]*[gGmM]", java_heap):
        raise ValueError("java_heap must be a positive size in M or G")
    if not arms or len(set(arms)) != len(arms) or not set(arms) <= set(ARMS):
        raise ValueError("arms must be distinct stage0/stage1_taxon names")
    if isinstance(timeout, bool) or not isinstance(timeout, int) or timeout <= 0:
        raise ValueError("timeout must be a positive integer number of seconds")
    scratch = Path(scratch or os.environ.get("TMPDIR") or Path.home() / "tmp").resolve()
    scratch.mkdir(parents=True, exist_ok=True)
    out.mkdir(parents=True, exist_ok=False)
    reports: dict[str, Any] = {}
    started = time.monotonic()
    try:
        with tempfile.TemporaryDirectory(prefix="gspa-taxon-cap-", dir=scratch) as work:
            result = _run_capacity_taxon(
                part,
                layout,
                jar,
                reference,
                taxon_lineage,
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
                "schema": "gspa.price-taxon-capacity/1",
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


def _verify_run_candidate_accounting(arm: str, report: dict[str, Any], expected_claims: int) -> None:
    if report.get("input_candidates") != expected_claims:
        raise ValueError(f"{arm}: input_candidates ({report.get('input_candidates')}) != {expected_claims}")
    if report.get("candidate_subset_verified") is not True:
        raise ValueError(f"{arm}: candidate_subset_verified is not True")
    pred_rows = report.get("prediction_rows")
    removed = report.get("removed_candidates", 0)
    if type(pred_rows) is not int or type(removed) is not int or pred_rows < 0 or removed < 0:
        raise ValueError(f"{arm}: invalid prediction_rows/removed_candidates in {report}")
    if arm == "stage0":
        if removed != 0 or pred_rows != expected_claims or report.get("candidate_universe_preserved") is not True:
            raise ValueError(f"stage0 candidate population differs from full matrix input: {report}")
    elif arm == "stage1_taxon":
        if pred_rows + removed != expected_claims:
            raise ValueError(
                f"stage1_taxon candidate accounting mismatch: "
                f"prediction_rows ({pred_rows}) + removed_candidates ({removed}) != {expected_claims}"
            )
        if report.get("candidate_universe_preserved") != (removed == 0):
            raise ValueError(f"stage1_taxon inconsistent candidate_universe_preserved flag: {report}")
        taxon_enf = report.get("taxon_enforcement", {})
        if taxon_enf.get("removed_candidates") != removed:
            raise ValueError(
                f"stage1_taxon removed_candidates ({removed}) != "
                f"taxon_enforcement.removed_candidates ({taxon_enf.get('removed_candidates')})"
            )
    else:
        raise ValueError(f"Unsupported arm {arm!r}")


def prepare(
    export_path: Path,
    audit_path: Path,
    jar: Path,
    reference: Path,
    taxon_lineage: Path,
    out: Path,
    *,
    allow_synthetic_subset: bool = False,
) -> dict[str, Any]:
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

    taxon_lineage_resolved = Path(taxon_lineage).resolve()
    layout = Path(export["layout"])
    layout_report = json.loads((layout / "report.json").read_text())
    tax_summary = verify_price22_taxon_lineages(
        layout_report,
        taxon_lineage_resolved,
        allow_synthetic_subset=allow_synthetic_subset,
    )

    files = [
        export_path,
        audit_path,
        jar,
        taxon_lineage_resolved,
        Path(shutil.which("java")),
        Path(sys.executable),
        Path(__file__).resolve(),
    ]
    files += [p for p in V3_DIR.glob("*.py")]
    files += V.streaming_parser_sources()
    files += [V3_DIR / "requirements-streaming.txt"]
    files += [reference / name for name in ("go.owl", "ec2go.txt", "kegg_pathways.tsv", "kegg_modules.tsv")]
    files += [p for p in (reference / "reasoner-cache").rglob("*") if p.is_file()]
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
        schema="gspa.price-taxon-production/1",
        cells=rows,
        arms=list(ARMS),
        jar=str(jar),
        layout=str(layout),
        reference=str(reference),
        taxon_lineage=str(taxon_lineage_resolved),
        taxon_lineage_sha256=tax_summary["taxon_lineage_sha256"],
        enforce_stage2=True,
        enforce_consistency_mode="remove",
        taxonomy_verification=tax_summary,
        java=str(Path(shutil.which("java")).resolve()),
        python=str(Path(sys.executable).resolve()),
        expected_claims=export["expected_claims"],
        scope=__doc__,
        model_fitting=False,
        labels_read=False,
        interpretation=(
            "stage0 is production integration without priors (--enable-priors none); "
            "stage1_taxon enables default Stage 1 priors plus Stage 2 taxon/EC consistency "
            "enforcement in remove mode (--enforce-stage2 --enforce-consistency-mode remove "
            "--taxon-lineage <lineage_tsv>)."
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
    return dict(
        cells=len(rows),
        expected_claims=plan["expected_claims"],
        arms=list(ARMS),
        taxon_lineage_sha256=plan["taxon_lineage_sha256"],
    )


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
        run_cell_capacity(
            Path(cell["part"]),
            Path(plan["layout"]),
            Path(plan["jar"]),
            Path(plan["reference"]),
            Path(plan["taxon_lineage"]),
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
        for arm, r in verified["runs"].items():
            _verify_run_candidate_accounting(arm, r, cell["claims"])
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
    """Recheck every input, command contract, stage ledger, and candidate accounting before certifying."""
    reg = json.loads((out / "registration.json").read_text())
    if reg["plan_sha256"] != V.sha256_file(str(out / "plan.json")):
        raise ValueError("production plan changed")
    plan = json.loads((out / "plan.json").read_text())
    if tuple(plan.get("arms", ())) != ARMS:
        raise ValueError(f"Expected arms {ARMS}, got {plan.get('arms')}")
    taxon_lineage_resolved = str(Path(plan["taxon_lineage"]).resolve())
    if reg["inputs"].get(taxon_lineage_resolved) != plan.get("taxon_lineage_sha256"):
        raise ValueError("taxon_lineage_sha256 in plan.json is not bound in registration.json inputs")

    layout_report = json.loads((Path(plan["layout"]) / "report.json").read_text())
    selected = completed_attempts(out, plan)
    inputs = dict(reg["inputs"])
    artifacts: dict[str, str] = {}
    reports = []
    total_removed_candidates = 0
    total_ec_violations = 0

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

        expected_taxon = "NCBITaxon_" + str(layout_report["assemblies"][cell["assembly"]]["taxid"])
        commands_by_arm: dict[str, list[str]] = {}
        for arm, report in result["runs"].items():
            _verify_run_candidate_accounting(arm, report, cell["claims"])
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
            if arm == "stage1_taxon" and spec["inputs"].get(taxon_lineage_resolved) != plan["taxon_lineage_sha256"]:
                raise ValueError("stage1_taxon receipt does not bind the registered taxon lineage TSV")
            commands_by_arm[arm] = list(spec["command"])
            if V.output_hashes(spec) != receipt["artifacts"]:
                raise ValueError("prediction or ledger changed")
            artifacts.update(receipt["artifacts"])

            output_tsv = Path(spec["outputs"][cell["assembly"]])
            taxon_meta = verify_taxon_outputs(
                arm,
                output_tsv,
                expected_taxon=expected_taxon,
            )
            if report.get("taxon_enforcement") != taxon_meta:
                raise ValueError("Stage 2 taxon enforcement ledger metadata changed since cell completion")
            if arm == "stage1_taxon":
                total_removed_candidates += taxon_meta["removed_candidates"]
                total_ec_violations += taxon_meta["ec_violations"]

            for in_path, digest in spec["inputs"].items():
                if in_path in inputs and inputs[in_path] != digest:
                    raise ValueError("inconsistent inference input hashes")
                inputs[in_path] = digest

        verify_arm_command_pair(
            commands_by_arm["stage0"],
            commands_by_arm["stage1_taxon"],
            taxon_lineage=taxon_lineage_resolved,
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
                schema="gspa.price-taxon-collection-verification/1",
                status="verified",
                cells=len(reports),
                arms=len(plan["arms"]),
                expected_claims=plan["expected_claims"],
                taxon_lineage=plan["taxon_lineage"],
                taxon_lineage_sha256=plan["taxon_lineage_sha256"],
                enforce_stage2=plan["enforce_stage2"],
                enforce_consistency_mode=plan["enforce_consistency_mode"],
                total_removed_candidates=total_removed_candidates,
                total_ec_violations=total_ec_violations,
                taxonomy_verification=plan["taxonomy_verification"],
                inputs=inputs,
                artifacts=artifacts,
                reports=reports,
                parser=V.parser_identity(),
                interpretation=(
                    "Complete Stage 1 + Stage 2 taxon/EC consistency enforcement (remove mode) "
                    "vs Stage 0 Price collection."
                ),
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
        "taxon-lineage",
        "scratch",
        "container",
        "launch-prefix",
    ):
        p.add_argument("--" + name, type=Path)
    p.add_argument("--allow-synthetic-subset", action="store_true")
    p.add_argument("--index", type=int)
    p.add_argument("--attempt")
    p.add_argument("--heap", default="24g")
    p.add_argument("--timeout", type=int, default=7200)
    p.add_argument("--returncode", type=int)
    args = p.parse_args(argv)

    if args.action == "prepare":
        if not all((args.export, args.audit, args.jar, args.reference)):
            p.error("prepare requires --export, --audit, --jar, --reference")
        lineage = (args.taxon_lineage or DEFAULT_LINEAGE_TSV).resolve()
        res = prepare(
            args.export.resolve(),
            args.audit.resolve(),
            args.jar.resolve(),
            args.reference.resolve(),
            lineage,
            args.out.resolve(),
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
