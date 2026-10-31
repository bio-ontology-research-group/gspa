#!/usr/bin/env python3
"""WP4 recurrence analysis: label-blind seal, labelled score, report tables.

Design: ``tasks/v3-wp4-recurrence-design-2026-10-29.md`` §5 (operational
definition), §6 (seal before labels), §7 (labels), §8 (T1-T9, power gate),
§11 (negative-result wording) and §13 (D1: joint recurrence = R2 - R0 inside
each relation; D3: no refitting; D4: R0 under the fenced reference is compared
with 10406 inputs and reported, not interpreted).

``seal``   reads only production ``gspa integrate --batch-manifest`` outputs
           (per-genome ``<genome>.integrated.tsv``, ``batch_stats.tsv``,
           ``verified_run.json``, ``features/<genome>/candidate_features.jsonl``),
           plan tables and label-blind sequence inputs, and writes the sealed
           table plus its sha256 in ``seal.json``.  It refuses label-file paths.
``score``  refuses a table whose hash differs from the seal, joins strict
           labels and computes T1-T9.  Unconfirmed is never negative: no false
           positive, precision, Brier or ECE is computed anywhere.
``report`` renders Markdown tables (numbers only in tables) and the §11
           negative-result sentences.

This module never imports ``match_swissprot_labels`` or ``homology_baseline``.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import random
import re
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from audit_wp4_leakage import label_pattern  # noqa: E402  (name matching only)

SEAL_SCHEMA = "gspa.wp4-recurrence-seal/1"
SCORE_SCHEMA = "gspa.wp4-recurrence-score/1"
TABLE_FILE = "sealed_table.tsv"
GAPS_FILE = "sealed_gaps.tsv"
SEAL_FILE = "seal.json"
METRICS_FILE = "metrics.json"

# Pre-registered analysis plan (§5, §8).  ``score`` refuses a seal whose
# recorded plan differs from these constants: no post-hoc change of delta.
ANALYSIS_PLAN = {
    "delta": 0.05,
    "power_gate_min_assertions": 20,
    "alpha": 0.05,
    "tolerance": 1e-12,
    "bootstrap_replicates": 2000,
    "bootstrap_seed": 20261029,
    "relations": ["PATHWAY", "UNIVERSAL", "ORGANISM"],
    "anchor": "SEQONLY",
    "feature_kind_preference": ["none", "logistic", "gbdt"],
    "learned_kinds": ["logistic", "gbdt"],
    "fragment_draft_tier": "25kb",
}
RELATIONS = tuple(ANALYSIS_PLAN["relations"])
CONTEXTS = RELATIONS + ("SEQONLY",)
LEARNED = tuple(ANALYSIS_PLAN["learned_kinds"])
KINDS = ("none",) + LEARNED
DELTA = ANALYSIS_PLAN["delta"]
TOL = ANALYSIS_PLAN["tolerance"]

# §2 table: target assembly -> 10406 split (binding design facts).
TARGET_SPLITS = {
    "GCA_000007065.1": "test", "GCA_000008765.1": "test", "GCA_000147055.1": "test",
    "GCA_000009985.1": "test", "GCA_000027305.1": "calibration",
    "GCA_000009085.1": "calibration", "GCA_000145295.1": "calibration",
    "GCA_000007005.1": "train", "GCA_000012285.1": "train", "GCA_000011205.1": "train",
    "GCA_000009965.1": "train", "GCA_000025685.1": "train", "GCA_000012525.1": "train",
    "GCA_000011585.1": "train", "GCA_000046845.1": "train", "GCA_000011125.1": "train",
    "GCA_000195915.1": "train", "GCA_000008485.1": "train",
}

ARM_REQUIRED = ("regime", "ctx", "kind", "inference_dir", "holdout")
FEATURE_MAP = {
    "cohort_available": "cohort.w.available",
    "complementary_available": "cohort.w.complementary_available",
    "joint_coverage": "cohort.w.joint_coverage",
    "complementarity": "cohort.w.complementarity",
    "complementary_unannotated": "cohort.w.complementary_unannotated",
    "support_fraction": "cohort.w.support_fraction",
    "local_available": "cohort.w.local_available",
    "local_window_truncated": "cohort.w.local_window_truncated",
    "local_window_unknown": "cohort.w.local_window_unknown",
    "family_known": "cohort.w.family_known",
}
SEAL_COLUMNS = [
    "arm_id", "block", "target", "genome_id", "split", "regime", "ctx", "kind", "variant",
    "holdout", "tier", "protein_id", "function_type", "function_id", "role", "candidate_key",
    "in_claims", "sequence_evidence", "in_integrated", "posterior_prob", "features_present",
    "features_sha256", "cohort_available", "complementary_available", "joint_coverage",
    "complementarity", "complementary_unannotated", "support_fraction", "local_available",
    "local_window_truncated", "local_window_unknown", "family_known", "model_score_status",
    "model_probability", "model_margin", "display_threshold", "above_display_threshold",
    "orthogroup_id", "c_peers_with_family",
]
GAP_COLUMNS = ["arm_id", "genome_id", "ec_number"]
NUMERIC_COMPARE = ("posterior_prob", "joint_coverage", "complementarity",
                   "complementary_unannotated", "support_fraction", "local_window_truncated",
                   "local_window_unknown", "model_probability")
EC_COMPLETE = re.compile(r"^EC:\d+\.\d+\.\d+\.n?\d+$")
NA = "NA"


class AnalysisError(ValueError):
    """Inputs cannot support a sealed, complete and label-blind analysis."""


def require(condition, message):
    if not condition:
        raise AnalysisError(message)


# ---------------------------------------------------------------- utilities
def sha256_file(path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def label_blind(path):
    """Seal-side guard: a label-shaped path is refused before it is opened."""
    path = Path(path)
    pattern = label_pattern(path.name)
    require(pattern is None, f"seal refuses a label file ({pattern}): {path}")
    return path


def read_tsv(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, fields, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t",
                                lineterminator="\n", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def save_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def fmt_value(value):
    if value is None:
        return NA
    if isinstance(value, bool):
        return "1" if value else "0"
    if isinstance(value, float):
        require(math.isfinite(value), "non-finite value in sealed table")
        return repr(value)
    return str(value)


def parse_float(text):
    if text in (None, "", NA):
        return None
    return float(text)


def parse_flag(text):
    if text in (None, "", NA):
        return None
    return text == "1"


def normalise_function(function_type, function_id, transfers=None):
    value = function_id.strip().replace(" ", "")
    if function_type == "EC":
        if not value.startswith("EC:"):
            value = "EC:" + value
        if transfers:
            value = transfers.get(value, value)
    return value


def ec_complete(function_id):
    return bool(EC_COMPLETE.match(function_id))


# -------------------------------------------------------------------- seal
def arm_targets(text):
    return [item.strip() for item in (text or "").split(",") if item.strip()]


def parse_sealed_genomes(text, targets):
    """``a,b`` or ``copy=target`` entries; plain ids map to themselves if a target."""
    result = []
    entries = arm_targets(text) if text else list(targets)
    for entry in entries:
        if "=" in entry:
            genome, target = (part.strip() for part in entry.split("=", 1))
        elif entry in targets:
            genome, target = entry, entry
        else:
            require(len(targets) == 1,
                    f"sealed genome {entry} needs copy=target with several block targets")
            genome, target = entry, targets[0]
        require(target in targets, f"sealed genome {genome} maps to non-target {target}")
        result.append((genome, target))
    require(len({g for g, _ in result}) == len(result), "duplicate sealed genome in arm")
    return result


def load_arms(arm_paths, runs, peer_tiers=None):
    arms = []
    for arm_path in arm_paths:
        arm_path = label_blind(arm_path).resolve()
        rows = read_tsv(arm_path)
        require(rows, f"empty arms table: {arm_path}")
        plan = {}
        plan_path = arm_path.parent / "plan_report.json"
        if plan_path.is_file():
            plan = json.loads(plan_path.read_text())
        for line, row in enumerate(rows, 2):
            for column in ARM_REQUIRED:
                require(row.get(column), f"{arm_path}:{line}: missing {column}")
            require(row["ctx"] in CONTEXTS, f"{arm_path}:{line}: unknown ctx {row['ctx']}")
            require(row["kind"] in KINDS, f"{arm_path}:{line}: unknown kind {row['kind']}")
            block = row.get("block") or plan.get("block")
            require(block, f"{arm_path}:{line}: no block column and no plan_report.json block")
            targets = arm_targets(row.get("target") or plan.get("target"))
            require(targets, f"{arm_path}:{line}: no target column and no plan_report.json target")
            variant = row.get("variant") or "base"
            tier = row.get("tier") or ""
            if not tier and variant.startswith("addone:") and peer_tiers:
                tier = peer_tiers.get((targets[0], variant.split(":", 1)[1]), "")
            if variant.startswith("fragment:") and not tier:
                tier = variant.split(":")[1] if variant.count(":") >= 2 else ""
            run_dir = resolve_run_dir(runs, row, block)
            genomes = parse_sealed_genomes(row.get("sealed_genomes"), targets)
            arms.append({
                "arms_tsv": str(arm_path), "arms_line": line, "block": block,
                "targets": targets, "genomes": genomes, "regime": row["regime"],
                "ctx": row["ctx"], "kind": row["kind"], "variant": variant, "tier": tier,
                "holdout": row["holdout"], "split_override": row.get("split") or "",
                "run_dir": str(run_dir), "manifest": row.get("manifest", ""),
                "model": row.get("model", "-"), "plan": plan,
                "arm_id": f"{block}/{row['regime']}/{row['ctx']}/{row['kind']}/{variant}",
            })
    ids = [arm["arm_id"] for arm in arms]
    require(len(ids) == len(set(ids)), "duplicate arm ids across arms tables")
    keys = set()
    for arm in arms:
        for genome, _ in arm["genomes"]:
            key = (genome, arm["regime"], arm["ctx"], arm["kind"], arm["variant"])
            require(key not in keys, f"genome/arm key sealed twice: {key}")
            keys.add(key)
    return arms


def resolve_run_dir(runs, row, block):
    runs = Path(runs)
    if row.get("run_dir"):
        path = Path(row["run_dir"])
        return path if path.is_absolute() else runs / path
    parts = list(Path(row["inference_dir"]).parts)
    if parts and parts[-1] == "inference":
        parts = parts[:-1]
    if block in parts:
        index = len(parts) - 1 - parts[::-1].index(block)
        return runs.joinpath(*parts[index:])
    return runs.joinpath(block, row["regime"], row["ctx"], row["kind"])


def read_integrated(path, function_types):
    path = label_blind(path)
    scores = {}
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        require(header is not None and header[:3] == ["protein_id", "type", "function_id"]
                and "posterior_prob" in header, f"invalid integrated header: {path}")
        column = header.index("posterior_prob")
        for line, row in enumerate(reader, 2):
            if not row:
                continue
            if row[1] not in function_types:
                continue
            value = float(row[column])
            require(math.isfinite(value), f"{path}:{line}: non-finite posterior")
            key = (row[0], row[1], normalise_function(row[1], row[2]))
            scores[key] = max(value, scores.get(key, -math.inf))
    return scores


def feature_value(features, name):
    value = features.get(name)
    return None if value is None else float(value)


def read_features(path, function_types, genome):
    path = label_blind(path)
    rows = {}
    with path.open() as handle:
        for line, text in enumerate(handle, 1):
            if not text.strip():
                continue
            record = json.loads(text)
            ftype = record.get("function_type")
            if ftype not in function_types:
                continue
            require(record.get("genome_id") in (None, genome),
                    f"{path}:{line}: feature row genome {record.get('genome_id')} != {genome}")
            features = record.get("features") or {}
            available = record.get("available_contexts") or {}
            values = {name: feature_value(features, column) for name, column in FEATURE_MAP.items()}
            cohort = values["cohort_available"] == 1.0
            complementary = cohort and values["complementary_available"] == 1.0
            local = values["local_available"] == 1.0
            # §5: unavailable is null, never 0 (the exporter writes 0 + indicator 0).
            for name in ("joint_coverage", "complementarity", "complementary_unannotated"):
                if not complementary:
                    values[name] = None
            if not cohort:
                values["support_fraction"] = None
            for name in ("local_window_truncated", "local_window_unknown"):
                if not local:
                    values[name] = None
            model = record.get("model") or {}
            status = model.get("score_status") if model else None
            above = model.get("above_display_threshold") if model else None
            fid = normalise_function(ftype, record["function_id"])
            key = (record["protein_id"], ftype, fid, record.get("role") or "")
            require(key not in rows, f"{path}:{line}: duplicate candidate {key}")
            rows[key] = {
                "candidate_key": record.get("candidate_key", ""),
                "sequence_evidence": available.get("sequence_evidence"),
                "features_sha256": hashlib.sha256(json.dumps(
                    features, sort_keys=True, separators=(",", ":")).encode()).hexdigest(),
                "values": values,
                "model_score_status": status,
                "model_probability": model.get("probability") if model else None,
                "model_margin": model.get("margin") if model else None,
                "display_threshold": model.get("display_threshold") if model else None,
                "above_display_threshold": above,
            }
    return rows


def read_claims(path, function_types):
    path = label_blind(path)
    claims = set()
    with Path(path).open() as handle:
        for text in handle:
            if not text.strip():
                continue
            record = json.loads(text)
            if record.get("function_type") in function_types:
                claims.add((record["protein_id"], record["function_type"],
                            normalise_function(record["function_type"], record["function_id"])))
    return claims


def verify_receipt(inference, genome, integrated_path, allow_missing):
    receipt_path = inference / "verified_run.json"
    if not receipt_path.is_file():
        require(allow_missing, f"missing verified-batch receipt: {receipt_path}")
        return {"status": "missing"}
    receipt = json.loads(label_blind(receipt_path).read_text())
    require(receipt.get("returncode") == 0, f"receipt returncode is not 0: {receipt_path}")
    artifacts = receipt.get("artifacts") or {}
    name = f"{genome}.integrated.tsv"
    matches = [value for key, value in artifacts.items() if Path(key).name == name]
    require(len(matches) == 1, f"receipt does not bind {name}: {receipt_path}")
    actual = sha256_file(integrated_path)
    require(matches[0] == actual, f"integrated table sha256 differs from receipt: {integrated_path}")
    return {"status": "verified", "receipt_sha256": sha256_file(receipt_path)}


def batch_stats_row(inference, genome):
    path = inference / "batch_stats.tsv"
    if not path.is_file():
        return None
    rows = [row for row in read_tsv(label_blind(path)) if row.get("genome_id") == genome]
    require(len(rows) == 1, f"batch_stats.tsv lacks exactly one row for {genome}: {path}")
    return {key: rows[0].get(key) for key in ("n_claims", "n_integrated")}


def local_manifest(arm):
    """The plan's local manifest (arms.tsv sits in the block root)."""
    arms_root = Path(arm["arms_tsv"]).parent
    candidates = []
    if arm["manifest"]:
        parts = list(Path(arm["manifest"]).parts)
        if arm["block"] in parts:
            index = len(parts) - 1 - parts[::-1].index(arm["block"])
            relative = Path(*parts[index + 1:])
            candidates.append(arms_root / relative.with_name("manifest.local.tsv"))
            candidates.append(arms_root / relative)
    candidates.append(arms_root / arm["regime"] / arm["ctx"] / arm["kind"] / "manifest.local.tsv")
    return next((path for path in candidates if path.is_file()), None)


def gap_ecs(arm, genome):
    manifest = local_manifest(arm)
    if manifest is None:
        return None
    rows = [row for row in read_tsv(label_blind(manifest)) if row.get("genome_id") == genome]
    if len(rows) != 1 or not rows[0].get("gaps_path"):
        return None
    path = Path(rows[0]["gaps_path"])
    if not path.is_file():
        return None
    ecs = set()
    with label_blind(path).open() as handle:
        for text in handle:
            if text.strip():
                value = json.loads(text).get("ec_number")
                if value:
                    ecs.add(normalise_function("EC", value))
    return ecs, sha256_file(path), str(path)


def read_peer_tiers(path):
    if path is None:
        return {}
    tiers = {}
    for row in read_tsv(label_blind(path)):
        target = row.get("target_accession") or row.get("target")
        peer = row.get("peer_accession") or row.get("peer")
        tier = row.get("tier")
        require(target and peer and tier in ("C", "F", "ND"),
                f"peer tier rows need target_accession, peer_accession, tier C/F/ND: {path}")
        tiers[(target, peer)] = tier
    return tiers


def family_counts(orthogroups, peer_tiers, targets):
    """Label-blind: orthogroup of each target protein, and C peers carrying it."""
    if orthogroups is None:
        return {}
    c_peers = defaultdict(set)
    for (target, peer), tier in peer_tiers.items():
        if tier == "C":
            c_peers[target].add(peer)
    wanted_genomes = set(targets) | set().union(*c_peers.values()) if c_peers else set(targets)
    membership, family_genomes = {}, defaultdict(set)
    with label_blind(orthogroups).open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            genome = row["genome_id"]
            if genome not in wanted_genomes:
                continue
            family_genomes[row["orthogroup_id"]].add(genome)
            if genome in targets:
                membership[(genome, row["protein_id"])] = row["orthogroup_id"]
    result = {}
    for (genome, protein), family in membership.items():
        count = len(family_genomes[family] & c_peers[genome]) if peer_tiers else None
        result[(genome, protein)] = (family, count)
    return result


def seal(arm_paths, runs, out, function_types=("EC",), orthogroups=None, peer_tiers_path=None,
         claims_dir=None, allow_missing_receipts=False):
    out = Path(out)
    require(not (out / SEAL_FILE).exists(), f"seal already exists: {out / SEAL_FILE}")
    peer_tiers = read_peer_tiers(peer_tiers_path)
    arms = load_arms(arm_paths, runs, peer_tiers)
    all_targets = sorted({t for arm in arms for t in arm["targets"]})
    families = family_counts(orthogroups, peer_tiers, all_targets)
    claims_cache, rows, gap_rows, arm_records = {}, [], [], []
    for arm in sorted(arms, key=lambda a: a["arm_id"]):
        run_dir = Path(arm["run_dir"])
        inference = run_dir / "inference" if (run_dir / "inference").is_dir() else run_dir
        record = {key: arm[key] for key in ("arm_id", "block", "targets", "regime", "ctx",
                                            "kind", "variant", "tier", "holdout", "run_dir",
                                            "arms_tsv", "model")}
        record["genomes"] = {}
        for genome, target in arm["genomes"]:
            integrated = inference / f"{genome}.integrated.tsv"
            features_path = run_dir / "features" / genome / "candidate_features.jsonl"
            require(integrated.is_file(), f"missing integrated table: {integrated}")
            require(features_path.is_file(), f"missing feature export: {features_path}")
            receipt = verify_receipt(inference, genome, integrated, allow_missing_receipts)
            posterior = read_integrated(integrated, function_types)
            features = read_features(features_path, function_types, genome)
            claims, claims_info = None, None
            if genome == target:
                claims, claims_info = target_claims(arm, genome, claims_dir, function_types,
                                                    claims_cache)
            split = arm["split_override"] or TARGET_SPLITS.get(target, "unknown")
            keys = set(features)
            role_free = {(p, t, f) for (p, t, f, _r) in features}
            keys.update((p, t, f, "") for (p, t, f) in posterior if (p, t, f) not in role_free)
            for key in sorted(keys):
                protein, ftype, fid, role = key
                feature = features.get(key)
                family = families.get((target, protein)) if genome == target else None
                values = feature["values"] if feature else {name: None for name in FEATURE_MAP}
                rows.append({
                    "arm_id": arm["arm_id"], "block": arm["block"], "target": target,
                    "genome_id": genome, "split": split, "regime": arm["regime"],
                    "ctx": arm["ctx"], "kind": arm["kind"], "variant": arm["variant"],
                    "holdout": arm["holdout"], "tier": arm["tier"], "protein_id": protein,
                    "function_type": ftype, "function_id": fid, "role": role,
                    "candidate_key": feature["candidate_key"] if feature else "",
                    "in_claims": fmt_value(None if claims is None else (protein, ftype, fid) in claims),
                    "sequence_evidence": fmt_value(feature["sequence_evidence"] if feature else None),
                    "in_integrated": fmt_value((protein, ftype, fid) in posterior),
                    "posterior_prob": fmt_value(posterior.get((protein, ftype, fid))),
                    "features_present": fmt_value(feature is not None),
                    "features_sha256": feature["features_sha256"] if feature else NA,
                    **{name: fmt_value(values[name]) for name in FEATURE_MAP},
                    "model_score_status": fmt_value(feature["model_score_status"] if feature else None),
                    "model_probability": fmt_value(feature["model_probability"] if feature else None),
                    "model_margin": fmt_value(feature["model_margin"] if feature else None),
                    "display_threshold": fmt_value(feature["display_threshold"] if feature else None),
                    "above_display_threshold": fmt_value(
                        feature["above_display_threshold"] if feature else None),
                    "orthogroup_id": family[0] if family else NA,
                    "c_peers_with_family": fmt_value(family[1] if family else None),
                })
            genome_record = {
                "target": target, "split": split,
                "integrated": {"path": str(integrated), "sha256": sha256_file(integrated)},
                "features": {"path": str(features_path), "sha256": sha256_file(features_path)},
                "receipt": receipt, "batch_stats": batch_stats_row(inference, genome),
                "claims": claims_info,
            }
            if arm["variant"] == "removal" or arm["variant"].startswith("fragment:"):
                gaps = gap_ecs(arm, genome)
                genome_record["gaps"] = None if gaps is None else {"path": gaps[2], "sha256": gaps[1]}
                if gaps is not None:
                    gap_rows.extend({"arm_id": arm["arm_id"], "genome_id": genome, "ec_number": ec}
                                    for ec in sorted(gaps[0]))
            record["genomes"][genome] = genome_record
        arm_records.append(record)
    rows.sort(key=lambda r: (r["arm_id"], r["genome_id"], r["protein_id"], r["function_type"],
                             r["function_id"], r["role"]))
    gap_rows.sort(key=lambda r: (r["arm_id"], r["genome_id"], r["ec_number"]))
    out.mkdir(parents=True, exist_ok=True)
    write_tsv(out / TABLE_FILE, SEAL_COLUMNS, rows)
    write_tsv(out / GAPS_FILE, GAP_COLUMNS, gap_rows)
    inputs = {str(Path(p).resolve()): sha256_file(p) for p in arm_paths}
    for path in (orthogroups, peer_tiers_path):
        if path is not None:
            inputs[str(Path(path).resolve())] = sha256_file(path)
    receipt = {
        "schema": SEAL_SCHEMA, "label_blind": True,
        "table": {"file": TABLE_FILE, "sha256": sha256_file(out / TABLE_FILE), "rows": len(rows),
                  "columns": SEAL_COLUMNS},
        "gaps": {"file": GAPS_FILE, "sha256": sha256_file(out / GAPS_FILE), "rows": len(gap_rows)},
        "function_types": list(function_types), "analysis_plan": ANALYSIS_PLAN,
        "arms": arm_records, "inputs": inputs, "analyzer_sha256": sha256_file(__file__),
        "receipts_required": not allow_missing_receipts,
    }
    save_json(out / SEAL_FILE, receipt)
    return receipt


def target_claims(arm, genome, claims_dir, function_types, cache):
    plan_claims = (arm["plan"] or {}).get("target_claims") or {}
    path = None
    if claims_dir is not None and (Path(claims_dir) / f"{genome}.claims.jsonl").is_file():
        path = Path(claims_dir) / f"{genome}.claims.jsonl"
    elif plan_claims.get("path") and Path(plan_claims["path"]).is_file() and \
            Path(plan_claims["path"]).name == f"{genome}.claims.jsonl":
        path = Path(plan_claims["path"])
    if path is None:
        return None, None
    digest = sha256_file(path)
    if plan_claims.get("sha256") and Path(plan_claims.get("path", "")).name == path.name:
        require(plan_claims["sha256"] == digest,
                f"target claims sha256 differs from plan_report: {path}")
    if (str(path), digest) not in cache:
        cache[(str(path), digest)] = read_claims(path, function_types)
    return cache[(str(path), digest)], {"path": str(path), "sha256": digest}


# --------------------------------------------------------------- scoring io
class Cand:
    """One (protein, function) candidate of one arm, aggregated over roles."""
    __slots__ = ("posterior", "in_claims", "seqev", "q", "comp", "unann", "support",
                 "trunc", "unknown", "local", "status", "prob", "above", "fsha", "numbers")

    def __init__(self):
        self.posterior = None
        self.in_claims = None
        self.seqev = None
        self.q = self.comp = self.unann = self.support = None
        self.trunc = self.unknown = None
        self.local = False
        self.status = None
        self.prob = None
        self.above = None
        self.fsha = []
        self.numbers = []


def _max(a, b):
    if a is None:
        return b
    if b is None:
        return a
    return max(a, b)


def _any(a, b):
    if a is None:
        return b
    if b is None:
        return a
    return a or b


def check_seal(seal_dir, expect_seal_sha256=None):
    seal_dir = Path(seal_dir)
    seal_path = seal_dir / SEAL_FILE
    require(seal_path.is_file(), f"missing {seal_path}")
    if expect_seal_sha256 is not None:
        require(sha256_file(seal_path) == expect_seal_sha256,
                "seal.json sha256 differs from the frozen value; refusing to score")
    receipt = json.loads(seal_path.read_text())
    require(receipt.get("schema") == SEAL_SCHEMA and receipt.get("label_blind") is True,
            "not a label-blind WP4 seal")
    for part, name in (("table", TABLE_FILE), ("gaps", GAPS_FILE)):
        path = seal_dir / name
        require(path.is_file(), f"missing sealed file: {path}")
        require(sha256_file(path) == receipt[part]["sha256"],
                f"sealed {part} sha256 differs from seal.json; refusing to score")
    require(receipt.get("analysis_plan") == ANALYSIS_PLAN,
            "sealed analysis plan differs from the pre-registered constants; refusing to score")
    return receipt


def load_sealed(seal_dir, transfers=None):
    """Return arm metadata and candidates keyed by (genome, regime, ctx, kind, variant)."""
    seal_dir = Path(seal_dir)
    arms, cands, c_peers = {}, defaultdict(dict), {}
    with (seal_dir / TABLE_FILE).open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require(reader.fieldnames == SEAL_COLUMNS, "sealed table columns differ from the schema")
        for row in reader:
            key = (row["genome_id"], row["regime"], row["ctx"], row["kind"], row["variant"])
            if key not in arms:
                arms[key] = {"block": row["block"], "target": row["target"], "split": row["split"],
                             "tier": row["tier"], "arm_id": row["arm_id"], "holdout": row["holdout"]}
            fid = normalise_function(row["function_type"], row["function_id"], transfers)
            ckey = (row["protein_id"], row["function_type"], fid)
            cand = cands[key].get(ckey)
            if cand is None:
                cand = cands[key][ckey] = Cand()
            cand.posterior = _max(cand.posterior, parse_float(row["posterior_prob"]))
            cand.in_claims = _any(cand.in_claims, parse_flag(row["in_claims"]))
            cand.seqev = _any(cand.seqev, parse_flag(row["sequence_evidence"]))
            for attr, column in (("q", "joint_coverage"), ("comp", "complementarity"),
                                 ("unann", "complementary_unannotated"),
                                 ("support", "support_fraction"),
                                 ("trunc", "local_window_truncated"),
                                 ("unknown", "local_window_unknown")):
                setattr(cand, attr, _max(getattr(cand, attr), parse_float(row[column])))
            cand.local = cand.local or parse_float(row["local_available"]) == 1.0
            if row["c_peers_with_family"] != NA and row["genome_id"] == row["target"]:
                c_peers[(row["genome_id"], row["protein_id"])] = int(row["c_peers_with_family"])
            status = row["model_score_status"]
            if status == "scored":
                cand.status = "scored"
            elif status == "unsupported" and cand.status is None:
                cand.status = "unsupported"
            prob = parse_float(row["model_probability"]) if status == "scored" else None
            cand.prob = _max(cand.prob, prob)
            cand.above = _any(cand.above, parse_flag(row["above_display_threshold"]))
            cand.fsha.append(row["role"] + ":" + row["features_sha256"])
            cand.numbers.append(tuple(parse_float(row[c]) for c in NUMERIC_COMPARE))
    # Arms with zero candidate rows are still arms.
    receipt = json.loads((seal_dir / SEAL_FILE).read_text())
    for arm in receipt["arms"]:
        for genome, info in arm["genomes"].items():
            key = (genome, arm["regime"], arm["ctx"], arm["kind"], arm["variant"])
            arms.setdefault(key, {"block": arm["block"], "target": info["target"],
                                  "split": info["split"], "tier": arm["tier"],
                                  "arm_id": arm["arm_id"], "holdout": arm["holdout"]})
            cands.setdefault(key, {})
    gaps = defaultdict(set)
    for row in read_tsv(seal_dir / GAPS_FILE):
        gaps[(row["arm_id"], row["genome_id"])].add(row["ec_number"])
    gap_available = set()
    for arm in receipt["arms"]:
        for genome, info in arm["genomes"].items():
            if info.get("gaps"):
                gap_available.add((arm["arm_id"], genome))
    return arms, cands, gaps, gap_available, c_peers


def read_transfers(path):
    if path is None:
        return None
    transfers = {}
    for row in read_tsv(path):
        old = normalise_function("EC", row["old"])
        new = normalise_function("EC", row["new"])
        require(old not in transfers or transfers[old] == new, f"ambiguous EC transfer {old}")
        transfers[old] = new
    return transfers


def load_labels(path, genomes, transfers=None):
    """Strict experimental labels (§7): complete EC with ECO:0000269; experimental GO."""
    ec, go = defaultdict(set), defaultdict(set)
    rows = read_tsv(path)
    require(rows and {"genome_id", "protein_id", "function_id", "evidence", "tier"} <= set(rows[0]),
            f"label table lacks the cohort label columns: {path}")
    for row in rows:
        if row["genome_id"] not in genomes or row["tier"] != "experimental":
            continue
        fid = row["function_id"].strip()
        if fid.startswith("EC:"):
            fid = normalise_function("EC", fid, transfers)
            if row["evidence"] == "ECO:0000269" and ec_complete(fid):
                ec[(row["genome_id"], row["protein_id"])].add(fid)
        elif fid.startswith("GO:"):
            # The cohort label tier is authoritative (EXP/IDA/IPI/IMP/IGI/IEP and the
            # high-throughput H* codes); the curated-inferential tier is already excluded.
            go[(row["genome_id"], row["protein_id"])].add(fid)
    return ec, go


def go_ancestors(obo_path):
    if obo_path is None:
        return None
    parents, term = defaultdict(set), None
    with Path(obo_path).open() as handle:
        for line in handle:
            line = line.strip()
            if line == "[Term]":
                term = None
            elif line.startswith("id: GO:"):
                term = line[4:]
            elif term and line.startswith("is_a: "):
                parents[term].add(line[6:].split()[0])
            elif term and line.startswith("relationship: part_of "):
                parents[term].add(line.split()[2])
    cache = {}

    def closure(node):
        if node not in cache:
            seen, stack = set(), [node]
            while stack:
                for parent in parents.get(stack.pop(), ()):
                    if parent not in seen:
                        seen.add(parent)
                        stack.append(parent)
            cache[node] = seen
        return cache[node]
    return closure


# --------------------------------------------------------------- statistics
def sign_test_one_sided(toward, away):
    """Exact one-sided binomial P(X >= toward | n = toward + away, 1/2)."""
    n = toward + away
    if n == 0:
        return None
    return sum(math.comb(n, k) for k in range(toward, n + 1)) / 2 ** n


def mcnemar_exact(gains, harms):
    n = gains + harms
    if n == 0:
        return None
    tail = sum(math.comb(n, k) for k in range(0, min(gains, harms) + 1)) / 2 ** n
    return min(1.0, 2 * tail)


def holm(pvalues):
    """Holm step-down adjustment for a {name: p or None} map."""
    items = sorted((p, name) for name, p in pvalues.items() if p is not None)
    adjusted, running, m = {}, 0.0, len(items)
    for index, (p, name) in enumerate(items):
        running = max(running, min(1.0, (m - index) * p))
        adjusted[name] = running
    return {name: adjusted.get(name) for name in pvalues}


def block_bootstrap(values_by_block, replicates, seed):
    blocks = sorted(values_by_block)
    if not blocks or not any(values_by_block.values()):
        return None
    rng = random.Random(seed)
    means = []
    for _ in range(replicates):
        sample = [value for block in (rng.choice(blocks) for _ in blocks)
                  for value in values_by_block[block]]
        if sample:
            means.append(sum(sample) / len(sample))
    means.sort()
    if not means:
        return None
    low = means[int(math.floor(0.025 * (len(means) - 1)))]
    high = means[int(math.ceil(0.975 * (len(means) - 1)))]
    return {"low": low, "high": high, "replicates": len(means), "blocks": len(blocks)}


def normalized_entropy(shares):
    if shares is None or len(shares) < 2:
        return None
    entropy = -sum(s * math.log(s) for s in shares.values() if s > 0)
    return entropy / math.log(len(shares))


def quantiles(values):
    if not values:
        return None
    ordered = sorted(values)

    def q(fraction):
        position = fraction * (len(ordered) - 1)
        low, high = math.floor(position), math.ceil(position)
        return ordered[low] + (ordered[high] - ordered[low]) * (position - low)
    return {"min": ordered[0], "q25": q(0.25), "median": q(0.5), "q75": q(0.75),
            "max": ordered[-1], "n": len(ordered)}


def mean(values):
    return sum(values) / len(values) if values else None


# -------------------------------------------------------------- the model
class Analysis:
    def __init__(self, arms, cands, gaps, gap_available, labels_ec, labels_go=None,
                 go_closure=None):
        self.arms, self.cands = arms, cands
        self.gaps, self.gap_available = gaps, gap_available
        self.labels_ec, self.labels_go = labels_ec, labels_go or {}
        self.go_closure = go_closure
        self.targets = sorted({key[0] for key, meta in arms.items()
                               if key[4] == "base" and key[0] == meta["target"]})
        self.candidate_sets, self.candidate_set_mismatch = self._candidate_sets()

    # -- arm access
    def arm(self, genome, regime, ctx, kind, variant="base"):
        key = (genome, regime, ctx, kind, variant)
        return self.cands[key] if key in self.arms else None

    def feature_arm(self, genome, regime, ctx, variant="base"):
        for kind in ANALYSIS_PLAN["feature_kind_preference"]:
            arm = self.arm(genome, regime, ctx, kind, variant)
            if arm is not None:
                return arm, (genome, regime, ctx, kind, variant)
        return None, None

    def variants(self, genome, prefix):
        return sorted({key[4] for key in self.arms if key[0] == genome and key[4].startswith(prefix)})

    def assertions(self, genome):
        result = []
        for (g, protein), functions in sorted(self.labels_ec.items()):
            if g == genome:
                result.extend((genome, protein, f) for f in sorted(functions))
        return result

    def _candidate_sets(self):
        """E(p): complete ECs in p's raw claims (claims file; else sequence evidence)."""
        sets, mismatch = {}, 0
        for genome in self.targets:
            per_arm = []
            for key, cands in self.cands.items():
                if key[0] != genome or key[4] != "base":
                    continue
                current = defaultdict(set)
                for (protein, ftype, fid), cand in cands.items():
                    if ftype != "EC" or not ec_complete(fid):
                        continue
                    member = cand.in_claims if cand.in_claims is not None else cand.seqev
                    if member:
                        current[protein].add(fid)
                if cands:
                    per_arm.append((key, dict(current)))
            reference = per_arm[0][1] if per_arm else {}
            for key, current in per_arm[1:]:
                if key[2] != "SEQONLY" and current != reference:
                    mismatch += 1
            sets[genome] = reference
        return sets, mismatch

    def E(self, genome, protein):
        return self.candidate_sets.get(genome, {}).get(protein, set())

    # -- §5 quantities
    @staticmethod
    def shares(arm, protein, E):
        if arm is None:
            return None
        available = {}
        for e in E:
            cand = arm.get((protein, "EC", e))
            if cand is not None and cand.q is not None:
                available[e] = cand.q
        total = sum(available.values())
        if total <= 0:
            return None
        return {e: q / total for e, q in available.items()}

    @staticmethod
    def rank_pessimistic(values, focus):
        if focus not in values or values[focus] is None:
            return None
        return 1 + sum(1 for e, v in values.items()
                       if e != focus and v is not None and v >= values[focus])

    def classify(self, genome, protein, F, ctx, g0_variant="base", gf_variant="base",
                 gf_regime="R2"):
        """§5 classification of one labelled assertion in one relation."""
        E = self.E(genome, protein)
        g0, _ = self.feature_arm(genome, "R0", ctx, g0_variant)
        gf, gf_key = self.feature_arm(genome, gf_regime, ctx, gf_variant)
        result = {"E": len(E), "block": self.arms[gf_key]["block"] if gf_key else None}
        if g0 is None or gf is None:
            result["category"] = "arm_missing"
            return result
        if F not in E:
            result["category"] = "no_candidate"
            return result
        if len(E) == 1:
            result["category"] = "no_alternative"
            return result
        s0, sf = self.shares(g0, protein, E), self.shares(gf, protein, E)
        if s0 is None or sf is None or F not in s0 or F not in sf:
            result["category"] = "f_unavailable"
            return result
        delta_s = sf[F] - s0[F]
        result["delta_s"] = delta_s
        result["category"] = ("toward" if delta_s > DELTA else
                              "away" if delta_s < -DELTA else "unchanged")
        c0, cf = g0.get((protein, "EC", F)), gf.get((protein, "EC", F))
        if c0.comp is None or cf.comp is None:
            result["complementary"] = None
        else:
            result["complementary"] = (cf.comp - c0.comp) > DELTA
        support0 = {e: g0[(protein, "EC", e)].support for e in E if (protein, "EC", e) in g0}
        supportf = {e: gf[(protein, "EC", e)].support for e in E if (protein, "EC", e) in gf}
        result["rank_g0"] = self.rank_pessimistic(support0, F)
        result["rank_gfull"] = self.rank_pessimistic(supportf, F)
        return result

    def protein_entropy_change(self, genome, protein, ctx, gf_variant="base"):
        E = self.E(genome, protein)
        g0, _ = self.feature_arm(genome, "R0", ctx)
        gf, _ = self.feature_arm(genome, "R2", ctx, gf_variant)
        if g0 is None or gf is None:
            return "arm_missing", None
        h0 = normalized_entropy(self.shares(g0, protein, E))
        hf = normalized_entropy(self.shares(gf, protein, E))
        if h0 is None or hf is None:
            return "h_undefined", None
        change = hf - h0
        return ("narrowed" if change < -DELTA else "broadened" if change > DELTA
                else "unchanged"), change


T2_CATEGORIES = ("toward", "away", "unchanged", "f_unavailable", "no_alternative",
                 "no_candidate", "arm_missing")


def partition(counts, categories, denominator):
    total = sum(counts.get(c, 0) for c in categories)
    require(total == denominator, f"partition {categories} sums to {total}, not {denominator}")
    return {"denominator": denominator, **{c: counts.get(c, 0) for c in categories},
            "partition_ok": True}


# ---------------------------------------------------------------- T1-T9
def t1_ledger(an):
    rows = []
    for genome in an.targets:
        assertions = an.assertions(genome)
        proteins = {p for _, p, _ in assertions}
        base = {"target": genome, "split": an_split(an, genome),
                "assertions": len(assertions), "proteins": len(proteins)}
        ledger = defaultdict(int)
        family_1 = family_3 = 0
        family_known = True
        for _, protein, F in assertions:
            E = an.E(genome, protein)
            ledger["no_candidate" if F not in E else "f_candidate_only" if len(E) == 1
                   else "multi_candidate"] += 1
            count = c_peer_count(an, genome, protein)
            if count is None:
                family_known = False
            else:
                family_1 += count >= 1
                family_3 += count >= 3
        row = dict(base, **partition(ledger, ("no_candidate", "f_candidate_only",
                                              "multi_candidate"), len(assertions)))
        row["family_in_ge1_C_peers"] = family_1 if family_known else None
        row["family_in_ge3_C_peers"] = family_3 if family_known else None
        for ctx in RELATIONS:
            sub = defaultdict(int)
            for _, protein, F in assertions:
                E = an.E(genome, protein)
                if F not in E or len(E) == 1:
                    continue
                g0, _ = an.feature_arm(genome, "R0", ctx)
                gf, _ = an.feature_arm(genome, "R2", ctx)
                s0 = an.shares(g0, protein, E)
                sf = an.shares(gf, protein, E) if gf is not None else None
                a0 = s0 is not None and F in s0
                af = sf is not None and F in sf
                if g0 is None or gf is None:
                    sub["arm_missing"] += 1
                elif not a0:
                    sub["f_unavailable_g0"] += 1
                elif not af:
                    sub["f_unavailable_gfull_only"] += 1
                else:
                    sub["f_available_both"] += 1
            row[f"multi_{ctx}"] = partition(sub, ("f_unavailable_g0", "f_unavailable_gfull_only",
                                                  "f_available_both", "arm_missing"),
                                            ledger["multi_candidate"])
        rows.append(row)
    total = {"target": "all", "assertions": sum(r["assertions"] for r in rows),
             "proteins": sum(r["proteins"] for r in rows)}
    for category in ("no_candidate", "f_candidate_only", "multi_candidate"):
        total[category] = sum(r[category] for r in rows)
    return {"rows": rows, "total": total}


def an_split(an, genome):
    for key, meta in an.arms.items():
        if key[0] == genome:
            return meta["split"]
    return "unknown"


def c_peer_count(an, genome, protein):
    # Stored in the sealed table; recover it from the raw TSV lazily.
    return an.c_peers.get((genome, protein)) if hasattr(an, "c_peers") else None


def t2_narrowing(an):
    result, pvalues = {}, {}
    for ctx in RELATIONS:
        pooled, by_split = defaultdict(int), defaultdict(lambda: defaultdict(int))
        deltas_by_block, ranks = defaultdict(list), defaultdict(int)
        complementary = defaultdict(int)
        denominator = 0
        per_assertion = []
        for genome in an.targets:
            split = an_split(an, genome)
            for _, protein, F in an.assertions(genome):
                denominator += 1
                state = an.classify(genome, protein, F, ctx)
                category = state["category"]
                pooled[category] += 1
                by_split[split][category] += 1
                per_assertion.append({"target": genome, "protein": protein, "function": F,
                                      "category": category, "delta_s": state.get("delta_s")})
                if "delta_s" in state:
                    deltas_by_block[state["block"]].append(state["delta_s"])
                    r0, rf = state["rank_g0"], state["rank_gfull"]
                    if r0 is None or rf is None:
                        ranks["rank_unavailable"] += 1
                    else:
                        ranks["rank_improved" if rf < r0 else "rank_worsened" if rf > r0
                              else "rank_same"] += 1
                if category == "toward":
                    flag = state.get("complementary")
                    complementary["complementarity_unavailable" if flag is None else
                                  "complementary" if flag else "not_complementary"] += 1
        table = partition(pooled, T2_CATEGORIES, denominator)
        testable_n = table["toward"] + table["away"] + table["unchanged"]
        p = sign_test_one_sided(table["toward"], table["away"])
        pvalues[ctx] = p
        all_deltas = [d for values in deltas_by_block.values() for d in values]
        result[ctx] = {
            "pooled": table,
            "by_split": {split: partition(counts, T2_CATEGORIES, sum(counts.values()))
                         for split, counts in sorted(by_split.items())},
            "power_gate": {"assertions_available_both": testable_n,
                           "minimum": ANALYSIS_PLAN["power_gate_min_assertions"],
                           "testable": testable_n >= ANALYSIS_PLAN["power_gate_min_assertions"]},
            "sign_test_one_sided_p": p,
            "mean_delta_s": mean(all_deltas),
            "mean_delta_s_block_bootstrap": block_bootstrap(
                deltas_by_block, ANALYSIS_PLAN["bootstrap_replicates"],
                ANALYSIS_PLAN["bootstrap_seed"]),
            "secondary_rank": partition(ranks, ("rank_improved", "rank_worsened", "rank_same",
                                                "rank_unavailable"), testable_n),
            "complementary_of_toward": partition(complementary, (
                "complementary", "not_complementary", "complementarity_unavailable"),
                table["toward"]),
            "assertions": per_assertion,
        }
    adjusted = holm(pvalues)
    for ctx in RELATIONS:
        entry = result[ctx]
        entry["holm_p"] = adjusted[ctx]
        if not entry["power_gate"]["testable"]:
            entry["verdict"] = "not_testable"
        elif (entry["pooled"]["toward"] > entry["pooled"]["away"] and adjusted[ctx] is not None
              and adjusted[ctx] < ANALYSIS_PLAN["alpha"]):
            entry["verdict"] = "narrowing"
        else:
            entry["verdict"] = "no_narrowing"
    return result


def relation_pairs(t2):
    """§11: do relations differ on paired assertions (toward status, exact McNemar)?"""
    status = {}
    for ctx in RELATIONS:
        for item in t2[ctx]["assertions"]:
            if item["category"] in ("toward", "away", "unchanged"):
                status.setdefault((item["target"], item["protein"], item["function"]), {})[ctx] = \
                    item["category"] == "toward"
    result = {}
    for i, a in enumerate(RELATIONS):
        for b in RELATIONS[i + 1:]:
            paired = [v for v in status.values() if a in v and b in v]
            only_a = sum(1 for v in paired if v[a] and not v[b])
            only_b = sum(1 for v in paired if v[b] and not v[a])
            result[f"{a}|{b}"] = {"paired": len(paired), "toward_only_first": only_a,
                                  "toward_only_second": only_b,
                                  "mcnemar_p": mcnemar_exact(only_a, only_b)}
    differ = any(v["mcnemar_p"] is not None and v["mcnemar_p"] < ANALYSIS_PLAN["alpha"]
                 for v in result.values())
    return {"pairs": result, "relations_differ": differ}


def shuffle_family(ctx):
    return ("S1:", "S2:") if ctx == "PATHWAY" else ("S1:",)


def t3_t4_nulls(an, t2):
    """Real vs each S1/S2 seed on the same genomes (the null is computed, never assumed 0)."""
    t3, t4 = {}, {}
    for ctx in RELATIONS:
        seeds = sorted({key[4] for key in an.arms if key[2] == ctx and key[1] == "R2"
                        and key[4].startswith(shuffle_family(ctx))})
        null3, null4 = [], []
        for seed in seeds:
            genomes = sorted({key[0] for key in an.arms if key[4] == seed and key[2] == ctx
                              and key[0] in an.targets})
            real = defaultdict(int)
            null = defaultdict(int)
            for genome in genomes:
                for _, protein, F in an.assertions(genome):
                    state = an.classify(genome, protein, F, ctx)
                    shuffled = an.classify(genome, protein, F, ctx, gf_variant=seed)
                    for counter, item in ((real, state), (null, shuffled)):
                        if item["category"] == "toward":
                            counter["toward"] += 1
                            counter["complementary"] += bool(item.get("complementary"))
                        if item["category"] == "away":
                            counter["away"] += 1
            null3.append({"seed": seed, "genomes": genomes, "real_toward": real["toward"],
                          "null_toward": null["toward"],
                          "real_complementary": real["complementary"],
                          "null_complementary": null["complementary"]})
            real_lf, null_lf = label_free_counts(an, genomes, ctx), \
                label_free_counts(an, genomes, ctx, seed)
            null4.append({"seed": seed, "genomes": genomes,
                          "real_narrowed": real_lf["narrowed"], "null_narrowed": null_lf["narrowed"],
                          "multi_candidate_proteins": real_lf["denominator"]})
        toward = t2[ctx]["pooled"]["toward"]
        complementary = t2[ctx]["complementary_of_toward"]["complementary"]
        t3[ctx] = {
            "toward": toward, "complementary": complementary,
            "not_complementary": t2[ctx]["complementary_of_toward"]["not_complementary"],
            "complementarity_unavailable":
                t2[ctx]["complementary_of_toward"]["complementarity_unavailable"],
            "null_by_seed": null3,
            "null_status": "computed" if null3 else "null_missing",
            "null_all_zero": bool(null3) and all(
                item["null_toward"] == 0 for item in null3),
            "toward_exceeds_all_seeds": (all(i["real_toward"] > i["null_toward"] for i in null3)
                                         if null3 else None),
            "complementary_exceeds_all_seeds": (
                all(i["real_complementary"] > i["null_complementary"] for i in null3)
                if null3 else None),
            "power_gate_testable": t2[ctx]["power_gate"]["testable"],
        }
        lf = label_free_counts(an, an.targets, ctx)
        t4[ctx] = dict(lf, null_by_seed=null4,
                       null_status="computed" if null4 else "null_missing",
                       null_all_zero=bool(null4) and all(i["null_narrowed"] == 0 for i in null4),
                       narrowed_exceeds_all_seeds=(
                           all(i["real_narrowed"] > i["null_narrowed"] for i in null4)
                           if null4 else None))
    return t3, t4


def label_free_counts(an, genomes, ctx, variant="base"):
    counts, changes, denominator = defaultdict(int), [], 0
    for genome in genomes:
        for protein in sorted(an.candidate_sets.get(genome, {})):
            if len(an.E(genome, protein)) < 2:
                continue
            denominator += 1
            category, change = an.protein_entropy_change(genome, protein, ctx, variant)
            counts[category] += 1
            if change is not None:
                changes.append(change)
    table = partition(counts, ("narrowed", "broadened", "unchanged", "h_undefined",
                               "arm_missing"), denominator)
    defined = table["narrowed"] + table["broadened"] + table["unchanged"]
    table["share_narrowed_of_defined"] = table["narrowed"] / defined if defined else None
    table["delta_h"] = quantiles(changes)
    return table


def t5_add_one(an):
    result = {}
    for ctx in RELATIONS:
        by_tier = defaultdict(lambda: {"counts": defaultdict(int), "dF": [], "dC": []})
        for genome in an.targets:
            g0, _ = an.feature_arm(genome, "R0", ctx)
            for variant in an.variants(genome, "addone:"):
                arm, key = an.feature_arm(genome, "ADD1", ctx, variant)
                if arm is None:
                    arm, key = next(((an.cands[k], k) for k in sorted(an.arms)
                                     if k[0] == genome and k[2] == ctx and k[4] == variant),
                                    (None, None))
                if arm is None:
                    continue
                tier = an.arms[key]["tier"] or "unknown"
                bucket = by_tier[tier]
                for _, protein, F in an.assertions(genome):
                    E = an.E(genome, protein)
                    if F not in E:
                        bucket["counts"]["no_candidate"] += 1
                        continue
                    if len(E) == 1:
                        bucket["counts"]["no_alternative"] += 1
                        continue
                    c0, cj = (g0 or {}).get((protein, "EC", F)), arm.get((protein, "EC", F))
                    if c0 is None or cj is None or c0.q is None or cj.q is None:
                        bucket["counts"]["f_unavailable"] += 1
                        continue
                    d_f = cj.q - c0.q
                    competitors = []
                    for e in E - {F}:
                        a, b = (g0 or {}).get((protein, "EC", e)), arm.get((protein, "EC", e))
                        if a is not None and b is not None and a.q is not None and b.q is not None:
                            competitors.append(b.q - a.q)
                    if not competitors:
                        bucket["counts"]["competitor_unavailable"] += 1
                        continue
                    d_c = max(competitors)  # pessimistic: the largest competitor gain
                    bucket["dF"].append(d_f)
                    bucket["dC"].append(d_c)
                    bucket["counts"]["f_gain_above_competitor" if d_f - d_c > TOL
                                     else "f_gain_not_above_competitor"] += 1
        result[ctx] = {}
        for tier, bucket in sorted(by_tier.items()):
            denominator = sum(bucket["counts"].values())
            entry = partition(bucket["counts"], (
                "f_gain_above_competitor", "f_gain_not_above_competitor", "f_unavailable",
                "competitor_unavailable", "no_alternative", "no_candidate"), denominator)
            entry["mean_delta_q_f"] = mean(bucket["dF"])
            entry["mean_delta_q_best_competitor"] = mean(bucket["dC"])
            result[ctx][tier] = entry
    return result


def protein_ec(arm):
    by_protein = defaultdict(dict)
    for (protein, ftype, fid), cand in arm.items():
        if ftype == "EC":
            by_protein[protein][fid] = cand
    return by_protein


def recovery_states(an, genome, arm, ftype="EC"):
    """Per labelled assertion: recovered / below / abstained / absent, top-1, RR."""
    labels = an.labels_ec if ftype == "EC" else an.labels_go
    by_protein = defaultdict(dict)
    for (protein, t, fid), cand in arm.items():
        if t == ftype:
            by_protein[protein][fid] = cand
    states = {}
    for (g, protein), functions in sorted(labels.items()):
        if g != genome:
            continue
        for F in sorted(functions):
            cand = by_protein.get(protein, {}).get(F)
            if cand is None:
                state = "candidate_absent"
            elif cand.status != "scored" or cand.prob is None:
                state = "abstained" if cand.status == "unsupported" else "not_scored"
            elif cand.above:
                state = "recovered"
            else:
                state = "below_threshold"
            top1, rr = False, 0.0
            if cand is not None and cand.prob is not None:
                others = [c.prob for f, c in by_protein[protein].items()
                          if f != F and c.prob is not None]
                rank = 1 + sum(1 for value in others if value >= cand.prob)  # ties pessimistic
                top1, rr = rank == 1, 1.0 / rank
            states[(protein, F)] = {"state": state, "top1": top1, "rr": rr}
    return states


def call_counts(an, genome, arm, ftype="EC"):
    labels = an.labels_ec if ftype == "EC" else an.labels_go
    counts = defaultdict(int)
    for (protein, t, fid), cand in arm.items():
        if t != ftype:
            continue
        counts[{"scored": "scored", "unsupported": "abstained_unchanged_unsupported"}.get(
            cand.status, "no_model_score")] += 1
        if cand.above:
            truth = labels.get((genome, protein), set())
            if ftype == "GO" and an.go_closure is not None:
                truth = set(truth) | {a for term in truth for a in an.go_closure(term)}
            # Unconfirmed is never negative: no false-positive category exists.
            counts["calls_confirmed" if fid in truth else "calls_unconfirmed"] += 1
    return dict(counts)


RECOVERY_CATEGORIES = ("recovered", "below_threshold", "abstained", "not_scored",
                       "candidate_absent")


def t6_learned(an, ftype="EC", splits=("test",)):
    per_arm, paired = [], []
    for genome in an.targets:
        if an_split(an, genome) not in splits:
            continue
        n = sum(len(v) for (g, _), v in (an.labels_ec if ftype == "EC" else an.labels_go).items()
                if g == genome)
        for kind in LEARNED:
            states = {}
            for regime, ctx in [("R0", "SEQONLY")] + [(r, c) for c in RELATIONS for r in ("R0", "R2")]:
                arm = an.arm(genome, regime, ctx, kind)
                if arm is None:
                    continue
                s = recovery_states(an, genome, arm, ftype)
                states[(regime, ctx)] = s
                counts = defaultdict(int)
                for item in s.values():
                    counts[item["state"]] += 1
                row = {"target": genome, "regime": regime, "ctx": ctx, "kind": kind,
                       **partition(counts, RECOVERY_CATEGORIES, n),
                       "strict_top1": sum(i["top1"] for i in s.values()),
                       "mrr": (sum(i["rr"] for i in s.values()) / n) if n else None,
                       "candidates": call_counts(an, genome, arm, ftype)}
                per_arm.append(row)
            for ctx in RELATIONS:
                for label, base in (("R2_vs_R0", ("R0", ctx)), ("R2_vs_SEQONLY", ("R0", "SEQONLY"))):
                    if ("R2", ctx) not in states or base not in states:
                        continue
                    after, before = states[("R2", ctx)], states[base]
                    gains = sum(1 for k in after if after[k]["state"] == "recovered"
                                and before[k]["state"] != "recovered")
                    harms = sum(1 for k in after if before[k]["state"] == "recovered"
                                and after[k]["state"] != "recovered")
                    paired.append({"target": genome, "ctx": ctx, "kind": kind, "contrast": label,
                                   "assertions": n, "gains": gains, "harms": harms,
                                   "mcnemar_p": mcnemar_exact(gains, harms)})
    pooled = defaultdict(lambda: {"gains": 0, "harms": 0, "assertions": 0})
    for row in paired:
        bucket = pooled[(row["contrast"], row["ctx"], row["kind"])]
        for key in ("gains", "harms", "assertions"):
            bucket[key] += row[key]
    pooled_rows = [{"contrast": c, "ctx": x, "kind": k, **v,
                    "mcnemar_p": mcnemar_exact(v["gains"], v["harms"])}
                   for (c, x, k), v in sorted(pooled.items())]
    abstain_all = defaultdict(lambda: True)
    for row in per_arm:
        if row["regime"] == "R2":
            scored = row["candidates"].get("scored", 0)
            abstain_all[(row["ctx"], row["kind"])] &= scored == 0
    return {"function_type": ftype, "splits": list(splits), "per_arm": per_arm,
            "paired": paired, "paired_pooled": pooled_rows,
            "r2_all_abstained": {f"{c}/{k}": v for (c, k), v in sorted(abstain_all.items())},
            "note": "underpowered by design; reported per genome; supports no positive claim alone"}


CONTROLS = (
    # name, regime, variant, baseline (regime, variant), ctx filter, expectation
    ("R0-SC", "R0-SC", "base", ("R0", "base"), None, "zero"),
    ("R1", "R1", "base", ("R0", "base"), None, "zero"),
    ("dup_i", None, "dup:i", ("R2", "base"), None, "zero"),
    ("dup_v", None, "dup:v", ("R2", "base"), None, "zero"),
    ("seqonly_parity", "R2", "base", ("R0", "base"), "SEQONLY", "zero"),
    ("dup_ii", None, "dup:ii", ("R2", "base"), None, "inflation"),
    ("dup_iii", None, "dup:iii", ("R2", "base"), None, "inflation"),
    ("dup_iv", None, "dup:iv", ("R2", "base"), None, "inflation"),
    ("d4_fenced_vs_10406", "R0", "ref10406", ("R0", "base"), None, "report_only"),
)


def compare_arms(after, before):
    keys = set(after) | set(before)
    changed, max_abs, increases, deltas_q, deltas_s = 0, 0.0, 0, [], []
    for key in keys:
        a, b = after.get(key), before.get(key)
        if a is None or b is None:
            changed += 1
            continue
        moved = sorted(a.fsha) != sorted(b.fsha)
        for x, y in zip(sorted(a.numbers, key=repr), sorted(b.numbers, key=repr)):
            for u, v in zip(x, y):
                if (u is None) != (v is None):
                    moved = True
                elif u is not None:
                    diff = abs(u - v)
                    max_abs = max(max_abs, diff)
                    moved = moved or diff > TOL
        if len(a.numbers) != len(b.numbers):
            moved = True
        changed += moved
        if a.q is not None and b.q is not None:
            deltas_q.append(a.q - b.q)
            increases += (a.q - b.q) > TOL
        if a.support is not None and b.support is not None:
            deltas_s.append(a.support - b.support)
    return {"candidates": len(keys), "changed": changed, "max_abs_delta": max_abs,
            "mean_delta_joint_coverage": mean(deltas_q),
            "mean_delta_support_fraction": mean(deltas_s), "joint_coverage_increased": increases}


def t7_controls(an):
    rows = []
    for genome in an.targets:
        for name, regime, variant, (b_regime, b_variant), ctx_filter, expect in CONTROLS:
            for key in sorted(an.arms):
                g, r, ctx, kind, v = key
                if g != genome or v != variant or (regime is not None and r != regime):
                    continue
                if ctx_filter is not None and ctx != ctx_filter:
                    continue
                base = an.arm(genome, b_regime, ctx, kind, b_variant)
                if base is None:
                    continue
                stats = compare_arms(an.cands[key], base)
                defect = expect == "zero" and (stats["changed"] > 0 or stats["max_abs_delta"] > TOL)
                rows.append({"control": name, "target": genome, "ctx": ctx, "kind": kind,
                             "expectation": expect, **stats, "defect": defect})
        # Copy layout: intact target inside removal/fragment batches equals the plain batch.
        for key in sorted(an.arms):
            g, r, ctx, kind, v = key
            if g == genome and (v == "removal" or v.startswith("fragment:")):
                base = an.arm(genome, r, ctx, kind, "base")
                if base is None:
                    continue
                stats = compare_arms(an.cands[key], base)
                rows.append({"control": f"copy_layout:{v}", "target": genome, "ctx": ctx,
                             "kind": kind, "regime": r, "expectation": "zero", **stats,
                             "defect": stats["changed"] > 0 or stats["max_abs_delta"] > TOL})
    return {"rows": rows, "defects": sum(r["defect"] for r in rows),
            "candidate_set_mismatches": an.candidate_set_mismatch}


def t8_removal(an, removals):
    if removals is None:
        return {"status": "removals_manifest_missing"}
    rows = defaultdict(lambda: defaultdict(int))
    unann = defaultdict(list)
    copy_arms = defaultdict(list)
    for key in an.arms:
        if key[4] == "removal":
            copy_arms[key[0]].append(key)
    denominators = defaultdict(int)
    for item in removals:
        copy = item["copy_genome"]
        F = normalise_function("EC", item["function_id"])
        for key in sorted(copy_arms.get(copy, [])):
            _, regime, ctx, kind, _ = key
            bucket_key = (item["class"], regime, ctx, kind)
            denominators[bucket_key] += 1
            arm_id = an.arms[key]["arm_id"]
            counts = rows[bucket_key]
            if (arm_id, copy) not in an.gap_available:
                counts["gaps_unavailable"] += 1
            elif F in an.gaps.get((arm_id, copy), set()):
                counts["f_in_observed_gaps"] += 1
            else:
                counts["f_not_in_observed_gaps"] += 1
            others = [c for (p, t, f), c in an.cands[key].items()
                      if t == "EC" and f == F and p != item["removed_protein"]]
            values = [c.unann for c in others if c.unann is not None]
            if values:
                unann[bucket_key].append(max(values))
            if any(c.above for c in others):
                counts["f_reassigned_above_threshold_unconfirmed"] += 1
    result = []
    for bucket_key in sorted(denominators):
        counts = rows[bucket_key]
        entry = partition(counts, ("f_in_observed_gaps", "f_not_in_observed_gaps",
                                   "gaps_unavailable"), denominators[bucket_key])
        observed = entry["f_in_observed_gaps"] + entry["f_not_in_observed_gaps"]
        entry.update({"class": bucket_key[0], "regime": bucket_key[1], "ctx": bucket_key[2],
                      "kind": bucket_key[3],
                      "gap_rate": entry["f_in_observed_gaps"] / observed if observed else None,
                      "mean_complementary_unannotated": mean(unann[bucket_key]),
                      "f_reassigned_above_threshold_unconfirmed":
                          counts["f_reassigned_above_threshold_unconfirmed"]})
        result.append(entry)
    comparisons = []
    index = {(r["class"], r["regime"], r["ctx"], r["kind"]): r for r in result}
    for (cls, regime, ctx, kind), row in sorted(index.items()):
        if cls != "labelled":
            continue
        base = index.get(("random", regime, ctx, kind))
        comparisons.append({
            "regime": regime, "ctx": ctx, "kind": kind,
            "labelled_rate": row["gap_rate"],
            "random_base_rate": base["gap_rate"] if base else None,
            "null_status": "computed" if base and base["gap_rate"] is not None else "null_missing",
            "labelled_above_base": (row["gap_rate"] > base["gap_rate"]
                                    if base and row["gap_rate"] is not None
                                    and base["gap_rate"] is not None else None)})
    harm = []
    for (cls, regime, ctx, kind), row in sorted(index.items()):
        if kind in LEARNED and ctx in RELATIONS and cls == "labelled":
            anchor = index.get((cls, regime, "SEQONLY", kind))
            if anchor is not None:
                harm.append({"regime": regime, "ctx": ctx, "kind": kind,
                             "reassigned": row["f_reassigned_above_threshold_unconfirmed"],
                             "reassigned_seqonly": anchor["f_reassigned_above_threshold_unconfirmed"],
                             "harm": row["f_reassigned_above_threshold_unconfirmed"] >
                             anchor["f_reassigned_above_threshold_unconfirmed"]})
    return {"rows": result, "labelled_vs_random": comparisons, "reassignment_vs_seqonly": harm}


def t9_fragmentation(an):
    rows = []
    variants = sorted({key[4] for key in an.arms if key[4].startswith("fragment:")})
    for variant in variants:
        parts = variant.split(":")
        tier, seed = (parts[1], parts[2]) if len(parts) >= 3 else (variant, "")
        for ctx in RELATIONS:
            counts, denominator = defaultdict(int), 0
            plain_toward = frag_toward = plain_away = frag_away = 0
            recovery = defaultdict(lambda: {"plain": 0, "fragment": 0})
            for genome in an.targets:
                copies = sorted({key[0] for key, meta in an.arms.items()
                                 if key[4] == variant and meta["target"] == genome
                                 and key[0] != genome})
                if not copies:
                    continue
                require(len(copies) == 1, f"several fragment copies of {genome} in {variant}")
                copy = copies[0]
                plain, _ = an.feature_arm(genome, "R0", ctx)
                frag, _ = an.feature_arm(copy, "R0", ctx, variant)
                for _, protein, F in an.assertions(genome):
                    denominator += 1
                    p = (plain or {}).get((protein, "EC", F))
                    f = (frag or {}).get((protein, "EC", F))
                    if p is None:
                        counts["no_candidate"] += 1
                    elif f is None:
                        counts["sequence_lost"] += 1
                    elif not p.local or p.trunc is None or p.unknown is None:
                        counts["window_unobserved_in_plain"] += 1
                    elif p.trunc == 0 and p.unknown == 0 and f.local and \
                            ((f.trunc or 0) > 0 or (f.unknown or 0) > 0):
                        counts["window_full_to_truncated_or_unknown"] += 1
                    else:
                        counts["window_unchanged"] += 1
                    plain_state = an.classify(genome, protein, F, ctx)["category"]
                    frag_state = classify_copy(an, copy, genome, protein, F, ctx, variant)
                    plain_toward += plain_state == "toward"
                    plain_away += plain_state == "away"
                    frag_toward += frag_state == "toward"
                    frag_away += frag_state == "away"
                for kind in LEARNED:
                    plain_arm = an.arm(genome, "R2", ctx, kind)
                    frag_arm = an.arm(copy, "R2", ctx, kind, variant)
                    if plain_arm is None or frag_arm is None:
                        continue
                    recovery[kind]["plain"] += sum(
                        s["state"] == "recovered" for s in recovery_states(an, genome, plain_arm).values())
                    recovery[kind]["fragment"] += sum(
                        s["state"] == "recovered"
                        for s in recovery_states_copy(an, genome, frag_arm).values())
            if denominator == 0:
                continue
            entry = partition(counts, ("sequence_lost", "window_full_to_truncated_or_unknown",
                                       "window_unchanged", "window_unobserved_in_plain",
                                       "no_candidate"), denominator)
            entry.update({"variant": variant, "tier": tier, "seed": seed, "ctx": ctx,
                          "plain_toward": plain_toward, "plain_away": plain_away,
                          "fragment_toward": frag_toward, "fragment_away": frag_away,
                          "delta_narrowing": (frag_toward - frag_away) - (plain_toward - plain_away),
                          "recovery": {k: dict(v, delta=v["fragment"] - v["plain"])
                                       for k, v in sorted(recovery.items())}})
            rows.append(entry)
    return {"rows": rows}


def classify_copy(an, copy, genome, protein, F, ctx, variant):
    """R2 - R0 inside fragment batches; E(p) is the intact target's claim set."""
    E = an.E(genome, protein)
    g0, _ = an.feature_arm(copy, "R0", ctx, variant)
    gf, _ = an.feature_arm(copy, "R2", ctx, variant)
    if g0 is None or gf is None:
        return "arm_missing"
    if F not in E:
        return "no_candidate"
    if len(E) == 1:
        return "no_alternative"
    s0, sf = an.shares(g0, protein, E), an.shares(gf, protein, E)
    if s0 is None or sf is None or F not in s0 or F not in sf:
        return "f_unavailable"
    d = sf[F] - s0[F]
    return "toward" if d > DELTA else "away" if d < -DELTA else "unchanged"


def recovery_states_copy(an, genome, arm):
    return recovery_states(an, genome, arm)


def read_removals(path):
    if path is None:
        return None
    rows = read_tsv(path)
    require(rows and {"copy_genome", "target_genome", "removed_protein", "class",
                      "function_id"} <= set(rows[0]),
            "removals manifest needs copy_genome, target_genome, removed_protein, class, function_id")
    for row in rows:
        require(row["class"] in ("labelled", "random"), f"unknown removal class {row['class']}")
    return rows


def score(labels, seal_dir, out, expect_seal_sha256=None, removals=None, ec_transfers=None,
          go_obo=None, label_policy=None):
    receipt = check_seal(seal_dir, expect_seal_sha256)
    transfers = read_transfers(ec_transfers)
    arms, cands, gaps, gap_available, c_peers = load_sealed(seal_dir, transfers)
    targets = {key[0] for key, meta in arms.items() if key[0] == meta["target"]}
    labels_ec, labels_go = load_labels(labels, targets, transfers)
    an = Analysis(arms, cands, gaps, gap_available, labels_ec, labels_go, go_ancestors(go_obo))
    an.c_peers = c_peers
    t2 = t2_narrowing(an)
    t3, t4 = t3_t4_nulls(an, t2)
    metrics = {
        "schema": SCORE_SCHEMA, "analysis_plan": ANALYSIS_PLAN,
        "seal": {"seal_sha256": sha256_file(Path(seal_dir) / SEAL_FILE),
                 "table_sha256": receipt["table"]["sha256"],
                 "gaps_sha256": receipt["gaps"]["sha256"]},
        "labels": {"path": str(Path(labels).resolve()), "sha256": sha256_file(labels),
                   "policy_sha256": sha256_file(label_policy) if label_policy else None,
                   "strict_ec_assertions": sum(len(v) for v in labels_ec.values()),
                   "experimental_go_assertions": sum(len(v) for v in labels_go.values())},
        "targets": sorted(an.targets),
        "T1": t1_ledger(an), "T2": t2, "relations": relation_pairs(t2),
        "T3": t3, "T4": t4, "T5": t5_add_one(an),
        "T6": t6_learned(an),
        "T6_GO": t6_learned(an, "GO") if "GO" in receipt["function_types"] else None,
        "T7": t7_controls(an), "T8": t8_removal(an, read_removals(removals)),
        "T9": t9_fragmentation(an),
        "unconfirmed_policy": "unconfirmed calls are counted, never false positives; "
                              "no precision, Brier or ECE",
    }
    out = Path(out)
    save_json(out / METRICS_FILE, metrics)
    save_json(out / "score_receipt.json", {
        "schema": SCORE_SCHEMA, "metrics_sha256": sha256_file(out / METRICS_FILE),
        "inputs": {str(Path(p).resolve()): sha256_file(p)
                   for p in [labels, Path(seal_dir) / SEAL_FILE, removals, ec_transfers, go_obo,
                             label_policy] if p is not None},
        "analyzer_sha256": sha256_file(__file__)})
    return metrics


# ------------------------------------------------------------------ report
NEGATIVE = {
    "T2_no_narrowing": "No narrowing toward labelled function from close relatives in the {ctx} relation.",
    "T2_not_testable": "Narrowing in the {ctx} relation is not testable: fewer labelled multi-candidate "
                       "assertions with the labelled function available in both genome sets than the "
                       "pre-registered power gate requires. This is neither a positive nor a negative result.",
    "T2T3_background": "In the {ctx} relation the signal cannot be told apart from density-matched background.",
    "T3_none": "No evidence that complementary pieces combine in the {ctx} relation.",
    "T4_shuffle": "Label-free narrowing in the {ctx} relation is no more frequent than under shuffles.",
    "T6_no_gain": "Close-relative recurrence does not improve the current {kind} models in the {ctx} "
                  "relation ({contrast}).",
    "T6_abstain": "The current {kind} models cannot use close-relative recurrence in the {ctx} relation: "
                  "they abstain on R2 features.",
    "relations_equal": "Relations do not differ on paired assertions: no relation advantage.",
    "T7_defect": "Engine or holdout defect: a target feature or probability changed beyond tolerance "
                 "in an invariance control. Stop and fix before interpreting anything.",
    "T8_no_localization": "Recurrence does not localize genuinely missing functions ({regime}, {ctx}).",
    "T8_harm": "Learned reassignment of the removed function is more frequent than under SEQONLY "
               "({regime}, {ctx}, {kind}): harm.",
    "T9_draft": "The benefit does not hold on draft assemblies ({ctx}, {tier} tier).",
}
SCOPE = ("Claims are limited to these genomes, this engine configuration and DIAMOND/Swiss-Prot "
         "candidates. Context concentration is not evidence that a function is correct. The "
         "joint-recurrence effect is the R2 - R0 contrast inside each relation (D1); models are "
         "frozen (D3); R0 is rerun under the fenced reference and its difference from the 10406 "
         "inputs is reported, not interpreted (D4).")


def cell(value, kind="g"):
    if value is None:
        return "NA"
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        return f"{value:.3g}" if kind == "p" else f"{value:.4g}"
    return str(value).replace("|", "\\|")


def md_table(header, rows):
    lines = ["| " + " | ".join(header) + " |",
             "|" + "|".join("---" if i == 0 else "---:" for i in range(len(header))) + "|"]
    lines.extend("| " + " | ".join(cell(v) for v in row) + " |" for row in rows)
    return "\n".join(lines)


def negative_results(metrics):
    sentences = []
    t2, t3, t4 = metrics["T2"], metrics["T3"], metrics["T4"]
    for ctx in RELATIONS:
        if t2[ctx]["verdict"] == "not_testable":
            sentences.append(NEGATIVE["T2_not_testable"].format(ctx=ctx))
        else:
            if t2[ctx]["verdict"] == "no_narrowing":
                sentences.append(NEGATIVE["T2_no_narrowing"].format(ctx=ctx))
            if t3[ctx]["toward_exceeds_all_seeds"] is False or \
                    t3[ctx]["complementary_exceeds_all_seeds"] is False:
                sentences.append(NEGATIVE["T2T3_background"].format(ctx=ctx))
            if t3[ctx]["complementary"] == 0:
                sentences.append(NEGATIVE["T3_none"].format(ctx=ctx))
        if t4[ctx]["narrowed_exceeds_all_seeds"] is False:
            sentences.append(NEGATIVE["T4_shuffle"].format(ctx=ctx))
    for row in metrics["T6"]["paired_pooled"]:
        if row["gains"] - row["harms"] <= 0:
            sentences.append(NEGATIVE["T6_no_gain"].format(kind=row["kind"], ctx=row["ctx"],
                                                           contrast=row["contrast"].replace("_", " ")))
    for key, value in metrics["T6"]["r2_all_abstained"].items():
        if value:
            ctx, kind = key.split("/")
            sentences.append(NEGATIVE["T6_abstain"].format(ctx=ctx, kind=kind))
    if not metrics["relations"]["relations_differ"]:
        sentences.append(NEGATIVE["relations_equal"])
    if metrics["T7"]["defects"]:
        sentences.insert(0, NEGATIVE["T7_defect"])
    t8 = metrics["T8"]
    for row in t8.get("labelled_vs_random", []):
        if row["labelled_above_base"] is False:
            sentences.append(NEGATIVE["T8_no_localization"].format(regime=row["regime"], ctx=row["ctx"]))
    for row in t8.get("reassignment_vs_seqonly", []):
        if row["harm"]:
            sentences.append(NEGATIVE["T8_harm"].format(**{k: row[k] for k in ("regime", "ctx", "kind")}))
    for row in metrics["T9"]["rows"]:
        if row["tier"] == ANALYSIS_PLAN["fragment_draft_tier"] and \
                row["fragment_toward"] - row["fragment_away"] <= 0 < row["plain_toward"] - row["plain_away"]:
            sentences.append(NEGATIVE["T9_draft"].format(ctx=row["ctx"], tier=row["tier"]))
    return list(dict.fromkeys(sentences))


def render(metrics):
    parts = ["# WP4 recurrence tables (generated)", "",
             f"Seal sha256 `{metrics['seal']['seal_sha256']}`, table sha256 "
             f"`{metrics['seal']['table_sha256']}`.", "", SCOPE, ""]
    t1 = metrics["T1"]
    parts += ["## T1 ledger", "", md_table(
        ["Target", "Split", "Assertions", "Proteins", "No candidate", "F candidate only",
         "Multi-candidate", "Family in >=1 C peer", "Family in >=3 C peers"],
        [[r["target"], r["split"], r["assertions"], r["proteins"], r["no_candidate"],
          r["f_candidate_only"], r["multi_candidate"], r["family_in_ge1_C_peers"],
          r["family_in_ge3_C_peers"]] for r in t1["rows"]]
        + [["all", "", t1["total"]["assertions"], t1["total"]["proteins"],
            t1["total"]["no_candidate"], t1["total"]["f_candidate_only"],
            t1["total"]["multi_candidate"], None, None]]), ""]
    parts += ["### T1 multi-candidate assertions by relation", "", md_table(
        ["Target", "Relation", "Multi-candidate", "F unavailable G0", "F unavailable Gfull only",
         "F available both", "Arm missing"],
        [[r["target"], ctx, r[f"multi_{ctx}"]["denominator"], r[f"multi_{ctx}"]["f_unavailable_g0"],
          r[f"multi_{ctx}"]["f_unavailable_gfull_only"], r[f"multi_{ctx}"]["f_available_both"],
          r[f"multi_{ctx}"]["arm_missing"]] for r in t1["rows"] for ctx in RELATIONS]), ""]
    t2 = metrics["T2"]
    parts += ["## T2 narrowing (R2 - R0 within each relation)", "", md_table(
        ["Relation", "Assertions", "Toward", "Away", "Unchanged", "F unavailable",
         "No alternative", "No candidate", "Arm missing", "Sign test p (one-sided)", "Holm p",
         "Mean delta s", "Bootstrap 95% low", "Bootstrap 95% high", "Power gate", "Verdict"],
        [[ctx, t2[ctx]["pooled"]["denominator"]] + [t2[ctx]["pooled"][c] for c in T2_CATEGORIES]
         + [t2[ctx]["sign_test_one_sided_p"], t2[ctx]["holm_p"], t2[ctx]["mean_delta_s"],
            (t2[ctx]["mean_delta_s_block_bootstrap"] or {}).get("low"),
            (t2[ctx]["mean_delta_s_block_bootstrap"] or {}).get("high"),
            "testable" if t2[ctx]["power_gate"]["testable"] else "not testable",
            t2[ctx]["verdict"].replace("_", " ")] for ctx in RELATIONS]), ""]
    parts += ["### T2 by split", "", md_table(
        ["Relation", "Split", "Assertions"] + [c.replace("_", " ") for c in T2_CATEGORIES],
        [[ctx, split, v["denominator"]] + [v[c] for c in T2_CATEGORIES]
         for ctx in RELATIONS for split, v in t2[ctx]["by_split"].items()]), ""]
    parts += ["### T2 secondary: rank of F by support fraction (ties pessimistic)", "", md_table(
        ["Relation", "F available both", "Rank improved", "Rank worsened", "Rank same",
         "Rank unavailable"],
        [[ctx, t2[ctx]["secondary_rank"]["denominator"], t2[ctx]["secondary_rank"]["rank_improved"],
          t2[ctx]["secondary_rank"]["rank_worsened"], t2[ctx]["secondary_rank"]["rank_same"],
          t2[ctx]["secondary_rank"]["rank_unavailable"]] for ctx in RELATIONS]), ""]
    pairs = metrics["relations"]["pairs"]
    parts += ["### Relation comparison on paired assertions", "", md_table(
        ["Relations", "Paired", "Toward only first", "Toward only second", "Exact McNemar p"],
        [[k.replace("|", " vs "), v["paired"], v["toward_only_first"], v["toward_only_second"],
          v["mcnemar_p"]] for k, v in pairs.items()]), ""]
    t3 = metrics["T3"]
    parts += ["## T3 complementarity (toward events)", "", md_table(
        ["Relation", "Toward", "Complementary", "Not complementary", "Complementarity unavailable",
         "Null", "Toward above every seed", "Complementary above every seed"],
        [[ctx, t3[ctx]["toward"], t3[ctx]["complementary"], t3[ctx]["not_complementary"],
          t3[ctx]["complementarity_unavailable"], t3[ctx]["null_status"].replace("_", " "),
          t3[ctx]["toward_exceeds_all_seeds"], t3[ctx]["complementary_exceeds_all_seeds"]]
         for ctx in RELATIONS]), ""]
    parts += ["### T3 real vs shuffle seeds (same genomes)", "", md_table(
        ["Relation", "Seed", "Genomes", "Real toward", "Null toward", "Real complementary",
         "Null complementary"],
        [[ctx, s["seed"], len(s["genomes"]), s["real_toward"], s["null_toward"],
          s["real_complementary"], s["null_complementary"]]
         for ctx in RELATIONS for s in t3[ctx]["null_by_seed"]]), ""]
    t4 = metrics["T4"]
    parts += ["## T4 label-free concentration (all multi-candidate target proteins)", "", md_table(
        ["Relation", "Proteins", "Narrowed", "Broadened", "Unchanged", "H undefined",
         "Arm missing", "Share narrowed", "Delta H median", "Delta H q25", "Delta H q75",
         "Narrowed above every seed"],
        [[ctx, t4[ctx]["denominator"], t4[ctx]["narrowed"], t4[ctx]["broadened"],
          t4[ctx]["unchanged"], t4[ctx]["h_undefined"], t4[ctx]["arm_missing"],
          t4[ctx]["share_narrowed_of_defined"], (t4[ctx]["delta_h"] or {}).get("median"),
          (t4[ctx]["delta_h"] or {}).get("q25"), (t4[ctx]["delta_h"] or {}).get("q75"),
          t4[ctx]["narrowed_exceeds_all_seeds"]] for ctx in RELATIONS]), ""]
    parts += ["### T4 real vs shuffle seeds", "", md_table(
        ["Relation", "Seed", "Genomes", "Multi-candidate proteins", "Real narrowed", "Null narrowed"],
        [[ctx, s["seed"], len(s["genomes"]), s["multi_candidate_proteins"], s["real_narrowed"],
          s["null_narrowed"]] for ctx in RELATIONS for s in t4[ctx]["null_by_seed"]]), ""]
    t5 = metrics["T5"]
    parts += ["## T5 add-one (assertion, peer) pairs", "", md_table(
        ["Relation", "Tier", "Pairs", "F gain above best competitor", "Not above (ties incl.)",
         "F unavailable", "Competitor unavailable", "No alternative", "No candidate",
         "Mean delta Q(F)", "Mean delta Q(best competitor)"],
        [[ctx, tier, v["denominator"], v["f_gain_above_competitor"], v["f_gain_not_above_competitor"],
          v["f_unavailable"], v["competitor_unavailable"], v["no_alternative"], v["no_candidate"],
          v["mean_delta_q_f"], v["mean_delta_q_best_competitor"]]
         for ctx in RELATIONS for tier, v in t5[ctx].items()]), ""]
    for key, title in (("T6", "T6 learned models, test targets (EC, primary)"),
                       ("T6_GO", "T6 learned models, test targets (GO, secondary)")):
        t6 = metrics[key]
        if t6 is None:
            parts += [f"## {title}", "", md_table(["Status"], [["GO not sealed"]]), ""]
            continue
        parts += [f"## {title}", "", md_table(
            ["Target", "Regime", "Context", "Kind", "Assertions", "Recovered", "Below threshold",
             "Abstained", "Not scored", "Candidate absent", "Strict top-1", "MRR",
             "Candidates scored", "Candidates abstained", "Calls confirmed", "Calls unconfirmed"],
            [[r["target"], r["regime"], r["ctx"], r["kind"], r["denominator"], r["recovered"],
              r["below_threshold"], r["abstained"], r["not_scored"], r["candidate_absent"],
              r["strict_top1"], r["mrr"], r["candidates"].get("scored", 0),
              r["candidates"].get("abstained_unchanged_unsupported", 0),
              r["candidates"].get("calls_confirmed", 0), r["candidates"].get("calls_unconfirmed", 0)]
             for r in t6["per_arm"]]), ""]
        parts += [f"### {key} paired gains and harms (pooled over test targets)", "", md_table(
            ["Contrast", "Context", "Kind", "Assertions", "Gains", "Harms", "Exact McNemar p"],
            [[r["contrast"].replace("_", " "), r["ctx"], r["kind"], r["assertions"], r["gains"],
              r["harms"], r["mcnemar_p"]] for r in t6["paired_pooled"]]), ""]
    t7 = metrics["T7"]
    parts += ["## T7 controls", "", md_table(
        ["Control", "Target", "Context", "Kind", "Expectation", "Candidates", "Changed",
         "Max abs delta", "Mean delta joint coverage", "Joint coverage increased", "Defect"],
        [[r["control"], r["target"], r["ctx"], r["kind"], r["expectation"].replace("_", " "),
          r["candidates"], r["changed"], r["max_abs_delta"], r["mean_delta_joint_coverage"],
          r["joint_coverage_increased"], r["defect"]] for r in t7["rows"]]), ""]
    t8 = metrics["T8"]
    if "rows" in t8:
        parts += ["## T8 true removal", "", md_table(
            ["Class", "Regime", "Context", "Kind", "Removals", "F in observed gaps",
             "F not in gaps", "Gaps unavailable", "Gap rate", "Mean complementary unannotated",
             "F reassigned above threshold (unconfirmed)"],
            [[r["class"], r["regime"], r["ctx"], r["kind"], r["denominator"], r["f_in_observed_gaps"],
              r["f_not_in_observed_gaps"], r["gaps_unavailable"], r["gap_rate"],
              r["mean_complementary_unannotated"], r["f_reassigned_above_threshold_unconfirmed"]]
             for r in t8["rows"]]), ""]
        parts += ["### T8 labelled rate vs random base rate", "", md_table(
            ["Regime", "Context", "Kind", "Labelled rate", "Random base rate", "Null",
             "Labelled above base"],
            [[r["regime"], r["ctx"], r["kind"], r["labelled_rate"], r["random_base_rate"],
              r["null_status"].replace("_", " "), r["labelled_above_base"]]
             for r in t8["labelled_vs_random"]]), ""]
    else:
        parts += ["## T8 true removal", "", md_table(["Status"], [[t8["status"].replace("_", " ")]]), ""]
    t9 = metrics["T9"]
    parts += ["## T9 fragmentation", "", md_table(
        ["Relation", "Tier", "Seed", "Assertions", "Sequence lost", "Window FULL to TRUNCATED/UNKNOWN",
         "Window unchanged", "Window unobserved in plain", "No candidate", "Plain toward",
         "Fragment toward", "Delta narrowing"],
        [[r["ctx"], r["tier"], r["seed"], r["denominator"], r["sequence_lost"],
          r["window_full_to_truncated_or_unknown"], r["window_unchanged"],
          r["window_unobserved_in_plain"], r["no_candidate"], r["plain_toward"],
          r["fragment_toward"], r["delta_narrowing"]] for r in t9["rows"]]), ""]
    parts += ["## Power gate", "", md_table(
        ["Relation", "Labelled multi-candidate assertions with F available in G0 and Gfull",
         "Minimum", "Status"],
        [[ctx, t2[ctx]["power_gate"]["assertions_available_both"], t2[ctx]["power_gate"]["minimum"],
          "testable" if t2[ctx]["power_gate"]["testable"] else "not testable"]
         for ctx in RELATIONS]), ""]
    return "\n".join(parts) + "\n"


def report(scores_dir, out):
    scores_dir, out = Path(scores_dir), Path(out)
    metrics_path = scores_dir / METRICS_FILE
    receipt_path = scores_dir / "score_receipt.json"
    metrics = json.loads(metrics_path.read_text())
    if receipt_path.is_file():
        require(json.loads(receipt_path.read_text()).get("metrics_sha256") == sha256_file(metrics_path),
                "metrics.json differs from score_receipt.json")
    out.mkdir(parents=True, exist_ok=True)
    (out / "v3-wp4-recurrence-tables.md").write_text(render(metrics))
    sentences = negative_results(metrics)
    lines = ["# WP4 negative results (§11), selected from the scored tables", ""]
    lines += [f"- {s}" for s in sentences] or ["- None of the §11 conditions is met."]
    lines += ["", "## All §11 templates (ready to use)", ""]
    lines += [f"- `{key}`: {text}" for key, text in NEGATIVE.items()]
    (out / "v3-wp4-negative-results.md").write_text("\n".join(lines) + "\n")
    return sentences


# --------------------------------------------------------------------- cli
def build_parser():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    steps = parser.add_subparsers(dest="step", required=True)
    s = steps.add_parser("seal", help="write the label-blind table and its sha256")
    s.add_argument("--arms", nargs="+", type=Path, required=True)
    s.add_argument("--runs", type=Path, required=True)
    s.add_argument("--out", type=Path, required=True)
    s.add_argument("--function-types", default="EC")
    s.add_argument("--orthogroups", type=Path)
    s.add_argument("--peer-tiers", type=Path)
    s.add_argument("--claims-dir", type=Path)
    s.add_argument("--allow-missing-receipts", action="store_true")
    c = steps.add_parser("score", help="join labels to a verified seal and compute T1-T9")
    c.add_argument("--labels", type=Path, required=True)
    c.add_argument("--seal", type=Path, required=True)
    c.add_argument("--out", type=Path, required=True)
    c.add_argument("--expect-seal-sha256")
    c.add_argument("--removals", type=Path)
    c.add_argument("--ec-transfers", type=Path)
    c.add_argument("--go-obo", type=Path)
    c.add_argument("--label-policy", type=Path)
    r = steps.add_parser("report", help="render the Markdown tables and §11 sentences")
    r.add_argument("--scores", type=Path, required=True)
    r.add_argument("--out", type=Path, required=True)
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    try:
        if args.step == "seal":
            types = tuple(t.strip() for t in args.function_types.split(",") if t.strip())
            require(set(types) <= {"EC", "GO"} and types, "--function-types must be EC and/or GO")
            result = seal(args.arms, args.runs, args.out, types, args.orthogroups,
                          args.peer_tiers, args.claims_dir, args.allow_missing_receipts)
            print(json.dumps({"table_sha256": result["table"]["sha256"],
                              "rows": result["table"]["rows"]}, indent=2))
        elif args.step == "score":
            metrics = score(args.labels, args.seal, args.out, args.expect_seal_sha256,
                            args.removals, args.ec_transfers, args.go_obo, args.label_policy)
            print(json.dumps({ctx: metrics["T2"][ctx]["verdict"] for ctx in RELATIONS}, indent=2))
        else:
            for sentence in report(args.scores, args.out):
                print(sentence)
        return 0
    except (AnalysisError, OSError, KeyError, ValueError, json.JSONDecodeError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
