#!/usr/bin/env python3
"""Model-genome stratum evaluation glue for the positive-unlabelled (PU) learning plan (row 4, §4.3, P5).

Evaluates frozen PN and PU model artifacts on the three model genomes:
  - E. coli K-12 MG1655 (GCF_000005845.2)
  - B. subtilis 168 (GCF_000009045.1)
  - M. tuberculosis H37Rv (GCF_000195955.2)

Stratum classification: "held-out, development-examined; never trained on".

Reuses ``prepare_training_data.GoView`` for ``is_a``/``part_of`` GO ancestor propagation of positives
and ``pu_learning.pu_evaluate`` for scoring, candidate metrics, negative specificity, pairwise contrasts,
and null-layer metrics. Negatives are exact assertions and never propagated to ancestors. Any candidate
row matching both a positive and an accepted negative is recorded as a conflict and excluded.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
V3 = HERE.parent
REPO = V3.parent.parent
sys.path.insert(0, str(V3))
sys.path.insert(0, str(REPO / "scripts"))
import prepare_training_data as PTD  # noqa: E402
import train_context_model as T  # noqa: E402
from pu_learning import pu_evaluate as PE  # noqa: E402
from pu_learning import pu_trainer as PT  # noqa: E402

SCHEMA = "gspa.pu-model-stratum-evaluation/1"
STRATUM_LABEL = "held-out, development-examined; never trained on"


# ----------------------------------------------------------------------------- normalisation & labels

def normalise_function_id(
    value: str | None,
    function_type: str | None = None,
    go: PTD.GoView | None = None,
) -> str:
    """Normalise function IDs for matching:
    - Strip 'EC:' prefix on EC IDs so 'EC:x.y.z.w' and 'x.y.z.w' map to 'x.y.z.w'.
    - Normalise GO IDs ('GO_xxxx' -> 'GO:xxxx') and resolve alt_ids via GoView.primary if provided.
    """
    if value is None:
        return ""
    v = str(value).strip().replace("_", ":")
    if not v:
        return ""
    ft = str(function_type or "").strip().upper()
    if v.upper().startswith("EC:") or ft == "EC":
        if v.upper().startswith("EC:"):
            v = v[3:].strip()
        return v
    if v.upper().startswith("GO:") or ft == "GO":
        if not v.upper().startswith("GO:"):
            v = "GO:" + v
        else:
            v = "GO:" + v[3:]
        if go is not None and not go.obsolete(v):
            v = go.primary(v)
        return v
    if v[0].isdigit() and "." in v:
        return v
    return v


def load_model_positives(
    positives_path: Path | str,
    go: PTD.GoView | None = None,
    target_set: str = "model",
) -> dict:
    """Load model-set positives from positives.jsonl and propagate GO terms to is_a/part_of ancestors via GoView.

    Returns dict with:
      - ``propagated_keys``: set of (genome_id, protein_id, norm_function_id) after GO ancestor closure
      - ``exact_keys``: set of (genome_id, protein_id, norm_function_id) before ancestor closure
      - ``provenance_by_key``: mapping (genome_id, protein_id, norm_function_id) -> list of provenance records
      - ``counts``: summary counts
    """
    p = Path(positives_path)
    if not p.is_file():
        raise FileNotFoundError(f"positives file not found: {p}")

    exact_keys: set[tuple[str, str, str]] = set()
    propagated_keys: set[tuple[str, str, str]] = set()
    provenance_by_key: dict[tuple[str, str, str], list] = defaultdict(list)
    counts = Counter()

    with p.open() as fh:
        for line in fh:
            if not line.strip():
                continue
            rec = json.loads(line)
            counts["total_records"] += 1
            if "set" in rec and str(rec["set"]).strip() != target_set:
                counts["skipped_other_set"] += 1
                continue
            counts["model_records"] += 1
            gid = str(rec["genome_id"]).strip()
            pid = str(rec["protein_id"]).strip()
            raw_fid = str(rec["function_id"]).strip().replace("_", ":")

            if raw_fid.upper().startswith("GO:"):
                fid = "GO:" + raw_fid[3:]
                if go is not None:
                    if go.obsolete(fid):
                        counts["dropped_obsolete_go"] += 1
                        continue
                    primary = go.primary(fid)
                    closure = set(go.closure(primary)) | ({primary} - PTD.GO_ROOTS)
                else:
                    primary = fid
                    closure = {primary}
                exact_keys.add((gid, pid, primary))
                for anc in closure:
                    propagated_keys.add((gid, pid, anc))
                    if rec.get("provenance"):
                        provenance_by_key[(gid, pid, anc)].extend(rec["provenance"])
            else:
                norm_ec = normalise_function_id(raw_fid, "EC", go=None)
                exact_keys.add((gid, pid, norm_ec))
                propagated_keys.add((gid, pid, norm_ec))
                if rec.get("provenance"):
                    provenance_by_key[(gid, pid, norm_ec)].extend(rec["provenance"])

    counts["unique_exact_keys"] = len(exact_keys)
    counts["unique_propagated_keys"] = len(propagated_keys)
    return {
        "exact_keys": exact_keys,
        "propagated_keys": propagated_keys,
        "provenance_by_key": dict(provenance_by_key),
        "counts": dict(counts),
    }


def load_model_negatives(
    negatives_path: Path | str,
    go: PTD.GoView | None = None,
    target_set: str = "model",
) -> dict:
    """Load accepted model-set negatives from negative_reviews.tsv.

    Negatives are exact assertions and are NEVER propagated to GO ancestors.
    """
    p = Path(negatives_path)
    if p.is_dir():
        if (p / "negative_reviews.tsv").is_file():
            p = p / "negative_reviews.tsv"
        elif (p / "registry" / "negative_reviews.tsv").is_file():
            p = p / "registry" / "negative_reviews.tsv"
    if not p.is_file():
        raise FileNotFoundError(f"negatives file not found: {negatives_path}")

    negative_keys: set[tuple[str, str, str]] = set()
    records_by_key: dict[tuple[str, str, str], dict] = {}
    counts = Counter()

    with p.open(newline="") as fh:
        for rec in csv.DictReader(fh, delimiter="\t"):
            counts["total_rows"] += 1
            row_set = str(rec.get("set", "")).strip()
            if row_set and row_set != target_set:
                counts["skipped_other_set"] += 1
                continue
            disp = str(rec.get("disposition", "")).strip()
            if disp != "accepted":
                counts[f"skipped_disposition_{disp or 'empty'}"] += 1
                continue
            counts["accepted_model_rows"] += 1
            gid = str(rec["genome_id"]).strip()
            pid = str(rec["protein_id"]).strip()
            raw_fid = str(rec["function_id"]).strip().replace("_", ":")
            if raw_fid.upper().startswith("GO:"):
                fid = "GO:" + raw_fid[3:]
                if go is not None:
                    if go.obsolete(fid):
                        counts["dropped_obsolete_go"] += 1
                        continue
                    fid = go.primary(fid)
                norm_fid = fid
            else:
                norm_fid = normalise_function_id(raw_fid, "EC", go=None)

            key = (gid, pid, norm_fid)
            negative_keys.add(key)
            records_by_key[key] = {
                "genome_id": gid,
                "protein_id": pid,
                "accession": str(rec.get("accession", "")).strip(),
                "function_id": norm_fid,
                "raw_function_id": str(rec["function_id"]).strip(),
                "evidence": str(rec.get("evidence", "")).strip(),
                "pmids": str(rec.get("pmids", "")).strip(),
            }

    counts["unique_negative_keys"] = len(negative_keys)
    return {
        "negative_keys": negative_keys,
        "records_by_key": records_by_key,
        "counts": dict(counts),
    }


# ----------------------------------------------------------------------------- export receipt verification

def _resolve_features_root(export_dir: Path) -> Path:
    export_dir = Path(export_dir)
    if not export_dir.is_dir():
        raise FileNotFoundError(f"export directory does not exist: {export_dir}")
    direct = [p for p in export_dir.iterdir() if p.is_dir() and (p / "candidate_features.jsonl").is_file()]
    if direct:
        return export_dir
    sub = export_dir / "features"
    if sub.is_dir():
        nested = [p for p in sub.iterdir() if p.is_dir() and (p / "candidate_features.jsonl").is_file()]
        if nested:
            return sub
    raise ValueError(f"no <genome>/candidate_features.jsonl found under {export_dir}")


def _parse_receipt_hashes(receipt_path: Path) -> dict[str, str]:
    p = Path(receipt_path)
    if not p.is_file():
        raise FileNotFoundError(f"missing export receipt: {p}")
    text = p.read_text().strip()
    if not text:
        raise ValueError(f"empty export receipt: {p}")

    by_path: dict[str, str] = {}
    if text.startswith("{"):
        obj = json.loads(text)
        files = obj.get("files")
        if isinstance(files, list):
            for entry in files:
                by_path[str(entry["path"]).strip()] = str(entry["sha256"]).strip()
        elif isinstance(files, dict):
            for k, v in files.items():
                by_path[str(k).strip()] = str(v["sha256"] if isinstance(v, dict) else v).strip()
        elif isinstance(obj.get("genomes"), dict):
            for gid, ginfo in obj["genomes"].items():
                if isinstance(ginfo, str):
                    by_path[f"{gid}/candidate_features.jsonl"] = ginfo.strip()
                elif isinstance(ginfo, dict):
                    sha = (
                        ginfo.get("sha256")
                        or ginfo.get("candidate_features_sha256")
                        or ginfo.get("rows_sha256")
                    )
                    if sha:
                        by_path[f"{gid}/candidate_features.jsonl"] = str(sha).strip()
        else:
            for k, v in obj.items():
                if isinstance(v, str) and len(v.strip()) == 64:
                    by_path[str(k).strip()] = v.strip()
    else:
        for line in text.splitlines():
            if not line.strip():
                continue
            parts = line.strip().split(None, 1)
            if len(parts) != 2:
                raise ValueError(f"malformed line in export receipt {p}: {line!r}")
            sha, rel = parts[0].strip(), parts[1].strip()
            by_path[rel.lstrip("./")] = sha

    if not by_path:
        raise ValueError(f"no file hashes parsed from receipt: {p}")
    return by_path


def verify_export_receipt(export_dir: Path, receipt_path: Path) -> tuple[Path, list[str], dict[str, str]]:
    """Verify every genome's candidate_features.jsonl (and any other referenced files) against receipt_path."""
    features_root = _resolve_features_root(export_dir)
    expected_map = _parse_receipt_hashes(receipt_path)
    genome_dirs = sorted(p for p in features_root.iterdir() if p.is_dir() and (p / "candidate_features.jsonl").is_file())
    if not genome_dirs:
        raise ValueError(f"no genome export directories found in {features_root}")

    verified: dict[str, str] = {}
    genomes: list[str] = []

    for gdir in genome_dirs:
        gid = gdir.name
        genomes.append(gid)
        cfile = gdir / "candidate_features.jsonl"
        actual_sha = T.sha256_file(cfile)

        candidate_keys = (
            f"{gid}/candidate_features.jsonl",
            f"features/{gid}/candidate_features.jsonl",
            gid,
            str(cfile),
        )
        exp_sha = None
        matched_key = None
        for ck in candidate_keys:
            if ck in expected_map:
                exp_sha = expected_map[ck]
                matched_key = ck
                break
        if exp_sha is None:
            raise ValueError(f"export file for genome {gid} not listed in receipt {receipt_path}")
        if actual_sha != exp_sha:
            raise ValueError(
                f"sha256 mismatch for {cfile}: expected {exp_sha} from {receipt_path}, got {actual_sha}"
            )
        verified[matched_key or f"{gid}/candidate_features.jsonl"] = actual_sha

        # Also check export_manifest.json if present in the genome directory
        manifest_file = gdir / "export_manifest.json"
        if manifest_file.is_file():
            mobj = json.loads(manifest_file.read_text())
            if mobj.get("rows_sha256") and mobj["rows_sha256"] != actual_sha:
                raise ValueError(
                    f"sha256 mismatch between {cfile} ({actual_sha}) and {manifest_file} ({mobj['rows_sha256']})"
                )

    # Verify any additional files listed in expected_map that exist under export_dir, features_root, or receipt_path.parent
    search_roots = [Path(export_dir), features_root, Path(receipt_path).parent]
    for rel_key, exp_sha in sorted(expected_map.items()):
        if rel_key in verified:
            continue
        resolved = None
        for sroot in search_roots:
            cand = sroot / rel_key
            if cand.is_file():
                resolved = cand
                break
        if resolved is not None:
            actual = T.sha256_file(resolved)
            if actual != exp_sha:
                raise ValueError(
                    f"sha256 mismatch for {resolved}: expected {exp_sha} from {receipt_path}, got {actual}"
                )
            verified[rel_key] = actual

    return features_root, genomes, verified


# ----------------------------------------------------------------------------- streaming & labelling

def stream_and_label_export(
    features_root: Path,
    genomes: list[str],
    feature_names: list[str],
    pos_info: dict,
    neg_info: dict,
    go: PTD.GoView | None = None,
) -> dict:
    """Stream verified candidate_features.jsonl files, label each row (positive / negative / unlabelled),
    and exclude and list any row that is both positive and negative (conflict).
    """
    pos_keys = pos_info["propagated_keys"]
    pos_exact_keys = pos_info["exact_keys"]
    neg_keys = neg_info["negative_keys"]

    rows: list[dict] = []
    conflicts: list[dict] = []
    matched_pos_keys: set[tuple[str, str, str]] = set()
    matched_neg_keys: set[tuple[str, str, str]] = set()
    conflict_neg_keys: set[tuple[str, str, str]] = set()

    per_genome_counts: dict[str, Counter] = {g: Counter() for g in genomes}
    pooled_counts = Counter()

    for gid in genomes:
        cfile = features_root / gid / "candidate_features.jsonl"
        with cfile.open() as fh:
            for line in fh:
                if not line.strip():
                    continue
                erow = json.loads(line)
                row_gid = str(erow.get("genome_id") or gid).strip()
                pid = str(erow["protein_id"]).strip()
                ft = str(erow["function_type"]).strip()
                raw_fid = str(erow["function_id"]).strip()
                cand_key = str(erow["candidate_key"]).strip()
                role = str(erow["role"]).strip()

                per_genome_counts[gid]["total_export_rows"] += 1
                pooled_counts["total_export_rows"] += 1

                norm_fid = normalise_function_id(raw_fid, ft, go=go)
                lookup_key = (row_gid, pid, norm_fid)

                is_pos = lookup_key in pos_keys
                is_pos_exact = lookup_key in pos_exact_keys
                is_neg = lookup_key in neg_keys

                if is_pos and is_neg:
                    conflict_rec = {
                        "example_id": f"{row_gid}|{cand_key}",
                        "genome_id": row_gid,
                        "protein_id": pid,
                        "function_type": ft,
                        "function_id": norm_fid,
                        "raw_function_id": raw_fid,
                        "role": role,
                        "candidate_key": cand_key,
                        "is_exact_positive": is_pos_exact,
                        "reason": "positive_and_negative_conflict",
                    }
                    conflicts.append(conflict_rec)
                    conflict_neg_keys.add(lookup_key)
                    per_genome_counts[gid]["conflicts_excluded"] += 1
                    pooled_counts["conflicts_excluded"] += 1
                    continue

                if is_pos:
                    label = 1
                    label_basis = "positive_truth"
                    is_exp_neg = False
                    stratum_label = "positive"
                    matched_pos_keys.add(lookup_key)
                    per_genome_counts[gid]["positive_rows"] += 1
                    pooled_counts["positive_rows"] += 1
                    if is_pos_exact:
                        per_genome_counts[gid]["positive_rows_exact"] += 1
                        pooled_counts["positive_rows_exact"] += 1
                    else:
                        per_genome_counts[gid]["positive_rows_propagated_ancestor"] += 1
                        pooled_counts["positive_rows_propagated_ancestor"] += 1
                elif is_neg:
                    label = 0
                    label_basis = "explicit_negative_truth"
                    is_exp_neg = True
                    stratum_label = "negative"
                    matched_neg_keys.add(lookup_key)
                    per_genome_counts[gid]["negative_rows"] += 1
                    pooled_counts["negative_rows"] += 1
                else:
                    label = 0
                    label_basis = "closed_world_unconfirmed"
                    is_exp_neg = False
                    stratum_label = "unlabelled"
                    per_genome_counts[gid]["unlabelled_rows"] += 1
                    pooled_counts["unlabelled_rows"] += 1

                feats = erow["features"]
                x = np.array([float(feats[n]) for n in feature_names], dtype=float)
                # Keep canonical function_id (with EC: prefix removed for EC; PE.norm_function_id handles both)
                rows.append(
                    {
                        "example_id": f"{row_gid}|{cand_key}",
                        "genome_id": row_gid,
                        "group": row_gid,
                        "split": "model_stratum",
                        "stratum": STRATUM_LABEL,
                        "stratum_label": stratum_label,
                        "label": label,
                        "label_basis": label_basis,
                        "is_explicit_negative": is_exp_neg,
                        "is_exact_positive": bool(is_pos_exact),
                        "x": x,
                        "locus_group_id": None,
                        "protein_id": pid,
                        "role": role,
                        "function_type": ft,
                        "function_id": norm_fid,
                        "raw_function_id": raw_fid,
                        "candidate_key": cand_key,
                        "sequence_sources": erow.get("sequence_sources"),
                        "sequence_contracts": erow.get("sequence_contracts"),
                        "quality_evidence": erow.get("quality_evidence"),
                        "cohort_context_mode": erow.get("cohort_context_mode", "UNIVERSAL"),
                    }
                )
                per_genome_counts[gid]["kept_rows"] += 1
                pooled_counts["kept_rows"] += 1

    unjoined_negatives = [
        neg_info["records_by_key"][k]
        for k in sorted(neg_keys - matched_neg_keys - conflict_neg_keys)
    ]

    return {
        "rows": rows,
        "conflicts": conflicts,
        "unjoined_negatives": unjoined_negatives,
        "matched_positive_keys": len(matched_pos_keys),
        "matched_negative_keys": len(matched_neg_keys),
        "pooled_counts": dict(pooled_counts),
        "per_genome_counts": {g: dict(c) for g, c in per_genome_counts.items()},
    }


# ----------------------------------------------------------------------------- negative-bearing competitions

def evaluate_negative_competitions(art: dict, scored_rows: list[dict]) -> dict:
    """Evaluate competitions containing at least one explicit negative member, plus all-negative competitions.

    Reports the rank of each negative member among the supported members of its competition and P(null)
    when the artifact has a trained null layer (plus naive null 1 - max p_hat for all artifacts).
    """
    by_comp: dict[tuple, list[dict]] = defaultdict(list)
    for sr in scored_rows:
        ck = PT.competition_key(sr)
        if ck is not None:
            by_comp[ck].append(sr)

    locus_trained = bool(art.get("locus_model", {}).get("trained"))
    thr = float(art.get("display_threshold", {}).get("value", 0.5))

    neg_competitions: list[dict] = []
    all_neg_competitions: list[dict] = []
    neg_ranks: list[float] = []
    neg_rel_ranks: list[float] = []
    neg_recip_ranks: list[float] = []
    neg_ranked_first: list[float] = []
    p_null_vals: list[float] = []
    naive_null_vals: list[float] = []

    for ck in sorted(by_comp):
        members = by_comp[ck]
        neg_members = [m for m in members if PE._is_explicit_negative(m)]
        if not neg_members:
            continue

        pos_members = [m for m in members if PE._is_positive(m)]
        unl_members = [m for m in members if PE._is_unlabelled(m)]
        sup_members = [m for m in members if m.get("supported") and m.get("probability") is not None]
        all_supported = len(sup_members) == len(members)
        is_all_negative = len(neg_members) == len(members)

        if all_supported and locus_trained and len(members) >= 1:
            p0, p_cands = PE._assign_locus_for_members(art, members)
            p_null = float(p0)
            p_null_vals.append(p_null)
        else:
            p_null = None
            p_cands = None

        if sup_members:
            naive_null = float(1.0 - max(float(m["probability"]) for m in sup_members))
            naive_null_vals.append(naive_null)
        else:
            naive_null = None

        sup_probs = [float(m["probability"]) for m in sup_members]
        neg_member_details: list[dict] = []
        for idx_m, m in enumerate(members):
            if not PE._is_explicit_negative(m):
                continue
            if m.get("supported") and m.get("probability") is not None:
                mp = float(m["probability"])
                higher = sum(1 for sp in sup_probs if sp > mp)
                tied = sum(1 for sp in sup_probs if sp == mp) - 1
                rk = 1.0 + higher + 0.5 * max(tied, 0)
                rel_rk = rk / len(sup_probs) if sup_probs else None
                neg_ranks.append(rk)
                if rel_rk is not None:
                    neg_rel_ranks.append(rel_rk)
                neg_recip_ranks.append(1.0 / rk)
                neg_ranked_first.append(1.0 if higher == 0 and tied <= 0 else 0.0)
            else:
                mp = None
                rk = None
                rel_rk = None

            neg_member_details.append(
                {
                    "example_id": m.get("example_id"),
                    "protein_id": m.get("protein_id"),
                    "candidate_key": m.get("candidate_key"),
                    "supported": bool(m.get("supported")),
                    "unsupported_reason": m.get("unsupported_reason"),
                    "probability": mp,
                    "below_threshold": (mp < thr) if mp is not None else None,
                    "rank_in_competition": rk,
                    "relative_rank": rel_rk,
                    "locus_posterior": float(p_cands[idx_m]) if p_cands is not None else None,
                }
            )

        comp_rec = {
            "competition": list(ck),
            "genome_id": ck[0],
            "function_type": ck[1],
            "function_id": ck[2],
            "role": ck[3],
            "n_members": len(members),
            "n_supported_members": len(sup_members),
            "n_positive_members": len(pos_members),
            "n_negative_members": len(neg_members),
            "n_unlabelled_members": len(unl_members),
            "is_all_negative": is_all_negative,
            "all_members_supported": all_supported,
            "p_null": p_null,
            "naive_p_null": naive_null,
            "negative_members": neg_member_details,
        }
        neg_competitions.append(comp_rec)
        if is_all_negative:
            all_neg_competitions.append(comp_rec)

    return {
        "locus_trained": locus_trained,
        "n_negative_bearing_competitions": len(neg_competitions),
        "n_all_negative_competitions": len(all_neg_competitions),
        "n_supported_negative_members": len(neg_ranks),
        "mean_negative_rank": float(np.mean(neg_ranks)) if neg_ranks else None,
        "mean_negative_relative_rank": float(np.mean(neg_rel_ranks)) if neg_rel_ranks else None,
        "mean_negative_reciprocal_rank": float(np.mean(neg_recip_ranks)) if neg_recip_ranks else None,
        "fraction_negative_ranked_first": float(np.mean(neg_ranked_first)) if neg_ranked_first else None,
        "mean_p_null": float(np.mean(p_null_vals)) if p_null_vals else None,
        "mean_naive_p_null": float(np.mean(naive_null_vals)) if naive_null_vals else None,
        "competitions": neg_competitions,
        "all_negative_competitions": all_neg_competitions,
    }


def _is_pu_artifact(arm_name: str, art: dict) -> bool:
    if art.get("training", {}).get("pu") is not None:
        return True
    if art.get("probability_target") == PT.PU_PROBABILITY_TARGET:
        return True
    lower = arm_name.lower()
    if "pu" in lower or arm_name.startswith(("A3", "A4", "A5", "S1", "S2")):
        return True
    return False


def _evaluate_subset(
    rows_sub: list[dict],
    scored_by_arm: dict[str, list[dict]],
    models_art: dict[str, dict],
    arm_c_info: dict[str, tuple[float, str]],
    pu_arms: list[str],
    pn_arms: list[str],
) -> dict:
    """Compute per-arm metrics, negative specificity, negative competitions, null layer, and pairwise contrasts."""
    arms_out: dict[str, dict] = {}
    for arm_name, art in models_art.items():
        c_val, c_note = arm_c_info[arm_name]
        thr = float(art.get("display_threshold", {}).get("value", 0.5))
        scored = scored_by_arm[arm_name]
        cand_m = PE.candidate_metrics(scored, c=c_val, threshold=thr)
        neg_spec = PE.negative_specificity(scored, threshold=thr)
        neg_comps = evaluate_negative_competitions(art, scored)
        null_m = PE.null_layer_metrics(art, rows_sub, split=None, c=c_val)
        arms_out[arm_name] = {
            "arm": arm_name,
            "kind": art.get("kind"),
            "threshold": thr,
            "c": c_val,
            "c_note": c_note,
            "pi": art.get("training", {}).get("pu", {}).get("pi"),
            "locus_trained": bool(art.get("locus_model", {}).get("trained")),
            "candidate_metrics": cand_m,
            "negative_specificity": neg_spec,
            "negative_competitions": neg_comps,
            "null_layer_metrics": null_m,
        }

    contrasts: dict[str, dict] = {}
    if pu_arms and pn_arms:
        pairs_to_run = [(pu_a, pn_a) for pu_a in pu_arms for pn_a in pn_arms]
    else:
        all_names = list(models_art.keys())
        pairs_to_run = [
            (all_names[i], all_names[j])
            for i in range(len(all_names))
            for j in range(len(all_names))
            if i != j
        ]

    for arm_a, arm_b in pairs_to_run:
        thr_a = float(models_art[arm_a].get("display_threshold", {}).get("value", 0.5))
        thr_b = float(models_art[arm_b].get("display_threshold", {}).get("value", 0.5))
        contrasts[f"{arm_a}_vs_{arm_b}"] = PE.pairwise_contrast(
            scored_by_arm[arm_a], scored_by_arm[arm_b], thr_a, thr_b
        )

    return {
        "n_rows": len(rows_sub),
        "arms": arms_out,
        "pairwise_contrasts": contrasts,
    }


def _fmt(v, nd: int = 4) -> str:
    if v is None:
        return "—"
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, int):
        return str(v)
    if isinstance(v, float):
        return f"{v:.{nd}f}"
    return str(v)


def write_markdown_tables(ev: dict, out_dir: Path) -> Path:
    """Render pooled and per-genome markdown tables to out_dir / 'report_tables.md'."""
    md_path = Path(out_dir) / "report_tables.md"
    if md_path.exists():
        raise FileExistsError(f"refusing to overwrite existing file: {md_path}")

    lines = [
        "# Model-Genome Stratum Evaluation Tables (Row 4, §4.3, P5)",
        "",
        f"- Stratum classification: **{ev['stratum']}**",
        f"- Genomes: {', '.join(ev['genomes'])}",
        f"- Total export rows: {ev['data_summary']['pooled_counts'].get('total_export_rows', 0)}",
        f"- Kept rows: {ev['data_summary']['pooled_counts'].get('kept_rows', 0)} "
        f"(positives: {ev['data_summary']['pooled_counts'].get('positive_rows', 0)} "
        f"[exact: {ev['data_summary']['pooled_counts'].get('positive_rows_exact', 0)}, "
        f"propagated ancestors: {ev['data_summary']['pooled_counts'].get('positive_rows_propagated_ancestor', 0)}], "
        f"negatives: {ev['data_summary']['pooled_counts'].get('negative_rows', 0)}, "
        f"unlabelled: {ev['data_summary']['pooled_counts'].get('unlabelled_rows', 0)})",
        f"- Excluded positive/negative conflicts: {len(ev['conflicts'])}",
        f"- Unjoined accepted model negatives: {len(ev['unjoined_negatives'])}",
        "",
        "## 1. Pooled Candidate-Layer & Negative-Specificity Metrics",
        "",
        "| Arm | Kind | c (Note) | Threshold | Supported | Abstained | Calls | Pos Called | Unconf Called | Neg Called | Supported Prec | PU Expected Prec | Recall | AUROC (Pos vs Unl) | PU-ECE | Neg Specificity |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]

    pooled_arms = ev["pooled"]["arms"]
    for arm_name, ares in pooled_arms.items():
        cm = ares["candidate_metrics"]
        ns = ares["negative_specificity"]
        c_str = f"{_fmt(ares['c'], 4)} ({'PU' if 'training.pu.c' in ares['c_note'] else 'PN c=1'})"
        lines.append(
            f"| {arm_name} | {ares['kind']} | {c_str} | {_fmt(ares['threshold'])} | "
            f"{cm['n']} | {cm['abstained']} | {cm['n_calls']} | {cm['labelled_positives_called']} | "
            f"{cm['unconfirmed_called']} | {cm['explicit_negatives_called']} | "
            f"{_fmt(cm['supported_precision'])} | {_fmt(cm['pu_expected_precision'])} | "
            f"{_fmt(cm['recall'])} | {_fmt(cm['auroc_positive_vs_unlabelled'])} | "
            f"{_fmt(cm['pu_ece'])} | {_fmt(ns['specificity'])} ({ns['n_below_threshold']}/{ns['n_supported']}) |"
        )

    lines.extend(
        [
            "",
            "## 2. Per-Genome Candidate-Layer & Negative-Specificity Metrics",
            "",
            "| Genome | Arm | Supported | Abstained | Positives (Sup/All) | Negatives (Sup/All) | Calls | Pos Called | Neg Called | Supported Prec | Recall | AUROC | Neg Specificity |",
            "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for gid in ev["genomes"]:
        garms = ev["per_genome"][gid]["arms"]
        for arm_name, ares in garms.items():
            cm = ares["candidate_metrics"]
            ns = ares["negative_specificity"]
            lines.append(
                f"| {gid} | {arm_name} | {cm['n']} | {cm['abstained']} | "
                f"{cm['positives']}/{cm['all_split_positives']} | {ns['n_supported']}/{ns['n_explicit_negatives']} | "
                f"{cm['n_calls']} | {cm['labelled_positives_called']} | {cm['explicit_negatives_called']} | "
                f"{_fmt(cm['supported_precision'])} | {_fmt(cm['recall'])} | "
                f"{_fmt(cm['auroc_positive_vs_unlabelled'])} | {_fmt(ns['specificity'])} |"
            )

    lines.extend(
        [
            "",
            "## 3. Pooled Pairwise Contrasts (PU vs PN)",
            "",
            "| Contrast (A vs B) | Both Supported | Corrections (A only) | Corrections (B only) | Harms on Negatives (A only) | Harms on Negatives (B only) | Unconfirmed (A only) | Unconfirmed (B only) |",
            "|---|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for cname, cval in ev["pooled"]["pairwise_contrasts"].items():
        lines.append(
            f"| {cname} | {cval['n_both_supported']} | {cval['corrections_a_over_b']} | "
            f"{cval['corrections_b_over_a']} | {cval['harms_a_over_b']} | {cval['harms_b_over_a']} | "
            f"{cval['unconfirmed_called_only_by_a']} | {cval['unconfirmed_called_only_by_b']} |"
        )

    lines.extend(
        [
            "",
            "## 4. Negative-Bearing Competitions & Null-Layer Summary (Pooled)",
            "",
            "| Arm | Locus Trained | Neg-Bearing Comps | All-Neg Comps | Sup Neg Members | Mean Neg Rank | Mean Neg Rel Rank | Neg Ranked First | Mean P(null) on Neg Comps | Mean Naive Null on Neg Comps | All Groups (>=2) Marginal NLL |",
            "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for arm_name, ares in pooled_arms.items():
        nc = ares["negative_competitions"]
        nl = ares["null_layer_metrics"]
        lines.append(
            f"| {arm_name} | {_fmt(ares['locus_trained'])} | {nc['n_negative_bearing_competitions']} | "
            f"{nc['n_all_negative_competitions']} | {nc['n_supported_negative_members']} | "
            f"{_fmt(nc['mean_negative_rank'])} | {_fmt(nc['mean_negative_relative_rank'])} | "
            f"{_fmt(nc['fraction_negative_ranked_first'])} | {_fmt(nc['mean_p_null'])} | "
            f"{_fmt(nc['mean_naive_p_null'])} | {_fmt(nl.get('marginal_nll'))} |"
        )

    lines.extend(
        [
            "",
            "## 5. Abstention Reasons by Arm (Pooled)",
            "",
            "| Arm | Total Rows | Supported | Abstained | Reasons |",
            "|---|---:|---:|---:|---|",
        ]
    )
    for arm_name, ares in pooled_arms.items():
        cm = ares["candidate_metrics"]
        r_str = ", ".join(f"{k}:{v}" for k, v in sorted(cm["unsupported_reasons"].items())) or "none"
        lines.append(f"| {arm_name} | {cm['total_rows']} | {cm['n']} | {cm['abstained']} | {r_str} |")

    md_path.write_text("\n".join(lines) + "\n")
    return md_path


# ----------------------------------------------------------------------------- main runner

def _parse_models_arg(model_specs: list[str]) -> dict[str, Path]:
    if not model_specs:
        raise ValueError("at least one --models <arm=path> specification is required")
    out: dict[str, Path] = {}
    for spec in model_specs:
        if "=" not in spec:
            raise ValueError(f"invalid --models entry {spec!r}; expected <arm>=<path>")
        arm, path_str = spec.split("=", 1)
        arm = arm.strip()
        p = Path(path_str.strip())
        if not arm:
            raise ValueError(f"empty arm name in --models entry {spec!r}")
        if arm in out:
            raise ValueError(f"duplicate arm name {arm!r} in --models")
        if not p.is_file():
            raise FileNotFoundError(f"model file for arm {arm!r} does not exist: {p}")
        out[arm] = p
    return out


def run(args) -> dict:
    out_dir = Path(args.out)
    if out_dir.exists():
        raise FileExistsError(f"refusing to overwrite existing output directory: {out_dir}")

    models_paths = _parse_models_arg(args.models)
    features_root, genomes, verified_export_hashes = verify_export_receipt(args.export_dir, args.export_receipt)

    models_art: dict[str, dict] = {}
    models_sha: dict[str, str] = {}
    feature_names: list[str] | None = None
    for arm_name, mpath in models_paths.items():
        art = T.load_artifact(mpath)
        models_art[arm_name] = art
        models_sha[arm_name] = T.sha256_file(mpath)
        art_names = list(art["contract"]["features"])
        if feature_names is None:
            feature_names = art_names
        elif art_names != feature_names:
            raise ValueError(f"model {arm_name} has a different feature contract column list than previous models")

    assert feature_names is not None

    go = PTD.GoView(args.go_obo)
    pos_info = load_model_positives(args.positives, go=go, target_set="model")
    neg_info = load_model_negatives(args.negatives, go=go, target_set="model")

    streamed = stream_and_label_export(
        features_root,
        genomes,
        feature_names,
        pos_info,
        neg_info,
        go=go,
    )
    all_rows = streamed["rows"]

    arm_c_info: dict[str, tuple[float, str]] = {}
    pu_arms: list[str] = []
    pn_arms: list[str] = []
    for arm_name, art in models_art.items():
        pu_block = art.get("training", {}).get("pu")
        if isinstance(pu_block, dict) and pu_block.get("c") is not None:
            arm_c_info[arm_name] = (float(pu_block["c"]), "from artifact training.pu.c")
        else:
            arm_c_info[arm_name] = (1.0, "c = 1.0 (PN artifact has no training.pu.c; closed-world reference)")
        if _is_pu_artifact(arm_name, art):
            pu_arms.append(arm_name)
        else:
            pn_arms.append(arm_name)

    scored_all_by_arm: dict[str, list[dict]] = {}
    for arm_name, art in models_art.items():
        scored_all_by_arm[arm_name] = PE.score(art, all_rows)

    pooled_eval = _evaluate_subset(all_rows, scored_all_by_arm, models_art, arm_c_info, pu_arms, pn_arms)

    per_genome_eval: dict[str, dict] = {}
    for gid in genomes:
        g_rows = [r for r in all_rows if r["genome_id"] == gid]
        g_scored = {
            arm_name: [sr for sr in scored_all_by_arm[arm_name] if sr["genome_id"] == gid]
            for arm_name in models_art
        }
        per_genome_eval[gid] = _evaluate_subset(g_rows, g_scored, models_art, arm_c_info, pu_arms, pn_arms)

    out_dir.mkdir(parents=True, exist_ok=False)
    evaluation = {
        "schema": SCHEMA,
        "plan": "tasks/v3-row4-pu-learning-2026-11-02.md",
        "stratum": STRATUM_LABEL,
        "genomes": genomes,
        "inputs": {
            "export_dir": str(args.export_dir),
            "features_root": str(features_root),
            "export_receipt": str(args.export_receipt),
            "export_receipt_sha256": T.sha256_file(args.export_receipt),
            "verified_export_hashes": verified_export_hashes,
            "positives_path": str(args.positives),
            "positives_sha256": T.sha256_file(args.positives),
            "negatives_path": str(args.negatives),
            "negatives_sha256": T.sha256_file(args.negatives),
            "go_obo_path": str(args.go_obo),
            "go_obo_sha256": T.sha256_file(args.go_obo),
            "models": {k: {"path": str(models_paths[k]), "sha256": models_sha[k]} for k in models_paths},
        },
        "data_summary": {
            "positives_loaded": pos_info["counts"],
            "negatives_loaded": neg_info["counts"],
            "matched_positive_keys": streamed["matched_positive_keys"],
            "matched_negative_keys": streamed["matched_negative_keys"],
            "pooled_counts": streamed["pooled_counts"],
            "per_genome_counts": streamed["per_genome_counts"],
        },
        "conflicts": streamed["conflicts"],
        "unjoined_negatives": streamed["unjoined_negatives"],
        "pooled": pooled_eval,
        "per_genome": per_genome_eval,
    }

    eval_path = out_dir / "evaluation.json"
    eval_path.write_text(T.canonical_json(evaluation))
    write_markdown_tables(evaluation, out_dir)
    return evaluation


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--export-dir", type=Path, required=True, help="Directory with <genome>/candidate_features.jsonl")
    p.add_argument("--export-receipt", type=Path, required=True, help="Receipt (sha256 list or JSON) for export files")
    p.add_argument("--positives", type=Path, required=True, help="positives.jsonl (set == 'model')")
    p.add_argument("--negatives", type=Path, required=True, help="negative_reviews.tsv (disposition == 'accepted', set == 'model')")
    p.add_argument("--go-obo", type=Path, required=True, help="Pinned go-basic.obo for positive ancestor propagation")
    p.add_argument("--models", action="extend", nargs="+", required=True, help="One or more <arm>=<model_json_path> pairs")
    p.add_argument("--out", type=Path, required=True, help="Output directory (must not exist)")
    return p


def main(argv=None) -> int:
    ap = build_parser()
    args = ap.parse_args(argv)
    try:
        run(args)
    except (ValueError, FileExistsError, FileNotFoundError, T.ProvenanceError, PTD.PreparationError) as exc:
        print(f"REFUSED: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
