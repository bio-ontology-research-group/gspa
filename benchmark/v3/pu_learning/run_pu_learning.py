#!/usr/bin/env python3
"""Runner for positive-unlabelled (PU) training and evaluation of the learned stage (row 4).

Stages (each refuses to overwrite its outputs; writes JSON with sorted keys):
  prepare   Verify input hashes against restore.json, check split disjointness, write registration.json.
  train     Fit A3 (pu logistic), A4 (pu_rank logistic), A5 (pu gbdt), S1 prior sweep,
            S2 unlabelled export pool (family_excluded with --export-dir), and S3 leave-one-genus-out;
            write frozen artifacts and SHA256SUMS.models before any test row is scored.
  evaluate  Verify SHA256SUMS.models, compute test metrics 1-7 for A1-A5 and S1-S3, write evaluation.json
            and per-arm evaluation/<arm>.json.
  report    Render markdown tables and P0-P8 checklist to <work>/<regime>/report_tables.md.
"""
from __future__ import annotations

import argparse
import copy
import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "benchmark" / "v3"))
sys.path.insert(0, str(REPO / "scripts"))
import train_context_model as T  # noqa: E402
from pu_learning import pu_evaluate as PE  # noqa: E402
from pu_learning import pu_trainer as PT  # noqa: E402

DEFAULT_WORK = Path("/mnt/data2/leechuck/gspa-v3-data/pu-learning-v1")
DEFAULT_SPLITS = REPO / "benchmark/v3/results/training-cohort/splits.tsv"
DEFAULT_NEG_REGISTRY = REPO / "benchmark/v3/results/complete-locus-labels-v1/registry"
DEFAULT_EXPORT_RECEIPT = REPO / "benchmark/v3/results/complete-locus-labels-v1/inputs/restore.json"
PLAN_SHA256_PATH = REPO / "benchmark/v3/results/pu-learning-v1/PLAN_SHA256"


# ----------------------------------------------------------------------------- paths & receipts

def regime_input_paths(work: Path, regime: str, contract_arg: Path | None) -> dict[str, Path]:
    base = work / "inputs" / "family-excluded-default"
    if contract_arg is None:
        contract_path = base / "feature_contract.json"
    elif contract_arg.is_absolute() or contract_arg.exists():
        contract_path = contract_arg
    else:
        contract_path = work / contract_arg

    if regime == "family_excluded":
        return {
            "base": base,
            "contract": contract_path,
            "table": base / "table" / "training_table.jsonl",
            "manifest": base / "table" / "training_table_manifest.json",
            "pn_logistic": base / "models" / "model_logistic.json",
            "pn_gbdt": base / "models" / "model_gbdt.json",
            "restore": base / "restore.json",
        }
    if regime == "homology":
        return {
            "base": base,
            "contract": contract_path,
            "table": base / "homology" / "table" / "training_table.jsonl",
            "manifest": base / "homology" / "table" / "training_table_manifest.json",
            "pn_logistic": base / "homology" / "models" / "model_logistic.json",
            "pn_gbdt": base / "homology" / "models" / "model_gbdt.json",
            "restore": base / "restore.json",
        }
    raise T.ProvenanceError(f"unknown regime {regime!r}")


def verify_restore_receipt(base_dir: Path, receipt_path: Path, only_rel_paths: set[str] | None = None) -> dict[str, str]:
    """Verify file sha256 values against a restore.json receipt (fail closed)."""
    if not receipt_path.exists():
        raise T.ProvenanceError(f"missing restore receipt: {receipt_path}")
    receipt = json.loads(receipt_path.read_text())
    files = receipt.get("files")
    if not isinstance(files, list) or not files:
        raise T.ProvenanceError(f"invalid restore receipt (no files list): {receipt_path}")

    by_rel = {str(entry["path"]): str(entry["sha256"]) for entry in files}
    targets = sorted(only_rel_paths) if only_rel_paths is not None else sorted(by_rel)
    verified: dict[str, str] = {}
    for rel in targets:
        if rel not in by_rel:
            raise T.ProvenanceError(f"file {rel!r} not listed in receipt {receipt_path}")
        full = base_dir / rel
        if not full.exists():
            raise T.ProvenanceError(f"missing input file {full} required by receipt {receipt_path}")
        actual = T.sha256_file(full)
        expected = by_rel[rel]
        if actual != expected:
            raise T.ProvenanceError(
                f"sha256 mismatch for {full}: expected {expected} from {receipt_path}, got {actual}"
            )
        verified[rel] = actual
    return verified


def read_splits_metadata(splits_path: Path, rows: list[dict]) -> dict[str, dict]:
    """Read genome metadata (split, group, species, genus, family) from splits.tsv, with row fallback for synthetic tests."""
    info: dict[str, dict] = {}
    if splits_path and splits_path.exists():
        with open(splits_path, newline="") as fh:
            for rec in csv.DictReader(fh, delimiter="\t"):
                gid = str(rec["genome_id"]).strip()
                info[gid] = {
                    "genome_id": gid,
                    "split": str(rec.get("split", "")).strip(),
                    "group": str(rec.get("group", "")).strip(),
                    "source": str(rec.get("source", "")).strip(),
                    "species": str(rec.get("species", "")).strip(),
                    "genus": str(rec.get("genus", "")).strip() or gid,
                    "family": str(rec.get("family", "")).strip() or str(rec.get("group", "")).strip(),
                }
    for r in rows:
        gid = str(r["genome_id"])
        if gid not in info:
            info[gid] = {
                "genome_id": gid,
                "split": str(r["split"]),
                "group": str(r["group"]),
                "source": str(r.get("dataset", "synthetic")),
                "species": gid,
                "genus": str(r.get("genus") or r["group"]),
                "family": str(r["group"]),
            }
    return info


def _write_json_no_overwrite(path: Path, obj: dict) -> str:
    if path.exists():
        raise FileExistsError(f"refusing to overwrite existing file: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    text = T.canonical_json(obj)
    path.write_text(text)
    return T.sha256_file(path)


# ----------------------------------------------------------------------------- Stage 1: prepare

def stage_prepare(args) -> dict:
    rdir = args.work / args.regime
    reg_path = rdir / "registration.json"
    if reg_path.exists():
        raise FileExistsError(f"refusing to overwrite existing registration: {reg_path}")

    paths = regime_input_paths(args.work, args.regime, args.contract)
    receipt_path = args.restore_receipt or paths["restore"]
    verified_hashes = verify_restore_receipt(paths["base"], receipt_path)

    contract = T.load_contract(paths["contract"])
    rows = T.read_table(paths["table"], contract)
    T.check_disjoint(rows, ["group", "genome_id"])
    identity = T.split_identity(rows, ["group", "genome_id"])
    splits_meta = read_splits_metadata(args.splits, rows)

    per_split = {}
    for s in T.SPLITS:
        sub = [r for r in rows if r["split"] == s]
        pos = sum(1 for r in sub if r["label"] == 1 and r.get("label_basis") != "explicit_negative_truth")
        unl = sum(1 for r in sub if r["label"] == 0 and r.get("label_basis") != "explicit_negative_truth")
        exp_neg = sum(1 for r in sub if r.get("label_basis") == "explicit_negative_truth")

        by_comp: dict[tuple, list[dict]] = defaultdict(list)
        for r in sub:
            if r.get("label_basis") == "explicit_negative_truth":
                continue
            ck = PT.competition_key(r)
            if ck is not None:
                by_comp[ck].append(r)

        ge2 = [ms for ms in by_comp.values() if len(ms) >= 2]
        ge2_pos1 = sum(1 for ms in ge2 if sum(1 for m in ms if m["label"] == 1) == 1)
        ge2_pos0 = sum(1 for ms in ge2 if sum(1 for m in ms if m["label"] == 1) == 0)
        ge2_multi = sum(1 for ms in ge2 if sum(1 for m in ms if m["label"] == 1) >= 2)

        per_split[s] = {
            "rows": int(len(sub)),
            "positives": int(pos),
            "unlabelled": int(unl),
            "explicit_negatives": int(exp_neg),
            "competitions_total": int(len(by_comp)),
            "competitions_ge2_members": int(len(ge2)),
            "competitions_ge2_single_positive": int(ge2_pos1),
            "competitions_ge2_no_label": int(ge2_pos0),
            "competitions_ge2_multi_positive": int(ge2_multi),
            "genomes": int(len({r["genome_id"] for r in sub})),
            "groups": int(len({r["group"] for r in sub})),
        }

    table_genomes = sorted({r["genome_id"] for r in rows})
    plan_sha_line = PLAN_SHA256_PATH.read_text().splitlines()[0].strip() if PLAN_SHA256_PATH.exists() else None

    reg = {
        "regime": args.regime,
        "plan_sha256_record": plan_sha_line,
        "restore_receipt": str(receipt_path),
        "verified_input_hashes": verified_hashes,
        "table_path": str(paths["table"]),
        "table_sha256": T.sha256_file(paths["table"]),
        "contract_path": str(paths["contract"]),
        "contract_sha256": T.sha256_file(paths["contract"]),
        "pn_logistic_path": str(paths["pn_logistic"]),
        "pn_logistic_sha256": T.sha256_file(paths["pn_logistic"]),
        "pn_gbdt_path": str(paths["pn_gbdt"]),
        "pn_gbdt_sha256": T.sha256_file(paths["pn_gbdt"]),
        "split_identity": identity,
        "per_split": per_split,
        "genus_by_genome": {g: splits_meta[g]["genus"] for g in table_genomes},
        "family_by_genome": {g: splits_meta[g]["family"] for g in table_genomes},
        "split_by_genome": {g: splits_meta[g]["split"] for g in table_genomes},
    }
    _write_json_no_overwrite(reg_path, reg)
    return reg


# ----------------------------------------------------------------------------- S2 export reader (streaming, verified)

def _collect_u2_extra_rows(
    train_calib_rows: list[dict],
    contract: dict,
    export_dir: Path,
    export_receipt: Path,
    max_extra_per_comp: int = 10,
) -> list[dict]:
    """Stream verified export files for train/calibration genomes and collect up to 10 extra unlabelled rows
    per table competition (genome, type, id, role), sorted deterministically by candidate_key.
    """
    names = contract["names"]
    comp_meta: dict[tuple, dict] = {}
    table_cand_keys_by_genome: dict[str, set[str]] = defaultdict(set)

    for r in train_calib_rows:
        gid = str(r["genome_id"])
        if r.get("candidate_key"):
            table_cand_keys_by_genome[gid].add(str(r["candidate_key"]))
        ck = PT.competition_key(r)
        if ck is not None and ck not in comp_meta:
            comp_meta[ck] = {"group": str(r["group"]), "split": str(r["split"])}

    target_genomes = sorted({ck[0] for ck in comp_meta})
    rel_paths_to_verify = {f"{gid}/candidate_features.jsonl" for gid in target_genomes}
    verify_restore_receipt(export_dir, export_receipt, only_rel_paths=rel_paths_to_verify)

    extra_rows: list[dict] = []
    for gid in target_genomes:
        genome_comps = {ck for ck in comp_meta if ck[0] == gid}
        existing_cks = table_cand_keys_by_genome[gid]
        pool_by_comp: dict[tuple, list[dict]] = defaultdict(list)

        feat_file = export_dir / gid / "candidate_features.jsonl"
        with open(feat_file) as fh:
            for line in fh:
                if not line.strip():
                    continue
                erow = json.loads(line)
                cand_key = str(erow["candidate_key"])
                if cand_key in existing_cks:
                    continue
                ft = str(erow["function_type"])
                fid = PE.norm_function_id(ft, erow["function_id"])
                role = str(erow["role"])
                ck = (gid, ft, fid, role)
                if ck in genome_comps:
                    pool_by_comp[ck].append(erow)

        for ck in sorted(pool_by_comp):
            meta = comp_meta[ck]
            candidates = sorted(pool_by_comp[ck], key=lambda item: str(item["candidate_key"]))[:max_extra_per_comp]
            for erow in candidates:
                cand_key = str(erow["candidate_key"])
                feats = erow["features"]
                x = np.array([float(feats[n]) for n in names], dtype=float)
                extra_rows.append(
                    {
                        "example_id": f"{gid}|{cand_key}",
                        "genome_id": gid,
                        "group": meta["group"],
                        "split": meta["split"],
                        "label": 0,
                        "label_basis": "no_truth_not_negative",
                        "x": x,
                        "locus_group_id": None,
                        "role": str(erow["role"]),
                        "function_type": str(erow["function_type"]),
                        "function_id": PE.norm_function_id(erow["function_type"], erow["function_id"]),
                        "candidate_key": cand_key,
                        "sequence_sources": erow.get("sequence_sources"),
                        "sequence_contracts": erow.get("sequence_contracts"),
                        "cohort_context_mode": erow.get("cohort_context_mode", "UNIVERSAL"),
                    }
                )
    return extra_rows


# ----------------------------------------------------------------------------- Stage 2: train

def _sweep_pi_tag(val: float) -> str:
    if abs(val - 0.20) < 1e-12:
        return "0.20"
    return f"{val:.4f}"


def stage_train(args) -> dict:
    rdir = args.work / args.regime
    mdir = rdir / "models"
    sums_path = rdir / "SHA256SUMS.models"
    if mdir.exists() or sums_path.exists():
        raise FileExistsError(f"refusing to overwrite existing models or SHA256SUMS.models in {rdir}")

    paths = regime_input_paths(args.work, args.regime, args.contract)
    receipt_path = args.restore_receipt or paths["restore"]
    verify_restore_receipt(paths["base"], receipt_path)

    contract = T.load_contract(paths["contract"])
    all_rows = T.read_table(paths["table"], contract)
    T.check_disjoint(all_rows, ["group", "genome_id"])
    identity = T.split_identity(all_rows, ["group", "genome_id"])
    splits_meta = read_splits_metadata(args.splits, all_rows)

    # Filter strictly to train and calibration rows before any fit.
    tr = [r for r in all_rows if r["split"] == "train"]
    ca = [r for r in all_rows if r["split"] == "calibration"]
    assert all(r["split"] in ("train", "calibration") for r in tr + ca)

    intended_use = "mechanism_test" if args.synthetic else "production"
    common_kw = dict(
        weighting="group",
        folds=args.folds,
        seed=args.seed,
        threshold_rule="max_f1",
        min_no_label_groups=args.min_no_label_groups,
        min_pos_groups=args.min_pos_groups,
        min_cal_no_label_groups=args.min_cal_no_label_groups,
        min_cal_pos_groups=args.min_cal_pos_groups,
        intended_use=intended_use,
        synthetic=bool(args.synthetic),
    )

    mdir.mkdir(parents=True, exist_ok=False)
    model_records: dict[str, dict] = {}

    def freeze_arm(arm_name: str, kind: str, art: dict):
        arm_dir = mdir / arm_name
        arm_dir.mkdir(parents=True, exist_ok=False)
        mpath = arm_dir / f"model_{kind}.json"
        sha = T.write_frozen(art, mpath)
        rel = str(mpath.relative_to(rdir))
        model_records[arm_name] = {
            "arm": arm_name,
            "kind": kind,
            "rel_path": rel,
            "sha256": sha,
            "pi": art.get("training", {}).get("pu", {}).get("pi"),
            "c": art.get("training", {}).get("pu", {}).get("c"),
            "locus_trained": bool(art.get("locus_model", {}).get("trained")),
        }

    # Main PU arms: A3 (pu, logistic), A4 (pu_rank, logistic), A5 (pu, gbdt)
    art_a3 = PT.train_pu(tr, ca, contract, "logistic", "pu", pi=None, identity=identity, **common_kw)
    freeze_arm("A3", "logistic", art_a3)

    art_a4 = PT.train_pu(tr, ca, contract, "logistic", "pu_rank", pi=None, identity=identity, **common_kw)
    freeze_arm("A4", "logistic", art_a4)

    art_a5 = PT.train_pu(tr, ca, contract, "gbdt", "pu", pi=None, identity=identity, **common_kw)
    freeze_arm("A5", "gbdt", art_a5)

    # S1 prior sweep for A3 at pi in {labelled_rate, 0.20, pi_hat, min(2*pi_hat, 0.8)}
    r_rate = float(art_a3["training"]["pu"]["labelled_rate"])
    pi_hat = float(art_a3["training"]["pu"]["pi"])
    pi_upper = float(min(2.0 * pi_hat, 0.80))
    s1_specs = [
        ("labelled_rate", r_rate),
        ("0.20", 0.20),
        ("pi_hat", pi_hat),
        ("min_2pi_hat_0.8", pi_upper),
    ]
    s1_manifest = []
    seen_s1_arms: dict[str, str] = {}
    for role_name, pi_val in s1_specs:
        tag = _sweep_pi_tag(pi_val)
        arm_name = f"S1_pi_{tag}"
        if arm_name in seen_s1_arms:
            arm_name = f"S1_pi_{pi_val:.6f}_{role_name}"
        seen_s1_arms[arm_name] = role_name
        art_s1 = PT.train_pu(tr, ca, contract, "logistic", "pu", pi=pi_val, identity=identity, **common_kw)
        freeze_arm(arm_name, "logistic", art_s1)
        s1_manifest.append({"role": role_name, "pi": float(pi_val), "arm": arm_name})

    # S2 unlabelled export pool (family_excluded only, when --export-dir is given)
    s2_info = {"trained": False, "reason": "not applicable or --export-dir not supplied"}
    if args.regime == "family_excluded" and args.export_dir is not None:
        exp_receipt = args.export_receipt or DEFAULT_EXPORT_RECEIPT
        extra_rows = _collect_u2_extra_rows(tr + ca, contract, Path(args.export_dir), Path(exp_receipt))
        extra_tr = [r for r in extra_rows if r["split"] == "train"]
        extra_ca = [r for r in extra_rows if r["split"] == "calibration"]
        tr_s2 = tr + extra_tr
        ca_s2 = ca + extra_ca
        assert all(r["split"] in ("train", "calibration") for r in tr_s2 + ca_s2)
        # Update split_identity for train/calibration while keeping test untouched from registration/initial identity
        ident_s2 = copy.deepcopy(identity)
        for s_name, s_rows in (("train", tr_s2), ("calibration", ca_s2)):
            ids = sorted(r["example_id"] for r in s_rows)
            ident_s2[s_name] = {
                "n": len(s_rows),
                "example_ids_sha256": T.hashlib.sha256("\n".join(ids).encode()).hexdigest(),
                "n_groups": len({r["group"] for r in s_rows}),
                "n_genome_ids": len({r["genome_id"] for r in s_rows}),
                "groups_sha256": T.hashlib.sha256("\n".join(sorted({r["group"] for r in s_rows})).encode()).hexdigest(),
                "genome_ids_sha256": T.hashlib.sha256(
                    "\n".join(sorted({r["genome_id"] for r in s_rows})).encode()
                ).hexdigest(),
            }
        art_s2 = PT.train_pu(tr_s2, ca_s2, contract, "logistic", "pu", pi=None, identity=ident_s2, **common_kw)
        freeze_arm("S2_U2", "logistic", art_s2)
        s2_info = {
            "trained": True,
            "arm": "S2_U2",
            "extra_train_rows": int(len(extra_tr)),
            "extra_calibration_rows": int(len(extra_ca)),
            "pi_reestimated": float(art_s2["training"]["pu"]["pi"]),
            "c_reestimated": float(art_s2["training"]["pu"]["c"]),
        }

    # S3 leave-one-genus-out (LOGO) for A1-style PN logistic and A3 PU logistic
    s3_folds = []
    if not args.skip_s3:
        genus_by_genome = {g: splits_meta[g]["genus"] for g in sorted({r["genome_id"] for r in all_rows})}
        frozen_cal_genomes = sorted({r["genome_id"] for r in ca})
        genera = sorted(set(genus_by_genome.values()))
        for genus in genera:
            held_out_genomes = sorted(g for g, gen in genus_by_genome.items() if gen == genus)
            rem_cal_genomes = [g for g in frozen_cal_genomes if g not in set(held_out_genomes)]
            if len(rem_cal_genomes) < 2:
                s3_folds.append(
                    {
                        "genus": genus,
                        "held_out_genomes": held_out_genomes,
                        "skipped": True,
                        "reason": f"excluding genus {genus} leaves {len(rem_cal_genomes)} calibration genome(s) (< 2 required)",
                    }
                )
                continue

            fold_tr: list[dict] = []
            fold_ca: list[dict] = []
            fold_te_ids: list[str] = []
            fold_te_groups: set[str] = set()
            fold_te_genomes: set[str] = set()
            for r in all_rows:
                gid = r["genome_id"]
                if gid in held_out_genomes:
                    fold_te_ids.append(r["example_id"])
                    fold_te_groups.add(r["group"])
                    fold_te_genomes.add(gid)
                elif gid in rem_cal_genomes:
                    rc = dict(r)
                    rc["split"] = "calibration"
                    fold_ca.append(rc)
                else:
                    rt = dict(r)
                    rt["split"] = "train"
                    fold_tr.append(rt)

            # Construct split identity without passing any fold test rows to train() or train_pu()
            fold_ident = {}
            for s_name, s_rows in (("train", fold_tr), ("calibration", fold_ca)):
                ids = sorted(r["example_id"] for r in s_rows)
                fold_ident[s_name] = {
                    "n": len(s_rows),
                    "example_ids_sha256": T.hashlib.sha256("\n".join(ids).encode()).hexdigest(),
                    "n_groups": len({r["group"] for r in s_rows}),
                    "n_genome_ids": len({r["genome_id"] for r in s_rows}),
                    "groups_sha256": T.hashlib.sha256(
                        "\n".join(sorted({r["group"] for r in s_rows})).encode()
                    ).hexdigest(),
                    "genome_ids_sha256": T.hashlib.sha256(
                        "\n".join(sorted({r["genome_id"] for r in s_rows})).encode()
                    ).hexdigest(),
                }
            fold_ident["test"] = {
                "n": len(fold_te_ids),
                "example_ids_sha256": T.hashlib.sha256("\n".join(sorted(fold_te_ids)).encode()).hexdigest(),
                "n_groups": len(fold_te_groups),
                "n_genome_ids": len(fold_te_genomes),
                "groups_sha256": T.hashlib.sha256("\n".join(sorted(fold_te_groups)).encode()).hexdigest(),
                "genome_ids_sha256": T.hashlib.sha256("\n".join(sorted(fold_te_genomes)).encode()).hexdigest(),
            }
            fold_ident["disjoint_keys_enforced"] = ["genome_id"]

            assert all(r["split"] in ("train", "calibration") for r in fold_tr + fold_ca)
            pn_arm = f"S3_logo_{genus}_pn"
            pu_arm = f"S3_logo_{genus}_pu"

            art_pn = T.train(
                fold_tr,
                fold_ca,
                contract,
                "logistic",
                weighting="group",
                folds=args.folds,
                seed=args.seed,
                threshold_rule="max_f1",
                min_null_groups=args.min_no_label_groups,
                min_pos_groups=args.min_pos_groups,
                identity=fold_ident,
                intended_use=intended_use,
                synthetic=bool(args.synthetic),
            )
            freeze_arm(pn_arm, "logistic", art_pn)

            art_pu = PT.train_pu(
                fold_tr,
                fold_ca,
                contract,
                "logistic",
                "pu",
                pi=None,
                identity=fold_ident,
                **common_kw,
            )
            freeze_arm(pu_arm, "logistic", art_pu)

            s3_folds.append(
                {
                    "genus": genus,
                    "held_out_genomes": held_out_genomes,
                    "calibration_genomes": rem_cal_genomes,
                    "skipped": False,
                    "pn_arm": pn_arm,
                    "pu_arm": pu_arm,
                }
            )

    manifest = {
        "regime": args.regime,
        "models": model_records,
        "s1_sweep": s1_manifest,
        "s2_u2": s2_info,
        "s3_logo_folds": s3_folds,
    }
    _write_json_no_overwrite(mdir / "training_manifest.json", manifest)

    # Write SHA256SUMS.models covering every model_*.json under models/
    lines = []
    for p in sorted(mdir.glob("*/model_*.json")):
        rel = str(p.relative_to(rdir))
        sha = T.sha256_file(p)
        lines.append(f"{sha}  {rel}")
    sums_path.write_text("\n".join(lines) + "\n")
    return manifest


# ----------------------------------------------------------------------------- Stage 3: evaluate

def _verify_models_sha256(rdir: Path) -> dict[str, str]:
    sums_path = rdir / "SHA256SUMS.models"
    if not sums_path.exists():
        raise T.ProvenanceError(f"missing {sums_path}; train stage must complete before evaluate")
    recorded: dict[str, str] = {}
    for line in sums_path.read_text().splitlines():
        if not line.strip():
            continue
        parts = line.strip().split(None, 1)
        if len(parts) != 2:
            raise T.ProvenanceError(f"malformed line in {sums_path}: {line!r}")
        sha, rel = parts[0].strip(), parts[1].strip()
        full = rdir / rel
        if not full.exists():
            raise T.ProvenanceError(f"frozen model missing: {full}")
        actual = T.sha256_file(full)
        if actual != sha:
            raise T.ProvenanceError(f"model file changed after freezing: {full} (expected {sha}, got {actual})")
        recorded[rel] = sha

    # Ensure no unrecorded model_*.json exists on disk
    on_disk = {str(p.relative_to(rdir)) for p in (rdir / "models").glob("*/model_*.json")}
    if on_disk != set(recorded):
        raise T.ProvenanceError(f"model files on disk differ from SHA256SUMS.models: {on_disk ^ set(recorded)}")
    return recorded


def _verify_sha256sum_file(root_dir: Path, sha256sum_path: Path, only_rel_paths: set[str] | None = None) -> dict[str, str]:
    """Verify a standard sha256sum file (<sha256>  <rel_path>) relative to root_dir."""
    if not sha256sum_path.exists():
        raise T.ProvenanceError(f"missing sha256 receipt: {sha256sum_path}")
    verified: dict[str, str] = {}
    for line in sha256sum_path.read_text().splitlines():
        if not line.strip():
            continue
        parts = line.strip().split(None, 1)
        if len(parts) != 2:
            raise T.ProvenanceError(f"malformed line in {sha256sum_path}: {line!r}")
        expected, rel = parts[0].strip(), parts[1].strip()
        if only_rel_paths is not None and rel not in only_rel_paths:
            continue
        full = root_dir / rel
        if not full.exists():
            raise T.ProvenanceError(f"missing file referenced in {sha256sum_path}: {full}")
        actual = T.sha256_file(full)
        if actual != expected:
            raise T.ProvenanceError(
                f"sha256 mismatch for {full}: expected {expected} from {sha256sum_path}, got {actual}"
            )
        verified[rel] = actual
    if only_rel_paths is not None:
        missing_in_receipt = only_rel_paths - set(verified)
        if missing_in_receipt:
            raise T.ProvenanceError(
                f"files not listed in {sha256sum_path}: {sorted(missing_in_receipt)}"
            )
    return verified


def _find_export_negative_rows(
    accepted_negatives: set[tuple[str, str, str]],
    table_rows: list[dict],
    splits_meta: dict[str, dict],
    contract: dict,
    export_dir: Path,
    export_receipt: Path,
) -> list[dict]:
    """Stream verified cohort export files to pull any accepted negative candidates not already in table_rows."""
    names = contract["names"]
    in_table = {
        (str(r["genome_id"]), PE.extract_protein_id(r), PE.norm_function_id(r.get("function_type"), r.get("function_id")))
        for r in table_rows
    }
    needed_by_genome: dict[str, set[tuple[str, str]]] = defaultdict(set)
    for gid, pid, fid in accepted_negatives - in_table:
        if (export_dir / gid / "candidate_features.jsonl").exists():
            needed_by_genome[gid].add((pid, fid))

    if not needed_by_genome:
        return []

    rel_paths = {f"{gid}/candidate_features.jsonl" for gid in sorted(needed_by_genome)}
    verify_restore_receipt(export_dir, export_receipt, only_rel_paths=rel_paths)

    found: list[dict] = []
    for gid in sorted(needed_by_genome):
        targets = needed_by_genome[gid]
        gmeta = splits_meta.get(gid, {"split": "train", "group": gid})
        feat_file = export_dir / gid / "candidate_features.jsonl"
        with open(feat_file) as fh:
            for line in fh:
                if not line.strip():
                    continue
                erow = json.loads(line)
                pid = str(erow["protein_id"]).strip()
                ft = str(erow["function_type"]).strip()
                fid = PE.norm_function_id(ft, erow["function_id"])
                if (pid, fid) in targets:
                    cand_key = str(erow["candidate_key"])
                    x = np.array([float(erow["features"][n]) for n in names], dtype=float)
                    found.append(
                        {
                            "example_id": f"{gid}|{cand_key}",
                            "genome_id": gid,
                            "group": gmeta["group"],
                            "split": gmeta["split"],
                            "label": 0,
                            "label_basis": "explicit_negative_truth",
                            "is_explicit_negative": True,
                            "x": x,
                            "locus_group_id": None,
                            "role": str(erow["role"]),
                            "function_type": ft,
                            "function_id": fid,
                            "candidate_key": cand_key,
                            "sequence_sources": erow.get("sequence_sources"),
                            "sequence_contracts": erow.get("sequence_contracts"),
                            "cohort_context_mode": erow.get("cohort_context_mode", "UNIVERSAL"),
                            "source_origin": "feature_export",
                        }
                    )
    return found


def _find_model_stratum_negative_rows(
    accepted_negatives: set[tuple[str, str, str]],
    existing_rows: list[dict],
    contract: dict,
    model_export_root: Path,
) -> list[dict]:
    """Stream verified model-stratum export files (<work>/model-stratum/export) for accepted negatives."""
    sha_file = model_export_root / "export.sha256"
    feat_root = model_export_root / "features"
    if not sha_file.exists() or not feat_root.exists():
        return []

    names = contract["names"]
    already = {
        (str(r["genome_id"]), PE.extract_protein_id(r), PE.norm_function_id(r.get("function_type"), r.get("function_id")))
        for r in existing_rows
    }
    needed_by_genome: dict[str, set[tuple[str, str]]] = defaultdict(set)
    for gid, pid, fid in accepted_negatives - already:
        if (feat_root / gid / "candidate_features.jsonl").exists():
            needed_by_genome[gid].add((pid, fid))

    if not needed_by_genome:
        return []

    rel_paths = {f"features/{gid}/candidate_features.jsonl" for gid in sorted(needed_by_genome)}
    _verify_sha256sum_file(model_export_root, sha_file, only_rel_paths=rel_paths)

    found: list[dict] = []
    for gid in sorted(needed_by_genome):
        targets = needed_by_genome[gid]
        feat_file = feat_root / gid / "candidate_features.jsonl"
        with open(feat_file) as fh:
            for line in fh:
                if not line.strip():
                    continue
                erow = json.loads(line)
                pid = str(erow["protein_id"]).strip()
                ft = str(erow["function_type"]).strip()
                fid = PE.norm_function_id(ft, erow["function_id"])
                if (pid, fid) in targets:
                    cand_key = str(erow["candidate_key"])
                    x = np.array([float(erow["features"][n]) for n in names], dtype=float)
                    found.append(
                        {
                            "example_id": f"{gid}|{cand_key}",
                            "genome_id": gid,
                            "group": gid,
                            "split": "heldout_development_examined",
                            "label": 0,
                            "label_basis": "explicit_negative_truth",
                            "is_explicit_negative": True,
                            "x": x,
                            "locus_group_id": None,
                            "role": str(erow["role"]),
                            "function_type": ft,
                            "function_id": fid,
                            "candidate_key": cand_key,
                            "sequence_sources": erow.get("sequence_sources"),
                            "sequence_contracts": erow.get("sequence_contracts"),
                            "cohort_context_mode": erow.get("cohort_context_mode", "UNIVERSAL"),
                            "source_origin": "model_stratum_export",
                        }
                    )
    return found


def stage_evaluate(args) -> dict:
    rdir = args.work / args.regime
    eval_path = rdir / "evaluation.json"
    eval_dir = rdir / "evaluation"
    if eval_path.exists() or eval_dir.exists():
        raise FileExistsError(f"refusing to overwrite existing evaluation in {rdir}")

    recorded_model_shas = _verify_models_sha256(rdir)
    paths = regime_input_paths(args.work, args.regime, args.contract)
    receipt_path = args.restore_receipt or paths["restore"]
    verify_restore_receipt(paths["base"], receipt_path)

    contract = T.load_contract(paths["contract"])
    rows = T.read_table(paths["table"], contract)
    splits_meta = read_splits_metadata(args.splits, rows)
    train_manifest = json.loads((rdir / "models" / "training_manifest.json").read_text())

    cal_rows = [r for r in rows if r["split"] == "calibration"]
    test_rows = [r for r in rows if r["split"] == "test"]

    # Load artifacts for A1, A2 (frozen PN) and trained PU arms
    arms_art: dict[str, dict] = {
        "A1": T.load_artifact(paths["pn_logistic"]),
        "A2": T.load_artifact(paths["pn_gbdt"]),
    }
    arms_sha: dict[str, str] = {
        "A1": T.sha256_file(paths["pn_logistic"]),
        "A2": T.sha256_file(paths["pn_gbdt"]),
    }
    for arm_name, rec in sorted(train_manifest["models"].items()):
        if arm_name.startswith("S3_logo_"):
            continue
        full_p = rdir / rec["rel_path"]
        arms_art[arm_name] = T.load_artifact(full_p)
        arms_sha[arm_name] = recorded_model_shas[rec["rel_path"]]

    # Reference label frequency c from A3's train OOF Elkan-Noto estimate
    c_ref = float(arms_art["A3"]["training"]["pu"]["c"])

    # Accepted negatives from registry + table + optional cohort export + optional model-stratum export
    accepted_neg = PE.load_accepted_negatives(args.negative_registry, table_rows=rows)
    neg_pool_rows = list(rows)
    if args.export_dir is not None:
        exp_receipt = args.export_receipt or DEFAULT_EXPORT_RECEIPT
        extra_neg = _find_export_negative_rows(
            accepted_neg, rows, splits_meta, contract, Path(args.export_dir), Path(exp_receipt)
        )
        neg_pool_rows.extend(extra_neg)
    model_stratum_export = args.work / "model-stratum" / "export"
    if model_stratum_export.exists():
        neg_pool_rows.extend(
            _find_model_stratum_negative_rows(accepted_neg, neg_pool_rows, contract, model_stratum_export)
        )
    joined_neg_rows, unjoined_negatives = PE.join_accepted_negatives(accepted_neg, neg_pool_rows) if hasattr(PE, "join_accepted_negatives") else PE.join_negatives_to_rows(accepted_neg, neg_pool_rows)

    eval_dir.mkdir(parents=True, exist_ok=False)
    scored_test_by_arm: dict[str, list[dict]] = {}
    per_arm_eval: dict[str, dict] = {}

    for arm_name, art in arms_art.items():
        c_arm = float(art.get("training", {}).get("pu", {}).get("c", c_ref))
        thr = float(art["display_threshold"]["value"])

        scored_cal = PE.score(art, cal_rows)
        scored_test = PE.score(art, test_rows)
        scored_neg = PE.score(art, joined_neg_rows)
        scored_test_by_arm[arm_name] = scored_test

        arm_res = {
            "arm": arm_name,
            "model_sha256": arms_sha[arm_name],
            "kind": art["kind"],
            "threshold": thr,
            "c": c_arm,
            "pi": art.get("training", {}).get("pu", {}).get("pi"),
            "pi_source": art.get("training", {}).get("pu", {}).get("pi_source"),
            "calibration_metrics": PE.candidate_metrics(scored_cal, c_arm, thr),
            "test_metrics": PE.candidate_metrics(scored_test, c_arm, thr),
            "deletion_test": PE.deletion_test(art, test_rows),
            "null_layer_calibration": PE.null_layer_metrics(art, rows, "calibration", c_arm),
            "null_layer_test": PE.null_layer_metrics(art, rows, "test", c_arm),
            "negative_specificity": PE.negative_specificity(scored_neg, thr),
            "unjoined_negatives": unjoined_negatives,
        }
        per_arm_eval[arm_name] = arm_res
        _write_json_no_overwrite(eval_dir / f"{arm_name}.json", arm_res)

    # Pairwise contrasts: A3, A4, A5 (and S1/S2) vs A1 and vs A2 on test
    contrasts: dict[str, dict] = {}
    for pu_arm in [k for k in arms_art if k not in ("A1", "A2")]:
        thr_pu = float(arms_art[pu_arm]["display_threshold"]["value"])
        for pn_arm in ("A1", "A2"):
            thr_pn = float(arms_art[pn_arm]["display_threshold"]["value"])
            # Merge joined explicit negatives on test if any
            contrasts[f"{pu_arm}_vs_{pn_arm}_test"] = PE.pairwise_contrast(
                scored_test_by_arm[pu_arm], scored_test_by_arm[pn_arm], thr_pu, thr_pn
            )
            # Also contrast on all joined negative rows
            contrasts[f"{pu_arm}_vs_{pn_arm}_negatives"] = PE.pairwise_contrast(
                PE.score(arms_art[pu_arm], joined_neg_rows),
                PE.score(arms_art[pn_arm], joined_neg_rows),
                thr_pu,
                thr_pn,
            )

    # S3 leave-one-genus-out evaluation on each fold's held-out genus test rows
    s3_results = []
    metric_keys = (
        "supported_precision",
        "pu_expected_precision",
        "recall",
        "auroc_positive_vs_unlabelled",
        "pu_ece",
    )
    diffs_by_metric: dict[str, list[float]] = defaultdict(list)

    for fold in train_manifest.get("s3_logo_folds", []):
        if fold.get("skipped"):
            s3_results.append(fold)
            continue
        genus = fold["genus"]
        held_genomes = set(fold["held_out_genomes"])
        fold_test_rows = [dict(r, split="test") for r in rows if r["genome_id"] in held_genomes]

        pn_rec = train_manifest["models"][fold["pn_arm"]]
        pu_rec = train_manifest["models"][fold["pu_arm"]]
        art_pn = T.load_artifact(rdir / pn_rec["rel_path"])
        art_pu = T.load_artifact(rdir / pu_rec["rel_path"])

        c_fold = float(art_pu["training"]["pu"]["c"])
        thr_pn = float(art_pn["display_threshold"]["value"])
        thr_pu = float(art_pu["display_threshold"]["value"])

        sc_pn = PE.score(art_pn, fold_test_rows)
        sc_pu = PE.score(art_pu, fold_test_rows)
        m_pn = PE.candidate_metrics(sc_pn, c_fold, thr_pn)
        m_pu = PE.candidate_metrics(sc_pu, c_fold, thr_pu)
        del_pn = PE.deletion_test(art_pn, fold_test_rows)
        del_pu = PE.deletion_test(art_pu, fold_test_rows)
        contrast = PE.pairwise_contrast(sc_pu, sc_pn, thr_pu, thr_pn)

        fold_diff = {}
        for mk in metric_keys:
            vp, vn = m_pu.get(mk), m_pn.get(mk)
            if vp is not None and vn is not None:
                d = float(vp - vn)
                fold_diff[mk] = d
                diffs_by_metric[mk].append(d)
            else:
                fold_diff[mk] = None

        fold_out = {
            "genus": genus,
            "held_out_genomes": fold["held_out_genomes"],
            "skipped": False,
            "n_test_rows": len(fold_test_rows),
            "c_fold": c_fold,
            "pi_fold": float(art_pu["training"]["pu"]["pi"]),
            "pn_metrics": m_pn,
            "pu_metrics": m_pu,
            "pu_minus_pn": fold_diff,
            "contrast_pu_vs_pn": contrast,
            "deletion_pn": del_pn,
            "deletion_pu": del_pu,
        }
        s3_results.append(fold_out)
        _write_json_no_overwrite(eval_dir / f"S3_logo_{genus}.json", fold_out)

    s3_summary = {}
    for mk in metric_keys:
        vals = diffs_by_metric.get(mk, [])
        if vals:
            s3_summary[mk] = {
                "n_folds": len(vals),
                "median": float(np.median(vals)),
                "min": float(np.min(vals)),
                "max": float(np.max(vals)),
                "n_positive_diff": int(sum(1 for v in vals if v > 0)),
                "n_zero_diff": int(sum(1 for v in vals if v == 0)),
                "n_negative_diff": int(sum(1 for v in vals if v < 0)),
            }
        else:
            s3_summary[mk] = {"n_folds": 0, "median": None, "min": None, "max": None}

    summary = {
        "regime": args.regime,
        "c_reference": c_ref,
        "arms": per_arm_eval,
        "pairwise_contrasts": contrasts,
        "negatives": {
            "n_accepted_registry_and_table": int(len(accepted_neg)),
            "n_joined_scored": int(len(joined_neg_rows)),
            "n_unjoined": int(len(unjoined_negatives)),
            "unjoined_negatives": unjoined_negatives,
            "per_arm_specificity": {k: v["negative_specificity"] for k, v in per_arm_eval.items()},
        },
        "s1_sweep": train_manifest.get("s1_sweep", []),
        "s2_u2": train_manifest.get("s2_u2", {}),
        "s3_logo": {
            "folds": s3_results,
            "spread_pu_minus_pn": s3_summary,
        },
    }
    _write_json_no_overwrite(eval_path, summary)
    return summary


# ----------------------------------------------------------------------------- Stage 4: report

def _fmt(v, digits: int = 4) -> str:
    if v is None:
        return "—"
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, int):
        return str(v)
    if isinstance(v, float):
        return f"{v:.{digits}f}"
    return str(v)


def _check_criterion_arm(arm_eval: dict, a1_eval: dict) -> dict[str, bool]:
    del_pu = arm_eval["deletion_test"]
    del_a1 = a1_eval["deletion_test"]
    auroc_del = del_pu.get("auroc_p_null")
    if auroc_del is None:
        auroc_del = del_pu.get("auroc_naive_null")
    naive_a1 = del_a1.get("auroc_naive_null")
    p3_ok = (
        auroc_del is not None
        and naive_a1 is not None
        and float(auroc_del) >= 0.70
        and float(auroc_del) > float(naive_a1)
    )

    tm_pu = arm_eval["test_metrics"]
    tm_a1 = a1_eval["test_metrics"]
    p4_ok = (
        tm_pu.get("supported_precision") is not None
        and tm_a1.get("supported_precision") is not None
        and tm_pu.get("recall") is not None
        and tm_a1.get("recall") is not None
        and float(tm_pu["supported_precision"]) >= float(tm_a1["supported_precision"]) - 1e-12
        and float(tm_pu["recall"]) >= float(tm_a1["recall"]) - 1e-12
    )

    neg_pu = arm_eval["negative_specificity"]
    neg_a1 = a1_eval["negative_specificity"]
    p5_ok = int(neg_pu["n_called"]) <= int(neg_a1["n_called"])

    p6_ok = (
        tm_pu.get("pu_ece") is not None
        and tm_a1.get("pu_ece") is not None
        and float(tm_pu["pu_ece"]) <= float(tm_a1["pu_ece"]) + 1e-12
    )
    return {"P3": p3_ok, "P4": p4_ok, "P5": p5_ok, "P6": p6_ok}


def stage_report(args) -> Path:
    rdir = args.work / args.regime
    report_path = rdir / "report_tables.md"
    if report_path.exists():
        raise FileExistsError(f"refusing to overwrite existing report: {report_path}")

    reg = json.loads((rdir / "registration.json").read_text())
    ev = json.loads((rdir / "evaluation.json").read_text())
    arms = ev["arms"]
    main_arms = [a for a in ("A1", "A2", "A3", "A4", "A5") if a in arms]

    lines = [
        f"# PU Learning Evaluation Tables — Regime `{args.regime}`",
        "",
        "## 1. Candidate-Layer Test Metrics (Main Arms)",
        "",
        "| Arm | Kind | Threshold | Supported Precision | PU-Adjusted Precision | Recall | AUROC (Pos vs Unl) | PU-ECE | PN-ECE (Ref) | Calls | Abstained |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for a in main_arms:
        ae = arms[a]
        tm = ae["test_metrics"]
        lines.append(
            f"| {a} | {ae['kind']} | {_fmt(ae['threshold'])} | {_fmt(tm['supported_precision'])} | "
            f"{_fmt(tm['pu_expected_precision'])} | {_fmt(tm['recall'])} | "
            f"{_fmt(tm['auroc_positive_vs_unlabelled'])} | {_fmt(tm['pu_ece'])} | "
            f"{_fmt(tm['pn_ece_reference'])} | {tm['n_calls']} | {tm['abstained']} |"
        )

    lines.extend(
        [
            "",
            "## 2. Pairwise Contrasts on Test Split",
            "",
            "| Contrast (A vs B) | Both Supported | Corrections (A only) | Corrections (B only) | Harms (A only) | Harms (B only) | Unconfirmed (A only) | Unconfirmed (B only) |",
            "|---|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for cname, cval in sorted(ev["pairwise_contrasts"].items()):
        if not cname.endswith("_test"):
            continue
        lines.append(
            f"| {cname} | {cval['n_both_supported']} | {cval['corrections_a_over_b']} | "
            f"{cval['corrections_b_over_a']} | {cval['harms_a_over_b']} | {cval['harms_b_over_a']} | "
            f"{cval['unconfirmed_called_only_by_a']} | {cval['unconfirmed_called_only_by_b']} |"
        )

    lines.extend(
        [
            "",
            "## 3. Experimental Negatives Specificity",
            "",
            f"- Accepted negatives in registry + table: {ev['negatives']['n_accepted_registry_and_table']}",
            f"- Joined and scored rows: {ev['negatives']['n_joined_scored']}",
            f"- Unjoined negatives: {ev['negatives']['n_unjoined']}",
            "",
            "| Arm | Threshold | Supported Negatives | Below Threshold | Called | Specificity | Probabilities |",
            "|---|---:|---:|---:|---:|---:|---|",
        ]
    )
    for a in main_arms:
        ns = arms[a]["negative_specificity"]
        p_str = ", ".join(_fmt(p) for p in ns["probabilities"]) if ns["probabilities"] else "—"
        lines.append(
            f"| {a} | {_fmt(ns['threshold'])} | {ns['n_supported']} | {ns['n_below_threshold']} | "
            f"{ns['n_called']} | {_fmt(ns['specificity'])} | {p_str} |"
        )

    lines.extend(
        [
            "",
            "## 4. Null Layer (Gates and Marginal NLL vs Baselines)",
            "",
            "| Arm | Split | Trained | Groups (>=2) | Single-Pos | No-Label | Marginal NLL | Uniform Marginal NLL | Mean log(n+1) | Mean P(null) Single-Pos | Mean P(null) No-Label | Spearman P(null) vs log1p(n) |",
            "|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for a in main_arms:
        for sp_key, sp_label in (("null_layer_calibration", "calibration"), ("null_layer_test", "test")):
            nl = arms[a][sp_key]
            cc = nl.get("candidate_count_check", {})
            lines.append(
                f"| {a} | {sp_label} | {_fmt(nl['trained'])} | {nl['n_groups_ge2']} | "
                f"{nl['n_single_positive_groups']} | {nl['n_no_label_groups']} | "
                f"{_fmt(nl['marginal_nll'])} | {_fmt(nl['uniform_marginal_nll'])} | "
                f"{_fmt(nl['mean_log_n_plus_1'])} | {_fmt(nl['mean_p_null_single_positive'])} | "
                f"{_fmt(nl['mean_p_null_no_label'])} | {_fmt(cc.get('spearman_p_null_vs_log1p_n'))} |"
            )

    lines.extend(
        [
            "",
            "## 5. Deletion Test (Single-Positive Test Competitions, >=2 Members)",
            "",
            "| Arm | Locus Trained | Competitions | Mean P(null) Intact | Mean P(null) Deleted | AUROC P(null) | Mean Naive Null Intact | Mean Naive Null Deleted | AUROC Naive Null |",
            "|---|---|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for a in main_arms:
        dt = arms[a]["deletion_test"]
        lines.append(
            f"| {a} | {_fmt(dt['locus_trained'])} | {dt['n_competitions']} | "
            f"{_fmt(dt['mean_p_null_intact'])} | {_fmt(dt['mean_p_null_deleted'])} | "
            f"{_fmt(dt['auroc_p_null'])} | {_fmt(dt['mean_naive_null_intact'])} | "
            f"{_fmt(dt['mean_naive_null_deleted'])} | {_fmt(dt['auroc_naive_null'])} |"
        )

    lines.extend(
        [
            "",
            "## 6. S1 Prior Sweep and S2 Unlabelled Pool",
            "",
            "| Arm | Role | pi | Supported Precision | PU-Adjusted Precision | Recall | AUROC (Pos vs Unl) | PU-ECE | Deletion AUROC P(null) | Negatives Called |",
            "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for s1 in ev.get("s1_sweep", []):
        aname = s1["arm"]
        if aname in arms:
            ae = arms[aname]
            tm = ae["test_metrics"]
            dt = ae["deletion_test"]
            ns = ae["negative_specificity"]
            lines.append(
                f"| {aname} | {s1['role']} | {_fmt(s1['pi'])} | {_fmt(tm['supported_precision'])} | "
                f"{_fmt(tm['pu_expected_precision'])} | {_fmt(tm['recall'])} | "
                f"{_fmt(tm['auroc_positive_vs_unlabelled'])} | {_fmt(tm['pu_ece'])} | "
                f"{_fmt(dt['auroc_p_null'])} | {ns['n_called']} |"
            )
    if "S2_U2" in arms:
        ae = arms["S2_U2"]
        tm = ae["test_metrics"]
        dt = ae["deletion_test"]
        ns = ae["negative_specificity"]
        lines.append(
            f"| S2_U2 | export_u2 | {_fmt(ae['pi'])} | {_fmt(tm['supported_precision'])} | "
            f"{_fmt(tm['pu_expected_precision'])} | {_fmt(tm['recall'])} | "
            f"{_fmt(tm['auroc_positive_vs_unlabelled'])} | {_fmt(tm['pu_ece'])} | "
            f"{_fmt(dt['auroc_p_null'])} | {ns['n_called']} |"
        )

    lines.extend(
        [
            "",
            "## 7. S3 Leave-One-Genus-Out (LOGO) Spread (PU minus PN)",
            "",
            "| Metric | Folds | Median (PU - PN) | Min (PU - PN) | Max (PU - PN) | Folds > 0 | Folds = 0 | Folds < 0 |",
            "|---|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for mk, mstat in sorted(ev["s3_logo"]["spread_pu_minus_pn"].items()):
        lines.append(
            f"| {mk} | {mstat['n_folds']} | {_fmt(mstat.get('median'))} | "
            f"{_fmt(mstat.get('min'))} | {_fmt(mstat.get('max'))} | "
            f"{mstat.get('n_positive_diff', 0)} | {mstat.get('n_zero_diff', 0)} | {mstat.get('n_negative_diff', 0)} |"
        )

    # Checklist P0 - P8
    a1 = arms["A1"]
    a3 = arms["A3"]
    c_a3 = _check_criterion_arm(a3, a1)
    nl_cal_a3 = a3["null_layer_calibration"]
    p0_ok = bool(reg.get("plan_sha256_record")) and (rdir / "SHA256SUMS.models").exists()
    p2_ok = (
        bool(nl_cal_a3.get("trained"))
        and nl_cal_a3.get("marginal_nll") is not None
        and nl_cal_a3.get("uniform_marginal_nll") is not None
        and float(nl_cal_a3["marginal_nll"]) < float(nl_cal_a3["uniform_marginal_nll"])
    )

    s1_by_role = {item["role"]: item["arm"] for item in ev.get("s1_sweep", [])}
    arm_020 = arms.get(s1_by_role.get("0.20", ""), a3)
    c_020 = _check_criterion_arm(arm_020, a1)
    p7_ok = (
        c_a3["P3"]
        and c_a3["P4"]
        and c_a3["P5"]
        and c_020["P3"]
        and c_020["P4"]
        and c_020["P5"]
    )

    lines.extend(
        [
            "",
            "## 8. Pass Criteria Checklist (P0–P8)",
            "",
            "| Criterion | Measured Value | Status |",
            "|---|---|---|",
            f"| P0 (Plan & artifact hashes frozen before test) | plan_sha256={bool(reg.get('plan_sha256_record'))}, SHA256SUMS.models=present | {'PASS' if p0_ok else 'FAIL'} |",
            "| P1 (JVM reader & production genome replay parity) | Python artifact schema/calibration verified; JVM replay external | PENDING_JVM_REPLAY |",
            f"| P2 (Null layer trained & cal marginal NLL < uniform) | trained={nl_cal_a3.get('trained')}, cal_nll={_fmt(nl_cal_a3.get('marginal_nll'))} vs uniform={_fmt(nl_cal_a3.get('uniform_marginal_nll'))} | {'PASS' if p2_ok else 'FAIL'} |",
            f"| P3 (Deletion test AUROC >= 0.70 and > A1 naive null) | A3_p_null={_fmt(a3['deletion_test'].get('auroc_p_null'))}, A3_naive={_fmt(a3['deletion_test'].get('auroc_naive_null'))}, A1_naive={_fmt(a1['deletion_test'].get('auroc_naive_null'))} | {'PASS' if c_a3['P3'] else 'FAIL'} |",
            f"| P4 (Supported precision A3 >= A1 at >= recall) | A3_prec={_fmt(a3['test_metrics']['supported_precision'])}, A1_prec={_fmt(a1['test_metrics']['supported_precision'])}, A3_rec={_fmt(a3['test_metrics']['recall'])}, A1_rec={_fmt(a1['test_metrics']['recall'])} | {'PASS' if c_a3['P4'] else 'FAIL'} |",
            f"| P5 (Harms on negatives A3 <= A1) | A3_called={a3['negative_specificity']['n_called']}, A1_called={a1['negative_specificity']['n_called']} | {'PASS' if c_a3['P5'] else 'FAIL'} |",
            f"| P6 (Test PU-ECE A3 <= A1) | A3_pu_ece={_fmt(a3['test_metrics']['pu_ece'])}, A1_pu_ece={_fmt(a1['test_metrics']['pu_ece'])} | {'PASS' if c_a3['P6'] else 'FAIL'} |",
            f"| P7 (Robustness at pi_hat and pi=0.20 + LOGO reported) | pi_hat(P3={c_a3['P3']},P4={c_a3['P4']},P5={c_a3['P5']}), pi_0.20(P3={c_020['P3']},P4={c_020['P4']},P5={c_020['P5']}) | {'PASS' if p7_ok else 'FAIL'} |",
            "| P8 (Price-148 evaluation under frozen Price split) | External to cohort table (awaits frozen Price Step C claims) | PENDING_PRICE |",
            "",
        ]
    )
    report_path.write_text("\n".join(lines))
    return report_path


# ----------------------------------------------------------------------------- CLI

def _add_common_args(parser: argparse.ArgumentParser, suppress_defaults: bool = False) -> None:
    def d(val):
        return argparse.SUPPRESS if suppress_defaults else val

    parser.add_argument("--work", type=Path, default=d(DEFAULT_WORK))
    parser.add_argument("--regime", choices=("family_excluded", "homology"), default=d("family_excluded"))
    parser.add_argument("--contract", type=Path, default=d(None))
    parser.add_argument("--restore-receipt", type=Path, default=d(None))
    parser.add_argument("--splits", type=Path, default=d(DEFAULT_SPLITS))
    parser.add_argument("--negative-registry", type=Path, default=d(DEFAULT_NEG_REGISTRY))
    parser.add_argument("--export-dir", type=Path, default=d(None))
    parser.add_argument("--export-receipt", type=Path, default=d(None))
    if suppress_defaults:
        parser.add_argument("--synthetic", action="store_true", default=argparse.SUPPRESS)
        parser.add_argument("--skip-s3", action="store_true", default=argparse.SUPPRESS)
    else:
        parser.add_argument("--synthetic", action="store_true", default=False)
        parser.add_argument("--skip-s3", action="store_true", default=False)
    parser.add_argument("--folds", type=int, default=d(5))
    parser.add_argument("--seed", type=int, default=d(42))
    parser.add_argument("--min-no-label-groups", type=int, default=d(20))
    parser.add_argument("--min-pos-groups", type=int, default=d(20))
    parser.add_argument("--min-cal-no-label-groups", type=int, default=d(5))
    parser.add_argument("--min-cal-pos-groups", type=int, default=d(5))


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__)
    _add_common_args(ap, suppress_defaults=False)
    sub = ap.add_subparsers(dest="stage", required=True)
    for st in ("prepare", "train", "evaluate", "report"):
        sp = sub.add_parser(st)
        _add_common_args(sp, suppress_defaults=True)
    return ap


def main(argv=None) -> int:
    # Support common flags either before or after the stage positional argument.
    ap = build_parser()
    args = ap.parse_args(argv)
    try:
        if args.stage == "prepare":
            stage_prepare(args)
        elif args.stage == "train":
            stage_train(args)
        elif args.stage == "evaluate":
            stage_evaluate(args)
        elif args.stage == "report":
            stage_report(args)
    except (T.ProvenanceError, FileExistsError) as e:
        print(f"REFUSED: {e}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
