#!/usr/bin/env python3
"""Price-148 evaluation glue for the positive-unlabelled (PU) learning plan (row 4, §4.2, P8).

Reuses ``benchmark/v3/price_contract_learned/evaluate.py`` (``collect``, ``read_raw``, ``learned_scored_only``,
``at_threshold``, ``learned_counts``) and ``benchmark/v3/price_calibration_context/row1_calibration.py``
(``pairs``, ``choose_threshold``, ``evaluate_arm``, ``calibration_metrics``, ``contrast``, ``load``, ``sha``,
``save``) by importing them directly without copying.

Order (plan E1): thresholds are chosen on the calibration genomes of the frozen row-1 split (primary and
swapped) and written to ``policy.json`` and hashed before any test-genome metric is computed. Unconfirmed
calls are never negatives; statistics that need a binary outcome treat them as unsupported and are labelled
reference-support statistics. Arms missing one or more cells are reported as missing, never scored as zero.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

if not hasattr(hashlib, "file_digest"):
    def _file_digest(fileobj, digest):
        h = hashlib.new(digest) if isinstance(digest, str) else digest()
        while chunk := fileobj.read(1 << 18):
            h.update(chunk)
        return h

    hashlib.file_digest = _file_digest  # type: ignore[attr-defined]

HERE = Path(__file__).resolve().parent
V3 = HERE.parent
REPO = V3.parent.parent
sys.path.insert(0, str(V3))
import ec_identity as E  # noqa: E402
from price_calibration_context import row1_calibration as R1  # noqa: E402
from price_contract_learned import evaluate as EV  # noqa: E402

SCHEMA = "gspa.pu-price-evaluation/1"
BASE_ARMS = ("raw", "stage0", "stage1", "learned_pn", "learned_pn_gbdt", "learned_pu", "learned_pu_gbdt")
RANK_ARM = "learned_pu_rank"
LEARNED_ARMS = ("learned_pn", "learned_pn_gbdt", "learned_pu", "learned_pu_gbdt", "learned_pu_rank")
PU_LEARNED_ARMS = ("learned_pu", "learned_pu_gbdt", "learned_pu_rank")
CONTRASTS = (
    ("stage1", "learned_pu"),
    ("learned_pn", "learned_pu"),
    ("learned_pn_gbdt", "learned_pu_gbdt"),
    ("learned_pu", "learned_pu_rank"),
)
SENSITIVITY_GENOME = "GCA_000006965.1"

DEFAULT_MODEL_SEARCH_DIRS = (
    Path("/mnt/data2/leechuck/gspa-v3-data/pu-learning-v1/homology/models"),
    Path("/mnt/data2/leechuck/gspa-v3-data/pu-learning-v1/family_excluded/models"),
    Path("/mnt/data2/leechuck/gspa-v3-data/pu-learning-v1/gate-deviation/homology/models"),
    Path("/mnt/data2/leechuck/gspa-v3-data/pu-learning-v1/gate-deviation/family_excluded/models"),
    Path("/mnt/data2/leechuck/gspa-v3-data/pu-learning-v1/inputs/family-excluded-default/homology/models"),
    Path("/mnt/data2/leechuck/gspa-v3-data/pu-learning-v1/inputs/family-excluded-default/models"),
)


def _verify_and_record_cell_hashes(cells_dir: Path) -> dict[str, dict]:
    """Verify every cell directory in ``cells_dir`` has matching receipt.json and retained.json hashes."""
    cells_path = Path(cells_dir)
    if not cells_path.is_dir():
        raise ValueError(f"cells directory does not exist: {cells_path}")
    cell_hashes: dict[str, dict] = {}
    for d in sorted(cells_path.iterdir(), key=lambda p: int(p.name) if p.name.isdigit() else -1):
        if not d.is_dir():
            continue
        rec_file = d / "receipt.json"
        ret_file = d / "retained.json"
        if not rec_file.is_file() and not ret_file.is_file():
            continue
        if not rec_file.is_file():
            raise ValueError(f"{d}: missing receipt.json")
        if not ret_file.is_file():
            raise ValueError(f"{d}: missing retained.json")
        rec = R1.load(rec_file)
        actual_ret_sha = R1.sha(ret_file)
        if actual_ret_sha != rec.get("retained_sha256"):
            raise ValueError(f"{d}: retained file changed (expected {rec.get('retained_sha256')}, got {actual_ret_sha})")
        assembly = rec.get("cell", {}).get("assembly", d.name)
        cell_hashes[str(assembly)] = {
            "cell_dir": str(d),
            "index": rec.get("cell", {}).get("index", d.name),
            "receipt_sha256": R1.sha(rec_file),
            "retained_sha256": actual_ret_sha,
        }
    return cell_hashes


def _count_claims(claims_dir: Path, genomes: list[str], pn_receipts: dict) -> tuple[dict, dict[str, str]]:
    """Count GO rows, complete four-level ECs, and partial ECs in claims_dir and record claims file hashes."""
    totals = Counter()
    per_genome = {}
    claims_hashes: dict[str, str] = {}
    for g in sorted(genomes):
        cfile = Path(claims_dir) / f"{g}.claims.jsonl"
        if not cfile.is_file():
            continue
        csha = R1.sha(cfile)
        claims_hashes[g] = csha
        exp_sha = pn_receipts.get(g, {}).get("cell", {}).get("claims_sha256")
        if exp_sha and exp_sha != csha:
            raise ValueError(f"{cfile}: claims sha256 mismatch against receipt (expected {exp_sha}, got {csha})")
        g_counts = Counter()
        with cfile.open() as fh:
            for line in fh:
                if not line.strip():
                    continue
                r = json.loads(line)
                fid = str(r.get("function_id", ""))
                g_counts["total_rows"] += 1
                if fid.startswith("GO:"):
                    g_counts["go_rows"] += 1
                elif EV.complete(fid):
                    g_counts["complete_ec_rows"] += 1
                else:
                    g_counts["partial_ec_rows"] += 1
        per_genome[g] = dict(g_counts)
        totals.update(g_counts)
    return {"totals": dict(totals), "per_genome": per_genome}, claims_hashes


def _remap_arm_collection(
    raw_data: dict,
    raw_learned: dict,
    raw_receipts: dict,
    raw_dropped: Counter,
    arm_map: dict[str, tuple[str, ...]],
) -> tuple[dict, dict, dict, Counter]:
    """Map source arm names in a cell collection (e.g. 'learned' -> 'learned_pu') into canonical arm names."""
    data_out: dict[str, dict] = defaultdict(dict)
    learned_out: dict[str, dict] = defaultdict(dict)
    receipts_out: dict[str, dict] = {}
    dropped_out: Counter = Counter()

    for g, rec in raw_receipts.items():
        rcopy = copy.deepcopy(rec)
        rcopy.setdefault("arms", {})
        for target_arm, src_candidates in arm_map.items():
            for src_arm in src_candidates:
                if src_arm in rec.get("arms", {}):
                    rcopy["arms"][target_arm] = copy.deepcopy(rec["arms"][src_arm])
                    break
        receipts_out[g] = rcopy

    for target_arm, src_candidates in arm_map.items():
        for src_arm in src_candidates:
            if src_arm in raw_data and raw_data[src_arm]:
                data_out[target_arm] = raw_data[src_arm]
                dropped_out[target_arm] += raw_dropped.get(src_arm, 0)
                break
        for src_arm in src_candidates:
            if src_arm in raw_learned and raw_learned[src_arm]:
                learned_out[target_arm] = raw_learned[src_arm]
                break

    return data_out, learned_out, receipts_out, dropped_out


def _index_models_by_sha(search_dirs: tuple[Path, ...] = DEFAULT_MODEL_SEARCH_DIRS) -> dict[str, tuple[Path, float]]:
    """Scan known model directories and map model file sha256 -> (path, display_threshold)."""
    by_sha: dict[str, tuple[Path, float]] = {}
    for sdir in search_dirs:
        if not sdir.is_dir():
            continue
        for mpath in sorted(sdir.rglob("model_*.json")):
            if not mpath.is_file():
                continue
            try:
                sha = R1.sha(mpath)
                obj = R1.load(mpath)
                dt = obj.get("display_threshold", {})
                val = dt.get("value") if isinstance(dt, dict) else dt
                if val is not None:
                    by_sha[sha] = (mpath, float(val))
            except Exception:  # noqa: BLE001
                continue
    return by_sha


def _extract_threshold_from_model_file(mpath: Path) -> tuple[float, str] | None:
    p = Path(mpath)
    if not p.is_file():
        return None
    obj = R1.load(p)
    dt = obj.get("display_threshold", {})
    val = dt.get("value") if isinstance(dt, dict) else dt
    if val is None:
        return None
    return float(val), R1.sha(p)


def _resolve_display_thresholds(
    learned_arms: list[str],
    receipts_by_arm: dict[str, dict[str, dict]],
    explicit_models: dict[str, Path],
    pn_cells_dir: Path,
) -> tuple[dict[str, float], dict[str, dict]]:
    """Resolve each learned arm's artifact display_threshold and provenance."""
    sha_index: dict[str, tuple[Path, float]] | None = None
    pn_policy_thresholds: dict[str, float] = {}
    pn_policy_path = Path(pn_cells_dir).parent / "evaluation" / "policy.json"
    if pn_policy_path.is_file():
        try:
            pol_obj = R1.load(pn_policy_path)
            vals = pol_obj.get("model_display_threshold", {}).get("values", {})
            if "learned" in vals:
                pn_policy_thresholds["learned_pn"] = float(vals["learned"])
            if "learned_gbdt" in vals:
                pn_policy_thresholds["learned_pn_gbdt"] = float(vals["learned_gbdt"])
        except Exception:  # noqa: BLE001
            pass

    values: dict[str, float] = {}
    sources: dict[str, dict] = {}

    for arm in learned_arms:
        # 1. Explicit model file CLI argument
        if arm in explicit_models and explicit_models[arm] is not None:
            res = _extract_threshold_from_model_file(explicit_models[arm])
            if res is None:
                raise ValueError(f"cannot read display_threshold from model file for {arm}: {explicit_models[arm]}")
            val, msha = res
            values[arm] = val
            sources[arm] = {"value": val, "source": str(explicit_models[arm]), "model_sha256": msha}
            continue

        # 2. Inspect receipt for the arm
        arm_recs = receipts_by_arm.get(arm, {})
        found = False
        for g in sorted(arm_recs):
            arec = arm_recs[g].get("arms", {}).get(arm, {})
            lrec = arec.get("learned") or {}
            mrec = lrec.get("model") or {}

            # 2a. Direct display_threshold in receipt learned/model block
            for cand_dt in (lrec.get("display_threshold"), mrec.get("display_threshold")):
                if cand_dt is not None:
                    val = float(cand_dt["value"] if isinstance(cand_dt, dict) else cand_dt)
                    values[arm] = val
                    sources[arm] = {"value": val, "source": "cell_receipt", "model_sha256": mrec.get("sha256")}
                    found = True
                    break
            if found:
                break

            # 2b. Command --learned-model path if accessible on local filesystem
            cmd = arec.get("command") or []
            if "--learned-model" in cmd:
                idx = cmd.index("--learned-model")
                if idx + 1 < len(cmd):
                    cand_path = Path(cmd[idx + 1])
                    res = _extract_threshold_from_model_file(cand_path)
                    if res is not None:
                        val, msha = res
                        values[arm] = val
                        sources[arm] = {"value": val, "source": str(cand_path), "model_sha256": msha}
                        found = True
                        break

            # 2c. Match model sha256 against local model directories
            msha = mrec.get("sha256")
            if msha and msha != "0" * 64:
                if sha_index is None:
                    sha_index = _index_models_by_sha()
                if msha in sha_index:
                    mpath, val = sha_index[msha]
                    values[arm] = val
                    sources[arm] = {"value": val, "source": str(mpath), "model_sha256": msha}
                    found = True
                    break

        if found:
            continue

        # 3. Fallback to frozen PN evaluation policy.json for PN arms
        if arm in pn_policy_thresholds:
            val = pn_policy_thresholds[arm]
            values[arm] = val
            sources[arm] = {"value": val, "source": str(pn_policy_path), "model_sha256": None}
            continue

        # 4. Fallback to 0.5
        values[arm] = 0.5
        sources[arm] = {"value": 0.5, "source": "default_0.5_fallback", "model_sha256": None}

    return values, sources


def _check_arm_completeness(
    arm_names: list[str],
    data: dict[str, dict],
    expected_genomes: list[str],
    reference_by_genome: dict[str, set[str]],
) -> tuple[list[str], dict[str, dict], dict[str, dict]]:
    """Separate complete arms from arms missing any cell or reference protein."""
    complete_arms: list[str] = []
    missing_arms: dict[str, dict] = {}
    arm_status: dict[str, dict] = {}
    exp_set = set(expected_genomes)

    for arm in arm_names:
        by_g = data.get(arm, {})
        valid_genomes: set[str] = set()
        incomplete_proteins: dict[str, list[str]] = {}
        for g, prot in by_g.items():
            exp_pids = reference_by_genome.get(g, set())
            missing_pids = sorted(exp_pids - set(prot.keys()))
            if missing_pids:
                incomplete_proteins[g] = missing_pids
            else:
                valid_genomes.add(g)

        missing_gs = sorted(exp_set - valid_genomes)
        present_gs = sorted(valid_genomes & exp_set)
        if not missing_gs and not incomplete_proteins and len(present_gs) == len(exp_set):
            complete_arms.append(arm)
            arm_status[arm] = {
                "status": "complete",
                "present_genomes": present_gs,
                "missing_genomes": [],
                "n_present": len(present_gs),
                "n_expected": len(exp_set),
            }
        else:
            rec = {
                "status": "missing",
                "reason": "missing_cells",
                "present_genomes": present_gs,
                "missing_genomes": missing_gs,
                "incomplete_reference_proteins": incomplete_proteins,
                "n_present": len(present_gs),
                "n_expected": len(exp_set),
            }
            missing_arms[arm] = rec
            arm_status[arm] = rec

    return complete_arms, missing_arms, arm_status


def _fmt(x, nd: int = 3) -> str:
    if x is None:
        return "—"
    if isinstance(x, float):
        return f"{x:.{nd}f}"
    return str(x)


def _ci(b, nd: int = 3) -> str:
    if not b or b.get("estimate") is None:
        return "—"
    return f"{_fmt(b['estimate'], nd)} [{_fmt(b.get('ci_low'), nd)}, {_fmt(b.get('ci_high'), nd)}]"


def write_markdown_tables(report: dict, out_dir: Path) -> Path:
    """Write markdown tables (arms x metrics; contrasts; learned counts) to out_dir / 'report_tables.md'."""
    md_path = Path(out_dir) / "report_tables.md"
    if md_path.exists():
        raise FileExistsError(f"refusing to overwrite existing file: {md_path}")

    lines = [
        "# Price-148 PU Evaluation Tables (Row 4, §4.2, P8)",
        "",
        f"- Convention: {report['convention']}",
        f"- Policy sha256: `{report['policy_sha256']}`",
        f"- Split sha256: `{report['split_sha256']}`",
        "",
        "## 1. Arms × Test Metrics (Calibrated WP5 Policy, Fixed 0.5, and Artifact Display Threshold)",
        "",
        "| Analysis | Arm | Status | Proteins | Cal Threshold | Cal Calls | Cal Supported | Cal Unconfirmed | Support Share [95% CI] | Ref Recall [95% CI] | Top-1 Match | Top-1 Share [95% CI] | Top-3 Match | Display Thr | Display Calls (Sup / Unconf) | Pos-Pair Brier | Top-Call ECE | LS Top-Call ECE |",
        "|---|---|---|---:|---:|---:|---:|---:|---|---|---:|---|---:|---:|---|---:|---:|---:|",
    ]

    all_listed_arms = list(report["arm_status"].keys())
    for key, an in report["analyses"].items():
        for arm in all_listed_arms:
            r = an["arms"].get(arm)
            if r is None or r.get("status") == "missing":
                m_gs = ",".join((r or {}).get("missing_genomes", [])) or "all"
                lines.append(
                    f"| {key} | {arm} | missing ({m_gs}) | — | — | — | — | — | — | — | — | — | — | — | — | — | — | — |"
                )
                continue
            cc = r.get("calibrated_calls") or {}
            cm = r.get("calibration_metrics") or {}
            cd = r.get("at_model_display_threshold") or {}
            ls = r.get("learned_scored_only") or {}
            lscm = ls.get("calibration_metrics") or {}
            disp_str = (
                f"{cd.get('calls')} ({cd.get('supported')} / {cd.get('unconfirmed')})"
                if cd.get("calls") is not None
                else "—"
            )
            lines.append(
                f"| {key} | {arm} | {r['status']} | {r['proteins']} | {_fmt(r['threshold'], 6)} | "
                f"{_fmt(cc.get('calls'))} | {_fmt(cc.get('supported'))} | {_fmt(cc.get('unconfirmed'))} | "
                f"{_ci(cc.get('support_share'))} | {_ci(cc.get('reference_recall'))} | "
                f"{r['top1_match']} | {_ci(r['top1']['share'])} | {r['top3_match']} | "
                f"{_fmt(cd.get('threshold'), 4)} | {disp_str} | "
                f"{_fmt(cm.get('positive_pair_brier'), 4)} | {_fmt(cm.get('top_call_ece'), 4)} | "
                f"{_fmt(lscm.get('top_call_ece'), 4)} |"
            )

    lines.extend(
        [
            "",
            "## 2. Pairwise Contrasts on Test Genomes",
            "",
            "| Analysis | Contrast | Thr Before | Thr After | Cal Corrected | Cal Harmed | Fixed 0.5 Corrected | Fixed 0.5 Harmed | Fixed 0.5 Net [95% CI] | Display Corrected | Display Harmed | Top-1 Gained | Top-1 Lost | Top-1 Net [95% CI] | Changed Top-1 |",
            "|---|---|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|---:|---|---:|",
        ]
    )
    for key, an in report["analyses"].items():
        for cname, c in an["contrasts"].items():
            if c.get("status") == "missing_arm":
                lines.append(
                    f"| {key} | {cname} | — | — | — | — | — | — | — | — | — | — | — | — | — |"
                )
                continue
            disp_c = c.get("at_model_display_threshold", {}).get("calibrated", {})
            lines.append(
                f"| {key} | {cname} | {_fmt(c['thresholds']['before'], 6)} | {_fmt(c['thresholds']['after'], 6)} | "
                f"{c['calibrated']['corrected']} | {c['calibrated']['harmed']} | "
                f"{c['fixed_0_5']['corrected']} | {c['fixed_0_5']['harmed']} | {_ci(c['fixed_0_5']['net'], 2)} | "
                f"{_fmt(disp_c.get('corrected'))} | {_fmt(disp_c.get('harmed'))} | "
                f"{c['top1']['gained']} | {c['top1']['lost']} | {_ci(c['top1']['net'], 2)} | "
                f"{c['paired']['changed_top1']} |"
            )

    lines.extend(
        [
            "",
            "## 3. Learned Arm Candidate & Abstention Counts (All Genomes and Without GCA_000006965.1)",
            "",
            "| Scope | Arm | Cells | Executed | Candidates | Scored | Unsupported | Rescored | Contract Unavail | Ref EC Scored | Ref EC Unsupported | Unsupported Reasons |",
            "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|",
        ]
    )
    for scope_label, scope_key in (("all", "learned"), ("without_GCA_000006965.1", "learned_without_GCA_000006965.1")):
        lmap = report.get(scope_key, {})
        for arm in all_listed_arms:
            if arm not in LEARNED_ARMS:
                continue
            if arm not in lmap or lmap[arm].get("status") == "missing":
                lines.append(f"| {scope_label} | {arm} | — | — | — | — | — | — | — | — | — | missing |")
                continue
            lc = lmap[arm]
            ref_ec = lc.get("reference_ec", {})
            ref_scored = ref_ec.get("scored", 0)
            ref_unsup = sum(v for k, v in ref_ec.items() if k != "scored")
            reasons_str = ", ".join(f"{k}:{v}" for k, v in sorted(lc.get("unsupported_reasons", {}).items())) or "none"
            lines.append(
                f"| {scope_label} | {arm} | {lc['cells']} | {lc['executed']} | {lc['candidates']} | "
                f"{lc['scored']} | {lc['unsupported']} | {lc['rescored']} | {lc['contract_unavailable']} | "
                f"{ref_scored} | {ref_unsup} | {reasons_str} |"
            )

    md_path.write_text("\n".join(lines) + "\n")
    return md_path


def run(args) -> dict:
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=False)

    ec_audit_obj = R1.load(args.ec_audit)
    entries = ec_audit_obj["entries"]
    E.validate_transfer_graph(entries)

    # 1. Verify every retained.json against its receipt.json and record all hashes
    pn_cell_hashes = _verify_and_record_cell_hashes(args.pn_cells)
    pu_cell_hashes = _verify_and_record_cell_hashes(args.pu_cells)
    rank_cell_hashes = _verify_and_record_cell_hashes(args.pu_rank_cells) if args.pu_rank_cells else {}

    pn_raw_data, pn_raw_learned, pn_raw_receipts, pn_raw_dropped = EV.collect(args.pn_cells, args.claims_dir)
    pu_raw_data, pu_raw_learned, pu_raw_receipts, pu_raw_dropped = EV.collect(args.pu_cells, args.claims_dir)
    if args.pu_rank_cells:
        rk_raw_data, rk_raw_learned, rk_raw_receipts, rk_raw_dropped = EV.collect(args.pu_rank_cells, args.claims_dir)
    else:
        rk_raw_data, rk_raw_learned, rk_raw_receipts, rk_raw_dropped = {}, {}, {}, Counter()

    pn_map = {
        "raw": ("raw",),
        "stage0": ("stage0",),
        "stage1": ("stage1",),
        "learned_pn": ("learned", "learned_pn"),
        "learned_pn_gbdt": ("learned_gbdt", "learned_pn_gbdt"),
    }
    pu_map = {
        "learned_pu": ("learned", "learned_pu"),
        "learned_pu_gbdt": ("learned_gbdt", "learned_pu_gbdt"),
    }
    rk_map = {
        "learned_pu_rank": ("learned", "learned_pu_rank"),
    }

    pn_data, pn_learned, pn_receipts, pn_dropped = _remap_arm_collection(
        pn_raw_data, pn_raw_learned, pn_raw_receipts, pn_raw_dropped, pn_map
    )
    pu_data, pu_learned, pu_receipts, pu_dropped = _remap_arm_collection(
        pu_raw_data, pu_raw_learned, pu_raw_receipts, pu_raw_dropped, pu_map
    )
    rk_data, rk_learned, rk_receipts, rk_dropped = _remap_arm_collection(
        rk_raw_data, rk_raw_learned, rk_raw_receipts, rk_raw_dropped, rk_map
    )

    data: dict[str, dict] = defaultdict(dict)
    learned: dict[str, dict] = defaultdict(dict)
    receipts_by_arm: dict[str, dict[str, dict]] = {}
    dropped = Counter()

    for arm in ("raw", "stage0", "stage1", "learned_pn", "learned_pn_gbdt"):
        if arm in pn_data:
            data[arm] = pn_data[arm]
        if arm in pn_learned:
            learned[arm] = pn_learned[arm]
        receipts_by_arm[arm] = pn_receipts
        dropped[arm] = pn_dropped.get(arm, 0)

    for arm in ("learned_pu", "learned_pu_gbdt"):
        if arm in pu_data:
            data[arm] = pu_data[arm]
        if arm in pu_learned:
            learned[arm] = pu_learned[arm]
        receipts_by_arm[arm] = pu_receipts
        dropped[arm] = pu_dropped.get(arm, 0)

    if args.pu_rank_cells:
        if "learned_pu_rank" in rk_data:
            data["learned_pu_rank"] = rk_data["learned_pu_rank"]
        if "learned_pu_rank" in rk_learned:
            learned["learned_pu_rank"] = rk_learned["learned_pu_rank"]
        receipts_by_arm["learned_pu_rank"] = rk_receipts
        dropped["learned_pu_rank"] = rk_dropped.get("learned_pu_rank", 0)

    split = R1.load(args.split)
    split_sha = R1.sha(args.split)
    ec_audit_sha = R1.sha(args.ec_audit)

    expected_genomes = sorted(set(pn_receipts.keys()) | set(split.get("assignment", {}).keys()))
    if not expected_genomes:
        raise ValueError("no genomes found in pn-cells or split")

    reference_by_genome = {
        g: set(data["stage0"].get(g, {}).keys()) or set(data["raw"].get(g, {}).keys())
        for g in expected_genomes
    }
    pop = {(g, pid) for g, pids in reference_by_genome.items() for pid in pids}

    requested_arms = list(BASE_ARMS) + ([RANK_ARM] if args.pu_rank_cells else [])
    requested_learned_arms = [a for a in requested_arms if a in LEARNED_ARMS]

    arms, missing_arms, arm_status = _check_arm_completeness(
        requested_arms, data, expected_genomes, reference_by_genome
    )

    claims_summary, claims_hashes = _count_claims(args.claims_dir, expected_genomes, pn_receipts)

    explicit_models = {
        "learned_pn": args.pn_model,
        "learned_pn_gbdt": args.pn_model_gbdt,
        "learned_pu": args.pu_model,
        "learned_pu_gbdt": args.pu_model_gbdt,
        "learned_pu_rank": args.pu_rank_model,
    }
    display_values, display_sources = _resolve_display_thresholds(
        requested_learned_arms, receipts_by_arm, explicit_models, args.pn_cells
    )

    cal = {g for g, a in split["assignment"].items() if a == "calibration"}
    test = {g for g, a in split["assignment"].items() if a == "test"}

    # ---- E1: thresholds on calibration genomes only (primary and swapped), written and hashed BEFORE any test metric
    policies: dict[str, dict] = {}
    for name, cal_set in (("primary", cal), ("swapped", test)):
        pol: dict[str, dict] = {}
        for arm in arms:
            rows = []
            for g, prot in EV.subset(data[arm], cal_set).items():
                rows += R1.pairs(prot, entries)[0]
            pol[arm] = R1.choose_threshold(rows)
        for m_arm, m_info in missing_arms.items():
            pol[m_arm] = {
                "threshold": None,
                "status": "missing",
                "missing_genomes": m_info["missing_genomes"],
            }
        policies[name] = pol

    fixed_display_policies = {
        arm: {
            "threshold": display_values[arm],
            "status": "fixed",
            "source": display_sources[arm]["source"],
            "model_sha256": display_sources[arm]["model_sha256"],
        }
        for arm in requested_learned_arms
    }

    policy = {
        "schema": SCHEMA + ":policy",
        "target_support": R1.PRECISION,
        "min_calls": R1.MIN_DECISIONS,
        "selection_rule": R1.SELECTION_RULE,
        "split_sha256": split_sha,
        "calibration_genomes": {"primary": sorted(cal), "swapped": sorted(test)},
        "arms": policies,
        "missing_arms": missing_arms,
        "model_display_threshold": {
            "values": display_values,
            "fixed_policies": fixed_display_policies,
            "sources": display_sources,
            "note": "artifact display_threshold chosen on its cohort training calibration split (fixed policy)",
        },
    }
    policy_sha = R1.save(out / "policy.json", policy)

    # ---- E2, E3, E4 & sensitivity without GCA_000006965.1 (computed strictly AFTER policy.json is on disk)
    assert (out / "policy.json").is_file() and R1.sha(out / "policy.json") == policy_sha

    populations = {
        "all": None,
        f"without_{SENSITIVITY_GENOME}": {SENSITIVITY_GENOME},
    }

    def drop_pop(by_genome: dict, popname: str) -> dict:
        rule = populations[popname]
        if rule is None:
            return by_genome
        return {g: p for g, p in by_genome.items() if g not in rule}

    recorded_hashes = {
        "split_sha256": split_sha,
        "ec_audit_sha256": ec_audit_sha,
        "policy_sha256": policy_sha,
        "pn_cells": pn_cell_hashes,
        "pu_cells": pu_cell_hashes,
        "pu_rank_cells": rank_cell_hashes,
        "claims_files": claims_hashes,
        "model_display_threshold_sources": display_sources,
    }

    report = {
        "schema": SCHEMA,
        "plan": "tasks/v3-row4-pu-learning-2026-11-02.md",
        "policy_sha256": policy_sha,
        "split_sha256": split_sha,
        "ec_audit_sha256": ec_audit_sha,
        "hashes": recorded_hashes,
        "arms": arms,
        "missing_arms": missing_arms,
        "arm_status": arm_status,
        "dropped_non_complete_ec_rows": dict(dropped),
        "claims_counts": claims_summary,
        "convention": (
            "unconfirmed calls are not negatives; support share counts them as unsupported (reference support)"
        ),
        "analyses": {},
        "all_genomes": {},
        "all_genomes_without_GCA_000006965.1": {},
        "learned": {},
        "learned_without_GCA_000006965.1": {},
    }

    for popname in populations:
        for name, test_set, cal_set in (("primary", test, cal), ("swapped", cal, test)):
            key = f"{name}/{popname}"
            effective_test = sorted(g for g in test_set if not (populations[popname] and g in populations[popname]))
            an: dict[str, dict] = {
                "test_genomes": effective_test,
                "calibration_genomes": sorted(cal_set),
                "arms": {},
                "contrasts": {},
            }
            for arm in requested_arms:
                if arm in missing_arms:
                    an["arms"][arm] = {
                        "status": "missing",
                        "missing_genomes": missing_arms[arm]["missing_genomes"],
                        "present_genomes": missing_arms[arm]["present_genomes"],
                        "threshold": None,
                        "genomes": None,
                        "proteins": None,
                        "top1": None,
                        "top1_match": None,
                        "top3_match": None,
                        "fixed_0_5": None,
                        "calibrated_calls": None,
                        "calibration_metrics": None,
                        "at_0_5": None,
                    }
                    continue

                by = drop_pop(EV.subset(data[arm], test_set), popname)
                pol = policies[name][arm]
                res, _, _ = R1.evaluate_arm(by, entries, pol)
                res["calibration_metrics"] = R1.calibration_metrics(by, entries)
                res["at_0_5"] = EV.at_threshold(by, entries, 0.5)
                if arm in LEARNED_ARMS:
                    display_t = display_values[arm]
                    at_disp = EV.at_threshold(by, entries, display_t) or {}
                    res["at_model_display_threshold"] = {"threshold": display_t, **at_disp}
                    ls = EV.learned_scored_only(by, learned.get(arm, {}))
                    ls_disp = EV.at_threshold(ls, entries, display_t) or {}
                    res["learned_scored_only"] = {
                        "calibration_metrics": R1.calibration_metrics(ls, entries),
                        "at_0_5": EV.at_threshold(ls, entries, 0.5),
                        "at_model_display_threshold": {"threshold": display_t, **ls_disp},
                    }
                an["arms"][arm] = res

            for before, after in CONTRASTS:
                if before not in requested_arms and after not in requested_arms:
                    continue
                if before in missing_arms or after in missing_arms or before not in arms or after not in arms:
                    missing_pair = [a for a in (before, after) if a not in arms]
                    an["contrasts"][f"{before}_to_{after}"] = {
                        "status": "missing_arm",
                        "missing_arms": missing_pair,
                    }
                    continue
                by_b = drop_pop(EV.subset(data[before], test_set), popname)
                by_a = drop_pop(EV.subset(data[after], test_set), popname)
                c_res = R1.contrast(
                    by_b,
                    by_a,
                    entries,
                    policies[name][before]["threshold"],
                    policies[name][after]["threshold"],
                )
                tb_disp = display_values.get(before, 0.5)
                ta_disp = display_values.get(after, 0.5)
                c_disp = R1.contrast(by_b, by_a, entries, tb_disp, ta_disp)
                c_res["at_model_display_threshold"] = {
                    "thresholds": {"before": tb_disp, "after": ta_disp},
                    "calibrated": c_disp["calibrated"],
                }
                an["contrasts"][f"{before}_to_{after}"] = c_res

            report["analyses"][key] = an

    # Descriptive all-genomes summaries (with and without GCA_000006965.1)
    for popname, target_key in (("all", "all_genomes"), (f"without_{SENSITIVITY_GENOME}", "all_genomes_without_GCA_000006965.1")):
        for arm in requested_arms:
            if arm in missing_arms:
                report[target_key][arm] = {
                    "status": "missing",
                    "missing_genomes": missing_arms[arm]["missing_genomes"],
                }
                continue
            by = drop_pop(data[arm], popname)
            if not by:
                continue
            pol = {"threshold": None, "status": "not_calibrated"}
            m, _, _ = R1.evaluate_arm(by, entries, pol)
            m["at_0_5"] = EV.at_threshold(by, entries, 0.5)
            if arm in LEARNED_ARMS:
                dt = display_values[arm]
                m["at_model_display_threshold"] = {"threshold": dt, **(EV.at_threshold(by, entries, dt) or {})}
            report[target_key][arm] = m

    # ---- E5: learned_counts per learned arm (all genomes and without GCA_000006965.1)
    pop_no_sino = {(g, pid) for (g, pid) in pop if g != SENSITIVITY_GENOME}
    for arm in requested_learned_arms:
        if arm in missing_arms:
            report["learned"][arm] = {
                "status": "missing",
                "missing_genomes": missing_arms[arm]["missing_genomes"],
            }
            report["learned_without_GCA_000006965.1"][arm] = {
                "status": "missing",
                "missing_genomes": missing_arms[arm]["missing_genomes"],
            }
            continue
        arm_recs = receipts_by_arm[arm]
        report["learned"][arm] = EV.learned_counts(arm_recs, learned, arm, pop)
        arm_recs_no_sino = {g: r for g, r in arm_recs.items() if g != SENSITIVITY_GENOME}
        learned_no_sino = {a: {g: v for g, v in by_g.items() if g != SENSITIVITY_GENOME} for a, by_g in learned.items()}
        report["learned_without_GCA_000006965.1"][arm] = EV.learned_counts(
            arm_recs_no_sino, learned_no_sino, arm, pop_no_sino
        )

    R1.save(out / "evaluation.json", report)
    write_markdown_tables(report, out)
    return report


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pn-cells", type=Path, required=True, help="Directory with PN cells/<index>/{receipt,retained}.json")
    p.add_argument("--pu-cells", type=Path, required=True, help="Directory with PU cells/<index>/{receipt,retained}.json")
    p.add_argument("--pu-rank-cells", type=Path, default=None, help="Optional directory with PU-rank cells/<index>/{receipt,retained}.json")
    p.add_argument("--claims-dir", type=Path, required=True, help="Directory with <genome>.claims.jsonl files")
    p.add_argument("--split", type=Path, required=True, help="Frozen Price row-1 split.json")
    p.add_argument("--ec-audit", type=Path, required=True, help="EC identity transfer audit JSON")
    p.add_argument("--out", type=Path, required=True, help="Output directory (must not exist)")
    p.add_argument("--pn-model", type=Path, default=None, help="Optional PN logistic artifact for display_threshold")
    p.add_argument("--pn-model-gbdt", type=Path, default=None, help="Optional PN GBDT artifact for display_threshold")
    p.add_argument("--pu-model", type=Path, default=None, help="Optional PU logistic artifact for display_threshold")
    p.add_argument("--pu-model-gbdt", type=Path, default=None, help="Optional PU GBDT artifact for display_threshold")
    p.add_argument("--pu-rank-model", type=Path, default=None, help="Optional PU-rank logistic artifact for display_threshold")
    return p


def main(argv=None) -> int:
    ap = build_parser()
    args = ap.parse_args(argv)
    try:
        run(args)
    except (ValueError, FileExistsError, FileNotFoundError) as exc:
        print(f"REFUSED: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
