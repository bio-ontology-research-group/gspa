#!/usr/bin/env python3
"""Join Price-148 and Swiss-Prot experimental EC truth to frozen Universe A and Universe B bundles.

Refuses to run unless every public file's SHA-256 equals the digest recorded in PUBLIC_SHA256
(the truth-separation gate). Opens label files only after the gate passes.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "benchmark" / "v3"))
import evaluate_stage3_agent as E  # noqa: E402
import v3_runner as runner  # noqa: E402
from price_contract_learned import row11_truth as R11  # noqa: E402

COMPLETE = R11.COMPLETE
ec2go = R11.ec2go
sha = R11.sha


def verify_separation_gate(public_sha256_path: Path, public_files: list[Path]) -> dict[str, str]:
    """Verify that PUBLIC_SHA256 exists and every recorded and requested public file matches its SHA-256."""
    sha_path = Path(public_sha256_path)
    if not sha_path.is_file():
        raise ValueError(f"Separation gate failed: PUBLIC_SHA256 file does not exist: {sha_path}")
    recorded_by_raw: dict[str, str] = {}
    recorded_by_resolved: dict[str, str] = {}
    for line in sha_path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        parts = line.split(None, 1)
        if len(parts) != 2:
            raise ValueError(f"Separation gate failed: malformed line in {sha_path}: {line!r}")
        digest_hex, p_str = parts[0].strip(), parts[1].strip()
        p = Path(p_str)
        if not p.is_file():
            raise ValueError(f"Separation gate failed: recorded public file missing on disk: {p_str}")
        actual = runner.sha256_file(str(p))
        if actual != digest_hex:
            raise ValueError(
                f"Separation gate failed: public file {p_str} sha256 ({actual}) does not match recorded ({digest_hex})"
            )
        recorded_by_raw[p_str] = digest_hex
        recorded_by_resolved[str(p.resolve())] = digest_hex

    if not recorded_by_resolved:
        raise ValueError(f"Separation gate failed: no public files recorded in {sha_path}")

    verified: dict[str, str] = {}
    for req in public_files:
        req_p = Path(req)
        if not req_p.is_file():
            raise ValueError(f"Separation gate failed: requested public file does not exist: {req}")
        resolved = str(req_p.resolve())
        if str(req_p) not in recorded_by_raw and resolved not in recorded_by_resolved:
            raise ValueError(
                f"Separation gate failed: public file {req} is not listed in {sha_path}"
            )
        expected = recorded_by_raw.get(str(req_p)) or recorded_by_resolved[resolved]
        actual = runner.sha256_file(str(req_p))
        if actual != expected:
            raise ValueError(
                f"Separation gate failed: public file {req} sha256 ({actual}) != recorded ({expected})"
            )
        verified[str(req_p)] = actual
    return verified


def load_all_labels(
    price_labels_dir: Path,
    swissprot_tsv: Path,
    leakage_json: Path,
) -> dict[tuple[str, str], set[tuple[str, str]]]:
    """Load Price-148 complete ECs and Swiss-Prot experimental ECs, tagging circular Swiss-Prot labels."""
    labels: dict[tuple[str, str], set[tuple[str, str]]] = defaultdict(set)
    for p in sorted(Path(price_labels_dir).glob("*.ec.tsv")):
        g = p.name.split(".ec.tsv")[0]
        with p.open("r", encoding="utf-8", newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                m = COMPLETE.match(r["ec"].strip())
                if m:
                    labels[(g, r["protein_id"].strip())].add((m.group(1), "price"))

    leak_data = json.loads(Path(leakage_json).read_text(encoding="utf-8"))
    circ = {
        (
            c["genome_id"].strip(),
            c["protein_id"].strip(),
            c["function_id"].strip().removeprefix("EC:"),
        )
        for c in leak_data.get("swissprot_price_labels_in_reference", [])
    }
    with Path(swissprot_tsv).open("r", encoding="utf-8", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            m = COMPLETE.match(r["function_id"].strip())
            if m:
                gid = r["genome_id"].strip()
                pid = r["protein_id"].strip()
                ec_val = m.group(1)
                src = "swissprot_circular" if (gid, pid, ec_val) in circ else "swissprot"
                labels[(gid, pid)].add((ec_val, src))
    return labels


def join_universe_a(
    public_a: dict,
    labels: dict[tuple[str, str], set[tuple[str, str]]],
    mapping: dict[str, set[str]],
    public_sha256: str = "",
) -> tuple[dict[str, list[dict]], dict]:
    """Apply the row11_truth rule to Universe A target-centred cases."""
    E.reject_truth(public_a)
    cases = public_a.get("cases", [])
    counts: Counter = Counter()
    by_split: dict[str, Counter] = {"calibration": Counter(), "test": Counter()}
    truth_rows: dict[str, list[dict]] = {"calibration": [], "test": []}
    seen_truth_keys: dict[str, set[tuple[str, str, str]]] = {"calibration": set(), "test": set()}
    detail: list[dict] = []
    target_kinds: Counter = Counter()
    cases_by_type_split: dict[str, Counter] = {"calibration": Counter(), "test": Counter()}

    for case in cases:
        an = case["anomaly"]
        g, split, target = case["genome_id"], case["split"], an["target_function_id"]
        kind = "GO" if target.startswith("GO:") else "EC"
        atype = an.get("anomaly_type") or an.get("type") or "UNKNOWN"
        target_kinds[(atype, kind)] += 1
        cases_by_type_split[split][f"{atype}/{kind}"] += 1
        for c in an.get("candidates", []):
            pid = c["protein_id"]
            lab = labels.get((g, pid), set())
            status, source, neg_rule = "unknown", "", False
            if lab:
                if kind == "EC":
                    t = target.removeprefix("EC:")
                    hit = [s for e, s in lab if e == t]
                    neg_rule = not hit and COMPLETE.match(target) is not None
                else:
                    hit = [s for e, s in lab if target in mapping.get(e, ())]
                if hit:
                    status = "positive"
                    non_circ = sorted(set(hit) - {"swissprot_circular"})
                    source = "swissprot_circular" if not non_circ else non_circ[0]
            key = status if status == "unknown" else f"positive/{source}"
            counts[key] += 1
            by_split[split][key] += 1
            if neg_rule:
                counts["negative_under_existing_rule_not_used"] += 1
                by_split[split]["negative_under_existing_rule_not_used"] += 1
            if lab:
                counts["candidate_pairs_on_labelled_proteins"] += 1
                by_split[split]["candidate_pairs_on_labelled_proteins"] += 1
                detail.append({
                    "case": an["anomaly_id"],
                    "split": split,
                    "genome": g,
                    "protein": pid,
                    "target": target,
                    "labels": [list(x) for x in sorted(lab)],
                    "status": status,
                    "source": source,
                })
            if status == "positive" and source != "swissprot_circular":
                tkey = (g, pid, target)
                if tkey not in seen_truth_keys[split]:
                    seen_truth_keys[split].add(tkey)
                    truth_rows[split].append({
                        "genome_id": g,
                        "protein_id": pid,
                        "function_id": target,
                        "status": "positive",
                        "source": source,
                    })

    report_a = {
        "public_sha256": public_sha256 or E.digest(public_a),
        "cases": len(cases),
        "cases_by_split": dict(Counter(c["split"] for c in cases)),
        "cases_by_type_and_split": {sp: dict(c) for sp, c in cases_by_type_split.items()},
        "candidate_pairs": sum(len(c["anomaly"].get("candidates", [])) for c in cases),
        "case_types": {f"{t}/{k}": n for (t, k), n in sorted(target_kinds.items(), key=str)},
        "pairs": dict(counts),
        "pairs_by_split": {k: dict(v) for k, v in by_split.items()},
        "labelled_pairs_detail": detail,
        "methods_applicable": all(bool(truth_rows[s]) for s in ("calibration", "test")),
    }
    return truth_rows, report_a


def join_universe_b(
    public_b: dict,
    private_b: dict,
    labels: dict[tuple[str, str], set[tuple[str, str]]],
    public_sha256: str = "",
) -> tuple[dict, dict]:
    """Join EC truth to Universe B locus cases, excluding circular Swiss-Prot labels."""
    E.reject_truth(public_b)
    if public_b.get("schema") != "gspa.locus-public/1":
        raise ValueError("expected Universe B public bundle v1")
    priv_cases = private_b.get("cases") or {k: v for k, v in private_b.items() if k.startswith("L_")}
    cases = public_b.get("cases", [])
    if set(c["anomaly_id"] for c in cases) != set(priv_cases):
        raise ValueError("private_B case IDs do not match public_B case IDs")

    truth_cases: dict[str, dict] = {}
    cov_by_split: dict[str, Counter] = {"calibration": Counter(), "test": Counter()}
    rank_hist_by_split: dict[str, Counter] = {"calibration": Counter(), "test": Counter()}

    for case in cases:
        aid = case["anomaly_id"]
        split = case["split"]
        gid = case["genome_id"]
        loc = case["locus"]
        priv = priv_cases[aid]
        role = priv["role"]
        pid = priv["protein_id"]
        if pid != loc["protein_id"]:
            raise ValueError(f"protein_id mismatch on {aid}")

        if role == "decoy":
            truth_cases[aid] = {
                "anomaly_id": aid,
                "split": split,
                "genome_id": gid,
                "role": "decoy",
                "positives": None,
                "labelled_ecs_not_among_candidates": [],
                "sources": [],
                "circular_excluded_ecs": [],
                "truth_ranks": [],
                "best_truth_rank": None,
            }
            cov_by_split[split]["n_decoys"] += 1
            continue

        cov_by_split[split]["n_labelled_cases"] += 1
        lab = labels.get((gid, pid), set())
        srcs_by_ec: dict[str, set[str]] = defaultdict(set)
        for ec_val, src in lab:
            srcs_by_ec[ec_val].add(src)

        valid_ecs: dict[str, list[str]] = {}
        circular_ecs: list[str] = []
        for ec_val, srcs in sorted(srcs_by_ec.items()):
            non_circ = sorted(srcs - {"swissprot_circular"})
            if non_circ:
                valid_ecs[ec_val] = non_circ
            else:
                circular_ecs.append(f"EC:{ec_val}")

        if circular_ecs:
            cov_by_split[split]["circular_excluded_cases"] += 1
            cov_by_split[split]["circular_excluded_pairs"] += len(circular_ecs)

        candidates = loc.get("candidates", [])
        positives: list[str] = []
        truth_ranks: list[int] = []
        cand_ec_set: set[str] = set()
        for idx, cand in enumerate(candidates, 1):
            fid = cand["function_id"]
            ec_bare = fid.removeprefix("EC:")
            cand_ec_set.add(ec_bare)
            if ec_bare in valid_ecs:
                positives.append(fid)
                truth_ranks.append(idx)

        not_among = sorted(f"EC:{ec_val}" for ec_val in valid_ecs if ec_val not in cand_ec_set)
        all_sources = sorted({s for s_list in valid_ecs.values() for s in s_list})
        best_rank = min(truth_ranks) if truth_ranks else None

        if best_rank is not None:
            cov_by_split[split]["truth_among_candidates"] += 1
            rank_hist_by_split[split][str(best_rank)] += 1
            if best_rank == 1:
                cov_by_split[split]["rank_1"] += 1
            elif 2 <= best_rank <= 10:
                cov_by_split[split]["rank_2_10"] += 1
            else:
                cov_by_split[split]["rank_gt_10"] += 1
        else:
            cov_by_split[split]["absent"] += 1

        truth_cases[aid] = {
            "anomaly_id": aid,
            "split": split,
            "genome_id": gid,
            "role": "labelled",
            "positives": positives,
            "labelled_ecs_not_among_candidates": not_among,
            "sources": all_sources,
            "sources_by_ec": {f"EC:{k}": v for k, v in valid_ecs.items()},
            "circular_excluded_ecs": circular_ecs,
            "truth_ranks": truth_ranks,
            "best_truth_rank": best_rank,
        }

    coverage_summary: dict[str, dict] = {}
    total_cov: Counter = Counter()
    total_rank_hist: Counter = Counter()
    for sp in ("calibration", "test"):
        c = cov_by_split[sp]
        coverage_summary[sp] = {
            "n_labelled_cases": c["n_labelled_cases"],
            "n_decoys": c["n_decoys"],
            "truth_among_candidates": c["truth_among_candidates"],
            "rank_1": c["rank_1"],
            "rank_2_10": c["rank_2_10"],
            "rank_gt_10": c["rank_gt_10"],
            "absent": c["absent"],
            "circular_excluded_cases": c["circular_excluded_cases"],
            "circular_excluded_pairs": c["circular_excluded_pairs"],
            "by_rank": dict(rank_hist_by_split[sp]),
        }
        for k in (
            "n_labelled_cases",
            "n_decoys",
            "truth_among_candidates",
            "rank_1",
            "rank_2_10",
            "rank_gt_10",
            "absent",
            "circular_excluded_cases",
            "circular_excluded_pairs",
        ):
            total_cov[k] += c[k]
        total_rank_hist.update(rank_hist_by_split[sp])
    coverage_summary["total"] = {
        **dict(total_cov),
        "by_rank": dict(total_rank_hist),
    }

    pub_digest = public_sha256 or E.digest(public_b)
    truth_b = {
        "schema": "gspa.locus-truth/1",
        "public_sha256": pub_digest,
        "cases": truth_cases,
        "candidate_truth_coverage": coverage_summary,
        "counts_by_split": coverage_summary,
        **truth_cases,
    }
    report_b = {
        "public_sha256": pub_digest,
        "n_cases": len(cases),
        "candidate_truth_coverage": coverage_summary,
        "counts_by_split": coverage_summary,
    }
    return truth_b, report_b


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--public-sha256", required=True, type=Path)
    parser.add_argument("--public-a", type=Path, default=None)
    parser.add_argument("--public-b", required=True, type=Path)
    parser.add_argument("--private-b", required=True, type=Path)
    parser.add_argument("--price-labels", required=True, type=Path)
    parser.add_argument("--swissprot", "--swissprot-ecs", dest="swissprot", required=True, type=Path)
    parser.add_argument("--leakage", required=True, type=Path)
    parser.add_argument("--ec2go", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args(argv)

    req_public = [args.public_b]
    if args.public_a is not None:
        req_public.insert(0, args.public_a)
    verified_shas = verify_separation_gate(args.public_sha256, req_public)

    labels = load_all_labels(args.price_labels, args.swissprot, args.leakage)
    mapping = ec2go(args.ec2go)

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    report_a = None
    if args.public_a is not None:
        pub_a = E.read(args.public_a)
        truth_rows_a, report_a = join_universe_a(
            pub_a,
            labels,
            mapping,
            public_sha256=verified_shas[str(args.public_a)],
        )
        for sp in ("calibration", "test"):
            E.write(out_dir / f"truth_A_{sp}.json", truth_rows_a[sp])

    pub_b = E.read(args.public_b)
    priv_b = E.read(args.private_b)
    truth_b, report_b = join_universe_b(
        pub_b,
        priv_b,
        labels,
        public_sha256=verified_shas[str(args.public_b)],
    )
    E.write(out_dir / "truth_B.json", truth_b)

    full_report = {
        "schema": "gspa.row11-truth-report/1",
        "separation_gate": {
            "public_sha256_file": str(args.public_sha256),
            "verified_files": verified_shas,
            "status": "verified",
        },
        "label_sources": {
            "price_tables": len(list(Path(args.price_labels).glob("*.ec.tsv"))),
            "swissprot_sha256": sha(args.swissprot),
            "leakage_sha256": sha(args.leakage),
            "ec2go_sha256": sha(args.ec2go),
        },
        "universe_a": report_a,
        "universe_b": report_b,
    }
    E.write(out_dir / "truth_report.json", full_report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
