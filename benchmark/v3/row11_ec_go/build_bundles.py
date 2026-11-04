#!/usr/bin/env python3
"""Build truth-separated public/private anomaly bundles for the row-11 EC-aware comparison.

Subcommands:
  universe-a: target-centred production cases via evaluate_stage3_agent.prepare (+ optional per-genome cap).
  universe-b: locus-anchored cases from *.anomalies.locus.jsonl + 2:1 matched decoys, shuffled with fixed seed.
  hash:       record SHA-256 of every frozen public bundle before any label file is opened.
"""
from __future__ import annotations

import argparse
import copy
import csv
import json
import random
import sys
from collections import Counter
from pathlib import Path
from unittest.mock import patch

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "benchmark" / "v3"))
import evaluate_stage3_agent as E  # noqa: E402
import v3_runner as runner  # noqa: E402

LOCUS_ANOMALY_SCHEMA = "gspa.locus-anomalies/1"
LOCUS_PUBLIC_SCHEMA = "gspa.locus-public/1"
LOCUS_PRIVATE_SCHEMA = "gspa.locus-private/1"
FORBIDDEN_PUBLIC_KEYS = E.FORBIDDEN | {"role", "decoy", "decoys", "decoy_of", "labelled", "is_decoy", "is_labelled"}

UNIVERSE_B_LIMITATIONS = [
    "Universe B is locus-anchored on proteins with experimental EC truth (locus identity only, never the EC label) "
    "together with truth-free matched decoys.",
    "Decoy rule: for each labelled locus case, select up to decoys_per_case locus cases from the same genome and "
    "anomaly_type whose protein is not in loci.tsv and shares no cohort orthogroup_id with any loci.tsv protein of "
    "that genome, nearest in top candidate current_posterior (tie-broken by protein_id), each decoy used at most once.",
]


def reject_public_leakage(value, _parent_key: str = "") -> None:
    """Enforce evaluate_stage3_agent.reject_truth and forbid role/decoy/loci-source keys in public bundles."""
    E.reject_truth(value)
    if isinstance(value, dict):
        bad = FORBIDDEN_PUBLIC_KEYS & set(value)
        if _parent_key != "supporting_claims" and "source" in value:
            bad = bad | {"source"}
        if bad:
            raise ValueError(f"forbidden role/truth fields in public bundle: {sorted(bad)}")
        for k, child in value.items():
            reject_public_leakage(child, k)
    elif isinstance(value, list):
        for child in value:
            reject_public_leakage(child, _parent_key)


def verify_and_load_runs(source: dict) -> dict[str, dict]:
    """Verify production runs as evaluate_stage3_agent.prepare does, supporting abstained learned arms."""
    if source.get("schema") != "gspa.agent-evidence-sources/2":
        raise ValueError("expected production evidence sources v2")
    runs = {name: E.verified_run(path) for name, path in source.get("runs", {}).items()}
    abstained = source.get("abstained_arms") or {}
    if set(runs) != {"stage1", "learned_logistic", "learned_gbdt"}:
        if "stage1" in runs and set(runs) | set(abstained) == {"stage1", "learned_logistic", "learned_gbdt"}:
            base = runs["stage1"]
            for name, info in abstained.items():
                rec_path = info.get("run_record")
                if not rec_path or not Path(rec_path).is_file():
                    raise ValueError(f"{name}: missing run_record for abstained arm")
                old = E.read(rec_path)
                identity = old.get("inference_identity")
                if not identity:
                    raise ValueError(f"{name}: abstained run lacks input identity")
                if old.get("returncode") != 0:
                    raise ValueError(f"{name}: abstained run had non-zero returncode")
                reasons = set((old.get("failed") or {}).values())
                if reasons != {"invalid execution evidence: required stage learned_integration did not execute"}:
                    raise ValueError(f"{name}: unexpected failure reasons in abstained arm: {reasons}")
                runs[name] = {**old, "_abstained": True, "genomes": {
                    gid: {"path": f"__abstained__:{name}:{gid}", "status": "abstained"}
                    for gid in base["genomes"]
                }}
        else:
            raise ValueError("require baseline plus actual logistic and GBDT production runs")
    base = runs["stage1"]
    for name, run in runs.items():
        if run["inference_identity"]["genomes"] != base["inference_identity"]["genomes"]:
            raise ValueError(f"{name}: unmatched claims or genome context")
        if run["inference_identity"]["launch"] != base["inference_identity"]["launch"]:
            raise ValueError(f"{name}: unmatched production configuration")
        if name != "stage1":
            args = run["inference_identity"]["arm"]["args"]
            if len(args) != 2 or args[0] != "--learned-model":
                raise ValueError("model arm must differ only by its frozen learned artifact")
            model_path = str(Path(args[1]).resolve())
            expected_files = {k: v for k, v in run["inference_identity"]["files"].items() if k != model_path}
            if expected_files != base["inference_identity"]["files"]:
                raise ValueError("model arm reference or executable differs from baseline")
    return runs


def apply_cap_universe_a(bundle: dict, cap_per_genome: int | None = None) -> dict:
    """Keep the first N cases per genome in file order (truth-free) and record kept/total in coverage."""
    if cap_per_genome is not None and cap_per_genome < 0:
        raise ValueError("cap_per_genome must be non-negative")
    out = copy.deepcopy(bundle)
    kept_cases = []
    kept_counts: Counter = Counter()
    total_counts: Counter = Counter()
    for case in out.get("cases", []):
        gid = case["genome_id"]
        total_counts[gid] += 1
        if cap_per_genome is None or kept_counts[gid] < cap_per_genome:
            kept_cases.append(case)
            kept_counts[gid] += 1
    out["cases"] = kept_cases
    out["cap_per_genome"] = cap_per_genome
    for cov in out.get("coverage", []):
        gid = cov["genome_id"]
        total = cov.get("n_anomalies", total_counts[gid])
        kept = kept_counts[gid] if cap_per_genome is not None else total
        cov["total"] = total
        cov["kept"] = kept
        cov["cap_per_genome"] = cap_per_genome
    return out


def build_universe_a(source: dict, cap_per_genome: int | None = None) -> dict:
    """Build Universe A public bundle via evaluate_stage3_agent.prepare with optional per-genome cap."""
    abstained = source.get("abstained_arms") or {}
    if set(source.get("runs", {})) == {"stage1"} and set(abstained) == {"learned_logistic", "learned_gbdt"}:
        runs = verify_and_load_runs(source)
        patched_source = copy.deepcopy(source)
        for name, info in abstained.items():
            patched_source["runs"][name] = info["run_record"]
        orig_verified_run = E.verified_run
        orig_parse = runner.parse_integrated

        def _verified_run(path):
            for name, info in abstained.items():
                if str(path) == str(info["run_record"]):
                    return runs[name]
            return orig_verified_run(path)

        def _parse_integrated(path, kind):
            if str(path).startswith("__abstained__:"):
                return {}
            return orig_parse(path, kind)

        with patch.object(E, "verified_run", side_effect=_verified_run), patch.object(
            runner, "parse_integrated", side_effect=_parse_integrated
        ):
            bundle = E.prepare(patched_source)
        bundle["source_sha256"] = E.digest(source)
    else:
        bundle = E.prepare(source)
    return apply_cap_universe_a(bundle, cap_per_genome=cap_per_genome)


def read_loci_tsv(loci_path: Path) -> list[dict]:
    """Read frozen loci.tsv (genome_id, protein_id, source) and verify LOCI_SHA256 if present next to it."""
    loci_path = Path(loci_path)
    sha_file = loci_path.parent / "LOCI_SHA256"
    if sha_file.is_file():
        expected = sha_file.read_text(encoding="utf-8").strip().split()[0]
        actual = runner.sha256_file(str(loci_path))
        if actual != expected:
            raise ValueError(f"loci.tsv sha256 mismatch against {sha_file}: {actual} != {expected}")
    rows = []
    seen = set()
    with loci_path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames or not {"genome_id", "protein_id", "source"}.issubset(set(reader.fieldnames)):
            raise ValueError(f"{loci_path}: expected columns genome_id, protein_id, source")
        if set(reader.fieldnames) & E.FORBIDDEN:
            raise ValueError(f"{loci_path}: forbidden truth column in loci file")
        for r in reader:
            gid = r["genome_id"].strip()
            pid = r["protein_id"].strip()
            src = r["source"].strip()
            if (gid, pid) in seen:
                raise ValueError(f"{loci_path}: duplicate (genome_id, protein_id) = {(gid, pid)}")
            seen.add((gid, pid))
            rows.append({"genome_id": gid, "protein_id": pid, "source": src})
    return rows


def read_orthogroups_tsv(orthogroups_path: Path) -> dict[tuple[str, str], set[str]]:
    """Read cohort_orthogroups.tsv into {(genome_id, protein_id): {orthogroup_id, ...}}."""
    out: dict[tuple[str, str], set[str]] = {}
    with Path(orthogroups_path).open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames or not {"genome_id", "protein_id", "orthogroup_id"}.issubset(set(reader.fieldnames)):
            raise ValueError(f"{orthogroups_path}: expected columns genome_id, protein_id, orthogroup_id")
        for r in reader:
            gid = r["genome_id"].strip()
            pid = r["protein_id"].strip()
            og = r["orthogroup_id"].strip()
            if og and og != "-":
                out.setdefault((gid, pid), set()).add(og)
    return out


def classify_missing_locus(stage1_ec_posteriors: dict[tuple[str, str], float], protein_id: str) -> str:
    """Classify a loci.tsv protein that has no locus anomaly case from its stage1 EC rows."""
    scores = [p for (pid, _ec), p in stage1_ec_posteriors.items() if pid == protein_id]
    n_ge_05 = sum(1 for p in scores if p >= 0.5)
    n_ge_02 = sum(1 for p in scores if p >= 0.2)
    if n_ge_05 == 1 and n_ge_02 == 1:
        return "resolved_single_call"
    if n_ge_02 == 0:
        return "no_candidate"
    return "unresolved_not_a_case"


def make_locus_case(
    gid: str,
    split: str,
    record: dict,
    predictions: dict[str, dict[tuple[str, str], float]],
) -> dict:
    """Sanitise one locus anomaly record, compute opaque L_ id, and attach baseline/model/raw scores."""
    E.reject_truth(record)
    if record.get("schema") != LOCUS_ANOMALY_SCHEMA:
        raise ValueError(f"expected schema {LOCUS_ANOMALY_SCHEMA!r}, got {record.get('schema')!r}")
    loc = copy.deepcopy(record)
    pid = loc["protein_id"]
    loc.pop("anomaly_id", None)
    aid = "L_" + E.digest([gid, pid, loc])[:24]
    loc["anomaly_id"] = aid

    candidates = loc.get("candidates", [])
    fids = [c["function_id"] for c in candidates]
    if len(set(fids)) != len(fids):
        raise ValueError(f"duplicate candidate function_id on {gid}/{pid}")
    posteriors = [float(c["current_posterior"]) for c in candidates]
    if posteriors != sorted(posteriors, reverse=True):
        raise ValueError(f"candidates not sorted by current_posterior descending on {gid}/{pid}")

    for c in candidates:
        ec_bare = c["function_id"].removeprefix("EC:")
        base_prob = predictions["stage1"].get((pid, ec_bare), 0.0)
        c["baseline_probability"] = E.probability(base_prob)
        c["model_probabilities"] = {}
        for mode in ("learned_logistic", "learned_gbdt"):
            val = predictions.get(mode, {}).get((pid, ec_bare))
            c["model_probabilities"][mode] = E.probability(val) if val is not None else None
        c["raw_sequence_score"] = max(
            (float(e["raw_score"]) for e in c.get("supporting_claims", []) if e.get("raw_score") is not None),
            default=None,
        )
    case_obj = {"anomaly_id": aid, "split": split, "genome_id": gid, "locus": loc}
    reject_public_leakage(case_obj)
    return case_obj


def assemble_universe_b(
    genome_data: dict[str, dict],
    loci_rows: list[dict],
    orthogroups: dict[tuple[str, str], set[str]],
    decoys_per_case: int = 2,
    seed: int = 20261102,
    source_sha256: str = "",
    source_receipts: dict | None = None,
    shared_production_wall_seconds: float = 0.0,
    limitations: list[str] | None = None,
) -> tuple[dict, dict]:
    """Assemble Universe B public and private bundles from per-genome locus records and EC posteriors.

    Pure function testable without production run receipts.
    """
    if decoys_per_case < 0:
        raise ValueError("decoys_per_case must be non-negative")

    loci_by_genome: dict[str, list[dict]] = {}
    all_loci_pairs: set[tuple[str, str]] = set()
    for row in loci_rows:
        gid, pid = row["genome_id"], row["protein_id"]
        loci_by_genome.setdefault(gid, []).append(row)
        all_loci_pairs.add((gid, pid))

    public_cases: list[dict] = []
    private_cases: dict[str, dict] = {}
    public_loci_by_genome: dict[str, dict] = {}
    public_loci_by_split: dict[str, Counter] = {
        "calibration": Counter(),
        "test": Counter(),
    }
    private_loci_by_genome: dict[str, dict] = {}
    public_shortfalls: list[dict] = []
    private_shortfalls: list[dict] = []

    for gid in sorted(genome_data):
        ginfo = genome_data[gid]
        split = ginfo["split"]
        if split not in ("calibration", "test"):
            raise ValueError(f"{gid}: only calibration and test splits are admitted")
        locus_records = ginfo["locus_records"]
        predictions = ginfo["predictions"]

        by_pid: dict[str, dict] = {}
        for rec in locus_records:
            E.reject_truth(rec)
            if rec.get("schema") != LOCUS_ANOMALY_SCHEMA:
                raise ValueError(f"{gid}: expected {LOCUS_ANOMALY_SCHEMA!r}")
            pid = rec["protein_id"]
            if pid in by_pid:
                raise ValueError(f"{gid}: duplicate protein_id {pid} in locus anomaly file")
            by_pid[pid] = rec

        g_loci = loci_by_genome.get(gid, [])
        loci_pids_in_g = {r["protein_id"] for r in g_loci}
        loci_ogs_in_g = {
            og
            for pid in loci_pids_in_g
            for og in orthogroups.get((gid, pid), set())
            if og and og != "-"
        }

        status_counts = Counter({
            "case": 0,
            "resolved_single_call": 0,
            "no_candidate": 0,
            "unresolved_not_a_case": 0,
        })
        identities_by_status: dict[str, list[dict]] = {
            "case": [],
            "resolved_single_call": [],
            "no_candidate": [],
            "unresolved_not_a_case": [],
        }

        labelled_for_g: list[tuple[dict, dict, dict]] = []
        for lrow in g_loci:
            pid = lrow["protein_id"]
            src = lrow["source"]
            if pid in by_pid:
                status = "case"
                case_obj = make_locus_case(gid, split, by_pid[pid], predictions)
                aid = case_obj["anomaly_id"]
                labelled_for_g.append((lrow, by_pid[pid], case_obj))
                identities_by_status[status].append(
                    {"genome_id": gid, "protein_id": pid, "source": src, "anomaly_id": aid}
                )
            else:
                status = classify_missing_locus(predictions["stage1"], pid)
                identities_by_status[status].append(
                    {"genome_id": gid, "protein_id": pid, "source": src}
                )
            status_counts[status] += 1

        used_decoy_pids: set[str] = set()
        n_decoys_for_g = 0
        shortfall_for_g = 0

        for lrow, raw_lab_rec, lab_case_obj in labelled_for_g:
            lab_aid = lab_case_obj["anomaly_id"]
            lab_pid = lrow["protein_id"]
            lab_src = lrow["source"]
            lab_type = raw_lab_rec["anomaly_type"]
            lab_top = float(raw_lab_rec["candidates"][0]["current_posterior"]) if raw_lab_rec.get("candidates") else 0.0

            public_cases.append(lab_case_obj)
            private_cases[lab_aid] = {
                "role": "labelled",
                "genome_id": gid,
                "split": split,
                "protein_id": lab_pid,
                "source": lab_src,
                "decoy_of": None,
            }

            eligible_decoys = []
            for cand_pid, cand_rec in by_pid.items():
                if cand_rec["anomaly_type"] != lab_type:
                    continue
                if cand_pid in loci_pids_in_g or (gid, cand_pid) in all_loci_pairs:
                    continue
                if cand_pid in used_decoy_pids:
                    continue
                cand_ogs = orthogroups.get((gid, cand_pid), set())
                if cand_ogs & loci_ogs_in_g:
                    continue
                cand_top = float(cand_rec["candidates"][0]["current_posterior"]) if cand_rec.get("candidates") else 0.0
                eligible_decoys.append((abs(cand_top - lab_top), cand_pid, cand_rec))

            eligible_decoys.sort(key=lambda item: (item[0], item[1]))
            chosen = eligible_decoys[:decoys_per_case]

            for _diff, dec_pid, dec_rec in chosen:
                used_decoy_pids.add(dec_pid)
                dec_case_obj = make_locus_case(gid, split, dec_rec, predictions)
                dec_aid = dec_case_obj["anomaly_id"]
                public_cases.append(dec_case_obj)
                private_cases[dec_aid] = {
                    "role": "decoy",
                    "genome_id": gid,
                    "split": split,
                    "protein_id": dec_pid,
                    "source": None,
                    "decoy_of": lab_aid,
                }
                n_decoys_for_g += 1

            if len(chosen) < decoys_per_case:
                miss = decoys_per_case - len(chosen)
                shortfall_for_g += miss
                public_shortfalls.append({
                    "genome_id": gid,
                    "split": split,
                    "anomaly_type": lab_type,
                    "requested": decoys_per_case,
                    "found": len(chosen),
                    "shortfall": miss,
                })
                private_shortfalls.append({
                    "genome_id": gid,
                    "split": split,
                    "labelled_anomaly_id": lab_aid,
                    "protein_id": lab_pid,
                    "anomaly_type": lab_type,
                    "requested": decoys_per_case,
                    "found": len(chosen),
                    "shortfall": miss,
                })

        g_summary_public = {
            "split": split,
            "case": status_counts["case"],
            "resolved_single_call": status_counts["resolved_single_call"],
            "no_candidate": status_counts["no_candidate"],
            "unresolved_not_a_case": status_counts["unresolved_not_a_case"],
            "n_bundle_cases": status_counts["case"] + n_decoys_for_g,
            "shortfalls": shortfall_for_g,
            "n_total_locus_anomalies": len(locus_records),
        }
        public_loci_by_genome[gid] = g_summary_public
        for k in ("case", "resolved_single_call", "no_candidate", "unresolved_not_a_case", "n_bundle_cases", "shortfalls"):
            public_loci_by_split[split][k] += g_summary_public[k]

        private_loci_by_genome[gid] = {
            **g_summary_public,
            "n_labelled": status_counts["case"],
            "n_decoys": n_decoys_for_g,
            "identities": identities_by_status,
        }

    ids = [c["anomaly_id"] for c in public_cases]
    if len(set(ids)) != len(ids):
        raise ValueError("duplicate opaque case ID in Universe B")

    random.Random(seed).shuffle(public_cases)

    by_split_public = {}
    total_counter: Counter = Counter()
    for sp in ("calibration", "test"):
        c = public_loci_by_split[sp]
        by_split_public[sp] = {
            "case": c["case"],
            "resolved_single_call": c["resolved_single_call"],
            "no_candidate": c["no_candidate"],
            "unresolved_not_a_case": c["unresolved_not_a_case"],
            "n_bundle_cases": c["n_bundle_cases"],
            "shortfalls": c["shortfalls"],
        }
        total_counter.update(by_split_public[sp])

    full_limitations = list(limitations or []) + UNIVERSE_B_LIMITATIONS
    public_bundle = {
        "schema": LOCUS_PUBLIC_SCHEMA,
        "cases": public_cases,
        "coverage": {
            "loci": {
                "by_genome": public_loci_by_genome,
                "by_split": by_split_public,
                "total": dict(total_counter),
            },
            "shortfalls": public_shortfalls,
        },
        "source_sha256": source_sha256,
        "source_receipts": source_receipts or {},
        "shared_production_wall_seconds": shared_production_wall_seconds,
        "limitations": full_limitations,
    }
    reject_public_leakage(public_bundle)

    private_by_split = {}
    for sp in ("calibration", "test"):
        n_lab = sum(1 for v in private_cases.values() if v["split"] == sp and v["role"] == "labelled")
        n_dec = sum(1 for v in private_cases.values() if v["split"] == sp and v["role"] == "decoy")
        private_by_split[sp] = {
            **by_split_public[sp],
            "labelled": n_lab,
            "decoy": n_dec,
        }

    private_bundle = {
        "schema": LOCUS_PRIVATE_SCHEMA,
        "public_sha256": E.digest(public_bundle),
        "seed": seed,
        "decoys_per_case": decoys_per_case,
        "cases": private_cases,
        "coverage": {
            "loci": {
                "by_genome": private_loci_by_genome,
                "by_split": private_by_split,
            },
            "shortfalls": private_shortfalls,
        },
        **private_cases,
    }
    return public_bundle, private_bundle


def build_universe_b(
    source: dict,
    loci_path: Path,
    orthogroups_path: Path,
    decoys_per_case: int = 2,
    seed: int = 20261102,
) -> tuple[dict, dict]:
    """Verify production runs and assemble Universe B public and private bundles."""
    runs = verify_and_load_runs(source)
    base = runs["stage1"]
    loci_rows = read_loci_tsv(loci_path)
    orthogroups = read_orthogroups_tsv(orthogroups_path)

    genome_data: dict[str, dict] = {}
    for gid in sorted(base["genomes"]):
        split = source["splits"][gid]["split"]
        if split not in ("calibration", "test"):
            raise ValueError(f"{gid}: only calibration and test splits are admitted")
        output = base["genomes"][gid]["path"]
        locus_file = output + ".anomalies.locus.jsonl"
        if locus_file not in base.get("artifacts", {}).get(gid, {}):
            raise ValueError(f"{gid}: locus anomaly export not covered by production receipt")
        records = [json.loads(line) for line in Path(locus_file).read_text(encoding="utf-8").splitlines() if line.strip()]
        predictions = {}
        for mode, r in runs.items():
            if r.get("_abstained"):
                predictions[mode] = {}
            else:
                predictions[mode] = runner.parse_integrated(r["genomes"][gid]["path"], "ec")
        genome_data[gid] = {
            "split": split,
            "locus_records": records,
            "predictions": predictions,
        }

    receipts = {name: runner.hash_record(path) for name, path in source.get("runs", {}).items()}
    for name, info in (source.get("abstained_arms") or {}).items():
        if info.get("run_record"):
            receipts[name] = runner.hash_record(info["run_record"])
    shared_wall = sum(float(r.get("runtime_sec", 0.0)) for r in runs.values())

    return assemble_universe_b(
        genome_data=genome_data,
        loci_rows=loci_rows,
        orthogroups=orthogroups,
        decoys_per_case=decoys_per_case,
        seed=seed,
        source_sha256=E.digest(source),
        source_receipts=receipts,
        shared_production_wall_seconds=shared_wall,
        limitations=source.get("limitations", []),
    )


def write_public_hashes(public_files: list[Path], out_path: Path) -> list[tuple[str, str]]:
    """Compute sha256 of each public file and write '<hex>  <path>' lines immutably."""
    if not public_files:
        raise ValueError("at least one public file is required")
    rows = []
    for p in public_files:
        if not Path(p).is_file():
            raise FileNotFoundError(f"public file not found: {p}")
        rows.append((runner.sha256_file(str(p)), str(p)))
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("x", encoding="utf-8") as fh:
        for digest_hex, path_str in rows:
            fh.write(f"{digest_hex}  {path_str}\n")
    return rows


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="cmd", required=True)

    ua = sub.add_parser("universe-a", help="Build Universe A target-centred public bundle")
    ua.add_argument("--sources", required=True, type=Path)
    ua.add_argument("--out", required=True, type=Path)
    ua.add_argument("--cap-per-genome", type=int, default=None)

    ub = sub.add_parser("universe-b", help="Build Universe B locus-anchored public and private bundles")
    ub.add_argument("--sources", required=True, type=Path)
    ub.add_argument("--loci", required=True, type=Path)
    ub.add_argument("--orthogroups", required=True, type=Path)
    ub.add_argument("--decoys-per-case", type=int, default=2)
    ub.add_argument("--seed", type=int, default=20261102)
    ub.add_argument("--out-public", required=True, type=Path)
    ub.add_argument("--out-private", required=True, type=Path)

    ha = sub.add_parser("hash", help="Record SHA-256 of frozen public bundles")
    ha.add_argument("--public", required=True, nargs="+", type=Path)
    ha.add_argument("--out", required=True, type=Path)

    args = parser.parse_args(argv)
    if args.cmd == "universe-a":
        bundle = build_universe_a(E.read(args.sources), cap_per_genome=args.cap_per_genome)
        E.write(args.out, bundle)
    elif args.cmd == "universe-b":
        pub, priv = build_universe_b(
            source=E.read(args.sources),
            loci_path=args.loci,
            orthogroups_path=args.orthogroups,
            decoys_per_case=args.decoys_per_case,
            seed=args.seed,
        )
        E.write(args.out_public, pub)
        E.write(args.out_private, priv)
    elif args.cmd == "hash":
        write_public_hashes(args.public, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
