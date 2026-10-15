#!/usr/bin/env python3
"""
Pillar 4 Benchmark Dataset Builder (De-identified Opaque-Identifier Edition).

Constructs:
1. results/pillar4_agent_benchmark_120.jsonl:
   - 120 structured genome-scale anomalies across 4 categories:
     - 30 PATHWAY_GAP (18 TP Promotable, 12 TN Decoys)
     - 30 INCOMPLETE_COMPLEX (18 TP Promotable, 12 TN Decoys)
     - 30 BORDERLINE_ESSENTIAL (18 TP Promotable, 12 TN Decoys)
     - 30 BORDERLINE_CONFLICT (15 Real Conflicts, 15 Cross-Clade Constraint Traps)
   - Spanning 4 prokaryotes (ecoli, bsubtilis, mtb, paeruginosa) and 3 eukaryotes (scerevisiae, celegans, athaliana).
   - Strictly de-identified and constructed uniformly from pre-integrated proteome candidate evidence:
     * All gene symbols are masked ("uncharacterized").
     * Every candidate in `candidate_pids` is populated by the same symmetric lookup
       `_build_candidate_from_proteome_store(...)` without any truth-conditioned branching.
     * Candidate order is sorted lexicographically by `protein_id`.

2. results/pillar4_lro_50_cases.jsonl:
   - 50 unique-protein Leave-Reaction-Out (LRO) metabolic gap cases from `cases_mg1655_all.tsv`
     and `runs_rlgc_all/protein/shard_*/case_*/suggestions.tsv` (matching `parts[3] == go_id`),
     where all 50 cases naturally contain `true_protein_id` in the top-5 RLGC suggestions
     without forced insertion, and all candidate features are populated symmetrically by
     `_build_candidate_from_proteome_store(...)`.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path
from typing import Any, Dict, List


REPO_ROOT = Path(__file__).resolve().parent
while REPO_ROOT != REPO_ROOT.parent and not (REPO_ROOT / "gspa-core").exists():
    REPO_ROOT = REPO_ROOT.parent
if not (REPO_ROOT / "gspa-core").exists():
    REPO_ROOT = Path("/tmp")
DEFAULT_RESULTS_DIR = REPO_ROOT / "results" if (REPO_ROOT / "gspa-core").exists() else Path("/tmp/results")
REF_DIR = REPO_ROOT / "benchmark" / "manifests" / "reference_data"
PROTEOMES_DIR = Path("/data/hohndor/gspa/proteomes")


def _prob_to_logodds(p: float) -> float:
    p_c = max(1e-4, min(1.0 - 1e-4, p))
    return round(math.log(p_c / (1.0 - p_c)), 4)


def _load_proteome_bundle() -> Dict[str, Any]:
    bundle_path = REF_DIR / "pillar4_proteome_candidates.json"
    if not bundle_path.exists():
        raise FileNotFoundError(f"Missing proteome candidate bundle: {bundle_path}")
    return json.loads(bundle_path.read_text(encoding="utf-8"))


def _build_candidate_from_proteome_store(
    store_key: str,
    pid: str,
    cand_store: Dict[str, Dict[str, Any]],
    fallback_q: float = 0.05,
) -> Dict[str, Any]:
    """
    Uniform, label-blind candidate record builder.
    Looks up pre-integrated proteome evidence for `pid` by `store_key` without inspecting ground truth.
    """
    raw_entry = cand_store.get(store_key)
    if raw_entry is not None:
        rec = dict(raw_entry)
        rec.pop("real_unimatrix_esm2_cos", None)
        rec.pop("real_co_strand_neighbors", None)
        rec["protein_id"] = pid
        rec["gene_symbol"] = "uncharacterized"
        return rec

    q_val = round(float(fallback_q), 4)
    return {
        "protein_id": pid,
        "gene_symbol": "uncharacterized",
        "description": f"Uncharacterized protein {pid}",
        "current_posterior": q_val,
        "likelihood_logodds": _prob_to_logodds(max(0.02, min(0.95, q_val * 5.0))),
        "final_logodds": _prob_to_logodds(max(0.02, min(0.95, q_val * 5.5))),
        "prior_contributions": {"rlgc_prior": q_val},
        "supporting_claims": [
            {
                "source": "rlgc_suggester",
                "evidence_type": "GENOMIC_CONTEXT",
                "raw_score": q_val,
                "calibrated_prob": q_val,
            }
        ],
        "operon_id": "op_unassigned",
        "in_pathway_operon": False,
        "operon_neighbors": {},
        "other_annotations": [],
        "rlgc_bayes_factor": 1.0,
        "suggestion_score": q_val,
    }


def build_lro_50_cases(bundle: Dict[str, Any]) -> List[Dict[str, Any]]:
    """
    Builds the 50 unique-protein challenging LRO metabolic gap testset.
    Reads `suggestions.tsv` (matching column 3 `function_id == go_id`) when mounted, or the
    extracted `lro_50_manifest` from `pillar4_proteome_candidates.json`.
    Verifies that `true_pid` is naturally present in `top5_pids` without insertion, and populates
    every candidate symmetrically via `_build_candidate_from_proteome_store`.
    """
    local_lro_dir = REF_DIR / "lro_data"
    lro_dir = PROTEOMES_DIR / "bench_ecoli" / "leave_reaction_out"
    if not lro_dir.exists() and local_lro_dir.exists():
        lro_dir = local_lro_dir
    cases_path = lro_dir / "cases_mg1655_all.tsv"
    runs_dir = (PROTEOMES_DIR / "bench_ecoli" / "leave_reaction_out") / "runs_rlgc_all" / "protein"

    cases_meta: Dict[int, Dict[str, str]] = {}
    if cases_path.exists():
        with open(cases_path, "r", encoding="utf-8") as fh:
            for idx, row in enumerate(csv.DictReader(fh, delimiter="\t")):
                c_idx = int(row.get("case_idx", idx))
                cases_meta[c_idx] = row

    cand_store = bundle["proteome_candidates"]
    lro_manifest = bundle["lro_50_manifest"]

    lro_anomalies: List[Dict[str, Any]] = []
    for seq_num, entry in enumerate(lro_manifest, start=1):
        c_idx = int(entry["case_idx"])
        true_pid = entry["true_protein_id"]
        rxn_id = entry["reaction_id"]
        go_id = entry["go_term"]
        meta = cases_meta.get(c_idx, {})

        cand_pids: List[str] = []
        cand_qs: List[float] = []
        sug_path = runs_dir / f"shard_{c_idx // 85}" / f"case_{c_idx:04d}" / "suggestions.tsv"
        if sug_path.exists():
            with open(sug_path, "r", encoding="utf-8") as fh:
                for line in fh:
                    if line.startswith("kind\t") or line.startswith("#") or not line.strip():
                        continue
                    parts = line.strip().split("\t")
                    if len(parts) >= 10 and parts[3].strip() == go_id and parts[8]:
                        cand_pids = [p.strip() for p in parts[8].split(",") if p.strip()]
                        cand_qs = [float(q.strip()) for q in parts[9].split(",") if q.strip()]
                        break
        if not cand_pids:
            cand_pids = list(entry["rlgc_ordered_pids"])
            cand_qs = list(entry["rlgc_ordered_qs"])

        top5_pids = list(cand_pids[:5])
        top5_qs = list(cand_qs[:5]) if len(cand_qs) >= len(top5_pids) else [0.05] * len(top5_pids)
        if true_pid not in top5_pids:
            raise ValueError(f"Expected {true_pid} in natural top-5 RLGC suggestions for case {c_idx}, got {top5_pids}")

        rlgc_pos = top5_pids.index(true_pid) + 1
        candidates: List[Dict[str, Any]] = []
        for pid, q_val in zip(top5_pids, top5_qs):
            candidates.append(
                _build_candidate_from_proteome_store(
                    store_key=f"LRO|{c_idx}|{pid}",
                    pid=pid,
                    cand_store=cand_store,
                    fallback_q=q_val,
                )
            )

        candidates_sorted = sorted(candidates, key=lambda x: x["protein_id"])
        lro_anomalies.append({
            "anomaly_id": entry.get("anomaly_id", f"LRO_GAP_{seq_num:03d}"),
            "anomaly_type": "PATHWAY_GAP",
            "organism_taxon": "NCBITaxon_562",
            "target_function_id": go_id,
            "target_function_type": "GO",
            "target_function_label": meta.get("reaction_id") or rxn_id,
            "go_aspect": "MF",
            "true_protein_id": true_pid,
            "rlgc_ordered_pids": list(top5_pids),
            "rlgc_baseline_rank": rlgc_pos - 1,
            "expected_action": "PROMOTE",
            "expected_protein_id": true_pid,
            "expected_function_id": go_id,
            "context": dict(entry["context"]),
            "candidates": candidates_sorted,
        })

    return lro_anomalies


def build_benchmark_120_anomalies(bundle: Dict[str, Any]) -> List[Dict[str, Any]]:
    """
    Builds the 120-anomaly benchmark across 4 categories by populating each anomaly's
    `candidate_pids` uniformly from `proteome_candidates` without truth-conditioned branching.
    """
    cand_store = bundle["proteome_candidates"]
    anom_manifest = bundle["anom_120_manifest"]

    anomalies: List[Dict[str, Any]] = []
    for spec in anom_manifest:
        aid = spec["anomaly_id"]
        cands = [
            _build_candidate_from_proteome_store(
                store_key=f"ANOM|{aid}|{pid}",
                pid=pid,
                cand_store=cand_store,
            )
            for pid in spec["candidate_pids"]
        ]
        cands_sorted = sorted(cands, key=lambda x: x["protein_id"])
        anom_obj = {k: v for k, v in spec.items() if k != "candidate_pids"}
        anom_obj["candidates"] = cands_sorted
        anomalies.append(anom_obj)

    return anomalies


def main() -> int:
    parser = argparse.ArgumentParser(description="Pillar 4 Dataset Builder")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_RESULTS_DIR, help="Output directory")
    args = parser.parse_args()

    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"[pillar4-builder] Generating de-identified Pillar 4 Benchmark Datasets in {out_dir}...")

    bundle = _load_proteome_bundle()
    anomalies_120 = build_benchmark_120_anomalies(bundle)
    out_120 = out_dir / "pillar4_agent_benchmark_120.jsonl"
    with open(out_120, "w", encoding="utf-8") as fh:
        for a in anomalies_120:
            fh.write(json.dumps(a) + "\n")
    print(f"  Wrote {len(anomalies_120)} anomalies to {out_120}")

    lro_50 = build_lro_50_cases(bundle)
    out_lro = out_dir / "pillar4_lro_50_cases.jsonl"
    with open(out_lro, "w", encoding="utf-8") as fh:
        for c in lro_50:
            fh.write(json.dumps(c) + "\n")
    print(f"  Wrote {len(lro_50)} unique-protein LRO gap cases to {out_lro}")
    unique_pids = len(set(c["true_protein_id"] for c in lro_50))
    print(f"  Verified unique true_protein_id count in LRO-50: {unique_pids}/{len(lro_50)}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
