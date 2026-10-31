#!/usr/bin/env python3
"""Experiments 4 and 5: re-attribute the stored Stage 2B cells with a synteny covariate and candidate limits (design §4, §5).

For each cell the stored sidecar payload (`fba_input.json`) and repair result (`fba_output.json`) are reused; only the
backward locus responsibility is recomputed with the production sidecar code (imported from
`scripts/stoichiometric_locus_refiner.py`) under the arms below. No linear program is solved. Each arm writes a
sparse update record (`<index>.<arm>.json`: the (protein, enzyme) pairs whose probability changes by more than
1e-5, as `LocusResponsibilityUpdater.applyUpdatedScores` writes them) plus diagnostics.

Arms:
  replay      empty covariates, original candidate floor      (must reproduce the stored stage 2B output)
  syn         synteny covariate from layout + KEGG modules     (experiment 4, S2B-syn)
  syn_k       synteny covariate + min score 0.05 + top 10     (experiment 4, S2B-syn-k)
  syn_shuf1   synteny covariate from a shuffled layout, seed 1 (control)
  K1          candidate min score 0.05                        (experiment 5)
  K2          candidate top 10 per enzyme number              (experiment 5)
  K3          both                                            (experiment 5)
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import importlib.util
import json
import os
from pathlib import Path
import random
import subprocess
import sys
import tempfile
import time
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

from .common import CANDIDATE_FLOOR, SYNTENY_WINDOW, NEIGHBOUR_MIN_SCORE, design_receipt, sha256_file, write_json
from .price_score import locate_cell

ARMS: Dict[str, Dict[str, Any]] = {
    "replay": {"synteny": False, "min_score": CANDIDATE_FLOOR, "top_k": 0},
    "syn": {"synteny": True, "min_score": CANDIDATE_FLOOR, "top_k": 0},
    "syn_k": {"synteny": True, "min_score": 0.05, "top_k": 10},
    "syn_shuf1": {"synteny": True, "min_score": CANDIDATE_FLOOR, "top_k": 0, "shuffle_layout_seed": 1},
    "K1": {"synteny": False, "min_score": 0.05, "top_k": 0},
    "K2": {"synteny": False, "min_score": CANDIDATE_FLOOR, "top_k": 10},
    "K3": {"synteny": False, "min_score": 0.05, "top_k": 10},
}
REPLAY_TOL = 1e-6


def load_sidecar(path: Path):
    spec = importlib.util.spec_from_file_location("stoichiometric_locus_refiner", str(path))
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def shuffle_layout(layout: Mapping[str, Tuple[str, int, int]], seed: int) -> Dict[str, Tuple[str, int, int]]:
    ids = sorted(layout)
    positions = [layout[i] for i in ids]
    rng = random.Random(seed)
    rng.shuffle(ids)
    return dict(zip(ids, positions))


def diff_scores(before: Mapping[str, Mapping[str, float]], after: Mapping[str, Mapping[str, float]], tol: float = 1e-5) -> Dict[str, Dict[str, float]]:
    out: Dict[str, Dict[str, float]] = {}
    for gid, ec_map in after.items():
        base = before.get(gid, {})
        for ec, v in ec_map.items():
            if abs(float(v) - float(base.get(ec, 0.0))) > tol:
                out.setdefault(gid, {})[ec] = float(v)
    return out


def top_calls(scores: Mapping[str, Mapping[str, float]]) -> Dict[str, Tuple[str, float]]:
    out: Dict[str, Tuple[str, float]] = {}
    for gid, ec_map in scores.items():
        if ec_map:
            best_score = max(float(s) for s in ec_map.values())
            # deterministic: highest score, then lexically smallest enzyme number
            best_ec = min(e for e, s in ec_map.items() if float(s) == best_score)
            out[gid] = (best_ec, best_score)
    return out


def top_call_changes(before: Mapping[str, Mapping[str, float]], after: Mapping[str, Mapping[str, float]], floor: float) -> Dict[str, int]:
    tb, ta = top_calls(before), top_calls(after)
    changed = tie_artefacts = 0
    for gid, (ec_b, s_b) in tb.items():
        ec_a, s_a = ta.get(gid, (ec_b, s_b))
        if ec_a != ec_b:
            changed += 1
            if s_b <= floor + 1e-9:
                tie_artefacts += 1
    return {"top_call_changes": changed, "tie_artefact_top_call_changes": tie_artefacts}


def run_cell(
    cell: Mapping[str, Any],
    located: Mapping[str, Any],
    sidecar: Path,
    fba_data_dir: Path,
    layout_tsv: Path,
    kegg_tsv: Path,
    out_dir: Path,
    arms: Sequence[str],
    cli_check: bool,
    python_bin: str,
) -> Dict[str, Any]:
    t0 = time.time()
    slr = load_sidecar(sidecar)
    met = Path(located["metabolism_dir"])
    payload = json.loads((met / "fba_input.json").read_text(encoding="utf-8"))
    stored = json.loads((met / "fba_output.json").read_text(encoding="utf-8"))
    stored_updated = stored.pop("updated_gene_ec_scores", {})
    repair = stored
    gene_ec_scores: Dict[str, Dict[str, float]] = payload["gene_ec_scores"]
    floor = min((float(s) for m in gene_ec_scores.values() for s in m.values()), default=0.0)
    n_scores = sum(len(m) for m in gene_ec_scores.values())
    refiner = slr.StoichiometricLocusRefiner(fba_data_dir)
    agg = slr.SparseRoleAggregator(top_m=int(payload.get("top_m", 3)), tau_sub=float(payload.get("tau_sub", 0.25)))
    rxn_probs = payload.get("rxn_probs")
    if not rxn_probs:
        ec_probs, _ = agg.aggregate_ec_scores(gene_ec_scores)
        rxn_probs = agg.compute_reaction_probabilities(ec_probs, refiner.r2ec, refiner.rxn_meta)
    layout = slr.load_layout_tsv(layout_tsv)
    ec_modules = slr.load_kegg_modules_tsv(kegg_tsv)
    diagnostics: Dict[str, Any] = {
        "cell": dict(cell), "n_scores": n_scores, "floor_score": floor, "n_repaired_rxns": len(repair.get("repaired_rxns", [])),
        "stored_changed_scores": len([1 for gid, m in stored_updated.items() for ec, v in m.items() if abs(float(v) - float(gene_ec_scores.get(gid, {}).get(ec, 0.0))) > 1e-5]),
        "arms": {},
    }
    for arm in arms:
        spec = ARMS[arm]
        operon_map: Dict[str, float] = {}
        if spec["synteny"]:
            lay = shuffle_layout(layout, spec["shuffle_layout_seed"]) if spec.get("shuffle_layout_seed") else layout
            operon_map = slr.compute_synteny_covariate(gene_ec_scores, lay, ec_modules, window=SYNTENY_WINDOW, min_score=NEIGHBOUR_MIN_SCORE)
        updated = refiner.attribute_locus_responsibility(
            repair, gene_ec_scores, rxn_probs,
            rlgc_bf_by_gene_ec=payload.get("rlgc_bf_by_gene_ec") or None,
            plm_sim_by_gene_ec=payload.get("plm_sim_by_gene_ec") or None,
            operon_partner_boost_by_gene_ec=operon_map or None,
            candidate_min_score=float(spec["min_score"]), candidate_top_k=int(spec["top_k"]),
        )
        updates = diff_scores(gene_ec_scores, updated)
        arm_diag: Dict[str, Any] = {
            "n_changed_scores": sum(len(m) for m in updates.values()),
            "n_proteins_changed": len(updates),
            "n_operon_covariates": len(operon_map),
            "candidate_min_score": float(spec["min_score"]), "candidate_top_k": int(spec["top_k"]),
            **top_call_changes(gene_ec_scores, updated, floor),
        }
        if arm == "replay":
            max_diff = 0.0
            for gid, m in updated.items():
                sm = stored_updated.get(gid, {})
                for ec, v in m.items():
                    max_diff = max(max_diff, abs(float(v) - float(sm.get(ec, gene_ec_scores[gid].get(ec, 0.0)))))
            arm_diag["replay_max_abs_diff"] = max_diff
            arm_diag["replay_ok"] = max_diff <= REPLAY_TOL
        write_json(out_dir / f"{cell['index']}.{arm}.json", {"cell": dict(cell), "arm": arm, "updates": updates, "diagnostics": arm_diag})
        diagnostics["arms"][arm] = arm_diag
        del updated
    if cli_check:
        with tempfile.TemporaryDirectory(prefix="slr-cli-") as td:
            out_json = Path(td) / "out.json"
            cmd = [python_bin, str(sidecar), "--fba-data-dir", str(fba_data_dir), "--input-json", str(met / "fba_input.json"), "--output-json", str(out_json),
                   "--repair-json", str(met / "fba_output.json"), "--layout-tsv", str(layout_tsv), "--kegg-modules-tsv", str(kegg_tsv),
                   "--candidate-min-score", "0.05", "--candidate-top-k", "10", "--gram", payload.get("gram", "neg"), "--medium", payload.get("medium", "minimal")]
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            cli_out = json.loads(out_json.read_text())
            cli_updates = diff_scores(gene_ec_scores, cli_out.get("updated_gene_ec_scores", {}))
            ref = json.loads((out_dir / f"{cell['index']}.syn_k.json").read_text())["updates"] if "syn_k" in arms else None
            same = (ref is not None and all(abs(cli_updates.get(g, {}).get(e, -1) - v) <= 1e-9 for g, m in ref.items() for e, v in m.items()) and sum(len(m) for m in cli_updates.values()) == sum(len(m) for m in ref.values()))
            diagnostics["cli_check"] = {"command": cmd, "n_changed_scores": sum(len(m) for m in cli_updates.values()), "matches_function_level_syn_k": bool(same),
                                        "attribution_options": cli_out.get("attribution_options"), "synteny_covariate": cli_out.get("synteny_covariate")}
    diagnostics["seconds"] = time.time() - t0
    write_json(out_dir / f"{cell['index']}.diagnostics.json", diagnostics)
    return diagnostics


def _worker(task: Tuple[Dict[str, Any], Dict[str, str]]) -> Dict[str, Any]:
    cell, a = task
    collection = Path(a["collection"])
    located = locate_cell(collection, cell)
    layout_tsv = Path(a["layout_dir"]) / f"{cell['assembly']}.layout.tsv"
    cli_check = a["cli_check_predictors"] and cell["predictor"] in a["cli_check_predictors"].split(",")
    try:
        d = run_cell(cell, located, Path(a["sidecar"]), Path(a["fba_data_dir"]), layout_tsv, Path(a["kegg_modules"]), Path(a["out_dir"]), a["arms"].split(","), bool(cli_check), a["python_bin"])
        print(json.dumps({"cell": cell["index"], "predictor": cell["predictor"], "assembly": cell["assembly"], "seconds": round(d["seconds"], 1),
                          "replay_ok": d["arms"].get("replay", {}).get("replay_ok"), "changed": {k: v["n_changed_scores"] for k, v in d["arms"].items()}}), flush=True)
        return {"cell": cell["index"], "status": "ok"}
    except Exception as exc:  # keep the failure visible, continue with other cells
        err = {"cell": cell["index"], "status": "failed", "error": repr(exc)}
        write_json(Path(a["out_dir"]) / f"{cell['index']}.failed.json", err)
        print(json.dumps(err), flush=True)
        return err


def parse_cells(spec: str, n: int = 132) -> List[int]:
    out: List[int] = []
    for part in spec.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            lo, hi = part.split("-")
            out.extend(range(int(lo), int(hi) + 1))
        else:
            out.append(int(part))
    return [i for i in out if 0 <= i < n]


def main(argv: Optional[Sequence[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--collection", type=Path, required=True)
    p.add_argument("--sidecar", type=Path, required=True)
    p.add_argument("--fba-data-dir", type=Path, required=True)
    p.add_argument("--layout-dir", type=Path, required=True)
    p.add_argument("--kegg-modules", type=Path, required=True)
    p.add_argument("--cells", default="0-131")
    p.add_argument("--arms", default=",".join(ARMS))
    p.add_argument("--cli-check-predictors", default="mapred")
    p.add_argument("--python-bin", default=sys.executable)
    p.add_argument("--workers", type=int, default=4)
    p.add_argument("--out-dir", type=Path, required=True)
    a = p.parse_args(argv)
    plan = json.loads((a.collection / "plan.json").read_text())
    cells = [c for c in plan["cells"] if c["index"] in set(parse_cells(a.cells))]
    a.out_dir.mkdir(parents=True, exist_ok=True)
    args = {"collection": str(a.collection), "sidecar": str(a.sidecar), "fba_data_dir": str(a.fba_data_dir), "layout_dir": str(a.layout_dir),
            "kegg_modules": str(a.kegg_modules), "out_dir": str(a.out_dir), "arms": a.arms, "cli_check_predictors": a.cli_check_predictors, "python_bin": a.python_bin}
    write_json(a.out_dir / f"receipt-{os.environ.get('SLURM_JOB_ID', 'local')}-{cells[0]['index'] if cells else 'none'}.json",
               {**design_receipt(), "sidecar_sha256": sha256_file(a.sidecar), "kegg_modules_sha256": sha256_file(a.kegg_modules), "arms": {k: ARMS[k] for k in a.arms.split(",")}, "cells": [c["index"] for c in cells]})
    tasks = [(c, args) for c in cells]
    if a.workers > 1 and len(tasks) > 1:
        with ProcessPoolExecutor(max_workers=a.workers) as pool:
            results = list(pool.map(_worker, tasks))
    else:
        results = [_worker(t) for t in tasks]
    failed = [r for r in results if r.get("status") != "ok"]
    print(json.dumps({"cells": len(results), "failed": len(failed)}))
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
