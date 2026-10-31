#!/usr/bin/env python3
"""Price scoring of sparse update records against stage 0 with the frozen metric functions (design §0).

A variant arm is stage 0 with the recorded updates applied (`posterior_prob` replaced for the listed
(protein, enzyme number) pairs). The Price metrics read only the 148 reference proteins' scores
(`evaluate_price_reference.read_scores` retains reference proteins only), so applying the sparse updates to the
retained stage-0 scores is identical to scoring a rewritten output file. Every metric function is imported from the
frozen engine (`price_reference_metrics`, `ec_identity`, `v3_metrics`) and the Stage 2B evaluator
(`cafa_stats`, `cafa_fmax`, `project_genome_pairs`, `summarize_corrections_and_harms`).
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
import traceback
from typing import Any, Dict, List, Mapping, Optional, Sequence

V3_DIR = Path(__file__).resolve().parents[1]
if str(V3_DIR) not in sys.path:
    sys.path.insert(0, str(V3_DIR))

import ec_identity as E  # noqa: E402
import evaluate_price_reference as R  # noqa: E402
import price_reference_metrics as M  # noqa: E402
import v3_metrics as V3M  # noqa: E402
from stage2b.evaluate_price_stage2b import (  # noqa: E402
    CAFA_VEC_LEN,
    cafa_fmax,
    cafa_stats,
    project_genome_pairs,
    sha,
    summarize_corrections_and_harms,
)
from stage2b.run_price_stage2b_collection import ARMS  # noqa: E402

# The scoring container has no cobra, so the design receipt is computed here rather than imported from .common
DESIGN_PATH = V3_DIR.parents[1] / "tasks" / "v3-stage2b-followups-design-2026-10-31.md"
DESIGN_SHA256 = "5a468ece2c892f3cc6b47622f8ef77f444d1e40f04bb11ca6b38ac5779ecb768"


DESIGN_BODY_SHA256 = "25a1a8676078e6e32d7c30120bb4024d7630c26509e18253ffcffdb612ef5f14"  # everything before "## 9. Deviations"


def design_receipt() -> Dict[str, Any]:
    current = sha(DESIGN_PATH) if DESIGN_PATH.exists() else None
    body = None
    if DESIGN_PATH.exists():
        raw = DESIGN_PATH.read_bytes()
        marker = b"## 9. Deviations"
        body = hashlib.sha256(raw[: raw.index(marker)] if marker in raw else raw).hexdigest()
    return {"design_path": str(DESIGN_PATH), "design_sha256_frozen": DESIGN_SHA256, "design_sha256_current": current, "design_unchanged": current == DESIGN_SHA256,
            "design_body_sha256_frozen": DESIGN_BODY_SHA256, "design_body_sha256_current": body, "design_body_unchanged": body == DESIGN_BODY_SHA256}

EC_AUDIT_SHA256 = "5302494b16ecd3f76fd96d24e48e27c298206b27608be72c5c66f1d197f4f79f"


def locate_cell(collection: Path, cell: Mapping[str, Any]) -> Dict[str, Any]:
    """The uniquely complete attempt of a cell: production stage output paths, their recorded digests and metabolism dir."""
    done_paths = list((collection / "attempts" / str(cell["index"])).glob("*/completed.json"))
    if len(done_paths) != 1:
        raise ValueError(f"cell {cell['index']} is not uniquely complete")
    directory = done_paths[0].parent
    verification = json.loads((directory / "production/verification.json").read_text())
    if verification["status"] != "verified" or set(verification["runs"]) != set(ARMS):
        raise ValueError(f"cell {cell['index']}: invalid verification")
    outputs: Dict[str, Dict[str, str]] = {}
    for arm in ARMS:
        receipt = json.loads((directory / "production" / arm / "verified_run.json").read_text())
        output = receipt["spec"]["outputs"][cell["assembly"]]
        outputs[arm] = {"path": output, "sha256": receipt["artifacts"][output]}
    return {
        "directory": str(directory),
        "outputs": outputs,
        "metabolism_dir": str(Path(outputs["stage2b"]["path"]).with_name(Path(outputs["stage2b"]["path"]).name + ".metabolism")),
    }


def load_updates(path: Path) -> Dict[str, Dict[str, float]]:
    obj = json.loads(Path(path).read_text())
    return {str(p): {M.ec(e): float(v) for e, v in ecs.items()} for p, ecs in obj.get("updates", {}).items()}


def apply_updates(retained: Mapping[str, Mapping[str, float]], updates: Mapping[str, Mapping[str, float]]) -> Dict[str, Dict[str, float]]:
    out = {pid: dict(scores) for pid, scores in retained.items()}
    n = 0
    for pid, ecs in updates.items():
        if pid not in out:
            continue
        for ec, v in ecs.items():
            out[pid][ec] = M.score(v)
            n += 1
    return out


def score_cells(
    collection: Path,
    manifest: Path,
    audit_path: Path,
    predictors: Sequence[str],
    updates_dir: Path,
    arms: Sequence[str],
    out: Path,
    *,
    n_boot: int = 2000,
    seed: int = 13,
    comparisons: Optional[Sequence[Sequence[str]]] = None,
) -> Dict[str, Any]:
    collection, manifest, audit_path, updates_dir, out = (Path(p).resolve() for p in (collection, manifest, audit_path, updates_dir, out))
    if out.exists():
        raise ValueError("refusing to overwrite evaluation")
    if sha(audit_path) != EC_AUDIT_SHA256:
        raise ValueError("EC identity audit differs from the registered audit")
    plan = json.loads((collection / "plan.json").read_text())
    if tuple(plan.get("arms", ())) != ARMS or len(plan["cells"]) != 132:
        raise ValueError("not the registered 132-cell Stage 2B collection")
    genomes = {c["assembly"] for c in plan["cells"]}
    refs, exclusions = R.references(manifest, genomes)
    audit = json.loads(audit_path.read_text())
    entries = audit["entries"]
    E.validate_transfer_graph(entries)
    selected = [c for c in plan["cells"] if c["predictor"] in predictors]
    if len(selected) != len(genomes) * len(predictors):
        raise ValueError("incomplete predictor/genome panel")
    eval_arms = ["stage0", "stage2b", *arms]
    modes = ("exact", "current_unambiguous")
    rows: Dict[str, Dict[str, Dict[str, Dict[str, Dict[str, Any]]]]] = {m: {p: {a: {} for a in eval_arms} for p in predictors} for m in modes}
    projected_pairs: Dict[str, Dict[str, Dict[str, Dict[str, Dict[Any, float]]]]] = {m: {p: {a: {} for a in eval_arms} for p in predictors} for m in modes}
    projected_truth: Dict[str, Dict[str, set]] = {m: {} for m in modes}
    out.mkdir(parents=True)
    applied_counts: Dict[str, Dict[str, int]] = {}
    bound: Dict[str, str] = {}
    try:
        for cell in selected:
            gid, pred, idx = cell["assembly"], cell["predictor"], cell["index"]
            located = locate_cell(collection, cell)
            targets = {pid: terms for (g, pid), terms in refs.items() if g == gid}
            retained: Dict[str, Dict[str, Dict[str, float]]] = {}
            for arm in ARMS:
                path = located["outputs"][arm]["path"]
                retained[arm] = R.read_scores(path, targets, raw=False, predictor=pred, expected_rows=cell["claims"], expected_sha=located["outputs"][arm]["sha256"])
                bound[path] = located["outputs"][arm]["sha256"]
            for arm in arms:
                upd_path = updates_dir / f"{idx}.{arm}.json"
                if not upd_path.exists():
                    raise FileNotFoundError(f"missing updates for cell {idx} arm {arm}: {upd_path}")
                upd = load_updates(upd_path)
                retained[arm] = apply_updates(retained["stage0"], upd)
                applied_counts.setdefault(arm, {})[str(idx)] = sum(1 for pid in upd if pid in targets)
                bound[str(upd_path)] = sha(upd_path)
            for arm in eval_arms:
                values = retained[arm]
                for mode in modes:
                    rows[mode][pred][arm][gid] = {pid: M.protein(terms, values[pid], entries, mode, integrated=True) for pid, terms in targets.items()}
                    t_pairs, p_pairs = project_genome_pairs(targets, values, entries, mode)
                    projected_truth[mode][gid] = t_pairs
                    projected_pairs[mode][pred][arm][gid] = p_pairs
            print(json.dumps(dict(cell=idx, predictor=pred, assembly=gid, status="scored")), flush=True)

        pairs_to_compare = list(comparisons) if comparisons else [("stage0", "stage2b"), ("stage0", "stage0")] + [("stage0", a) for a in arms] + [("stage2b", a) for a in arms]
        summaries: Dict[str, Dict[str, Any]] = {}
        for mode, panel in rows.items():
            summaries[mode] = {}
            for pred, arms_map in panel.items():
                results: Dict[str, Any] = {}
                top1_credit_vectors = {a: {g: [r["top1_credit"] for _, r in sorted(records.items())] for g, records in by.items()} for a, by in arms_map.items()}
                micro_vecs: Dict[str, Dict[str, List[float]]] = {a: {} for a in eval_arms}
                cafa_vecs: Dict[str, Dict[str, List[float]]] = {a: {} for a in eval_arms}
                for g in sorted(genomes):
                    t_g = projected_truth[mode][g]
                    univ_g = V3M.candidate_universe(t_g, *(projected_pairs[mode][pred][a][g] for a in eval_arms))
                    for a in eval_arms:
                        p_g = projected_pairs[mode][pred][a][g]
                        micro_vecs[a][g] = V3M.pair_stats(p_g, t_g, universe=univ_g)
                        cafa_vecs[a][g] = cafa_stats(p_g, t_g)
                for arm, by in arms_map.items():
                    s = M.summary([r for records in by.values() for r in records.values()])
                    s["top1_genome_bootstrap"] = M.bootstrap_top1(top1_credit_vectors[arm], repetitions=n_boot, seed=seed)
                    tm = V3M.sum_vectors(micro_vecs[arm].values())
                    s["micro_fmax"] = {**V3M.fmax(tm), "genome_bootstrap": V3M.genome_bootstrap(micro_vecs[arm], lambda v: V3M.fmax(v)["fmax"], n_boot=n_boot, seed=seed)}
                    tc = V3M.sum_vectors(cafa_vecs[arm].values())
                    s["cafa_fmax"] = {**cafa_fmax(tc), "genome_bootstrap": V3M.genome_bootstrap(cafa_vecs[arm], lambda v: cafa_fmax(v)["fmax"], n_boot=n_boot, seed=seed)}
                    results[arm] = s
                paired: Dict[str, Any] = {}
                for before, after in pairs_to_compare:
                    flat = lambda a: {(g, p): r for g, records in arms_map[a].items() for p, r in records.items()}
                    pm = {g: micro_vecs[after][g] + micro_vecs[before][g] for g in sorted(genomes)}
                    pc = {g: cafa_vecs[after][g] + cafa_vecs[before][g] for g in sorted(genomes)}
                    ch = summarize_corrections_and_harms(projected_pairs[mode][pred][before], projected_pairs[mode][pred][after], projected_truth[mode], arms_map[before], arms_map[after], tau=0.5, n_boot=n_boot, seed=seed)
                    paired[f"{before}_to_{after}"] = dict(
                        **M.paired(flat(before), flat(after)),
                        top1_delta_bootstrap=M.bootstrap_top1(top1_credit_vectors[after], baseline=top1_credit_vectors[before], repetitions=n_boot, seed=seed),
                        micro_fmax_delta_bootstrap=V3M.genome_bootstrap(pm, V3M.paired_stat(lambda v: V3M.fmax(v)["fmax"], V3M.PAIR_VEC_LEN), n_boot=n_boot, seed=seed),
                        cafa_fmax_delta_bootstrap=V3M.genome_bootstrap(pc, V3M.paired_stat(lambda v: cafa_fmax(v)["fmax"], CAFA_VEC_LEN), n_boot=n_boot, seed=seed),
                        corrections_and_harms=ch,
                    )
                summaries[mode][pred] = dict(arms=results, paired=paired)
        result = dict(
            schema="gspa.stage2b-followups-price-score/1",
            status="scored",
            policy=M.POLICY,
            predictors=list(predictors),
            arms=eval_arms,
            genomes=sorted(genomes),
            reference_proteins=len(refs),
            exclusions=exclusions,
            ec_audit_sha256=EC_AUDIT_SHA256,
            applied_update_counts_on_reference_proteins=applied_counts,
            bound_inputs=bound,
            summaries=summaries,
            records=rows,
            **design_receipt(),
        )
        (out / "report.json").write_text(json.dumps(result, indent=1, sort_keys=True) + "\n")
        write_summary_tsv(result, out / "summary_exact.tsv")
        return result
    except Exception:
        (out / "failure.log").write_text(traceback.format_exc())
        raise


def write_summary_tsv(result: Mapping[str, Any], path: Path) -> None:
    cols = ["predictor", "arm", "top1", "top1_credit", "top3", "top5", "mrr", "calls_matching", "calls_other", "micro_fmax", "cafa_fmax", "cafa_ci_low", "cafa_ci_high",
            "vs_stage0_changed_top1", "vs_stage0_gained", "vs_stage0_lost", "vs_stage0_pair_corrected", "vs_stage0_pair_harmed", "vs_stage0_cafa_delta", "vs_stage0_cafa_delta_low", "vs_stage0_cafa_delta_high"]
    lines = ["\t".join(cols)]
    for pred, d in result["summaries"]["exact"].items():
        for arm, s in d["arms"].items():
            p = d["paired"].get(f"stage0_to_{arm}", {})
            ch = p.get("corrections_and_harms", {}).get("pair_at_0_5", {})
            cd = p.get("cafa_fmax_delta_bootstrap", {})
            f5 = s.get("fixed_0_5", {})
            lines.append("\t".join(str(x) for x in [
                pred, arm, s["top1_match"], round(s["top1_reference_accuracy"] * s["proteins"], 2), s["top3_match"], s["top5_match"], round(s["lexical_mean_reciprocal_rank"], 4),
                f5.get("matched", ""), f5.get("unconfirmed", ""), round(s["micro_fmax"]["fmax"], 4), round(s["cafa_fmax"]["fmax"], 4),
                round(s["cafa_fmax"]["genome_bootstrap"].get("ci_low", float("nan")), 4), round(s["cafa_fmax"]["genome_bootstrap"].get("ci_high", float("nan")), 4),
                p.get("changed_top1", ""), p.get("top1_gained_reference_match", ""), p.get("top1_lost_reference_match", ""), ch.get("corrected", ""), ch.get("harmed", ""),
                round(cd.get("estimate", float("nan")), 4) if cd else "", round(cd.get("ci_low", float("nan")), 4) if cd else "", round(cd.get("ci_high", float("nan")), 4) if cd else "",
            ]))
    path.write_text("\n".join(lines) + "\n")


def main(argv: Optional[Sequence[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("collection", "manifest", "ec-audit", "updates-dir", "out"):
        p.add_argument("--" + name, type=Path, required=True)
    p.add_argument("--predictors", nargs="+", default=["topec", "mapred", "enzbert", "dpz", "graphec", "clean"])
    p.add_argument("--arms", nargs="+", required=True, help="update arm labels; file <index>.<arm>.json per cell")
    p.add_argument("--n-boot", type=int, default=2000)
    p.add_argument("--seed", type=int, default=13)
    a = p.parse_args(argv)
    score_cells(a.collection, a.manifest, a.ec_audit, a.predictors, a.updates_dir, a.arms, a.out, n_boot=a.n_boot, seed=a.seed)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
