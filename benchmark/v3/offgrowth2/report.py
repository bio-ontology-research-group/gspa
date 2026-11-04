"""Aggregate per-genome results into summary JSON, TSV tables, and Markdown tables (spec §6).

CLI usage::

    python -m offgrowth2.report \\
        --results-dir DIR --out-dir OUT --tier A|B|C --candidate-table six|clean \\
        --labels L0.tsv,L1.tsv --threshold-labels L0.tsv --test-labels L1.tsv \\
        --q-fixed-threshold T [--group-tsv groups.tsv]
"""
from __future__ import annotations

import argparse
import csv
import gc
import json
import statistics
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple, Union

V3_DIR = Path(__file__).resolve().parents[1]
if str(V3_DIR) not in sys.path:
    sys.path.insert(0, str(V3_DIR))

from offgrowth.score import holm_bonferroni_adjust  # noqa: E402

from .panels import load_org_id_map  # noqa: E402
from .score import (  # noqa: E402
    BOOTSTRAP_RESAMPLES,
    BOOTSTRAP_SEED,
    FIXED_PRECISION_MIN_DECISIONS,
    FIXED_PRECISION_TARGET,
    PRIMARY_ARMS,
    ArmRunEvaluation,
    DistinctDecision,
    GenomeScoreResult,
    GenomeScoringContext,
    OrphanStep,
    apply_q_fixed_threshold_to_result,
    build_genome_context,
    clear_genome_caches,
    find_q_fixed_precision_threshold,
    fixed_precision_rule,
    genome_block_bootstrap_precision,
    genome_block_bootstrap_recovery,
    labels_for_genome,
    load_genome_results_json,
    load_genome_score_result,
    load_label_tsv,
    paired_mcnemar_comparison,
    save_genome_score_result,
    score_genome,
)


def discover_genome_result_files(results_dir: Union[Path, str]) -> List[Path]:
    """Find all ``genome_results2.json.gz`` (or ``genome_results2.json``) under ``results_dir``."""
    root = Path(results_dir)
    if root.is_file():
        return [root]
    if not root.exists():
        raise FileNotFoundError(f"Results directory does not exist: {root}")

    gz_files = sorted(
        p for p in root.rglob("genome_results2.json.gz")
        if p.is_file() and "per-genome" not in p.parts
    )
    if gz_files:
        return gz_files
    json_files = sorted(
        p for p in root.rglob("genome_results2.json")
        if p.is_file() and "per-genome" not in p.parts
    )
    if json_files:
        return json_files
    any_gz = sorted(
        p for p in root.rglob("*.json.gz")
        if p.is_file() and "per-genome" not in p.parts
    )
    if any_gz:
        return any_gz
    any_json = sorted(
        p for p in root.rglob("*.json")
        if p.is_file()
        and "per-genome" not in p.parts
        and p.name != "summary.json"
        and not p.name.endswith(".receipt.json")
    )
    return any_json


def load_group_tsv(group_tsv: Optional[Union[Path, str]]) -> Dict[str, str]:
    """Load ``genome_id -> group`` mapping from ``--group-tsv``."""
    if group_tsv is None:
        return {}
    path = Path(group_tsv)
    if not path.exists():
        return {}
    out: Dict[str, str] = {}
    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            gid = (row.get("genome_id") or row.get("org_id") or "").strip()
            grp = (row.get("group") or "").strip()
            if gid and grp:
                out[gid] = grp
    return out


def _fmt_float(val: Optional[float], digits: int = 4) -> str:
    if val is None:
        return ""
    return f"{float(val):.{digits}f}"


def _write_tsv(path: Path, fieldnames: Sequence[str], rows: Sequence[Mapping[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=list(fieldnames),
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _markdown_table(
    caption: str,
    headers: Sequence[str],
    rows: Sequence[Sequence[Any]],
) -> str:
    """Render a one-line caption followed by a GitHub-flavored Markdown table (no prose)."""
    lines = [
        caption.strip(),
        "| " + " | ".join(str(h) for h in headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    for r in rows:
        lines.append("| " + " | ".join(str(cell) for cell in r) + " |")
    return "\n".join(lines)


def _aggregate_group(
    group_name: str,
    genome_results: Sequence[GenomeScoreResult],
    arms_order: Sequence[str],
    n_resamples: int = BOOTSTRAP_RESAMPLES,
    seed: int = BOOTSTRAP_SEED,
) -> Tuple[Dict[str, Any], List[Dict[str, Any]]]:
    """Aggregate recovery, precision, abstention, volume, McNemar, and fixed-precision for a genome group."""
    genome_ids = sorted(gr.genome_id for gr in genome_results)
    n_genomes = len(genome_ids)

    # Collect off-growth masked evaluations per arm (fallback to all masked runs if none has stratum=='off_growth')
    offgrowth_masked_by_arm: Dict[str, List[ArmRunEvaluation]] = {}
    for arm in arms_order:
        if arm in ("Q", "Q_fixed"):
            offgrowth_masked_by_arm[arm] = []
            continue
        all_masked = [
            ev
            for gr in genome_results
            for ev in gr.masked_evals_by_arm.get(arm, ())
        ]
        og_only = [
            ev for ev in all_masked
            if ev.stratum in (None, "", "off_growth")
        ]
        offgrowth_masked_by_arm[arm] = og_only if og_only else all_masked

    # McNemar comparisons vs N1 on gene_recovered and reaction_recovered
    n1_masked = offgrowth_masked_by_arm.get("N1", [])
    mcnemar_gene_raw: Dict[str, Dict[str, Any]] = {}
    mcnemar_rxn_raw: Dict[str, Dict[str, Any]] = {}
    comp_arms = [a for a in arms_order if a not in ("N1", "Q", "Q_fixed")]
    for a in comp_arms:
        mcnemar_gene_raw[a] = paired_mcnemar_comparison(
            n1_masked, offgrowth_masked_by_arm.get(a, []), metric="gene_recovered"
        )
        mcnemar_rxn_raw[a] = paired_mcnemar_comparison(
            n1_masked, offgrowth_masked_by_arm.get(a, []), metric="reaction_recovered"
        )

    gene_holm = holm_bonferroni_adjust([mcnemar_gene_raw[a]["p_value"] for a in comp_arms]) if comp_arms else []
    rxn_holm = holm_bonferroni_adjust([mcnemar_rxn_raw[a]["p_value"] for a in comp_arms]) if comp_arms else []
    for idx, a in enumerate(comp_arms):
        mcnemar_gene_raw[a]["p_holm"] = gene_holm[idx]
        mcnemar_rxn_raw[a]["p_holm"] = rxn_holm[idx]

    arms_summary: Dict[str, Any] = {}
    arms_tsv_rows: List[Dict[str, Any]] = []

    for arm in arms_order:
        masked_evs = offgrowth_masked_by_arm.get(arm, [])
        unmasked_evs = [
            gr.unmasked_eval_by_arm[arm]
            for gr in genome_results
            if arm in gr.unmasked_eval_by_arm
        ]

        gene_rec_boot = genome_block_bootstrap_recovery(
            masked_evs, metric="gene_recovered", n_resamples=n_resamples, seed=seed
        )
        rxn_rec_boot = genome_block_bootstrap_recovery(
            masked_evs, metric="reaction_recovered", n_resamples=n_resamples, seed=seed
        )
        pw_rec_boot = genome_block_bootstrap_recovery(
            masked_evs, metric="pathway_recovered", n_resamples=n_resamples, seed=seed
        )

        # Decisions by label set: test (L1), threshold (L0), pooled (L0+L1)
        dec_by_set = {
            "test": [d for gr in genome_results for d in gr.test_decisions_by_arm.get(arm, ())],
            "threshold": [d for gr in genome_results for d in gr.threshold_decisions_by_arm.get(arm, ())],
            "pooled": [d for gr in genome_results for d in gr.pooled_decisions_by_arm.get(arm, ())],
        }
        prec_by_set: Dict[str, Dict[str, Any]] = {}
        fp_by_set: Dict[str, Dict[str, Any]] = {}
        for lset_name, d_list in dec_by_set.items():
            prec_by_set[lset_name] = genome_block_bootstrap_precision(
                d_list, all_genome_ids=genome_ids, n_resamples=n_resamples, seed=seed
            )
            fp_by_set[lset_name] = fixed_precision_rule(
                d_list,
                min_precision=FIXED_PRECISION_TARGET,
                min_decisions=FIXED_PRECISION_MIN_DECISIONS,
            )

        # Abstention and volume across unmasked runs
        tot_cand_rxns = sum(ev.n_candidate_reactions for ev in unmasked_evs)
        tot_called_rxns = sum(ev.n_called_reactions for ev in unmasked_evs)
        tot_rxn_only = sum(ev.n_reaction_only_calls for ev in unmasked_evs)
        tot_gene_calls = sum(ev.n_distinct_gene_calls for ev in unmasked_evs)
        tot_pw_calls = sum(ev.n_pathway_calls for ev in unmasked_evs)
        tot_lps = sum(ev.n_lps for ev in unmasked_evs)
        tot_cpu = sum(ev.cpu_seconds for ev in unmasked_evs)

        gene_abst = float(tot_rxn_only) / float(tot_called_rxns) if tot_called_rxns > 0 else 0.0
        rxn_abst = (
            float(tot_cand_rxns - tot_called_rxns) / float(tot_cand_rxns)
            if tot_cand_rxns > 0
            else 0.0
        )

        n_unm = len(unmasked_evs)
        mean_gene_vol = float(tot_gene_calls) / float(n_unm) if n_unm > 0 else 0.0
        mean_rxn_only_vol = float(tot_rxn_only) / float(n_unm) if n_unm > 0 else 0.0
        mean_pw_vol = float(tot_pw_calls) / float(n_unm) if n_unm > 0 else 0.0
        mean_lps = float(tot_lps) / float(n_unm) if n_unm > 0 else 0.0
        mean_cpu = float(tot_cpu) / float(n_unm) if n_unm > 0 else 0.0

        mc_g = mcnemar_gene_raw.get(arm)
        mc_r = mcnemar_rxn_raw.get(arm)

        # Record rows in arms.tsv
        if arm not in ("Q", "Q_fixed"):
            arms_tsv_rows.append({
                "group": group_name,
                "arm": arm,
                "metric": "gene_recovery",
                "label_set": "masked_offgrowth",
                "point": _fmt_float(gene_rec_boot["point"]),
                "ci_low": _fmt_float(gene_rec_boot["ci_low"]),
                "ci_high": _fmt_float(gene_rec_boot["ci_high"]),
                "numerator": gene_rec_boot["numerator"],
                "denominator": gene_rec_boot["denominator"],
                "n_genomes": gene_rec_boot["n_genomes"],
                "mcnemar_p": _fmt_float(mc_g["p_value"], 6) if mc_g else "",
                "mcnemar_p_holm": _fmt_float(mc_g["p_holm"], 6) if mc_g else "",
            })
            arms_tsv_rows.append({
                "group": group_name,
                "arm": arm,
                "metric": "reaction_recovery",
                "label_set": "masked_offgrowth",
                "point": _fmt_float(rxn_rec_boot["point"]),
                "ci_low": _fmt_float(rxn_rec_boot["ci_low"]),
                "ci_high": _fmt_float(rxn_rec_boot["ci_high"]),
                "numerator": rxn_rec_boot["numerator"],
                "denominator": rxn_rec_boot["denominator"],
                "n_genomes": rxn_rec_boot["n_genomes"],
                "mcnemar_p": _fmt_float(mc_r["p_value"], 6) if mc_r else "",
                "mcnemar_p_holm": _fmt_float(mc_r["p_holm"], 6) if mc_r else "",
            })
            arms_tsv_rows.append({
                "group": group_name,
                "arm": arm,
                "metric": "pathway_recovery",
                "label_set": "masked_offgrowth",
                "point": _fmt_float(pw_rec_boot["point"]),
                "ci_low": _fmt_float(pw_rec_boot["ci_low"]),
                "ci_high": _fmt_float(pw_rec_boot["ci_high"]),
                "numerator": pw_rec_boot["numerator"],
                "denominator": pw_rec_boot["denominator"],
                "n_genomes": pw_rec_boot["n_genomes"],
                "mcnemar_p": "",
                "mcnemar_p_holm": "",
            })

        for lset_name in ("test", "threshold", "pooled"):
            pb = prec_by_set[lset_name]
            fp = fp_by_set[lset_name]
            arms_tsv_rows.append({
                "group": group_name,
                "arm": arm,
                "metric": "precision",
                "label_set": lset_name,
                "point": _fmt_float(pb["precision"]),
                "ci_low": _fmt_float(pb["precision_ci_low"]),
                "ci_high": _fmt_float(pb["precision_ci_high"]),
                "numerator": pb["n_supported"],
                "denominator": pb["n_labelled"],
                "n_genomes": pb["n_genomes"],
                "mcnemar_p": "",
                "mcnemar_p_holm": "",
            })
            arms_tsv_rows.append({
                "group": group_name,
                "arm": arm,
                "metric": "harm",
                "label_set": lset_name,
                "point": _fmt_float(pb["harm"]),
                "ci_low": _fmt_float(pb["harm_ci_low"]),
                "ci_high": _fmt_float(pb["harm_ci_high"]),
                "numerator": pb["n_harmful"],
                "denominator": pb["n_labelled"],
                "n_genomes": pb["n_genomes"],
                "mcnemar_p": "",
                "mcnemar_p_holm": "",
            })
            rel_frac = (
                float(fp["n_related_substrate"]) / float(fp["n_harmful"])
                if fp["n_harmful"] > 0
                else 0.0
            )
            arms_tsv_rows.append({
                "group": group_name,
                "arm": arm,
                "metric": "related_substrate_harm_share",
                "label_set": lset_name,
                "point": _fmt_float(rel_frac),
                "ci_low": "",
                "ci_high": "",
                "numerator": fp["n_related_substrate"],
                "denominator": fp["n_harmful"],
                "n_genomes": pb["n_genomes"],
                "mcnemar_p": "",
                "mcnemar_p_holm": "",
            })

        arms_tsv_rows.append({
            "group": group_name,
            "arm": arm,
            "metric": "gene_abstention_rate",
            "label_set": "unmasked",
            "point": _fmt_float(gene_abst),
            "ci_low": "",
            "ci_high": "",
            "numerator": tot_rxn_only,
            "denominator": tot_called_rxns,
            "n_genomes": n_unm,
            "mcnemar_p": "",
            "mcnemar_p_holm": "",
        })
        arms_tsv_rows.append({
            "group": group_name,
            "arm": arm,
            "metric": "reaction_abstention_rate",
            "label_set": "unmasked",
            "point": _fmt_float(rxn_abst),
            "ci_low": "",
            "ci_high": "",
            "numerator": tot_cand_rxns - tot_called_rxns,
            "denominator": tot_cand_rxns,
            "n_genomes": n_unm,
            "mcnemar_p": "",
            "mcnemar_p_holm": "",
        })

        arms_summary[arm] = {
            "gene_recovery": gene_rec_boot,
            "reaction_recovery": rxn_rec_boot,
            "pathway_recovery": pw_rec_boot,
            "precision_by_label_set": prec_by_set,
            "fixed_precision_by_label_set": fp_by_set,
            "gene_abstention_rate": gene_abst,
            "reaction_abstention_rate": rxn_abst,
            "volume": {
                "total_called_reactions": tot_called_rxns,
                "total_distinct_gene_calls": tot_gene_calls,
                "total_reaction_only_calls": tot_rxn_only,
                "total_pathway_calls": tot_pw_calls,
                "mean_gene_calls_per_genome": mean_gene_vol,
                "mean_reaction_only_calls_per_genome": mean_rxn_only_vol,
                "mean_pathway_calls_per_genome": mean_pw_vol,
            },
            "cost": {
                "total_lps": tot_lps,
                "total_cpu_seconds": tot_cpu,
                "mean_lps_per_genome": mean_lps,
                "mean_cpu_seconds_per_genome": mean_cpu,
            },
            "mcnemar_vs_N1_gene": mc_g,
            "mcnemar_vs_N1_reaction": mc_r,
        }

    group_summary = {
        "group": group_name,
        "n_genomes": n_genomes,
        "genome_ids": genome_ids,
        "arms": arms_summary,
    }
    return group_summary, arms_tsv_rows


_PHASE1_CONFIG: Tuple[
    Path,
    Dict[str, str],
    Optional[Dict[str, Dict[str, List[str]]]],
    Optional[Dict[str, Dict[str, List[str]]]],
    Optional[Union[Path, str]],
    List[str],
    Optional[float],
    bool,
] = (Path("."), {}, None, None, None, [], None, False)


def _score_genome_phase1_worker(rf_str: str) -> str:
    """Phase-1 worker: load one genome payload, build context, score all arms, write `<out-dir>/per-genome/<genome>.json.gz`, free memory."""
    (
        per_genome_dir,
        org_id_map,
        thresh_label_store,
        test_label_store,
        matches_tsv,
        arms,
        q_thresh,
        resume,
    ) = _PHASE1_CONFIG
    rf = Path(rf_str)

    if resume and rf.name.startswith("genome_results2.json") and rf.parent.name:
        fast_out = per_genome_dir / f"{rf.parent.name}.json.gz"
        if fast_out.is_file() and fast_out.stat().st_size > 0:
            return str(fast_out)

    payload = load_genome_results_json(rf)
    gid = str(payload.get("genome_id") or rf.stem.split(".")[0])
    out_file = per_genome_dir / f"{gid}.json.gz"
    if resume and out_file.is_file() and out_file.stat().st_size > 0:
        del payload
        return str(out_file)

    org_id = org_id_map.get(gid, str(payload.get("org_id") or gid))

    inputs_dict = payload.get("inputs") or {}
    cand_path: Optional[Path] = None
    trans_path: Optional[Path] = None
    if "candidates" not in payload:
        for p_str in inputs_dict:
            p_obj = Path(p_str)
            if p_obj.exists():
                if "candidate" in p_obj.name or p_obj.name.endswith(".clean.tsv"):
                    cand_path = p_obj
                elif p_obj.name.endswith("-Transporter.tbl"):
                    trans_path = p_obj

    g_known = payload.get("known_labels") or {}
    g_thresh = labels_for_genome(gid, thresh_label_store, g_known)
    g_test = labels_for_genome(gid, test_label_store, g_known)

    ctx = build_genome_context(
        payload,
        candidates=cand_path,
        transporters=trans_path,
        matches_tsv=matches_tsv,
        org_id=org_id,
        threshold_labels=g_thresh,
        test_labels=g_test,
    )
    gr = score_genome(
        payload,
        ctx=ctx,
        arms=arms,
        q_fixed_threshold=q_thresh,
    )
    save_genome_score_result(gr, out_file)

    clear_genome_caches(gid)
    del payload, ctx, gr
    gc.collect()
    return str(out_file)


def generate_report(
    results_dir: Union[Path, str],
    out_dir: Union[Path, str],
    tier: str = "A",
    candidate_table: str = "six",
    labels: Optional[Union[Path, str, Sequence[Union[Path, str]]]] = None,
    threshold_labels: Optional[Union[Path, str, Sequence[Union[Path, str]]]] = None,
    test_labels: Optional[Union[Path, str, Sequence[Union[Path, str]]]] = None,
    q_fixed_threshold: Optional[float] = None,
    group_tsv: Optional[Union[Path, str]] = None,
    matches_tsv: Optional[Union[Path, str]] = None,
    org_map_tsv: Optional[Union[Path, str]] = None,
    arms: Sequence[str] = PRIMARY_ARMS,
    n_resamples: int = BOOTSTRAP_RESAMPLES,
    seed: int = BOOTSTRAP_SEED,
    resume: bool = False,
) -> Dict[str, Any]:
    """Read all ``genome_results2.json.gz`` under ``results_dir`` and write all report artifacts to ``out_dir``."""
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    per_genome_dir = out_path / "per-genome"
    per_genome_dir.mkdir(parents=True, exist_ok=True)

    result_files = discover_genome_result_files(results_dir)
    if not result_files:
        raise FileNotFoundError(f"No genome result files found under {results_dir}")

    all_label_store = load_label_tsv(labels) if labels is not None else None
    thresh_label_store = load_label_tsv(threshold_labels) if threshold_labels is not None else all_label_store
    test_label_store = load_label_tsv(test_labels) if test_labels is not None else all_label_store

    group_map = load_group_tsv(group_tsv)
    org_id_map = load_org_id_map(org_map_tsv)

    # Phase 1: score each genome in a worker from its file path, write <out-dir>/per-genome/<genome>.json.gz, free everything.
    import os as _os
    import multiprocessing as _mp

    global _PHASE1_CONFIG
    _PHASE1_CONFIG = (
        per_genome_dir,
        org_id_map,
        thresh_label_store,
        test_label_store,
        matches_tsv,
        list(arms),
        q_fixed_threshold,
        bool(resume),
    )

    _jobs = max(1, int(_os.environ.get("OFFGROWTH2_SCORE_JOBS", "1") or "1"))
    rf_strs = [str(rf) for rf in result_files]
    if _jobs > 1 and len(rf_strs) > 1 and "fork" in _mp.get_all_start_methods():
        with _mp.get_context("fork").Pool(
            processes=min(_jobs, len(rf_strs)),
            maxtasksperchild=1,
        ) as pool:
            per_genome_paths = list(pool.map(_score_genome_phase1_worker, rf_strs))
    else:
        per_genome_paths = [_score_genome_phase1_worker(rf_s) for rf_s in rf_strs]

    # Phase 2 in the parent: read compact per-genome results and write summary/tables.
    genome_results: List[GenomeScoreResult] = [
        load_genome_score_result(p_str) for p_str in per_genome_paths
    ]
    genome_results.sort(key=lambda g: g.genome_id)

    # Determine Q fixed-precision threshold on L0 from phase-1 results if not explicitly passed
    resolved_q_thresh = q_fixed_threshold
    if resolved_q_thresh is None:
        resolved_q_thresh = find_q_fixed_precision_threshold(
            genome_results,
            min_precision=FIXED_PRECISION_TARGET,
            min_decisions=FIXED_PRECISION_MIN_DECISIONS,
        )

    for gr in genome_results:
        apply_q_fixed_threshold_to_result(gr, resolved_q_thresh)
        gr.q_best_per_protein = []

    arms_with_q: List[str] = list(arms) + ["Q"]
    if resolved_q_thresh is not None:
        arms_with_q.append("Q_fixed")

    # Group aggregations: always "pooled", plus each distinct group in group_map
    pooled_summary, all_arms_rows = _aggregate_group(
        "pooled",
        genome_results,
        arms_order=arms_with_q,
        n_resamples=n_resamples,
        seed=seed,
    )
    groups_summary: Dict[str, Any] = {"pooled": pooled_summary}

    if group_map:
        by_grp: Dict[str, List[GenomeScoreResult]] = {}
        for gr in genome_results:
            grp = group_map.get(gr.genome_id, "other")
            by_grp.setdefault(grp, []).append(gr)
        for grp in sorted(by_grp):
            g_sum, g_rows = _aggregate_group(
                grp,
                by_grp[grp],
                arms_order=arms_with_q,
                n_resamples=n_resamples,
                seed=seed,
            )
            groups_summary[grp] = g_sum
            all_arms_rows.extend(g_rows)

    # 1. Write arms.tsv
    arms_cols = [
        "group",
        "arm",
        "metric",
        "label_set",
        "point",
        "ci_low",
        "ci_high",
        "numerator",
        "denominator",
        "n_genomes",
        "mcnemar_p",
        "mcnemar_p_holm",
    ]
    _write_tsv(out_path / "arms.tsv", arms_cols, all_arms_rows)

    # 2. Write cases.tsv
    cases_rows: List[Dict[str, Any]] = []
    for gr in genome_results:
        gid = gr.genome_id
        grp = group_map.get(gid, "pooled")
        # Metadata per case from genome payload
        for arm in arms:
            for ev in sorted(gr.masked_evals_by_arm.get(arm, ()), key=lambda e: str(e.case_id or "")):
                # Look up case metadata from the run if available
                cases_rows.append({
                    "genome_id": gid,
                    "group": grp,
                    "case_id": ev.case_id or "",
                    "stratum": ev.stratum or "off_growth",
                    "arm": arm,
                    "reaction_recovered": int(ev.reaction_recovered),
                    "gene_recovered": int(ev.gene_recovered),
                    "pathway_recovered": int(ev.pathway_recovered),
                    "recovering_fill_size": "" if ev.recovering_fill_size is None else ev.recovering_fill_size,
                    "hidden_protein_rank": "" if ev.hidden_protein_rank is None else ev.hidden_protein_rank,
                    "n_called_reactions": ev.n_called_reactions,
                    "n_gene_calls": ev.n_gene_called_reactions,
                    "n_pathway_calls": ev.n_pathway_calls,
                    "n_lps": ev.n_lps,
                    "cpu_seconds": _fmt_float(ev.cpu_seconds, 4),
                })
    cases_rows.sort(key=lambda r: (str(r["genome_id"]), str(r["case_id"]), str(r["arm"])))
    cases_cols = [
        "genome_id",
        "group",
        "case_id",
        "stratum",
        "arm",
        "reaction_recovered",
        "gene_recovered",
        "pathway_recovered",
        "recovering_fill_size",
        "hidden_protein_rank",
        "n_called_reactions",
        "n_gene_calls",
        "n_pathway_calls",
        "n_lps",
        "cpu_seconds",
    ]
    _write_tsv(out_path / "cases.tsv", cases_cols, cases_rows)

    # 3. Write decisions.tsv (distinct decisions of unmasked runs)
    decisions_rows: List[Dict[str, Any]] = []
    for gr in genome_results:
        gid = gr.genome_id
        grp = group_map.get(gid, "pooled")
        for lset_name, dec_dict in (
            ("test", gr.test_decisions_by_arm),
            ("threshold", gr.threshold_decisions_by_arm),
            ("pooled", gr.pooled_decisions_by_arm),
        ):
            for arm in arms_with_q:
                for d in dec_dict.get(arm, ()):
                    decisions_rows.append({
                        "genome_id": gid,
                        "group": grp,
                        "arm": arm,
                        "label_set": lset_name,
                        "protein_id": d.protein_id,
                        "ec": d.ec,
                        "score": _fmt_float(d.score, 6),
                        "reactions": ",".join(d.reactions),
                        "classification": d.classification,
                        "related_substrate": int(d.related_substrate),
                        "truth_ecs": ",".join(d.truth_ecs),
                    })
    decisions_rows.sort(
        key=lambda r: (
            str(r["label_set"]),
            str(r["genome_id"]),
            str(r["arm"]),
            str(r["protein_id"]),
            str(r["ec"]),
        )
    )
    decisions_cols = [
        "genome_id",
        "group",
        "arm",
        "label_set",
        "protein_id",
        "ec",
        "score",
        "reactions",
        "classification",
        "related_substrate",
        "truth_ecs",
    ]
    _write_tsv(out_path / "decisions.tsv", decisions_cols, decisions_rows)

    # 4. Write pathway_calls.tsv
    pw_rows: List[Dict[str, Any]] = []
    for gr in genome_results:
        gid = gr.genome_id
        grp = group_map.get(gid, "pooled")
        for arm in arms:
            unm = gr.unmasked_eval_by_arm.get(arm)
            if unm is not None:
                for pc in unm.pathway_calls:
                    pw_rows.append({
                        "genome_id": gid,
                        "group": grp,
                        "run_kind": "unmasked",
                        "case_id": "",
                        "arm": arm,
                        "test_id": pc.test_id,
                        "contig": pc.contig,
                        "cluster_span": pc.cluster_span,
                        "reactions": ",".join(pc.reactions),
                        "genes": ",".join(pc.genes),
                    })
            for mev in gr.masked_evals_by_arm.get(arm, ()):
                for pc in mev.pathway_calls:
                    pw_rows.append({
                        "genome_id": gid,
                        "group": grp,
                        "run_kind": "masked",
                        "case_id": mev.case_id or "",
                        "arm": arm,
                        "test_id": pc.test_id,
                        "contig": pc.contig,
                        "cluster_span": pc.cluster_span,
                        "reactions": ",".join(pc.reactions),
                        "genes": ",".join(pc.genes),
                    })
    pw_rows.sort(
        key=lambda r: (
            str(r["genome_id"]),
            str(r["run_kind"]),
            str(r["case_id"]),
            str(r["arm"]),
            str(r["test_id"]),
            str(r["reactions"]),
        )
    )
    pw_cols = [
        "genome_id",
        "group",
        "run_kind",
        "case_id",
        "arm",
        "test_id",
        "contig",
        "cluster_span",
        "reactions",
        "genes",
    ]
    _write_tsv(out_path / "pathway_calls.tsv", pw_cols, pw_rows)

    # 5. Write orphan_steps.tsv
    orphan_rows: List[Dict[str, Any]] = []
    for gr in genome_results:
        gid = gr.genome_id
        grp = group_map.get(gid, "pooled")
        for os_entry in gr.orphan_steps:
            orphan_rows.append({
                "genome_id": gid,
                "group": grp,
                "rxn_id": os_entry.rxn_id,
                "ecs": ",".join(os_entry.ecs),
                "max_genome_score": _fmt_float(os_entry.max_genome_score, 6),
                "smallest_fill": os_entry.smallest_fill,
                "n_tests": os_entry.n_tests,
                "test_types": ",".join(os_entry.test_types),
                "calling_tests": ",".join(os_entry.calling_tests),
            })
    orphan_rows.sort(key=lambda r: (str(r["genome_id"]), str(r["rxn_id"])))
    orphan_cols = [
        "genome_id",
        "group",
        "rxn_id",
        "ecs",
        "max_genome_score",
        "smallest_fill",
        "n_tests",
        "test_types",
        "calling_tests",
    ]
    _write_tsv(out_path / "orphan_steps.tsv", orphan_cols, orphan_rows)

    # 6. Write cost.tsv
    cost_rows: List[Dict[str, Any]] = []
    all_t1_fills_n1: List[int] = []
    all_t1_stable_n1: List[int] = []
    for gr in genome_results:
        gid = gr.genome_id
        grp = group_map.get(gid, "pooled")
        for arm in arms:
            unm = gr.unmasked_eval_by_arm.get(arm)
            if unm is not None:
                if arm == "N1":
                    all_t1_fills_n1.extend(unm.t1_fill_sizes)
                    all_t1_stable_n1.extend(unm.t1_stable_fill_sizes)
                med_fill = statistics.median(unm.t1_fill_sizes) if unm.t1_fill_sizes else 0.0
                med_stab = statistics.median(unm.t1_stable_fill_sizes) if unm.t1_stable_fill_sizes else 0.0
                cost_rows.append({
                    "genome_id": gid,
                    "group": grp,
                    "run_kind": "unmasked",
                    "case_id": "",
                    "arm": arm,
                    "n_admitted_tests": unm.n_admitted_tests,
                    "n_lps": unm.n_lps,
                    "cpu_seconds": _fmt_float(unm.cpu_seconds, 4),
                    "stability_rate": _fmt_float(unm.stability_rate, 4),
                    "t1_median_fill_size": _fmt_float(med_fill, 2),
                    "t1_median_stable_fill_size": _fmt_float(med_stab, 2),
                    "measured_panel_fallback": int(unm.measured_panel_fallback),
                })
            for mev in gr.masked_evals_by_arm.get(arm, ()):
                med_fill = statistics.median(mev.t1_fill_sizes) if mev.t1_fill_sizes else 0.0
                med_stab = statistics.median(mev.t1_stable_fill_sizes) if mev.t1_stable_fill_sizes else 0.0
                cost_rows.append({
                    "genome_id": gid,
                    "group": grp,
                    "run_kind": "masked",
                    "case_id": mev.case_id or "",
                    "arm": arm,
                    "n_admitted_tests": mev.n_admitted_tests,
                    "n_lps": mev.n_lps,
                    "cpu_seconds": _fmt_float(mev.cpu_seconds, 4),
                    "stability_rate": _fmt_float(mev.stability_rate, 4),
                    "t1_median_fill_size": _fmt_float(med_fill, 2),
                    "t1_median_stable_fill_size": _fmt_float(med_stab, 2),
                    "measured_panel_fallback": int(mev.measured_panel_fallback),
                })
    cost_rows.sort(
        key=lambda r: (
            str(r["genome_id"]),
            str(r["run_kind"]),
            str(r["case_id"]),
            str(r["arm"]),
        )
    )
    cost_cols = [
        "genome_id",
        "group",
        "run_kind",
        "case_id",
        "arm",
        "n_admitted_tests",
        "n_lps",
        "cpu_seconds",
        "stability_rate",
        "t1_median_fill_size",
        "t1_median_stable_fill_size",
        "measured_panel_fallback",
    ]
    _write_tsv(out_path / "cost.tsv", cost_cols, cost_rows)

    # 7. Build report.md (tables with one-line captions only, no prose)
    md_blocks: List[str] = []
    table_idx = 1

    for grp_key in ["pooled"] + sorted(k for k in groups_summary if k != "pooled"):
        g_sum = groups_summary[grp_key]
        t_rows: List[List[str]] = []
        for arm in arms_with_q:
            a_info = g_sum["arms"][arm]
            gr_b = a_info["gene_recovery"]
            rr_b = a_info["reaction_recovery"]
            pr_b = a_info["pathway_recovery"]
            p_test = a_info["precision_by_label_set"]["test"]
            vol = a_info["volume"]
            cst = a_info["cost"]
            t_rows.append([
                arm,
                f"{gr_b['numerator']}/{gr_b['denominator']} ({_fmt_float(gr_b['point'])})",
                f"[{_fmt_float(gr_b['ci_low'])}, {_fmt_float(gr_b['ci_high'])}]",
                f"{rr_b['numerator']}/{rr_b['denominator']} ({_fmt_float(rr_b['point'])})",
                f"{pr_b['numerator']}/{pr_b['denominator']} ({_fmt_float(pr_b['point'])})",
                f"{p_test['n_supported']}/{p_test['n_labelled']} ({_fmt_float(p_test['precision'])})",
                f"[{_fmt_float(p_test['precision_ci_low'])}, {_fmt_float(p_test['precision_ci_high'])}]",
                _fmt_float(p_test["harm"]),
                _fmt_float(a_info["gene_abstention_rate"]),
                _fmt_float(a_info["reaction_abstention_rate"]),
                _fmt_float(vol["mean_gene_calls_per_genome"], 2),
                _fmt_float(vol["mean_reaction_only_calls_per_genome"], 2),
                _fmt_float(vol["mean_pathway_calls_per_genome"], 2),
                _fmt_float(cst["mean_lps_per_genome"], 1),
            ])
        md_blocks.append(
            _markdown_table(
                caption=f"Table {table_idx}: Recovery, precision, abstention, volume, and LP cost by arm ({grp_key}, n={g_sum['n_genomes']} genomes).",
                headers=[
                    "Arm",
                    "Gene recovery",
                    "Gene rec 95% CI",
                    "Rxn recovery",
                    "Pathway recovery",
                    "Test precision",
                    "Precision 95% CI",
                    "Test harm",
                    "Gene abstention",
                    "Rxn abstention",
                    "Gene calls/genome",
                    "Rxn-only/genome",
                    "Pathway calls/genome",
                    "LPs/genome",
                ],
                rows=t_rows,
            )
        )
        table_idx += 1

    # McNemar table (pooled)
    mc_rows: List[List[str]] = []
    for arm in arms:
        if arm == "N1":
            continue
        a_info = pooled_summary["arms"][arm]
        mg = a_info["mcnemar_vs_N1_gene"] or {}
        mr = a_info["mcnemar_vs_N1_reaction"] or {}
        mc_rows.append([
            arm,
            str(mg.get("n_cases", 0)),
            str(mg.get("b_better_only", 0)),
            str(mg.get("c_baseline_only", 0)),
            _fmt_float(mg.get("p_value"), 6),
            _fmt_float(mg.get("p_holm"), 6),
            str(mr.get("b_better_only", 0)),
            str(mr.get("c_baseline_only", 0)),
            _fmt_float(mr.get("p_value"), 6),
            _fmt_float(mr.get("p_holm"), 6),
        ])
    md_blocks.append(
        _markdown_table(
            caption=f"Table {table_idx}: Paired exact one-sided McNemar comparisons of N1 against comparator arms on off-growth masked cases.",
            headers=[
                "Comparator arm",
                "Cases",
                "Gene N1-only",
                "Gene arm-only",
                "Gene p",
                "Gene p (Holm)",
                "Rxn N1-only",
                "Rxn arm-only",
                "Rxn p",
                "Rxn p (Holm)",
            ],
            rows=mc_rows,
        )
    )
    table_idx += 1

    # Fixed-precision rule table across label sets
    fp_rows: List[List[str]] = []
    for arm in arms_with_q:
        a_info = pooled_summary["arms"][arm]
        for lset_name in ("threshold", "test", "pooled"):
            fp = a_info["fixed_precision_by_label_set"][lset_name]
            fp_rows.append([
                arm,
                lset_name,
                str(fp["n_supported"]),
                str(fp["n_harmful"]),
                str(fp["n_related_substrate"]),
                str(fp["n_unconfirmed"]),
                str(fp["n_labelled"]),
                _fmt_float(fp["precision"]),
                "PASS" if fp["passes"] else "FAIL",
            ])
    md_blocks.append(
        _markdown_table(
            caption=f"Table {table_idx}: Fixed-precision rule evaluation (precision >= {FIXED_PRECISION_TARGET:.2f} on >= {FIXED_PRECISION_MIN_DECISIONS} labelled decisions) by arm and label set.",
            headers=[
                "Arm",
                "Label set",
                "Supported",
                "Harmful",
                "Related-substrate",
                "Unconfirmed",
                "Labelled",
                "Precision",
                "Rule status",
            ],
            rows=fp_rows,
        )
    )
    table_idx += 1

    # Orphan steps table
    orphan_md_rows = [
        [
            r["genome_id"],
            r["rxn_id"],
            r["ecs"] or "-",
            str(r["max_genome_score"]),
            str(r["smallest_fill"]),
            str(r["n_tests"]),
            r["test_types"] or "-",
        ]
        for r in orphan_rows
    ]
    if not orphan_md_rows:
        orphan_md_rows = [["-", "-", "-", "0.000000", "0", "0", "-"]]
    md_blocks.append(
        _markdown_table(
            caption=f"Table {table_idx}: Orphan reaction steps called by N1 with no genome candidate >= { FIXED_PRECISION_TARGET * 0 + 0.30:.2f}.",
            headers=[
                "Genome",
                "Reaction",
                "ECs",
                "Max genome score",
                "Smallest fill",
                "Calling tests",
                "Test types",
            ],
            rows=orphan_md_rows,
        )
    )

    (out_path / "report.md").write_text("\n\n".join(md_blocks) + "\n", encoding="utf-8")

    # 8. Write summary.json
    additional_funcs_all = [
        {
            "genome_id": d.genome_id,
            "protein_id": d.protein_id,
            "ec": d.ec,
            "score": d.score,
            "classification": d.classification,
            "reactions": list(d.reactions),
        }
        for gr in genome_results
        for d in gr.additional_functions
    ]
    t1_fill_summary = {
        "n_tests": len(all_t1_fills_n1),
        "median_fill": float(statistics.median(all_t1_fills_n1)) if all_t1_fills_n1 else 0.0,
        "max_fill": max(all_t1_fills_n1) if all_t1_fills_n1 else 0,
        "median_stable_fill": float(statistics.median(all_t1_stable_n1)) if all_t1_stable_n1 else 0.0,
        "max_stable_fill": max(all_t1_stable_n1) if all_t1_stable_n1 else 0,
    }

    summary_payload: Dict[str, Any] = {
        "schema": "gspa.offgrowth2.summary/1",
        "tier": str(tier),
        "candidate_table": str(candidate_table),
        "q_fixed_threshold": resolved_q_thresh,
        "n_genomes": len(genome_results),
        "genome_ids": [gr.genome_id for gr in genome_results],
        "groups": groups_summary,
        "arms": pooled_summary["arms"],
        "t1_fill_distribution": t1_fill_summary,
        "n_orphan_steps": len(orphan_rows),
        "n_additional_functions": len(additional_funcs_all),
        "additional_functions": additional_funcs_all,
    }
    (out_path / "summary.json").write_text(
        json.dumps(summary_payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return summary_payload


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="Aggregate offgrowth2 genome results into summary JSON, TSVs, and Markdown tables.")
    ap.add_argument("--results-dir", type=Path, required=True, help="Directory with genome_results2.json.gz files")
    ap.add_argument("--out-dir", type=Path, required=True, help="Output directory for summary.json, *.tsv, and report.md")
    ap.add_argument("--tier", type=str, default="A", choices=["A", "B", "C"], help="Evaluation tier (A, B, or C)")
    ap.add_argument("--candidate-table", type=str, default="six", choices=["six", "clean"], help="Candidate table kind (six or clean)")
    ap.add_argument("--labels", type=str, default=None, help="Comma-separated label TSV path(s)")
    ap.add_argument("--threshold-labels", "--labels-for-threshold", dest="threshold_labels", type=str, default=None, help="Threshold label TSV (L0)")
    ap.add_argument("--test-labels", "--labels-for-test", dest="test_labels", type=str, default=None, help="Test label TSV (L1)")
    ap.add_argument("--q-fixed-threshold", type=float, default=None, help="Fixed threshold T for comparator Q")
    ap.add_argument("--group-tsv", type=Path, default=None, help="Optional genome_id, group TSV (e.g. Pseudomonas split)")
    ap.add_argument("--matches-tsv", type=Path, default=None, help="Optional measured experiment_panel_matches.tsv")
    ap.add_argument("--org-map-tsv", type=Path, default=None, help="Optional price22_fitness_browser_organisms.tsv")
    ap.add_argument("--bootstrap-resamples", type=int, default=BOOTSTRAP_RESAMPLES, help="Bootstrap resamples (default 2000)")
    ap.add_argument("--resume", action="store_true", help="Skip genomes whose per-genome result file already exists")
    args = ap.parse_args(argv)

    generate_report(
        results_dir=args.results_dir,
        out_dir=args.out_dir,
        tier=args.tier,
        candidate_table=args.candidate_table,
        labels=args.labels,
        threshold_labels=args.threshold_labels,
        test_labels=args.test_labels,
        q_fixed_threshold=args.q_fixed_threshold,
        group_tsv=args.group_tsv,
        matches_tsv=args.matches_tsv,
        org_map_tsv=args.org_map_tsv,
        n_resamples=args.bootstrap_resamples,
        resume=args.resume,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
