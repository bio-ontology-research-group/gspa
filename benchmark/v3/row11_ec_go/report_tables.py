#!/usr/bin/env python3
"""Render Markdown summary tables from row-11 Universe A and Universe B JSON artifacts."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "benchmark" / "v3"))
import evaluate_stage3_agent as E  # noqa: E402

METHOD_ORDER = (
    "abstain_only",
    "rank_by_raw_seq",
    "rank_by_stage1_posterior",
    "rank_by_context_support",
    "learned_logistic",
    "learned_gbdt",
    "qwen_agent",
)


def _fmt_float(val, digits: int = 3) -> str:
    if val is None:
        return "n/a"
    try:
        return f"{float(val):.{digits}f}"
    except (TypeError, ValueError):
        return str(val)


def _fmt_ci(rate, ci, digits: int = 3) -> str:
    if rate is None:
        return "n/a"
    if not ci or len(ci) != 2 or ci[0] is None or ci[1] is None:
        return _fmt_float(rate, digits)
    return f"{float(rate):.{digits}f} ({float(ci[0]):.{digits}f}--{float(ci[1]):.{digits}f})"


def _infer_set_name(path: Path, root: Path) -> str:
    rel_parts = path.relative_to(root).parts
    for token in rel_parts:
        low = token.lower()
        if "ec6c" in low or "ec6+contract" in low or "contract" in low:
            return "ec6c"
        if "ec6" in low:
            return "ec6"
    if len(rel_parts) > 1:
        return rel_parts[0]
    return "default"


def _discover_artifacts(results_dir: Path) -> dict[str, dict]:
    sets: dict[str, dict] = {}

    def bucket(name: str) -> dict:
        return sets.setdefault(
            name,
            {
                "public_a": None,
                "public_b": None,
                "private_b": None,
                "truth_b": None,
                "truth_report": None,
                "proposals_b": {},
                "policies_b": {},
                "evals_b": {},
                "leak_b": None,
                "compare_b": None,
                "policies_a": {},
                "evals_a": {},
            },
        )

    for p in sorted(Path(results_dir).rglob("*.json")):
        try:
            data = json.loads(p.read_text(encoding="utf-8"))
        except Exception:
            continue
        if not isinstance(data, dict):
            continue
        schema = data.get("schema", "")
        sname = _infer_set_name(p, Path(results_dir))
        b = bucket(sname)
        if schema == "gspa.agent-public/2":
            b["public_a"] = data
        elif schema == "gspa.locus-public/1":
            b["public_b"] = data
        elif schema == "gspa.locus-private/1":
            b["private_b"] = data
        elif schema == "gspa.locus-truth/1":
            b["truth_b"] = data
        elif schema in ("gspa.row11-truth-report/1", "gspa.price-contract-row11-truth/1"):
            b["truth_report"] = data
        elif schema == "gspa.locus-proposals/1":
            b["proposals_b"][data.get("mode", p.stem)] = data
        elif schema == "gspa.locus-policy/1":
            b["policies_b"][data.get("mode", p.stem)] = data
        elif schema == "gspa.locus-evaluation/1":
            b["evals_b"][data.get("mode", p.stem)] = data
        elif schema == "gspa.locus-leak-check/1":
            b["leak_b"] = data
        elif schema == "gspa.locus-compare/1":
            b["compare_b"] = data
        elif schema == "gspa.agent-policy/2":
            b["policies_a"][data.get("mode", p.stem)] = data
        elif schema == "gspa.agent-evaluation/2":
            b["evals_a"][data.get("mode", p.stem)] = data

    return sets


def render_tables(results_dir: Path) -> str:
    sets = _discover_artifacts(Path(results_dir))
    lines: list[str] = ["# Row-11 EC-Aware Comparison Tables", ""]

    # (a) Loci coverage per split
    lines.append("## (a) Universe B Loci Coverage per Split")
    lines.append("")
    lines.append("| set | split | case (labelled) | resolved_single_call | no_candidate | unresolved_not_a_case | decoy | bundle_total | shortfalls |")
    lines.append("|---|---|---:|---:|---:|---:|---:|---:|---:|")
    for sname, b in sorted(sets.items()):
        priv = b["private_b"]
        pub = b["public_b"]
        if priv and "coverage" in priv and "loci" in priv["coverage"]:
            by_sp = priv["coverage"]["loci"].get("by_split", {})
            for sp in ("calibration", "test"):
                if sp not in by_sp:
                    continue
                row = by_sp[sp]
                n_case = row.get("case", row.get("labelled", 0))
                n_dec = row.get("decoy", row.get("n_bundle_cases", 0) - n_case)
                lines.append(
                    f"| {sname} | {sp} | {n_case} | {row.get('resolved_single_call', 0)} | "
                    f"{row.get('no_candidate', 0)} | {row.get('unresolved_not_a_case', 0)} | "
                    f"{n_dec} | {row.get('n_bundle_cases', n_case + n_dec)} | {row.get('shortfalls', 0)} |"
                )
        elif pub and "coverage" in pub and "loci" in pub["coverage"]:
            by_sp = pub["coverage"]["loci"].get("by_split", {})
            for sp in ("calibration", "test"):
                if sp not in by_sp:
                    continue
                row = by_sp[sp]
                n_case = row.get("case", 0)
                n_tot = row.get("n_bundle_cases", n_case)
                n_dec = n_tot - n_case
                lines.append(
                    f"| {sname} | {sp} | {n_case} | {row.get('resolved_single_call', 0)} | "
                    f"{row.get('no_candidate', 0)} | {row.get('unresolved_not_a_case', 0)} | "
                    f"{n_dec} | {n_tot} | {row.get('shortfalls', 0)} |"
                )
    lines.append("")

    # (b) Truth among candidates per split
    lines.append("## (b) Universe B Truth Among Candidates per Split")
    lines.append("")
    lines.append("| set | split | labelled_cases | truth_among_candidates | rank_1 | rank_2_10 | absent | circular_excluded_cases |")
    lines.append("|---|---|---:|---:|---:|---:|---:|---:|")
    for sname, b in sorted(sets.items()):
        tb = b["truth_b"] or (b["truth_report"] or {}).get("universe_b")
        if not tb:
            continue
        cov = tb.get("candidate_truth_coverage") or tb.get("counts_by_split") or {}
        for sp in ("calibration", "test", "total"):
            if sp not in cov:
                continue
            r = cov[sp]
            lines.append(
                f"| {sname} | {sp} | {r.get('n_labelled_cases', 0)} | {r.get('truth_among_candidates', 0)} | "
                f"{r.get('rank_1', 0)} | {r.get('rank_2_10', 0)} | {r.get('absent', 0)} | "
                f"{r.get('circular_excluded_cases', 0)} |"
            )
    lines.append("")

    # (c) Per method and set on Universe B
    lines.append("## (c) Universe B Method Comparison")
    lines.append("")
    lines.append(
        "| set | method | threshold | cal_known | cal_precision | test_correct | test_harmful | "
        "test_unconfirmed | test_abstained | correct_rate (Wilson 95%) | harm_rate (Wilson 95%) | "
        "recall_truth_ec | unknown_call_rate | leak_verdict | decision_wall_s | prompt_tokens | completion_tokens |"
    )
    lines.append("|---|---|---:|---:|---:|---:|---:|---:|---:|---|---|---:|---:|---|---:|---:|---:|")
    for sname, b in sorted(sets.items()):
        leak_by_mode = (b["leak_b"] or {}).get("by_mode", {})
        present_modes = set(b["evals_b"]) | set(b["policies_b"]) | set(b["proposals_b"]) | set(leak_by_mode)
        modes = [m for m in METHOD_ORDER if m in present_modes]
        for m in sorted(present_modes):
            if m not in modes:
                modes.append(m)
        for mode in modes:
            pol = b["policies_b"].get(mode, {})
            ev = b["evals_b"].get(mode, {})
            prop = b["proposals_b"].get(mode, {})
            lk = leak_by_mode.get(mode, {})
            lc = ev.get("labelled_counts") or ev.get("counts") or {}
            usage = ev.get("api_usage") or (ev.get("cost") or {}).get("api_usage") or prop.get("api_usage") or {}
            p_tok = usage.get("prompt_tokens", "n/a") if usage else "n/a"
            c_tok = usage.get("completion_tokens", "n/a") if usage else "n/a"
            wall_s = ev.get("decision_wall_seconds", prop.get("decision_wall_seconds"))
            lines.append(
                f"| {sname} | {mode} | {_fmt_float(pol.get('threshold'), 4)} | "
                f"{pol.get('calibration_known_decisions', 0)} | {_fmt_float(pol.get('calibration_precision'), 3)} | "
                f"{lc.get('correct', 0)} | {lc.get('harmful', 0)} | {lc.get('unconfirmed', 0)} | {lc.get('abstained', 0)} | "
                f"{_fmt_ci(ev.get('correct_rate'), ev.get('correct_rate_ci95'))} | "
                f"{_fmt_ci(ev.get('harm_rate'), ev.get('harm_rate_ci95'))} | "
                f"{_fmt_float(ev.get('recall_of_truth_ec'), 3)} | "
                f"{_fmt_float(ev.get('unknown_call_rate'), 3)} | "
                f"{lk.get('verdict', 'n/a')} (p={_fmt_float(lk.get('fisher_exact_p'), 3)}, diff={_fmt_float(lk.get('abs_diff'), 3)}) | "
                f"{_fmt_float(wall_s, 2)} | {p_tok} | {c_tok} |"
            )
    lines.append("")

    # (d) Universe A
    lines.append("## (d) Universe A Production Bundle and Method Counts")
    lines.append("")
    lines.append("### Universe A Cases and Truth Intersection")
    lines.append("")
    lines.append("| set | split | cases | case_types | candidate_pairs_on_labelled | positive_price | positive_swissprot | circular_excluded |")
    lines.append("|---|---|---:|---|---:|---:|---:|---:|")
    for sname, b in sorted(sets.items()):
        tr = b["truth_report"] or {}
        ua = tr.get("universe_a") if "universe_a" in tr else (tr if tr.get("schema") == "gspa.price-contract-row11-truth/1" else None)
        if not ua:
            continue
        by_sp_cases = ua.get("cases_by_split", {})
        by_sp_types = ua.get("cases_by_type_and_split", {})
        by_sp_pairs = ua.get("pairs_by_split", {})
        for sp in ("calibration", "test"):
            p_sp = by_sp_pairs.get(sp, {})
            t_sp = by_sp_types.get(sp) or ua.get("case_types", {})
            t_str = ", ".join(f"{k}:{v}" for k, v in sorted(t_sp.items())) or "none"
            lines.append(
                f"| {sname} | {sp} | {by_sp_cases.get(sp, 0)} | {t_str} | "
                f"{p_sp.get('candidate_pairs_on_labelled_proteins', 0)} | "
                f"{p_sp.get('positive/price', 0)} | {p_sp.get('positive/swissprot', 0)} | "
                f"{p_sp.get('positive/swissprot_circular', 0)} |"
            )
    lines.append("")
    lines.append("### Universe A Evaluation Counts per Method")
    lines.append("")
    lines.append("| set | method | threshold | test_cases | correct | harmful | unknown | abstained | unchanged | error | decision_wall_s |")
    lines.append("|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for sname, b in sorted(sets.items()):
        modes = [m for m in METHOD_ORDER if m in b["evals_a"] or m in b["policies_a"]]
        for m in sorted(set(b["evals_a"]) | set(b["policies_a"])):
            if m not in modes:
                modes.append(m)
        for mode in modes:
            pol = b["policies_a"].get(mode, {})
            ev = b["evals_a"].get(mode, {})
            cnt = ev.get("counts", {})
            lines.append(
                f"| {sname} | {mode} | {_fmt_float(pol.get('threshold'), 4)} | {ev.get('n_cases', 0)} | "
                f"{cnt.get('correct', 0)} | {cnt.get('harmful', 0)} | {cnt.get('unknown', 0)} | "
                f"{cnt.get('abstained', 0)} | {cnt.get('unchanged', 0)} | {cnt.get('error', 0)} | "
                f"{_fmt_float(ev.get('decision_wall_seconds'), 2)} |"
            )
    lines.append("")

    # (e) Paired contrasts
    lines.append("## (e) Universe B Paired Contrasts")
    lines.append("")
    lines.append(
        "| set | reference | comparator | n_test_loci | n_genomes | "
        "correct_diff (case bootstrap 95%) | correct_diff (genome bootstrap 95%) | "
        "harm_diff (case bootstrap 95%) | harm_diff (genome bootstrap 95%) |"
    )
    lines.append("|---|---|---|---:|---:|---|---|---|---|")
    for sname, b in sorted(sets.items()):
        cmp_b = b["compare_b"]
        if not cmp_b:
            continue
        for c in cmp_b.get("contrasts", []):
            cd_case = c["correct_rate_diff"]["paired_bootstrap"]
            cd_gen = c["correct_rate_diff"]["genome_blocked_bootstrap"]
            hd_case = c["harm_rate_diff"]["paired_bootstrap"]
            hd_gen = c["harm_rate_diff"]["genome_blocked_bootstrap"]
            lines.append(
                f"| {sname} | {c['reference_mode']} | {c['comparator_mode']} | "
                f"{c['n_labelled_test_cases']} | {c['n_test_genomes']} | "
                f"{_fmt_ci(cd_case['estimate'], [cd_case['ci_low'], cd_case['ci_high']])} | "
                f"{_fmt_ci(cd_gen['estimate'], [cd_gen['ci_low'], cd_gen['ci_high']])} | "
                f"{_fmt_ci(hd_case['estimate'], [hd_case['ci_low'], hd_case['ci_high']])} | "
                f"{_fmt_ci(hd_gen['estimate'], [hd_gen['ci_low'], hd_gen['ci_high']])} |"
            )
    lines.append("")
    return "\n".join(lines)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args(argv)

    md = render_tables(args.results)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("x", encoding="utf-8") as fh:
        fh.write(md)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
