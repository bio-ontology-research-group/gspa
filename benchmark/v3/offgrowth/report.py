#!/usr/bin/env python3
"""Plain-word, table-first Markdown report builder for off-growth enzyme discovery (§6-§9).

Follows the project reporting rules:
- All quantitative results live in Markdown tables.
- Prose uses plain words to explain what the tables show and what their limits are,
  without repeating numbers that the tables already carry.
- Negative findings and unmapped/unrouted panel entries are stated plainly.
"""

from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence

from .arms import SolverInfo
from .panel import UnmappedPanelEntry
from .score import (
    ArmSummaryMetrics,
    CaseArmEvaluation,
    H3EfficiencyResult,
    McNemarResult,
    ShuffledControlSummary,
    SolverReplaySummary,
)


ARM_DISPLAY_ORDER: Sequence[str] = (
    "A0",
    "A1",
    "A2",
    "A3",
    "A1_A2_A3",
    "A4@10%",
    "A4@25%",
)


def _fmt_rate_ci(rate: float, ci: Sequence[float]) -> str:
    lo, hi = (ci[0], ci[1]) if len(ci) >= 2 else (0.0, 0.0)
    return f"{rate:.3f} [{lo:.3f}, {hi:.3f}]"


def _fmt_opt_float(val: Optional[float], digits: int = 1) -> str:
    if val is None:
        return "NA"
    return f"{val:.{digits}f}"


def render_experiment_markdown_report(
    arm_summaries: Mapping[str, ArmSummaryMetrics],
    h1_result: McNemarResult,
    h2_results: Sequence[McNemarResult],
    h3_result: H3EfficiencyResult,
    shuffled_control: ShuffledControlSummary,
    replay_summary: SolverReplaySummary,
    solver_info: SolverInfo,
    unmapped_panel_entries: Sequence[UnmappedPanelEntry] = (),
    unmasked_summary: Optional[Mapping[str, ArmSummaryMetrics]] = None,
    masked_candidate_summary: Optional[Mapping[str, ArmSummaryMetrics]] = None,
) -> str:
    """Render the complete table-first Markdown report in plain words."""
    ordered_arms = [a for a in ARM_DISPLAY_ORDER if a in arm_summaries] + [
        a for a in sorted(arm_summaries) if a not in ARM_DISPLAY_ORDER
    ]

    lines: List[str] = []
    lines.append("# Off-Growth Enzyme Discovery Report")
    lines.append("")
    lines.append(
        "This report summarizes the pre-registered experiment testing whether targeted "
        "gap-filling beyond single-medium biomass growth recovers experimentally verified "
        "bacterial enzymes after per-case masking."
    )
    lines.append("")

    # Section 1: Recovery by Stratum
    lines.append("## 1. Enzyme Recovery by Stratum")
    lines.append("")
    lines.append(
        "Each truth case is classified before evaluation into the on-growth stratum "
        "(deleting the truth reaction from the unmasked model drops default-medium "
        "biomass growth) or the off-growth stratum (default-medium growth persists "
        "without the truth reaction). Confidence intervals come from the genome-block "
        "bootstrap."
    )
    lines.append("")
    lines.append(
        "| Arm | Total Cases | Total Recovered | Total Recovery [95% CI] | "
        "On-Growth Cases | On-Growth Recovered | On-Growth Recovery [95% CI] | "
        "Off-Growth Cases | Off-Growth Recovered | Off-Growth Recovery [95% CI] |"
    )
    lines.append("| --- | ---: | ---: | --- | ---: | ---: | --- | ---: | ---: | --- |")
    for arm in ordered_arms:
        m = arm_summaries[arm]
        lines.append(
            f"| `{m.arm}` | {m.n_cases_total} | {m.n_recovered_total} | "
            f"{_fmt_rate_ci(m.recovery_rate_total, m.recovery_ci_total)} | "
            f"{m.n_cases_on_growth} | {m.n_recovered_on_growth} | "
            f"{_fmt_rate_ci(m.recovery_rate_on_growth, m.recovery_ci_on_growth)} | "
            f"{m.n_cases_off_growth} | {m.n_recovered_off_growth} | "
            f"{_fmt_rate_ci(m.recovery_rate_off_growth, m.recovery_ci_off_growth)} |"
        )
    lines.append("")

    has_p = any(m.n_cases_price148 > 0 for m in arm_summaries.values())
    has_sp = any(m.n_cases_swissprot > 0 for m in arm_summaries.values())
    if has_p and has_sp:
        lines.append("### Recovery by Truth Case Source Stratum")
        lines.append("")
        lines.append(
            "Stratification of recovery across the Price-148 truth set and the "
            "Swiss-Prot experimental enzyme truth set."
        )
        lines.append("")
        lines.append(
            "| Arm | Case Source | Cases | Recovered | Recovery [95% CI] | "
            "Off-Growth Cases | Off-Growth Recovered | Off-Growth Recovery [95% CI] |"
        )
        lines.append("| --- | --- | ---: | ---: | --- | ---: | ---: | --- |")
        for arm in ordered_arms:
            m = arm_summaries[arm]
            lines.append(
                f"| `{m.arm}` | `price148` | {m.n_cases_price148} | {m.n_recovered_price148} | "
                f"{_fmt_rate_ci(m.recovery_rate_price148, m.recovery_ci_price148)} | "
                f"{m.n_cases_price148_off_growth} | {m.n_recovered_price148_off_growth} | "
                f"{_fmt_rate_ci(m.recovery_rate_price148_off_growth, m.recovery_ci_price148_off_growth)} |"
            )
            lines.append(
                f"| `{m.arm}` | `swissprot` | {m.n_cases_swissprot} | {m.n_recovered_swissprot} | "
                f"{_fmt_rate_ci(m.recovery_rate_swissprot, m.recovery_ci_swissprot)} | "
                f"{m.n_cases_swissprot_off_growth} | {m.n_recovered_swissprot_off_growth} | "
                f"{_fmt_rate_ci(m.recovery_rate_swissprot_off_growth, m.recovery_ci_swissprot_off_growth)} |"
            )
        lines.append("")
    elif has_sp and not has_p:
        lines.append("Truth case source: `swissprot` (Swiss-Prot experimental enzymes).")
        lines.append("")
    elif has_p and not has_sp:
        lines.append("Truth case source: `price148` (Price et al. 2018).")
        lines.append("")

    # Section 2: Call Quality, Harm, and Unconfirmed Calls
    lines.append("## 2. Call Quality, Harm, and Unconfirmed Calls")
    lines.append("")
    lines.append(
        "Added reactions are assigned to the top-scoring masked candidate gene for a "
        "matching enzyme class. Calls on proteins with known experimental labels are "
        "scored as supported when the reaction matches the experimental enzyme class "
        "and harmful when it disagrees. Calls on proteins without experimental labels "
        "are listed as unconfirmed and are never counted as negatives."
    )
    lines.append("")
    lines.append(
        "| Arm | Gene Calls | Reaction-Only Calls | Supported Calls | Harmful Calls | "
        "Unconfirmed Calls | Harm Rate [95% CI] | Exceeds Baseline Harm Gate |"
    )
    lines.append("| --- | ---: | ---: | ---: | ---: | ---: | --- | --- |")
    for arm in ordered_arms:
        m = arm_summaries[arm]
        gate_flag = "yes" if m.exceeds_a0_harm_gate else "no"
        lines.append(
            f"| `{m.arm}` | {m.n_gene_calls} | {m.n_reaction_only_calls} | "
            f"{m.n_supported_calls} | {m.n_harmful_calls} | "
            f"{m.n_unconfirmed_calls} | {_fmt_rate_ci(m.harm_rate, m.harm_ci)} | "
            f"{gate_flag} |"
        )
    lines.append("")

    # Section 3: Simulation and Linear Programming Cost
    lines.append("## 3. Simulation and Linear Programming Cost")
    lines.append("")
    lines.append(
        "Every linear program and weighted gap-fill solve is counted together with "
        "processor time, wall-clock time, and peak resident memory."
    )
    lines.append("")
    lines.append(
        "| Arm | Tests Considered | Tests Executed | Total LPs | Total Gap-Fills | "
        "CPU Seconds | Wall Seconds | Peak RSS (KiB) | LPs per Recovered Enzyme |"
    )
    lines.append("| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
    for arm in ordered_arms:
        m = arm_summaries[arm]
        lines.append(
            f"| `{m.arm}` | {m.total_tests_considered} | {m.total_tests_executed} | "
            f"{m.total_lps} | {m.total_gapfills} | {m.total_cpu_seconds:.3f} | "
            f"{m.total_wall_seconds:.3f} | {m.max_peak_rss_kb} | "
            f"{_fmt_opt_float(m.lps_per_recovered_enzyme, 1)} |"
        )
    lines.append("")

    # Section 4: Pre-Registered Statistical Tests
    lines.append("## 4. Pre-Registered Statistical Tests")
    lines.append("")
    lines.append(
        "The paired comparisons test whether the targeted gap-filling arms recover "
        "more off-growth enzymes than default-medium biomass gap-filling alone, and "
        "whether the deterministic rule ordering retains most of the off-growth "
        "recovery within a reduced test budget."
    )
    lines.append("")
    lines.append(
        "| Hypothesis | Treatment Arm | Baseline Arm | Stratum | Paired Cases | "
        "Both Recovered | Treatment Only | Baseline Only | Neither | Raw p-value | "
        "Holm p-value | Significant |"
    )
    lines.append("| --- | --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |")
    for res in [h1_result, *h2_results]:
        sig = "yes" if res.significant_at_0_05 else "no"
        lines.append(
            f"| `{res.hypothesis_id}` | `{res.arm_treatment}` | `{res.arm_baseline}` | "
            f"`{res.stratum}` | {res.n_pairs} | {res.both_recovered} | "
            f"{res.treatment_only} | {res.baseline_only} | {res.neither_recovered} | "
            f"{res.p_value_raw:.4g} | {res.p_value_holm:.4g} | {sig} |"
        )
    lines.append("")
    lines.append(
        "| Budget Arm | Full Arm | Off-Growth Recovered (Budget) | "
        "Off-Growth Recovered (Full) | Retention Ratio [95% CI] | "
        "Tests Executed / Considered | Budget Fraction | Meets Retention Target |"
    )
    lines.append("| --- | --- | ---: | ---: | --- | ---: | ---: | --- |")
    h3_pass = "yes" if h3_result.passes_80pct_retention else "no"
    lines.append(
        f"| `{h3_result.arm_budget}` | `{h3_result.arm_full}` | "
        f"{h3_result.off_growth_recovered_budget} | {h3_result.off_growth_recovered_full} | "
        f"{_fmt_rate_ci(h3_result.retention_ratio, h3_result.retention_ci)} | "
        f"{h3_result.tests_executed_budget} / {h3_result.tests_considered_full} | "
        f"{h3_result.test_budget_fraction:.3f} | {h3_pass} |"
    )
    lines.append("")

    # Section 5: Controls
    lines.append("## 5. Controls and Solver Replay")
    lines.append("")
    lines.append(
        "The shuffled-candidate control permutes candidate protein identifiers within "
        "each genome across fixed random seeds and reassigns calls on the added "
        "reactions. The solver replay control re-solves a deterministic sample of "
        "weighted gap-fill linear programs with the secondary solver."
    )
    lines.append("")
    lines.append(
        "| Shuffled Seed | Arm | Total Cases | Total Recovered | Total Recovery Rate | "
        "Off-Growth Cases | Off-Growth Recovered | Off-Growth Recovery Rate |"
    )
    lines.append("| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |")
    for item in shuffled_control.per_seed:
        lines.append(
            f"| {item.seed} | `{item.arm}` | {item.n_cases} | {item.n_recovered} | "
            f"{item.recovery_rate:.3f} | {item.n_off_growth_cases} | "
            f"{item.n_off_growth_recovered} | {item.off_growth_recovery_rate:.3f} |"
        )
    lines.append("")
    lines.append(
        "| Primary Solver | Replay Solver | Highspy Installed | Replayed Gap-Fills | "
        "Agreed Gap-Fills | Agreement Rate | Max Relative Objective Diff |"
    )
    lines.append("| --- | --- | --- | ---: | ---: | ---: | ---: |")
    hi_flag = "yes" if solver_info.highspy_available else "no"
    lines.append(
        f"| `{solver_info.primary_solver}` | `{solver_info.replay_solver}` | {hi_flag} | "
        f"{replay_summary.n_replayed} | {replay_summary.n_agreed} | "
        f"{replay_summary.agreement_rate:.3f} | {replay_summary.max_relative_obj_diff:.3e} |"
    )
    lines.append("")

    if unmasked_summary:
        lines.append("### Unmasked Reachability Control")
        lines.append("")
        lines.append(
            "Running every arm on the unmasked draft model with unmasked candidates "
            "and scoring the same way separates reaction reachability in the universal "
            "database from masking loss."
        )
        lines.append("")
        lines.append(
            "| Arm | Total Cases | Unmasked Recovered | Unmasked Recovery Rate | "
            "Off-Growth Cases | Unmasked Off-Growth Recovered | Unmasked Off-Growth Rate |"
        )
        lines.append("| --- | ---: | ---: | ---: | ---: | ---: | ---: |")
        for arm in ordered_arms:
            if arm not in unmasked_summary:
                continue
            um = unmasked_summary[arm]
            lines.append(
                f"| `{um.arm}` | {um.n_cases_total} | {um.n_recovered_total} | "
                f"{um.recovery_rate_total:.3f} | {um.n_cases_off_growth} | "
                f"{um.n_recovered_off_growth} | {um.recovery_rate_off_growth:.3f} |"
            )
        lines.append("")

    if masked_candidate_summary:
        lines.append("### Masked-Candidate Attribution Control")
        lines.append("")
        lines.append(
            "Scoring the masked arm runs against the masked candidate table (where the "
            "target protein's truth enzyme claim has been deleted) yields zero recovery "
            "by construction."
        )
        lines.append("")
        lines.append(
            "| Arm | Total Cases | Masked-Candidate Recovered | Masked-Candidate Recovery Rate | "
            "Off-Growth Cases | Masked-Candidate Off-Growth Recovered | Masked-Candidate Off-Growth Rate |"
        )
        lines.append("| --- | ---: | ---: | ---: | ---: | ---: | ---: |")
        for arm in ordered_arms:
            if arm not in masked_candidate_summary:
                continue
            mc = masked_candidate_summary[arm]
            lines.append(
                f"| `{mc.arm}` | {mc.n_cases_total} | {mc.n_recovered_total} | "
                f"{mc.recovery_rate_total:.3f} | {mc.n_cases_off_growth} | "
                f"{mc.n_recovered_off_growth} | {mc.recovery_rate_off_growth:.3f} |"
            )
        lines.append("")

    # Section 6: Unmapped or Unrouted Source Panel Compounds
    lines.append("## 6. Unmapped or Unrouted Source Panel Compounds")
    lines.append("")
    lines.append(
        "Every compound from the input source panel that could not be mapped to a "
        "ModelSEED identifier with a verified elemental formula and universal "
        "transport route is listed below rather than dropped silently."
    )
    lines.append("")
    lines.append("| Plate | Well | Element | Compound Name | SEED ID | Reason | Detail |")
    lines.append("| --- | --- | --- | --- | --- | --- | --- |")
    if unmapped_panel_entries:
        for u in unmapped_panel_entries:
            lines.append(
                f"| `{u.plate}` | `{u.well}` | `{u.element}` | {u.compound_name} | "
                f"`{u.seed_id or 'NA'}` | `{u.reason}` | {u.detail} |"
            )
    else:
        lines.append("| NA | NA | NA | None | NA | none | All panel rows mapped and routed |")
    lines.append("")

    return "\n".join(lines)


def write_experiment_outputs(
    out_dir: Path,
    arm_summaries: Mapping[str, ArmSummaryMetrics],
    evaluations_by_arm: Mapping[str, Sequence[CaseArmEvaluation]],
    h1_result: McNemarResult,
    h2_results: Sequence[McNemarResult],
    h3_result: H3EfficiencyResult,
    shuffled_control: ShuffledControlSummary,
    replay_summary: SolverReplaySummary,
    solver_info: SolverInfo,
    unmapped_panel_entries: Sequence[UnmappedPanelEntry] = (),
    unmasked_summary: Optional[Mapping[str, ArmSummaryMetrics]] = None,
    masked_candidate_summary: Optional[Mapping[str, ArmSummaryMetrics]] = None,
    markdown_path: Optional[Path] = None,
) -> Path:
    """Write JSON summary and Markdown report to `out_dir`."""
    out_dir.mkdir(parents=True, exist_ok=True)
    md_text = render_experiment_markdown_report(
        arm_summaries=arm_summaries,
        h1_result=h1_result,
        h2_results=h2_results,
        h3_result=h3_result,
        shuffled_control=shuffled_control,
        replay_summary=replay_summary,
        solver_info=solver_info,
        unmapped_panel_entries=unmapped_panel_entries,
        unmasked_summary=unmasked_summary,
        masked_candidate_summary=masked_candidate_summary,
    )
    report_dest = markdown_path if markdown_path is not None else (out_dir / "report.md")
    report_dest.parent.mkdir(parents=True, exist_ok=True)
    report_dest.write_text(md_text, encoding="utf-8")

    payload = {
        "solver_info": asdict(solver_info),
        "arm_summaries": {k: asdict(v) for k, v in arm_summaries.items()},
        "h1": asdict(h1_result),
        "h2": [asdict(x) for x in h2_results],
        "h3": asdict(h3_result),
        "shuffled_control": asdict(shuffled_control),
        "solver_replay": asdict(replay_summary),
        "unmapped_panel": [asdict(u) for u in unmapped_panel_entries],
        "evaluations_by_arm": {
            arm: [asdict(ev) for ev in evals]
            for arm, evals in evaluations_by_arm.items()
        },
    }
    if unmasked_summary is not None:
        payload["unmasked_summary"] = {k: asdict(v) for k, v in unmasked_summary.items()}
    if masked_candidate_summary is not None:
        payload["masked_candidate_summary"] = {
            k: asdict(v) for k, v in masked_candidate_summary.items()
        }

    summary_json = out_dir / "summary.json"
    summary_json.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return report_dest
