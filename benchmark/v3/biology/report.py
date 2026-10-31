"""Render Markdown tables and plain-word summaries for GSPA biology checks.

All quantitative results live strictly inside Markdown tables (genome x predictor
x arm, plus metabolic stratum breakdown for experimental truth agreement).
Prose paragraphs state the biological rule, qualitative findings, and scope
limits in plain words without repeating numbers from the tables.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from typing import Any, Dict, List, Optional

V3_DIR = Path(__file__).resolve().parents[1]
if str(V3_DIR) not in sys.path:
    sys.path.insert(0, str(V3_DIR))

from biology.checks import (  # noqa: E402
    STRATA_ORDER,
    evaluate_collection_summaries,
    load_summaries_from_path,
)


def _fmt_list(items: List[str]) -> str:
    if not items:
        return "none"
    return ", ".join(items)


def _fmt_bool(val: Optional[bool]) -> str:
    if val is True:
        return "yes"
    if val is False:
        return "no"
    return "unknown"


def render_translation_table(cells: List[Dict[str, Any]]) -> str:
    lines = [
        "| Genome | Organism | Exempt | Predictor | Arm | Call-Threshold Routes | Call-Threshold Gat Routes | Call-Threshold Missing | Retained Routes | Retained Gat Routes | Retained Missing | Pass (Call) |",
        "| --- | --- | --- | --- | --- | ---: | ---: | --- | ---: | ---: | --- | --- |",
    ]
    for cell in cells:
        asm = cell["assembly"]
        pred = cell["predictor"]
        ex = cell.get("exemption", {})
        org = ex.get("organism") or "unknown"
        exempt_str = "yes" if ex.get("exempt") else "no"
        for arm_name in ("stage0", "stage1"):
            tr = cell["arms"][arm_name]["translation"]
            cal = tr["called"]
            ret = tr["retained"]
            cal_ratio = f"{cal['satisfied_amino_acids']}/{cal['total_amino_acids']}"
            ret_ratio = f"{ret['satisfied_amino_acids']}/{ret['total_amino_acids']}"
            lines.append(
                f"| {asm} | {org} | {exempt_str} | {pred} | {arm_name} | "
                f"{cal_ratio} | {cal['indirect_gat_count']} | {_fmt_list(cal['missing_amino_acids'])} | "
                f"{ret_ratio} | {ret['indirect_gat_count']} | {_fmt_list(ret['missing_amino_acids'])} | "
                f"{_fmt_bool(cal['passed'])} |"
            )
    return "\n".join(lines)


def render_core_metabolism_table(cells: List[Dict[str, Any]]) -> str:
    lines = [
        "| Genome | Organism | Exempt | Predictor | Arm | Call-Threshold Present | Call-Threshold Missing | Retained Present | Retained Missing | Pass (Call) |",
        "| --- | --- | --- | --- | --- | ---: | --- | ---: | --- | --- |",
    ]
    for cell in cells:
        asm = cell["assembly"]
        pred = cell["predictor"]
        ex = cell.get("exemption", {})
        org = ex.get("organism") or "unknown"
        exempt_str = "yes" if ex.get("exempt") else "no"
        for arm_name in ("stage0", "stage1"):
            cm = cell["arms"][arm_name]["core_metabolism"]
            cal = cm["called"]
            ret = cm["retained"]
            cal_ratio = f"{cal['present_roles']}/{cal['total_roles']}"
            ret_ratio = f"{ret['present_roles']}/{ret['total_roles']}"
            lines.append(
                f"| {asm} | {org} | {exempt_str} | {pred} | {arm_name} | "
                f"{cal_ratio} | {_fmt_list(cal['missing_roles'])} | "
                f"{ret_ratio} | {_fmt_list(ret['missing_roles'])} | "
                f"{_fmt_bool(cal['passed'])} |"
            )
    return "\n".join(lines)


def render_taxon_consistency_table(cells: List[Dict[str, Any]]) -> str:
    lines = [
        "| Genome | Predictor | Arm | Enforcement On | Enforce Stage2 | Symbolic Stage | Taxon In Hierarchy | Priors Fired Rows | Call Violating Rows | Call Violating ECs | Retained Violating Rows | Retained Violating ECs |",
        "| --- | --- | --- | --- | --- | --- | --- | ---: | ---: | ---: | ---: | ---: |",
    ]
    for cell in cells:
        asm = cell["assembly"]
        pred = cell["predictor"]
        for arm_name in ("stage0", "stage1"):
            tc = cell["arms"][arm_name]["taxon_consistency"]
            enf_on = _fmt_bool(tc.get("enforcement_on"))
            enf_s2 = _fmt_bool(tc.get("enforce_stage2"))
            sym_st = tc.get("symbolic_enforcement_status") or "unknown"
            in_hier = _fmt_bool(tc.get("genome_taxon_in_bundled_hierarchy"))
            lines.append(
                f"| {asm} | {pred} | {arm_name} | {enf_on} | {enf_s2} | {sym_st} | {in_hier} | "
                f"{tc['priors_fired_rows']} | {tc['violation_rows_called']} | "
                f"{tc['violating_distinct_ecs_called']} | {tc['violation_rows_retained']} | "
                f"{tc['violating_distinct_ecs_retained']} |"
            )
    return "\n".join(lines)


def render_truth_agreement_cell_table(cells: List[Dict[str, Any]]) -> str:
    lines = [
        "| Genome | Predictor | Arm | Truth Proteins | Retained Four-Level Exact | Retained Four-Level Canonical | Retained Three-Level | Call Four-Level Exact | Call Four-Level Canonical | Call Three-Level |",
        "| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for cell in cells:
        asm = cell["assembly"]
        pred = cell["predictor"]
        for arm_name in ("stage0", "stage1"):
            ta = cell["arms"][arm_name]["truth_agreement"]
            ret = ta["retained"]["overall"]
            cal = ta["called"]["overall"]
            n_tot = ret["n_truth_proteins"]
            lines.append(
                f"| {asm} | {pred} | {arm_name} | {n_tot} | "
                f"{ret['ec4_exact_matches']}/{n_tot} | {ret['ec4_canonical_matches']}/{n_tot} | {ret['ec3_matches']}/{n_tot} | "
                f"{cal['ec4_exact_matches']}/{n_tot} | {cal['ec4_canonical_matches']}/{n_tot} | {cal['ec3_matches']}/{n_tot} |"
            )
    return "\n".join(lines)


def render_truth_agreement_stratum_table(
    truth_by_predictor: Dict[str, Any]
) -> str:
    lines = [
        "| Predictor | Arm | Metabolic Stratum | Truth Proteins | Retained Four-Level Exact | Retained Four-Level Canonical | Retained Three-Level | Call Four-Level Exact | Call Four-Level Canonical | Call Three-Level |",
        "| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for pred in sorted(truth_by_predictor.keys()):
        for arm_name in ("stage0", "stage1"):
            ret_block = truth_by_predictor[pred][arm_name]["retained"]
            cal_block = truth_by_predictor[pred][arm_name]["called"]
            for stratum in STRATA_ORDER:
                r_st = ret_block["by_stratum"][stratum]
                c_st = cal_block["by_stratum"][stratum]
                n_tot = r_st["n_truth_proteins"]
                lines.append(
                    f"| {pred} | {arm_name} | {stratum} | {n_tot} | "
                    f"{r_st['ec4_exact_matches']}/{n_tot} | {r_st['ec4_canonical_matches']}/{n_tot} | {r_st['ec3_matches']}/{n_tot} | "
                    f"{c_st['ec4_exact_matches']}/{n_tot} | {c_st['ec4_canonical_matches']}/{n_tot} | {c_st['ec3_matches']}/{n_tot} |"
                )
            r_all = ret_block["overall"]
            c_all = cal_block["overall"]
            n_all = r_all["n_truth_proteins"]
            lines.append(
                f"| {pred} | {arm_name} | All_Strata | {n_all} | "
                f"{r_all['ec4_exact_matches']}/{n_all} | {r_all['ec4_canonical_matches']}/{n_all} | {r_all['ec3_matches']}/{n_all} | "
                f"{c_all['ec4_exact_matches']}/{n_all} | {c_all['ec4_canonical_matches']}/{n_all} | {c_all['ec3_matches']}/{n_all} |"
            )
    return "\n".join(lines)


def _qualitative_overview(evaluation: Dict[str, Any]) -> List[str]:
    """Generate plain-word qualitative summary lines with zero numbers."""
    cells = evaluation.get("cells", [])
    if not cells:
        return [
            "No extracted cell summaries were supplied for evaluation.",
        ]

    all_identical = all(c.get("stage0_equals_stage1_sha256", False) for c in cells)
    any_exempt = any(c.get("exemption", {}).get("exempt", False) for c in cells)
    any_enforcement_on = any(
        c["arms"]["stage1"]["taxon_consistency"].get("enforcement_on", False)
        for c in cells
    )

    paragraphs: List[str] = []
    if all_identical:
        paragraphs.append(
            "Across every evaluated cell, the full integration arm produces the exact same output digest as the raw evidence arm because default prior adjustments and symbolic stage enforcement do not modify enzyme commission candidates in this configuration."
        )
    else:
        paragraphs.append(
            "The tables below compare the raw predictor evidence arm and the full integration arm side by side across each evaluated genome and predictor."
        )

    if any_exempt:
        paragraphs.append(
            "Genomes belonging to obligate intracellular endosymbionts or host-restricted parasites with extreme genome reduction are marked exempt in the completeness checks, whereas free-living bacteria are required to carry the complete translation and core metabolic repertoire."
        )
    else:
        paragraphs.append(
            "Every evaluated genome in this panel is a free-living heterotrophic bacterium; none is exempt from the translation or core metabolism completeness rules."
        )

    if not any_enforcement_on:
        paragraphs.append(
            "Cell configuration records show that symbolic consistency enforcement was disabled and skipped during integration, so taxonomic constraint violations present in the raw predictor candidates remain unchanged in the full integration output."
        )
    else:
        paragraphs.append(
            "Cell configuration records indicate that symbolic consistency enforcement was active in at least part of the evaluated collection."
        )

    return paragraphs


def render_markdown_report(evaluation: Dict[str, Any]) -> str:
    """Render the complete biology check report in Markdown."""
    cells = evaluation.get("cells", [])
    truth_by_pred = evaluation.get("truth_agreement_by_predictor", {})

    sections: List[str] = [
        "# Biology Check Report for GSPA Predictions",
        "",
        "## Overview",
        "",
    ]
    for para in _qualitative_overview(evaluation):
        sections.append(para)
        sections.append("")

    sections.extend(
        [
            "## Translation Charging Routes",
            "",
            "Free-living bacteria must carry an aminoacyl-tRNA charging route for every standard proteinogenic amino acid, accepting the indirect amidotransferase route for glutamine and asparagine when the direct synthetase is absent.",
            "",
            render_translation_table(cells),
            "",
            "## Core Information Processing and Central Metabolism",
            "",
            "Free-living bacteria are checked for universal enzymes spanning transcription, chromosomal replication, DNA topology and repair, the lower triose-phosphate glycolytic trunk, phosphoribosyl diphosphate synthesis, and proton-translocating ATP synthesis.",
            "",
            render_core_metabolism_table(cells),
            "",
            "## Taxon Consistency Against Bundled Bacterial Constraints",
            "",
            "Predictions are checked against the bundled taxonomic constraints for bacteria by mapping enzyme identifiers to ontology terms and propagating taxonomic ancestor and disjointness relations.",
            "",
            render_taxon_consistency_table(cells),
            "",
            "## Experimental Truth Agreement by Genome and Predictor",
            "",
            "For experimentally characterized benchmark proteins, the top predicted enzyme call is compared against the assigned experimental function at full specificity and at subclass prefix level; proteins outside the experimental truth set are unconfirmed and are never scored as wrong.",
            "",
            render_truth_agreement_cell_table(cells),
            "",
            "## Experimental Truth Agreement by Metabolic Stratum",
            "",
            "Experimental truth agreement is stratified into core biomass pathways, condition-specific catabolic assays, and peripheral secondary metabolism.",
            "",
            render_truth_agreement_stratum_table(truth_by_pred),
            "",
            "## Mutual Exclusivity Rule Disposition",
            "",
            "In international enzyme nomenclature a single identifier covers both forward and reverse catalysis of a reversible reaction, while distinct identifiers for opposing physiological steps represent chemically different reactions that can occur together on bifunctional polypeptides; therefore no direction-reversal mutual exclusivity failure rule is applied.",
            "",
        ]
    )
    return "\n".join(sections)


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="Path to an evaluated biology JSON file or a directory of extracted cell summaries",
    )
    parser.add_argument(
        "--out",
        type=Path,
        required=True,
        help="Output Markdown report path",
    )
    args = parser.parse_args(argv)

    in_path = Path(args.input)
    if in_path.is_file():
        raw = json.loads(in_path.read_text(encoding="utf-8"))
        if raw.get("schema") == "gspa.biology_collection_evaluation/1":
            evaluation = raw
        else:
            summaries = load_summaries_from_path(in_path)
            evaluation = evaluate_collection_summaries(summaries)
    else:
        summaries = load_summaries_from_path(in_path)
        evaluation = evaluate_collection_summaries(summaries)

    md = render_markdown_report(evaluation)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(md, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
