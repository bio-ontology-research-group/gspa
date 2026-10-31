"""Plain-unittest suite for benchmark/v3/biology (extract, checks, report, sbatch)."""
from __future__ import annotations

import os
import hashlib
import io
import json
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest

V3_DIR = Path(__file__).resolve().parents[2]
REPO_ROOT = V3_DIR.parents[1]
if str(V3_DIR) not in sys.path:
    sys.path.insert(0, str(V3_DIR))

from biology.checks import (  # noqa: E402
    CORE_BACTERIAL_EC_RULES,
    TRANSLATION_CHARGING_RULES,
    check_genome_exemption,
    evaluate_cell_summary,
    evaluate_collection_summaries,
    evaluate_core_metabolism_arm,
    evaluate_taxon_consistency_arm,
    evaluate_translation_arm,
    evaluate_truth_agreement_arm,
)
from biology.extract import (  # noqa: E402
    INTEGRATED_TSV_HEADER,
    extract_cell,
    extract_collection,
    load_bacterial_taxon_constraints,
    load_ec_audit_entries,
    load_price148_targets,
    stream_integrated_handle,
)
from biology.report import render_markdown_report  # noqa: E402


def _make_tsv_bytes(rows: list[tuple[str, str, str, str, str, str, str, str, str, str]]) -> bytes:
    lines = ["\t".join(INTEGRATED_TSV_HEADER)]
    for r in rows:
        lines.append("\t".join(r))
    return ("\n".join(lines) + "\n").encode("utf-8")


def _ec_row(
    protein_id: str,
    ec_id: str,
    prob: float,
    priors_fired: str = "",
) -> tuple[str, str, str, str, str, str, str, str, str, str]:
    return (
        protein_id,
        "EC",
        ec_id,
        "",
        f"{prob:.6f}",
        "0.100000",
        "0.100000",
        "1",
        priors_fired,
        "1",
    )


class TestBiologySuite(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.ec_entries = load_ec_audit_entries()
        cls.taxon_constraints = load_bacterial_taxon_constraints(
            ec_entries=cls.ec_entries
        )
        cls.price_targets = load_price148_targets()

    def test_rules_have_explicit_sources(self) -> None:
        self.assertEqual(len(TRANSLATION_CHARGING_RULES), 20)
        seen_aa = set()
        for rule in TRANSLATION_CHARGING_RULES:
            self.assertTrue(rule["amino_acid"])
            self.assertNotIn(rule["amino_acid"], seen_aa)
            seen_aa.add(rule["amino_acid"])
            self.assertTrue(rule["direct_ecs"])
            self.assertTrue(rule["source"] and len(rule["source"]) > 10)

        self.assertEqual(len(CORE_BACTERIAL_EC_RULES), 10)
        for rule in CORE_BACTERIAL_EC_RULES:
            self.assertTrue(rule["role_id"])
            self.assertTrue(rule["accepted_ecs"])
            self.assertTrue(rule["source"] and len(rule["source"]) > 10)

    def test_streaming_extraction_tie_break_and_prefix(self) -> None:
        rows = [
            _ec_row("WP_B", "EC:1.1.1.1", 0.80),
            _ec_row("WP_A", "1.1.1.1", 0.80),  # tie at 0.80 -> WP_A wins deterministically
            _ec_row("WP_C", "EC:1.1.1.2", 0.90),  # 3-level 1.1.1 best is 0.90 on WP_C
            _ec_row("WP_D", "EC:2.7.7.6", 0.40),  # retained only (< 0.5 call threshold)
            _ec_row("WP_E", "EC:5.99.1.2", 0.75),  # transfers to 5.6.2.1
            _ec_row("WP_F", "EC:1.13.12.7", 0.65),  # violates NCBITaxon_2 (firefly luciferase)
        ]
        raw_bytes = _make_tsv_bytes(rows)
        expected_sha = hashlib.sha256(raw_bytes).hexdigest()

        summary = stream_integrated_handle(
            io.BytesIO(raw_bytes),
            source_name="fixture.integrated.tsv",
            call_threshold=0.5,
            taxon_constraints=self.taxon_constraints,
            ec_entries=self.ec_entries,
            genome_taxon_id="211586",
        )

        self.assertEqual(summary["input_sha256"], expected_sha)
        self.assertEqual(summary["input_bytes"], len(raw_bytes))
        self.assertEqual(summary["total_rows"], 6)
        self.assertEqual(summary["called_rows"], 5)

        # Deterministic tie break on 1.1.1.1
        self.assertEqual(summary["ec4"]["called"]["1.1.1.1"], [0.8, "WP_A", 2])
        self.assertEqual(summary["ec4"]["retained"]["1.1.1.1"], [0.8, "WP_A", 2])

        # Three-level prefix 1.1.1 takes max across 1.1.1.1 and 1.1.1.2
        self.assertEqual(summary["ec3"]["called"]["1.1.1"], [0.9, "WP_C", 3])

        # Below-threshold 2.7.7.6 is in retained but not in called
        self.assertIn("2.7.7.6", summary["ec4"]["retained"])
        self.assertNotIn("2.7.7.6", summary["ec4"]["called"])

        # Transferred EC 5.99.1.2 populates canonical 5.6.2.1
        self.assertIn("5.6.2.1", summary["ec4"]["called"])
        self.assertEqual(summary["ec4"]["called"]["5.6.2.1"][1], "WP_E")

        # Taxon violation recorded for 1.13.12.7
        tc = summary["taxon_consistency"]
        self.assertTrue(tc["genome_taxon_in_bundled_hierarchy"])
        self.assertEqual(tc["violation_rows_called"], 1)
        self.assertEqual(tc["violation_rows_retained"], 1)
        self.assertIn("1.13.12.7", tc["violating_ecs_called"])

    def test_translation_gatcab_indirect_alternative_accepted(self) -> None:
        # Build a genome that lacks direct GlnRS (6.1.1.18) and direct AsnRS (6.1.1.22)
        # but carries indirect GatCAB/GatDE routes 6.3.5.7 (Gln) and 6.3.5.6 (Asn),
        # plus non-discriminating AspRS (6.1.1.23) and GluRS (6.1.1.24).
        direct_16 = [
            "6.1.1.1",
            "6.1.1.2",
            "6.1.1.3",
            "6.1.1.4",
            "6.1.1.5",
            "6.1.1.6",
            "6.1.1.7",
            "6.1.1.9",
            "6.1.1.10",
            "6.1.1.11",
            "6.1.1.14",
            "6.1.1.15",
            "6.1.1.16",
            "6.1.1.19",
            "6.1.1.20",
            "6.1.1.21",
        ]
        rows = [
            _ec_row(f"PROT_{idx}", ec, 0.85)
            for idx, ec in enumerate(direct_16)
        ]
        rows.append(_ec_row("PROT_ASP_ND", "6.1.1.23", 0.88))
        rows.append(_ec_row("PROT_GLU_ND", "6.1.1.24", 0.87))
        rows.append(_ec_row("PROT_GAT_ASN", "6.3.5.6", 0.91))
        rows.append(_ec_row("PROT_GAT_GLN", "6.3.5.7", 0.92))

        arm_summary = stream_integrated_handle(
            io.BytesIO(_make_tsv_bytes(rows)),
            call_threshold=0.5,
            ec_entries=self.ec_entries,
        )
        res = evaluate_translation_arm(
            arm_summary, level="called", ec_entries=self.ec_entries
        )

        self.assertTrue(res["passed"])
        self.assertEqual(res["satisfied_amino_acids"], 20)
        self.assertEqual(res["missing_count"], 0)
        self.assertEqual(res["indirect_gat_count"], 2)
        self.assertEqual(res["routes"]["Gln"]["route_mode"], "indirect_gat")
        self.assertEqual(res["routes"]["Gln"]["best_ec"], "6.3.5.7")
        self.assertEqual(res["routes"]["Asn"]["route_mode"], "indirect_gat")
        self.assertEqual(res["routes"]["Asn"]["best_ec"], "6.3.5.6")

        # Removing both 6.1.1.18 and 6.3.5.7 causes Gln to fail on a free-living genome
        rows_missing_gln = [r for r in rows if r[2] != "6.3.5.7"]
        arm_missing = stream_integrated_handle(
            io.BytesIO(_make_tsv_bytes(rows_missing_gln)),
            call_threshold=0.5,
            ec_entries=self.ec_entries,
        )
        res_missing = evaluate_translation_arm(
            arm_missing, level="called", ec_entries=self.ec_entries
        )
        self.assertFalse(res_missing["passed"])
        self.assertEqual(res_missing["missing_amino_acids"], ["Gln"])

    def test_exempt_genome_handling(self) -> None:
        # Free-living Price genome is not exempt
        ex_free = check_genome_exemption("GCA_000006965.1")
        self.assertFalse(ex_free["exempt"])

        # Obligate endosymbiont Buchnera (taxid 9 or organism genus Buchnera) is exempt
        ex_buchnera = check_genome_exemption(
            "GCF_000009605.1", taxon_id="9", organism="Buchnera aphidicola"
        )
        self.assertTrue(ex_buchnera["exempt"])
        self.assertIn("Buchnera", ex_buchnera["reason"])

        # Explicit caller exemption list is also honored
        ex_custom = check_genome_exemption(
            "GCF_TEST_SYMBIONT.1", exempt_assemblies=["GCF_TEST_SYMBIONT.1"]
        )
        self.assertTrue(ex_custom["exempt"])

        # An incomplete proteome passes completeness when exempt=True while still recording missing items
        empty_arm = stream_integrated_handle(
            io.BytesIO(_make_tsv_bytes([_ec_row("P1", "2.7.7.6", 0.90)])),
            call_threshold=0.5,
            ec_entries=self.ec_entries,
        )
        tr_exempt = evaluate_translation_arm(
            empty_arm, level="called", ec_entries=self.ec_entries, exemption=ex_buchnera
        )
        self.assertTrue(tr_exempt["exempt"])
        self.assertTrue(tr_exempt["passed"])
        self.assertEqual(tr_exempt["missing_count"], 20)

        cm_exempt = evaluate_core_metabolism_arm(
            empty_arm, level="called", ec_entries=self.ec_entries, exemption=ex_buchnera
        )
        self.assertTrue(cm_exempt["exempt"])
        self.assertTrue(cm_exempt["passed"])
        self.assertEqual(cm_exempt["present_roles"], 1)
        self.assertEqual(cm_exempt["missing_count"], 9)

    def test_core_metabolism_accepts_transferred_ecs(self) -> None:
        # Supply old pre-2018 ECs 5.99.1.2, 5.99.1.3, 3.6.3.14 alongside the other 7 core ECs
        core_ecs = [
            "2.7.7.6",
            "2.7.7.7",
            "5.99.1.2",
            "5.99.1.3",
            "6.5.1.2",
            "5.3.1.1",
            "2.7.2.3",
            "4.2.1.11",
            "2.7.6.1",
            "3.6.3.14",
        ]
        rows = [_ec_row(f"CORE_{i}", ec, 0.80) for i, ec in enumerate(core_ecs)]
        arm = stream_integrated_handle(
            io.BytesIO(_make_tsv_bytes(rows)),
            call_threshold=0.5,
            ec_entries=self.ec_entries,
        )
        res = evaluate_core_metabolism_arm(
            arm, level="called", ec_entries=self.ec_entries
        )
        self.assertTrue(res["passed"])
        self.assertEqual(res["present_roles"], 10)
        self.assertEqual(res["missing_roles"], [])

    def test_unconfirmed_labels_never_negative_and_stratum_split(self) -> None:
        targets = {
            "PROT_STRAT_A": {
                "price_entry": "PE_A",
                "host_assembly": "GCF_TEST.1",
                "proteome_protein_id": "PROT_STRAT_A",
                "ground_truth_ecs": ["5.3.1.7"],
                "ground_truth_ec3s": ["5.3.1"],
                "metabolic_stratum": "Stratum_A_Core_Biomass",
            },
            "PROT_STRAT_B": {
                "price_entry": "PE_B",
                "host_assembly": "GCF_TEST.1",
                "proteome_protein_id": "PROT_STRAT_B",
                "ground_truth_ecs": ["4.2.1.43"],
                "ground_truth_ec3s": ["4.2.1"],
                "metabolic_stratum": "Stratum_B_Catabolic_Assay_Specific",
            },
        }
        # PROT_STRAT_A gets exact 4-level match 5.3.1.7 at 0.90.
        # PROT_STRAT_B gets 4.2.1.25 at 0.70 (wrong 4-level, right 3-level prefix 4.2.1!).
        # Plus 50 unconfirmed proteins with high-confidence predictions that must NEVER count as negatives.
        rows = [
            _ec_row("PROT_STRAT_A", "5.3.1.7", 0.90),
            _ec_row("PROT_STRAT_A", "1.1.1.1", 0.20),
            _ec_row("PROT_STRAT_B", "4.2.1.25", 0.70),
            _ec_row("PROT_STRAT_B", "4.2.1.43", 0.40),
        ]
        for idx in range(50):
            rows.append(_ec_row(f"UNCONFIRMED_PROT_{idx}", "3.2.1.1", 0.99))

        arm = stream_integrated_handle(
            io.BytesIO(_make_tsv_bytes(rows)),
            call_threshold=0.5,
            target_proteins=targets,
            ec_entries=self.ec_entries,
        )
        truth_eval = evaluate_truth_agreement_arm(
            arm,
            manifest_targets_for_assembly=targets,
            level="called",
            ec_entries=self.ec_entries,
        )

        self.assertTrue(truth_eval["unconfirmed_never_negative"])
        overall = truth_eval["overall"]
        # Only the 2 experimental truth proteins are in the denominator; the 50 unconfirmed proteins are ignored
        self.assertEqual(overall["n_truth_proteins"], 2)
        self.assertEqual(overall["ec4_exact_matches"], 1)
        self.assertEqual(overall["ec3_matches"], 2)

        strat_a = truth_eval["by_stratum"]["Stratum_A_Core_Biomass"]
        self.assertEqual(strat_a["n_truth_proteins"], 1)
        self.assertEqual(strat_a["ec4_exact_matches"], 1)
        self.assertEqual(strat_a["ec3_matches"], 1)

        strat_b = truth_eval["by_stratum"]["Stratum_B_Catabolic_Assay_Specific"]
        self.assertEqual(strat_b["n_truth_proteins"], 1)
        self.assertEqual(strat_b["ec4_exact_matches"], 0)
        self.assertEqual(strat_b["ec3_matches"], 1)

    def test_collection_extraction_and_report_prose_has_no_numbers(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            coll_dir = root / "collection"
            att_dir = coll_dir / "attempts" / "0" / "job-100" / "production"
            (att_dir / "stage0").mkdir(parents=True)
            (att_dir / "stage1").mkdir(parents=True)
            (coll_dir / "attempts" / "0" / "job-100" / "completed.json").write_text(
                '{"status": "completed"}\n', encoding="utf-8"
            )

            plan = {
                "cells": [
                    {"index": 0, "predictor": "topec", "assembly": "GCF_000236665.1"},
                    {"index": 1, "predictor": "topec", "assembly": "GCF_000015245.1"},
                ]
            }
            (coll_dir / "plan.json").write_text(json.dumps(plan), encoding="utf-8")

            tsv_bytes = _make_tsv_bytes(
                [
                    _ec_row("WP_014235449.1", "EC:4.2.1.113", 0.72),
                    _ec_row("WP_014235128.1", "EC:2.7.7.6", 0.81),
                    _ec_row("WP_014235129.1", "EC:6.3.5.7", 0.64),
                    _ec_row("WP_014235130.1", "EC:1.13.12.7", 0.55),
                ]
            )
            for arm_name, priors_flag in (("stage0", "--no-priors"), ("stage1", "")):
                arm_dir = att_dir / arm_name
                (arm_dir / "GCF_000236665.1.integrated.tsv").write_bytes(tsv_bytes)
                (arm_dir / "batch_stats.tsv").write_text(
                    "genome_id\tenforce_stage2\tenforce_mode\tstatus\n"
                    "GCF_000236665.1\tfalse\trelax\tcompleted\n",
                    encoding="utf-8",
                )
                (arm_dir / "manifest.tsv").write_text(
                    "genome_id\ttaxon_id\tgaps_path\nGCF_000236665.1\t640081\t\n",
                    encoding="utf-8",
                )
                (arm_dir / "GCF_000236665.1.integrated.tsv.stages.json").write_text(
                    json.dumps(
                        {
                            "stages": [
                                {"name": "symbolic_enforcement", "status": "skipped"}
                            ]
                        }
                    ),
                    encoding="utf-8",
                )
                cmd = ["java", "-jar", "gspa.jar"]
                if priors_flag:
                    cmd.append(priors_flag)
                (arm_dir / "command.json").write_text(
                    json.dumps({"command": cmd}), encoding="utf-8"
                )

            out_dir = root / "out"
            receipts_idx = extract_collection(
                collection_dir=coll_dir, out_dir=out_dir, workers=1
            )
            self.assertEqual(receipts_idx["extracted_cells"], 1)
            self.assertEqual(receipts_idx["pending_cells"], 1)

            summary_file = out_dir / "cells" / "000_topec_GCF_000236665.1.json"
            self.assertTrue(summary_file.exists())
            cell_summary = json.loads(summary_file.read_text(encoding="utf-8"))

            evaluation = evaluate_collection_summaries([cell_summary])
            cell_eval = evaluation["cells"][0]
            tc_stage1 = cell_eval["arms"]["stage1"]["taxon_consistency"]
            self.assertFalse(tc_stage1["enforcement_on"])
            self.assertFalse(tc_stage1["enforce_stage2"])
            self.assertEqual(tc_stage1["symbolic_enforcement_status"], "skipped")
            self.assertEqual(tc_stage1["violation_rows_called"], 1)

            md = render_markdown_report(evaluation)
            self.assertIn("| Genome | Organism | Exempt | Predictor | Arm |", md)
            self.assertIn("| GCF_000236665.1 |", md)

            # Verify no prose line contains any digit (0-9): all numbers live strictly in tables
            for line in md.splitlines():
                stripped = line.strip()
                if not stripped or stripped.startswith("|") or stripped.startswith("#"):
                    continue
                self.assertIsNone(
                    re.search(r"[0-9]", stripped),
                    f"Prose line unexpectedly contains a number: {stripped!r}",
                )

    def test_sbatch_dry_run_prints_commands_without_submitting(self) -> None:
        sbatch_script = V3_DIR / "shaheen" / "biology-extract.sbatch"
        # Inside a Slurm job the script would switch to worker mode; test the launcher with Slurm variables removed.
        env = {k: v for k, v in os.environ.items() if not k.startswith("SLURM_")}
        proc = subprocess.run(
            ["bash", str(sbatch_script)],
            capture_output=True,
            text=True,
            check=True,
            env=env,
        )
        out = proc.stdout
        self.assertIn("--account=k10018", out)
        self.assertIn("--partition=workq", out)
        self.assertIn("nid00022", out)
        self.assertIn("nid00497", out)
        self.assertIn("nid01033", out)
        self.assertIn("/scratch/project/k10018/gspa-v3/out/biology-rc1", out)

    def test_real_archived_price_cell_end_to_end(self) -> None:
        archived_cell = V3_DIR / "results" / "price-ec-compatibility-v1" / "after"
        s0 = archived_cell / "stage0" / "GCF_000236665.1.integrated.tsv"
        s1 = archived_cell / "stage1" / "GCF_000236665.1.integrated.tsv"
        if not (s0.exists() and s1.exists()):
            self.skipTest("Archived compatibility cell not present")

        summary, receipt = extract_cell(
            stage0_tsv=s0,
            stage1_tsv=s1,
            predictor="topec",
            assembly="GCF_000236665.1",
            cell_index=0,
            taxon_constraints=self.taxon_constraints,
            ec_entries=self.ec_entries,
            price_targets_by_assembly=self.price_targets,
        )
        self.assertTrue(summary["stage0_equals_stage1_sha256"])
        self.assertEqual(
            receipt["inputs"]["stage0_integrated_tsv"]["rows"],
            3481,
        )
        ev = evaluate_collection_summaries(
            [summary],
            taxon_constraints=self.taxon_constraints,
            ec_entries=self.ec_entries,
            price_targets_by_assembly=self.price_targets,
        )
        cell_ev = ev["cells"][0]
        self.assertFalse(cell_ev["exemption"]["exempt"])
        self.assertEqual(
            cell_ev["arms"]["stage0"]["translation"]["called"]["satisfied_amino_acids"],
            9,
        )
        self.assertEqual(
            cell_ev["arms"]["stage0"]["translation"]["called"]["indirect_gat_count"],
            1,
        )
        self.assertEqual(
            cell_ev["arms"]["stage0"]["core_metabolism"]["called"]["present_roles"],
            8,
        )
        self.assertEqual(
            cell_ev["arms"]["stage0"]["core_metabolism"]["retained"]["present_roles"],
            9,
        )
        self.assertEqual(
            cell_ev["arms"]["stage1"]["taxon_consistency"]["violation_rows_called"],
            14,
        )
        self.assertFalse(
            cell_ev["arms"]["stage1"]["taxon_consistency"]["enforcement_on"]
        )


if __name__ == "__main__":
    unittest.main()
